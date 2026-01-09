package seqmer

import (
	"bufio"
	"encoding/csv"
	"encoding/json"
	"fmt"
	"globals"
	"image/color"
	"io/fs"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"time"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
	"gonum.org/v1/plot/vg/draw"
)

//import bitsy "github.com/yourbasic/bit"

//
// ## sequence-related structures and functions
//

// these nucs are currently only used in RC
var nucs = []byte{'G', 'A', 'C', 'T'}
var cnucs = []byte{'C', 'T', 'G', 'A'}

// note, sequence names may not be unique within a given file
// might be useful to have supersets
type SeqSets struct {
	setType   *string    // file name e.g a sequence, a contig, a genome, a sample, a geographic location
	name      string     // name (as opposed to path) may not be unique
	path      string     // path to get to this node, although redundant to sets, it is a unique name though
	subSets   []*SeqSets // name of the sequence following break character
	superSets []*SeqSets
	// for each SubSets in read order, pointer to (kmer, [positions]) map
	filetype *string
	kmers    *Oligos //[]*map[string][]int // kmers just a pointer to Oligos which stores kmers and their positions
	// record total number of strings in the file
	numSubSets int

	klen       int
	linelimit  int
	dofilter   bool
	entrystart string
}

// Init for SeqSets creates new root structure of the sequence sets
func (seqs *SeqSets) Init(filetype string, dofilter bool, llimit int, klen int) {
	// seqs.setType = seqfile  // there is not necessarily a sequence file associated with a set
	seqs.filetype = &filetype // for now, this is copied to each new SeqSet in a hierarchy
	seqs.dofilter = dofilter

	seqs.linelimit = llimit
	seqs.klen = klen

	seqs.subSets = make([]*SeqSets, 0)
	// seqs.superSets = make([]*SeqSets, 0) // to begin, we will define the root as nil supersets
}

// AddTo for SeqSets adds
func (seqs *SeqSets) AddTo(filetype *string, dofilter bool, llimit int, klen int) {
	// seqs.setType = seqfile  // there is not necessarily a sequence file associated with a set
	seqs.filetype = filetype // for now, this is copied to each new SeqSet in a hierarchy
	seqs.dofilter = dofilter

	seqs.linelimit = llimit
	seqs.klen = klen

	seqs.subSets = make([]*SeqSets, 0)
	// seqs.superSets = make([]*SeqSets, 0) // to begin, we will define the root as nil supersets
}

// Sequences holds sequence info
// this may be depracated in lieu of SeqSets
type Sequences struct {
	Name      string
	seqfile   string
	seqmap    map[string]string
	seqfilter map[string]*filter
	Filterdir string
	minlength int
	record    bool
	dofilter  bool
	linelimit int
	linemin   int
	count     int
	isalign   bool
	firstlen  int
	Totallen  int
	// could add length map
	outfile    string
	filetype   string
	entrystart string
	// fancy reading and printing selection
	offset      int
	interval    int
	keepcompany int
}

// Init creates new parameter structure of hash types
func (seqs *Sequences) Init(seqfile string, minlen int, llimit int, minline int, dorecord bool, filetype string, dofilter bool) {
	seqs.Name = seqfile
	seqs.seqfile = seqfile
	seqs.seqmap = make(map[string]string)
	seqs.seqfilter = make(map[string]*filter)
	seqs.Filterdir = "none"

	seqs.minlength = minlen
	seqs.linelimit = llimit
	seqs.linemin = minline
	seqs.record = dorecord
	seqs.count = 0
	seqs.isalign = false
	seqs.dofilter = dofilter
	// fmt.Println("input and recorded filter states", dofilter, seqs.dofilter)
	seqs.firstlen = 0
	seqs.Totallen = 0
	seqs.outfile = "sequences.txt"
	seqs.filetype = filetype
	if filetype == "fastq" {
		seqs.entrystart = "@"
	} else {
		seqs.entrystart = ">"
	}
}

// Newseqfile updates new sequence file name
func (seqs *Sequences) Newseqfile(seqfile string) {
	seqs.seqfile = seqfile
}

// InitFancy creates selective sequence reading
func (seqs *Sequences) InitFancy(offset int, interval int, keepcompany int) {
	seqs.offset = offset
	seqs.interval = interval
	seqs.keepcompany = keepcompany
}

// Print prints out sequences
func (seqs *Sequences) Print(printmode string, headers bool) {
	fmt.Println("Printing sequences to ", seqs.outfile)
	fkout, _ := os.Create(seqs.outfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output

	fullmode := "long"
	if headers {
		fmt.Fprintf(kwriter, "%s\t%s\t", "Name", "Length")
		if printmode == fullmode {
			fmt.Fprintf(kwriter, "%s\t", "sequence")
		}
		fmt.Fprintf(kwriter, "\n")
	}
	for name, seq := range seqs.seqmap {
		seqlen := len(seq)
		fmt.Fprintf(kwriter, "%s\t%d", name, seqlen)
		if printmode == fullmode {
			fmt.Fprintf(kwriter, "%s", seq)
		}
		fmt.Fprintf(kwriter, "\n")
	}
}

// Print prints out sequences, parsing genbank name
func (seqs *Sequences) Printparse(printmode string, headers bool) {
	fmt.Println("Printing sequences to ", seqs.outfile)
	fkout, _ := os.Create(seqs.outfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output

	parsename := "shortparsename"
	const firstsplitter = "|"
	const splitter = "="
	if printmode == parsename {
		if headers {
			for name, _ := range seqs.seqmap {
				fmt.Fprintf(kwriter, "%s", name)
				fmt.Fprintf(kwriter, "\n")
			}
			fmt.Fprintf(kwriter, "\n")
		}
		count := 0
		for name, _ := range seqs.seqmap {
			fields := strings.Fields(name)
			count++
			fmt.Fprintf(kwriter, "%d\t", count)
			var namemap map[string]string
			namemap = make(map[string]string)
			lastkey := "empty"
			namemap[lastkey] = lastkey
			// for i, j := 0, 1; i < 10; i, j = i+1, j+1
			bitcount := 0
			for _, bit := range fields {
				bitcount++
				bit := strings.Trim(bit, "[]")
				bitpair := strings.SplitN(bit, splitter, 2)
				if bitcount == 1 {
					bitpair = strings.SplitN(bit, firstsplitter, 2)
				}
				if len(bitpair) > 1 {
					namemap[bitpair[0]] = bitpair[1]
					lastkey = bitpair[0]
				} else {
					bitcount--
					namemap[lastkey] = namemap[lastkey] + " " + bit
				}
			}
			for key, value := range namemap {
				if key != "empty" {
					fmt.Fprintf(kwriter, " %s = %s ", key, value)
				}
			}
			fmt.Fprintf(kwriter, "\t%d\n", len(namemap))
		}
	}
}

// filter holds info about read filters
type filter struct {
	name      string
	readlen   int
	pstart    int
	pend      int
	Wuleft    int
	Wuright   int
	fp        int
	lp        int
	direction string

	olist []string // stable ordered list pointing to kmer info
	klist []*oligo // stable ordered list pointing to kmer info
} //

// reads holds info about the reads (but not the sequences)
type reads struct {
	name   string
	length int
	start  int
	end    int
} //

// recordseq adds seq, name, and count total information to seqs
// this is not used much because files too big
func (seqs *Sequences) recordseq(seq string, name string) {
	seqlen := len(seq)
	if (seqlen > seqs.minlength) && seqs.record {
		blurb := "Error 73, sequence name already exists => "
		if _, exists := seqs.seqmap[name]; exists {
			panic(blurb + name)
		}
		// a fancier version would add a number to the name or something instead of panicking
		seqs.seqmap[name] = seq
		seqs.count += 1
		seqs.Totallen += len(seq)
		fmt.Println("Recorded sequence ", name, "length", seqlen)
	}
}

// Printseq prints sequence information in controllable fashion
// not used much because files so big
func (seqs *Sequences) Printseq(printmode string) {
	var name, seq string
	var count, lcount int

	// reading stuff
	fmt.Println("File to open is ", seqs.seqfile)
	fpin, err := os.Open(seqs.seqfile)
	globals.Check(err)
	defer fpin.Close()
	scanner := bufio.NewScanner(fpin)

	// read, record, count kmers
	noisy := true
	for scanner.Scan() {
		lcount += 1
		line := scanner.Text()              // should not include eol
		trimline := strings.TrimSpace(line) // trim off leading and lagging whitespace
		if strings.HasPrefix(line, ">") {
			fmt.Println("header line", line)
			seqs.recordseq(seq, name)
			seq = ""
			name = strings.TrimPrefix(trimline, ">")
			count += 1
			if noisy {
				fmt.Println("New seq", name, "number", count)
			}
		} else {
			seq = seq + trimline
			keepcompany(lcount, len(seq), 10000, 2000, 100000)
		}
	}
	seqs.recordseq(seq, name)
	fmt.Println("Lines counted1\n", count, lcount)
	fmt.Println("Lines counted\n", count, lcount)
}

//
// //  oligo-related structures and functions // //
//

// oligo holds info for individual kmer
type oligo struct {
	name    string
	revcomp string
	poses   []int
	kcount  int
} // will make global seqs

// Init creates new parameter structure of hash types
func (k *oligo) Init(kmer string, kcount int) {
	k.name = kmer
	k.revcomp = rc(kmer)
	k.kcount = kcount
	k.poses = make([]int, 0)
}

func (k *oligo) Getname() string {
	return k.name
}

func (k *oligo) Getkcount() int {
	return k.kcount
}

// oligos holds kmer map of all kmers, tracks length and total kmers
type Oligos struct {
	name       string
	kmap       map[string]*oligo
	rmap       map[string]*oligo
	kcount     map[string]int
	locprimer  map[int]int
	readcounts map[int]int
	frags      map[string]*fragment
	alligos    map[string]*kneighbors // might convert this to a structure to include klen
	offmin     int
	total      int
	klen       int
	minprint   int
	Outfile    string
	POutfile   string
	kfile      string
	printNs    bool
	remnant    string
} // will make global kmers

// Fragment holds individual fragment info so that positions are assigned properly
// when there are multiple fragments
// continued use of kmap using kmer string rather than kper *Kmer is suboptimal
type fragment struct {
	name   string
	kcount map[*kneighbors]int
	total  int
	kmap   map[*kneighbors]*oligo
	klocs  []*kneighbors
}

// addfragment checks for existence of named sequence fragment
// return if exists, or initialize if not
// continued use of kmer string rather than kper *Kmer is suboptimal
func (kmers *Oligos) addfragment(fragname string) *fragment {
	if kmers.frags[fragname] == nil {
		kmers.frags[fragname] = new(fragment)
		frag := kmers.frags[fragname]
		frag.name = fragname
		frag.kcount = make(map[*kneighbors]int)
		frag.kmap = make(map[*kneighbors]*oligo)
		frag.klocs = make([]*kneighbors, 0)
		return frag // return newly created fragment
	}
	return kmers.frags[fragname] // if fragment already existed, return that
}

// addpose creates the storm and the PCloud list adds to tempest and returns the storm
func (f *fragment) addpose(kptr *kneighbors, location int) {
	o := f.kmap[kptr]
	if o == nil {
		o = new(oligo)
		o.poses = make([]int, 0)
		f.kmap[kptr] = o
	}
	if kptr != nil {
		o.poses = append(o.poses, location)
	}
	f.kcount[kptr]++                // this might be bad form if kptr is nil, but for now...
	f.klocs = append(f.klocs, kptr) // this might be bad form if kptr is nil, but for now...
	f.total++
}

// Init creates new parameter structure of hash types
func (kmers *Oligos) Init(klen int, koutfile string, poutfile string, doNs bool, kminprint int, dofind bool, alligos *alloligos) {
	// fmt.Println("Running Init kmers (Oligos) ", koutfile)
	kmers.name = "generic_seqfile"
	kmers.kmap = make(map[string]*oligo)
	kmers.rmap = make(map[string]*oligo)
	kmers.kcount = make(map[string]int)
	kmers.locprimer = make(map[int]int)
	kmers.readcounts = make(map[int]int)
	kmers.frags = make(map[string]*fragment)
	kmers.alligos = alligos.alligos

	kmers.total = 0
	kmers.minprint = kminprint
	kmers.klen = klen
	kmers.Outfile = koutfile
	kmers.POutfile = poutfile
	kmers.printNs = doNs
}

// kneighbors holds the neighborhood of a kmer, held in alligos string map
type kneighbors struct {
	// KmerString string // KmerString deprecated by use of alligos hash that is created once with all the kmer strings
	neighbors []*kneighbors // neighbors must be initialized by FindAlligos()
}

// alloligos holds a kmer string map to set of kneighbors
type alloligos struct {
	alligos map[string]*kneighbors
}

// reverseoligos holds a *kneighbors kptr map of kmer strings for reverse lookup
type reverseoligos struct {
	revoligos map[*kneighbors]string
	klen      int
}

// Newalligos creates a new alligos structure
func Newalligos(klen int, dofind bool, nofind bool) *alloligos {
	alloligos := new(alloligos)
	alloligos.alligos = make(map[string]*kneighbors)
	if nofind {
		fmt.Println("Newalligos returning empty map because nofind")
	} else {
		alloligos.findAlligos(klen, dofind)
	}
	return alloligos
}

// Newalligos creates a new alligos structure
func NewRevAlligos(klen int) *reverseoligos {
	ro := new(reverseoligos)
	ro.revoligos = make(map[*kneighbors]string)
	ro.klen = klen
	return ro
}

func (r *reverseoligos) Reversefill(alloligos *alloligos) {
	alligos := alloligos.alligos
	for kmer, kptr := range alligos {
		if r.revoligos[kptr] == "" {
			r.revoligos[kptr] = kmer
		}
	}
}

// findAlligos finds all oligos and puts them into list
func (a *alloligos) findAlligos(klen int, dofind bool) {
	now := time.Now()
	fmt.Println("starting findAlligos", now)
	alligos := a.alligos
	var kmer string
	for i := 0; i < klen; i++ {
		kmer += string(nucs[0])
	}
	findneighbors(alligos, kmer, klen, 0, dofind)
	now = time.Now()
	fmt.Println("ending findAlligos", now)
}

// findneighbors recursively finds all neighbors that don't already exist, then adds them to alligos
// returns *Kmer for given kmer string
// dofind is intended to allow choice on finding neighbors or not
// if not, should return all oligos pointers but without list
func findneighbors(alligos map[string]*kneighbors, kmer string, klen int, loc int, dofind bool) *kneighbors {
	var mutant string
	if loc < klen {
		if alligos[kmer] == nil {
			alligos[kmer] = new(kneighbors)
		}
		kbits := []byte(kmer)
		mutbits := []byte(kmer)
		for n := range nucs {
			if nucs[n] != kbits[loc] {
				mutbits[loc] = nucs[n]
				mutant = string(mutbits)
				if alligos[mutant] == nil {
					findneighbors(alligos, strings.Clone(mutant), klen, loc, dofind) // send clone so it has new string
				}
				if dofind {
					alligos[kmer].neighbors = append(alligos[kmer].neighbors, alligos[mutant])
				}
			}
		}
		// fmt.Println("find neighbors loc => +1 , kmer ", loc, kmer)
		findneighbors(alligos, kmer, klen, loc+1, dofind) // go to next position
	} // end if location not bigger than klen
	return alligos[kmer]
}

// clone := strings.Clone(s)

// Storm contains information on the set of Pclouds that make up the storm
// formerly kmerstorm
type storm struct {
	parentseq  string      // parent sequence fragment name
	parentfrag string      // parent sequence fragment name
	seedmer    *kneighbors // focal or seed kmer, at the center of the storm
	seedrefpos int         // position of seed in starting reference sequence
	seedpos    int         // seed position in storm
	endrefpos  int         // end position in storm
	size       int         // number of clouds in storm
	PClouds    []*PCloud   // list of PClouds
}

// Print2 storm
func (s *storm) Print2(revigos *reverseoligos, stormid int) {
	fmt.Println("\t\tstorm", stormid, revigos.revoligos[s.seedmer], s.seedrefpos, s.endrefpos, s.size, len(s.PClouds), s.parentseq, s.parentfrag)
	for pid, _ := range s.PClouds {
		fmt.Printf("%d\t", pid)
	}
	fmt.Println()
	for _, pcloud := range s.PClouds {
		fmt.Printf("%s\t", revigos.revoligos[pcloud.seedmer])
	}
	fmt.Println()
}

// Print storm
func (s *storm) Print(revigos *reverseoligos, id int) {
	if id < 10 || id > 3800 {
		fmt.Println("\t\tstorm", id, revigos.revoligos[s.seedmer], s.seedrefpos, s.endrefpos, s.size, len(s.PClouds))
	}
}

// Pcloud holds a cloud of kmers that are probably homologous, formerly kmercloud
type PCloud struct {
	seedmer   *kneighbors   // original kmer in cloud
	count     int           // number of kmers in cloud
	loc       int           // location of pcloud in storm
	cloudmers []*kneighbors // list of kmers in the cloud
	// should this be a list or a hash?
}

// newcloud creates the cloud and the cloudmers list adds to storm and returns the cloud
func (s *storm) newcloud(kptr *kneighbors, cloudloc int) *PCloud {
	cloud := new(PCloud)
	cloud.cloudmers = make([]*kneighbors, 0)
	cloud.loc = cloudloc
	// fmt.Println("adding cloud to storm", cloudloc, kptr)
	cloud.seedmer = kptr
	cloud.cloudmers = append(cloud.cloudmers, cloud.seedmer) //	redundant, but sometimes handy to have seed here too
	cloud.count = 1
	s.PClouds = append(s.PClouds, cloud)
	return cloud
}

// newstormfinish closes out the initial storm creation
// fails if cloudcount is 0
func (s *storm) newstormfinish(cloudcount int, seedloc int, pos int, seedrefpos int, parentseq string, parentfrag string) {
	s.size = cloudcount + 1
	if seedloc >= cloudcount {
		fmt.Println("fixing storm finish length")
		seedloc = cloudcount - 1
		seedrefpos = pos
		if seedloc < 0 { // a hack, but seems to work
			seedloc = 0
		}
	}
	s.seedmer = s.PClouds[seedloc].seedmer // storm seedmer same as cloud seedmer at seed location
	s.seedpos = seedloc
	s.seedrefpos = seedrefpos
	s.endrefpos = pos
	s.parentseq = parentseq
	s.parentfrag = parentfrag

	// if pos < 1000 {
	// 	fmt.Println("seedstorms finish ccount pos seedloc", cloudcount, seedloc, pos)
	// }

}

// Tempest contains a map of kmer storms, meant to represent complete set of storms for homologous sequences
type Tempest struct {
	SeedName     string
	count        int // umber of storms in the tempest
	stormseedloc int
	stormsize    int
	kstorms      []*storm
}

// Print prints ths storm
func (t *Tempest) Print(revigos *reverseoligos, loudness int) {
	fmt.Println("printing storms at loudness", loudness)
	fmt.Println("\tstorm size ", t.stormsize)
	fmt.Println("\tstorm seed location ", t.stormseedloc)
	fmt.Println("\tkstorms built ", t.count)
	if loudness > 1 {
		for i, storm := range t.kstorms {
			storm.Print(revigos, i)
		}
	}
}

// Print prints ths storm
func (t *Tempest) Print2(revigos *reverseoligos) {
	fmt.Println("\tstorm size ", t.stormsize)
	fmt.Println("\tstorm seed location ", t.stormseedloc)
	fmt.Println("\tkstorms built ", t.count)
	for i, storm := range t.kstorms {
		fmt.Println("i", i)
		storm.Print2(revigos, i)
	}
}

// newstorm creates the storm and the PCloud list adds to tempest and returns the storm
func (t *Tempest) newstorm() *storm {
	s := new(storm)
	s.PClouds = make([]*PCloud, 0)
	t.kstorms = append(t.kstorms, s)
	t.count++
	return s
}

// Init creates the PCloud list for structure storm
func (t *Tempest) Init(ssize int) {
	t.kstorms = make([]*storm, 0)
	t.stormsize = ssize
	t.stormseedloc = int(ssize / 2)
}

func (t *Tempest) GetNumberOfStorms() int {
	return t.count
}

// stormhit contains info on a hit to a storm
type stormhit struct {
	matches    []int // list of match positions
	matchcount int
	frag       *fragment
	storm      *storm
	boundary   int
	location   int // first match position locates the hit
	// location is redundant to first int in matches list
}

// addpos
func (sh *stormhit) addpos(pos int) {
	sh.matches = append(sh.matches, pos)
	// sh.cloudlocs = append(sh.cloudlocs, cloud.location)
	sh.matchcount += 1
}

// Print prints out the identified storm homologs
func (sh *stormhit) Print(kwriter *bufio.Writer, revigos *reverseoligos, stormID int) {
	if stormID < 500 && sh.matchcount < 19 {
		fmt.Fprintf(kwriter, "%d\t%d\n", stormID, sh.matchcount)
	}
}

// if stormID < 500 || stormID > 5000 && sh.matchcount < 19 {
// fmt.Println("hit at storm", stormID, "details", sh.matchcount)
// for _, match := range sh.matches {
// 	fmt.Printf("\t%d", match)
// 	fmt.Println()
// }

// stormologs is a list of potentially homologous hits to a storm
type stormologs struct {
	goodhits       []*stormhit // map hits to given storm ID
	goodhitmin     int         // minimum used to qualify as good hit
	hitcount       int
	uniquehitcount int
	matchtotal     int
	frag           string
	Seqname        string
}

// newstorm creates the storm and the PCloud list adds to tempest and returns the storm
func newstormologs(goodmin int, name string) *stormologs {
	hitlst := new(stormologs)
	hitlst.goodhits = make([]*stormhit, 0) // should be a call to init, but so it later
	hitlst.goodhitmin = goodmin
	hitlst.Seqname = name
	return hitlst
}

// addnewhit adds a new hit to stormologs and increments counter
func (slogs *stormologs) addnewhit(pos int) {
	newhit := new(stormhit)
	newhit.addpos(pos)
	newhit.location = pos //	want first time to record stormhit.location for convenience
	slogs.goodhits = append(slogs.goodhits, newhit)
	slogs.hitcount++
}

// Print prints out the identified storm homologs
func (slogs *stormologs) Print(slogs_outfile string, revigos *reverseoligos, goodmin int) {
	fmt.Println("Printing stormologs to ", slogs_outfile)
	fkout, _ := os.Create(slogs_outfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output

	// fmt.Println("\n\t\tPrinting stormologs\n", goodmin)
	//	fmt.Fprintf(kwriter, "%s\t%s\t", "Name", "Length")
	// fmt.Printf("\n\t\tPrinting stormologs %d\n", goodmin)
	fmt.Fprintf(kwriter, "refloc\ttargetloc\tmatchcount\n")
	for _, goodhit := range slogs.goodhits {
		refpos := goodhit.storm.seedrefpos
		fmt.Fprintf(kwriter, "%d\t", refpos)
		fmt.Fprintf(kwriter, "%d\t", goodhit.location)
		fmt.Fprintf(kwriter, "%d\n", goodhit.matchcount)
	}
	fmt.Fprintf(kwriter, "\n")
}

// func rollingStats(data plotter.XYs, w int) (mean, stddev []float64) {
// 	mean = make([]float64, len(data))
// 	stddev = make([]float64, len(data))

// 	var sum, sumSq float64
// 	queue := make([]float64, 0, w)

// 	for i, pt := range data {
// 		// push new sample
// 		queue = append(queue, pt.Y)
// 		sum += pt.Y
// 		sumSq += pt.Y * pt.Y

// 		// pop oldest sample if window is full
// 		if len(queue) > w {
// 			old := queue[0]
// 			queue = queue[1:]
// 			sum -= old
// 			sumSq -= old * old
// 		}

// 		n := float64(len(queue))
// 		m := sum / n
// 		mean[i] = m
// 		variance := sumSq/n - m*m
// 		if variance < 0 {
// 			variance = 0 // numeric safety
// 		}
// 		stddev[i] = math.Sqrt(variance)
// 	}
// 	return
// }

// type chunk struct {
// 	xStart, xEnd float64
// 	mean, sigma  float64
// }

// func makeChunks(data plotter.XYs, chunkSize int) []chunk {
// 	if chunkSize <= 0 {
// 		panic("chunkSize must be > 0")
// 	}
// 	var chunks []chunk

// 	for i := 0; i < len(data); i += chunkSize {
// 		end := i + chunkSize
// 		if end > len(data) {
// 			end = len(data)
// 		}
// 		c := chunk{xStart: data[i].X, xEnd: data[end-1].X}

// 		// mean & variance
// 		var sum, sumSq float64
// 		for _, pt := range data[i:end] {
// 			sum += pt.Y
// 			sumSq += pt.Y * pt.Y
// 		}
// 		n := float64(end - i)
// 		c.mean = sum / n
// 		vari := sumSq/n - c.mean*c.mean
// 		if vari < 0 {
// 			vari = 0 // numeric safety
// 		}
// 		c.sigma = math.Sqrt(vari)
// 		chunks = append(chunks, c)
// 	}
// 	return chunks
// }

// func realignBlocks(pts plotter.XYs, tol float64) (adjusted plotter.XYs) {
// 	if len(pts) == 0 {
// 		return nil
// 	}

// 	// sort by X (already sorted in demo data; sort if yours is not)

// 	blockΔ := pts[0].Y - pts[0].X // current offset
// 	for _, p := range pts {
// 		δ := p.Y - p.X
// 		if math.Abs(δ-blockΔ) > tol {
// 			// new block starts here
// 			blockΔ = δ
// 			// fmt.Println(δ)
// 		}
// 		adjusted = append(adjusted, plotter.XY{X: p.X, Y: p.Y - blockΔ})
// 	}
// 	return
// }

// type contig struct {
// 	x0, x1 float64 // span on reference (x-axis)
// 	delta  float64 // offset
// 	fwd    bool    // orientation
// }

// // find contiguous blocks on the Y-axis that keep the same orientation
// func splitIntoContigs(pts plotter.XYs, maxGap float64) []contig {
// 	if len(pts) == 0 {
// 		return nil
// 	}
// 	// sort by Y (keeps concatenated contigs together)
// 	sort.Slice(pts, func(i, j int) bool { return pts[i].Y < pts[j].Y })

// 	var cs []contig
// 	start := 0
// 	for i := 1; i < len(pts); i++ {
// 		dx := pts[i].X - pts[i-1].X
// 		dy := pts[i].Y - pts[i-1].Y
// 		stepFwd := dx*dy >= 0
// 		currFwd := pts[i-1].Y-pts[start].Y >= 0
// 		if stepFwd != currFwd || math.Abs(dy/dx) > 1.1 { //math.Abs(dy) > maxGap {
// 			cs = append(cs, summariseBlock(pts[start:i], currFwd))
// 			start = i
// 		}
// 	}
// 	cs = append(cs, summariseBlock(pts[start:], pts[len(pts)-1].Y-pts[start].Y >= 0))
// 	fmt.Println("cs")
// 	fmt.Println(len(cs))
// 	return cs
// }

// // Assign each point its contig index
// func pointContigIndex(pts plotter.XYs, cs []contig) []int {
// 	idx := make([]int, len(pts))
// 	for i, p := range pts {
// 		for k, c := range cs {
// 			if p.X >= c.x0 && p.X <= c.x1 {
// 				idx[i] = k
// 				break
// 			}
// 		}
// 	}
// 	return idx
// }

// // reduce one block to its orientation and delta
// func summariseBlock(b plotter.XYs, fwd bool) contig {
// 	var sum float64
// 	for _, p := range b {
// 		if fwd {
// 			sum += p.Y - p.X
// 		} else {
// 			sum += p.Y + p.X
// 		}
// 	}
// 	delta := sum / float64(len(b))
// 	return contig{x0: b[0].X, x1: b[len(b)-1].X, delta: delta, fwd: fwd}
// }

// // apply the correction
// func realign(pts plotter.XYs, cs []contig) plotter.XYs {
// 	out := make(plotter.XYs, 0, len(pts))
// 	for _, p := range pts {
// 		for _, c := range cs {
// 			if p.X >= c.x0 && p.X <= c.x1 {
// 				if c.fwd {
// 					// fmt.Println("fwd")
// 					out = append(out, plotter.XY{X: p.X, Y: p.Y - c.delta})
// 				} else {
// 					// fmt.Println("rev")
// 					out = append(out, plotter.XY{X: p.X, Y: -p.Y + c.delta})
// 				}
// 				break
// 			}
// 		}
// 	}
// 	fmt.Println(len(pts))
// 	return out
// }

// func distinctColours(N int) []color.Color {
// 	colours := make([]color.Color, N)
// 	for i := 0; i < N; i++ {
// 		h := float64(i) / float64(N)       // 0 … <1
// 		r, g, b := HSVtoRGB(h, 0.65, 0.95) // sat=0.65, val=0.95
// 		colours[i] = color.NRGBA{
// 			R: uint8(r * 255),
// 			G: uint8(g * 255),
// 			B: uint8(b * 255),
// 			A: 255,
// 		}
// 	}
// 	return colours
// }

// func HSVtoRGB(h, s, v float64) (r, g, b float64) {
// 	if s == 0 {
// 		return v, v, v
// 	}
// 	h = math.Mod(h*6, 6)
// 	i := math.Floor(h)
// 	f := h - i
// 	p := v * (1 - s)
// 	q := v * (1 - s*f)
// 	t := v * (1 - s*(1-f))
// 	switch int(i) {
// 	case 0:
// 		r, g, b = v, t, p
// 	case 1:
// 		r, g, b = q, v, p
// 	case 2:
// 		r, g, b = p, v, t
// 	case 3:
// 		r, g, b = p, q, v
// 	case 4:
// 		r, g, b = t, p, v
// 	default: // case 5:
// 		r, g, b = v, p, q
// 	}
// 	return
// }

// PlotStorms plots the identified storm homologs
func (slogs *stormologs) PlotStorms(slogs_outfile string, revigos *reverseoligos, goodmin int, scattersize float64, plotxdim int, plotydim int, fontsize int, titleFontSize int, averageWindow int, lineWidth int, spacing float64) (float64, float64, float64, int, plotter.XYs, plotter.XYs) {

	if len(slogs.goodhits) == 0 {

		return 0.0, 0.0, 0.0, 0, plotter.XYs{}, plotter.XYs{}
	}
	fmt.Println("Plotting stormologs to ", slogs_outfile)

	titleSize := vg.Points(float64(titleFontSize))
	integerLabels := func(min, max float64) []plot.Tick {
		// Start from the default ticks so we keep the same positions.
		ticks := plot.DefaultTicks{}.Ticks(min, max)

		for i := range ticks {
			ticks[i].Label = fmt.Sprintf("%.0f", ticks[i].Value) // no decimals
		}
		return ticks
	}
	//
	// process data plots
	//
	// scatter plot, values need to be float64 type
	points_positions := make(plotter.XYs, len(slogs.goodhits))
	points_score := make(plotter.XYs, len(slogs.goodhits))
	stormRefPositionMap := make(map[int]int)

	for i, goodhit := range slogs.goodhits {
		refpos := goodhit.storm.seedrefpos
		stormRefPositionMap[refpos]++
		points_positions[i].X = float64(refpos)
		points_positions[i].Y = float64(goodhit.location)
		points_score[i].X = float64(refpos)
		points_score[i].Y = float64(goodhit.matchcount)
	}

	// plot the BED diagram

	p_bed := plot.New()

	var refPositions []int
	for k := range stormRefPositionMap {
		refPositions = append(refPositions, k)
	}
	sort.Ints(refPositions)

	storm_counts_points := make(plotter.XYs, len(refPositions))
	for i, k := range refPositions {
		storm_counts_points[i].X = float64(k)
		storm_counts_points[i].Y = float64(stormRefPositionMap[k])
	}

	p_bed.Title.Text = "Storm Counts"
	p_bed.X.Label.Text = "Reference" + " Sequence Position"
	p_bed.Y.Label.Text = "Number Found in Target Sequence"
	p_bed.Title.TextStyle.Font.Size = titleSize
	p_bed.X.Label.TextStyle.Font.Size = titleSize
	p_bed.Y.Label.TextStyle.Font.Size = titleSize

	p_bed.X.Tick.Marker = plot.TickerFunc(integerLabels)
	p_bed.Y.Tick.Marker = plot.TickerFunc(integerLabels)
	p_bed.X.Tick.Label.Font.Size = vg.Points(float64(fontsize))
	p_bed.Y.Tick.Label.Font.Size = vg.Points(float64(fontsize))

	// p_bed.Add(plotter.NewGrid())

	bed_positions, err := plotter.NewScatter(storm_counts_points)
	if err != nil {
		fmt.Printf("Failed to create scatter plot: %v", err)
	}
	bed_positions.GlyphStyle.Shape = draw.CircleGlyph{}
	bed_positions.GlyphStyle.Color = plotter.DefaultLineStyle.Color
	bed_positions.GlyphStyle.Radius = vg.Points(scattersize)

	p_bed.Add(bed_positions)

	p_bed.X.Min = 0
	p_bed.Y.Min = 0
	p_bed.X.Tick.Marker = plot.TickerFunc(func(min, max float64) []plot.Tick {
		// start at the first multiple of spacing ≥ min
		start := math.Ceil(min/spacing) * spacing
		var ticks []plot.Tick
		for v := start; v <= max; v += spacing {
			ticks = append(ticks, plot.Tick{
				Value: v,
				Label: fmt.Sprintf("%.0f", v),
			})
		}
		return ticks
	})

	bedFileName := strings.TrimSuffix(slogs_outfile, ".png") + "_BED.png"
	if err := p_bed.Save(vg.Length(plotxdim)*vg.Inch, vg.Length(plotydim)*vg.Inch, bedFileName); err != nil {
		fmt.Printf("Failed to save plot: %v", err)
	}

	// mean, _ := rollingStats(points_score, averageWindow)

	// meanLineXY := make(plotter.XYs, len(slogs.goodhits))
	// upperXY := make(plotter.XYs, len(slogs.goodhits))
	// lowerXY := make(plotter.XYs, len(slogs.goodhits))
	// for i, pt := range points_score {
	// 	meanLineXY[i].X = pt.X
	// 	meanLineXY[i].Y = mean[i]

	// upperXY[i].X = pt.X
	// upperXY[i].Y = mean[i] + std[i]

	// lowerXY[i].X = pt.X
	// lowerXY[i].Y = mean[i] - std[i]
	// }
	// // build the filled polygon: upper -> reverse(lower)
	// bandXY := make(plotter.XYs, 0, 2*len(slogs.goodhits))
	// bandXY = append(bandXY, upperXY...)
	// for i := len(lowerXY) - 1; i >= 0; i-- {
	// 	bandXY = append(bandXY, lowerXY[i])
	// }

	// chunks_score := makeChunks(points_score, averageWindow)
	// var meanLineXY_score plotter.XYs
	// var upperXY_score, lowerXY_score plotter.XYs
	// for _, c := range chunks_score {
	// 	// duplicate start value to create steps
	// 	meanLineXY_score = append(meanLineXY_score, plotter.XY{X: c.xStart, Y: c.mean})
	// 	meanLineXY_score = append(meanLineXY_score, plotter.XY{X: c.xEnd, Y: c.mean})

	// 	upperXY_score = append(upperXY_score, plotter.XY{X: c.xStart, Y: c.mean + c.sigma})
	// 	upperXY_score = append(upperXY_score, plotter.XY{X: c.xEnd, Y: c.mean + c.sigma})

	// 	lowerXY_score = append(lowerXY_score, plotter.XY{X: c.xStart, Y: c.mean - c.sigma})
	// 	lowerXY_score = append(lowerXY_score, plotter.XY{X: c.xEnd, Y: c.mean - c.sigma})
	// }
	// var bandXY_score plotter.XYs
	// bandXY_score = append(bandXY_score, upperXY_score...)
	// for i := len(lowerXY_score) - 1; i >= 0; i-- {
	// 	bandXY_score = append(bandXY_score, lowerXY_score[i])
	// }

	// positions
	// chunks_positions := makeChunks(points_positions, averageWindow)
	// var meanLineXY_positions plotter.XYs
	// var upperXY_positions, lowerXY_positions plotter.XYs
	// for _, c := range chunks_positions {
	// 	// duplicate start value to create steps
	// 	meanLineXY_positions = append(meanLineXY_positions, plotter.XY{X: c.xStart, Y: c.mean})
	// 	meanLineXY_positions = append(meanLineXY_positions, plotter.XY{X: c.xEnd, Y: c.mean})

	// 	upperXY_positions = append(upperXY_positions, plotter.XY{X: c.xStart, Y: c.mean + c.sigma})
	// 	upperXY_positions = append(upperXY_positions, plotter.XY{X: c.xEnd, Y: c.mean + c.sigma})

	// 	lowerXY_positions = append(lowerXY_positions, plotter.XY{X: c.xStart, Y: c.mean - c.sigma})
	// 	lowerXY_positions = append(lowerXY_positions, plotter.XY{X: c.xEnd, Y: c.mean - c.sigma})
	// }
	// var bandXY_positions plotter.XYs
	// bandXY_positions = append(bandXY_positions, upperXY_positions...)
	// for i := len(lowerXY_positions) - 1; i >= 0; i-- {
	// 	bandXY_positions = append(bandXY_positions, lowerXY_positions[i])
	// }

	// uv := make(plotter.XYs, len(points_positions))
	// for i, pt := range points_positions {
	// 	u := (pt.X + pt.Y) / 2 // along the diagonal
	// 	v := pt.Y - pt.X       // deviation from diagonal
	// 	uv[i].X, uv[i].Y = u, v
	// }

	// maxGap := 1000.0 // choose gap that triggers new contig
	// contigs := splitIntoContigs(points_positions, maxGap)
	// aligned := realign(points_positions, contigs)
	// contigIdx := pointContigIndex(points_positions, contigs)

	// N := len(contigs)
	// colours := distinctColours(N)

	// sc, _ := plotter.NewScatter(points_positions)
	// sc.GlyphStyle.Radius = vg.Points(scattersize)
	// sc.GlyphStyleFunc = func(i int) draw.GlyphStyle {
	// 	gs := sc.GlyphStyle
	// 	gs.Shape = draw.CircleGlyph{}
	// 	gs.Color = colours[contigIdx[i]]
	// 	return gs
	// }

	// p_overlap := plot.New()

	// p_overlap.Title.Text = "Original scatter coloured by contig"
	// p_overlap.X.Label.Text = "Sequence-A position"
	// p_overlap.Y.Label.Text = "Sequence-B position"

	// minY, maxY := points_positions[0].Y, points_positions[0].Y
	// for _, pt := range points_positions {
	// 	if pt.Y < minY {
	// 		minY = pt.Y
	// 	}
	// 	if pt.Y > maxY {
	// 		maxY = pt.Y
	// 	}
	// }

	// for k, c := range contigs {
	// 	polyXY := plotter.XYs{
	// 		{X: c.x0, Y: minY},
	// 		{X: c.x1, Y: minY},
	// 		{X: c.x1, Y: maxY},
	// 		{X: c.x0, Y: maxY},
	// 	}
	// 	poly, _ := plotter.NewPolygon(polyXY)
	// 	// same hue as points, but translucent
	// 	col := colours[k]
	// 	nrgb := color.NRGBAModel.Convert(col).(color.NRGBA)
	// 	poly.Color = color.NRGBA{R: nrgb.R, G: nrgb.G, B: nrgb.B, A: 30}
	// 	poly.LineStyle.Width = 0 // no outline
	// 	p_overlap.Add(poly)
	// }

	// p_overlap.Add(sc)

	// if err := p_overlap.Save(vg.Length(plotxdim)*vg.Inch, vg.Length(plotydim)*vg.Inch, slogs_outfile+"contig_overlap.png"); err != nil {
	// 	fmt.Printf("Failed to create scatter plot: %v", err)
	// }

	//
	// contig plot
	//

	// p_contig := plot.New()

	// p_contig.Title.Text = "Aligned contigs (colour = original Y) Number of Contigs: " + strconv.Itoa(len(contigs))
	// p_contig.X.Label.Text = "Sequence-A position"
	// p_contig.Y.Label.Text = "Sequence-B position (shifted)"

	// // reference diagonal
	// ref_contig, _ := plotter.NewLine(plotter.XYs{{aligned[0].X, aligned[0].X},
	// 	{aligned[len(aligned)-1].X, aligned[len(aligned)-1].X}})
	// ref_contig.Color = color.Gray{Y: 180}
	// ref_contig.Width = vg.Points(5)
	// p_contig.Add(ref_contig)

	// // realigned points
	// sc_contig, _ := plotter.NewScatter(aligned)
	// sc_contig.GlyphStyle.Radius = vg.Points(scattersize)
	// sc_contig.GlyphStyle.Color = color.RGBA{0, 90, 0, 255}

	// minY_contig, maxY_contig := points_positions[0].Y, points_positions[0].Y
	// for _, p := range points_positions {
	// 	if p.Y < minY_contig {
	// 		minY = p.Y
	// 	}
	// 	if p.Y > maxY_contig {
	// 		maxY_contig = p.Y
	// 	}
	// }
	// fmt.Println(len(slogs.goodhits))
	// fmt.Println("--------------Max y: ")
	// fmt.Println(maxY_contig)
	// fmt.Println("--------------Min y: ")
	// fmt.Println(minY_contig)

	// // colour-mapping function  (blue → red)
	// sc_contig.GlyphStyleFunc = func(i int) draw.GlyphStyle {
	// 	yOrig := points_positions[i].Y
	// 	t := (yOrig - minY_contig) / (maxY_contig - minY_contig) // 0..1
	// 	col := color.NRGBA{
	// 		R: uint8(255 * t),
	// 		G: 0,
	// 		B: uint8(255 * (1 - t)),
	// 		A: 200,
	// 	}
	// 	gs := sc_contig.GlyphStyle // copy default
	// 	gs.Color = col
	// 	gs.Shape = draw.CircleGlyph{}
	// 	return gs
	// }

	// p_contig.Add(sc_contig)

	// p_contig.Title.TextStyle.Font.Size = titleSize
	// p_contig.X.Label.TextStyle.Font.Size = titleSize
	// p_contig.Y.Label.TextStyle.Font.Size = titleSize
	// p_contig.Add(plotter.NewGrid())
	// p_contig.X.Tick.Marker = plot.TickerFunc(integerLabels)
	// p_contig.Y.Tick.Marker = plot.TickerFunc(integerLabels)
	// p_contig.X.Tick.Label.Font.Size = vg.Points(float64(fontsize))
	// p_contig.Y.Tick.Label.Font.Size = vg.Points(float64(fontsize))

	// p_contig.X.Min = 0
	// p_contig.Y.Min = 0
	// p_contig.X.Tick.Marker = plot.TickerFunc(func(min, max float64) []plot.Tick {
	// 	// start at the first multiple of spacing ≥ min
	// 	start := math.Ceil(min/spacing) * spacing
	// 	var ticks []plot.Tick
	// 	for v := start; v <= max; v += spacing {
	// 		ticks = append(ticks, plot.Tick{
	// 			Value: v,
	// 			Label: fmt.Sprintf("%.0f", v),
	// 		})
	// 	}
	// 	return ticks
	// })
	// p_contig.Y.Tick.Marker = plot.TickerFunc(func(min, max float64) []plot.Tick {
	// 	// start at the first multiple of spacing ≥ min
	// 	start := math.Ceil(min/spacing) * spacing
	// 	var ticks []plot.Tick
	// 	for v := start; v <= max; v += spacing {
	// 		ticks = append(ticks, plot.Tick{
	// 			Value: v,
	// 			Label: fmt.Sprintf("%.0f", v),
	// 		})
	// 	}
	// 	return ticks
	// })

	// if err := p_contig.Save(vg.Length(plotxdim)*vg.Inch, vg.Length(plotydim)*vg.Inch, slogs_outfile+"contig_aligned.png"); err != nil {
	// 	fmt.Printf("Failed to create scatter plot: %v", err)
	// }

	// //
	// // align position plot
	// //

	// adj := realignBlocks(points_positions, 1000.0)

	// p_align := plot.New()

	// p_align.Title.Text = "Blocks realigned on diagonal (y←y-Δblock)"
	// p_align.X.Label.Text = "Sequence-A position"
	// p_align.Y.Label.Text = "Sequence-B position (shifted)"

	// // reference diagonal y = x
	// ref, _ := plotter.NewLine(plotter.XYs{{adj[0].X, adj[0].X}, {adj[len(adj)-1].X, adj[len(adj)-1].X}})
	// ref.Color = color.Gray{200}
	// ref.Width = vg.Points(0.6)
	// p_align.Add(ref)

	// // realigned points
	// sc_align, _ := plotter.NewScatter(adj)
	// sc_align.GlyphStyle.Radius = vg.Points(scattersize)
	// sc_align.GlyphStyle.Color = color.RGBA{0, 90, 0, 255}
	// p_align.Add(sc_align)
	// p_align.Title.TextStyle.Font.Size = titleSize
	// p_align.X.Label.TextStyle.Font.Size = titleSize
	// p_align.Y.Label.TextStyle.Font.Size = titleSize
	// p_align.Add(plotter.NewGrid())
	// p_align.X.Tick.Marker = plot.TickerFunc(integerLabels)
	// p_align.Y.Tick.Marker = plot.TickerFunc(integerLabels)
	// p_align.X.Tick.Label.Font.Size = vg.Points(float64(fontsize))
	// p_align.Y.Tick.Label.Font.Size = vg.Points(float64(fontsize))

	// p_align.X.Min = 0
	// p_align.Y.Min = 0
	// p_align.X.Tick.Marker = plot.TickerFunc(func(min, max float64) []plot.Tick {
	// 	// start at the first multiple of spacing ≥ min
	// 	start := math.Ceil(min/spacing) * spacing
	// 	var ticks []plot.Tick
	// 	for v := start; v <= max; v += spacing {
	// 		ticks = append(ticks, plot.Tick{
	// 			Value: v,
	// 			Label: fmt.Sprintf("%.0f", v),
	// 		})
	// 	}
	// 	return ticks
	// })
	// p_align.Y.Tick.Marker = plot.TickerFunc(func(min, max float64) []plot.Tick {
	// 	// start at the first multiple of spacing ≥ min
	// 	start := math.Ceil(min/spacing) * spacing
	// 	var ticks []plot.Tick
	// 	for v := start; v <= max; v += spacing {
	// 		ticks = append(ticks, plot.Tick{
	// 			Value: v,
	// 			Label: fmt.Sprintf("%.0f", v),
	// 		})
	// 	}
	// 	return ticks
	// })

	// if err := p_align.Save(vg.Length(plotxdim)*vg.Inch, vg.Length(plotydim)*vg.Inch, slogs_outfile+"aligned_blocks.png"); err != nil {
	// 	fmt.Printf("Failed to create scatter plot: %v", err)
	// }

	//
	// match position plot
	//
	// Create a new plot for positons
	p_positions := plot.New()

	p_positions.Title.Text = "Storm Positions Matches"
	p_positions.X.Label.Text = "Reference" + " Sequence Position"
	p_positions.Y.Label.Text = "Target" + " Sequence Position"
	p_positions.Title.TextStyle.Font.Size = titleSize
	p_positions.X.Label.TextStyle.Font.Size = titleSize
	p_positions.Y.Label.TextStyle.Font.Size = titleSize

	p_positions.X.Tick.Marker = plot.TickerFunc(integerLabels)
	p_positions.Y.Tick.Marker = plot.TickerFunc(integerLabels)
	p_positions.X.Tick.Label.Font.Size = vg.Points(float64(fontsize))
	p_positions.Y.Tick.Label.Font.Size = vg.Points(float64(fontsize))

	// p_positions.Add(plotter.NewGrid())

	scatter_positions, err := plotter.NewScatter(points_positions)
	if err != nil {
		fmt.Printf("Failed to create scatter plot: %v", err)
	}
	scatter_positions.GlyphStyle.Shape = draw.CircleGlyph{}
	scatter_positions.GlyphStyle.Color = plotter.DefaultLineStyle.Color
	scatter_positions.GlyphStyle.Radius = vg.Points(scattersize)

	p_positions.Add(scatter_positions)

	p_positions.X.Min = 0
	p_positions.Y.Min = 0
	p_positions.X.Tick.Marker = plot.TickerFunc(func(min, max float64) []plot.Tick {
		// start at the first multiple of spacing ≥ min
		start := math.Ceil(min/spacing) * spacing
		var ticks []plot.Tick
		for v := start; v <= max; v += spacing {
			ticks = append(ticks, plot.Tick{
				Value: v,
				Label: fmt.Sprintf("%.0f", v),
			})
		}
		return ticks
	})
	p_positions.Y.Tick.Marker = plot.TickerFunc(func(min, max float64) []plot.Tick {
		// start at the first multiple of spacing ≥ min
		start := math.Ceil(min/spacing) * spacing
		var ticks []plot.Tick
		for v := start; v <= max; v += spacing {
			ticks = append(ticks, plot.Tick{
				Value: v,
				Label: fmt.Sprintf("%.0f", v),
			})
		}
		return ticks
	})

	// meanLine_positions, _ := plotter.NewLine(meanLineXY_positions)
	// meanLine_positions.LineStyle.Width = vg.Points(float64(lineWidth))
	// meanLine_positions.LineStyle.Color = color.RGBA{R: 0, G: 90, B: 0, A: 255}
	// p_positions.Add(meanLine_positions)
	// // p_score.Legend.Add("rolling mean", meanLine)

	// //shaded variance band
	// band_positions, _ := plotter.NewPolygon(bandXY_positions)
	// band_positions.Color = color.NRGBA{R: 0, G: 90, B: 0, A: 20} // translucent fill
	// band_positions.LineStyle.Width = 0
	// band_positions.LineStyle.Color = color.NRGBA{R: 0, G: 90, B: 0, A: 255}
	// band_positions.LineStyle.Width = vg.Points(0)
	// p_positions.Add(band_positions)

	//
	// match score plot
	//

	ptsLine := make(plotter.XYs, len(points_score))
	copy(ptsLine, points_score)

	sort.Slice(ptsLine, func(i, j int) bool { return ptsLine[i].X < ptsLine[j].X })

	ln_score, _ := plotter.NewLine(ptsLine)
	ln_score.Color = color.RGBA{R: 0, G: 0, B: 90, A: 255} // same hue as points
	ln_score.Width = vg.Points(float64(lineWidth))

	// Create a new plot for match score
	p_score := plot.New()

	p_score.Title.Text = "Storm Match Count"
	p_score.X.Label.Text = "Reference" + " Sequence Position"
	p_score.Y.Label.Text = "Match Count"
	p_score.Title.TextStyle.Font.Size = titleSize
	p_score.X.Label.TextStyle.Font.Size = titleSize
	p_score.Y.Label.TextStyle.Font.Size = titleSize
	// p_score.Add(plotter.NewGrid())
	p_score.X.Tick.Marker = plot.TickerFunc(integerLabels)
	p_score.Y.Tick.Marker = plot.TickerFunc(integerLabels)
	p_score.X.Tick.Label.Font.Size = vg.Points(float64(fontsize))
	p_score.Y.Tick.Label.Font.Size = vg.Points(float64(fontsize))

	scatter_score, err := plotter.NewScatter(points_score)
	if err != nil {
		fmt.Printf("Failed to create scatter plot: %v", err)
	}
	scatter_score.GlyphStyle.Shape = draw.CircleGlyph{}
	scatter_score.GlyphStyle.Color = plotter.DefaultLineStyle.Color
	scatter_score.GlyphStyle.Radius = vg.Points(scattersize)

	var sum, sumSq float64
	numberOfPoints := float64(len(points_score))
	// if len(points_score) == 0 {
	// 	fmt.Println("no storm-hits to plot")
	// 	return 0.0, 0.0, 0.0, len(points_score), ptsLine
	// }
	maxMatch := points_score[0].Y
	for _, pt := range points_score {
		sum += pt.Y
		sumSq += pt.Y * pt.Y
		if pt.Y > maxMatch {
			maxMatch = pt.Y
		}
	}
	meanY := sum / numberOfPoints
	varPop := (sumSq / numberOfPoints) - meanY*meanY
	varianceMatch := varPop * numberOfPoints / (numberOfPoints - 1) // unbiased sample variance

	meanLineXY_score := plotter.XYs{
		{X: points_score[0].X, Y: meanY},                   // leftmost x
		{X: points_score[len(points_score)-1].X, Y: meanY}, // rightmost x
	}
	meanLine_score, _ := plotter.NewLine(meanLineXY_score)
	meanLine_score.LineStyle.Width = vg.Points(float64(lineWidth))
	meanLine_score.LineStyle.Color = color.NRGBA{R: 255 * 0.8, G: 165 * 0.8, B: 0, A: 255}
	meanLine_score.Dashes = []vg.Length{vg.Points(8), vg.Points(4)}

	p_score.Legend.Top = true           // top    (false → bottom)
	p_score.Legend.Left = false         // right  (true  → left)
	p_score.Legend.XOffs = vg.Points(0) // 5-pt gap from frame (optional)

	// register the line with the legend
	p_score.Legend.Add("Match Count", ln_score)
	p_score.Legend.Add("Mean Match Count", meanLine_score)

	p_score.Add(ln_score)
	p_score.Add(meanLine_score)
	p_score.Add(scatter_score)
	// p_score.Legend.Add("rolling mean", meanLine)

	//shaded variance band
	// band, _ := plotter.NewPolygon(bandXY_score)
	// band.Color = color.NRGBA{R: 0, G: 90, B: 0, A: 20} // translucent fill
	// band.LineStyle.Width = 0
	// band.LineStyle.Color = color.NRGBA{R: 0, G: 90, B: 0, A: 255}
	// band.LineStyle.Width = vg.Points(0)
	// p_score.Add(band)

	p_score.X.Min = 0
	p_score.Y.Min = 0
	p_score.X.Tick.Marker = plot.TickerFunc(func(min, max float64) []plot.Tick {
		// start at the first multiple of spacing ≥ min
		start := math.Ceil(min/spacing) * spacing
		var ticks []plot.Tick
		for v := start; v <= max; v += spacing {
			ticks = append(ticks, plot.Tick{
				Value: v,
				Label: fmt.Sprintf("%.0f", v),
			})
		}
		return ticks
	})

	// Save the plot to a PNG file
	// os.MkdirAll(kmerPlotDir, os.ModePerm)
	if err := p_positions.Save(vg.Length(plotxdim)*vg.Inch, vg.Length(plotydim)*vg.Inch, slogs_outfile); err != nil {
		fmt.Printf("Failed to save plot: %v", err)
	}

	// Save the plot to a PNG file
	// os.MkdirAll(kmerPlotDir, os.ModePerm)
	kmerPlotFileName := strings.TrimSuffix(slogs_outfile, ".png") + "_matchscore.png"
	fmt.Println("Plotting storm match count to ", kmerPlotFileName)
	if err := p_score.Save(vg.Length(plotxdim)*vg.Inch, vg.Length(plotydim)*vg.Inch, kmerPlotFileName); err != nil {
		fmt.Printf("Failed to save plot: %v", err)
	}

	return meanY, varianceMatch, maxMatch, len(points_score), points_score, storm_counts_points

}

type weatherSystem struct {
	runtime            float64
	referenceLength    float64
	targetLength       float64
	numberOfStorms     float64
	meanDivergence     float64
	varianceDivergence float64
	matchPlotPoints    plotter.XYs
}

type WeatherSystem struct {
	RefName                     string
	TargetName                  string
	MatchPlotPoints             plotter.XYs // points that get plotted in Match Score Plot (x,y)
	BedPlotPoints               plotter.XYs // points that get plotted in Match Score Plot (x,y)
	TimeKmerize                 float64     // use time.Duration? variable.Seconds() gives a float64
	TimeStormMatch              float64     // variable.Seconds() gives a float64
	TimeStormPlot               float64     // variable.Seconds() gives a float64
	MeanMatchScore              float64
	VarianceMatchScore          float64
	CorrectedMeanMatchScore     float64
	CorrectedVarianceMatchScore float64
	CorrectedStormHits          float64
	StormsFound                 int
	MaxMatchScore               float64
	NumberOfMatchPoints         int
	SequenceLength              int
	Stormsize                   int
	StormCount                  int
}

func (ws *WeatherSystem) Init(
	refName, targetName string,
	matchPlotPoints, bedPlotPoints plotter.XYs,
	timeKmerize, timeStormMatch, timeStormPlot float64,
	meanMatchScore, varianceMatchScore, maxMatchScore float64,
	numberOfMatchPoints, sequenceLength, stormsize, stormcount int,
) {

	ws.RefName = refName
	ws.TargetName = targetName
	ws.MatchPlotPoints = matchPlotPoints
	ws.BedPlotPoints = bedPlotPoints
	ws.TimeKmerize = timeKmerize
	ws.TimeStormMatch = timeStormMatch
	ws.TimeStormPlot = timeStormPlot
	ws.MeanMatchScore = meanMatchScore
	ws.VarianceMatchScore = varianceMatchScore
	ws.CorrectedMeanMatchScore = 0.0
	ws.CorrectedVarianceMatchScore = 0.0
	ws.MaxMatchScore = maxMatchScore
	ws.CorrectedStormHits = 0
	ws.StormsFound = 0
	ws.NumberOfMatchPoints = numberOfMatchPoints
	ws.SequenceLength = sequenceLength
	ws.Stormsize = stormsize
	ws.StormCount = stormcount
}

func (ws *WeatherSystem) Print() {
	fmt.Println("------------------------------------------------------------")
	fmt.Printf(" reference sequence   : %s\n", ws.RefName)
	fmt.Printf(" target sequence      : %s\n", ws.TargetName)
	fmt.Printf(" sequence length      : %d bp\n", ws.SequenceLength)
	fmt.Printf(" match points         : %d\n", ws.NumberOfMatchPoints)
	fmt.Printf(" mean match score     : %.4f\n", ws.MeanMatchScore)
	fmt.Printf(" variance match score : %.4f\n", ws.VarianceMatchScore)
	fmt.Printf(" max match score      : %.4f\n", ws.MaxMatchScore)
	fmt.Printf(" time k-merize        : %v\n", ws.TimeKmerize)
	fmt.Printf(" time storm-match     : %v\n", ws.TimeStormMatch)
	fmt.Printf(" time storm-plot      : %v\n", ws.TimeStormPlot)
	fmt.Printf(" Sequence Length      : %d\n", ws.SequenceLength)
	fmt.Printf(" Storm size           : %d\n", ws.Stormsize)
	fmt.Printf(" Storm count          : %d\n", ws.StormCount)
	fmt.Println("------------------------------------------------------------")
}

type WeatherForecast struct {
	Systems []*WeatherSystem
	Count   int
}

// Add attaches one system to the forecast
func (wf *WeatherForecast) Add(ws *WeatherSystem) {
	wf.Systems = append(wf.Systems, ws)
	wf.Count += 1
}

func (wf *WeatherForecast) PlotWeatherForecast(outfile string, scattersize float64, plotxdim int, plotydim int, fontsize int, titleFontSize int, averageWindow int, lineWidth int, spacing float64) {

	wf.PlotKmerizeVsLength(outfile, scattersize, plotxdim, plotydim, fontsize, titleFontSize, averageWindow, lineWidth, spacing)
	wf.PlotStormTimeVsLength(outfile, scattersize, plotxdim, plotydim, fontsize, titleFontSize, averageWindow, lineWidth, spacing)
	wf.PlotMeanScoreHeatmap(outfile, scattersize, plotxdim, plotydim, fontsize, titleFontSize, averageWindow, lineWidth, spacing)
	wf.PlotOverlayMatchPoints(outfile, scattersize, plotxdim, plotydim, fontsize, titleFontSize, averageWindow, lineWidth, spacing)
}

func (wf *WeatherForecast) PlotKmerizeVsLength(outfile string, scattersize float64, plotxdim int, plotydim int, fontsize int, titleFontSize int, averageWindow int, lineWidth int, spacing float64) {
	// Set up style

	titleSize := vg.Points(float64(titleFontSize))
	integerLabels := func(min, max float64) []plot.Tick {
		// Start from the default ticks so we keep the same positions.
		ticks := plot.DefaultTicks{}.Ticks(min, max)

		for i := range ticks {
			ticks[i].Label = fmt.Sprintf("%.0f", ticks[i].Value) // no decimals
		}
		return ticks
	}

	//
	// Plot the sequence lengths against Time to Kmerize
	var pts plotter.XYs
	for _, ws := range wf.Systems {
		if ws == nil {
			continue
		}
		x := float64(ws.SequenceLength)
		y := ws.TimeKmerize
		pts = append(pts, plotter.XY{X: x, Y: y})
	}

	sort.Slice(pts, func(i, j int) bool { return pts[i].X < pts[j].X })

	sc, _ := plotter.NewScatter(pts)

	sc.GlyphStyle.Shape = draw.CircleGlyph{}
	sc.GlyphStyle.Radius = vg.Points(scattersize)
	sc.GlyphStyle.Color = plotter.DefaultLineStyle.Color

	p_kmerize := plot.New()

	p_kmerize.Title.Text = "Time to k-merize vs sequence length"
	p_kmerize.X.Label.Text = "Sequence length (bp)"
	p_kmerize.Y.Label.Text = "Time k-merize (s)"
	p_kmerize.Title.TextStyle.Font.Size = titleSize
	p_kmerize.X.Label.TextStyle.Font.Size = titleSize
	p_kmerize.Y.Label.TextStyle.Font.Size = titleSize
	p_kmerize.X.Tick.Marker = plot.TickerFunc(integerLabels)
	p_kmerize.Y.Tick.Marker = plot.TickerFunc(integerLabels)
	p_kmerize.X.Tick.Label.Font.Size = vg.Points(float64(fontsize))
	p_kmerize.Y.Tick.Label.Font.Size = vg.Points(float64(fontsize))

	// plot a diagonal y=x line
	maxX := pts[0].X
	for _, pt := range pts[1:] {
		if pt.X > maxX {
			maxX = pt.X
		}
	}

	// draw y = ax over the same span
	var sumXY, sumX2 float64
	minX, maxX := pts[0].X, pts[0].X
	for _, p := range pts {
		sumXY += p.X * p.Y
		sumX2 += p.X * p.X
		if p.X < minX {
			minX = p.X
		}
		if p.X > maxX {
			maxX = p.X
		}
	}

	coeff := sumXY / sumX2

	lineXY := plotter.XYs{
		{X: 0.0, Y: 0.0},
		{X: maxX, Y: coeff * maxX},
	}
	ref, _ := plotter.NewLine(lineXY)
	ref.Color = color.Gray{160} // light-grey
	ref.Width = vg.Points(1.4)
	ref.Dashes = []vg.Length{vg.Points(6), vg.Points(4)}

	p_kmerize.Add(ref)
	p_kmerize.Legend.Add(fmt.Sprintf("y = %.3g·x", coeff), ref)

	p_kmerize.Add(sc)

	p_kmerize.X.Min = 0
	p_kmerize.Y.Min = 0
	p_kmerize.X.Tick.Marker = plot.TickerFunc(func(min, max float64) []plot.Tick {
		// start at the first multiple of spacing ≥ min
		start := math.Ceil(min/spacing) * spacing
		var ticks []plot.Tick
		for v := start; v <= max; v += spacing {
			ticks = append(ticks, plot.Tick{
				Value: v,
				Label: fmt.Sprintf("%.0f", v),
			})
		}
		return ticks
	})

	p_kmerize.Y.Tick.Marker = plot.TickerFunc(func(min, max float64) []plot.Tick {

		ticks := plot.DefaultTicks{}.Ticks(min, max)

		for i := range ticks {
			ticks[i].Label = fmt.Sprintf("%.3f", ticks[i].Value) // 3 decimals
		}
		return ticks
	})

	if err := p_kmerize.Save(vg.Length(plotxdim)*vg.Inch, vg.Length(plotydim)*vg.Inch, outfile+"kmerizeTime.PNG"); err != nil {
		fmt.Printf("Failed to save plot: %v", err)
	}

}

func (wf *WeatherForecast) PlotStormTimeVsLength(outfile string, scattersize float64, plotxdim int, plotydim int, fontsize int, titleFontSize int, averageWindow int, lineWidth int, spacing float64) {
	// Set up style

	titleSize := vg.Points(float64(titleFontSize))
	integerLabels := func(min, max float64) []plot.Tick {
		// Start from the default ticks so we keep the same positions.
		ticks := plot.DefaultTicks{}.Ticks(min, max)

		for i := range ticks {
			ticks[i].Label = fmt.Sprintf("%.0f", ticks[i].Value) // no decimals
		}
		return ticks
	}

	//
	// Plot the sequence lengths against Time to Kmerize
	var pts plotter.XYs
	for _, ws := range wf.Systems {
		if ws == nil {
			continue
		}
		x := float64(ws.SequenceLength)
		y := ws.TimeStormMatch
		pts = append(pts, plotter.XY{X: x, Y: y})
	}

	sort.Slice(pts, func(i, j int) bool { return pts[i].X < pts[j].X })

	sc, _ := plotter.NewScatter(pts)

	sc.GlyphStyle.Shape = draw.CircleGlyph{}
	sc.GlyphStyle.Radius = vg.Points(scattersize)
	sc.GlyphStyle.Color = plotter.DefaultLineStyle.Color

	p_stormtime := plot.New()

	p_stormtime.Title.Text = "Time to Storm Match vs sequence length"
	p_stormtime.X.Label.Text = "Sequence length (bp)"
	p_stormtime.Y.Label.Text = "Time Storm Match (s)"
	p_stormtime.Title.TextStyle.Font.Size = titleSize
	p_stormtime.X.Label.TextStyle.Font.Size = titleSize
	p_stormtime.Y.Label.TextStyle.Font.Size = titleSize
	p_stormtime.X.Tick.Marker = plot.TickerFunc(integerLabels)
	p_stormtime.Y.Tick.Marker = plot.TickerFunc(integerLabels)
	p_stormtime.X.Tick.Label.Font.Size = vg.Points(float64(fontsize))
	p_stormtime.Y.Tick.Label.Font.Size = vg.Points(float64(fontsize))

	// plot a diagonal y=x line
	maxX := pts[0].X
	for _, pt := range pts[1:] {
		if pt.X > maxX {
			maxX = pt.X
		}
	}

	// draw y = ax over the same span
	var sumXY, sumX2 float64
	minX, maxX := pts[0].X, pts[0].X
	for _, p := range pts {
		sumXY += p.X * p.Y
		sumX2 += p.X * p.X
		if p.X < minX {
			minX = p.X
		}
		if p.X > maxX {
			maxX = p.X
		}
	}

	coeff := sumXY / sumX2

	lineXY := plotter.XYs{
		{X: 0.0, Y: 0.0},
		{X: maxX, Y: coeff * maxX},
	}
	ref, _ := plotter.NewLine(lineXY)
	ref.Color = color.Gray{160} // light-grey
	ref.Width = vg.Points(1.4)
	ref.Dashes = []vg.Length{vg.Points(6), vg.Points(4)}

	p_stormtime.Add(ref)
	p_stormtime.Legend.Add(fmt.Sprintf("y = %.3g·x", coeff), ref)

	p_stormtime.Add(sc)

	p_stormtime.X.Min = 0
	p_stormtime.Y.Min = 0
	p_stormtime.X.Tick.Marker = plot.TickerFunc(func(min, max float64) []plot.Tick {
		// start at the first multiple of spacing ≥ min
		start := math.Ceil(min/spacing) * spacing
		var ticks []plot.Tick
		for v := start; v <= max; v += spacing {
			ticks = append(ticks, plot.Tick{
				Value: v,
				Label: fmt.Sprintf("%.0f", v),
			})
		}
		return ticks
	})
	p_stormtime.Y.Tick.Marker = plot.TickerFunc(func(min, max float64) []plot.Tick {

		ticks := plot.DefaultTicks{}.Ticks(min, max)

		for i := range ticks {
			ticks[i].Label = fmt.Sprintf("%.3f", ticks[i].Value) // 3 decimals
		}
		return ticks
	})

	if err := p_stormtime.Save(vg.Length(plotxdim)*vg.Inch, vg.Length(plotydim)*vg.Inch, outfile+"stormMatchTime.PNG"); err != nil {
		fmt.Printf("Failed to save plot: %v", err)
	}

}

type grayScale struct{ n int } // n = how many colours you want

func (p grayScale) Colors() []color.Color {
	if p.n <= 0 {
		return nil
	}
	clrs := make([]color.Color, p.n)
	for i := 0; i < p.n; i++ {
		v := uint8(float64(p.n-1-i) * 255 / float64(p.n-1))
		clrs[i] = color.NRGBA{R: v, G: v, B: v, A: 255}
	}
	return clrs
}

type matchGrid struct {
	refs, targets []string
	score         [][]float64 // score[col][row]
}

func (g *matchGrid) Dims() (c, r int)   { return len(g.targets), len(g.refs) }
func (g *matchGrid) Z(c, r int) float64 { return g.score[c][r] }
func (g *matchGrid) X(c int) float64    { return float64(c) }
func (g *matchGrid) Y(r int) float64    { return float64(r) }

func buildGrid(wf *WeatherForecast) *matchGrid {
	if wf == nil {
		return nil
	}

	// collect unique names
	refSet := make(map[string]bool)
	tgtSet := make(map[string]bool)
	for _, ws := range wf.Systems {
		if ws == nil {
			continue
		}
		refSet[ws.RefName] = true
		tgtSet[ws.TargetName] = true
	}

	// sort for deterministic order
	var refs, tgts []string
	for n := range refSet {
		refs = append(refs, n)
	}
	for n := range tgtSet {
		tgts = append(tgts, n)
	}
	sort.Strings(refs)
	sort.Strings(tgts)

	// index helpers
	refIdx := make(map[string]int)
	tgtIdx := make(map[string]int)
	for i, n := range refs {
		refIdx[n] = i
	}
	for j, n := range tgts {
		tgtIdx[n] = j
	}

	// fill matrix with NaN (so missing combinations are blank)
	score := make([][]float64, len(tgts))
	for j := range score {
		score[j] = make([]float64, len(refs))
		for i := range score[j] {
			score[j][i] = math.NaN()
		}
	}

	// insert real scores
	for _, ws := range wf.Systems {
		if ws == nil {
			continue
		}
		r := refIdx[ws.RefName]
		c := tgtIdx[ws.TargetName]
		score[c][r] = ws.MeanMatchScore
	}

	return &matchGrid{refs: refs, targets: tgts, score: score}
}

func (wf *WeatherForecast) PlotMeanScoreHeatmap(outfile string, scattersize float64, plotxdim int, plotydim int, fontsize int, titleFontSize int, averageWindow int, lineWidth int, spacing float64) {
	grid := buildGrid(wf)
	if grid == nil || len(grid.refs) == 0 || len(grid.targets) == 0 {
		return
	}

	// colour map
	pal := grayScale{n: 256}
	hm := plotter.NewHeatMap(grid, pal)

	titleSize := vg.Points(float64(titleFontSize))
	integerLabels := func(min, max float64) []plot.Tick {
		// Start from the default ticks so we keep the same positions.
		ticks := plot.DefaultTicks{}.Ticks(min, max)

		for i := range ticks {
			ticks[i].Label = fmt.Sprintf("%.0f", ticks[i].Value) // no decimals
		}
		return ticks
	}

	// build plot
	p := plot.New()
	p.Title.Text = "Mean match score (ref rows × target cols)"
	p.X.Label.Text = "Target Sequence"
	p.Y.Label.Text = "Reference Storm"
	p.Title.TextStyle.Font.Size = titleSize
	p.X.Label.TextStyle.Font.Size = titleSize
	p.Y.Label.TextStyle.Font.Size = titleSize
	p.X.Tick.Marker = plot.TickerFunc(integerLabels)
	p.Y.Tick.Marker = plot.TickerFunc(integerLabels)
	p.X.Tick.Label.Font.Size = vg.Points(float64(fontsize))
	p.Y.Tick.Label.Font.Size = vg.Points(float64(fontsize))

	// custom tickers that use sequence names
	p.X.Tick.Marker = plot.TickerFunc(func(min, max float64) []plot.Tick {
		var ticks []plot.Tick
		for c, name := range grid.targets {
			ticks = append(ticks, plot.Tick{Value: float64(c), Label: name})
		}
		return ticks
	})
	p.Y.Tick.Marker = plot.TickerFunc(func(min, max float64) []plot.Tick {
		var ticks []plot.Tick
		for r, name := range grid.refs {
			ticks = append(ticks, plot.Tick{Value: float64(r), Label: name})
		}
		return ticks
	})

	p.Add(hm)

	p.X.Tick.Label.Rotation = math.Pi / 2
	p.X.Tick.Label.XAlign = draw.XRight
	p.X.Tick.Label.YAlign = draw.YCenter

	// colour-bar
	// cb := plotter.NewColorBar(hm.Palette)
	// cb.Vertical = true
	// p.Y.Right = true // enable right axis so the color bar can sit there
	// p.Add(cb)

	if err := p.Save(vg.Length(plotxdim)*vg.Inch, vg.Length(plotydim)*vg.Inch, outfile+"heatmap.PNG"); err != nil {
		fmt.Printf("Failed to save plot: %v", err)
	}
	return
}

func makeColours(N int) []color.Color {
	colours := make([]color.Color, N)
	for i := 0; i < N; i++ {
		h := float64(i) / float64(N) // 0 .. <1
		r, g, b := HSVtoRGB(h, 0.75, 0.95)
		colours[i] = color.NRGBA{
			R: uint8(r * 255),
			G: uint8(g * 255),
			B: uint8(b * 255),
			A: 255,
		}
	}
	return colours
}

func HSVtoRGB(h, s, v float64) (r, g, b float64) {
	if s == 0 {
		return v, v, v // grey
	}
	h = math.Mod(h*6, 6)
	i := math.Floor(h)
	f := h - i
	p := v * (1 - s)
	q := v * (1 - s*f)
	t := v * (1 - s*(1-f))
	switch int(i) {
	case 0:
		r, g, b = v, t, p
	case 1:
		r, g, b = q, v, p
	case 2:
		r, g, b = p, v, t
	case 3:
		r, g, b = p, q, v
	case 4:
		r, g, b = t, p, v
	default: // 5
		r, g, b = v, p, q
	}
	return
}

func (wf *WeatherForecast) PlotOverlayMatchPoints(outfile string, scattersize float64, plotxdim int, plotydim int, fontsize int, titleFontSize int, averageWindow int, lineWidth int, spacing float64) {

	refMap := make(map[string][]*WeatherSystem)
	for _, ws := range wf.Systems {
		if ws == nil {
			continue
		}
		refMap[ws.RefName] = append(refMap[ws.RefName], ws)
	}

	for ref, _ := range refMap {

		var systems []*WeatherSystem
		for _, ws := range wf.Systems {
			if ws == nil {
				continue // no data
			}
			if ws.RefName == ws.TargetName { // ← skip ref==target
				continue
			}
			if len(ws.MatchPlotPoints) == 0 {
				continue
			}
			systems = append(systems, ws)
		}

		// one colour per WeatherSystem
		colours := makeColours(len(systems))

		titleSize := vg.Points(float64(titleFontSize))
		integerLabels := func(min, max float64) []plot.Tick {
			// Start from the default ticks so we keep the same positions.
			ticks := plot.DefaultTicks{}.Ticks(min, max)

			for i := range ticks {
				ticks[i].Label = fmt.Sprintf("%.0f", ticks[i].Value) // no decimals
			}
			return ticks
		}

		// build plot
		p := plot.New()

		p.Title.Text = fmt.Sprintf("Matches vs %s", ref)
		p.X.Label.Text = "Reference sequence position"
		p.Y.Label.Text = "Match Score"
		// p.Legend.Top = true
		// p.Legend.Left = false
		// p.Legend.XOffs = vg.Points(5)
		p.Title.TextStyle.Font.Size = titleSize
		p.X.Label.TextStyle.Font.Size = titleSize
		p.Y.Label.TextStyle.Font.Size = titleSize
		p.X.Tick.Marker = plot.TickerFunc(integerLabels)
		p.Y.Tick.Marker = plot.TickerFunc(integerLabels)
		p.X.Tick.Label.Font.Size = vg.Points(float64(fontsize))
		p.Y.Tick.Label.Font.Size = vg.Points(float64(fontsize))

		var allPoints plotter.XYs

		for i, ws := range systems {
			if ws == nil || len(ws.MatchPlotPoints) == 0 {
				continue
			}
			allPoints = append(allPoints, ws.MatchPlotPoints...)
			pts := make(plotter.XYs, len(ws.MatchPlotPoints))
			copy(pts, ws.MatchPlotPoints) // keep original order intact

			sort.Slice(pts, func(a, b int) bool { return pts[a].X < pts[b].X })

			// scatter
			sc, _ := plotter.NewScatter(pts)
			sc.GlyphStyle.Shape = draw.CircleGlyph{}
			sc.GlyphStyle.Radius = vg.Points(2.5)
			sc.GlyphStyle.Color = colours[i]

			// thin line through the points (so cloud shape is visible)
			ln, _ := plotter.NewLine(pts)
			ln.Color = colours[i]
			ln.Width = vg.Points(1)

			// register in legend using the target-sequence name
			// p.Legend.Add(ws.TargetName, sc)

			// add to plot (line first, then scatter on top)
			p.Add(ln, sc)
		}

		if err := p.Save(vg.Length(plotxdim)*vg.Inch, vg.Length(plotydim)*vg.Inch, outfile+"match_overlay"+ref+".PNG"); err != nil {
			fmt.Printf("Failed to save plot: %v", err)
		}

		// compute ensemble statistics
		binSize := float64(1000)
		meanXY, upperXY, lowerXY := binStats(allPoints, binSize)

		// build mean line
		meanLine, _ := plotter.NewLine(meanXY)
		meanLine.Color = color.RGBA{R: 255, G: 140, B: 0, A: 255} // dark orange
		meanLine.Width = vg.Points(1.5)

		// variance band: upper -> reversed(lower)
		var bandXY plotter.XYs
		bandXY = append(bandXY, upperXY...)
		for i := len(lowerXY) - 1; i >= 0; i-- {
			bandXY = append(bandXY, lowerXY[i])
		}
		band, _ := plotter.NewPolygon(bandXY)
		band.Color = color.NRGBA{R: 255, G: 165, B: 0, A: 60} // translucent orange
		band.LineStyle.Width = 0

		// create plot
		p_mean := plot.New()

		p_mean.Title.Text = fmt.Sprintf("Mean Matches vs %s", ref)
		p_mean.X.Label.Text = "Reference sequence position"
		p_mean.Y.Label.Text = "Match Score"
		p_mean.Title.TextStyle.Font.Size = titleSize
		p_mean.X.Label.TextStyle.Font.Size = titleSize
		p_mean.Y.Label.TextStyle.Font.Size = titleSize
		p_mean.X.Tick.Marker = plot.TickerFunc(integerLabels)
		p_mean.Y.Tick.Marker = plot.TickerFunc(integerLabels)
		p_mean.X.Tick.Label.Font.Size = vg.Points(float64(fontsize))
		p_mean.Y.Tick.Label.Font.Size = vg.Points(float64(fontsize))

		p_mean.Add(band)
		p_mean.Add(meanLine)

		if err := p_mean.Save(vg.Length(plotxdim)*vg.Inch, vg.Length(plotydim)*vg.Inch, outfile+"mean_match_overlay"+ref+".PNG"); err != nil {
			fmt.Printf("Failed to save plot: %v", err)
		}

	}

}

func binStats(pts plotter.XYs, binSize float64) (meanXY, upperXY, lowerXY plotter.XYs) {
	if len(pts) == 0 {
		return
	}
	// make sure points are sorted by X so bins are contiguous
	sort.Slice(pts, func(i, j int) bool { return pts[i].X < pts[j].X })

	binStart := math.Floor(pts[0].X/binSize) * binSize
	var sum, sumSq float64
	var n int
	for _, p := range pts {
		if p.X-binStart >= binSize {
			if n > 0 {
				m := sum / float64(n)
				variance := sumSq/float64(n) - m*m
				if variance < 0 {
					variance = 0
				}
				sigma := math.Sqrt(variance)

				mid := binStart + binSize/2
				meanXY = append(meanXY, plotter.XY{X: mid, Y: m})
				upperXY = append(upperXY, plotter.XY{X: mid, Y: m + sigma})
				lowerXY = append(lowerXY, plotter.XY{X: mid, Y: m - sigma})
			}
			// reset for next bin
			binStart = math.Floor(p.X/binSize) * binSize
			sum, sumSq, n = 0, 0, 0
		}
		sum += p.Y
		sumSq += p.Y * p.Y
		n++
	}
	// flush last bin
	if n > 0 {
		m := sum / float64(n)
		variance := sumSq/float64(n) - m*m
		if variance < 0 {
			variance = 0
		}
		sigma := math.Sqrt(variance)
		mid := binStart + binSize/2
		meanXY = append(meanXY, plotter.XY{X: mid, Y: m})
		upperXY = append(upperXY, plotter.XY{X: mid, Y: m + sigma})
		lowerXY = append(lowerXY, plotter.XY{X: mid, Y: m - sigma})
	}
	return
}

func refNameOf(fc *WeatherForecast) string {
	for _, ws := range fc.Systems {
		if ws != nil && ws.RefName != "" {
			return ws.RefName
		}
	}
	return "(unknown)"
}

func (wf *WeatherForecast) Print() {
	fmt.Println()
	fmt.Printf("WeatherForecast  —  %d system(s)\n", len(wf.Systems))
	fmt.Println("============================================================")
	for i, ws := range wf.Systems {
		if ws == nil {
			fmt.Printf("[%d] <nil entry>\n", i)
			continue
		}
		fmt.Printf("[%d] ", i)
		ws.Print()
	}
}

func (wf *WeatherForecast) Save(outfile string) {
	f, _ := os.Create(outfile + "weather_forecast.json")

	defer f.Close()

	enc := json.NewEncoder(f)
	enc.SetIndent("", "  ") // pretty-print
	enc.Encode(wf)          // Write `wf` as JSON

}

func Load(infile string) (*WeatherForecast, error) {
	f, err := os.Open(infile)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	var wf WeatherForecast
	if err := json.NewDecoder(f).Decode(&wf); err != nil {
		return nil, err
	}
	return &wf, nil
}

func correctedStats(ptsMatch plotter.XYs, stormSize float64) (mean, variance, n float64) {
	var sum, sumSq float64
	for _, p := range ptsMatch {
		if p.Y <= stormSize {
			sum += p.Y
			sumSq += p.Y * p.Y
			n += 1.0
		}
	}
	if n < 1.0 {
		return 0.0, 0.0, 0.0
	}
	mean = sum / float64(n)
	variance = sumSq/float64(n) - mean*mean // population variance
	if variance < 0 {
		variance = 0
	}
	return
}

func UpdateForecastsInOutput(outputDir string, stormSize float64) error {
	dirEntries, err := os.ReadDir(outputDir)
	if err != nil {
		return err
	}

	for _, entry := range dirEntries {
		if !entry.IsDir() {
			continue
		}

		jsonPath := filepath.Join(outputDir, entry.Name(), "weather_forecast.json")
		fmt.Println(jsonPath)
		if _, err := os.Stat(jsonPath); err != nil {
			continue // file not found – skip silently
		}

		fc, err := Load(jsonPath) // ← your existing decoder
		if err != nil {
			continue
			// return fmt.Errorf("decode %s: %w", jsonPath, err)
		}
		fmt.Println("here 0")
		for _, ws := range fc.Systems {
			if ws == nil {
				continue
			}

			mean, varp, corrected_hits := correctedStats(ws.MatchPlotPoints, stormSize)

			ws.CorrectedMeanMatchScore = mean
			ws.CorrectedVarianceMatchScore = varp
			ws.CorrectedStormHits = corrected_hits
			ws.StormsFound = len(ws.BedPlotPoints)
			// fmt.Println("here out")
			// fmt.Println(ws.CorrectedMeanMatchScore)
			// fmt.Println(ws.CorrectedVarianceMatchScore)
			// fmt.Println(ws.CorrectedStormHits)
		}
		fmt.Println("here 1")
		// overwrite the JSON file (pretty)
		tmp := jsonPath + ".tmp"
		f, err := os.Create(tmp)
		if err != nil {
			continue
		}
		fmt.Println("here 2")
		enc := json.NewEncoder(f)
		enc.SetIndent("", "  ")
		fmt.Println("here 4")
		if err := enc.Encode(fc); err != nil {
			fmt.Println("here 5")
			f.Close()
			continue
		}
		f.Close()
		fmt.Println("here 3")
		if err := os.Rename(tmp, jsonPath); err != nil {
			continue
		}
		fmt.Printf("updated %s\n", jsonPath)
	}
	return nil
}

type WeatherPortfolio struct {
	Forecasts []*WeatherForecast
	Count     int
}

func NewWeatherPortfolio() *WeatherPortfolio {
	return &WeatherPortfolio{
		Forecasts: make([]*WeatherForecast, 0),
		Count:     0,
	}
}

func (wp *WeatherPortfolio) Add(fc *WeatherForecast) {
	if fc == nil {
		return
	}
	wp.Forecasts = append(wp.Forecasts, fc)
	wp.Count = len(wp.Forecasts)
}

func LoadForecastsFromDir(rootDir string) (*WeatherPortfolio, error) {
	var jsonFiles []string
	err := filepath.WalkDir(rootDir, func(path string, d fs.DirEntry, walkErr error) error {
		if walkErr != nil { // abort traversal on I/O error
			return walkErr
		}
		if d.IsDir() || filepath.Ext(path) != ".json" {
			return nil // ignore directories and non-json files
		}
		jsonFiles = append(jsonFiles, path)
		return nil
	})
	if err != nil {
		return nil, err
	}

	wp := NewWeatherPortfolio()

	for _, file := range jsonFiles {
		fc, err := Load(file) // your existing JSON-decoder
		if err != nil {
			fmt.Println("error with file: %w ", err)
			fmt.Println(file)
			// return nil, fmt.Errorf("load %s: %w", file, err)
		} else {
			// skip file if there is an issue loading weather
			wp.Add(fc)
			fmt.Println(file)
		}

	}
	return wp, nil
}

type gridMany struct {
	rows   []string    // one label per forecast
	cols   []string    // all distinct target names
	values [][]float64 // values[col][row] == score   (NaN = missing)
}

func referenceName(fc *WeatherForecast) string {
	for _, ws := range fc.Systems {
		if ws != nil && ws.RefName != "" {
			return ws.RefName
		}
	}
	return "(unknown)"
}

func (g *gridMany) Dims() (c, r int)   { return len(g.cols), len(g.rows) }
func (g *gridMany) Z(c, r int) float64 { return g.values[c][r] }
func (g *gridMany) X(c int) float64    { return float64(c) } // column index
func (g *gridMany) Y(r int) float64    { return float64(r) } // row index

func buildManyGrid(forecasts []*WeatherForecast) *gridMany {
	if len(forecasts) == 0 {
		return nil
	}

	// collect unique target names
	colSet := make(map[string]bool)
	for _, fc := range forecasts {
		for _, ws := range fc.Systems {
			if ws != nil {
				colSet[ws.TargetName] = true
			}
		}
	}
	if len(colSet) == 0 {
		return nil
	}
	cols := make([]string, 0, len(colSet))
	for n := range colSet {
		cols = append(cols, n)
	}
	sort.Strings(cols)

	// index helpers
	colIdx := make(map[string]int, len(cols))
	for i, n := range cols {
		colIdx[n] = i
	}

	// build matrix initialised with NaN
	values := make([][]float64, len(cols))
	for i := range values {
		values[i] = make([]float64, len(forecasts))
		for j := range values[i] {
			values[i][j] = math.NaN()
		}
	}

	// fill matrix
	for r, fc := range forecasts {
		for _, ws := range fc.Systems {
			if ws == nil {
				continue
			}
			c := colIdx[ws.TargetName]
			values[c][r] = ws.MeanMatchScore
		}
	}

	// rows are now reference names (duplicates get suffix _2, _3, …)
	rows := make([]string, len(forecasts))
	seen := make(map[string]int)

	for i, fc := range forecasts {
		name := refNameOf(fc)
		if cnt, ok := seen[name]; ok { // duplicate → make unique
			cnt++
			name = fmt.Sprintf("%s_%d", name, cnt)
			seen[refNameOf(fc)] = cnt
		} else {
			seen[name] = 0
		}
		rows[i] = name
	}
	return &gridMany{rows: rows, cols: cols, values: values}
}

func buildManyGridHitFrequency(forecasts []*WeatherForecast) *gridMany {

	fmt.Println("++++++++++++ buildManyGridHitFrequency ++++++++++++++")

	if len(forecasts) == 0 {
		return nil
	}

	// collect unique target names
	colSet := make(map[string]bool)
	for _, fc := range forecasts {
		for _, ws := range fc.Systems {
			if ws != nil {
				colSet[ws.TargetName] = true
			}
		}
	}
	if len(colSet) == 0 {
		return nil
	}
	cols := make([]string, 0, len(colSet))
	for n := range colSet {
		cols = append(cols, n)
	}
	sort.Strings(cols)

	// index helpers
	colIdx := make(map[string]int, len(cols))
	for i, n := range cols {
		colIdx[n] = i
	}

	// build matrix initialised with NaN
	values := make([][]float64, len(cols))
	for i := range values {
		values[i] = make([]float64, len(forecasts))
		for j := range values[i] {
			values[i][j] = math.NaN()
		}
	}

	// fill matrix
	for r, fc := range forecasts {
		for _, ws := range fc.Systems {
			if ws == nil {
				continue
			}
			c := colIdx[ws.TargetName]
			fmt.Println(float64(ws.StormsFound) / float64(ws.StormCount))
			values[c][r] = float64(ws.StormsFound) / float64(ws.StormCount)
		}
	}

	// rows are now reference names (duplicates get suffix _2, _3, …)
	rows := make([]string, len(forecasts))
	seen := make(map[string]int)

	for i, fc := range forecasts {
		name := refNameOf(fc)
		if cnt, ok := seen[name]; ok { // duplicate → make unique
			cnt++
			name = fmt.Sprintf("%s_%d", name, cnt)
			seen[refNameOf(fc)] = cnt
		} else {
			seen[name] = 0
		}
		rows[i] = name
	}
	return &gridMany{rows: rows, cols: cols, values: values}
}

func buildManyGridCorrected(forecasts []*WeatherForecast) *gridMany {
	if len(forecasts) == 0 {
		return nil
	}

	// collect unique target names
	colSet := make(map[string]bool)
	for _, fc := range forecasts {
		for _, ws := range fc.Systems {
			if ws != nil {
				colSet[ws.TargetName] = true
			}
		}
	}
	if len(colSet) == 0 {
		return nil
	}
	cols := make([]string, 0, len(colSet))
	for n := range colSet {
		cols = append(cols, n)
	}
	sort.Strings(cols)

	// index helpers
	colIdx := make(map[string]int, len(cols))
	for i, n := range cols {
		colIdx[n] = i
	}

	// build matrix initialised with NaN
	values := make([][]float64, len(cols))
	for i := range values {
		values[i] = make([]float64, len(forecasts))
		for j := range values[i] {
			values[i][j] = math.NaN()
		}
	}

	// fill matrix
	for r, fc := range forecasts {
		for _, ws := range fc.Systems {
			if ws == nil {
				continue
			}
			c := colIdx[ws.TargetName]
			values[c][r] = ws.CorrectedMeanMatchScore
		}
	}

	// rows are now reference names (duplicates get suffix _2, _3, …)
	rows := make([]string, len(forecasts))
	seen := make(map[string]int)

	for i, fc := range forecasts {
		name := refNameOf(fc)
		if cnt, ok := seen[name]; ok { // duplicate → make unique
			cnt++
			name = fmt.Sprintf("%s_%d", name, cnt)
			seen[refNameOf(fc)] = cnt
		} else {
			seen[name] = 0
		}
		rows[i] = name
	}
	return &gridMany{rows: rows, cols: cols, values: values}
}

func PlotTimeKmerizeVsRefLength(
	forecasts []*WeatherForecast,
	outPNG string,
	scattersize float64, plotxdim int, plotydim int, fontsize int,
	titleFontSize int, averageWindow int, lineWidth int, spacing float64) {

	boxWidth := 20.0

	// var pts plotter.XYs
	group := make(map[int][]float64)
	for _, fc := range forecasts {
		if fc == nil {
			continue
		}

		// find the self row = reference length
		refLen := -1
		for _, ws := range fc.Systems {
			if ws != nil && ws.RefName == ws.TargetName {
				refLen = ws.SequenceLength
				break
			}
		}
		if refLen < 0 {
			continue // skip forecast with no self-row
		}

		// add every non-self system as a point
		for _, ws := range fc.Systems {
			if ws == nil || ws.RefName == ws.TargetName {
				continue
			}
			y := ws.TimeKmerize
			group[refLen] = append(group[refLen], y)
			// pts = append(pts, plotter.XY{
			// 	X: float64(refLen),
			// 	Y: y,
			// })
		}
	}

	var lengths []int
	for l := range group {
		lengths = append(lengths, l)
	}
	sort.Ints(lengths)

	p := plot.New()
	p.Title.Text = "Time to k-merize vs reference length"
	p.Title.TextStyle.Font.Size = vg.Points(float64(titleFontSize))
	p.X.Label.Text, p.Y.Label.Text = "Reference length (bp)", "Time k-merize (s)"
	p.X.Label.TextStyle.Font.Size = vg.Points(float64(titleFontSize))
	p.Y.Label.TextStyle.Font.Size = vg.Points(float64(titleFontSize))
	p.X.Tick.Label.Font.Size = vg.Points(float64(fontsize))
	p.Y.Tick.Label.Font.Size = vg.Points(float64(fontsize))
	p.X.Tick.Label.Rotation = math.Pi / 2
	p.X.Tick.Label.XAlign = draw.XRight
	p.X.Tick.Label.YAlign = draw.YCenter

	/* colours */
	boxCol := color.RGBA{R: 0, G: 90, B: 0, A: 255}

	/* add one box per length */
	for idx, l := range lengths {
		values := make(plotter.Values, len(group[l]))
		for i, v := range group[l] {
			values[i] = v
		}
		box, _ := plotter.NewBoxPlot(vg.Points(boxWidth), float64(idx), values)

		box.FillColor = color.NRGBA{R: 0, G: 90, B: 0, A: 60}
		box.BoxStyle.Color = boxCol     // border of the box
		box.MedianStyle.Color = boxCol  // median line
		box.WhiskerStyle.Color = boxCol // whiskers
		// box.OutlierStyle.Color = boxCol
		p.Add(box)
	}

	/* x-axis: label the tick at each integer position with the length */
	if len(lengths) == 0 {
		fmt.Println("no reference-length groups – nothing to plot")
	}

	labels := make([]string, len(lengths))
	for i, l := range lengths {
		labels[i] = fmt.Sprintf("%d", l)
	}
	p.NominalX(labels...)

	if err := p.Save(1.5*vg.Length(plotxdim)*vg.Inch, 1.5*vg.Length(plotydim)*vg.Inch, outPNG); err != nil {
		fmt.Printf("Failed to save plot: %v", err)
	}

}

func PlotTimeStormVsRefLength(
	forecasts []*WeatherForecast,
	outPNG string,
	scattersize float64, plotxdim int, plotydim int, fontsize int,
	titleFontSize int, averageWindow int, lineWidth int, spacing float64) {

	boxWidth := 20.0

	// var pts plotter.XYs
	group := make(map[int][]float64)
	for _, fc := range forecasts {
		if fc == nil {
			continue
		}

		// find the self row = reference length
		refLen := -1
		for _, ws := range fc.Systems {
			if ws != nil && ws.RefName == ws.TargetName {
				refLen = ws.SequenceLength
				break
			}
		}
		if refLen < 0 {
			continue // skip forecast with no self-row
		}

		// add every non-self system as a point
		for _, ws := range fc.Systems {
			if ws == nil || ws.RefName == ws.TargetName {
				continue
			}
			y := ws.TimeStormMatch
			group[refLen] = append(group[refLen], y)
			// pts = append(pts, plotter.XY{
			// 	X: float64(refLen),
			// 	Y: y,
			// })
		}
	}

	var lengths []int
	for l := range group {
		lengths = append(lengths, l)
	}
	sort.Ints(lengths)

	p := plot.New()
	p.Title.Text = "Time to Storm Match vs reference length"
	p.Title.TextStyle.Font.Size = vg.Points(float64(titleFontSize))
	p.X.Label.Text, p.Y.Label.Text = "Reference length (bp)", "Time Storm Match (s)"
	p.X.Label.TextStyle.Font.Size = vg.Points(float64(titleFontSize))
	p.Y.Label.TextStyle.Font.Size = vg.Points(float64(titleFontSize))
	p.X.Tick.Label.Font.Size = vg.Points(float64(fontsize))
	p.Y.Tick.Label.Font.Size = vg.Points(float64(fontsize))
	p.X.Tick.Label.Rotation = math.Pi / 2
	p.X.Tick.Label.XAlign = draw.XRight
	p.X.Tick.Label.YAlign = draw.YCenter

	/* colours */
	boxCol := color.RGBA{R: 0, G: 90, B: 0, A: 255}

	/* add one box per length */
	for idx, l := range lengths {
		values := make(plotter.Values, len(group[l]))
		for i, v := range group[l] {
			values[i] = v
		}
		box, _ := plotter.NewBoxPlot(vg.Points(boxWidth), float64(idx), values)

		box.FillColor = color.NRGBA{R: 0, G: 90, B: 0, A: 60}
		box.BoxStyle.Color = boxCol     // border of the box
		box.MedianStyle.Color = boxCol  // median line
		box.WhiskerStyle.Color = boxCol // whiskers
		// box.OutlierStyle.Color = boxCol
		p.Add(box)
	}

	/* x-axis: label the tick at each integer position with the length */
	if len(lengths) == 0 {
		fmt.Println("no reference-length groups – nothing to plot")
	}

	labels := make([]string, len(lengths))
	for i, l := range lengths {
		labels[i] = fmt.Sprintf("%d", l)
	}
	p.NominalX(labels...)

	if err := p.Save(1.5*vg.Length(plotxdim)*vg.Inch, 1.5*vg.Length(plotydim)*vg.Inch, outPNG); err != nil {
		fmt.Printf("Failed to save plot: %v", err)
	}

}

func PlotAllForecastsHeatmap(
	forecasts []*WeatherForecast,
	outPNG string,
	scattersize float64, plotxdim int, plotydim int, fontsize int,
	titleFontSize int, averageWindow int, lineWidth int, spacing float64,
	maxMatchValue float64) {

	fmt.Println("1")
	grid := buildManyGrid(forecasts)
	if grid == nil || len(grid.rows) == 0 || len(grid.cols) == 0 {
		fmt.Println("nothing to plot – grid is empty")
		return
	}

	csvPath := strings.TrimSuffix(outPNG, filepath.Ext(outPNG)) + ".csv"

	file, err := os.Create(csvPath)
	if err != nil {
		fmt.Printf("cannot create CSV: %v\n", err)
	} else {
		w := csv.NewWriter(file)

		// header row: empty cell + column labels
		header := make([]string, 1+len(grid.cols))
		header[0] = "Reference \\ Target"
		copy(header[1:], grid.cols)
		_ = w.Write(header)

		// data rows
		for r, rowName := range grid.rows {
			record := make([]string, 1+len(grid.cols))
			record[0] = rowName
			for c := range grid.cols {
				v := grid.values[c][r]
				if math.IsNaN(v) {
					record[c+1] = "" // missing combination
				} else {
					record[c+1] = strconv.FormatFloat(v, 'f', 4, 64)
				}
			}
			_ = w.Write(record)
		}
		w.Flush()
		_ = file.Close()
		fmt.Println("CSV written to", csvPath)
	}

	fmt.Println("2")

	for c := range grid.values {
		for r, v := range grid.values[c] {
			if !math.IsNaN(v) && v > maxMatchValue {
				grid.values[c][r] = maxMatchValue
			}
		}
	}

	// colour map (same fallback as earlier)
	pal := grayScale{n: 256}
	hm := plotter.NewHeatMap(grid, pal)
	fmt.Println("3")
	titleSize := vg.Points(float64(titleFontSize))

	fmt.Println("4")

	p := plot.New()
	p.Title.Text = "Mean match score (ref rows × target cols)"
	p.X.Label.Text = "Target Sequence"
	p.Y.Label.Text = "Reference Sequence"
	p.Title.TextStyle.Font.Size = titleSize
	p.X.Label.TextStyle.Font.Size = titleSize
	p.Y.Label.TextStyle.Font.Size = titleSize

	p.X.Tick.Marker = plot.TickerFunc(func(min, max float64) []plot.Tick {
		var ticks []plot.Tick
		for c, name := range grid.cols {
			ticks = append(ticks, plot.Tick{Value: float64(c), Label: name})
		}
		return ticks
	})
	p.Y.Tick.Marker = plot.TickerFunc(func(min, max float64) []plot.Tick {
		var ticks []plot.Tick
		for r, name := range grid.rows {
			ticks = append(ticks, plot.Tick{Value: float64(r), Label: name})
		}
		return ticks
	})

	p.X.Tick.Label.Rotation = math.Pi / 2
	p.X.Tick.Label.XAlign = draw.XRight
	p.X.Tick.Label.YAlign = draw.YCenter
	p.X.Tick.Label.Font.Size = vg.Points(14) // column names
	p.Y.Tick.Label.Font.Size = vg.Points(14) // row names
	p.Add(hm)
	p.Add(plotter.NewGrid())

	if err := p.Save(1.5*vg.Length(plotxdim)*vg.Inch, 1.5*vg.Length(plotydim)*vg.Inch, outPNG); err != nil {
		fmt.Printf("Failed to save plot: %v", err)
	}
}

func PlotAllForecastsHitFrequencyHeatmap(
	forecasts []*WeatherForecast,
	outPNG string,
	scattersize float64, plotxdim int, plotydim int, fontsize int,
	titleFontSize int, averageWindow int, lineWidth int, spacing float64,
	maxMatchValue float64) {

	fmt.Println("1")
	grid := buildManyGridHitFrequency(forecasts)
	if grid == nil || len(grid.rows) == 0 || len(grid.cols) == 0 {
		fmt.Println("nothing to plot – grid is empty")
		return
	}

	csvPath := strings.TrimSuffix(outPNG, filepath.Ext(outPNG)) + ".csv"

	file, err := os.Create(csvPath)
	if err != nil {
		fmt.Printf("cannot create CSV: %v\n", err)
	} else {
		w := csv.NewWriter(file)

		// header row: empty cell + column labels
		header := make([]string, 1+len(grid.cols))
		header[0] = "Reference \\ Target"
		copy(header[1:], grid.cols)
		_ = w.Write(header)

		// data rows
		for r, rowName := range grid.rows {
			record := make([]string, 1+len(grid.cols))
			record[0] = rowName
			for c := range grid.cols {
				v := grid.values[c][r]
				if math.IsNaN(v) {
					record[c+1] = "" // missing combination
				} else {
					record[c+1] = strconv.FormatFloat(v, 'f', 4, 64)
				}
			}
			_ = w.Write(record)
		}
		w.Flush()
		_ = file.Close()
		fmt.Println("CSV written to", csvPath)
	}

	fmt.Println("2")

	for c := range grid.values {
		for r, v := range grid.values[c] {
			if !math.IsNaN(v) && v > maxMatchValue {
				grid.values[c][r] = maxMatchValue
			}
		}
	}

	// colour map (same fallback as earlier)
	pal := grayScale{n: 256}
	hm := plotter.NewHeatMap(grid, pal)
	fmt.Println("3")
	titleSize := vg.Points(float64(titleFontSize))

	fmt.Println("4")

	p := plot.New()
	p.Title.Text = "Hit Frequency (ref rows × target cols)"
	p.X.Label.Text = "Target Sequence"
	p.Y.Label.Text = "Reference Sequence"
	p.Title.TextStyle.Font.Size = titleSize
	p.X.Label.TextStyle.Font.Size = titleSize
	p.Y.Label.TextStyle.Font.Size = titleSize

	p.X.Tick.Marker = plot.TickerFunc(func(min, max float64) []plot.Tick {
		var ticks []plot.Tick
		for c, name := range grid.cols {
			ticks = append(ticks, plot.Tick{Value: float64(c), Label: name})
		}
		return ticks
	})
	p.Y.Tick.Marker = plot.TickerFunc(func(min, max float64) []plot.Tick {
		var ticks []plot.Tick
		for r, name := range grid.rows {
			ticks = append(ticks, plot.Tick{Value: float64(r), Label: name})
		}
		return ticks
	})

	p.X.Tick.Label.Rotation = math.Pi / 2
	p.X.Tick.Label.XAlign = draw.XRight
	p.X.Tick.Label.YAlign = draw.YCenter
	p.X.Tick.Label.Font.Size = vg.Points(14) // column names
	p.Y.Tick.Label.Font.Size = vg.Points(14) // row names
	p.Add(hm)
	p.Add(plotter.NewGrid())

	if err := p.Save(1.5*vg.Length(plotxdim)*vg.Inch, 1.5*vg.Length(plotydim)*vg.Inch, outPNG); err != nil {
		fmt.Printf("Failed to save plot: %v", err)
	}
}

func buildManyGridBed(forecasts []*WeatherForecast) *gridMany {
	if len(forecasts) == 0 {
		return nil
	}

	// --- collect all unique X positions ---------------------------
	posSet := make(map[int]struct{})
	for _, fc := range forecasts {
		for _, ws := range fc.Systems {
			if ws == nil {
				continue
			}
			for _, p := range ws.BedPlotPoints {
				x := int(p.X)
				posSet[x] = struct{}{}
			}
		}
	}
	if len(posSet) == 0 {
		return nil
	}

	// --- turn into sorted slice -----------------------------------
	var posList []int
	for x := range posSet {
		posList = append(posList, x)
	}
	sort.Ints(posList)

	// columns as strings (e.g. "123")
	cols := make([]string, len(posList))
	for i, x := range posList {
		cols[i] = strconv.Itoa(x)
	}

	// index helper
	colIdx := make(map[int]int, len(posList))
	for i, x := range posList {
		colIdx[x] = i
	}

	// --- allocate matrix ------------------------------------------
	values := make([][]float64, len(cols))
	for c := range values {
		values[c] = make([]float64, len(forecasts))
		for r := range values[c] {
			values[c][r] = 0 // mark "no data"
		}
	}

	// --- fill matrix ----------------------------------------------
	rows := make([]string, len(forecasts))
	for r, fc := range forecasts {
		rows[r] = fmt.Sprintf("forecast_%d", r+1)
		for _, ws := range fc.Systems {
			if ws == nil {
				continue
			}
			for _, p := range ws.BedPlotPoints {
				x := int(p.X)
				c := colIdx[x]
				values[c][r] = p.Y
			}
		}
	}
	return &gridMany{rows: rows, cols: cols, values: values}
}

func PlotAllForecastsRefCoverageHeatmap(
	forecasts []*WeatherForecast,
	outPNG string,
	scattersize float64, plotxdim int, plotydim int, fontsize int,
	titleFontSize int, averageWindow int, lineWidth int, spacing float64,
	maxMatchValue float64) {

	fmt.Println("PlotAllForecastsRefCoverageHeatmap")
	fmt.Println("1")
	grid := buildManyGridBed(forecasts)
	if grid == nil || len(grid.rows) == 0 || len(grid.cols) == 0 {
		fmt.Println("nothing to plot – grid is empty")
		return
	}

	csvPath := strings.TrimSuffix(outPNG, filepath.Ext(outPNG)) + ".csv"

	file, err := os.Create(csvPath)
	if err != nil {
		fmt.Printf("cannot create CSV: %v\n", err)
	} else {
		w := csv.NewWriter(file)

		// header row: empty cell + column labels
		header := make([]string, 1+len(grid.cols))
		header[0] = "Target \\ Reference Storm Position"
		copy(header[1:], grid.cols)
		_ = w.Write(header)

		// data rows
		for r, rowName := range grid.rows {
			record := make([]string, 1+len(grid.cols))
			record[0] = rowName
			for c := range grid.cols {
				v := grid.values[c][r]
				if math.IsNaN(v) {
					record[c+1] = "" // missing combination
				} else {
					record[c+1] = strconv.FormatFloat(v, 'f', 4, 64)
				}
			}
			_ = w.Write(record)
		}
		w.Flush()
		_ = file.Close()
		fmt.Println("CSV written to", csvPath)
	}

	fmt.Println("2")

	pal := grayScale{n: 256}
	hm := plotter.NewHeatMap(grid, pal)
	fmt.Println("3")
	titleSize := vg.Points(float64(titleFontSize))

	fmt.Println("4")

	p := plot.New()
	p.Title.Text = "Reference Coverage"
	p.X.Label.Text = "Target Sequence"
	p.Y.Label.Text = "Reference Storm Positions"
	p.Title.TextStyle.Font.Size = titleSize
	p.X.Label.TextStyle.Font.Size = titleSize
	p.Y.Label.TextStyle.Font.Size = titleSize

	p.X.Tick.Marker = plot.TickerFunc(func(min, max float64) []plot.Tick {
		var ticks []plot.Tick
		for c, name := range grid.cols {
			ticks = append(ticks, plot.Tick{Value: float64(c), Label: name})
		}
		return ticks
	})
	p.Y.Tick.Marker = plot.TickerFunc(func(min, max float64) []plot.Tick {
		var ticks []plot.Tick
		for r, name := range grid.rows {
			ticks = append(ticks, plot.Tick{Value: float64(r), Label: name})
		}
		return ticks
	})

	p.X.Tick.Label.Rotation = math.Pi / 2
	p.X.Tick.Label.XAlign = draw.XRight
	p.X.Tick.Label.YAlign = draw.YCenter
	p.X.Tick.Label.Font.Size = vg.Points(14) // column names
	p.Y.Tick.Label.Font.Size = vg.Points(14) // row names
	p.Add(hm)
	p.Add(plotter.NewGrid())

	if err := p.Save(1.5*vg.Length(plotxdim)*vg.Inch, 1.5*vg.Length(plotydim)*vg.Inch, outPNG); err != nil {
		fmt.Printf("Failed to save plot: %v", err)
	}
}

func PlotAllForecastsCorrectedHeatmap(
	forecasts []*WeatherForecast,
	outPNG string,
	scattersize float64, plotxdim int, plotydim int, fontsize int,
	titleFontSize int, averageWindow int, lineWidth int, spacing float64,
	maxMatchValue float64) {

	fmt.Println("1")
	grid := buildManyGridCorrected(forecasts)
	if grid == nil || len(grid.rows) == 0 || len(grid.cols) == 0 {
		fmt.Println("nothing to plot – grid is empty")
		return
	}

	csvPath := strings.TrimSuffix(outPNG, filepath.Ext(outPNG)) + ".csv"

	file, err := os.Create(csvPath)
	if err != nil {
		fmt.Printf("cannot create CSV: %v\n", err)
	} else {
		w := csv.NewWriter(file)

		// header row: empty cell + column labels
		header := make([]string, 1+len(grid.cols))
		header[0] = "Reference \\ Target"
		copy(header[1:], grid.cols)
		_ = w.Write(header)

		// data rows
		for r, rowName := range grid.rows {
			record := make([]string, 1+len(grid.cols))
			record[0] = rowName
			for c := range grid.cols {
				v := grid.values[c][r]
				if math.IsNaN(v) {
					record[c+1] = "" // missing combination
				} else {
					record[c+1] = strconv.FormatFloat(v, 'f', 4, 64)
				}
			}
			_ = w.Write(record)
		}
		w.Flush()
		_ = file.Close()
		fmt.Println("CSV written to", csvPath)
	}

	fmt.Println("2")

	for c := range grid.values {
		for r, v := range grid.values[c] {
			if !math.IsNaN(v) && v > maxMatchValue {
				grid.values[c][r] = maxMatchValue
			}
		}
	}

	// colour map (same fallback as earlier)
	pal := grayScale{n: 256}
	hm := plotter.NewHeatMap(grid, pal)
	fmt.Println("3")
	titleSize := vg.Points(float64(titleFontSize))

	fmt.Println("4")

	p := plot.New()
	p.Title.Text = "Mean match score (ref rows × target cols)"
	p.X.Label.Text = "Target Sequence"
	p.Y.Label.Text = "Reference Sequence"
	p.Title.TextStyle.Font.Size = titleSize
	p.X.Label.TextStyle.Font.Size = titleSize
	p.Y.Label.TextStyle.Font.Size = titleSize

	p.X.Tick.Marker = plot.TickerFunc(func(min, max float64) []plot.Tick {
		var ticks []plot.Tick
		for c, name := range grid.cols {
			ticks = append(ticks, plot.Tick{Value: float64(c), Label: name})
		}
		return ticks
	})
	p.Y.Tick.Marker = plot.TickerFunc(func(min, max float64) []plot.Tick {
		var ticks []plot.Tick
		for r, name := range grid.rows {
			ticks = append(ticks, plot.Tick{Value: float64(r), Label: name})
		}
		return ticks
	})

	p.X.Tick.Label.Rotation = math.Pi / 2
	p.X.Tick.Label.XAlign = draw.XRight
	p.X.Tick.Label.YAlign = draw.YCenter
	p.X.Tick.Label.Font.Size = vg.Points(14) // column names
	p.Y.Tick.Label.Font.Size = vg.Points(14) // row names
	p.Add(hm)
	p.Add(plotter.NewGrid())

	if err := p.Save(1.5*vg.Length(plotxdim)*vg.Inch, 1.5*vg.Length(plotydim)*vg.Inch, outPNG); err != nil {
		fmt.Printf("Failed to save plot: %v", err)
	}
}

// Relative contains a set of stormologs versus a particular genome (set of fragments)
type Relatives struct {
	count  int
	pieces []*stormologs
}

// Newrel creates the relative and the pieces list returns the relative
func Newrel() *Relatives {
	r := new(Relatives)
	r.pieces = make([]*stormologs, 0)
	return r
}

// Addrel adds the set of stormolog hits to the relative
func (r *Relatives) Addrel(setofhits *stormologs) {
	r.pieces = append(r.pieces, setofhits)
}

// Print prints out the identified storm homologs
func (r *Relatives) Print(rels_outfile string, revigos *reverseoligos, goodmin int) {
	fmt.Println("Printing relatives to ", rels_outfile)
	fkout, _ := os.Create(rels_outfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output

	fmt.Println("\n\t\tPrinting relatives\n", len(r.pieces), goodmin)
	//	fmt.Fprintf(kwriter, "%s\t%s\t", "Name", "Length")
	fmt.Printf("\n\t\tPrinting rels %d\n", goodmin)
	fmt.Fprintf(kwriter, "There were %d relatives identified\n", len(r.pieces))
	// fmt.Fprintf(kwriter, "refloc\ttargetloc\tmatchcount\n")
	//	goodhits   []*stormhit // map hits to given storm ID

	//type stormhit struct {
	//	matches    []int // list of match positions
	//	matchcount int
	//	frag       *fragment
	//	storm      *storm
	//	boundary   int
	//	location   int // first match position locates the hit

	piececount := 0
	for _, piece := range r.pieces {
		piececount++
		hitcount := 0
		matchtotal := 0
		for _, gh := range piece.goodhits {
			hitcount++
			matchtotal += gh.matchcount
		}
		piece.matchtotal = matchtotal
		piece.hitcount = hitcount
	}

	fmt.Fprintf(kwriter, "ID\thits\thits2\ttotal\tmean\tmin\n")
	for i, piece := range r.pieces {
		hits := len(piece.goodhits)
		min := piece.goodhitmin
		total := piece.matchtotal
		mean := float64(total) / float64(piece.hitcount)
		fmt.Fprintf(kwriter, "%d\t%s\t%d\t%d\t%d\t%f\t%d\n", i, piece.Seqname, hits, piece.hitcount, total, mean, min)
		//		fmt.Fprintf(kwriter, "ID %d had %d hits to storms (hits %d) with min %d\n", i, len(piece.goodhits), piece.hitcount, piece.goodhitmin)
	}

	fmt.Fprintf(kwriter, "\n")
}

// func (r *Relatives) PlotRelatives(rels_plotfile string, revigos *reverseoligos, goodmin int) {

// 	for i, piece := range r.pieces {
// 		hits := len(piece.goodhits)
// 		min := piece.goodhitmin
// 		total := piece.matchtotal
// 		mean := float64(total) / float64(piece.hitcount)
// 		fmt.Fprintf(kwriter, "ID\thits\thits2\tmin\n")
// 		fmt.Fprintf(kwriter, "%d\t%d\t%d\t%lf\t%d\n", i, hits, piece.hitcount, total, mean, min)
// 		//		fmt.Fprintf(kwriter, "ID %d had %d hits to storms (hits %d) with min %d\n", i, len(piece.goodhits), piece.hitcount, piece.goodhitmin)
// 	}

// 	// Create a new plot
// 	p := plot.New()

// 	p.Title.Text = "Kmer Positions: " + querygenome + " vs " + testgenome
// 	p.X.Label.Text = querygenome + " Sequence Position"
// 	p.Y.Label.Text = testgenome + " Sequence Position"
// 	p.Add(plotter.NewGrid())

// 	// scatter plot, values need to be float64 type
// 	points := make(plotter.XYs, len(xSequencePositions))
// 	for i := range points {
// 		points[i].X = float64(xSequencePositions[i])
// 		points[i].Y = float64(yGenomePositions[i])
// 	}

// 	scatter, err := plotter.NewScatter(points)
// 	if err != nil {
// 		fmt.Printf("Failed to create scatter plot: %v", err)
// 	}
// 	scatter.GlyphStyle.Color = plotter.DefaultLineStyle.Color
// 	scatter.GlyphStyle.Radius = vg.Points(scattersize)

// 	p.Add(scatter)

// 	// Save the plot to a PNG file
// 	os.MkdirAll(kmerPlotDir, os.ModePerm)
// 	kmerPlotFileName = kmerPlotDir + kmerPlotFileName
// 	if err := p.Save(vg.Length(plotxdim)*vg.Inch, vg.Length(plotydim)*vg.Inch, kmerPlotFileName); err != nil {
// 		fmt.Printf("Failed to save plot: %v", err)
// 	}

// 	// Save the x and y values to a CSV file
// 	csvFileName := kmerPlotFileName[:strings.LastIndex(kmerPlotFileName, ".")] + "_data.csv"
// 	csvFile, err := os.Create(csvFileName)
// 	if err != nil {
// 		fmt.Printf("Failed to create CSV file: %v", err)
// 	}
// 	defer csvFile.Close()

// 	writer := csv.NewWriter(csvFile)
// 	defer writer.Flush()

// 	// Write CSV headers
// 	if err := writer.Write([]string{querygenome + " Sequence Position", testgenome + " Sequence Position"}); err != nil {
// 		fmt.Printf("Failed to write CSV headers: %v", err)
// 	}

// 	// Write CSV data
// 	for i := range xSequencePositions {
// 		record := []string{
// 			strconv.Itoa(xSequencePositions[i]),
// 			strconv.Itoa(yGenomePositions[i]),
// 		}
// 		if err := writer.Write(record); err != nil {
// 			fmt.Printf("Failed to write CSV record: %v", err)
// 		}
// 	}

// }

// goodhit.Print(kwriter, revigos, i)
// fmt.Fprintf(kwriter, "%d\t%d\n", stormID, sh.matchcount)

// SeedStorms uses kloc list of kmers.frags to build seed storms from
// non-overlapping kmers as seeds for the storm PClouds
// Each storm contains stormsize PClouds. SeedStorms expects empty (0 count) but not nil storm pointers
func (kmers *Oligos) SeedStorms(stormsize int, storms *Tempest) {
	fmt.Println("Beginning SeedStorms")
	if kmers != nil && storms != nil && storms.count == 0 {
		storms.SeedName = kmers.name
		// storms.Init(stormsize) // this could have been done when storms was created
		klen := kmers.klen
		klen2 := klen + klen
		// fmt.Println("SeedStorms initiating new storms ", klen, storms.stormsize)
		var storm *storm
		var seedrefpos int
		for _, frag := range kmers.frags { // for each fragment in kmers
			// go through each fragment of the original sequence file in non-overlapping chunks of klen
			// fmt.Println("SeedStorms going through kmers in fragment ", frag.name, klen, storms.stormsize)
			for pos, cloudloc := 0, 0; (pos + klen) < frag.total; pos += kmers.klen {
				// if pos < 1000 {
				// 	fmt.Println("seedstorms in", pos, klen, cloudloc)
				// }
				if cloudloc == 0 { // will rezero every time a storm is completed
					storm = storms.newstorm()
				}
				// for each klen section (a kmer)
				storm.newcloud(frag.klocs[pos], cloudloc)
				if cloudloc == storms.stormseedloc {
					seedrefpos = pos
				}
				cloudloc++
				if cloudloc >= stormsize || (pos+klen2) >= frag.total {
					storm.newstormfinish(cloudloc-1, storms.stormseedloc, pos, seedrefpos, kmers.name, frag.name)
					seedrefpos = 0
					cloudloc = 0
				}
			}
			// fmt.Println("SeedStorms finished frag,  ", frag.name, storms.count, storms.stormseedloc, storms.stormsize)
		}
	} else {
		fmt.Println("SeedStorms didn't do anything (storms, stormsize): ", stormsize)
	}
}

// FindHits uses kmer position information to find matches to clouds
func (kmers *Oligos) FindHits(storms *Tempest, boundary int, extendtoside int, goodmin int, stormSize int) *stormologs {
	sh := newstormologs(goodmin, kmers.name)
	stormparse := 0
	for _, fragment := range kmers.frags { // fragments in the test sequences
		// For every fragment iterate over every reference storm
		// fmt.Println("---------- length fragment ", fragment.total)
		for _, storm := range storms.kstorms { // for each storm in list of storms in tempest
			// fmt.Println("Storm ref  location: ", storm.seedrefpos)
			// fmt.Println("Storm test location: ", storm.seedpos)
			// fmt.Println("Storm size: ", storm.size)
			stormparse++
			// for the current fragment / storm pair
			hitlst := newstormologs(goodmin, kmers.name)
			locinstorm := storm.seedpos

			// if stormid < 5 {
			// 	fmt.Println("first seed", locinstorm, stormid, storm.seedpos)
			// }
			// add locations of PCloud members to hitlist, starting from seed and extend outwards
			hitlst.addkmerlocs(storm.PClouds, locinstorm, storm.size, fragment.kmap, boundary, extendtoside, stormparse)

			for i := 1; i <= (storm.size / 2); i++ { // make sure this works under all conditions to extend out the storm
				// fmt.Println("i ", i)
				locinstorm = storm.seedpos + i
				// if stormid < 5 && (extendtoside-i) > -2 {
				// 	// fmt.Println("first i", i, locinstorm, stormid, extendtoside, boundary, stormparse)
				// }
				hitlst.addkmerlocs(storm.PClouds, locinstorm, storm.size, fragment.kmap, boundary, extendtoside-i, stormparse)
				locinstorm = storm.seedpos - i
				hitlst.addkmerlocs(storm.PClouds, locinstorm, storm.size, fragment.kmap, boundary, extendtoside-i, stormparse)
			}
			sh.addgood_stormologs(hitlst, fragment, storm, boundary, goodmin, stormSize)
			// fmt.Println("				Count ", count)
		}
	}
	return sh
}

// addkmerlocs adds fragment locations of kmers in a PCloud to a hit list
// iff a particular PCloud kmer position is close enough to an existing hit list
// else generates a new hit and adds pos as location
// expects that kmapr uses `*kmer`
// extendstate is used to limit how many pclouds get added as possible new hits
// called as hitlist.addkmerlocs(storm.PClouds[storm.seedpos-i],fragment.kmap, boundary)
func (hitlst *stormologs) addkmerlocs(clouds []*PCloud, loc int, size int, kmapr map[*kneighbors]*oligo, boundary int, extendstate int, count int) {

	if count < 0 { // turned off
		fmt.Println("## addkmerlocs ", boundary, extendstate, count)
		fmt.Println("# storm kmer locs add ", count, extendstate)
	}
	if loc >= 0 && loc < size { // make sure that loc is in range
		cloud := clouds[loc]
		for _, kptr := range cloud.cloudmers {
			if kmapr[kptr] == nil {
				// fmt.Println("this would throw a nil pointer dereference error if allowed to pass")
			} else {
				poses := kmapr[kptr].poses  // positions of given kmer in the given sequence fragment
				for _, pos := range poses { // if near another pos, add to that hit, otherwise create new hit
					matched := false
					for _, oldhit := range hitlst.goodhits {
						match_distance := abs(oldhit.location - pos)
						if match_distance < boundary { // near a previous match
							oldhit.addpos(pos)
							matched = true
						}
					}
					if matched == false && extendstate >= 0 { // add new hit to hitlst and set location as well as first pos
						hitlst.addnewhit(pos)
					}
				}
			}
		}
	}
}

// a hit is a potentially homologous site, and a match is a kmer match between storm and fragment

// addgoodhits adds a hitlist to a set of hits for all storms
// called as stormhits.addhitlist(hitlst, fragment, storm, boundary)
func (tempesthits *stormologs) addgood_stormologs(hitlst *stormologs, frag *fragment, storm *storm, boundary int, minhit int, stormSize int) {

	// ifStormHit := false
	for _, hit := range hitlst.goodhits {
		if hit.matchcount >= minhit {
			hit.frag = frag
			hit.storm = storm
			hit.boundary = boundary
			if hit.matchcount > stormSize {
				hit.matchcount = stormSize
			}
			tempesthits.goodhits = append(tempesthits.goodhits, hit)
			// ifStormHit = true
		}
	}
	tempesthits.goodhitmin = minhit
	// if ifStormHit {
	// 	tempesthits.uniquehitcount++
	// }
}

func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

// Getoutfile outputs kmer counts and positions
func (kmers *Oligos) Getoutfile(infile string) {
	base := strings.SplitN(infile, ".", 2)
	strlen := strconv.Itoa(kmers.klen)
	kmers.Outfile = kmers.Outfile + strlen + "_" + base[0] + ".xls"
}

// GetoutfileMulti outputs kmer counts and positions
func (kmers *Oligos) GetoutfileMulti(kcountpre string, basename string, directory string) {
	strlen := strconv.Itoa(kmers.klen)
	kmers.Outfile = directory + kcountpre + strlen + "_" + basename + ".xls"
}

// Kprint outputs kmer counts
func (kmers *Oligos) Kprint() {
	// fmt.Println("Opening Kmer Count Output File", kmers.Outfile)

	// Get the directory part of the filePath
	dir := filepath.Dir(kmers.Outfile)

	// Create the directory structure
	if err := os.MkdirAll(dir, os.ModePerm); err != nil {
		fmt.Printf("error creating directories: %v", err)
	}

	fkout, _ := os.Create(kmers.Outfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output

	fmt.Fprintln(kwriter, "kmer\tcount")
	// fmt.Println("size of kmers.kcount", len(kmers.kcount))
	// fmt.Println("minprint, printNs", kmers.minprint, kmers.printNs)
	for kmer, kcount := range kmers.kcount {
		if (kcount >= kmers.minprint) && (kmers.printNs || (!strings.Contains(kmer, "N"))) {
			fmt.Fprintf(kwriter, "%s\t%d\n", kmer, kcount)
		}
	}
}

// Kposprint outputs kmer positions
func (kmers *Oligos) Kposprint(kmin int, kmax int) {
	// fmt.Println("Opening Kmer Position Output File", kmers.POutfile)

	// Get the directory part of the filePath
	dir := filepath.Dir(kmers.POutfile)

	// Create the directory structure
	if err := os.MkdirAll(dir, os.ModePerm); err != nil {
		fmt.Printf("error creating directories: %v", err)
	}

	fkout, _ := os.Create(kmers.POutfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output

	fmt.Fprintln(kwriter, "kmer\tcount\tpositions")
	for kmer, kcount := range kmers.kcount { // fix this so works with structure
		if (kcount >= kmin) && (kmers.printNs || (!strings.Contains(kmer, "N"))) {
			if kmers.kmap[kmer] != nil {
				poses := kmers.kmap[kmer].poses
				if (poses != nil) && (kcount < kmax) && (kcount >= kmin) {
					fmt.Fprintf(kwriter, "%s\t%d", kmer, kcount)
					for i := range poses {
						fmt.Fprintf(kwriter, "\t%d", poses[i])
					}
					fmt.Fprintf(kwriter, "\n")
				}
			}
		}
	}
}

// kadd adds 1 to the kinfo for given token
// see if this works, but maybe add simple map to kmers, then create kinfo if count is 1
func (kmers *Oligos) kadd(kmer string, kpos int) {
	var kinfo *oligo
	if _, ok := kmers.kmap[kmer]; !ok {
		kmers.kmap[kmer] = new(oligo)
	}
	kinfo = kmers.kmap[kmer] // have to reassign because out of loop
	kinfo.kcount += 1
}

func (kmers *Oligos) GetTotal() int {

	return kmers.total
}

func (kmers *Oligos) Getkmap() map[string]*oligo {

	return kmers.kmap
}

// Countref counts kmers in string, adds to stored counts in kmers
func (kmers *Oligos) Countref(seqs *Sequences, seq string, name string) {
	var kmer string
	for i := 0; i < (len(seq) - kmers.klen + 1); i++ {
		kmer = seq[i : i+kmers.klen]
		if kmers.kmap[kmer] == nil {
			kmers.kmap[kmer] = new(oligo)
			kmers.kmap[kmer].Init(kmer, kmers.kcount[kmer])
		}
		kmers.kcount[kmer]++
		kmers.total++
	}
	// create remnant in case it needs to be pre-pendend to next sequence fragment
	nextpos := len(seq) - kmers.klen + 1
	kmers.remnant = seq[nextpos : nextpos+kmers.klen-1]
}

// getmersInFragment counts kmers in sequence fragment string, adds to stored counts in kmers
// positions of kmers are stored in kinfo.poses, or kinfo=kmers.kmap[kmer]
func (kmers *Oligos) getmersInFragment(seqs *Sequences, seq string, name string, linepos int) {
	var kmer string
	// var firsttimefrag int
	// if kmers.frags == nil {
	// 	firsttimefrag = 1
	// }
	fragment := kmers.addfragment(name) // create and initialize if not already there
	for i := 0; i < (len(seq) - kmers.klen + 1); i++ {
		kmer = seq[i : i+kmers.klen]
		if kmers.alligos[kmer] == nil {
			kmers.alligos[kmer] = new(kneighbors)
		}
		fragment.addpose(kmers.alligos[kmer], linepos+i)
		// if firsttimefrag == 1 { // this is total hack for backwards compatibility
		// 	kmers.kcount[kmer]++
		// 	kmers.total++
		// 	kmers.kmap[kmer].poses = append(kmers.kmap[kmer].poses, (linepos + i))
		// } // note that poses will be incorrect if there is more than 1 fragment per kmers structure
	}

	// create remnant in case it needs to be pre-pendend to next sequence fragment
	nextpos := len(seq) - kmers.klen + 1
	if nextpos < 0 {
		// kmer is longer than the length of the line
		kmers.remnant = kmers.remnant + seq
	} else {
		kmers.remnant = seq[nextpos : nextpos+kmers.klen-1]
	}
}

// addmer adds kmer to kmers count
func (kmers *Oligos) addmer(kmer string) {
	if kmers.kmap[kmer] == nil {
		kmers.kmap[kmer] = new(oligo)
		kmers.kmap[kmer].Init(kmer, kmers.kcount[kmer])
	}
	kmers.kcount[kmer]++
	kmers.total++
}

// combine the identified kmer and counts from a second kmer set (kmer_add) with an original set (kmers_ori)
// TODO: how to combine concurrently to prevent race conditions?
func (kmers_ori *Oligos) Combiner(kmer_add *Oligos) {

	for kmer_key, kmer_entry := range kmer_add.kmap {

		if previous_entry, exists := kmers_ori.kmap[kmer_key]; exists {
			// kmer has already been found, merge data
			previous_entry.poses = append(previous_entry.poses, kmer_entry.poses...)
			previous_entry.kcount += previous_entry.kcount
		} else {
			// new kmer
			kmers_ori.kmap[kmer_key] = new(oligo)
			kmers_ori.kmap[kmer_key].name = kmer_key
			kmers_ori.kmap[kmer_key].poses = append(kmers_ori.kmap[kmer_key].poses, kmer_entry.poses...)
			kmers_ori.kmap[kmer_key].kcount = kmer_entry.kcount
		}
	}
}

// CopyStruct creates a copy of the given struct and returns a new pointer
func CopyOligos(original *Oligos) *Oligos {
	return &Oligos{
		name:       original.name,
		kmap:       original.kmap,
		rmap:       original.rmap,
		kcount:     original.kcount,
		locprimer:  original.locprimer,
		readcounts: original.readcounts,
		offmin:     original.offmin,
		total:      original.total,
		klen:       original.klen,
		minprint:   original.minprint,
		Outfile:    original.Outfile,
		POutfile:   original.POutfile,
		kfile:      original.kfile,
		printNs:    original.printNs,
		remnant:    original.remnant,
	}
}

//
// large database sequence readers and kmerizers //
//

func checkfilter(seqfilter map[string]*filter, name string) bool {
	var passfilter bool
	if seqfilter[name] != nil {
		passfilter = true
	} else {
		passfilter = false
	}
	return passfilter
}

func splitOnFileExtension(fileName string) (string, string) {
	// Find the index of the last dot
	extensionIndex := strings.LastIndex(fileName, ".")

	// If there's no period in the file name, return the entire file name as the base
	// and an empty string for the extension
	if extensionIndex == -1 {
		return fileName, ""
	}

	// Split the file name into base and extension
	baseName := fileName[:extensionIndex]
	extension := fileName[extensionIndex+1:]

	return baseName, extension
}

// Function to remove a substring from the beginning of a string
func removePrefix(str, prefix string) string {
	if strings.HasPrefix(str, prefix) {
		return str[len(prefix):]
	}
	return str
}

// function to take read in the parent sequence file and split into a set of new files of length ksplit.
// Output into a new sequence file name altered from the
// parent file name with the position of the child sequence within the parent sequence.
func (seqs *Sequences) Split(parentBaseDir string, ksplit int, splitSequenceDir string) map[string]*Sequences {

	// store new (split) sequences for kmerize
	seqSets := make(map[string]*Sequences)

	var name, seq string
	var count, lcount int
	lcount = 0
	count = 0

	// reading stuff and setup
	fmt.Println("File to open is ", seqs.seqfile)
	fpin, err := os.Open(seqs.seqfile)
	globals.Check(err)
	defer fpin.Close()
	scanner := bufio.NewScanner(fpin)
	entrycount := 0 // track line in each entry
	entrylimit := 0 // default every line in entry is counted
	parentLocation := 0
	// hackish fastq, just reading the first line after name line
	if seqs.filetype == "fastq" {
		entrylimit = 1
	}

	passfilter := false // flag to see if name is in filter list

	parentBaseName, parentExtension := splitOnFileExtension(seqs.Name)
	baseChildName := splitSequenceDir + removePrefix(parentBaseName, parentBaseDir)

	remnant := ""
	lengthParentLine := 0
	ifDetermineParentLine := true
	// read, record, count kmers
	for scanner.Scan() {
		lcount += 1
		line := scanner.Text()              // should not include eol
		trimline := strings.TrimSpace(line) // trim off leading and lagging whitespace
		nofilter := !seqs.dofilter
		if strings.HasPrefix(line, seqs.entrystart) {
			name = strings.TrimPrefix(trimline, seqs.entrystart)
			count += 1
			entrycount = 1
			remnant = ""

			if seqs.dofilter {
				passfilter = checkfilter(seqs.seqfilter, name)
			}
		} else {
			if (nofilter || passfilter) && (lcount < seqs.linelimit) && (entrylimit < 1 || entrycount <= entrylimit) {
				seq = remnant + trimline
				entrycount++
				if ifDetermineParentLine {
					lengthParentLine = len(seq)
					ifDetermineParentLine = false
				}
				if lcount <= seqs.linelimit {

					if len(seq) > ksplit {
						// cut off a new child sequence
						split_seq := seq[:ksplit]
						seq = seq[ksplit:]

						childSeqFileName := baseChildName + "_" + strconv.Itoa(ksplit) + "_" + strconv.Itoa(parentLocation) + "." + parentExtension

						split_seq_name := childSeqFileName

						// save split_seq sequence to file
						first_line_child := seqs.entrystart + name + " split into lengths of " + strconv.Itoa(ksplit) + " from position " + strconv.Itoa(parentLocation)

						var seqChunks []string
						seqChunks = append(seqChunks, first_line_child)
						// Loop through the string in steps of size lengthParentLine
						for i := 0; i < len(split_seq); i += lengthParentLine {
							// Calculate the end index for the current chunk
							end := i + lengthParentLine

							// Ensure that the end index does not exceed the string length
							if end > len(split_seq) {
								end = len(split_seq)
							}

							// Append the substring to the chunks slice
							seqChunks = append(seqChunks, split_seq[i:end])
						}

						// Create or open the file for writing
						if err := writeSequenceToFile(seqChunks, childSeqFileName); err != nil {
							fmt.Println("Error writing to file:", err)
						}

						// create new seqs and refmers, Init() them
						seqs_split := new(Sequences)
						seqs_split.Init(split_seq_name, seqs.minlength, seqs.linelimit, seqs.linemin, seqs.record, seqs.filetype, seqs.dofilter)
						seqSets[split_seq_name] = seqs_split
						parentLocation += ksplit
					}
					// else add a new line

				}
				remnant = seq
			}
		}
		// if at end of file, discard the final sequence if len(seq) < ksplit to keep consistent split chunks
	}
	fmt.Println("Seqs and Lines counted\n", count, lcount)

	return seqSets
}

// // function to take read in the parent sequence file and split into a set of new files of length ksplit.
// // Output into a new sequence file name altered from the
// // parent file name with the position of the child sequence within the parent sequence.
// func (seqs *Sequences) ReadSplit(parentBaseDir string, readSplitLength int, randomWidthInt int, randomSeedInt int, splitSequenceDir string) map[string]*Sequences {

// 	// store new (split) sequences for kmerize
// 	seqSets := make(map[string]*Sequences)

// 	var name, seq string
// 	var count, lcount int
// 	var split_position int
// 	lcount = 0
// 	count = 0

// 	rng := rand.New(rand.NewSource(int64(randomSeedInt)))
// 	split_position = rng.Intn(randomWidthInt)

// 	// reading stuff and setup
// 	fmt.Println("File to open is ", seqs.seqfile)
// 	fpin, err := os.Open(seqs.seqfile)
// 	globals.Check(err)
// 	defer fpin.Close()
// 	scanner := bufio.NewScanner(fpin)
// 	entrycount := 0 // track line in each entry
// 	entrylimit := 0 // default every line in entry is counted
// 	parentLocation := 0
// 	// hackish fastq, just reading the first line after name line
// 	if seqs.filetype == "fastq" {
// 		entrylimit = 1
// 	}

// 	passfilter := false // flag to see if name is in filter list

// 	parentBaseName, parentExtension := splitOnFileExtension(seqs.Name)
// 	baseChildName := splitSequenceDir + removePrefix(parentBaseName, parentBaseDir)

// 	remnant := ""
// 	lengthParentLine := 0
// 	ifDetermineParentLine := true
// 	// read, record, count kmers
// 	for scanner.Scan() {
// 		lcount += 1
// 		line := scanner.Text()              // should not include eol
// 		trimline := strings.TrimSpace(line) // trim off leading and lagging whitespace
// 		nofilter := !seqs.dofilter
// 		if strings.HasPrefix(line, seqs.entrystart) {
// 			name = strings.TrimPrefix(trimline, seqs.entrystart)
// 			count += 1
// 			entrycount = 1
// 			remnant = ""

// 			if seqs.dofilter {
// 				passfilter = checkfilter(seqs.seqfilter, name)
// 			}
// 		} else {
// 			if (nofilter || passfilter) && (lcount < seqs.linelimit) && (entrylimit < 1 || entrycount <= entrylimit) {
// 				seq = remnant + trimline
// 				entrycount++
// 				if ifDetermineParentLine {
// 					lengthParentLine = len(seq)
// 					ifDetermineParentLine = false
// 				}
// 				if lcount <= seqs.linelimit {

// 					if len(seq) > split_position+readSplitLength {

// 						fmt.Println("New read start position: " + strconv.Itoa(parentLocation+split_position) + " with parent start pos " + strconv.Itoa(parentLocation) + " and random start " + strconv.Itoa(split_position))
// 						// cut off a new child sequence
// 						split_seq := seq[split_position : split_position+readSplitLength]
// 						seq = seq[split_position+readSplitLength:]

// 						childSeqFileName := baseChildName + "_read" + "_" + strconv.Itoa(readSplitLength) + "_" + strconv.Itoa(parentLocation) + "." + parentExtension

// 						split_seq_name := childSeqFileName

// 						// save split_seq sequence to file
// 						first_line_child := seqs.entrystart + name + " split into lengths of " + strconv.Itoa(readSplitLength) + " from position " + strconv.Itoa(parentLocation+split_position)

// 						var seqChunks []string
// 						seqChunks = append(seqChunks, first_line_child)
// 						// Loop through the string in steps of size lengthParentLine
// 						for i := 0; i < len(split_seq); i += lengthParentLine {
// 							// Calculate the end index for the current chunk
// 							end := i + lengthParentLine

// 							// Ensure that the end index does not exceed the string length
// 							if end > len(split_seq) {
// 								end = len(split_seq)
// 							}

// 							// Append the substring to the chunks slice
// 							seqChunks = append(seqChunks, split_seq[i:end])
// 						}

// 						// Create or open the file for writing
// 						if err := writeSequenceToFile(seqChunks, childSeqFileName); err != nil {
// 							fmt.Println("Error writing to file:", err)
// 						}

// 						// create new seqs and refmers, Init() them
// 						seqs_split := new(Sequences)
// 						seqs_split.Init(split_seq_name, seqs.minlength, seqs.linelimit, seqs.linemin, seqs.record, seqs.filetype, seqs.dofilter)
// 						seqSets[split_seq_name] = seqs_split
// 						parentLocation += split_position + readSplitLength
// 						split_position = rng.Intn(randomWidthInt)
// 					}
// 					// else add a new line

// 				}
// 				remnant = seq
// 			}
// 		}
// 		// if at end of file, discard the final sequence if len(seq) < readSplitLength to keep consistent split chunks
// 	}
// 	fmt.Println("Seqs and Lines counted\n", count, lcount)

// 	return seqSets
// }

// Function to write a slice of strings to a file
func writeSequenceToFile(lines []string, fileName string) error {
	// Create or open the file for writing
	// Get the directory part of the filePath
	dir := filepath.Dir(fileName)

	// Create the directory structure
	if err := os.MkdirAll(dir, os.ModePerm); err != nil {
		return fmt.Errorf("error creating directories: %v", err)
	}

	file, err := os.Create(fileName)
	if err != nil {
		return fmt.Errorf("error creating file: %w", err)
	}
	defer file.Close()

	// Create a new buffered writer
	writer := bufio.NewWriter(file)

	// Write each line to the file
	for _, line := range lines {
		_, err := writer.WriteString(line + "\n")
		if err != nil {
			return fmt.Errorf("error writing to file: %w", err)
		}
	}

	// Flush buffered data to the file
	if err := writer.Flush(); err != nil {
		return fmt.Errorf("error flushing buffer: %w", err)
	}

	return nil
}

// Kmerize reads fasta or fastq file, turns into kmers and counts them
// filtering is implemented but would require the seqs.seqfilter list to be non-nil, needs re-testing
func (seqs *Sequences) Kmerize(kmers *Oligos) {
	var name, seq string
	var count, lcount int

	// reading stuff and setup
	// now := time.Now()
	// fmt.Println("opening file at time", now)
	// fmt.Println("File to open is ", seqs.seqfile)
	fpin, err := os.Open(seqs.seqfile)
	globals.Check(err)
	defer fpin.Close()
	scanner := bufio.NewScanner(fpin)
	entrycount := 0 // track line in each entry
	entrylimit := 0 // default every line in entry is counted
	position := 0
	// hackish fastq, just reading the first line after name line
	if seqs.filetype == "fastq" {
		entrylimit = 1
	}
	kmers.name = seqs.seqfile
	// fmt.Println("In Kmerize, dofilter ", seqs.dofilter)
	// fmt.Println("File type is ", seqs.filetype, "and entry limit is", entrylimit)
	passfilter := false // flag to see if name is in filter list

	// read, record, count kmers
	for scanner.Scan() {
		lcount += 1
		line := scanner.Text()              // should not include eol
		trimline := strings.TrimSpace(line) // trim off leading and lagging whitespace
		nofilter := !seqs.dofilter
		if strings.HasPrefix(line, seqs.entrystart) {
			name = strings.TrimPrefix(trimline, seqs.entrystart)
			count += 1
			entrycount = 1
			kmers.remnant = ""
			//fmt.Println("New seq", name, "number", count)
			if seqs.dofilter {
				passfilter = checkfilter(seqs.seqfilter, name)
			}
		} else {
			if (nofilter || passfilter) && (lcount < seqs.linelimit) && (entrylimit < 1 || entrycount <= entrylimit) {
				seqs.Totallen += len(trimline)
				seq = kmers.remnant + trimline
				entrycount++
				if (lcount <= seqs.linelimit) && (lcount > seqs.linemin) {
					// kmers.Countref(seqs, seq, name)
					kmers.getmersInFragment(seqs, seq, name, position)
					position += (len(seq) - kmers.klen + 1)
				}
				keepcompany(lcount, len(seq), 50000, 10000, 50000)
			}
		}
	}

	// fmt.Println("Seqs and Lines counted\n", count, lcount)
}

// keepcompany outputs count and seqlen at intervals
// we could control this through globals, but...
func keepcompany(count int, seqlen int, early int, earlyinterval int, interval int) {
	if count < early {
		if (count % earlyinterval) == 0 {
			fmt.Println(count, seqlen)
		}
	}
	if (count % interval) == 0 {
		fmt.Println(count, seqlen)
	}
}

//
// Basic Tools
//

// Index is from go by example to find the index positio of a string match in a list
func Index(vs []string, t string) int {
	for i, v := range vs {
		if v == t {
			return i
		}
	}
	return -1
}

// strmatch checks that strings match or exits
func strmatch(query string, match string, words string) {
	if query != match {
		fmt.Println("Exiting string comparison failure upon read", query, match, words)
		os.Exit(123)
	}
}

// rc returns reverse complement
func rc(kmer string) string {
	kbits := []byte(kmer)
	rcbits := []byte(kmer)
	klen := len(kmer)
	for i := range kbits {
		for n := range nucs {
			if nucs[n] == kbits[i] {
				rcbits[klen-i-1] = cnucs[n]
			}
		}
	}
	return string(rcbits)
}

func counthits2(intslice []int) int {
	var hits int
	if intslice != nil {
		for i := range intslice {
			if intslice[i] != 0 {
				hits += 1
			}
		}
	}
	return hits
}

func counthits(intslice []int) int {
	var hits int
	if intslice != nil {
		for i := range intslice {
			if intslice[i] > 0 {
				hits += 1
			}
		}
	}
	return hits
}

func getsum(intslice []int) int {
	var sum int
	if intslice != nil {
		for i := range intslice {
			if intslice[i] > 0 {
				sum += intslice[i]
			}
		}
	}
	return (sum)
}

func sumdiffsqr(intslice []int, mean float64) float64 {
	var sumsqr, realval, diff float64
	if intslice != nil {
		for i := range intslice {
			if intslice[i] > 0 {
				realval = float64(intslice[i])
				diff = realval - mean
				sumsqr += diff * diff
			}
		}
	}
	return (sumsqr)
}
