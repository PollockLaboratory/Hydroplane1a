package main

import (
	"bufio"
	"flag"
	"fmt"
	"globals"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"time"

	"seqmer"

	"gonum.org/v1/plot/plotter"
)

var prog program // program structure defined below

type program struct {
	id       int
	name     string
	nickname string
	authors  string
	began    string
	modified string
	uses     string
	runrec   string
	computer string
}

// these are the default files that hold parameters and control the program settings
const factfile = "factory"
const controlfile = "control"
const modefile = "mode"

// Set adds and prints program information
func (p *program) Set(writer io.Writer) {
	p.name = "hydro_plane" // count reference kmers from seqs file
	fmt.Fprintln(writer, "\n\tRunning Program", p.name)
	p.authors = "David Pollock, Connah Johnson"
	p.began = "January 20, 2025"
	p.modified = "May 31, 2025"
	fmt.Fprintln(writer, "\tAuthors:", p.authors, "Last Modified:", p.modified)
}

func main() {

	//-------------------------------------------------------------------------//
	//						Process the control file						   //
	//-------------------------------------------------------------------------//

	now := time.Now()
	fmt.Println("Beginning main program", now)
	fmt.Println("Setting access.log File to append")
	fappend, _ := os.OpenFile("access.log", os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)
	defer fappend.Close()
	writer := bufio.NewWriter(fappend)
	defer writer.Flush() // need this to get output

	// set up the program and globals by reading controls or flag line input
	// kcount3 -mode testmode -control testcontrol
	prog.Set(writer)
	modeptr := flag.String("mode", modefile, "source for mode file")
	contptr := flag.String("control", controlfile, "source for control file")
	flag.Parse()
	fmt.Fprintln(writer, "\n\tMode file is: ", *modeptr)
	fmt.Fprintln(writer, "\n\tControl file is: ", *contptr)
	var globs = globals.New()
	globs.ProgSetUp(*contptr, *modeptr) // should change through factory and command line
	fmt.Fprintln(writer, "The program was run on", globs.Runstart)
	globs.Print(os.Stdout, "\nStatus after Setup")

	fmt.Println("Program mode:")
	fmt.Println("Kmerize mode:	", globs.Getb("kmerizemode"))
	fmt.Println("Read split mode:	", globs.Getb("readsplitmode"))
	fmt.Println("Match sequence mode:	", globs.Getb("matchsequencemode"))
	fmt.Println("Target gene mode:	", globs.Getb("targetgenemode"))
	fmt.Println("First sequence as a reference mode:	", globs.Getb("firstsequencereference"))
	fmt.Println("Start program.")

	kcountdir := globs.Getf("kcountdir")
	pcountdir := globs.Getf("pcountdir")
	kcountfile_cntrl := globs.Getf("kcountfile")
	pcountfile_cntrl := globs.Getf("pcountfile")
	klen_cntrl := globs.Geti("klen")

	// get set of query sequence file names
	fmt.Println("Pre seq route and filepath", globs.Gets("refdir"), globs.Gets("seqdir"), globs.Gets("seqprefix"), globs.Gets("seqext"))
	refseqroute := getallext(globs.Gets("seqprefix"), globs.Gets("seqext"), globs.Gets("refdir")+globs.Gets("seqdir"))
	refseqfiles, err := filepath.Glob(refseqroute) // returns names of all files that match route, including the path
	printerr(err)

	// the firstseq file name needs to be integrated with the file pah if this is to work
	firstseq := globs.Gets("firstseq")
	if firstseq != "no" {
		refseqfiles = append([]string{firstseq}, refseqfiles...)
	}
	fmt.Println("seq route and filepath", refseqroute, refseqfiles)
	for count, afilepath := range refseqfiles {
		fmt.Println("\n\nchecking file number (filepath)", count+1, afilepath)
	}

	// storm parameters
	extendside := globs.Geti("extendtoside")
	stormsize := globs.Geti("stormsize") // return kmers.storms if already exists

	// os.Exit(0)

	// ifFirstPass := true

	// process the reference genome to populate the storms and act as the sequence to compare all others against
	// referenceSeqsSets := new(seqmer.SeqSets) // create global;
	// referenceSeqsSets.Init(globs.Gets("filetype"), globs.Getb("dofilter"), globs.Geti("linelimit"), klen_cntrl)

	ifRunSim := true

	if ifRunSim {
		// for each reference sequence
		for count_outer, afilepath_outer := range refseqfiles {
			fmt.Println(count_outer)
			fmt.Println(afilepath_outer)
			dir2, file2_outer := filepath.Split(afilepath_outer)
			ext2_outer := filepath.Ext(afilepath_outer)
			parts := strings.SplitN(dir2, globs.Gets("seqdir"), 2)
			// sequence_dir := strings.TrimSuffix(parts[1], "/")
			// newOutDir := filepath.Join(globs.Gets("outdir"), sequence_dir) // "output/GCA_949511955.1"
			newOutDir := filepath.Join(globs.Gets("outdir"), parts[1]) // "output/GCA_949511955.1"
			newOutDir = filepath.Clean(newOutDir)
			basename2_outer := strings.TrimSuffix(file2_outer, ext2_outer)
			newOutDir = newOutDir + "/"
			if err := os.MkdirAll(newOutDir, 0o755); err != nil {
				fmt.Printf("mkdir %s", newOutDir)
			}
			fmt.Println("directory ready:", newOutDir)

			//-------------------------------------------------------------------------//
			//					Set up global structures to store data				   //
			//-------------------------------------------------------------------------//

			// set up storms global
			storms := new(seqmer.Tempest)
			storms.Init(stormsize) // this could have been done when storms was created

			rels := seqmer.Newrel()

			now = time.Now()
			// fmt.Println("About to create all oligos for klen ", klen_cntrl, now)

			alligos := seqmer.Newalligos(klen_cntrl, globs.Getb("dofindalligos"), globs.Getb("nofindalligos")) // nofind controls if cum
			revigos := seqmer.NewRevAlligos(klen_cntrl)                                                        // reverse oligo map to find kmers given kptr
			goodmin := globs.Geti("goodmin")

			//-------------------------------------------------------------------------//
			//					Process reference file to Storms					   //
			//-------------------------------------------------------------------------//
			seqs := new(seqmer.Sequences)   // create global;
			seq_kmers := new(seqmer.Oligos) // create global;
			seqs.Init(afilepath_outer, globs.Geti("minseqlen"), globs.Geti("linelimit"), globs.Geti("minline"), globs.Getb("recordseq"), globs.Gets("filetype"), globs.Getb("dofilter"))

			kcountfile := getkoutname(kcountfile_cntrl, klen_cntrl, basename2_outer)
			pcountfile := getpoutname(pcountfile_cntrl, klen_cntrl, basename2_outer)

			seq_kmers.Init(klen_cntrl, kcountdir+kcountfile, pcountdir+pcountfile, globs.Getb("printNs"), globs.Geti("kminprint"), globs.Getb("dofindalligos"), alligos)
			// if taking in reference sequences then kmerize, else use the already kmerized sequences
			if globs.Getb("kmerizemode") {
				// if (globs.Getb("kmerizemode") == true) && (ifFirstPass == true) {
				// create new seqs and seq_kmers, Init() them
				seqs.Kmerize(seq_kmers)
				seq_kmers.Kprint()
				seq_kmers.Kposprint(globs.Geti("kmin"), globs.Geti("kmax"))
			}
			seq_kmers.SeedStorms(globs.Geti("stormsize"), storms) // only triggers for the first sequence call

			weatherForecast := new(seqmer.WeatherForecast)

			//-------------------------------------------------------------------------//
			//			Process each sequence file									   //
			//-------------------------------------------------------------------------//

			now = time.Now()
			globs.Delta()
			// fmt.Println("About to go through files", now)
			// if count_outer <= 10 {
			// 	continue
			// }
			fmt.Println("===========================================================")
			fmt.Println("===========================================================")
			fmt.Println("===========================================================")
			fmt.Println("===========================================================")
			fmt.Println("===========================================================")

			// for each reference sequence
			for count, afilepath := range refseqfiles {

				// set up timeing
				timeKmerize := 0.0
				timeStormMatch := 0.0
				timeStormPlot := 0.0
				meanMatchScore := 0.0
				varianceMatchScore := 0.0
				maxMatchScore := 0.0
				numberOfMatchPoints := 0
				sequenceLength := 0
				matchPlotPoints := plotter.XYs{}
				bedPlotPoints := plotter.XYs{}

				//-------------------------------------------------------------------------//
				//			Sequence I/O and convert to Kmers and positions				   //
				//-------------------------------------------------------------------------//

				// set up file name and path info
				// fmt.Println("\n\nreading Xfile number (filepath)", count+1, afilepath)
				// fmt.Println()
				_, file2 := filepath.Split(afilepath)
				ext2 := filepath.Ext(afilepath)
				basename2 := strings.TrimSuffix(file2, ext2)
				// fmt.Println("\n\n Printing Xfile 2 path, dir", dir2, file2, ext2, basename2)
				// set up seqs and seq_kmers structures for current file and kmerize or read in kmer counts and positions
				seqs := new(seqmer.Sequences)   // create global;
				seq_kmers := new(seqmer.Oligos) // create global;
				seqs.Init(afilepath, globs.Geti("minseqlen"), globs.Geti("linelimit"), globs.Geti("minline"), globs.Getb("recordseq"), globs.Gets("filetype"), globs.Getb("dofilter"))

				kcountfile := getkoutname(kcountfile_cntrl, klen_cntrl, basename2)
				pcountfile := getpoutname(pcountfile_cntrl, klen_cntrl, basename2)

				seq_kmers.Init(klen_cntrl, kcountdir+kcountfile, pcountdir+pcountfile, globs.Getb("printNs"), globs.Geti("kminprint"), globs.Getb("dofindalligos"), alligos)

				// if taking in reference sequences then kmerize, else use the already kmerized sequences
				if globs.Getb("kmerizemode") {
					// create new seqs and seq_kmers, Init() them
					preKmerize := time.Now()
					seqs.Kmerize(seq_kmers)
					timeKmerize = time.Since(preKmerize).Seconds()
					seq_kmers.Kprint()
					seq_kmers.Kposprint(globs.Geti("kmin"), globs.Geti("kmax"))

					sequenceLength = seqs.Totallen
				}

				// fmt.Println("Finished Kmers x for file:", basename2)

				//-------------------------------------------------------------------------//
				//					Build storms by comparing to global storms			   //
				//-------------------------------------------------------------------------//

				if globs.Getb("build_storms") {
					// fmt.Println("building storms now ", now, globs.Getb("build_storms"))
					// now = time.Now()
					globs.Delta()

					// compare storms with the new sequences, determine hit matches
					preStormMatch := time.Now()
					stormhits := seq_kmers.FindHits(storms, globs.Geti("boundary"), extendside, goodmin, globs.Geti("stormsize")) //
					timeStormMatch = time.Since(preStormMatch).Seconds()

					// fmt.Println("\n\t\ttrying to be Printing stormologs", strconv.Itoa(klen_cntrl), strconv.Itoa(goodmin))
					stormoutfile := getsoutname(goodmin, stormsize, extendside, klen_cntrl, basename2)
					stormhits.Print(newOutDir+stormoutfile, revigos, goodmin) // don't do this if need for speed

					//storms.Build(stormhits, seq_kmers, buildminhits)               // build up storms from gaps in matches
					//			observigos := seq_kmers.make_observedkmers_hash()
					revigos.Reversefill(alligos)
					if count < 1 {
						storms.Print(revigos, 2)
					}
					rels.Addrel(stormhits)
					now = time.Now()
					globs.Delta()

					// fmt.Println("done building storms now ", now, globs.Getb("build_storms"))
					if globs.Getb("plot_storms_mode") {
						// provide plots for storm matches
						now = time.Now()
						globs.Delta()
						// fmt.Println("plotting storms now ", now, globs.Getb("plot_storms_mode"))
						stormplotfile := getsplotname(goodmin, stormsize, extendside, klen_cntrl, basename2)
						preStormPlot := time.Now()

						meanMatchScore, varianceMatchScore, maxMatchScore, numberOfMatchPoints, matchPlotPoints, bedPlotPoints = stormhits.PlotStorms(newOutDir+stormplotfile, revigos, goodmin, globs.Getn("scattersize"), globs.Geti("plotxdim"), globs.Geti("plotydim"), globs.Geti("plottickfontsize"), globs.Geti("plotaxisfontsize"), globs.Geti("movingaveragewindow"), globs.Geti("plotlinewidth"), globs.Getn("plottickspacing"))

						timeStormPlot = time.Since(preStormPlot).Seconds()
						now = time.Now()
						globs.Delta()
						// fmt.Println("done plotting storms now ", now, globs.Getb("plot_storms_mode"))
					}
				}

				// send data to weatherSystem
				weatherSystem := new(seqmer.WeatherSystem)
				_, file2 = filepath.Split(storms.SeedName)
				ext2 = filepath.Ext(storms.SeedName)
				basename_ref := strings.TrimSuffix(file2, ext2)
				weatherSystem.Init(basename_ref, basename2, matchPlotPoints, bedPlotPoints, timeKmerize, timeStormMatch, timeStormPlot, meanMatchScore, varianceMatchScore, maxMatchScore, numberOfMatchPoints, sequenceLength, stormsize, storms.GetNumberOfStorms())

				// send weather system to weather forecast
				weatherForecast.Add(weatherSystem)

				// plot combined figures, heatmap, layered matchplots

				// fixoutputs to be into directories

				//we should probably make sure that seq_kmers and seqs are disposed of here
				fmt.Println("Finished Kmers for file:", basename2)
			}

			// ifFirstPass = false

			// weatherForecast.Print()

			weatherForecast.PlotWeatherForecast(newOutDir, globs.Getn("scattersize"), globs.Geti("plotxdim"), globs.Geti("plotydim"), globs.Geti("plottickfontsize"), globs.Geti("plotaxisfontsize"), globs.Geti("movingaveragewindow"), globs.Geti("plotlinewidth"), globs.Getn("plottickspacing"))

			weatherForecast.Save(newOutDir)

		}

	} else {

		if err := seqmer.UpdateForecastsInOutput(globs.Gets("outdir"), globs.Getn("stormsize")); err != nil {
			fmt.Println("Error correcting forecasts: %v", err)
		}
		// display
		portfolio, err := seqmer.LoadForecastsFromDir(globs.Gets("outdir"))
		if err != nil {
			fmt.Println("cannot load forecasts: %v", err)
		}

		seqmer.PlotTimeKmerizeVsRefLength(
			portfolio.Forecasts,
			globs.Gets("outdir")+"all_referenceLength_kmerizeTime.png",
			globs.Getn("scattersize"), globs.Geti("plotxdim"), globs.Geti("plotydim"), globs.Geti("plottickfontsize"), globs.Geti("plotaxisfontsize"), globs.Geti("movingaveragewindow"), globs.Geti("plotlinewidth"), globs.Getn("plottickspacing"))
		seqmer.PlotTimeStormVsRefLength(
			portfolio.Forecasts,
			globs.Gets("outdir")+"all_referenceLength_stormTime.png",
			globs.Getn("scattersize"), globs.Geti("plotxdim"), globs.Geti("plotydim"), globs.Geti("plottickfontsize"), globs.Geti("plotaxisfontsize"), globs.Geti("movingaveragewindow"), globs.Geti("plotlinewidth"), globs.Getn("plottickspacing"))

		fmt.Println("PlotAllForecastsHeatmap")
		seqmer.PlotAllForecastsHeatmap(
			portfolio.Forecasts,
			globs.Gets("outdir")+"all_forecasts_heatmap.png",
			globs.Getn("scattersize"), globs.Geti("plotxdim"), globs.Geti("plotydim"), globs.Geti("plottickfontsize"), globs.Geti("plotaxisfontsize"), globs.Geti("movingaveragewindow"), globs.Geti("plotlinewidth"), globs.Getn("plottickspacing"), globs.Getn("stormsize"))
		seqmer.PlotAllForecastsCorrectedHeatmap(
			portfolio.Forecasts,
			globs.Gets("outdir")+"all_forecasts_corrected_heatmap.png",
			globs.Getn("scattersize"), globs.Geti("plotxdim"), globs.Geti("plotydim"), globs.Geti("plottickfontsize"), globs.Geti("plotaxisfontsize"), globs.Geti("movingaveragewindow"), globs.Geti("plotlinewidth"), globs.Getn("plottickspacing"), globs.Getn("stormsize"))

		seqmer.PlotAllForecastsHitFrequencyHeatmap(
			portfolio.Forecasts,
			globs.Gets("outdir")+"all_forecasts_hit_frequency_heatmap.png",
			globs.Getn("scattersize"), globs.Geti("plotxdim"), globs.Geti("plotydim"), globs.Geti("plottickfontsize"), globs.Geti("plotaxisfontsize"), globs.Geti("movingaveragewindow"), globs.Geti("plotlinewidth"), globs.Getn("plottickspacing"), globs.Getn("stormsize"))

		// seqmer.PlotAllForecastsRefCoverageHeatmap(
		// 	portfolio.Forecasts,
		// 	globs.Gets("outdir")+"all_forecasts_reference_coverage_heatmap.png",
		// 	globs.Getn("scattersize"), globs.Geti("plotxdim"), globs.Geti("plotydim"), globs.Geti("plottickfontsize"), globs.Geti("plotaxisfontsize"), globs.Geti("movingaveragewindow"), globs.Geti("plotlinewidth"), globs.Getn("plottickspacing"), globs.Getn("stormsize"))

	}
	// save storms to file
	// storms.Print2(revigos)
	// rels.Print(globs.Gets("outdir")+globs.Gets("reloutfile"), revigos, goodmin)

	// globs.Print(writer, "\nStatus after Program Completion")
	// now = time.Now()
	// globs.Delta()
	// fmt.Fprintln(writer, "\nThe program finished", now)

	// if globs.Getb("plot_relatives_mode") {

	// 	// provide plots for storm matches
	// 	rels.PlotRelatives(globs.Gets("outdir")+globs.Gets("relplotfile"), revigos, goodmin)
	// }

}

func processFile(count int, afilepath string) error {
	// --- Everything below is verbatim from your loop -------------------------

	fmt.Println("\n\nreading Xfile number (filepath)", count+1, afilepath)
	dir2, file2 := filepath.Split(afilepath)
	ext2 := filepath.Ext(afilepath)
	basename2 := strings.TrimSuffix(file2, ext2)
	fmt.Println("\n\n Printing Xfile 2 path, dir", dir2, file2, ext2, basename2)

	/*  ... the rest of your original code ...  */

	time.Sleep(200 * time.Millisecond) // pretend work
	fmt.Println("Finished Kmers for file:", basename2)
	return nil
}

// func PlotScatterKmerMatch(matchData []ReadMatch, kmerPlotDir string, kmerPlotFileName string, querygenome string, testgenome string, plotxdim int, plotydim int, scattersize float64) {

// 	// Prepare data for plotting
// 	var xSequencePositions []int
// 	var yGenomePositions []int

// 	for _, read := range matchData {
// 		xPosition := extractXPosition(read.ReadName)

// 		if len(read.KmerDetails) > 0 {
// 			for _, kmerDetail := range read.KmerDetails {
// 				// Calculate global positions by adding x_position
// 				globalPositions := make([]int, len(kmerDetail.ReadPositions))
// 				for i, readPos := range kmerDetail.ReadPositions {
// 					globalPositions[i] = xPosition + readPos
// 				}
// 				genomePositions := kmerDetail.GenomePositions

// 				// For plotting, we need direct pairing of global positions with genome positions
// 				for j := 0; j < len(globalPositions) && j < len(genomePositions); j++ {
// 					xSequencePositions = append(xSequencePositions, globalPositions[j])
// 					yGenomePositions = append(yGenomePositions, genomePositions[j])
// 				}
// 			}
// 		}
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

// func PlotScatterReadMatch(matchData []ReadMatch, readmatchPlotDir string, readmatchPlotFileName string, querygenome string, testgenome string, plotxdim int, plotydim int, scattersize float64) {

// 	// Process data for plotting
// 	var xPositions []int
// 	var matchScores []int

// 	for _, read := range matchData {
// 		xPosition := extractXPosition(read.ReadName)
// 		xPositions = append(xPositions, xPosition)
// 		matchScores = append(matchScores, read.MatchScore)
// 	}

// 	// Create a new plot
// 	p := plot.New()

// 	p.Title.Text = "Read Match Positions: " + querygenome + " vs " + testgenome
// 	p.X.Label.Text = querygenome + " Sequence Position"
// 	p.Y.Label.Text = "Match Score"
// 	p.Add(plotter.NewGrid())

// 	// Create a scatter plot
// 	points := make(plotter.XYs, len(xPositions))
// 	for i := range points {
// 		points[i].X = float64(xPositions[i])
// 		points[i].Y = float64(matchScores[i])
// 	}

// 	scatter, err := plotter.NewScatter(points)
// 	if err != nil {
// 		fmt.Printf("Failed to create scatter plot: %v", err)
// 	}
// 	scatter.GlyphStyle.Color = plotter.DefaultGlyphStyle.Color
// 	scatter.GlyphStyle.Radius = vg.Points(scattersize)

// 	p.Add(scatter)

// 	// Save the plot to a PNG file
// 	os.MkdirAll(readmatchPlotDir, os.ModePerm)
// 	readmatchPlotFileName = readmatchPlotDir + readmatchPlotFileName
// 	if err := p.Save(vg.Length(plotxdim)*vg.Inch, vg.Length(plotydim)*vg.Inch, readmatchPlotFileName); err != nil {
// 		fmt.Printf("Failed to save plot: %v", err)
// 	}

// 	// Save the x and y values to a CSV file
// 	csvFileName := readmatchPlotFileName[:strings.LastIndex(readmatchPlotFileName, ".")] + "_data.csv"
// 	csvFile, err := os.Create(csvFileName)
// 	if err != nil {
// 		fmt.Printf("Failed to create CSV file: %v", err)
// 	}
// 	defer csvFile.Close()

// 	writer := csv.NewWriter(csvFile)
// 	defer writer.Flush()

// 	// Write CSV headers
// 	if err := writer.Write([]string{querygenome + " Sequence Position", "Match Score"}); err != nil {
// 		fmt.Printf("Failed to write CSV headers: %v", err)
// 	}

// 	// Write CSV data
// 	for i := range xPositions {
// 		record := []string{
// 			strconv.Itoa(xPositions[i]),
// 			strconv.Itoa(matchScores[i]),
// 		}
// 		if err := writer.Write(record); err != nil {
// 			fmt.Printf("Failed to write CSV record: %v", err)
// 		}
// 	}

// }

func getpos(s []string, element string) int {
	for i, e := range s {
		if e == element {
			return i
		}
	}
	return -1
}

// getallext gets all files with input extension from the input directory
func getallext(prefix string, extension string, directory string) string {
	//extension := "fq"
	localdir := "./"
	//directory := "seqfiles/"
	dirpath := localdir + directory
	allfiles := "*/*"
	// allfiles := "/" + beginning + "*/*"
	pattern := prefix + allfiles + "." + extension
	route := dirpath + pattern
	return route
}

func getallFullext(extension string, beginning string, directory string) string {
	//extension := "fq"
	localdir := "./"
	//directory := "seqfiles/"
	dirpath := localdir + directory
	allfiles := "*"
	// allfiles := "/" + beginning + "*/*"
	pattern := allfiles + "." + extension
	route := dirpath + pattern
	return route
}

func printerr(err error) {
	if err != nil {
		fmt.Println("File reading error", err)
	}
}

func getkoutname(base string, klen int, seqsource string) string {
	outname := base
	klenstr := strconv.Itoa(klen)
	outname = outname + klenstr + "_" + seqsource + ".xls"
	return outname
}

func getsoutname(goodmin int, stormsize int, extend int, klen int, seqsource string) string {
	goodminstr := strconv.Itoa(goodmin)
	sizestr := strconv.Itoa(stormsize)
	extstr := strconv.Itoa(extend)
	klenstr := strconv.Itoa(klen)
	outname := seqsource + "_" + klenstr + "_" + sizestr + "_" + extstr + "_" + goodminstr + "_" + "Stormhits.xls"
	return outname
}

func getsplotname(goodmin int, stormsize int, extend int, klen int, seqsource string) string {
	goodminstr := strconv.Itoa(goodmin)
	sizestr := strconv.Itoa(stormsize)
	extstr := strconv.Itoa(extend)
	klenstr := strconv.Itoa(klen)
	outname := seqsource + "_" + klenstr + "_" + sizestr + "_" + extstr + "_" + goodminstr + "_" + "stormplot.png"
	return outname
}

func getpoutname(base string, klen int, seqsource string) string {
	outname := base
	klenstr := strconv.Itoa(klen)
	outname = outname + klenstr + "_" + seqsource + ".xls"
	return outname
}

func getkoutnamesplit(base string, klen int, seqsource string, splitname string) string {
	outname := base
	klenstr := strconv.Itoa(klen)
	outname = outname + klenstr + "_" + splitname + "_" + seqsource + ".xls"
	return outname
}

func getpoutnamesplit(base string, klen int, seqsource string, splitname string) string {
	outname := base
	klenstr := strconv.Itoa(klen)
	outname = outname + klenstr + "_" + splitname + "_" + seqsource + ".xls"
	return outname
}

// func printGlobalKmerRepo() {
// 	for kmerName, _ := range globalKmerRepo {
// 		fmt.Println(kmerName)
// 	}
// }
