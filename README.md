# Hydroplane1a
Hydroplane is a software platform designed to analyze evolutionary relationships among microbial and viral genomes using metagenomics data. The framework clusters homologous (matching) kmers to identify genome-wide homologous regions and model evolutionary divergence. By assessing the spatial proximity of kmers, Hydroplane detects both identical and divergent sequences, providing detailed insights into genome clusters at multiple hierarchial levels. Hydroplane achieves efficient genome-wide analysis with O(L) complexity, where L is the total length of sequence reads, enabling scalability on standard desktop systems. The implementation in the Go programming language ensures ease of use, minimal reliance on external libraries, and cross-platform compatibility across MacOS, Unix, and Windows. 

see also https://github.com/NWBRaVE/Hydroplane

## Operation and Description
Code is written in the language Go $^1$. The main program is in folder hydroplane, while the globals and seqmer packages are in globals_hydro and seqmer_hydro respectively, along with their respective required go.mod and go.sum files. The hydroplane folder also contains factory, mode, and control files, which set up the desired directory structure and parameters to run the program. These files should exist in the folder from which you run the program. The reference directory and sequence directory in the example control file are set to ../references/109PromarinusGenomes18Nov2023/ so those directories should exist when running or change the refdir and seqdir parameters as desired to point to reference sequences on your computer from where you will run the program. In the control files, everything after a `#` symbol is a comment, and further description of parameters are contained in the control file and/or as comments in the code. The folder 109PromarinusGenomes18Nov2023/ contains two of the genomes (not all 109) mentioned in the paper for example running purposes, and the program will recursively search for files starting with "GCA" and ending with "fna", which are expected to be in "fasta" format (editable parameters seqprefix, seqext, and geneseqext, respectively). Further genomes or other reference sequences can be downloaded from public repositories (e.g. https://www.ncbi.nlm.nih.gov/datasets/genome/ ) and placed in this folder or accessed by renaming the seqdir parameter. 

Because this is a reasearch program, all program parameters are editable, providing maximum user power and flexibility, but there are inherent dangers in this so you should change them at your own risk. 

The hydroplane program is under continued development, and we aim to provide more examples of running modes with fewer parameters changed in the control file, and further examples and explanation on how to run it. If you need further information or want to encourage us to hurry up, please contact the authors. 

$^1$ Please follow instructions at https://go.dev/doc/install to install Go if not already installed. 

## References
Hydroplane I: one-shot probabilistic evolutionary analysis for scalable organizational identification

Connah G. M. Johnson $^1$, David D. Pollock $^2$ $^3$ $^*$ 

$^1$ Physical and Computational Sciences Directorate, Pacific Northwest National Laboratory, Richland, Washington, 99352, USA  
$^2$ Department of Biochemistry and Molecular Genetics, University of Colorado School of Medicine, Aurora, CO, USA.  
$^3$ Biological Sciences Division, Earth and Biological Sciences Directorate, Pacific Northwest National Laboratory, Richland, WA, USA.  

$^*$ Corresponding author
David D. Pollock, david.pollock@cuanschutz.edu, david.pollock@pnnl.gov (I check the second one less often)

_In Review_, Jan 2026

### endnote ref format

@article {Johnson2025.10.14.682387,
	author = {Johnson $^1$, Connah G. M. and Pollock $^2$ $^3$ $^*$, David D.},
	title = {Hydroplane I: one-shot probabilistic evolutionary analysis for scalable organizational identification},
	elocation-id = {2025.10.14.682387},
	year = {2025},
	doi = {10.1101/2025.10.14.682387},
	publisher = {Cold Spring Harbor Laboratory},
	journal = {bioRxiv}
}

## Acknowledments
Go code based from https://github.com/PollockLaboratory/AnCovMulti (https://pubmed.ncbi.nlm.nih.gov/33545711/) and further work from Pollock Laboratory, owes a debt of gratitude to many discussions with Kristen Wade as part of her thesis work (Homology unleashed: nonparametric and alignment-free methods for understanding genome evolution, https://digitalcollections.cuanschutz.edu/work/ns/7ccff632-145a-4f36-99f5-31f9b73bbd25), focused on understanding evolutionary patterns of the non-coding and repetitive content of the genome, which undergo significantly different selective pressures than coding regions. For theoretical background see also https://pubmed.ncbi.nlm.nih.gov/22144907/ https://pubmed.ncbi.nlm.nih.gov/32095164/ https://pubmed.ncbi.nlm.nih.gov/28107347/ https://pubmed.ncbi.nlm.nih.gov/18541131/ and https://pubmed.ncbi.nlm.nih.gov/19783593/ . 

This work was most recently supported by the University of Colorado School of Medicine at CUAnschutz and the NW-BRaVE for Biopreparedness project funded by the U. S. Department of Energy (DOE), Office of Science, Office of Biological and Environmental Research, under FWP 81832. 










