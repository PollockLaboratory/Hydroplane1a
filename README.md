# Hydroplane1a
preliminary version from the hydroplane project

Hydroplane is a software platform designed to analyze evolutionary relationships among microbial and viral genomes using metagenomics data. The framework clusters homologous (matching) kmers to identify genome-wide homologous regions and model evolutionary divergence. By assessing the spatial proximity of kmers, Hydroplane detects both identical and divergent sequences, providing detailed insights into genome clusters at multiple hierarchial levels. Hydroplane achieves efficient genome-wide analysis with O(L) complexity, where L is the total length of sequence reads, enabling scalability on standard desktop systems. The implementation in the Go programming language ensures ease of use, minimal reliance on external libraries, and cross-platform compatibility across MacOS, Unix, and Windows. 

see also https://github.com/NWBRaVE/Hydroplane

## References
Hydroplane I: one-shot probabilistic evolutionary analysis for scalable organizational identification

Connah G. M. Johnson $^1$, David D. Pollock $^2$ $^3$ $^*$ 

$^1$ Physical and Computational Sciences Directorate, Pacific Northwest National Laboratory, Richland, Washington, 99352, USA  
$^2$ Department of Biochemistry and Molecular Genetics, University of Colorado School of Medicine, Aurora, CO, USA.  
$^3$ Biological Sciences Division, Earth and Biological Sciences Directorate, Pacific Northwest National Laboratory, Richland, WA, USA.  

$^*$ Corresponding author
David D. Pollock, david.pollock@cuanschutz.edu, david.pollock@pnnl.gov (I check the second one less often)

_In Review_, Jan 2026

### ref format

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










