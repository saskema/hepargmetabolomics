# Supplementary Material
The code provided by this repository is part of the supplementary material for the scientific article "Toxicometabolomics of the new psychoactive substances α-PBP and α-PEP studied in HepaRG cell incubates by means of untargeted metabolomics revealed unexpected amino acid adducts " by Sascha K. Manier et al., Arch Toxicol, 2020 (DOI:10.1007/s00204-020-02742-1 ). See the following sections for detailed information about it. Corresponding files can be found at MetaboLights using the study identifier MTBLS1481. Please note that the corresponding study at MetaboLights might still be in curation and thus not yet accessible.

# Spectra
These are the spectra in mzXML format recorded during the study for each significant feature (if not annotated as an adduct or isotope by CAMERA) using a CE of 10, 20, and 40.

# Source Code
## Dependencies
The code deployed by this repository requires the following packages:

- tidyverse (CRAN)
- ggrepel (CRAN)
- gplots (CRAN)
- MASS (CRAN)
- caret (CRAN)
- XCMS (Bioconductor)
- CAMERA (Bioconductor)

 If not installed you can use the following lines to install them:
 
	 `install.packages(c("tidyverse", "ggrepel", "gplots", "MASS", "caret")
	 source("https://bioconductor.org/biocLite.R")
	 biocLite(c("xcms", "CAMERA"))`
	
Additionally, evaluation script require functions provided by the metaboLib.R library. The corresponding library can be found [here](https://github.com/saskema/metaboLib). Please note that metabolib.R has additional dependencies.

## aPBP\_cells\_neg.R etc.
### General
These scripts were used to evaluate the sample analyses by performing peak picking and subsequent staistical analysis. Scripts that contain "cells" in their name were used for those samples after preparation of the HepaRG cells, scripts that contain "media" were used for those samples that were obtained from the corresponding surrounding cell media. The suffixes "pos" and "neg" describe their usage for each ionization mode.

### Usage
To use the scripts as deployed in this repository, place the matching script in the same folder were the files are grouped in folders corresponding to their study group ("Blank", "Low", "High") and run it.

### Result
You will obtain the following graphics as PDF files:

- evaluateResults/evaluate.pdf
- evaluateResults/evaluate.RData

## centWaveOpt.R
### General
This script was used to perform peak picking optimization according to DOI:10.1002/dta.2552.
