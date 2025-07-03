# rabr-biofilm-analysis

## Overview
What the project does
Why the project is useful
## Installation
How users can get started with the project
Where users can get help with your project
Who maintains and contributes to the project




# POPULATION ANALYSIS OF CYANOBACTERIA IN ROTATING ALGAL BIOFILM REACTORS

## List of files within Registry:
### 16S Analysis Files
16S.sh – Executable file that sends the 16S batch file to the high performance computing center.\
16S.oligos – File that contains the 16S primers. Used by the 16S batch file.\
16S.batch – Classifies the 16S sequences using mothur and Silva.\
A_DADA2.R R – Classifies the lab scale, CVWRF, greenhouse RABRs, and the trickling filter 16S sequences using DADA2.\
E_DADA2.R – Classifies the pilot and control RABR 16S sequences using DADA2.\
tax_versus.R – Compares the mothur and DADA2 analysis results for the 16S sequences.\
16S_abund.R – Generates relative abundance graphs of 16S sequences from the mothur analysis. Was manually edited and rerun to select for each taxonomic level (Phylum, Class, Order, Family, Genus). Resulting figures were color adjusted to match with relative abundance graphs from the 16S DADA2 analysis.\
16S_alphaD.R – Generates alpha diversity and productivity graphs of 16S sequences from the mothur analysis.\
16S_betaD.R – Generates beta diversity graphs of 16S sequences from the mothur analysis.\
16S_lefse.R – Generates LEfSe graphs for 16S sequences from the mothur analysis. Has commented sections that contain code that must be run through mothur to proceed.\
16S_rarefaction.R – Generates rarefaction graphs of 16S sequences from the mothur analysis. Included a figure that is colored by section, as well as individual graphs for each section.\
16S_StrainProductivity.R - Generates figures that compare specific 16S strain relative abundance to sample productivity from the mothur analysis.\
16S_metadata.xlsx - Contains metadata used for mothur 16S R files.

16S_D_abund.R – Generates relative abundance graphs of 16S sequences from the DADA2 analysis. Was manually edited and rerun to select for each taxonomic level (Phylum, Class, Order, Family, Genus). Resulting figures were color adjusted to match with relative abundance graphs from the 16S mothur analysis.\
16S_ D_alphaD.R – Generates alpha diversity and productivity graphs of 16S sequences from the DADA2 analysis.\
16S_ D_betaD.R – Generates beta diversity graphs of 16S sequences from the DADA2 analysis.\
16S_ D_lefse.R – Generates LEfSe graphs for 16S sequences from the DADA2 analysis. Has commented sections that contain code that must be run through mothur to proceed.\
16S_ D_rarefaction.R – Generates rarefaction graphs of 16S sequences from the DADA2 analysis. Included a figure that is colored by section, as well as individual graphs for each section.\
16S_D_StrainProductivity.R - Generates figures that compare specific 16S strain relative abundance to sample productivity from the DADA2 analysis.\
16S_D_metadata.xlsx - Contains metadata used for the DADA2 16S R files.

### 18S Analysis Files
18S.sh – Executable file that sends the 18S batch file to the high performance computing center.\
18S.oligos – File that contains the 18S primers. Used by the 18S batch file.\
18S.batch – Classifies the 18S sequences using mothur and Silva.\
18S_abund.R – Generates relative abundance graphs of 18S sequences. Automatically returns figures for each taxonomic level (Supergroup, Division, Subdivision, Class, Order, Family, Genus).\
18S_alphaD.R – Generates alpha diversity and productivity graphs of 18S sequences.\
18S_betaD.R – Generates beta diversity graphs of 18S sequences.\
18S_lefse.R – Generates LEfSe graphs for 18S sequences. Has commented sections that contain code that must be run through mothur to proceed.\
18S_rarefaction.R – Generates rarefaction graphs of 18S sequences. Included a figure that is colored by section, as well as individual graphs for each section.\
18S_StrainProductivity.R - Generates figures that compare specific 18S strain relative abundance to sample productivity.\
18S_metadata.xlsx - Contains metadata used for 18S R files.

### 23S Analysis Files
23S_DADA2.R – Classifies the pilot and control RABR 23S sequences using DADA2.\
A23S_DADA2.R – Classifies all remaining biofilm 23S sequences using DADA2. This includes the lab scale, CVWRF, greenhouse RABRs, and the trickling filter.\
23S_D_abund.R – Generates relative abundance graphs of 23S sequences. Was manually edited and rerun to select for each taxonomic level (Phylum, Class, Order, Family, Genus).\
23S_D_alphaD.R – Generates alpha diversity and productivity graphs of 23S sequences.\
23S_D_betaD.R – Generates beta diversity graphs of 23S sequences.\
23S_D_lefse.R – Generates LEfSe graphs for 23S sequences. Has commented sections that contain code that must be run through mothur to proceed.\
23S_D_rarefaction.R – Generates rarefaction graphs of 23S sequences. Included a figure that is colored by section, as well as individual graphs for each section.\
23S_StrainProductivity.R - Generates figures that compare specific 23S strain relative abundance to sample productivity.\
23S_metadata.xlsx - Contains metadata used for 23S R files.

### ITS Analysis Files
ITS.sh – Executable file that sends the ITS_fungi and ITS_euk batch files to the high performance computing center.\
ITS.oligos – File that contains the ITS primers. Used by the ITS_fungi and ITS_euk batch files.\
ITS_fungi.batch – Classifies the ITS sequences using mothur and Silva and PR2.\
ITS_euk.batch – Classifies the ITS sequences using mothur and Silva and PR2.\
ITS_abund.R – Generates relative abundance graphs of ITS. Automatically returns figures for each taxonomic level (Supergroup, Division, Subdivision, Class, Order, Family, Genus).\
ITS_alphaD.R – Generates alpha diversity and productivity graphs of ITS sequences.\
ITS_betaD.R – Generates beta diversity graphs of ITS sequences.\
ITS_lefse.R – Generates LEfSe graphs for ITS sequences. Has commented sections that contain code that must be run through mothur to proceed.\
ITS_rarefaction.R – Generates rarefaction graphs of ITS sequences. Included a figure that is colored by section, as well as individual graphs for each section.
ITS_StrainProductivity.R - Generates figures that compare specific ITS strain relative abundance to sample productivity.\
ITS_metadata.xlsx - Contains metadata used for ITS R files\
