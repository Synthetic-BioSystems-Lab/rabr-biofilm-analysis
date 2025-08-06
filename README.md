# rabr-biofilm-analysis

## Overview
Population analysis of rotating algal bioreactor biofilm 16S, 18S, 23S, ITS sequences. Comparison of mothur and DADA2 methodologies for such analysis.

## Installation
Download all files.
### Running Mothur
Download batch, oligos, and sh files to the same folder as desired fastq files.\
Download taxonomic database.\
Run sh file using mothur.\
Use files starting with 'final' for R Studio analysis.
### Running DADA2
Download DADA2 files to folder with access to desired fastq files and taxonomic database.\
Download any packages the files require.\
Run DADA2 using RStudiio.\
Use seqtab and tax files for R Studio Analysis.
### Running R Studio
Download RStudio files with access to cooresponding analysis files.\
Download any packages the files require.\
Create plots, csvs, and processed data directories for each analysis type.\
Set directory to access proper folders files in each RStudio code.\
Run RStudio files. Lefse R files have commented sections that contain code that must be run through mothur to proceed.\
RStudio files produce figures and csv files.

# POPULATION ANALYSIS OF CYANOBACTERIA IN ROTATING ALGAL BIOFILM REACTORS

## List of files within Registry:
### 16SDADA2 Files
#### DADA2
A_DADA2.R R – Classifies the lab scale, CVWRF, greenhouse RABRs, and the trickling filter 16S sequences using DADA2.\
E_DADA2.R – Classifies the pilot and control RABR 16S sequences using DADA2.
#### RStudio
16S_D_abund.R – Generates relative abundance graphs of 16S sequences from the DADA2 analysis. Automatically returns figures for each taxonomic level (Phylum, Class, Order, Family, Genus).\
16S_ D_alphaD.R – Generates alpha diversity and productivity graphs of 16S sequences from the DADA2 analysis.\
16S_ D_betaD.R – Generates beta diversity graphs of 16S sequences from the DADA2 analysis.\
16S_ D_lefse.R – Generates LEfSe graphs for 16S sequences from the DADA2 analysis. Has commented sections that contain code that must be run through mothur to proceed.\
16S_ D_rarefaction.R – Generates rarefaction graphs of 16S sequences from the DADA2 analysis. Included a figure that is colored by section, as well as individual graphs for each section.\
16S_D_StrainProductivity.R - Generates figures that compare specific 16S strain relative abundance to sample productivity from the DADA2 analysis.\
16S_D_metadata.xlsx - Contains metadata used for the DADA2 16S R files.
### 16Smothur
#### RStudio
16S_abund.R – Generates relative abundance graphs of 16S sequences from the mothur analysis. Automatically returns figures for each taxonomic level (Phylum, Class, Order, Family, Genus).\
16S_alphaD.R – Generates alpha diversity and productivity graphs of 16S sequences from the mothur analysis.\
16S_betaD.R – Generates beta diversity graphs of 16S sequences from the mothur analysis.\
16S_lefse.R – Generates LEfSe graphs for 16S sequences from the mothur analysis. Has commented sections that contain code that must be run through mothur to proceed.\
16S_rarefaction.R – Generates rarefaction graphs of 16S sequences from the mothur analysis. Included a figure that is colored by section, as well as individual graphs for each section.\
16S_StrainProductivity.R - Generates figures that compare specific 16S strain relative abundance to sample productivity from the mothur analysis.\
16S_metadata.xlsx - Contains metadata used for mothur 16S R files.
#### mothur
16S.sh – Executable file that sends the 16S batch file to the high performance computing center.\
16S.oligos – File that contains the 16S primers. Used by the 16S batch file.\
16S.batch – Classifies the 16S sequences using mothur and Silva.
### 16Smothur_vs_DADA2
tax_versus.R – Compares the mothur and DADA2 analysis results for the 16S sequences.
### 18Smothur
#### RStudio
18S_abund.R – Generates relative abundance graphs of 18S sequences. Automatically returns figures for each taxonomic level (Supergroup, Division, Subdivision, Class, Order, Family, Genus).\
18S_alphaD.R – Generates alpha diversity and productivity graphs of 18S sequences.\
18S_betaD.R – Generates beta diversity graphs of 18S sequences.\
18S_lefse.R – Generates LEfSe graphs for 18S sequences. Has commented sections that contain code that must be run through mothur to proceed.\
18S_rarefaction.R – Generates rarefaction graphs of 18S sequences. Included a figure that is colored by section, as well as individual graphs for each section.\
18S_StrainProductivity.R - Generates figures that compare specific 18S strain relative abundance to sample productivity.\
18S_metadata.xlsx - Contains metadata used for 18S R files.
#### mothur
18S.sh – Executable file that sends the 18S batch file to the high performance computing center.\
18S.oligos – File that contains the 18S primers. Used by the 18S batch file.\
18S.batch – Classifies the 18S sequences using mothur and Silva.
### 23SDADA2
#### DADA2
23S_DADA2.R – Classifies the pilot and control RABR 23S sequences using DADA2.\
A23S_DADA2.R – Classifies all remaining biofilm 23S sequences using DADA2. This includes the lab scale, CVWRF, greenhouse RABRs, and the trickling filter.
#### RStudio
23S_D_abund.R – Generates relative abundance graphs of 23S sequences. Automatically returns figures for each taxonomic level (Phylum, Class, Order, Family, Genus).\
23S_D_alphaD.R – Generates alpha diversity and productivity graphs of 23S sequences.\
23S_D_betaD.R – Generates beta diversity graphs of 23S sequences.\
23S_D_lefse.R – Generates LEfSe graphs for 23S sequences. Has commented sections that contain code that must be run through mothur to proceed.\
23S_D_rarefaction.R – Generates rarefaction graphs of 23S sequences. Included a figure that is colored by section, as well as individual graphs for each section.\
23S_StrainProductivity.R - Generates figures that compare specific 23S strain relative abundance to sample productivity.\
23S_metadata.xlsx - Contains metadata used for 23S R files.
### ITSmothur
#### RStudio
ITS_abund.R – Generates relative abundance graphs of ITS. Automatically returns figures for each taxonomic level (Supergroup, Division, Subdivision, Class, Order, Family, Genus).\
ITS_alphaD.R – Generates alpha diversity and productivity graphs of ITS sequences.\
ITS_betaD.R – Generates beta diversity graphs of ITS sequences.\
ITS_lefse.R – Generates LEfSe graphs for ITS sequences. Has commented sections that contain code that must be run through mothur to proceed.\
ITS_rarefaction.R – Generates rarefaction graphs of ITS sequences. Included a figure that is colored by section, as well as individual graphs for each section.\
ITS_StrainProductivity.R - Generates figures that compare specific ITS strain relative abundance to sample productivity.\
ITS_metadata.xlsx - Contains metadata used for ITS R files.
#### mothur
ITS.sh – Executable file that sends the ITS_fungi and ITS_euk batch files to the high performance computing center.\
ITS.oligos – File that contains the ITS primers. Used by the ITS_fungi and ITS_euk batch files.\
ITS_fungi.batch – Classifies the ITS sequences using mothur and Silva and PR2.\
ITS_euk.batch – Classifies the ITS sequences using mothur and Silva and PR2.
