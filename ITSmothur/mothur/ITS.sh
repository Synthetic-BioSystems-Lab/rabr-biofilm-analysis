#!/bin/bash 
#SBATCH --partition=kingspeak 
#SBATCH --account=miller 
#SBATCH --nodes=4 
#SBATCH --ntasks=4 
#SBATCH -o ITS_output.txt 
#SBATCH -e ITS_error.txt 
#SBTACH --time=05:00:00 
#SBATCH --job-name=ITSmatthews 
 
# fungi database 
mothur ITS_fungi.batch 
 
# eukaryotic database 
mothur ITS_euk.batch 
