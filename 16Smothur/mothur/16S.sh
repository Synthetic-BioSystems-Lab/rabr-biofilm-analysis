#!/bin/bash 
#SBATCH --partition=kingspeak 
#SBATCH --account=miller 
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH -o 16S_output.txt 
#SBATCH -e 16S_error.txt 
#SBTACH --time=03:00:00: 
#SBATCH --job-name=16Smatthews 
 
mothur 16S.batch 
