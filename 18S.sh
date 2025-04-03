#!/bin/bash 
#SBATCH --partition=kingspeak 
#SBATCH --account=miller 
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH -o 18S_output.txt 
#SBATCH -e 18S_error.txt 
#SBTACH --time=03:00:00 
#SBATCH --job-name=18Smatthews 
 
# load module 
mothur 18S.batch 
