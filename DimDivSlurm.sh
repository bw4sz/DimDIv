#!/bin/bash

#SBATCH -J DimDiv
#SBATCH -o DimDiv.out
#SBATCH -p normal
#SBATCH -N 1 # Total nodes requested
#SBATCH -t 03:00:00  

#SBATCH --mail-user=benweinstein2010@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load R_mkl

R --no-save /home1/02443/bw4sz/DimDiv/DimDivCluster.R < out.R
