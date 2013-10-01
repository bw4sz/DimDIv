#!/bin/bash

#SBATCH -J DimDiv2
#SBATCH -o /home1/02443/bw4sz/DimDiv/DimDiv123.out
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH -A TG-TRA120007   
#SBATCH -n 100 # Total cores

#SBATCH --mail-user=benweinstein2010@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load R_mkl

source /work/01125/yye00/ParallelR/sourceme.sh


echo "Dim Div 1"

ibrun RMPISNOW < /home1/02443/bw4sz/DimDiv/DimDivCluster3.R > /home1/02443/bw4sz/DimDiv/DimDiv1.out


echo "Dim Div 2"

ibrun RMPISNOW < /home1/02443/bw4sz/DimDiv/DimDivCluster3.R > /home1/02443/bw4sz/DimDiv/DimDiv2.out



echo "Dim Div 3"

ibrun RMPISNOW < /home1/02443/bw4sz/DimDiv/DimDivCluster3.R > /home1/02443/bw4sz/DimDiv/DimDiv3.out

echo "done"

