#!/bin/bash

#SBATCH -J cytoSPACE
#SBATCH -p fat01
#SBATCH --mem=15G
#SBATCH -N 1
#SBATCH -n 5
#SBATCH -t 0-12:00
#SBATCH -o %j.out
#SBATCH -e %j.err

module load conda
source activate cytospace

sample=$1

inputDir="~/project/lung_10xST/data/cytoSPACE_input"
scDir=$inputDir/"SC_exp"
stDir=$inputDir/$sample"_ST"
outDir="~/project/lung_10xST/results/cytoSPACE_origQC"

mkdir $outDir/$sample"_cytospace"

cytospace \
    -sp $scDir/subscRNA_data.txt \
    -ctp $scDir/subcell_type_labels.txt \
    -stp $stDir/ST_data.txt \
    -cp $stDir/Coordinates.txt \
    -o $outDir/$sample"_cytospace" \
    -sm lap_CSPR 
   # --downsample-off
