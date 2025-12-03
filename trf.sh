#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=2-00:00:00
#SBATCH --account=phd-ssalimi42
#SBATCH --partition=batch-impulse

spack load trf

input_1=$1 #specific region sequence like FoTeR element as fasta format
#output=$2 #output name and as gff
#input_2=$3 # full name of input_1


#trf $input_1 2 7 10 80 10 50 100 
trf $input_1.fasta 2 7 7 80 10 50 2000 -h

#~/.local/bin/trf2gff  -o ${input_1}_trf.gff < ${input_1}.fasta.2.7.7.80.10.50.2000.dat #set up this #.2.7.10.80.10.50.100.dat# based on trf's input
#Final_FoTeR_5176.fasta.2.7.10.80.10.50.100.dat

python3 trf2bed.py \
       --dat ${input_1}.fasta.2.7.7.80.10.50.2000.dat \
       --bed ${input_1}_trf.bed \
       --tool repeatseq
