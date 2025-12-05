#!/bin/bash
#SBATCH --account=phd-ssalimi42
#SBATCH --partition=batch-impulse
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00

spack load singularityce
spack load samtools/f6zgopc cdhit  gffread bedops



input=$1 #genome name

gffread -x ${input}_cds.fasta -g $input.fasta $input.gff3
gff2bed < $input.gff3 > $input.bed



input=$1 #genome name


singularity run \
-B /usr/bin/which \
--bind /work \
--bind $(spack location -i ncbi-rmblastn):/rmblast \
--env PATH=/usr/local/bin:$PATH \
~/work/edta2.sif EDTA.pl \
--genome $input.fasta --cds ${input}_cds.fasta --exclude $input.bed \
--species others --overwrite 1 --sensitive 1 --anno 1 --threads ${SLURM_CPUS_PER_TASK} --force 1


#--curatedlib $PWD/EDTA/database/rice7.0.0.liban


awk 'BEGIN{OFS="\t"}
     !/^#/ {
        match($9,/classification=([^;]+)/,a);
        class=a[1];
        if (class=="") class="NA";
        print $1, $4-1, $5, $3, class, $6, $7
     }' ${input}.fasta.mod.EDTA.TEanno.gff3 > ${input}_edta.bed
