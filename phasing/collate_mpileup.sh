#!/bin/bash

#Script to callate all mpileup outputs from each subdivided bed file

source activate /nemo/lab/swantonc/working/ruizc/apps/anaconda/envs/nanopore_variant_calling

ml purge
module load parallel
module load R/4.0.0-foss-2020a

SCRATCH_DIR='outputs/phasing/20230627/scratch/'

echo "Mpileup collate for sample "${1}

#Get all pileups together
cat ${SCRATCH_DIR}/${T_SAMPLE}.pileup.split*.txt > outputs/phasing/20230627/${T_SAMPLE}.pileup.txt

#Merge together mutation groups present in the same reads
Rscript scripts/phasing/mpileup_same_reads_snvs.R ${1}

#Remove temprary bam
rm -f ${SCRATCH_DIR}/${1}_filtered.bam*

