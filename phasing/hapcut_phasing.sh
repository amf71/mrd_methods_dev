#!/bin/bash

## Run Hapcut2 to phase events on long read ONT BAM

## Libraries
module load Anaconda3/2023.03
source activate /nemo/lab/swantonc/working/ruizc/apps/anaconda/envs/nanopore_variant_calling

# Following arg Inputs
# 1 = reference.fa
# 2 = VCF of events to phase together where possible
# 3 = ONT BAM
# 4 = Intermeidate output fragment file
# 5 = output haplotype file

extractHAIRS --ONT 1 --ref ${1} --indels 1 --realign_variants 1 --VCF ${2} --bam ${3} --out ${4}
HAPCUT2 --verbose 1 --threshold 1 --fragments ${4} --VCF ${2} --output ${5}





