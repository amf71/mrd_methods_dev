#!/bin/bash

## Run Beagle to phase SNPs statistically rather than using the ONT BAM

## Libraries

# Following arg Inputs
# 1 = reference.fa
# 2 = ONT BAM
# 3 = VCF of events to phase together where possible
# 4 = Intermeidate output fragment file
# 5 = output haplotype file

extractHAIRS --ONT 1 --ref ${1} --indels 1 --realign_variants 1 --bam ${2} --VCF ${3} --out ${4}
HAPCUT2 --verbose 1 --fragments ${4} --VCF ${3} --output ${5}
HAPCUT2 --verbose 1 --threshold 0.5 --fragments ${4} --VCF ${3} --output ${5}.thresholded





