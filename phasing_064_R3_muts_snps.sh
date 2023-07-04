#!/bin/bash
#SBATCH --job-name=hapcut_phasing  # Job name
#SBATCH --ntasks=1                   # Run one task
#SBATCH --mem=32G                   # Job Memory
#SBATCH --time=10:00:00               # Time limit hrs:min:sec
#SBATCH --out='%x.%j.out'    ### may need to modify if moving script
#SBATCH --error='%x.%j.err'  ### may need to modify if moving script
#SBATCH --array=1

REF="/camp/project/proj-tracerx-lung/tracerx/_PIPELINE/TRACERx-RNAseq-assets/v1/reference/hg38/hg38.fa"
VCF="/camp/project/proj-tracerx-lung/tctProjects/frankella/ctDNA_methods_dev/inputs/phasing/HaplotypeCaller_LTX0064_GL.vcf"
ONT_BAM="/nemo/lab/swantonc/working/ruizc/2022-07-14-dna_long_read_mutation_phasing/input/2022-07-27-mpileup_cohort/LTX0064_SU_N.sorted.bam"
FRAG_OUT="/camp/project/proj-tracerx-lung/tctProjects/frankella/ctDNA_methods_dev/outputs/mut_cn_tracking/20230418/LTX0064_SU_N_snps_muts_phased.frag.file" 
FINAL_OUT="/camp/project/proj-tracerx-lung/tctProjects/frankella/ctDNA_methods_dev/outputs/mut_cn_tracking/20230418/LTX0064_SU_N_snps_muts_phased.haplotype"
FINAL_2_OUT="/camp/project/proj-tracerx-lung/tctProjects/frankella/ctDNA_methods_dev/outputs/mut_cn_tracking/20230418/LTX0064_SU_N_snps_lowerthres_muts_phased.haplotype"

extractHAIRS --ONT 1 --ref ${REF} --indels 1 --realign_variants 1 --VCF ${VCF} --bam ${ONT_BAM}  --out ${FRAG_OUT}
HAPCUT2 --verbose 1 --fragments ${FRAG_OUT} --VCF ${VCF} --output ${FINAL_OUT}
HAPCUT2 --verbose 1 --threshold 0.5 --fragments ${FRAG_OUT} --VCF ${VCF} --output ${FINAL_2_OUT}






