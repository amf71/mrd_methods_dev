#!/bin/bash
#SBATCH --job-name=mpileup_ont  # Job name
#SBATCH --ntasks=1                   # Run one task
#SBATCH --mem=4G                   # Job Memory
#SBATCH --time=00:30:00               # Time limit hrs:min:sec
#SBATCH --output=outputs/phasing/20230627/scratch/logs/mpileup_%j%t.log   # Standard output and error log

###        Script to phase mutations (from Illumina WGS) together and onto         ###
### parental alleles (with SNPs calledfrom Illumina WGS) using Nanopore sequencing ###

## Libraries
source activate /nemo/lab/swantonc/working/ruizc/apps/anaconda/envs/nanopore_variant_calling

ml purge
module load parallel
module load R/4.0.0-foss-2020a
module load BCFtools/1.9-foss-2018b

## Variables
# DIRs
INPUT_DIR='inputs/phasing/'
TUM_BAM_DIR='/nemo/lab/swantonc/working/ruizc/2022-07-14-dna_long_read_mutation_phasing/results/2022-07-27-mpileup_cohort/input'
OUT_DIR='outputs/phasing/20230627/'
SCRATCH_DIR='outputs/phasing/20230627/scratch/'

# BAMs
ONT_BAM_T="/nemo/project/proj-tracerx-lung/datasets/tracerx/sequence_data/nanopore/dna/2_minimap2Res/ds849_5_P1.sorted.bam"
ONT_BAM_N="/nemo/lab/swantonc/working/ruizc/2022-07-14-dna_long_read_mutation_phasing/input/2022-07-27-mpileup_cohort/LTX0064_SU_N.sorted.bam"

# VCFs
MUT_VCF="${INPUT_DIR}Mutect2_filtered_LTX0064_SU_T1-R3_vs_LTX0064_BS_GL.vcf"
SNP_VCF_N="${INPUT_DIR}HaplotypeCaller_LTX0064_GL.vcf"

# Sample names
T_SAMPLE="LTX0064_SU_T1-R3"
N_SAMPLE="LTX0064_SU_N"

# Ref files
REF="/nemo/project/proj-tracerx-lung/tracerx/_PIPELINE/TRACERx-Nanopore-assets/v1/reference/hg38/hg38.fa"

# Parameters
# Select the number of het SNPs around the mutation to extract with mpileup
# Please keep this even or I'll have to write unnecessary code to deal with this
NUMSNPS='6'

# Select the number of het SNPs around the mutation to exract with mpileup
NUM_PARALLEL='100'

## Filter / Select Mutations and SNPs we want to phase

#Generate a bed file from the VCF (keep only single nucleotide positions for now) for muts and snps
#Keep only het SNPs with Allele freq > 0.25 and <0.75
#bcftools view -v snps --samples ${T_SAMPLE} -i 'FILTER="PASS"' ${MUT_VCF} > ${SCRATCH_DIR}/${T_SAMPLE}_filt_snvs_mpileup.vcf
bcftools view --min-af 0.25 --max-af 0.75 ${SNP_VCF_N} > ${SCRATCH_DIR}/${N_SAMPLE}_filt_hetsnps_mpileup.vcf

# Run haplotype phasing on SNPs on normal ONT BAM to get parental alleles
jobid1=$(sbatch --parsable \
                --time=10:00:00 \
                --mem 32G \
                -J hapcut_phasing \
                --output=outputs/phasing/20230627/scratch/logs/hapcut_phasing_%j%t.log \
                scripts/phasing/hapcut_phasing.sh ${REF} \
                                                  ${SCRATCH_DIR}/${N_SAMPLE}_filt_hetsnps_mpileup.vcf \
                                                  ${ONT_BAM_N} \
                                                  ${OUT_DIR}${N_SAMPLE}.snps_phased.frag.het.only.file \
                                                  ${OUT_DIR}${N_SAMPLE}.snps_phased.frag.het.only.haplotype)



#Get the positions of PASS SNVs/Het SNPs
#grep -v '#' ${SCRATCH_DIR}/${T_SAMPLE}_filt_snvs_mpileup.vcf | awk -v OFS="\t" ' {print $1, $2}' > ${SCRATCH_DIR}/${T_SAMPLE}.snvs.bed
#grep -v '#' ${SCRATCH_DIR}/${N_SAMPLE}_filt_hetsnps_mpileup.vcf | awk -v OFS="\t" ' {print $1, $2}' > ${SCRATCH_DIR}/${N_SAMPLE}.hetsnps.bed
#
## Select only het SNPs around mutations (how many determined by 'NUMSNPS')
## 4th Arg = output 
## A better coder than me could do this in bash but this will only be like 2% of total run time anyway
#Rscript scripts/phasing/select_hetsnps.R \
#      ${SCRATCH_DIR}/${T_SAMPLE}.snvs.bed \
#      ${SCRATCH_DIR}/${N_SAMPLE}.hetsnps.bed \
#      ${NUMSNPS} \
#      ${SCRATCH_DIR}/${N_SAMPLE}.hetsnps.selected.bed
#
#
## Combine selected het SNPs and mutations together 
#cat ${SCRATCH_DIR}/${T_SAMPLE}.snvs.bed ${SCRATCH_DIR}/${N_SAMPLE}.hetsnps.selected.bed > ${SCRATCH_DIR}/${T_SAMPLE}.hetsnp.snvs.bed
#
### Run mpileup in Parallel for all selected SNPs + SNVs
#
##Subdivide the large bed file into multiple smaller bed files, ensure no lines are split
#split --number l/${NUM_PARALLEL} ${SCRATCH_DIR}/${T_SAMPLE}.hetsnp.snvs.bed ${SCRATCH_DIR}/${T_SAMPLE}.hetsnp.snvs.split -da 4
#
##ONT_BAM=$(ls ${BAM_DIR}/${T_SAMPLE}*bam)
#
## Run in parallel through all split bams
#ALL_BED_SPLIT=($(ls ${SCRATCH_DIR}${T_SAMPLE}.hetsnp.snvs.split* | cut -d '.' -f 4))
#
## Run mpileup in an array job for all split bed files
#jobid2=$(sbatch --parsable \
#                --time=10:00:00 \
#                --array=1-${#ALL_BED_SPLIT[@]} \
#                --mem 32G \
#                -J mpileup_split \
#                --output=outputs/phasing/20230627/scratch/logs/mpileup_split_%j%t.log \
#                scripts/phasing/mpileup_module.sh ${T_SAMPLE} ${ONT_BAM_T} ${REF})
#
### Collate all mpileup subfiles and run R script to generate phased blocks + assign parental alleles
##sbatch --parsable \
##       --dependency="afterok:$jobid1:$jobid2" \
##       --time=03-00:00:00 \
##       --mem 32G \
##       -J collate_mpileup \
##       --output=outputs/phasing/20230627/scratch/logs/collate_mpileup_%j%t.log \
##       scripts/phasing/collate_mpileup_generate_blocks.sh ${T_SAMPLE} ${N_SAMPLE}
#
### END ###