#!/bin/bash
#Script to run SAMtools mpileup in subsampled bed files

#Libraries
source activate /nemo/lab/swantonc/working/ruizc/apps/anaconda/envs/nanopore_variant_calling

ml purge
module load parallel
module load R/4.0.0-foss-2020a
module load SAMtools/1.16.1-GCC-11.3.0

SCRATCH_DIR='outputs/phasing/20230627/scratch/'

#Select the split bed
ALL_BED_SPLIT=($(ls ${SCRATCH_DIR}/${1}.hetsnp.snvs.split* | cut -d '.' -f 4))
NB_BED=$((${SLURM_ARRAY_TASK_ID} - 1))

SPLIT=${ALL_BED_SPLIT[${NB_BED}]}

#Subset BAM for region of interest, take first-last position of the bed per chromosome
CHRS=($(cut -f 1 ${SCRATCH_DIR}/${1}.hetsnp.snvs.${SPLIT} | sort -u))

touch ${SCRATCH_DIR}/${1}.regions.${SPLIT}
for THIS_CHR in "${CHRS[@]}"
do
  STARTPOS=$(grep ${THIS_CHR} ${SCRATCH_DIR}/${1}.hetsnp.snvs.${SPLIT} | head -n 1 | cut -f 2)
  ENDPOS=$(grep ${THIS_CHR} ${SCRATCH_DIR}/${1}.hetsnp.snvs.${SPLIT} | tail -n 1 | cut -f 2)
  echo -e ${THIS_CHR}"\t"${STARTPOS}"\t"${ENDPOS} >> ${SCRATCH_DIR}/${1}.regions.${SPLIT}
done

#Subset bam
echo "Sarting subset for" ${2}
samtools view -bh -L ${SCRATCH_DIR}/${1}.regions.${SPLIT} ${2} > ${SCRATCH_DIR}/${1}.${SPLIT}.bam
echo "BAM subset finished"

#Subset reference
echo "Sarting subset reference for" ${SPLIT}
egrep -A 1 ${CHRS_PATTERN} ${3}
samtools faidx ${3} ${CHRS[@]} > ${SCRATCH_DIR}/${1}.${SPLIT}.fa
samtools faidx ${SCRATCH_DIR}/${1}.${SPLIT}.fa

echo "Mpileup run in sample" ${1}", split" ${SPLIT}
#Run mpileup, remove secondary reads (-B option removes the BAQ step, but vastly reduces memory usage)
# See https://github.com/samtools/bcftools/issues/1584
# See https://www.biostars.org/p/9466154/
samtools mpileup -f ${SCRATCH_DIR}/${1}.${SPLIT}.fa \
                 -l ${SCRATCH_DIR}/${1}.hetsnp.snvs.${SPLIT} \
                 -B \
                 --output-extra QNAME,MAPQ \
                 --excl-flags SECONDARY \
                 ${SCRATCH_DIR}/${1}.${SPLIT}.bam > ${SCRATCH_DIR}/${1}.pileup.${SPLIT}.txt

rm -f ${SCRATCH_DIR}/${1}.${SPLIT}.bam
rm -f  ${SCRATCH_DIR}/${1}.${SPLIT}.fa

