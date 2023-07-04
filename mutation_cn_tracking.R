#!/usr/bin/env Rscript

#SBATCH --job-name=mut_phas
#SBATCH --time=10:00:00 
#SBATCH --mem=16G
#SBATCH --qos=normal
#SBATCH --out='outputs/mut_cn_tracking/logs/%j%t.out'    ### may need to modify if moving script
#SBATCH --error='outputs/mut_cn_tracking/logs/%j%t.err'  ### may need to modify if moving script

# print date in output file & error file so it can be matched back to output #

print( paste0(Sys.time(), "\n\n\n") )
warning( paste0(Sys.time(), "\n\n\n") )

#==============================================================================#
#==============================================================================#
######                                                                    ######
######  Script test whether you can exact good CN from Phased mutations   ######
######                                                                    ######
#==============================================================================#
#==============================================================================#

# Author: Alex Frankell
# Date: 2023-06-29
setwd("/camp/project/proj-tracerx-lung/tctProjects/frankella/ctDNA_methods_dev/")
setwd("/Volumes/proj-tracerx-lung/tctProjects/frankella/ctDNA_methods_dev/")

####################################################
#### Source required functions & load libraries ####
####################################################

# suppress warning on R build version #
library(fst)
library(data.table) 
library(vcfR) 
library(GenomicRanges) 
library(dplyr) 

#############################################
#### Make a folder for this analysis run ####
#############################################

date <- gsub("-","",Sys.Date())

analysis_name <- 'mut_cn_tracking'
out_name <- 'outputs'
out_dir_general <- paste(out_name, analysis_name, sep='/')
if( !file.exists(out_dir_general) ) file.create( out_dir_general )

out_dir_logs <- paste(out_dir_general, 'logs', sep='/')
if( !file.exists(out_dir_logs) ) file.create( out_dir_logs )

outputs.folder <- paste0( out_dir_general, "/", date, "/" )

if( !file.exists(outputs.folder) ) file.create( outputs.folder )


##############################################
#### Get Inputs required for all analyses ####
##############################################

# Phased muts short read vcf
muts_phas_path  <- "inputs/phasing/Mutect2_filtered_LTX0064_SU_T1-R3_vs_LTX0064_BS_GL.vcf"

# Short read Phasing
phasing_path  <- "outputs/phasing/20230627/phased_muts_long.tsv"

# CN from phased region
cn_path  <- "/camp/lab/swantonc/working/opich/shared/afrankell/wgs/LTX0064_SU_T1-R3_subclones.txt"
#cn_path  <- "/Volumes/lab-swantonc/working/opich/shared/afrankell/wgs/LTX0064_SU_T1-R3_subclones.txt"

# purity from pahsed region (from battenberg)
purity_path <- "/camp/lab/swantonc/working/opich/shared/afrankell/wgs/LTX0064_SU_T1-R3_purity_ploidy.txt"
#purity_path <- "/Volumes/lab-swantonc/working/opich/shared/afrankell/wgs/LTX0064_SU_T1-R3_purity_ploidy.txt"

#het SNPs from phased tumour region to get allele identities to match phasing
#snps_path <- 'outputs/phasing/LTX0064_haplotype.phased.VCF'
bat_snps_path <- '/camp/project/proj-tracerx-lung/tctProjects/frankella/ctDNA_methods_dev/inputs/LTX0064_R3_Battenberg/LTX0064_SU_T1-R3.BAFsegmented.txt'
#bat_snps_path <- '/Volumes/proj-tracerx-lung/tctProjects/frankella/ctDNA_methods_dev/inputs/LTX0064_R3_Battenberg/LTX0064_SU_T1-R3.BAFsegmented.txt'

#mutations from lower purity region of same tumour to test
muts_unphas_path <- '/camp/lab/swantonc/working/opich/shared/afrankell/wgs_temp/Mutect2_filtered_LTX0064_SU_T1-R6_vs_LTX0064_BS_GL.vcf.gz'
#muts_unphas_path <- '/Volumes/lab-swantonc/working/opich/shared/afrankell/wgs_temp/Mutect2_filtered_LTX0064_SU_T1-R6_vs_LTX0064_BS_GL.vcf.gz'

# Read in
muts_vcf  <-  vcfR::read.vcfR( muts_phas_path )
cn    <-  data.table::fread( cn_path )
bat_snps <- data.table::fread( snps_path )
setnames(bat_snps,  c('chr', 'pos', 'baf', 'baf_phased', 'baf_seg'))
muts_test_vcf <- vcfR::read.vcfR( muts_unphas_path )
phas <- data.table::fread( phasing_path )
purity_pliody <- data.table::fread( purity_path )

####################################
#### bring needed info together ####
####################################

### Get average allele specific cn accross all the clones in this sample ###

#remove unwanted columns 
cn <- cn[, 1:13]

# Melt by clone
setnames(cn, gsub('_A', '', names(cn)))
cn <- data.table::melt(cn, id.vars = c('chr', 'startpos', 'endpos'),
           measure = patterns('^nMaj', '^nMin', '^frac'),
           value.name = c("nMaj", "nMin", "frac"),
           variable.name = 'clone_id')

# remove NAs where fewer than max number of clones for a segment
cn <- cn[ !is.na(nMaj) ]

# average across clones 
cn <- cn[, .(nMaj = sum(nMaj * frac),
             nMin = sum(nMin * frac)),
         by = .(chr, startpos, endpos) ]
cn[, ntot := nMaj + nMin ]

## merge phasing status with short read VAFs ##
# Filter for PASS
muts_vcf <- muts_vcf[ getFILTER(muts_vcf) == 'PASS' ]
AD <- extract.gt(muts_vcf, element='AD')

muts <- data.table(  chr =  getCHROM(muts_vcf),
                     pos =  getPOS(muts_vcf),
                     ref =  getREF(muts_vcf),
                     alt =  getALT(muts_vcf),
                     varcount_tumour  =  as.numeric( tstrsplit(AD[, 2], split = ',')[[2]] ),
                     refcount_tumour  =  as.numeric( tstrsplit(AD[, 2], split = ',')[[1]] ), 
                     varcount_norm  =  as.numeric( tstrsplit(AD[, 1], split = ',')[[2]] ),
                     refcount_norm  =  as.numeric( tstrsplit(AD[, 1], split = ',')[[1]] ) )


phased_muts <- as.data.table( inner_join( muts, phas ) )

## merge in allele specific CN data ##
# get phsaed mut granges obj
phased_muts[, `:=`(start=pos, end=pos) ]
phased_muts[, chr := gsub('chr', '', chr) ]
ph_granges <- GenomicRanges::makeGRangesFromDataFrame(phased_muts)

#get phsed cn object
setnames(cn, c('startpos', 'endpos'), c('start', 'end'))
cn_granges <- GenomicRanges::makeGRangesFromDataFrame(cn)

# get indices for mutations in cn obj and then copy over information
match <- GenomicRanges::findOverlaps(ph_granges, cn_granges, select = 'arbitrary')
phased_muts[, maj_cn := cn[ match, nMaj ] ]
phased_muts[, min_cn := cn[ match, nMin ] ]
phased_muts[, tot_cn := cn[ match, ntot ] ]

# copy over the segment borders of the mutation so we resist to this area when calculating
# the hapbaf
phased_muts[, seg_end := cn[ match, end ] ]
phased_muts[, seg_start := cn[ match, start ] ]


#####################################################
#### Match alleles and determine early mutations ####
#####################################################

## Match the phase to the A/B/Maj/Min allele ##

## Function to extract BAF of each haplotype for each mutation
extract_hap_baf <- function( snps_by_chr, mut_chr, hap_id, SNP_positions ){

 # get SNPs & calculate average BAF for correct phase
 #no y in Battenberg data
 if(mut_chr == 'Y') return( NA )
 snps_chr_sub <- snps_by_chr[ names(snps_by_chr) == mut_chr ][[1]]
 hapbafs <- snps_chr_sub[ match(SNP_positions, pos), mean( ifelse(hap_id == 'A', baf, 1-baf), na.rm=TRUE)  ]
 
 return( hapbafs )
}

rowi = 3
mut_chr = phased_muts[rowi, chr]
hap_id = split_hapids[[rowi]]
SNP_positions = split_hapids_pos[[rowi]]


split_hapids_pos <- phased_muts[, tstrsplit( hap_id_SNP_pos, split = ':')]
split_hapids_pos <- lapply(1:nrow(split_hapids_pos), function(rowi){
    out <- as.numeric( na.omit( as.numeric(split_hapids_pos[rowi,])) )
    # limit to within segment borders
    out <- out[ out > phased_muts[ rowi, seg_start] & out < phased_muts[ rowi, seg_end] ]
    return(out)
})

split_hapids     <- phased_muts[, tstrsplit( hap1_id, split = ':')]
split_hapids <- lapply(1:nrow(split_hapids), function(rowi) as.character( na.omit( as.character(split_hapids[rowi,]) ) ) )

snps_by_chr <- lapply( bat_snps[, unique(chr)], function(snpchr) bat_snps[ chr == snpchr ])
names(snps_by_chr) <- bat_snps[, unique(chr)]


phased_muts[, hapbaf := sapply(1:.N, function(rowi){
    extract_hap_baf(snps_by_chr, phased_muts[rowi, chr], split_hapids[[rowi]], split_hapids_pos[[rowi]])
} ) ]

phased_muts[, AI := maj_cn/(maj_cn + min_cn) ]

# add purity
phased_muts[, purity := purity_pliody$purity ]

# Calculate mut CPN
phased_muts[, vaf := varcount_tumour / (varcount_tumour + refcount_tumour)]
phased_muts[, mut_cpn := (vaf * 1 / purity) * ((purity * tot_cn) + 2 * (1 - purity))]

## Assign whether mutation is on all copies of an allele or not (early / late but done properly!) ie if mutCPN = CN of allele
# If there is imbalence which allele
phased_muts[, block_hapbaf := mean(hapbaf, na.rm=T), by = block_id ]
phased_muts[, phased_allele_cn := ifelse( ifelse( phase == 0, block_hapbaf, 1 - block_hapbaf) > 0.5, maj_cn, min_cn) ]

phased_muts[, early_mut := round(mut_cpn) >= phased_allele_cn ]


############################################################
#### extract pure BAF/LogR from early mutations in bins ####
############################################################

# split genome to 5mb bins
# if not at least 2 mutations from each allele merge bins until this is the case for all
# Some analysis to look for any need for normalisation for coverage etc. Normalise by coverage in the normal sample (data in VCF)
# Get mutation based LogR and BAF average per bin (with CIs) & plot it


#############
#### END ####
#############
