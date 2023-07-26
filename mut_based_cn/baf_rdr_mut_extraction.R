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
library(Rsamtools)

## Function to combine p values with fisher method from metap package (can't install..?) ##
sumlog <- function(p, log.p = FALSE) {
    
    p <- p[ !is.na(p) ] ## added line
    if(length(p)==1) return(list(p=p)) ## added line
    if(length(p)==0) return(list(p=NA)) ## added line

   keep <- (p > 0) & (p <= 1)
   invalid <- sum(1L * keep) < 2
   if(invalid) {
      warning("Must have at least two valid p values")
      res <- list(chisq = NA_real_, df = NA_integer_,
         p = NA_real_, validp = p[keep])
   } else {
      lnp <- log(p[keep])
      chisq <- (-2) * sum(lnp)
      df <- 2 * length(lnp)
      if(length(lnp) != length(p)) {
         warning("Some studies omitted")
      }
      res <- list(chisq = chisq, df = df,
         p = pchisq(chisq, df, lower.tail = FALSE,
            log.p = log.p), validp = p[keep])
    }
   class(res) <- c("sumlog", "metap")
   res
}


#############################################
#### Make a folder for this analysis run ####
#############################################

date <- gsub("-","",Sys.Date())

analysis_name <- 'mut_cn_processing'
out_name <- 'outputs'
out_dir_general <- paste(out_name, analysis_name, sep='/')
if( !file.exists(out_dir_general) ) dir.create( out_dir_general )

out_dir_logs <- paste(out_dir_general, 'logs', sep='/')
if( !file.exists(out_dir_logs) ) dir.create( out_dir_logs )

outputs.folder <- paste0( out_dir_general, "/", date, "/" )

if( !file.exists(outputs.folder) ) dir.create( outputs.folder )


##############################################
#### Get Inputs required for all analyses ####
##############################################

# muts called form short read vcf R3
muts_phas_path  <- "inputs/phasing/Mutect2_filtered_LTX0064_SU_T1-R3_vs_LTX0064_BS_GL.vcf"

# mutations phased using Nanopore from R3
phasing_path  <- "outputs/phasing/20230627/phased_mutations.tsv"

# HAPCUT2 output from phasing SNP with Nanopore Normal sample
phasing_vcf_path <- 'outputs/phasing/20230627/LTX0064_SU_N.snps_phased.frag.het.only.haplotype.phased.VCF'

# CN from phased region 3 / 6
cn_path  <- "inputs/LTX0064_R3_Battenberg/LTX0064_SU_T1-R3_subclones.txt"

# purity from regions (from battenberg)
purity_path <- "inputs/LTX0064_R3_Battenberg/LTX0064_SU_T1-R3_purity_ploidy.txt"

# het SNPs from phased tumour region 3 
snps_path <- 'outputs/phasing/20230627/scratch/LTX0064_SU_N_filt_hetsnps_mpileup.vcf'
bat_snps_dir <- 'inputs/LTX0064_R3_Battenberg/'
bat_snps_pattern <- 'LTX0064_SU_T1-R3_alleleFrequencies*'

## BAM file for R6 to extract mutations from R3
R6_bam_path <- 'inputs/LTX0064_R3_Battenberg/LTX0064_SU_T1-R6.recal.bam'
R6_index_path <- 'inputs/LTX0064_R3_Battenberg/LTX0064_SU_T1-R6.recal.bam.bai'

### Read in ###
muts_vcf  <-  vcfR::read.vcfR( muts_phas_path )
phas <- data.table::fread( phasing_path )
haplotypes_vcf <- read.vcfR( phasing_vcf_path )
purity_pliody <- data.table::fread( purity_path )
cn    <-  data.table::fread( cn_path )
cn_test    <-  data.table::fread( cn_test_path )
snps_vcf  <-  vcfR::read.vcfR( snps_path )

# collate tumour SNP read counts from battenberg
files <- list.files( bat_snps_dir, pattern = bat_snps_pattern, full.names = TRUE )
bat_snps <- rbindlist( lapply(files, fread ))


#############################################
#### Combine SNP, mut, phasing & cn data ####
#############################################

### annotate the battenberg SNPs baf ###
snps <- data.table(  chr =  gsub( 'chr', '', getCHROM(snps_vcf)),
                     pos =  getPOS(snps_vcf),
                     ref =  getREF(snps_vcf),
                     alt =  getALT(snps_vcf) )
setnames(bat_snps,  c('chr', 'pos', 'countA', 'countC', 'countG', 'countT', 'depth'))
bat_snps <- as.data.table( inner_join( bat_snps, snps ) )
bat_snps <- bat_snps[ nchar(ref) == 1 & nchar(alt) == 1 ]
bat_snps[, `:=`( var_count = get(paste0('count',alt)),  
                 ref_count = get(paste0('count',ref)) ), 
         1:nrow(bat_snps)]
bat_snps[, baf := var_count/(var_count + ref_count)]

# Extract full haplotype blocks
hap_blocks <- data.table( chr =      gsub( 'chr', '', getCHROM(haplotypes_vcf)),
                          pos =      getPOS(haplotypes_vcf),
                          ref =      getREF(haplotypes_vcf),
                          alt =      getALT(haplotypes_vcf),
                          block_id = as.character(extract.gt(haplotypes_vcf, element='PS')),
                          gt =       as.character(extract.gt(haplotypes_vcf, element='GT')) )

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
phased_muts[, seg_start := cn[ match, start ] ]
phased_muts[, seg_end := cn[ match, end ] ]

#####################################################
#### Match alleles and determine early mutations ####
#####################################################

## In cases where blocks are very long split them up to smaller bins ##
# Should allow finer CN resollution and also minimises
# prob of mutation getting incorrectly phased with another

#phased_muts[, block_id := block_id_orig ]
phased_muts[, bin_id := as.character(block_id) ]

#phased_muts[, block_nmuts := .N, by = block_id ]
#unique(phased_muts[, .(block_id,block_nmuts)])[, sum(block_nmuts > 15)]
#block <- 158536185
#nmuts <- phased_muts[ block_id_orig == block, .N ]

split_block <- function( block, nmuts, mut_length_lim = 15, mut_target_length = 10 ){
  if( nmuts > mut_length_lim ){
    n_new_blocks <- round( nmuts / mut_target_length )
    length_new_blocks <- ceiling( nmuts / n_new_blocks)
    excess <- length_new_blocks * n_new_blocks - nmuts
    length_new_blocks <- rep(length_new_blocks, n_new_blocks)
    if( excess > 0 ) length_new_blocks[ 1:excess ] <- length_new_blocks[ 1:excess ] - 1 
    suffixes <- unlist( lapply(1:n_new_blocks, function(i) rep(i, length_new_blocks[i])) )
    return( paste(block, suffixes, sep = '_'))
  } else {
    return( rep(block, nmuts) )
  }
}

phased_muts[ !is.na(block_id), 
             bin_id := split_block(unique(bin_id), .N, mut_length_lim = 15, mut_target_length = 10 ), 
             by = block_id ]



#phased_muts[ grepl( '_', bin_id), .(block_id, bin_id, pos)]
#phased_muts[ !block_id_orig == block_id, .(block_id, block_id_orig, pos)][1:20]
#phased_muts[ block_id_orig == 16675390, .(block_id, block_id_orig, pos)]


## copy new blocks onto baf data ##
#First work out where the edges of the blocks/new blocks are
# edges of old blocks from snp data
bat_snps[, pre_mut_pos := data.table::shift(pos, n = 1), by = chr ]
bat_snps[, post_mut_pos := data.table::shift(pos, n = -1), by = chr ]
bat_snps[, mid_before_pos := ifelse(!is.na(pre_mut_pos), round(pos - (pos - pre_mut_pos) / 2), pos) ]
bat_snps[, mid_after_pos  := ifelse(!is.na(post_mut_pos), round(pos + (post_mut_pos - pos) / 2), pos)  ]
bat_snps[, block_start := min(mid_before_pos), by = block ]
bat_snps[, block_end := min(mid_before_pos), by = block ]
bat_snps[, n_snps_block := , ]


# edges of get edges of bins from mut data
phased_muts[, pre_mut_pos := data.table::shift(pos, n = 1), by = chr ]
phased_muts[, post_mut_pos := data.table::shift(pos, n = -1), by = chr ]
phased_muts[, mid_before_pos := ifelse(!is.na(pre_mut_pos), round(pos - (pos - pre_mut_pos) / 2), pos) ]
phased_muts[, mid_after_pos  := ifelse(!is.na(post_mut_pos), round(pos + (post_mut_pos - pos) / 2), pos)  ]
phased_muts[, block_start := bat_snps[ match(block_id, block), block_start ]]
phased_muts[, block_end := bat_snps[ match(block_id, block), block_end ]]

block_sum <- unique( phased_muts[, .(mid_before_bin = min(mid_before_pos),
                             mid_after_bin = max(mid_after_pos),
                             block_start, block_end, block_id, 
                             bin_start = block_start,
                             bin_end = block_end,
                             block_pos = mean(pos), chr ),
                         by = bin_id] )

block_sum[ !bin_id == block_id, 
           `:=`( bin_start = ifelse(block_pos == min(block_pos), block_start, mid_before_bin ),
                 bin_end   = ifelse(block_pos == max(block_pos), block_end,   mid_after_bin  ) ),
           by = block_id ]

# overlay new bins on mutations and snps
phased_muts[, bin_start := block_sum[ match(phased_muts$bin_id, block_sum$bin_id), bin_start ]]
phased_muts[, bin_end := block_sum[ match(phased_muts$bin_id, block_sum$bin_id), bin_end ]]

bat_snps_granges <- copy(bat_snps)
setnames(bat_snps_granges, 'pos', 'start')
bat_snps_granges[, end := start ]
bat_snps_granges <- GenomicRanges::makeGRangesFromDataFrame(bat_snps_granges)
block_sum_granges <- copy(block_sum[, .(chr, start = bin_start, end = bin_end)])
block_sum_granges <- GenomicRanges::makeGRangesFromDataFrame(block_sum_granges, na.rm=TRUE)
match <- GenomicRanges::findOverlaps(bat_snps_granges, block_sum_granges, select = 'arbitrary')
bat_snps[, bin_start := block_sum[ match, bin_start ] ]
bat_snps[, bin_end   := block_sum[ match, bin_end ] ]
bat_snps[, bin_id   := block_sum[ match, bin_id ] ]

bat_snps[, table(is.na(bin_start))]

match <- GenomicRanges::findOverlaps(bat_snps_granges, cn_granges, select = 'arbitrary')
bat_snps[, AI := cn[ match, AI ] ]
bat_snps[, seg_name := cn[ match, paste(chr, start, end, sep = ':') ] ]



## Match the phase to the A/B/Maj/Min allele ##
## Function to extract BAF of each haplotype for each mutation
extract_hap_baf <- function( snps_by_chr, mut_chr, block, hap_by_chr, seg_start, seg_end ){ 

 #no y in Battenberg data
 if(mut_chr == 'Y') return( NA )

 # get SNPs & calculate average BAF for correct phase
 snps_chr_sub <- snps_by_chr[ names(snps_by_chr) == mut_chr ][[1]]
 hap_chr_sub <- hap_by_chr[ names(hap_by_chr) == mut_chr ][[1]]

 block_indices <- hap_chr_sub[, block_id == block ]
 block_positions <- hap_chr_sub[ block_indices, pos ]
 pos_seg_indices <- block_positions > seg_start & block_positions < seg_end
 block_phases <- hap_chr_sub[ block_indices, gt ]
 block_positions <- block_positions[ pos_seg_indices ]
 block_phases <- block_phases[ pos_seg_indices ]

 phas_baf <- snps_chr_sub[ match(block_positions, pos), 
                          mean( ifelse(block_phases == '0|1', baf, ifelse(block_phases == '1|0', 1-baf, NA)), na.rm=TRUE)  ]
 return( phas_baf )
}

# to speed up seperate by chrs at the start
snps_by_chr <- lapply( bat_snps[, unique(chr)], function(snpchr) bat_snps[ chr == snpchr ])
names(snps_by_chr) <- bat_snps[, unique(chr)]
hap_by_chr <- lapply( hap_blocks[, unique(chr)], function(snpchr) hap_blocks[ chr == snpchr ])
names(hap_by_chr) <- hap_blocks[, unique(chr)]

# run hapbaf calculation
phased_muts[, phas_baf := extract_hap_baf(snps_by_chr, 
                                          unique(chr),
                                          unique(block_id), 
                                          hap_by_chr,
                                          unique(seg_start),  
                                          unique(seg_end) ), 
            by = .(block_id, seg_start, seg_end)  ]

# add purity
phased_muts[, purity := purity_pliody$purity ]

# Calculate mut CPN
phased_muts[, vaf := varcount_tumour / (varcount_tumour + refcount_tumour)]
phased_muts[, mut_cpn := (vaf * 1 / purity) * ((purity * tot_cn) + 2 * (1 - purity))]

## Assign whether mutation is on all copies of an allele or not (early / late but done properly!) ie if mutCPN = CN of allele
# If there is imbalence which allele
phased_muts[, block_phas_baf := mean(phas_baf, na.rm=T), by = .(block_id, seg_start, seg_end) ]
phased_muts[, phased_allele_cn := ifelse( ifelse( phase == 'hap1', block_phas_baf, ifelse(phase == 'hap2', 1-block_phas_baf, NA)) > 0.5, maj_cn, min_cn) ]

# Parameter to limit to pretty clonal mutations
clonal_mut_cpn_lim <- 0.75
phased_muts[, early_mut := round(mut_cpn) >= phased_allele_cn & mut_cpn > clonal_mut_cpn_lim ]


############################################################
#### extract pure BAF/LogR from early mutations in bins ####
############################################################

# Do this first just from the same sample as phased

## GC, Rep timing and matched normal depth correction
# First extract GC content of local areas
library(BSgenome.Hsapiens.NCBI.GRCh38)

calc_surrounding_gc <- function( chr, pos, window, ref = Hsapiens ){

    if( grepl('chr', chr) ) chr <- gsub('chr', '', chr)
    seq <- BSgenome::getSeq(ref, chr, pos - window/2, pos + window/2)
    seq <- lapply(1:length(seq), function(i) strsplit(as.character(seq[[i]]), '')[[1]])
    gc  <- lapply(seq, function(bases) bases == 'G'| bases == 'C')  
    return( sapply(gc, function(x) sum(x)/length(x) ))

}


phased_muts[, gc := calc_surrounding_gc(unique(chr), pos, window = 300 ), by = chr ]

## also replication timing correction
# Code from ASCAT git repo
UWrepliseqbigwigs <- c("wgEncodeUwRepliSeqBg02esWaveSignalRep1.bigWig","wgEncodeUwRepliSeqBjWaveSignalRep1.bigWig",
                       "wgEncodeUwRepliSeqGm06990WaveSignalRep1.bigWig","wgEncodeUwRepliSeqGm12801WaveSignalRep1.bigWig",
                       "wgEncodeUwRepliSeqGm12812WaveSignalRep1.bigWig","wgEncodeUwRepliSeqGm12813WaveSignalRep1.bigWig",
                       "wgEncodeUwRepliSeqGm12878WaveSignalRep1.bigWig","wgEncodeUwRepliSeqHelas3WaveSignalRep1.bigWig",
                       "wgEncodeUwRepliSeqHepg2WaveSignalRep1.bigWig","wgEncodeUwRepliSeqHuvecWaveSignalRep1.bigWig",
                       "wgEncodeUwRepliSeqImr90WaveSignalRep1.bigWig","wgEncodeUwRepliSeqK562WaveSignalRep1.bigWig",
                       "wgEncodeUwRepliSeqMcf7WaveSignalRep1.bigWig","wgEncodeUwRepliSeqNhekWaveSignalRep1.bigWig",
                       "wgEncodeUwRepliSeqSknshWaveSignalRep1.bigWig")
UCSCdownloadpath <- "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/"
print("Downloading UW Repli-seq data from UCSC")
invisible(lapply(X = UWrepliseqbigwigs, FUN = function(x) download.file(url = paste0(UCSCdownloadpath, x), destfile = x)))
wavelets <- lapply(UWrepliseqbigwigs, rtracklayer::import.bw)

## match each SNP locus with the nearest replication timing estimate from each cell line
seqlevelsStyle(wavelets[[1]]) <- "Ensembl"
seqlevels(wavelets[[1]]) <- gsub('chr', '', seqlevels(wavelets[[1]]))
nearestidxs <- nearest(x = ph_granges, subject = wavelets[[1]], select = "arbitrary")

mcols(ph_granges) <- do.call(cbind, lapply(X = wavelets, FUN = function(x) mcols(x)[nearestidxs, "score"]))
colnames(mcols(ph_granges)) <- sub(pattern = "^wgEncodeUwRepliSeq", replacement = "", x = sub(pattern = "WaveSignalRep1.bigWig$", replacement = "", x = UWrepliseqbigwigs))

locidf <- as.data.frame(ph_granges)[, -c(3:5)]
locidf <- as.data.table( locidf )
setnames(locidf, 1:2,  c('chr', 'pos'))

phased_muts <- as.data.table( inner_join(phased_muts, locidf) )
cols <- names(locidf)[3:17]
correlation <- cor(phased_muts[!is.na(Hepg2), varcount_tumour],
           phased_muts[ !is.na(Hepg2), ..cols])
# should use the rep timing that best corrlates - probably doesn't make loads of diff
best <- colnames(correlation)[ correlation == min(correlation) ]

# calculate overall depth in match normal at each locus
phased_muts[, normal_depth := varcount_norm + refcount_norm ]

# Run model to account for these covariates - as per ascat
model = lm(varcount_tumour ~ splines::ns(x = gc, df = 5, intercept = T) + 
                             splines::ns(x = get(best), df = 5, intercept = T) + 
                             splines::ns(x = normal_depth, df = 5, intercept = T), y=F, model = F, data = phased_muts, na.action="na.exclude")
phased_muts[, varcount_tumour_corrected := residuals(model) + mean(varcount_tumour) ]
phased_muts[ varcount_tumour_corrected < 0, varcount_tumour_corrected := 0]

## Calcualate mutation based BAF and logR per block ##
# summarise haplotype counts and number of mutations per block
mut_baf_logR <- phased_muts[ (early_mut), 
                                .( mean_hap1_varcount = mean(varcount_tumour[ phase == 'hap1' ],na.rm=T),
                                   mean_hap2_varcount = mean(varcount_tumour[ phase == 'hap2' ],na.rm=T),
                                   mean_hap1_varnorm = mean((varcount_tumour/normal_depth)[ phase == 'hap1' ],na.rm=T),
                                   mean_hap2_varnorm = mean((varcount_tumour/normal_depth)[ phase == 'hap2' ],na.rm=T),
                                   mean_hap1_varcount_cor = mean(varcount_tumour_corrected[ phase == 'hap1' ],na.rm=T),
                                   mean_hap2_varcount_cor = mean(varcount_tumour_corrected[ phase == 'hap2' ],na.rm=T),
                                   n_hap1_muts = sum(phase == 'hap1'), n_hap2_muts = sum(phase == 'hap2'),
                                   block_start = min(pos),
                                   block_end = max(pos) ),
                             by = .(block_id, seg_start, seg_end, chr, maj_cn, min_cn, tot_cn) ]

# when you have LOH one allele will have no mutations. In this case set Varcount to 0. 
mut_baf_logR[ min_cn == 0, 
                `:=`(mean_hap1_varcount = ifelse( is.na(mean_hap1_varcount), 0, mean_hap1_varcount),
                     mean_hap2_varcount = ifelse( is.na(mean_hap2_varcount), 0, mean_hap2_varcount),
                     mean_hap1_varcount_cor = ifelse( is.na(mean_hap1_varcount_cor), 0, mean_hap1_varcount_cor),
                     mean_hap2_varcount_cor = ifelse( is.na(mean_hap2_varcount_cor), 0, mean_hap2_varcount_cor))]

# summarise baf/RDR across each block
mut_baf_logR <- mut_baf_logR[ !is.na(mean_hap1_varcount) & !is.na(mean_hap2_varcount), 
                                .( chr, block_start, block_end, seg_start, seg_end, block_id, maj_cn, min_cn, tot_cn, 
                                   n_hap1_muts, n_hap2_muts, 
                                   baf = mean_hap1_varcount / (mean_hap1_varcount + mean_hap2_varcount),
                                   rdr = mean_hap1_varcount + mean_hap2_varcount,
                                   baf_norm = mean_hap1_varnorm / (mean_hap1_varnorm + mean_hap2_varnorm),
                                   rdr_norm = mean_hap1_varnorm + mean_hap2_varnorm,
                                   baf_cor = mean_hap1_varcount_cor / (mean_hap1_varcount_cor + mean_hap2_varcount_cor),
                                   rdr_cor = mean_hap1_varcount_cor + mean_hap2_varcount_cor ) ]

# caluclate different ways of looking at baf (eg hatchet)
mut_baf_logR[, baf_inv := 1 - baf]
mut_baf_logR[, baf_consistent := abs(baf - 0.5)]
mut_baf_logR[, baf_norm_consistent := abs(baf_norm - 0.5)]
mut_baf_logR[, baf_cor_consistent := abs(baf_cor - 0.5)]

# write out 
fwrite(phased_muts, paste0(outputs.folder, date, '_phased_counts_by_mut_R3.tsv'), sep = '\t')
fwrite(mut_baf_logR, paste0(outputs.folder, date, '_baf_rdr_by_block_R3.tsv'), sep = '\t')

###########################################
#### Try applying to low purity region ####
###########################################

### First extract read counts for all mutations called in high purity region in low purity region
extr_muts <- phased_muts[ (early_mut), .(chr = paste0('chr', chr), start = pos, end = pos, ref, alt) ]
extr_muts_gr <- GenomicRanges::makeGRangesFromDataFrame(extr_muts)
sbp = ScanBamParam( which = extr_muts_gr )

pup = PileupParam(max_depth=10000,
                  min_base_quality=13,
                  min_mapq=0,
                  min_nucleotide_depth=1,
                  min_minor_allele_depth=0,
                  distinguish_strands=FALSE,
                  distinguish_nucleotides=TRUE,
                  ignore_query_Ns=TRUE,
                  include_deletions=FALSE,
                  include_insertions=FALSE)

read_counts <- as.data.table( pileup(file = R6_bam_path, index = R6_index_path, scanBamParam=sbp, pileupParam=pup) )

#above takes 20 minutes - so save
fst::write_fst(read_counts, paste0(outputs.folder, date, '_R6_early_phased_muts_piledup_in_R3.fst'))
# read_counts <- read_fst('outputs/mut_cn_processing/20230725/20230719_R6_early_phased_muts_piledup_in_R3.fst',as.data.table=T)

read_counts_gr <- copy(read_counts)
setnames(read_counts_gr, c('seqnames', 'pos'), c('chr', 'start'))
read_counts_gr[, end := start ]
read_counts_gr <- GenomicRanges::makeGRangesFromDataFrame(read_counts_gr)

match <- GenomicRanges::findOverlaps(read_counts_gr, extr_muts_gr, select = 'arbitrary')
read_counts[, ref := extr_muts[ match, ref ] ]
read_counts[, alt := extr_muts[ match, alt ] ]

# warning where varcount 0 as no row
read_counts <- unique(read_counts[, .(ref, alt, 
                                refcount_tumour = count[nucleotide == ref],
                                varcount_tumour = count[nucleotide == alt]), 
                            by = .(chr = gsub('chr', '', seqnames), pos) ])
read_counts[ is.na(refcount_tumour), refcount_tumour := 0 ]
read_counts[ is.na(varcount_tumour), varcount_tumour := 0 ]


test_data <- as.data.table( inner_join(read_counts, 
                                       phased_muts[ (early_mut), 
                                                    .(chr, pos, ref, alt, block_id, phase, 
                                                      normal_depth, gc, get(best), 
                                                      ref_seg_start = seg_start, ref_seg_end = seg_end,
                                                      loh_ref_sample = min_cn == 0 ) ]) )
setnames(test_data, 'V9', best)

model = lm(varcount_tumour ~ splines::ns(x = gc, df = 5, intercept = T) + 
                             splines::ns(x = get(best), df = 5, intercept = T) + 
                             splines::ns(x = normal_depth, df = 5, intercept = T), y=F, model = F, data = test_data, na.action="na.exclude")
test_data[, varcount_tumour_corrected := residuals(model) + mean(varcount_tumour) ]
test_data[ varcount_tumour_corrected < 0, varcount_tumour_corrected := 0]

quasar_test <- test_data[, .( ref_seg_start, ref_seg_end, chr, loh_ref_sample,
                             mean_hap1_varcount = mean(varcount_tumour[ phase == 'hap1' ], na.rm = TRUE),
                             mean_hap2_varcount = mean(varcount_tumour[ phase == 'hap2' ], na.rm = TRUE),
                             mean_hap1_varcount_cor = mean(varcount_tumour_corrected[ phase == 'hap1' ], na.rm = TRUE),
                             mean_hap2_varcount_cor = mean(varcount_tumour_corrected[ phase == 'hap2' ], na.rm = TRUE),
                             n_hap1_muts = sum(phase == 'hap1'), n_hap2_muts = sum(phase == 'hap2'),
                             n_hap1_muts_present = sum(phase == 'hap1' & varcount_tumour > 4), 
                             n_hap2_muts_present = sum(phase == 'hap2' & varcount_tumour > 4),
                             block_start = min(pos), block_end = max(pos) ),
                         by = block_id ]

# when you have LOH one allele will have no mutations. In this case set Varcount to 0. 
quasar_test[ (loh_ref_sample), 
                `:=`(mean_hap1_varcount = ifelse( is.na(mean_hap1_varcount), 0, mean_hap1_varcount),
                     mean_hap2_varcount = ifelse( is.na(mean_hap2_varcount), 0, mean_hap2_varcount),
                     mean_hap1_varcount_cor = ifelse( is.na(mean_hap1_varcount_cor), 0, mean_hap1_varcount_cor),
                     mean_hap2_varcount_cor = ifelse( is.na(mean_hap2_varcount_cor), 0, mean_hap2_varcount_cor))]

# Need to deal with clonal lllusion mutations - this case actually has a lot and messing up everything
# Going to have to presume all mutations at ~0 varcount are clonal illusion (will include new LOHs)
# Then re-infer new LOH by looking for blocks where all mutations are from one haplotype - only works for long segments with >7-8 muts
# Then look for refererence segments where all blocks of sufficent length have this propety and if so call all blocks where only 1 hap present as loh
# Think about resegmentation for this later!

# <5 varcount = absent for this purity fine 
quasar_test[, `:=`(mean_hap1_varcount = ifelse( mean_hap1_varcount < 5, NA, mean_hap1_varcount),
                   mean_hap2_varcount = ifelse( mean_hap2_varcount < 5, NA, mean_hap2_varcount),
                   mean_hap1_varcount_cor = ifelse( mean_hap1_varcount_cor < 5, NA, mean_hap1_varcount_cor),
                   mean_hap2_varcount_cor = ifelse( mean_hap2_varcount_cor < 5, NA, mean_hap2_varcount_cor))]

quasar_test[, `:=`( loh_test = binom.test( c( sum(!is.na(mean_hap1_varcount)), sum(!is.na(mean_hap2_varcount)) ) )$p.value,
                    nmuts_present = sum(!is.na(mean_hap1_varcount)) + sum(!is.na(mean_hap2_varcount)) ), 
            by = block_id ]

quasar_test[, sum_varcount := sum(mean_hap1_varcount, mean_hap2_varcount, na.rm=T), 1:nrow(quasar_test)]

# very unlikely to ever get all 10 mutations as clonal illusion
quasar_test[,`:=`(seg_loh_p = as.numeric(sumlog(loh_test[ nmuts_present > 4 ])$p),
                  is_mirror_vs_loh_ref = mean(sum_varcount) < 5 & (sum(n_hap1_muts) + sum(n_hap2_muts)) > 10), ## also look for cases with mirrored LOH
            by = .(chr, ref_seg_start, ref_seg_end)]

quasar_test[ (seg_loh_p < 0.0001 & n_hap1_muts > 0 & n_hap2_muts > 0) | is_mirror_vs_loh_ref, 
                `:=`(mean_hap1_varcount = ifelse( is.na(mean_hap1_varcount), 0, mean_hap1_varcount),
                     mean_hap2_varcount = ifelse( is.na(mean_hap2_varcount), 0, mean_hap2_varcount),
                     mean_hap1_varcount_cor = ifelse( is.na(mean_hap1_varcount_cor), 0, mean_hap1_varcount_cor),
                     mean_hap2_varcount_cor = ifelse( is.na(mean_hap2_varcount_cor), 0, mean_hap2_varcount_cor))]

quasar_test <- quasar_test[ !is.na(mean_hap1_varcount_cor) & !is.na(mean_hap1_varcount_cor), 
                                .( chr, block_start, block_end, ref_seg_start, ref_seg_end, block_id, loh_ref_sample, 
                                   seg_loh_p, is_mirror_vs_loh_ref, loh_test,
                                   n_hap1_muts, n_hap2_muts, nmuts = n_hap1_muts + n_hap2_muts,
                                   mean_hap1_varcount, mean_hap2_varcount, 
                                   nmuts_present = n_hap1_muts + n_hap2_muts,
                                   baf = mean_hap1_varcount / (mean_hap1_varcount + mean_hap2_varcount),
                                   rdr = mean_hap1_varcount + mean_hap2_varcount,
                                   baf_cor = mean_hap1_varcount_cor / (mean_hap1_varcount_cor + mean_hap2_varcount_cor),
                                   rdr_cor = mean_hap1_varcount_cor + mean_hap2_varcount_cor ) ]
quasar_test[, baf_inv := 1 - baf]
quasar_test[, baf_consistent := abs(baf - 0.5)]
quasar_test[, baf_cor_inv := 1 - baf_cor]
quasar_test[, baf_cor_consistent := abs(baf_cor - 0.5)]

# write out
fwrite(quasar_test, paste0(outputs.folder, date, '_quasar_baf_rdr_by_block_R6.tsv'), sep = '\t')
fwrite(test_data, paste0(outputs.folder, date, '_quasar_phased_counts_by_mut_R6.tsv'), sep = '\t')

#############
#### END ####
#############
