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
library(biomaRt)
library(ggplot2)

# my across genome plotting functions
source('../repos/my_R_packages/amfFunctions/R/FrankellA.functions.R')

#############################################
#### Make a folder for this analysis run ####
#############################################

date <- gsub("-","",Sys.Date())

analysis_name <- 'mut_cn_tracking'
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

# Phased muts short read vcf
muts_phas_path  <- "inputs/phasing/Mutect2_filtered_LTX0064_SU_T1-R3_vs_LTX0064_BS_GL.vcf"

# data from baf_rdr_mut_extraction.R processing script
phased_muts_path <- 'outputs/mut_cn_processing/20230726/20230726_phased_counts_by_mut_R3.tsv'
mut_baf_logR_path <- "outputs/mut_cn_processing/20230726/20230726_baf_rdr_by_block_R3.tsv"
test_data_path <- 'outputs/mut_cn_processing/20230726/20230726_quasar_phased_counts_by_mut_R6.tsv'
quasar_test_path <- 'outputs/mut_cn_processing/20230726/20230726_quasar_baf_rdr_by_block_R6.tsv'

# CN from phased region
cn_path  <- "inputs/LTX0064_R3_Battenberg/LTX0064_SU_T1-R3_subclones.txt"
cn_test_path  <- "inputs/LTX0064_R6_Battenberg//LTX0064_SU_T1-R6_subclones.txt"

#het SNPs from phased tumour region to get allele identities to match phasing
snps_path <- 'outputs/phasing/20230627/scratch/LTX0064_SU_N_filt_hetsnps_mpileup.vcf'
bat_snps_path <- 'outputs/mut_cn_processing/20230725/20230725_bat_snps.fst'

# mutation phasing data
phasing_long_path  <- "outputs/phasing/20230627/phased_muts_long.tsv"

# HAPCUT2 output
hapcut_blocks_path <- 'outputs/mut_cn_processing/20230725/20230725_hapcut2_full_blocks.fst'


### Read in ###

# Phased mutations data processed R3
phased_muts <- fread( phased_muts_path )

# Data processed by block baf/logR R3
mut_baf_logR <- fread( mut_baf_logR_path )

# Phased mutations extracted R6
test_data <- fread( test_data_path )

# Data processed by block baf/logR R6 quasar
quasar_test <- fread( quasar_test_path )

# cn data
cn    <-  data.table::fread( cn_path )
cn_test    <-  data.table::fread( cn_test_path )

# Long format phasing output
phasing_long  <- fread( phasing_long_path )

# tumour R3 SNP read counts / bafs from battenberg
bat_snps <- fst::read_fst( bat_snps_path, as.data.table=T )

# R3 muts vcf
muts_vcf  <-  vcfR::read.vcfR( muts_phas_path )

# full phasing data
hap_blocks <- fst::read_fst( hapcut_blocks_path, as.data.table=T )


###########################################################################
#### Make some QC plots on which mutations aren't phasable and why not ####
###########################################################################

phased_muts[, hapdepth_plot := ifelse(hapdepth > 50, 50, hapdepth )] ## phas to phased_muts

pdf( paste0( outputs.folder, date, '_hapdepth_vs_hapvaf.pdf') )
ggplot(phased_muts, aes(x = hapdepth_plot, y = hap1_vaf )) +
    geom_density_2d_filled( bins = 20 ) +
    geom_vline( xintercept = 2.5 ) +
    geom_hline( yintercept = 0.25 ) +
    geom_hline( yintercept = 0.75 ) +
    scale_x_continuous( breaks = seq(0, 50, 5), labels = c(seq(0, 45, 5), '>50')) +
    #geom_text( x = 25, y = 0.5, size = 4, colour = 'white',
    #           label = '73% of mutations >=3 hapdepth & >0.75 or <0.25 hapAF') +
    theme_classic() +
    labs( x = 'Number of hetSNP + Mut overlapping ONT reads',
          y = 'Mean haplotype AF of hetSNP + Mutation reads') +
    theme( text = element_text( size = 20),
           legend.position = "none" )
dev.off()

# Above doesn't include mutations with no reads overlapping a hetSNP
plot_data <- phasing_long[(gt == '0|1' | gt == '1|0'), 
                            .(n_hetsnps = .N, mean_dist = mean(dist_mut_to_snp), 
                              hapdepth = sum(mut_varcount), total_depth = sum(shared_depth),
                              hap1_vaf), 
                          by = .(chr, pos) ]
plot_data[, hapdepth_plot := ifelse(hapdepth > 50, 50, hapdepth )]

pdf( paste0( outputs.folder, date, '_hapdepth_hist.pdf') )
ggplot(plot_data, aes(x = hapdepth_plot)) +
    geom_histogram( bins = 100 ) +
    scale_x_continuous( breaks = seq(0, 50, 5), labels = c(seq(0, 45, 5), '>50') ) +
    theme_classic() +
    geom_vline( xintercept = 2.5 ) +
    geom_text( x = 25, y = 5000, size = 4, 
               label = '81% mutations > 0 hetSNP overlapping reads\n76% at least 3 overlapping reads (phasing threshold)') +
    labs( x = 'Number of hetSNP + Mut overlapping ONT reads',
          y = 'Freq') +
    theme( text = element_text( size = 20),
           legend.position = "none" )
dev.off()

## Check if the mutations we can't phase are subclonal / low VAF
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
rm(muts_vcf)

plot_data <- as.data.table( inner_join( muts[, .( chr = gsub('chr', '', chr), pos, 
                                                  vaf = varcount_tumour / (refcount_tumour + varcount_tumour) )], 
                                        plot_data[, .(chr, pos, hapdepth, hap1_vaf) ] ) )

pdf( paste0( outputs.folder, date, '_shortread_vaf_vs_0hapdepth.pdf') )
ggplot(plot_data, aes(x = hapdepth == 0, y = vaf )) +
    geom_boxplot() +
    theme_classic() +
    labs( y = 'VAF in Illumina WGS') +
    theme( text = element_text( size = 20),
           legend.position = "none" )
dev.off()

plot_data <- unique(plot_data[, .(chr, pos, vaf, hap1_vaf, hapdepth)])
plot_data[ vaf>0.20 & !is.na(hapdepth), table( (hap1_vaf < 0.25 | hap1_vaf > 0.75) &  hapdepth > 2) ]


#########################################################
####   check to see how good phasing is using areas  ####
####        of imbalence in the tumour samples       ####
#########################################################

# get to fractional cn from bat output
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

#copy over AI for each segment to Bat SNPs to pull out ares of imbalence and thier haplotypes
cn[, AI := nMaj/(nMin+nMaj) ]
setnames(cn, c('startpos', 'endpos'), c('start', 'end'))
cn_granges <- GenomicRanges::makeGRangesFromDataFrame(cn)
bat_snps_granges <- copy(bat_snps)
setnames(bat_snps_granges, 'pos', 'start')
bat_snps_granges[, end := start ]
bat_snps_granges <- GenomicRanges::makeGRangesFromDataFrame(bat_snps_granges)
match <- GenomicRanges::findOverlaps(bat_snps_granges, cn_granges, select = 'arbitrary')
bat_snps[, AI := cn[ match, AI ] ]
bat_snps[, seg_name := cn[ match, paste(chr, start, end, sep = ':') ] ]

# just look at high imbalence and SNPs that are very clearly on different alleles
bat_snps[ AI > 0.6 , phase := ifelse(baf > 0.58, 'hap1', ifelse( baf < 0.4, 'hap2', NA ))]

# Overlap hapcut phasing
hap_blocks_granges <- copy(hap_blocks)
setnames(hap_blocks_granges, 'pos', 'start')
hap_blocks_granges[, end := start ]
hap_blocks_granges <- GenomicRanges::makeGRangesFromDataFrame(hap_blocks_granges)
match <- GenomicRanges::findOverlaps(bat_snps_granges, hap_blocks_granges, select = 'arbitrary')
hap_blocks[, block_length := max(pos, na.rm=T) - min(pos, na.rm=T), by = block_id ]
bat_snps[, block := hap_blocks[ match, block_id ] ]
bat_snps[, block_length := hap_blocks[ match, block_length ] ]
bat_snps[, gt := hap_blocks[ match, gt ] ]
bat_snps[, N_seg := length(unique(seg_name)), by = block ]

# Overlap bins from mut phasing
phased_muts_granges <- GenomicRanges::makeGRangesFromDataFrame(phased_muts[, .(chr, start = bin_start, end = bin_end)])
match <- GenomicRanges::findOverlaps(bat_snps_granges, phased_muts_granges, select = 'arbitrary')
bat_snps[, bin_id := phased_muts[ match, bin_id ] ]


# Overlay onto phasing whether there is consistency with BAF seperation based 
# phasing in each block

good_phas_fun <- function( phase, gt, n_snps_limit = 5 ){
    
    if(all(is.na(phase))) return( NA )
    if(n_snps_limit > min(table(phase)) - 1) n_snps_limit <- min(table(phase)) - 1

    phas1_correct <- table( phase[ gt == '1|0' ] )
    phas1_correct <- phas1_correct[ phas1_correct > n_snps_limit ]
    phas1_correct <- length(phas1_correct) < 2
    phas2_correct <- table( phase[ gt == '0|1' ] )
    phas2_correct <- phas2_correct[ phas2_correct > n_snps_limit ]
    phas2_correct <- length(phas2_correct) < 2

    return( phas2_correct & phas1_correct )

}

bat_snps[ !is.na(bin_id), quality_bin := good_phas_fun(phase, gt, n_snps_limit = 0), by = bin_id  ]
phased_muts[, quality_bin := bat_snps[ match(phased_muts$bin_id, bat_snps$bin_id), quality_bin ]]
mut_baf_logR[, quality_bin := phased_muts[ match(mut_baf_logR$bin_id, phased_muts$bin_id), quality_bin ]]


phased_muts[, poss_misphased := phased_allele_cn - mut_cpn < -0.5 ]
phased_muts[, table(poss_misphased, quality_bin)]
sum <- phased_muts[ !phase=='', .(nmut = .N, block_length = max(pos) - min(pos)), by = block_id ]
sum[, gn_frac := block_length / 3000000000 ]
sum[, sum(gn_frac) * 100, by = nmut > 6 ][ order(nmut) ]
# need to get these numbers up with longer phasing (+ also more accurate!)


############################################################
#### Plot RDR / BAF from mutations in Reference region  ####
############################################################

# add cummulative for plotting across chromosomes

mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart)
a <- getBM(c("end_position","chromosome_name"),"chromosome_name",c(1:22,"X"), mart=mart, useCache = FALSE)
chrlength.human <- sapply( unique(a$chromosome_name), function(x) max(a[a$chromosome_name ==x, "end_position" ]))
names(chrlength.human)[23] <- 23
chrlength.human <- chrlength.human[order(as.numeric(names(chrlength.human)))]
chrlength.human_before_cum <- sapply(1:length(chrlength.human), function(i){
    if( i == 1 ) return( 0 )
    return( sum(chrlength.human[1:(i-1)]) )
})
names(chrlength.human_before_cum) <- names(chrlength.human)

setnames(mut_baf_logR, c('seg_start', 'seg_end'), c('start', 'end'))
mut_baf_logR[, `:=`(Start.cummulitive = start + chrlength.human_before_cum[ names(chrlength.human_before_cum) == chr ], 
                    End.cummulitive = end + chrlength.human_before_cum[ names(chrlength.human_before_cum) == chr ]),
               1:nrow(mut_baf_logR) ]

plot_cum <- make.axis.adjust.cummulative(chrlength=chrlength.human)
plot_cum$chr <- factor(plot_cum$chr, levels = 1:23 )
plot_cum <- plot_cum[order( as.numeric(plot_cum$chr) ),]

### Plot just the CN data from Battenberg

plot_cn_data <- mut_baf_logR[, .(Start.cummulitive, End.cummulitive, maj_cn, min_cn, tot_cn) ]
plot_cn_data <- melt(plot_cn_data, id.vars = c('Start.cummulitive', 'End.cummulitive'), 
                     measure = patterns( '*_cn$' ),
                     variable.name = 'allele',
                     value.name = 'cn' )
plot_cn_data[, cn := as.numeric(cn) ]
plot_cn_data[, cn := ifelse( cn > 5, 5, cn) ]

plot_cn_data[, cn_plot := ifelse(allele == 'maj_cn', cn + 0.01, ifelse( allele == 'min_cn',  cn - 0.01, cn))]

pdf(paste0(outputs.folder, date, '_R3_ascn.pdf'), width = 15, height = 5)
ggplot() +
    geom_rect(data = plot_cn_data[ !allele == 'tot_cn' ], aes(xmin = Start.cummulitive, xmax = End.cummulitive , ymin = cn_plot, ymax = cn_plot, colour = allele)) +
    scale_x_continuous(label = plot_cum$chr, breaks = plot_cum$center) +
    scale_y_continuous( breaks = 0:5, labels = c(0:4, '>5'))+
    scale_colour_brewer( palette = 'Set1' ) +
    labs( x = 'Chromosome', 
          y = 'cn' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

pdf(paste0(outputs.folder, date, 'totcn.pdf'), width = 15, height = 5)
    ggplot() +
    geom_rect(data = plot_cn_data[ allele == 'tot_cn' ], 
               aes(xmin = Start.cummulitive, xmax = End.cummulitive , ymin = cn, ymax = cn), colour = 'black') +
    scale_x_continuous(label = plot_cum$chr, breaks = plot_cum$center) +
    scale_y_continuous( breaks = 0:5, labels = c(0:4, '>5'), limits = c(0,5))+
    labs( x = 'Chromosome', 
          y = 'cn' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()


# Now plot the RDR and BAF across the genome

mut_baf_logR[, block_centre := (block_end - block_start)/2 + block_start ]
mut_baf_logR[, bin_centre_cum := block_centre + chrlength.human_before_cum[ names(chrlength.human_before_cum) == chr ], 
               1:nrow(mut_baf_logR) ]

# rdr & baf across genome
mut_baf_logR[, plot_rdr := ifelse( rdr > 1200, 1200, rdr) ]
plot_data <- mut_baf_logR[ , .(mean_rdr = mean(plot_rdr, na.rm=T)), by = .(Start.cummulitive, End.cummulitive)]

pdf(paste0(outputs.folder, date, '_rdr_accross_genome.pdf'), width = 15, height = 5)
    ggplot() +
    geom_point(data = mut_baf_logR, shape = 21,
                aes(x = bin_centre_cum, y = plot_rdr, fill = baf_consistent)) +
    geom_rect(data = plot_data, aes(xmin = Start.cummulitive, 
                                    xmax = End.cummulitive , 
                                    ymin = mean_rdr, 
                                    ymax = mean_rdr), 
              colour = 'red', size = 0.5) +
    scale_fill_distiller( palette = 'YlGnBu', direction = 1 ) +
    scale_x_continuous(label = plot_cum$chr, breaks = plot_cum$center) +
    scale_y_continuous(breaks = seq(0, 1200, 200), labels = c(seq(0, 1000, 200), '>1200'), limits =) +
    labs( x = 'Chromosome', 
          y = 'Sum mean SNV varcount across alleles',
          fill = 'hapBAF' ) +
    theme_bw() +
    theme( text = element_text( size = 16 ))
dev.off()

phased_muts[, pos_cum := pos + chrlength.human_before_cum[ names(chrlength.human_before_cum) == chr ], 
             1: nrow(phased_muts)]

pdf(paste0(outputs.folder, date, 'raw_varcounts_across_genome_phased.pdf'), width = 15, height = 5)
    ggplot(phased_muts[ !phase == '' & (early_mut) ]) +
    geom_point(aes(x = pos_cum, y = varcount_tumour, colour = phase ), alpha = 0.5 ) +
    scale_x_continuous(label = plot_cum$chr, breaks = plot_cum$center) +
    labs( x = 'Chromosome', 
          y = 'Mutation Varcount' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()


## make a hatchet style plot
pdf(paste0(outputs.folder, date, '_hatchet.pdf'), width = 9)
    ggplot(mut_baf_logR) +
    geom_point(aes(x = baf_consistent, y = log2(rdr), colour = chr ), alpha = 0.5 ) +
    labs( x = '0.5-BAF', 
          y = 'log2 sum allele readcounts' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

#look at hist of LogR with LOH 
pdf(paste0(outputs.folder, date, '_loh_BAF_hist.pdf'), width = 9)
    ggplot(mut_baf_logR[ baf_consistent == 0.5]) +
    geom_histogram(aes(x = log2(rdr), fill = chr ) ) +
    labs( x = 'log2 sum allele readcounts' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

#look at hist of all LogR  
pdf(paste0(outputs.folder, date, '_BAF_hist.pdf'), width = 9)
    ggplot(mut_baf_logR) +
    geom_histogram(aes(x = baf_consistent, fill = chr ) ) +
    labs( x = '0.5-BAF' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

#look at hist of all BAF  
pdf(paste0(outputs.folder, date, '_rdr_hist.pdf'), width = 9)
    ggplot(mut_baf_logR) +
    geom_histogram(aes(x = log2(rdr), fill = chr ) ) +
    labs( x = 'log2 sum allele readcounts' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()


#look at hist of all LogR  -qual
pdf(paste0(outputs.folder, date, '_BAF_hist_qual.pdf'), width = 9)
    ggplot(mut_baf_logR[ (quality_bin) ]) +
    geom_histogram(aes(x = baf_consistent, fill = chr ) ) +
    labs( x = '0.5-BAF' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

#look at hist of all BAF  -qual
pdf(paste0(outputs.folder, date, '_rdr_hist_qual.pdf'), width = 9)
    ggplot(mut_baf_logR[ (quality_bin) ]) +
    geom_histogram(aes(x = log2(rdr), fill = chr ) ) +
    labs( x = 'log2 sum allele readcounts' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

## make a hatchet style plot
pdf(paste0(outputs.folder, date, '_hatchet_quality2.pdf'), width = 9)
    ggplot(mut_baf_logR[ (quality_bin) & n_hap1_muts > 5 & n_hap2_muts > 5 ]) +
    geom_point(aes(x = baf_consistent, y = log2(rdr), colour = chr ), alpha = 0.5 ) +
    labs( x = '0.5-BAF', 
          y = 'log2 sum allele readcounts' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

# compare baf and AI in quality bins
pdf(paste0(outputs.folder, date, '_rdr_vs_totcn_quality.pdf'), width = 9)
    mut_baf_logR[ (quality_bin), plot(log2(rdr), tot_cn) ]
dev.off()

# compare baf and AI in quality bins
pdf(paste0(outputs.folder, date, '_baf_vs_AI_quality.pdf'), width = 9)
    mut_baf_logR[ (quality_bin) & n_hap1_muts > 1 & n_hap2_muts > 1, plot(AI, baf_consistent) ]
dev.off()

# drdr across genome for quality bins
plot_data <- mut_baf_logR[ (quality_bin) , .(mean_rdr = mean(rdr, na.rm=T)), by = .(Start.cummulitive, End.cummulitive)]

pdf(paste0(outputs.folder, date, '_rdr_accross_genome_quality.pdf'), width = 15, height = 5)
    ggplot() +
    geom_point(data = mut_baf_logR[(quality_bin)], shape = 21,
                aes(x = bin_centre_cum, y = rdr, fill = baf_consistent)) +
    geom_rect(data = plot_data, aes(xmin = Start.cummulitive, 
                                    xmax = End.cummulitive , 
                                    ymin = mean_rdr, 
                                    ymax = mean_rdr), 
              colour = 'red', size = 0.5) +
    scale_fill_distiller( palette = 'YlGnBu', direction = 1 ) +
    scale_x_continuous(label = plot_cum$chr, breaks = plot_cum$center) +
    labs( x = 'Chromosome', 
          y = 'Sum mean SNV varcount across alleles',
          fill = 'hapBAF' ) +
    theme_bw() +
    theme( text = element_text( size = 16 ))
dev.off()

# restrict to blocks with more mutations
mut_baf_logR[, nmuts := n_hap1_muts + n_hap2_muts]
mut_baf_logR[, nmuts_perAl := ifelse( min_cn == 0, nmuts, nmuts/2 )]

plot_data <- mut_baf_logR[ (quality_bin) & nmuts_perAl > 2 , .(mean_rdr = mean(rdr, na.rm=T)), by = .(Start.cummulitive, End.cummulitive)]

pdf(paste0(outputs.folder, date, '_rdr_accross_genome_quality_3muts_per_al.pdf'), width = 15, height = 5)
    ggplot() +
    geom_point(data = mut_baf_logR[(quality_bin)  & nmuts_perAl > 2], shape = 21,
                aes(x = bin_centre_cum, y = rdr, fill = baf_consistent)) +
    geom_rect(data = plot_data, aes(xmin = Start.cummulitive, 
                                    xmax = End.cummulitive , 
                                    ymin = mean_rdr, 
                                    ymax = mean_rdr), 
              colour = 'red', size = 0.5) +
    scale_fill_distiller( palette = 'YlGnBu', direction = 1 ) +
    scale_x_continuous(label = plot_cum$chr, breaks = plot_cum$center) +
    labs( x = 'Chromosome', 
          y = 'Sum mean SNV varcount across alleles',
          fill = 'hapBAF' ) +
    theme_bw() +
    theme( text = element_text( size = 16 ))
dev.off()
# doesn't really seem to decrease noise

# quantity variability in RDR/BAF within blocks in the same CN segment (ie noise)
mut_baf_logR[, nmuts_perAl := ifelse( min_cn == 0, nmuts, nmuts/2 )]
mut_baf_logR[, AI := maj_cn / (maj_cn + min_cn) ]

seg_sum <- mut_baf_logR[ n_hap1_muts > 1 & n_hap2_muts > 1, 
                         .(AI, tot_cn, mean_rdr = mean(rdr), 
                          norm_sd_rdr = sd(rdr)/mean(rdr), 
                          mean_rdr_percn = mean(rdr/tot_cn), 
                          norm_sd_rdr_perc = sd(rdr/tot_cn)/mean(rdr/tot_cn), 
                          mean_baf = mean(baf_consistent), 
                          norm_sd_baf = sd(baf_consistent)/mean(baf_consistent),
                          mean_baf_perAI = mean(baf_consistent/AI),
                          norm_sd_baf_perAI = sd(baf_consistent/AI)/mean(baf_consistent/AI),
                          is_all_loh = all(min_cn == 0),
                          is_any_loh = any(min_cn == 0),
                          nmuts_perAl = sum(nmuts_perAl)), 
                        by = .(start, end,  quality_bin)]

pdf( paste0(outputs.folder, date, '_rdr_sd_by_quality_bin.pdf') )
ggplot(seg_sum, aes(x = nmuts_perAl, y = norm_sd_baf)) +
    geom_point() 
dev.off()


## from the var counts looks like some of the phasing is still wrong - means averaging over mutation isn't helping as much as it should
test_muts <- phased_muts[ !phase == '' & (early_mut) & chr == '17' & tot_cn > 2.1 ]
test_muts[, hap_count := ifelse(varcount_tumour_corrected > 50, 1, 0)]
test_muts[, hap_status := paste(phase, hap_count)]
test_muts[, num_hapstates := length(unique(hap_status)), by = block_id ]
test_muts[ num_hapstates > 2, .(chr, pos, hap_status, varcount_tumour_corrected, hap1_vaf, hapdepth, block_id) ]

# variants that look like they are incorrectly phased on chr 17:
# both very clear from the hapbaf/hapdepth
incor_phase_mut_positions <- c(70793230, 72473349)
incor_phase_mut_blocks <- c(70760312, 71806347)

plot_data <- bat_snps[ block == incor_phase_mut_blocks[1] ]

pdf( paste0(outputs.folder, date, '_SNPphasingtesteg1.pdf') )
ggplot(plot_data, aes(x = pos, y = baf, colour = gt)) +
    geom_point() +
    geom_vline( xintercept = incor_phase_mut_positions[1] )
dev.off()

plot_data <- bat_snps[ block == incor_phase_mut_blocks[2] ]

pdf( paste0(outputs.folder, date, '_SNPphasingtesteg2.pdf') )
ggplot(plot_data, aes(x = pos, y = baf, colour = gt)) +
    geom_point() +
    geom_vline( xintercept = incor_phase_mut_positions[2] )
dev.off()

test_muts[ pos %in% incor_phase_mut_positions, 
            .(chr, pos, ref, alt, hap1_vaf, hapdepth, hap1_id, hap1_bases, 
              hap2_bases, maj_cn, min_cn, mut_cpn, phased_allele_cn, block_id, phas_baf) ]
phasing_long[ pos_mut %in% incor_phase_mut_positions,
                .(mutid, pos_hetsnp, ref_hetsnp, alt_hetsnp, SNV_MT_SNP_MT_N, 
                  SNV_MT_SNP_WT_N, snp_vaf_on_mut, gt, block_id, hap1_vaf,
                  hap1_bases, hap2_bases, phase) ]
phasing_long[ block_id %in% incor_phase_mut_blocks,
                .(mutid, pos_hetsnp, ref_hetsnp, alt_hetsnp, SNV_MT_SNP_MT_N, 
                  SNV_MT_SNP_WT_N, snp_vaf_on_mut, gt, block_id, hap1_vaf,
                  hap1_bases, hap2_bases, phase) ]


# rdr & baf across genome -corrected
mut_baf_logR[, plot_rdr := ifelse( rdr_cor > 1200, 1200, rdr_cor) ]
plot_data <- mut_baf_logR[ , .(mean_rdr = mean(plot_rdr, na.rm=T)), by = .(Start.cummulitive, End.cummulitive)]

pdf(paste0(outputs.folder, date, '_rdr_cor_accross_genome.pdf'), width = 15, height = 5)
    ggplot() +
    geom_point(data = mut_baf_logR, shape = 21,
                aes(x = bin_centre_cum, y = plot_rdr, fill = baf_consistent)) +
    geom_rect(data = plot_data, aes(xmin = Start.cummulitive, 
                                    xmax = End.cummulitive , 
                                    ymin = mean_rdr, 
                                    ymax = mean_rdr), 
              colour = 'red', size = 0.5) +
    scale_fill_distiller( palette = 'YlGnBu', direction = 1 ) +
    scale_x_continuous(label = plot_cum$chr, breaks = plot_cum$center) +
    scale_y_continuous(breaks = seq(0, 1200, 200), labels = c(seq(0, 1000, 200), '>1200'), limits =) +
    labs( x = 'Chromosome', 
          y = 'Sum mean SNV varcount across alleles',
          fill = 'hapBAF' ) +
    theme_bw() +
    theme( text = element_text( size = 16 ))
dev.off()

phased_muts[, pos_cum := pos + chrlength.human_before_cum[ names(chrlength.human_before_cum) == chr ], 
             1: nrow(phased_muts)]

pdf(paste0(outputs.folder, date, '_corrected_varcounts_across_genome_phased.pdf'), width = 15, height = 5)
    ggplot(phased_muts[ !phase == '' & (early_mut) ]) +
    geom_point(aes(x = pos_cum, y = varcount_tumour_corrected, colour = phase ), alpha = 0.5 ) +
    scale_x_continuous(label = plot_cum$chr, breaks = plot_cum$center) +
    labs( x = 'Chromosome', 
          y = 'Mutation Varcount Corrected' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()


## make a hatchet style plot
pdf(paste0(outputs.folder, date, '_hatchet_cor.pdf'), width = 9)
    ggplot(mut_baf_logR) +
    geom_point(aes(x = baf_cor_consistent, y = log2(rdr_cor), colour = chr ), alpha = 0.5 ) +
    labs( x = '0.5-BAF', 
          y = 'log2 sum allele readcounts' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()


#look at hist of all LogR  
pdf(paste0(outputs.folder, date, '_BAF_cor_hist.pdf'), width = 9)
    ggplot(mut_baf_logR) +
    geom_histogram(aes(x = baf_cor_consistent, fill = chr ) ) +
    labs( x = '0.5-BAF' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

#look at hist of all BAF  
pdf(paste0(outputs.folder, date, '_rd_corr_hist.pdf'), width = 9)
    ggplot(mut_baf_logR) +
    geom_histogram(aes(x = log2(rdr_cor), fill = chr ) ) +
    labs( x = 'log2 sum allele readcounts' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()


#look at hist of all LogR  -qual
pdf(paste0(outputs.folder, date, '_BAF_cor_hist_qual.pdf'), width = 9)
    ggplot(mut_baf_logR[ (quality_bin) ]) +
    geom_histogram(aes(x = baf_cor_consistent, fill = chr ) ) +
    labs( x = '0.5-BAF' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

#look at hist of all BAF  -qual
pdf(paste0(outputs.folder, date, '_rdr_cor_hist_qual.pdf'), width = 9)
    ggplot(mut_baf_logR[ (quality_bin) ]) +
    geom_histogram(aes(x = log2(rdr_cor), fill = chr ) ) +
    labs( x = 'log2 sum allele readcounts' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

## make a hatchet style plot
pdf(paste0(outputs.folder, date, '_hatchet_cor_quality.pdf'), width = 9)
    ggplot(mut_baf_logR[ (quality_bin | (maj_cn / tot_cn) < 0.6) ]) +
    geom_point(aes(x = baf_cor_consistent, y = rdr_cor, colour = as.factor(chr) ) ) +
    labs( x = '0.5-BAF', 
          y = 'sum allele readcounts',
          colour = '' ) +
    scale_y_continuous(limits = c(0,1000)) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

## make a hatchet style plot + density
pdf(paste0(outputs.folder, date, '_hatchet_cor_quality_density.pdf'), width = 9)
    ggplot(mut_baf_logR[ (quality_bin) ], aes(x = baf_consistent, y = log2(rdr))) +
    geom_point(aes( colour = chr ), alpha = 0.5 ) +
    stat_density_2d( alpha = 0.25,geom = "raster", contour = FALSE ) +
    labs( x = '0.5-BAF', 
          y = 'log2 sum allele readcounts' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

## make a hatchet style plot
pdf(paste0(outputs.folder, date, '_hatchet_cor_quality_3muts_perAl.pdf'), width = 9)
    ggplot(mut_baf_logR[ (quality_bin) & nmuts_perAl > 3 ]) +
    geom_point(aes(x = baf_consistent, y = log2(rdr), colour = chr ), alpha = 0.5 ) +
    labs( x = '0.5-BAF', 
          y = 'log2 sum allele readcounts' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

# compare baf and AI in quality bins
pdf(paste0(outputs.folder, date, '_rdr_cor_vs_totcn_quality.pdf'), width = 9)
    mut_baf_logR[ (quality_bin), plot(log2(rdr_cor), tot_cn) ]
dev.off()

# compare baf and AI in quality bins
pdf(paste0(outputs.folder, date, '_baf_cor_vs_AI_quality.pdf'), width = 9)
    mut_baf_logR[ (quality_bin), plot(AI, baf_cor_consistent) ]
dev.off()

# drdr across genome for quality bins
plot_data <- mut_baf_logR[ (quality_bin | (maj_cn / tot_cn) < 0.6 ) , .(mean_rdr = mean(rdr_cor, na.rm=T)), by = .(Start.cummulitive, End.cummulitive)]

pdf(paste0(outputs.folder, date, '_rdr_cor_accross_genome_quality_.pdf'), width = 15, height = 5)
    ggplot() +
    geom_point(data = mut_baf_logR[(quality_bin | (maj_cn / tot_cn) < 0.6 )], shape = 21,
                aes(x = bin_centre_cum, y = rdr_cor, fill = baf_cor_consistent)) +
    geom_rect(data = plot_data, aes(xmin = Start.cummulitive, 
                                    xmax = End.cummulitive , 
                                    ymin = mean_rdr, 
                                    ymax = mean_rdr), 
              colour = 'red', size = 0.5) +
    scale_fill_distiller( palette = 'YlGnBu', direction = 1 ) +
    scale_x_continuous(label = plot_cum$chr, breaks = plot_cum$center) +
    scale_y_continuous(limits = c(0,1000)) +
    labs( x = 'Chromosome', 
          y = 'Sum mean SNV varcount across alleles',
          fill = 'hapBAF' ) +
    theme_bw() +
    theme( text = element_text( size = 16 ))
dev.off()


# drdr across genome for quality bins using refn control
plot_data <- mut_baf_logR[ (quality_bin | (maj_cn / tot_cn) < 0.6 ) , .(mean_rdr = mean(rdr_cor_rc, na.rm=T)), by = .(Start.cummulitive, End.cummulitive)]

pdf(paste0(outputs.folder, date, '_rdr_cor_accross_genome_refcontrol.pdf'), width = 15, height = 5)
    ggplot() +
    geom_point(data = mut_baf_logR, shape = 21,
                aes(x = bin_centre_cum, y = rdr_cor_rc)) +
    geom_rect(data = plot_data, aes(xmin = Start.cummulitive, 
                                    xmax = End.cummulitive , 
                                    ymin = mean_rdr, 
                                    ymax = mean_rdr), 
              colour = 'red', size = 0.5) +
    scale_fill_distiller( palette = 'YlGnBu', direction = 1 ) +
    scale_x_continuous(label = plot_cum$chr, breaks = plot_cum$center) +
    #scale_y_continuous(limits = c(0,1000)) +
    labs( x = 'Chromosome', 
          y = 'Sum mean SNV varcount across alleles',
          fill = 'hapBAF' ) +
    theme_bw() +
    theme( text = element_text( size = 16 ))
dev.off()

unique( mut_baf_logR[ chr == '3' & bin_start > 150000000 & rdr_cor < 600 ] )

## make a plot of BAFs in each bin to check cn segmenation phasing
bat_snps[, pos_cum := pos + chrlength.human_before_cum[ names(chrlength.human_before_cum) == chr ], 
             1:nrow(bat_snps)]

pdf(paste0(outputs.folder, date, '_raw_varcounts_SNPs_phased_chr3.pdf'), width = 15, height = 5)
    ggplot(bat_snps[ !is.na(phase) & chr == '3' & pos > 150000000 ]) +
    geom_point(aes(x = pos_cum, y = var_count, colour = phase ), alpha = 0.5 ) +
    scale_x_continuous(label = plot_cum$chr, breaks = plot_cum$center) +
    scale_y_continuous(limits=c(0, 1200)) +
    labs( x = 'Chromosome', 
          y = 'SNP Varcount' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()


# drdr across genome for quality bins & limited baf
plot_data <- mut_baf_logR[ (quality_bin) , .(mean_rdr = mean(rdr_cor, na.rm=T)), by = .(Start.cummulitive, End.cummulitive)]
mut_baf_logR[, baf_cor_consistent_lim := ifelse( baf_cor_consistent > 0.2, 0.2, baf_cor_consistent) ]

pdf(paste0(outputs.folder, date, '_rdr_cor_accross_genome_quality_baf_lim0.2.pdf'), width = 15, height = 5)
    ggplot() +
    geom_point(data = mut_baf_logR[(quality_bin)], shape = 21,
                aes(x = bin_centre_cum, y = rdr_cor, fill = baf_cor_consistent_lim)) +
    geom_rect(data = plot_data, aes(xmin = Start.cummulitive, 
                                    xmax = End.cummulitive , 
                                    ymin = mean_rdr, 
                                    ymax = mean_rdr), 
              colour = 'red', size = 0.5) +
    scale_fill_distiller( palette = 'YlGnBu', direction = 1 ) +
    scale_x_continuous(label = plot_cum$chr, breaks = plot_cum$center) +
    labs( x = 'Chromosome', 
          y = 'Sum mean SNV varcount across alleles',
          fill = 'hapBAF' ) +
    theme_bw() +
    theme( text = element_text( size = 16 ))
dev.off()

# restrict to blocks with more mutations
mut_baf_logR[, nmuts := n_hap1_muts + n_hap2_muts]
mut_baf_logR[, nmuts_perAl := ifelse( min_cn == 0, nmuts, nmuts/2 )]

plot_data <- mut_baf_logR[ (quality_bin) & nmuts_perAl > 2 , .(mean_rdr = mean(rdr_cor, na.rm=T)), by = .(Start.cummulitive, End.cummulitive)]

pdf(paste0(outputs.folder, date, '_rdr_cor_accross_genome_quality_3muts_per_al.pdf'), width = 15, height = 5)
    ggplot() +
    geom_point(data = mut_baf_logR[(quality_bin)  & nmuts_perAl > 2], shape = 21,
                aes(x = bin_centre_cum, y = rdr_cor, fill = baf_consistent)) +
    geom_rect(data = plot_data, aes(xmin = Start.cummulitive, 
                                    xmax = End.cummulitive , 
                                    ymin = mean_rdr, 
                                    ymax = mean_rdr), 
              colour = 'red', size = 0.5) +
    scale_fill_distiller( palette = 'YlGnBu', direction = 1 ) +
    scale_x_continuous(label = plot_cum$chr, breaks = plot_cum$center) +
    labs( x = 'Chromosome', 
          y = 'Sum mean SNV varcount across alleles',
          fill = 'hapBAF' ) +
    theme_bw() +
    theme( text = element_text( size = 16 ))
dev.off()
# doesn't really seem to decrease noise

###########################################
#### Plot for test lower purity region ####
###########################################

test_data[, quality_bin := bat_snps[ match( test_data$bin_id, bat_snps$bin_id ), quality_bin ]]
quasar_test[, quality_bin := bat_snps[ match( quasar_test$bin_id, bat_snps$bin_id ), quality_bin ]]


quasar_test[, bin_centre := (bin_end - bin_start)/2 + bin_start ]
quasar_test[, bin_centre_cum := bin_centre + chrlength.human_before_cum[ names(chrlength.human_before_cum) == chr ], 
               1:nrow(quasar_test) ]
setnames(quasar_test, c('ref_seg_start', 'ref_seg_end'), c('start', 'end'))
quasar_test[, `:=`(Start.cummulitive = start + chrlength.human_before_cum[ names(chrlength.human_before_cum) == chr ], 
                    End.cummulitive = end + chrlength.human_before_cum[ names(chrlength.human_before_cum) == chr ]),
               1:nrow(quasar_test) ]

# rdr & baf across genome -cor
plot_data <- quasar_test[ , .(mean_rdr = mean(rdr_cor, na.rm=T)), by = .(Start.cummulitive, End.cummulitive)]

pdf(paste0(outputs.folder, date, '_rdr_cor_accross_genome_test.pdf'), width = 15, height = 5)
    ggplot() +
    geom_point(data = quasar_test, shape = 21,
                aes(x = bin_centre_cum, y = rdr_cor, fill = baf_cor_consistent)) +
    geom_rect(data = plot_data, aes(xmin = Start.cummulitive, 
                                    xmax = End.cummulitive , 
                                    ymin = mean_rdr, 
                                    ymax = mean_rdr), 
              colour = 'red', size = 0.5) +
    scale_fill_distiller( palette = 'YlGnBu', direction = 1 ) +
    scale_x_continuous(label = plot_cum$chr, breaks = plot_cum$center) +
    labs( x = 'Chromosome', 
          y = 'Sum mean SNV varcount across alleles',
          fill = 'hapBAF' ) +
    theme_bw() +
    theme( text = element_text( size = 16 ))
dev.off()

# rdr & baf across genome - quality cor
plot_data <- quasar_test[ (quality_bin & (!is.na(baf_cor_consistent) | is_mirror_vs_loh_ref)), .(mean_rdr = mean(rdr_cor, na.rm=T)), by = .(Start.cummulitive, End.cummulitive)]

pdf(paste0(outputs.folder, date, '_rdr_cor_accross_genome_quality_test.pdf'), width = 15, height = 5)
    ggplot() +
    geom_point(data = quasar_test[ (quality_bin & (!is.na(baf_cor_consistent) | is_mirror_vs_loh_ref)) ], shape = 21,
                aes(x = bin_centre_cum, y = rdr_cor, fill = baf_cor_consistent)) +
    geom_rect(data = plot_data, aes(xmin = Start.cummulitive, 
                                    xmax = End.cummulitive , 
                                    ymin = mean_rdr, 
                                    ymax = mean_rdr), 
              colour = 'red', size = 0.5) +
    scale_fill_distiller( palette = 'YlGnBu', direction = 1 ) +
    scale_x_continuous(label = plot_cum$chr, breaks = plot_cum$center) +
    labs( x = 'Chromosome', 
          y = 'Sum mean SNV varcount across alleles',
          fill = 'hapBAF' ) +
    theme_bw() +
    theme( text = element_text( size = 16 ))
dev.off()

pdf(paste0(outputs.folder, date, '_rdr_cor_accross_genome_quality_longblocks_test.pdf'), width = 15, height = 5)
    ggplot() +
    geom_point(data = quasar_test[ (quality_bin) & n_hap1_muts > 1 & n_hap2_muts > 1 ], shape = 21,
                aes(x = bin_centre_cum, y = rdr_cor, fill = baf_cor_consistent)) +
    geom_rect(data = plot_data, aes(xmin = Start.cummulitive, 
                                    xmax = End.cummulitive , 
                                    ymin = mean_rdr, 
                                    ymax = mean_rdr), 
              colour = 'red', size = 0.5) +
    scale_fill_distiller( palette = 'YlGnBu', direction = 1 ) +
    scale_x_continuous(label = plot_cum$chr, breaks = plot_cum$center) +
    labs( x = 'Chromosome', 
          y = 'Sum mean SNV varcount across alleles',
          fill = 'hapBAF' ) +
    theme_bw() +
    theme( text = element_text( size = 16 ))
dev.off()


# rdr & baf across genome - cor
plot_data <- quasar_test[ , .(mean_rdr = mean(rdr, na.rm=T)), by = .(Start.cummulitive, End.cummulitive)]

pdf(paste0(outputs.folder, date, '_rdr_accross_genome_quality_test.pdf'), width = 15, height = 5)
    ggplot() +
    geom_point(data = quasar_test[ (quality_bin) ], shape = 21,
                aes(x = bin_centre_cum, y = rdr, fill = baf_consistent)) +
    geom_rect(data = plot_data, aes(xmin = Start.cummulitive, 
                                    xmax = End.cummulitive , 
                                    ymin = mean_rdr, 
                                    ymax = mean_rdr), 
              colour = 'red', size = 0.5) +
    scale_fill_distiller( palette = 'YlGnBu', direction = 1 ) +
    scale_x_continuous(label = plot_cum$chr, breaks = plot_cum$center) +
    labs( x = 'Chromosome', 
          y = 'Sum mean SNV varcount across alleles',
          fill = 'hapBAF' ) +
    theme_bw() +
    theme( text = element_text( size = 16 ))
dev.off()

# Look at raw counts per mutation


test_data[, pos_cum := pos + chrlength.human_before_cum[ names(chrlength.human_before_cum) == chr ], 
             1:nrow(test_data)]

pdf(paste0(outputs.folder, date, '_raw_varcounts_across_genome_phased_test.pdf'), width = 15, height = 5)
    ggplot(test_data) +
    geom_point(aes(x = pos_cum, y = varcount_tumour, colour = phase ), alpha = 0.5 ) +
    scale_x_continuous(label = plot_cum$chr, breaks = plot_cum$center) +
    labs( x = 'Chromosome', 
          y = 'Mutation Varcount' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

pdf(paste0(outputs.folder, date, '_cor_varcounts_across_genome_phased_test.pdf'), width = 15, height = 5)
    ggplot(test_data) +
    geom_point(aes(x = pos_cum, y = varcount_tumour, colour = phase ), alpha = 0.5 ) +
    scale_x_continuous(label = plot_cum$chr, breaks = plot_cum$center) +
    labs( x = 'Chromosome', 
          y = 'Mutation Varcount' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

pdf(paste0(outputs.folder, date, '_varcounts_hist_test.pdf'), width = 15, height = 5)
 test_data[ varcount_tumour < 100, hist(varcount_tumour, breaks = 100)]
dev.off()

#add on CN data

#remove unwanted columns 
cn_test <- cn_test[, 1:13]

# Melt by clone
setnames(cn_test, gsub('_A', '', names(cn_test)))
cn_test <- data.table::melt(cn_test, id.vars = c('chr', 'startpos', 'endpos'),
           measure = patterns('^nMaj', '^nMin', '^frac'),
           value.name = c("nMaj", "nMin", "frac"),
           variable.name = 'clone_id')

# remove NAs where fewer than max number of clones for a segment
cn_test <- cn_test[ !is.na(nMaj) ]

# average across clones 
cn_test <- cn_test[, .(nMaj = sum(nMaj * frac),
             nMin = sum(nMin * frac)),
         by = .(chr, startpos, endpos) ]
cn_test[, ntot := nMaj + nMin ]

ph_granges <- GenomicRanges::makeGRangesFromDataFrame(quasar_test)

#get phsed cn object
setnames(cn_test, c('startpos', 'endpos'), c('start', 'end'))
cn_granges <- GenomicRanges::makeGRangesFromDataFrame(cn_test)

# get indices for mutations in cn obj and then copy over information
match <- GenomicRanges::findOverlaps(ph_granges, cn_granges, select = 'arbitrary')
quasar_test[, maj_cn := cn_test[ match, nMaj ] ]
quasar_test[, min_cn := cn_test[ match, nMin ] ]
quasar_test[, tot_cn := cn_test[ match, ntot ] ]

# copy over the segment borders of the mutation so we resist to this area when calculating
# the hapbaf
quasar_test[, seg_start := cn_test[ match, start ] ]
quasar_test[, seg_end := cn_test[ match, end ] ]

plot_cn_data <- quasar_test[, .(Start.cummulitive, End.cummulitive, maj_cn, min_cn, tot_cn) ]
plot_cn_data <- melt(plot_cn_data, id.vars = c('Start.cummulitive', 'End.cummulitive'), 
                     measure = patterns( '*_cn$' ),
                     variable.name = 'allele',
                     value.name = 'cn' )
plot_cn_data[, cn := as.numeric(cn) ]
plot_cn_data[, cn := ifelse( cn > 5, 5, cn) ]

plot_cn_data[, cn_plot := ifelse(allele == 'maj_cn', cn + 0.01, ifelse( allele == 'min_cn',  cn - 0.01, cn))]

pdf(paste0(outputs.folder, date, '_ascn_test.pdf'), width = 15, height = 5)
ggplot() +
    geom_rect(data = plot_cn_data[ !allele == 'tot_cn' ], aes(xmin = Start.cummulitive, xmax = End.cummulitive , ymin = cn_plot, ymax = cn_plot, colour = allele)) +
    scale_x_continuous(label = plot_cum$chr, breaks = plot_cum$center) +
    scale_y_continuous( breaks = 0:5, labels = c(0:4, '>5'))+
    scale_colour_brewer( palette = 'Set1' ) +
    labs( x = 'Chromosome', 
          y = 'cn' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()


## make a hatchet style plot 
pdf(paste0(outputs.folder, date, '_hatchet_cor_quality_test.pdf'), width = 9)
    ggplot(mut_baf_logR[ (quality_bin) ]) +
    geom_point(aes(x = baf_cor_consistent, y = log2(rdr_cor), colour = chr ), alpha = 0.5 ) +
    labs( x = '0.5-BAF', 
          y = 'log2 sum allele readcounts' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()


### look at hap varcount per block
plot_data <- melt(quasar_test[ (quality_bin), .(block_id, bin_centre_cum, hap1=mean_hap1_varcount, hap2=mean_hap2_varcount)], 
                    id.vars = c('block_id', 'bin_centre_cum'), measure = patterns( '^hap') )

pdf(paste0(outputs.folder, date, '_hapcounts_per_block_test_quality.pdf'), width = 15, height = 5)
    ggplot() +
    geom_point(data = plot_data, shape = 21,
                aes(x = bin_centre_cum, y = value, fill = variable)) +
    scale_x_continuous(label = plot_cum$chr, breaks = plot_cum$center) +
    labs( x = 'Chromosome', 
          y = 'SNV varcount in each haplotype',
          fill = 'Phase' ) +
    theme_bw() +
    theme( text = element_text( size = 16 ))
dev.off()


# try plotting all the blocks in one chromosome - nice seperation in 15
blocks <- test_data[ chr == '15' & pos > 60000000 , 
                     .(block_start = min(pos)), by = block_id ]

pdf(paste0(outputs.folder, date, '_chr_15_over8E7_raw_block_data_test.pdf'), width = 15, height = 5)
    ggplot() +
    geom_point(data = test_data[ chr == '15' & pos > 60000000 & (quality_bin) ], shape = 21,
                aes(x = pos, y = varcount_tumour, fill = phase)) +
    geom_segment( data = blocks, aes(x = block_start, xend = block_start, y = 0, yend = 100) ) +
    labs( x = 'Chromosome', 
          y = 'SNV varcount in each haplotype',
          fill = 'Phase' ) +
    theme_bw() +
    theme( text = element_text( size = 16 ))
dev.off()

# loaddss of apparent mis-phasing - is this the same for the refernce sample?

chr_test <- '7'
pos_start_test <- 60000000
pos_end_test <- 160000000

blocks <- phased_muts[ !phase == '' & (early_mut)  & chr == chr_test & pos > pos_start_test & pos < pos_end_test, 
                     .(bin_start = min(pos)), by = bin_id ]


pdf(paste0(outputs.folder, date, '_chr_7_raw_block_data.pdf'), width = 15, height = 5)
    ggplot() +
    geom_point(data = phased_muts[ !phase == '' & chr == chr_test & pos > pos_start_test & 
                                    pos < pos_end_test & (early_mut) ], shape = 21,
                aes(x = pos, y = varcount_tumour_corrected , fill = phase)) +
    #geom_segment( data = blocks, aes(x = bin_start, xend = bin_start, y = 0, yend = 1), alpha = 0.5, colour = 'grey' ) +
    #scale_y_continuous( limits = c(0.25,1)) +
    labs( x = 'Chromosome', 
          y = 'SNV varcount corrected in each haplotype',
          fill = 'Phase' ) +
    theme_bw() +
    theme( text = element_text( size = 16 ))
dev.off()

pdf(paste0(outputs.folder, date, '_chr_7_raw_block_data_vaf.pdf'), width = 15, height = 5)
    ggplot() +
    geom_point(data = phased_muts[ !phase == '' & chr == chr_test & pos > pos_start_test & 
                                    pos < pos_end_test & (early_mut) ], shape = 21,
                aes(x = pos, y = varcount_tumour_corrected/ (varcount_tumour_corrected + refcount_tumour) , fill = phase)) +
    #geom_segment( data = blocks, aes(x = bin_start, xend = bin_start, y = 0, yend = 1), alpha = 0.5, colour = 'grey' ) +
    #scale_y_continuous( limits = c(0.25,1)) +
    labs( x = 'Chromosome', 
          y = 'SNV vaf corrected in each haplotype',
          fill = 'Phase' ) +
    theme_bw() +
    theme( text = element_text( size = 16 ))
dev.off()


plot_cn_data <- mut_baf_logR[, .(Start.cummulitive, End.cummulitive, maj_cn, min_cn, tot_cn, chr, start, end) ]
plot_cn_data <- unique( melt(plot_cn_data, id.vars = c('Start.cummulitive', 'End.cummulitive', 'chr', 'start', 'end'), 
                     measure = patterns( '*_cn$' ),
                     variable.name = 'allele',
                     value.name = 'cn' ) )
plot_cn_data[, cn := as.numeric(cn) ]
plot_cn_data[, cn := ifelse( cn > 5, 5, cn) ]

plot_cn_data[, cn_plot := ifelse(allele == 'maj_cn', cn + 0.01, ifelse( allele == 'min_cn',  cn - 0.01, cn))]

pdf(paste0(outputs.folder, date, '_R3_ascn_test.pdf'), width = 15, height = 5)
ggplot() +
    geom_rect(data = plot_cn_data[ !allele == 'tot_cn' & chr == chr_test & end > pos_start_test & start < pos_end_test ], 
              aes(xmin = Start.cummulitive, xmax = End.cummulitive , ymin = cn_plot, ymax = cn_plot, colour = allele)) +
    scale_x_continuous(label = plot_cum$chr, breaks = plot_cum$center) +
    scale_y_continuous( breaks = 0:5, labels = c(0:4, '>5'), limits = c(0,5))+
    scale_colour_brewer( palette = 'Set1' ) +
    labs( x = 'Chromosome', 
          y = 'cn' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

# check snp phasing where looks wrong

plot_data <- bat_snps[ block == 83825235 ]

pdf( paste0(outputs.folder, date, '_SNPphasingtesteg3.pdf') )
ggplot(plot_data, aes(x = pos, y = baf, colour = gt)) +
    geom_point() #+
    #geom_vline( xintercept = incor_phase_mut_positions[1] )
dev.off()

phased_muts[ block_id == 83825235 ,.(pos, phase, block_id, quality_bin)]
# still called quality!

## make plots for reference control

plot_data <- melt( mut_baf_logR[, .(chr, bin_centre_cum, hap1=hap1_vaf_cor, hap2=hap2_vaf_cor) ],
                   id.vars = c('chr', 'bin_centre_cum'))

pdf(paste0(mut_baf_logR, date, '_raw_vafs_across_genome_phased.pdf'), width = 15, height = 5)
    ggplot(plot_data) +
    geom_point(aes(x = bin_centre_cum, y = value, colour = variable ), alpha = 0.5 ) +
    scale_x_continuous(label = plot_cum$chr, breaks = plot_cum$center) +
    labs( x = 'Chromosome', 
          y = 'Mutation hap VAF' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

#############
#### END ####
#############

