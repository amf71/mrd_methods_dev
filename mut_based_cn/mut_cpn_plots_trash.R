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

# Phased mutations data processed R3
fwrite(phased_muts, paste0(outputs.folder, date, '_phased_counts_by_mut_R3.tsv'), sep = '\t')

# Data processed by block baf/logR R3
fwrite(mut_baf_logR, paste0(outputs.folder, date, '_baf_rdr_by_block_R3.tsv'), sep = '\t')

# Phased mutations extracted R6
fwrite(test_data, paste0(outputs.folder, date, '_quasar_phased_counts_by_mut_R6.tsv'), sep = '\t')

# Data processed by block baf/logR R6 quasar
fwrite(quasar_test, paste0(outputs.folder, date, '_quasar_baf_rdr_by_block_R6.tsv'), sep = '\t')


# CN from phased region
cn_path  <- "/camp/lab/swantonc/working/opich/shared/afrankell/wgs/LTX0064_SU_T1-R3_subclones.txt"
cn_test_path  <- "inputs/LTX0064_R6_Battenberg//LTX0064_SU_T1-R6_subclones.txt"

# purity from pahsed region (from battenberg)
purity_path <- "/camp/lab/swantonc/working/opich/shared/afrankell/wgs/LTX0064_SU_T1-R3_purity_ploidy.txt"

#het SNPs from phased tumour region to get allele identities to match phasing
snps_path <- 'outputs/phasing/20230627/scratch/LTX0064_SU_N_filt_hetsnps_mpileup.vcf'
bat_snps_dir <- '/camp/project/proj-tracerx-lung/tctProjects/frankella/ctDNA_methods_dev/inputs/LTX0064_R3_Battenberg/'
bat_snps_pattern <- 'LTX0064_SU_T1-R3_alleleFrequencies*'

# Read in
cn    <-  data.table::fread( cn_path )
cn_test    <-  data.table::fread( cn_test_path )
phasing_long  <- fread( phasing_long_path )

# collate tumour SNP read counts from battenberg
files <- list.files( bat_snps_dir, pattern = bat_snps_pattern, full.names = TRUE )
bat_snps <- rbindlist( lapply(files, fread ))


############################################################
#### Make some QC plots on which mutations are phasable ####
############################################################
phas[, hapdepth_plot := ifelse(hapdepth > 50, 50, hapdepth )] ## phas to phased_muts

# Above doesn't include mutations with no reads overlapping a hetSNP
plot_data <- phasing_long[(gt == '0|1' | gt == '1|0'), 
                            .(n_hetsnps = .N, mean_dist = mean(dist_mut_to_snp), 
                              hapdepth = sum(mut_varcount), total_depth = sum(shared_depth)), 
                          by = mutid ]
plot_data[, hapdepth_plot := ifelse(hapdepth > 50, 50, hapdepth )]

plot_data[, mean_dist_plot := ifelse(mean_dist > 10000, 10000, mean_dist )]

## Look at how many SNPs associated with each mutations and average distances
pdf( paste0( outputs.folder, date, '_N_snps_vs_snp_dist.pdf') )
ggplot(plot_data, aes(x = n_hetsnps, y = mean_dist_plot )) +
    geom_density_2d_filled( bins = 20 ) +
    scale_y_continuous( breaks = seq(0, 10000, 2000), labels = c(seq(0, 8000, 2000), '>10000') ) +
    scale_x_continuous( breaks = seq(0, 10, 1) ) +
    theme_classic() +
    labs( x = 'Number of hetSNPs associated with each mut',
          y = 'Mean distance from hetSNPs to mut') +
    theme( text = element_text( size = 20),
           legend.position = "none" )
dev.off()

## Look at how many SNPs associated with each mutations and average distances in those without any depth
pdf( paste0( outputs.folder, date, '_N_snps_vs_snp_dist_hapdepth0.pdf') )
ggplot(plot_data[ hapdepth == 0 ], aes(x = n_hetsnps, y = mean_dist_plot )) +
    geom_density_2d_filled( bins = 20 ) +
    scale_y_continuous( breaks = seq(0, 10000, 2000), labels = c(seq(0, 8000, 2000), '>10000') ) +
    scale_x_continuous( breaks = seq(0, 10, 1) ) +
    theme_classic() +
    labs( x = 'Number of hetSNPs associated with each mut',
          y = 'Mean distance from hetSNPs to mut') +
    theme( text = element_text( size = 20),
           legend.position = "none" )
dev.off()

pdf( paste0( outputs.folder, date, '_total_depth_vs_0hapdepth.pdf') )
ggplot(plot_data, aes(x = hapdepth == 0, y = total_depth )) +
    geom_boxplot() +
    scale_y_continuous(  trans = 'log10' ) +
    #scale_x_continuous( breaks = seq(0, 10, 1) ) +
    theme_classic() +
    #labs( x = 'Number of hetSNPs associated with each mut',
    #      y = 'Mean distance from hetSNPs to mut') +
    theme( text = element_text( size = 20),
           legend.position = "none" )
dev.off()

#####################################################
#### Match alleles and determine early mutations ####
#####################################################

### Lets do some checks to see how good phasing is using areas of imbalence in the tumour samples

#copy over AI for each segment to Bat SNPs to pull out ares of imbalence and thier haplotypes
cn[, AI := nMaj/(nMin+nMaj) ]
cn_granges <- GenomicRanges::makeGRangesFromDataFrame(cn)
bat_snps_granges <- copy(bat_snps)
setnames(bat_snps_granges, 'pos', 'start')
bat_snps_granges[, end := start ]
bat_snps_granges <- GenomicRanges::makeGRangesFromDataFrame(bat_snps_granges)
match <- GenomicRanges::findOverlaps(bat_snps_granges, cn_granges, select = 'arbitrary')
bat_snps[, AI := cn[ match, AI ] ]
bat_snps[, seg_name := cn[ match, paste(chr, start, end, sep = ':') ] ]

# just look at complete imbalence and SNPs that are very clear
bat_snps[ AI > 0.75 , phase := ifelse(baf > 0.6, 'hap1', ifelse( baf < 0.4, 'hap2', NA ))]

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
bat_snps[, N_haps := length(unique(paste(phase[!is.na(phase)], gt))), by = block ]
bat_snps[, N_seg := length(unique(seg_name)), by = block ]
#bat_snps[ N_seg==1, table(N_haps)]
#bat_snps[ block == 70760312 ]

plot_data <- bat_snps[ block == bat_snps[ AI > 0.75, unique(block)[2]] ]

pdf( paste0(outputs.folder, date, '_SNPphasingtesteg.pdf') )
ggplot(plot_data, aes(x = start, y = baf, colour = gt)) +
    geom_point()
dev.off()

# Overlay onto phasing whether there is consistency with BAF seperation based 
# phasing in each block
block_summary <- unique(bat_snps[ !is.na(block), .(block, block_length, N_haps, N_seg)])
phased_muts[, quality_block := block_summary[ match( block_id, block ), N_haps == 2 ]]


phased_muts[ , poss_misphased := phased_allele_cn - mut_cpn < -0.5 ]
phased_muts[, table(poss_misphased, quality_block)]
sum <- phased_muts[!phase=='', .(nmut = .N, block_length=max(pos)-min(pos)), by = block_id]
sum[, gn_frac := block_length/3000000000]
sum[, sum(gn_frac)*100, by = nmut>6 ][ order(nmut)]
# need to get these numbers up with longer phasing (+ also more accurate!)

# add on the SNP based BAF/RDR



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
rowi = 6
mut_chr = phased_muts[rowi, chr]
seg_start = phased_muts[rowi, seg_start]
seg_end = phased_muts[rowi, seg_end]
block = phased_muts[rowi, block_id]

snps_by_chr <- lapply( bat_snps[, unique(chr)], function(snpchr) bat_snps[ chr == snpchr ])
names(snps_by_chr) <- bat_snps[, unique(chr)]
hap_by_chr <- lapply( hap_blocks[, unique(chr)], function(snpchr) hap_blocks[ chr == snpchr ])
names(hap_by_chr) <- hap_blocks[, unique(chr)]

phased_muts[, phas_baf := extract_hap_baf(snps_by_chr, 
                                          unique(chr),
                                          unique(block_id), 
                                          hap_by_chr,
                                          unique(seg_start),  
                                          unique(seg_end) ), 
            by = .(block_id, seg_start, seg_end)  ]

phased_muts[, AI := maj_cn/(maj_cn + min_cn) ]

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

# Do this first just from the same sample

# first try without any normalisation
phased_muts[, normal_depth := varcount_norm + refcount_norm ]

mut_baf_logR <- phased_muts[ (early_mut), 
                                .( mean_hap1_varcount = mean(varcount_tumour[ phase == 'hap1' ],na.rm=T),
                                   mean_hap2_varcount = mean(varcount_tumour[ phase == 'hap2' ],na.rm=T),
                                   mean_hap1_varnorm = mean((varcount_tumour/normal_depth)[ phase == 'hap1' ],na.rm=T),
                                   mean_hap2_varnorm = mean((varcount_tumour/normal_depth)[ phase == 'hap2' ],na.rm=T),
                                   n_hap1_muts = sum(phase == 'hap1'), n_hap2_muts = sum(phase == 'hap2'),
                                   block_start = min(pos),
                                   block_end = max(pos) ),
                             by = .(block_id, seg_start, seg_end, chr, maj_cn, min_cn, tot_cn, quality_block) ]

# when you have LOH one allele will have no mutations. In this case set Varcount to 0. 
mut_baf_logR[ min_cn == 0, 
                `:=`(mean_hap1_varcount = ifelse( is.na(mean_hap1_varcount), 0, mean_hap1_varcount),
                     mean_hap2_varcount = ifelse( is.na(mean_hap2_varcount), 0, mean_hap2_varcount))]


mut_baf_logR <- mut_baf_logR[ !is.na(mean_hap1_varcount) & !is.na(mean_hap2_varcount), 
                                .( chr, block_start, block_end, seg_start, seg_end, block_id, maj_cn, min_cn, tot_cn, 
                                   n_hap1_muts, n_hap2_muts, quality_block, 
                                   baf = mean_hap1_varcount / (mean_hap1_varcount + mean_hap2_varcount),
                                   rdr = mean_hap1_varcount + mean_hap2_varcount,
                                   baf_norm = mean_hap1_varnorm / (mean_hap1_varnorm + mean_hap2_varnorm),
                                   rdr_norm = mean_hap1_varnorm + mean_hap2_varnorm ) ]
mut_baf_logR[, baf_inv := 1 - baf]
mut_baf_logR[, baf_consistent := abs(baf - 0.5)]
mut_baf_logR[, AI := abs((maj_cn - min_cn)/tot_cn)]
mut_baf_logR[, baf_norm_consistent := abs(baf - 0.5)]

# plot across the genome

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
home <- '/camp/project/'
mut_baf_logR <- as.data.table(add.cummulative.start.end(as.data.frame(mut_baf_logR)))
mut_baf_logR[, `:=`(Start.cummulitive = start + chrlength.human_before_cum[ names(chrlength.human_before_cum) == chr ], 
                    End.cummulitive = end + chrlength.human_before_cum[ names(chrlength.human_before_cum) == chr ]),
               1:nrow(mut_baf_logR) ]

plot_cum <- make.axis.adjust.cummulative(chrlength=chrlength.human)
plot_cum$chr <- factor(plot_cum$chr, levels = 1:23 )
plot_cum <- plot_cum[order( as.numeric(plot_cum$chr) ),]

plot_cn_data <- mut_baf_logR[, .(Start.cummulitive, End.cummulitive, maj_cn, min_cn, tot_cn) ]
plot_cn_data <- melt(plot_cn_data, id.vars = c('Start.cummulitive', 'End.cummulitive'), 
                     measure = patterns( '*_cn$' ),
                     variable.name = 'allele',
                     value.name = 'cn' )
plot_cn_data[, cn := as.numeric(cn) ]
plot_cn_data[, cn := ifelse( cn > 5, 5, cn) ]

plot_cn_data[, cn_plot := ifelse(allele == 'maj_cn', cn + 0.01, ifelse( allele == 'min_cn',  cn - 0.01, cn))]

pdf(paste0(outputs.folder, date, '_ascn.pdf'), width = 15, height = 5)
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

mut_baf_logR[, block_centre := (block_end - block_start)/2 + block_start ]
mut_baf_logR[, block_centre_cum := block_centre + chrlength.human_before_cum[ names(chrlength.human_before_cum) == chr ], 
               1:nrow(mut_baf_logR) ]

# rdr & baf across genome
mut_baf_logR[, plot_rdr := ifelse( rdr > 1200, 1200, rdr) ]
plot_data <- mut_baf_logR[ , .(mean_rdr = mean(plot_rdr, na.rm=T)), by = .(Start.cummulitive, End.cummulitive)]

pdf(paste0(outputs.folder, date, '_rdr_accross_genome.pdf'), width = 15, height = 5)
    ggplot() +
    geom_point(data = mut_baf_logR, shape = 21,
                aes(x = block_centre_cum, y = plot_rdr, fill = baf_consistent)) +
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
    ggplot(mut_baf_logR[ (quality_block) ]) +
    geom_histogram(aes(x = baf_consistent, fill = chr ) ) +
    labs( x = '0.5-BAF' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

#look at hist of all BAF  -qual
pdf(paste0(outputs.folder, date, '_rdr_hist_qual.pdf'), width = 9)
    ggplot(mut_baf_logR[ (quality_block) ]) +
    geom_histogram(aes(x = log2(rdr), fill = chr ) ) +
    labs( x = 'log2 sum allele readcounts' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

## make a hatchet style plot
pdf(paste0(outputs.folder, date, '_hatchet_quality2.pdf'), width = 9)
    ggplot(mut_baf_logR[ (quality_block) & n_hap1_muts > 5 & n_hap2_muts > 5 ]) +
    geom_point(aes(x = baf_consistent, y = log2(rdr), colour = chr ), alpha = 0.5 ) +
    labs( x = '0.5-BAF', 
          y = 'log2 sum allele readcounts' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

# compare baf and AI in quality bins
pdf(paste0(outputs.folder, date, '_rdr_vs_totcn_quality.pdf'), width = 9)
    mut_baf_logR[ (quality_block), plot(log2(rdr), tot_cn) ]
dev.off()

# compare baf and AI in quality bins
pdf(paste0(outputs.folder, date, '_baf_vs_AI_quality.pdf'), width = 9)
    mut_baf_logR[ (quality_block) & n_hap1_muts > 1 & n_hap2_muts > 1, plot(AI, baf_consistent) ]
dev.off()

# drdr across genome for quality bins
plot_data <- mut_baf_logR[ (quality_block) , .(mean_rdr = mean(rdr, na.rm=T)), by = .(Start.cummulitive, End.cummulitive)]

pdf(paste0(outputs.folder, date, '_rdr_accross_genome_quality.pdf'), width = 15, height = 5)
    ggplot() +
    geom_point(data = mut_baf_logR[(quality_block)], shape = 21,
                aes(x = block_centre_cum, y = rdr, fill = baf_consistent)) +
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

plot_data <- mut_baf_logR[ (quality_block) & nmuts_perAl > 2 , .(mean_rdr = mean(rdr, na.rm=T)), by = .(Start.cummulitive, End.cummulitive)]

pdf(paste0(outputs.folder, date, '_rdr_accross_genome_quality_3muts_per_al.pdf'), width = 15, height = 5)
    ggplot() +
    geom_point(data = mut_baf_logR[(quality_block)  & nmuts_perAl > 2], shape = 21,
                aes(x = block_centre_cum, y = rdr, fill = baf_consistent)) +
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

seg_sum <- mut_baf_logR[ n_hap1_muts > 1 & n_hap2_muts > 1, 
                         .(tot_cn, mean_rdr = mean(rdr), 
                          norm_sd_rdr = sd(rdr)/mean(rdr), 
                          mean_rdr_percn = mean(rdr/tot_cn), 
                          norm_sd_rdr_perc = sd(rdr/tot_cn)/mean(rdr/tot_cn),
                          AI, mean_baf = mean(baf_consistent), 
                          norm_sd_baf = sd(baf_consistent)/mean(baf_consistent),
                          mean_baf_perAI = mean(baf_consistent/AI),
                          norm_sd_baf_perAI = sd(baf_consistent/AI)/mean(baf_consistent/AI),
                          is_all_loh = all(min_cn == 0),
                          is_any_loh = any(min_cn == 0),
                          nmuts_perAl = sum(nmuts_perAl)), 
                        by = .(seg_start, seg_end,  quality_block)]

pdf( paste0(outputs.folder, date, '_rdr_sd_by_quality_block.pdf') )
ggplot(seg_sum, aes(x = start, y = baf, colour = gt)) +
    geom_point() +
    geom_vline( xintercept = incor_phase_mut_positions[1] )
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
ggplot(plot_data, aes(x = start, y = baf, colour = gt)) +
    geom_point() +
    geom_vline( xintercept = incor_phase_mut_positions[1] )
dev.off()

plot_data <- bat_snps[ block == incor_phase_mut_blocks[2] ]

pdf( paste0(outputs.folder, date, '_SNPphasingtesteg2.pdf') )
ggplot(plot_data, aes(x = start, y = baf, colour = gt)) +
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



## make a hatchet style plot
pdf(paste0(outputs.folder, date, 'plot9.pdf', width = 9)
    ggplot(mut_baf_logR) +
    geom_point(aes(x = baf_consistent, y = log2(rdr), colour = chr ), alpha = 0.5 ) +
    labs( x = '0.5-BAF', 
          y = 'log2 sum allele readcounts' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()


mut_baf_logR1 <- phased_muts_test[ (early_mut), 
                                .( mean_hap1_varcount = mean(varcount_tumour[ phase == 'hap1' ],na.rm=T),
                                   mean_hap2_varcount = mean(varcount_tumour[ phase == 'hap2' ],na.rm=T),
                                   n_hap1_muts = sum(phase == 'hap1'), n_hap2_muts = sum(phase == 'hap2'),
                                   block_start = min(pos),
                                   block_end = max(pos) ),
                             by = .(block_id, chr) ]

mut_baf_logR <- mut_baf_logR1[ !is.na(mean_hap1_varcount) & !is.na(mean_hap2_varcount), 
                                .( chr, block_start, block_end, block_id, n_hap1_muts, n_hap2_muts, 
                                   n_hap1_muts, n_hap2_muts, quality_block,
                                   baf = mean_hap1_varcount / (mean_hap1_varcount + mean_hap2_varcount),
                                   rdr = mean_hap1_varcount + mean_hap2_varcount ) ]
mut_baf_logR[, baf_inv := 1 - baf]
mut_baf_logR[, baf_consistent := abs(baf - 0.5)]

pdf(paste0(outputs.folder, date, 'plot10.pdf', width = 9)
    ggplot(mut_baf_logR[n_hap1_muts>1 & n_hap2_muts>1 ]) +
    geom_point(aes(x = baf_consistent, y = log2(rdr), colour = chr ), alpha = 0.5 ) +
    labs( x = '0.5-BAF', 
          y = 'log2 sum allele readcounts' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

## Lets try GC correction
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
colnames(ph_granges) <- c("seqnames", "start", colnames(ph_granges)[-c(1,2)])
locidf[,3:ncol(ph_granges)]=round(ph_granges[,3:ncol(ph_granges)],6)
locidf <- as.data.table( locidf )
setnames(locidf, 1:2,  c('chr', 'pos'))

phased_muts <- as.data.table( inner_join(locidf, phased_muts) )
cols <- names(locidf)[3:17]
out <- cor(phased_muts[!is.na(Hepg2), varcount_tumour],
           phased_muts[ !is.na(Hepg2), ..cols])

model = lm(varcount_tumour ~ splines::ns(x = gc, df = 5, intercept = T) + 
                             splines::ns(x = Mcf7, df = 5, intercept = T) + 
                             splines::ns(x = normal_depth, df = 5, intercept = T), y=F, model = F, data = phased_muts, na.action="na.exclude")
phased_muts[, varcount_tumour_corrected := residuals(model) + mean(varcount_tumour) ]
phased_muts[ varcount_tumour_corrected < 0, varcount_tumour_corrected := 0]


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
                             by = .(block_id, seg_start, seg_end, chr, maj_cn, min_cn, tot_cn, quality_block) ]

# when you have LOH one allele will have no mutations. In this case set Varcount to 0. 
mut_baf_logR[ min_cn == 0, 
                `:=`(mean_hap1_varcount = ifelse( is.na(mean_hap1_varcount), 0, mean_hap1_varcount),
                     mean_hap2_varcount = ifelse( is.na(mean_hap2_varcount), 0, mean_hap2_varcount),
                     mean_hap1_varcount_cor = ifelse( is.na(mean_hap1_varcount_cor), 0, mean_hap1_varcount_cor),
                     mean_hap2_varcount_cor = ifelse( is.na(mean_hap2_varcount_cor), 0, mean_hap2_varcount_cor))]


mut_baf_logR <- mut_baf_logR[ !is.na(mean_hap1_varcount) & !is.na(mean_hap2_varcount), 
                                .( chr, block_start, block_end, seg_start, seg_end, block_id, maj_cn, min_cn, tot_cn, 
                                   n_hap1_muts, n_hap2_muts, quality_block, 
                                   baf = mean_hap1_varcount / (mean_hap1_varcount + mean_hap2_varcount),
                                   rdr = mean_hap1_varcount + mean_hap2_varcount,
                                   baf_norm = mean_hap1_varnorm / (mean_hap1_varnorm + mean_hap2_varnorm),
                                   rdr_norm = mean_hap1_varnorm + mean_hap2_varnorm,
                                   baf_cor = mean_hap1_varcount_cor / (mean_hap1_varcount_cor + mean_hap2_varcount_cor),
                                   rdr_cor = mean_hap1_varcount_cor + mean_hap2_varcount_cor ) ]
mut_baf_logR[, baf_inv := 1 - baf]
mut_baf_logR[, baf_consistent := abs(baf - 0.5)]
mut_baf_logR[, AI := abs((maj_cn - min_cn)/tot_cn)]
mut_baf_logR[, baf_norm_consistent := abs(baf_norm - 0.5)]
mut_baf_logR[, baf_cor_consistent := abs(baf_cor - 0.5)]


mut_baf_logR[, block_centre := (block_end - block_start)/2 + block_start ]
mut_baf_logR[, block_centre_cum := block_centre + chrlength.human_before_cum[ names(chrlength.human_before_cum) == chr ], 
               1:nrow(mut_baf_logR) ]
setnames(mut_baf_logR, c('seg_start', 'seg_end'), c('start', 'end'))
mut_baf_logR[, `:=`(Start.cummulitive = start + chrlength.human_before_cum[ names(chrlength.human_before_cum) == chr ], 
                    End.cummulitive = end + chrlength.human_before_cum[ names(chrlength.human_before_cum) == chr ]),
               1:nrow(mut_baf_logR) ]

# rdr & baf across genome
mut_baf_logR[, plot_rdr := ifelse( rdr_cor > 1200, 1200, rdr_cor) ]
plot_data <- mut_baf_logR[ , .(mean_rdr = mean(plot_rdr, na.rm=T)), by = .(Start.cummulitive, End.cummulitive)]

pdf(paste0(outputs.folder, date, '_rdr_cor_accross_genome.pdf'), width = 15, height = 5)
    ggplot() +
    geom_point(data = mut_baf_logR, shape = 21,
                aes(x = block_centre_cum, y = plot_rdr, fill = baf_consistent)) +
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
    ggplot(mut_baf_logR[ (quality_block) ]) +
    geom_histogram(aes(x = baf_cor_consistent, fill = chr ) ) +
    labs( x = '0.5-BAF' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

#look at hist of all BAF  -qual
pdf(paste0(outputs.folder, date, '_rdr_cor_hist_qual.pdf'), width = 9)
    ggplot(mut_baf_logR[ (quality_block) ]) +
    geom_histogram(aes(x = log2(rdr_cor), fill = chr ) ) +
    labs( x = 'log2 sum allele readcounts' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

## make a hatchet style plot
pdf(paste0(outputs.folder, date, '_hatchet_cor_quality.pdf'), width = 9)
    ggplot(mut_baf_logR[ (quality_block) ]) +
    geom_point(aes(x = baf_cor_consistent, y = log2(rdr_cor), colour = chr ), alpha = 0.5 ) +
    labs( x = '0.5-BAF', 
          y = 'log2 sum allele readcounts' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

## make a hatchet style plot + density
pdf(paste0(outputs.folder, date, '_hatchet_cor_quality_density.pdf'), width = 9)
    ggplot(mut_baf_logR[ (quality_block) ], aes(x = baf_consistent, y = log2(rdr))) +
    geom_point(aes( colour = chr ), alpha = 0.5 ) +
    stat_density_2d( alpha = 0.25,geom = "raster", contour = FALSE ) +
    labs( x = '0.5-BAF', 
          y = 'log2 sum allele readcounts' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

## make a hatchet style plot
pdf(paste0(outputs.folder, date, '_hatchet_cor_quality_3muts_perAl.pdf'), width = 9)
    ggplot(mut_baf_logR[ (quality_block) & nmuts_perAl > 3 ]) +
    geom_point(aes(x = baf_consistent, y = log2(rdr), colour = chr ), alpha = 0.5 ) +
    labs( x = '0.5-BAF', 
          y = 'log2 sum allele readcounts' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()

# compare baf and AI in quality bins
pdf(paste0(outputs.folder, date, '_rdr_cor_vs_totcn_quality.pdf'), width = 9)
    mut_baf_logR[ (quality_block), plot(log2(rdr_cor), tot_cn) ]
dev.off()

# compare baf and AI in quality bins
pdf(paste0(outputs.folder, date, '_baf_cor_vs_AI_quality.pdf'), width = 9)
    mut_baf_logR[ (quality_block), plot(AI, baf_cor_consistent) ]
dev.off()

# drdr across genome for quality bins
plot_data <- mut_baf_logR[ (quality_block) , .(mean_rdr = mean(rdr_cor, na.rm=T)), by = .(Start.cummulitive, End.cummulitive)]

pdf(paste0(outputs.folder, date, '_rdr_cor_accross_genome_quality.pdf'), width = 15, height = 5)
    ggplot() +
    geom_point(data = mut_baf_logR[(quality_block)], shape = 21,
                aes(x = block_centre_cum, y = rdr_cor, fill = baf_cor_consistent)) +
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

# drdr across genome for quality bins & limited baf
plot_data <- mut_baf_logR[ (quality_block) , .(mean_rdr = mean(rdr_cor, na.rm=T)), by = .(Start.cummulitive, End.cummulitive)]
mut_baf_logR[, baf_cor_consistent_lim := ifelse( baf_cor_consistent > 0.2, 0.2, baf_cor_consistent) ]

pdf(paste0(outputs.folder, date, '_rdr_cor_accross_genome_quality_baf_lim0.2.pdf'), width = 15, height = 5)
    ggplot() +
    geom_point(data = mut_baf_logR[(quality_block)], shape = 21,
                aes(x = block_centre_cum, y = rdr_cor, fill = baf_cor_consistent_lim)) +
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

plot_data <- mut_baf_logR[ (quality_block) & nmuts_perAl > 2 , .(mean_rdr = mean(rdr_cor, na.rm=T)), by = .(Start.cummulitive, End.cummulitive)]

pdf(paste0(outputs.folder, date, '_rdr_cor_accross_genome_quality_3muts_per_al.pdf'), width = 15, height = 5)
    ggplot() +
    geom_point(data = mut_baf_logR[(quality_block)  & nmuts_perAl > 2], shape = 21,
                aes(x = block_centre_cum, y = rdr_cor, fill = baf_consistent)) +
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
#### Try applying to low purity region ####
###########################################



### First extract read counts for all mutations called in high purity region in low purity region
extr_muts <- phased_muts[ (early_mut), .(chr = paste0('chr', chr), start = pos, end = pos, ref, alt) ]
extr_muts_gr <- GenomicRanges::makeGRangesFromDataFrame(extr_muts)
sbp = ScanBamParam(which=extr_muts_gr, isDuplicate=)

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
fst::write_fst(read_counts, paste(outputs.folder, date, '_R6_early_phased_muts_piledup_in_R3.fst'))


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
                                                    .(chr, pos, ref, alt, quality_block, 
                                                      block_id, phase, normal_depth, gc, Mcf7, 
                                                      ref_seg_start = seg_start, ref_seg_end = seg_end,
                                                      loh_ref_sample = min_cn == 0 ) ]) )

model = lm(varcount_tumour ~ splines::ns(x = gc, df = 5, intercept = T) + 
                             splines::ns(x = Mcf7, df = 5, intercept = T) + 
                             splines::ns(x = normal_depth, df = 5, intercept = T), y=F, model = F, data = test_data, na.action="na.exclude")
test_data[, varcount_tumour_corrected := residuals(model) + mean(varcount_tumour) ]
test_data[ varcount_tumour_corrected < 0, varcount_tumour_corrected := 0]

quasar_test <- test_data[, .( ref_seg_start, ref_seg_end, chr, loh_ref_sample, quality_block,
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
# Going to have to presume all mutations at ~0 varcount are clonal illusion (will include new LOH)
# The re-infer new LOH by looking for blocks where all mutations are from one haplotype - only works for long segments with >7-8 muts
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
                                   n_hap1_muts, n_hap2_muts, quality_block, nmuts = n_hap1_muts + n_hap2_muts,
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



quasar_test[, block_centre := (block_end - block_start)/2 + block_start ]
quasar_test[, block_centre_cum := block_centre + chrlength.human_before_cum[ names(chrlength.human_before_cum) == chr ], 
               1:nrow(quasar_test) ]
setnames(quasar_test, c('ref_seg_start', 'ref_seg_end'), c('start', 'end'))
quasar_test[, `:=`(Start.cummulitive = start + chrlength.human_before_cum[ names(chrlength.human_before_cum) == chr ], 
                    End.cummulitive = end + chrlength.human_before_cum[ names(chrlength.human_before_cum) == chr ]),
               1:nrow(quasar_test) ]

## why still getting random 0 blocks eg on chr3
quasar_test[ (quality_block) & chr == 1 & rdr_cor < 5]


# rdr & baf across genome -cor
plot_data <- quasar_test[ , .(mean_rdr = mean(rdr_cor, na.rm=T)), by = .(Start.cummulitive, End.cummulitive)]

pdf(paste0(outputs.folder, date, '_rdr_cor_accross_genome_test.pdf'), width = 15, height = 5)
    ggplot() +
    geom_point(data = quasar_test, shape = 21,
                aes(x = block_centre_cum, y = rdr_cor, fill = baf_cor_consistent)) +
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
plot_data <- quasar_test[ (quality_block & (!is.na(baf_cor_consistent) | is_mirror_vs_loh_ref)), .(mean_rdr = mean(rdr_cor, na.rm=T)), by = .(Start.cummulitive, End.cummulitive)]

pdf(paste0(outputs.folder, date, '_rdr_cor_accross_genome_quality_test.pdf'), width = 15, height = 5)
    ggplot() +
    geom_point(data = quasar_test[ (quality_block & (!is.na(baf_cor_consistent) | is_mirror_vs_loh_ref)) ], shape = 21,
                aes(x = block_centre_cum, y = rdr_cor, fill = baf_cor_consistent)) +
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
    geom_point(data = quasar_test[ (quality_block) & n_hap1_muts > 1 & n_hap2_muts > 1 ], shape = 21,
                aes(x = block_centre_cum, y = rdr_cor, fill = baf_cor_consistent)) +
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
    geom_point(data = quasar_test[ (quality_block) ], shape = 21,
                aes(x = block_centre_cum, y = rdr, fill = baf_consistent)) +
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
    ggplot(mut_baf_logR[ (quality_block) ]) +
    geom_point(aes(x = baf_cor_consistent, y = log2(rdr_cor), colour = chr ), alpha = 0.5 ) +
    labs( x = '0.5-BAF', 
          y = 'log2 sum allele readcounts' ) +
    theme_bw() +
    theme( text = element_text( size = 20 ))
dev.off()






#############
#### END ####
#############
