#==============================================================================#
#==============================================================================#
######                                                                    ######
######  Script plot LOD for samples and relationships with other metrics  ######
######                                                                    ######
#==============================================================================#
#==============================================================================#

# Author: Alex Frankell
# Date: 2022-09-27

setwd("/Volumes/proj-tracerx-lung/tctProjects/frankella/ctDNA_methods_dev//")

####################################################
#### Source required functions & load libraries ####
####################################################

# suppress warning on R build version #
library(fst)
library(data.table) 
library(ggplot2) 
library(RColorBrewer) 


#############################################
#### Make a folder for this analysis run ####
#############################################

date <- gsub("-","",Sys.Date())

# the output folder paths will already exist if specified in following args, if not use default specified here #

outputs.folder <- "outputs/fixed_panel/"  ## this folder must already exist  

outputs.folder <- paste0( outputs.folder, date, "/" )

system( paste0('mkdir ', outputs.folder ) )


##############################################
#### Get Inputs required for all analyses ####
##############################################

SVcalls_path        <- "outputs/PCAWG_calls/20220927/20220927_PCAWG_SVs.fst"
pacwg_samples_path  <- "inputs/pcawg_sample_sheet.tsv"

# read in variant level table
SVcalls <- fst::read_fst( SVcalls_path, as.data.table=T )

# read in sample level table
pacwg_samples <- data.table::fread( pacwg_samples_path )
# overlay tumours types
SVcalls[, tumour_type := pacwg_samples[ match(sample, aliquot_id), dcc_project_code ] ]


#########################################################
#### In which tumour types is this at all realistic? ####
#########################################################

# Need probably at least 40 SVs in >80% of tumours
SV_counts <- SVcalls[, .(TSVB = .N), by = .(sample, tumour_type) ]
tumour_types <- SV_counts[, .(frac_over_25 = sum(TSVB>30)/.N*100), by = tumour_type ]
tumour_types[, tumour_type := factor(tumour_type, levels = tumour_types[ order(frac_over_25, decreasing = T), tumour_type])]

cols <- rep(brewer.pal(12, "Paired"), 3)[1:tumour_types[, .N]]

pdf( paste(outputs.folder, date, '_tumour_types_with_enough_SVs'))

ggplot(tumour_types, aes(x = tumour_type, y = frac_over_25) ) +
  geom_col( aes(fill = tumour_type )) +
  scale_fill_manual( values = cols ) +
  geom_hline( yintercept = 80, linetype = 2 ) +
  scale_y_continuous( breaks = seq(0, 100, 20)) +
  labs( x = '',
        y = '% of tumours with > 40 SVs') +
  theme_classic() +
  theme( text = element_text( size = 20),
         axis.text.x = element_text( angle = 45, hjust = 1 ),
         legend.position='none')

dev.off()

###############################################################
#### Look at the top fragile site areas - frequency of SVs ####
###############################################################

# Lets try a 500Kb panel around the three biggest Fragile sites
SVcalls[, FRA3B_centre := (chrom1 == 3 & start1 > 60400000 & end1 < 60550000) | (chrom2 == 3 & start2 > 64000000 & end2 < 65500000) ]
SVcalls[, FRA320B_centre := (chrom1 == 20 & start1 > 14850000 & end1 < 15050000) | (chrom2 == 20 & start2 > 14850000 & end2 < 15050000) ]
SVcalls[, FRA16D_centre := (chrom1 == 16 & start1 > 78675000 & end1 < 78825000) | (chrom2 == 16 & start2 > 78675000 & end2 < 78825000) ]

60550000 - 60400000 + 15050000 - 14850000 + 78825000 - 78675000
# [1] 500000

# 4 Mb panel - top 8 FSs (500Kb each)
SVcalls[, FRA3B_centre := (chrom1 == 3 & start1 > 60200000 & end1 < 60700000) | (chrom2 == 3 & start2 > 60200000 & end2 < 60700000) ]
SVcalls[, FRA320B_centre := (chrom1 == 20 & start1 > 14750000 & end1 < 15250000) | (chrom2 == 20 & start2 > 14750000 & end2 < 15250000) ]
SVcalls[, FRA16D_centre := (chrom1 == 16 & start1 > 78350000 & end1 < 788500000) | (chrom2 == 16 & start2 > 78350000 & end2 < 788500000) ]
SVcalls[, IMMP2L_centre := (chrom1 == 7 & start1 > 110800000 & end1 < 111300000) | (chrom2 == 16 & start2 > 110800000 & end2 < 111300000) ]
SVcalls[, NAALADL2_centre := (chrom1 == 3 & start1 > 174750000 & end1 < 175250000) | (chrom2 == 16 & start2 > 174750000 & end2 < 175250000) ]
SVcalls[, LRP1B_centre := (chrom1 == 2 & start1 > 140750000 & end1 < 141250000) | (chrom2 == 16 & start2 > 140750000 & end2 < 141250000) ]
SVcalls[, PDE4D_centre := (chrom1 == 5 & start1 > 59500000 & end1 < 60000000) | (chrom2 == 16 & start2 > 59500000 & end2 < 60000000) ]
SVcalls[, CCSER1_centre := (chrom1 == 4 & start1 > 90600000 & end1 < 91100000) | (chrom2 == 16 & start2 > 90600000 & end2 < 91100000) ]

#do this properly...
FSs <- as.data.table( readxl::read_excel('~/Downloads/FSs.xlsx') )
FSs[, chr := gsub('Chr', '', tstrsplit(`Genomic co-ordinates`, split = ':')[[1]]) ]
FSs[, pos := gsub('Chr', '', tstrsplit(`Genomic co-ordinates`, split = ':')[[2]]) ]
FSs[, start := as.numeric(tstrsplit(pos, split = '-')[[1]]) ]
FSs[, end := as.numeric(tstrsplit(pos, split = '-')[[2]]) ]
FSs[, size := end - start ]
FSs[, centre :=  start + size/2 ]

# annotate number of SVs in total
FSs[, PCAWG_SVs_total := SVcalls[ (chrom1 == chr & start1 > start & end1 < end) | (chrom2 == chr & start2 > start & end2 < end), .N ], 1:nrow(FSs) ]
FSs[, PCAWG_SVs_EAC := SVcalls[ tumour_type == 'ESAD-UK' & (chrom1 == chr & start1 > start & end1 < end) | (chrom2 == chr & start2 > start & end2 < end), .N ], 1:nrow(FSs) ]


window <- 1000000
window <- 1000000
SVcalls[, in_panel := FALSE ]
freq_threshold <- 220
# panel size, Mb
FSs[ PCAWG_SVs_total > freq_threshold, .N * window/1000000]
# number sites
FSs[ PCAWG_SVs_total > freq_threshold, .N ]

for( rowi in which(FSs$PCAWG_SVs_total > freq_threshold) ){
  SVcalls[, in_panel := (chrom1 == FSs$chr[rowi] & start1 > FSs$centre[rowi] - window/2  & end1 < FSs$centre[rowi] + window/2) | (chrom2 == FSs$chr[rowi] & start2 > FSs$centre[rowi] - window/2 & end2 < FSs$centre[rowi] + window/2) | in_panel ]
}

#SVcalls[, in_panel := FRA3B_centre | FRA320B_centre | FRA16D_centre | IMMP2L_centre | NAALADL2_centre | PDE4D_centre | LRP1B_centre | CCSER1_centre ]

panel <- SVcalls[, .(panel_SVs = sum(in_panel)), by = .(sample, tumour_type) ]
type_ranks <- panel[, .(mean_svs_in_panel = mean(panel_SVs)), by = tumour_type]
panel[, tumour_type := factor(tumour_type, levels = type_ranks[ order(mean_svs_in_panel, decreasing = T), tumour_type])]

panel[tumour_type== 'ESAD-UK', median(panel_SVs)]

#pdf( paste(outputs.folder, date, '_SV_freq_in_top_3_FRA_hotspots_500kb.pdf'), width = 8)

ggplot(panel, aes(x = panel_SVs)) +
  geom_histogram( binwidth = 0.5, aes( fill = tumour_type ) ) +
  scale_fill_manual( values = cols ) +
  labs( x = ' # of SVs in 500Kb fragile site panel (FHIT, WWOX & MACROD2)',
        y = 'Frequency' ) +
  theme_classic() +
  theme( text = element_text( size = 13),
        legend.position = 'none') +
  facet_wrap( ~ tumour_type, scales = 'free')

dev.off()

### This is close to working for EAC but even then not really (<50% of patients 
### have at least 3 SVs which could get to ~0.02% high sen LOD) and not close to 
### working at all for other cancer types

##############################
#### Look at larger panels ####
##############################

# This might be possible if we deplete for references sequences leaving/enriching for Chymeric fragments  
# probably can do 100 Mb

# Lets try a ~8Mb panel around the three biggest Fragile sites
SVcalls[, FRA3B_centre := (chrom1 == 3 & start1 > 60000000 & end1 < 61000000) | (chrom2 == 3 & start2 > 60000000 & end2 < 61000000) ]
SVcalls[, FRA320B_centre := (chrom1 == 20 & start1 > 145000000 & end1 < 15500000) | (chrom2 == 20 & start2 > 145000000 & end2 < 15500000) ]
SVcalls[, FRA16D_centre := (chrom1 == 16 & start1 > 78300000 & end1 < 79300000) | (chrom2 == 16 & start2 > 78300000 & end2 < 79300000) ]

SVcalls[, in_panel := FRA3B_centre | FRA320B_centre | FRA16D_centre ]

panel <- SVcalls[, .(panel_SVs = sum(in_panel)), by = .(sample, tumour_type) ]
type_ranks <- panel[, .(mean_svs_in_panel = mean(panel_SVs)), by = tumour_type]
panel[, tumour_type := factor(tumour_type, levels = type_ranks[ order(mean_svs_in_panel, decreasing = T), tumour_type])]

pdf( paste(outputs.folder, date, '_SV_freq_in_top_3_FRA_hotspots_500kb.pdf'), width = 8)

ggplot(panel, aes(x = panel_SVs)) +
  geom_histogram( binwidth = 0.5, aes( fill = tumour_type ) ) +
  scale_fill_manual( values = cols ) +
  labs( x = ' # of SVs in 8Mb fragile site panel (FHIT, WWOX & MACROD2)',
        y = 'Frequency' ) +
  theme_classic() +
  theme( text = element_text( size = 13),
         legend.position = 'none') +
  facet_wrap( ~ tumour_type, scales = 'free')

dev.off()


### How many SVs could be captured using sequenctail enrichment?



window <- 10000000
SVcalls[, enrichable := FALSE ]
freq_threshold <- 275
# panel size, Mb
FSs[ PCAWG_SVs_total > freq_threshold, .N * window/1000000]
# number sites
FSs[ PCAWG_SVs_total > freq_threshold, .N ]

for( rowi in which(FSs$PCAWG_SVs_total > freq_threshold) ){
  SVcalls[, enrichable := ((chrom1 == FSs$chr[rowi] & start1 > FSs$centre[rowi] - window/2  & end1 < FSs$centre[rowi]) & (chrom2 == FSs$chr[rowi] & start2 > FSs$centre[rowi] & end2 < FSs$centre[rowi] + window/2)) | enrichable ]
}

panel <- SVcalls[, .(panel_SVs = sum(enrichable)), by = .(sample, tumour_type) ]
type_ranks <- panel[, .(mean_svs_in_panel = mean(panel_SVs)), by = tumour_type]
panel[, tumour_type := factor(tumour_type, levels = type_ranks[ order(mean_svs_in_panel, decreasing = T), tumour_type])]

panel[tumour_type== 'ESAD-UK', median(panel_SVs)]

#pdf( paste(outputs.folder, date, '_SV_freq_in_top_3_FRA_hotspots_500kb.pdf'), width = 8)

ggplot(panel, aes(x = panel_SVs)) +
  geom_histogram( binwidth = 0.5, aes( fill = tumour_type ) ) +
  scale_fill_manual( values = cols ) +
  labs( x = ' # of SVs in 500Kb fragile site panel (FHIT, WWOX & MACROD2)',
        y = 'Frequency' ) +
  theme_classic() +
  theme( text = element_text( size = 13),
         legend.position = 'none') +
  facet_wrap( ~ tumour_type, scales = 'free')

dev.off()





#########################################
#### Design 500Kb panel from scratch ####
#########################################






#############
#### END ####
#############
