#==============================================================================#
#==============================================================================#
######                                                                    ######
######  Script plot LOD for samples and relationships with other metrics  ######
######                                                                    ######
#==============================================================================#
#==============================================================================#

# Author: Alex Frankell
# Date: 2021-06-22

setwd("/Volumes/proj-tracerx-lung/tctProjects/frankella/archer_ctdna/")

####################################################
#### Source required functions & load libraries ####
####################################################

# suppress warning on R build version #
library(fst)
library(data.table) 
library(dplyr) 
library(amfFunctions) 
library(ggplot2) 
library(cowplot) 
library(cloneMap) 
library(fishplot) 
library(RColorBrewer) 

#############################################
#### Make a folder for this analysis run ####
#############################################

date <- gsub("-","",Sys.Date())

# the output folder paths will already exist if specified in following args, if not use default specified here #

outputs.folder <- "plots/template/"  ## this folder must already exist  

outputs.folder <- paste0( outputs.folder, date, "/" )

system( paste0('mkdir ', outputs.folder ) )


##############################################
#### Get Inputs required for all analyses ####
##############################################

# these input file paths will already exist if speciified in following args, if not use defaults specified here #

ctDNA_data_path          <- "../../lungTx/Tx421/release/ctDNA_archer/20220202/20220308_tracked_mutations_primary_and_met_data_eclipse_annotated_short.fst"
ctDNA_samples_path       <- "../../lungTx/Tx421/release/ctDNA_archer/20220202/20220308_output_main_sample_table_primary__onlyannotated.tsv"

# read in variant level table
ctDNA_data <- fread( ctDNA_data_path )

# read in sample level table
ctDNA_samples <- fread( ctDNA_samples_path )


####################
#### Analysis 1 ####
####################

### Look at clone detection

# Calculate the limit of detection for each subclone 

all_patient_ids <- c('LTX474', 'LTX208')

for( patient in all_patient_ids ){

pdf( paste0(outputs.folder, date, '_', patient, '_name_of_file.pdf') )

ggplot( ctDNA_data[ patient_name == patient ], aes() )...

dev.off()

}

# cowplot
plot1 <- ggplot( ctDNA_data[ patient_name == patient ], aes() )...
plot2 <- ggplot( ctDNA_data[ patient_name == patient ], aes() )...
combined_plot <- combine_plot( plot1, plot2, ncol = 1) ##not the right function name!!! but from cowplot


####################
#### Analysis 2 ####
####################


pdf( paste0(outputs.folder, date, '_name_of_file2.pdf') )

ggplot()...

dev.off()


#############
#### END ####
#############
