#!/usr/bin/env Rscript

#SBATCH --job-name=ctDNA_overlay
#SBATCH --time=2:00:00 
#SBATCH --mem=10G
#SBATCH --qos=normal
#SBATCH --out='plots/ctDNA_overlay/command_history/%x.%j.out'    ### may need to modify if moving script
#SBATCH --error='plots/ctDNA_overlay/command_history/%x.%j.err'  ### may need to modify if moving script

# print date in output file & error file so it can be matched back to output #

#=========================================================#
#  ml R/3.6.0-foss-2019 before lunching batch submission  #
#=========================================================#

print( paste0(Sys.time(), "\n\n\n") )
warning( paste0(Sys.time(), "\n\n\n") )


#========================================================================================#
#========================================================================================#
######                                                                              ######
######  Script bind together all the PCAWG data into manageable cohort-level files  ######
######                                                                              ######
#========================================================================================#
#========================================================================================#

# Author: Alex Frankell
# Date: 2021-06-22

setwd("/Volumes/proj-tracerx-lung/tctProjects/frankella/mrd_methods_dev/")

####################################################
#### Source required functions & load libraries ####
####################################################

# suppress warning on R build version #
library(vcfR)
library(data.table) 
library(fst) 

#############################################
#### Make a folder for this analysis run ####
#############################################

date <- gsub("-","",Sys.Date())

# the output folder paths will already exist if specified in following args, if not use default specified here #

outputs.folder <- "outputs/PCAWG_calls/"  ## this folder must already exist  

outputs.folder <- paste0( outputs.folder, date, "/" )

system( paste0('mkdir ', outputs.folder ) )


#############################
#### Read in all the SVs ####
#############################

files.dir <- 'inputs/PCAWG_sv_calls/'
paths <- paste0(files.dir, list.files( files.dir ))
paths <- paths[ !grepl('.md5$', paths) ]

SVs <- lapply(paths[1:100], function(path){
  print(which(paths==path))
  SV <- fread(path)
  id <- gsub( files.dir, '', tstrsplit(path, split = '.pcawg')[[1]] )
  SV[, sample := id ]
  return(SV)
})
SVs <- rbindlist( SVs )

SVs2 <- lapply(paths[101:800], function(path){
  print(which(paths==path))
  SV <- fread(path)
  id <- gsub( files.dir, '', tstrsplit(path, split = '.pcawg')[[1]] )
  SV[, sample := id ]
  return(SV)
})
SVs2 <- rbindlist( SVs2 )

SVs3 <- lapply(paths[801:1400], function(path){
  print(which(paths==path))
  SV <- fread(path)
  id <- gsub( files.dir, '', tstrsplit(path, split = '.pcawg')[[1]] )
  SV[, sample := id ]
  return(SV)
})
SVs3 <- rbindlist( SVs3 )

SVs4 <- lapply(paths[1401:length(paths)], function(path){
  print(which(paths==path))
  SV <- fread(path)
  id <- gsub( files.dir, '', tstrsplit(path, split = '.pcawg')[[1]] )
  SV[, sample := id ]
  return(SV)
})
SVs4 <- rbindlist( SVs4 )

SVs <- rbindlist( list(SVs, SVs2, SVs3, SVs4) )

write_fst(SVs, paste0( outputs.folder, date, '_PCAWG_SVs.fst'))


#############
#### END ####
#############
