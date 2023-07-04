#==============================================================================#
#==============================================================================#
######                                                                    ######
######    Script select heterozygous SNPs nearby mutation to perform      ######
######             mpileup and therefore phase the mutations              ######
######                                                                    ######
#==============================================================================#
#==============================================================================#

# Author: Alex Frankell
# Date: 2023-06-29

####################################################
#### Source required functions & load libraries ####
####################################################

library(data.table) 
library(GenomicRanges) 

##################################
#### Interpret following Args ####
##################################

cmdArgs    <- commandArgs(trailingOnly = TRUE)
mut_bed_path  <- cmdArgs[1]
snp_bed_path  <- cmdArgs[2]
num_snps      <- as.numeric( cmdArgs[3] )
outfile       <- cmdArgs[4]

## Testing 
#mut_bed_path  <- 'outputs/phasing/20230627/scratch/LTX0064_SU_T1-R3.snvs.bed'
#snp_bed_path  <- 'outputs/phasing/20230627/scratch/LTX0064.hetsnps.bed'
#num_snps      <- as.numeric( '6' )
#outfile       <- 'outputs/phasing/20230627/scratch/LTX0064.hetsnps.selected.bed'

#############################
#### Get inputs required ####
#############################

muts <- data.table::fread( mut_bed_path )
setnames(muts, c('chr', 'start'))
muts[, end := start ]
muts <- GenomicRanges::makeGRangesFromDataFrame(muts)

snps <- data.table::fread( snp_bed_path )
setnames(snps, c('chr', 'start'))
snps[, end := start ]
snps <- GenomicRanges::makeGRangesFromDataFrame(snps)

####################################
#### Extract the nearest X SNPs ####
####################################

# Ensure both ordered by chr and position
muts <- muts[ order(muts) ]
snps <- snps[ order(snps) ]

# speed up by doing per chromosome
chrs <- unique(seqnames(snps))

# extract the cummulative max index per chromosome (sum of previous chrs)
# to add to per chromosome indices 
chr_max_indices_cum <- sapply(1:length(chrs), function(i){
    if(i == 1) return(0)
    chrs_sel <- chrs[1:(i-1)]
    return( length( snps[ seqnames(snps) %in% chrs_sel, ] ) )
} )
names( chr_max_indices_cum ) <- chrs

# extract indicies per chr, correct to full table indices and bind together
selected_snps <- unlist( lapply( chrs, function(chr, num_sel = num_snps){
        
        # Subset to Chrs
        muts_chr <- muts[ seqnames(muts) == chr, ]
        if( length( muts_chr ) == 0 ) next
        snps_chr <- snps[ seqnames(snps) == chr, ]

        # Exract index of nearest hetSNP to each mut in SNP file
        SNPindices <- GenomicRanges::nearest(muts_chr, snps_chr)

        # Exract specified number of SNPs up and downstream  of nearest hetSNP
        # already ordered by pos 
        max_index <- length(snps_chr)
        SNPindices <- unlist( lapply(SNPindices, function(index){
            
            out <- (index-round(num_sel/2)):(index+round(num_sel/2))
            
            #check this is within reasonable boundaries
            # if not ensure you still select the correct number 
            if( any(out < 1) ){
                len_removed <- sum(out < 1)
                out <- out[ !out < 1 ]
                out <- c(out, (max(out)+1):(max(out)+1+len_removed))
            }
            if( any(out > max_index) ){
                len_removed <- sum(out > max_index)
                out <- out[ !out > max_index ]
                out <- c((max(out)-1-len_removed):(max(out)-1), out)
            }
            
            return(out)
        }))

        # remove dups if any mutations clones enough together
        SNPindices <- unique( SNPindices )

        #convert indices to those in full table before subsetting by chr
        SNPindices <- SNPindices + chr_max_indices_cum[ names(chr_max_indices_cum) == chr ]

        return(SNPindices)
}))

out_hetSNPs <- snps[ selected_snps ]

###################
#### Write out ####
###################

out <- as.data.table(out_hetSNPs)

data.table::fwrite(out[, .(seqnames, start)], file = outfile,
                   col.names = FALSE, sep = '\t')

#############
#### END ####
#############
