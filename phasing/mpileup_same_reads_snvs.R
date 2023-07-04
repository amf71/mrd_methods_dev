library(vcfR)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(Biostrings)
library(gtools)
library(data.table)

#Script to extract mpileup same read somatic mutation phasing
#Load sample of interest
Args <- commandArgs(trailingOnly=TRUE)[1]
T_sample <- Args[1]
N_sample <- Args[2]

# parameter (could specify in parent script)
# The max distance in bp to bother trying to phase a hetsnp to a mutation
# Shouldn't make much diff as usually hetSNP every Kb anyway
max_read_length <- 20000

# testing
T_sample="LTX0064_SU_T1-R3"
N_sample="LTX0064_SU_N"

print(T_sample)
print(N_sample)


#Get functions---------------------------------------------
#Get a function to group Indels in reads from the mpileup output:
#Base before the +/- symbol followed by n bases  
group_indels <- function(read_bases_vec, pos){
  #Read how many positions need to be taken into account
  #Add following element in the vector if the integer is more than 1 digit
  if (suppressWarnings(!is.na(as.numeric(read_bases_vec[pos + 2])))){
    nb_bases <- as.numeric(paste(read_bases_vec[c(pos + 1, pos + 2)], collapse = ""))
    nb_bases <- nb_bases + 1
  }else {
    nb_bases <- as.numeric(read_bases_vec[pos + 1])
  }
  
  #Get the relevant positions into a single string (base before the +- and number of bases affected)
  string_indel <- read_bases_vec[(pos - 1):(pos + nb_bases + 1)]
  string_indel <- paste(string_indel, collapse = "")
  
  #Replace the affected positions 
  read_bases_out <- read_bases_vec[pos-1] <- string_indel
  read_bases_out <- read_bases_vec[-(pos:(pos + nb_bases + 1))]
  return(read_bases_out)
}

#Get a function to keep only reads with the mutant allele and good mapping quality
keep_reads <- function(mpileup_row, mapq_thres = 10, mode = "alt"){
  #Skip if no read counts were found
  if (as.numeric(mpileup_row["counts"]) > 0){
    #Get a vector with the base reads
    read_bases <- toupper(as.character(mpileup_row["read_bases"]))
    read_bases <- unlist(strsplit(read_bases, ""))
    
    #Deal with "^", they mean the base is the first base in the read, merge into a single string:
    #Previous base, followed by the ^ symbol, followed by ASCII phred score symbol
    if (any(read_bases == "^")){
      first_base_nb <- sum(read_bases == "^")
      for (this_first_base in 1:first_base_nb){
        first_base_pos <- which(read_bases == "^")
        this_pos_first_base <- first_base_pos[1]
        string_first_base <- paste(read_bases[(this_pos_first_base - 1):(this_pos_first_base + 1)], collapse = "")
        #Replace the affected positions 
        read_bases[this_pos_first_base-1] <- string_first_base
        read_bases <- read_bases[-((this_pos_first_base):(this_pos_first_base + 1))]
      }
    }
    
    #Deal with "$", they mean the base is the lastt base in the read, merge into a single string:
    #Previous base, followed by the $ symbol
    #Any $ symbol part of ASCII phred scores will have been taken care of in the previous step
    if (any(read_bases == "$")){
      last_base_nb <- sum(read_bases == "$")
      for (this_last_base in 1:last_base_nb){
        last_base_pos <- which(read_bases == "$")
        this_pos_last_base <- last_base_pos[1]
        string_last_base <- paste(read_bases[(this_pos_last_base - 1):this_pos_last_base], collapse = "")
        #Replace the affected positions 
        read_bases[this_pos_last_base-1] <- string_last_base
        read_bases <- read_bases[-((this_pos_last_base):(this_pos_last_base))]
      }
    }
    
    #Any + or - symbol part of ASCII phred scores will have been taken care of in the previous steps
    if (any(read_bases %in% c("+", "-"))){
      #Group the base before the +- symbol, and the following n bases into a single string
      indel_nb <- sum(read_bases %in% c("+", "-"))
      for (this_indel in 1:indel_nb){
        indel_pos <- which(read_bases %in% c("+", "-"))
        this_pos <- indel_pos[1]
        read_bases <- group_indels(read_bases, this_pos)
      }
    }
    
    #Get the mapping qualities
    map_q <- as.character(mpileup_row["map_qual"])
    map_q <- as.numeric(unlist(strsplit(map_q, ",")))
    
    #Get the positions of the reads covering the allele of interest
    if (mode == "alt"){
      match_allele <- as.character(mpileup_row["alt"])
    } else if (mode == "ref"){
      match_allele <- c(".", ",", "^.", "^,", "$.", "$,") 
    } else {
      stop("Mode has to be either 'ref' or 'alt'")
    }
    
    read_ids_pos <- which(read_bases %in% match_allele & map_q >= mapq_thres)
    
    if (length(read_ids_pos) > 0){
      #Get the reads of interest (if any)
      read_ids <- as.character(mpileup_row["reads_id"])
      read_ids <- unlist(strsplit(read_ids, ","))
      read_ids <- read_ids[read_ids_pos]
    } else {
      read_ids <- ""
    } # Close if position lengt == 0 else/if
  } else {
    read_ids <- ""
  } # Close if no coverage, skip else/if
  return(paste(read_ids, collapse = ";"))
}

#Turn read identification results into long format
set_reads_ids_long <- function(df_mpileup, col){
  read_ids_list <- strsplit(df_mpileup[, col], split = ";")
  
  mpileup_results_long <- data.frame(chr = rep(df_mpileup$chr, sapply(read_ids_list, length)),
                                     pos = rep(df_mpileup$pos, sapply(read_ids_list, length)),
                                     ref = rep(df_mpileup$ref, sapply(read_ids_list, length)),
                                     alt = rep(df_mpileup$alt, sapply(read_ids_list, length)),
                                     counts = rep(df_mpileup$counts, sapply(read_ids_list, length)),
                                     read_ids = unlist(read_ids_list),
                                     stringsAsFactors = FALSE)
  
  mpileup_results_long$mutid <- paste(mpileup_results_long$chr,
                                      mpileup_results_long$pos,
                                      mpileup_results_long$ref,
                                      mpileup_results_long$alt,
                                      sep = ":")
  return(mpileup_results_long)
}


#Get a function to extract the mutation patterns per read
get_mut_profile <- function(read_mpileup_df){
  #Add a column with whether or not the mutation is present with its id
  read_mpileup_df$mutid_present <- paste(read_mpileup_df$mutid,
                                         ifelse(read_mpileup_df$mut_present, 1, 0),
                                         sep = ":")
  joint_muts <- paste(mixedsort(read_mpileup_df$mutid_present), collapse = ";")
  return(joint_muts)
}

#Function to merge reads with missing mutations
replace_missing_mut <- function(groupid_mut, all_mutids_df){
  
  #Get the number of reads supporting the focal groupid
  initial_depth <- all_mutids_df$coverage[all_mutids_df$mutids_read == groupid_mut]
  #Remove focal  group
  all_mutids_df <- all_mutids_df[all_mutids_df$mutids_read != groupid_mut, ]
  
  #Get group ids only
  all_mutids <- all_mutids_df$mutids_read
  
  
  #Get the position of the first and last mutation in the group
  first_mut <- gsub(pattern = ":[TCGA]:[TCGA]:[01];.+",
                    replacement = "",
                    x = groupid_mut)
  last_mut <- gsub(pattern = ".+;",
                   replacement = "",
                   x = groupid_mut)
  
  last_mut <- gsub(pattern = ":[TCGA].+",
                   replacement = "",
                   x = last_mut)
  
  #Get the chromosome for the first mutation, keep only mutation groups from that chromosome
  #Exclude rare, but occurring instances of reads spanning multiple chromosmes)
  chr_first <- gsub(pattern = ":.+",
                    replacement = "",
                    x = first_mut)
  
  chr_last <- gsub(pattern = ":.+",
                   replacement = "",
                   x = last_mut)
  
  if(chr_last == chr_first){
    mut_groups_chr_pattern <- gsub(pattern = ":.+",
                                   replacement = "",
                                   x = first_mut)
    mut_groups_chr_pattern <- paste0("^", mut_groups_chr_pattern,
                                     ".+;", mut_groups_chr_pattern, ":[0-9]+:[ATCG]:[ATCG]:[01]$")
    
    mut_groups_chr <- grep(mut_groups_chr_pattern, all_mutids, value = TRUE)
    
    #Get groups contained within the current group
    start_pos_group <- as.numeric(gsub(pattern = ".+:",
                                       replacement = "",
                                       x = first_mut))
    
    end_pos_group <- as.numeric(gsub(pattern = ".+:",
                                     replacement = "",
                                     x = last_mut))
    
    mut_groups_start_pos <- as.numeric(gsub(pattern = "(^chr([0-9]+|[XY]):)([0-9]+)(:[ATCG]:.+;.+$)",
                                            replacement = "\\3",
                                            x = mut_groups_chr))
    
    mut_groups_end_pos <- as.numeric(gsub(pattern = "(.+chr([0-9]+|[XY]):)([0-9]+)(:[ATCG]:.+$)",
                                          replacement = "\\3",
                                          x = mut_groups_chr))
    
    #Get the groups that are contained within the focal group
    overlap_groups <- mut_groups_chr[start_pos_group <= mut_groups_start_pos &
                                       end_pos_group >= mut_groups_end_pos]
    
    #For each overlapping group, select the ones that can be merged into the larger one
    groupid_mut_split <- unlist(strsplit(groupid_mut, ";"))
    
    total_reads_add <- 0
    merged_groups <- character()
    
    #Skip if no groups overlap
    if(length(overlap_groups) > 0){
      
      for(this_group in overlap_groups){
        this_group_split <-  unlist(strsplit(this_group, ";"))
        
        #To merge if all combinations of mutations are contained within the focal group
        to_merge <- all(this_group_split %in% groupid_mut_split)
        
        #If merged, add the proportion equivalent to the proportion of overlapping mutations in common
        if (to_merge){
          reads_this_group <- all_mutids_df$coverage[all_mutids_df$mutids_read == this_group]
          add_reads <- reads_this_group * sum(this_group_split %in% groupid_mut_split) / length(groupid_mut_split)
          total_reads_add <- total_reads_add + add_reads
          
          #Keep track of the groups merged
          merged_groups <- c(merged_groups, this_group)
        } #End if merging
        
      } #End merging loop
      
      #Get the total number of "reads" per merged group
      total_reads <- total_reads_add + initial_depth
      
      #Get output together
      merged_groups_char <- paste(merged_groups, collapse = "-")
      merged_out <- c(groupid_mut, merged_groups_char, total_reads)
      
    } else { #End empty overlapping groups if 
      merged_out <-  c(groupid_mut, "", initial_depth)
      
    } # End empty overlapping groups else
    
  } else { #End different chromosomes if 
    merged_out <-  c(groupid_mut, "", initial_depth)
    
  } #End different chromosomes else
  
  return(merged_out)
} 

#Load the mpileup results showing reads shared across mutations-----------------------------------------------------------------------
mpileup_file_out <- paste0("outputs/phasing/20230627/", T_sample, ".pileup.txt")
mpileup_results <- fread(mpileup_file_out,
                         sep = "\t",
                         na.strings = NULL,
                         quote = "",
                         data.table = FALSE)

colnames(mpileup_results) <- c("chr", "pos", "ref", "counts", "read_bases", "qual", "reads_id", "map_qual")

#Add the alternative allele from the original VCFs
#snvs
input_vcf_mpileup <- paste0("outputs/phasing/20230627/scratch/", T_sample, "_filt_snvs_mpileup.vcf")
somatic_vcf <- read.vcfR(input_vcf_mpileup)

#snps
patient <- strsplit(T_sample, split = '_')[[1]][1]
input_vcf_mpileup <- paste0("outputs/phasing/20230627/scratch/", patient, "_filt_hetsnps_mpileup.vcf")
hetsnps_vcf <- read.vcfR(input_vcf_mpileup)

#Get the reference and alternative allele for each mutation
all_muts <- data.frame(chr = getCHROM(somatic_vcf),
                       pos = getPOS(somatic_vcf),
                       ref = getREF(somatic_vcf),
                       alt = getALT(somatic_vcf),
                       type = 'mut',
                       stringsAsFactors = FALSE)

#Get the reference and alternative allele for each het snp
all_hetsnps <- data.frame(chr = getCHROM(hetsnps_vcf),
                          pos = getPOS(hetsnps_vcf),
                          ref = getREF(hetsnps_vcf),
                          alt = getALT(hetsnps_vcf),
                          type = 'het_snp',
                          stringsAsFactors = FALSE)

#Merge
mpileup_results$ref <- toupper(mpileup_results$ref)
mpileup_results <- inner_join(mpileup_results, rbind(all_muts, all_hetsnps) )

#Get reads containing variant allele
print("Get mutations per read, alt")
mpileup_results$mut_read <- apply(mpileup_results, 1, function(x) keep_reads(x, 10, "alt"))

#Get reads containing the reference allele
print("Get mutations per read, wt")
mpileup_results$wt_read <- apply(mpileup_results, 1, function(x) keep_reads(x, 10, "ref"))

# Only keep muts for de novo phasing of blocks
mpileup_results_muts <- mpileup_results[ type == 'mut' ]

#Keep only variants with coverage
mpileup_results_filt <- mpileup_results_muts[mpileup_results_muts$mut_read != "" &
                                             mpileup_results_muts$wt_read != "", ]

#Get into long format, one row per read for wt and mut reads
print("Turn into long format mut")
mpileup_long_mut <- set_reads_ids_long(mpileup_results_filt, "mut_read")
#Set mutation detected
mpileup_long_mut$mut_present <- TRUE

print("Turn into long format wt")
mpileup_long_wt <- set_reads_ids_long(mpileup_results_filt, "wt_read")
mpileup_long_wt$mut_present <- FALSE

#Join
mpileup_results_long <- rbind(mpileup_long_mut, mpileup_long_wt)

#Remove duplicated reads (for some reason)
mpileup_results_long <- unique(mpileup_results_long)

#Keep raw mpileup output for QC purposes
mpileup_results_raw <- mpileup_results_long

#Keep only reads that appear more than once (covering more than one mutation)
mpileup_results_long <- mpileup_results_long %>% group_by(read_ids) %>%
  filter(n_distinct(mutid) > 1)

#Count reads supporting different mutation combinations------------------------------
#Get unique reads
read_ids <- unique(mpileup_results_long$read_ids)

#Per read, get the mutations present/absent
print("Get present/absent")
muts_per_read <- mpileup_results_long %>%
  group_by(read_ids) %>%
  do(summarise(., mutids_read = get_mut_profile(.)))

#Merge identical reads
muts_per_read <- muts_per_read %>% group_by(mutids_read) %>%
  summarise(coverage = n())

#Remove reads with no mutations, as those are not informative
all_zero_reads <- strsplit(muts_per_read$mutids_read, ";")
all_zero_reads <- unlist(lapply(all_zero_reads, function(x) sum(grepl("[ATCG]:0", x)) == length(x)))

muts_per_read <- muts_per_read[!all_zero_reads, ]

#Merge mutation groups contained in each other-----------------------------------------------------------
print("Merge muts")
merged_mut_groups <- do.call("rbind", lapply(muts_per_read$mutids_read, function(x) replace_missing_mut(x, muts_per_read)))
  
#Add colnames
colnames(merged_mut_groups) <- c("group_id", "merged_groups", "coverage")
merged_mut_groups <- data.frame(merged_mut_groups, stringsAsFactors = FALSE)

#Keep only mutation groups that have the most groups merged (i.e. they are never merged themselves)
merged_groups_id <- unique(merged_mut_groups$merged_groups)
merged_groups_id <- unique(unlist(strsplit(merged_groups_id, "-")))

#Remove groups that have been merged
merged_mut_groups_filt <- merged_mut_groups[!merged_mut_groups$group_id %in% merged_groups_id, ]

merged_mut_groups_filt$coverage <- as.numeric(merged_mut_groups_filt$coverage)
#Output
out_file <- paste0("outputs/phasing/20230627/", T_sample, "_merged_mut_groups.csv")
write.csv(merged_mut_groups_filt, out_file,
          row.names = FALSE)

out_file_raw <- paste0("outputs/phasing/20230627/", T_sample, "_merged_mut_raw.csv")
write.csv(mpileup_results_raw, out_file_raw,
          row.names = FALSE)

out_file_all_raw <- paste0("outputs/phasing/20230627/", T_sample, '_', T_sample, "_merged_raw.csv")
write.csv(mpileup_results, out_file_all_raw, row.names = FALSE)

# Now merge together the HAPCUT phasing to get parental alleles phased to mutations ------------------------------

phasing_vcf_path <- '/nemo/project/proj-tracerx-lung/tctProjects/frankella/ctDNA_methods_dev/outputs/mut_cn_tracking/20230418/LTX0064_SU_N_snps_muts_phased.haplotype.phased.VCF'
phasing_vcf <- read.vcfR( phasing_vcf_path )

# Just extract info we need
phased <- data.table(chr =      getCHROM(phasing_vcf),
                     pos =      getPOS(phasing_vcf),
                     ref =      getREF(phasing_vcf),
                     alt =      getALT(phasing_vcf),
                     block_id = as.character(extract.gt(phasing_vcf, element='PS')),
                     gt =       as.character(extract.gt(phasing_vcf, element='GT')), 
                     stringsAsFactors = FALSE)

phased[, block_size_N_SNPs := .N, by = block_id ]
phased[, block_size_Mb := signif(max(pos)-min(pos)/1000000, 2), by = block_id ]

## use the mpileup to manually phase the mutations to nearby SNPs
mpileup_results <- as.data.table(mpileup_results)

# associate the hetSNPs with a nearby mutation
mpileup_results[ type == 'mut', mutid := paste(chr, pos, ref, sep = ':') ]

# there is probably a faster/cleaner way to do this... 

collate_nearest_snps <- function(mut_id, mpileup, max_phasing_dist, n_snps){

  # Collate require info on mut
  mutinfo <- mpileup_results[ mutid == mut_id, .(chr, pos, ref, alt, counts, mut_read, wt_read, mutid) ]
  setnames(mutinfo, paste0( names(mutinfo), '_mut'))
  setnames(mutinfo, 'mutid_mut', 'mutid')


  # Collate require info on nearby snps
  # only try to phase if within XKb (20 at the mo) 
  nearbysnps <- mpileup_results[ type == 'het_snp' &
                                 chr == mutinfo$chr_mut & 
                                 pos > mutinfo$pos - max_phasing_dist &
                                 pos < mutinfo$pos + max_phasing_dist, 
                                  .( chr, pos, ref, counts, mut_read, wt_read,
                                     dist = abs(pos - mutinfo$pos), id = mutid ) ]
  setnames(nearbysnps, paste0( names(nearbysnps), '_hetsnp'))

  
  if( nrow(nearbysnps) == 0 ) return( 'No hetSNPs within 20Kb' )

  # Take nearest 10 hetSNPs to mutation
  # supresswarnings when less than 10 hetSNPs within 20Kb
  nearbysnps <- nearbysnps[ order(dist_hetsnp) ][ 1:(ifelse(.N<n_snps, .N, n_snps)) ]
  setnames(nearbysnps, 'dist_hetsnp', 'dist_mut_to_snp')

  return( cbind(mutinfo, nearbysnps) )

}

muts <- mpileup_results[ type == 'mut', unique(mutid) ]
mut_phasing <- lapply(muts, function(mutid){
  collate_nearest_snps(mutid, mpileup_results, max_phasing_dist = max_read_length, n_snps = 10)
}) 

is_phasable <- sapply(mut_phasing, function(x) !all(x == 'No hetSNPs within 20Kb'))

unphased_muts <- muts[ !is_phasable ]
#write this out to record this - DO THIS LATER WHERE NOT ENOUGH OVERLAPPING READS
#fwrite(as.data.table(unphased_muts), 'outputs/phasing/20230627/unphased_muts.tsv', sep='\t')

phased_muts <- rbindlist( mut_phasing[ is_phasable ] )

### Now assess whether any read IDs are shared between the mutation and each SNP
phased_muts[, `:=`(SNV_MT_SNP_MT_N = sum(strsplit(mut_read_mut, split = ';')[[1]] %in% strsplit(mut_read_hetsnp, split = ';')[[1]]),
                   SNV_WT_SNP_MT_N = sum(strsplit(wt_read_mut, split = ';')[[1]]  %in% strsplit(mut_read_hetsnp, split = ';')[[1]]),
                   SNV_MT_SNP_WT_N = sum(strsplit(mut_read_mut, split = ';')[[1]] %in% strsplit(wt_read_hetsnp,  split = ';')[[1]]),
                   SNV_WT_SNP_WT_N = sum(strsplit(wt_read_mut, split = ';')[[1]]  %in% strsplit(wt_read_hetsnp,  split = ';')[[1]])),
            , 1:nrow(phased_muts) ]

# Get some stats
phased_muts[, `:=`(shared_depth   = sum( SNV_MT_SNP_MT_N, SNV_WT_SNP_MT_N, SNV_MT_SNP_WT_N, SNV_WT_SNP_WT_N ),
                   mut_varcount = sum(SNV_MT_SNP_MT_N, SNV_MT_SNP_WT_N, na.rm=T), 
                   snp_vaf_on_mut = SNV_MT_SNP_MT_N / sum(SNV_MT_SNP_MT_N, SNV_MT_SNP_WT_N, na.rm=T),
                   snp_vaf_on_wt  = SNV_WT_SNP_MT_N / sum(SNV_WT_SNP_MT_N, SNV_WT_SNP_WT_N, na.rm=T)),
             , 1:nrow(phased_muts) ]

## save this while testing - takes 10 minutes to run the above
fwrite(phased_muts, 'outputs/phasing/20230627/tmp_phased_muts.tsv', sep='\t')
#phased_muts <- fread( 'outputs/phasing/20230627/tmp_phased_muts.tsv' )

# overlay block and genotype for each snp from HAPCUT2
phased[, id_hetsnp := paste(chr, pos, ref, sep = ':') ]
phased_muts <- as.data.table( inner_join( phased_muts, 
                                          phased[, .(id_hetsnp, gt, block_id )] ) )

# get counts from each haplotype per SNP
phased_muts[, gt := gsub('\\/', '|', gt) ]
phased_muts[, Mut_hap1_varcount := ifelse(gt == '0|1', SNV_MT_SNP_MT_N, ifelse(gt == '1|0', SNV_MT_SNP_WT_N, NA)) ]
phased_muts[, Mut_hap2_varcount := ifelse(gt == '0|1', SNV_MT_SNP_WT_N, ifelse(gt == '1|0', SNV_MT_SNP_MT_N, NA)) ]

# Get Haplotype VAF per mutation
# assign each mutation a block
phased_muts[, mut_block := block_id[ dist_mut_to_snp == min(dist_mut_to_snp) ][1], by = mutid ]
#remove SNPs which aren't in the same block as the mutation
phased_muts <- phased_muts[ block_id == mut_block ]
# summarise per mutation
muts_only_phased <- unique( phased_muts[mut_varcount > 0, 
                                         .(chr = chr_mut, pos = pos_mut, ref = ref_mut, #alt = alt_mut,
                                           block_id, hapvaf = sum(Mut_hap1_varcount)/(sum(Mut_hap1_varcount) + sum(Mut_hap2_varcount)), 
                                           hapdepth = sum(Mut_hap1_varcount) + sum(Mut_hap2_varcount), 
                                           hap1_id = paste(ifelse(gt == '0|1', 'A', 'R'), collapse = ':'), 
                                           hap_id_SNP_pos = paste(pos_hetsnp, collapse = ':')), 
                                        by = mutid ] )

# categorise phases
muts_only_phased[, phase := ifelse( hapvaf > 0.75 & hapdepth > 2, 1, 
                                 ifelse(hapvaf < 0.25 & hapdepth > 2, 0, NA) ) ]

fwrite(muts_only_phased, 'outputs/phasing/20230627/phased_mutations.tsv', sep = '\t')


# summarise blocks
muts_only_phased[, `:=`(num_muts_block = sum(!is.na(phase)), 
                        sum_phase1_block = sum(phase == 1,na.rm=T),
                        sum_phase0_block = sum(phase == 0,na.rm=T)),
                 by = block_id ]

muts_only_phased[ hapdepth>5, sum(hapvaf>0.9)/.N ]
muts_only_phased[ hapdepth>5, sum(hapvaf<0.1)/.N ]

muts_only_phased[ hapdepth==2, table(hapvaf) ]
muts_only_phased[ hapdepth==3, table(hapvaf) ]

muts_only_phased[, table(table(block_id))]

muts_only_phased[, block_muts_N := .N, by = block_id ]

pdf( 'outputs/phasing/20230627/plot.pdf' )
muts_only_phased[, hist(hapvaf)]
dev.off()

# # get overall score for both alleles based on SNP var on mutant alleles (WT alleles - aren't super informative unsurprisingly w purity CN etc)
# if its abigious then limit to closer variants







phased_muts[, id_hetsnp := paste(chr_hetsnp, pos_hetsnp, ref_hetsnp)]

plot_data <- data.table::melt(phased_muts, id.vars = c('mutid','id_hetsnp', 'shared_depth', 'mut_varcount'),
                    measure.vars = c('snp_vaf_on_mut', 'snp_vaf_on_wt'), 
                    value.name = 'vaf',
                  variable.name = 'type')

plot_data[, id := paste(mutid, id_hetsnp)]
ids <- plot_data[, unique(id)]

pdf( 'outputs/phasing/20230627/plot.pdf' )
ggplot(plot_data[ mut_varcount > 1 & id %in% ids[10000:15000] ], 
          aes(x = type, fill = type, y = vaf)) +
    geom_point() +
    geom_violin() +
    geom_line( aes(group =  id ), alpha = 0.2 ) +
    labs( x = '',
          y = 'vaf of snp on mut loci overlapping reads') +
    theme_classic() +
    theme( text = element_text(size=20))
dev.off()
pdf( 'outputs/phasing/20230627/plot1.pdf' )
phased_muts[ shared_depth > 3, hist(snp_vaf_on_mut, snp_vaf_on_wt)]
dev.off()

# Make a longer format where each mutation has 1 row per nearby SNP (any SNP within 5Kb) 
# or nearest 4 SNPs (can't be more than 50Kb away or mutation unphasable) - need at least 1 to get phase
# record distance to mutation, fraction/number of mutant reads that phase to VAR/REF - 3 cols. 
# mutation pos depth, mutation vaf, snp depth & vaf, a call whether mut phases to ref/var/uncertain
# check that this SNP is also a het snp in the phasing - if not don't use it

# Summarise per mutation and get a string indicating whether the mutation phases to the REF or VAR
# for each SNP also keep a string of positions for the SNPs used in the ID, how many SNPs used, how 
# close they were on average to mutation

# Look at mutation position in phasing and extract het SNP ref/vars from above for each allele present

# Create a score for how closely both alleles match and pick the best score. If score is poor for both
# then try rescoring with a smaller number of SNPs close to the mutation. Keep the score, the ID from the
# milepileup & ID from hapcut, the Phased block from hapcut and the 0/1 allele status (consistent accross 
# block)

# output this & long version with all SNPs as rows


tmp <- mpileup_results[type == 'mut'][1000]
tmp[, .(chr, pos)]
mpileup_results[type == 'het_snp' & chr == tmp$chr & pos > tmp$pos-1000 &  pos < tmp$pos+1000, .(chr, pos)]
phased_muts[ pos_hetsnp == 113248716 ]





