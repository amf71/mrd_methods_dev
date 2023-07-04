library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(intervals)
library(gtools)
library(GenomicRanges)
library(fst)
library(tidyr)
library(RColorBrewer)
library(GenomeInfoDb)
library(vcfR)
#Script to analyse the output from the mutation merging using Mpilueup, including overlaying phased mutations with long read data with tree building information

#Functions-----------------------
#Function to load the merged mutation info, including sample name
get_merged_mut <- function(path){
  merged_muts_df <- read.csv(path, stringsAsFactors = FALSE)
  merged_muts_df$sample_name <- gsub(pattern = "(.+/)(.+)(_merged.+)",
                                     replacement = "\\2",
                                     x = path)
  return(merged_muts_df)
}

#Function to load the variant positions, including sample name
get_vars_df <- function(path){
  vars_df <- read.table(path, stringsAsFactors = FALSE,
                        sep = "\t", header = FALSE)
  colnames(vars_df) <- c("chr", "start")
  vars_df$sample_name <- gsub(pattern = "(.+/)(.+)(\\.snvs.bed)",
                              replacement = "\\2",
                              x = path)
  return(vars_df)
}

#Function to count the number of mutations detected by mpileup
get_mpileup_muts <- function(path){
  raw_mpileup <- read.csv(path,
                          stringsAsFactors = FALSE)

  #Count number of mutations detected by mpileup
  nb_muts_detected <- length(unique(raw_mpileup$mutid[raw_mpileup$mut_present]))
  
  #Get sample name
  sample_name <- gsub(pattern = "(results/)(.+)(_merged_.+)",
                      replacement = "\\2",
                      x = path)
  return(c(sample_name, nb_muts_detected))
}

#Long read BAM QC------------------------------------------------------------------------------
bamqc_lr <- read.table("input/all_bamqc_stats_coverage.txt",
                       sep = "\t",
                       header = TRUE)

#Load read length per sample
cohort_sample_full <- read.csv("input/cohort_summary_length.csv",
                               stringsAsFactors = FALSE)

#Keep only median, ensure compatible sample names
median_length <- cohort_sample_full[, c("fif_perc", "patient", "sample_name")]
median_length$sample_name <- gsub(pattern = "DNA[-_]",
                                  replacement = "",
                                  x = median_length$sample_name)

median_length$sample_name <- gsub(pattern = "(SU)(-)(T[0-9])_",
                                  replacement = "\\1_\\3-",
                                  x = median_length$sample_name)

median_length$sample_name <- paste0(median_length$patient, "_", median_length$sample_name)

#Get all merged mut files------------------------------------------------
merged_mut_files <- list.files("results", pattern = "merged_mut_groups.csv", full.names = TRUE)
merged_mut_files <- merged_mut_files[!merged_mut_files %in% c("results/20230620/LTX0388_SU_T1-R3_merged_mut_groups.csv",
                                                             "results/test/LTX0038_SU_T1-R3_merged_mut_groups.csv")]

all_merged_muts <- do.call("rbind", lapply(merged_mut_files, get_merged_mut))

#Plot distributions
ggplot(all_merged_muts, aes(x = coverage)) + geom_histogram(binwidth = 1) +
  facet_grid(sample_name ~ .) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12))
  
#Get number of mutations per group as well as distance covered in the genome
all_merged_muts$nb_muts <- str_count(all_merged_muts$group_id, ";") +1

#Get first and last chromosome of the group, flag groups with different start-ends
all_merged_muts$chr <- gsub(pattern = ":.+",
                            replacement = "",
                            x = all_merged_muts$group_id)
last_chr <- gsub(pattern = "(.+;)(chr.+)(:[0-9]+:[ATCG]:[ATCG]:[01]$)",
                 replacement = "\\2",
                 x = all_merged_muts$group_id)

all_merged_muts$diff_chrom <- all_merged_muts$chr != last_chr

#Get distances covered (only when on the same chromosome)
all_merged_muts$start <- as.numeric(gsub(pattern = "(^chr..{0,1}:)([0-9]+)(:[ATCG]:.+)",
                                    replacement = "\\2",
                                    x = all_merged_muts$group_id) )

all_merged_muts$end <- as.numeric(gsub(pattern = "(.+;chr..{0,1}:)([0-9]+)(:[ATCG]:[ATCG]:[01]$)",
                                       replacement = "\\2",
                                       x = all_merged_muts$group_id))

all_merged_muts$dist <- ifelse(all_merged_muts$diff_chrom,
                               NA, all_merged_muts$end - all_merged_muts$start)

#Plot number of mutations and distance by coverage
ggplot(all_merged_muts, aes(x = coverage, y = nb_muts)) +
  geom_point() + geom_smooth(method = "lm") +
  stat_cor() +
  facet_grid(sample_name ~ .) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12))

ggplot(all_merged_muts, aes(x = coverage, y = log10(dist))) +
  geom_point() + geom_smooth(method = "lm") +
  stat_cor() +
  ylab("log10(Size)") +
  facet_grid(sample_name ~ .) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12))

ggplot(all_merged_muts, aes(x = dist)) +
  facet_grid(sample_name ~ .) +
  geom_histogram(fill = "#1C9099") +
  scale_x_log10() +
  xlab("log10(Size -bp-)") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12))

ggplot(all_merged_muts, aes(x = nb_muts)) + geom_histogram() +
  facet_grid(sample_name ~ .) +
  xlab("Mutations per group") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12))

#Get number of variants per chromosome, compare number of groups per chromosome------------------------------------------------------
all_chroms <- unique(all_merged_muts$chr)

all_bed <- list.files("tmp/", pattern = ".+.snvs.bed", full.names = TRUE)

#Keep only bed files from samples present in the grouped mutations file
samples_merged_muts <- unique(all_merged_muts$sample_name)
samples_merged_muts_pattern <- paste(samples_merged_muts,
                                     collapse = "|")

all_bed <- grep(pattern = samples_merged_muts_pattern, x = all_bed, value = TRUE)

all_variants <-  do.call("rbind", lapply(all_bed, get_vars_df))

#Get the number of variants detected by mpileup only
all_mpileup_raw <- list.files("results", full.names = TRUE,
                              pattern = "merged_mut_raw.csv")

nb_vars_mpileup <-  do.call("rbind", lapply(all_mpileup_raw, get_mpileup_muts))

#Get the number of variants per chromosome and per sample
nb_vars_chr <- all_variants %>% group_by(chr, sample_name) %>% tally()

#Add number of groups and median length of the groups per sample and per chromosome
nb_mut_groups <- all_merged_muts[!all_merged_muts$diff_chrom, ] %>% group_by(chr, sample_name) %>% 
  summarise(nb_groups = n_distinct(group_id),
            median_size = median(dist))

#Merge and plot
nb_mut_groups <- inner_join(nb_mut_groups, nb_vars_chr)
nb_mut_groups$chr <- factor(nb_mut_groups$chr,
                            levels = mixedsort(unique(nb_mut_groups$chr)))

ggplot(nb_mut_groups, aes(x = chr, y = median_size)) + geom_bar(stat = "identity") +
  facet_grid(sample_name ~ .) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 16),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12),
        legend.position = "none")

ggplot(nb_mut_groups, aes(x = chr, y = nb_groups / n)) + geom_bar(stat = "identity") +
  facet_grid(sample_name ~ .) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 16),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12),
        legend.position = "none")

ggplot(nb_mut_groups, aes(x = chr, y = nb_groups)) + geom_bar(stat = "identity") +
  facet_grid(sample_name ~ .) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 16),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12),
        legend.position = "none")

#Plot the proportion of mutations grouped per chromosome---------------------------------------------------------
nb_grouped_muts <- all_merged_muts[!all_merged_muts$diff_chrom, ]
nb_grouped_muts$mut_ids <- gsub(pattern = "(:[01];)|(:[01]$)",
                                replacement = "",
                                x = nb_grouped_muts$group_id)

nb_grouped_muts$mut_ids <- gsub(pattern = "([ATCG])(chr)",
                                replacement = "\\1;\\2",
                                x = nb_grouped_muts$mut_ids)

#Get unique list of grouped mutations
list_muts_id <- strsplit(nb_grouped_muts$mut_ids, ";")
names(list_muts_id) <- nb_grouped_muts$sample_name

nb_grouped_muts <- data.frame(mut_id = unlist(list_muts_id),
                              sample_name = rep(names(list_muts_id), sapply(list_muts_id, length)))

#Ensure only unique mutations
nb_grouped_muts <- unique(nb_grouped_muts)

nb_grouped_muts$chr <- gsub(pattern = ":.+",
                            replacement = "",
                            nb_grouped_muts$mut_id)

nb_grouped_muts <- nb_grouped_muts %>% group_by(sample_name, chr) %>% summarise(nb_grouped_muts = n_distinct(mut_id))

#Add number of mutations
nb_grouped_muts <- inner_join(nb_grouped_muts, nb_vars_chr)
nb_grouped_muts$ungrouped_muts <- nb_grouped_muts$n - nb_grouped_muts$nb_grouped_muts
nb_grouped_muts$prop_grouped <- nb_grouped_muts$nb_grouped_muts / nb_grouped_muts$n

#Plot proportion
nb_grouped_muts$chr <- factor(nb_grouped_muts$chr,
                            levels = mixedsort(unique(nb_grouped_muts$chr)))

ggplot(nb_grouped_muts,
       aes(x = chr, y = prop_grouped)) + geom_bar(stat = "identity") +
  facet_grid(sample_name ~ .) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 16),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12))

#Turn into long format
nb_grouped_muts_long <- nb_grouped_muts
nb_grouped_muts_long$n <- NULL
nb_grouped_muts_long$prop_grouped <- NULL

nb_grouped_muts_long <- reshape2::melt(nb_grouped_muts_long, c("sample_name", "chr"),
                                       variable.name = "grouping", value.name = "nb_muts")

nb_grouped_muts_long$chr <- factor(nb_grouped_muts_long$chr,
                              levels = mixedsort(unique(nb_grouped_muts_long$chr)))

ggplot(nb_grouped_muts_long, aes(x = chr, y = nb_muts, fill = grouping)) + geom_bar(stat = "identity") +
  facet_grid(sample_name ~ .) +
  scale_fill_manual(values = brewer.pal(2, "Set2"),
                    labels = c("In a group", "Ungrouped")) +
  theme_bw() +
  ylab("Number of mutations") +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_blank())

#Get the total proportion grouped per sample, total TMB per Mb as well as coverage and average read length------------------------------------------------------
hg38_genome <- Seqinfo(genome="hg38")
total_size <- sum(seqlengths(hg38_genome[paste0("chr", c(1:22, "X", "Y"))]))
total_size_mb <- total_size / 10^6

total_muts_sample <- nb_grouped_muts %>% group_by(sample_name) %>%
  dplyr::summarise(total_muts = sum(n),
                   total_grouped = sum(nb_grouped_muts))

#Get mpileup mutation counts
colnames(nb_vars_mpileup) <- c("sample_name", "total_muts_mpileup")
nb_vars_mpileup <- data.frame(nb_vars_mpileup)
nb_vars_mpileup$total_muts_mpileup <- as.numeric(as.character(nb_vars_mpileup$total_muts_mpileup))
nb_vars_mpileup$sample_name <- as.character(nb_vars_mpileup$sample_name)

#Use total number of mutations from mpileup when available
total_muts_sample <- left_join(total_muts_sample, nb_vars_mpileup)

total_muts_sample$prop_grouped <- ifelse(!is.na(total_muts_sample$total_muts_mpileup),
                                         total_muts_sample$total_grouped / total_muts_sample$total_muts_mpileup,
                                         total_muts_sample$total_grouped / total_muts_sample$total_muts)

colnames(bamqc_lr)[colnames(bamqc_lr) == "Sample"] <- "sample_name"
total_muts_sample_qc <- inner_join(total_muts_sample, bamqc_lr)
total_muts_sample_qc$muts_mb <- ifelse(!is.na(total_muts_sample$total_muts_mpileup),
                                       total_muts_sample$total_muts_mpileup / total_size_mb,
                                       total_muts_sample$total_muts / total_size_mb)

total_muts_sample_qc <- inner_join(total_muts_sample_qc, median_length)

ggplot(total_muts_sample_qc, aes(x = average_length, y = muts_mb,
                                color = prop_grouped, size = coverage)) +
  
  ylim(c(5, NA)) + xlim(c(2500, NA)) +
  geom_point() +
  scale_colour_gradient(low = "white", high = "darkblue",
                        limits = c(0, NA),
                        name = "Proportion\ngrouped") +
  xlab("Average read length") + ylab("Mutations (per Mb)") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

ggplot(total_muts_sample_qc, aes(x = prop_grouped, y = muts_mb,
                                 color = average_length, size = coverage)) +
  
  ylim(c(5, NA)) +
  geom_point() +
  scale_colour_gradient(low = "white", high = "darkblue",
                        limits = c(0, NA),
                        name = "Read length") +
  xlab("Proportion grouped") + ylab("Mutations (per Mb)") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

#Clean table and output
out_prop_table <- total_muts_sample_qc[, c("sample_name", "muts_mb", "average_length", "coverage", "prop_grouped")]

write.csv(out_prop_table, row.names = FALSE,
          file = "~/Downloads/out_prop_table.csv")
