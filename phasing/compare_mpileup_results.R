library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(intervals)
library(gtools)
library(GenomicRanges)
library(liftOver)
library(fst)
library(tidyr)
library(RColorBrewer)
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

extract_ccf <- function(pyclone_full, sample){
  
  sample_id <- gsub(pattern = ".+:",
                    replacement = "",
                    x = sample)
  
  #Divide the ccf vector
  split_ccf <- strsplit(pyclone_full, ";")
  ccf_sample <- unlist(lapply(split_ccf, function(x) grep(pattern = sample_id,
                                                          x = x,
                                                          value = TRUE)))
  ccf_sample <- as.numeric(gsub(pattern = ".+:",
                                replacement = "",
                                x = ccf_sample))
  
}


#Get all merged mut files------------------------------------------------
merged_mut_files <- list.files("results", pattern = "merged_mut_groups.csv", full.names = TRUE)

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

ggplot(all_merged_muts[all_merged_muts$sample_name == "LTX0388_SU_T1-R1", ], aes(x = coverage, y = log10(dist))) +
  geom_point() + geom_smooth(method = "lm") +
  stat_cor() +
  ylab("log10(Size)") +
  facet_grid(sample_name ~ .) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12))

ggplot(all_merged_muts[all_merged_muts$sample_name == "LTX0388_SU_T1-R1", ], aes(x = dist)) +
  geom_histogram(fill = "#1C9099") +
  scale_x_log10() +
  xlab("log10(Size -bp-)") +
  theme_bw() +
  theme(axis.text = element_text(size = 24),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 24))

ggplot(all_merged_muts[all_merged_muts$sample_name == "LTX0388_SU_T1-R1", ], aes(x = nb_muts)) + geom_histogram() +
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

ggplot(nb_grouped_muts[nb_grouped_muts$sample_name == "LTX0388_SU_T1-R1", ],
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

ggplot(nb_grouped_muts_long[nb_grouped_muts_long$sample_name == "LTX0388_SU_T1-R1", ], aes(x = chr, y = nb_muts, fill = grouping)) + geom_bar(stat = "identity") +
  # facet_grid(sample_name ~ .) +
  scale_fill_manual(values = brewer.pal(2, "Set2"),
                    labels = c("In a group", "Ungrouped")) +
  theme_bw() +
  ylab("Number of mutations") +
  theme(axis.text.y = element_text(size = 24),
        axis.text.x = element_text(size = 24, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 28),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 24),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 28),
        legend.title = element_blank())

#Compare to WES trees----------------------------------------------------------------------------------------------
#Lift over variants to hg19
#Turn into a granges object
all_variants_uniq <- all_variants
all_variants_uniq$sample_name <- NULL
all_variants_uniq <- unique(all_variants_uniq)

all_variants_gr <- GRanges(Rle(all_variants_uniq$chr), IRanges(start = all_variants_uniq$start, width = 1),
                            mutid_hg38 = paste(all_variants_uniq$chr, all_variants_uniq$start, sep = "_"))

#Load the chain file
chain_file <- import.chain("input/hg38ToHg19.over.chain")

#Lift over
seqlevelsStyle(all_variants_gr) = "UCSC"
h19_vars_gr <- liftOver(all_variants_gr, chain_file)
h19_vars_df <- as.data.frame(unlist(h19_vars_gr))
h19_vars_df <- h19_vars_df[, c("seqnames", "start", "mutid_hg38")]
colnames(h19_vars_df) <- c("chr_hg19", "pos_hg19", "mutid_hg38") 

#Check whether mutations within the same haplo block tend to also be in the same clone/branch of the tree
#Load muttable
muttable_df <- read_fst("input/tx421_RegionMutTable.fst")
muttable_df <- muttable_df[muttable_df$PASS, ]
muttable_df <- muttable_df[as.logical(muttable_df$isSNV), ]
muttable_df <- muttable_df[muttable_df <- muttable_df$Is.present.region, ]

#Add hg19 positions to compare with Tx421 WES mutations
mut_groups_id <- strsplit(all_merged_muts$group_id, split = ";")
mpileup_groups_long <- data.frame(mutation_id = unlist(mut_groups_id),
                                  group_id = rep(all_merged_muts$group_id, sapply(mut_groups_id, length)),
                                  sample_name = rep(all_merged_muts$sample_name, sapply(mut_groups_id, length)),
                                  nb_muts = rep(all_merged_muts$nb_muts, sapply(mut_groups_id, length)),
                                  stringsAsFactors = FALSE)

mpileup_groups_long$mutid_hg38 <- gsub(pattern = "(.+)(:)(.+)(:)(.+)(:)(.+)(:)(.+)",
                                       replacement = "\\1_\\3",
                                       x = mpileup_groups_long$mutation_id)

mpileup_groups <-  inner_join(mpileup_groups_long, h19_vars_df)

mpileup_groups$mutid_hg19 <- paste(mpileup_groups$chr_hg19,
                                   mpileup_groups$pos_hg19,
                                   sep = ":")

#Get the columns of interest for the relevant samples
muttable_short <- muttable_df[, c("mutation_id", "RegionID", "PyCloneCluster_SC", "PyCloneCCF_SC", "mut_cpn_SC")]

#Clean CCF to keep only CCF of interest for the region
muttable_short$ccf_clean <- apply(muttable_short, 1,
                                  function(x) extract_ccf(x["PyCloneCCF_SC"], x["RegionID"]))

#Keep only samples of interest
muttable_short$sample_name <- gsub(pattern = "(LTX)([0-9]+)(:)(.+)(\\.)(.+)",
                                   replacement = "\\10\\2_\\4-\\6",
                                   x = muttable_short$RegionID)

muttable_short <- muttable_short[muttable_short$sample_name %in% samples_merged_muts, ]

muttable_short <- muttable_short[!is.na(muttable_short$PyCloneCluster_SC), ]
muttable_short$ccf_clean <- unlist(muttable_short$ccf_clean)

muttable_short$mutid_hg19 <- gsub(pattern = "(LTX[0-9]+:)(.+:[0-9]+)(:.+)",
                                  replacement = "chr\\2",
                                  x = muttable_short$mutation_id)

muttable_short$mutation_id <- NULL

mpileup_groups_clusters <- inner_join(muttable_short, mpileup_groups)

#Plot the number of variants per haplo block shared across multiple clusters
mpileup_groups_clusters <- mpileup_groups_clusters %>% group_by(group_id, sample_name) %>% mutate(same_cluster = all(PyCloneCluster_SC[1] == PyCloneCluster_SC)) %>%
  group_by(group_id, PyCloneCCF_SC, sample_name) %>% mutate(nb_vars_cluster = n())

mpileup_groups_clusters$PyCloneCluster_SC <- as.factor(mpileup_groups_clusters$PyCloneCluster_SC)
#Plot
ggplot(mpileup_groups_clusters[mpileup_groups_clusters$sample_name == "LTX0388_SU_T1-R1", ],
       aes(x = PyCloneCluster_SC, y = ccf_clean)) +
  #geom_boxplot(outlier.colour = NA) +
  geom_point(aes(size = nb_muts), alpha = 0.6, position = position_jitter(width = 0.2, seed = 123)) + 
  geom_path(aes(group = group_id),  position = position_jitter(width = 0.2, seed = 123)) +
  #facet_grid(. ~ sample_name, scales = "free") +
  scale_x_discrete(labels = c("1", "2", "3", "4")) +
  scale_size_continuous(range = c(4, 8),
                        name = "# mutations\ngrouped") +
  theme_bw() +
  ylab("Cancer cell fraction") + xlab("Clone ID") +
  theme(axis.text = element_text(size = 24),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 24),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 28))

ggplot(mpileup_groups_clusters, aes(x = PyCloneCluster_SC, y = mut_cpn_SC)) + geom_boxplot(outlier.colour = NA) +
  geom_point(aes(color = same_cluster, size = nb_muts), position = position_jitter(width = 0.2, seed = 123)) + 
  geom_path(aes(group = group_id, color = same_cluster),  position = position_jitter(width = 0.2, seed = 123)) +
  facet_grid(. ~ sample_name, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))
