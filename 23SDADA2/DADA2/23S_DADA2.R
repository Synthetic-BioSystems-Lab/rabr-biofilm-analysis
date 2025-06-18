library(dada2) 
library(tidyverse) 
library(tidyr) 
library(ggplot2) 
library(ggtext) 

setwd("~/Miller Lab/Rscripts_PilotRABR/DADA2")

#Functions# 
# Create Directories Function # 
create_directory <- function(dir_path) { 
  #This code will check if your directory already  
  #exists, and will create that directory if it 
  #can't find it.) 
  #check if directory exists 
  
  if(!dir.exists(dir_path)) { 
    dir.create(dir_path, recursive = TRUE) 
    message("Directory created at ", dir_path) 
  } else { 
    message("Directory already exists at ", dir_path) 
  } 
} 

# Set up paths and directories # 
create_directory("./23SFastq") 
create_directory("./23SFastq/23SF") 
create_directory("./23SFastq/23SR") 
create_directory("./23SFastq/23SF/filtered") 
create_directory("./23SFastq/23SR/filtered") 
path <- "./23SFastq" 
forward_data_path <- "./23SFastq/23SF" 
reverse_data_path <- "./23SFastq/23SR" 
forward_data_path_filtered_23S <- file.path(forward_data_path, "filtered") #filtered forward files go into pathF/filtered subdirectory 
reverse_data_path_filtered_23S <- file.path(reverse_data_path, "filtered") #filtered rev files -> pathR filtered subdirectory 

# Checking quality of data # 
forward_read_quality_23S <- sort(list.files(path, pattern="_R1_001.fastq.gz",  
                                             full.names = TRUE)) 
reverse_read_quality_23S <- sort(list.files(path, pattern="_R2_001.fastq.gz",  
                                             full.names = TRUE)) 
message("Here is a graph of the quality of the forward reads.") 
plotQualityProfile(forward_read_quality_23S[1:2]) 
ggsave("23SFastq/23SDADA2_quality_forward.tiff", 
       width=9, height=4) 
message("Here is a graph of the quality of the reverse reads.") 
plotQualityProfile(reverse_read_quality_23S[1:2]) 
ggsave("23SFastq/23SDADA2_quality_reverse.tiff", 
       width=9, height=4) 

# Filtering Forward and Reverse Files # 
#filters files based on R1 and R2. R1 means forward, R2 means reverse 
#dir() selects files, file.copy moves them 
#This code is based on: https://stackoverflow.com/questions/49593823/r-file-copy-function 
forward_data_files <- dir(path, pattern = "*R1_001.fastq.gz",  
                          all.files = FALSE) 
file.copy(file.path(path, forward_data_files), forward_data_path, overwrite=TRUE) 

#Reverse Reads 
reverse_data_path <- "./23SFastq/23SR" 
data_filesR <- dir(path, pattern = "*R2_001.fastq.gz",  
                   all.files = FALSE) 
file.copy(file.path(path, data_filesR), reverse_data_path, overwrite=TRUE) 
#removing intermediate files 

# File Dissecting for Filtering # 
forward_data_path_filtered_23S <- file.path(forward_data_path, "filtered")  
reverse_data_path_filtered_23S <- file.path(reverse_data_path, "filtered")  
forward_fastq_list_23S <- sort(list.files(forward_data_path, pattern = "fastq")) 
reverse_fastq_list_23S <- sort(list.files(reverse_data_path, pattern = "fastq")) 
if(length(forward_fastq_list_23S) != length(reverse_fastq_list_23S)) 
  stop("Forward and reverse reads do not match. :(") 

# Filtering # 
#Note: These Parameters aren't optimal for all data sets 
out <- filterAndTrim(fwd=file.path(forward_data_path, forward_fastq_list_23S),
                     filt=file.path(forward_data_path_filtered_23S,forward_fastq_list_23S),
                     trimLeft = 20, truncLen = 200, 
                     maxEE = 2, truncQ = 2,
                     maxN=0, rm.phix = TRUE, compress = TRUE, verbose = TRUE, 
                     multithread = FALSE) 
head(out) 

forward_data_path_filtered_23S <- "./23SFastq/23SF/filtered" #forward_data_path_filtered_23S has now been changed to the "filtered" directory 
forward_fastq_filesE <- list.files(forward_data_path_filtered_23S, pattern="fastq.gz", full.names = TRUE) 
sample_names_23S <- sapply(strsplit(basename(forward_fastq_filesE), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz 
# if(!identical(sample_names_23S, samplenamesRE))  
#   stop("Forward and reverse reads don't match") 

# Infer Sequence Variants # 
names(forward_fastq_filesE) <- sample_names_23S 
set.seed(100) 
# Learn forward & reverse error rates 
#This one takes a while 
forward_error_rate <- learnErrors(forward_fastq_filesE, nbases = 1e8, multithread = TRUE) 

# Sample Inference & Merger of Paired-End Reads # 
mergers <- vector("list", length(sample_names_23S)) 
names(mergers) <- sample_names_23S
for(sam in sample_names_23S) { 
  cat("Processing:", sam, "\n") 
  derepF <- derepFastq(forward_fastq_filesE[[sam]]) 
  ddF <- dada(derepF, err = forward_error_rate, multithread = TRUE) 
  #merger <- mergePairs(ddF, derepF) #argument "dadaR" is missing, with no default
  mergers[[sam]] <- ddF #merger 
} 
rm(derepF) 

# construct sequence table and remove chimeras # 
seqtabE <- makeSequenceTable(mergers) 
create_directory("./23SFastq/results") 
saveRDS(seqtabE, "./23SFastq/results/seqtabE.rds")



# Merge Runs, Remove Chimeras, Assign Taxonomy # 
pathRDS <- "./23SFastq/results/seqtabE.rds" # Where the RDS files were saved 
#Remove chimeras 
seqtabE <- removeBimeraDenovo(readRDS(pathRDS), method="consensus", multithread = TRUE) 
#Assign taxonomy 

#Note: Silva file comes from https://zenodo.org/records/1172783 
#silva_nr_v132_train_set.fa is the name of the file 


dim(seqtabE)
summary(colSums(seqtabE))
sample(getSequences(seqtabE), 10)

taxE <- assignTaxonomy(seqtabE, "./microgreen_algaebase_dada2.fa",  
                       multithread = TRUE)
# taxE <- assignTaxonomy(seqtabE, "./silva_nr_v132_train_set.fa.gz",  
#                        multithread = TRUE) DIDN'T WORK

# Write to disk 

saveRDS(seqtabE, "./23SFastq/results/seqtabE.rds") # This is where the sequencing table file is saved 
saveRDS(taxE, "./23SFastq/results/taxE.rds") # This is where the taxonomy file is saved 

# Visualization # 
#We used ggplot2 and stuff to make this graph, not phyloseq 
#Code based on 23S's code that he used to visualize his data. 
#Note that this code also calculates relative abundace in order to  
# visualize this data. 

#Rename first col of taxE > "ASVs"  
asv_taxE <- cbind(rownames(taxE), data.frame(taxE), row.names=NULL) %>% 
  rename("ASVs" = "rownames(taxE)") 

asv_seqtabE <- cbind(rownames(seqtabE), data.frame(seqtabE), row.names=NULL) %>% 
  rename("sample_id" = "rownames(seqtabE)") %>% 
  pivot_longer(!sample_id, names_to = "ASVs", values_to = "count") 
#merge otu and taxon 
composite <- inner_join(asv_seqtabE, asv_taxE, by="ASVs") 

#remove empty rows, add relative abundance, pivot for readability 
trimmed_composite <- composite[!(composite$count=="0"),] %>% 
  group_by(sample_id) %>% 
  #relative abundance calculated here 
  mutate(rel_abund = count / sum(count)) %>% 
  ungroup() %>% 
  pivot_longer(cols = c("Kingdom", "Phylum", "Class",  
                        "Order", "Family", "Genus"),  
               names_to = "level", values_to = "taxon") 

taxon_rel_abund <- trimmed_composite %>%
  filter(level=="Phylum") %>%
  group_by(sample_id, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups = "drop") %>%
  group_by(sample_id, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups = "drop") %>%
  mutate(taxon = str_replace(taxon, "(.*)_unclassified", "Unclassified *\\1*"), 
         taxon = str_replace(taxon, "^(\\S*)$", "*\\1*")) 

#set aside other 
taxon_pool <- taxon_rel_abund %>% 
  group_by(taxon) %>% 
  summarize(pool = max(mean_rel_abund) < 5, .groups="drop") 

#assemble others and make RA stacked bar plot 
taxon_assembled <- inner_join(taxon_rel_abund, taxon_pool, by="taxon") %>% 
  mutate(taxon = if_else(pool, "Other", taxon)) %>% 
  group_by(sample_id, taxon) %>% 
  summarize(mean_rel_abund = sum(mean_rel_abund), .groups = "drop") %>% 
  ggplot(aes(x = sample_id, y = mean_rel_abund, fill = taxon)) + 
  geom_col() + 
  #scale_fill_manual(name = NULL, values = c(brewer.pal(6, "Dark2"), "gray")) + 
  scale_fill_discrete(name=NULL) + 
  scale_y_continuous(expand=c(0,0)) + 
  labs(x = NULL, 
       y = "Mean Relative Abundance (%)") + 
  theme_classic() + 
  theme(legend.text = element_markdown(), 
        axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0), 
        legend.key.size = unit(10, "pt")) 

ggsave("23SFastq/results/23SDADA2_23S_stacked_bar.tiff", 
       width=9, height=4) 
