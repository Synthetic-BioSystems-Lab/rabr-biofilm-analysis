#R Analysis is done in RStudio, 4.3.2, “Eye Holes” 
## Amanda files, DADA2 analysis ## 
#https://benjjneb.github.io/dada2/bigdata_paired.html 
#don't do the basic DADA2 protocol cuz we have too much data 
#Note: assume that primers have been removed 
#Note: assume we have DADA2 up to date (even if I have doubts) 
#Note: The filename parsing in the workflow I'm using 
#expects the paired forward and reverse fastq files to have the 
#same samplename_ prefix. Might cause problems. 
#Note: The workflow expects demultiplexed, per-sample, gzipped fastq files  
#for the primer-free forward reads of a Hiseq run to be in the directory pathF. 
#I made a homemade algorithm to accomplish this. 

# libraries # 
library(dada2)
library(tidyverse) 
library(tidyr) 
library(ggplot2) 
library(ggtext) 

setwd("~/Miller Lab/Rscripts_PilotRABR/DADA2")

# Objects that I renamed # 
#filtFs/fnFs=>"forward_fastq_filesA", filtRs=>"reverse_fastq_filesA" 
#sample.names => sample_names_Amanda 
#filtpathF => forward_data_path_filtered_Amanda 
#pathF=>forward_data_path, pathR=>reverse_data_path 

# Functions # 
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
create_directory("./amandaFastq") 
create_directory("./amandaFastq/amandaF") 
create_directory("./amandaFastq/amandaR") 
create_directory("./amandaFastq/amandaF/filtered") 
create_directory("./amandaFastq/amandaR/filtered") 
path <- "./amandaFastq" 
forward_data_path <- "./amandaFastq/amandaF" 
reverse_data_path <- "./amandaFastq/amandaR" 
forward_data_path_filtered_Amanda <- file.path(forward_data_path, "filtered") #filtered forward files go into pathF/filtered subdirectory 
reverse_data_path_filtered_Amanda <- file.path(reverse_data_path, "filtered") #filtered rev files -> pathR filtered subdirectory 


# Checking quality of data # 
forward_read_quality_Amanda <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE)) 
reverse_read_quality_Amanda <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE)) 
message("Here is a graph of the quality of the reverse reads.") 
plotQualityProfile(reverse_read_quality_Amanda[1:3]) 
message("Here is a graph of the quality of the forward reads.") 
plotQualityProfile(forward_read_quality_Amanda[1:3]) 

# Filtering Forward and Reverse Files # 
#filters files based on R1 and R2. R1 means forward, R2 means reverse 
#dir() selects files, file.copy moves them 
#This code is based on: https://stackoverflow.com/questions/49593823/r-file-copy-function 
#Forward Reads 
forward_data_path <- "./amandaFastq/AmandaF" 
data_filesF <- dir(path, pattern = "*R1_001.fastq.gz",  
                   all.files = FALSE) 
file.copy(file.path(path, data_filesF), forward_data_path, overwrite=TRUE) 
#Reverse Reads 
reverse_data_path <- "./amandaFastq/AmandaR" 
data_filesR <- dir(path, pattern = "*R2_001.fastq.gz",  
                   all.files = FALSE) 
file.copy(file.path(path, data_filesR), reverse_data_path, overwrite=TRUE) 

# File Dissecting for Filtering # 
forward_data_path_filtered_Amanda <- file.path(forward_data_path,  
                                               "filtered") #filtered forward files go into pathF/filtered subdirectory 
reverse_data_path_filtered_Amanda <- file.path(reverse_data_path,  
                                               "filtered") #filtered rev files -> pathR filtered subdirectory 
forward_fastq_list_Amanda <- sort(list.files(forward_data_path,  
                                             pattern = "fastq.gz")) 
reverse_fastq_list_Amanda <- sort(list.files(reverse_data_path,  
                                             pattern = "fastq.gz")) 
if(length(forward_fastq_list_Amanda) != length(reverse_fastq_list_Amanda)) 
  stop("Forward and reverse reads do not match. :(") 

# Filtering # 
#Note: These Parameters aren't optimal for all data sets
out <- filterAndTrim(fwd=file.path(forward_data_path, forward_fastq_list_Amanda),
                     filt=file.path(forward_data_path_filtered_Amanda,forward_fastq_list_Amanda),
                     rev = file.path(reverse_data_path, reverse_fastq_list_Amanda),
                     filt.rev = file.path(reverse_data_path_filtered_Amanda, reverse_fastq_list_Amanda),
                     trimLeft = c(20,30), truncLen = c(240,200),
                     maxEE = 2, truncQ = 2,
                     maxN=0, rm.phix = TRUE, compress = TRUE, verbose = TRUE,
                     multithread = TRUE)
head(out) 

# File Dissecting for Inferring Sequence Variants # 
#Note: The workflow has two different lines of code that start with "File parsing". 
#I've changed that to File Dissecting, and I've included what it's doing it for 
#in the title. 
forward_fastq_filesA <- list.files(forward_data_path_filtered_Amanda, pattern="fastq.gz", full.names = TRUE) 
reverse_fastq_filesA <- list.files(reverse_data_path_filtered_Amanda, pattern = "fastq.gz", full.names = TRUE) 
sample_names_Amanda <- sapply(strsplit(basename(forward_fastq_filesA), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz 
sample_names_Amanda <- sapply(strsplit(basename(reverse_fastq_filesA), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz 
if(!identical(sample_names_Amanda, sample_names_Amanda))  
  stop("Forward and reverse reads don't match") 

# Infer Sequence Variants # 
names(forward_fastq_filesA) <- sample_names_Amanda 
names(reverse_fastq_filesA) <- sample_names_Amanda 
set.seed(100) 
# Learn forward & reverse error rates 
forward_error_rateA <- learnErrors(forward_fastq_filesA, nbases = 1e8, multithread = TRUE) 
reverse_error_rateA <- learnErrors(reverse_fastq_filesA, nbases = 1e8, multithread = TRUE) 

# Sample Inference & Merger of Paired-End Reads # 
mergers <- vector("list", length(sample_names_Amanda)) 
names(mergers) <- sample_names_Amanda
for(sam in sample_names_Amanda) { 
  cat("Processing:", sam, "\n") 
  derepF <- derepFastq(forward_fastq_filesA[[sam]]) 
  ddF <- dada(derepF, err = forward_error_rateA, multithread = TRUE) 
  derepR <- derepFastq(reverse_fastq_filesA[[sam]]) 
  ddR <- dada(derepR, err = reverse_error_rateA, multithread = TRUE) 
  merger <- mergePairs(ddF, derepF, ddR, derepR) 
  mergers[[sam]] <- merger 
} 
rm(derepF, derepR) 

# construct sequence table and remove chimeras # 
seqtabA <- makeSequenceTable(mergers) 
create_directory("./amandaFastq/results") 
saveRDS(seqtabA, "./amandaFastq/results/seqtabA.rds") 

# Merge Runs, Remove Chimeras, Assign Taxonomy # 
pathRDS <- "./amandaFastq/results/seqtabA.rds" # Where the RDS files were saved 
#Remove chimeras 
seqtabA <- removeBimeraDenovo(readRDS(pathRDS),  
                              method="consensus", multithread = TRUE) 
#Assign taxonomy 
#Note: Silva file comes from https://zenodo.org/records/1172783 
#silva_nr_v132_train_set.fa is the name of the file 
taxA <- assignTaxonomy(seqtabA, "./silva_nr_v132_train_set.fa.gz",
                       multithread = TRUE) 

# Write to disk 
saveRDS(seqtabA, "./amandaFastq/results/seqtabA.rds") # This is where the sequencing table file is saved 
saveRDS(taxA, "./amandaFastq/results/taxA.rds") # This is where the taxonomy file is saved 

# Visualization # 
#We used ggplot2 and stuff to make this graph, not phyloseq 
#Code based on Eric's code used to visualize his data 
#note that this code also calculates relative abundace in order to  
# visualize this data. 

#Rename first col of taxA > "ASVs"  
asv_taxA <- cbind(rownames(taxA), data.frame(taxA), row.names=NULL) %>% 
  rename("ASVs" = "rownames(taxA)") 

asv_seqtabA <- cbind(rownames(seqtabA), data.frame(seqtabA), row.names=NULL) %>% 
  rename("sample_id" = "rownames(seqtabA)") %>% 
  pivot_longer(!sample_id, names_to = "ASVs", values_to = "count") 
#merge otu and taxon 
composite <- inner_join(asv_seqtabA, asv_taxA, by="ASVs") 

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
inner_join(taxon_rel_abund, taxon_pool, by="taxon") %>% 
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

ggsave("amandaFastq/results/amandaDADA2_16S_stacked_bar.tiff", 
       width=9, height=4) 
