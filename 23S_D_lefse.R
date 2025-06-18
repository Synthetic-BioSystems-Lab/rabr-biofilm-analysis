library(tidyr)
library(ggtext)
library(readxl)
library(glue)

setwd("~/Miller Lab/Rscripts_PilotRABR")
#DADA2/23SDADA2visualization/")

gather_metadata <- function(target, page) {
  read_excel("DADA2/23SDADA2visualization/23S_D_metadata.xlsx", sheet=page) %>%
    rename_all(tolower) %>%
    mutate(sample = str_replace_all(sample, "-", "_")) %>%
    select(sample, target)
}

#import ASV counts
seqtabE <- readRDS("DADA2/23SDADA2visualization/seqtabE.rds")
seqtabA <- readRDS("DADA2/23SDADA2visualization/seqtabA.rds")
taxE <- readRDS("DADA2/23SDADA2visualization/taxE.rds")
taxA <- readRDS("DADA2/23SDADA2visualization/taxA.rds")

#Rename first col of taxE > "ASVs"  
asv_taxE <- cbind(rownames(taxE), data.frame(taxE), row.names=NULL) %>% 
  rename("ASVs" = "rownames(taxE)") 

asv_seqtabE <- cbind(rownames(seqtabE), data.frame(seqtabE), row.names=NULL) %>% 
  rename("sample_id" = "rownames(seqtabE)") %>% 
  pivot_longer(!sample_id, names_to = "ASVs", values_to = "count") 

asv_taxA <- cbind(rownames(taxA), data.frame(taxA), row.names=NULL) %>% 
  rename("ASVs" = "rownames(taxA)") 

asv_seqtabA <- cbind(rownames(seqtabA), data.frame(seqtabA), row.names=NULL) %>% 
  rename("sample_id" = "rownames(seqtabA)") %>% 
  pivot_longer(!sample_id, names_to = "ASVs", values_to = "count") 

#merge otu and taxon 
compositeA <- inner_join(asv_seqtabA, asv_taxA, by="ASVs")
compositeE <- inner_join(asv_seqtabE, asv_taxE, by="ASVs")
composite <- rbind(compositeE, compositeA)

#tax file
taxonomy <- composite %>%
  select("ASVs", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus") %>%
  rename_all(tolower)

taxonomy <- taxonomy[!(taxonomy$genus %in% NA),]
taxonomy <- taxonomy[!duplicated(taxonomy), ]
taxonomy <- taxonomy %>% mutate(rownumber = 1:nrow(taxonomy),
         taxon=glue("{genus} (ASV {rownumber})")) %>%
  select(asvs, taxon)


## dli analysis -> significantly different
dli_metadata <- gather_metadata("dli_level", 1)

# shared <- rbind(asv_seqtabA, asv_seqtabE)
# shared <- shared %>%
#   pivot_wider(names_from = ASVs, values_from = count, values_fill = 0) %>%
#   mutate(sample_id = str_replace_all(sample_id, "-", "_"))
shared <- asv_seqtabE %>%
  pivot_wider(names_from = ASVs, values_from = count, values_fill = 0) %>%
  mutate(sample_id = str_replace_all(sample_id, "-", "_"))%>%
  rename("Group"="sample_id")

# for(i in 2:ncol(shared)) {
#   name <- as.character(i-1)
#   while(nchar(name)<5){
#     name <- paste("0", name, sep="")
#   }
#   name <- paste("Otu", name, sep="") #Otu + 000s + (i-1)
#   names(shared)[i]<-paste(name)
# }

dli_shared_design <- inner_join(shared, dli_metadata, by=c("Group"="sample"))
numOtus <- ncol(dli_shared_design)-2
label <- 0.03
dli_shared_design <- cbind(label,dli_shared_design[,1],numOtus,dli_shared_design[,2:ncol(dli_shared_design)])

# comparisons

# high vs low
high_low <- dli_shared_design %>%
  filter(dli_level == "high" | dli_level == "low")

high_low %>%
  select(-dli_level) %>%
  write_tsv("DADA2/23SDADA2visualization/processed_data/23S.high_low.shared")

high_low %>%
  select(Group, dli_level) %>%
  write_tsv("DADA2/23SDADA2visualization/processed_data/23S.high_low.design")

# run mothur lefse from desktop
#in terminal
# mothur\mothur
#lefse(shared=23S.high_low.shared, design=23S.high_low.design, inputdir=DADA2/23SDADA2visualization/processed_data)
      


# make plots

# high vs low
read_tsv("DADA2/23SDADA2visualization/processed_data/23S.high_low.0.03.lefse_summary") %>%
  drop_na(LDA) %>%
  filter(LDA > 3) %>%
  inner_join(., taxonomy, by=c("OTU" = "asvs")) %>%
  mutate(LDA = if_else(Class == "low", -1 * LDA, LDA),
         taxon = fct_reorder(taxon, LDA)) %>%
  ggplot(aes(x=LDA, y=taxon, fill=Class)) +
  geom_col() +
  labs(y=NULL, x="LDA Score (log 10)", title="Discriminant Genera between the\nhigh and low DLI") +
  scale_x_continuous(limits = c(-6, 6), breaks=seq(-6, 6, by=2)) +
  theme_classic() +
  theme(axis.text.y = element_markdown(size=6),
        plot.title=element_text(hjust=0.5))

ggsave("DADA2/23SDADA2visualization/23S_D_plots/23SD_lefse_pilot_dli_high_low.tiff", width=6, height=5)


##LabRABR vs Pilot
## dli analysis -> significantly different

pilot_metadata <- gather_metadata("section", 1) 

temp_metadata <- gather_metadata("section", 2)

pilot_metadata <- dplyr::bind_rows(pilot_metadata, temp_metadata)

pilot_metadata %>%
  mutate(section = factor(section,
                          levels=c("pilot", "81RABR")))

shared <- dplyr::bind_rows(asv_seqtabA, asv_seqtabE) %>%
  pivot_wider(names_from = ASVs, values_from = count, values_fill = 0) %>%
  mutate(sample_id = str_replace_all(sample_id, "-", "_"))%>%
  rename("Group"="sample_id") 

shared[is.na(shared)] <- 0


for(i in 1:nrow(shared)) {
  row <- shared[i,]
  row[1] <- gsub("_23S", "", row[1])
  shared[i,1] <- row[1]
}
#shared <- shared[-1]

for(i in 1:nrow(pilot_metadata)) {
  row <- pilot_metadata[i,]
  row[1] <- gsub("_23S", "", row[1])
  pilot_metadata[i,1] <- row[1]
}

# # test
# for(i in 2:ncol(shared)) {
#   name <- as.character(i-1)
#   while(nchar(name)<5){
#     name <- paste("0", name, sep="")
#   }
#   name <- paste("Otu", name, sep="") #Otu + 000s + (i-1)
#   names(shared)[i]<-paste(name)
# }

pilot_shared_design <- inner_join(shared, pilot_metadata, by=c("Group"="sample"))
numOtus <- ncol(pilot_shared_design)-2
label <- 0.03
pilot_shared_design <- cbind(label,pilot_shared_design[,1],numOtus,pilot_shared_design[,2:ncol(pilot_shared_design)])

# comparisons

# LabRABR v Pilot
LabvPilot <- pilot_shared_design %>%
  filter(section == "pilot" | section == "81RABR")

LabvPilot %>%
  select(-section) %>%
  write_tsv("DADA2/23SDADA2visualization/processed_data/23S.labvPilot.shared")

LabvPilot %>%
  select(Group, section) %>%
  write_tsv("DADA2/23SDADA2visualization/processed_data/23S.labvPilot.design")

# run mothur lefse from desktop
#in terminal
# mothur\mothur
#lefse(shared=23S.labvPilot.shared, design=23S.labvPilot.design, inputdir=DADA2\23SDADA2visualization\processed_data)
 


# make plots
bob <- read_tsv("DADA2/23SDADA2visualization/processed_data/23S.labvPilot.0.03.lefse_summary") %>%
  drop_na(LDA) %>%
  filter(LDA > 4) %>%
  inner_join(., taxonomy, by=c("OTU" = "asvs")) %>%
  mutate(LDA = if_else(Class == "81RABR", -1 * LDA, LDA),
         taxon = fct_reorder(taxon, LDA))
bob <-bob[!duplicated(bob), ]
# Lab v Pilot
read_tsv("DADA2/23SDADA2visualization/processed_data/23S.labvPilot.0.03.lefse_summary") %>%
  drop_na(LDA) %>%
  filter(LDA > 4) %>%
  inner_join(., taxonomy, by=c("OTU" = "asvs")) %>%
  mutate(LDA = if_else(Class == "81RABR", -1 * LDA, LDA),
         taxon = fct_reorder(taxon, LDA)) %>%
  ggplot(aes(x=LDA, y=taxon, fill=Class)) +
  geom_col() +
  labs(y=NULL, x="LDA Score (log 10)", title="Discriminant Genera between the\nLab and Pilot RABRs") +
  scale_x_continuous(limits = c(-6, 6), breaks=seq(-6, 6, by=2)) +
  theme_classic() +
  theme(axis.text.y = element_markdown(size=12),
        plot.title=element_text(hjust=0.5))

ggsave("DADA2/23SDADA2visualization/23S_D_plots/23SD_lefse_labvpilot.tiff", width=6, height=5)

##TF vs Pilot
pilot_metadata <- gather_metadata("section", 1) 

temp_metadata <- gather_metadata("section", 4)

pilot_metadata <- dplyr::bind_rows(pilot_metadata, temp_metadata)

pilot_metadata %>%
  mutate(section = factor(section,
                          levels=c("pilot", "TF")))


for(i in 1:nrow(pilot_metadata)) {
  row <- pilot_metadata[i,]
  row[1] <- gsub("_23S", "", row[1])
  pilot_metadata[i,1] <- row[1]
}

pilot_shared_design <- inner_join(shared, pilot_metadata, by=c("Group"="sample"))
numOtus <- ncol(pilot_shared_design)-2
label <- 0.03
pilot_shared_design <- cbind(label,pilot_shared_design[,1],numOtus,pilot_shared_design[,2:ncol(pilot_shared_design)])


# comparisons

# TF v Pilot
TFvPilot <- pilot_shared_design %>%
  filter(section == "pilot" | section == "TF")

TFvPilot %>%
  select(-section) %>%
  write_tsv("DADA2/23SDADA2visualization/processed_data/23S.TFvPilot.shared")

TFvPilot %>%
  select(Group, section) %>%
  write_tsv("DADA2/23SDADA2visualization/processed_data/23S.TFvPilot.design")

# run mothur lefse from desktop
#in terminal
# mothur\mothur
#lefse(shared=23S.TFvPilot.shared, design=23S.TFvPilot.design, inputdir=DADA2\23SDADA2visualization\processed_data)
 


# make plots
# TF v Pilot
read_tsv("DADA2/23SDADA2visualization/processed_data/23S.TFvPilot.0.03.lefse_summary") %>%
  drop_na(LDA) %>%
  filter(LDA > 3.6) %>%
  inner_join(., taxonomy, by=c("OTU" = "asvs")) %>%
  mutate(LDA = if_else(Class == "TF", -1 * LDA, LDA),
         taxon = fct_reorder(taxon, LDA))  %>%
  ggplot(aes(x=LDA, y=taxon, fill=Class)) +
  geom_col() +
  labs(y=NULL, x="LDA Score (log 10)", title="Discriminant Genera between the\nTrickling Filter and Pilot RABRs") +
  scale_x_continuous(limits = c(-6, 6), breaks=seq(-6, 6, by=2)) +
  theme_classic() +
  theme(axis.text.y = element_markdown(size=6),
        plot.title=element_text(hjust=0.5))

ggsave("DADA2/23SDADA2visualization/23S_D_plots/23SD_lefse_TFvpilot.tiff", width=6, height=5)

##CVWRF vs Pilot
pilot_metadata <- gather_metadata("section", 1) 

temp_metadata <- gather_metadata("section", 3)

pilot_metadata <- dplyr::bind_rows(pilot_metadata, temp_metadata)

pilot_metadata %>%
  mutate(section = factor(section,
                          levels=c("pilot", "CVWRF")))
#shared <- shared[-1]

for(i in 1:nrow(pilot_metadata)) {
  row <- pilot_metadata[i,]
  row[1] <- gsub("_23S", "", row[1])
  pilot_metadata[i,1] <- row[1]
}

pilot_shared_design <- inner_join(shared, pilot_metadata, by=c("Group"="sample"))
numOtus <- ncol(pilot_shared_design)-2
label <- 0.03
pilot_shared_design <- cbind(label,pilot_shared_design[,1],numOtus,pilot_shared_design[,2:ncol(pilot_shared_design)])

# comparisons

# CVWRF v Pilot
CVWRFvPilot <- pilot_shared_design %>%
  filter(section == "pilot" | section == "CVWRF")

CVWRFvPilot %>%
  select(-section) %>%
  write_tsv("DADA2/23SDADA2visualization/processed_data/23S.CVWRFvPilot.shared")

CVWRFvPilot %>%
  select(Group, section) %>%
  write_tsv("DADA2/23SDADA2visualization/processed_data/23S.CVWRFvPilot.design")

# run mothur lefse from desktop
#in terminal
# mothur\mothur
#lefse(shared=23S.CVWRFvPilot.shared, design=23S.CVWRFvPilot.design, inputdir=DADA2\23SDADA2visualization\processed_data)
 


# make plots

# CVWRF v Pilot
read_tsv("DADA2/23SDADA2visualization/processed_data/23S.CVWRFvPilot.0.03.lefse_summary") %>%
  drop_na(LDA) %>%
  filter(LDA > 1) %>%
  inner_join(., taxonomy, by=c("OTU" = "asvs")) %>%
  mutate(LDA = if_else(Class == "CVWRF", -1 * LDA, LDA),
         taxon = fct_reorder(taxon, LDA)) %>%
  ggplot(aes(x=LDA, y=taxon, fill=Class)) +
  geom_col() +
  labs(y=NULL, x="LDA Score (log 10)", title="Discriminant Genera between the\nCVWRF and Pilot RABRs") +
  scale_x_continuous(limits = c(-6, 6), breaks=seq(-6, 6, by=2)) +
  theme_classic() +
  theme(axis.text.y = element_markdown(size=6),
        plot.title=element_text(hjust=0.5))

ggsave("DADA2/23SDADA2visualization/23S_D_plots/23SD_lefse_CVWRFvpilot.tiff", width=6, height=5)

#other comparisons?