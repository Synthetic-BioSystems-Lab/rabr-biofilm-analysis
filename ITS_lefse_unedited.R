library(tidyr)
library(ggtext)
library(readxl)
library(glue)
library(dplyr)
library(tidyverse)

setwd("~/Miller Lab/Rscripts_PilotRABR")

gather_metadata <- function(target, page) {
  read_excel("ITSpilot/ITS_metadata.xlsx", sheet=page) %>%
    rename_all(tolower) %>%
    mutate(sample = str_replace_all(sample, "-", "_")) %>%
    select(sample, target)
}

#name <- "pr2_euk" #done
#name <-  "pr2_fungi" #done
#name <-  "silv_euk" #done
name <-  "silv_fungi" #done #run for each lineup
# taxonomy file
taxonomy <- read_tsv(paste("ITSpilot/final_", name, ".agc.0.05.cons.taxonomy", sep="")) %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus"),
           sep=";") %>%
  mutate(pretty_otu = str_replace(string=otu,
                                  pattern="tu0*",
                                  replacement = "TU "),
         genus = str_replace(string=genus,
                             pattern="(.*)",
                             replacement="*\\1*"),
         genus = str_replace(string=genus,
                             pattern="\\*(.*)_unclassified\\*",
                             replacement="Unclassified<br>*\\1*"),
         taxon = glue("{genus}<br>({pretty_otu})")) %>%
  select(otu, taxon)

## dli analysis -> significantly different
dli_metadata <- gather_metadata("dli_level", 1)

shared <- read_tsv(paste("ITSpilot/final_",name,".agc.shared", sep="")) 

dli_shared_design <- inner_join(shared, dli_metadata, by=c("Group"="sample"))

# comparisons

# high vs low
high_low <- dli_shared_design %>%
  filter(dli_level == "high" | dli_level == "low")

high_low %>%
  select(-dli_level) %>%
  write_tsv("ITSpilot/processed_data/ITS.high_low.shared")

high_low %>%
  select(Group, dli_level) %>%
  write_tsv("ITSpilot/processed_data/ITS.high_low.design")

# run mothur lefse from desktop
#in terminal
# mothur\mothur
#lefse(shared=ITS.high_low.shared, design=ITS.high_low.design, inputdir=ITSpilot\processed_data)
     #NOT WORKING


# make plots

# high vs low
read_tsv("ITSpilot/processed_data/ITS.high_low.0.05.lefse_summary") %>%
  drop_na(LDA) %>%
  filter(LDA > 3) %>%
  inner_join(., taxonomy, by=c("OTU" = "otu")) %>%
  mutate(LDA = if_else(Class == "low", -1 * LDA, LDA),
         taxon = fct_reorder(taxon, LDA)) %>%
  ggplot(aes(x=LDA, y=taxon, fill=Class)) +
  geom_col() +
  labs(y=NULL, x="LDA Score (log 10)", title="Discriminant Genera between the\nhigh and low DLI") +
  scale_x_continuous(limits = c(-6, 6), breaks=seq(-6, 6, by=2)) +
  theme_classic() +
  theme(axis.text.y = element_markdown(size=6),
        plot.title=element_text(hjust=0.5))

ggsave(paste("ITSpilot/ITSplots/ITS_", name, "lefse_pilot_dli_high_low.tiff", sep=""), width=6, height=5)


##LabRABR vs Pilot
## dli analysis -> significantly different

pilot_metadata <- gather_metadata("section", 1) 

temp_metadata <- gather_metadata("section", 2)

pilot_metadata <- dplyr::bind_rows(pilot_metadata, temp_metadata)

pilot_metadata %>%
  mutate(section = factor(section,
                          levels=c("pilot", "81RABR")))

shared <- read_tsv(paste("ITSpilot/final_",name,".agc.shared", sep="")) 



for(i in 1:nrow(shared)) {
  row <- shared[i,]
  row[2] <- gsub("_ITS", "", row[2])
  shared[i,2] <- row[2]
}
#shared <- shared[-1]

for(i in 1:nrow(pilot_metadata)) {
  row <- pilot_metadata[i,]
  row[1] <- gsub("_ITS", "", row[1])
  pilot_metadata[i,1] <- row[1]
}

pilot_shared_design <- inner_join(shared, pilot_metadata, by=c("Group"="sample"))

# comparisons

# LabRABR v Pilot
LabvPilot <- pilot_shared_design %>%
  filter(section == "pilot" | section == "81RABR")

LabvPilot %>%
  select(-section) %>%
  write_tsv("ITSpilot/processed_data/ITS.labvPilot.shared")

LabvPilot %>%
  select(Group, section) %>%
  write_tsv("ITSpilot/processed_data/ITS.labvPilot.design")

# run mothur lefse from desktop
#in terminal
# mothur\mothur
#lefse(shared=ITS.labvPilot.shared, design=ITS.labvPilot.design, inputdir=ITSpilot\processed_data)
#NOT WORKING


# make plots

# Lab v Pilot
read_tsv("ITSpilot/processed_data/ITS.labvPilot.0.05.lefse_summary") %>%
  drop_na(LDA) %>%
  filter(LDA > 3.3) %>%
  inner_join(., taxonomy, by=c("OTU" = "otu")) %>%
  mutate(LDA = if_else(Class == "81RABR", -1 * LDA, LDA),
         taxon = fct_reorder(taxon, LDA)) %>%
  ggplot(aes(x=LDA, y=taxon, fill=Class)) +
  geom_col() +
  labs(y=NULL, x="LDA Score (log 10)", title="Discriminant Genera between the\nLab and Pilot RABRs") +
  scale_x_continuous(limits = c(-6, 6), breaks=seq(-6, 6, by=2)) +
  theme_classic() +
  theme(axis.text.y = element_markdown(size=6),
        plot.title=element_text(hjust=0.5))

ggsave(paste("ITSpilot/ITSplots/ITS_", name, "lefse_labvpilot.tiff", sep=""), width=6, height=5)

##TF vs Pilot

gather_metadata <- function(target, page) {
  read_excel("ITSpilot/ITS_metadata.xlsx", sheet=page) %>%
    rename_all(tolower) %>%
    mutate(sample = str_replace_all(sample, "-", "_")) %>%
    select(sample, target)
}

pilot_metadata <- gather_metadata("section", 1) 

temp_metadata <- gather_metadata("section", 4)

pilot_metadata <- dplyr::bind_rows(pilot_metadata, temp_metadata)

pilot_metadata %>%
  mutate(section = factor(section,
                          levels=c("pilot", "TF")))

shared <- read_tsv(paste("ITSpilot/final_",name,".agc.shared", sep="")) 



for(i in 1:nrow(shared)) {
  row <- shared[i,]
  row[2] <- gsub("_ITS", "", row[2])
  shared[i,2] <- row[2]
}
#shared <- shared[-1]

for(i in 1:nrow(pilot_metadata)) {
  row <- pilot_metadata[i,]
  row[1] <- gsub("_ITS", "", row[1])
  pilot_metadata[i,1] <- row[1]
}

pilot_shared_design <- inner_join(shared, pilot_metadata, by=c("Group"="sample"))

# comparisons

# TF v Pilot
TFvPilot <- pilot_shared_design %>%
  filter(section == "pilot" | section == "TF")

TFvPilot %>%
  select(-section) %>%
  write_tsv("ITSpilot/processed_data/ITS.TFvPilot.shared")

TFvPilot %>%
  select(Group, section) %>%
  write_tsv("ITSpilot/processed_data/ITS.TFvPilot.design")

# run mothur lefse from desktop
#in terminal
# mothur\mothur
#lefse(shared=ITS.TFvPilot.shared, design=ITS.TFvPilot.design, inputdir=ITSpilot\processed_data)
#NOT WORKING


# make plots

# TF v Pilot
read_tsv("ITSpilot/processed_data/ITS.TFvPilot.0.05.lefse_summary") %>%
  drop_na(LDA) %>%
  filter(LDA > 3.5) %>%
  inner_join(., taxonomy, by=c("OTU" = "otu")) %>%
  mutate(LDA = if_else(Class == "TF", -1 * LDA, LDA),
         taxon = fct_reorder(taxon, LDA)) %>%
  ggplot(aes(x=LDA, y=taxon, fill=Class)) +
  geom_col() +
  labs(y=NULL, x="LDA Score (log 10)", title="Discriminant Genera between the\nTrickling Filter and Pilot RABRs") +
  scale_x_continuous(limits = c(-6, 6), breaks=seq(-6, 6, by=2)) +
  theme_classic() +
  theme(axis.text.y = element_markdown(size=6),
        plot.title=element_text(hjust=0.5))

ggsave(paste("ITSpilot/ITSplots/ITS_", name, "lefse_TFvpilot.tiff", sep=""), width=6, height=5)

##CVWRF vs Pilot

gather_metadata <- function(target, page) {
  read_excel("ITSpilot/ITS_metadata.xlsx", sheet=page) %>%
    rename_all(tolower) %>%
    mutate(sample = str_replace_all(sample, "-", "_")) %>%
    select(sample, target)
}

pilot_metadata <- gather_metadata("section", 1) 

temp_metadata <- gather_metadata("section", 3)

pilot_metadata <- dplyr::bind_rows(pilot_metadata, temp_metadata)

pilot_metadata %>%
  mutate(section = factor(section,
                          levels=c("pilot", "CVWRF")))

shared <- read_tsv(paste("ITSpilot/final_",name,".agc.shared", sep="")) 



for(i in 1:nrow(shared)) {
  row <- shared[i,]
  row[2] <- gsub("_ITS", "", row[2])
  shared[i,2] <- row[2]
}
#shared <- shared[-1]

for(i in 1:nrow(pilot_metadata)) {
  row <- pilot_metadata[i,]
  row[1] <- gsub("_ITS", "", row[1])
  pilot_metadata[i,1] <- row[1]
}

pilot_shared_design <- inner_join(shared, pilot_metadata, by=c("Group"="sample"))

# comparisons

# CVWRF v Pilot
CVWRFvPilot <- pilot_shared_design %>%
  filter(section == "pilot" | section == "CVWRF")

CVWRFvPilot %>%
  select(-section) %>%
  write_tsv("ITSpilot/processed_data/ITS.CVWRFvPilot.shared")

CVWRFvPilot %>%
  select(Group, section) %>%
  write_tsv("ITSpilot/processed_data/ITS.CVWRFvPilot.design")

# run mothur lefse from desktop
#in terminal
# mothur\mothur
#lefse(shared=ITS.CVWRFvPilot.shared, design=ITS.CVWRFvPilot.design, inputdir=ITSpilot\processed_data)
#NOT WORKING


# make plots

# CVWRF v Pilot
read_tsv("ITSpilot/processed_data/ITS.CVWRFvPilot.0.05.lefse_summary") %>%
  drop_na(LDA) %>%
  filter(LDA > .1) %>%
  inner_join(., taxonomy, by=c("OTU" = "otu")) %>%
  mutate(LDA = if_else(Class == "CVWRF", -1 * LDA, LDA),
         taxon = fct_reorder(taxon, LDA)) %>%
  ggplot(aes(x=LDA, y=taxon, fill=Class)) +
  geom_col() +
  labs(y=NULL, x="LDA Score (log 10)", title="Discriminant Genera between the\nCVWRF and Pilot RABRs") +
  scale_x_continuous(limits = c(-6, 6), breaks=seq(-6, 6, by=2)) +
  theme_classic() +
  theme(axis.text.y = element_markdown(size=6),
        plot.title=element_text(hjust=0.5))

ggsave(paste("ITSpilot/ITSplots/ITS_", name, "lefse_CVWRFvpilot.tiff", sep=""), width=6, height=5)

#other comparisons?

