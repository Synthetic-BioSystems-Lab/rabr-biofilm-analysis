library(tidyr)
library(ggtext)
library(tidyverse)
library(readxl)
library(glue)

setwd("~/Miller Lab/Rscripts_PilotRABR")

# taxonomy file
taxonomy <- read_tsv("16Spilot/final.opti_mcc.0.03.cons.taxonomy") %>%
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

gather_metadata <- function(target, page) {
  read_excel("16Spilot/16S_metadata.xlsx", sheet=page) %>%
    rename_all(tolower) %>%
    mutate(sample = str_replace_all(sample, "-", "_")) %>%
    select(sample, target)
}

## dli analysis -> significantly different
dli_metadata <- gather_metadata("dli_level", 1)

shared <- read_tsv("16Spilot/final.opti_mcc.0.03.subsample.shared") 

for(i in 1:nrow(shared)) {
  row <- shared[i,]
  row[2] <- gsub("_16S", "", row[2])
  shared[i,2] <- row[2]
}
#shared <- shared[-1]

for(i in 1:nrow(dli_metadata)) {
  row <- dli_metadata[i,]
  row[1] <- gsub("_16S", "", row[1])
  dli_metadata[i,1] <- row[1]
}

dli_shared_design <- inner_join(shared, dli_metadata, by=c("Group"="sample"))

# comparisons

# high vs low
high_low <- dli_shared_design %>%
  filter(dli_level == "high" | dli_level == "low")

high_low %>%
  select(-dli_level) %>%
  write_tsv("16Spilot/processed_data/16S.high_low.shared")

high_low %>%
  select(Group, dli_level) %>%
  write_tsv("16Spilot/processed_data/16S.high_low.design")

# run mothur lefse from desktop
#in terminal
# mothur\mothur
#lefse(shared=16S.high_low.shared, design=16S.high_low.design, inputdir=16Spilot\processed_data)


# make plots

# high vs low
read_tsv("mothur/16S.high_low.0.03.lefse_summary") %>%
  drop_na(LDA) %>%
  filter(LDA > 3.2) %>%
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

ggsave("16Spilot/16Splots/16S_lefse_pilot_dli_high_low.tiff", width=6, height=5)


##LabRABR vs Pilot
## dli analysis -> significantly different

gather_metadata <- function(target, page) {
  read_excel("16Spilot/16S_metadata.xlsx", sheet=page) %>%
    rename_all(tolower) %>%
    mutate(sample = str_replace_all(sample, "-", "_")) %>%
    select(sample, target)
}

pilot_metadata <- gather_metadata("section", 1) 

temp_metadata <- gather_metadata("section", 2)

pilot_metadata <- dplyr::bind_rows(pilot_metadata, temp_metadata)

pilot_metadata %>%
  mutate(section = factor(section,
                          levels=c("pilot", "81RABR")))

shared <- read_tsv("16Spilot/final.opti_mcc.0.03.subsample.shared") 



for(i in 1:nrow(shared)) {
  row <- shared[i,]
  row[2] <- gsub("_16S", "", row[2])
  shared[i,2] <- row[2]
}
#shared <- shared[-1]

for(i in 1:nrow(pilot_metadata)) {
  row <- pilot_metadata[i,]
  row[1] <- gsub("_16S", "", row[1])
  pilot_metadata[i,1] <- row[1]
}

pilot_shared_design <- inner_join(shared, pilot_metadata, by=c("Group"="sample"))

# comparisons

# LabRABR v Pilot
LabvPilot <- pilot_shared_design %>%
  filter(section == "pilot" | section == "81RABR")

LabvPilot %>%
  select(-section) %>%
  write_tsv("16Spilot/processed_data/16S.labvPilot.shared")

LabvPilot %>%
  select(Group, section) %>%
  write_tsv("16Spilot/processed_data/16S.labvPilot.design")

# run mothur lefse from desktop
#in terminal
# mothur\mothur
#lefse(shared=16S.labvPilot.shared, design=16S.labvPilot.design, inputdir=16Spilot\processed_data)
#NOT WORKING


# make plots

# Lab v Pilot
read_tsv("16Spilot/processed_data/16S.labvPilot.0.03.lefse_summary") %>%
  drop_na(LDA) %>%
  filter(LDA > 4.1) %>%
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

ggsave("16Spilot/16Splots/16S_lefse_labvpilot.tiff", width=6, height=5)

##TF vs Pilot

gather_metadata <- function(target, page) {
  read_excel("16Spilot/16S_metadata.xlsx", sheet=page) %>%
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

shared <- read_tsv("16Spilot/final.opti_mcc.0.03.subsample.shared") 



for(i in 1:nrow(shared)) {
  row <- shared[i,]
  row[2] <- gsub("_16S", "", row[2])
  shared[i,2] <- row[2]
}
#shared <- shared[-1]

for(i in 1:nrow(pilot_metadata)) {
  row <- pilot_metadata[i,]
  row[1] <- gsub("_16S", "", row[1])
  pilot_metadata[i,1] <- row[1]
}

pilot_shared_design <- inner_join(shared, pilot_metadata, by=c("Group"="sample"))

# comparisons

# TF v Pilot
TFvPilot <- pilot_shared_design %>%
  filter(section == "pilot" | section == "TF")

TFvPilot %>%
  select(-section) %>%
  write_tsv("16Spilot/processed_data/16S.TFvPilot.shared")

TFvPilot %>%
  select(Group, section) %>%
  write_tsv("16Spilot/processed_data/16S.TFvPilot.design")

# run mothur lefse from desktop
#in terminal
# mothur\mothur
#lefse(shared=16S.TFvPilot.shared, design=16S.TFvPilot.design, inputdir=16Spilot\processed_data)
#NOT WORKING


# make plots

# TF v Pilot
read_tsv("16Spilot/processed_data/16S.TFvPilot.0.03.lefse_summary") %>%
  drop_na(LDA) %>%
  filter(LDA > 4.13) %>%
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

ggsave("16Spilot/16Splots/16S_lefse_TFvpilot.tiff", width=6, height=5)

##CVWRF vs Pilot

gather_metadata <- function(target, page) {
  read_excel("16Spilot/16S_metadata.xlsx", sheet=page) %>%
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

shared <- read_tsv("16Spilot/final.opti_mcc.0.03.subsample.shared") 



for(i in 1:nrow(shared)) {
  row <- shared[i,]
  row[2] <- gsub("_16S", "", row[2])
  shared[i,2] <- row[2]
}
#shared <- shared[-1]

for(i in 1:nrow(pilot_metadata)) {
  row <- pilot_metadata[i,]
  row[1] <- gsub("_16S", "", row[1])
  pilot_metadata[i,1] <- row[1]
}

pilot_shared_design <- inner_join(shared, pilot_metadata, by=c("Group"="sample"))

# comparisons

# CVWRF v Pilot
CVWRFvPilot <- pilot_shared_design %>%
  filter(section == "pilot" | section == "CVWRF")

CVWRFvPilot %>%
  select(-section) %>%
  write_tsv("16Spilot/processed_data/16S.CVWRFvPilot.shared")

CVWRFvPilot %>%
  select(Group, section) %>%
  write_tsv("16Spilot/processed_data/16S.CVWRFvPilot.design")

# run mothur lefse from desktop
#in terminal
# mothur\mothur
#lefse(shared=16S.CVWRFvPilot.shared, design=16S.CVWRFvPilot.design, inputdir=16Spilot\processed_data)
#NOT WORKING


# make plots

# CVWRF v Pilot
read_tsv("16Spilot/processed_data/16S.CVWRFvPilot.0.03.lefse_summary") %>%
  drop_na(LDA) %>%
  filter(LDA > 4.2) %>%
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

ggsave("16Spilot/16Splots/16S_lefse_CVWRFvpilot.tiff", width=6, height=5)

#other comparisons?