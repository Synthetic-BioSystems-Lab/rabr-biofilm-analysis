library(tidyverse) 
library(ggplot2)
library(ggtext)
library(readxl)
library(NLP)
library(ggpubr)

setwd("~/Miller Lab/Rscripts_PilotRABR")

final_shared <- "18Spilotv2/final_p0.asv.shared"
final_tax <- "18Spilotv2/final_p0.asv.0.03.cons.taxonomy"

gather_metadata <- function(target, t2, t3, page) {
  read_excel("18Spilot/18S_metadata.xlsx", sheet=page) %>%
    rename_all(tolower) %>%
    mutate(sample = str_replace_all(sample, "-", "_")) %>%
    select(sample, target, t2, t3)
}

#import otu counts
otu_counts <- read_tsv(final_shared) %>%
  select(-label, -numASVs) %>%
  pivot_longer(-Group, names_to = "otu", values_to = "count") %>%
  rename(sample_id = Group)

# import taxonomy
taxonomy <- read_tsv(final_tax) %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""), taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy, into=c("domain", "supergroup", "division", "subdivision", "class", "order", "family", "genus", "species"), sep = ";")

#merge otu and taxon
composite <- inner_join(otu_counts, taxonomy, by="otu")

# import metadata
all_metadata <- gather_metadata("section", "dry_productivity_substratum", "date", 1) %>%
  rename(productivity = dry_productivity_substratum)

all_metadata <- dplyr::bind_rows(all_metadata, gather_metadata("section", "productivity", "date", 2))

all_metadata %>%
  mutate(section = factor(section,
                          levels=c("pilot", "81RABR")))

#remove empty rows, add relative abundance, pivot for readability
trimmed_composite <- composite[!(composite$count=="0"),] %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols = c("domain", "supergroup", "division", "subdivision", "class", "order", "family", "genus", "species"), names_to = "level", values_to = "taxon")

otu_rel_abund <- inner_join(trimmed_composite, all_metadata,  by=c('sample_id'='sample'))

#function
tax = c("Agaricomycotina", "Chlorellales", "Hypotrichia")
lvl = "order"
RA_Prod <- function(lvl, tax) {
  taxon_rel_abund <- otu_rel_abund %>%
    filter(level==lvl) %>%
    group_by(section, sample_id, taxon, date, productivity) %>%
    summarize(xrel_abund = 100*sum(rel_abund), .groups = "drop") %>%
    # group_by(sample_id, taxon, section) %>%
    # summarize(mean_rel_abund = 100*mean(rel_abund), .groups = "drop") %>%
    mutate(taxon = str_replace(taxon, "(.*)_unclassified", "Unclassified \\1"),
           taxon = str_replace(taxon, "^(\\S*)$", "\\1"),
           rel_abund = xrel_abund) %>%
    select(section, sample_id, taxon, date, productivity, rel_abund)
  
  #select for tax
  tax_count = length(tax)
  
  filtered_rel_abund <- taxon_rel_abund %>%
    filter(taxon==tax[1])
  for(i in 2:tax_count) {
    save <- taxon_rel_abund %>%
      filter(taxon==tax[i])
    filtered_rel_abund <- dplyr::bind_rows(filtered_rel_abund, save)
  }
  write.csv(filtered_rel_abund,paste("18Spilotv2/18Scsvs/18S_StrainProd_", lvl, ".csv", sep=""), row.names = FALSE)
  #graph
  filtered_rel_abund %>%
    ggplot(aes(x = productivity, y = rel_abund, color = taxon)) +
    geom_point(size=3) + 
    #scale_fill_manual(name = NULL, values = c(brewer.pal(6, "Dark2"), "gray")) +
    scale_fill_discrete(name="Organism") +
    scale_y_continuous(expand=c(0,0)) +
    # facet_grid(~section, scale="free_x", space="free", 
    #            labeller = labeller(section=pretty)) +
    labs(title=paste("Organism Relative Abundance vs Productivity for ", str_to_title(lvl), sep=""),
         x = "Productivity (g/m2/day)",
         y = "Relative Abundance") +
    theme_classic() +
    ylim(0, NA) +
    xlim(0, NA) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.text = element_markdown(), 
          axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0),
          legend.key.size = unit(10, "pt"),
          strip.background = element_blank(),
          strip.text = element_markdown()) +
    stat_cor() +
    geom_smooth(method=lm, se=FALSE)
  
  ggsave(paste("18Spilotv2/18Splots/18S_RAvsProd_", lvl, ".tiff", sep=""), width=7, height=5)
  
}

tax = c("Agaricomycotina", "Chlorellales", "Hypotrichia")
lvl = "order"
rabr = 'pilot'
x_spot = 3.2
RA_Prod_RABR <- function(lvl, tax, rabr, x_spot) {
  taxon_rel_abund <- otu_rel_abund %>%
    filter(level==lvl) %>%
    filter(section==rabr) %>%
    group_by(section, sample_id, taxon, date, productivity) %>%
    summarize(xrel_abund = 100*sum(rel_abund), .groups = "drop") %>%
    # group_by(sample_id, taxon, section) %>%
    # summarize(mean_rel_abund = 100*mean(rel_abund), .groups = "drop") %>%
    mutate(taxon = str_replace(taxon, "(.*)_unclassified", "Unclassified \\1"),
           taxon = str_replace(taxon, "^(\\S*)$", "\\1"),
           rel_abund = xrel_abund) %>%
    select(section, sample_id, taxon, date, productivity, rel_abund)
  
  #select for tax
  tax_count = length(tax)
  
  filtered_rel_abund <- taxon_rel_abund %>%
    filter(taxon==tax[1])
  for(i in 2:tax_count) {
    save <- taxon_rel_abund %>%
      filter(taxon==tax[i])
    filtered_rel_abund <- dplyr::bind_rows(filtered_rel_abund, save)
  }
  #write.csv(filtered_rel_abund,paste("18Spilotv2/18Scsvs/18S_StrainProd_", lvl, ".csv", sep=""), row.names = FALSE)
  #graph
  filtered_rel_abund %>%
    ggplot(aes(x = productivity, y = rel_abund, color = taxon)) +
    geom_point(size=3) + 
    #scale_fill_manual(name = NULL, values = c(brewer.pal(6, "Dark2"), "gray")) +
    scale_fill_discrete(name="Organism") +
    scale_y_continuous(expand=c(0,0)) +
    # facet_grid(~section, scale="free_x", space="free", 
    #            labeller = labeller(section=pretty)) +
    labs(title=paste("Organism Relative Abundance vs Productivity for ", str_to_title(lvl), sep=""),
         x = "Productivity (g/m2/day)",
         y = "Relative Abundance") +
    theme_classic() +
    ylim(0, NA) +
    xlim(0, NA) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.text = element_markdown(), 
          axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0),
          legend.key.size = unit(10, "pt"),
          strip.background = element_blank(),
          strip.text = element_markdown()) +
    stat_cor() +
    geom_smooth(method=lm, se=FALSE)
  
  ggsave(paste("18Spilotv2/18Splots/18S_RAvsProd_", lvl, "_", rabr, ".tiff", sep=""), width=7, height=5)
  
}


#"domain", "supergroup", "division", "subdivision", "class", "order", "family", "genus", "species"
RA_Prod("genus", c("Kahliella", "Micractinium", "Psychoda", "Trichosporon"))
RA_Prod("family", c("Chlorellales_X", "Insecta", "Oxytrichidae", "Pezizomycetes", "Platyophryida", "Sphaeropleales_X", "Tremellomycetes"))
RA_Prod("order", c("Agaricomycotina", "Chlorellales", "Hypotrichia"))
RA_Prod("class", c("Ascomycota", "Basidiomycota", "Chlorophyceae", "Colpodea", "Spirotrichea", "Trebouxiophyceae"))
RA_Prod("subdivision", c("Chlorophyta_X", "Ciliophora", "Fungi", "Gyrista", "Metazoa"))
RA_Prod("division", c("Alveolata", "Chlorophyta", "Opisthokonta", "Stramenopiles"))
RA_Prod("supergroup", c("Archaeplastida", "Obazoa", "TSAR"))

RA_Prod_RABR("genus", c("Kahliella", "Micractinium", "Psychoda", "Trichosporon"), 'pilot', 3.2)
RA_Prod_RABR("family", c("Chlorellales_X", "Insecta", "Oxytrichidae", "Pezizomycetes", "Platyophryida", "Sphaeropleales_X", "Tremellomycetes"), 'pilot', 3.2)
RA_Prod_RABR("order", c("Agaricomycotina", "Chlorellales", "Hypotrichia"), 'pilot', 3.2)
RA_Prod_RABR("class", c("Ascomycota", "Basidiomycota", "Chlorophyceae", "Colpodea", "Spirotrichea", "Trebouxiophyceae"), 'pilot', 3.2)
RA_Prod_RABR("subdivision", c("Chlorophyta_X", "Ciliophora", "Fungi", "Gyrista", "Metazoa"), 'pilot', 3.2)
RA_Prod_RABR("division", c("Alveolata", "Chlorophyta", "Opisthokonta", "Stramenopiles"), 'pilot', 3.2)
RA_Prod_RABR("supergroup", c("Archaeplastida", "Obazoa", "TSAR"), 'pilot', 3.2)
RA_Prod_RABR("genus", c("Kahliella", "Micractinium", "Psychoda", "Trichosporon"), '81RABR', 16)
RA_Prod_RABR("family", c("Chlorellales_X", "Insecta", "Oxytrichidae", "Pezizomycetes", "Platyophryida", "Sphaeropleales_X", "Tremellomycetes"), '81RABR', 16)
RA_Prod_RABR("order", c("Agaricomycotina", "Chlorellales", "Hypotrichia"), '81RABR', 16)
RA_Prod_RABR("class", c("Ascomycota", "Basidiomycota", "Chlorophyceae", "Colpodea", "Spirotrichea", "Trebouxiophyceae"), '81RABR', 16)
RA_Prod_RABR("subdivision", c("Chlorophyta_X", "Ciliophora", "Fungi", "Gyrista", "Metazoa"), '81RABR', 16)
RA_Prod_RABR("division", c("Alveolata", "Chlorophyta", "Opisthokonta", "Stramenopiles"), '81RABR', 16)
RA_Prod_RABR("supergroup", c("Archaeplastida", "Obazoa", "TSAR"), '81RABR', 16)
