library(tidyverse) 
library(ggplot2)
library(ggtext)
library(readxl)
library(NLP)
library(ggpubr)

setwd("~/Miller Lab/Rscripts_PilotRABR")

name <- "pr2_euk"

gather_metadata <- function(target, t2, t3, page) {
  read_excel("ITSfinal2/ITS_metadata.xlsx", sheet=page) %>%
    rename_all(tolower) %>%
    mutate(sample = str_replace_all(sample, "-", "_")) %>%
    select(sample, target, t2, t3)
}

#import otu counts
otu_counts <- read_tsv(paste("ITSfinal2/final_",name,".agc.shared", sep="")) %>%
  select(-label, -numOtus) %>%
  pivot_longer(-Group, names_to = "otu", values_to = "count") %>%
  rename(sample_id = Group)

# import taxonomy
taxonomy <- read_tsv(paste("ITSfinal2/final_", name, ".agc.0.05.cons.taxonomy", sep="")) %>%
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
tax = c("Chlamydomonadales", "Chlorellales", "Sphaeropleales", "Unclassified Chlorophyceae", "Unclassified Eukaryota")
lvl = "order"
xcor = 1
ycor = 'top'
RA_Prod <- function(lvl, tax, xcor, ycor) {
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
  
  write.csv(filtered_rel_abund,paste("ITSfinal2/ITScsvs/IT_pr2_euk_StrainProd_", lvl, ".csv", sep=""), row.names = FALSE)
  
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
    geom_smooth(method=lm, se=FALSE) +
    stat_cor(label.x=xcor, label.y.npc=ycor)
  
  ggsave(paste("ITSfinal2/ITSplots/ITS_RAvsProd_", lvl, ".tiff", sep=""), width=7, height=5)
  
}

tax = c("Chlamydomonadales", "Chlorellales", "Sphaeropleales", "Unclassified Eukaryota")
lvl = "order"
rabr = '81RABR'
x_spot = 16
RA_Prod_RABR <- function(lvl, tax, rabr, xcor, ycor) {
  taxon_rel_abund <- otu_rel_abund %>%
    filter(level==lvl) %>%
    filter(section==rabr)%>%
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
  
  #write.csv(filtered_rel_abund,paste("ITSfinal2/ITScsvs/IT_pr2_euk_StrainProd_", lvl, ".csv", sep=""), row.names = FALSE)
  
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
          strip.text = element_markdown())  +
    geom_smooth(method=lm, se=FALSE) +
    stat_cor(label.x=xcor, label.y.npc=ycor)
  
  ggsave(paste("ITSfinal2/ITSplots/ITS_RAvsProd_", lvl, "_", rabr, ".tiff", sep=""), width=7, height=5)
  
}


#"domain", "supergroup", "division", "subdivision", "class", "order", "family", "genus", "species"
RA_Prod("genus", c("Chlorella", "Ettlia", "Unclassified Eukaryota"), 15, 'center')
RA_Prod("family", c("Chlamydomonadales_X", "Chlorellales_X", "Sphaeropleales_X", "Unclassified Eukaryota"), 15, 'center')
RA_Prod("order", c("Chlamydomonadales", "Chlorellales", "Sphaeropleales", "Unclassified Chlorophyceae", "Unclassified Eukaryota"), 15, 'center')
RA_Prod("class", c("Chlorophyceae", "Trebouxiophyceae", "Unclassified Eukaryota"), 15, 'center')
RA_Prod("subdivision", c("Chlorophyta_X", "Proteobacteria_X", "Unclassified Eukaryota"), 15, 'center')
RA_Prod("division", c("Chlorophyta", "Proteobacteria", "Unclassified Eukaryota"), 15, 'center')
RA_Prod("supergroup", c("Archaeplastida", "PANNAM", "Unclassified Eukaryota"), 15, 'center')

RA_Prod_RABR("genus", c("Chlorella", "Ettlia", "Unclassified Eukaryota"), 'pilot', 1.5, 'top')
RA_Prod_RABR("family", c("Chlamydomonadales_X", "Chlorellales_X", "Sphaeropleales_X", "Unclassified Eukaryota"), 'pilot', 1, 'center')
RA_Prod_RABR("order", c("Chlamydomonadales", "Chlorellales", "Sphaeropleales", "Unclassified Chlorophyceae", "Unclassified Eukaryota"), 'pilot', 1, 'center')
RA_Prod_RABR("class", c("Chlorophyceae", "Trebouxiophyceae", "Unclassified Eukaryota"), 'pilot', 1.5, 'top')
RA_Prod_RABR("subdivision", c("Chlorophyta_X", "Proteobacteria_X", "Unclassified Eukaryota"), 'pilot', 1, 'center')
RA_Prod_RABR("division", c("Chlorophyta", "Proteobacteria", "Unclassified Eukaryota"), 'pilot', 0.5, 'center')
RA_Prod_RABR("supergroup", c("Archaeplastida", "PANNAM", "Unclassified Eukaryota"), 'pilot', 1, 'center')
RA_Prod_RABR("genus", c("Chlorella", "Ettlia", "Unclassified Eukaryota"), "81RABR", 1, 'center')
RA_Prod_RABR("family", c("Chlamydomonadales_X", "Chlorellales_X", "Sphaeropleales_X", "Unclassified Eukaryota"), "81RABR", 1, 'center')
RA_Prod_RABR("order", c("Chlamydomonadales", "Chlorellales", "Sphaeropleales", "Unclassified Eukaryota"), "81RABR", 1, 'center')
RA_Prod_RABR("class", c("Chlorophyceae", "Trebouxiophyceae", "Unclassified Eukaryota"), "81RABR", 1, 'center')
RA_Prod_RABR("subdivision", c("Chlorophyta_X", "Proteobacteria_X", "Unclassified Eukaryota"), "81RABR", 1, 'center')
RA_Prod_RABR("division", c("Chlorophyta", "Proteobacteria", "Unclassified Eukaryota"), "81RABR", 1, 'center')
RA_Prod_RABR("supergroup", c("Archaeplastida", "PANNAM", "Unclassified Eukaryota"), "81RABR", 1, 'center')
