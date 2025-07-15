library(tidyverse) 
library(ggplot2)
library(ggtext)
library(readxl)
library(NLP)
library(ggpubr)

setwd("~/Miller Lab/Rscripts_PilotRABR")

final_shared <- "16Spilot/final.opti_mcc.shared"
final_tax <- "16Spilot/final.opti_mcc.0.03.cons.taxonomy"

gather_metadata <- function(target, t2, t3, page) {
  read_excel("16Spilot/16S_metadata.xlsx", sheet=page) %>%
    rename_all(tolower) %>%
    mutate(sample = str_replace_all(sample, "-", "_")) %>%
    select(sample, target, t2, t3)
}

#import otu counts
otu_counts <- read_tsv(final_shared) %>%
  select(-label, -numOtus) %>%
  pivot_longer(-Group, names_to = "otu", values_to = "count") %>%
  rename(sample_id = Group)

# import taxonomy
taxonomy <- read_tsv(final_tax) %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""), taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), sep = ";")

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
  pivot_longer(cols = c("kingdom", "phylum", "class", "order", "family", "genus"), names_to = "level", values_to = "taxon")

otu_rel_abund <- inner_join(trimmed_composite, all_metadata,  by=c('sample_id'='sample'))

#function
tax = c("Synechococcus_PCC-7942", "Tychonema_CCAP_1459-11B", "Symphothece_PCC-7002")
lvl = "genus"
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
  
  write.csv(filtered_rel_abund,paste("16Spilot/16Scsvs/16S_StrainProd_", lvl, ".csv", sep=""), row.names = FALSE)
  
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
    stat_cor(label.x=xcor, label.y.npc=ycor) +  
    geom_smooth(method=lm, se=FALSE)
  
  ggsave(paste("16Spilot/16Splots/16S_RAvsProd_", lvl, ".tiff", sep=""), width=7, height=5)
  
}
tax = c("Synechococcus_PCC-7942",  "Symphothece_PCC-7002")
lvl = "genus"
rabr = 'pilot'

RA_Prod_RABR <- function(lvl, tax, rabr, xcor, ycor) {
  taxon_rel_abund <- otu_rel_abund %>%
    filter(section==rabr) %>%
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
  
  #write.csv(filtered_rel_abund,paste("16Spilot/16Scsvs/16S_StrainProd_", lvl, ".csv", sep=""), row.names = FALSE)
  
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
    stat_cor(label.x=xcor, label.y.npc=ycor) +
    geom_smooth(method=lm, se=FALSE)
  
  ggsave(paste("16Spilot/16Splots/16S_RAvsProd_", lvl, "_", rabr, ".tiff", sep=""), width=7, height=5)
  
}

RA_Prod("genus", c("Synechococcus_PCC-7942", "Tychonema_CCAP_1459-11B", "Symphothece_PCC-7002"), 15, 'top')
RA_Prod("family", c("Burkholderiaceae", "Weeksellaceae", "Synechococcaceae"), 15, 'top')
RA_Prod("order", c("Nostocales", "Betaproteobacteriales", "Flavobacteriales", "Synechococcales"), 15, 'center')
RA_Prod("class", c("Alphaproteobacteria", "Bacteroidia", "Gammaproteobacteria", "Oxyphotobacteria"), 15, 'center')
RA_Prod("phylum", c("Bacteroidetes", "Proteobacteria", "Cyanobacteria"), 15, 'center')

RA_Prod_RABR("genus", c("Synechococcus_PCC-7942", "Symphothece_PCC-7002"), 'pilot', 0, 'top')
RA_Prod_RABR("family", c("Burkholderiaceae", "Weeksellaceae", "Synechococcaceae"), 'pilot', 3.2, 'bottom')
RA_Prod_RABR("order", c("Nostocales", "Betaproteobacteriales", "Flavobacteriales", "Synechococcales"), 'pilot', 0, 'center')
RA_Prod_RABR("class", c("Alphaproteobacteria", "Bacteroidia", "Gammaproteobacteria", "Oxyphotobacteria"), 'pilot', 3.2, 'center')
RA_Prod_RABR("phylum", c("Bacteroidetes", "Proteobacteria", "Cyanobacteria"), 'pilot', 3.2, 'center')
RA_Prod_RABR("genus", c("Synechococcus_PCC-7942", "Tychonema_CCAP_1459-11B", "Symphothece_PCC-7002"), '81RABR', 16, 'top')
RA_Prod_RABR("family", c("Burkholderiaceae", "Weeksellaceae", "Synechococcaceae"), '81RABR', 16, 'top')
RA_Prod_RABR("order", c("Nostocales", "Betaproteobacteriales", "Flavobacteriales", "Synechococcales"), '81RABR', 16, 'center')
RA_Prod_RABR("class", c("Alphaproteobacteria", "Bacteroidia", "Gammaproteobacteria", "Oxyphotobacteria"), '81RABR', 16, 'center')
RA_Prod_RABR("phylum", c("Bacteroidetes", "Proteobacteria", "Cyanobacteria"), '81RABR', 16, 'center')

