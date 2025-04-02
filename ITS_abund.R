library(tidyverse) 
library(ggplot2)
library(ggtext)
library(readxl)


#relative abund stacked bar graph
#relative abund heat map
#relative abund bar graph

setwd("~/Miller Lab/Rscripts_PilotRABR")

gather_metadata <- function(target, t2, t3, page) {
  read_excel("ITSpilot/ITS_metadata.xlsx", sheet=page) %>%
    rename_all(tolower) %>%
    mutate(sample = str_replace_all(sample, "-", "_")) %>%
    select(sample, target, t2, t3)
}
name <- "pr2_euk"
Tlevel <- "phylum"


abund <- function(name, Tlevel) {
  #pr2 euk
  #import otu counts
  otu_counts <- read_tsv(paste("ITSpilot/final_",name,".agc.shared", sep="")) %>%
    select(-label, -numOtus) %>%
    pivot_longer(-Group, names_to = "otu", values_to = "count") %>%
    rename(sample_id = Group)
  
  # import taxonomy
  taxonomy <- read_tsv(paste("ITSpilot/final_", name, ".agc.0.05.cons.taxonomy", sep="")) %>%
    select("OTU", "Taxonomy") %>%
    rename_all(tolower) %>%
    mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""), taxonomy = str_replace(taxonomy, ";$", "")) %>%
    separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), sep = ";")
  
  # import metadata
  all_metadata <- gather_metadata("section", "date", "label", 1) 
  
  for(i in 2:6) {
    all_metadata <- dplyr::bind_rows(all_metadata, gather_metadata("section", "date", "label", i))
  }
  
  all_metadata %>%
    mutate(section = factor(section,
                            levels=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control")))
  
  #merge otu and taxon
  composite <- inner_join(otu_counts, taxonomy, by="otu")
  
  #remove empty rows, add relative abundance, pivot for readability
  trimmed_composite <- composite[!(composite$count=="0"),] %>%
    group_by(sample_id) %>%
    mutate(rel_abund = count / sum(count)) %>%
    ungroup() %>%
    pivot_longer(cols = c("kingdom", "phylum", "class", "order", "family", "genus"), names_to = "level", values_to = "taxon")
  
  sample_order <- all_metadata %>%
    arrange(date) %>%
    mutate(order = 1:nrow(.)) %>%
    select(sample, order, date)
  
  otu_rel_abund <- inner_join(trimmed_composite, all_metadata,  by=c('sample_id'='sample'))
  
  taxon_rel_abund <- otu_rel_abund %>%
    filter(level==Tlevel) %>%
    group_by(section, sample_id, taxon, date, label) %>%
    summarize(rel_abund = 100*sum(rel_abund), .groups = "drop") %>%
    # group_by(sample_id, taxon) %>%
    # summarize(mean_rel_abund = 100*mean(rel_abund), .groups = "drop") %>%
    mutate(taxon = str_replace(taxon, "(.*)_unclassified", "Unclassified *\\1*"),
           taxon = str_replace(taxon, "^(\\S*)$", "*\\1*"))
  
  #set aside other
  taxon_pool <- taxon_rel_abund %>%
    group_by(section, taxon, rel_abund) %>%
    summarize(mean=mean(rel_abund), .groups="drop") %>%
    group_by(taxon) %>%
    summarize(pool = max(mean) < 5, 
              mean = mean(mean),
              .groups="drop")
  # # metadata
  # metadata <- read_excel("16Spilot/16S_Metadata.xlsx") %>%
  #   select(Sample, DLI_Level) 
  
  # metadata <- metadata[3:10,]
  pretty <- c("pilot" = "Pilot-scale RABRs",
              "81RABR" = "Lab-scale RABRs",
              "control" = "Control",
              "CVWRF" = "CVWRF",
              "GHR" = "GHR",
              "TF" = "Trickling Filter")
  
  #assemble others and make RA stacked bar plot
  prep <- inner_join(taxon_rel_abund, taxon_pool, by="taxon") %>%
    mutate(taxon = if_else(pool, "Other", taxon)) %>%
    group_by(sample_id, label, section, taxon) %>%
    summarize(rel_abund = sum(rel_abund), .groups = "drop") %>%
    inner_join(., sample_order, by=c("sample_id"="sample")) %>%
    mutate(label = factor(label),
           label= fct_reorder(label, order))
  
  prep %>%
    ggplot(aes(x = label, y = rel_abund, fill = taxon)) +
    geom_col() + 
    #scale_fill_manual(name = NULL, values = c(brewer.pal(6, "Dark2"), "gray")) +
    scale_fill_discrete(name=NULL) +
    scale_y_continuous(expand=c(0,0)) +
    facet_grid(~section, scale="free_x", space="free", 
               labeller = labeller(section=pretty)) +
    labs(title=paste(str_to_title(Tlevel), " Relative Abundance of RABRs"),
         x = NULL,
         y = "Mean Relative Abundance (%)") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.text = element_markdown(), 
          axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0),
          legend.key.size = unit(10, "pt"),
          strip.background = element_blank(),
          strip.text = element_markdown())
  
  ggsave(paste("ITSpilot/ITSplots/ITS_", name, "_stacked_bar_", Tlevel, ".tiff", sep=""), width=9, height=4)
  
  # Relative Abundance
  
  # Heat Map
  prep %>%
    ggplot(aes(x=label, fill = rel_abund, y = taxon)) +
    geom_tile() +
    geom_text(aes(label = format(round(rel_abund, 0), nsmal=0))) +
    #scale_fill_manual(name=NULL, values = c(brewer.pal(6, "Dark2"), "gray")) +
    scale_x_discrete(expand= c(0,0),
                     position = "top") +
    scale_fill_gradient(name="Mean<br>Relative<br>Abund (%)",
                        low = "#FFFFFF", high = "#FF0000",
                        expand = c(0,0)) +
    #scale_y_continuous(expand=c(0,0)) + 
    labs(title=paste(str_to_title(Tlevel), " Relative Abundance of RABRs"),
         x=NULL,
         y=NULL) +
    facet_grid(~section, scale="free_x", space="free",  switch="x",
               labeller = labeller(section=pretty)) +
    #      y="Mean Relative Abundance (%)") +
    theme_classic() +
    theme(#axis.text.x = element_markdown(),
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_markdown(),
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0),
      legend.title = element_markdown(size=8),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      legend.text = element_text(size=7),
      legend.key.height = unit(10, "pt"),
      strip.background = element_blank(),
      strip.text = element_markdown())
    
  ggsave(paste("ITSpilot/ITSplots/ITS_", name, "_heat_map_", Tlevel, ".tiff", sep=""), width=9, height=4)

}

for (i in c("pr2_euk")) {
  for (j in c("phylum", "class", "order", "family", "genus")) {
    abund(i, j)
  }
}
# , "pr2_fungi", "silv_euk", "silv_fungi"
# abund("pr2_euk", "phylum")
# abund("pr2_fungi", "phylum")
# abund("silv_euk", "phylum")
# abund("silv_fungi", "phylum")
