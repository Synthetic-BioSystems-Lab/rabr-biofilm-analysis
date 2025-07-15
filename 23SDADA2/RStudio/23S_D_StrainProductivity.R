library(tidyverse) 
library(ggplot2)
library(ggtext)
library(readxl)
library(NLP)
library(ggpubr)

setwd("~/Miller Lab/Rscripts_PilotRABR/DADA2/23SDADA2visualization")

gather_metadata <- function(target, t2, t3, page) {
  read_excel("23S_D_metadata.xlsx", sheet=page) %>%
    rename_all(tolower) %>%
    mutate(sample = str_replace_all(sample, "-", "_")) %>%
    select(sample, target, t2, t3)
}

#import ASV counts
seqtabE <- readRDS("seqtabE.rds")
seqtabA <- readRDS("seqtabA.rds")
taxE <- readRDS("taxE.rds")
taxA <- readRDS("taxA.rds")

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

# import metadata
all_metadata <- gather_metadata("section", "dry_productivity_substratum", "date", 1) %>%
  rename(productivity = dry_productivity_substratum)

all_metadata <- dplyr::bind_rows(all_metadata, gather_metadata("section", "productivity", "date", 2))

all_metadata %>%
  mutate(section = factor(section,
                          levels=c("pilot", "81RABR")))

#remove empty rows, add relative abundance, pivot for readability 
mid_composite <- composite[!(composite$count=="0"),] %>% 
  group_by(sample_id) %>% 
  #relative abundance calculated here 
  mutate(rel_abund = count / sum(count)) %>% 
  mutate(sample_id = str_replace_all(sample_id, "-", "_")) %>%
  ungroup() 

col <- 10
i = 1
j = 5
for(i in 1:nrow(mid_composite)) {
  if(is.na(mid_composite[i, 4])) {
    
  } else{
    for(j in 5:10){
      if(is.na(mid_composite[i, j])){
        prior <- mid_composite[i, j-1]
        if(grepl('Unclassified', prior, fixed=TRUE)){
          mid_composite[i,j] <- mid_composite[i, j-1]
        } else{
          prior <- paste("Unclassified ",prior, sep="")
          mid_composite[i,j] <- prior
        }
      }
    }
  }
}

trimmed_composite <- mid_composite %>% 
  pivot_longer(cols = c("Kingdom", "Phylum", "Class",  
                        "Order", "Family", "Genus"),  
               names_to = "level", values_to = "taxon")

for(i in 1:nrow(trimmed_composite)) {
  if(is.na(trimmed_composite[i, 7])) {
    trimmed_composite[i, 7] <- "Unclassified"
  }
}

otu_rel_abund <- inner_join(trimmed_composite, all_metadata,  by=c('sample_id'='sample'))

#function
tax = c("Chlorella", "Leptolyngbya", "Nannochloropsis", "Synechococcus")
lvl = "Genus"
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
  write.csv(filtered_rel_abund,paste("23Scsvs/23S_StrainProd_",lvl,".csv", sep=""), row.names = FALSE)
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
  
  ggsave(paste("23S_D_plots/23SD_RAvsProd_", lvl, ".tiff", sep=""), width=7, height=5)
  
}

tax = c("Chlorella", "Leptolyngbya", "Nannochloropsis", "Synechococcus")
lvl = "Genus"
rabr = 'pilot'
x_spot = 3.2
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
  #write.csv(filtered_rel_abund,paste("23Scsvs/23S_StrainProd_",lvl,".csv", sep=""), row.names = FALSE)
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
  
  ggsave(paste("23S_D_plots/23SD_RAvsProd_", lvl, "_", rabr, ".tiff", sep=""), width=7, height=5)
  
}


RA_Prod("Genus", c("Chlorella", "Leptolyngbya", "Nannochloropsis", "Synechococcus"), 15, 'top')
RA_Prod("Family", c("Chlorellaceae", "Leptolyngbyaceae", "Microcoleaceae", "Monodopsidaceae", "Synechococcaceae"), 15, 'top')
RA_Prod("Order", c("Chlorellales", "Eustigmatales", "Oscillatoriales", "Prasiolales", "Synechococcales"), 15, 'top')
RA_Prod("Class", c("Chlorophyceae", "Cyanophyceae", "Eustigmatophyceae", "Trebouxiophyceae"), 15, 'top')
RA_Prod("Phylum", c("Chlorophyta", "Ochrophyta", "Cyanobacteria"), 15, 'center')

RA_Prod_RABR("Genus", c("Chlorella", "Leptolyngbya", "Nannochloropsis", "Synechococcus"), 'pilot', 1.5, 'top')
RA_Prod_RABR("Family", c("Chlorellaceae", "Leptolyngbyaceae", "Microcoleaceae", "Monodopsidaceae", "Synechococcaceae"), 'pilot', 1.2, 'top')
RA_Prod_RABR("Order", c("Chlorellales", "Eustigmatales", "Oscillatoriales", "Prasiolales", "Synechococcales"), 'pilot', 1, 'top')
RA_Prod_RABR("Class", c("Chlorophyceae", "Cyanophyceae", "Eustigmatophyceae", "Trebouxiophyceae"), 'pilot', 1.5, 'center')
RA_Prod_RABR("Phylum", c("Chlorophyta", "Ochrophyta", "Cyanobacteria"), 'pilot', 1, 'center')
RA_Prod_RABR("Genus", c("Chlorella", "Leptolyngbya", "Nannochloropsis", "Synechococcus"), '81RABR', 16, 'top')
RA_Prod_RABR("Family", c("Chlorellaceae", "Leptolyngbyaceae", "Microcoleaceae", "Monodopsidaceae", "Synechococcaceae"), '81RABR', 16, 'top')
RA_Prod_RABR("Order", c("Chlorellales", "Eustigmatales", "Oscillatoriales", "Prasiolales", "Synechococcales"), '81RABR', 16, 'top')
RA_Prod_RABR("Class", c("Chlorophyceae", "Cyanophyceae", "Eustigmatophyceae", "Trebouxiophyceae"), '81RABR', 16, 'top')
RA_Prod_RABR("Phylum", c("Chlorophyta", "Ochrophyta", "Cyanobacteria"), '81RABR', 16, 'center')

