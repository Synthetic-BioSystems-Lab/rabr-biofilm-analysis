library(tidyverse) 
library(ggplot2)
library(ggtext)
library(readxl)


#relative abund stacked bar graph
#relative abund heat map
#relative abund bar graph

setwd("~/Miller Lab/Rscripts_PilotRABR/DADA2/23SDADA2visualization")

# final_shared <- "final.opti_mcc.shared"
# final_tax <- "final.opti_mcc.0.03.cons.taxonomy"

gather_metadata <- function(target, t2, t3, t4, page) {
  read_excel("23S_D_metadata.xlsx", sheet=page) %>%
    rename_all(tolower) %>%
    mutate(sample = str_replace_all(sample, "-", "_")) %>%
    select(sample, target, t2, t3, t4)
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
all_metadata <- gather_metadata("section", "date", "label","order", 1) 

for(i in 2:6) {
  all_metadata <- dplyr::bind_rows(all_metadata, gather_metadata("section", "date", "label","order", i))
}

all_metadata %>%
  mutate(section = factor(section,
                          levels=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control")))

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

# sample_order <- all_metadata %>%
#   arrange(date) %>%
#   mutate(order = 1:nrow(.)) %>%
#   select(sample, order, date)

tax <- "Order"
poollvl <- 5
abund <- function(tax, poollvl) {
  #"Kingdom", "Phylum", "Class", "Order", "Family", "Genus"
  taxon_rel_abund <- otu_rel_abund %>%
    filter(level==tax) %>%
    group_by(section, sample_id, taxon, date, label, order) %>%
    summarize(rel_abund = 100*sum(rel_abund), .groups = "drop") %>%
    # group_by(sample_id, taxon, section) %>%
    # summarize(mean_rel_abund = 100*mean(rel_abund), .groups = "drop") %>%
    mutate(taxon = str_replace(taxon, "(.*)_unclassified", "Unclassified *\\1*"),
           taxon = str_replace(taxon, "^(\\S*)$", "*\\1*"))
    #filter(sample_id != "C1_23S" & sample_id != "C2_23S")
  
  taxon_rel_abund <- taxon_rel_abund %>%
    mutate(section = str_replace_all(section, "control", "1_control")) %>%
    mutate(section = str_replace_all(section, "81RABR", "2_81RABR")) %>%
    mutate(section = str_replace_all(section, "GHR", "3_GHR")) %>%
    mutate(section = str_replace_all(section, "pilot", "4_pilot")) %>%
    mutate(section = str_replace_all(section, "CVWRF", "5_CVWRF")) %>%
    mutate(section = str_replace_all(section, "TF", "6_TF"))
  
  #set aside other
  # taxon_pool <- taxon_rel_abund %>%
  #   group_by(section, taxon, rel_abund) %>%
  #   summarize(mean=mean(rel_abund), .groups="drop") %>%
  #   group_by(taxon) %>%
  #   summarize(pool = max(mean) < poollvl, 
  #             mean = mean(mean),
  #             .groups="drop")
  
  # RvsS <- taxon_rel_abund %>%
  #   filter(sample_id %in% c("11S_23S", "11R_23S")) %>%
  #   group_by(sample_id) %>%
  #   pivot_wider(names_from = sample_id, values_from = mean_rel_abund) %>%
  #   data.table::setnames("11S_23S",'Srel_abund') %>%
  #   data.table::setnames("11R_23S",'Rrel_abund') %>%
  #   mutate(diff = (Srel_abund - Rrel_abund) / Rrel_abund) %>%
  #   ungroup() 
  
  pretty <- c("4_pilot" = "Pilot-scale RABRs",
              "2_81RABR" = "Lab-scale RABRs",
              "1_control" = "Control",
              "5_CVWRF" = "CVWRF",
              "3_GHR" = "GHR",
              "6_TF" = "TF")
  #<br>
  write.csv(taxon_rel_abund,paste("23Scsvs/23S_abund_", tolower(tax), ".csv", sep=""), row.names = FALSE)
  
  
  #assemble others and make RA stacked bar plot
  prep <- taxon_rel_abund %>% # inner_join(taxon_rel_abund, taxon_pool, by="taxon") %>%
    #filter(sample_id %in% c("S1_23S", "S2_23S", "S3_23S")) %>%
    mutate(taxon = if_else(rel_abund < poollvl, "Other", taxon)) %>%
    group_by(sample_id, label, section, taxon, order) %>%
    summarize(rel_abund = sum(rel_abund), .groups = "drop") %>%
    # mutate(taxon = factor(taxon),
    #        taxon = fct_reorder(taxon, mean, .desc=TRUE),
    #        taxon = fct_shift(taxon, n=1)) %>%
    #inner_join(., sample_order, by=c('sample_id'='sample')) %>%
    mutate(label = factor(label),
           label= fct_reorder(label, order)) %>% #%>% select !other
    filter(taxon != "Other")
  
  prep %>%
    ggplot(aes(x = label, y = rel_abund, fill = taxon)) +
    geom_col() + 
    #scale_fill_manual(name = NULL, values = c(brewer.pal(6, "Dark2"), "gray")) +
    scale_fill_discrete(name=NULL) +
    scale_y_continuous(expand=c(0,0)) +
    facet_grid(~section, scale="free_x", space="free", 
               labeller = labeller(section=pretty)) +
    labs(title=paste(tax, " Relative Abundance of RABRs", sep=""),
         x = NULL,
         y = "Mean Relative Abundance (%)") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.text = element_markdown(), 
          axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0),
          legend.key.size = unit(10, "pt"),
          strip.background = element_blank(),
          strip.text = element_markdown()) +
    ylim(0, 100)
  
  ggsave(paste("23S_D_plots/23SD_stacked_bar_", tolower(tax),".tiff", sep=""), width=13, height=4)
  
  # Relative Abundance
  
  # Heat Map
  # inner_join(taxon_rel_abund, taxon_pool, by="taxon") %>%
  #   mutate(taxon = if_else(pool, "Other", taxon)) %>%
  #   group_by(sample_id, section, taxon) %>%
  #   summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  #   inner_join(., sample_order, by="sample_id") %>%
  #   mutate(sample_id = factor(sample_id),
  #          sample_id= fct_reorder(sample_id, order)) %>%
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
    labs(title=paste(tax, " Relative Abundance of RABRs", sep=""),
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
  #coord_fixed(ratio = 4)
  
  ggsave(paste("23S_D_plots/23SD_heat_map_", tolower(tax), ".tiff", sep=""), width=13, height=4)
}

abund("Phylum", 5)
abund("Class", 5)
abund("Order", 5)
abund("Family", 5)
abund("Genus", 5)

# abund("Phylum", 5)
# abund("Class", 5)
# abund("Order", 5)
# abund("Family", 10)
# abund("Genus", 10)
