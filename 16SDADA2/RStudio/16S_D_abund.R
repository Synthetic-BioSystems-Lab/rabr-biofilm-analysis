library(tidyverse) 
library(ggplot2)
library(ggtext)
library(readxl)


#relative abund stacked bar graph
#relative abund heat map
#relative abund bar graph

setwd("~/Miller Lab/Rscripts_PilotRABR/DADA2/16SDADA2visualization")

# final_shared <- "final.opti_mcc.shared"
# final_tax <- "final.opti_mcc.0.03.cons.taxonomy"

gather_metadata <- function(target, t2, t3, t4, page) {
  read_excel("16S_D_metadata.xlsx", sheet=page) %>%
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
all_metadata <- gather_metadata("section", "date", "label", "order", 1) 

for(i in 2:6) {
  all_metadata <- dplyr::bind_rows(all_metadata, gather_metadata("section", "date", "label", "order", i))
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
    for(j in 5:9){
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
  if(is.na(trimmed_composite[i, 6])) {
    trimmed_composite[i, 6] <- "Unclassified"
  }
}

otu_rel_abund <- inner_join(trimmed_composite, all_metadata,  by=c('sample_id'='sample'))

lvl <- "Phylum"
pl <- 5
#color <- genusColor
abund <- function(lvl, pl, color){
  #"Kingdom", "Phylum", "Class", "Order", "Family", "Genus"
  taxon_rel_abund <- otu_rel_abund %>%
    filter(level==lvl) %>%
    group_by(section, sample_id, taxon, date, label, order) %>%
    summarize(rel_abund = 100*sum(rel_abund), .groups = "drop") %>%
    # group_by(sample_id, taxon, section) %>%
    # summarize(mean_rel_abund = 100*mean(rel_abund), .groups = "drop") %>%
    mutate(taxon = str_replace(taxon, "(.*)_unclassified", "Unclassified *\\1*"),
           taxon = str_replace(taxon, "^(\\S*)$", "*\\1*"))
    #filter(sample_id != "C1_16S" & sample_id != "C2_16S")
  
  taxon_rel_abund <- taxon_rel_abund %>%
    mutate(section = str_replace_all(section, "control", "1_control")) %>%
    mutate(section = str_replace_all(section, "81RABR", "2_81RABR")) %>%
    mutate(section = str_replace_all(section, "GHR", "3_GHR")) %>%
    mutate(section = str_replace_all(section, "pilot", "4_pilot")) %>%
    mutate(section = str_replace_all(section, "CVWRF", "5_CVWRF")) %>%
    mutate(section = str_replace_all(section, "TF", "6_TF"))
  
  write.csv(taxon_rel_abund, paste("16SDcsvs/16SD_abund_", tolower(lvl), ".csv", sep=""), row.names = FALSE)
  
  #set aside other
  # taxon_pool <- taxon_rel_abund %>%
  #   group_by(section, taxon, rel_abund) %>%
  #   summarize(mean=mean(rel_abund), .groups="drop") %>%
  #   group_by(taxon) %>%
  #   summarize(pool = max(mean) < pl, 
  #             mean = mean(mean),
  #             .groups="drop")
  
  # sample_order <- all_metadata %>%
  #   arrange(date) %>%
  #   mutate(order = 1:nrow(.)) %>%
  #   select(sample, order, date)
  
  # RvsS <- taxon_rel_abund %>%
  #   filter(sample_id %in% c("11S_16S", "11R_16S")) %>%
  #   group_by(sample_id) %>%
  #   pivot_wider(names_from = sample_id, values_from = mean_rel_abund) %>%
  #   data.table::setnames("11S_16S",'Srel_abund') %>%
  #   data.table::setnames("11R_16S",'Rrel_abund') %>%
  #   mutate(diff = (Srel_abund - Rrel_abund) / Rrel_abund) %>%
  #   ungroup() 
  
  pretty <- c("4_pilot" = "Pilot-scale RABRs",
              "2_81RABR" = "Lab-scale RABRs",
              "1_control" = "Control",
              "5_CVWRF" = "CVWRF",
              "3_GHR" = "GHR",
              "6_TF" = "TF")
  #<br>
  
  #assemble others and make RA stacked bar plot
  prep <- taxon_rel_abund %>% # inner_join(taxon_rel_abund, taxon_pool, by="taxon") %>%
    #filter(sample_id %in% c("S1_16S", "S2_16S", "S3_16S")) %>%
    mutate(taxon = if_else(rel_abund < pl, "Other", taxon)) %>%
    group_by(sample_id, label, section, taxon, order) %>%
    summarize(rel_abund = sum(rel_abund), .groups = "drop") %>%
    # mutate(taxon = factor(taxon),
    #        taxon = fct_reorder(taxon, mean, .desc=TRUE),
    #        taxon = fct_shift(taxon, n=1)) %>%
    #inner_join(., sample_order, by=c("sample_id"="sample")) %>%
    mutate(label = factor(label),
           label= fct_reorder(label, order)) %>% #%>% select !other
    filter(taxon != "Other")
  
  write.csv(prep,paste("16S_D_plots/16SD_abund_", tolower(lvl), ".csv", sep=""), row.names = FALSE)
  
  prep %>%
    ggplot(aes(x = label, y = rel_abund, fill = taxon)) +
    geom_col() + 
    scale_fill_manual(name = NULL, values = color) +
    #scale_fill_discrete(name=NULL) +
    scale_y_continuous(expand=c(0,0)) +
    facet_grid(~section, scale="free_x", space="free", 
               labeller = labeller(section=pretty)) +
    labs(title=paste(lvl, " Relative Abundance of RABRs", sep=""),
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
  
  ggsave(paste("16S_D_plots/16SD_stacked_bar_", tolower(lvl), ".tiff", sep=""), width=13, height=4)
  
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
    labs(title=paste(lvl, " Relative Abundance of RABRs", sep=""),
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
  
  ggsave(paste("16S_D_plots/16SD_heat_map_", tolower(lvl), ".tiff", sep=""), width=13, height=4)
}

genusColor <- c("*Chryseobacterium*"="#FA5B3D", "*Clostridium_sensu_stricto_13*"="#474440", 
                "*Comamonas*"="#C04000", "*Pseudomonas*"="#00008B", "Unclassified *Burkholderiaceae*", "#8B008B", 
                "Unclassified *Xanthomonadaceae*" = "#808000", "*Subsaxibacter*"= "#C86558",
                "*uncultured*"="#0000CD", "Unclassified *Chitinophagaceae*"	="#AB3DA9",
                "*Algoriphagus*"="#2E8B57", "*Symphothece_PCC-7002*"=	"#008000", "Unclassified *Bacteria*"="#465DAA", 
                "Unclassified *Cyclobacteriaceae*"="#663399", "*Thauera*"="#3B3734", "*Acinetobacter*"="#6B506B", 
                "*Fibrobacter*"="#DE25DA", "*Macellibacteroides*"="#0000B3", "*Psychrobacter*"="#000080", 
                "*Succinispira*"="#A0522D", "Unclassified *Bacteroidia*"="#FF1493", "*Devosia*"="#0010D9", 
                "*Hydrogenophaga*"="#D2691E", "*Rhodanobacter*"="#A020F0", "Unclassified *Rhodanobacteraceae*"="#8B8D8B", 
                "Unclassified *Rhodocyclaceae*"="#4682B4", "*Pseudoxanthomonas*"="#0000FF", 
                "*Flavobacterium*"="#C0C0C0", "*Simplicispira*"="#B8860B", "*Giesbergeria*"="#2E2B28", 
                "*Bacteroidetes_vadinHA17_ge*"="#FF80FF", "*Lacihabitans*"="#00FFFF", "*Bergeyella*"="#EB44E8", 
                "*Nitrospira*"="#DF8879", "*Pleurocapsa_PCC-7319*"="#FF00FF", "Unclassified *Weeksellaceae*"="#B03060", 
                "*Arenimonas*"="#C71585", "Unclassified *Enterobacteriaceae*"="#483D8B", "*Aeromonas*"="#54504C", 
                "Unclassified Eukaryota"="#696969", "*Brevundimonas*"="#800000", "Unclassified Chitinophagaceae"="#008080", 
                "*Phormidesmis_ANT.LACV5.1*"="#FFFF00", "*Thermomonas*"="#B04238", "*Synechococcus_PCC-7942*"="#8B0000", 
                "Unclassified Bacteria"="#2F4F4F", "Unclassified Chloroplast"="#9370DB", "*Arcobacter*"="#808080", 
                "Unclassified Burkholderiaceae"="#9932CC", "*Oscillochloris*"="#CBD6E4", 
                "*Tychonema_CCAP_1459-11B*"="#DD5FB3", "*Acidovorax*"="#A4A2A8", "*Limnohabitans*"="#4169E1", 
                "*Caenimonas*"="#FF0000", "*Candidatus_Nitrosopumilus*"="#B3BFD1", "*Clade_Ia*"="#D7E1EE", 
                "*Synechococcus_CC9902*"="#DC143C", "Unclassified SAR86_clade"="#BFCDDB", "Unclassified Rhodocyclaceae"="#991F17")

familyColor <- c("*Burkholderiaceae*"="#FA5B3D", "*Chitinophagaceae*"="#474440", "*Clostridiaceae_1*"="#C04000", 
                 "*Flavobacteriaceae*"="#00008B", "*Pseudomonadaceae*"="#8B008B", "*Weeksellaceae*"="#808000", 
                 "*Xanthomonadaceae*"="#C86558", "*Cyanobacteriaceae*"="#0000CD", "*Cyclobacteriaceae*"="#AB3DA9", 
                 "Unclassified *Bacteria*"="#2E8B57", "*Rhodocyclaceae*"="#008000", "*Acidaminococcaceae*"="#465DAA", 
                 "*Fibrobacteraceae*"="#663399", "*Moraxellaceae*"="#3B3734", "*Rikenellaceae*"="#6B506B", 
                 "*Tannerellaceae*"="#DE25DA", "Unclassified *Bacteroidia*"="#0000B3", "*Devosiaceae*"="#000080", 
                 "*Rhodanobacteraceae*"="#A0522D", "*Rhodobacteraceae*"="#FF1493", "*Bacteroidetes_vadinHA17*"="#0010D9", 
                 "*Spirosomaceae*"="#D2691E", "*Nitrospiraceae*"="#A020F0", "*Xenococcaceae*"="#8B8D8B", "*Enterobacteriaceae*"="#4682B4", 
                 "*Saprospiraceae*"="#0000FF", "*Aeromonadaceae*"="#C0C0C0", "Unclassified Eukaryota"="#B8860B", 
                 "*Caulobacteraceae*"="#2E2B28", "*Phormidesmiaceae*"="#FF80FF", "*Synechococcaceae*"="#00FFFF", 
                 "Unclassified Bacteria"="#EB44E8", "Unclassified Chloroplast"="#DF8879", "*Arcobacteraceae*"="#FF00FF", 
                 "*Sphingomonadaceae*"="#B03060", "*Chloroflexaceae*"="#C71585", "*Phormidiaceae*"="#483D8B", "*Clade_I*"="#54504C", 
                 "*Cyanobiaceae*"="#696969", "*Nitrosopumilaceae*"="#800000", "Unclassified SAR86_clade"="#008080")

orderColor <- c("*Bacteroidales*"="#FA5B3D", "*Betaproteobacteriales*"="#474440", "*Chitinophagales*"="#C04000", 
                "*Clostridiales*"="#00008B", "*Flavobacteriales*"="#8B008B", "*Pseudomonadales*"="#808000", "*Xanthomonadales*"="#C86558", 
                "*Cytophagales*"="#0000CD", "*Nostocales*"="#AB3DA9", "Unclassified *Bacteria*"="#2E8B57", "*Desulfovibrionales*"="#008000", 
                "*Fibrobacterales*"="#465DAA", "*Selenomonadales*"="#663399", "Unclassified *Bacteroidia*"="#3B3734", 
                "*Rhizobiales*"="#6B506B", "*Rhodobacterales*"="#DE25DA", "*Nitrospirales*"="#0000B3", 
                "*Enterobacteriales*"="#000080", "*Aeromonadales*"="#A0522D", "Unclassified Eukaryota"="#FF1493", 
                "*Caulobacterales*"="#0010D9", "*Phormidesmiales*"="#D2691E", "*Chloroplast*"="#A020F0", "*Synechococcales*"="#8B8D8B", 
                "Unclassified Bacteria"="#4682B4", "*Campylobacterales*"="#0000FF", "*Sphingomonadales*"="#C0C0C0", 
                "*Chloroflexales*"="#B8860B", "*Nitrosopumilales*"="#2E2B28", "*SAR11_clade*"="#FF80FF", "*SAR86_clade*"="#00FFFF")

classColor <- c("*Bacteroidia*"="#FA5B3D", "*Clostridia*"="#474440", "*Gammaproteobacteria*"="#C04000", 
                "*Oxyphotobacteria*"="#00008B", "Unclassified *Bacteria*"="#8B008B", "*Deltaproteobacteria*"="#808000", 
                "*Fibrobacteria*"="#C86558", "*Negativicutes*"="#0000CD", "*Alphaproteobacteria*"="#AB3DA9", "*Nitrospira*"="#2E8B57", 
                "Unclassified Eukaryota"="#008000", "Unclassified Bacteria"="#465DAA", "*Campylobacteria*"="#663399", 
                "*Chloroflexia*"="#3B3734", "*Nitrososphaeria*"="#6B506B")

phylumColor <- c("*Bacteroidetes*"="#FA5B3D", "*Firmicutes*"="#474440", "*Proteobacteria*"="#C04000", 
                 "*Cyanobacteria*"="#00008B", "Unclassified *Bacteria*"="#8B008B", "*Fibrobacteres*"="#808000", 
                 "*Nitrospirae*"="#C86558", "Unclassified Eukaryota"="#0000CD", "Unclassified Bacteria"="#AB3DA9", 
                 "*Epsilonbacteraeota*"="#2E8B57", "*Chloroflexi*"="#008000", "*Actinobacteria*"="#465DAA", 
                 "*Thaumarchaeota*"="#663399")

abund("Phylum", 5, phylumColor)
abund("Class", 5, classColor)
abund("Order", 5, orderColor)
abund("Family", 5, familyColor)
abund("Genus", 5, genusColor)

# abund("Phylum", 10)
# abund("Class", 15)
# abund("Order", 15)
# abund("Family", 15)
# abund("Genus", 17)

