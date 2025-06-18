#libraries
library(tidyverse) 
library(ggplot2)
library(ggpubr)
library(ggtext)
library(NLP)
library(readxl)

#prep (take from abund files)
setwd("~/Miller Lab/Rscripts_PilotRABR")

gather_metadata <- function(target, t2, t3, page) {
  read_excel("16Spilot/16S_metadata.xlsx", sheet=page) %>%
    rename_all(tolower) %>%
    mutate(sample = str_replace_all(sample, "-", "_")) %>%
    select(sample, target, t2, t3)
}

# import metadata
all_metadata <- gather_metadata("section", "date", "label", 1) 

for(i in 2:6) {
  all_metadata <- dplyr::bind_rows(all_metadata, gather_metadata("section", "date", "label", i))
}

all_metadata %>%
  mutate(section = factor(section,
                          levels=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control")))

pretty <- c("pilot" = "Pilot-scale RABRs",
            "81RABR" = "Lab-scale RABRs",
            "control" = "Control",
            "CVWRF" = "CVWRF",
            "GHR" = "GHR",
            "TF" = "Trickling Filter")

#mothur prep
final_shared <- "16Spilot/final.opti_mcc.shared"
final_tax <- "16Spilot/final.opti_mcc.0.03.cons.taxonomy"

#import otu counts
Motu_counts <- read_tsv(final_shared) %>%
  select(-label, -numOtus) %>%
  pivot_longer(-Group, names_to = "otu", values_to = "count") %>%
  rename(sample_id = Group)

# import taxonomy
Mtaxonomy <- read_tsv(final_tax) %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""), taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";")

#merge otu and taxon
Mcomposite <- inner_join(Motu_counts, Mtaxonomy, by="otu")

#remove empty rows, add relative abundance, pivot for readability
Mtrimmed_composite <- Mcomposite[!(Mcomposite$count=="0"),] %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), names_to = "level", values_to = "taxon")

Motu_rel_abund <- inner_join(Mtrimmed_composite, all_metadata,  by=c('sample_id'='sample'))

#DADA2 prep
#import ASV counts
seqtabE <- readRDS("DADA2/16SDADA2visualization/seqtabE.rds")
seqtabA <- readRDS("DADA2/16SDADA2visualization/seqtabA.rds")
taxE <- readRDS("DADA2/16SDADA2visualization/taxE.rds")
taxA <- readRDS("DADA2/16SDADA2visualization/taxA.rds")

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
DcompositeA <- inner_join(asv_seqtabA, asv_taxA, by="ASVs")
DcompositeE <- inner_join(asv_seqtabE, asv_taxE, by="ASVs")
Dcomposite <- rbind(DcompositeE, DcompositeA)

#remove empty rows, add relative abundance, pivot for readability 
Dtrimmed_composite <- Dcomposite[!(Dcomposite$count=="0"),] %>% 
  group_by(sample_id) %>% 
  #relative abundance calculated here 
  mutate(rel_abund = count / sum(count)) %>% 
  mutate(sample_id = str_replace_all(sample_id, "-", "_")) %>%
  ungroup() %>% 
  pivot_longer(cols = c("Kingdom", "Phylum", "Class",  
                        "Order", "Family", "Genus"),  
               names_to = "level", values_to = "taxon")

Dotu_rel_abund <- inner_join(Dtrimmed_composite, all_metadata,  by=c('sample_id'='sample'))




#command to make plots
tax = "Synechococcus_PCC-7942"
lvl = "Genus"
pool_lvl = 10
versus <- function(tax, lvl, pool_lvl) {
  Mtaxon_rel_abund <- Motu_rel_abund %>%
    filter(level==lvl) %>%
    group_by(section, sample_id, taxon, date, label) %>%
    summarize(rel_abund = 100*sum(rel_abund), .groups = "drop") %>%
    # group_by(sample_id, taxon, section) %>%
    # summarize(mean_rel_abund = 100*mean(rel_abund), .groups = "drop") %>%
    mutate(taxon = str_replace(taxon, "(.*)_unclassified", "Unclassified \\1"),
           taxon = str_replace(taxon, "^(\\S*)$", "\\1"),
           M_rel_abund = rel_abund) %>%
    select(section, sample_id, taxon, date, label, M_rel_abund)
  
  Dtaxon_rel_abund <- Dotu_rel_abund %>%
    filter(level==lvl) %>%
    group_by(section, sample_id, taxon, date, label) %>%
    summarize(rel_abund = 100*sum(rel_abund), .groups = "drop") %>%
    # group_by(sample_id, taxon, section) %>%
    # summarize(mean_rel_abund = 100*mean(rel_abund), .groups = "drop") %>%
    mutate(taxon = str_replace(taxon, "(.*)_unclassified", "Unclassified \\1"),
           taxon = str_replace(taxon, "^(\\S*)$", "\\1"),
           D_rel_abund = rel_abund) %>%
    select(section, sample_id, taxon, date, label, D_rel_abund)
  
  #merge by inner_join
  joined_rel_abund <- inner_join(Mtaxon_rel_abund, Dtaxon_rel_abund, by=c("section", "sample_id", "taxon", "date", "label"))
  if(tax==""){
    taxon_pool <- joined_rel_abund %>%
      #group_by(section, taxon, M_rel_abund, D_rel_abund) %>%
      #summarize(M_mean=mean(M_rel_abund), D_mean=mean(D_rel_abund), .groups="drop") %>%
      group_by(taxon) %>%
      summarize(pool = max(M_rel_abund) < pool_lvl|max(D_rel_abund) < pool_lvl, 
                #mean = mean(mean),
                .groups="drop")
    #plot?
    prep <- inner_join(joined_rel_abund, taxon_pool, by="taxon") %>%
      mutate(taxon = if_else(pool, "Other", taxon))
    prep <- prep[!(prep$taxon %in% "Other"),]
    
    
    reg<-lm(D_rel_abund ~ M_rel_abund, data=prep)                       
    #get intercept and slope value 
    coeff<-coefficients(reg)           
    intercept<-round(coeff[1], digits=3)
    slope<- round(coeff[2], digits=3) 
    write.csv(prep,paste("mothurvDADA2/vcsvs/16S_vs_", lvl, ".csv", sep=""), row.names = FALSE)
    
    prep %>%
      ggplot(aes(x=M_rel_abund, y=D_rel_abund)) +
      geom_point() + 
      labs(x="mothur Relative Abundance (%)", y="DADA2 Relative Abundance (%)") +
      ggtitle(paste(lvl, " Alignment", sep = ""), subtitle = paste("Slope=", slope, ", Intercept=", intercept, sep = "")) +
      theme_classic() +
      ylim(0, 100) +
      xlim(0, 100) +
      theme(plot.title=element_text(hjust=0.5),
            axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0)) +
      stat_cor()+
      geom_smooth(method=lm, se=FALSE)
    
    ggsave(paste("mothurvDADA2/vplots/16S_vs_", lvl, ".tiff", sep=""), width=4, height=4.05)
    
    prep %>%
      ggplot(aes(x=M_rel_abund, y=D_rel_abund, color = section)) +
      geom_point() + 
      labs(x="mothur Relative Abundance (%)", y="DADA2 Relative Abundance (%)") +
      ggtitle(paste(lvl, " Alignment", sep = "")) +
      theme_classic() +
      ylim(0, 100) +
      xlim(0, 100) +
      theme(plot.title=element_text(hjust=0.5),
            axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0)) +
      stat_cor()
    
    ggsave(paste("mothurvDADA2/vplots/16S_vs_", lvl, "_sectioned.tiff", sep=""), width=4, height=4.05)
    
  } else {
    # select for "strain"
    strain_rel_abund <- joined_rel_abund %>%
      filter(taxon==tax)
    
    # order by "strain" RA for mothur
    prep_rel_abund <- strain_rel_abund %>%
      arrange(M_rel_abund) %>%
      mutate(order = 1:nrow(.))%>%
      mutate(label = factor(label),
             label= fct_reorder(label, order),
             percent_diff = 100*abs(M_rel_abund - D_rel_abund) / ((M_rel_abund+D_rel_abund)/2)) %>%
      pivot_longer(c(M_rel_abund, D_rel_abund), names_to = "method", values_to = "rel_abund")
    
    write.csv(prep_rel_abund,paste("mothurvDADA2/vcsvs/16S_vs_", lvl, "_", tax, ".csv", sep=""), row.names = FALSE)
    
    # plot RA vs Samples with DADA2 and mothur separate colors
    prep_rel_abund %>%
      ggplot(aes(x = label, y = rel_abund, color = method)) +
      geom_point(size=3) + 
      #scale_fill_manual(name = NULL, values = c(brewer.pal(6, "Dark2"), "gray")) +
      scale_fill_discrete(name=NULL) +
      scale_y_continuous(expand=c(0,0)) +
      # facet_grid(~section, scale="free_x", space="free", 
      #            labeller = labeller(section=pretty)) +
      labs(title=paste("RA Comparison of ", tax, " levels using Mothur and DADA2", sep=""),
           x = NULL,
           y = "Mean Relative Abundance (%)") +
      theme_classic() +
      ylim(0, NA) +
      theme(plot.title = element_text(hjust = 0.5),
            legend.text = element_markdown(), 
            axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0),
            legend.key.size = unit(10, "pt"),
            strip.background = element_blank(),
            strip.text = element_markdown())
    
    ggsave(paste("mothurvDADA2/vplots/strain_vs/16S_vs_", lvl, "_", tax, ".tiff", sep=""), width=7, height=5)
    
    # calc RA diff? or % diff?
    prep_rel_abund %>%
      ggplot(aes(x = label, y = percent_diff, color = section)) +
      geom_point(size=3) + 
      #scale_fill_manual(name = NULL, values = c(brewer.pal(6, "Dark2"), "gray")) +
      scale_fill_discrete(name=NULL) +
      scale_y_continuous(expand=c(0,0)) +
      # facet_grid(~section, scale="free_x", space="free", 
      #            labeller = labeller(section=pretty)) +
      labs(title=paste("RA % Difference of ", tax, " levels using Mothur and DADA2", sep=""),
           x = NULL,
           y = "Relative Abundance Percent Difference") +
      theme_classic() +
      ylim(0, NA) +
      theme(plot.title = element_text(hjust = 0.5),
            legend.text = element_markdown(), 
            axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0),
            legend.key.size = unit(10, "pt"),
            strip.background = element_blank(),
            strip.text = element_markdown())
    
    ggsave(paste("mothurvDADA2/vplots/strain_vs/16S_vs_", lvl, "_", tax, "percdiff.tiff", sep=""), width=7, height=5)
    
  }
}

versus("", "Phylum", 10)
versus("", "Class", 10)
versus("", "Order", 10)
versus("", "Family", 10)
versus("", "Genus", 10)

versus("Phormidiaceae", "Family", 10)
versus("Tychonema_CCAP_1459-11B", "Genus", 10)
versus("Symphothece_PCC-7002", "Genus", 10)
versus("Nostocales", "Order", 10)
versus("Bacteroidetes", "Phylum", 10)
versus("Proteobacteria", "Phylum", 10)
versus("Cyanobacteria", "Phylum", 10)
versus("Alphaproteobacteria", "Class", 10)
versus("Bacteroidia", "Class", 10)
versus("Gammaproteobacteria", "Class", 10)
versus("Oxyphotobacteria", "Class", 10)
versus("Betaproteobacteriales", "Order", 10)
versus("Flavobacteriales", "Order", 10)
versus("Synechococcales", "Order", 10)
versus("Burkholderiaceae", "Family", 10)
versus("Weeksellaceae", "Family", 10)
versus("Synechococcaceae", "Family", 10)
versus("Synechococcus_PCC-7942", "Genus", 10)
versus("Symphothece_PCC-7002", "Genus", 10)

