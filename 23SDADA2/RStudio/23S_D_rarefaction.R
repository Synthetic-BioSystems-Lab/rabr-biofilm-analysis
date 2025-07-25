library(tidyverse) 
library(vegan)
library(ggtext)
library(readxl)

set.seed(1111)
setwd("~/Miller Lab/Rscripts_PilotRABR/DADA2/23SDADA2visualization")

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
composite <- rbind(compositeE, compositeA) %>%
  mutate(sample_id = str_replace_all(sample_id, "-", "_"))

# Rarefaction
shared <- rbind(asv_seqtabA, asv_seqtabE)%>%
  pivot_wider(names_from = ASVs, values_from = count, values_fill = 0) %>%
  mutate(sample_id = str_replace_all(sample_id, "-", "_"))
loc_shared_df <- shared %>%
  as.data.frame()

rownames(loc_shared_df) <- loc_shared_df$sample_id
loc_shared_df <- loc_shared_df[, -1]

# plot
rarecurve_data <- rarecurve(loc_shared_df, step=100)

map_dfr(rarecurve_data, bind_rows) %>%
  bind_cols(Group = rownames(loc_shared_df), .) %>%
  pivot_longer(-Group) %>%
  mutate(n_seqs = as.numeric(str_replace(name, "N", ""))) %>%
  select(-name) %>%
  ggplot(aes(x=n_seqs, y=value, group=Group)) +
  geom_line() +
  theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  labs(x="Number of Sequences", y="Number of OTUs", 
       title="Rarefaction Curves for All Samples") +
  theme(plot.title=element_text(hjust=0.5),
        legend.text = element_markdown(),
        legend.key.size = unit(10, "pt"),
        strip.background = element_blank(),
        strip.placement="outside",
        strip.text.x = element_markdown())

ggsave("23S_D_plots/23SD_rarefaction.tiff", width=5, height=4)

#test grouping
gather_metadata <- function(target, t2, t3, page) {
  read_excel("23S_D_metadata.xlsx", sheet=page) %>%
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

prep_rarecurve <- map_dfr(rarecurve_data, bind_rows) %>%
  bind_cols(Group = rownames(loc_shared_df), .) %>%
  inner_join(all_metadata, by=c("Group"="sample")) %>%
  pivot_longer(-c(Group, section, date, label)) %>%
  mutate(n_seqs = as.numeric(str_replace(name, "N", ""))) %>%
  select(-name)

Pilot_color <- "#BEBEBE"
CVWRF_color <- "#0000FF"
Labrabr_color <- "#FF0000"
TF_color <- "#00FF00"
GHR_color <- "cyan"
Control_color <- "#BE00FF"

write.csv(prep_rarecurve,"23Scsvs/23S_rare_allcolor.csv", row.names = FALSE)

breaks <- c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control")
labels <- c("Pilot RABR", "CVWRF", "Lab-scale RABRs", "Trickling Filter", "GHR", "Control")
prep_rarecurve %>%
  ggplot(aes(x=n_seqs, y=value, group=Group, colour = section)) +
  geom_line() + geom_point(size=1) +
  theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  labs(x="Number of Sequences", y="Number of ASVs", 
       title="Rarefaction Curves for All Samples") +
  theme(plot.title=element_text(hjust=0.5),
        legend.text = element_markdown(),
        legend.key.size = unit(10, "pt"),
        strip.background = element_blank(),
        strip.placement="outside",
        strip.text.x = element_markdown())+
  scale_color_manual(name=NULL, breaks=breaks,
                     labels=labels,
                     values=c(Pilot_color, CVWRF_color, Labrabr_color, TF_color, GHR_color, Control_color)) +
  scale_fill_manual(name=NULL,
                    breaks=breaks,
                    labels=labels,
                    values=c(Pilot_color, CVWRF_color, Labrabr_color, TF_color, GHR_color, Control_color))

ggsave("23S_D_plots/23SD_rarefaction_color.tiff", width=7, height=4) 

# pilot vs labRABR vs TF vs GH vs 
# Rarefaction

# Lab RABRs
# "23S_R27_11_3_21", "17_R36_10_18_22_23S", "23S_R43_11_15_21",
# "11_R45_10_18_21_23S", "16_R45_11_15_22_23S", "23S_R46_11_5_21",
# "23S_R58_10_28_21", "18_R60_11_1_21_23S", "12_R60_11_15_22_23S",
# "23S_R7_11_15_21", "21_R72_11_15_21_23S", "22_R75_11_15_21_23S"

# Pilot
# "10_5_23S", "19_23S", "26_23S", "11S_23S", "11R_23S", "S1_23S", "S2_23S", "S3_23S"

# TF
# "3_TF_5_25_22_23S", "23S_TF_6_9_22", "10_TF_6_22_22_23S", "TF_7_6_21", "TF_9_11_21", "TF_11_9_21_R1"

# CVWRF
# "1_CVWRF_PR_6_22_22_23S", "4_CVWRF_PSR_6_22_22_23S"

# GH
# "23S_GHR_6_15_22", "6_GHR_5_1_22_23S"

# Control
# "C1_23S", "C2_23S"

rare <- function(section, samples) {
  loc_shared_df <- shared %>%
    filter(sample_id %in% samples) %>%
    as.data.frame()
  
  rownames(loc_shared_df) <- loc_shared_df$Group
  loc_shared_df <- loc_shared_df[, -1]
  
  # plot
  rarecurve_data <- rarecurve(loc_shared_df, step=100)
  
  DF <- map_dfr(rarecurve_data, bind_rows) %>%
    bind_cols(Group = rownames(loc_shared_df), .)  %>%
    pivot_longer(-Group) %>%
    mutate(n_seqs = as.numeric(str_replace(name, "N", ""))) %>%
    select(-name)
  
  ggplot(data = DF, aes(x=n_seqs, y=value, group=Group)) +
    geom_line() +
    #geom_text(aes(label=Group),
              #data = DF %>% filter(n_seqs == median(n_seqs)),
              #nudge_x = 0.35, 
              #size = 4) +
    theme_classic() +
    scale_y_continuous(expand=c(0,0)) +
    labs(x="Number of Sequences", y="Number of OTUs", 
         title=paste("Rarefaction Curves for ", section, sep="")) +
    theme(plot.title=element_text(hjust=0.5),
          legend.text = element_markdown(),
          legend.key.size = unit(10, "pt"),
          strip.background = element_blank(),
          strip.placement="outside",
          strip.text.x = element_markdown())
  
  ggsave(paste("23S_D_plots/23SD_rarefaction_", tolower(section), ".tiff", sep=""), width=5, height=4)
  
}

rare("Lab RABRs", c("R27_11_3_21_23S", "R36_10_18_21_23S", "R43_11_15_21_23S", 
                    "R45_10_18_21_23S", "R45_11_15_21_23S", "R46_11_5_21_23S", 
                    "R58_10_28_21_23S", "R60_11_1_21_23S", "R60_11_15_21_23S", 
                    "R7_11_15_21_23S", "R72_11_15_21_23S", "R75_11_15_21_23S"))
rare("Pilot", c("10_5_23S", "19_23S", "26_23S", "11S_23S", "11R_23S", "S1_23S", "S2_23S", "S3_23S"))
rare("TF", c("TF_5_25_22_23S", "TF_6_9_22_23S", "TF_6_22_22_23S", "TF_7_6_21_23S", "TF_9_11_21_23S", "TF_11_9_21_23S"))
rare("CVWRF", c("CVWRF_PR_6_9_22_23S", "CVWRF_PSR_6_22_22_23S"))
rare("GH", c("GHR_6_15_22_23S", "GHR_5_1_22_23S"))
rare("Control", c("C1_23S", "C2_23S"))
