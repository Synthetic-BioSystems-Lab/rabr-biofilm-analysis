library(tidyverse) 
library(vegan)
library(ggtext)
library(readxl)

set.seed(1111)
setwd("~/Miller Lab/Rscripts_PilotRABR")

# Rarefaction
loc_shared <- read_tsv("16Spilot/final.opti_mcc.0.03.subsample.shared") %>%
  select(Group, starts_with("Otu")) %>%
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  mutate(total=sum(value)) %>%
  arrange(total) %>%
  group_by(name) %>%
  mutate(total=sum(value)) %>% 
  filter(total != 0) %>%
  ungroup() %>%
  select(-total) 
  #filter(Group != "C1_16S" & Group != "C2_16S")

loc_shared_df <- loc_shared %>%
  pivot_wider(names_from="name", values_from="value", values_fill = 0) %>%
  as.data.frame()

rownames(loc_shared_df) <- loc_shared_df$Group
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

ggsave("16Spilot/16Splots/16S_rarefaction.tiff", width=5, height=4)

#test grouping
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

breaks <- c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control")
labels <- c("Pilot RABR", "CVWRF", "Lab-scale RABRs", "Trickling Filter", "GHR", "Control")

write.csv(prep_rarecurve,"16Spilot/16Scsvs/16S_rare_allcolor.csv", row.names = FALSE)

prep_rarecurve %>%
  ggplot(aes(x=n_seqs, y=value, group=Group, colour = section)) +
  geom_line() + geom_point(size=1) +
  theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  labs(x="Number of Sequences", y="Number of OTUs", 
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
ggsave("16Spilot/16Splots/16S_rarefaction_color.tiff", width=7, height=4) 

# pilot vs labRABR vs TF vs GH vs 
# Rarefaction
rare <-function(section, samples){
  loc_shared <- read_tsv("16Spilot/final.opti_mcc.0.03.subsample.shared") %>%
    select(Group, starts_with("Otu")) %>%
    pivot_longer(-Group) %>%
    group_by(Group) %>%
    mutate(total=sum(value)) %>%
    arrange(total) %>%
    group_by(name) %>%
    mutate(total=sum(value)) %>% 
    filter(total != 0) %>%
    ungroup() %>%
    select(-total) %>%
    filter(Group %in% samples)
  
  loc_shared_df <- loc_shared %>%
    pivot_wider(names_from="name", values_from="value", values_fill = 0) %>%
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
  
  ggsave(paste("16Spilot/16Splots/16S_rarefaction_", section, ".tiff", sep=""), width=5, height=4)
}
rare("labRABR", c("R27_11_3_21_16S", "R36_10_28_21_16S", "R43_11_15_21_16S", 
                  "R45_10_18_21_16S", "R45_11_15_21_16S", "R46_11_5_21_16S", 
                  "R58_10_28_21_16S", "R60_11_1_21_16S", "R60_11_15_21_16S", 
                  "R7_11_15_21_16S", "R72_11_15_21_16S", "R75_11_15_21_16S"))
rare("pilot", c("10_5_16S", "19_16S", "26_16S", "11S_16S", "11R_16S", "S1_16S", "S2_16S", "S3_16S"))
rare("TF", c("TF_5_25_22_16S", "TF_6_9_22_16S", "TF_6_22_22_16S", 
             "TF_7_6_21_16S", "TF_9_11_21_16S", "TF_11_9_21_16S"))
rare("CVWRF", c("CVWRF_PR_6_9_22_16S", "CVWRF_PSR_6_22_22_16S"))
rare("GHR", c("GHR_6_15_22_16S", "GHR_5_1_22_16S"))
rare("Control", c("C1_16S", "C2_16S"))

