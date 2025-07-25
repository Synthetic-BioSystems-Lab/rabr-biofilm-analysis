library(tidyverse) 
library(vegan)
library(readxl)
library(ggtext)

set.seed(1111)
setwd("~/Miller Lab/Rscripts_PilotRABR")
name <- "pr2_euk"
loc <- "Lab_RABRs"
string <- Lab_RABRs


rarefaction_all <- function(name) {
  # Rarefaction
  loc_shared <- read_tsv(paste("ITSfinal2/final_",name,".agc.shared", sep="")) %>%
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
  ggsave(paste("ITSfinal2/ITSplots/ITS_", name, "_rarefaction.tiff", sep=""), width=5, height=4)

  #test grouping
  gather_metadata <- function(target, t2, t3, page) {
    read_excel("ITSfinal2/ITS_metadata.xlsx", sheet=page) %>%
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
  
  write.csv(prep_rarecurve,paste("ITSfinal2/ITScsvs/IT_", name, "_rare_allcolor.csv", sep=""), row.names = FALSE)
  
  breaks <- c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control")
  labels <- c("Pilot RABR", "CVWRF", "Lab-scale RABRs", "Trickling Filter", "GHR", "Control")
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
  ggsave(paste("ITSfinal2/ITSplots/ITS_", name, "_rarefaction_color.tiff", sep=""), width=7, height=4) 
  
}




rarefaction <- function(name, loc, string) {
  # pilot vs labRABR vs TF vs GH vs 
  # Rarefaction
  loc_shared <- read_tsv(paste("ITSfinal2/final_",name,".agc.shared", sep="")) %>%
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
    filter(Group %in% string)
  
  
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
         title=paste("Rarefaction Curves for ", loc, " Samples", sep="")) +
    theme(plot.title=element_text(hjust=0.5),
          legend.text = element_markdown(),
          legend.key.size = unit(10, "pt"),
          strip.background = element_blank(),
          strip.placement="outside",
          strip.text.x = element_markdown())
  
  ggsave(paste("ITSfinal2/ITSplots/ITS_", name, "_rarefaction_", loc, ".tiff", sep=""), width=5, height=4)
}

#rarefaction <- function(name, loc, string)

Lab_RABRs <- c("R27_11_3_21_ITS", "R36_10_28_21_ITS", "R43_11_15_21_ITS", 
"R45_10_18_21_ITS", "R45_11_15_21_ITS", "R46_11_5_21_ITS", "R58_10_28_21_ITS", 
"R60_11_1_21_ITS", "R60_11_15_21_ITS", "R7_11_15_21_ITS", "R72_11_15_21_ITS", "R75_11_15_21_ITS")
Pilot_RABRs <- c("10_5_ITS", "19_ITS", "26_ITS", "11S_ITS", "11R_ITS", "S1_ITS", "S2_ITS", "S3_ITS")
Trickling_Filter <- c("TF_5_25_22_ITS", "TF_6_9_22_ITS", "TF_6_22_22_ITS", "TF_7_6_21_ITS", "TF_9_11_21_ITS", "TF_11_9_21_ITS")
CVWRF <- c("CVWRF_PR_6_9_22_ITS", "CVWRF_PSR_6_22_22_ITS")
GH <- c("GHR_6_15_22_ITS", "GHR_5_1_22_ITS")
Control <- c("C1_ITS", "C2_ITS")

for (i in c("pr2_euk")) {
  rarefaction_all(i)
  rarefaction(i, "Lab_RABRs", Lab_RABRs)
  rarefaction(i, "Pilot_RABRs", Pilot_RABRs)
  rarefaction(i, "Trickling_Filter", Trickling_Filter)
  rarefaction(i, "CVWRF", CVWRF)
  rarefaction(i, "GH", GH)
  rarefaction(i, "Control", Control)
}
