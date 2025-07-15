library(vegan)
library(ggtext)
library(tidyverse)
library(readxl)

setwd("~/Miller Lab/Rscripts_PilotRABR")

gather_metadata <- function(target, t2, t3, page) {
  read_excel("ITSfinal2/ITS_metadata.xlsx", sheet=page) %>%
    rename_all(tolower) %>%
    mutate(sample = str_replace_all(sample, "-", "_")) %>%
    select(sample, target, t2, t3)
}

name <- "pr2_euk"

beta_div <- function(name) {
  #Beta Diversity 2
  all_shared <- read_tsv(paste("ITSfinal2/final_",name,".agc.shared", sep="")) %>%
    select(Group, starts_with("Otu")) %>%
    pivot_longer(-Group) %>%
    group_by(Group) %>%
    mutate(total=sum(value)) %>%
    # arrange(total) %>%
    # group_by(name) %>%
    # summarize(total=sum(value)) %>% filter(total ==0) %>%
    ungroup() %>%
    select(-total) %>%
    pivot_wider(Group) %>%
    filter(Group %in% c("10_5_ITS", "19_ITS", "26_ITS", "11S_ITS", "11R_ITS",
                        "S1_ITS", "S2_ITS", "S3_ITS"))
  
  # for(i in 1:nrow(all_shared)) {
  #   row <- all_shared[i,]
  #   row[1] <- gsub("_ITS", "", row[1])
  #   all_shared[i,1] <- row[1]
  # }
  
  
  all_shared <- as.data.frame(all_shared)
  
  rownames(all_shared) <- all_shared$Group
  all_shared <- all_shared[, -1]
  all_shared <- as.matrix(all_shared)
  
  
  all_dist <- avgdist(all_shared, dmethod="bray", iterations=1000, sample=2073)
  all_nmds <- metaMDS(all_dist) #low stress > insufficient data?
  
  all_scores <- scores(all_nmds) %>%
    as_tibble(rownames = "Group")
  
  all_metadata <- gather_metadata("section", "date", "label", 1)
  
  all_metadata_nmds <- inner_join(all_metadata, all_scores, by=c('sample'='Group')) 
  
  Five_color <- "#BEBEBE"
  Twelve_color <- "#0000FF"
  Nineteen_color <- "#FF0000"
  TwentySix_color <- "#00FF00"
  Two_color <- "cyan"
  
  date_centroid <- all_metadata_nmds %>%
    group_by(date) %>%
    summarize(NMDS1 = mean(NMDS1),
              NMDS2 = mean(NMDS2), .groups="drop")
  
  date_star <- all_metadata_nmds %>%
    group_by(date) %>%
    mutate(centroid1 = mean(NMDS1),
           centroid2 = mean(NMDS2)) %>%
    ungroup()
  
  breaks=c("2023_10_05", "2023_10_12", "2023_10_19", "2023_10_26", "2023_11_02")
  labels=c("2023_10_05", "2023_10_12", "2023_10_19", "2023_10_26", "2023_11_02")
  
  write.csv(date_star,paste("ITSfinal2/ITScsvs/IT_", name, "_beta_date.csv", sep=""), row.names = FALSE)
  
  ggplot(date_star, aes(x=NMDS1, xend=centroid1,
                        y=NMDS2, yend=centroid2, color=date)) +
    geom_segment() +
    geom_point() +
    geom_point(data=date_centroid,
               mapping=aes(x=NMDS1, fill=date, y=NMDS2),
               shape=22, size=5, show.legend=FALSE,
               #color="black",
               inherit.aes=FALSE) +
    #coord_fixed(xlim=c(-1, 1), ylim=c(1, 1)) +
    labs(x="NMDS Axis 1", y="NMDS Axis 2", title="All Samples NMDS by Date") +
    scale_color_manual(name=NULL, breaks=breaks,
                       labels=labels,
                       values=c(Five_color, Twelve_color, Nineteen_color, TwentySix_color, Two_color)) +
    scale_fill_manual(name=NULL, breaks=breaks,
                      labels=labels,
                      values=c(Five_color, Twelve_color, Nineteen_color, TwentySix_color, Two_color)) +
    theme_classic() + xlim(-0.4, 0.4) + ylim(-0.4, 0.4) +
    theme(legend.key.size = unit(0.25, "cm"),
          legend.position = c(0.85, 0.9),
          legend.background = element_rect(fill="NA",
                                           color="black"),
          legend.margin = margin(t=2, r=3, b=3, l=3),
          plot.title=element_text(hjust=0.5))
  
  ggsave(paste("ITSfinal2/ITSplots/ITS_", name, "_nmds_all_date.tiff", sep=""), width=5, height=4)
  
  #Beta Diversity PilotvLabRABR
  samples_wanted <- c("R27_11_3_21_ITS", "R36_10_28_21_ITS", "R43_11_15_21_ITS", 
                      "R45_10_18_21_ITS", "R45_11_15_21_ITS", "R46_11_5_21_ITS", "R58_10_28_21_ITS", 
                      "R60_11_1_21_ITS", "R60_11_15_21_ITS", "11_15_21_ITS", "R72_11_15_21_ITS", 
                      "R75_11_15_21_ITS", "10_5_ITS", 
                      "19_ITS", "26_ITS", "11S_ITS", "11R_ITS", "S1_ITS", 
                      "S2_ITS", "S3_ITS")
  
  
  all_shared <- read_tsv(paste("ITSfinal2/final_",name,".agc.shared", sep="")) %>%
    select(Group, starts_with("Otu")) %>%
    filter(Group %in% samples_wanted) %>%
    pivot_longer(-Group) %>%
    group_by(Group) %>%
    # mutate(total=sum(value)) %>%
    # arrange(total) %>%
    group_by(name) %>%
    mutate(total=sum(value)) %>%
    filter(total != 0) %>%
    # summarize(total=sum(value)) %>% filter(total ==0) %>%
    ungroup() %>%
    select(-total) %>%
    pivot_wider(Group) %>%
    as.data.frame
  
  #all_shared <- as.data.frame(all_shared)
  
  rownames(all_shared) <- all_shared$Group
  all_shared <- all_shared[, -1]
  all_shared <- as.matrix(all_shared)
  
  
  all_dist <- avgdist(all_shared, dmethod="bray", iterations=1000, sample=2073)
  all_nmds <- metaMDS(all_dist) #low stress > insufficient data?
  
  all_scores <- scores(all_nmds) %>%
    as_tibble(rownames = "Group")
  
  all_metadata <- gather_metadata("section", "date", "label", 1)
  
  for(i in 2:6) {
    all_metadata <- dplyr::bind_rows(all_metadata, gather_metadata("section", "date", "label", i))
  }
  
  all_metadata_nmds <- inner_join(all_metadata, all_scores, by=c('sample'='Group')) %>%
    mutate(section = factor(section,
                         levels=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control")))
  
  Pilot_color <- "#BEBEBE"
  CVWRF_color <- "#0000FF"
  Labrabr_color <- "#FF0000"
  TF_color <- "#00FF00"
  GHR_color <- "cyan"
  Control_color <- "#BE00FF"
  
  all_centroid <- all_metadata_nmds %>%
    group_by(section) %>%
    summarize(NMDS1 = mean(NMDS1),
              NMDS2 = mean(NMDS2), .groups="drop")
  
  all_star <- all_metadata_nmds %>%
    group_by(section) %>%
    mutate(centroid1 = mean(NMDS1),
           centroid2 = mean(NMDS2)) %>%
    ungroup()
  
  ggplot(all_star, aes(x=NMDS1, xend=centroid1,
                        y=NMDS2, yend=centroid2, color=section)) +
    geom_segment() +
    geom_point() +
    geom_point(data=all_centroid,
               mapping=aes(x=NMDS1, fill=section, y=NMDS2),
               shape=22, size=5, show.legend=FALSE,
               #color="black",
               inherit.aes=FALSE) +
    #coord_fixed(xlim=c(-1, 1), ylim=c(1, 1)) +
    labs(x="NMDS Axis 1", y="NMDS Axis 2", title="Lab vs Pilot RABR NMDS") +
    scale_color_manual(name=NULL, breaks=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control"),
                       labels=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control"),
                       values=c(Pilot_color, CVWRF_color, Labrabr_color, TF_color, GHR_color, Control_color)) +
    scale_fill_manual(name=NULL, breaks=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control"),
                      labels=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control"),
                      values=c(Pilot_color, CVWRF_color, Labrabr_color, TF_color, GHR_color, Control_color)) +
    theme_classic() + xlim(-0.75, 0.75) + ylim(-0.75, 0.75) +
    theme(legend.key.size = unit(0.25, "cm"),
          legend.position = c(0.85, 0.9),
          legend.background = element_rect(fill="NA",
                                           color="black"),
          legend.margin = margin(t=-2, r=3, b=3, l=3),
          plot.title=element_text(hjust=0.5))
  
  ggsave(paste("ITSfinal2/ITSplots/ITS_", name, "_nmds_pilotv81RABR.tiff", sep=""), width=5, height=4)
  
  #Beta Diversity PilotvLabRABRvTF
  samples_wanted <- c("R27_11_3_21_ITS", "R36_10_28_21_ITS", "R43_11_15_21_ITS", 
                      "R45_10_18_21_ITS", "R45_11_15_21_ITS", "R46_11_5_21_ITS", "R58_10_28_21_ITS", 
                      "R60_11_1_21_ITS", "R60_11_15_21_ITS", "11_15_21_ITS", "R72_11_15_21_ITS", 
                      "R75_11_15_21_ITS", "10_5_ITS", 
                      "19_ITS", "26_ITS", "11S_ITS", "11R_ITS", "S1_ITS", 
                      "S2_ITS", "S3_ITS", "TF_5_25_22_ITS", "TF_6_9_22_ITS", 
                      "TF_6_22_22_ITS", "TF_7_6_21_ITS", "TF_9_11_21_ITS", "TF_11_9_21_ITS")
  
  
  all_shared <- read_tsv(paste("ITSfinal2/final_",name,".agc.shared", sep="")) %>%
    select(Group, starts_with("Otu")) %>%
    filter(Group %in% samples_wanted) %>%
    pivot_longer(-Group) %>%
    group_by(Group) %>%
    # mutate(total=sum(value)) %>%
    # arrange(total) %>%
    group_by(name) %>%
    mutate(total=sum(value)) %>%
    filter(total != 0) %>%
    # summarize(total=sum(value)) %>% filter(total ==0) %>%
    ungroup() %>%
    select(-total) %>%
    pivot_wider(Group) %>%
    as.data.frame
  
  #all_shared <- as.data.frame(all_shared)
  
  rownames(all_shared) <- all_shared$Group
  all_shared <- all_shared[, -1]
  all_shared <- as.matrix(all_shared)
  
  
  all_dist <- avgdist(all_shared, dmethod="bray", iterations=1000, sample=2073)
  all_nmds <- metaMDS(all_dist) #low stress > insufficient data?
  
  all_scores <- scores(all_nmds) %>%
    as_tibble(rownames = "Group") 
  
  all_metadata <- gather_metadata("section", "date", "label", 1)
  
  for(i in 2:6) {
    all_metadata <- dplyr::bind_rows(all_metadata, gather_metadata("section", "date", "label", i))
  }
  
  all_metadata_nmds <- inner_join(all_metadata, all_scores, by=c('sample'='Group')) %>%
    mutate(section = factor(section,
                            levels=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control")))
  
  Pilot_color <- "#BEBEBE"
  CVWRF_color <- "#0000FF"
  Labrabr_color <- "#FF0000"
  TF_color <- "#00FF00"
  GHR_color <- "cyan"
  Control_color <- "#BE00FF"
  
  all_centroid <- all_metadata_nmds %>%
    group_by(section) %>%
    summarize(NMDS1 = mean(NMDS1),
              NMDS2 = mean(NMDS2), .groups="drop")
  
  all_star <- all_metadata_nmds %>%
    group_by(section) %>%
    mutate(centroid1 = mean(NMDS1),
           centroid2 = mean(NMDS2)) %>%
    ungroup()
  
  ggplot(all_star, aes(x=NMDS1, xend=centroid1,
                       y=NMDS2, yend=centroid2, color=section)) +
    geom_segment() +
    geom_point() +
    geom_point(data=all_centroid,
               mapping=aes(x=NMDS1, fill=section, y=NMDS2),
               shape=22, size=5, show.legend=FALSE,
               #color="black",
               inherit.aes=FALSE) +
    #coord_fixed(xlim=c(-1, 1), ylim=c(1, 1)) +
    labs(x="NMDS Axis 1", y="NMDS Axis 2", title="NMDS by RABR") +
    scale_color_manual(name=NULL, breaks=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control"),
                       labels=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control"),
                       values=c(Pilot_color, CVWRF_color, Labrabr_color, TF_color, GHR_color, Control_color)) +
    scale_fill_manual(name=NULL, breaks=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control"),
                      labels=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control"),
                      values=c(Pilot_color, CVWRF_color, Labrabr_color, TF_color, GHR_color, Control_color)) +
    theme_classic() + xlim(-0.75, 0.75) + ylim(-0.75, 0.75) +
    theme(legend.key.size = unit(0.25, "cm"),
          legend.position = c(0.85, 0.9),
          legend.background = element_rect(fill="NA",
                                           color="black"),
          legend.margin = margin(t=-2, r=3, b=3, l=3),
          plot.title=element_text(hjust=0.5))
  
  ggsave(paste("ITSfinal2/ITSplots/ITS_", name, "_nmds_pilotv81RABRvTF.tiff", sep=""), width=5, height=4)
  
  #Beta Diversity PilotvLabRABRvTFvCVWRF
  samples_wanted <- c("R27_11_3_21_ITS", "R36_10_28_21_ITS", "R43_11_15_21_ITS", 
                      "R45_10_18_21_ITS", "R45_11_15_21_ITS", "R46_11_5_21_ITS", "R58_10_28_21_ITS", 
                      "R60_11_1_21_ITS", "R60_11_15_21_ITS", "11_15_21_ITS", "R72_11_15_21_ITS", 
                      "R75_11_15_21_ITS", "10_5_ITS", 
                      "19_ITS", "26_ITS", "11S_ITS", "11R_ITS", "S1_ITS", 
                      "S2_ITS", "S3_ITS", "TF_5_25_22_ITS", "TF_6_9_22_ITS", 
                      "TF_6_22_22_ITS", "TF_7_6_21_ITS", "TF_9_11_21_ITS", "TF_11_9_21_ITS",
                      "CVWRF_PR_6_9_22_ITS", "CVWRF_PSR_6_22_22_ITS")
  
  
  all_shared <- read_tsv(paste("ITSfinal2/final_",name,".agc.shared", sep="")) %>%
    select(Group, starts_with("Otu")) %>%
    filter(Group %in% samples_wanted) %>%
    pivot_longer(-Group) %>%
    group_by(Group) %>%
    # mutate(total=sum(value)) %>%
    # arrange(total) %>%
    group_by(name) %>%
    mutate(total=sum(value)) %>%
    filter(total != 0) %>%
    # summarize(total=sum(value)) %>% filter(total ==0) %>%
    ungroup() %>%
    select(-total) %>%
    pivot_wider(Group) %>%
    as.data.frame
  
  #all_shared <- as.data.frame(all_shared)
  
  rownames(all_shared) <- all_shared$Group
  all_shared <- all_shared[, -1]
  all_shared <- as.matrix(all_shared)
  
  
  all_dist <- avgdist(all_shared, dmethod="bray", iterations=1000, sample=2073)
  all_nmds <- metaMDS(all_dist) #low stress > insufficient data?
  
  all_scores <- scores(all_nmds) %>%
    as_tibble(rownames = "Group") 
  
  all_metadata <- gather_metadata("section", "date", "label", 1)
  
  for(i in 2:6) {
    all_metadata <- dplyr::bind_rows(all_metadata, gather_metadata("section", "date", "label", i))
  }
  
  all_metadata_nmds <- inner_join(all_metadata, all_scores, by=c('sample'='Group')) %>%
    mutate(section = factor(section,
                            levels=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control")))
  
  Pilot_color <- "#BEBEBE"
  CVWRF_color <- "#0000FF"
  Labrabr_color <- "#FF0000"
  TF_color <- "#00FF00"
  GHR_color <- "cyan"
  Control_color <- "#BE00FF"
  
  all_centroid <- all_metadata_nmds %>%
    group_by(section) %>%
    summarize(NMDS1 = mean(NMDS1),
              NMDS2 = mean(NMDS2), .groups="drop")
  
  all_star <- all_metadata_nmds %>%
    group_by(section) %>%
    mutate(centroid1 = mean(NMDS1),
           centroid2 = mean(NMDS2)) %>%
    ungroup()
  
  ggplot(all_star, aes(x=NMDS1, xend=centroid1,
                       y=NMDS2, yend=centroid2, color=section)) +
    geom_segment() +
    geom_point() +
    geom_point(data=all_centroid,
               mapping=aes(x=NMDS1, fill=section, y=NMDS2),
               shape=22, size=5, show.legend=FALSE,
               #color="black",
               inherit.aes=FALSE) +
    #coord_fixed(xlim=c(-1, 1), ylim=c(1, 1)) +
    labs(x="NMDS Axis 1", y="NMDS Axis 2", title="NMDS by RABR") +
    scale_color_manual(name=NULL, breaks=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control"),
                       labels=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control"),
                       values=c(Pilot_color, CVWRF_color, Labrabr_color, TF_color, GHR_color, Control_color)) +
    scale_fill_manual(name=NULL, breaks=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control"),
                      labels=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control"),
                      values=c(Pilot_color, CVWRF_color, Labrabr_color, TF_color, GHR_color, Control_color)) +
    theme_classic() + xlim(-0.75, 0.75) + ylim(-0.75, 0.75) +
    theme(legend.key.size = unit(0.25, "cm"),
          legend.position = c(0.85, 0.9),
          legend.background = element_rect(fill="NA",
                                           color="black"),
          legend.margin = margin(t=-2, r=3, b=3, l=3),
          plot.title=element_text(hjust=0.5))
  
  ggsave(paste("ITSfinal2/ITSplots/ITS_", name, "_nmds_pilotv81RABRvTFvCVWRF.tiff", sep=""), width=5, height=4)
  
  #Beta Diversity PilotvLabRABRvTFvCVWRF
  samples_wanted <- c("R27_11_3_21_ITS", "R36_10_28_21_ITS", "R43_11_15_21_ITS", 
                      "R45_10_18_21_ITS", "R45_11_15_21_ITS", "R46_11_5_21_ITS", "R58_10_28_21_ITS", 
                      "R60_11_1_21_ITS", "R60_11_15_21_ITS", "11_15_21_ITS", "R72_11_15_21_ITS", 
                      "R75_11_15_21_ITS", "10_5_ITS", 
                      "19_ITS", "26_ITS", "11S_ITS", "11R_ITS", "S1_ITS", 
                      "S2_ITS", "S3_ITS", "TF_5_25_22_ITS", "TF_6_9_22_ITS", 
                      "TF_6_22_22_ITS", "TF_7_6_21_ITS", "TF_9_11_21_ITS", "TF_11_9_21_ITS",
                      "CVWRF_PR_6_9_22_ITS", "CVWRF_PSR_6_22_22_ITS",
                      "GHR_6_15_22_ITS", "GHR_5_1_22_ITS")
  
  
  all_shared <- read_tsv(paste("ITSfinal2/final_",name,".agc.shared", sep="")) %>%
    select(Group, starts_with("Otu")) %>%
    filter(Group %in% samples_wanted) %>%
    pivot_longer(-Group) %>%
    group_by(Group) %>%
    # mutate(total=sum(value)) %>%
    # arrange(total) %>%
    group_by(name) %>%
    mutate(total=sum(value)) %>%
    filter(total != 0) %>%
    # summarize(total=sum(value)) %>% filter(total ==0) %>%
    ungroup() %>%
    select(-total) %>%
    pivot_wider(Group) %>%
    as.data.frame
  
  #all_shared <- as.data.frame(all_shared)
  
  rownames(all_shared) <- all_shared$Group
  all_shared <- all_shared[, -1]
  all_shared <- as.matrix(all_shared)
  
  
  all_dist <- avgdist(all_shared, dmethod="bray", iterations=1000, sample=2073)
  all_nmds <- metaMDS(all_dist) #low stress > insufficient data?
  
  all_scores <- scores(all_nmds) %>%
    as_tibble(rownames = "Group") 
  
  all_metadata <- gather_metadata("section", "date", "label", 1)
  
  for(i in 2:6) {
    all_metadata <- dplyr::bind_rows(all_metadata, gather_metadata("section", "date", "label", i))
  }
  
  all_metadata_nmds <- inner_join(all_metadata, all_scores, by=c('sample'='Group')) %>%
    mutate(section = factor(section,
                            levels=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control")))
  
  Pilot_color <- "#BEBEBE"
  CVWRF_color <- "#0000FF"
  Labrabr_color <- "#FF0000"
  TF_color <- "#00FF00"
  GHR_color <- "cyan"
  Control_color <- "#BE00FF"
  
  all_centroid <- all_metadata_nmds %>%
    group_by(section) %>%
    summarize(NMDS1 = mean(NMDS1),
              NMDS2 = mean(NMDS2), .groups="drop")
  
  all_star <- all_metadata_nmds %>%
    group_by(section) %>%
    mutate(centroid1 = mean(NMDS1),
           centroid2 = mean(NMDS2)) %>%
    ungroup()
  
  ggplot(all_star, aes(x=NMDS1, xend=centroid1,
                       y=NMDS2, yend=centroid2, color=section)) +
    geom_segment() +
    geom_point() +
    geom_point(data=all_centroid,
               mapping=aes(x=NMDS1, fill=section, y=NMDS2),
               shape=22, size=5, show.legend=FALSE,
               #color="black",
               inherit.aes=FALSE) +
    #coord_fixed(xlim=c(-1, 1), ylim=c(1, 1)) +
    labs(x="NMDS Axis 1", y="NMDS Axis 2", title="NMDS by RABR") +
    scale_color_manual(name=NULL, breaks=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control"),
                       labels=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control"),
                       values=c(Pilot_color, CVWRF_color, Labrabr_color, TF_color, GHR_color, Control_color)) +
    scale_fill_manual(name=NULL, breaks=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control"),
                      labels=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control"),
                      values=c(Pilot_color, CVWRF_color, Labrabr_color, TF_color, GHR_color, Control_color)) +
    theme_classic() + xlim(-0.75, 0.75) + ylim(-0.75, 0.75) +
    theme(legend.key.size = unit(0.25, "cm"),
          legend.position = c(0.85, 0.9),
          legend.background = element_rect(fill="NA",
                                           color="black"),
          legend.margin = margin(t=-2, r=3, b=3, l=3),
          plot.title=element_text(hjust=0.5))
  
  ggsave(paste("ITSfinal2/ITSplots/ITS_", name, "_nmds_pilotv81RABRvTFvCVWRFvGH.tiff", sep=""), width=5, height=4)
  
  # All samples nmds
  #Beta Diversity 
  all_shared <- read_tsv(paste("ITSfinal2/final_",name,".agc.shared", sep="")) %>%
    select(Group, starts_with("Otu")) %>%
    #filter(Group %in% samples_wanted) %>%
    pivot_longer(-Group) %>%
    group_by(Group) %>%
    # mutate(total=sum(value)) %>%
    # arrange(total) %>%
    group_by(name) %>%
    mutate(total=sum(value)) %>%
    filter(total != 0) %>%
    # summarize(total=sum(value)) %>% filter(total ==0) %>%
    ungroup() %>%
    select(-total) %>%
    pivot_wider(Group) %>%
    as.data.frame
  
  #all_shared <- as.data.frame(all_shared)
  
  rownames(all_shared) <- all_shared$Group
  all_shared <- all_shared[, -1]
  all_shared <- as.matrix(all_shared)
  
  
  all_dist <- avgdist(all_shared, dmethod="bray", iterations=1000, sample=2073)
  all_nmds <- metaMDS(all_dist) #low stress > insufficient data?
  
  all_scores <- scores(all_nmds) %>%
    as_tibble(rownames = "Group")
  
  all_metadata <- gather_metadata("section", "date", "label", 1)
  
  for(i in 2:6) {
    all_metadata <- dplyr::bind_rows(all_metadata, gather_metadata("section", "date", "label", i))
  }
  
  all_metadata_nmds <- inner_join(all_metadata, all_scores, by=c('sample'='Group')) %>%
    mutate(section = factor(section,
                            levels=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control")))
  
  Pilot_color <- "#BEBEBE"
  CVWRF_color <- "#0000FF"
  Labrabr_color <- "#FF0000"
  TF_color <- "#00FF00"
  GHR_color <- "cyan"
  Control_color <- "#BE00FF"
  
  all_centroid <- all_metadata_nmds %>%
    group_by(section) %>%
    summarize(NMDS1 = mean(NMDS1),
              NMDS2 = mean(NMDS2), .groups="drop")
  
  all_star <- all_metadata_nmds %>%
    group_by(section) %>%
    mutate(centroid1 = mean(NMDS1),
           centroid2 = mean(NMDS2)) %>%
    ungroup()
  
  write.csv(all_star,paste("ITSfinal2/ITScsvs/IT_", name, "_beta_all.csv", sep=""), row.names = FALSE)
  
  
  ggplot(all_star, aes(x=NMDS1, xend=centroid1,
                       y=NMDS2, yend=centroid2, color=section)) +
    geom_segment() +
    geom_point() +
    geom_point(data=all_centroid,
               mapping=aes(x=NMDS1, fill=section, y=NMDS2),
               shape=22, size=5, show.legend=FALSE,
               #color="black",
               inherit.aes=FALSE) +
    #coord_fixed(xlim=c(-1, 1), ylim=c(1, 1)) +
    labs(x="NMDS Axis 1", y="NMDS Axis 2", title="All Samples NMDS by RABR") +
    scale_color_manual(name=NULL, breaks=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control"),
                       labels=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control"),
                       values=c(Pilot_color, CVWRF_color, Labrabr_color, TF_color, GHR_color, Control_color)) +
    scale_fill_manual(name=NULL, breaks=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control"),
                      labels=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control"),
                      values=c(Pilot_color, CVWRF_color, Labrabr_color, TF_color, GHR_color, Control_color)) +
    theme_classic() + xlim(-0.75, 0.75) + ylim(-0.75, 0.75) +
    theme(legend.key.size = unit(0.25, "cm"),
          legend.position = c(0.85, 0.9),
          legend.background = element_rect(fill="NA",
                                           color="black"),
          legend.margin = margin(t=2, r=3, b=3, l=3),
          plot.title=element_text(hjust=0.5))
  
  ggsave(paste("ITSfinal2/ITSplots/ITS_", name, "_nmds_all.tiff", sep=""), width=5, height=4)
  
  # Temp comparison
  # select for Lab samples
  samples_wanted <- c("R27_11_3_21_ITS", "R36_10_28_21_ITS", "R43_11_15_21_ITS", 
                      "R45_10_18_21_ITS", "R45_11_15_21_ITS", "R46_11_5_21_ITS", "R58_10_28_21_ITS", 
                      "R60_11_1_21_ITS", "R60_11_15_21_ITS", "R7_11_15_21_ITS", "R72_11_15_21_ITS", "R75_11_15_21_ITS")
  
  
  temp_shared <- read_tsv(paste("ITSfinal2/final_",name,".agc.shared", sep="")) %>%
    select(Group, starts_with("Otu")) %>%
    filter(Group %in% samples_wanted) %>%
    pivot_longer(-Group) %>%
    group_by(Group) %>%
    # mutate(total=sum(value)) %>%
    # arrange(total)
    # filter(total > 29272) %>%
    group_by(name) %>%
    mutate(total=sum(value)) %>% 
    filter(total != 0) %>%
    ungroup() %>%
    select(-total) %>%
    pivot_wider(Group) %>%
    as.data.frame
  
  rownames(temp_shared) <- temp_shared$Group
  temp_shared <- temp_shared[, -1]
  temp_shared <- as.matrix(temp_shared)
  
  temp_dist <- avgdist(temp_shared, dmethod="bray", iterations=1000, sample=2073)
  temp_nmds <- metaMDS(temp_dist)
  
  temp_scores <- scores(temp_nmds) %>%
    as_tibble(rownames = "Group") 
  
  temp_metadata <- gather_metadata("temp", "date", "label", 2)
  
  temp_metadata_nmds <- inner_join(temp_metadata, temp_scores, by=c('sample'='Group')) %>%
    mutate(temp = factor(temp,
                                levels=c("10",
                                         "20",
                                         "25")))
  
  T1_color <- "#BEBEBE"
  T2_color <- "#0000FF"
  T3_color <- "#FF0000"
  
  temp_centroid <- temp_metadata_nmds %>%
    group_by(temp) %>%
    summarize(NMDS1 = mean(NMDS1),
              NMDS2 = mean(NMDS2), .groups="drop")
  
  temp_star <- temp_metadata_nmds %>%
    group_by(temp) %>%
    mutate(centroid1 = mean(NMDS1),
           centroid2 = mean(NMDS2)) %>%
    ungroup()
  
  ggplot(temp_star, aes(x=NMDS1, xend=centroid1,
                        y=NMDS2, yend=centroid2, color=temp)) + 
    geom_segment() +
    geom_point() +
    geom_point(data=temp_centroid, 
               mapping=aes(x=NMDS1, fill=temp, y=NMDS2),
               shape=22, size=5, show.legend=FALSE,
               inherit.aes=FALSE) +
    #coord_fixed(xlim=c(-1, 1), ylim=c(1, 1)) +
    labs(x="NMDS Axis 1", y="NMDS Axis 2", title="Lab-scale RABRs NMDS
         by Temperature") +
    scale_color_manual(name=NULL, 
                       breaks=c("10",
                                "20",
                                "25"),
                       labels=c("10 C",
                                "20 C",
                                "25 C"),
                       values=c(T1_color, T2_color, T3_color)) +
    scale_fill_manual(name=NULL, 
                      breaks=c("10",
                               "20",
                               "25"),
                      labels=c("10 C",
                               "20 C",
                               "25 C"),
                      values=c(T1_color, T2_color, T3_color)) +
    theme_classic() + xlim(-0.75, 0.75) + ylim(-0.75, 0.75) +
    theme(legend.key.size = unit(0.35, "cm"),
          legend.position = c(0.85, 0.9),
          legend.background = element_rect(fill="NA",
                                           color="black"),
          legend.margin = margin(t=3, r=3, b=3, l=3),
          plot.title=element_text(hjust=0.5))
  
  ggsave(paste("ITSfinal2/ITSplots/ITS_", name, "_nmds_lab_temp.tiff", sep=""), width=5, height=4)
}

beta_div("pr2_euk")
#beta_div("pr2_fungi")
#beta_div("silv_euk")
#beta_div("silv_fungi")
