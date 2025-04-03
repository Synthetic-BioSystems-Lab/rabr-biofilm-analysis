library(tidyverse)
library(ggplot2)
library(ggtext)
library(readxl)
library(NLP)
library(ggpubr)


setwd("~/Miller Lab/Rscripts_PilotRABR")

gather_metadata <- function(target, t2, t3, page) {
  read_excel("ITSpilot/ITS_metadata.xlsx", sheet=page) %>%
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

name <- "pr2_euk"

slow_alpha_div <- function(name) {
  #import otu counts
  otu_counts <- read_tsv(paste("ITSpilot/final_",name,".agc.shared", sep="")) %>%
    select(-label, -numOtus) %>%
    pivot_longer(-Group, names_to = "otu", values_to = "count") %>%
    rename(sample_id = Group)
  
  # Alpha Diversity
  
  richness <- function(x) {
    sum(x > 0)
  }
  
  shannon <- function(x) {
    rabund <-x[x>0]/sum(x)
    -sum(rabund * log(rabund))
  }
  simpson <- function(x){
    n <- sum(x)
    sum(x*(x-1) / (n*(n-1)))
  }
  alpha_count <- inner_join(otu_counts, all_metadata, by=c('sample_id'='sample')) %>%
    group_by(sample_id, section, date, label) %>%
    summarize(sobs = richness(count),
              shannon = shannon(count),
              simpson=simpson(count),
              invsimpson = 1/simpson,
              n = num(count)) %>% #n is number of individuals in each sample
    pivot_longer(cols = c(sobs, shannon, invsimpson, simpson),
                 names_to = "metric")
  
  alpha_count %>%
    ggplot(aes(x=n, y=value)) +
    geom_point() +
    geom_smooth() +
    facet_wrap(~metric, nrow=4, scales="free_y")
  
  ggsave(paste("ITSpilot/ITSplots/ITS_", name, "_alpha_div0_all.tiff", sep=""), width=8, height=8)
  
  alpha_count %>%
    filter(section == "pilot") %>%
    ggplot(aes(x=n, y=value)) +
    geom_point() +
    geom_smooth() +
    facet_wrap(~metric, nrow=4, scales="free_y")
  
  ggsave(paste("ITSpilot/ITSplots/ITS_", name, "_alpha_div0_pilot.tiff", sep=""), width=8, height=8)
  
  alpha_count %>%
    filter(section == "81RABR") %>%
    ggplot(aes(x=n, y=value)) +
    geom_point() +
    geom_smooth() +
    facet_wrap(~metric, nrow=4, scales="free_y")
  
  ggsave(paste("ITSpilot/ITSplots/ITS_", name, "_alpha_div0_81RABR.tiff", sep=""), width=8, height=8)
  
  alpha_count %>%
    filter(section == "TF") %>%
    ggplot(aes(x=n, y=value)) +
    geom_point() +
    geom_smooth() +
    facet_wrap(~metric, nrow=4, scales="free_y")
  
  ggsave(paste("ITSpilot/ITSplots/ITS_", name, "_alpha_div0_TF.tiff", sep=""), width=8, height=8)
  
  alpha_count %>%
    filter(section == "control") %>%
    ggplot(aes(x=n, y=value)) +
    geom_point() +
    geom_smooth() +
    facet_wrap(~metric, nrow=4, scales="free_y")
  
  ggsave(paste("ITSpilot/ITSplots/ITS_", name, "_alpha_div0_control.tiff", sep=""), width=8, height=8)
}

alpha_div <- function(name) {
  # Amanda Method
  #import otu counts
  otu_counts <- read_tsv(paste("ITSpilot/final_",name,".agc.shared", sep="")) %>%
    select(-label, -numOtus) %>%
    pivot_longer(-Group, names_to = "otu", values_to = "count") %>%
    rename(sample_id = Group)
  
  # Alpha Diversity
  
  richness <- function(x) {
    sum(x > 0)
  }
  
  shannon <- function(x) {
    rabund <-x[x>0]/sum(x)
    -sum(rabund * log(rabund))
  }
  simpson <- function(x){
    n <- sum(x)
    sum(x*(x-1) / (n*(n-1)))
  }
  
  dli_pilot_metadata <- gather_metadata("dli_level","date", "label", 1) %>%
    mutate(dli_level = factor(dli_level,
                              levels=c("high",
                                       "low"))) %>%
    select(sample, dli_level)
  
  alpha_diversity <- inner_join(otu_counts, all_metadata, by=c('sample_id'='sample')) %>%
    group_by(sample_id, section, date, label) %>%
    summarize(sobs = richness(count),
              shannon = shannon(count),
              simpson=simpson(count),
              invsimpson = 1/simpson)
  
  dli_metadata_alpha <- inner_join(dli_pilot_metadata, alpha_diversity,
                                    by=c('sample'='sample_id'))
  
  #Five_color <- "#BEBEBE"
  high_color <- "#0000FF"
  low_color <- "#FF0000"
  #TwentySix_color <- "#00FF00"
  #Two_color <- "cyan"
  
  dli_count <- dli_metadata_alpha %>%
    #specify package
    dplyr::count(dli_level)
  
  high_n <- dli_count %>%
    filter(dli_level == "high") %>%
    pull(n)
  
  low_n <- dli_count %>%
    filter(dli_level == "low") %>%
    pull(n)
  
  
  # box plot with jitters
  dli_metadata_alpha %>%
    ggplot(aes(x=dli_level, y=invsimpson, fill=dli_level)) +
    # geom_boxplot(show.legend=FALSE, outlier.shape=NA, alpha=0.25, 
    #              width=0.6, coef=0) +
    stat_summary(fun.data=median_hilow, fun.args=0.5, show.legend=FALSE, 
                 geom="crossbar", alpha=0.25, width=0.6) +
    geom_jitter(show.legend=FALSE, width=0.25, shape=21, color="black") +
    labs(x=NULL, y="Inverse Simpson Index") +
    ggtitle("All Samples Alpha Diversity by DLI") +
    scale_x_discrete(breaks=c("high", "low"),
                     labels=c("high", "low")) +
    scale_fill_manual(name=NULL,
                      breaks=c("high", "low"),
                      labels=c("high", "low"),
                      values=c(high_color, low_color)) +
    theme_classic() +
    ylim(0, 10) +
    theme(axis.text.x = element_markdown(), plot.title=element_text(hjust=0.5))
  
  ggsave(paste("ITSpilot/ITSplots/ITS_", name, "_alpha_div_pilot_dli.tiff", sep=""), width=5, height=4)
  
  # Amanda Method all Samples
  #EXPAND METADATA
  
  all_alpha_diversity <- inner_join(otu_counts, all_metadata, by=c('sample_id'='sample')) %>%
    group_by(sample_id, section, date, label) %>%
    summarize(sobs = richness(count),
              shannon = shannon(count),
              simpson=simpson(count),
              invsimpson = 1/simpson)
  
  all_metadata_alpha <- all_alpha_diversity
  
  Pilot_color <- "#BEBEBE"
  CVWRF_color <- "#0000FF"
  Labrabr_color <- "#FF0000"
  TF_color <- "#00FF00"
  GHR_color <- "cyan"
  Control_color <- "#BE00FF"
  
  all_count <- all_metadata_alpha %>%
    #specify package
    dplyr::count(section)
  
  pilot_n <- all_count %>%
    filter(section == "pilot") %>%
    pull(n)
  
  CVWRF_n <- all_count %>%
    filter(section == "CVWRF") %>%
    pull(n)
  
  Labrabr_n <- all_count %>%
    filter(section == "81RABR") %>%
    pull(n)
  
  TF_n <- all_count %>%
    filter(section == "TF") %>%
    pull(n)
  
  GHR_n <- all_count %>%
    filter(section == "GHR") %>%
    pull(n)
  
  control_n <- all_count %>%
    filter(section == "control") %>%
    pull(n)
  
  breaks <- c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control")
  labels <- c("Pilot RABR", "CVWRF", "Lab-scale RABRs", "Trickling Filter", "GHR", "Control")
  # box plot with jitters
  all_metadata_alpha %>%
    ggplot(aes(x=section, y=invsimpson, fill=section)) +
    # geom_boxplot(show.legend=FALSE, outlier.shape=NA, alpha=0.25, 
    #              width=0.6, coef=0) +
    stat_summary(fun.data=median_hilow, fun.args=0.5, show.legend=FALSE, 
                 geom="crossbar", alpha=0.25, width=0.6) +
    geom_jitter(show.legend=FALSE, width=0.25, shape=21, color="black") +
    labs(x=NULL, y="Inverse Simpson Index") +
    ggtitle("All Samples Alpha Diversity by Section") +
    scale_x_discrete(breaks=breaks,
                     labels=labels) +
    scale_fill_manual(name=NULL,
                      breaks=breaks,
                      labels=labels,
                      values=c(Pilot_color, CVWRF_color, Labrabr_color, TF_color, GHR_color, Control_color)) +
    theme_classic() +
    ylim(0, 15) +
    theme(axis.text.x = element_markdown(), plot.title=element_text(hjust=0.5))
  
  ggsave(paste("ITSpilot/ITSplots/ITS_", name, "_alpha_div_all_sections.tiff", sep=""), width=5, height=4)
  
  ## Amanda Method pilot vs 81RABR
  
  #EXPAND METADATA
  
  p_metadata_alpha <- all_metadata %>%
    filter(section == "pilot")
  pvl_metadata_alpha <- all_metadata %>%
    filter(section == "81RABR") %>%
    rbind(.,p_metadata_alpha) %>%
    select(sample) %>%
    inner_join(., all_alpha_diversity, by=c("sample"="sample_id"))
  
  Pilot_color <- "#BEBEBE"
  labRABR_color <- "#0000FF"
  
  all_count <- pvl_metadata_alpha %>%
    #specify package
    dplyr::count(section)
  
  pilot_n <- all_count %>%
    filter(section == "pilot") %>%
    pull(n)
  
  labRABR_n <- all_count %>%
    filter(section == "81RABR") %>%
    pull(n)
  
  # box plot with jitters
  pvl_metadata_alpha %>%
    ggplot(aes(x=section, y=invsimpson, fill=section)) +
    # geom_boxplot(show.legend=FALSE, outlier.shape=NA, alpha=0.25, 
    #              width=0.6, coef=0) +
    stat_summary(fun.data=median_hilow, fun.args=0.5, show.legend=FALSE, 
                 geom="crossbar", alpha=0.25, width=0.6) +
    geom_jitter(show.legend=FALSE, width=0.25, shape=21, color="black") +
    labs(x=NULL, y="Inverse Simpson Index") +
    ggtitle("Standard vs Lab RABR Alpha Diversity by Section") +
    scale_x_discrete(breaks=c("pilot", "81RABR"),
                     labels=c("pilot", "81RABR")) +
    scale_fill_manual(name=NULL,
                      breaks=c("pilot", "81RABR"),
                      labels=c("pilot", "81RABR"),
                      values=c(Pilot_color, labRABR_color)) +
    theme_classic() +
    ylim(0, 10) +
    theme(axis.text.x = element_markdown(), plot.title=element_text(hjust=0.5))
  
  ggsave(paste("ITSpilot/ITSplots/ITS_", name, "_alpha_div_pilotv81.tiff", sep=""), width=5, height=4)
  
  # Alpha Div: Productivity v Invsimpson
  prod_metadata1 <- gather_metadata("section", "dry_productivity_substratum", "date", 1) %>%
    rename(productivity = dry_productivity_substratum)
  prod_meta_alpha1 <- inner_join(prod_metadata1, alpha_diversity,
                                 by=c('sample'='sample_id'))
  prod_meta_alpha1 %>%
    ggplot(aes(x=productivity, y=invsimpson)) +
    geom_point() + 
    labs(x="Pilot RABR Productivity (g/m2/day)", y="Inverse Simpson") +
    ggtitle("Productivity vs Inverse Simpson") +
    theme_classic() +
    ylim(0, 20) +
    theme(axis.text.x = element_markdown(), plot.title=element_text(hjust=0.5)) +
    stat_cor()+
    geom_smooth(method=lm, se=FALSE)
  ggsave(paste("ITSpilot/ITSplots/ITS_", name, "_alpha_div_vPilot_Productivity.tiff", sep=""), width=3.5, height=2.3)
  
  #lab rabr productivity v invsimpson
  prod_metadata2 <- gather_metadata("section", "productivity", "date", 2)
  prod_meta_alpha2 <- inner_join(prod_metadata2, alpha_diversity,
                                 by=c('sample'='sample_id'))
  prod_meta_alpha2 %>%
    ggplot(aes(x=productivity, y=invsimpson)) +
    geom_point() + 
    labs(x="Lab RABR Productivity (g/m2/day)", y="Inverse Simpson") +
    ggtitle("Productivity vs Inverse Simpson") +
    theme_classic() +
    ylim(0, 15) +
    xlim(0, 15) +
    theme(axis.text.x = element_markdown(), plot.title=element_text(hjust=0.5)) +
    stat_cor()+
    geom_smooth(method=lm, se=FALSE)
  ggsave(paste("ITSpilot/ITSplots/ITS_", name, "_alpha_div_v81RABR_Productivity.tiff", sep=""), width=3.5, height=2.3)
  
  # pilot RABR date productivity
  date_meta_alpha <- all_metadata_alpha %>%
    filter(section == "pilot") %>%
    select(sample_id, date, invsimpson) %>%
    rename(sample=sample_id)
  
  RnS <- date_meta_alpha %>%
    filter(sample %in% c("11S_ITS", "11R_ITS")) %>%
    group_by(sample) %>%
    pivot_wider(names_from = sample, values_from = invsimpson) %>%
    data.table::setnames("11S_ITS",'S') %>%
    data.table::setnames("11R_ITS",'R') %>%
    mutate(sample = "11_2_ITS", invsimpson = ((S + R)/2)) %>%
    select(sample, date, invsimpson)
  
  S123 <- date_meta_alpha %>%
    filter(sample %in% c("S1_ITS", "S2_ITS", "S3_ITS")) %>%
    group_by(sample) %>%
    pivot_wider(names_from = sample, values_from = invsimpson) %>%
    mutate(sample = "10_12_ITS", invsimpson = ((S1_ITS + S2_ITS + S3_ITS)/3)) %>%
    select(sample, date, invsimpson)
  
  mod_date_meta_alpha <- date_meta_alpha %>%
    filter(sample %in% c("10_5_ITS", "19_ITS", "26_ITS")) %>%
    rbind(.,RnS) %>%
    rbind(.,S123)
  
  
  mod_date_meta_alpha %>%
    ggplot(aes(x=date, y=invsimpson)) +
    geom_point() + 
    labs(x="Date", y="Inverse Simpson") +
    ggtitle("Timeline of Inverse Simpson") +
    theme_classic() +
    ylim(0, 20) +
    theme(plot.title=element_text(hjust=0.5),
          axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
  ggsave(paste("ITSpilot/ITSplots/ITS_", name, "_alpha_div_timeline.tiff", sep=""), width=4, height=2.3)
  
}
# Div vs Light?
# Div vs Temp?

# 11R v 11S smooth is more diverse than rough, though rough has much higher productivity
# Sample  InvSimpson
# 11R_ITS_14.69596
# 11S_ITS_19.28912

slow_alpha_div("pr2_euk") #done
slow_alpha_div("pr2_fungi") #done
slow_alpha_div("silv_euk") #done
slow_alpha_div("silv_fungi") #done
alpha_div("pr2_euk")  #done
alpha_div("pr2_fungi") #done
alpha_div("silv_euk") #done
alpha_div("silv_fungi") #done
