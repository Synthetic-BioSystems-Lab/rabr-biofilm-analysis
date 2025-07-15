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

# import metadata
all_metadata <- gather_metadata("section", "date", "label", 1) 

for(i in 2:6) {
  all_metadata <- dplyr::bind_rows(all_metadata, gather_metadata("section", "date", "label", i))
}

all_metadata %>%
  mutate(section = factor(section,
                          levels=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control")))

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
alpha_count <- inner_join(composite, all_metadata, by=c('sample_id'='sample')) %>%
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

ggsave("23S_D_plots/23SD_alpha_div0_all.tiff", width=8, height=8)

alpha_count %>%
  filter(section == "pilot") %>%
  ggplot(aes(x=n, y=value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~metric, nrow=4, scales="free_y")

ggsave("23S_D_plots/23SD_alpha_div0_pilot.tiff", width=8, height=8)

alpha_count %>%
  filter(section == "81RABR") %>%
  ggplot(aes(x=n, y=value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~metric, nrow=4, scales="free_y")

ggsave("23S_D_plots/23SD_alpha_div0_81RABR.tiff", width=8, height=8)

alpha_count %>%
  filter(section == "TF") %>%
  ggplot(aes(x=n, y=value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~metric, nrow=4, scales="free_y")

ggsave("23S_D_plots/23SD_alpha_div0_TF.tiff", width=8, height=8)

alpha_count %>%
  filter(section == "control") %>%
  ggplot(aes(x=n, y=value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~metric, nrow=4, scales="free_y")

ggsave("23S_D_plots/23SD_alpha_div0_control.tiff", width=8, height=8)

# Amanda Method
dli_pilot_metadata <- gather_metadata("dli_level","date", "label", 1) %>%
  mutate(dli_level = factor(dli_level,
                       levels=c("high",
                                "low")))

alpha_diversity <- inner_join(composite, all_metadata, by=c('sample_id'='sample')) %>%
  group_by(sample_id, section, date, label) %>%
  summarize(sobs = richness(count),
            shannon = shannon(count),
            simpson=simpson(count),
            invsimpson = 1/simpson)

dli_metadata_alpha <- dli_pilot_metadata %>%
  select(sample, dli_level) %>%
  inner_join(., alpha_diversity, by=c('sample'='sample_id'))

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

write.csv(dli_metadata_alpha,"23Scsvs/23S_alpha_dli.csv", row.names = FALSE)

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
  ylim(0, 45) +
  theme(axis.text.x = element_markdown(), plot.title=element_text(hjust=0.5))

ggsave("23S_D_plots/23SD_alpha_div_pilot_dli.tiff", width=5, height=4)

# Amanda Method all Samples
#EXPAND METADATA

all_metadata_alpha <- alpha_diversity

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

breaks <- c("4_pilot", "5_CVWRF", "2_81RABR", "6_TF", "3_GHR", "1_control")
labels <- c("Pilot RABR", "CVWRF", "Lab-scale RABRs", "Trickling Filter", "GHR", "Control")

all_metadata_alpha <- all_metadata_alpha %>%
  mutate(section = str_replace_all(section, "control", "1_control")) %>%
  mutate(section = str_replace_all(section, "81RABR", "2_81RABR")) %>%
  mutate(section = str_replace_all(section, "GHR", "3_GHR")) %>%
  mutate(section = str_replace_all(section, "pilot", "4_pilot")) %>%
  mutate(section = str_replace_all(section, "CVWRF", "5_CVWRF")) %>%
  mutate(section = str_replace_all(section, "TF", "6_TF"))

write.csv(all_metadata_alpha,"23Scsvs/23S_alpha_all.csv", row.names = FALSE)

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
  ylim(0, 45) +
  theme(axis.text.x = element_markdown(), plot.title=element_text(hjust=0.5))

ggsave("23S_D_plots/23SD_alpha_div_all_sections.tiff", width=5, height=4)

## Amanda Method pilot vs 81RABR

#EXPAND METADATA
p_metadata_alpha <- all_metadata_alpha %>%
  filter(section == "4_pilot")
pvl_metadata_alpha <- all_metadata_alpha %>%
  filter(section == "2_81RABR") %>%
  rbind(.,p_metadata_alpha)

Pilot_color <- "#BEBEBE"
labRABR_color <- "#0000FF"

all_count <- pvl_metadata_alpha %>%
  #specify package
  dplyr::count(section)

pilot_n <- all_count %>%
  filter(section == "4_pilot") %>%
  pull(n)

labRABR_n <- all_count %>%
  filter(section == "2_81RABR") %>%
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
  scale_x_discrete(breaks=c("4_pilot", "2_81RABR"),
                   labels=c("Pilot", "Lab RABR")) +
  scale_fill_manual(name=NULL,
                    breaks=c("4_pilot", "2_81RABR"),
                    labels=c("Pilot", "Lab RABR"),
                    values=c(Pilot_color, labRABR_color)) +
  theme_classic() +
  ylim(0, 45) +
  theme(axis.text.x = element_markdown(), plot.title=element_text(hjust=0.5))

ggsave("23S_D_plots/23SD_alpha_div_pilotv81.tiff", width=5, height=4)

# Alpha Div: Productivity v Invsimpson
prod_metadata1 <- gather_metadata("section", "dry_productivity_substratum", "date", 1) %>%
  rename(productivity = dry_productivity_substratum) %>%
  select("sample", "productivity")
prod_meta_alpha1 <- inner_join(prod_metadata1, alpha_diversity,
                              by=c('sample'='sample_id'))
prod_meta_alpha1 %>%
  ggplot(aes(x=productivity, y=invsimpson)) +
  geom_point() + 
  labs(x="Pilot RABR Productivity (g/m2/day)", y="Inverse Simpson") +
  ggtitle("Productivity vs Inverse Simpson") +
  theme_classic() +
  ylim(0, 45) +
  theme(axis.text.x = element_markdown(), plot.title=element_text(hjust=0.5)) +
  stat_cor()+
  geom_smooth(method=lm, se=FALSE)
ggsave("23S_D_plots/23SD_alpha_div_vPilot_Productivity.tiff", width=3.5, height=2.3)

#lab rabr productivity v invsimpson
prod_metadata2 <- gather_metadata("section", "productivity", "date", 2) %>%
  select("sample", "productivity")
prod_meta_alpha2 <- inner_join(prod_metadata2, alpha_diversity,
                              by=c('sample'='sample_id'))
prod_meta_alpha2 %>%
  ggplot(aes(x=productivity, y=invsimpson)) +
  geom_point() + 
  labs(x="Lab RABR Productivity (g/m2/day)", y="Inverse Simpson") +
  ggtitle("Productivity vs Inverse Simpson") +
  theme_classic() +
  ylim(0, 30) +
  xlim(0, NA) +
  theme(axis.text.x = element_markdown(), plot.title=element_text(hjust=0.5)) +
  stat_cor()+
  geom_smooth(method=lm, se=FALSE)
ggsave("23S_D_plots/23SD_alpha_div_v81RABR_Productivity.tiff", width=3.5, height=2.3)

# pilot RABR date productivity
date_meta_alpha <- all_metadata_alpha %>%
  filter(section == "4_pilot") %>%
  select(sample_id, date, invsimpson)

RnS <- date_meta_alpha %>%
  filter(sample_id %in% c("11S_23S", "11R_23S")) %>%
  group_by(sample_id) %>%
  pivot_wider(names_from = sample_id, values_from = invsimpson) %>%
  data.table::setnames("11S_23S",'S') %>%
  data.table::setnames("11R_23S",'R') %>%
  mutate(sample_id = "11_2_23S", invsimpson = ((S + R)/2)) %>%
  select(sample_id, date, invsimpson)

S123 <- date_meta_alpha %>%
  filter(sample_id %in% c("S1_23S", "S2_23S", "S3_23S")) %>%
  group_by(sample_id) %>%
  pivot_wider(names_from = sample_id, values_from = invsimpson) %>%
  mutate(sample_id = "10_12_23S", invsimpson = ((S1_23S + S2_23S + S3_23S)/3)) %>%
  select(sample_id, date, invsimpson)

mod_date_meta_alpha <- date_meta_alpha %>%
  filter(sample_id %in% c("10_5_23S", "19_23S", "26_23S")) %>%
  rbind(.,RnS) %>%
  rbind(.,S123) %>%
  ungroup() %>%
  select(sample_id, date, invsimpson)


mod_date_meta_alpha %>%
  ggplot(aes(x=date, y=invsimpson)) +
  geom_point() + 
  labs(x="Date", y="Inverse Simpson") +
  ggtitle("Timeline of Inverse Simpson") +
  theme_classic() +
  ylim(0, 45) +
  theme(plot.title=element_text(hjust=0.5),
        axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
ggsave("23S_D_plots/23SD_alpha_div_timeline.tiff", width=4, height=2.3)



# Div vs Light?
# Div vs Temp?

# 11R v 11S smooth is more diverse than rough, though rough has much higher productivity
# Sample  InvSimpson
# 11R_23S_14.69596
# 11S_23S_19.28912

