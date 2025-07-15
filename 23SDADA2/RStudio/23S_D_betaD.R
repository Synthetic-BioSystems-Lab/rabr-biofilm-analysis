library(vegan)
library(ggtext)
library(tidyverse)
library(readxl)

setwd("~/Miller Lab/Rscripts_PilotRABR/DADA2/23SDADA2visualization")

gather_metadata <- function(target, t2, t3, page) {
  read_excel("23S_D_metadata.xlsx", sheet=page) %>%
    rename_all(tolower) %>%
    mutate(sample = str_replace_all(sample, "-", "_")) %>%
    select(sample, target, t2, t3)
}

#Beta Diversity 2
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

shared <- rbind(asv_seqtabA, asv_seqtabE)%>%
  pivot_wider(names_from = ASVs, values_from = count, values_fill = 0) %>%
  mutate(sample_id = str_replace_all(sample_id, "-", "_"))
 

# for(i in 1:nrow(all_shared)) {
#   row <- all_shared[i,]
#   row[1] <- gsub("_23S", "", row[1])
#   all_shared[i,1] <- row[1]
# }
samples_wanted <- c("10_5_23S", "19_23S","26_23S", "11S_23S", "11R_23S", "S1_23S", "S2_23S", "S3_23S")

all_shared <- shared %>%
  filter(sample_id %in% samples_wanted) %>%
  as.data.frame

rownames(all_shared) <- all_shared$sample_id
all_shared <- all_shared[, -1]
all_shared <- as.matrix(all_shared)


all_dist <- avgdist(all_shared, dmethod="bray", iterations=1000, sample=3800)
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

write.csv(date_star,"23Scsvs/23S_beta_date.csv", row.names = FALSE)

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
  theme_classic() + xlim(-0.75, 0.75) + ylim(-0.75, 0.75) +
  theme(legend.key.size = unit(0.25, "cm"),
        legend.position = c(0.85, 0.9),
        legend.background = element_rect(fill="NA",
                                         color="black"),
        legend.margin = margin(t=2, r=3, b=3, l=3),
        plot.title=element_text(hjust=0.5))

ggsave("23S_D_plots/23SD_nmds_all_date.tiff", width=5, height=4)

#Beta Diversity PilotvLabRABR
samples_wanted <- c("R27_11_3_21_23S", "R36_10_18_21_23S", "R43_11_15_21_23S", 
                    "R45_10_18_21_23S", "R45_11_15_21_23S", "R46_11_5_21_23S", 
                    "R58_10_28_21_23S", "R60_11_1_21_23S", "R60_11_15_21_23S", 
                    "R7_11_15_21_23S", "R72_11_15_21_23S", "R75_11_15_21_23S",
                    "10_5_23S", "19_23S",
                    "26_23S", "11S_23S", "11R_23S", "S1_23S", "S2_23S", "S3_23S")

all_shared <- shared %>%
  filter(sample_id %in% samples_wanted) %>%
  as.data.frame

#all_shared <- as.data.frame(all_shared)

rownames(all_shared) <- all_shared$sample_id
all_shared <- all_shared[, -1]
all_shared <- as.matrix(all_shared)


all_dist <- avgdist(all_shared, dmethod="bray", iterations=1000, sample=20073)
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

ggsave("23S_D_plots/23SD_nmds_pilotv81RABR.tiff", width=5, height=4)

#Beta Diversity PilotvLabRABRvTF
samples_wanted <- c("R27_11_3_21_23S", "R36_10_18_21_23S", "R43_11_15_21_23S", 
                    "R45_10_18_21_23S", "R45_11_15_21_23S", "R46_11_5_21_23S", 
                    "R58_10_28_21_23S", "R60_11_1_21_23S", "R60_11_15_21_23S", 
                    "R7_11_15_21_23S", "R72_11_15_21_23S", "R75_11_15_21_23S",
                    "10_5_23S", "19_23S",
                    "26_23S", "11S_23S", "11R_23S", "S1_23S", "S2_23S", "S3_23S",
                    "TF_5_25_22_23S", "TF_6_9_22_23S", "TF_6_22_22_23S", 
                    "TF_7_6_21_23S", "TF_9_11_21_23S", "TF_11_9_21_23S")


all_shared <- shared %>%
  filter(sample_id %in% samples_wanted) %>%
  as.data.frame

#all_shared <- as.data.frame(all_shared)

rownames(all_shared) <- all_shared$sample_id
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
  theme_classic() + xlim(-1, 1) + ylim(-1, 1) +
  theme(legend.key.size = unit(0.25, "cm"),
        legend.position = c(0.85, 0.9),
        legend.background = element_rect(fill="NA",
                                         color="black"),
        legend.margin = margin(t=2, r=3, b=3, l=3),
        plot.title=element_text(hjust=0.5))

ggsave("23S_D_plots/23SD_nmds_pilotv81RABRvTF.tiff", width=5, height=4)

#Beta Diversity PilotvLabRABRvTFvCVWRF
samples_wanted <- c("R27_11_3_21_23S", "R36_10_18_21_23S", "R43_11_15_21_23S", 
                    "R45_10_18_21_23S", "R45_11_15_21_23S", "R46_11_5_21_23S", 
                    "R58_10_28_21_23S", "R60_11_1_21_23S", "R60_11_15_21_23S", 
                    "R7_11_15_21_23S", "R72_11_15_21_23S", "R75_11_15_21_23S",
                    "10_5_23S", "19_23S",
                    "26_23S", "11S_23S", "11R_23S", "S1_23S", "S2_23S", "S3_23S",
                    "TF_5_25_22_23S", "TF_6_9_22_23S", "TF_6_22_22_23S", 
                    "TF_7_6_21_23S", "TF_9_11_21_23S", "TF_11_9_21_23S",
                    "CVWRF_PR_6_9_22_23S", "CVWRF_PSR_6_22_22_23S")


all_shared <- shared %>%
  filter(sample_id %in% samples_wanted) %>%
  as.data.frame

#all_shared <- as.data.frame(all_shared)

rownames(all_shared) <- all_shared$sample_id
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
  theme_classic() + xlim(-1, 1) + ylim(-1, 1) +
  theme(legend.key.size = unit(0.25, "cm"),
        legend.position = c(0.85, 0.9),
        legend.background = element_rect(fill="NA",
                                         color="black"),
        legend.margin = margin(t=2, r=3, b=3, l=3),
        plot.title=element_text(hjust=0.5))

ggsave("23S_D_plots/23SD_nmds_pilotv81RABRvTFvCVWRF.tiff", width=5, height=4)

#Beta Diversity PilotvLabRABRvTFvCVWRF
samples_wanted <- c("R27_11_3_21_23S", "R36_10_18_21_23S", "R43_11_15_21_23S", 
                    "R45_10_18_21_23S", "R45_11_15_21_23S", "R46_11_5_21_23S", 
                    "R58_10_28_21_23S", "R60_11_1_21_23S", "R60_11_15_21_23S", 
                    "R7_11_15_21_23S", "R72_11_15_21_23S", "R75_11_15_21_23S",
                    "10_5_23S", "19_23S",
                    "26_23S", "11S_23S", "11R_23S", "S1_23S", "S2_23S", "S3_23S",
                    "TF_5_25_22_23S", "TF_6_9_22_23S", "TF_6_22_22_23S", 
                    "TF_7_6_21_23S", "TF_9_11_21_23S", "TF_11_9_21_23S",
                    "CVWRF_PR_6_9_22_23S", "CVWRF_PSR_6_22_22_23S",
                    "GHR_6_15_22_23S", "GHR_5_1_22_23S")


all_shared <- shared %>%
  filter(sample_id %in% samples_wanted) %>%
  as.data.frame

#all_shared <- as.data.frame(all_shared)

rownames(all_shared) <- all_shared$sample_id
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
        legend.margin = margin(t=2, r=3, b=3, l=3),
        plot.title=element_text(hjust=0.5))

ggsave("23S_D_plots/23SD_nmds_pilotv81RABRvTFvCVWRFvGH.tiff", width=5, height=4)

# All samples nmds
#Beta Diversity 
all_shared <- shared %>%
  #filter(sample_id %in% samples_wanted) %>%
  as.data.frame

#all_shared <- as.data.frame(all_shared)

rownames(all_shared) <- all_shared$sample_id
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

write.csv(all_star,"23Scsvs/23S_beta_all.csv", row.names = FALSE)

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
        legend.margin = margin(t=1, r=3, b=3, l=3),
        plot.title=element_text(hjust=0.5))

ggsave("23S_D_plots/23SD_nmds_all.tiff", width=5, height=4)

# Temp comparison
# select for Lab samples
samples_wanted <- c("R27_11_3_21_23S", "R36_10_18_21_23S", "R43_11_15_21_23S", 
                    "R45_10_18_21_23S", "R45_11_15_21_23S", "R46_11_5_21_23S", 
                    "R58_10_28_21_23S", "R60_11_1_21_23S", "R60_11_15_21_23S", 
                    "R7_11_15_21_23S", "R72_11_15_21_23S", "R75_11_15_21_23S")

temp_shared <- shared %>%
  filter(sample_id %in% samples_wanted) %>%
  as.data.frame

#all_shared <- as.data.frame(all_shared)

rownames(temp_shared) <- temp_shared$sample_id
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

ggsave("23S_D_plots/23SD_nmds_lab_temp.tiff", width=5, height=4)

