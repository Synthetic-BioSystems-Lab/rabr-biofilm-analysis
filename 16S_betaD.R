library(vegan)
library(ggtext)
library(tidyverse)
library(readxl)

setwd("~/Miller Lab/Rscripts_PilotRABR")

gather_metadata <- function(target, t2, t3, page) {
  read_excel("16Spilot/16S_metadata.xlsx", sheet=page) %>%
    rename_all(tolower) %>%
    mutate(sample = str_replace_all(sample, "-", "_")) %>%
    select(sample, target, t2, t3)
}

#Beta Diversity 1

dist_matrix <- read_tsv("16Spilot/final.opti_mcc.braycurtis.0.03.square.ave.dist", 
                        skip=1, col_names = FALSE) %>%
  as.data.frame

rownames(dist_matrix) <- dist_matrix$X1
dist_matrix <- dist_matrix[, -1]
dist_matrix <- as.matrix(dist_matrix)

dist_nmds <- metaMDS(dist_matrix) # also insufficient data :(

#Beta Diversity 2

all_shared <- read_tsv("16Spilot/final.opti_mcc.shared") %>%
  select(Group, starts_with("Otu")) %>%
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  mutate(total=sum(value)) %>%
  # arrange(total) %>%
  # group_by(name) %>%
  # summarize(total=sum(value)) %>% filter(total ==0) %>%
  ungroup() %>%
  select(-total) %>%
  pivot_wider(Group)
  #filter(Group != "C1_16S" & Group != "C2_16S")

for(i in 1:nrow(all_shared)) {
  row <- all_shared[i,]
  row[1] <- gsub("_16S", "", row[1])
  all_shared[i,1] <- row[1]
}


all_shared <- as.data.frame(all_shared)

rownames(all_shared) <- all_shared$Group
all_shared <- all_shared[, -1]
all_shared <- as.matrix(all_shared)


all_dist <- avgdist(all_shared, dmethod="bray", iterations=1000, sample=20073)
all_nmds <- metaMDS(all_dist) #low stress > insufficient data?

all_scores <- scores(all_nmds) %>%
  as_tibble(rownames = "Group")

all_metadata <- read_excel("metadata_start.xlsx", sheet=1) %>%
  rename_all(tolower) %>%
  mutate(sample = str_replace_all(sample, "-", "_")) %>%
  select(sample, date)

all_metadata_nmds <- inner_join(all_metadata, all_scores, by=c('sample'='Group')) %>%
  mutate(date = factor(date,
                       levels=c("10_5_23",
                                "10_12_23",
                                "10_19_23",
                                "10_26_23",
                                "11_2_23")))

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
  scale_color_manual(name=NULL, breaks=c("10_5_23", "10_12_23", "10_19_23", "10_26_23", "11_2_23"),
                     labels=c("10_5_23", "10_12_23", "10_19_23", "10_26_23", "11_2_23"),
                     values=c(Five_color, Twelve_color, Nineteen_color, TwentySix_color, Two_color)) +
  scale_fill_manual(name=NULL, breaks=c("10_5_23", "10_12_23", "10_19_23", "10_26_23", "11_2_23"),
                    labels=c("10_5_23", "10_12_23", "10_19_23", "10_26_23", "11_2_23"),
                    values=c(Five_color, Twelve_color, Nineteen_color, TwentySix_color, Two_color)) +
  theme_classic() + xlim(-0.75, 0.75) + ylim(-0.75, 0.75) +
  theme(legend.key.size = unit(0.25, "cm"),
        legend.position = c(0.85, 0.9),
        legend.background = element_rect(fill="NA",
                                         color="black"),
        legend.margin = margin(t=-2, r=3, b=3, l=3),
        plot.title=element_text(hjust=0.5))

ggsave("16Spilot/16Splots/16S_nmds_all_date.tiff", width=5, height=4)

#Beta Diversity PilotvLabRABR
samples_wanted <- c("11_R45_10_18_21_16S", "12_R60_11_15_21_16S", "13_R27_11_3_21_16S",
                    "14_R46_11_5_21_16S", "15_R7_11_15_21_16S", "16_R45_11_15_21",
                    "17_R36_10_28_21_16S", "18_R60_11_21_16S", "19_R43_11_15_21_16S",
                    "20_R58_10_28_21_16S", "21_R72_11_15_21_16S", "R1_S2_15_10_4_21",
                    "R1_S2_15_10_12_21", "R3_S1_58_10_4_21", "10_5_16S", "19_16S",
                    "26_16S", "11S_16S", "11R_16S", "S1_16S", "S2_16S", "S3_16S")


all_shared <- read_tsv("16Spilot/final.opti_mcc.shared") %>%
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
        legend.margin = margin(t=-2, r=3, b=3, l=3),
        plot.title=element_text(hjust=0.5))

ggsave("16Spilot/16Splots/16S_nmds_pilotv81RABR.tiff", width=5, height=4)

#Beta Diversity PilotvLabRABRvTF
samples_wanted <- c("11_R45_10_18_21_16S", "12_R60_11_15_21_16S", "13_R27_11_3_21_16S",
                    "14_R46_11_5_21_16S", "15_R7_11_15_21_16S", "16_R45_11_15_21",
                    "17_R36_10_28_21_16S", "18_R60_11_21_16S", "19_R43_11_15_21_16S",
                    "20_R58_10_28_21_16S", "21_R72_11_15_21_16S", "R1_S2_15_10_4_21",
                    "R1_S2_15_10_12_21", "R3_S1_58_10_4_21", "10_5_16S", "19_16S",
                    "26_16S", "11S_16S", "11R_16S", "S1_16S", "S2_16S", "S3_16S",
                    "3_TF_5_25_22_16S", "7_TF_6_9_22_16S", "10_TF_6_22_22_16S",
                    "TF_7_6_21", "TF_9_11_21", "TF_11_9_21_R1")


all_shared <- read_tsv("16Spilot/final.opti_mcc.shared") %>%
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

ggsave("16Spilot/16Splots/16S_nmds_pilotv81RABRvTF.tiff", width=5, height=4)

#Beta Diversity PilotvLabRABRvTFvCVWRF
samples_wanted <- c("11_R45_10_18_21_16S", "12_R60_11_15_21_16S", "13_R27_11_3_21_16S",
                    "14_R46_11_5_21_16S", "15_R7_11_15_21_16S", "16_R45_11_15_21",
                    "17_R36_10_28_21_16S", "18_R60_11_21_16S", "19_R43_11_15_21_16S",
                    "20_R58_10_28_21_16S", "21_R72_11_15_21_16S", "R1_S2_15_10_4_21",
                    "R1_S2_15_10_12_21", "R3_S1_58_10_4_21", "10_5_16S", "19_16S",
                    "26_16S", "11S_16S", "11R_16S", "S1_16S", "S2_16S", "S3_16S",
                    "3_TF_5_25_22_16S", "7_TF_6_9_22_16S", "10_TF_6_22_22_16S",
                    "TF_7_6_21", "TF_9_11_21", "TF_11_9_21_R1",
                    "1_CVWRF_PR_6_22_22_16S", "4_CVWRF_PSR_2_22_22_16S")


all_shared <- read_tsv("16Spilot/final.opti_mcc.shared") %>%
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

ggsave("16Spilot/16Splots/16S_nmds_pilotv81RABRvTFvCVWRF.tiff", width=5, height=4)

#Beta Diversity PilotvLabRABRvTFvCVWRF
samples_wanted <- c("11_R45_10_18_21_16S", "12_R60_11_15_21_16S", "13_R27_11_3_21_16S",
                    "14_R46_11_5_21_16S", "15_R7_11_15_21_16S", "16_R45_11_15_21",
                    "17_R36_10_28_21_16S", "18_R60_11_21_16S", "19_R43_11_15_21_16S",
                    "20_R58_10_28_21_16S", "21_R72_11_15_21_16S", "R1_S2_15_10_4_21",
                    "R1_S2_15_10_12_21", "R3_S1_58_10_4_21", "10_5_16S", "19_16S",
                    "26_16S", "11S_16S", "11R_16S", "S1_16S", "S2_16S", "S3_16S",
                    "3_TF_5_25_22_16S", "7_TF_6_9_22_16S", "10_TF_6_22_22_16S",
                    "TF_7_6_21", "TF_9_11_21", "TF_11_9_21_R1",
                    "1_CVWRF_PR_6_22_22_16S", "4_CVWRF_PSR_2_22_22_16S",
                    "2_GHR_6_15_22_16S", "6_GHR_5_1_22_16S")


all_shared <- read_tsv("16Spilot/final.opti_mcc.shared") %>%
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

ggsave("16Spilot/16Splots/16S_nmds_pilotv81RABRvTFvCVWRFvGH.tiff", width=5, height=4)

# All samples nmds
#Beta Diversity PilotvLabRABR
all_shared <- read_tsv("16Spilot/final.opti_mcc.shared") %>%
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
        legend.margin = margin(t=-2, r=3, b=3, l=3),
        plot.title=element_text(hjust=0.5))

ggsave("16Spilot/16Splots/16S_nmds_all.tiff", width=5, height=4)

# Temp comparison
# select for Lab samples
samples_wanted <- c("11_R45_10_18_21_16S", "12_R60_11_15_21_16S", "13_R27_11_3_21_16S",
                    "14_R46_11_5_21_16S", "15_R7_11_15_21_16S", "16_R45_11_15_21",
                    "17_R36_10_28_21_16S", "18_R60_11_21_16S", "19_R43_11_15_21_16S",
                    "20_R58_10_28_21_16S", "21_R72_11_15_21_16S", "R1_S2_15_10_4_21",
                    "R1_S2_15_10_12_21", "R3_S1_58_10_4_21")

temp_shared <- read_tsv("16Spilot/final.opti_mcc.shared") %>%
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

temp_dist <- avgdist(temp_shared, dmethod="bray", iterations=1000, sample=20073)
temp_nmds <- metaMDS(temp_dist)

temp_scores <- scores(temp_nmds) %>%
  as_tibble(rownames = "Group") 

temp_metadata <- gather_metadata("temp", "date", "label", 2)

temp_metadata_nmds <- inner_join(temp_metadata, temp_scores, by=c('sample'='Group')) %>%
  mutate(temp = factor(temp,
                              levels=c("10",
                                       "20",
                                       "30")))

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
                              "30"),
                     labels=c("10 C",
                              "20 C",
                              "30 C"),
                     values=c(T1_color, T2_color, T3_color)) +
  scale_fill_manual(name=NULL, 
                    breaks=c("10",
                             "20",
                             "30"),
                    labels=c("10 C",
                             "20 C",
                             "30 C"),
                    values=c(T1_color, T2_color, T3_color)) +
  theme_classic() + xlim(-0.75, 0.75) + ylim(-0.75, 0.75) +
  theme(legend.key.size = unit(0.35, "cm"),
        legend.position = c(0.85, 0.9),
        legend.background = element_rect(fill="NA",
                                         color="black"),
        legend.margin = margin(t=3, r=3, b=3, l=3),
        plot.title=element_text(hjust=0.5))

ggsave("16Spilot/16Splots/16S_nmds_lab_temp.tiff", width=5, height=4)
