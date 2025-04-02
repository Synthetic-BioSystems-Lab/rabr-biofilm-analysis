library(tidyverse) 
library(ggplot2)
library(ggtext)
library(readxl)


#relative abund stacked bar graph
#relative abund heat map
#relative abund bar graph

setwd("~/Miller Lab/Rscripts_PilotRABR")

final_shared <- "18Spilotv2/final_p0.asv.shared"
final_tax <- "18Spilotv2/final_p0.asv.0.03.cons.taxonomy"

gather_metadata <- function(target, t2, t3, page) {
  read_excel("18Spilotv2/18S_metadata.xlsx", sheet=page) %>%
    rename_all(tolower) %>%
    mutate(sample = str_replace_all(sample, "-", "_")) %>%
    select(sample, target, t2, t3)
}

#import asv counts
asv_counts <- read_tsv(final_shared) %>%
  select(-label, -numASVs) %>%
  pivot_longer(-Group, names_to = "otu", values_to = "count") %>%
  rename(sample_id = Group)

# import taxonomy
taxonomy <- read_tsv(final_tax) %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""), taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), sep = ";")

# import metadata
all_metadata <- gather_metadata("section", "date", "label", 1) 

for(i in 2:6) {
  all_metadata <- dplyr::bind_rows(all_metadata, gather_metadata("section", "date", "label", i))
}

all_metadata %>%
  mutate(section = factor(section,
                          levels=c("pilot", "CVWRF", "81RABR", "TF", "GHR", "control")))

#merge asv and taxon
composite <- inner_join(asv_counts, taxonomy, by="otu")

#remove empty rows, add relative abundance, pivot for readability
trimmed_composite <- composite[!(composite$count=="0"),] %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols = c("kingdom", "phylum", "class", "order", "family", "genus"), names_to = "level", values_to = "taxon")

asv_rel_abund <- inner_join(trimmed_composite, all_metadata,  by=c('sample_id'='sample'))

taxon_rel_abund <- asv_rel_abund %>%
  filter(level=="class") %>%
  group_by(section, sample_id, taxon, date, label) %>%
  summarize(rel_abund = 100*sum(rel_abund), .groups = "drop") %>%
  # group_by(sample_id, taxon, section) %>%
  # summarize(mean_rel_abund = 100*mean(rel_abund), .groups = "drop") %>%
  mutate(taxon = str_replace(taxon, "(.*)_unclassified", "Unclassified *\\1*"),
         taxon = str_replace(taxon, "^(\\S*)$", "*\\1*")) 

#set aside other
taxon_pool <- taxon_rel_abund %>%
  group_by(section, taxon, rel_abund) %>%
  summarize(mean=mean(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarize(pool = max(mean) < 18.8, 
            mean = mean(mean),
            .groups="drop")

sample_order <- all_metadata %>%
  #filter(taxon == "*Brevundimonas*") %>%
  arrange(date) %>%
  mutate(order = 1:nrow(.)) %>%
  select(sample, order, date)

# RvsS <- taxon_rel_abund %>%
#   filter(sample_id %in% c("11S_18S", "11R_18S")) %>%
#   group_by(sample_id) %>%
#   pivot_wider(names_from = sample_id, values_from = mean_rel_abund) %>%
#   data.table::setnames("11S_18S",'Srel_abund') %>%
#   data.table::setnames("11R_18S",'Rrel_abund') %>%
#   mutate(diff = (Srel_abund - Rrel_abund) / Rrel_abund) %>%
#   ungroup() 

pretty <- c("pilot" = "Pilot-scale RABRs",
            "81RABR" = "Lab-scale RABRs",
            "control" = "Control",
            "CVWRF" = "CVWRF",
            "GHR" = "GHR",
            "TF" = "Trickling Filter")
#<br>

#assemble others and make RA stacked bar plot
prep <- inner_join(taxon_rel_abund, taxon_pool, by="taxon") %>%
  #filter(sample_id %in% c("S1_18S", "S2_18S", "S3_18S")) %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(sample_id, label, section, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups = "drop") %>%
  # mutate(taxon = factor(taxon),
  #        taxon = fct_reorder(taxon, mean, .desc=TRUE),
  #        taxon = fct_shift(taxon, n=1)) %>%
  inner_join(., sample_order, by=c("sample_id"="sample")) %>%
  mutate(label = factor(label),
         label= fct_reorder(label, order))

prep %>%
  ggplot(aes(x = label, y = rel_abund, fill = taxon)) +
  geom_col() + 
  #scale_fill_manual(name = NULL, values = c(brewer.pal(6, "Dark2"), "gray")) +
  scale_fill_discrete(name=NULL) +
  scale_y_continuous(expand=c(0,0)) +
  facet_grid(~section, scale="free_x", space="free", 
             labeller = labeller(section=pretty)) +
  labs(title="Class Relative Abundance of RABRs",
       x = NULL,
       y = "Mean Relative Abundance (%)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_markdown(), 
        axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0),
        legend.key.size = unit(10, "pt"),
        strip.background = element_blank(),
        strip.text = element_markdown())

ggsave("18Spilotv2/18Splots/18S_stacked_bar_class.tiff", width=9, height=4)

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
  labs(title="Class Relative Abundance of RABRs",
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

ggsave("18Spilotv2/18Splots/18S_heat_map_class.tiff", width=9, height=4)

