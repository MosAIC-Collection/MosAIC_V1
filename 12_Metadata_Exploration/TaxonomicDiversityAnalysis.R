#### Script to Explore Diversity of the Collection 
library(tidyverse)

MosAIC_clean <- read_tsv("MosAIC_V1/12_Metadata_Exploration/Table_S1_MosAIC_Metadata_Revised2_Final.txt") 

colnames(MosAIC_clean)

#### How many isolates identified in males vs females 
### Function to group by and summarise MosAIC variables
GroupByAndSummarise <- function(table, var1, var2){
  table %>% 
    select(internal_id, {{var1}}, {{var2}}) %>%
    group_by({{var2}}) %>%
    summarise(number_isolates = n())
}

# Sex
sex <- GroupByAndSummarise(MosAIC_clean, gtdb_species, mosquito_sex)
# Whole Body 
tissue <- GroupByAndSummarise(MosAIC_clean, gtdb_genus, mosquito_tissue)
# Life Stage
life_stage <- GroupByAndSummarise(MosAIC_clean, gtdb_genus, mosquito_life_stage)
# Lab Field Derived
lab_field <- GroupByAndSummarise(MosAIC_clean, gtdb_genus, lab_field_derived)
# Feeding Status
feed_status <- GroupByAndSummarise(MosAIC_clean, gtdb_genus, mosquito_female_feeding_status)
# Source 
source <- GroupByAndSummarise(MosAIC_clean, gtdb_species, sample_source)
# Mosquito species 
species <- GroupByAndSummarise(MosAIC_clean, gtdb_genus, mosquito_species)
# Lab 
lab <- GroupByAndSummarise(MosAIC_clean, gtdb_genus, source_lab)
# location 
location <- GroupByAndSummarise(MosAIC_clean, gtdb_genus, sample_location)

create_bar_plot <- function(data, facet_var) {
  
  data_filt <- data %>%
    filter({{facet_var}} != "NA") %>%
    filter({{facet_var}} != "unk") 
  
  data_wrangle_filt <- data_filt %>%
    select(internal_id, gtdb_family, facet_var) %>%
    group_by(gtdb_family, facet_var) %>%
    summarise(number_isolates = n())
  
  data_wrangle_filt %>%
    ggplot() + 
    aes(x = gtdb_family, y = number_isolates) + 
    geom_bar(stat = "identity") + 
    theme_bw(base_size = 20) + 
    facet_wrap(~facet_var, scales = "fixed") + 
    coord_flip() + 
    ylab("Number of Isolates") + 
    xlab("Family")
}



#### How many isolates in male / female 
MosAIC_clean %>%
  select(internal_id, gtdb_family, mosquito_female_feeding_status) %>%
  group_by(gtdb_family, mosquito_female_feeding_status) %>%
  summarise(number_isolates = n()) %>%
  filter(mosquito_female_feeding_status != "NA") %>%
  filter(mosquito_female_feeding_status != "unk") %>%
  ggplot() + 
  aes(x = gtdb_family, y = number_isolates) + 
  geom_bar(stat = "identity") + 
  theme_bw(base_size = 20) + 
  facet_wrap(~mosquito_female_feeding_status, scales = "fixed") + 
  coord_flip() + 
  ylab("Number of Isolates") + 
  xlab("Family") 

ggsave("270924_BloodVsSugar_Family_Distributions.pdf", width = 10, height = 10)

MosAIC_clean %>%
  select(internal_id, gtdb_family, mosquito_sex) %>%
  group_by(gtdb_family, mosquito_sex) %>%
  summarise(number_isolates = n()) %>%
  filter(mosquito_sex != "NA") %>%
  filter(mosquito_sex != "unk") %>%
  ggplot() + 
  aes(x = gtdb_family, y = number_isolates) + 
  geom_bar(stat = "identity") + 
  theme_bw(base_size = 20) + 
  facet_wrap(~mosquito_sex, scales = "fixed") + 
  coord_flip() + 
  ylab("Number of Isolates") + 
  xlab("Family") 

ggsave("270924_Sex_Family_Distributions.pdf", width = 10, height = 10)

MosAIC_clean %>%
  select(internal_id, gtdb_family, mosquito_tissue) %>%
  group_by(gtdb_family, mosquito_tissue) %>%
  summarise(number_isolates = n()) %>%
  filter(mosquito_tissue != "NA") %>%
  filter(mosquito_tissue != "unk") %>%
  ggplot() + 
  aes(x = gtdb_family, y = number_isolates) + 
  geom_bar(stat = "identity") + 
  theme_bw(base_size = 20) + 
  facet_wrap(~mosquito_tissue, scales = "fixed", nrow = 1) + 
  coord_flip() + 
  ylab("Number of Isolates") + 
  xlab("Family") 

ggsave("270924_Tissue_Family_Distributions.pdf", width = 10, height = 10)

MosAIC_clean %>%
  select(internal_id, gtdb_family, mosquito_life_stage) %>%
  group_by(gtdb_family, mosquito_life_stage) %>%
  summarise(number_isolates = n()) %>%
  filter(mosquito_life_stage != "NA") %>%
  filter(mosquito_life_stage != "unknown") %>%
  ggplot() + 
  aes(x = gtdb_family, y = number_isolates) + 
  geom_bar(stat = "identity") + 
  theme_bw(base_size = 20) + 
  facet_wrap(~mosquito_life_stage, scales = "fixed", nrow = 1) + 
  coord_flip() + 
  ylab("Number of Isolates") + 
  xlab("Family") 

ggsave("270924_LifeStage_Family_Distributions.pdf", width = 10, height = 10)

MosAIC_clean %>%
  select(internal_id, gtdb_family, lab_field_derived) %>%
  group_by(gtdb_family, lab_field_derived) %>%
  summarise(number_isolates = n()) %>%
  filter(lab_field_derived != "NA") %>%
  filter(lab_field_derived != "unknown") %>%
  ggplot() + 
  aes(x = gtdb_family, y = number_isolates) + 
  geom_bar(stat = "identity") + 
  theme_bw(base_size = 20) + 
  facet_wrap(~lab_field_derived, scales = "fixed", nrow = 1) + 
  coord_flip() + 
  ylab("Number of Isolates") + 
  xlab("Family") 

ggsave("270924_lab_field_Family_Distributions.pdf", width = 10, height = 10)

MosAIC_clean %>%
  select(internal_id, gtdb_family, sample_source) %>%
  group_by(gtdb_family, sample_source) %>%
  summarise(number_isolates = n()) %>%
  #filter(sample_source != "NA") %>%
  #filter(sample_source != "unknown") %>%
  ggplot() + 
  aes(x = gtdb_family, y = number_isolates) + 
  geom_bar(stat = "identity") + 
  theme_bw(base_size = 20) + 
  facet_wrap(~sample_source, scales = "fixed", nrow = 1) + 
  coord_flip() + 
  ylab("Number of Isolates") + 
  xlab("Family") 


ggsave("270924_Source_Family_Distributions.pdf", width = 10, height = 10)

### Percentages
MosAIC_clean %>%
  select(internal_id, gtdb_family, lab_field_derived) %>%
  group_by(gtdb_family, lab_field_derived) %>%
  summarise(number_isolates = n()) %>%
  filter(gtdb_family == "Bacillaceae" | gtdb_family == "Bacillaceae_G" | gtdb_family == "Bacillaceae_H" ) %>%
  ungroup() %>%
  group_by(lab_field_derived) %>%
  summarise(sum_isolates = sum(number_isolates)) %>%
  mutate(percent = sum_isolates / sum(sum_isolates))

MosAIC_clean %>%
  select(internal_id, gtdb_family, lab_field_derived) %>%
  group_by(gtdb_family, lab_field_derived) %>%
  summarise(number_isolates = n()) %>%
  filter(gtdb_family == "Enterobacteriaceae") %>%
  ungroup() %>%
  group_by(lab_field_derived) %>%
  summarise(sum_isolates = sum(number_isolates)) %>%
  mutate(percent = sum_isolates / sum(sum_isolates))

MosAIC_clean %>%
  select(internal_id, mosquito_tissue) %>%
  group_by(mosquito_tissue) %>%
  summarise(number_isolates = n()) %>%
  mutate(percent = number_isolates / sum(number_isolates))

### Looking at Enterobactericeae 
MosAIC_clean %>%
  filter(gtdb_family == "Enterobacteriaceae") %>%
  select(internal_id, gtdb_genus, sample_source) %>%
  group_by(gtdb_genus, sample_source) %>%
  summarise(number_isolates = n()) %>%
  filter(sample_source != "NA") %>%
  filter(sample_source != "unknown") %>%
  ggplot() + 
  aes(x = gtdb_genus, y = number_isolates) + 
  geom_bar(stat = "identity") + 
  theme_bw(base_size = 25) + 
  facet_wrap(~sample_source, scales = "fixed", nrow = 1) + 
  coord_flip() + 
  ylab("Number of Isolates") + 
  xlab("Genus") + 
  theme(axis.text.y = element_text(face = "italic"))



MosAIC_clean %>%
  filter(gtdb_family == "Enterobacteriaceae") %>%
  select(internal_id, gtdb_genus, lab_field_derived) %>%
  group_by(gtdb_genus, lab_field_derived) %>%
  summarise(number_isolates = n()) %>%
  filter(lab_field_derived != "NA") %>%
  filter(lab_field_derived != "unknown") %>%
  ggplot() + 
  aes(x = gtdb_genus, y = number_isolates) + 
  geom_bar(stat = "identity") + 
  theme_bw(base_size = 25) + 
  facet_wrap(~lab_field_derived, scales = "fixed", nrow = 1) + 
  coord_flip() + 
  ylab("Number of Isolates") + 
  xlab("Genus") + 
  theme(axis.text.y = element_text(face = "italic"))

ggsave("270924_LabField_Enterobactericeae.pdf", width = 10, height = 10)


MosAIC_clean %>%
  filter(gtdb_family == "Enterobacteriaceae") %>%
  select(internal_id, gtdb_genus, mosquito_sex) %>%
  group_by(gtdb_genus, mosquito_sex) %>%
  summarise(number_isolates = n()) %>%
  filter(mosquito_sex != "NA") %>%
  filter(mosquito_sex != "unknown") %>%
  ggplot() + 
  aes(x = gtdb_genus, y = number_isolates) + 
  geom_bar(stat = "identity") + 
  theme_bw(base_size = 25) + 
  facet_wrap(~mosquito_sex, scales = "fixed", nrow = 1) + 
  coord_flip() + 
  ylab("Number of Isolates") + 
  xlab("Genus") + 
  theme(axis.text.y = element_text(face = "italic"))

ggsave("270924_Sex_Enterobactericeae.pdf", width = 10, height = 10)


MosAIC_clean %>%
  filter(gtdb_family == "Enterobacteriaceae") %>%
  select(internal_id, gtdb_genus, tissue) %>%
  group_by(gtdb_genus, tissue) %>%
  summarise(number_isolates = n()) %>%
  filter(tissue != "NA") %>%
  filter(tissue != "unknown") %>%
  ggplot() + 
  aes(x = gtdb_genus, y = number_isolates) + 
  geom_bar(stat = "identity") + 
  theme_bw(base_size = 20) + 
  facet_wrap(~tissue, scales = "fixed", nrow = 1) + 
  coord_flip() + 
  ylab("Number of Isolates") + 
  xlab("Family")

MosAIC_clean %>%
  filter(gtdb_family == "Enterobacteriaceae") %>%
  select(internal_id, gtdb_genus, feeding_status) %>%
  group_by(gtdb_genus, feeding_status) %>%
  summarise(number_isolates = n()) %>%
  filter(feeding_status != "NA") %>%
  filter(feeding_status != "unknown") %>%
  ggplot() + 
  aes(x = gtdb_genus, y = number_isolates) + 
  geom_bar(stat = "identity") + 
  theme_bw(base_size = 20) + 
  facet_wrap(~feeding_status, scales = "fixed", nrow = 1) + 
  coord_flip() + 
  ylab("Number of Isolates") + 
  xlab("Family")

MosAIC_clean %>%
  filter(gtdb_family == "Enterobacteriaceae") %>%
  select(internal_id, gtdb_genus, lab_field_derived) %>%
  group_by(gtdb_genus, lab_field_derived) %>%
  summarise(number_isolates = n()) %>%
  filter(gtdb_genus == "Pantoea") %>%
  mutate(percent = number_isolates / sum(number_isolates))

MosAIC_clean %>%
  filter(gtdb_family == "Enterobacteriaceae") %>%
  select(internal_id, gtdb_genus, lab_field_derived) %>%
  group_by(gtdb_genus, lab_field_derived) %>%
  summarise(number_isolates = n()) %>%
  filter(gtdb_genus == "Enterobacter") %>%
  mutate(percent = number_isolates / sum(number_isolates))

MosAIC_clean %>%
  filter(gtdb_family == "Enterobacteriaceae") %>%
  select(internal_id, gtdb_genus, mosquito_sex) %>%
  group_by(mosquito_sex) %>%
  summarise(number_isolates = n()) %>%
  #filter(gtdb_genus == "Enterobacter") %>%
  mutate(percent = number_isolates / sum(number_isolates))
