### Load Packages 
library(tidyverse)
library(ggsankey)
library(janitor)

### Make Sankey Plots
MosAIC_tbl <- read_tsv("MosAIC_V1/04_Sankey_Diagram/Table S1 MosAIC_Metadata_Revised3_LB.txt") %>%
  clean_names()

sources <- c("unknown", "water", "cow manure", "non-mosquito Diptera", "mosquito", "larval water")
species <- c("unknown", "Toxorhynchites amboinensis", "Ochlerotatus trivittatus", "2_not applicable", 
             "Culex pipiens", "Culex erraticus", "Anopheles quadtimaculatus", "Anopheles punctipennis", 
             "Anopheles gambiae", "Aedes triseriatus", "Aedes atropalpus", "Aedes albopictus", "Aedes aegypti")

colours = c('1_Unknown' = "light grey", 
            'not applicable' = "light grey", 
            'non-mosquito diptera' = "blue", 
            '4_mosquito' = 'black', 
            '3_MLW' = "#0F5DA7", 
            '2_NMDE' = '#67993E', 
            '1_NMDH' = '#67993E', 
            'lab' = '#C9C89A',
            'field' = '#1A4926',
            '3_field' = '#1A4926',
            '2_field' = '#1A4926',
            '1_field' = '#1A4926',
            "Aedes aegypti" = "#867730", 
            "Aedes albopictus" = "#867730", 
            "Aedes atropalpus" = "#867730",
            "Aedes triseriatus" = "#867730",
            "Aedes trivittatus" = "#867730", 
            "Anopheles gambiae" = "#7C8562", 
            "Anopheles punctipennis" = "#7C8562",
            "Anopheles quadrimaculatus" = "#7C8562", 
            "Culex pipiens" = "#95748B", 
            "Culex erraticus" = "#95748B", 
            "Ochlerotatus trivittatus" = "#B59CDE", 
            "Toxorhynchites amboinensis" = "#4B4F36", 
            "larvae" = "#E99860", 
            "eggs" = "#987628", 
            "adults" = "#DFAE3D", 
            "male" = "#E1C0BC", 
            "female" = "#CAB994", 
            "non blood-fed" = "#DCEBD1", 
            "blood-fed" = "#F2B994", 
            "whole body" = "#ADC899", 
            "midgut" = "#CBC2A8")

### Clean Data 
MosAIC_tbl_clean <- MosAIC_tbl %>%
  mutate(sample_source = if_else(sample_source == "mosquito_larval_water", "3_MLW", sample_source)) %>%
  mutate(sample_source = if_else(sample_source == "non_mosquito_Diptera_host", "1_NMDH", sample_source)) %>%
  mutate(sample_source = if_else(sample_source == "non_mosquito_Diptera_environment", "2_NMDE", sample_source)) %>%
  mutate(sample_source = if_else(sample_source == "mosquito", "4_mosquito", sample_source)) %>%
  mutate(mosquito_species = if_else(mosquito_species == "Aedes_aegypti", "Aedes aegypti", mosquito_species)) %>%
  mutate(mosquito_species = if_else(mosquito_species == "Aedes_albopictus", "Aedes albopictus", mosquito_species)) %>%
  mutate(mosquito_species = if_else(mosquito_species == "Aedes_atropalpus", "Aedes atropalpus", mosquito_species)) %>%
  mutate(mosquito_species = if_else(mosquito_species == "Aedes_triseriatus", "Aedes triseriatus", mosquito_species)) %>%
  mutate(mosquito_species = if_else(mosquito_species == "Aedes_trivittatus", "Aedes trivittatus", mosquito_species)) %>%
  mutate(mosquito_species = if_else(mosquito_species == "Anopheles_gambiae", "Anopheles gambiae", mosquito_species)) %>%
  mutate(mosquito_species = if_else(mosquito_species == "Anopheles_punctipennis", "Anopheles punctipennis", mosquito_species)) %>%
  mutate(mosquito_species = if_else(mosquito_species == "Anopheles_quadrimaculatus", "Anopheles quadrimaculatus", mosquito_species)) %>%
  mutate(mosquito_species = if_else(mosquito_species == "Culex_erraticus", "Culex erraticus", mosquito_species)) %>%
  mutate(mosquito_species = if_else(mosquito_species == "Culex_pipiens", "Culex pipiens", mosquito_species)) %>%
  mutate(mosquito_species = if_else(mosquito_species == "Toxorhynchites_amboinensis", "Toxorhynchites amboinensis", mosquito_species)) %>%
  #mutate(lab_field_derived = if_else(sample_source == "3_MLW", "3_lab", lab_field_derived)) %>%
  #mutate(lab_field_derived = if_else(sample_source == "3_MLW", "3_field", lab_field_derived)) %>%
  mutate(lab_field_derived = if_else(sample_source == "2_NMDE", "2_field", lab_field_derived)) %>%
  mutate(lab_field_derived = if_else(sample_source == "2_NMDE", "2_field", lab_field_derived)) %>%
  mutate(lab_field_derived = if_else(sample_source == "1_NMDH", "1_field", lab_field_derived)) %>%
  mutate(lab_field_derived = if_else(sample_source == "1_NMDH", "1_field", lab_field_derived)) %>%
  mutate(mosquito_species = if_else(mosquito_species == "unk", "1_Unknown", mosquito_species)) %>%
  mutate(mosquito_sex = if_else(mosquito_sex == "unk", "1_Unknown", mosquito_sex)) %>%
  mutate(mosquito_life_stage = if_else(mosquito_life_stage == "unk", "1_Unknown", mosquito_life_stage)) %>%
  mutate(mosquito_female_feeding_status = if_else(mosquito_species == "unk", "1_Unknown", mosquito_female_feeding_status)) %>%
  mutate(mosquito_tissue = if_else(mosquito_tissue == "unk", "1_Unknown", mosquito_tissue)) %>%
  mutate(mosquito_female_feeding_status = if_else(mosquito_female_feeding_status == "unk", "1_Unknown", mosquito_female_feeding_status)) %>%
  mutate(mosquito_tissue = if_else(mosquito_tissue == "whole_body", "whole body", mosquito_tissue)) %>%
  mutate(mosquito_female_feeding_status = if_else(mosquito_female_feeding_status == "non_blood_fed", "non blood-fed", mosquito_female_feeding_status)) %>%
  mutate(mosquito_female_feeding_status = if_else(mosquito_female_feeding_status == "blood_fed", "blood-fed", mosquito_female_feeding_status)) %>%
  rename("Source" = sample_source, 
         "Mosquito species" = mosquito_species, 
         "Sex" = mosquito_sex, 
         "Life stage" = mosquito_life_stage, 
         "Feeding status" = mosquito_female_feeding_status, 
         "Lab / Field Derived" = lab_field_derived, 
         "Mosquito Tissue" = mosquito_tissue)

MosAIC_Sankey_1 <- MosAIC_tbl_clean %>%
  select(Source, `Lab / Field Derived`, `Mosquito species`, `Life stage`, `Mosquito Tissue`, Sex, `Feeding status`) %>%
  arrange(Source, ) %>%
  make_long(Source, `Lab / Field Derived`, `Mosquito species`, `Life stage`, `Mosquito Tissue`, Sex, `Feeding status`) 


# Function to make Sankey, use on different DFs
MakeItSankey <- function(table){
  pl <- ggplot(table, aes(x = x,                        
                          next_x = next_x,                                     
                          node = node,
                          next_node = next_node,        
                          fill = factor(node),
                          label = node))             
  
  pl <- pl + geom_sankey(flow.alpha = 0.7,         
                         node.color = "black",     
                         show.legend = TRUE, 
                         na.rm = TRUE, 
                         #smooth = 25,  
                         width =  0.1)       
  
  pl <- pl + geom_sankey_label(size = 4.5, 
                               color = "black", 
                               fill = "white", 
                               type = "sankey", 
                               na.rm = TRUE)
  
  pl <- pl + theme_bw(base_size = 20)
  pl <- pl + theme(legend.position = 'none')
  pl <- pl + theme(axis.title = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks = element_blank())
  
  pl <- pl + scale_fill_manual(values = colours)
  pl
}

MakeItSankey(MosAIC_Sankey_1)
ggsave(file = "Figure1_Mosquito_Isolation_sources_201023.pdf", width = 10, height = 11)

# MISC
MosAIC_tbl_clean2 <- MosAIC_tbl %>%
  select(source, mosquito_species, sex, life_stage, feeding_status) %>% 
  mutate(source = str_replace(source, "non-mosquito diptera", "2_non-mosquito Diptera")) %>%
  mutate(life_stage = str_replace(life_stage, "larva", "larvae")) %>%
  mutate(life_stage = str_replace(life_stage, "egg", "eggs")) %>%
  mutate(life_stage = str_replace(life_stage, "adult", "adults")) %>%
  mutate(feeding_status = str_replace(feeding_status, "sugar", "sugar-fed")) %>%
  mutate(feeding_status = str_replace(feeding_status, "blood", "blood-fed")) %>%
  replace(is.na(.), "1_unknown") %>%
  #mutate(mosquito_species = if_else(mosquito_species == "Aedes aegypti" | 
  #mosquito_species == "Aedes albopictus" |
  #mosquito_species == "Aedes atropalpus" |
  #mosquito_species == "Aedes triseriatus", "Aedes", mosquito_species)) %>%
  #mutate(mosquito_species = if_else(mosquito_species == "Anopheles gambiae" | 
  #mosquito_species == "Anopheles punctipennis" |
  #mosquito_species == "Anopheles quadrimaculatus", "Anopheles", mosquito_species)) %>%
  #mutate(mosquito_species = if_else(mosquito_species == "Culex erraticus" | 
  #mosquito_species == "Culex pipiens" |
  #mosquito_species == "Culex pippiens", "Culex",  mosquito_species)) %>%
  mutate(source = if_else(source == "stable fly", "2_non-mosquito Diptera", source)) %>%
  mutate(source = if_else(life_stage == "larvae", "mosquito", source)) %>%
  mutate(source = if_else(life_stage == "adults", "mosquito", source)) %>%
  mutate(source = if_else(life_stage == "eggs", "mosquito", source)) %>%
  mutate(source = if_else(sex == "male" | sex == "female", "mosquito", source)) %>%
  mutate(source = if_else(source == "water" | source == "cow manure", "environmental", source)) %>%
  mutate(source = if_else(mosquito_species == "Aedes aegypti" 
                          | mosquito_species == "Aedes albopictus"
                          | mosquito_species == "Aedes atropalpus"
                          | mosquito_species == "Aedes triseriatus"
                          | mosquito_species == "Anopheles gambiae"
                          | mosquito_species == "Anopheles punctipennis"
                          | mosquito_species == "Anopheles quadrimaculatus"
                          | mosquito_species == "Culex erraticus"
                          | mosquito_species == "Culex pipiens"
                          | mosquito_species == "Culex pippiens"
                          | mosquito_species == "Ochlerotatus trivittatus"
                          | mosquito_species == "Toxorhynchites amboinensis", "mosquito", source)) %>%
  mutate(mosquito_species = if_else(source == "2_non-mosquito Diptera", "2_not applicable", mosquito_species)) %>%
  mutate(mosquito_species = if_else(source == "water", "1_unknown", mosquito_species)) %>%
  mutate(mosquito_species = if_else(mosquito_species == "Culex pippiens", "Culex pipiens", mosquito_species)) %>%
  mutate(mosquito_species = if_else(source == "environmental", "2_not applicable", mosquito_species)) %>%
  mutate(life_stage = if_else(source == "larval water", "2_not applicable", life_stage)) %>%
  mutate(life_stage = if_else(source == "2_non-mosquito Diptera", "2_not applicable", life_stage)) %>%
  mutate(life_stage = if_else(mosquito_species == "2_not applicable", NA, life_stage)) %>%
  mutate(sex = if_else(life_stage == "eggs", "1_unknown", sex)) %>%
  mutate(sex = if_else(feeding_status == "blood-fed", "female", sex)) %>%
  #mutate(sex = if_else(mosquito_species == "Not Applicable", NA, sex)) %>%
  mutate(sex = if_else(life_stage == "2_not applicable", NA, sex)) %>%
  mutate(feeding_status = if_else(sex == "male", "sugar-fed", feeding_status)) %>%
  mutate(feeding_status = if_else(sex == "female" & feeding_status != "blood-fed", "sugar-fed", feeding_status)) %>%
  mutate(feeding_status = if_else(sex == "2_not applicable", NA, feeding_status)) %>%
  #mutate(feeding_status = if_else(life_stage == "egg", NA, feeding_status)) %>%
  #mutate(feeding_status = if_else(life_stage == "larva", NA, feeding_status)) %>%
  rename("Source" = source, 
         "Mosquito species" = mosquito_species, 
         "Sex" = sex, 
         "Life stage" = life_stage, 
         "Feeding status" = feeding_status)

MosAIC_formatted <- MosAIC_tbl_clean %>%
  arrange(Source, ) %>%
  make_long(Source, `Mosquito species`, `Life stage`, Sex , `Feeding status`)




pl <- ggplot(MosAIC_formatted, aes(x = x,                        
                                   next_x = next_x,                                     
                                   node = node,
                                   next_node = next_node,        
                                   fill = factor(node),
                                   label = node))             

pl <- pl + geom_sankey(flow.alpha = 0.7,         
                       node.color = "black",     
                       show.legend = TRUE, 
                       na.rm = TRUE, 
                       #smooth = 25,  
                       width =  0.1)       

pl <- pl + geom_sankey_label(size = 4.5, 
                             color = "black", 
                             fill = "white", 
                             type = "sankey", 
                             na.rm = TRUE)

pl <- pl + theme_bw(base_size = 20)
pl <- pl + theme(legend.position = 'none')
pl <- pl + theme(axis.title = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks = element_blank())

pl <- pl + scale_fill_manual(values = c('1_unknown' = "light grey", 
                                        'not applicable' = "light grey", 
                                        'non-mosquito diptera' = "blue", 
                                        'mosquito' = 'black', 
                                        'larval water' = "#0F5DA7", 
                                        'environmental' = '#67993E', 
                                        "Aedes aegypti" = "#867730", 
                                        "Aedes albopictus" = "#867730", 
                                        "Aedes atropalpus" = "#867730",
                                        "Aedes triseriatus" = "#867730", 
                                        "Anopheles gambiae" = "#7C8562", 
                                        "Anopheles punctipennis" = "#7C8562",
                                        "Anopheles quadrimaculatus" = "#7C8562", 
                                        "Culex pipiens" = "#95748B", 
                                        "Culex erraticus" = "#95748B", 
                                        "Ochlerotatus trivittatus" = "#B59CDE", 
                                        "Toxorhynchites amboinensis" = "#4B4F36", 
                                        "larvae" = "#E99860", 
                                        "eggs" = "#987628", 
                                        "adults" = "#DFAE3D", 
                                        "male" = "#E1C0BC", 
                                        "female" = "#CAB994", 
                                        "sugar-fed" = "#DCEBD1", 
                                        "blood-fed" = "#F2B994"))
pl

ggsave(pl, file = "Figure1_Mosquito_Isolation_sources.pdf", width = 10, height = 10)
ggsave(pl, file = "Figure1_Mosquito_Isolation_sources.png", width = 10, height = 10)

MosAIC_Sankey_2 <- MosAIC_tbl_clean %>%
  filter(Source == "mosquito larval water") %>%
  select(Source, lab_field_derived) %>%
  arrange(Source, ) %>%
  make_long(Source, lab_field_derived)

MosAIC_Sankey_3 <- MosAIC_tbl_clean %>%
  filter(Source == "non_mosquito_Diptera_host") %>%
  select(Source, lab_field_derived, `Mosquito species`, `Life stage`, mosquito_tissue, Sex, `Feeding status`) %>%
  arrange(Source, ) %>%
  make_long(Source, lab_field_derived, `Mosquito species`, `Life stage`, mosquito_tissue, Sex, `Feeding status`) 

