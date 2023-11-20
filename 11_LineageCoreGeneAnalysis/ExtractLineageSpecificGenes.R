###
# Script to analyse the lineage-specific mosquito-associated gene annotations from Bakta
###

library(tidyverse)
library(packcircles)

### Funtion to join and wrangle all relevent tables
JoinAndWrangleData <- function(panaroo_data, panaroo_annotations, lineage_data, 
                               population_data, species){
  
  ### Read Panaroo Gene Presence Absence File
  panaroo_data_tib <- as_tibble(read.table(panaroo_data, 
                                             header = T,
                                             sep="\t",
                                             quote = "",
                                             comment.char = "",
                                             stringsAsFactors = F,
                                             check.names = F))
  
  ### Read Annotations
  panaroo_annotation <- read_csv(panaroo_annotations) %>%
    select(Gene, Annotation)
  
  ### Read PopPUNK Lineage Clusters 
  lineage_data <- read_tsv(lineage_data, col_names = c("Sample", "Cluster")) #col_names = c("Sample", "Cluster"))  
  
  ### Read Twilight Classifications + wrangle
  pop_class_full <- read_tsv(population_data) %>%
    rename("Gene" = gene_name) %>%
    mutate(Core = ifelse(str_detect(details, "Core:"), str_extract(details, "(?<=Core: )[\\d+\\-]+"), NA_character_),
         Intermediate = ifelse(str_detect(details, "Inter:"), str_extract(details, "(?<=Inter: )[\\d+\\-]+"), NA_character_),
         Rare = ifelse(str_detect(details, "Rare:"), str_extract(details, "(?<=Rare: )[\\w+\\-]+"), NA_character_)) %>%
    select(Gene, core, inter, rare, total, general_class, specific_class, details, Core, Intermediate, Rare)
  
  ### Pivot the Panaroo Data Long
  panaroo_data_tib_lon <- panaroo_data_tib %>%
    pivot_longer(!Gene, values_to = "Pres_Abs", names_to = "sample")
  
  ### Join Panaroo Data and Gene Annotations 
  panaroo_data_tib_lon_anot <- panaroo_data_tib_lon %>%
    left_join(panaroo_annotation, by = c("Gene" = "Gene"))
  
  ### Join Sample and PopPUNK Defined Cluster 
  panaroo_data_tib_lon_anot_clust <- panaroo_data_tib_lon_anot %>%
    left_join(lineage_data, by = c("sample" = "Sample"))
  
  ### join with Twilgith Classifications 
  panaroo_data_tib_lon_anot_clust_twilight <- panaroo_data_tib_lon_anot_clust %>%
    left_join(pop_class_full, by = c("Gene" = "Gene")) %>%
    mutate(dataset = species)
  
  ### Return complete dataframe 
  return(panaroo_data_tib_lon_anot_clust_twilight)
  
}

# Enterobacter asburiae
Enterobacter_asburiae_gene_table <- JoinAndWrangleData("EnterobacterModel/EnterobacterAsburiae/090423_gene_presence_absence_clean.Rtab", 
                   "EnterobacterModel/EnterobacterAsburiae/gene_presence_absence.csv", 
                   "EnterobacterModel/EnterobacterAsburiae/090423_PopPUNK_clusters_refine_clean.tsv", 
                   "EnterobacterModel/EnterobacterAsburiae/output_classification_table.tab", 
                   "Enterobacter_asburiae")

# Serratia marscecens
Serratia_marscecens_gene_table <- JoinAndWrangleData("SerratiaModel/Serratia_marscecens/090423_gene_presence_absence_clean.Rtab", 
                                                     "SerratiaModel/Serratia_marscecens/gene_presence_absence.csv", 
                                                     "SerratiaModel/Serratia_marscecens/090423_PopPUNK_clusters_refine_clean.tsv", 
                                                     "SerratiaModel/Serratia_marscecens/104023_output_classification_table.tsv", 
                                                     "Serratia_marscecens")

# Elizabethkingia anophelis
Elizabethkingia_anphelis_gene_table <- JoinAndWrangleData("ElizabethkingiaModel/260423_gene_presence_absence_clean.Rtab", 
                                                          "ElizabethkingiaModel/260423_gene_presence_absence.csv", 
                                                          "ElizabethkingiaModel/260423_PopPUNK_clusters_refine_clean.tsv", 
                                                          "ElizabethkingiaModel/260423_Eliz_output_classification_table.tsv", 
                                                          "Elizabethkingia_anophelis")


### Which Lineages are Mosquito Associated? 
SummariseLineage <- function(data, lineage, title, subtitle){
  Cluster <- data %>%
    filter(details == lineage) %>%
    filter(Pres_Abs != 0)
  
  Cluster %>%
    group_by(sample, Annotation) %>%
    summarise(number_annotations = n()) %>%
    ggplot() +
    aes(x = sample, y = Annotation, fill = number_annotations) +
    geom_tile() + 
    theme_bw(base_size = 15) + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
    labs(title = title, subtitle = subtitle) + 
    xlab("Sample")
  
}

## First variable of function is the DF, second is the Twilight classification, third is the Classification class, fourth gives title to the plot
SummariseLineage(Enterobacter_asburiae_gene_table, "Core: 1 Inter:  Rare:", "Lineage specific core Genes", "Mosquito-Associated Lineage 3")
ggsave(filename = "EnterobacterModel/EnterobacterAsburiae/190423_MosquitoCluster3.pdf")

SummariseLineage(Enterobacter_asburiae_gene_table, "Core: 7 Inter:  Rare:", "Lineage specific core Genes", "Mosquito-Associated Lineage 1")
ggsave(filename = "EnterobacterModel/EnterobacterAsburiae/190423_MosquitoCluster1.pdf")

SummariseLineage(Enterobacter_asburiae_gene_table, "Core: 8 Inter:  Rare:", "Lineage specific core Genes", "Mosquito-Associated Lineage 2")
ggsave(filename = "EnterobacterModel/EnterobacterAsburiae/190423_MosquitoCluster2.pdf")

SummariseLineage(Enterobacter_asburiae_gene_table, "Core: 1+8 Inter:  Rare:", "Multi lineage core Genes", "Mosquito-Associated Lineages 1 + 2")
ggsave(filename = "EnterobacterModel/EnterobacterAsburiae/190423_MosquitoCluster1_2.pdf")

SummariseLineage(Enterobacter_asburiae_gene_table, "Core: 7+8 Inter:  Rare:", "Multi lineage core Genes", "Mosquito-Associated Lineages 1 + 2")

SummariseLineage(Enterobacter_asburiae_gene_table, "Core:  Inter: 8 Rare:", "Lineage specific core Genes", "Mosquito-Associated Lineage 4")
ggsave(filename = "EnterobacterModel/EnterobacterAsburiae/190423_MosquitoCluster4.pdf")

SummariseLineage(Serratia_marscecens_gene_table, "Core: 12 Inter:  Rare:", "Lineage specific core Genes", "Mosquito-Associated Lineage 12")
ggsave(filename = "SerratiaModel/Serratia_marscecens/190423_MosquitoCluster12.pdf")

SummariseLineage(Serratia_marscecens_gene_table, "Core: 24 Inter:  Rare:", "Lineage specific core Genes", "Mosquito-Associated Lineage 12")
ggsave(filename = "SerratiaModel/Serratia_marscecens/190423_MosquitoCluster24.pdf")

SummariseLineage(Serratia_marscecens_gene_table, "Core: 24 Inter:  Rare:", "Lineage specific core Genes", "Mosquito-Associated Lineage 12")

SummariseLineage(Serratia_marscecens_gene_table, "Core: 49 Inter:  Rare:", "Lineage specific core Genes", "Mosquito-Associated Lineage 12")
ggsave(filename = "SerratiaModel/Serratia_marscecens/190423_MosquitoCluster49.pdf")

SummariseLineage(Serratia_marscecens_gene_table, "Core: 68 Inter:  Rare:", "Lineage specific core Genes", "Mosquito-Associated Lineage 12")
ggsave(filename = "SerratiaModel/Serratia_marscecens/190423_MosquitoCluster68.pdf")

SummariseLineage(Serratia_marscecens_gene_table, "Core: 69 Inter:  Rare:", "Lineage specific core Genes", "Mosquito-Associated Lineage 12")
ggsave(filename = "SerratiaModel/Serratia_marscecens/190423_MosquitoCluster69.pdf")

SummariseLineage(Serratia_marscecens_gene_table, "Core: 12+24 Inter:  Rare:", "Multi lineage core Genes", "Mosquito-Associated Lineage 12")
ggsave(filename = "SerratiaModel/Serratia_marscecens/190423_MosquitoCluster12_68.pdf")

SummariseLineage(Serratia_marscecens_gene_table, "Core: 148 Inter:  Rare:", "Lineage specific core Genes", "Mosquito-Associated Lineage 12")
ggsave(filename = "SerratiaModel/Serratia_marscecens/190423_MosquitoCluster148.pdf")

SummariseLineage(Elizabethkingia_anphelis_gene_table, "Core: 10+7 Inter:  Rare:", "Multi lineage core Genes", "Mosquito-Associated Lineage 7+10")

### Join all the tables together
full_table <- rbind(Enterobacter_asburiae_gene_table, Serratia_marscecens_gene_table, Elizabethkingia_anphelis_gene_table)

### Plot genes specific to all mosquito associated lineages 
full_table %>%
  filter(
    details == "Core: 1 Inter:  Rare:" & dataset == "Enterobacter_asburiae" | 
      details == "Core: 7 Inter:  Rare:" & dataset == "Enterobacter_asburiae" | 
      details == "Core: 8 Inter:  Rare:" & dataset == "Enterobacter_asburiae" | 
      details == "Core: 12 Inter:  Rare:" & dataset == "Serratia_marscecens" | 
      details == "Core: 24 Inter:  Rare:" & dataset == "Serratia_marscecens" |
      details == "Core: 7 Inter:  Rare:" & dataset == "Elizabethkingia_model" |
      details == "Core: 10 Inter:  Rare:" & dataset == "Elizabethkingia_model") %>%
  filter(Pres_Abs != 0) %>%
  mutate(Annotation = str_remove(Annotation, "domain-containing protein"), 
         Annotation = str_remove(Annotation, "protein")) %>%
  group_by(sample, Annotation, dataset) %>%
  summarise(number_annotations = n()) %>%
  ggplot() + 
  aes(x = sample, y = Annotation, fill = number_annotations) + 
  geom_tile() + 
  facet_wrap(~dataset, scales = "free_x") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
  
ggsave(filename = "plots/Supplementary_Mosquito_Associated_Core_genes.pdf", width = 15, height = 20)
ggsave(filename = "plots/Supplementary_Mosquito_Associated_Core_genes.png", width = 15, height = 20)


### Convert table into a summarised dataframe with annotation and # of annotations to summarise
GetDataframe <- function(table, taget){
  full_table_simp <- table %>%
    filter(
      details == "Core: 1 Inter:  Rare:" & dataset == "Enterobacter_asburiae" | 
        details == "Core: 7 Inter:  Rare:" & dataset == "Enterobacter_asburiae" | 
        details == "Core: 8 Inter:  Rare:" & dataset == "Enterobacter_asburiae" | 
        details == "Core: 12 Inter:  Rare:" & dataset == "Serratia_marscecens" | 
        details == "Core: 24 Inter:  Rare:" & dataset == "Serratia_marscecens" |
        details == "Core: 7 Inter:  Rare:" & dataset == "Elizabethkingia_model") %>%
    filter(dataset == taget) %>%
    filter(Pres_Abs != 0) %>%
    filter(Annotation != "hypothetical protein") %>%
    mutate(Annotation = str_remove(Annotation, "domain-containing protein"), 
           Annotation = str_remove(Annotation, "Outer membrane protein"),
           Annotation = str_remove(Annotation, "domain involved in"), 
           Annotation = str_remove(Annotation, "site-specific")) %>%
    #Annotation = str_remove(Annotation, "protein")) %>%
    group_by(Annotation) %>%
    summarise(number_annotations = n()) %>%
    select(Annotation, number_annotations)
  #mutate(new_column = ifelse(n_distinct(dataset) == 1, as.character(unique(dataset)), "Overlapping_twice")) %>%
  #group_by(Annotation, new_column) %>%
  #summarise(total_annotations = sum(number_annotations))
  
  return(full_table_simp)
}

Data_EA <- GetDataframe(full_table, "Enterobacter_asburiae")
Data_SM <- GetDataframe(full_table, "Serratia_marscecens")
Data_Eliz <- GetDataframe(full_table, "Elizabethkingia_model")

# Most common annotation? 
Data_Eliz %>%
  arrange(desc(number_annotations))

# how many unqiue phage annotations? 
patterns <- c("phage", "Phage")

Data_Eliz %>% 
  filter(str_detect(Annotation, str_c(patterns, collapse = "|")))

Data_SM %>%
  filter(str_detect(Annotation, str_c(patterns, collapse = "|")))

Data_EA %>%
  filter(str_detect(Annotation, str_c(patterns, collapse = "|")))

### Make bubble plots
MakeCirclePlot <- function(data){
  
  packing <- circleProgressiveLayout(data$number_annotations, sizetype='area')
  
  packing$radius <- 0.9*packing$radius
  
  data <- cbind(data, packing)
  
  dat.gg <- circleLayoutVertices(packing, npoints=50)
  dat.gg$number_annotations <- rep(data$number_annotations, each=51)
  
  ggplot() + 
    
    # Make the bubbles
    geom_polygon(data = dat.gg, aes(x, y, group = id, fill=number_annotations), colour = "black", alpha = 0.6) +
    scale_fill_distiller(palette = "BuPu", direction = 1 ) +
    
    # Add text in the center of each bubble + control its size
    geom_text(data = data, aes(x, y, size=2, label = Annotation)) +
    scale_size_continuous(range = c(1,4)) +
    
    # General theme:
    theme_void() + 
    theme(legend.position="right") +
    coord_equal() 
  #facet_wrap(~dataset, scales = "free_y") 
}

MakeCirclePlot(Data_EA)
ggsave(filename = "plots/Enterobacter_asburiae_bubble_plot.pdf")

MakeCirclePlot(Data_SM)
ggsave(filename = "plots/Serratia_marscecens_bubble_plot.pdf")

MakeCirclePlot(Data_Eliz)
ggsave(filename = "plots/Elizabethkingia_anophelis_bubble_plot.pdf")

### Make a table with overlapping genes across mosquito associated lineages
table_with_overlap <- full_table %>%
  filter(
    details == "Core: 1 Inter:  Rare:" & dataset == "Enterobacter_asburiae" | 
      details == "Core: 7 Inter:  Rare:" & dataset == "Enterobacter_asburiae" | 
      details == "Core: 8 Inter:  Rare:" & dataset == "Enterobacter_asburiae" | 
      details == "Core: 12 Inter:  Rare:" & dataset == "Serratia_marscecens" | 
      details == "Core: 24 Inter:  Rare:" & dataset == "Serratia_marscecens" |
      details == "Core: 7 Inter:  Rare:" & dataset == "Elizabethkingia_model") %>%
  filter(Pres_Abs != 0) %>%
  mutate(Annotation = str_remove(Annotation, "-containing protein")) %>%
  #Annotation = str_remove(Annotation, "protein")) %>%
  group_by(Annotation, dataset) %>%
  summarise(number_annotations = n()) %>%
  group_by(Annotation) %>%
  mutate(new_column = case_when(
    n_distinct(dataset) == 1 ~ paste0("unqique_", unique(dataset)),
    n_distinct(dataset) == 2 ~ paste0("overlap_", paste(unique(dataset), collapse = "_")),
    n_distinct(dataset) >= 3 ~ paste0("overlap_", paste(unique(dataset), collapse = "_")),
    TRUE ~ ""
  )) %>%
  ungroup() 

### Graph with overlapping genes across the bacteria
table_with_overlap %>%
  filter(new_column != "unqique_Serratia_marscecens" &
           new_column != "unqique_Enterobacter_asburiae" &
           new_column != "unqique_Elizabethkingia_model" &
           new_column != "overlap_Elizabethkingia_model_Enterobacter_asburiae_Serratia_marscecens") %>%
  ggplot() + 
  aes(x = Annotation, y = new_column, fill = number_annotations) + 
  geom_tile() +
  coord_flip() + 
  theme_bw(base_size = 20) + 
  theme(axis.text.x = element_text(angle = 10, hjust = 1, size = 10)) + 
  xlab("Overlap") + 
  ylab("Annotation")

ggsave(filename = "plots/Overlapping_Mosquito_Associated_Genes.pdf")
ggsave(filename = "plots/Overlapping_Mosquito_Associated_Genes.png")
