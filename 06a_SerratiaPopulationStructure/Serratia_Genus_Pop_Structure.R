##
# Visualise the Population Structure of Serratia
##

# Load Packages
library(ape)
library(ggtree)
library(ggtreeExtra)
library(treeio)
library(phytools)
library(tidyverse)
library(ggnewscale)
library(janitor)

# Read Tree
Serratia_pop <- read.tree("SerratiaModel/280323_core_gene_alignment_SNP_filt.aln.contree")

# Root Tree
rooted_serratia_tree <- phytools::midpoint.root(Serratia_pop)

# Rename Tip Labels
rooted_serratia_tree$tip.label <- gsub("#", "_", rooted_serratia_tree$tip.label)

## colours:
cols <- c("Serratia entomophila" = "#4e79a7",
          "Serratia ficaria" = "#f28e2b",
          "Serratia fonticola" = "#e15759",
          "Serratia liquefaciens" = "#76b7b2",
          "Serratia marcescens" = "#59a14f",
          "Serratia plymuthica" = "#edc948",
          "Serratia proteamaculans" = "#b07aa1",
          "Serratia rubidaea" = "#ff9da7",
          "Serratia grimesii" = "#9c755f",
          "Serratia proteamaculans subs. quinovorans" = "#bab0ac",
          "Serratia proteamaculans\nsubs. quinovorans" = "#bab0ac",
          "Serratia sp." = "#d3d3d3",
          "Unknown" = "white",
          "uncultured Serratia" = "#a0a0a0",
          "Rhodotorula mucilaginosa/Serratia"="#BF5B17",
          "Rhodotorula mucilaginosa/\nSerratia"="#BF5B17",
          "Serratia aquatilis"="#FF7F00",
          "Serratia myotis"="#E5C494",
          "Serratia nematodiphila"="#FFED6F",
          "Enterobacter liquefaciens"="#8D9130",
          "Serratia quinivorans"="#D95F02",
          "Serratia Serratia"="#CAB2D6",
          "Serratia sp.XT-30"="#E6F5C9",
          "Serratia symbiont"="#386CB0",
          "Serratia oryzae" = "#a0a0a0",
          "Serratia odorifera"="#F781BF",
          "Serratia symbiotica"="#B3E2CD",
          "Serratia ureilytica"="#A6CEE3",
          "Serratia vespertilionis"="#FFFFCC",
          "Serratia nevei" = "#9604B3",
          "Serratia bockelmannii" = "#C5CC00",
          "Animal"="#FFEB61",
          "Clinical"="#C22F06",
          "Environmental"="#67993E",
          "Homo sapiens"="#DE9875",
          "Insect"="#67717F",
          "Missing"="white",
          "Mosquito"="black",
          "Plant"="#BAD24E",
          "Williams"="#7D9CA9",
          "Curated"="#667B8C",
          "MosAIC"="#00003D",
          "Food" = "#2A336C", 
          "Aedes aegypti" = "#0E4D64", 
          "Aedes triseriatus" = "#188977", 
          "Aedes atropalpus" = "#39A96B", 
          "Anopheles stephensi" = "#6B4850", 
          "Toxorhynchites amboinensis" = "#81B393")

# Initial Tree Drawing
S1 <- ggtree(rooted_serratia_tree, layout="circular", size=0.2, right = T, ladderize = T) + 
  layout_fan(angle = 90)

ggsave(S1, filename = "SerratiaModel/GenusPop/130423_InitialTree.pdf", height = 40)

# Read Metadata
serratia_metadata <- read.csv("SerratiaModel/serratia_metadata.csv", header = T,
                              quote = "",
                              stringsAsFactors = F,
                              comment.char = "")

# Edit Metadata to fit with ggtree
serratia_metadata_edit <- serratia_metadata %>%
  select(File_prefix,Short.species, Species, Paper,Generic_host,Country, Mosquito, Collection, host_simple) %>%
  rename(`Labelled species` = Short.species, 
         Dataset = Paper, Source = Generic_host, 
         `MOSaic Collection` = Mosquito, 
         Classification = "Species", 
         Host = "host_simple") %>%
  mutate(GTDB_Classification = NA, 
         NCBI_Classification = NA) %>%
  mutate(GTDB_Classification = if_else(Collection == "MosAIC", Classification, GTDB_Classification)) %>%
  mutate(NCBI_Classification = if_else(Collection == "Williams", Classification, NCBI_Classification))


# Function to extract relevent metadata for each ring in the tree
ExtractMetdata <- function(data, col){
  table2 <- data %>%
    dplyr::select(File_prefix, !!col) %>%
    distinct(File_prefix, .keep_all = T) %>%
    filter(!is.na(File_prefix)) %>%
    tibble::column_to_rownames(var = "File_prefix")
} 


GTDB_Classification <- ExtractMetdata(serratia_metadata_edit, "GTDB_Classification")
NCBI_Classification <- ExtractMetdata(serratia_metadata_edit, "NCBI_Classification")
Collection <- ExtractMetdata(serratia_metadata_edit, "Collection")
Host <- ExtractMetdata(serratia_metadata_edit, "Host")

# Make Trees with Appended Metadata
S2 <- S1 %>% ggtree::gheatmap(NCBI_Classification, color = NULL,
                                             colnames_offset_y = 3,
                                             colnames_angle = 0,
                                             hjust = 0,
                                             colnames_position = "top",
                                             colnames = T,
                                             width = 0.1,
                                             offset = 0, 
                                             font.size = 3) + 
  scale_fill_manual(values = c(cols), na.value = "white") +
  theme(legend.position = "right", 
        legend.key.size = unit(0.9, "cm"), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10)) + 
  guides(fill = guide_legend(ncol = 1, title = "Classification")) + 
  geom_treescale(x = min(S1$data$x) + 0.1, y = min(S1$data$y) + 10, offset = 10, fontsize = 3) #+ 
  #ggplot2::ylim(0, max(S1$data$y) + 50) 

S3 <- S2 + new_scale_fill()

S4 <- S3 %>% ggtree::gheatmap(GTDB_Classification, color = NULL,
                              colnames_offset_y = 3,
                              colnames_angle = 0,
                              hjust = 0,
                              colnames_position = "top",
                              colnames = T,
                              width = 0.1,
                              offset = 0.04, 
                              font.size = 3) + 
  scale_fill_manual(values = c(cols), na.value = "white") +
  theme(legend.position = "right", 
        legend.key.size = unit(0.9, "cm"), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10)) + 
  guides(fill = guide_legend(ncol = 1, title = "Classification")) 

S5 <- S4 + new_scale_fill()

S6 <- gheatmap(S5, Collection, colnames_offset_y = 3,
                colnames_angle = 0,
                hjust = 0,
                colnames_position = "top",
                colnames = T,
                width = 0.1,
                offset = 0.08, 
                font.size = 3) + 
  scale_fill_manual(values = c(cols)) + 
  guides(fill = guide_legend(ncol = 2, title = "Collection")) 

S7 <- S6 + new_scale_fill()

# Final Tree 
S8 <- gheatmap(S7, Host, colnames_offset_y = 3,
                colnames_angle = 0,
                hjust = 0,
                colnames_position = "top",
                colnames = T,
                width = 0.1,
                offset = 0.12, 
                font.size = 3) + 
  scale_fill_manual(values = c(cols)) + 
  guides(fill = guide_legend(nrow = 2, title = "Host", position = "bottom")) 

# Save Files
ggsave(filename = "SerratiaModel/PlotsFinal/Figure_3_Williamns_MosAIC_PopulationStructure_Circular.pdf", height = 10, width = 10)
ggsave(filename = "SerratiaModel/PlotsFinal/Figure_3_Williamns_MosAIC_PopulationStructure_Circular_Fixed.pdf",S8, height = 13, width = 10)

# MAKE SUBSETS
# Read in Serratia marscecens metadata - clean + filter for isolation source 
S_metadata <- read_csv("SerratiaModel/serratia_metadata.csv") %>%
  clean_names() 

# Substitute all "." with "_"  
#E_metadata_final[["accession"]] <- gsub("\\.", "_", E_metadata_final[["accession"]])

S_metadata_edit2 <- as.data.frame(S_metadata %>% 
                                    select(file_prefix, host_simple, mos_aic_source_specific, 
                                           mos_aic_source_lab, location_77, mos_aic_lab_field, 
                                           mosquito_spp) %>%
                                    rename(ID = "file_prefix") %>%
                                    mutate(host = 1))


SM_subset <- tree_subset(rooted_serratia_tree, "HN118_contigs.fa", levels_back = 10)
SM_subset_1 <- SM_subset %>%
  ggtree(size = 0.2, layout = "fan", open.angle = 185) 

SM_subset_2 <- SM_subset_1 %>% gheatmap(Host, colnames_offset_y = 20,
                                           colnames_angle = 60,
                                           hjust = 0.15,
                                           colnames_position = "top",
                                           colnames = T,
                                           width = 0.05,
                                           offset = 0.0, 
                                           font.size = 5) + 
  scale_fill_manual(values = c(cols)) + 
  guides(fill = guide_legend(ncol = 2, title = "Source")) 


SM_subset_3 <- SM_subset_2 + geom_fruit(
  data = subset(S_metadata_edit2, !is.na(mosquito_spp)),
  geom = geom_point, 
  mapping = aes(x = host, y = ID, col=mosquito_spp),
  offset = 0,
  grid.params=list(
    linetype=3,
    size=0.3
  )) + 
  scale_color_manual(values = cols, 
                     name = "Mosquito Isolation Specific",
                     guide=guide_legend(
                       keywidth=0.3, 
                       keyheight=0.3, 
                       order=1), 
                     na.translate=FALSE) + 
  new_scale_fill() 
  

SM_subset_4 <- SM_subset_3 + geom_fruit(
  data = S_metadata_edit2, 
  geom = geom_star, 
  mapping = aes(x = host, y = ID, starshape=mos_aic_source_lab),
  offset = -0.1, 
  grid.params=list(
    linetype=3,
    size=0.5
  )
) + 
  scale_fill_manual(values = MosAIC_cols,
                    name = "Source Lab",
                    guide=guide_legend(
                      keywidth=0.3, 
                      keyheight=0.3, 
                      order=3
                    ),
                    na.translate=FALSE) + 
  new_scale_fill()

# Complete Subset of the Phylogeny with 1 Species, 2 Lab 3 Environment
SM_subset_5 <- SM_subset_4 + geom_fruit(
  data = S_metadata_edit2, 
  geom = geom_star, 
  mapping = aes(x = host, y = ID, starshape=mos_aic_lab_field),
  offset = -0.1, 
  grid.params=list(
    linetype=3,
    size=0.5
  )
) +
  scale_starshape_manual(values = c("Coon" = "square diamond", "lab" = "circle", "Caragata, Eric" = "triangle", 
                                    "Brackney,Doug" = "regular triangle down", "Chen;S. and Walker;E." = "ellipse", 
                                    "Jacobs-Lorena,Marcelo" = "rhombus", "Wang, S" = "left-triangle1")) + 
  new_scale_fill()

ggsave(SM_subset_5, filename = "SerratiaModel/PlotsFinal/030523_SM_Subset_SourceMetadata.pdf", height = 5, width = 10)

# Second subset but for Serratia fonticola clades
SF_Subset <- tree_subset(rooted_serratia_tree, "HN157_contigs.fa", levels_back = 5)
SF_Subset_tree <- SF_Subset %>%
  ggtree(size = 0.2, layout = "fan", open.angle = 185) 

SF_Subset_1 <- SF_Subset_tree %>% gheatmap(Host, colnames_offset_y = 100,
                                           colnames_angle = 0,
                                           hjust = 0,
                                           colnames_position = "bottom",
                                           colnames = T,
                                           width = 0.05,
                                           offset = 0.0, 
                                           font.size = 5) + 
  scale_fill_manual(values = c(cols)) + 
  guides(fill = guide_legend(ncol = 2, title = "Source")) + 
  new_scale_fill()


SF_Subset_2 <- SF_Subset_1 + geom_fruit(
  data = subset(S_metadata_edit2, !is.na(mosquito_spp)),
  geom = geom_point, 
  mapping = aes(x = host, y = ID, col=mosquito_spp),
  offset = 0,
  grid.params=list(
    linetype=3,
    size=0.3
  )) + 
  scale_color_manual(values = cols, 
                     name = "Mosquito Isolation Specific",
                     guide=guide_legend(
                       keywidth=0.3, 
                       keyheight=0.3, 
                       order=1), 
                     na.translate=FALSE) + 
  new_scale_fill() 

SF_Subset_3 <- SF_Subset_2 + geom_fruit(
  data = S_metadata_edit2, 
  geom = geom_star, 
  mapping = aes(x = host, y = ID, starshape=mos_aic_source_lab),
  offset = -0.1, 
  grid.params=list(
    linetype=3,
    size=0.5
  )
) + 
  scale_fill_manual(values = MosAIC_cols,
                    name = "Source Lab",
                    guide=guide_legend(
                      keywidth=0.3, 
                      keyheight=0.3, 
                      order=3
                    ),
                    na.translate=FALSE) + 
  new_scale_fill()

SF_Subset_4 <- SF_Subset_3 + geom_fruit(
  data = S_metadata_edit2, 
  geom = geom_star, 
  mapping = aes(x = host, y = ID, starshape=mos_aic_lab_field),
  offset = -0.1, 
  grid.params=list(
    linetype=3,
    size=0.5
  )
) + 
  scale_fill_manual(values = cols) + 
  scale_starshape_manual(values = c("Coon" = "square diamond", "lab" = "circle", "Pei D" = "plus filled")) + 
  new_scale_fill()

ggsave(SF_Subset_4, filename = "SerratiaModel/PlotsFinal/030523_SF_Subset_SourceMetadata.pdf", height = 5, width = 10)

# Clades corresponding to Serratia marscecens
#Serratia_marscecens_clade <- viewClade(p, MRCA(p, "GCA_001902635.1_ASM190263v1_genomic", "GCA_001539745.1_12082_2_93_genomic"))
#MRCA(rooted_serratia_tree, "GCA_001902635.1_ASM190263v1_genomic", "GCA_001539745.1_12082_2_93_genomic")

#### Get accession for Just S. marscecens 
#Serratia_marscecens_accessions <- as_tibble(get_taxa_name(p)) %>%
  #rownames_to_column() %>%
  #filter(between(rowname, 265, 683)) %>%
  #select(value) 

#Serratia_marscecens_accessions %>%
  #filter(value == "GCA_001902635.1_ASM190263v1_genomic" | value == "GCA_001539745.1_12082_2_93_genomic")
  
#Serratia_marscecens_accessions$value <- ifelse(grepl("contigs", Serratia_marscecens_accessions$value), 
                                    #paste0(Serratia_marscecens_accessions$value, ".gff3"), 
                                    #paste0(Serratia_marscecens_accessions$value, ".gff"))

# Use these accessions to construct the species pangenome 
#write.table(x = Serratia_marscecens_accessions, file = "SerratiaModel/050323_SerratiaMAccesions.tsv", quote = F, sep = "\t", row.names = F, col.names = F)




