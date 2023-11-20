# Load Packages
library(ape)
library(tidyverse)
library(ggtree)
library(phytools)
library(janitor)
library(treeio)
library(ggnewscale)

# Define column Colours
Chavda_cols <- c("Citrobacter braakii" = "#4e79a7",
          "E. aerogenes" = "#f28e2b",
          "E. asburiae" = "#e15759",
          "E. bugandensis" = "#C3D673",
          "E. chengduensis" = "#76b7b2",
          "E. cloacae" = "#394B18",
          "E. cloacae complex" = "#336220",
          "Enterobacter hormaechei_A" = "#edc948", 
          "E. cloacae complex 'Hoffmann cluster III'" = "#297A31",
          "E. cloacae complex 'Hoffmann cluster IV'" = "#32915E",
          "E. cloacae subsp. cloacae" = "#1BC200",
          "E. cloacae subsp. dissolvens" = "#1FA918",
          "E. hormaechei" = "#344FAA",
          "E. hormaechei subsp. hormaechei" = "#2760AA",
          "E. hormaechei subsp. hoffmannii" = "#1B77A9",
          "E. hormaechei subsp. oharae" =  "#1191A7",
          "E. hormaechei subsp. steigerwaltii" = "#07A49A",
          "E. kobei" = "#913431",
          "Mosquito" = "Black", 
          "NonMosquito" = "grey",
          "E. ludwigii" = "#90318C",
          "E. mori" = "#827A87",
          "E. roggenkampii" = "#76b7b2",
          "E. xiangfangensis" = "#B88694",
          "Enterobacter sp" = "black",
          "NonMosquito" = "white", 
          "Mosquito" = "#F781BF",
          "Animal"="#FFEB61",
          "Clinical"="#C22F06",
          "Environmental"="#67993E",
          "Homo sapiens"="#DE9875",
          "Insect"="#67717F",
          "Missing"="white",
          "Mosquito"="#F781BF",
          "Plant"="#BAD24E",
          "Chavda"="#91B98B",
          "Curated"="#667B8C",
          "MosAIC"="#00003D", 
          "Aedes aegypti" = "#0E4D64", 
          "Aedes albopictus" = "#188977", 
          "Aedes atropalpus" = "#39A96B", 
          "Anopheles gambiae" = "#6B6850", 
          "Culex pipiens" = "#815A93")

# Chavda's Metadata
chavda_metadata <- read_tsv("EnterobacterModel/110423_Chavda_Mos_Curated_Metadata.tsv") %>%
  clean_names() %>%
  select(gan_bank_accession_s, proposed_name, host, assembly, host_alt, paper_name, host_simple, collection)
  
# Read Tree
Enterobacter_tree <- read.tree("EnterobacterModel/Genus_Pop/060423_core_gene_alignment_snps_tree.treefile")

# Sort out the Tip Labels
Enterobacter_tree_tib <- as_tibble(Enterobacter_tree) %>%
  separate(col = label, into = c("label2", "junk"), sep = "_ASM", remove = F) %>%
  separate(col = label2, into = c("label2", "junk"), sep = "_contigs", remove = F) %>%
  separate(col = label2, into = c("label2", "junk"), sep = "_SPADES", remove = F) %>%
  separate(col = label2, into = c("label2", "junk"), sep = "_SMART", remove = F) %>%
  separate(col = label2, into = c("label2", "junk"), sep = "_SAMN", remove = F) %>%
  separate(col = label2, into = c("label2", "junk"), sep = "_Assembly", remove = F) %>%
  separate(col = label2, into = c("label2", "junk"), sep = "_DC1", remove = F) %>%
  separate(col = label2, into = c("label2", "junk"), sep = "_L1", remove = F) %>%
  separate(col = label2, into = c("label2", "junk"), sep = "_UHGG", remove = F) %>%
  separate(col = label2, into = c("label2", "junk"), sep = "_Ente", remove = F) %>%
  separate(col = label2, into = c("label2", "junk"), sep = "11983", remove = F) %>%
  select(label, label2) 

# Rename tree tips
EA_tree_rename <- rename_taxa(Enterobacter_tree, Enterobacter_tree_tib, label, label2)

# Rooting at Enterobacter aerogenes (Now Klebsiella aerogenes)
EA_tree_root <- phytools::midpoint.root(EA_tree_rename)

# Reformat the Metadata for ggTree
chavda_meta_paper <- chavda_metadata %>%
  select(assembly, paper_name, host_simple, collection) %>%
  rename(c(Classification = "paper_name", 
           Host = "host_simple", 
           Collection = "collection")) %>%
  mutate(GTDB_Classification = NA, 
         NCBI_Classification = NA) %>%
  mutate(GTDB_Classification = if_else(Collection == "MosAIC", Classification, GTDB_Classification)) %>%
  mutate(NCBI_Classification = if_else(Collection == "Chavda" | Collection == "Curated", Classification, NCBI_Classification))

# Function to extract columns of interest
ExtractMetdata <- function(data, col){
  table2 <- data %>%
    dplyr::select(assembly, !!col) %>%
    distinct(assembly, .keep_all = T) %>%
    filter(!is.na(assembly)) %>%
    tibble::column_to_rownames(var = "assembly")
} 

GTDB_Classification <- ExtractMetdata(chavda_meta_paper, "GTDB_Classification")
NCBI_Classification <- ExtractMetdata(chavda_meta_paper, "NCBI_Classification")
Host <- ExtractMetdata(chavda_meta_paper, "Host")
Collection <- ExtractMetdata(chavda_meta_paper, "Collection")

# Make the Tree 
Enterobacter_Tree <- ggtree(EA_tree_root, layout =  "circular", size = 0.2) + 
  layout_fan(angle = 90)

#ggsave(plot = Enterobacter_Tree, filename = "EnterobacterModel/PlotsFinal/240423_Chavda_Mos_PopulationStructure_Tip.pdf")

# Append metadata using gheatmap
EP <- Enterobacter_Tree %>% ggtree::gheatmap(NCBI_Classification, color = NULL,
                        colnames_offset_y = 3,
                        colnames_angle = 0,
                        hjust = 0,
                        colnames_position = "top",
                        colnames = T,
                        width = 0.1,
                        offset = 0, 
                        font.size = 5) + 
  scale_fill_manual(values = c(Chavda_cols), na.value = "white") +
  theme(legend.position = "right", 
        legend.key.size = unit(0.9, "cm"), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10)) + 
  guides(fill = guide_legend(ncol = 1, title = "Classification")) + 
  geom_treescale(x = min(Enterobacter_Tree$data$x) + 0.1, y = min(Enterobacter_Tree$data$y) + 10) #+ 
  #ggplot2::ylim(0, max(Enterobacter_Tree$data$y) + 50) 

EP2 <- EP + new_scale_fill()

EP3 <- EP2 %>% ggtree::gheatmap(GTDB_Classification, color = NULL,
                                             colnames_offset_y = 3,
                                             colnames_angle = 0,
                                             hjust = 0,
                                             colnames_position = "top",
                                             colnames = T,
                                             width = 0.1,
                                             offset = 0.07, 
                                             font.size = 5) + 
  scale_fill_manual(values = c(Chavda_cols), na.value = "white") +
  theme(legend.position = "right", 
        legend.key.size = unit(0.9, "cm"), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10)) + 
  guides(fill = guide_legend(ncol = 1, title = "Classification")) + 
  geom_treescale(x = min(Enterobacter_Tree$data$x) + 0.1, y = min(Enterobacter_Tree$data$y) + 10) #+ 
#ggplot2::ylim(0, max(Enterobacter_Tree$data$y) + 50) 

EP4 <- EP3 + new_scale_fill()


EP5 <- gheatmap(EP4, Collection, colnames_offset_y = 3,
         colnames_angle = 0,
         hjust = 0,
         colnames_position = "top",
         colnames = T,
         width = 0.1,
         offset = 0.14, 
         font.size = 5) + 
  scale_fill_manual(values = c(Chavda_cols)) + 
  guides(fill = guide_legend(ncol = 1, title = "Collection")) + 
  new_scale_fill()
  
# Final phylogeny with Classification, collection and host
EP6 <- gheatmap(EP5, Host, colnames_offset_y = 3,
                colnames_angle = 0,
                hjust = 0,
                colnames_position = "top",
                colnames = T,
                width = 0.1,
                offset = 0.21, 
                font.size = 5) + 
  scale_fill_manual(values = c(Chavda_cols), na.value = "white") + 
  guides(fill = guide_legend(ncol = 2, title = "Source")) 
  
# Write to File
ggsave(plot = EP6, filename = "EnterobacterModel/PlotsFinal/Figure_3_Chavda_Mos_Curated_PopulationStructure_Circular.pdf", height = 10, width = 10)
ggsave(plot = EP6, filename = "EnterobacterModel/PlotsFinal/Figure_3_Chavda_Mos_Curated_PopulationStructure_Circular_Fixed.pdf", height = 13, width = 10)

# Subset Asburiaes 

# Read in Enterobacter asrbuiae metadata - clean + filter for isolation source 
E_metadata <- read_tsv("EnterobacterModel/110423_Chavda_Mos_Curated_Metadata.tsv") %>%
  clean_names() 

# Substitute all "." with "_"  
#E_metadata_final[["accession"]] <- gsub("\\.", "_", E_metadata_final[["accession"]])
E_metadata[["assembly"]] <- gsub("\\.", "_", E_metadata[["assembly"]])

E_metadata_edit <- E_metadata %>%
  select(assembly, host_simple) %>%
  distinct(assembly, .keep_all = T) %>%
  filter(!is.na(assembly)) %>%
  column_to_rownames(var = "assembly")

E_metadata_edit2 <- as.data.frame(E_metadata %>% 
                                    select(!gan_bank_accession_s) %>% 
                                    select(assembly, host_simple, mos_aic_source_specific, 
                                           mos_aic_source_lab, location, mos_aic_lab_field, 
                                           mosquito_spp) %>%
                                    rename(ID = "assembly") %>%
                                    mutate(host = 1))
# Make Subsets of the Tree
EA_Subset <- tree_subset(EA_tree_root, "USMM051", levels_back = 19)

EA_Subset_Tree <- EA_Subset %>%
  ggtree(size = 0.2, layout = "fan", open.angle = 185) 

EA_subset_1 <- EA_Subset_Tree %>% gheatmap(Host, colnames_offset_y = 20,
                            colnames_angle = 60,
                            hjust = 0.15,
                            colnames_position = "top",
                            colnames = T,
                            width = 0.05,
                            offset = 0.0, 
                            font.size = 5) + 
  scale_fill_manual(values = c(Chavda_cols)) + 
  guides(fill = guide_legend(ncol = 2, title = "Source")) + 
  new_scale_fill()


EA_subset_2 <- EA_subset_1 + geom_fruit(
  data = subset(E_metadata_edit2, !is.na(mosquito_spp)),
  geom = geom_point, 
  mapping = aes(x = host, y = ID, col=mosquito_spp),
  offset = 0,
  grid.params=list(
    linetype=3,
    size=0.5
  )) + scale_color_manual(values = Chavda_cols,
                         name = "Mosquito Isolation Specific",
                         guide=guide_legend(
                           keywidth=0.3, 
                           keyheight=0.3, 
                           order=1
                         ),
                         na.translate=FALSE) + 
  new_scale_fill() 

EA_subset_3 <- EA_subset_2 + geom_fruit(
  data = E_metadata_edit2, 
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

# Final subset with host, species, lab, environment
EA_subset_4 <- EA_subset_3 + geom_fruit(
  data = E_metadata_edit2, 
  geom = geom_star, 
  mapping = aes(x = host, y = ID, starshape=mos_aic_lab_field),
  offset = -0.1, 
  grid.params=list(
    linetype=3,
    size=0.5
  )
) + 
  scale_fill_manual(values = Chavda_cols) + 
  scale_starshape_manual(values = c("Coon" = "square diamond", "lab" = "circle", "field" = "hexagonal star", 
                                    "Povelones, Michael" = "regular pentagon", "Valiente Moro, Claire" = "antiparallelogram")) + 
  new_scale_fill()

ggsave(EA_subset_4, filename = "EnterobacterModel/PlotsFinal/030523_EA_Subset_SourceMetadata.pdf", height = 5, width = 10)

# Same again but with portion focussing on Enterobacter roggenkampii
ER_Subset <- tree_subset(EA_tree_root, "GCF_000958015.1", levels_back = 10)
ER_Subset_Tree <- ER_Subset %>%
  ggtree(size = 0.2, layout = "fan", open.angle = 185) 

ER_subset_1 <- ER_Subset_Tree %>% gheatmap(Host, colnames_offset_y = 20,
                                           colnames_angle = 60,
                                           hjust = 0.15,
                                           colnames_position = "top",
                                           colnames = T,
                                           width = 0.05,
                                           offset = 0.0, 
                                           font.size = 5) + 
  scale_fill_manual(values = c(Chavda_cols)) + 
  guides(fill = guide_legend(ncol = 2, title = "Source")) + 
  new_scale_fill()


ER_subset_2 <- ER_subset_1 + geom_fruit(
  data = E_metadata_edit2, 
  geom = geom_point, 
  mapping = aes(x = host, y = ID, col=mosquito_spp),
  offset = 0,
  grid.params=list(
    linetype=3,
    size=0.5
  )) + scale_color_manual(values = Chavda_cols,
                        name = "Mosquito Isolation Specific",
                        guide=guide_legend(
                          keywidth=0.3, 
                          keyheight=0.3, 
                          order=1
                        ),
                        na.translate=FALSE) + 
  new_scale_fill() 

ER_subset_3 <- ER_subset_2 + geom_fruit(
  data = E_metadata_edit2, 
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

ER_subset_4 <- ER_subset_3 + geom_fruit(
  data = E_metadata_edit2, 
  geom = geom_star, 
  mapping = aes(x = host, y = ID, starshape=mos_aic_lab_field),
  offset = -0.1, 
  grid.params=list(
    linetype=3,
    size=0.5
  )
) + 
  scale_fill_manual(values = MosAIC_cols) + 
  scale_starshape_manual(values = c("Coon" = "square diamond", "lab" = "circle", "field" = "hexagonal star")) + 
  new_scale_fill()

ggsave(ER_subset_4, filename = "EnterobacterModel/PlotsFinal/00523_ER_Subset_SourceMetadata.pdf", height = 5, width = 10)


#### Subset and Write Accessions of Enterobacter asburiae to File 
# Which nodes are the first and last E. asburiae located? 
EP2 <- ggtree(EA_tree, size = 0.2) 

Enterobacter_asburiae_accessions <- as_tibble(get_taxa_name(EP2)) %>%
  rownames_to_column() %>%
  filter(between(rowname, 1, 106)) %>%
  select(value)  %>%
  print(n = 106)

#Enterobacter_asburiae_accessions %>%
  #filter(value == "USMM053" | value == "GCF_900075425.1_")
  
Enterobacter_asburiae_accessions$value <- paste0(Enterobacter_asburiae_accessions$value, ".gff3")
                           
write.table(x = Enterobacter_asburiae_accessions, file = "EnterobacterModel/050323_EnterobacterAAccessions.tsv", quote = F, sep = "\t", row.names = F, col.names = F)


#### Stats 
#### How many species in this population? 
#Enterobacter_tree_tib %>%
  #head(451) %>%
  #distinct(label2) %>%
  #print(n = 500) # 451

#### How many different species did we put in this population?
#E_metadata %>% 
  #filter(host_simple == "Mosquito") %>%
  #group_by(proposed_name) %>%
  #summarise(n())

