# Load Packages
library(tidyverse)
library(ggtreeExtra)
library(ggtree)
library(ggnewscale)
library(tidyverse)
library(ape)
library(treeio)
library(MetBrewer)
library(janitor)
library(ggstar)

# Read Tree
Elizabethkingia_Pop <- read.tree("ElizabethkingiaModel/100423_core_gene_alignment_snp_tree.treefile")
# Convert to Tibble
Elizabethkingia_Pop_tib <- as_tibble(Elizabethkingia_Pop)
# Rename 
Elizabethkingia_Pop_tib$label2 <- substring(Elizabethkingia_Pop_tib$label, 1, 15)
Elizabethkingia_Pop_tib$label2 <- gsub("_contig.*", "", Elizabethkingia_Pop_tib$label2)
Elizabethkingia_Pop_tib[["label2"]] <- gsub("\\.", "_", Elizabethkingia_Pop_tib[["label2"]])

Elizabethkingia_Pop_tib <- Elizabethkingia_Pop_tib %>%
  select(label, label2)

# Fix Tip Labels
Elizabethkingia_Pop_rename <- rename_taxa(Elizabethkingia_Pop, Elizabethkingia_Pop_tib, label, label2)

# Make rooted tree
Rooted_Elizabethkingia <- phytools::midpoint.root(Elizabethkingia_Pop_rename) 
Rooted_Elizabethkingia$tip.label <- gsub("#", "_", Rooted_Elizabethkingia$tip.label)

# Define Colours
Eliz_cols <- c("Elizabethkingia anophelis" = "#4e79a7",
          "Animal"="#FFEB61",
          "Clinical"="#C22F06",
          "Environmental"="#67993E",
          "Homo sapiens"="#DE9875",
          "Insect"="#67717F",
          "Missing"="white",
          "Mosquito"="black",
          "Plant"="#BAD24E",
          "Hu"="#EAD4FB",
          "Curated"="#667B8C",
          "MosAIC"="#00003D",
          "Food" = "#2A336C", 
          "Toxorhynchites amboinensis" = "#81B393", 
          "Aedes aegypti" = "#0E4D64", 
          "Aedes albopictus" = "#188977")

# Make Initial Tree
Eliz_Tree <- ggtree(Rooted_Elizabethkingia, layout="circular", size=0.2, right = T, ladderize = T, open.angle = 90) + 
  layout_fan(angle = 90) #+ geom_tiplab(size = 0.5)

# Read Metadata
Elizabethkingia_metadata <- read_tsv("ElizabethkingiaModel/Elizabethkingiaanophelis_metadata.tsv")

# Fix File Extensions
Elizabethkingia_metadata$File_prefix <- gsub("\\.", "_", Elizabethkingia_metadata$File_prefix)

# Edit Metadata Table
Elizabethkingia_metadata_edit <- Elizabethkingia_metadata %>%
  rename(Classification = "Species", 
         Host = "host_simple") %>%
  select(File_prefix, Host, Classification, Collection) %>%
  mutate(GTDB_Classification = NA, 
         NCBI_Classification = NA) %>%
  mutate(GTDB_Classification = if_else(Collection == "MosAIC", Classification, GTDB_Classification)) %>%
  mutate(NCBI_Classification = if_else(Collection == "Hu", Classification, NCBI_Classification))
  
# Function to Extract Specific Columns
ExtractMetdata <- function(data, col){
  table2 <- data %>%
    dplyr::select(File_prefix, !!col) %>%
    tibble::column_to_rownames(var = "File_prefix")
} 

E_Classification_Hu <- ExtractMetdata(Elizabethkingia_metadata_edit, "NCBI_Classification")
E_Classification_MosAIC <- ExtractMetdata(Elizabethkingia_metadata_edit, "GTDB_Classification")
E_Collection <- ExtractMetdata(Elizabethkingia_metadata_edit, "Collection")
E_Host <- ExtractMetdata(Elizabethkingia_metadata_edit, "Host")

# Highlighting Mosquito Associated Lineages
Eliz1 <- Eliz_Tree %<+% Elizabethkingia_metadata_edit + 
  # geom_tippoint(aes(col = Host), size = 0.5) + 
  scale_color_manual(values = c(Eliz_cols)) + 
  #geom_text(aes(label=node, size = 0.01), size = 0.6,nudge_x = 0.01) + 
  geom_highlight(node = c(131, 475, 132, 474, 133, 476, 134, 471, 129, 
                          473, 130, 472, 128, 294, 291, 293, 
                          292, 288, 563, 290, 289), fill = "red", alpha = 0.2, extendto = max(Eliz_Tree$data$x)) 

#ggsave("ElizabethkingiaModel/130423_ElizabethkingiaModelInitial.pdf")

# Tree with Classifications
Eliz2 <- Eliz_Tree %>% ggtree::gheatmap(E_Classification_Hu, color = NULL,
                        colnames_offset_y = 3,
                        colnames_angle = 0,
                        hjust = 0,
                        colnames_position = "top",
                        colnames = T,
                        width = 0.1,
                        offset = 0, 
                        font.size = 4) + 
  scale_fill_manual(values = c(Eliz_cols), na.value = "white") +
  theme(legend.position = "right", 
        legend.key.size = unit(0.9, "cm"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 15)) + 
  guides(fill = guide_legend(ncol = 2, title = "Classification")) + 
  geom_treescale(x = min(Eliz1$data$x) + 0.1, y = min(Eliz1$data$y) + 10,  offset = 10, fontsize = 3) #+ 
  #ggplot2::ylim(0, max(Eliz1$data$y) + 50) 

# Tree with GTDB Classifications
Eliz3 <- Eliz2 %>% ggtree::gheatmap(E_Classification_MosAIC, color = NULL,
                                        colnames_offset_y = 3,
                                        colnames_angle = 0,
                                        hjust = 0,
                                        colnames_position = "top",
                                        colnames = T,
                                        width = 0.1,
                                        offset = 0.03, 
                                        font.size = 4) + 
  scale_fill_manual(values = c(Eliz_cols), na.value = "white") +
  theme(legend.position = "right", 
        legend.key.size = unit(0.9, "cm"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 15)) + 
  guides(fill = guide_legend(ncol = 2, title = "Classification")) 
  #geom_treescale(x = min(Eliz1$data$x) + 0.1, y = min(Eliz1$data$y) + 10) #+ 
  #ggplot2::ylim(0, max(Eliz1$data$y) + 50) 

Eliz4 <- Eliz3 + new_scale_fill()

Eliz5 <- gheatmap(Eliz4, E_Collection, color = NULL,
               colnames_offset_y = 3,
               colnames_angle = 0,
               hjust = 0,
               colnames_position = "top",
               colnames = T,
               width = 0.1,
               offset = 0.06, 
               font.size = 4) + 
  scale_fill_manual(values = c(Eliz_cols), na.value = "white") + 
  guides(fill = guide_legend(ncol = 2, title = "Hu Collection")) + 
  new_scale_fill()

# Full tree with Classifications, collection source and host
Eliz6 <- gheatmap(Eliz5, E_Host, color = NULL,
                  colnames_offset_y = 3,
                  colnames_angle = 0,
                  hjust = 0,
                  colnames_position = "top",
                  colnames = T,
                  width = 0.1,
                  offset = 0.09, 
                  font.size = 4) + 
  scale_fill_manual(values = c(Eliz_cols)) + 
  guides(fill = guide_legend(ncol = 2, title = "Host")) 

#ggsave(filename = "ElizabethkingiaModel/030523_ElizabethkingiaPopStructure.pdf", height = 10, width = 10)
#ggsave(Eliz6, filename = "ElizabethkingiaModel/030523_ElizabethkingiaPopStructure_Collection_Legend.pdf", height = 10, width = 10)

# Read in Elizabethkingia anophelis metadata - clean + filter for isolation source 
Eliz_metadata <- read_tsv("ElizabethkingiaModel/Elizabethkingiaanophelis_metadata.tsv") %>%
  clean_names() 

Eliz_metadata$file_prefix <- gsub("\\.", "_", Eliz_metadata$file_prefix)

# Substitute all "." with "_"  
#E_metadata_final[["accession"]] <- gsub("\\.", "_", E_metadata_final[["accession"]])

Eliz_metadata_edit2 <- as.data.frame(Eliz_metadata %>% 
                                    select(file_prefix, host_simple, mos_aic_source_specific, 
                                           mos_aic_source_lab, location, mos_aic_lab_field, 
                                           mosquito_spp) %>%
                                    rename(ID = "file_prefix") %>%
                                    mutate(host = 1))

Eliz_subset <- tree_subset(Rooted_Elizabethkingia, "MMO−105", levels_back = 6)
Eliz_subset_1 <- Eliz_subset %>%
  ggtree(size = 0.2, layout = "fan", open.angle = 185)

Eliz_subset_2 <- Eliz_subset_1 %>% gheatmap(E_Host, colnames_offset_y = 20,
                                        colnames_angle = 60,
                                        hjust = 0.15,
                                        colnames_position = "top",
                                        colnames = T,
                                        width = 0.05,
                                        offset = 0.0, 
                                        font.size = 5) + 
  scale_fill_manual(values = c(Eliz_cols)) + 
  guides(fill = guide_legend(ncol = 2, title = "Source")) 


Eliz_subset_3 <- Eliz_subset_2 + geom_fruit(
  data = subset(Eliz_metadata_edit2, !is.na(mosquito_spp)),
  geom = geom_point, 
  mapping = aes(x = host, y = ID, col=mosquito_spp),
  offset = 0,
  grid.params=list(
    linetype=3,
    size=0.3
  )) + 
  scale_color_manual(values = Eliz_cols, 
                     name = "Mosquito Isolation Specific",
                     guide=guide_legend(
                       keywidth=0.3, 
                       keyheight=0.3, 
                       order=1), 
                     na.translate=FALSE) + 
  new_scale_fill() 


Eliz_subset_4 <- Eliz_subset_3 + geom_fruit(
  data = Eliz_metadata_edit2, 
  geom = geom_star, 
  mapping = aes(x = host, y = ID, starshape=mos_aic_source_lab),
  offset = -0.1, 
  grid.params=list(
    linetype=3,
    size=0.5
  )
) + 
  scale_fill_manual(values = Eliz_cols,
                    name = "Source Lab",
                    guide=guide_legend(
                      keywidth=0.3, 
                      keyheight=0.3, 
                      order=3
                    ),
                    na.translate=FALSE) + 
  new_scale_fill()

# Subset with host, lab and environment
Eliz_subset_5 <- Eliz_subset_4 + geom_fruit(
  data = Eliz_metadata_edit2, 
  geom = geom_star, 
  mapping = aes(x = host, y = ID, starshape=mos_aic_lab_field),
  offset = -0.1, 
  grid.params=list(
    linetype=3,
    size=0.5
  )
) +
  scale_starshape_manual(values = c("Coon" = "square diamond", "lab" = "circle", "Caragata, Eric" = "triangle")) + 
  new_scale_fill()

ggsave(Eliz_subset_5, filename = "ElizabethkingiaModel/030523_ElizAnophelis_Subset_SourceMetadata.pdf", height = 5, width = 10)