# Load Packages 
library(tidyverse)
library(ggtreeExtra)
library(ggtree)
library(ggnewscale)
library(ape)
library(treeio)
library(MetBrewer)

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


# Read in the Enterobacter Tree 
serratia_tree <- ape::read.tree("MosAIC_V1/07a_SerratiaPangenome/060323_core_gene_alignment_snp_tree.treefile")
# Root - rooting at the split defined by the genera tree
serratia_tree_root <- phytools::midpoint.root(serratia_tree)

# Read Panaroo Gene Presence Absence File
S_panaroo_data <- read.table("MosAIC_V1/07a_SerratiaPangenome/090423_gene_presence_absence_clean.Rtab", 
                             header = T,
                             sep="\t",
                             quote = "",
                             comment.char = "",
                             stringsAsFactors = F,
                             check.names = F)

# Read the Twilight Gene Classification 
S_pop_class <- read_tsv("MosAIC_V1/07a_SerratiaPangenome/104023_output_classification_table.tsv") %>%
  rename("Gene" = gene_name) %>%
  select(Gene, specific_class)

# Join Panaroo Gene Presence Absence + Twilight Classifications 
S_pangenome_w_classification <- S_panaroo_data %>%
  pivot_longer(!Gene, names_to =  "strain", values_to = "presence") %>%
  left_join(S_pop_class, by =  "Gene") %>%
  print(n = 50000)

#pangenome_w_classification[["strain"]] <- gsub("\\.", "_", pangenome_w_classification[["strain"]])

# Read in PopPUNK Lineage Clusters 
S_lineage_data <- read_tsv("MosAIC_V1/07a_SerratiaPangenome/090423_PopPUNK_clusters_refine_clean.tsv") #col_names = c("Sample", "Cluster"))  
# Convert the Cluster column from dbl to fct
S_lineage_data$Cluster <- as.factor(S_lineage_data$Cluster)

# Join with Panaroo Gene Presence + Absence + Twilight Classifications
S_pangenome_classfication_lineage <- S_pangenome_w_classification %>%
  left_join(S_lineage_data, by = c("strain" =  "Sample"))

S_pangenome_classfication_lineage %>% 
  group_by(Cluster) %>% summarise(n()) %>% 
  print(n = 400)

# Read in Enterobacter asrbuiae metadata - clean + filter for isolation source 
S_metadata <- read_csv("MosAIC_V1/07a_SerratiaPangenome/serratia_metadata.csv") %>%
  select(File_prefix, Collection, host_simple, Species)

S_metadata$File_prefix <- substring(S_metadata$File_prefix, 1, 15)
S_metadata$File_prefix <- gsub("_contig.*", "", S_metadata$File_prefix)
S_metadata[["File_prefix"]] <- gsub("\\.", "_", S_metadata[["File_prefix"]])

colnames(S_metadata)
S_metadata$File_prefix
S_metadata %>% group_by(host_simple) %>%
  summarise(n())

# Establish Gals Twilight Collection Categories
twilight_order <- c("Collection core", 
                    "Multi-lineage core",
                    "Lineage specific core", 
                    #"Collection intermediate",
                    "Multi-lineage intermediate",
                    "Lineage specific intermediate",
                    #"Collection rare",
                    "Multi-lineage rare",
                    "Lineage specific rare",
                    "Core and intermediate",
                    "Core and rare",
                    "Core, intermediate and rare",
                    "Intermediate and rare")


### Change the Tip Labels of the Tree to Match PopPUNK Custers and Panaroo Gene Presence / Absence
#### Bit messy - is there a cleaner way to do this? 
serratia_tree_tib <- as_tibble(serratia_tree)
serratia_tree_tib$label2 <- substring(serratia_tree_tib$label, 1, 15)
serratia_tree_tib$label2 <- gsub("_contig.*", "", serratia_tree_tib$label2)

### Substitute all "." for "_"
serratia_tree_tib[["label2"]] <- gsub("\\.", "_", serratia_tree_tib[["label2"]])

serratia_tree_tib <- serratia_tree_tib %>%
  select(label, label2)
### Rename the tip labels in the tree file 
serratia_tree_rename <- rename_taxa(serratia_tree_root, serratia_tree_tib, label, label2)

serratia_tree_rename$tip.label

### Draw the renamed, rooted tree
S <- ggtree(serratia_tree_rename, layout="rectangular", size=0.2) 

### Annotate and Highlight with Mosquito Associated Lineages
S_source <- S %<+% S_metadata + 
  geom_tippoint(aes(col = host_simple), size = 0.5) + 
  scale_color_manual(values = c(cols)) + 
  #geom_text(aes(label=node, size = 0.01), size = 0.6,nudge_x = 0.01) + 
  geom_highlight(node = c(384, 467, 383, 465, 382, 466, 381, 464, 387, 
                          469, 386, 468, 385, 400, 403, 483, 402, 305, 
                          306, 309, 307), fill = "red", alpha = 0.2, extendto = max(S$data$x)) + 
  guides(col = guide_legend(ncol = 1, title = "Host", override.aes = list(size = 4))) + 
  theme(legend.position = "right", 
        legend.key.size = unit(1, "line"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 15)) + 
  geom_treescale(x = min(S$data$x), y = min(S$data$y) + 190) 

S_metadata_edit <- S_metadata %>%
  select(File_prefix, host_simple) %>%
  column_to_rownames(var = "File_prefix")

S_source2 <- S %>% gheatmap(S_metadata_edit, width = 0.05, colnames = F) + 
  scale_fill_manual(values = c(cols)) + 
  #geom_highlight(node = c(384, 467, 383, 465, 382, 466, 381, 464, 387, 
  #469, 386, 468, 385, 400, 403, 483, 402, 305, 
  #306, 309, 307), fill = "grey", alpha = 0.5, extendto = max(S$data$x)) + 
  guides(fill = guide_legend(ncol = 1, title = "Host", override.aes = list(size = 4)))
  #guides(fill = "none")
  
  
ggsave(plot = S_source2, "SerratiaModel/Serratia_marscecens/030523_InitialTreeLabels.pdf", height = 7.7*3, width = (5 + (5*plot_width))*3, units = "cm", limitsize = F)

# Useful function when plotting the heatmaps 
get_plot_width <- function(target_class){
  plot_width <- (ncol(target_class) / 50000) * 3
}

### For loop to generate separate heatmaps
heatmaps <- list()
for (i in 1:length(twilight_order)){
  class <- twilight_order[i]
  print(class)
  
  temp_pangenome <- S_pangenome_w_classification %>%
    filter(specific_class == class) %>%
    mutate(presence = as.factor(presence)) %>%
    select(!specific_class) %>%
    tidyr::pivot_wider(names_from = Gene, values_from = presence) %>%
    tibble::column_to_rownames(var = "strain")
  
  heatmaps[[i]] <- temp_pangenome
  
  plot_width <- (ncol(temp_pangenome) / 50000) * 3
  
  S_source2 %>% gheatmap(temp_pangenome,  
                        colnames = F, 
                        offset = 0.01, 
                        width = plot_width) + 
    scale_fill_manual(values = c("0" = "white", "1" = "black")) + 
    theme(legend.position = "none") 
  #theme(panel.background = element_rect(fill = "#bbbbbb", colour = NA))
  
  ggsave(glue::glue("SerratiaModel/Serratia_marscecens/140423_plot_{class}.pdf"), height = 7.7*3, width = (5 + (5*plot_width))*3, units = "cm", limitsize = F)

} # Would you like to create a new directory? - say Yes

##### Statistics
### How many genes in the collection core? (core genome)
S_pangenome_w_classification %>%
  filter(presence == 1) %>%
  group_by(specific_class, strain) %>%
  summarise(number_genes_per_strain = n()) %>%
  group_by(specific_class) %>%
  summarise(mean(number_genes_per_strain))


### Boxplot showing gene distributions
S_pangenome_w_classification %>%
  filter(presence == 1) %>%
  group_by(specific_class, strain) %>%
  summarise(number_genes_per_strain = n()) %>%
  ggplot() + 
  aes(x = specific_class, y = number_genes_per_strain) + 
  geom_boxplot() + 
  theme_bw(base_size = 20) + 
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) + 
  xlab("Sub-Classification") + 
  ylab("Number of Genes")

### Function to get gene hits per lineage
PlotPopPUNKAndGeneClassification <- function(data, category){
  data %>%
    filter(presence == 1) %>%
    group_by(specific_class, strain, Cluster) %>%
    summarise(number_genes_per_strain = n()) %>%
    filter(specific_class == category) %>%
    ggplot() + 
    aes(x = Cluster, y = number_genes_per_strain) + 
    geom_boxplot() + 
    theme_bw(base_size = 20) + 
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8)) + 
    xlab("PopPUNK Lineage") + 
    ylab("Number of Genes") + 
    geom_vline(xintercept = c(12, 24, 49, 68, 69, 144, 145, 146, 147, 148), linetype = 2, colour = "green")
}

PlotPopPUNKAndGeneClassification(S_pangenome_classfication_lineage, "Collection core")
PlotPopPUNKAndGeneClassification(S_pangenome_classfication_lineage, "Lineage specific core")
ggsave("SerratiaModel/PlotsFinal/SupplementaryFigure_LineageSpecificCoreGenes_PopPUNKCluster.pdf")


### Misc
S_pangenome_classfication_lineage %>%
  filter(presence == 1) %>%
  group_by(specific_class, strain, Cluster) %>%
  summarise(number_genes_per_strain = n()) %>%
  filter(specific_class == "Lineage specific intermediate") %>%
  filter(Cluster == 12 | Cluster == 24 | Cluster == 49 |
           Cluster == 68 | Cluster == 69 | Cluster == 144 |
           Cluster == 144 | Cluster == 145 | Cluster == 145 |
           Cluster == 146 | Cluster == 147 | Cluster == 148) %>%
  group_by(Cluster) %>%
  summarise(sum(number_genes_per_strain))



