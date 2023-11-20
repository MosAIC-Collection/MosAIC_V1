# Load Packages 
library(tidyverse)
library(ggtreeExtra)
library(ggtree)
library(ggnewscale)
library(tidyverse)
library(ape)
library(treeio)
library(MetBrewer)

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

# Read in the Enterobacter Tree 
Eliz_tree <- ape::read.tree("ElizabethkingiaModel/100423_core_gene_alignment_snp_tree.treefile")
# Root - rooting at the split defined by the genera tree
Eliz_tree_root <- phytools::midpoint.root(Eliz_tree)

# Read Panaroo Gene Presence Absence File
Eliz_panaroo_data <- read.table("ElizabethkingiaModel/260423_gene_presence_absence_clean.Rtab", 
                             header = T,
                             sep="\t",
                             quote = "",
                             comment.char = "",
                             stringsAsFactors = F,
                             check.names = F)
colnames(Eliz_panaroo_data)

# Read the Twilight Gene Classification 
Eliz_pop_class <- read_tsv("ElizabethkingiaModel/260423_Eliz_output_classification_table.tsv") %>%
  rename("Gene" = gene_name) %>%
  select(Gene, specific_class)

# Join Panaroo Gene Presence Absence + Twilight Classifications 
Eliz_pangenome_w_classification <- Eliz_panaroo_data %>%
  pivot_longer(!Gene, names_to =  "strain", values_to = "presence") %>%
  left_join(Eliz_pop_class, by =  "Gene") %>%
  print(n = 50000)

#pangenome_w_classification[["strain"]] <- gsub("\\.", "_", pangenome_w_classification[["strain"]])

# Read in PopPUNK Lineage Clusters 
Eliz_lineage_data <- read_tsv("ElizabethkingiaModel/260423_PopPUNK_clusters_refine_clean.tsv") #col_names = c("Sample", "Cluster"))  
# Convert the Cluster column from dbl to fct
Eliz_lineage_data$Cluster <- as.factor(Eliz_lineage_data$Cluster)

# Join with Panaroo Gene Presence + Absence + Twilight Classifications
Eliz_pangenome_classfication_lineage <- Eliz_pangenome_w_classification %>%
  left_join(Eliz_lineage_data, by = c("strain" =  "Taxon"))

Eliz_pangenome_classfication_lineage %>% 
  group_by(Cluster) %>% summarise(n()) %>% 
  print(n = 400)

# Read in Enterobacter asrbuiae metadata - clean + filter for isolation source 
Eliz_metadata <- read_tsv("ElizabethkingiaModel/Elizabethkingiaanophelis_metadata.tsv") %>%
  select(File_prefix, Collection, host_simple, Species)

Eliz_metadata$File_prefix <- substring(Eliz_metadata$File_prefix, 1, 15)
Eliz_metadata$File_prefix <- gsub("_contig.*", "", Eliz_metadata$File_prefix)
Eliz_metadata[["File_prefix"]] <- gsub("\\.", "_", Eliz_metadata[["File_prefix"]])

colnames(Eliz_metadata)
Eliz_metadata$File_prefix

# Establish Gals Twilight Collection Categories
twilight_order <- c("Collection core", 
                    "Multi-lineage core",
                    "Lineage specific core", 
                    #"Collection intermediate",
                    "Multi-lineage intermediate",
                    "Lineage specific intermediate",
                    #"Collection rare",
                    #"Multi-lineage rare",
                    "Lineage specific rare",
                    "Core and intermediate",
                    "Core and rare",
                    "Core, intermediate and rare",
                    "Intermediate and rare")


### Change the Tip Labels of the Tree to Match PopPUNK Custers and Panaroo Gene Presence / Absence
#### Bit messy - is there a cleaner way to do this? 
Eliz_tree_root_tib <- as_tibble(Eliz_tree_root)
Eliz_tree_root_tib$label2 <- substring(Eliz_tree_root_tib$label, 1, 15)
Eliz_tree_root_tib$label2 <- gsub("_contig.*", "", Eliz_tree_root_tib$label2)

### Substitute all "." for "_"
Eliz_tree_root_tib[["label2"]] <- gsub("\\.", "_", Eliz_tree_root_tib[["label2"]])

Eliz_tree_root_tib <- Eliz_tree_root_tib %>%
  select(label, label2)
### Rename the tip labels in the tree file 
Eliz_tree_rename <- rename_taxa(Eliz_tree_root, Eliz_tree_root_tib, label, label2)

Eliz_tree_rename$tip.label

### Draw the renamed, rooted tree
Eliz <- ggtree(Eliz_tree_rename, layout="rectangular", size=0.2) 

### Annotate and Highlight with Mosquito Associated Lineages
Eliz_source <- Eliz %<+% Eliz_metadata + 
  geom_tippoint(aes(col = host_simple), size = 0.5) + 
  scale_color_manual(values = c(Eliz_cols)) + 
  #geom_text(aes(label=node, size = 0.01), size = 0.6,nudge_x = 0.01) + 
  #geom_highlight(node = c(384, 467, 383, 465, 382, 466, 381, 464, 387, 
                          #469, 386, 468, 385, 400, 403, 483, 402, 305, 
                          #306, 309, 307), fill = "red", alpha = 0.2, extendto = max(S$data$x)) + 
  guides(col = guide_legend(ncol = 1, title = "Host", override.aes = list(size = 4))) + 
  theme(legend.position = "right", 
        legend.key.size = unit(1, "line"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 15)) + 
  geom_treescale(x = min(S$data$x), y = min(S$data$y) + 190) 

Eliz_metadata_edit <- Eliz_metadata %>%
  select(File_prefix, host_simple) %>%
  column_to_rownames(var = "File_prefix")

Eliz_source2 <- Eliz %>% gheatmap(Eliz_metadata_edit, width = 0.05, colnames = F) + 
  scale_fill_manual(values = c(Eliz_cols)) + 
  #geom_highlight(node = c(384, 467, 383, 465, 382, 466, 381, 464, 387, 
  #469, 386, 468, 385, 400, 403, 483, 402, 305, 
  #306, 309, 307), fill = "grey", alpha = 0.5, extendto = max(S$data$x)) + 
  #guides(fill = guide_legend(ncol = 1, title = "Host", override.aes = list(size = 4))) + 
  guides(fill = "none") + 
  new_scale_fill() 

ggsave(plot = Eliz_source2, "ElizabethkingiaModel/030523_InitialTreeLabels.pdf", height = 7.7*3, width = (5 + (5*plot_width))*3, units = "cm", limitsize = F)

# Useful function when plotting the heatmaps 
get_plot_width <- function(target_class){
  plot_width <- (ncol(target_class) / 50000) * 3
}

### For loop to generate separate heatmaps
heatmaps <- list()
for (i in 1:length(twilight_order)){
  class <- twilight_order[i]
  print(class)
  
  temp_pangenome <- Eliz_pangenome_w_classification %>%
    filter(specific_class == class) %>%
    mutate(presence = as.factor(presence)) %>%
    select(!c(specific_class)) %>%
    pivot_wider(names_from = Gene, values_from = presence) %>%
    column_to_rownames(var = "strain")
  
  heatmaps[[i]] <- temp_pangenome
  
  plot_width <- (ncol(temp_pangenome) / 50000) * 3
  
  Eliz_source2 %>% gheatmap(temp_pangenome, color = NULL, 
                        colnames = F, 
                        offset = 0.02, 
                        width = plot_width) + 
    scale_fill_manual(values = c("0" = "white", "1" = "black")) + 
    theme(legend.position = "none") 
  #theme(panel.background = element_rect(fill = "#bbbbbb", colour = NA))
  
  ggsave(glue::glue("ElizabethkingiaModel/260423_plot_{class}.pdf"), height = 7.7*3, width = (5 + (5*plot_width))*3, units = "cm", limitsize = F)
  
}

##### Statistics
### How many genes in the collection core? (core genome)
Eliz_pangenome_w_classification %>%
  filter(presence == 1) %>%
  group_by(specific_class, strain) %>%
  summarise(number_genes_per_strain = n()) %>%
  group_by(specific_class) %>%
  summarise(mean_genes = mean(number_genes_per_strain)) %>%
  summarise(sum(mean_genes))

Eliz_pangenome_w_classification %>%
  filter(presence == 1) %>%
  group_by(specific_class, strain) %>%
  summarise(number_genes_per_strain = n()) %>%
  group_by(specific_class) %>%
  summarise(mean(number_genes_per_strain))

write.table()

Eliz_pangenome_w_classification %>%
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
    geom_vline(xintercept = c(10, 7), linetype = 2, colour = "green") 
}

PlotPopPUNKAndGeneClassification(Eliz_pangenome_classfication_lineage, "Collection core")
PlotPopPUNKAndGeneClassification(Eliz_pangenome_classfication_lineage, "Lineage specific core")
ggsave("ElizabethkingiaModel/SupplementaryFigure_LineageSpecificCoreGenes_Elizabethkingia_anophelis_PopPUNKCluster.pdf")
PlotPopPUNKAndGeneClassification(Eliz_pangenome_classfication_lineage, "Multi-lineage core")

Eliz_pangenome_classfication_lineage %>%
  filter(presence == 1) %>%
  group_by(specific_class, strain, Cluster) %>%
  summarise(number_genes_per_strain = n()) %>%
  filter(specific_class == "Lineage specific intermediate") %>%
  filter(Cluster == 7) %>%
  group_by(Cluster) %>%
  summarise(mean(number_genes_per_strain))

