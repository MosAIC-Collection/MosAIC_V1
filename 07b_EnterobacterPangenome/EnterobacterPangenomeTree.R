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

# Read in the Enterobacter Tree 
enterobacter_tree <- ape::read.tree("MosAIC_V1/07b_EnterobacterPangenome/050423_core_gene_alignment_snp_tree.treefile")

# Root - rooting at the split defined by the genera tree
enterobacter_tree_root <- phytools::midpoint.root(enterobacter_tree)

# Read Panaroo Gene Presence Absence File
E_panaroo_data <- read.table("MosAIC_V1/07b_EnterobacterPangenome/090423_gene_presence_absence_clean.Rtab", 
                           header = T,
                           sep="\t",
                           quote = "",
                           comment.char = "",
                           stringsAsFactors = F,
                           check.names = F)

# Read the Twilight Gene Classification 
E_pop_class <- read_tsv("MosAIC_V1/07b_EnterobacterPangenome/output_classification_table.tab") %>%
  rename("Gene" = gene_name) %>%
  select(Gene, specific_class)

# Join Panaroo Gene Presence Absence + Twilight Classifications 
E_pangenome_w_classification <- E_panaroo_data %>%
  pivot_longer(!Gene, names_to =  "strain", values_to = "presence") %>% 
  left_join(E_pop_class, by =  "Gene") 

#pangenome_w_classification[["strain"]] <- gsub("\\.", "_", pangenome_w_classification[["strain"]])

# Read in PopPUNK Lineage Clusters 
E_lineage_data <- read_tsv("MosAIC_V1/07b_EnterobacterPangenome/090423_PopPUNK_clusters_refine_clean.tsv") #col_names = c("Sample", "Cluster"))  

# Convert the Cluster column from dbl to fct
E_lineage_data$Cluster <- as.factor(E_lineage_data$Cluster)

# Join with Panaroo Gene Presence + Absence + Twilight Classifications
E_pangenome_classfication_lineage <- E_pangenome_w_classification %>%
  left_join(E_lineage_data, by = c("strain" =  "Sample"))

# Read in Enterobacter asrbuiae metadata - clean + filter for isolation source 
E_metadata <- read_tsv("MosAIC_V1/07b_EnterobacterPangenome/Chavda_Mos_Curated_Metadata_2.txt") %>%
  clean_names() 

# Substitute all "." with "_"  
#E_metadata_final[["accession"]] <- gsub("\\.", "_", E_metadata_final[["accession"]])
E_metadata[["assembly"]] <- gsub("\\.", "_", E_metadata[["assembly"]])

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


# Look at Tip labels - they don't match the Panaroo + Twilight Labels - will cause problems when annotating 
enterobacter_tree$tip.label

### Change the Tip Labels of the Tree to Match PopPUNK Custers and Panaroo Gene Presence / Absence
#### Bit messy - is there a cleaner way to do this? 
enterobacter_tree_tib <- as_tibble(enterobacter_tree) %>%
  separate(col = label, into = c("label2", "junk"), sep = "_ASM", remove = F) %>%
  separate(col = label2, into = c("label2", "junk"), sep = "_contigs", remove = F) %>%
  separate(col = label2, into = c("label2", "junk"), sep = "_120", remove = F) %>%
  separate(col = label2, into = c("label2", "junk"), sep = "_254", remove = F) %>%
  separate(col = label2, into = c("label2", "junk"), sep = "_119", remove = F) %>%
  separate(col = label2, into = c("label2", "junk"), sep = "_L1_gen", remove = F) %>%
  separate(col = label2, into = c("label2", "junk"), sep = "_De_novo", remove = F) %>%
  separate(col = label2, into = c("label2", "junk"), sep = "_DC1", remove = F) %>%
  separate(col = label2, into = c("label2", "junk"), sep = "_UHGG", remove = F) %>%
  separate(col = label2, into = c("label2", "junk"), sep = "_Ent", remove = F) %>%
  select(label, label2)
  
### Substitute all "." for "_"
enterobacter_tree_tib[["label2"]] <- gsub("\\.", "_", enterobacter_tree_tib[["label2"]])

### Rename the tip labels in the tree file 
enterobacter_tree_rename <- rename_taxa(enterobacter_tree_root, enterobacter_tree_tib, label, label2)

enterobacter_tree_rename$tip.label

### Draw the renamed, rooted tree
E <- ggtree(enterobacter_tree_rename, layout="rectangular", size=0.2)

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


MosAIC_cols <- c("female_mosquito" = "#4e79a7",
                 "male_mosquito" = "#f28e2b",
                 "larval_water" = "#e15759",
                 "mosquito" = "#C3D673",
                 "mosquito_egg" = "#76b7b2",
                 "stable fly" = "#394B18",
                 "water" = "#336220",
                 "Coon" = "#edc948", 
                 "Hughes/Heinz" = "#297A31",
                 "Lampe, David" = "#32915E",
                 "Povelones, Michael" = "#1BC200",
                 "Valiente Moro, Claire" = "#1FA918",
                 "field" = "#344FAA",
                 "Ankazobe, Madagascar" = "#2760AA",
                 "E. hormaechei subsp. hoffmannii" = "#1B77A9",
                 "E. hormaechei subsp. oharae" =  "#1191A7",
                 "E. hormaechei subsp. steigerwaltii" = "#07A49A",
                 "lab" = "#913431",
                 "field" = "Black", 
                 "lab" = "grey",
                 "Arlington, WI" = "#90318C",
                 "Cook Co, IL" = "#827A87",
                 "LSTM" = "#76b7b2",
                 "Tamatave, Madagascar" = "#B88694",
                 "UGA" = "black",
                 "UW Madison" = "white", 
                 "UW Madison" = "#F781BF")


E_metadata %>%
  group_by(location) %>%
  summarise(n())

# Plot Initial Tree with Host Metadata + Highlight the Mosquito-Associated Clades
E_source <- E %<+% E_metadata + 
  geom_tippoint(aes(col = host_simple), size = 1) +
  geom_text(aes(label=node, size = 0.2), size = 2,nudge_x = 0.01) 
  #geom_highlight(node = c(100, 132, 99, 131, 101, 130, 98, 128, 
  #97, 129, 96, 123, 100, 126, 105, 125, 104, 124, 102, 92, 140, 
  #91, 139, 122, 98, 188, 96, 23, 157, 156, 
  #62, 211, 61, 210, 60), fill = "red", alpha = 0.2) + 
  scale_color_manual(values = c(Chavda_cols)) + 
  theme(legend.position = "right", 
        legend.key.size = unit(1, "line"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 15)) + 
  guides(col = guide_legend(ncol = 1, title = "Host", override.aes = list(size = 4))) + 
  geom_treescale(x = min(E$data$x), y = min(E$data$y) + 30) 
  
E_metadata_edit <- E_metadata %>%
  select(assembly, host_simple) %>%
  distinct(assembly, .keep_all = T) %>%
  filter(!is.na(assembly)) %>%
  column_to_rownames(var = "assembly")

E_source2 <- E %>% gheatmap(E_metadata_edit, width = 0.05, colnames = F) + 
  scale_fill_manual(values = c(Chavda_cols)) + 
  #geom_highlight(node = c(100, 132, 99, 131, 101, 130, 98, 128, 
  #97, 129, 96, 123, 100, 126, 105, 125, 104, 124, 102, 92, 140, 
  #91, 139, 122, 98, 188, 96, 23, 157, 156, 
  #62, 211, 61, 60), fill = "grey", alpha = 0.25, extendto = max(E$data$x)) +  
  #guides(fill = guide_legend(ncol = 1, title = "Host", override.aes = list(size = 4)))
  guides(fill = "none")
  

ggsave(plot = E_source2, "EnterobacterModel/EnterobacterAsburiae/030523_InitialTreeLabels.pdf", height = 7.7*3, width = (5 + (5*plot_width))*3, units = "cm", limitsize = F)


# Useful function when plotting the heatmaps 
get_plot_width <- function(target_class){
  plot_width <- (ncol(target_class) / 50000) * 3
}

### For loop to generate separate heatmaps
heatmaps <- list()
for (i in 1:length(twilight_order)){
  class <- twilight_order[i]
  print(class)
  
  temp_pangenome <- E_pangenome_w_classification %>%
    filter(specific_class == class) %>%
    mutate(presence = as.factor(presence)) %>%
    select(!specific_class) %>%
    tidyr::pivot_wider(names_from = Gene, values_from = presence) %>%
    tibble::column_to_rownames(var = "strain")
  
  heatmaps[[i]] <- temp_pangenome
  
  plot_width <- (ncol(temp_pangenome) / 50000) * 3
  
  E %>% gheatmap(temp_pangenome, color = NULL, 
                  colnames = F, 
                  offset = 0.01, 
                  width = plot_width) + 
    scale_fill_manual(values = c("0" = "white", "1" = "black")) + 
    theme(legend.position = "none") 
  #theme(panel.background = element_rect(fill = "#bbbbbb", colour = NA))
  
  ggsave(glue::glue("EnterobacterModel/EnterobacterAsburiae/140423_plot_{class}.pdf"), height = 7.7*3, width = (5 + (5*plot_width))*3, units = "cm", limitsize = F)
  
  
}

#### Make a more "honed" in version appending more of the metadata 
E_sub <- ggtree(enterobacter_tree_rename, size = 0.2, layout = "fan", open.angle = 230) 

E_metadata_edit2 <- as.data.frame(E_metadata %>% 
                                    select(!gan_bank_accession_s) %>% 
                                    select(assembly, host_simple, mos_aic_source_specific, 
                                           mos_aic_source_lab, location, mos_aic_lab_field, 
                                           mosquito_spp) %>%
                                    rename(ID = "assembly") %>%
                                    mutate(host = 1))

E_sub2 <- E_sub + geom_fruit(data = E_metadata_edit2, geom=geom_tile, 
                   mapping = aes(y = ID, x= host, fill = host_simple), 
                   size = 0.02, width = 0.01, offset = -0.15) + 
  scale_fill_manual(values = cols) + new_scale_fill() + 
  geom_highlight(node = c(100, 132, 99, 131, 101, 130, 98, 128, 
                          97, 129, 96, 123, 100, 126, 105, 125, 104, 124, 102, 92, 140, 
                          91, 139, 122, 98, 96, 23, 157, 156, 
                          62, 211, 61, 210, 60), fill = "grey",
                 size = 0.05, alpha = 0.1, extendto = 0.35)

E_sub3 <- E_sub2 + geom_fruit(
  data = E_metadata_edit2, 
  geom = geom_star, 
  mapping = aes(x = host, y = ID, starshape=mos_aic_source_specific),
  offset = -0.1,
  grid.params=list(
    linetype=3,
    size=0.3
  )) + scale_fill_manual(values = MosAIC_cols,
                         name = "Mosquito Isolation Specific",
                         guide=guide_legend(
                           keywidth=0.3, 
                           keyheight=0.3, 
                           order=3
  ),
  na.translate=FALSE) + 
  new_scale_fill() 

E_sub4 <- E_sub3 + geom_fruit(
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

E_sub5 <- E_sub4 + geom_fruit(
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
  new_scale_fill()

ggsave(plot = E_sub5,filename =  "EnterobacterModel/EnterobacterAsburiae/204023_Enterobacter_asburiae_subset_metadata.pdf")


# Alternatively - put on one plot 
# The Y-axis alignment might not be right however 
# Make Tree

#### Functiont to extract the Twilight classifications of each category
extract_gene_classification <- function(data, column){
  target <- data %>%
    filter(specific_class == column) %>%
    mutate(presence = as.factor(presence)) %>%
    pivot_wider(names_from = Gene, values_from = presence) %>%
    column_to_rownames(var = "strain")
  
  return(target)
}

collection_core <- extract_gene_classification(E_pangenome_w_classification, "Collection core")
multi_lineage_core <- extract_gene_classification(E_pangenome_w_classification, "Multi-lineage core")
lineage_specific_core <- extract_gene_classification(E_pangenome_w_classification, "Lineage specific core")

## Nothing in Collection intermediate - ignore
#collection_intermediate <- extract_gene_classification(E_pangenome_w_classification, "Collection intermediate")
multi_lineage_intermediate <- extract_gene_classification(E_pangenome_w_classification, "Multi-lineage intermediate")
lineage_specific_intermediate <- extract_gene_classification(E_pangenome_w_classification, "Lineage specific intermediate")

#collection_rare <- extract_gene_classification(E_pangenome_w_classification, "Collection rare")
multie_lineage_rare <- extract_gene_classification(E_pangenome_w_classification, "Multi-lineage rare")
lineage_specific_rare <- extract_gene_classification(E_pangenome_w_classification, "Lineage specific rare")

core_and_intermediate <- extract_gene_classification(E_pangenome_w_classification, "Core and intermediate")
core_and_rare <- extract_gene_classification(E_pangenome_w_classification, "Core and rare")
core_intermediate_rare <- extract_gene_classification(E_pangenome_w_classification, "Core, intermediate and rare")
intermediate_rare <- extract_gene_classification(E_pangenome_w_classification, "Intermediate and rare")

E2 <- gheatmap(E_Source, collection_core, width = get_plot_width(collection_core), 
               colnames = FALSE, color = NULL, 
               colnames_offset_y = -100, low = "white", high = "black", offset = 0.01) + 
  scale_fill_manual(values = c("0" = "white", "1" = "#3B7844"))

E3 <- E2 + new_scale_fill()

E4 <- gheatmap(E3, multi_lineage_core, width = get_plot_width(multi_lineage_core), 
               colnames = FALSE, color = NULL, 
               colnames_offset_y = -100, low = "white", high = "black", offset = 0.03) + 
  scale_fill_manual(values = c("0" = "white", "1" = "#9ACD6C"))

E5 <- E4 + new_scale_fill()

E6 <- gheatmap(E5, lineage_specific_core, width = get_plot_width(lineage_specific_core), 
               colnames = FALSE, color = NULL, 
               colnames_offset_y = -10, low = "white", high = "black", offset = 0.04) + 
  scale_fill_manual(values = c("0" = "white", "1" = "#267764"))

E7 <- E6 + new_scale_fill()

collection_intermediate

E8 <- gheatmap(E7, multi_lineage_intermediate, width = get_plot_width(multi_lineage_intermediate), 
               colnames = FALSE, color = NULL, 
               colnames_offset_y = -10, low = "white", high = "black", offset = 0.05) + 
  scale_fill_manual(values = c("0" = "white", "1" = "#6A1A00"))

E9 <- E8 + new_scale_fill()

E10 <- gheatmap(E9, lineage_specific_intermediate, width = get_plot_width(lineage_specific_intermediate), 
               colnames = FALSE, color = NULL, 
               colnames_offset_y = -10, low = "white", high = "black", offset = 0.06) + 
  scale_fill_manual(values = c("0" = "white", "1" = "#D12AAE"))

E11 <- E10 + new_scale_fill()

E12 <- gheatmap(E11, multie_lineage_rare, width = get_plot_width(multie_lineage_rare), 
  colnames = FALSE, color = NULL, 
  colnames_offset_y = -10, low = "white", high = "black", offset = 0.09) + 
  scale_fill_manual(values = c("0" = "white", "1" = "#7573FF"))

E13 <- E12 + new_scale_fill()

E14 <- gheatmap(E13, lineage_specific_rare, width = get_plot_width(lineage_specific_rare), 
                colnames = FALSE, color = NULL, 
                colnames_offset_y = -10, low = "white", high = "black", offset = 0.095) + 
  scale_fill_manual(values = c("0" = "white", "1" = "#A7D9FF"))

E15 <- E14 + new_scale_fill()

E16 <- gheatmap(E15, core_and_intermediate, width = get_plot_width(core_and_intermediate), 
                colnames = FALSE, color = NULL, 
                colnames_offset_y = -10, low = "white", high = "black", offset = 0.11) + 
  scale_fill_manual(values = c("0" = "white", "1" = "#6B512D"))

E17 <- E16 + new_scale_fill()

E18 <- gheatmap(E17, core_and_rare, width = get_plot_width(core_and_rare), 
                colnames = FALSE, color = NULL, 
                colnames_offset_y = -10, low = "white", high = "black", offset = 0.13) + 
  scale_fill_manual(values = c("0" = "white", "1" = "#3A853F"))

E19 <- E18 + new_scale_fill()

E20 <- gheatmap(E19, core_intermediate_rare, width = get_plot_width(core_intermediate_rare), 
                colnames = FALSE, color = NULL, 
                colnames_offset_y = -10, low = "white", high = "black", offset = 0.14) + 
  scale_fill_manual(values = c("0" = "white", "1" = "#476C9E"))

E21 <- E20 + new_scale_fill()

E22 <- gheatmap(E21, intermediate_rare, width = get_plot_width(intermediate_rare), 
                colnames = FALSE, color = NULL, 
                colnames_offset_y = -10, low = "white", high = "black", offset = 0.15) + 
  scale_fill_manual(values = c("0" = "white", "1" = "#B754B0"))

ggsave(filename = "EnterobacterModel/170223_EnterobacterPangenomeStructure.pdf", E22)

#### Statistics 
E_pangenome_w_classification %>%
  filter(presence == 1) %>%
  group_by(specific_class, strain) %>%
  summarise(number_genes_per_strain = n()) %>%
  group_by(specific_class) %>%
  summarise(mean_genes = mean(number_genes_per_strain)) %>%
  summarise(sum(mean_genes))

write.table()

E_pangenome_w_classification %>%
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
    theme(axis.text.x = element_text(angle = 40, hjust = 1)) + 
    xlab("PopPUNK Lineage") + 
    ylab("Number of Genes") + 
    geom_vline(xintercept = c(1, 7, 8), linetype = 2, colour = "green")
}

PlotPopPUNKAndGeneClassification(E_pangenome_classfication_lineage, "Collection core")
PlotPopPUNKAndGeneClassification(E_pangenome_classfication_lineage, "Lineage specific core")
ggsave("EnterobacterModel/PlotsFinal/SupplementaryFigure_LineageSpecificCoreGenes_PopPUNKCluster.pdf")
PlotPopPUNKAndGeneClassification(E_pangenome_classfication_lineage, "Multi-lineage core")

# How many lineage-specific core genes in mosquito-associated lineages? 
E_pangenome_classfication_lineage %>%
  filter(presence == 1) %>%
  group_by(specific_class, strain, Cluster) %>%
  summarise(number_genes_per_strain = n()) %>%
  filter(specific_class == "Lineage specific core") %>%
  filter(Cluster == 1 | Cluster == 7 | Cluster == 8) %>%
  group_by(Cluster) %>%
  summarise(mean(number_genes_per_strain))

# How many lineage-specific core genes in mosquito-associated lineages? 
E_pangenome_classfication_lineage %>%
  filter(presence == 1) %>%
  group_by(specific_class, strain, Cluster) %>%
  summarise(number_genes_per_strain = n()) %>%
  filter(specific_class == "Lineage specific intermediate") %>%
  filter(Cluster == 1 | Cluster == 7 | Cluster == 8) %>%
  group_by(Cluster) %>%
  summarise(mean(number_genes_per_strain))

E_Full_Class <- read_table("EnterobacterModel/120223_Enterobacter_classification.tab") %>%
  select(gene_name, total, specific_class)

E_Full_Class %>%
  group_by(total, core, inter, rare, total) %>%
  summarise(n())

E_pangenome_w_classification %>%
  filter(presence >= 1) 
  
E_pangenome_classfication_lineage %>%
  filter(presence >= 1)

E_pangenome_classfication_lineage %>%
  filter(presence >= 1) %>%
  filter(general_class != "Absent in large lineages") %>%
  group_by(Cluster, specific_class, general_class) %>%
  summarise(number_genes = n()) %>%
  ggplot() +
  aes(x = Cluster, y = number_genes, fill = specific_class) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(values = met.brewer("Renoir", 12)) + 
  theme_bw(base_size = 10) + 
  xlab("Cluster") + 
  ylab("Number of Genes") + 
  facet_wrap(~general_class, scales = "free_x") 
  
### 
# MISC
### group the lineages and check if Collection Core does actually represent >95% Frequency
E_pangenome_classfication_lineage %>%
  filter(presence >= 1) %>%
  filter(specific_class == "Collection core") %>%
  group_by(Cluster, strain) %>%
  summarise(n())


E_pangenome_classfication_lineage %>%
  filter(Cluster == 1) %>%
  group_by(specific_class) %>%
  summarise(number_genes = n()) %>%
  ggplot() +
  aes(x = specific_class, y = number_genes) + 
  geom_bar(stat = "identity", position = "dodge") + 
  theme_bw(base_size = 20) + 
  xlab("Cluster") + 
  ylab("Number of Genes") + 
  theme(axis.text.x = element_text(angle = 40, hjust = 1))


ggplot(E_pangenome_classfication_lineage, aes(x = Cluster, fill = specific_class)) + 
  geom_bar(color = "black", lwd = 0.2) +
  facet_grid(general_class~., scales = "free") + 
  xlab("") + theme_classic(base_size = 12) + ylab("Number of genes") +
  scale_fill_manual(values = met.brewer("Renoir", 12)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Number of lineages in which gene is present") 

ggsave(A, filename = file.path(plots_out, "barplot_per_lineage_count.pdf"), height = 10, width= 10)

## Distributions per Cluster 

gene_per_isoalte <- read_tsv("EnterobacterModel/genes_per_isolate.tsv") #col_names = c("cluster", "class", "collection", "count", "count2"))
gene_per_isoalte$class <- as.factor(gene_per_isoalte$class)

gene_per_isoalte %>%
  filter(specific_classification == "Multi-lineage intermediate") %>%
  ggplot() + 
  aes(x = class, y = count, group = class) + 
  geom_boxplot() + 
  theme_bw(base_size = 20) + 
  xlab("PopPUNK Cluster") + 
  ylab("Number of Genes") + 
  ggtitle("Lineage Specific Core Genes Per PopPUNK Cluster")


## Loop and Generate PDFs across all classes 
#plots <- list()
#for (i in 1:length(twilight_order)){
 # target <- twilight_order[i]
  #print(target)
  
  #filter <- gene_per_isoalte %>%
   # filter(specific_classification == target)
  
  #plot <- filter[[i]]
  
 # filter %>%
   # ggplot() + 
    #aes(x = class, y = count, group = class) + 
    #geom_boxplot() + 
    #theme_bw(base_size = 20) + 
    #xlab("PopPUNK Cluster") + 
    #ylab("Number of Genes") + 
    #ggtitle("{class} Per PopPUNK Cluster")

  #ggsave(glue::glue("EnterobacterModel/{class}_per_Cluster_Distribution.pdf"), units = "cm", limitsize = F)
#}

gene_per_isoalte %>%
  filter(specific_classification == "Lineage specific core") %>%
  ggplot() + 
  aes(x = class, y = count, group = class) + 
  geom_boxplot() + 
  theme_bw(base_size = 20) + 
  xlab("PopPUNK Cluster") + 
  ylab("Number of Genes") + 
  ggtitle("Lineage Specific Genes Per PopPUNK Cluster")

ggsave("EnterobacterModel/140223_lineage_specific_genes_per_cluster.pdf")

gene_per_isoalte %>%
  filter(specific_classification == "Multi-lineage core") %>%
  ggplot() + 
  aes(x = class, y = count, group = class) + 
  geom_boxplot() + 
  theme_bw(base_size = 20) + 
  xlab("PopPUNK Cluster") + 
  ylab("Number of Genes") + 
  ggtitle("Lineage Specific Genes Per PopPUNK Cluster")

ggsave("EnterobacterModel/140223_lineage_specific_genes_per_cluster.pdf")



plot_gene_distribution <- function(table, filter){
  table %>%
    filter(specific_classification == {{filter}}) %>%
    ggplot() + 
    aes(x = class, y = count, group = class) + 
    geom_boxplot() + 
    theme_bw(base_size = 20) + 
    xlab("PopPUNK Cluster") + 
    ylab("Number of Genes") + 
    ggtitle(filter + "Per PopPUNK Cluster")
}

plot_gene_distribution(gene_per_isoalte, "Collection core")

gene_per_isoalte %>%
  group_by(specific_classification) %>%
  summarise(n())
