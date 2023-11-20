#### Script to plot PopPUNK Clusters 
library(ape)
library(tidyverse)
library(ggtree)
library(phytools)
library(janitor)
library(treeio)

# Load Tree
EA_tree <- read.tree("EnterobacterModel/EnterobacterAsburiae/050423_core_gene_alignment_snp_tree.treefile")
EA_tree_root <- phytools::midpoint.root(EA_tree)
EA_tree_root$tip.label <- gsub("\\.", "_", EA_tree_root$tip.label)

SM_tree <- read.tree("SerratiaModel/Serratia_marscecens/060323_core_gene_alignment_snp_tree.treefile")
SM_tree_root <- phytools::midpoint.root(SM_tree)
SM_tree_root$tip.label <- gsub("\\.", "_", SM_tree_root$tip.label)

Eliz_tree <- read.tree("ElizabethkingiaModel/100423_core_gene_alignment_snp_tree.treefile")
Eliz_tree_root <- root(Eliz_tree, outgroup = "GCA_025192605.1_ASM2519260v1_genomic.fna") 
Eliz_tree_root$tip.label <- gsub("\\.", "_", Eliz_tree_root$tip.label)

# Load PopPUNK Data
EA_PopPunkClusters <- read_tsv("EnterobacterModel/EnterobacterAsburiae/090423_PopPUNK_clusters_refine_clean.tsv")
SM_PopPunkClusters <- read_tsv("SerratiaModel/Serratia_marscecens/090423_PopPUNK_clusters_refine_clean.tsv") #col_names = c("Sample", "Cluster"))  
Eliz_PopPunkClusters <- read_tsv("ElizabethkingiaModel/260423_PopPUNK_clusters_refine_clean.tsv") %>%
  rename(Sample = "Taxon")

## Function to Make a Renamed Tree based on PopPUNK Data 
RenameTreeTipsWithPopPUNK <- function(tree_data, PopPUNK_data){
  # Fix tip Labels
  tree_data_tib <- as_tibble(tree_data) %>%
    mutate(label2 = str_sub(label, 1, 15)) %>%
    separate(col = label2, into = c("label2", "junk"), sep = "_contig", remove = F) %>%
    select(label, label2) 
  
  # Fix PopPUNK Labels
  PopPUNK_data_clean <- PopPUNK_data %>%
    mutate(Sample = str_sub(Sample, 1, 15)) %>%
    separate(col = Sample, into = c("Taxon", "junk"), sep = "_contig", remove = F) %>%
    select(Taxon, Cluster)
  
  # Join the Datasets
  tree_tib_join <- tree_data_tib %>%
    left_join(PopPUNK_data_clean, by = c("label2" = "Taxon")) %>%
    select(label, Cluster)
  
  # Rename Tips with PopPUNK Clusters
  tree_rename <- rename_taxa(tree_data, tree_tib_join, label, Cluster)  
  
  return(tree_rename)
  
}

## Function to Make a Renamed Tree based on PopPUNK Data 
EA_tree_rename <- RenameTreeTipsWithPopPUNK(EA_tree_root, EA_PopPunkClusters)
SM_tree_rename <- RenameTreeTipsWithPopPUNK(SM_tree_root, SM_PopPunkClusters)
Eliz_tree_rename <- RenameTreeTipsWithPopPUNK(Eliz_tree_root, Eliz_PopPunkClusters)

# plot + save
ggtree(EA_tree_rename, size = 0.2, layout = "circular") + geom_tiplab(offset = 0.001, size = 2)
ggsave(filename = "EnterobacterModel/EnterobacterAsburiae/260423_EA_PopPUNKCluster_Tree.pdf")

ggtree(SM_tree_rename, size = 0.2, layout = "circular") + geom_tiplab(offset = 0.001, size = 1)
ggsave(filename = "SerratiaModel/Serratia_marscecens/260423_SM_PopPUNKCluster_Tree.pdf")

ggtree(Eliz_tree_rename, size = 0.2, layout = "circular") + geom_tiplab(offset = 0.001, size = 1)
ggsave(filename = "ElizabethkingiaModel/260423_Eliz_PopPUNKCluster_Tree.pdf")

#### Stats 
#### How many lineages defined by popPUNK? 
EA_PopPunkClusters %>%
  distinct(Cluster) %>%
  summarise(n()) # 45 Clusters 

SM_PopPunkClusters %>%
  distinct(Cluster) %>%
  summarise(n()) # 153 Clusters 

Eliz_PopPunkClusters %>%
  distinct(Cluster) %>%
  summarise(n()) # 141 Clusters

