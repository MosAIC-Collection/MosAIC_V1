### 
# Script to show that Elizabethkingia meningoseptica do not fall within the Elizabethkingia anophelis clades
###

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
Elizabethkingia_Pop <- read.tree("ElizabethkingiaModel/250823_core_gene_alignment_snp_filt.aln.treefile")
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
#Rooted_Elizabethkingia <- phytools::midpoint.root(Elizabethkingia_Pop_rename) 
Rooted_Elizabethkingia <- root(Elizabethkingia_Pop_rename, outgroup = "GCA_900460995_1") # root this at Chryseobacterium outgroup 
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
  layout_fan(angle = 90) + geom_tiplab(size = 0.5)
## Chryseobacterium forms an outgroup and the 7 E. meningoseptica form a deep branching clade 

# Save tree with tip labels - use this to root the tree without the Chrysoebacterium outgroup for display purposes
ggsave(filename = "ElizabethkingiaModel/280823_ElizabethkingiaModel_tree.pdf")
