### Packages 
library(TDbook)
library(ggstar)
library(ggtreeExtra)
library(ggtree)
library(ggnewscale)
library(ggstar)
library(tidyverse)
library(treeio)
library(MetBrewer)

### VFs
MosAIC_VF_clean <- read_tsv("210423_MosAIC_VF_clean.tsv")
MosAIC_VF_clean$number_file <- gsub("_contig.*", "", MosAIC_VF_clean$number_file)

### dRep Clusters
dRep_clusters <- readr::read_tsv("Cdb.txt")
dRep_clusters$genome <- gsub("_contig.*", "", dRep_clusters$genome)

dRep_winning_clusters <- readr::read_csv("Wdb.csv")
dRep_winning_clusters$genome <- gsub("_contig.*", "", dRep_winning_clusters$genome)

### GTDB 
GTDB_Class <- read_tsv("291122_gtdbtk_isolates_drep.bac120.summary.tsv")
GTDB_Class$user_genome <- gsub("_contig.*", "", GTDB_Class$user_genome)

GTDB_Full <- read_tsv("190123_gtdb.bac120.summary_nodrep.tsv")
GTDB_Full$user_genome <- gsub("_contig.*", "", GTDB_Full$user_genome)

GTDB_Class_Clean <- GTDB_Class %>%
  select(user_genome, classification)

GTDB_Class_Clean <- GTDB_Class %>%
  separate(classification, into = c("domain", "parent", "class", "order", "family", "genus", "species"), sep = ";") %>% 
  separate(species, into = c("junk", "species"), sep = "s__") %>%
  separate(family, into = c("junk2", "family"), sep = "f__") %>%
  separate(genus, into = c("junk3", "genus"), sep = "g__") %>%
  separate(domain, into = c("junk4", "domain"), sep = "d__") %>%
  separate(parent, into = c("junk5", "parent"), sep = "p__") %>%
  separate(class, into = c("junk6", "class"), sep = "c__") %>%
  separate(order, into = c("junk7", "order"), sep = "o__") %>%
  select(!c(junk, junk2, junk3, junk4, junk5, junk6, junk7)) %>%
  select(user_genome, family)

GTDB_Full_Clean <- GTDB_Full %>%
  separate(classification, into = c("domain", "parent", "class", "order", "family", "genus", "species"), sep = ";") %>% 
  separate(species, into = c("junk", "species"), sep = "s__") %>%
  separate(family, into = c("junk2", "family"), sep = "f__") %>%
  separate(genus, into = c("junk3", "genus"), sep = "g__") %>%
  separate(domain, into = c("junk4", "domain"), sep = "d__") %>%
  separate(parent, into = c("junk5", "parent"), sep = "p__") %>%
  separate(class, into = c("junk6", "class"), sep = "c__") %>%
  separate(order, into = c("junk7", "order"), sep = "o__") %>%
  select(!c(junk, junk2, junk3, junk4, junk5, junk6, junk7)) %>%
  select(user_genome, genus)

GTDB_Full_Clean <- GTDB_Full %>%
  separate(classification, into = c("domain", "parent", "class", "order", "family", "genus", "species"), sep = ";") %>% 
  separate(species, into = c("junk", "species"), sep = "s__") %>%
  separate(family, into = c("junk2", "family"), sep = "f__") %>%
  separate(genus, into = c("junk3", "genus"), sep = "g__") %>%
  separate(domain, into = c("junk4", "domain"), sep = "d__") %>%
  separate(parent, into = c("junk5", "parent"), sep = "p__") %>%
  separate(class, into = c("junk6", "class"), sep = "c__") %>%
  separate(order, into = c("junk7", "order"), sep = "o__") %>%
  select(!c(junk, junk2, junk3, junk4, junk5, junk6, junk7))

### Tree 
MosAIC_Tree <- ape::read.tree("MMO142ADJ_240223_concatenated_16s_fileappend_final_aln_trimal.fa.contree")
MosAIC_Tree <- phytools::midpoint.root(MosAIC_Tree)

MosAIC_Tree_Tib <- as_tibble(MosAIC_Tree)
MosAIC_Tree_Tib$label2 <- gsub("_contig.*", "", MosAIC_Tree_Tib$label)

MosAIC_Tree_Tib <- MosAIC_Tree_Tib %>%
  select(label, label2)

### Rename the tip labels in the tree file 
MosAIC_Tree_rename <- rename_taxa(MosAIC_Tree, MosAIC_Tree_Tib, label, label2)

### Join 
MosAIC_VF_clean_grouped <- MosAIC_VF_clean %>% 
  group_by(number_file, category) %>%
  summarise(num_associated_genes = n())

### Join this with dRep clusters 
MosAIC_VF_dRep <- MosAIC_VF_clean_grouped %>%
  left_join(dRep_clusters, by = c("number_file" = "genome")) %>%
  select(number_file, category, num_associated_genes, secondary_cluster, primary_cluster)

MosAIC_VF_Winning_Cluster <- MosAIC_VF_dRep %>%
  left_join(dRep_winning_clusters, by = c("secondary_cluster" = "cluster"), keep = T) %>%
  select(number_file, category, num_associated_genes, cluster)

MosAIC_VF_Winning_Cluster_Mean_VF <- MosAIC_VF_Winning_Cluster %>%
  group_by(cluster, category) %>%
  summarise(mean_VF = mean(num_associated_genes))

MosAIC_VF_Winning_Cluster_Mean_VF_Samples <- MosAIC_VF_Winning_Cluster_Mean_VF %>%
  left_join(dRep_winning_clusters) %>%
  select(genome, category, mean_VF, cluster)

extract_column_and_convert <- function(data, column){
  target <- data %>%
    ungroup() %>%
    filter(category == "Adherence") %>%
    select(!c(category, cluster)) %>%
    column_to_rownames(var = "genome")
  
  return(target)
}

MosAIC_VF_Winning_Cluster_Mean_VF_Samples_wide <- MosAIC_VF_Winning_Cluster_Mean_VF_Samples %>%
  pivot_wider(names_from = category, values_from = mean_VF) %>%
  ungroup() %>%
  select(!cluster) %>%
  column_to_rownames(var = "genome")

GTDB_Class_Clean_Data <- GTDB_Class_Clean %>%
  column_to_rownames(var = "user_genome")

MosAIC_VF_sum_per_clust <- MosAIC_VF_Winning_Cluster_Mean_VF_Samples %>%
  group_by(genome) %>%
  summarise(mean_VF_per_clust = sum(mean_VF)) %>%
  rename(ID = genome)

MosAIC_VF_Winning_Cluster_Mean_VF_Samples_wide

## Plot Clusters in Tree 
p <- ggtree(MosAIC_Tree_rename, layout="rectangular", size=0.2, open.angle=270, xlim = 1)

p2 <- gheatmap(p, GTDB_Class_Clean_Data, width = 0.05, offset = 0.0, colnames_position = "top", colnames = FALSE, legend_title = "Family") + 
  scale_fill_manual(values = c("#005C66", "#BD5129", "#9FBAE7", "#1000D2", "#0070C8", "#52CEFF", 
                               "#009464", "#E1B698", "#06BF34", "#9D49AB", "#94E4FF", "#E3F8FF", 
                               "#A2581B", "#00B9C8", "#9AC7D9", "#50E4BC", "#6E2A4D", "#D582D2", 
                               "#F0AC7A", "#F7D7EF", "#97ADFF", "#7500FF", "#FFF4B1", "#52E685", 
                               "#FF9146", "#E0FFF1", "#A3A3FF", "#FFAAB7", "#FFE7CC")) + 
  guides(fill = guide_legend(ncol = 3, title = "Family", override.aes = list(size = 4))) + 
  new_scale_fill()

p3 <- gheatmap(p2, MosAIC_VF_Winning_Cluster_Mean_VF_Samples_wide, offset = 0.02,
         colnames_position = "top",
         colnames_angle = 90,
         font.size = 5, 
         hjust = 0,
         colnames_offset_y = 3,
         legend_title = "Family") + 
  scale_fill_gradient2(
    low = "#84a98c",
    mid = "#52796f",
    high = "#ffd166",
    midpoint = 30,
    space = "Lab",
    na.value = "#e9ecef",
    guide = "colourbar",
    aesthetics = "fill"
  ) + 
  new_scale_fill() 
  

p4 <- p3 + ylim(0, 300)



p5 <- p4 +
  geom_fruit(data=MosAIC_VF_sum_per_clust, geom=geom_bar,
             mapping=aes(y=ID, x=mean_VF_per_clust),
             pwidth=0.5, 
             orientation="y", 
             stat="identity",
             offset = 1.09,
             axis.params=list(
               axis       = "x",
               text.size  = 5,
               hjust      = 1,
               vjust      = 1.5,
               nbreak     = 5,
             ),grid.params = list(size = 0.5),
  ) +
  #scale_fill_manual(values=c("#0000FF","#FFA500","#FF0000",
                             #"#800000", "#006400","#800080","#696969"),
                    #guide=guide_legend(keywidth = 0.3, 
                                      # keyheight = 0.3, order=4, ncol = 1))+
  geom_treescale(fontsize=2, linesize=0.3) +
  theme(legend.background=element_rect(fill=NA),
        legend.title=element_text(size=6.5),
        legend.text=element_text(size=10),
        legend.spacing.y = unit(0.02, "cm"),
  )


ggsave(p5, filename = "180523_MosAIC_VF_Heatmap.pdf", width = 20, height = 15)


#### Stats 
MosAIC_VF_clean %>%
  select(Gene) # 11,774 genes in the collection 

MosAIC_VF_GTDB <- MosAIC_VF_clean %>%
  left_join(GTDB_Full_Clean, by = c("number_file" = "user_genome")) %>%
  select(number_file, category, Gene, class, family, genus)

#### Total number of genes per class 
MosAIC_VF_GTDB %>%
  group_by(class) %>%
  summarise(nunber_virulence_factor_genes = n()) 
  
MosAIC_VF_GTDB %>%
  group_by(class, category, number_file) %>%
  summarise(nunber_virulence_factor_genes = n()) %>%
  ungroup() %>%
  group_by(class, category) %>%
  summarise(mean_num_VF_genes = mean(nunber_virulence_factor_genes)) %>%
  ggplot() + 
  aes(x = category, y = mean_num_VF_genes, fill = class) + 
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(values=c("#B999CC","#1D9A6C","#002ACD","#B2925C",
                               "#B64521", "#9ACD32","#D15FEE","#FFC0CB",
                               "#EE6A50","#8DEEEE", "#006400","#800000",
                               "#B0171F","#191970")) + 
  theme_bw(base_size = 20) + 
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) + 
  coord_flip()

MosAIC_VF_GTDB %>%
  group_by(class, category, number_file) %>%
  summarise(nunber_virulence_factor_genes = n()) %>%
  ungroup() %>%
  group_by(class, category) %>%
  summarise(mean_num_VF_genes = mean(nunber_virulence_factor_genes)) %>%
  print(n = 27)

#### Which class has the most virulence factors 
MosAIC_VF_GTDB %>%
  group_by(class, number_file) %>%
  summarise(number_genes = n()) %>%
  group_by(class) %>%
  summarise(mean_number_genes = mean(number_genes))

#### Uniquely present genes? 
MosAIC_VF_GTDB %>%
  group_by(Gene) %>%
  summarise(n()) 

#### Percentage of each category in each class (based on mean numbers of genes)
MosAIC_VF_GTDB %>%
  group_by(class, category, number_file) %>%
  summarise(number_genes = n()) %>%
  group_by(class, category) %>%
  summarise(mean_genes = mean(number_genes)) %>%
  group_by(category) %>%
  mutate(count = sum(mean_genes)) %>%
  mutate(per = paste0(round(100*mean_genes / count,2),'%')) %>%
  print(n = 40)

MosAIC_VF_GTDB %>%
  group_by(family, category, number_file) %>%
  summarise(number_genes = n()) %>%
  group_by(family, category) %>%
  summarise(mean_genes = mean(number_genes)) %>%
  group_by(category) %>%
  mutate(count = sum(mean_genes)) %>%
  mutate(per = paste0(round(100*mean_genes / count,2),'%')) %>%
  print(n = 100) %>%
  ggplot() + 
  aes(x = category, y = mean_genes, fill = family) + 
  geom_bar(stat = "identity", position = "fill") + 
  theme_bw(base_size = 20) + 
  scale_fill_manual(values = c("#005C66", "#BD5129", "#9FBAE7", "#1000D2", "#0070C8", "#52CEFF", 
                               "#009464", "#E1B698", "#06BF34", "#9D49AB", "#94E4FF", "#E3F8FF", 
                               "#A2581B", "#00B9C8", "#9AC7D9", "#50E4BC", "#6E2A4D", "#D582D2", 
                               "#F0AC7A", "#F7D7EF", "#97ADFF", "#7500FF", "#FACCA3", "#52E685", 
                               "#58600B", "#E0FFF1", "#A3A3FF", "#58600B", "#FFE7CC")) + 
  coord_flip() + 
  ylab("Proportion of Genes") + 
  xlab("Virulence Factor Category")

MosAIC_VF_GTDB %>%
  group_by(family, category, number_file) %>%
  summarise(number_genes = n()) %>%
  group_by(family, category) %>%
  summarise(mean_genes = mean(number_genes)) %>%
  group_by(category) %>%
  mutate(count = sum(mean_genes)) %>%
  mutate(per = paste0(round(100*mean_genes / count,2),'%')) %>%
  print(n = 100) %>%
  ggplot() 

### Number of virulence factors in Enterobacteriaceae
MosAIC_VF_GTDB %>%
  filter(family == "Enterobacteriaceae") %>%
  group_by(genus, category) %>%
  summarise(number_genes = n()) %>%
  ggplot() + 
  aes(x = genus, y = category, size = number_genes) + 
  geom_point() +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

#### Boxplots showing distributions of VF in Each Cluster
MosAIC_VF_clean %>%
  left_join(GTDB_Full_Clean, by = c("number_file" = "user_genome")) %>%
  select(number_file, category, Gene, class) %>%
  group_by(category, class, number_file) %>%
  summarise(number_genes = n()) %>%
  ggplot() + 
  aes(x = class, y = number_genes) + 
  geom_boxplot() + 
  theme_bw(base_size = 30) + 
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) + 
  xlab("") + 
  ylab("Number of VF Genes")

ggsave(filename = "230423_VF_Per_Class_Boxplot.pdf", width = 10, height = 10)

#### Boxplots showign distribution of genes in each category
MosAIC_VF_clean %>%
  left_join(GTDB_Full_Clean, by = c("number_file" = "user_genome")) %>%
  select(number_file, category, Gene, class) %>%
  group_by(category, class, number_file) %>%
  summarise(number_genes = n()) %>%
  ggplot() + 
  aes(x = class, y = number_genes) + 
  geom_boxplot() + 
  facet_wrap(facets = "category") + 
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) + 
  xlab("Bacterial Class") + 
  ylab("Number of VF Genes")

ggsave(filename = "230423_VF_Per_Category_Boxplot.pdf")


MosAIC_VF_Winning_Cluster_Order %>%
  ggplot(aes(x = genus, y = category)) +
  geom_boxplot() + 
  coord_flip() + 
  theme_bw(base_size = 20) + 
  ylab("Number of Identified Virulence Factors (VFDB)") + 
  xlab("Genera")

#Family_VF_Distribution <- MosAIC_VF_Winning_Cluster_Order %>%
  #ggplot() + 
  #aes(x = category, y = num_associated_genes, fill = family) + 
  #geom_bar(stat = "identity", position = "fill") + 
  #theme_bw(base_size = 20) + 
  #scale_fill_manual(values = c("#005C66", "#BD5129", "#9FBAE7", "#1000D2", "#0070C8", "#52CEFF", 
   #                            "#009464", "#E1B698", "#06BF34", "#9D49AB", "#94E4FF", "#E3F8FF", 
    #                           "#A2581B", "#00B9C8", "#9AC7D9", "#50E4BC", "#6E2A4D", "#D582D2", 
     #                          "#F0AC7A", "#F7D7EF", "#97ADFF", "#7500FF", "#FACCA3", "#52E685", 
      #                         "#58600B", "#E0FFF1", "#A3A3FF", "#58600B", "#FFE7CC")) + 
  #coord_flip() + 
  #ylab("Proportion of Genes") + 
  #xlab("Virulence Factor Category")

#ggsave(Family_VF_Distribution, filename = "230423_Family_VF_Distribution.pdf", width = 12)



