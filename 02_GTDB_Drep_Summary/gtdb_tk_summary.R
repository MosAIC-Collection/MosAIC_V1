##
# Summarise GTDB Classications from non dereplicated dataset 
##
#### Load Packages
library(tidyverse)
library(MetBrewer)

#### Read Data 
gtdb_nodrep <- read_tsv("MosAIC_V1/02_GTDB_Drep_Summary/190123_gtdb.bac120.summary_nodrep.tsv")

#### Clean
gtdb_class_clean <- gtdb_nodrep %>%
  separate(classification, into = c("domain", "parent", "class", "order", "family", "genus", "species"), sep = ";") %>% 
  separate(species, into = c("junk", "species"), sep = "s__") %>%
  separate(family, into = c("junk2", "family"), sep = "f__") %>%
  separate(genus, into = c("junk3", "genus"), sep = "g__") %>%
  separate(domain, into = c("junk4", "domain"), sep = "d__") %>%
  separate(parent, into = c("junk5", "parent"), sep = "p__") %>%
  separate(class, into = c("junk6", "class"), sep = "c__") %>%
  separate(order, into = c("junk7", "order"), sep = "o__") %>%
  select(!c(junk, junk2, junk3, junk4, junk5, junk6, junk7))

# same plot produced in gtdb_tk_summary.R
bar_gtdb_class <- gtdb_class_clean %>% 
  ggplot() + 
  aes(x = family, group = family) + 
  geom_histogram(bins = 70, stat = "count") + 
  theme_bw(base_size = 20) + 
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1), 
        legend.title = element_blank(), 
        legend.position = "none") + 
  ylab("Number of Isolates") + 
  xlab("GTDB-Tk Family Classification") + 
  facet_wrap(~class, scales = "free_x", ncol = 4)
  
ggsave(filename = "plots/030122_bar_gtdb_class.pdf", plot = bar_gtdb_class, width = 20, height = 8)

# How many genera, facetted by family
genera_class_bar <- gtdb_class_clean %>%
  ggplot() + 
  aes(x = genus) + 
  geom_histogram(stat = "count") + 
  theme_bw(base_size = 8) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1), 
        legend.title = element_blank(), 
        legend.position = "none") + 
  ylab("Number of Dereplicated Genera") + 
  xlab("Genus Classification") + 
  facet_wrap(~family, scales = "free_x", ncol = 4)


ggsave(filename = "plots/220223_22_genera_class_bar.pdf", plot = genera_class_bar, width = 8, height = 12)

## Top Genera
gtdb_class_clean %>%
  group_by(genus) %>%
  summarise(number_species = n()) %>%
  filter(number_species > 4) %>%
  ggplot() + 
  aes(x = reorder(genus, number_species), y = number_species) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() + 
  theme_bw(base_size = 20) + 
  ylab("Number of Genera") + 
  xlab("Genus")

ggsave(filename = "plots/220223_TopGeneraRepresentatives.pdf", width = 10, height = 10)

## Top Species - produces a different set of bacteria due to dRep clustering species >95% 
### Example is Serratia - serratia nevei is >95% similar to Serratia marscecens, therefore clusters with dRep
gtdb_class_clean %>%
  group_by(species) %>%
  summarise(number_species = n()) %>%
  filter(number_species > 4) %>%
  filter(species != "") %>%
  ggplot() + 
  aes(x = reorder(species, number_species), y = number_species) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() + 
  theme_bw(base_size = 20) + 
  ylab("Number of Species") + 
  xlab("Species")

ggsave(filename = "plots/220223_TopSpeciesRepresentatives.pdf", width = 10, height = 10)

#gtdb_class_clean %>%
  #filter(parent == "Proteobacteria") %>%
  #ggplot() + 
  #aes(x = user_genome, y = genus) + 
  #geom_point() + 
  #theme_bw(base_size = 15) + 
  #theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1)) + 
  #ylab("Genera") + 
  #xlab("Sample") + 
  #ggtitle(label = "Cultured Proteobacteria from the Mosquito Microbiome")
  
#ggsave(filename = "301122_Enterobacteriaceae_tile.pdf", plot = Enterobacteriaceae_tile, width = 14, height = 10)

#Enterobacteriaceae_liverpool <- gtdb_class_clean %>%
  #filter(family == "Enterobacteriaceae") %>%
  #filter(!grepl("AS", user_genome)) %>%
  #filter(!grepl("HN", user_genome)) %>%
  #filter(!grepl("US", user_genome)) %>%
  #filter(!grepl("M50", user_genome))

#Enterobacteriaceae_groups <- Enterobacteriaceae_liverpool %>%
  #group_by(species) %>%
  #summarise(number = n())
  #ggplot() + 
  #aes(x = species, y = number) + 
  #geom_col(width = 0.5) + 
  #theme_bw(base_size = 20) + 
  #theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) + 
  #ylab("Species") + 
  #xlab("Sample")

#ggsave(filename = "301122_Enterobacteriaceae_bar_liverpool.pdf", plot = Enterobacteriaceae_groups, width = 10, height = 8)

