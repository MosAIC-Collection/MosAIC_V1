##
# Summarise data from dreplicated data set 
##
#### Load Packages
library(tidyverse)
library(MetBrewer)
library(ggtreeExtra)
library(grid)

#### Read Data 
gtdb_class <- read_tsv("MosAIC_V1/02_GTDB_Drep_Summary/291122_gtdbtk_isolates_drep.bac120.summary.tsv")
drep <- read_csv("MosAIC_V1/02_GTDB_Drep_Summary/291122_Widb.csv")

#### Clean GTDB Classifications
gtdb_class_clean <- gtdb_class %>%
  separate(classification, into = c("domain", "parent", "class", "order", "family", "genus", "species"), sep = ";") %>% 
  separate(species, into = c("junk", "species"), sep = "s__") %>%
  separate(family, into = c("junk2", "family"), sep = "f__") %>%
  separate(genus, into = c("junk3", "genus"), sep = "g__") %>%
  separate(domain, into = c("junk4", "domain"), sep = "d__") %>%
  separate(parent, into = c("junk5", "parent"), sep = "p__") %>%
  separate(class, into = c("junk6", "class"), sep = "c__") %>%
  separate(order, into = c("junk7", "order"), sep = "o__") %>%
  select(!c(junk, junk2, junk3, junk4, junk5, junk6, junk7))

#### Clean Drep 
drep_clean <- drep %>%
  separate(genome, into = c("user_genome", "junk"), sep = ".fa") %>%
  select(!junk)
  
#### Get dRep Names
drep_clean_names <- drep_clean %>%
  select(user_genome)

#### Use this for dRep Names
write.table(x = drep_clean_names, file = "plots/PlotsFinal/Supplementary_Table_3_drep_Species_Representative_Sample_Names.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#### Join 
gtdb_drep_join <- gtdb_class_clean %>%
  left_join(drep_clean)

#### Check join is okay
is.na(gtdb_drep_join$user_genome) ### all FALSE - good

#### How many isolates in each species 
gtdb_drep_join$cluster_members

#### Function to get Facets of family, with each species present + number of isolates 
plot_bacterial_family_facets <- function(data, target){
  facet <- data %>%
    filter(class == target) %>% 
    ggplot() + 
    aes(x = reorder(genus, cluster_members), y = cluster_members) + 
    geom_col() + 
    theme_bw(base_size = 20) + 
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) + 
    facet_wrap(~family, scales = "free_x") + 
    ylab("Number of Strains") + 
    xlab("GTDB-Tk Species Classification")
  
  return(facet)
}

actinomycetia <- plot_bacterial_family_facets(gtdb_drep_join, "Actinomycetia")
alphaproteobacteria <- plot_bacterial_family_facets(gtdb_drep_join, "Alphaproteobacteria")
bacteroidia <- plot_bacterial_family_facets(gtdb_drep_join, "Bacteroidia")
bacilli <- plot_bacterial_family_facets(gtdb_drep_join, "Bacilli")
gammaproteobacteria <- plot_bacterial_family_facets(gtdb_drep_join, "Gammaproteobacteria")

ggsave(filename = "plots/PlotsFinal/Supplementary_Figure_4_Actinomycetia_Species_Representatives.pdf", plot = actinomycetia, width = 15, height = 10)
ggsave(filename = "plots/PlotsFinal/Supplementary_Figure_5_Alphaproteobacteria_Species_Representatives.pdf", plot = alphaproteobacteria, width = 15, height = 10)
ggsave(filename = "plots/PlotsFinal/Supplementary_Figure_6_Bacteroidia_Species_Representatives.pdf", plot = bacteroidia, width = 15, height = 10)
ggsave(filename = "plots/PlotsFinal/Supplementary_Figure_7_Bacilli_Species_Representatives.pdf", plot = bacilli, width = 15, height = 10)
ggsave(filename = "plots/PlotsFinal/Supplementary_Figure_8_Gammaproteobacteria_Species_Representatives.pdf", plot = gammaproteobacteria, width = 15, height = 10)

gtdb_drep_join %>%
  ggplot() + 
  aes(x = genus, y = cluster_members) + 
  geom_bar(stat = "identity") + 
  facet_grid(rows = vars(family), cols = vars(class), scales = "free", space = "free_x") + 
  coord_flip() + 
  theme_bw() + 
  theme(strip.text.y.right = element_text(angle = 0)) 
  

ggsave(filename = "plots/PlotsFinal/Supplementary_Figure3_MosAIC_Collection_FacetGrid.pdf", width = 20, height = 50, limitsize = F)

## Basic Plots
### What are the top species? 
top_species <- gtdb_drep_join %>%
  filter(cluster_members > 5) %>%
  ggplot() + 
  aes(reorder(x = species, cluster_members), y = cluster_members, fill = family) + 
  geom_bar(stat = "identity", position = "dodge") + 
  theme_bw(base_size = 35) + 
  ylab("Number of Isolates") + 
  xlab("Species") + 
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) + 
  labs(fill = "Family") + 
  scale_fill_manual(values = c("#BD5129",
                               "#0070C8",
                               "#A2581B",
                               "#D582D2",
                               "#F0AC7A",
                               "#FFF4B1",
                               "#FFAAB7", "red"), ) +
  theme(legend.position = "none", 
        plot.margin = unit(c(1, 1, 1, 5), "lines")) + 
  guides(fill = guide_legend(ncol = 2)) + 
  coord_flip()

ggsave(plot = top_species, filename = "plots/PlotsFinal/Figure_2_Abundant_Species_Representatives.pdf", width = 15, height = 15)

#### Distribution of Families subdivided by Class - 
Class_Family_Plot <- gtdb_drep_join %>%
  ggplot() + 
  aes(x = reorder(family, cluster_members), y = cluster_members) + 
  geom_col() + 
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1)) + 
  facet_wrap(~class, scales = "free_x") + 
  ylab("Number of Strains") + 
  xlab("GTDB-Tk Family Classification")

ggsave(plot = Class_Family_Plot, filename = "plots/PlotsFinal/Supplementary_Figure_9_Family_Representatives.pdf", width = 11, height = 8)

# A Whole List of the Isolates Collection at the species level 
# Blanks are those without a species classification
EverythingEverywhereAllAtOnce <- gtdb_drep_join %>%
  #filter(species != '') %>%
  ggplot() + 
  aes(reorder(x = species, cluster_members), y = cluster_members) + 
  geom_bar(stat = "identity", position = "dodge") + 
  theme_bw(base_size = 10) + 
  ylab("Number of Isolates") + 
  xlab("Species") + 
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 15)) + 
  labs(fill = "Family") + 
  scale_fill_manual(values = met.brewer("Redon")) + 
  theme(legend.position = "bottom") + 
  guides(fill = guide_legend(ncol = 2)) + 
  coord_flip()
  
ggsave(plot = EverythingEverywhereAllAtOnce, filename = "plots/PlotsFinal/Supplementary_Figure_10_MosAIC_Collection_Species_List.pdf", width = 8, height = 15)


#### Statistics 
### How many single species representatives? 
gtdb_drep_join %>%
  summarise(n()) # 142 

### How many representatives in each class? 
gtdb_drep_join %>%
  group_by(class) %>%
  summarise(n())

### How many different bacterial families
gtdb_drep_join %>%
  group_by(family) %>%
  summarise(number_SSR = n()) # 29 
  
### Top representatives in each class? 
gtdb_drep_join %>%
  group_by(class, family) %>%
  summarise(number_SSR = n()) %>%
  print(n = 400)

### Top Isolate Representative in each family (those big bars)
gtdb_drep_join %>%
  group_by(class, family) %>%
  summarise(number_isolates = sum(cluster_members)) %>%
  print(n = 400)

### What are these in percentages? 
gtdb_drep_join %>%
  group_by(class, family) %>%
  summarise(number_isolates = sum(cluster_members)) %>%
  mutate(number_isolates_perc = (number_isolates/392)*100) %>%
  print(n = 400) # enterobactericae are 40% of the collection

### Top species representatives 
gtdb_drep_join %>%
  select(family, species, cluster_members) %>%
  arrange(desc(cluster_members)) %>%
  print(n = 400)

### Genus classifications - potentially novel species
gtdb_drep_join %>%
  select(genus, species, cluster_members) %>%
  filter(species == "") %>%
  ggplot() + 
  aes(x = reorder(genus, cluster_members), y = cluster_members) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 20) + 
  coord_flip() + 
  xlab("Genus") + 
  ylab("Number of Isolates")
  
ggsave(filename = "plots/PlotsFinal/Supplementary_Figure_3_Genus_Classifications.pdf", width = 10, height = 12)
ggsave(filename = "plots/PlotsFinal/Supplementary_Figure_3_Genus_Classifications.png", width = 10, height = 12)

  
