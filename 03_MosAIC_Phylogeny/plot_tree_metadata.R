##
# Make MosAIC's phylogeny with Drep and GTDB classifications
##
### Load the Packages
library(ggstar)
library(ggtreeExtra)
library(ggtree)
library(ggnewscale)
library(ggstar)
library(tidyverse)
library(treeio)
library(MetBrewer)

### Load the tree - phylogeny based on the 16s rRNA region of Mosquito Associated Bacteria
tree_mos <- read.tree("MosAIC_V1/03_MosAIC_Phylogeny/270423_concatenated_16s_fileappend_final_with_outgroup_aln_trimal.fa.treefile")
tree_mos_root <- tidytree::root(phy = tree_mos, outgroup = "16S_rRNA")

### Load the meta data of names and sample IDs corresponding to the nodes of the tree 
dat1_mos <- read.csv("MosAIC_V1/03_MosAIC_Phylogeny/MMO142ADJ_240223_tree_tibble_join_td_select.csv")
dat1_mos <- dat1_mos[1:142,]

### Load the metadata of the full table 
dat2_mos <- read.csv("MosAIC_V1/03_MosAIC_Phylogeny/MMO142ADJ_240223_071222_tree_tibble_drep_join_full_table.csv")
dat2_mos <- dat2_mos[1:142,]


### Make the first tree - fan layout and include bootstrap values
ggtree(tree_mos_root, layout="fan", size=0.2, open.angle=20) 
#geom_point2(aes(subset=!isTip & UFboot > 80, x = branch, col = UFboot), size = 0.7, na.rm = T) + 
#scale_color_gradientn(colors = met.brewer("Pillement", n = 3))
#geom_hilight(data=nodedf, mapping=aes(node=node),
#extendto=8, alpha=0.5, fill="grey", color="grey50",
#size=0.01) 
#geom_cladelab(data=labdf, 
#mapping=aes(node=node, 
#label=label,
#offset.text=pos),
#hjust=0.5,
#angle="auto",
#barsize=NA,
#horizontal=FALSE, 
#fontsize=1.4,
#fontface="italic"
#)


p1 <- p + new_scale_color() 


### Add first layer which are points at the end of each node corresponding to the order of bacteria
p1 <- p1 %<+% dat1_mos + geom_point2(
  mapping=aes(col=class),
  position="identity") +
  scale_colour_manual(values=c("#FFC125","#87CEFA","#7B68EE","#808080",
                               "#800080", "#9ACD32","#D15FEE","#FFC0CB",
                               "#EE6A50","#8DEEEE", "#006400","#800000",
                               "#B0171F","#191970"),
                      guide=guide_legend(keywidth = 0.5, 
                                         keyheight = 0.5, order=1,
                                         override.aes=list(starshape=15)),
                      na.translate=FALSE) +
  scale_starshape_manual(values=c(15, 1),
                         guide=guide_legend(keywidth = 0.5, 
                                            keyheight = 0.5, order=2),
                         na.translate=FALSE) +
  scale_size_continuous(range = c(1, 2.5),
                        guide = guide_legend(keywidth = 0.5, 
                                             keyheight = 0.5, order=3,
                                             override.aes=list(starshape=15)))




### Function to extract specific columsna and convert in specific way for the rest of ggtree gheatmap function
extract_column_and_convert <- function(data, column){
  target <- data[, c(1, column)]
  target2 <- target %>%
    remove_rownames %>%
    column_to_rownames(var = "label")
  
  return(target2)
}

order_data <- extract_column_and_convert(dat2_mos, 5)
source_data <- extract_column_and_convert(dat2_mos, 16)
mosquito_data <- extract_column_and_convert(dat2_mos, 17)
lab_field <- extract_column_and_convert(dat2_mos, 18)
family_data <- extract_column_and_convert(dat2_mos, 6)

### Cleaning metadata - Replace Lab with lab in lab_field df 
lab_field$lab_field <- gsub("Lab", "lab", lab_field$lab_field)

### Second layer showing which species the bacteria were isolated from
p3 <- gheatmap(p1, mosquito_data, width = 0.05, offset = 0.02, colnames = FALSE, legend_title = title) + 
  scale_fill_manual(values = c("#1D9A6C", "#43BC62", "#7AD76F", "#C5EBA6", 
                               "#C080B2", "#CB96CC", "#CCADD8", "#7F85C1", "#96A8CC", "#ADC5D8", "#1852A8", "#040300"))

p4 <- p3 + new_scale_fill() 


### Third layer showing where the bacteria were isolated from 
p5 <- gheatmap(p4, source_data, width = 0.05, offset = 0.06, colnames_position = "top", colnames = FALSE, legend_title = c("Isolation Source")) +  
  scale_fill_manual(values = c("#A88100", "#5B7888", "#FF4D38", "#B06B7B", "#C998BB", "#5A86C6"))


### Potentially going to use this to highlight specific clade labels
p4 + geom_cladelab(node = 31, label = "fsdf")

p6 <- p5 + new_scale_fill()


###  Fourth layer to show where (if lab or field) the bacteria were isolated from
p7 <- gheatmap(p6, lab_field, width = 0.05, offset = 0.1, colnames = FALSE) + 
  scale_fill_manual(values = c("#67E05D", "#DCEA31"))


### Fifth Layer to get number of isolates within each cluster
p8 <- p7  + new_scale_fill() +
  #geom_fruit(data=dat2_mos, geom=geom_tile,
  #mapping=aes(y=label, x=domain, fill=class),
  #color = "grey50", size = 0.5) +
  #scale_alpha_continuous(range=c(0, 1),
  #guide=guide_legend(keywidth = 0.3, 
  #keyheight = 0.3, order=5)) +
  geom_fruit(
    data=dat2_mos, 
    geom=geom_bar, 
    mapping=aes(
      y=label, x=cluster_members
    ),
    pwidth=0.5,
    offset = 0.25, 
    orientation="y", 
    stat="identity",
    axis.params = list(
      axis = "x", 
      text.size = 1.8,
      hjust = 1, 
      vjust = 0.5, 
    ),
    grid.params = list()
  ) +
  #scale_fill_manual(values=c("#0000FF","#FFA500","#FF0000",
  #"#800000", "#006400","#800080","#696969"),
  #guide=guide_legend(keywidth = 0.3, 
  #keyheight = 0.3, order=4))
  geom_treescale(fontsize=2, linesize=0.1, x=1, y=0.05) +
  theme(legend.position=c(1.05, 0.5),
        legend.background=element_rect(fill=NA),
        legend.title=element_text(size=12),
        legend.text=element_text(size=10),
        legend.spacing.y = unit(0.2, "cm"),
  ) #+ 
#geom_fruit(data=dat2_mos, geom=geom_tile, 
#mapping = aes(y = label, x = source, fill = source), 
#offset = 0.0001)



p8

p8 + geom_highlight(node = c(21:86), to.bottom = TRUE, extendto = 0.74, alpha = 0.2)
p8 + geom_highlight(node = c(9:20), to.bottom = TRUE, extendto = 0.74, alpha = 0.2, size = 0.001)

p8 + geom_highlight(node = 39, to.bottom = TRUE, extendto = 0.74, alpha = 0.2, size = 0.001)

ggsave(filename = "plots/PlotsFinal/Supplementary_Figure_11_Extended_Metadata_MosAIC_Phylogeny_Circular.pdf", plot = p8)


#### Create another tree showing class classifications, bootstrap, family, and cluster number 
#### More representative of the dereplicated genomes

pv2 <- ggtree(tree_mos_root, layout="fan", size=0.2, open.angle=20) 
#geom_point2(aes(subset=!isTip & UFboot > 80, x = branch, col = UFboot), size = 0.7, na.rm = T) + 
#scale_color_gradientn(colors = met.brewer("Pillement", n = 3))
#geom_point2(aes(subset=!isTip & UFboot > 80, x = branch, col = UFboot), size = 0.7, na.rm = T) + 
#scale_color_gradientn(colors = met.brewer("Pillement", n = 3))

pv2 <- pv2 + new_scale_color() 


pv2_1 <- pv2 %<+% dat1_mos + geom_point2(
  mapping=aes(col=class),
  position="identity") +
  scale_colour_manual(values=c("#B999CC","#1D9A6C","#002ACD","#FFC60B",
                               "#B64521", "#9ACD32","#D15FEE","#FFC0CB",
                               "#EE6A50","#8DEEEE", "#006400","#800000",
                               "#B0171F","#191970"),
                      guide=guide_legend(keywidth = 0.5, 
                                         keyheight = 0.5, order=1,
                                         override.aes=list(starshape=15)),
                      na.translate=FALSE) +
  scale_starshape_manual(values=c(15, 1),
                         guide=guide_legend(keywidth = 0.5, 
                                            keyheight = 0.5, order=2),
                         na.translate=FALSE) +
  scale_size_continuous(range = c(1, 2.5),
                        guide = guide_legend(keywidth = 0.5, 
                                             keyheight = 0.5, order=3,
                                             override.aes=list(starshape=15)))

pv2_2 <- gheatmap(pv2_1, family_data, width = 0.05, offset = 0.0, colnames_position = "top", colnames = FALSE, legend_title = "Family") + 
  scale_fill_manual(values = c("#005C66", "#BD5129", "#9FBAE7", "#1000D2", "#0070C8", "#52CEFF", 
                               "#009464", "#E1B698", "#06BF34", "#9D49AB", "#94E4FF", "#E3F8FF", 
                               "#A2581B", "#00B9C8", "#9AC7D9", "#50E4BC", "#6E2A4D", "#D582D2", 
                               "#F0AC7A", "#F7D7EF", "#176E93", "#7500FF", "#FFF4B1", "#52E685", 
                               "#FF9146", "#E0FFF1", "#A3A3FF", "#FFAAB7", "#FFE7CC")) + 
  guides(fill = guide_legend(ncol = 1, title = "Classification", label.position = "bottom")) 

ggsave(pv2_2, filename = "plots/PlotsFinal/Figure2_ALT_LegendKey.pdf", width = 14, heigh = 35)

pv2_3 <- pv2_2  + new_scale_fill() +
  #geom_fruit(data=dat2_mos, geom=geom_tile,
  #mapping=aes(y=label, x=domain, fill=class),
  #color = "grey50", size = 0.5) +
  #scale_alpha_continuous(range=c(0, 1),
  #guide=guide_legend(keywidth = 0.3, 
  #keyheight = 0.3, order=5)) +
  geom_fruit(
    data=dat2_mos, 
    geom=geom_bar, 
    mapping=aes(
      y=label, x=cluster_members, fill = family
    ),
    pwidth=0.5,
    offset = 0.1, 
    orientation="y", 
    stat="identity",
    axis.params = list(
      axis = "x", 
      text.size = 3,
      hjust = 0.5, 
      vjust = 1.5, 
    ),
    grid.params = list(size = 0.5)
  ) +
  scale_fill_manual(values = c("#005C66", "#BD5129", "#9FBAE7", "#1000D2", "#0070C8", "#52CEFF", 
                               "#009464", "#E1B698", "#06BF34", "#9D49AB", "#94E4FF", "#E3F8FF", 
                               "#A2581B", "#00B9C8", "#9AC7D9", "#50E4BC", "#6E2A4D", "#D582D2", 
                               "#F0AC7A", "#F7D7EF", "#97ADFF", "#7500FF", "#FFF4B1", "#52E685", 
                               "#FF9146", "#E0FFF1", "#A3A3FF", "#FFAAB7", "#FFE7CC")) + 
  #scale_fill_manual(values=c("#0000FF","#FFA500","#FF0000",
  #"#800000", "#006400","#800080","#696969"),
  #guide=guide_legend(keywidth = 0.3, 
  #keyheight = 0.3, order=4))
  geom_treescale(fontsize=5, linesize=0.1, x=0.5, y=-2) +
  theme(legend.position=c(1.3, 0.5),
        legend.background=element_rect(fill=NA),
        legend.title=element_blank(),
        legend.text=element_text(size=18),
        legend.spacing.y = unit(0.2, "cm"),
  ) #+ 
#geom_fruit(data=dat2_mos, geom=geom_tile, 
#mapping = aes(y = label, x = source, fill = source), 
#offset = 0.0001)

ggsave(filename = "plots/PlotsFinal/Figure_2_MosAIC_dRep_Cluster_Phylogeny.pdf", plot = pv2_3, height = 10, width = 22)



