library(vegan)
library(micropan)

#### Gene Accumulation Curves 
# Load Dataframe
df = data.frame(cluster = character(0), richness = numeric(0),genomes = numeric(0),sd = numeric(0))
# Files in WD
files = list.files(path = "MosAIC_V1/08_GeneAccumulationCurve/", full.names = T, pattern = "*Rtab")
# Loop and split the files by "_" 
for (f in files){
  cluster = strsplit(basename(f), split = "_", fixed = T)[[1]][1]
  print(cluster)
  mydata <- data.frame(t(read.table(f, header = T, row.names = 1, comment.char = "", quote = )))
  sp <- specaccum(mydata, "random", permutations=100)
  df = rbind(df, data.frame(cluster = rep(cluster, length(sp$sites)),
                            richness = sp$richness,
                            genomes = sp$sites,
                            sd = sp$sd))
}

write.table(df, "MosAIC_V1/08_GeneAccumulationCurve/260423_Accumilation_Curve.csv", col.names = T, row.names = F, quote = F, sep = ',')

# Read new dataframe 
df = read.table("MosAIC_V1/08_GeneAccumulationCurve/260423_Accumilation_Curve.csv", header = T, comment.char = "", sep = ",", stringsAsFactors = F)
df = cbind(df, min= df$richness-df$sd, max = df$richness+df$sd )
#df$cluster = factor(df$cluster,df$Cluster)
df = df[which(df$genomes %% 10 == 0),] ## to visualise more clearly

# Plot
ggplot(df, aes(x = genomes, y = richness, col = cluster)) + geom_line(alpha = 1, linewidth = 0.4) + 
  geom_point(size = 2, alpha = 0.8) +
  theme_bw(base_size = 25) + xlab("Number of Genomes") +
  ylab("Predicted Number of Genes") +
  geom_errorbar(aes(ymin=min, ymax=max), width=.2,alpha = 0.5,
                position=position_dodge(0.05)) +
  #scale_color_manual(values = graphics$Color, guide = F) + scale_shape_manual(values = graphics$Shape, guide = F) + ggtitle("C") +
  scale_x_continuous(limits = c(0,400)) + 
  scale_color_manual(values = c("#AC6700", "#229EB2", "#D1CB67"), labels = c("Enterobacter asburiae", "Elizabethkingia anophelis", "Serratia marscecens")) + 
  guides(col = guide_legend(ncol = 1, title = "Target Species", override.aes = list(size = 4)))


ggsave("gene_accumulation_curves/MosAIC_Gene_Accumulation_Curves.pdf", height = 10, width = 15)
