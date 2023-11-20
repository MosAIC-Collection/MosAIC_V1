# Load Packages
library(tidyverse)
library("wesanderson")
library(RColorBrewer)
library(ggsci)
library(ggpubr)
library(MetBrewer)
library("cowplot")
library(janitor)


# Replace All MMO_ with MMO- for correct consistency
quast <- read_tsv("five_runs_combined_quast_report2.tsv") %>%
  filter(assembly != "Assembly")

# Replace All MMO_ with MMO- for correct consistency
quast_clean <- quast %>%
  separate(assembly, into = c("sample", "junk"), sep = "_contigs") %>%
  select(!junk) %>%
  mutate_at(2:29, as.numeric) 


# Replace All MMO_ with MMO- for correct consistency
quast_clean$sample <- str_replace(quast_clean$sample, 'MMO_', "MMO-")

# Histogram of Percentage of Reads Mapped and Properly Paired
perc_paired <- ggplot(quast_clean, aes(x = perc_properly_paired)) + 
  geom_histogram(bins = 150) + 
  theme_bw(base_size = 25) + 
  ylab("Count") + 
  xlab("Reads Mapped and Properly Paired (%)")

ggsave("plot/PlotsFinal/Supplementary_Figure3_Percentage_Reads_Mapped_And_Paired.pdf", perc_paired, height = 7, width = 8)

# Histogram of Average Read Coverage
avg_coverage <- ggplot(quast_clean, aes(x = avg_coverage)) + 
  geom_histogram(bins = 150) + 
  theme_bw(base_size = 35) + 
  geom_vline(xintercept = 10, colour = "red") + 
  ylab("Count") + 
  xlab("Average Coverage")

ggsave("plot/PlotsFinal/Figure_1_Avg_Read_Coverage.pdf", avg_coverage, height = 8, width = 8)

# Histogram of Reads Mapping to Contigs
perc_mapped <- ggplot(quast_clean, aes(x = perc_mapped)) + 
  geom_histogram(bins = 150) + 
  theme_bw(base_size = 25) + 
  ylab("Count") + 
  xlab("Reads Mapped (%)")

ggsave("plot/PlotsFinal/Supplementary_Figure3_Percentage_Reads_Mapped.pdf", perc_mapped, height = 7, width = 8)

# Histogram of Total Reads used in Assemblies
total_reads <- ggplot(quast_clean, aes(x = total_reads)) + 
  geom_histogram(bins = 150) + 
  theme_bw(base_size = 25) + 
  ylab("Count") + 
  xlab("Read Count")

ggsave("plot/PlotsFinal/Supplementary_Figure_3_Total_Reads_Distribution.pdf", total_reads, height = 7, width = 8)

# Scatter plot showing correlation between Read coverage and total reads used for assembly
total_reads_vs_avg_coverage <- ggplot(quast_clean, aes(x = total_reads, y = avg_coverage)) + 
  geom_point() + 
  geom_smooth(method = "lm", col = "black") + 
  stat_cor(label.x = 1, label.y = 1300, size = 6) + 
  theme_bw(base_size = 25) + 
  ylab("Average Assembly Coverage") + 
  xlab("Total Reads Used for Assembly") + 
  theme(legend.key = element_blank())

ggsave("plot/PlotsFinal/Supplementary_Figure_3_Total_Reads_vs_Avg_Coverage_ScatterPlot.pdf", total_reads_vs_avg_coverage, height = 7, width = 8)


quast_clean_long_contig <- quast_clean %>%
  select(sample, contigs_0bp, contigs_1000bp, 
         contigs_5000bp, contigs_10000bp, 
         contigs_25000bp, contigs_50000bp, total_reads) %>%
  pivot_longer(cols = !c(sample, total_reads), names_to = "category", values_to = "length") 

# Faceted Histograms of Contig Size Range 
contig_labels <- c(
  contigs_0bp = "Contigs > 0bp",
  contigs_1000bp = "Contigs > 1000bp",
  contigs_5000bp = "Contigs > 5000bp",
  contigs_10000bp = "Contigs > 10000bp",
  contigs_25000bp = "Contigs > 25000bp",
  contigs_50000bp = "Contigs > 50000bp"
)

contig_distribution <- quast_clean_long_contig %>%
  ggplot() + 
  aes(x = length) + 
  geom_histogram(binwidth = 1, bins = 70) + 
  theme_bw(base_size = 12) + 
  xlim(0, 500) + 
  facet_wrap(category ~ ., labeller = as_labeller(contig_labels), strip.position = "top", ncol = 1) + 
  xlab("Number of Contigs") + 
  ylab("Count") + 
  theme(axis.title.x = element_text(size = 20), 
        axis.text.x = element_text(size = 18), 
        axis.title.y = element_text(size = 20), 
        axis.text.y = element_text(size = 18))

ggsave("111112_contig_distribution.pdf", contig_distribution, height = 7, width = 8)

#### Summary Statistics 
quast_clean %>% 
  select(sample, total_reads) %>%
  filter(total_reads != 0) %>%
  summarise(mean_reads = mean(total_reads / 1000000), 
            range_reads = range(total_reads / 1000000))
