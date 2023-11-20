# Load The Packages
library(tidyverse)
library(janitor)

# Read Data
checkm <- read_tsv("five_runs_combined_checkM_report2.tsv") %>%
  clean_names()

QUAST <- read_tsv("five_runs_combined_quast_report2.tsv") %>%
  filter(assembly != "Assembly")

# Replace All MMO_ with MMO- for consistency with other sample names
QUAST$assembly <- str_replace(QUAST$assembly, 'MMO_', "MMO-")

# Add columns referencing high and low quality genomes based on CheckM metrics
## Our high Quality Threshold is >= 98 % complete and <= 5% contamination
checkm_mutated <- checkm %>%
  mutate(GSC_Completeness = 
           case_when(completeness >= 95 ~ "High_Completeness", 
                     completeness < 95 ~ "Low_Completeness")) %>%
  mutate(GSC_Contamination = 
           case_when(contamination > 1 ~ "High_Contamination", 
                     contamination <= 1 ~ "Low_Contamination")) %>%
  mutate(GSC_Genome = 
           case_when(contamination <= 5 & completeness >= 98 ~ "High Quality Genomes", 
                     contamination > 5 | completeness < 98 ~ "Low Quality Genomes"))

# How many assemblies do we have in total? - 416 genomes
checkm_mutated %>%
  summarise(n())

# Function to plot number of predicted genes against genome size
bin_dot_plot <- function(clean_bin_stats){
  
  compvscontam_dot_plot <- clean_bin_stats %>%
    ggplot() + 
    aes(x = completeness, y = contamination, col = GSC_Genome) %>%
    geom_point(alpha = 0.8, size = 1) + 
    scale_x_continuous(limits = c(0, 100 )) +
    geom_vline(xintercept = 98, 
                color = "red", size = 0.5) + 
    geom_hline(yintercept = 5, 
               color = "red", size = 0.5) + 
    theme_bw(base_size = 35) + 
    scale_color_manual(values = c("#188977", "grey"), guide = "none") + 
    xlab("Completeness (%)") + 
    ylab("Contamination (%)") + 
    ylim(0, 15) + 
    xlim(90,100) + 
    labs(colour = "Completeness (%)") + 
    theme(legend.position = "right", 
          legend.direction = "vertical", 
          legend.text = element_text(size = 18))
  
  return(compvscontam_dot_plot)
  
}

# Generate Dot plot = x-axis limited to > 90% complete
dot_plot <- bin_dot_plot(checkm_mutated)
ggsave("plot/PlotsFinal/Figure_1_CheckM_Dotplot_NoLINE.pdf", dot_plot, height = 8, width = 8)
ggsave("plot/PlotsFinal/Figure_1_CheckM_Dotplot_NoLINE.png", dot_plot, height = 8, width = 8)


# Convert our checkm dataframe to long format
bin_stats_bar <- checkm_mutated %>%
  pivot_longer(c(completeness, contamination), 
               names_to = "bin_rank", 
               values_to = "value")


## Bar Chart Showing Completess and Contamination Scores
Bar_Chart <- bin_stats_bar %>%
  filter(value < 100) %>%
  ggplot() + 
  aes(reorder(x = bin_id, value), y = value, fill = bin_rank) + 
  geom_bar(stat = "identity", position = "dodge") + 
  theme_bw(base_size = 25) + 
  scale_fill_manual(values = c("#188977", "red"), guide = "none") + 
  geom_hline(yintercept = 99, 
             color = "green", size = 0.3) + 
  #geom_hline(yintercept = 50, linetype = "dotted", 
  # color = "orange", size = 1) + 
  geom_hline(yintercept = 1, 
             color = "red", size = 0.3) + 
  xlab("Genome") + 
  ylab("CheckM Contamination and Completeness Score (%)") + 
  labs(fill = "CheckM Completion vs Contamination Score") +
  theme(legend.position = "bottom", 
        legend.justification = "right", 
        legend.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "mm"),
        legend.direction = "vertical", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 15), 
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 15)) + 
  guides(fill = guide_legend(nrow = 2))

ggsave("11122_bar_chart.pdf", Bar_Chart, height = 7, width = 8)

## dataframe with genomes passing CheckM completeness and contamination 
checkm_mutated_GSC_pass <- checkm_mutated %>%
  filter(completeness >= 98 & contamination <= 5)

## 396 left over
checkm_mutated_GSC_pass %>%
  summarise(n())

# plot contig number
contig_number_plot <- ggplot(checkm_mutated, aes(x = number_contigs, fill = GSC_Genome)) + 
  geom_histogram(bins = 70) + 
  theme_bw(base_size = 25) + 
  scale_fill_manual(values = c("#188977", "grey"), guide = "none") + 
  facet_wrap(facets = "GSC_Genome") + 
  ylab("Count") + 
  xlab("Number of Contigs")

ggsave("plot/PlotsFinal/Supplementary_Figure1_Contig_Number_Distribution.pdf", contig_number_plot, height = 7, width = 8)


# Plot N50
N50 <- ggplot(checkm_mutated, aes(x = n50_contigs, fill = GSC_Genome)) + 
  geom_histogram(bins = 100) + 
  theme_bw(base_size = 25) + 
  scale_fill_manual(values = c("#188977", "grey"), guide = "none") + 
  facet_wrap(facets = "GSC_Genome") + 
  ylab("Count") + 
  xlab("N50")

ggsave("plot/PlotsFinal/Supplementary_Figure1_N50_Distribution.pdf", N50, height = 7, width = 12)

# N50 But without the Facets + low Quality
N50_NoFacet <- checkm_mutated %>%
  filter(GSC_Genome == "High Quality Genomes") %>%
  ggplot(aes(x = n50_contigs, fill = GSC_Genome)) + 
  geom_histogram(bins = 100) + 
  theme_bw(base_size = 35) + 
  scale_fill_manual(values = c("#188977", "grey"), guide = "none") + 
  #facet_wrap(facets = "GSC_Genome") + 
  ylab("Count") + 
  xlab("N50 (Bp)")

ggsave(N50_NoFacet, filename = "plot/PlotsFinal/Figure_1_N50_Distribution.pdf", height = 8, width = 8)

# plot Genome Size
genome_size_plot <- ggplot(checkm_mutated, aes(x = genome_size_bp/1e+06, fill = GSC_Genome)) + 
  geom_histogram(bins = 100) + 
  theme_bw(base_size = 25) + 
  scale_fill_manual(values = c("#188977", "grey"), guide = "none") + 
  facet_wrap(facets = "GSC_Genome") + 
  ylab("Count") + 
  xlab("Genome Size (MBp)")

ggsave(genome_size_plot, filename = "plot/PlotsFinal/Figure_1_Genome_Size_Distribution.pdf", height = 8, width = 8)

genome_size_plot_NoFacet <- checkm_mutated %>%
  filter(GSC_Genome == "High Quality Genomes") %>%
  ggplot(aes(x = genome_size_bp/1e+06, fill = GSC_Genome)) + 
  geom_histogram(bins = 100) + 
  theme_bw(base_size = 35) + 
  scale_fill_manual(values = c("#188977", "grey"), guide = "none") + 
  #facet_wrap(facets = "GSC_Genome") + 
  ylab("Count") + 
  xlab("Genome Size (Mbp)")

ggsave("plot/PlotsFinal/Figure_1_Avg_Genome_Size_Distribution.pdf", genome_size_plot_NoFacet, height = 8, width = 8)

# Plot Mean Contig Size
mean_contigs_size <- ggplot(checkm_mutated, aes(x = mean_contig_length_bp, fill = GSC_Genome)) + 
  geom_histogram(bins = 70) + 
  scale_fill_manual(values = c("#188977", "grey"), guide = "none") + 
  theme_bw(base_size = 25) + 
  facet_wrap(facets = "GSC_Genome") + 
  ylab("Count") + 
  xlab("Mean Contig Length")

ggsave("plot/PlotsFinal/Supplementary_Figure1_Mean_Contig_Size_Distribution.pdf", mean_contigs_size, height = 7, width = 12)


# Plot Genome size vs number predicted genes
NoGenes_vs_GenomeSize_LM <- checkm_mutated %>%
  filter(GSC_Genome != "NA") %>%
  ggplot() + 
  aes(x = genome_size_bp/1e+06, y = number_predicted_genes, col = GSC_Genome) + 
  geom_point(size = 1) +
  geom_smooth(method = "lm") + 
  theme_bw(base_size = 25) + 
  scale_color_manual(values = c("dark green", "grey", "dark orange", "black")) + 
  facet_wrap(facets = "GSC_Genome") + 
  ylab("Number of Predicted Genes (Prodigal)") + 
  xlab("Genome Size (Mbp)") + 
  labs(col = "Genome Quality")

ggsave("plot/PlotsFinal/Supplementary_Figure2_NoGenes_vs_GenomeSizeLM.pdf", NoGenes_vs_GenomeSize_LM, height = 7, width = 12)


# Function to fit a linear model to Genomes size vs Number of predicted genes
fit_lin_model_get_residuals <- function(dataset, threshold){
  # Fit a linear model to num predicted genes against genome size
  lm = glm(formula = number_predicted_genes~genome_size_bp, data = dataset)
  # Get residual values to the fitted line - how much does each genome deviate from the line?
  res = abs(residuals(lm))
  # Assign values based on residual - here using 400 i.e if genome of x number of contigs deviates from the 
  # model by more than 400 CDS, it doesn't fit
  res[which(res > threshold)] = "no_fit"
  res[which(res != "no_fit")] = "fit"
  test = data.frame(dataset, col = res, stringsAsFactors = F) 
  
  return(test)
}

# Fit linear model with residual coloured in red based on threshold
checkm_mutated_90_lm <- as_tibble(fit_lin_model_get_residuals(checkm_mutated_90, threshold = 0))

# Plot Linear model with residuals over threshold in red (threshold specified in previous step)
lm_high_quality_genomes <- ggplot(checkm_mutated_90_lm, aes(x = genome_size_bp/1e+06, y = number_predicted_genes, color = col)) + 
  geom_point(size = 2, alpha = 0.8) + 
  geom_smooth(se = F, method = "glm") + 
  theme_bw(base_size = 25) +
  scale_color_manual(values = c("black", "#ff6666"), guide = "none") + 
  ylab("Number of Predicted Genes (Prodigal)") +
  xlab("Genome Size (MBP)") + 
  theme(legend.text = element_blank()) 

ggsave("plot/PlotsFinal/Supplementary_Figure2_lm_High_Quality_Genomes_Scatterplot.pdf", lm_high_quality_genomes, height = 7, width = 8)
ggsave("plot/PlotsFinal/Supplementary_Figure2_lm_High_Quality_Genomes_Scatterplot.png", lm_high_quality_genomes, height = 7, width = 8)

# filter! 
## join checkm with QUAST
joined_genome_table <- checkm_mutated_GSC_pass %>%
  left_join(y = QUAST, by = c("bin_id" = "assembly"))

joined_genome_table$avg_coverage <- as.numeric(joined_genome_table$avg_coverage)

## plot AVG Read Coverage
joined_genome_table %>%
  ggplot(aes(x = avg_coverage, fill = GSC_Genome)) + 
  geom_histogram(bins = 150) + 
  theme_bw(base_size = 35) + 
  scale_fill_manual(values = c("#188977", "grey"), guide = "none") + 
  geom_vline(xintercept = 10, colour = "red") + 
  ylab("Count") + 
  xlab("Average Read Coverage")

ggsave("plot/PlotsFinal/Figure_1_Avg_Read_Coverage.pdf", height = 8, width = 8)
ggsave("plot/PlotsFinal/Figure_1_Avg_Read_Coverage.png", height = 8, width = 8)

## Same number as before? - Sanity Check 
joined_genome_table %>%
  summarise(n()) # 396 - phew

## Our final table of high-quality isolates!
filtered_table <- joined_genome_table %>%
  mutate_at("avg_coverage", as.numeric) %>%
  #filter(number_contigs < 400) %>%
  filter(avg_coverage > 10)
  #filter(col == "fit") 

## How many final? 
filtered_table %>%
  summarise(n()) # 394

write.table(x = filtered_table, file = "plot/PlotsFinal/Supplementary_Table1_MosAIC_Filtered_Genome_Table.tsv", sep = "\t", row.names = FALSE)

## Get just the names so you can filter with bash - use these sample names to filter the genomes for the remainder of the analysis
filtered_table_names <- filtered_table %>%
  select(bin_id)

write.table(x = filtered_table_names, file = "plot/PlotsFinal/Supplementary_Table2_MosAIC_Filtered_Genome_Table_Names.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


