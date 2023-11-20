##
# Wrangle Abricate Results
##

## Load Packages
library(tidyverse)
library(janitor)

MosAIC_VF <- read_tsv("210423_MosAIC_VF_Search.tsv") %>%
  clean_names()

VFDB_Higher_Categories <- c("Adherence", "Invasion", "Effector delivery system", 
                            "Motility", "Exotoxin", "Exoenzyme", "Immune modulation", 
                            "Biofilm", "Nutritional/Metabolic factor", "Stress survival", 
                            "Post-translational modification", "Antimicrobial activity/Competitive advantage", 
                            "Regulation", "Others")

# Separate the PRODUCT column into three columns
MosAIC_VF_clean <- MosAIC_VF %>%
  mutate(category = str_extract_all(product, paste(VFDB_Higher_Categories, collapse = "|"))) %>%
  mutate(product = str_replace_all(product, "VFG\\S*?\\(", "(")) %>%
  mutate(product = gsub("\\(gb.*?\\)", "", product), 
         product = str_remove(product, "VFG048552")) %>%
  mutate(Gene = stringr::str_extract(product, "(?<=\\()[^)]+(?=\\))")) %>%
  mutate(product = stringr::str_replace(product, "\\([^)]+\\)", "")) %>%
  select(number_file, percent_coverage, percent_identity, category, Gene, product) %>%
  mutate(category = sapply(category, toString))
  
MosAIC_VF_clean$product <- sub("^\\((.*?)\\).*", "\\1", MosAIC_VF_clean$product)
MosAIC_VF_clean$Gene <- gsub("[()]", "", MosAIC_VF_clean$Gene)
MosAIC_VF_clean$product <- gsub("^\\(.*?\\)\\s*", "", MosAIC_VF_clean$product)

write.table(x = MosAIC_VF_clean, file = "210423_MosAIC_VF_clean.tsv", quote = F, sep = "\t", row.names = F)

