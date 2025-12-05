library(glmnet)
library(tidyverse)
library(phyloseq)
library(vegan)

asv_table <- readRDS("C:/Project/5_16s_thai_population/seqtab_final.rds")
# row = ASV (sequence in this case)
# col = sample ID
asv_table <- t(asv_table)

taxa_table <- readRDS("C:/Project/5_16s_thai_population/tax_final.rds")

meta_table <- read_tsv("C:/Project/5_16s_thai_population/thai_16s/metaData/metadata.tsv")
meta_table <- meta_table %>% column_to_rownames(var = "donor_id") 

## Cleaning metadata ----

fil_meta_data <- meta_table

fil_meta_data[fil_meta_data == "na"] <- NA

fil_meta_data <-fil_meta_data %>% 
  select(-contains("Dim")) %>%                                 # remove transformed columns (contain "Dim")
  select(where(~ !all(. %in% c(0, NA)) | !is.numeric(.))) %>%  # remove column that all zero
  select(-contains("PC")) %>%                                  # remove column that contain "PC" for now 
  select(where(~ !is.numeric(.) | mean(. == 0, na.rm = TRUE) <= 0.5)) %>%   # Remove numeric columns with >50% zeros
  select(where(~ !is.factor(.) | mean(. == NA, na.rm = TRUE) <= 0.5)) %>% # Remove factor columns with >50% na
  mutate(across(where(is.character), as.factor)) %>%           # convert character column to factor
  mutate(across(where(is.factor), as.numeric))

glimpse(fil_meta_data)

cor_matrix <- cor(fil_meta_data, method = "pearson")

remov_column <- cor_matrix %>% 
  as.data.frame() %>% 
  select(locality) %>% 
  rownames_to_column(var="feature") %>%
  filter(abs(locality) > 0.5 | is.na(locality)) %>% 
  filter(!feature %in% c("locality", "sex")) %>% 
  pull(feature)


fil2_meta_data <- meta_table %>% 
  select(-contains("Dim")) %>%                                 # remove transformed columns (contain "Dim")
  select(where(~ !all(. %in% c(0, NA)) | !is.numeric(.))) %>%  # remove column that all zero
  select(-contains("PC")) %>%                                  # remove column that contain "PC" for now -- this transformation is not bad though
  mutate(across(where(is.character), as.factor)) %>%           # convert character column to factor %>% 
  select(-c(remov_column)) %>%                                 # remove selected features
  select(where(~ !is.numeric(.) | mean(. == 0, na.rm = TRUE) <= 0.5)) %>%   # Remove numeric columns with >50% zeros
  select(where(~ !is.factor(.) | mean(. == "na", na.rm = TRUE) <= 0.5)) %>% # Remove factor columns with >50% na
  rownames_to_column(var="Sample")

glimpse(fil2_meta_data)

# Convert to phyloseq components
ASV <- otu_table(asv_table, taxa_are_rows = TRUE)
TAX <- tax_table(taxa_table)
SAMP <- sample_data(meta_table)

physeq <- phyloseq(ASV, TAX, SAMP)

taxonomic_df_all <- psmelt(physeq)


## Phylogenetic diversity ----

## Merge phylogenetic tree ----

# Extract original taxa names (sequences)
# Create mapping table for later use
# Convert to DNAStringSet

sequences <- taxa_names(physeq)  # These are the actual sequences
asv_ids <- paste0("ASV", seq_along(sequences))  # Create short IDs

mapping <- data.frame(ASV_ID = asv_ids, Sequence = sequences)

dna <- DNAStringSet(sequences)
names(dna) <- asv_ids  # Assign short ASV IDs as names

# Save to FASTA
writeXStringSet(dna, "C:/Project/5_16s_thai_population/asv_sequences_short.fasta")
# Save mapping table for later reconstruction
write.csv(mapping, "C:/Project/5_16s_thai_population/asv_mapping.csv", row.names = FALSE)

# Getting tree from bash (muscle > Gblocks > FastTree)

library(ggtree)

tree_fasttree <- read.tree("C:/Project/5_16s_thai_population/Phylogenetic tree/tree_asv")
ggtree(tree_fasttree, layout = "circular") +
  geom_tiplab(size = 3)

# Load mapping table
# Restore original sequence names
# Attach tree to phyloseq

mapping <- read.csv("C:/Project/5_16s_thai_population/Phylogenetic tree/asv_mapping.csv", stringsAsFactors = FALSE)
tree_fasttree$tip.label <- mapping$Sequence[match(tree_fasttree$tip.label, mapping$ASV_ID)]
physeq_tree <- merge_phyloseq(physeq, tree_fasttree)

# Get the phylogenetic tree
# Root the tree at the midpoint
# Assign the rooted tree back to phyloseq
library(phangorn)

tree_rooted <- midpoint(tree_fasttree)
physeq_treeroot <- merge_phyloseq(physeq, tree_rooted)

## Phylogenetic diversity analysis ----
library(picante)

# Alpha PD
meta_table_1 <- meta_table %>% rownames_to_column(var = "donor_id")
faith_pd <- pd(t(otu_table(physeq_treeroot)), phy_tree(physeq_treeroot), include.root = TRUE)
psv_value <- psv(t(otu_table(physeq_treeroot)), phy_tree(physeq_treeroot))

#Higher Faith’s PD → More evolutionary history is preserved in the sample.
#Lower Faith’s PD → The species in the community are closely related (evolutionarily constrained).


# Calculate alpha diversities and merge phylogenetic diversity

alpha_diversity <- taxonomic_df_all %>%
  group_by(Sample) %>%
  summarise(
    Richness = n_distinct(OTU[Abundance > 0]),
    Shannon = vegan::diversity(Abundance, index = "shannon"),
    Simpson = vegan::diversity(Abundance, index = "simpson"),
    Evenness = Shannon / log(Richness),
    .groups = "drop"
  ) %>% 
  left_join(., fil2_meta_data, by = "Sample")

faith_pd_df <- faith_pd %>% as.data.frame() %>% rownames_to_column(var = "Sample") %>% select(Sample, PD)
psv_df <- psv_value %>% as.data.frame() %>% rownames_to_column(var = "Sample") %>% select(Sample, PSVs)

alpha_diversity_2 <- alpha_diversity %>% 
  #select(-c(locality, latitude, country)) %>% 
  left_join(., faith_pd_df, by = "Sample") %>% 
  left_join(., psv_df, by = "Sample")

alpha_diversity_2[alpha_diversity_2 == "na"] <- NA

write.csv(alpha_diversity_2, "C:/Project/5_16s_thai_population/22052005_diversityMetadataSampleDF.csv", row.names = FALSE)
