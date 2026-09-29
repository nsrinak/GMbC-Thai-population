library(glmnet)
library(tidyverse)
library(phyloseq)
library(vegan)

## Load taxa data and metadata ----

asv_table <- readRDS("seqtab_final.rds")
asv_table <- t(asv_table)

taxa_table <- readRDS("tax_final.rds")

meta_table <- read_tsv("metadata.tsv")
meta_table <- meta_table %>% column_to_rownames(var = "donor_id") 

## Cleaning metadata ----

# Prepare numeric data frame for correlate with locality

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

# Correlation with NA aware argument use = "pairwise.complete.obs"

cor_matrix <- cor(fil_meta_data, method = "pearson", use = "pairwise.complete.obs")

# Here feature that show no variance (NA) and show high correlation with locality (>0.5) will be filtered out
# They are considered confounding factors

remov_column <- cor_matrix %>% 
  as.data.frame() %>% 
  select(locality) %>% 
  rownames_to_column(var="feature") %>%
  filter(abs(locality) > 0.50 | is.na(locality)) %>% 
  filter(!feature %in% c("locality")) %>% 
  pull(feature)

#cor_matrix %>% 
#  as.data.frame() %>% 
#  select(locality) %>% 
#  rownames_to_column(var="feature") %>%
#  filter(abs(locality) > 0.50 | is.na(locality)) %>% 
#  filter(!feature %in% c("locality")) %>% filter(is.na(locality))

# 21 features have 0 variance

#cor_matrix %>% 
#  as.data.frame() %>% 
#  select(locality) %>% 
#  rownames_to_column(var="feature") %>%
#  filter(abs(locality) > 0.50 | is.na(locality)) %>% 
#  filter(!feature %in% c("locality")) %>% filter(!is.na(locality))

# 31 features have PCC > 0.5

# Features correlated with locality are used to filter the metadata

fil2_meta_data <- meta_table
fil2_meta_data[fil2_meta_data == "na"] <- NA
fil2_meta_data <- fil2_meta_data %>% 
  select(-contains("Dim")) %>%                                 # remove transformed columns (contain "Dim")
  select(where(~ !all(. %in% c(0, NA)) | !is.numeric(.))) %>%  # remove column that all zero
  select(-contains("PC")) %>%                                  # remove column that contain "PC" for now -- this transformation is not bad though
  mutate(across(where(is.character), as.factor)) %>%           # convert character column to factor %>% 
  select(-all_of(remov_column)) %>%                            # remove selected features from correlation analysis
  select(where(~ !is.numeric(.) | mean(. == 0, na.rm = TRUE) <= 0.5)) %>%   # Remove numeric columns with >50% zeros
  select(where(~ !is.factor(.) | mean(. == "na", na.rm = TRUE) <= 0.5)) %>% # Remove factor columns with >50% na
  rownames_to_column(var="Sample")

# All zero removed: 17
# >50% zero removed: 11
# >50% NA removed: 0
# Confounded with locality removed: 52
  # 31 features have PCC > 0.5
    # within 31, 10 were strongly correlated to location
  # 21 features have 0 variance

# Manually correct variable type

to_factor <- c("sex", "c_section", "breast_fed", "tobacco", "admixture")
to_numeric <- c("Yogurt", "Grapes", "Oranges", "Banana", "OtherLeaves", "Eggplant", "Yams", "Manioc", "SweetPotatoes", "SweetHerbs", "Coffee", "DietSoda", "Wine", "Liquor", "Beer", "Cider")  
to_remove <- c("PrecipitationSeasonality", "latitude") 

fil2_meta_data <- fil2_meta_data %>%
  mutate(across(all_of(to_factor), as.factor),
         across(all_of(to_numeric), as.numeric)) %>%
  select(-all_of(to_remove))

# Convert to phyloseq components

ASV <- otu_table(asv_table, taxa_are_rows = TRUE)
TAX <- tax_table(taxa_table)

physeq <- phyloseq(ASV, TAX)

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

writeXStringSet(dna, "asv_sequences_short.fasta")

# Save mapping table for later reconstruction

write.csv(mapping, "asv_mapping.csv", row.names = FALSE)

# Getting tree from bash (muscle > Gblocks > FastTree)

library(ggtree)

tree_fasttree <- read.tree("tree_asv")
ggtree(tree_fasttree, layout = "circular") +
  geom_tiplab(size = 3)

# Load mapping table
# Restore original sequence names
# Attach tree to phyloseq

mapping <- read.csv("asv_mapping.csv", stringsAsFactors = FALSE)
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
#Lower Faith’s PD → The species in the community are closely related (evolutionary constrained).


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
  left_join(., faith_pd_df, by = "Sample") %>% 
  left_join(., psv_df, by = "Sample")

write.csv(alpha_diversity_2, "diversityMetadataSampleDF_pairwise.complete.csv", row.names = FALSE)
write_rds(alpha_diversity_2, "diversityMetadataSampleDF_pairwise.complete.rds")
