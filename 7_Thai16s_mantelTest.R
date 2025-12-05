library(tidyverse)
library(vegan)
library(phyloseq)
library(cowplot)

asv_table <- readRDS("C:/Project/5_16s_thai_population/seqtab_final.rds")
# row = ASV (sequence in this case)
# col = sample ID
asv_table <- t(asv_table)

taxa_table <- readRDS("C:/Project/5_16s_thai_population/tax_final.rds")


meta_table <- read_tsv("C:/Project/5_16s_thai_population/thai_16s/metaData/metadata.tsv")
meta_table <- meta_table %>% column_to_rownames(var = "donor_id") 
## Phyloseq object ----

# Convert to phyloseq components
ASV <- otu_table(asv_table, taxa_are_rows = TRUE)
TAX <- tax_table(taxa_table)
SAMP <- sample_data(meta_table)

physeq <- phyloseq(ASV, TAX, SAMP)

taxonomic_df_all <- psmelt(physeq)

taxonomic_df_all <- taxonomic_df_all %>%
  group_by(Sample) %>%
  mutate(Relative_Abundance = Abundance / sum(Abundance))

## Cleaning metadata ----

fil_meta_data <- meta_table

fil_meta_data[fil_meta_data == "na"] <- NA

fil_meta_data <-fil_meta_data %>% 
  select(-contains("Dim")) %>%                                 # remove transformed columns (contain "Dim")
  select(where(~ !all(. %in% c(0, NA)) | !is.numeric(.))) %>%  # remove column that all zero
  select(-contains("PC")) %>%                                  # remove column that contain "PC" for now -- this transformation is not bad though
  select(where(~ !is.numeric(.) | mean(. == 0, na.rm = TRUE) <= 0.5)) %>%   # Remove numeric columns with >50% zeros
  select(where(~ !is.factor(.) | mean(. == NA, na.rm = TRUE) <= 0.5)) %>% # Remove factor columns with >50% na
  mutate(across(where(is.character), as.factor)) %>%           # convert character column to factor
  mutate(across(where(is.factor), as.numeric))

glimpse(fil_meta_data)

cor_matrix <- cor(fil_meta_data, method = "pearson")
cor
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

## Mantel test ----
library(cluster)

fil2_meta_data %>% glimpse()

# Set row names to match sample IDs
rownames(fil2_meta_data) <- fil2_meta_data$Sample
fil2_meta_data2 <- fil2_meta_data[, -1]  # Remove sample column after setting row names

# Compute Gower’s distance
env_gower <- daisy(fil2_meta_data2, metric = "gower")

# Convert to matrix for visualization
env_gower<- as.matrix(env_gower)

# Run UMAP on bray distance matrix ----
umap_gower <- uwot::umap(env_gower, seed = 123)

# Convert UMAP results into a dataframe
umap_gower_df <- data.frame(
  SampleID = rownames(env_gower),
  UMAP1 = umap_gower[,1],
  UMAP2 = umap_gower[,2],
  Locality = fil2_meta_data$locality  # Use locality for coloring
)

ord_umap_gower <- ggplot(umap_gower_df, aes(x = UMAP1, y = UMAP2, color = Locality, fill = Locality)) +
  geom_point(size = 2, alpha = 0.7) +  # Points
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  #stat_ellipse(geom = "polygon", alpha = 0.07) +# Filled ellipses
  theme_bw()+
  theme(legend.position = "none")+
  labs(title = "")


# PCoA using ape ----
set.seed(123)
pcoa_gower <- ape::pcoa(env_gower)

# Variance explained
var_gower <- pcoa_gower$values$Relative_eig[1:2] * 100

# Data frame for plotting
pcoa_gower_df <- data.frame(
  SampleID = rownames(env_gower),
  PC1 = pcoa_gower$vectors[, 1],
  PC2 = pcoa_gower$vectors[, 2],
  Locality = fil2_meta_data$locality
)

# Plot Bray-Curtis PCoA ----
ord_PCoA_gower <- ggplot(pcoa_gower_df, aes(x = PC1, y = PC2, color = Locality, fill = Locality)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    title = "",
    x = paste0("PC1 (", round(var_gower[1], 2), "%)"),
    y = paste0("PC2 (", round(var_gower[2], 2), "%)")
  )

ggsave(plot = ord_PCoA_gower, 
       filename = "C:/Project/5_16s_thai_population/figure/03092025_ord_PCoA_gower.png", 
       width = 6.5/2, height = 4.5/2, units = "in", dpi = 300)

p_comb_gower_lifestyle <- plot_grid(ord_PCoA_gower, ord_umap_gower, ncol = 2)


ggsave(plot = p_comb_gower_lifestyle,
       filename = "C:/Project/5_16s_thai_population/figure/29072025_p_comb_gower_lifestyle.png",
       width = 4.33, height = 2.5, units = "in", dpi = 300)

ggsave(plot = ord_umap_gower,
       filename = "C:/Project/5_16s_thai_population/figure/23072025_UMAP_Gower_lifestyle.png",
       height = 10, width = 12, units = "cm")

set.seed(123)

relative_matrix <- taxonomic_df_all %>%
  select(Sample, OTU, Relative_Abundance, locality) %>%
  pivot_wider(names_from = OTU, values_from = Relative_Abundance, values_fill = 0) %>% 
  column_to_rownames(var = "Sample")

bray_curtis_2 <- vegdist(as.matrix(relative_matrix %>% select(-locality)), method = "bray")
bray_curtis_2 <- as.matrix(bray_curtis_2)

sample_order <- row.names(fil2_meta_data2)

bray_curtis_2 <- bray_curtis_2[sample_order, sample_order]
env_gower <- env_gower[sample_order, sample_order]

mantel_result <- mantel(bray_curtis_2, env_gower, method = "pearson", permutations = 10000)
print(mantel_result)

mantel
#> mantel_result

#Mantel statistic based on Pearson's product-moment correlation 

#Call:
#mantel(xdis = bray_curtis_2, ydis = env_gower, method = "pearson",      permutations = 10000) 

#Mantel statistic r: 0.1616 
#      Significance: 0.00059994 

#Upper quantiles of permutations (null model):
#   90%    95%  97.5%    99% 
#0.0584 0.0764 0.0906 0.1076 
#Permutation: free
#Number of permutations: 10000

# Convert distance matrices to vectors
bray_vector <- bray_curtis_2[upper.tri(bray_curtis_2)]
env_vector <- env_gower[upper.tri(env_gower)]

# Create a data frame for plotting
data <- data.frame(BrayCurtis = bray_vector, EnvDistance = env_vector)

# Create the scatter plot
p_mantel_bray <- ggplot(data, aes(x = BrayCurtis, y = EnvDistance)) +
  geom_point(color= "#6699cc", alpha=0.15, size = 2) +
  geom_smooth(method = "lm", color ="#999999", se = TRUE) +
  labs(title = " ", x = "Bray Curtis Distance", y = "Gower's Distance") +
  theme_bw()

## Beta PD --

physeq_treeroot <- readRDS("C:/Project/5_16s_thai_population/Phylogenetic tree/physeq_treeroot.rds")

unifrac_unweighted <- phyloseq::distance(physeq_treeroot, method = "unifrac", weighted = FALSE)
unifrac_weighted <- phyloseq::distance(physeq_treeroot, method = "unifrac", weighted = TRUE)

unifrac_unweighted_2 <- as.matrix(unifrac_unweighted)
unifrac_unweighted_2 <- unifrac_unweighted_2[sample_order, sample_order]

unifrac_weighted_2 <- as.matrix(unifrac_weighted)
unifrac_weighted_2 <- unifrac_weighted_2[sample_order, sample_order]

mantel_result_unifrac_unweighted <- mantel(unifrac_unweighted_2, env_gower, method = "pearson", permutations = 10000)
mantel_result_unifrac_weighted <- mantel(unifrac_weighted_2, env_gower, method = "pearson", permutations = 10000)

#> mantel_result_unifrac_unweighted
#
#Mantel statistic based on Pearson's product-moment correlation 
#
#Call:
#mantel(xdis = unifrac_unweighted_2, ydis = env_gower, method = "pearson",      permutations = 10000) 
#
#Mantel statistic r: 0.1993 
#      Significance: 9.999e-05 
#
#Upper quantiles of permutations (null model):
#   90%    95%  97.5%    99% 
#0.0604 0.0791 0.0961 0.1146 
#Permutation: free
#Number of permutations: 10000

#> mantel_result_unifrac_weighted
#
#Mantel statistic based on Pearson's product-moment correlation 
#
#Call:
#  mantel(xdis = unifrac_weighted_2, ydis = env_gower, method = "pearson",      permutations = 10000) 
#
#Mantel statistic r: 0.07565 
#Significance: 0.027297 
#
#Upper quantiles of permutations (null model):
#  90%    95%  97.5%    99% 
#  0.0478 0.0642 0.0774 0.0957 
#Permutation: free
#Number of permutations: 10000

unifrac_unweighted_vector <- unifrac_unweighted_2[upper.tri(unifrac_unweighted_2)]
unifrac_weighted_vector <- unifrac_weighted_2[upper.tri(unifrac_weighted_2)]

data_unifrac_unweighted <- data.frame(unifrac_unweighted = unifrac_unweighted_vector, EnvDistance = env_vector)
data_unifrac_weighted <- data.frame(unifrac_weighted = unifrac_weighted_vector, EnvDistance = env_vector)

p_mantel_unweighted <- ggplot(data_unifrac_unweighted, aes(x = unifrac_unweighted, y = EnvDistance)) +
  geom_point(color= "#336600", alpha=0.15, size = 2) +
  geom_smooth(method = "lm", color ="#999999", se = TRUE) +
  labs(title = " ", x = "Unweighted UniFrac Distance", y = " ") +
  theme_bw()
  #xlim(0.28,0.81)

p_mantel_weighted <- ggplot(data_unifrac_weighted, aes(x = unifrac_weighted, y = EnvDistance)) +
  geom_point(color= "#993300", alpha=0.15, size = 2) +
  geom_smooth(method = "lm", color ="#999999", se = TRUE) +
  labs(title = " ", x = "Weighted UniFrac Distance", y = "Gower's Distance") +
  theme_bw()

p_comb_mantel <- plot_grid(p_mantel_bray, p_mantel_unweighted, p_mantel_weighted, ncol = 2)

ggsave(plot = p_comb_mantel, 
       filename = "C:/Project/5_16s_thai_population/figure/29072025_corr_mantel_all.png", 
       width = 6.5, height = 4.5, units = "in", dpi = 300)
