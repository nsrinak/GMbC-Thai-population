library(tidyverse)
library(vegan)
library(phyloseq)
library(cowplot)

asv_table <- readRDS("seqtab_final.rds")
asv_table <- t(asv_table)

taxa_table <- readRDS("tax_final.rds")


meta <- readRDS("diversityMetadataSampleDF_pairwise.complete.rds")
meta <- meta %>% select(-PD, -PSVs, -Richness, -Shannon, -Simpson, -Evenness) %>% as.data.frame()

meta <- meta %>% select(-c("sex", "age", "height_cm", "bmi", "calprotectin", 
                           "chromogranin", "igm", "iga", "admixture", "PRS_IBD", "weight_kg" ))

## Phyloseq object ----

# Convert to phyloseq components
ASV <- otu_table(asv_table, taxa_are_rows = TRUE)
TAX <- tax_table(taxa_table)

physeq <- phyloseq(ASV, TAX)

taxonomic_df_all <- psmelt(physeq)

taxonomic_df_all <- taxonomic_df_all %>%
  group_by(Sample) %>%
  mutate(Relative_Abundance = Abundance / sum(Abundance))

## Mantel test ----
library(cluster)

# Set row names to match sample IDs
rownames(meta) <- meta$Sample
meta_1 <- meta[, -1] %>%   # Remove sample column after setting row names
  select(-locality) 

# Compute Gower’s distance
env_gower <- daisy(meta_1, metric = "gower")
env_gower <- as.matrix(env_gower)


set.seed(123)

loc_sample <- meta %>% select(Sample, locality) %>% unique()

relative_matrix <- taxonomic_df_all %>%
  left_join(.,loc_sample, by = "Sample") %>% 
  select(Sample, OTU, Relative_Abundance, locality) %>%
  pivot_wider(names_from = OTU, values_from = Relative_Abundance, values_fill = 0) %>% 
  column_to_rownames(var = "Sample")

bray_curtis_2 <- vegdist(as.matrix(relative_matrix %>% select(-locality)), method = "bray")
bray_curtis_2 <- as.matrix(bray_curtis_2)

sample_order <- row.names(meta_1)

bray_curtis_2 <- bray_curtis_2[sample_order, sample_order]
env_gower <- env_gower[sample_order, sample_order]

mantel_result <- mantel(bray_curtis_2, env_gower, method = "pearson", permutations = 10000)
print(mantel_result)

#Mantel statistic based on Pearson's product-moment correlation 
#
#Call:
#mantel(xdis = bray_curtis_2, ydis = env_gower, method = "pearson",      permutations = 10000) 
#
#Mantel statistic r:  0.15 
#      Significance: 0.0055994 
#
#Upper quantiles of permutations (null model):
#   90%    95%  97.5%    99% 
#0.0752 0.0963 0.1140 0.1370 
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
  labs(title = "All", x = "Bray Curtis Distance", y = "Gower's Distance") +
  theme_bw(base_size = 8)

## Beta PD --

physeq_treeroot <- readRDS("physeq_treeroot.rds")

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
#Mantel statistic r: 0.1335 
#      Significance: 0.013099 
#
#Upper quantiles of permutations (null model):
#   90%    95%  97.5%    99% 
#0.0748 0.0976 0.1168 0.1392 
#Permutation: free
#Number of permutations: 10000
#
#> mantel_result_unifrac_weighted
#
#Mantel statistic based on Pearson's product-moment correlation 
#
#Call:
#mantel(xdis = unifrac_weighted_2, ydis = env_gower, method = "pearson",      permutations = 10000) 
#
#Mantel statistic r: 0.05751 
#      Significance: 0.09779 
#
#Upper quantiles of permutations (null model):
#   90%    95%  97.5%    99% 
#0.0566 0.0768 0.0946 0.1174 
#Permutation: free
#Number of permutations: 10000

unifrac_unweighted_vector <- unifrac_unweighted_2[upper.tri(unifrac_unweighted_2)]
unifrac_weighted_vector <- unifrac_weighted_2[upper.tri(unifrac_weighted_2)]

data_unifrac_unweighted <- data.frame(unifrac_unweighted = unifrac_unweighted_vector, EnvDistance = env_vector)
data_unifrac_weighted <- data.frame(unifrac_weighted = unifrac_weighted_vector, EnvDistance = env_vector)

p_mantel_unweighted <- ggplot(data_unifrac_unweighted, aes(x = unifrac_unweighted, y = EnvDistance)) +
  geom_point(color= "#336600", alpha=0.15, size = 2) +
  geom_smooth(method = "lm", color ="#999999", se = TRUE) +
  labs(x = "Unweighted UniFrac Distance", y = "Gower's Distance") +
  theme_bw(base_size = 8)
  #xlim(0.28,0.81)

p_mantel_weighted <- ggplot(data_unifrac_weighted, aes(x = unifrac_weighted, y = EnvDistance)) +
  geom_point(color= "#993300", alpha=0.15, size = 2) +
  geom_smooth(method = "lm", color ="#999999", se = TRUE) +
  labs(x = "Weighted UniFrac Distance", y = "Gower's Distance") +
  theme_bw(base_size = 8)

p_comb_mantel <- plot_grid(p_mantel_bray, p_mantel_unweighted, p_mantel_weighted, ncol = 2)

ggsave(plot = p_comb_mantel, 
       filename = "corr_mantel_all.png", 
       width = 6.5, height = 4.5, units = "in", dpi = 300)
