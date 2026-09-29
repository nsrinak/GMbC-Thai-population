library(tidyverse)
library(vegan)
library(phyloseq)
library(cowplot)

asv_table <- readRDS("C:/Project/5_16s_thai_population/seqtab_final.rds")
# row = ASV (sequence in this case)
# col = sample ID
asv_table <- t(asv_table)

taxa_table <- readRDS("C:/Project/5_16s_thai_population/tax_final.rds")


meta <- readRDS("C:/Project/5_16s_thai_population/revise_aftermSystems/16082026_diversityMetadataSampleDF_pairwise.complete.rds")
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



## Separate location --
loc_sample <- meta %>% select(Sample, locality) %>% unique()

relative_matrix_bkk <- taxonomic_df_all %>%
  left_join(.,loc_sample, by = "Sample") %>% 
  select(Sample, OTU, Relative_Abundance, locality) %>%
  pivot_wider(names_from = OTU, values_from = Relative_Abundance, values_fill = 0) %>% 
  column_to_rownames(var = "Sample") %>%
  filter(locality =="bangkok") %>% 
  select(-locality)

relative_matrix_phat <- taxonomic_df_all %>%
  left_join(.,loc_sample, by = "Sample") %>% 
  select(Sample, OTU, Relative_Abundance, locality) %>%
  pivot_wider(names_from = OTU, values_from = Relative_Abundance, values_fill = 0) %>% 
  column_to_rownames(var = "Sample") %>%
  filter(locality =="phatthalung") %>% 
  select(-locality)

relative_matrix_tak <- taxonomic_df_all %>%
  left_join(.,loc_sample, by = "Sample") %>% 
  select(Sample, OTU, Relative_Abundance, locality) %>%
  pivot_wider(names_from = OTU, values_from = Relative_Abundance, values_fill = 0) %>% 
  column_to_rownames(var = "Sample") %>%
  filter(locality =="tak") %>% 
  select(-locality)

set.seed(123)
bray_curtis_bkk <- vegdist(as.matrix(relative_matrix_bkk), method = "bray")
bray_curtis_bkk <- as.matrix(bray_curtis_bkk)

bray_curtis_phat <- vegdist(as.matrix(relative_matrix_phat), method = "bray")
bray_curtis_phat <- as.matrix(bray_curtis_phat)

bray_curtis_tak <- vegdist(as.matrix(relative_matrix_tak), method = "bray")
bray_curtis_tak <- as.matrix(bray_curtis_tak)


sample_order_bkk <- row.names(meta %>% filter(locality == "bangkok"))

sample_order_phat <- row.names(meta %>% filter(locality == "phatthalung"))

sample_order_tak <- row.names(meta %>% filter(locality == "tak"))


bray_curtis_bkk <- bray_curtis_bkk[sample_order_bkk, sample_order_bkk]
env_gower_bkk <- env_gower[sample_order_bkk, sample_order_bkk]

bray_curtis_phat <- bray_curtis_phat[sample_order_phat, sample_order_phat]
env_gower_phat <- env_gower[sample_order_phat, sample_order_phat]

bray_curtis_tak <- bray_curtis_tak[sample_order_tak, sample_order_tak]
env_gower_tak <- env_gower[sample_order_tak, sample_order_tak]


mantel_result_bkk <- mantel(bray_curtis_bkk, env_gower_bkk, method = "pearson", permutations = 10000)

mantel_result_phat <- mantel(bray_curtis_phat, env_gower_phat, method = "pearson", permutations = 10000)

mantel_result_tak <- mantel(bray_curtis_tak, env_gower_tak, method = "pearson", permutations = 10000)


print(mantel_result_bkk)
print(mantel_result_phat)
print(mantel_result_tak)

#> print(mantel_result_bkk)
#
#Mantel statistic based on Pearson's product-moment correlation 
#
#Call:
#mantel(xdis = bray_curtis_bkk, ydis = env_gower_bkk, method = "pearson",      permutations = 10000) 
#
#Mantel statistic r: 0.08228 
#      Significance: 0.13119 
#
#Upper quantiles of permutations (null model):
#   90%    95%  97.5%    99% 
#0.0944 0.1191 0.1430 0.1703 
#Permutation: free
#Number of permutations: 10000
#
#> print(mantel_result_phat)
#
#Mantel statistic based on Pearson's product-moment correlation 
#
#Call:
#  mantel(xdis = bray_curtis_phat, ydis = env_gower_phat, method = "pearson",      permutations = 10000) 
#
#Mantel statistic r: 0.03482 
#Significance: 0.29717 
#
#Upper quantiles of permutations (null model):
#  90%    95%  97.5%    99% 
#  0.0927 0.1231 0.1473 0.1753 
#Permutation: free
#Number of permutations: 10000
#
#> print(mantel_result_tak)
#
#Mantel statistic based on Pearson's product-moment correlation 
#
#Call:
#mantel(xdis = bray_curtis_tak, ydis = env_gower_tak, method = "pearson",      permutations = 10000) 
#
#Mantel statistic r: -0.1815 
#      Significance: 0.94881 
#
#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#0.147 0.188 0.224 0.259 
#Permutation: free
#Number of permutations: 10000
#

# Convert distance matrices to vectors
bray_vector_bkk <- bray_curtis_bkk[upper.tri(bray_curtis_bkk)]
env_vector_bkk <- env_gower_bkk[upper.tri(env_gower_bkk)]

bray_vector_phat <- bray_curtis_phat[upper.tri(bray_curtis_phat)]
env_vector_phat <- env_gower_phat[upper.tri(env_gower_phat)]

bray_vector_tak <- bray_curtis_tak[upper.tri(bray_curtis_tak)]
env_vector_tak <- env_gower_tak[upper.tri(env_gower_tak)]

# Create a data frame for plotting
data_bkk <- data.frame(BrayCurtis = bray_vector_bkk, EnvDistance = env_vector_bkk)

data_phat <- data.frame(BrayCurtis = bray_vector_phat, EnvDistance = env_vector_phat)

data_tak <- data.frame(BrayCurtis = bray_vector_tak, EnvDistance = env_vector_tak)

# Create the scatter plot
p_mantel_bray_bkk <- ggplot(data_bkk, aes(x = BrayCurtis, y = EnvDistance)) +
  geom_point(color= "#6699cc", alpha=0.15, size = 1.5) +
  geom_smooth(method = "lm", color ="#999999", se = TRUE) +
  labs(x = "Bray-Curtis Distance", y = " ", title = "Bangkok") +
  theme_bw(base_size = 8)

p_mantel_bray_phat <- ggplot(data_phat, aes(x = BrayCurtis, y = EnvDistance)) +
  geom_point(color= "#6699cc", alpha=0.15, size = 1.5) +
  geom_smooth(method = "lm", color ="#999999", se = TRUE) +
  labs(x = "Bray-Curtis Distance", y = "", title = "Phatthalung") +
  theme_bw(base_size = 8)

p_mantel_bray_tak <- ggplot(data_tak, aes(x = BrayCurtis, y = EnvDistance)) +
  geom_point(color= "#6699cc", alpha=0.15, size = 1.5) +
  geom_smooth(method = "lm", color ="#999999", se = TRUE) +
  labs(x = "Bray-Curtis Distance", y = "", title = "Tak") +
  theme_bw(base_size = 8)

p_metal_bray_all <- plot_grid(p_mantel_bray, p_mantel_bray_bkk, p_mantel_bray_phat, p_mantel_bray_tak, ncol = 4)


# Beta PD ----

physeq_treeroot <- readRDS("C:/Project/5_16s_thai_population/Phylogenetic tree/physeq_treeroot.rds")

unifrac_unweighted <- phyloseq::distance(physeq_treeroot, method = "unifrac", weighted = FALSE)
unifrac_weighted <- phyloseq::distance(physeq_treeroot, method = "unifrac", weighted = TRUE)

set.seed(123)

unifrac_unweighted_bkk <- as.matrix(unifrac_unweighted)
unifrac_unweighted_bkk <- unifrac_unweighted_bkk[sample_order_bkk, sample_order_bkk]

unifrac_unweighted_phat <- as.matrix(unifrac_unweighted)
unifrac_unweighted_phat <- unifrac_unweighted_phat[sample_order_phat, sample_order_phat]

unifrac_unweighted_tak <- as.matrix(unifrac_unweighted)
unifrac_unweighted_tak <- unifrac_unweighted_tak[sample_order_tak, sample_order_tak]


unifrac_weighted_bkk <- as.matrix(unifrac_weighted)
unifrac_weighted_bkk <- unifrac_weighted_bkk[sample_order_bkk, sample_order_bkk]

unifrac_weighted_phat <- as.matrix(unifrac_weighted)
unifrac_weighted_phat <- unifrac_weighted_phat[sample_order_phat, sample_order_phat]

unifrac_weighted_tak <- as.matrix(unifrac_weighted)
unifrac_weighted_tak <- unifrac_weighted_tak[sample_order_tak, sample_order_tak]


mantel_result_unifrac_unweighted_bkk <- mantel(unifrac_unweighted_bkk, env_gower_bkk, method = "pearson", permutations = 10000)
mantel_result_unifrac_unweighted_phat<- mantel(unifrac_unweighted_phat, env_gower_phat, method = "pearson", permutations = 10000)
mantel_result_unifrac_unweighted_tak <- mantel(unifrac_unweighted_tak, env_gower_tak, method = "pearson", permutations = 10000)

print(mantel_result_unifrac_unweighted_bkk)
print(mantel_result_unifrac_unweighted_phat)
print(mantel_result_unifrac_unweighted_tak)

#> print(mantel_result_unifrac_unweighted_bkk)
#
#Mantel statistic based on Pearson's product-moment correlation 
#
#Call:
#mantel(xdis = unifrac_unweighted_bkk, ydis = env_gower_bkk, method = "pearson",      permutations = 10000) 
#
#Mantel statistic r: 0.04319 
#      Significance: 0.30207 
#
#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#0.127 0.174 0.211 0.254 
#Permutation: free
#Number of permutations: 10000
#
#> print(mantel_result_unifrac_unweighted_phat)
#
#Mantel statistic based on Pearson's product-moment correlation 
#
#Call:
#  mantel(xdis = unifrac_unweighted_phat, ydis = env_gower_phat,      method = "pearson", permutations = 10000) 
#
#Mantel statistic r: 0.05297 
#Significance: 0.25277 
#
#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#  0.113 0.150 0.179 0.215 
#Permutation: free
#Number of permutations: 10000
#
#> print(mantel_result_unifrac_unweighted_tak)
#
#Mantel statistic based on Pearson's product-moment correlation 
#
#Call:
#mantel(xdis = unifrac_unweighted_tak, ydis = env_gower_tak, method = "pearson",      permutations = 10000) 
#
#Mantel statistic r: -0.09034 
#      Significance: 0.79872 
#
#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#0.136 0.172 0.201 0.239 
#Permutation: free
#Number of permutations: 10000

mantel_result_unifrac_weighted_bkk <- mantel(unifrac_weighted_bkk, env_gower_bkk, method = "pearson", permutations = 10000)
mantel_result_unifrac_weighted_phat <- mantel(unifrac_weighted_phat, env_gower_phat, method = "pearson", permutations = 10000)
mantel_result_unifrac_weighted_tak <- mantel(unifrac_weighted_tak, env_gower_tak, method = "pearson", permutations = 10000)

print(mantel_result_unifrac_weighted_bkk)
print(mantel_result_unifrac_weighted_phat)
print(mantel_result_unifrac_weighted_tak)

#> print(mantel_result_unifrac_weighted_bkk)
#
#Mantel statistic based on Pearson's product-moment correlation 
#
#Call:
#mantel(xdis = unifrac_weighted_bkk, ydis = env_gower_bkk, method = "pearson",      permutations = 10000) 
#
#Mantel statistic r: -0.04895 
#      Significance: 0.73693 
#
#Upper quantiles of permutations (null model):
#   90%    95%  97.5%    99% 
#0.0938 0.1261 0.1541 0.1846 
#Permutation: free
#Number of permutations: 10000
#
#> print(mantel_result_unifrac_weighted_phat)
#
#Mantel statistic based on Pearson's product-moment correlation 
#
#Call:
#  mantel(xdis = unifrac_weighted_phat, ydis = env_gower_phat, method = "pearson",      permutations = 10000) 
#
#Mantel statistic r: 0.02647 
#Significance: 0.32997 
#
#Upper quantiles of permutations (null model):
#  90%    95%  97.5%    99% 
#  0.0972 0.1279 0.1554 0.1831 
#Permutation: free
#Number of permutations: 10000
#
#> print(mantel_result_unifrac_weighted_tak)
#
#Mantel statistic based on Pearson's product-moment correlation 
#
#Call:
#mantel(xdis = unifrac_weighted_tak, ydis = env_gower_tak, method = "pearson",      permutations = 10000) 
#
#Mantel statistic r: -0.07463 
#      Significance: 0.76512 
#
#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#0.132 0.168 0.197 0.234 
#Permutation: free
#Number of permutations: 10000

unifrac_unweighted_vector_bkk <- unifrac_unweighted_bkk[upper.tri(unifrac_unweighted_bkk)]
unifrac_weighted_vector_bkk <- unifrac_weighted_bkk[upper.tri(unifrac_weighted_bkk)]

unifrac_unweighted_vector_phat <- unifrac_unweighted_phat[upper.tri(unifrac_unweighted_phat)]
unifrac_weighted_vector_phat <- unifrac_weighted_phat[upper.tri(unifrac_weighted_phat)]

unifrac_unweighted_vector_tak <- unifrac_unweighted_tak[upper.tri(unifrac_unweighted_tak)]
unifrac_weighted_vector_tak <- unifrac_weighted_tak[upper.tri(unifrac_weighted_tak)]


data_unifrac_unweighted_bkk <- data.frame(unifrac_unweighted = unifrac_unweighted_vector_bkk, EnvDistance = env_vector_bkk)
data_unifrac_weighted_bkk <- data.frame(unifrac_weighted = unifrac_weighted_vector_bkk, EnvDistance = env_vector_bkk)

data_unifrac_unweighted_phat <- data.frame(unifrac_unweighted = unifrac_unweighted_vector_phat, EnvDistance = env_vector_phat)
data_unifrac_weighted_phat <- data.frame(unifrac_weighted = unifrac_weighted_vector_phat, EnvDistance = env_vector_phat)

data_unifrac_unweighted_tak <- data.frame(unifrac_unweighted = unifrac_unweighted_vector_tak, EnvDistance = env_vector_tak)
data_unifrac_weighted_tak <- data.frame(unifrac_weighted = unifrac_weighted_vector_tak, EnvDistance = env_vector_tak)

p_mantel_unweighted_bkk <- ggplot(data_unifrac_unweighted_bkk, aes(x = unifrac_unweighted, y = EnvDistance)) +
  geom_point(color= "#336600", alpha=0.15, size = 1.5) +
  geom_smooth(method = "lm", color ="#999999", se = TRUE) +
  labs(x = "Unweighted UniFrac Distance", y = " ") +
  theme_bw(base_size = 8)
  #xlim(0.28,0.81)

p_mantel_unweighted_phat <- ggplot(data_unifrac_unweighted_phat, aes(x = unifrac_unweighted, y = EnvDistance)) +
  geom_point(color= "#336600", alpha=0.15, size = 1.5) +
  geom_smooth(method = "lm", color ="#999999", se = TRUE) +
  labs(x = "Unweighted UniFrac Distance", y = "") +
  theme_bw(base_size = 8)
  #xlim(0.28,0.81)

p_mantel_unweighted_tak <- ggplot(data_unifrac_unweighted_tak, aes(x = unifrac_unweighted, y = EnvDistance)) +
  geom_point(color= "#336600", alpha=0.15, size = 1.5) +
  geom_smooth(method = "lm", color ="#999999", se = TRUE) +
  labs(x = "Unweighted UniFrac Distance", y = "") +
  theme_bw(base_size = 8)
  #xlim(0.28,0.81)

p_mantel_weighted_bkk <- ggplot(data_unifrac_weighted_bkk, aes(x = unifrac_weighted, y = EnvDistance)) +
  geom_point(color= "#993300", alpha=0.15, size = 1.5) +
  geom_smooth(method = "lm", color ="#999999", se = TRUE) +
  labs(x = "Weighted UniFrac Distance", y = " ") +
  theme_bw(base_size = 8)

p_mantel_weighted_phat <- ggplot(data_unifrac_weighted_phat, aes(x = unifrac_weighted, y = EnvDistance)) +
  geom_point(color= "#993300", alpha=0.15, size = 1.5) +
  geom_smooth(method = "lm", color ="#999999", se = TRUE) +
  labs(x = "Weighted UniFrac Distance", y = "") +
  theme_bw(base_size = 8)

p_mantel_weighted_tak <- ggplot(data_unifrac_weighted_tak, aes(x = unifrac_weighted, y = EnvDistance)) +
  geom_point(color= "#993300", alpha=0.15, size = 1.5) +
  geom_smooth(method = "lm", color ="#999999", se = TRUE) +
  labs(x = "Weighted UniFrac Distance", y = "") +
  theme_bw(base_size = 8)

p_mantel_weighted_unifrac_all <- plot_grid(p_mantel_weighted, p_mantel_weighted_bkk, p_mantel_weighted_phat, p_mantel_weighted_tak, ncol = 4)

p_mantel_unweighted_unifrac_all <- plot_grid(p_mantel_unweighted, p_mantel_unweighted_bkk, p_mantel_unweighted_phat, p_mantel_unweighted_tak, ncol = 4)


ggsave(plot = p_mantel_weighted_unifrac_all, 
       filename = "C:/Project/5_16s_thai_population/figure/revision_figures/16082026_p_mantel_weighted_unifrac_all.png", 
       width = 15, height = 4)

ggsave(plot = p_mantel_unweighted_unifrac_all, 
       filename = "C:/Project/5_16s_thai_population/figure/revision_figures/16082026_p_mantel_unweighted_unifrac_all.png", 
       width = 15, height = 4)

ggsave(plot = p_metal_bray_all, 
       filename = "C:/Project/5_16s_thai_population/figure/revision_figures/16082026_p_metal_bray_all.png", 
       width = 15, height = 4)


p_comb_separateLoc_mantel <- plot_grid(p_mantel_bray_bkk, p_mantel_bray_phat, p_mantel_bray_tak,
                                       p_mantel_unweighted_bkk, p_mantel_unweighted_phat, p_mantel_unweighted_tak,
                                       p_mantel_weighted_bkk, p_mantel_weighted_phat, p_mantel_weighted_tak,
                                       ncol = 3, rel_heights = c(1.1,1,1))

ggsave(plot = p_comb_separateLoc_mantel,
       filename = "C:/Project/5_16s_thai_population/figure/revision_figures/16082026_p_comb_separateLoc_mantel.png",
       width = 7, height = 7.5, units = "in", dpi = 300)


p_maltel_all_andLocation <- plot_grid(p_metal_bray_all,
                                      p_mantel_weighted_unifrac_all,
                                      p_mantel_unweighted_unifrac_all,
                                      nrow = 3)

ggsave(plot = p_maltel_all_andLocation,
       filename = "C:/Project/5_16s_thai_population/figure/revision_figures/09092026_p_maltel_all_andLocation.png",
       width = 7, height = 5.25, units = "in", dpi = 600)
