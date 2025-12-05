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



## Separate location --

relative_matrix_bkk <- taxonomic_df_all %>%
  select(Sample, OTU, Relative_Abundance, locality) %>%
  pivot_wider(names_from = OTU, values_from = Relative_Abundance, values_fill = 0) %>% 
  column_to_rownames(var = "Sample") %>%
  filter(locality =="bangkok") %>% 
  select(-locality)

relative_matrix_phat <- taxonomic_df_all %>%
  select(Sample, OTU, Relative_Abundance, locality) %>%
  pivot_wider(names_from = OTU, values_from = Relative_Abundance, values_fill = 0) %>% 
  column_to_rownames(var = "Sample") %>%
  filter(locality =="phatthalung") %>% 
  select(-locality)

relative_matrix_tak <- taxonomic_df_all %>%
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


sample_order_bkk <- row.names(fil2_meta_data2 %>% filter(locality == "bangkok"))

sample_order_phat <- row.names(fil2_meta_data2 %>% filter(locality == "phatthalung"))

sample_order_tak <- row.names(fil2_meta_data2 %>% filter(locality == "tak"))


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
#Mantel statistic r: 0.06644 
#      Significance: 0.17508 
#
#Upper quantiles of permutations (null model):
#   90%    95%  97.5%    99% 
#0.0914 0.1162 0.1373 0.1602 
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
#Mantel statistic r: 0.08164 
#Significance: 0.13529 
#
#Upper quantiles of permutations (null model):
#  90%    95%  97.5%    99% 
#  0.0963 0.1260 0.1501 0.1770 
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
#Mantel statistic r: -0.1469 
#      Significance: 0.91741 
#
#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#0.139 0.176 0.210 0.243 
#Permutation: free
#Number of permutations: 10000

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
  labs(x = "Bray-Curtis Distance", y = "Gower's Distance", title = "Bangkok") +
  theme_bw()

p_mantel_bray_phat <- ggplot(data_phat, aes(x = BrayCurtis, y = EnvDistance)) +
  geom_point(color= "#6699cc", alpha=0.15, size = 1.5) +
  geom_smooth(method = "lm", color ="#999999", se = TRUE) +
  labs(x = "Bray-Curtis Distance", y = "", title = "Phatthalung") +
  theme_bw()

p_mantel_bray_tak <- ggplot(data_tak, aes(x = BrayCurtis, y = EnvDistance)) +
  geom_point(color= "#6699cc", alpha=0.15, size = 1.5) +
  geom_smooth(method = "lm", color ="#999999", se = TRUE) +
  labs(x = "Bray-Curtis Distance", y = "", title = "Tak") +
  theme_bw()

p_metal_bray_all <- plot_grid(p_mantel_bray_bkk, p_mantel_bray_phat, p_mantel_bray_tak, ncol = 3)


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

#Mantel statistic based on Pearson's product-moment correlation 
#
#Call:
#mantel(xdis = unifrac_unweighted_bkk, ydis = env_gower_bkk, method = "pearson",      permutations = 10000) 
#
#Mantel statistic r: 0.0483 
#      Significance: 0.28717 
#
#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#0.122 0.161 0.197 0.241 
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
#Mantel statistic r: 0.1086 
#Significance: 0.11469 
#
#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#  0.116 0.149 0.177 0.212 
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
#Mantel statistic r: -0.06226 
#      Significance: 0.73033 
#
#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#0.129 0.164 0.192 0.228 
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
#Mantel statistic r: -0.0634 
#      Significance: 0.82442 
#
#Upper quantiles of permutations (null model):
#   90%    95%  97.5%    99% 
#0.0918 0.1214 0.1479 0.1803 
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
#Mantel statistic r: 0.05289 
#Significance: 0.23548 
#
#Upper quantiles of permutations (null model):
#  90%    95%  97.5%    99% 
#  0.0993 0.1281 0.1531 0.1849 
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
#Mantel statistic r: -0.04399 
#      Significance: 0.67243 
#
#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#0.123 0.160 0.188 0.219 
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
  labs(x = "Unweighted UniFrac Distance", y = "Gower's Distance") +
  theme_bw()
  #xlim(0.28,0.81)

p_mantel_unweighted_phat <- ggplot(data_unifrac_unweighted_phat, aes(x = unifrac_unweighted, y = EnvDistance)) +
  geom_point(color= "#336600", alpha=0.15, size = 1.5) +
  geom_smooth(method = "lm", color ="#999999", se = TRUE) +
  labs(x = "Unweighted UniFrac Distance", y = "") +
  theme_bw()
  #xlim(0.28,0.81)

p_mantel_unweighted_tak <- ggplot(data_unifrac_unweighted_tak, aes(x = unifrac_unweighted, y = EnvDistance)) +
  geom_point(color= "#336600", alpha=0.15, size = 1.5) +
  geom_smooth(method = "lm", color ="#999999", se = TRUE) +
  labs(x = "Unweighted UniFrac Distance", y = "") +
  theme_bw()
  #xlim(0.28,0.81)

p_mantel_weighted_bkk <- ggplot(data_unifrac_weighted_bkk, aes(x = unifrac_weighted, y = EnvDistance)) +
  geom_point(color= "#993300", alpha=0.15, size = 1.5) +
  geom_smooth(method = "lm", color ="#999999", se = TRUE) +
  labs(x = "Weighted UniFrac Distance", y = "Gower's Distance") +
  theme_bw()

p_mantel_weighted_phat <- ggplot(data_unifrac_weighted_phat, aes(x = unifrac_weighted, y = EnvDistance)) +
  geom_point(color= "#993300", alpha=0.15, size = 1.5) +
  geom_smooth(method = "lm", color ="#999999", se = TRUE) +
  labs(x = "Weighted UniFrac Distance", y = "") +
  theme_bw()

p_mantel_weighted_tak <- ggplot(data_unifrac_weighted_tak, aes(x = unifrac_weighted, y = EnvDistance)) +
  geom_point(color= "#993300", alpha=0.15, size = 1.5) +
  geom_smooth(method = "lm", color ="#999999", se = TRUE) +
  labs(x = "Weighted UniFrac Distance", y = "") +
  theme_bw()

p_mantel_weighted_unifrac_all <- plot_grid(p_mantel_weighted_bkk, p_mantel_weighted_phat, p_mantel_weighted_tak, ncol = 3)

p_mantel_unweighted_unifrac_all <- plot_grid(p_mantel_unweighted_bkk, p_mantel_unweighted_phat, p_mantel_unweighted_tak, ncol = 3)


ggsave(plot = p_mantel_weighted_unifrac_all, 
       filename = "C:/Project/5_16s_thai_population/figure/23072025_p_mantel_weighted_unifrac_all.png", 
       width = 15, height = 4)

ggsave(plot = p_mantel_unweighted_unifrac_all, 
       filename = "C:/Project/5_16s_thai_population/figure/23072025_p_mantel_unweighted_unifrac_all.png", 
       width = 15, height = 4)

ggsave(plot = p_metal_bray_all, 
       filename = "C:/Project/5_16s_thai_population/figure/23072025_p_metal_bray_all.png", 
       width = 15, height = 4)


p_comb_separateLoc_mantel <- plot_grid(p_mantel_bray_bkk, p_mantel_bray_phat, p_mantel_bray_tak,
                                       p_mantel_unweighted_bkk, p_mantel_unweighted_phat, p_mantel_unweighted_tak,
                                       p_mantel_weighted_bkk, p_mantel_weighted_phat, p_mantel_weighted_tak,
                                       ncol = 3, rel_heights = c(1.1,1,1))

ggsave(plot = p_comb_separateLoc_mantel,
       filename = "C:/Project/5_16s_thai_population/figure/03092025_p_comb_separateLoc_mantel.png",
       width = 7, height = 7.5, units = "in", dpi = 300)
