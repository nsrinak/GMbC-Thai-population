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

# Beta diversity - new - with PCoA ---
relative_matrix <- taxonomic_df_all %>%
  select(Sample, OTU, Relative_Abundance, locality) %>%
  pivot_wider(names_from = OTU, values_from = Relative_Abundance, values_fill = 0) %>% 
  column_to_rownames(var = "Sample")

set.seed(123)
bray_curtis <- vegdist(as.matrix(relative_matrix %>% select(-locality)), method = "bray")

# PCoA using ape ----
pcoa_bray <- ape::pcoa(bray_curtis)

# Variance explained
var_bray <- pcoa_bray$values$Relative_eig[1:2] * 100

# Data frame for plotting
pcoa_bray_df <- data.frame(
  SampleID = row.names(relative_matrix),
  PC1 = pcoa_bray$vectors[, 1],
  PC2 = pcoa_bray$vectors[, 2],
  Locality = relative_matrix$locality
)

# Plot Bray-Curtis PCoA ----
ord_PCoA_bray <- ggplot(pcoa_bray_df, aes(x = PC1, y = PC2, color = Locality, fill = Locality)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    title = "Bray-Curtis",
    x = paste0("PC1 (", round(var_bray[1], 2), "%)"),
    y = paste0("PC2 (", round(var_bray[2], 2), "%)")
  )

physeq_treeroot <- readRDS("C:/Project/5_16s_thai_population/Phylogenetic tree/physeq_treeroot.rds")
unifrac_unweighted <- phyloseq::distance(physeq_treeroot, method = "unifrac", weighted = FALSE)
unifrac_weighted <- phyloseq::distance(physeq_treeroot, method = "unifrac", weighted = TRUE)


# Unweighted UniFrac PCoA ----
pcoa_unifrac_unweighted <- ape::pcoa(unifrac_unweighted)
pcoa_unifrac_unweighted_df <- data.frame(
  SampleID = row.names(sample_data(physeq_treeroot)),
  PC1 = pcoa_unifrac_unweighted$vectors[, 1],
  PC2 = pcoa_unifrac_unweighted$vectors[, 2],
  Locality = sample_data(physeq_treeroot)$locality
)
var_unifrac_unweighted <- pcoa_unifrac_unweighted$values$Relative_eig[1:2] * 100

ord_PCoA_unifrac_unweighted <- ggplot(pcoa_unifrac_unweighted_df, aes(x = PC1, y = PC2, color = Locality, fill = Locality)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    title = "Unweighted UniFrac",
    x = paste0("PC1 (", round(var_unifrac_unweighted[1], 2), "%)"),
    y = paste0("PC2 (", round(var_unifrac_unweighted[2], 2), "%)")
  )

# Weighted UniFrac PCoA ----
pcoa_unifrac_weighted <- ape::pcoa(unifrac_weighted)
pcoa_unifrac_weighted_df <- data.frame(
  SampleID = row.names(sample_data(physeq_treeroot)),
  PC1 = pcoa_unifrac_weighted$vectors[, 1],
  PC2 = pcoa_unifrac_weighted$vectors[, 2],
  Locality = sample_data(physeq_treeroot)$locality
)
var_unifrac_weighted <- pcoa_unifrac_weighted$values$Relative_eig[1:2] * 100

ord_PCoA_unifrac_weighted <- ggplot(pcoa_unifrac_weighted_df, aes(x = PC1, y = PC2, color = Locality, fill = Locality)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    title = "Weighted UniFrac",
    x = paste0("PC1 (", round(var_unifrac_weighted[1], 2), "%)"),
    y = paste0("PC2 (", round(var_unifrac_weighted[2], 2), "%)")
  )


# aitchison distance PCA ----

library(compositions)  # For CLR transformation

# Extract the abundance matrix only
rel_abund <- as.matrix(relative_matrix %>% select(-locality))

# Replace 0s with a small value to avoid log(0)
rel_abund[rel_abund == 0] <- 1e-6

# Apply CLR transformation
clr_abund <- clr(rel_abund)

aitchison_dist <- dist(clr_abund, method = "euclidean")

pca_aitchison <- prcomp(clr_abund, scale. = FALSE)

pca_aitchison_df <- data.frame(
  SampleID = row.names(relative_matrix),
  PC1 = pca_aitchison$x[, 1],
  PC2 = pca_aitchison$x[, 2],
  Locality = relative_matrix$locality  # Use locality for coloring
)

# Calculate variance explained
variance_explained <- pca_aitchison$sdev^2 / sum(pca_aitchison$sdev^2) * 100

# Update the plot with percentage variance on axes
ord_PCA_aitchison <- ggplot(pca_aitchison_df, aes(x = PC1, y = PC2, color = Locality, fill = Locality)) +
  geom_point(size = 2, alpha = 0.7) +  # Points
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  #stat_ellipse(geom = "polygon", alpha = 0.07) +# Filled ellipses
  theme_bw()+
  theme(legend.position = "none")+
  labs(
    title = "Aitchison (PCA)",
    x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 2), "%)")
  )

# Combine all three PCoA plots ----
p_comb_PCoA <- plot_grid(ord_PCoA_bray, ord_PCA_aitchison, 
                         ord_PCoA_unifrac_unweighted, ord_PCoA_unifrac_weighted, 
                         ncol = 2)

ggsave(plot = p_comb_PCoA, 
       filename = "C:/Project/5_16s_thai_population/figure/02092025_PCoA_all_betadivers_axis.png", 
       width = 6.5, height = 4.5, units = "in", dpi = 300)

ggsave(plot = ord_PCA_aitchison,
       filename = "C:/Project/5_16s_thai_population/figure/23072025_PCA_aitchison_betadivers.png",
       width = 7, height = 6)


# Run UMAP on bray distance matrix ----
library(uwot)

relative_matrix <- taxonomic_df_all %>%
  select(Sample, OTU, Relative_Abundance, locality) %>%
  pivot_wider(names_from = OTU, values_from = Relative_Abundance, values_fill = 0) %>% 
  column_to_rownames(var = "Sample")

set.seed(123)
bray_curtis <- vegdist(as.matrix(relative_matrix %>% select(-locality)), method = "bray")

umap_bray_res <- uwot::umap(as.matrix(bray_curtis), seed = 123)

# Convert UMAP results into a dataframe
umap_bray_df <- data.frame(
  SampleID = rownames(relative_matrix),
  UMAP1 = umap_bray_res[,1],
  UMAP2 = umap_bray_res[,2],
  Locality = relative_matrix$locality  # Use locality for coloring
)

ord_umap_bray <- ggplot(umap_bray_df, aes(x = UMAP1, y = UMAP2, color = Locality, fill = Locality)) +
  geom_point(size = 2, alpha = 0.7) +  # Points
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  #stat_ellipse(geom = "polygon", alpha = 0.07) +# Filled ellipses
  theme_bw()+
  theme(legend.position = "none")+
  labs(title = "Bray Curtis")

# Run UMAP on aitchison (CLR norm) ----
umap_aitchison_res <- uwot::umap(as.matrix((clr_abund)), seed = 123)


# Convert UMAP results into a dataframe
umap_aitchison_df <- data.frame(
  SampleID = rownames(as.matrix(clr_abund)),
  UMAP1 = umap_aitchison_res[,1],
  UMAP2 = umap_aitchison_res[,2],
  Locality = relative_matrix$locality # Use locality for coloring
)

ord_umap_aitchison <- ggplot(umap_aitchison_df, aes(x = UMAP1, y = UMAP2, color = Locality, fill = Locality)) +
  geom_point(size = 2, alpha = 0.7) +  # Points
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  #stat_ellipse(geom = "polygon", alpha = 0.07) +# Filled ellipses
  theme_bw()+
  theme(legend.position = "none")+
  labs(title = "CLR normalization")

# Run UMAP on unweighted UniFrac distance matrix ----
umap_unifrac_unweight_res <- uwot::umap(as.matrix(unifrac_unweighted), seed = 123)

# Convert UMAP results into a dataframe
umap_unifrac_unweight_df <- data.frame(
  SampleID = rownames(sample_data(physeq_treeroot)),
  UMAP1 = umap_unifrac_unweight_res[,1],
  UMAP2 = umap_unifrac_unweight_res[,2],
  Locality = sample_data(physeq_treeroot)$locality  # Use locality for coloring
)

ord_umap_unifrac_unweight <- ggplot(umap_unifrac_unweight_df, aes(x = UMAP1, y = UMAP2, color = Locality, fill = Locality)) +
  geom_point(size = 2, alpha = 0.7) +  # Points
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  #stat_ellipse(geom = "polygon", alpha = 0.07) +# Filled ellipses
  theme_bw()+
  theme(legend.position = "none")+
  labs(title = "Unweighted UniFrac")


# Run UMAP on weighted UniFrac distance matrix ----
umap_unifrac_weight_res <- uwot::umap(as.matrix(unifrac_weighted), seed = 123)

# Convert UMAP results into a dataframe
umap_unifrac_weight_df <- data.frame(
  SampleID = rownames(sample_data(physeq_treeroot)),
  UMAP1 = umap_unifrac_weight_res[,1],
  UMAP2 = umap_unifrac_weight_res[,2],
  Locality = sample_data(physeq_treeroot)$locality  # Use locality for coloring
)

ord_umap_unifrac_weight <- ggplot(umap_unifrac_weight_df, aes(x = UMAP1, y = UMAP2, color = Locality, fill = Locality)) +
  geom_point(size = 2, alpha = 0.7) +  # Points
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  #stat_ellipse(geom = "polygon", alpha = 0.07) +  # Filled ellipses
  theme_bw()+
  theme(legend.position = "none")+
  labs(title = "Weighted UniFrac")


p_comb_umap <- plot_grid(ord_umap_bray, ord_umap_unifrac_unweight, ord_umap_unifrac_weight, ncol = 3)

ggsave(plot = p_comb_umap, 
       filename = "C:/Project/5_16s_thai_population/figure/29072025_UMAP_all_betadivers.png", 
       width = 6.5, height = 2.5, units = "in", dpi = 300)


## Brta diversity statistic ----

betadis_aitchison <- betadisper(aitchison_dist, group = relative_matrix$locality)
betadis_bray <- betadisper(bray_curtis, group = relative_matrix$locality)
betadis_unifrac_unweighted <- betadisper(unifrac_unweighted, group = sample_data(physeq_treeroot)$locality)
betadis_unifrac_weighted <- betadisper(unifrac_weighted, group = sample_data(physeq_treeroot)$locality)

perm_betadis_aitchison <- anova(betadis_aitchison, permutations = 10000)
perm_betadis_bray <- anova(betadis_bray, permutations = 10000)
perm_betadis_unifrac_unweighted <- anova(betadis_unifrac_unweighted, permutations = 10000)
perm_betadis_unifrac_weighted <- anova(betadis_unifrac_weighted, permutations = 10000)

permanova_aitchison <- adonis2(aitchison_dist ~ locality, data = relative_matrix, permutations = 10000)
permanova_bray <- adonis2(bray_curtis ~ locality, data = relative_matrix, permutations = 10000)
permanova_unifrac_unweighted <- adonis2(unifrac_unweighted ~ locality, data = a, permutations = 10000)
permanova_unifrac_weighted <- adonis2(unifrac_weighted ~ locality, data = a, permutations = 10000)

library(pairwiseAdonis)

pairwise_aitchison <- pairwise.adonis(aitchison_dist, factors = relative_matrix$locality, perm = 10000)
pairwise_bray <- pairwise.adonis(bray_curtis, factors = relative_matrix$locality, perm = 10000)
pairwise_unifrac_unweighted <- pairwise.adonis(unifrac_unweighted, factors = a$locality, perm = 10000)
pairwise_unifrac_weighted <- pairwise.adonis(unifrac_weighted, factors = a$locality, perm = 10000)
