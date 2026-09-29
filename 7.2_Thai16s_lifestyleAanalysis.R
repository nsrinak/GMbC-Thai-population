library(tidyverse)
library(vegan)
library(phyloseq)
library(cowplot)
library(cluster)
library(ggrepel)

meta <- readRDS("diversityMetadataSampleDF_pairwise.complete.rds")

meta <- meta %>% select(-PD, -PSVs, -Richness, -Shannon, -Simpson, -Evenness) %>% as.data.frame()

# Remove non-lifestyle feature 
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
# Set row names to match sample IDs
rownames(meta) <- meta$Sample
meta_1 <- meta[, -1] %>%   # Remove sample column after setting row names
        select(-locality) 
  
# Compute Gower’s distance
env_gower <- daisy(meta_1, metric = "gower")
env_gower <- as.matrix(env_gower)

# PCoA using ape ----
set.seed(123)
pcoa_gower <- ape::pcoa(env_gower)

# Variance explained
var_gower <- pcoa_gower$values$Relative_eig[1:2] * 100

# Data frame for plotting points
pcoa_gower_df <- data.frame(
  SampleID = rownames(env_gower),
  PC1 = pcoa_gower$vectors[, 1],
  PC2 = pcoa_gower$vectors[, 2],
  Locality = meta$locality
)

pcoa_gower_df <- pcoa_gower_df %>% 
  mutate(Locality = case_when(
    Locality == "bangkok" ~ "Bangkok",
    Locality == "phatthalung" ~ "Phatthalung",
    Locality == "tak" ~ "Tak"
  ))
  
# Biplot Arrow ----

# Isolate numeric variables
# numeric_cols <- sapply(meta_1, is.numeric)
# env_numeric <- meta_1[, numeric_cols, drop = FALSE]

# Calculate correlation between numeric variables and PCoA axes
arrow_cor <- cor(meta_1 %>% 
                   mutate(across(where(is.factor), as.numeric)), 
                 pcoa_gower$vectors[, 1:2], 
                 use = "pairwise.complete.obs")

# Create a data frame for the arrows
arrow_scale <- max(abs(pcoa_gower_df$PC1)) * 0.7  # Adjust 0.7 to make arrows longer or shorter

arrow_df <- data.frame(
  Variable = rownames(arrow_cor),
  xstart = 0,
  ystart = 0,
  xend = arrow_cor[, 1] * arrow_scale,
  yend = arrow_cor[, 2] * arrow_scale
)

# calculate the unscaled correlation length to find true contribution
arrow_df <- arrow_df %>%
  # Filter locality related features
  filter(!Variable %in% c("latitude", "PrecipitationSeasonality")) %>% 
  mutate(
    # Contribution is the Euclidean distance from the origin (0,0)
    contribution = sqrt(xend^2 + yend^2)
  )

write.csv(x = arrow_df, file = "Correlation_arrow.csv", row.names = FALSE)

# Apply filter

# Keep only the top N contributors 
arrow_df_filtered <- arrow_df %>%
  slice_max(order_by = contribution, n = 15)

# Plot PCoA with Biplot Arrows ----
ord_PCoA_gower <- ggplot() +
  # Points layer
  geom_point(data = pcoa_gower_df, aes(x = PC1, y = PC2, color = Locality, fill = Locality), size = 2, alpha = 0.7) +
  
  # Filtered arrows
  geom_segment(data = arrow_df_filtered, aes(x = xstart, y = ystart, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 0.4, alpha = 0.5) +
  
  # REPELLED LABELS
  geom_text_repel(data = arrow_df_filtered, aes(x = xend, y = yend, label = Variable),
                  color = "black", size = 2.25,
                  box.padding = 0,       # Padding around the text box
                  point.padding = 0.5,     # Padding around the xend/yend point
                  segment.color = NA,      # Set to NA so it doesn't draw extra lines from the arrow tip to the text
                  force = 2) +             # Higher numbers push labels further apart
  
  # Rest of your formatting
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  theme_bw(base_size = 9) +
  #coord_equal()+
  labs(
    title = "",
    x = paste0("PC1 (", round(var_gower[1], 2), "%)"),
    y = paste0("PC2 (", round(var_gower[2], 2), "%)")
  )+
  theme(legend.position = "none")

#print(ord_PCoA_gower)

ggsave(plot = ord_PCoA_gower, 
       filename = "p_ord_PCoA_gower_removeGenetic.png", 
       width = 2.3, height = 2.15, units = "in", dpi = 300)


#-------------------------------------------------------------------------------#

asv_table <- readRDS("seqtab_final.rds")
asv_table <- t(asv_table)

taxa_table <- readRDS("tax_final.rds")

meta_table <- read_tsv("metadata.tsv")
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

relative_matrix <- taxonomic_df_all %>%
  select(Sample, OTU, Relative_Abundance, locality) %>%
  pivot_wider(names_from = OTU, values_from = Relative_Abundance, values_fill = 0) %>% 
  column_to_rownames(var = "Sample")

## Aitchison and Gower PC 1 correlation

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

pca_aitchison_df_1 <- pca_aitchison_df %>% rename(PC1_aitchison = PC1,
                                                  PC2_aitchison = PC2)

gower_aitchison <- pcoa_gower_df %>% rename(PC1_gower = PC1,
                                            PC2_gower = PC2) %>% 
  select(-Locality) %>% 
  left_join(., pca_aitchison_df_1, by = "SampleID")

summary(lm(gower_aitchison$PC1_gower ~ gower_aitchison$PC1_aitchison))

#
#Call:
#  lm(formula = gower_aitchison$PC1_gower ~ gower_aitchison$PC1_aitchison)
#
#Residuals:
#  Min        1Q    Median        3Q       Max 
#-0.279221 -0.041672  0.002292  0.038684  0.140678 
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                   -2.010e-17  6.207e-03    0.00        1    
#gower_aitchison$PC1_aitchison  2.753e-03  2.717e-04   10.13   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.0639 on 104 degrees of freedom
#Multiple R-squared:  0.4969,	Adjusted R-squared:  0.492 
#F-statistic: 102.7 on 1 and 104 DF,  p-value: < 2.2e-16

p_gower_aitchison <- ggplot(data = gower_aitchison, aes(x = PC1_gower, y = PC1_aitchison, fill = Locality, colour = Locality))+
  geom_point(shape = 21, size = 2, alpha = 0.7)+
  geom_smooth(aes(group = 1), method = "lm", color ="#999999", se = TRUE)+
  theme_bw(base_size = 9)+
  theme(legend.position = "none")

ggsave(plot = p_gower_aitchison, 
       filename = "p_gower_aitchison_removeGenetic.png", 
       width = 2.25, height = 1.95, units = "in", dpi = 300)