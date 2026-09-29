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

taxaTable_combine <- read.csv(file = "C:/Project/5_16s_thai_population/13082026_metaAnalysis_combinedTaxaTable.csv")

remove_sample <- taxaTable_combine %>% 
  group_by(Sample, Study, locality) %>% 
  summarise(Abundance = sum(Abundance),
            count = n()) %>% 
  filter(Abundance < 1000) %>% pull(Sample)

## Phyloseq object ----

# Convert to phyloseq components
ASV <- otu_table(asv_table, taxa_are_rows = TRUE)
TAX <- tax_table(taxa_table)
SAMP <- sample_data(meta_table)

physeq <- phyloseq(ASV, TAX, SAMP)

sample_data(physeq)

taxonomic_df_all <- psmelt(physeq)

## Get representative ASV for a each genus and extract sequence for later building phylogenetic tree ----

taxa_table_df <- taxa_table %>% as.data.frame() %>% rownames_to_column(var = "ASV") %>% glimpse()

# There are multiple ASVs in a genus due to limitation of 16s rDNA sequencing technique 
# Across all samples, we selected the representative AVS that has the highest sum relative abundance
asv_table <- asv_table %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ASV") %>% 
  pivot_longer(cols = -ASV, names_to = "sample", values_to = "count") %>% 
  left_join(., taxa_table_df, by = "ASV") %>% 
  group_by(sample) %>% 
  mutate(ra = count / sum(count)) %>% 
  ungroup() %>% 
  group_by(ASV) %>%
  mutate(sumRA = sum(ra)) %>% 
  ungroup() %>% 
  group_by(Genus) %>% 
  mutate(max_sumRA = if_else(sumRA == max(sumRA), 1, 0)) %>% 
  filter(max_sumRA == 1) %>% 
  filter(!is.na(Genus)) %>% 
  select(Genus, ASV) %>% unique()

# Extract the sequences
library(Biostrings)
dna <- DNAStringSet(asv_table$ASV)
names(dna) <- asv_table$Genus 

writeXStringSet(dna, "C:/Project/5_16s_thai_population/30072025_asv_genus_sequences_short.fasta")

## Plots of prevalence analysis ----

# STILL UNSURE IF I AM GOING TO PUT THIS IN MS

thaiDataset <- taxaTable_combine %>% 
  filter(!Sample %in% remove_sample, 
         Study != "This study", Genus %in% unique(taxonomic_df_all$Genus)) %>% 
  select(OTU, Sample, Abundance, locality, Genus) %>% 
  mutate(locality = "Thai datasets")

thaiDataset %>% glimpse()

# Count appearance (1, 2, or all location) of microbial genera 
p_appearance <- taxonomic_df_all %>% 
  select(OTU, Sample, Abundance, locality, Genus) %>% 
  rbind(thaiDataset) %>% 
  filter(Abundance != 0 & !is.na(Genus)) %>% 
  select(locality, Genus) %>% 
  unique() %>% 
  group_by(Genus) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  mutate(count = as.character(count)) %>% 
  group_by(count) %>% 
  summarise(count_count = n()) %>% 
  ggplot(aes(x = count, y = count_count))+
  geom_bar(stat = "identity", fill = "grey")+
  theme_bw()+
  labs(x = "Appearance",
       y = "Number of genus")+
  theme(plot.margin = margin(t = 0, r = 25, b = 0, l = 10, unit = "pt"))

# Prevalence of each microbial genus in each location
p_prevalance <- taxonomic_df_all %>% 
  select(OTU, Sample, Abundance, locality, Genus) %>% 
  rbind(thaiDataset) %>% 
  filter(Abundance != 0 & !is.na(Genus)) %>%
  select(Sample, Genus, locality) %>% 
  unique() %>% 
  group_by(Genus, locality) %>% 
  summarise(count_genus = n()) %>% 
  ungroup() %>% 
  mutate(count_found_loc = if_else(locality == "bangkok", count_genus/47,
                                   if_else(locality == "phatthalung", count_genus/29, 
                                           if_else(locality == "tak", count_genus/30, count_genus/340)))) %>% 
  ggplot(aes(x = locality, y = count_found_loc, fill = locality))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(width = 0.25, alpha=0.2)+
  theme_bw()+
  labs(y = "Prevalence")+
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 310, hjust = 0),
        axis.title.x = element_blank(),
        plot.margin = margin(t = 5, r = 25, b = 0, l = 10, unit = "pt"))+
  scale_fill_manual(values = c("bangkok"= "#F8766D", "tak" = "#619CFF", "phatthalung" = "#00BA38", "Thai datasets" = "#FFB90F"))+
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))

# Plot prevalence against and appearance separated by locations
p_app_pre <- taxonomic_df_all %>% 
  select(OTU, Sample, Abundance, locality, Genus) %>% 
  rbind(thaiDataset) %>% 
  filter(Abundance != 0 & !is.na(Genus)) %>%
  select(Sample, Genus, locality) %>% 
  unique() %>% 
  group_by(Genus, locality) %>% 
  summarise(count_genus = n()) %>% 
  ungroup() %>% 
  mutate(count_found_loc = if_else(locality == "bangkok", count_genus/47,
                                   if_else(locality == "phatthalung", count_genus/29, 
                                           if_else(locality == "tak", count_genus/30, count_genus/340)))) %>%
  group_by(Genus) %>% 
  mutate(appear = n()) %>% 
  ggplot(aes(x = as.factor(appear), y = count_found_loc, fill = locality))+
  geom_boxplot()+
  theme_bw()+
  labs(x = "Appearance",
       y = "Prevalence")+
  scale_fill_manual(values = c("bangkok"= "#F8766D", "tak" = "#619CFF", "phatthalung" = "#00BA38", "Thai datasets" = "#FFB90F"))+
  theme(legend.position = "none",
        plot.margin = margin(t = 0, r = 25, b = 0, l = 10, unit = "pt"))

p_comb_prevalence <- plot_grid(p_pearson_prev, p_appearance, p_prevalance, p_app_pre, ncol = 1, rel_heights = c(1, 0.85,1,0.85))

ggsave(plot = p_comb_prevalence,
       filename = "C:/Project/5_16s_thai_population/figure/revision_figures/08092026_p_comb_prevalence_limitThaiAllGenus.png", 
       width = 2.15, height = 7, units = "in")

## Phylo genetic tree plot ----

library(ggtree)
library(ggtreeExtra)
library(stringr)

# Phylogenetic tree was generated in Bash.
# Multiple sequence alignment was conducted using MUSCLE v5.1 and Gblocks v0.91b
# was used to select high conserve aligned blocks.
# Finally, FastTree v2.2 was used to construct the phylogenetic tree.
# The script for this was attached separately.

tree_fasttree <- read.tree("C:/Project/5_16s_thai_population/Phylogenetic tree/29082025_tree_asv_genus_nogap-gb_boot")
tree_fasttree$tip.label

a <- taxonomic_df_all %>% 
  filter(Abundance != 0 & !is.na(Genus)) %>%
  select(Sample, Genus, locality, Phylum) %>% 
  unique() %>% 
  group_by(Phylum, Genus, locality) %>% 
  summarise(count_genus = n(), .groups = "drop") %>% 
  ungroup() %>% 
  mutate(count_found_loc = if_else(locality == "bangkok", count_genus/47,
                                   if_else(locality == "phatthalung", count_genus/29, count_genus/30)))%>% 
  mutate(genus_rmvspace = str_remove_all(Genus, "\\[.*?\\]")) %>% 
  mutate(genus_rmvspace_1 = str_replace_all(genus_rmvspace, " ", "_")) %>% 
  mutate(genus_rmvspace_1 = factor(genus_rmvspace_1, levels = tree_fasttree$tip.label))%>% 
  arrange(genus_rmvspace_1)

metadata_prevalence <- taxaTable_combine %>% 
  filter(!Sample %in% remove_sample) %>% 
  filter(Abundance != 0 & !is.na(Genus)) %>%
  select(Sample, Genus, Phylum) %>% 
  unique() %>% 
  group_by(Phylum, Genus) %>% 
  summarise(count_genus = n(), .groups = "drop") %>% 
  mutate(prev_all = count_genus/446) %>% 
  mutate(genus_rmvspace = str_remove_all(Genus, "\\[.*?\\]")) %>% 
  mutate(genus_rmvspace_1 = str_replace_all(genus_rmvspace, " ", "_")) %>% 
  mutate(genus_rmvspace_1 = factor(genus_rmvspace_1, levels = tree_fasttree$tip.label))%>% 
  arrange(genus_rmvspace_1)

phylo_data <- a %>% 
  mutate(label = genus_rmvspace_1) %>% 
  select(label, Genus, Phylum, count_found_loc, locality, genus_rmvspace_1) %>% 
  unique() 

df_prev_ThaiDataset <- metadata_prevalence %>% 
  select(Phylum, Genus, prev_all) %>% 
  mutate(locality = "Thai datasets") %>% 
  rename(Prevalence = prev_all)
 
df_prev_thisStudy <- phylo_data %>% 
  select(locality, Phylum, Genus, count_found_loc) %>% 
  rename(Prevalence = count_found_loc) 

df_prev_comb <- rbind(df_prev_thisStudy, df_prev_ThaiDataset) %>% 
  mutate(locality = case_when(
    locality == "bangkok" ~ "Bangkok",
    locality == "phatthalung" ~ "Phatthalung",
    locality == "tak" ~ "Tak",
    .default = locality
  ))

write.csv(x = df_prev_comb, 
          "C:/Project/5_16s_thai_population/revise_aftermSystems/21092026_prevalence_phylogenetic.csv", row.names = FALSE)

# 1. Get the unique 15 phyla in your data
phylum_list <- unique(phylo_data$Phylum)

# 2. Manually assign 15 colors 
my_colors <- c(
  "#DA362A", "#548EC0", "#ACCEE3", "#E49048", "#6EA7A0",
  "#4CA330", "#F0AA63", "#E6D27A", "#C1AAD2", "#7B55A5",
  "#B49B74", "#9ED176", "#E66968", "#B6A499", "#F4992D"
)
my_colors <- c("Actinobacteriota"  = "#ACCEE3",
               "Bacteroidota"      = "#548EC0",
               "Desulfobacterota"  = "#4CA330",
               "Euryarchaeota"     = "#E66968",
               "Firmicutes"        = "#DA362A",
               "Fusobacteriota"    = "#F0AA63",
               "Proteobacteria"    = "#E49048",
               "Verrucomicrobiota" = "#E6D27A",
               "Campylobacterota"  = "#6EA7A0",
               "Spirochaetota"     = "#C1AAD2",
               "Synergistota"      = "#7B55A5",
               "Elusimicrobiota"   = "#B49B74",
               "Thermoplasmatota"  = "#B6A499",
               "Halobacterota"     = "#F4992D",
               "Cyanobacteria"     = "#9ED176")

# 3. Name the color vector with Phylum names
phylum_colors <- setNames(my_colors, phylum_list)


# Apply to each dataset
bkk_PA <- phylo_data %>%
  filter(locality == "bangkok") %>% 
  mutate(count_bin = cut(count_found_loc,
                         breaks = c(-Inf, 0, 0.2, 0.4, 0.6, 0.8, Inf),
                         labels = c("0", "0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1"),
                         right = TRUE)) %>% 
  select(-locality, -Phylum, -Genus, -label, -count_found_loc)

phat_PA <- phylo_data %>%
  filter(locality == "phatthalung") %>% 
  mutate(count_bin = cut(count_found_loc,
                         breaks = c(-Inf, 0, 0.2, 0.4, 0.6, 0.8, Inf),
                         labels = c("0", "0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1"),
                         right = TRUE))%>% 
  select(-locality, -Phylum, -Genus, -label, -count_found_loc)

tak_PA <- phylo_data %>%
  filter(locality == "tak") %>% 
  mutate(count_bin = cut(count_found_loc,
                         breaks = c(-Inf, 0, 0.2, 0.4, 0.6, 0.8, Inf),
                         labels = c("0", "0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1"),
                         right = TRUE))%>% 
  select(-locality, -Phylum, -Genus, -label, -count_found_loc)

all_PA <- metadata_prevalence %>% 
  mutate(count_bin = cut(prev_all,
                        breaks = c(-Inf, 0, 0.2, 0.4, 0.6, 0.8, Inf),
                        labels = c("0", "0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1"),
                        right = TRUE))%>% 
  select(genus_rmvspace_1, count_bin)

count_colors <- c(
  "0" = "white",       # light gray
  "0-0.2" = "#a4f4a1",   # yellow
  "0.2-0.4" = "#7fc194", # orange
  "0.4-0.6" = "#5a8e87", # red
  "0.6-0.8" = "#355a7a", # purple
  "0.8-1" = "#10276e"    # blue
)


p_tree <-ggtree(tree_fasttree, layout = "circular", branch.length = "none") %<+% phylo_data +
  geom_tippoint(size = 0.8, aes(colour = Phylum)) +
  scale_colour_manual(values = my_colors)+
  theme_tree() +
  geom_fruit(data = bkk_PA, geom = geom_tile, mapping = aes(y = genus_rmvspace_1, fill = count_bin), width = 0.6, color = "black")+
  geom_fruit(data = phat_PA, geom = geom_tile, mapping = aes(y = genus_rmvspace_1, fill = count_bin), width = 0.6, color = "black")+
  geom_fruit(data = tak_PA, geom = geom_tile, mapping = aes(y = genus_rmvspace_1, fill = count_bin), width = 0.6, color = "black")+
  geom_fruit(data = all_PA, geom = geom_tile, mapping = aes(y = genus_rmvspace_1, fill = count_bin), width = 0.6, color = "black")+
  scale_fill_manual(name = "Prevalence fraction", values = count_colors) +
  theme(legend.position = "right")
  
ggsave(p_tree, file="C:/Project/5_16s_thai_population/figure/revision_figures/07092026_tree_addAlldata_fixlabel.png", dpi = 600, width = 7, height = 6, units = "in")

p_tree_nl <-ggtree(tree_fasttree, layout = "circular", branch.length = "none") %<+% phylo_data +
  geom_tippoint(size = 0.7, aes(colour = Phylum)) +
  scale_colour_manual(values = my_colors)+
  theme_tree() +
  geom_fruit(data = bkk_PA, geom = geom_tile, mapping = aes(y = genus_rmvspace_1, fill = count_bin), width = 0.6, color = "black")+
  geom_fruit(data = phat_PA, geom = geom_tile, mapping = aes(y = genus_rmvspace_1, fill = count_bin), width = 0.6, color = "black")+
  geom_fruit(data = tak_PA, geom = geom_tile, mapping = aes(y = genus_rmvspace_1, fill = count_bin), width = 0.6, color = "black")+
  geom_fruit(data = all_PA, geom = geom_tile, mapping = aes(y = genus_rmvspace_1, fill = count_bin), width = 0.6, color = "black")+
  scale_fill_manual(name = "Prevalence fraction", values = count_colors) +
  theme(legend.position = "none")

ggsave(p_tree_nl, file="C:/Project/5_16s_thai_population/figure/revision_figures/07092026_tree_nolegend_addAlldata_fixlabel.png",dpi = 600, width = 5, height = 4, units = "in")


# All prevalence that exclude this study

metadata_prevalence_excludeThisStudy <- taxaTable_combine %>% 
  filter(!Sample %in% remove_sample, 
         Study != "This study") %>% 
  filter(Abundance != 0 & !is.na(Genus)) %>%
  select(Sample, Genus, Phylum) %>% 
  unique() %>% 
  group_by(Phylum, Genus) %>% 
  summarise(count_genus = n(), .groups = "drop") %>% 
  mutate(prev_all = count_genus/340) %>% 
  mutate(genus_rmvspace = str_remove_all(Genus, "\\[.*?\\]")) %>% 
  mutate(genus_rmvspace_1 = str_replace_all(genus_rmvspace, " ", "_")) %>% 
  mutate(genus_rmvspace_1 = factor(genus_rmvspace_1, levels = tree_fasttree$tip.label))%>% 
  arrange(genus_rmvspace_1)

all_PA_excludeThisStudy <- metadata_prevalence_excludeThisStudy %>% 
  mutate(count_bin = cut(prev_all,
                         breaks = c(-Inf, 0, 0.2, 0.4, 0.6, 0.8, Inf),
                         labels = c("0", "0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1"),
                         right = TRUE))%>% 
  select(genus_rmvspace_1, count_bin)

p_tree_nl_excludeThisStudy <-ggtree(tree_fasttree, layout = "circular", branch.length = "none") %<+% phylo_data +
  geom_tippoint(size = 1.1, aes(colour = Phylum)) +
  scale_colour_manual(values = my_colors)+
  theme_tree() +
  geom_fruit(data = bkk_PA, geom = geom_tile, mapping = aes(y = genus_rmvspace_1, fill = count_bin), width = 0.6, color = "black")+
  geom_fruit(data = phat_PA, geom = geom_tile, mapping = aes(y = genus_rmvspace_1, fill = count_bin), width = 0.6, color = "black")+
  geom_fruit(data = tak_PA, geom = geom_tile, mapping = aes(y = genus_rmvspace_1, fill = count_bin), width = 0.6, color = "black")+
  geom_fruit(data = all_PA_excludeThisStudy, geom = geom_tile, mapping = aes(y = genus_rmvspace_1, fill = count_bin), width = 0.6, color = "black")+
  scale_fill_manual(name = "Prevalence fraction", values = count_colors) +
  theme(legend.position = "none")

ggsave(p_tree_nl_excludeThisStudy, file="C:/Project/5_16s_thai_population/figure/revision_figures/07092026_tree_nolegend_addAlldata_excludeThisStudy_fixLabel.png",
       dpi = 600, width = 7, height = 6, units = "in")



# ---------------------------------------------------------------------------- #
prev_thisStudy <- phylo_data %>% select(genus_rmvspace_1, count_found_loc, locality)
prev_thaiDataset <- metadata_prevalence_excludeThisStudy %>% rename(count_found_loc = prev_all) %>% 
  filter(genus_rmvspace_1 %in% unique(prev_thisStudy$genus_rmvspace_1)) %>% 
  mutate(locality = "Thai datasets") %>% 
  select(genus_rmvspace_1, count_found_loc, locality)

combin_prev_all <- rbind(prev_thisStudy, prev_thaiDataset)

mat_prev_all <- combin_prev_all %>%
  pivot_wider(
    names_from = locality, 
    values_from = count_found_loc
  ) %>% 
  unique() %>% 
  column_to_rownames(var = "genus_rmvspace_1") %>% 
  as.matrix()

cor_loc <- cor(mat_prev_all, use = "pairwise.complete.obs", method = "pearson")
cor_loc[upper.tri(cor_loc, diag = FALSE)] <- NA
cor_loc <- cor_loc %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Locality_1") %>% 
  pivot_longer(cols = -Locality_1, names_to = "Locality_2", values_to = "PCC") %>% 
  mutate(PCC = if_else(Locality_2 == "tak" & Locality_1 == "phatthalung", NA,
                       if_else(Locality_1 == "tak" & Locality_2 == "phatthalung", 0.88, PCC)))

p_pearson_prev <- ggplot(cor_loc, aes(x = Locality_1, y = Locality_2, fill = PCC)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(PCC, 2)), color = "black", size = 2.5) +
  scale_fill_gradient2(
    low = "dodgerblue3", 
    mid = "white", 
    high = "#CD2626", 
    midpoint = 0, 
    limit = c(-1,1), na.value = "lightgrey"
  ) +
  coord_fixed() +
  theme_bw(base_size = 9) +
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))+
  scale_y_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))+
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 310, hjust = 0),
    panel.grid = element_blank(),
    legend.position = "none",
    #legend.title = element_text(size = 7), 
    #legend.text = element_text(size = 6),
    #legend.key.size = unit(0.3, "cm"),
    plot.margin = margin(t = 10, r = 25, b = 0, l = 10, unit = "pt"))

p_comb_prevalence <- plot_grid(p_pearson_prev, p_appearance, p_prevalance, p_app_pre, ncol = 1, rel_heights = c(1, 0.85,1,0.85))
ggsave(plot = p_comb_prevalence,
       filename = "C:/Project/5_16s_thai_population/figure/revision_figures/08092026_p_comb_prevalence_limitThaiAllGenus.png", 
       width = 2.1, height = 7, units = "in", dpi = 600)




mat <- as.matrix(mat_prev_all)
n <- ncol(mat)
vars <- colnames(mat)

# 1. Initialize an empty list to collect pairwise results
results <- list()
k <- 1

# 2. Loop through unique pairs (upper triangle only)
for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    # Compute correlation test with pairwise handling
    test <- cor.test(mat[, i], mat[, j], use = "pairwise.complete.obs", method = "pearson")
    
    results[[k]] <- data.frame(
      Var1 = vars[i],
      Var2 = vars[j],
      r = test$estimate,
      p_value = test$p.value,
      stringsAsFactors = FALSE
    )
    k <- k + 1
  }
}

# 3. Combine list into a single data frame
df_cor <- do.call(rbind, results)

# Optional: Add significance stars or FDR-adjusted p-values
df_cor$p_adj <- p.adjust(df_cor$p_value, method = "fdr")
df_cor$sig_star <- ifelse(df_cor$p_value < 0.001, "***",
                          ifelse(df_cor$p_value < 0.01,  "**",
                                 ifelse(df_cor$p_value < 0.05,  "*", "ns")))

# View output
as.data.frame(df_cor)


write.csv(x = as.data.frame(df_cor), file = "C:/Project/5_16s_thai_population/revise_aftermSystems/21092026_Prevalence_correlation.csv", row.names = FALSE)
