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

## Rarefying ----

# Rarefaction curve plot all ----
spnumber_all <- specnumber(t(asv_table))
Srare_all <- rarefy(t(asv_table), sample = min(rowSums(t(asv_table))))

png("C:/Project/5_16s_thai_population/figure/12112025_all_species_rarefaction.png", width = 3.5, height = 3.5, res = 120, units = "in")
plot(spnumber_all, Srare_all,
     xlab = "Observed No. of Species",
     ylab = "Rarefied No. of Species",
     col = "black", pch = 19, 
     cex.lab = 0.7,     # axis label size
     cex.axis = 0.7,)
abline(0, 1, col = "black", lwd = 1.5)
dev.off()

rarecurve_all <- rarecurve(t(asv_table),
                            sample = min(rowSums(t(asv_table))),
                            step = 10,
                            tidy = TRUE)

p_all_rarecurv <- ggplot(rarecurve_all, aes(x = Sample, y = Species))+
  geom_point(colour = "black", size = 0.4, alpha = 0.7)+
  geom_vline(xintercept = min(rowSums(t(asv_table))))+
  theme_bw()+
  labs(x = "Sequencing read", y = "ASV")+
  ylim(0, 300)+xlim(0, 62000)

ggsave(plot = p_all_rarecurv,
       filename = "C:/Project/5_16s_thai_population/figure/12112025_p_all_rarecurv.png", 
       width = 3.5, height = 3.5, units = "in")

# Rarefaction curve plots for each location ----
sampleID_bkk <- meta_table %>% filter(locality == "bangkok") %>% rownames()
sampleID_tak <- meta_table %>% filter(locality == "tak") %>% rownames()
sampleID_phat <- meta_table %>% filter(locality == "phatthalung") %>% rownames()

matrix_bkk <- t(asv_table[, sampleID_bkk])
spnumber_bkk <- specnumber(matrix_bkk)
Srare_bkk <- rarefy(matrix_bkk, sample = min(rowSums(matrix_bkk)))

matrix_tak <- t(asv_table[, sampleID_tak])
spnumber_tak <- specnumber(matrix_tak)
Srare_tak <- rarefy(matrix_tak, sample = min(rowSums(matrix_tak)))

matrix_phat <- t(asv_table[, sampleID_phat])
spnumber_phat <- specnumber(matrix_phat)
Srare_phat <- rarefy(matrix_phat, sample = min(rowSums(matrix_phat)))

png("C:/Project/5_16s_thai_population/figure/12112025_bkk_species_rarefaction.png", width = 3.5, height = 3.5, res = 120, units = "in")
plot(spnumber_bkk, Srare_bkk,
     xlab = "Observed No. of Species",
     ylab = "Rarefied No. of Species",
     col = "#F8766D", pch = 19, 
     cex.lab = 0.7,     # axis label size
     cex.axis = 0.7,)
abline(0, 1, col = "black", lwd = 1.5)
dev.off()

png("C:/Project/5_16s_thai_population/figure/12112025_tak_species_rarefaction.png", width = 3.5, height = 3.5, res = 120, units = "in")
plot(spnumber_tak, Srare_tak,
     xlab = "Observed No. of Species",
     ylab = "Rarefied No. of Species",
     col = "#619CFF", pch = 19, 
     cex.lab = 0.7,     # axis label size
     cex.axis = 0.7,)
abline(0, 1, col = "black", lwd = 1.5)
dev.off()

png("C:/Project/5_16s_thai_population/figure/12112025_phat_species_rarefaction.png", width = 3.5, height = 3.5, res = 120, units = "in")
plot(spnumber_phat, Srare_phat,
     xlab = "Observed No. of Species",
     ylab = "Rarefied No. of Species",
     col = "#00BA38", pch = 19, 
     cex.lab = 0.7,     # axis label size
     cex.axis = 0.7,)
abline(0, 1, col = "black", lwd = 1.5)
dev.off()

rarecurve_bkk <- rarecurve(matrix_bkk,
                           sample = min(rowSums(matrix_bkk)),
                           step = 10,
                           tidy = TRUE)

rarecurve_tak <- rarecurve(matrix_tak,
                           sample = min(rowSums(matrix_tak)),
                           step = 10,
                           tidy = TRUE)

rarecurve_phat <- rarecurve(matrix_phat,
                           sample = min(rowSums(matrix_phat)),
                           step = 10,
                           tidy = TRUE)

p_bkk_rarecurv <- ggplot(rarecurve_bkk, aes(x = Sample, y = Species))+
  geom_point(colour = "#F8766D", size = 0.4, alpha = 0.7)+
  geom_vline(xintercept = min(rowSums(matrix_bkk)))+
  theme_bw()+
  labs(x = "Sequencing read", y = "ASV")+
  ylim(0, 300)+xlim(0, 62000)
  
p_tak_rarecurv <- ggplot(rarecurve_tak, aes(x = Sample, y = Species))+
  geom_point(colour = "#619CFF", size = 0.4, alpha = 0.7)+
  geom_vline(xintercept = min(rowSums(matrix_tak)))+
  theme_bw()+
  labs(x = "Sequencing read", y = "ASV")+
  ylim(0, 300)+xlim(0, 62000)

p_phat_rarecurv <- ggplot(rarecurve_phat, aes(x = Sample, y = Species))+
  geom_point(colour = "#00BA38", size = 0.4, alpha = 0.7)+
  geom_vline(xintercept = min(rowSums(matrix_phat)))+
  theme_bw()+
  labs(x = "Sequencing read", y = "ASV")+
  ylim(0, 300)+xlim(0, 62000)

ggsave(plot = p_bkk_rarecurv,
       filename = "C:/Project/5_16s_thai_population/figure/12112025_p_bkk_rarecurv.png", 
       width = 3.5, height = 3.5, units = "in")

ggsave(plot = p_phat_rarecurv,
       filename = "C:/Project/5_16s_thai_population/figure/12112025_p_phat_rarecurv.png", 
       width = 3.5, height = 3.5, units = "in")

ggsave(plot = p_tak_rarecurv,
       filename = "C:/Project/5_16s_thai_population/figure/12112025_p_tak_rarecurv.png", 
       width = 3.5, height = 3.5, units = "in")



# Rarefying samples to min sequencing dept ----
set.seed(12112025)
asv_table_rarefied <- rrarefy(t(asv_table), sample = min(rowSums(t(asv_table))))
#rowSums(asv_table_rarefied)
#col_sums=colSums(asv_table_rarefied)
#col_sums[col_sums == 0]

## Phyloseq object ----

# Convert to phyloseq components
ASV <- otu_table(t(asv_table_rarefied), taxa_are_rows = TRUE)
TAX <- tax_table(taxa_table)
SAMP <- sample_data(meta_table)

physeq <- phyloseq(ASV, TAX, SAMP)

#sample_data(physeq)

taxonomic_df_all <- psmelt(physeq)


## Plots of prevalence analysis ----

# Count appearance (1, 2, or all location) of microbial genera 
p_appearance <- taxonomic_df_all %>% 
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
       y = "Number of genus")

# Prevalence of each microbial genus in each location
p_prevalance <- taxonomic_df_all %>% 
  filter(Abundance != 0 & !is.na(Genus)) %>%
  select(Sample, Genus, locality) %>% 
  unique() %>% 
  group_by(Genus, locality) %>% 
  summarise(count_genus = n()) %>% 
  ungroup() %>% 
  mutate(count_found_loc = if_else(locality == "bangkok", count_genus/47,
                                   if_else(locality == "phatthalung", count_genus/29, count_genus/30))) %>% 
  ggplot(aes(x = locality, y = count_found_loc, fill = locality))+
  geom_boxplot()+
  geom_jitter(width = 0.25, alpha=0.4)+
  theme_bw()+
  labs(x = "",
       y = "Prevalence")+
  theme(legend.position = "none", axis.text.x = element_text(angle = 315, hjust = 0))+
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))

# Plot prevalence against and appearance separated by locations
p_app_pre <- taxonomic_df_all %>% 
  filter(Abundance != 0 & !is.na(Genus)) %>%
  select(Sample, Genus, locality) %>% 
  unique() %>% 
  group_by(Genus, locality) %>% 
  summarise(count_genus = n()) %>% 
  ungroup() %>% 
  mutate(count_found_loc = if_else(locality == "bangkok", count_genus/47,
                                   if_else(locality == "phatthalung", count_genus/29, count_genus/30))) %>%
  group_by(Genus) %>% 
  mutate(appear = n()) %>% 
  ggplot(aes(x = as.factor(appear), y = count_found_loc, fill = locality))+
  geom_boxplot()+
  theme_bw()+
  labs(x = "Appearance",
       y = "Prevalence")+
  theme(legend.position = "none")

p_comb_prevalence <- plot_grid(p_appearance, p_prevalance, p_app_pre, ncol = 1)

ggsave(plot = p_comb_prevalence,
       filename = "C:/Project/5_16s_thai_population/figure/13112025_p_comb_prevalence_rarefied.png", 
       width = 2, height = 7, units = "in")


## Get representative ASV for a each genus and extract sequence for later building phylogenetic tree ----
taxa_table_df <- taxa_table %>% as.data.frame() %>% rownames_to_column(var = "ASV") %>% glimpse()

# There are multiple ASVs in a genus due to limitation of 16s rDNA sequencing technique 
# Across all samples, we selected the representative AVS that has the highest sum relative abundance

asv_table_rarefied_extract <- t(asv_table_rarefied) %>% 
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
  #slice_max(order_by = sumRA, n = 1, with_ties = FALSE) %>% 
  mutate(max_sumRA = if_else(sumRA == max(sumRA), 1, 0)) %>% 
  filter(max_sumRA == 1) %>% 
  ungroup() %>% 
  filter(sumRA != 0 & !is.na(Genus)) %>% 
  select(Genus, ASV, sumRA) %>% unique() 

# There are duplicated genera when using rarefied samples in this analysis
# Here, I just simply selected the first sequence that show up

asv_table_rarefied_extract <- asv_table_rarefied_extract[!duplicated(asv_table_rarefied_extract$Genus), ]

# Extract the sequences
library(Biostrings)
dna <- DNAStringSet(asv_table_rarefied_extract$ASV)
names(dna) <- gsub(" ", "_", asv_table_rarefied_extract$Genus) 
names(dna)

writeXStringSet(dna, "C:/Project/5_16s_thai_population/13112025_asv_genus_sequences_short_filter0_rarefied.fasta")


## Phylogenetic tree plot ----

library(ggtree)
library(ggtreeExtra)
library(stringr)

# Phylogenetic tree was generated in Bash.
# Multiple sequence alignment was conducted using MUSCLE v5.1 and Gblocks v0.91b
# was used to select high conserve aligned blocks.
# Finally, FastTree v2.2 was used to construct the phylogenetic tree.
# The script for this was attached separately.

tree_fasttree <- read.tree("C:/Project/5_16s_thai_population/Phylogenetic tree/13112025_asv_genus_sequences_short_rarefied_tree")
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

phylo_data <- a %>% mutate(label = genus_rmvspace_1) %>% select(label, Genus, Phylum, count_found_loc, locality, genus_rmvspace_1) %>% unique() 
#view(phylo_data)
# 1. Get the unique 15 phyla in your data
phylum_list <- unique(phylo_data$Phylum)

# 2. Manually assign 14 colors 
# Since order change from previous tree, I brought oerder below, and manually re-order 

#Firmicutes      Bacteroidota  Actinobacteriota    Proteobacteria  Campylobacterota  Desulfobacterota    Fusobacteriota Verrucomicrobiota     Spirochaetota 
#"#DA362A"         "#548EC0"         "#ACCEE3"         "#E49048"         "#6EA7A0"         "#4CA330"         "#F0AA63"         "#E6D27A"         "#C1AAD2" 
#Synergistota   Elusimicrobiota     Euryarchaeota  Thermoplasmatota     Halobacterota     Cyanobacteria 
#"#7B55A5"         "#B49B74"         "#9ED176"         "#E66968"         "#B6A499"         "#F4992D" 

my_colors <- c(
  "#DA362A", "#7B55A5", "#E6D27A", "#B49B74", "#ACCEE3",
  "#C1AAD2", "#F0AA63", "#E49048", "#4CA330", "#6EA7A0", 
  "#9ED176", "#E66968", "#B6A499", "#548EC0"
)

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


count_colors <- c(
  "0" = "white",       # light gray
  "0-0.2" = "#a4f4a1",   # yellow
  "0.2-0.4" = "#7fc194", # orange
  "0.4-0.6" = "#5a8e87", # red
  "0.6-0.8" = "#355a7a", # purple
  "0.8-1" = "#10276e"    # blue
)

# 

p_tree <-ggtree(tree_fasttree, layout = "circular", branch.length = "none") %<+% phylo_data +
  geom_tippoint(size = 1.2, aes(colour = Phylum)) +
  scale_colour_manual(values = phylum_colors)+
  theme_tree() +
  geom_fruit(data = bkk_PA, geom = geom_tile, mapping = aes(y = genus_rmvspace_1, fill = count_bin), width = 0.6, color = "black")+
  geom_fruit(data = phat_PA, geom = geom_tile, mapping = aes(y = genus_rmvspace_1, fill = count_bin), width = 0.6, color = "black")+
  geom_fruit(data = tak_PA, geom = geom_tile, mapping = aes(y = genus_rmvspace_1, fill = count_bin), width = 0.6, color = "black")+
  scale_fill_manual(name = "Prevalence fraction", values = count_colors) +
  theme(legend.position = "right")

ggsave(p_tree, file="C:/Project/5_16s_thai_population/figure/13112025_tree_rarefied.png", dpi = 600, width = 7, height = 6, units = "in")

p_tree_nl <-ggtree(tree_fasttree, layout = "circular", branch.length = "none") %<+% phylo_data +
  geom_tippoint(size = 0.5, aes(colour = Phylum)) +
  scale_colour_manual(values = phylum_colors)+
  theme_tree() +
  geom_fruit(data = bkk_PA, geom = geom_tile, mapping = aes(y = genus_rmvspace_1, fill = count_bin), width = 0.6, color = "black")+
  geom_fruit(data = phat_PA, geom = geom_tile, mapping = aes(y = genus_rmvspace_1, fill = count_bin), width = 0.6, color = "black")+
  geom_fruit(data = tak_PA, geom = geom_tile, mapping = aes(y = genus_rmvspace_1, fill = count_bin), width = 0.6, color = "black")+
  scale_fill_manual(name = "Prevalence fraction", values = count_colors) +
  theme(legend.position = "none")

ggsave(p_tree_nl, file="C:/Project/5_16s_thai_population/figure/13112025_tree_rarefied_nolegend.png", 
       dpi = 600, width = 5, height = 4, units = "in")
