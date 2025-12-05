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
       filename = "C:/Project/5_16s_thai_population/figure/01092025_p_comb_prevalence.png", 
       width = 2, height = 7, units = "in")

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

phylo_data <- a %>% mutate(label = genus_rmvspace_1) %>% select(label, Genus, Phylum, count_found_loc, locality, genus_rmvspace_1) %>% unique() 

# 1. Get the unique 15 phyla in your data
phylum_list <- unique(phylo_data$Phylum)

# 2. Manually assign 15 colors 
my_colors <- c(
  "#DA362A", "#548EC0", "#ACCEE3", "#E49048", "#6EA7A0",
  "#4CA330", "#F0AA63", "#E6D27A", "#C1AAD2", "#7B55A5",
  "#B49B74", "#9ED176", "#E66968", "#B6A499", "#F4992D"
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


p_tree <-ggtree(tree_fasttree, layout = "circular", branch.length = "none") %<+% phylo_data +
  geom_tippoint(size = 1.2, aes(colour = Phylum)) +
  scale_colour_manual(values = phylum_colors)+
  theme_tree() +
  geom_fruit(data = bkk_PA, geom = geom_tile, mapping = aes(y = genus_rmvspace_1, fill = count_bin), width = 0.6, color = "black")+
  geom_fruit(data = phat_PA, geom = geom_tile, mapping = aes(y = genus_rmvspace_1, fill = count_bin), width = 0.6, color = "black")+
  geom_fruit(data = tak_PA, geom = geom_tile, mapping = aes(y = genus_rmvspace_1, fill = count_bin), width = 0.6, color = "black")+
  scale_fill_manual(name = "Prevalence fraction", values = count_colors) +
  theme(legend.position = "right")
  
#ggsave(p_tree, file="C:/Project/5_16s_thai_population/figure/01092025_tree.png", dpi = 600, width = 7, height = 6, units = "in")

p_tree_nl <-ggtree(tree_fasttree, layout = "circular", branch.length = "none") %<+% phylo_data +
  geom_tippoint(size = 0.9, aes(colour = Phylum)) +
  scale_colour_manual(values = phylum_colors)+
  theme_tree() +
  geom_fruit(data = bkk_PA, geom = geom_tile, mapping = aes(y = genus_rmvspace_1, fill = count_bin), width = 0.6, color = "black")+
  geom_fruit(data = phat_PA, geom = geom_tile, mapping = aes(y = genus_rmvspace_1, fill = count_bin), width = 0.6, color = "black")+
  geom_fruit(data = tak_PA, geom = geom_tile, mapping = aes(y = genus_rmvspace_1, fill = count_bin), width = 0.6, color = "black")+
  scale_fill_manual(name = "Prevalence fraction", values = count_colors) +
  theme(legend.position = "none")

#ggsave(p_tree_nl, file="C:/Project/5_16s_thai_population/figure/01092025_tree_nolegend.png",dpi = 600, width = 5, height = 4, units = "in")

