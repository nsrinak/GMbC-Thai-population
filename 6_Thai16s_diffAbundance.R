library(tidyverse)
library(ggplot2)
library(cowplot)
library(phyloseq)
library(ANCOMBC)

# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #

asv_table <- readRDS("C:/Project/5_16s_thai_population/seqtab_final.rds")
asv_table <- t(asv_table)

taxa_table <- readRDS("C:/Project/5_16s_thai_population/tax_final.rds")

taxa_table <- taxa_table %>% 
  as.data.frame() %>% 
  mutate(
    # Resolve missing Family → Order
    Family_resolved = coalesce(Family, Order),
    
    # Resolve missing Genus → Family
    Genus_resolved = coalesce(Genus, Family_resolved),
    
    # Build taxon label
    Taxon_label = case_when(
      !is.na(Species) ~ paste0("s__", Genus_resolved, "_", Species),
      !is.na(Genus)   ~ paste0("g__", Genus_resolved),
      !is.na(Family)  ~ paste0("f__", Family_resolved),
      !is.na(Order)   ~ paste0("o__", Order),
      TRUE            ~ "unclassified"
    ),
    
    # Make names unique
    Taxon_label = make.unique(Taxon_label, sep = "_")
  ) %>% 
  as.matrix()

meta_table <- read_tsv("C:/Project/5_16s_thai_population/thai_16s/metaData/metadata.tsv")
meta_table <- meta_table %>% column_to_rownames(var = "donor_id") 

## Phyloseq object ----

# Convert to phyloseq components
ASV <- otu_table(asv_table, taxa_are_rows = TRUE)
TAX <- tax_table(taxa_table)
SAMP <- sample_data(meta_table)

physeq <- phyloseq(ASV, TAX, SAMP)

# Define formula for fixed effects and group
formula = "locality + age + sex + bmi"               # Fixed effects
group = "locality"                             # Main grouping variable
random = "Subject"                             # Random effect for repeated measures (optional)

# Run the ANCOM-BC function
set.seed(123)
ancombc_genus <- ancombc2(
  data = physeq,              
  tax_level = "Genus",
  fix_formula = "locality + age + sex + bmi",
  rand_formula = NULL,
  prv_cut = 0.1, 
  lib_cut = 1000, 
  s0_perc = 0.05,
  group = "locality", 
  p_adj_method = "BH",
  pseudo_sens = TRUE,
  struc_zero = TRUE, 
  neg_lb = TRUE,
  alpha = 0.05, 
  n_cl = 8, 
  verbose = TRUE,
  global = FALSE, 
  pairwise = TRUE, 
  dunnet = FALSE, 
  trend = FALSE,
  iter_control = list(tol = 1e-2, max_iter = 20, 
                      verbose = TRUE),
  em_control = list(tol = 1e-5, max_iter = 100),
  lme_control = lme4::lmerControl(),
  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100))

set.seed(123)
ancombc_species <- ancombc2(
  data = physeq,              
  tax_level = NULL,
  fix_formula = "locality + age + sex + bmi",
  rand_formula = NULL,
  prv_cut = 0.1, 
  lib_cut = 1000, 
  s0_perc = 0.05,
  group = "locality", 
  p_adj_method = "BH",
  pseudo_sens = TRUE,
  struc_zero = TRUE, 
  neg_lb = TRUE,
  alpha = 0.05, 
  n_cl = 8, 
  verbose = TRUE,
  global = FALSE, 
  pairwise = TRUE, 
  dunnet = FALSE, 
  trend = FALSE,
  iter_control = list(tol = 1e-2, max_iter = 20, 
                      verbose = TRUE),
  em_control = list(tol = 1e-5, max_iter = 100),
  lme_control = lme4::lmerControl(),
  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100))

# --------------------------------------------------------------------------- #
ancombc_genus$res_pair %>% glimpse()

ancombc_genus$res_pair %>% 
  filter(q_localityphatthalung < 0.05|q_localitytak < 0.05|q_localitytak_localityphatthalung < 0.05) %>% 
  pull(taxon) %>% unique()

ancombc_species$res_pair %>% glimpse()

ancombc_species$res_pair %>% 
  filter(q_localityphatthalung < 0.05|q_localitytak < 0.05|q_localitytak_localityphatthalung < 0.05) %>% 
  pull(taxon) %>% unique()

asvToLabel <-taxa_table %>% 
  as.data.frame %>% 
  rownames_to_column(var = "taxon") %>%
  select(taxon, Taxon_label) %>% 
  unique()

ancombc_species$res_pair <- ancombc_species$res_pair %>% 
  left_join(., asvToLabel, by = "taxon") %>% 
  glimpse()

write.csv(x = ancombc_species$res_pair, file = "C:/Project/5_16s_thai_population/revise_aftermSystems/21092026_ANCOMBC2_ASV.csv", row.names = FALSE)
write.csv(x = ancombc_genus$res_pair, file = "C:/Project/5_16s_thai_population/revise_aftermSystems/21092026_ANCOMBC2_genus.csv", row.names = FALSE)
# --------------------------------------------------------------------------- #

ancombc_res <- taxa_table %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "asv") %>%
  left_join(., ancombc_genus$res_pair,   by = c("Genus" = "taxon")) %>% 
  left_join(., ancombc_species$res_pair, by = c("asv"   = "taxon")) %>% 
  mutate(lfc_localityphatthalung.x             = if_else(q_localityphatthalung.x             < 0.05, lfc_localityphatthalung.x, 0),
         lfc_localitytak.x                     = if_else(q_localitytak.x                     < 0.05, lfc_localitytak.x, 0),
         lfc_localitytak_localityphatthalung.x = if_else(q_localitytak_localityphatthalung.x < 0.05, lfc_localitytak_localityphatthalung.x, 0),
         lfc_localityphatthalung.y             = if_else(q_localityphatthalung.y             < 0.05, lfc_localityphatthalung.y, 0),
         lfc_localitytak.y                     = if_else(q_localitytak.y                     < 0.05, lfc_localitytak.y, 0),
         lfc_localitytak_localityphatthalung.y = if_else(q_localitytak_localityphatthalung.y < 0.05, lfc_localitytak_localityphatthalung.y, 0)) %>%
  select(Phylum, Genus, Taxon_label, 
         lfc_localityphatthalung.x, lfc_localitytak.x, lfc_localitytak_localityphatthalung.x,
         lfc_localityphatthalung.y, lfc_localitytak.y, lfc_localitytak_localityphatthalung.y) %>% 
  pivot_longer(cols = -c(Phylum, Genus, Taxon_label), names_to = "name", values_to = "value") %>% 
  mutate(level = if_else(grepl("\\.x$", name), "Genus", "Species")) %>% 
  filter(value != 0, 
         !is.na(value))

# Prepare data and explicit factor levels
plot_data <- ancombc_res %>% 
  arrange(Genus, Taxon_label) %>% 
  mutate(
    Genus_unique = paste(Genus, name, sep = "."),
    Genus_unique = make.unique(Genus_unique, sep = "_x_"),
    value = if_else(grepl("_x_", Genus_unique) & level == "Genus", NA, value),
    Genus = if_else(is.na(Genus), "Unknown", Genus),
    Genus = sub("\\[", replacement = "", x = Genus),
    Genus = sub("\\]", replacement = "", x = Genus)) %>% 
  filter(!is.na(value)) %>% 
  arrange(Phylum, Genus, Taxon_label) %>% 
  mutate(number_phylum = as.numeric(as.factor(Phylum)),
         number_genus = as.numeric(as.factor(Genus)),
         number_taxon = as.numeric(as.factor(Taxon_label)),
         number_genus = if_else(Genus == "Unknown", number_genus * 20, number_genus)) %>% glimpse()

# Prepare plot_data: Sort and set explicit factor levels
plot_data_prepared <- plot_data %>%
  mutate(
    # Compute sorting rank
    sort_score = ((number_phylum * 100) + number_genus + (number_taxon / 1e6)) * -1,
    # Convert Taxon_label to a factor ordered by sort_score
    Taxon_label = fct_reorder(Taxon_label, sort_score)
  )

# Calculate horizontal line positions (now safe from NAs)
hline_positions <- plot_data_prepared %>%
  select(Taxon_label, Genus) %>%
  distinct() %>%
  mutate(y_pos = as.numeric(Taxon_label)) %>% # Returns 1, 2, 3... corresponding to y-axis index
  group_by(Genus) %>%
  summarise(max_pos = max(y_pos), .groups = "drop") %>%
  filter(max_pos < max(as.numeric(plot_data_prepared$Taxon_label))) %>%
  transmute(intercept = max_pos + 0.5)

# Plot with horizontal separator lines
# Define custom colors for each facet column (matching the unique values in 'name')
facet_colors <- c(
  "lfc_localityphatthalung.x"            = "palevioletred1",
  "lfc_localityphatthalung.y"            = "palevioletred1",
  "lfc_localitytak_localityphatthalung.x"= "lightsteelblue3",
  "lfc_localitytak_localityphatthalung.y"= "lightsteelblue3",
  "lfc_localitytak.x"                    = "paleturquoise3",
  "lfc_localitytak.y"                    = "paleturquoise3"
  # Add or adjust key-value pairs matching your 'name' column levels
)

p_main <- plot_data %>% 
  ggplot(aes(x = value, 
             y = reorder(Taxon_label, ((number_phylum * 100) + number_genus + (number_taxon / 1e6)) * -1),
             fill = name)) +  # Map fill to facet column variable
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey", linewidth = 0.25) +
  geom_hline(data = hline_positions, aes(yintercept = intercept), 
             color = "grey", linetype = "dashed", linewidth = 0.25) + # Thinner hline, matching vline
  geom_col() +
  scale_fill_manual(values = facet_colors) + # Set manual bar colors per facet
  facet_wrap(~ name, nrow = 1) +
  labs(x = "Log2 fold change") + # Set x-axis label (Log2 fold change)
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 7),
    axis.text.y = element_text(size = 7, vjust = 0.25),
    axis.text.x = element_text(size = 7),
    strip.text = element_blank(),     # Remove facet title labels
    strip.background = element_blank(), # Remove background box of facet headers
    legend.position = "none"          # Hide fill legend since facets already separate them
  )

p_main

ggsave(plot = p_main,
       filename = "C:/Project/5_16s_thai_population/figure/revision_figures/07092026_p_main_LFC_genusSpecies.png",
       height = 7.25, width = 7.25, units = "in", dpi = 600)

my_colors <- c("Actinobacteriota" = "#ACCEE3",
  "Bacteroidota" = "#548EC0",
  "Desulfobacterota" = "#4CA330",
  "Euryarchaeota" = "#E66968",
  "Firmicutes" = "#DA362A",
  "Fusobacteriota" = "#F0AA63",
  "Proteobacteria" = "#E49048",
  "Verrucomicrobiota" = "#E6D27A")

p_strip <- plot_data_prepared %>%
  select(Taxon_label, Genus, Phylum) %>%
  distinct() %>%
  ggplot(aes(x = "", y = Taxon_label, fill = Phylum)) +
  geom_tile(color = "white", linewidth = 0.3) + # Thin white borders (ComplexHeatmap look)
  scale_x_discrete(expand = c(0, 0), position = "top") + # Remove x-padding & place header at top
  scale_y_discrete(expand = c(0, 0)) + # Remove y-padding so tiles touch edge-to-edge
  scale_fill_manual(values = my_colors)+
  theme_minimal(base_size = 10) +
  theme(
    # Clean grid/panel elements
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.5),
    
    # Align axes and titles like ComplexHeatmap headers
    axis.text.x = element_blank(),
    axis.text.y = element_blank(), # Keep labels on far left
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    
    # Clean, compact legend layout
    legend.position = "left",
    legend.key.size = unit(3.5, "mm"),
    legend.margin = margin(r = 5)
  )

ggsave(plot = p_strip,
       filename = "C:/Project/5_16s_thai_population/figure/revision_figures/07092026_p_strip_phylum.png",
       height = 6.95, width = 1.55, units = "in", dpi = 600)

