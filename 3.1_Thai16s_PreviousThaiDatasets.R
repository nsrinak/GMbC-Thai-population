library(tidyverse)
library(phyloseq)
library(vegan)
library(cowplot)
library(ggplot2)

## Combine all taxonomic information from every studies into the same table ----

path <- "meta_analysis/"

studies <-list.dirs(path = path, recursive = FALSE, full.names = FALSE)

taxonomic_df_all <- data.frame("Study" = character(),
                               "OTU" = character(),
                               "Sample" = character(),
                               "Abundance" = numeric(),
                               "Kingdom" = character(),
                               "Phylum" = character(),
                               "Class" = character(),
                               "Order" = character(),
                               "Family" = character(),
                               "Genus" = character(),
                               "Species" = character())

for (i in studies) {
  asv_table <- readRDS(paste(path, i, "/results/seqtab_final.rds", sep = ""))
  asv_table <- t(asv_table)

  taxa_table <- readRDS(paste(path, i, "/results/tax_final.rds", sep = ""))

  ASV <- otu_table(asv_table, taxa_are_rows = TRUE)
  TAX <- tax_table(taxa_table)
  
  physeq <- phyloseq(ASV, TAX)

  taxonomic_df <- psmelt(physeq)
  taxonomic_df$Study <- i
  
  taxonomic_df_all <- rbind(taxonomic_df_all, taxonomic_df)
}


glimpse(taxonomic_df_all)

# Get metadata from all studies 

metadata_df_all <- data.frame("Study" = character(), 
                              "Sample" = character(),
                              "sample_name" = character(),
                              "locality" = character(),
                              "age" = character(),
                              "sex" = character())

for (i in studies) {
  metadata_df <- read.csv(file = paste(path, i, "/metadata_", i, ".csv", sep = ""))
  metadata_df$Study <- i
  metadata_df$Sample <- metadata_df$sample
  metadata_df$sample <- NULL
  
  metadata_df_all <- rbind(metadata_df_all, metadata_df)
}

glimpse(metadata_df_all)

# Joint to the taxonomic table

taxonomic_df_all <- taxonomic_df_all %>% 
  left_join(., metadata_df_all, by = c("Sample", "Study"))

## Load the taxonomic information of my study ----

asv_table <- readRDS("seqtab_final.rds")
asv_table <- t(asv_table)

taxa_table <- readRDS("tax_final.rds")

meta_table <- read_tsv("metaData/metadata.tsv")
meta_table <- meta_table %>% column_to_rownames(var = "donor_id") 

ASV <- otu_table(asv_table, taxa_are_rows = TRUE)
TAX <- tax_table(taxa_table)
SAMP <- sample_data(meta_table)

physeq <- phyloseq(ASV, TAX, SAMP)

taxonomic_df_thisstudy <- psmelt(physeq)

selected_df <- taxonomic_df_thisstudy[,c("OTU", "Sample", "Abundance", "Kingdom", "Phylum", 
                          "Class", "Order", "Family", "Genus", "Species", 
                          "age", "sex", "locality")] %>% 
  mutate(locality = case_when(locality == "bangkok" ~ "Bangkok",
                              locality == "phatthalung" ~ "Phatthalung",
                              locality == "tak" ~ "Tak"))

selected_df$Study <- "This study"

glimpse(selected_df)

taxaTable_combine <- rbind(taxonomic_df_all %>% select(-sample_name), selected_df)

write.csv(x = taxaTable_combine, file = "metaAnalysis_combinedTaxaTable.csv")

remove_sample <- taxaTable_combine %>% 
  group_by(Sample, Study, locality) %>% 
  summarise(Abundance = sum(Abundance),
            count = n()) %>% 
  filter(Abundance < 1000) %>% pull(Sample)

my_colors <- c(
  "Bangkok"="#0072B2", 
  "Buriram"="#E69F00", 
  "Chachoengsao"="#009E73", 
  "Chiang Mai"="#D55E00",
  "Chiang Rai"="#CC79A7", 
  "Maha Sarakham"="#56B4E9", 
  "Nakhon Pathom"="#F0E442", 
  "Phatthalung"="#000000",
  "Tak"="#999999"
)

p_totalLocality <-taxaTable_combine %>%
  select(locality, Study, Sample) %>% 
  unique() %>% 
  group_by(Study, locality) %>% 
  summarise(number_sample = n(), .groups = "drop") %>% 
  group_by(Study) %>% 
  mutate(total_sample = sum(number_sample)) %>% 
  ggplot(aes(x = reorder(Study, total_sample, FUN = "max"), y = number_sample, fill = locality))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = my_colors)+
  guides(fill = guide_legend(nrow = 3))+
  theme_bw()+
  theme(axis.text.x = element_text(angle =45, hjust = 1, size = 7),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "bottom")+
  labs(x = "", y="Sample", fill = "Locality")

ggsave(plot = p_totalLocality_legend,
       filename = "p_totalLocality_legend.png", 
       width = 4.5, height = 2.75, units = "in", dpi = 300)

## Plot map of Thailand ----

cities <- data.frame(
  name = c(
    "Maha Sarakham",
    "Chiang Rai",
    "Bangkok",
    "Buriram",
    "Nakhon Pathom",
    "Chachoengsao",
    "Phatthalung",
    "Tak",
    "Chiang Mai"
  ),
  
  lon = c(
    103.3000,
    99.8325,
    100.4993,
    103.1030,
    100.0590,
    101.0677,
    100.0320,
    99.1656,
    98.9853
  ),
  
  lat = c(
    16.1851,
    19.9105,
    13.7564,
    14.9943,
    13.8199,
    13.6904,
    7.5927,
    16.8816,
    18.7883
  )
)

library(sf)
library(rnaturalearth)
library(ggrepel)

# Thailand province map
th_provinces <- ne_states(country = "Thailand", returnclass = "sf")

# Map
Yesp_map <- ggplot() +
  geom_sf(data = th_provinces,fill = "grey90",color = "white",linewidth = 0.3) +
  geom_sf(data = st_union(th_provinces),fill = NA,color = "black",linewidth = 0.4) +
  geom_point(data = cities,aes(x = lon, y = lat),color = "steelblue",size = 2.5) +
  geom_text_repel(data = cities,aes(x = lon, y = lat, label = name),vjust = -1,size = 2.8) +
  theme_void()

ggsave(plot = Yesp_map,
       filename = "Yesp_map_meta.png", 
       width = 3.5, height = 5, units = "in")

## Dominant Phylum ----

taxaTable_combine %>% 
  filter(!Sample %in% remove_sample) %>% 
  select(Sample, Study, locality) %>% 
  unique() %>% 
  summarise(co = n())

#Number of samples
#Study              count
#<chr>           <int>
#1 Grant_2019          6
#2 Gruneck_2020        8
#3 Kisuse_2018        56
#4 Senaprom_2024      40
#5 Somnuk_2023        30
#6 Therdtatha_2025    30
#7 Therdtatha_2026    11
#8 This study        106
#9 Vangay_2019       183

#Except for this study
# 364

#Number of sample with < 5000 counts
# 24

#Number of sample after filtering
# A tibble: 9 × 2
#Study              co
#<chr>           <int>
#1 Grant_2019          6
#2 Gruneck_2020        8
#3 Kisuse_2018        56
#4 Senaprom_2024      40
#5 Somnuk_2023        30
#6 Therdtatha_2025    30
#7 Therdtatha_2026    11
#8 This study        106
#9 Vangay_2019       159

#All sample after filtering
# 446

p_domTaxa_number <-taxaTable_combine %>% 
  filter(!Sample %in% remove_sample) %>% 
  group_by(Sample) %>% 
  mutate(RA = Abundance/sum(Abundance),
         max_RA = max(RA)) %>% ungroup() %>% 
  filter(RA == max_RA) %>%
  group_by(Phylum) %>% 
  mutate(med_RA = median(RA)) %>% 
  ggplot(aes(x = reorder(Phylum, med_RA, FUN = "max")))+
  geom_bar(fill = "steelblue")+
  theme_bw()+
  theme(axis.text.x = element_text(angle =45, hjust = 1))+
  labs(x = "Dominant phylum",
      y = "Count")

ggsave(plot = p_domTaxa_number,
       filename = "p_domTaxa_number.png", 
       width = 3, height = 2, units = "in", dpi = 300)

p_domTaxa <- taxaTable_combine %>% 
  filter(!Sample %in% remove_sample) %>% 
  group_by(Sample) %>% 
  mutate(RA = Abundance/sum(Abundance),
         max_RA = max(RA)) %>% ungroup() %>% 
  filter(RA == max_RA) %>% 
  group_by(Phylum) %>% 
  mutate(med_RA = median(RA)) %>% 
  ggplot(aes(x = reorder(Phylum, med_RA, FUN = "max"), y = RA))+
  geom_jitter(shape = 21, fill = "steelblue", alpha = 0.6)+
  geom_boxplot(outliers = FALSE)+
  theme_bw()+
  theme(axis.text.x = element_text(angle =45, hjust = 1))+
  labs(x = "Dominant phylum",
       y = "Relative abundance")

ggsave(plot = p_domTaxa,
       filename = "p_domTaxa.png", 
       width = 3, height = 3, units = "in", dpi = 300)

n_sample <- length(unique(taxaTable_combine[!taxaTable_combine$Sample %in% remove_sample,]$Sample))

## Beta diversity 
# Species level

relative_matrix <- taxaTable_combine %>% 
  filter(!Sample %in% remove_sample) %>%
  mutate(genus_species = paste(Genus, Species, sep = "_")) %>% 
  group_by(locality, Sample, genus_species) %>% 
  mutate(Abundance = sum(Abundance)) %>% 
  group_by(Sample) %>% 
  mutate(Relative_Abundance = Abundance/sum(Abundance)) %>% ungroup() %>% 
  select(Sample, genus_species, Relative_Abundance, locality, Study) %>%
  unique() %>% 
  pivot_wider(names_from = genus_species, values_from = Relative_Abundance, values_fill = 0) %>% 
  column_to_rownames(var = "Sample")

set.seed(123)
bray_curtis <- vegdist(as.matrix(relative_matrix %>% select(-locality, -Study)), method = "bray")

# PCoA using ape ----
pcoa_bray <- ape::pcoa(bray_curtis)

# Variance explained
var_bray <- pcoa_bray$values$Relative_eig[1:2] * 100

# Data frame for plotting
pcoa_bray_df <- data.frame(
  SampleID = row.names(relative_matrix),
  PC1 = pcoa_bray$vectors[, 1],
  PC2 = pcoa_bray$vectors[, 2],
  Locality = relative_matrix$locality,
  Study = relative_matrix$Study
)

# Plot Bray-Curtis PCoA Species level
ord_PCoA_bray_Species <- ggplot(pcoa_bray_df, aes(x = PC1, y = PC2, color = Locality, fill = Locality)) +
  geom_point(size = 1.25, alpha = 0.7, shape = 21) +
  geom_point(
    data = subset(pcoa_bray_df, Study == "This study"),
    shape = 21,
    fill = NA,
    color = "black",
    stroke = 0.5,
    size = 1.25
  )+
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = my_colors)+
  scale_color_manual(values = my_colors)+
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 6)) +
  guides(color = guide_legend(nrow = 3),
         fill = guide_legend(nrow = 3))+
  labs(
    x = paste0("PC1 (", round(var_bray[1], 2), "%)"),
    y = paste0("PC2 (", round(var_bray[2], 2), "%)")
  )

ggsave(plot = ord_PCoA_bray_Species,
       filename = "C:/Project/5_16s_thai_population/figure/metaAnalysis/23082026_ord_PCoA_bray_Species_borderedGMbC.png", 
       width = 2.15, height = 1.65, units = "in", dpi = 600)

## Prevalence all ----

tax_level <- c("Kingdom", "Phylum", "Class", "Order",
                "Family", "Genus", "Species")

modf_species <- taxaTable_combine %>% 
  filter(!Sample %in% remove_sample) %>% 
  mutate(Species = paste(Genus, Species, sep = "_")) %>% 
  filter(Abundance > 0) %>% 
  group_by(Sample) %>% 
  mutate(RA = Abundance/sum(Abundance)) %>% ungroup() 

sens_prev_all <- data.frame("percent" = numeric(),
                           "number" = numeric(),
                           "taxonomic_level" = character())

prev_abund_locality <- data.frame("Sample" = character(),
                                  "locality" = character(),
                                  "Taxa" = character(),
                                  "Abundance" = numeric(),
                                  "RA" = numeric(),
                                  "Prevalence" = numeric(),
                                  "taxonomic_level" = character())

for (t in tax_level) {
  prev_table <- modf_species %>% 
    select(all_of(t), Sample) %>% distinct() %>% 
    filter(!is.na(.data[[t]]), 
           .data[[t]] != "NA_NA") %>% 
    group_by(.data[[t]]) %>% 
    summarise(
      prevalence = n() / n_sample,
      .groups = "drop")
  
  fil_df <- modf_species %>% 
    select(all_of(t), Sample, RA, Abundance, locality) %>% distinct() %>% 
    filter(!is.na(.data[[t]]), 
           .data[[t]] != "NA_NA") %>% 
    group_by(.data[[t]], Sample, locality) %>% 
    summarise(RA = sum(RA, na.rm = TRUE), 
              Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>% 
    left_join(., prev_table, by = t) %>% 
    rename("Taxa" = t)
  
  fil_df$taxonomic_level <- t
  
  prev_abund_locality <- rbind(prev_abund_locality, fil_df)
  
  for (i in seq(1, 100, 1)) {
    len <- dim(prev_table[prev_table$prevalence * 100 > i,])[1]
    sens_prev_df <- data.frame("percent" = i/100,
                             "number" = len,
                             "taxonomic_level" = t)
    sens_prev_all <- rbind(sens_prev_all, sens_prev_df)
  }
}

sens_prev_all <- sens_prev_all[sens_prev_all$number != 0,] 

sens_prev_all$taxonomic_level <- factor(
  sens_prev_all$taxonomic_level,
  levels = c(
    "Kingdom", "Phylum", "Class", "Order",
    "Family", "Genus", "Species"
  )
)

prev_abund_locality$taxonomic_level <- factor(
  prev_abund_locality$taxonomic_level,
  levels = c(
    "Kingdom", "Phylum", "Class", "Order",
    "Family", "Genus", "Species"
  )
)


p_prevalence_all_GS <- ggplot(sens_prev_all[sens_prev_all$taxonomic_level %in% c("Genus", "Species"),], aes(x = percent, y = number, fill = taxonomic_level))+
  geom_line(linetype = "dashed", color = "grey")+
  geom_point(shape = 21, alpha = 0.6)+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x = "Prevalence",
       y = "Number of taxa")+
  facet_wrap(~ taxonomic_level, nrow = 1)


ggsave(plot = p_prevalence_all_GS,
       filename = "p_prevalence_all_GS.png", 
       width = 2.5, height = 1.5, units = "in", dpi = 300)

p_prevalence_abundance_GS <- ggplot(prev_abund_locality[prev_abund_locality$taxonomic_level %in% c("Genus", "Species"), ], aes(x = prevalence, y = RA, fill = taxonomic_level))+
  geom_point(shape = 21, alpha = 0.4)+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x = "Prevalence",
       y = "RA")+
  facet_wrap(~ taxonomic_level)

ggsave(plot = p_prevalence_abundance_GS,
       filename = "p_prevalence_abundance_GS.png", 
       width = 2.5, height = 1.5, units = "in", dpi = 300)

## Prevalence per locality ----

locality <- unique(taxaTable_combine$locality)

prev_per_loc <- data.frame("locality" = character(),
                           "Taxa" = character(),
                           "prevalence" = numeric(),
                           "taxonomic_level" = character())

for (i in locality) {
  n_sample_loc <- modf_species %>% 
    filter(locality == i) %>% 
    summarise(n = n_distinct(Sample)) %>% 
    pull(n)
  for (t in tax_level) {
    prev_table <- modf_species %>%
      filter(locality == i) %>% 
      select(all_of(t), Sample) %>% distinct() %>% 
      filter(!is.na(.data[[t]]), 
             .data[[t]] != "NA_NA") %>% 
      group_by(.data[[t]]) %>% 
      summarise(
        prevalence = n() / n_sample_loc,
        .groups = "drop")
    
    fil_df <- modf_species %>%
      filter(locality == i) %>% 
      select(all_of(t), Sample, RA, Abundance) %>% distinct() %>% 
      filter(!is.na(.data[[t]]), 
             .data[[t]] != "NA_NA") %>% 
      group_by(.data[[t]], Sample) %>% 
      summarise(RA = sum(RA, na.rm = TRUE), 
                Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>% 
      left_join(., prev_table, by = t) %>% 
      rename("Taxa" = t)
    
    fil_df$taxonomic_level <- t
    fil_df$locality <- i
    
    prev_per_loc <- rbind(prev_per_loc, fil_df)
  }
}

prev_per_loc$taxonomic_level <- factor(
  prev_per_loc$taxonomic_level,
  levels = c(
    "Kingdom", "Phylum", "Class", "Order",
    "Family", "Genus", "Species"
  )
)


write.csv(x = prev_per_loc %>% 
            filter(taxonomic_level == "Species") %>% 
            select(-Sample,-RA,-Abundance, -taxonomic_level) %>% 
            unique(), 
          file = "Prevalence_eachProvince_species_ThaiMicrobiomes.csv", 
          row.names = FALSE)

# SPECIES level prevalence ----

bubble_df <- prev_per_loc %>% 
  select(-c(Sample, RA, Abundance)) %>% 
  filter(
    taxonomic_level == "Species",
    prevalence > 0.5
  ) %>%
  mutate(Taxa = sub("_NA", "", Taxa)) %>% 
  distinct() %>% 
  group_by(Taxa) %>% 
  mutate(count = n()) %>% 
  filter(count > 3) %>% 
  ungroup() %>% 
  mutate(
    Taxa = forcats::fct_reorder(Taxa, prevalence)
  )

p_bubble <- ggplot(
  bubble_df,
  aes(
    x = locality,
    y = Taxa,
    #size = prevalence,
    fill = prevalence
  )
) +
  geom_point(
    shape = 21,
    colour = "black",
    alpha = 0.7,
    stroke = 0.3,
    size = 2
  ) +
  scale_size_continuous(
    name = "Prevalence",
    range = c(1, 4)
  ) +
  scale_fill_gradientn(
    name = "Prevalence",
    colours = c(
      "#D0F0E8",
      "#66C2A4",
      "#2C7FB8",
      "#5E3C99"
    ),
    limits = c(0.5, 1),
    breaks = c(0.5, 0.7, 0.9, 1)
  ) +
  labs(
    x = "Locality",
    y = "Species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 9
    ),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.position = "right"
  )

ggsave(plot = p_bubble,
       filename = "p_bubble_prev50count3_nosize.png", 
       width = 5, height = 5.5, units = "in", dpi = 300)