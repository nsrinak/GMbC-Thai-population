library(tidyverse)
library(stringr)
library(cowplot)
library(ggplot2)
library(ggvenn)

## LM alpha diversity results ----

lm_alpha <- read.csv("C:/Project/5_16s_thai_population/revise_aftermSystems/24082026_LMfromLasso_AlphaDivesr_withFDR_80cutoff.csv")

lm_alpha_modelFDR <- data.frame()

for (i in c("bangkok", "phatthalung", "tak")) {
  res_loc <- lm_alpha %>% filter(Locality == i) %>% 
    select(Diversity, adj.r.squared, Locality, model_pvalue) %>% unique()
  fdr <- p.adjust(p = res_loc$model_pvalue, method = "BH")
  
  res_loc$modelFDR <- fdr
  lm_alpha_modelFDR <- bind_rows(lm_alpha_modelFDR, res_loc)
}

p_bar_adjRSqure_alpha <- lm_alpha_modelFDR %>% 
  filter(modelFDR < 0.05) %>% 
  ggplot(aes(y = reorder(Diversity, adj.r.squared, FUN = "max"), x = adj.r.squared, fill = Locality))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_bw(base_size = 7)+
  scale_x_continuous(
    labels = scales::label_number(accuracy = 0.1)
  )+
  theme(legend.position = "none")+
  labs(y="", x = "Adjusted R-Squared")+
  facet_wrap(
    ~ Locality, ncol = 1,
    labeller = as_labeller(c(
      "bangkok" = "Bangkok",
      "tak" = "Tak",
      "phatthalung" = "Phatthalung"
    ))
  )

ggsave(plot = p_bar_adjRSqure_alpha,
       filename = "C:/Project/5_16s_thai_population/figure/revision_figures/07092026_p_bar_adjRSqure_alpha_1col.png", 
       width = 1.4, height = 4, units = "in", dpi = 300)

library(ggtext)

host_genetic <- c("sexM")

p_heatmap_alpha_LM <- lm_alpha %>% 
  filter(FDR < 0.05) %>% 
  mutate(Locality = case_when(
           Locality == "bangkok" ~ "Bangkok",
           Locality == "phatthalung" ~ "Phatthalung",
           Locality == "tak" ~ "Tak"
         )) %>% 
  mutate(term = as.character(term),
         term = if_else(term %in% host_genetic, paste0("**", term, "**"), term)) %>% 
  ggplot(aes(y = Diversity, x = term, fill = t_value))+
  geom_tile()+
  scale_fill_gradientn(
    colors = c("#723194", "white", "#8dce6b"),
    limits = c(-12, 12),
    na.value = "white")+
  facet_grid(
      . ~ Locality,
      scales = "free_x",
      space = "free_x"
    )+
  labs(x = "", y = "", fill = "t value")+
  theme_bw() +
  theme(
    axis.text.x = element_markdown(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 6.5
    ),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "bottom"
  )

p_heatmap_alpha_LM

ggsave(plot = p_heatmap_alpha_LM,
       filename = "C:/Project/5_16s_thai_population/figure/revision_figures/07092026_p_heatmap_alpha_LM_unlimitScale.png", 
       width = 3.75, height = 2.25, units = "in", dpi = 300)

## LM abundance result ----

lm_abundance <- read.csv("C:/Project/5_16s_thai_population/revise_aftermSystems/24082026_LMfromLasso_Abundance_withFDR_80cutoff_asv_pairwiseComplete_CLR.csv")

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
    Taxon_label = make.unique(Taxon_label, sep = "_")) %>% 
  rownames_to_column(var = "OTU")

lm_abundance <- lm_abundance %>% 
  left_join(., taxa_table, by = "OTU")

write.csv(x = lm_abundance, file = "C:/Project/5_16s_thai_population/revise_aftermSystems/21092026_lm_abundance_with_GenusName.csv", row.names = FALSE)


lm_abundance_modelFDR <- data.frame()

for (i in c("bangkok", "phatthalung", "tak")) {
  res_loc <- lm_abundance %>% 
  filter(Locality == i) %>% 
  mutate(Genus = Taxon_label) #%>% 
  #select(Genus, adj.r.squared, Locality, model_pvalue) %>% unique()
  fdr <- p.adjust(p = res_loc$model_pvalue, method = "BH")
  
  res_loc$modelFDR <- fdr
  lm_abundance_modelFDR <- bind_rows(lm_abundance_modelFDR, res_loc)
}

sele_ge <- lm_abundance_modelFDR %>% 
  filter(modelFDR < 0.05) %>% 
  group_by(Locality) %>% 
  slice_max(n = 20, order_by =adj.r.squared) %>% 
  pull(Genus) %>% unique()

p_bar_adjRSqure_genus <- lm_abundance_modelFDR %>% 
  filter(modelFDR < 0.05, !is.na(Genus),
         Genus %in% sele_ge) %>% 
  ggplot(aes(y = reorder(Genus, adj.r.squared, FUN = "max"), x = adj.r.squared, fill = Locality))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_bw(base_size = 8.5)+
  theme(legend.position = "none")+
  labs(y="", x = "Adjusted R-Squared")+
  facet_wrap(
    ~ Locality, ncol = 3,
    labeller = as_labeller(c(
      "bangkok" = "Bangkok",
      "tak" = "Tak",
      "phatthalung" = "Phatthalung"
    ))
  )

ggsave(plot = p_bar_adjRSqure_genus,
       filename = "C:/Project/5_16s_thai_population/figure/revision_figures/07092026_p_bar_adjRSqure_genus_asv_complete.pairwise_AllLabel_top20.png", 
       width = 5, height = 5, units = "in", dpi = 600)

#------
abun_bkk <- lm_abundance_modelFDR %>% filter(modelFDR < 0.05 & Locality == "bangkok", !is.na(Genus)) %>% pull(Genus)
abun_phat <- lm_abundance_modelFDR %>% filter(modelFDR < 0.05 & Locality == "phatthalung", !is.na(Genus)) %>% pull(Genus)
abun_tak <- lm_abundance_modelFDR %>% filter(modelFDR < 0.05 & Locality == "tak", !is.na(Genus)) %>% pull(Genus)

unique_bkk <- setdiff(abun_bkk, union(abun_phat, abun_tak))
unique_phat <- setdiff(abun_phat, union(abun_bkk, abun_tak))
unique_tak <- setdiff(abun_tak, union(abun_bkk, abun_phat))
all_unique <- unique(c(
  setdiff(abun_bkk,  union(abun_phat, abun_tak)),
  setdiff(abun_phat, union(abun_bkk, abun_tak)),
  setdiff(abun_tak,  union(abun_bkk, abun_phat))
))

all_intersect_taxa <- intersect(intersect(abun_bkk, abun_phat), abun_tak)
  
venn_list_all <- list("Bangkok"=abun_bkk,
                      "Phatthalung"=abun_phat,
                      "Tak"=abun_tak)

p_venn_genus_sigLM <-ggvenn(
  venn_list_all,
  fill_color = c("#F8766D", "#00BA38", "#619CFF"), # Pastel colors
  stroke_color = "black", stroke_size = 0.5,        # Remove border lines
  #set_name_size = 2,        # Adjust set label size
  text_size = 2.7,             # Adjust text inside the regions
  show_percentage = FALSE
)

ggsave(plot = p_venn_genus_sigLM ,
       filename = "C:/Project/5_16s_thai_population/figure/revision_figures/07092026_p_venn_genus_sigLM_asv_complete.pairwise_allLabel.png",
       width = 1.8, height = 1.8)


lifestyle_bkk <- lm_abundance %>% filter(FDR < 0.05 & Locality == "bangkok") %>% pull(term) %>% unique()
lifestyle_phat <- lm_abundance %>% filter(FDR < 0.05 & Locality == "phatthalung") %>% pull(term) %>% unique()
lifestyle_tak <- lm_abundance %>% filter(FDR < 0.05 & Locality == "tak") %>% pull(term) %>% unique()

unique_bkk <- setdiff(lifestyle_bkk, union(lifestyle_phat, lifestyle_tak))
unique_phat <- setdiff(lifestyle_phat, union(lifestyle_bkk, lifestyle_tak))
unique_tak <- setdiff(lifestyle_tak, union(lifestyle_bkk, lifestyle_phat))
all_unique <- unique(c(
  setdiff(lifestyle_bkk,  union(lifestyle_phat, lifestyle_tak)),
  setdiff(lifestyle_phat, union(lifestyle_bkk, lifestyle_tak)),
  setdiff(lifestyle_tak,  union(lifestyle_bkk, lifestyle_phat))
))

all_union <- intersect(intersect(lifestyle_bkk, lifestyle_phat), lifestyle_tak)

venn_list_all <- list("Bangkok"=lifestyle_bkk,
                      "Phatthalung"=lifestyle_phat,
                      "Tak"=lifestyle_tak)

p_venn_lifestyle_sigLM <-ggvenn(
  venn_list_all,
  fill_color = c("#F8766D", "#00BA38", "#619CFF"), # Pastel colors
  stroke_color = "black", stroke_size = 0.5,        # Remove border lines
  #set_name_size = 2,        # Adjust set label size
  text_size = 2.7,             # Adjust text inside the regions
  show_percentage = FALSE
)

ggsave(plot = p_venn_lifestyle_sigLM ,
       filename = "C:/Project/5_16s_thai_population/figure/revision_figures/24082026_p_venn_lifestyle_sigLM_asv_complete.pairwise.png",
       width = 1.8, height = 1.8)




sele_ge <- lm_abundance_modelFDR %>% 
  filter(modelFDR < 0.05) %>% 
  group_by(Locality) %>% 
  slice_max(n = 20, order_by =adj.r.squared) %>% 
  pull(Genus) %>% unique()


host_genetic <- c("sexM", "height_cm", "age")

plot_main_df <- lm_abundance %>% 
  mutate(Locality = case_when(
    Locality == "bangkok" ~ "Bangkok",
    Locality == "phatthalung" ~ "Phatthalung",
    Locality == "tak" ~ "Tak"
  )) %>%
  mutate(
    term = as.character(term),
    term = if_else(term %in% host_genetic, paste0("**", term, "**"), term),
    Genus = if_else(is.na(Genus), "Unknown", Genus)
  ) %>% 
  filter(
    !is.na(Taxon_label),
    FDR < 0.05,
    Taxon_label %in% sele_ge
  ) %>% 
  group_by(Taxon_label) %>% 
  mutate(count_genus = n()) %>% 
  group_by(term) %>% 
  mutate(count_feature = n()) %>% 
  ungroup() %>% 
  filter(count_feature > 3) %>% 
  arrange(Phylum, Genus, Taxon_label) %>% 
  mutate(
    number_phylum = as.numeric(as.factor(Phylum)),
    number_genus = as.numeric(as.factor(Genus)),
    number_taxon = as.numeric(as.factor(Taxon_label)),
    number_genus = if_else(Genus == "Unknown", number_genus * 20, number_genus),
    # 1. Create a sort_key column using your formula
    sort_key = ((number_phylum * 100) + number_genus + (number_taxon / 1e6)) * -1,
    # 2. Lock Taxon_label as a factor with this explicit order
    Taxon_label = reorder(Taxon_label, sort_key)
  )

# Calculate hline positions based on the newly locked factor levels
hline_positions <- plot_main_df %>%
  select(Taxon_label, Genus) %>%
  distinct() %>%
  # 3. Extract the actual 1..N index matching the plot y-axis
  mutate(y_pos = as.numeric(Taxon_label)) %>% 
  group_by(Genus) %>%
  summarise(max_pos = max(y_pos), .groups = "drop") %>%
  filter(max_pos < max(as.numeric(plot_main_df$Taxon_label))) %>%
  transmute(intercept = max_pos + 0.5)

p_lm_abundance <- plot_main_df %>% 
  ggplot(aes(
    x = reorder(term, count_feature, FUN = "max"),
    y = Taxon_label, # Factor levels are already set above
    fill = t_value
  )) +
  geom_tile() +
  geom_hline(data = hline_positions, aes(yintercept = intercept), color = "grey", linetype = "dashed", linewidth = 0.25) +
  scale_fill_gradientn(
    colours = c("#723194", "white", "#8dce6b"),
    limits = c(-12, 12),
    oob = scales::squish,
    name = "t-value"
  ) +
  facet_grid(
    . ~ Locality,
    scales = "free_x",
    space = "free_x"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_markdown(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 6.5
    ),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "bottom"
  ) +
  labs(x = "Lifestyle feature", y = "Taxa")

p_lm_abundance

ggsave(plot = p_lm_abundance ,
       filename = "C:/Project/5_16s_thai_population/figure/revision_figures/07092026_p_lm_abundance_allLabel.png",
       width = 7.5, height = 6.5, dpi = 600)


my_colors <- c("Actinobacteriota" = "#ACCEE3",
               "Bacteroidota" = "#548EC0",
               "Desulfobacterota" = "#4CA330",
               "Euryarchaeota" = "#E66968",
               "Firmicutes" = "#DA362A",
               "Fusobacteriota" = "#F0AA63",
               "Proteobacteria" = "#E49048",
               "Verrucomicrobiota" = "#E6D27A")

p_strip <- plot_main_df %>%
  select(Taxon_label, Genus, Phylum, number_phylum, number_genus, number_taxon) %>%
  distinct() %>%
  ggplot(aes(x = "", 
             y = reorder(Taxon_label, ((number_phylum * 100) + number_genus + (number_taxon / 1e6)) * -1), 
             fill = Phylum)) +
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
    legend.key.size = unit(3.5, "mm")
  )

p_strip

ggsave(plot = p_strip ,
       filename = "C:/Project/5_16s_thai_population/figure/revision_figures/09092026_p_stripp_lm_abundance_allLabel.png",
       width = 1.62, height = 4.78, dpi = 600)
