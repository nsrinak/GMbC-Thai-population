library(tidyverse)
library(stringr)
library(cowplot)
library(ggplot2)
library(ggvenn)

## LM alpha diversity results --

lm_alpha <- read.csv("C:/Project/5_16s_thai_population/26052025_LMfromLasso_AlphaDivesr_withFDR_80cutoff.csv")

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
  theme_bw(base_size = 9)+
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

ggsave(plot = p_bar_adjRSqure_alpha,
       filename = "C:/Project/5_16s_thai_population/figure/16112025_p_bar_adjRSqure_alpha.png", 
       width = 2.5, height = 2, units = "in", dpi = 300)

p_heatmap_alpha_LM <- lm_alpha %>% 
  #mutate(sig_loc_genus = paste0(Genus, "_", Locality)) %>% 
  #filter(sig_loc_genus %in% sigLocGenus) %>% 
  filter(FDR < 0.05) %>% 
  mutate(label = ifelse(FDR < 0.05, round(t_value, 2), NA)) %>% 
  ggplot(aes(y = reorder(Diversity, adj.r.squared, FUN = "max"), x = term, fill = t_value))+
  geom_tile()+
  #geom_text(aes(label = label), size = 3, na.rm = TRUE) +
  scale_fill_gradientn(
    colors = c("#723194", "white", "#8dce6b"),
    limits = c(-4, 4),
    na.value = "white")+
  facet_grid(~ Locality, scales = "free_x", space = "free_x",
             labeller = as_labeller(c(
               "bangkok" = "Bangkok",
               "tak" = "Tak",
               "phatthalung" = "Phatthalung"
             )))+
  labs(x = "", y = "", fill = "t value")+
  theme_bw(base_size = 9)+
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 315, hjust =0))

ggsave(plot = p_heatmap_alpha_LM,
       filename = "C:/Project/5_16s_thai_population/figure/16112025_p_heatmap_alpha_LM.png", 
       width = 5, height = 2.5, units = "in", dpi = 300)

## LM abundance result --

lm_abundance <- read.csv("C:/Project/5_16s_thai_population/25052025_LMfromLasso_Abundance_withFDR_80cutoff.csv")

lm_abundance_modelFDR <- data.frame()

for (i in c("bangkok", "phatthalung", "tak")) {
  res_loc <- lm_abundance %>% filter(Locality == i) %>% 
  select(Genus, adj.r.squared, Locality, model_pvalue) %>% unique()
  fdr <- p.adjust(p = res_loc$model_pvalue, method = "BH")
  
  res_loc$modelFDR <- fdr
  lm_abundance_modelFDR <- bind_rows(lm_abundance_modelFDR, res_loc)
}

p_bar_adjRSqure_genus <- lm_abundance_modelFDR %>% 
  filter(modelFDR < 0.05) %>% 
  ggplot(aes(y = reorder(Genus, adj.r.squared, FUN = "max"), x = adj.r.squared, fill = Locality))+
  geom_bar(stat = "identity", position = "dodge")+
  theme_bw(base_size = 10)+
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
       filename = "C:/Project/5_16s_thai_population/figure/02092025_p_bar_adjRSqure_genus.png", 
       width = 5, height = 7, units = "in", dpi = 600)

abun_bkk <- lm_abundance_modelFDR %>% filter(modelFDR < 0.05 & Locality == "bangkok") %>% pull(Genus)
abun_phat <- lm_abundance_modelFDR %>% filter(modelFDR < 0.05 & Locality == "phatthalung") %>% pull(Genus)
abun_tak <- lm_abundance_modelFDR %>% filter(modelFDR < 0.05 & Locality == "tak") %>% pull(Genus)

unique_bkk <- setdiff(abun_bkk, union(abun_phat, abun_tak))
unique_phat <- setdiff(abun_phat, union(abun_bkk, abun_tak))
unique_tak <- setdiff(abun_tak, union(abun_bkk, abun_phat))
all_unique <- unique(c(
  setdiff(abun_bkk,  union(abun_phat, abun_tak)),
  setdiff(abun_phat, union(abun_bkk, abun_tak)),
  setdiff(abun_tak,  union(abun_bkk, abun_phat))
))

venn_list_all <- list("Bangkok"=abun_bkk,
                      "Phatthalung"=abun_phat,
                      "Tak"=abun_tak)

p_venn_genus_sigLM <-ggvenn(
  venn_list_all,
  fill_color = c("#F8766D", "#00BA38", "#619CFF"), # Pastel colors
  stroke_color = "black",        # Remove border lines
  #set_name_size = 2,        # Adjust set label size
  #text_size = 2,             # Adjust text inside the regions
  show_percentage = FALSE
)

ggsave(plot = p_venn_genus_sigLM ,
       filename = "C:/Project/5_16s_thai_population/figure/02092025_p_venn_genus_sigLM .png",
       width = 3.5, height = 3.5)


p_heatmap_abundance_LM <- lm_abundance %>% 
  #mutate(sig_loc_genus = paste0(Genus, "_", Locality)) %>% 
  #filter(sig_loc_genus %in% sigLocGenus) %>% 
  filter(FDR < 0.05) %>% 
  mutate(label = ifelse(FDR < 0.05, round(t_value, 2), NA)) %>% 
  ggplot(aes(y = reorder(Genus, adj.r.squared, FUN = "max"), x = term, fill = t_value))+
  geom_tile()+
  #geom_text(aes(label = label), size = 3, na.rm = TRUE) +
  scale_fill_gradientn(
    colors = c("#723194", "white", "#8dce6b"),
    limits = c(-6, 6),
    na.value = "white")+
  facet_grid(~ Locality, scales = "free_x", space = "free_x",
    labeller = as_labeller(c(
    "bangkok" = "Bangkok",
    "tak" = "Tak",
    "phatthalung" = "Phatthalung"
  )))+
  labs(x = "", y = "")+
  theme_bw(base_size = 9)+
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(plot = p_heatmap_abundance_LM,
       filename = "C:/Project/5_16s_thai_population/figure/02092025_p_heatmap_abundance_LM_label.png", 
       width = 7.5, height = 7.5, units = "in", dpi = 600)
