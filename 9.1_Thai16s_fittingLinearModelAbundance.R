library(tidyverse)
library(ggplot2)
library(dplyr)
library(broom)

alpha_diversity_2 <- read.csv("C:/Project/5_16s_thai_population/22052005_diversityMetadataSampleDF.csv")

ancombc_res <- readRDS("C:/Project/5_16s_thai_population/221224_ancombc_results.rds")

abundance_df <- ancombc_res$feature_table %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Genus") %>%
  pivot_longer(cols = -Genus,names_to = "Sample", values_to = "bc_abundance") %>% 
  left_join(., alpha_diversity_2, by =c("Sample"))

consensus_features_all <- read.csv("C:/Project/5_16s_thai_population/25052025_Lassoconsensus_features_all_abundance.csv")
# Optionally group features for downstream use
consensus_feature_list <- consensus_features_all %>%
  filter(Frequency >= 0.8) %>% 
  group_by(Locality, Variable) %>%
  summarise(Features = list(Feature), .groups = "drop")

lm_locality <- consensus_feature_list %>% pull(Locality)
lm_variable <- consensus_feature_list %>%  pull(Variable)
lm_feature <- consensus_feature_list %>% pull(Features)

results_all <- data.frame()

for (i in 1:(length(lm_locality))) {
  abun_df <- abundance_df %>% filter(locality==lm_locality[i], Genus == lm_variable[i])
  fomu <- paste("bc_abundance", " ~ ", paste(lm_feature[i][[1]], collapse = " + "), sep = "")
  res_lm <- lm(fomu, abun_df)
  
  # Extract overall model p-value using ANOVA
  f_stat <- summary(res_lm)$fstatistic
  model_pvalue <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
  
  # Coefficients summary with t-values
  tidy_res <- tidy(res_lm) %>%
    rename(t_value = statistic) %>%
    mutate(Locality = lm_locality[i],
           Genus = lm_variable[i],
           Formula = fomu)
  
  # Model-level summary (R², adj R²)
  model_info <- glance(res_lm) %>%
    select(r.squared, adj.r.squared, sigma, AIC, BIC) %>%
    mutate(Locality = lm_locality[i],
           Genus = lm_variable[i],
           Formula = fomu,
           model_pvalue = model_pvalue)
  
  # Join model-level info to each coefficient row
  tidy_res <- left_join(tidy_res, model_info, by = c("Locality", "Genus", "Formula"))
  
  # Save
  results_all <- bind_rows(results_all, tidy_res)
}

result_modelFDR <- data.frame()

for (i in unique(lm_locality)) {
  res_loc <- results_all %>% filter(Locality == i, term == "(Intercept)")
  fdr <- p.adjust(p = res_loc$p.value, method = "BH")
  
  res_loc$modelFDR <- fdr
  result_modelFDR <- bind_rows(result_modelFDR, res_loc)
}


result_FDR <- data.frame()

for (i  in unique(lm_locality)) {
  res_loc <- results_all %>% filter(Locality == i, term != "(Intercept)")
  fdr <- p.adjust(p = res_loc$p.value, method = "BH")
  
  res_loc$FDR <- fdr
  result_FDR <- bind_rows(result_FDR, res_loc)
}

write.csv(result_FDR, file = "C:/Project/5_16s_thai_population/25052025_LMfromLasso_Abundance_withFDR_80cutoff.csv")

