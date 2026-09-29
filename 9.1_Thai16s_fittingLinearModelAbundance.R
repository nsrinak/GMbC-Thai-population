library(tidyverse)
library(ggplot2)
library(dplyr)
library(broom)

abundance_df <- readRDS("C:/Project/5_16s_thai_population/revise_aftermSystems/24082026_CLR_abundace_table.rds")

consensus_features_all <- read.csv("C:/Project/5_16s_thai_population/revise_aftermSystems/24082026_Lassoconsensus_features_all_abundance_asv_pairwise.complete_CLR.csv")

# Optionally group features for downstream use
consensus_feature_list <- consensus_features_all %>%
  filter(Frequency >= 0.8) %>% 
  group_by(Locality, Variable) %>%
  summarise(Features = list(Feature), .groups = "drop")

# 1. Clean the features directly inside the dataframe pipeline
consensus_feature_list_clean <- consensus_features_all %>%
  filter(Frequency >= 0.8) %>% 
  # Clean the individual feature names FIRST
  mutate(Feature_Clean = case_when(
    str_detect(Feature, "^tobacco") ~ "tobacco",
    str_detect(Feature, "^admixture") ~ "admixture",
    str_detect(Feature, "^sex") ~ "sex",
    str_detect(Feature, "^breast_fed") ~ "breast_fed",
    str_detect(Feature, "^c_section") ~ "c_section",
    TRUE ~ Feature
  )) %>%
  # IMPORTANT: Remove duplicates if multiple sub-features collapsed into the same name
  distinct(Locality, Variable, Feature_Clean) %>%
  # Now group and collapse into a formula string
  group_by(Locality, Variable) %>%
  summarise(
    Formula_String = paste("CLR_val ~ ", paste(Feature_Clean, collapse = " + ")),
    .groups = "drop"
  )

# 2. Run the loop safely using the clean dataframe rows
results_all <- data.frame()

for (i in 1:nrow(consensus_feature_list_clean)) {
  
  # Pull variables safely from the current row
  current_locality <- consensus_feature_list_clean$Locality[i]
  current_otu    <- consensus_feature_list_clean$Variable[i]
  fomu             <- consensus_feature_list_clean$Formula_String[i]
  
  # Filter data
  abun_df <- abundance_df %>% 
    filter(locality == current_locality, OTU == current_otu)
  
  # Skip if no data matches (prevents lm errors)
  if(nrow(abun_df) == 0) next
  
  # Fit the linear model
  res_lm <- lm(as.formula(fomu), data = abun_df)
  
  # Extract overall model p-value using ANOVA
  f_stat <- summary(res_lm)$fstatistic
  model_pvalue <- if(!is.null(f_stat)) pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE) else NA
  
  # Model-level summary (R², adj R²)
  model_info <- glance(res_lm) %>%
    select(r.squared, adj.r.squared, sigma, AIC, BIC) %>%
    mutate(model_pvalue = model_pvalue)
  
  # Coefficients summary with t-values
  tidy_res <- tidy(res_lm) %>%
    rename(t_value = statistic) %>%
    mutate(Locality = current_locality,
           OTU = current_otu,
           Formula = fomu)
  
  # Bind model info side-by-side safely without a join key risk
  tidy_res <- bind_cols(tidy_res, model_info)
  
  # Save results
  results_all <- bind_rows(results_all, tidy_res)
}

result_FDR <- data.frame()

for (i  in unique(consensus_feature_list_clean$Locality)) {
  res_loc <- results_all %>% filter(Locality == i, term != "(Intercept)")
  fdr <- p.adjust(p = res_loc$p.value, method = "BH")
  
  res_loc$FDR <- fdr
  result_FDR <- bind_rows(result_FDR, res_loc)
}


asv_name <- read.csv("C:/Project/5_16s_thai_population/03062026_linked_ASV_combinedName.csv")
asv_name <- asv_name %>% 
  mutate(OTU = taxon) %>% 
  select(OTU, unique_combined_1)

result_FDR <- result_FDR %>% 
  left_join(., asv_name, by = "OTU")

write.csv(result_FDR, file = "C:/Project/5_16s_thai_population/revise_aftermSystems/24082026_LMfromLasso_Abundance_withFDR_80cutoff_asv_pairwiseComplete_CLR.csv")

