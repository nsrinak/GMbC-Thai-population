library(tidyverse)
library(glmnet)
library(ggplot2)

alpha_diversity_2 <- read.csv("C:/Project/5_16s_thai_population/22052005_diversityMetadataSampleDF.csv")

ancombc_res <- readRDS("C:/Project/5_16s_thai_population/221224_ancombc_results.rds")

abundance_df <- ancombc_res$feature_table %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Genus") %>%
  pivot_longer(cols = -Genus,names_to = "Sample", values_to = "bc_abundance") %>% 
  left_join(., alpha_diversity_2, by =c("Sample"))

unq_genus <- abundance_df %>% filter(Genus != "Unknown") %>% pull(Genus) %>% unique()
loc <- unique(alpha_diversity_2$locality)

all_results <- list()

for (i in loc) {
  for (j in unq_genus) {
    
    message("Processing ", i, " - ", j)
    
    df <- abundance_df %>% drop_na() %>% filter(locality == i, Genus == j)
    y <- df %>% select(-c("Genus", "Sample", "bc_abundance", "Richness", "Shannon",
                          "Simpson", "Evenness", "locality", "PD", "PSVs", "FaithPD_Diet")) %>% 
      select(where(~ n_distinct(.) >= 3)) %>% 
      model.matrix(~ . - 1, data = .)
    
    selection_matrix <- matrix(0, nrow = length(features), ncol = 300)
    rownames(selection_matrix) <- features
    
    discovered <- c()
    cumulative <- c()
    
    for (k in 1:300) {
      cv_fit <- cv.glmnet(y, df$bc_abundance, alpha = 1)
      final_model <- glmnet(y, df$bc_abundance, lambda = cv_fit$lambda.min, alpha = 1)
      selected <- setdiff(rownames(coef(final_model))[coef(final_model)[, 1] != 0], "(Intercept)")
      
      selection_matrix[, k] <- as.integer(features %in% selected)
      discovered <- union(discovered, selected)
      cumulative[k] <- length(discovered)
    }
    # Feature selection frequency
    freq <- rowSums(selection_matrix) / 300
    
    # Cutoff summary
    cutoff_df <- data.frame(
      Cutoff = seq(0.1, 1, 0.1) * 100,
      NumFeatures = sapply(seq(0.1, 1, 0.1), function(c) sum(freq >= c)),
      Locality = i,
      Variable = j
    )
    
    # Cumulative features found over time
    cumulative_df <- data.frame(
      Run = 1:300,
      CumulativeFeatures = cumulative,
      Locality = i,
      Variable = j
    )
    
    # All features and their selection frequency
    feature_df <- data.frame(
      Feature = rownames(selection_matrix),
      Frequency = freq,
      Locality = i,
      Variable = j
    )
    
    # Consensus set: features found in ≥ 80% of runs
    threshold <- 0.80
    consensus_df <- feature_df %>%
      filter(Frequency >= threshold)
    
    # Store results
    all_results[[paste(i, j, sep = "_")]] <- list(
      cutoff = cutoff_df,
      cumulative = cumulative_df,
      features = feature_df,
      consensus = consensus_df
    )
  }
}

# Combine results across all runs
cutoff_all <- bind_rows(lapply(all_results, `[[`, "cutoff"))
cumulative_all <- bind_rows(lapply(all_results, `[[`, "cumulative"))
features_all <- bind_rows(lapply(all_results, `[[`, "features"))
consensus_features_all <- bind_rows(lapply(all_results, `[[`, "consensus"))

# Optionally group features for downstream use
consensus_feature_list <- consensus_features_all %>%
  group_by(Locality, Variable) %>%
  summarise(Features = list(Feature), .groups = "drop")

write.csv(cutoff_all, file = "C:/Project/5_16s_thai_population/25052025_Lassocutoff_all_abundance.csv")
write.csv(cumulative_all, file = "C:/Project/5_16s_thai_population/25052025_Lassocumulative_all_abundance.csv")
write.csv(features_all, file = "C:/Project/5_16s_thai_population/25052025_Lassofeatures_all_abundance.csv")
write.csv(consensus_features_all, file = "C:/Project/5_16s_thai_population/25052025_Lassoconsensus_features_all_abundance.csv")

p_cutoff <- ggplot(cutoff_all, aes(x= Cutoff, y = NumFeatures, colour = Variable))+
  geom_line()+
  facet_wrap(~ Locality)+
  theme_bw()

p_cumulative <- ggplot(cumulative_all, aes(x = Run, y=CumulativeFeatures, colour = Variable))+
  geom_line()+
  facet_wrap(~ Locality)+
  theme_bw()

ggsave(plot = p_cutoff,
       filename = "C:/Project/5_16s_thai_population/figure/25052025_p_cutoff_featureLasso_abundance.png", 
       width = 44, height = 20, units = "cm")

ggsave(plot = p_cumulative,
       filename = "C:/Project/5_16s_thai_population/figure/25052025_p_cumulative_featureLasso_abundance.png", 
       width = 44, height = 20, units = "cm")
