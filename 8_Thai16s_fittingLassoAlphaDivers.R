library(tidyverse)
library(glmnet)
library(ggplot2)

alpha_diversity_2 <- readRDS("C:/Project/5_16s_thai_population/revise_aftermSystems/16082026_diversityMetadataSampleDF_pairwise.complete.rds")

# List of response variables
loc <- unique(alpha_diversity_2$locality)

y_vars <- c("Richness", "Evenness", "Shannon", "Simpson", "PD", "PSVs")

# Initialize storage
all_results <- list()

# Loop over locality and diversity variables
for (l in loc) {
  for (y_name in y_vars) {
    
    message("Processing ", l, " - ", y_name)
    
    # Subset and prepare inputs
    loc_alpha_div <- alpha_diversity_2 %>% 
      filter(locality == l) %>% 
      #Remove columns containing any NA
      select(where(~ !any(is.na(.))))
    
    loc_metaData <- loc_alpha_div %>% 
      select(-c("Richness", "Evenness", "Shannon", "Simpson", "PD", "PSVs", "FaithPD_Diet")) %>% 
      unique() %>% 
      column_to_rownames(var = "Sample") %>% 
      # Keep a column if it is non-numeric OR if it is numeric and has at least 3 distinct values.
      select(where(~ !is.numeric(.) || n_distinct(.) >= 3)) %>% 
      # Remove variables with only 1 level
      select(where(~ n_distinct(., na.rm = TRUE) >= 2)) %>% 
      model.matrix(~ . - 1, data = .)
    
    features <- colnames(loc_metaData)
    
    y <- loc_alpha_div[[y_name]]
    
    if (length(unique(y)) <= 1 || any(is.na(y))) next
    
    selection_matrix <- matrix(0, nrow = length(features), ncol = 400)
    rownames(selection_matrix) <- features
    
    discovered <- c()
    cumulative <- c()
    
    for (i in 1:400) {
      cv_fit <- cv.glmnet(loc_metaData, y, alpha = 1)
      final_model <- glmnet(loc_metaData, y, lambda = cv_fit$lambda.min, alpha = 1)
      selected <- setdiff(rownames(coef(final_model))[coef(final_model)[, 1] != 0], "(Intercept)")
      
      selection_matrix[, i] <- as.integer(features %in% selected)
      discovered <- union(discovered, selected)
      cumulative[i] <- length(discovered)
    }
    
    # Feature selection frequency
    freq <- rowSums(selection_matrix) / 400
    
    # Cutoff summary
    cutoff_df <- data.frame(
      Cutoff = seq(0.1, 1, 0.1) * 100,
      NumFeatures = sapply(seq(0.1, 1, 0.1), function(c) sum(freq >= c)),
      Locality = l,
      Variable = y_name
    )
    
    # Cumulative features found over time
    cumulative_df <- data.frame(
      Run = 1:400,
      CumulativeFeatures = cumulative,
      Locality = l,
      Variable = y_name
    )
    
    # All features and their selection frequency
    feature_df <- data.frame(
      Feature = rownames(selection_matrix),
      Frequency = freq,
      Locality = l,
      Variable = y_name
    )
    
    # Consensus set: features found in ≥ 80% of runs
    threshold <- 0.80
    consensus_df <- feature_df %>%
      filter(Frequency >= threshold)
    
    # Store results
    all_results[[paste(l, y_name, sep = "_")]] <- list(
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

write.csv(cutoff_all, file = "C:/Project/5_16s_thai_population/revise_aftermSystems/23082026_Lassocutoff_all_alphaDivesr_pairwise.complete.csv")
write.csv(cumulative_all, file = "C:/Project/5_16s_thai_population/revise_aftermSystems/23082026_Lassocumulative_all_alphaDivesr_pairwise.complete.csv")
write.csv(features_all, file = "C:/Project/5_16s_thai_population/revise_aftermSystems/23082026_Lassofeatures_all_alphaDivesr_pairwise.complete.csv")
write.csv(consensus_features_all, file = "C:/Project/5_16s_thai_population/revise_aftermSystems/23082026_Lassoconsensus_features_all_alphaDivesr_pairwise.complete.csv")

p_cutoff <- ggplot(cutoff_all, aes(x= Cutoff, y = NumFeatures, colour = Variable))+
  geom_line()+
  facet_wrap(~ Locality)+
  theme_bw()

p_cumulative <- ggplot(cumulative_all, aes(x = Run, y=CumulativeFeatures, colour = Variable))+
  geom_line()+
  facet_wrap(~ Locality)+
  theme_bw()

ggsave(plot = p_cutoff,
       filename = "C:/Project/5_16s_thai_population/figure/revision_figures/23082026_p_cutoff_featureLasso.png", 
       width = 16, height = 12, units = "cm")

ggsave(plot = p_cumulative,
       filename = "C:/Project/5_16s_thai_population/figure/revision_figures/23082026_p_cumulative_featureLasso.png", 
       width = 16, height = 12, units = "cm")


