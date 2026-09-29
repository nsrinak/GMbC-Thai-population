#install.packages(c("randomForest", "missForest"))

library(tidyverse)
library(missForest)
library(randomForest)
library(ggplot2)

# Load data 

meta <- readRDS("C:/Project/5_16s_thai_population/revise_aftermSystems/16082026_diversityMetadataSampleDF_pairwise.complete.rds")

meta <- meta %>% select(-PD, -PSVs, -Richness, -Shannon, -Simpson, -Evenness) %>% column_to_rownames(var = "Sample")

# Remove non-lifestyle feature 

meta <- meta %>% select(-c("sex", "age", "height_cm", "bmi", "calprotectin", 
                           "chromogranin", "igm", "iga", "admixture", "PRS_IBD", "weight_kg"))
  
# Clean up and select predictors
#rf_prep <- meta %>%
#  # Ensure all character columns are recognized as proper categorical factors
#  mutate(across(where(is.character), as.factor)) %>% 
#  select(-c(country, latitude, longitude, locality_density, AnnualMeanTemp,
#            PrecipitationSeasonality, TempSeasonality, AnnualPrecipitation,
#            locality_density_category, majority_subsistence_strategy, self_declared_ethnicity,
#            admixture, urbanism, industrialization, water_source))

rf_prep <- meta

# Impute missing values using Random Forest imputation
set.seed(123) # For reproducibility
imputation_res <- missForest(meta)

# Extract the completed dataset
rf_clean <- imputation_res$ximp

set.seed(123)
rf_model <- randomForest(
  as.factor(locality) ~ ., 
  data = rf_clean, 
  importance = TRUE, 
  ntree = 500 # Number of trees in the forest
)

# Print model summary and the confusion matrix to check performance
print(rf_model)


# 1. Extract variable importance scores
importance_df <- as.data.frame(importance(rf_model)) %>%
  # Move row names (variables) into a column
  tibble::rownames_to_column(var = "Variable") %>%
  # MeanDecreaseAccuracy measures total global predictive power
  select(Variable, MeanDecreaseAccuracy) 

write.csv(x = as.data.frame(importance(rf_model)) %>% tibble::rownames_to_column(var = "Variable") , file = "C:/Project/5_16s_thai_population/revise_aftermSystems/21092026_Important_randomForest.csv", row.names = FALSE)

# 2. Filter for the top 15 most important variables
top_importance <- importance_df %>%
  slice_max(order_by = MeanDecreaseAccuracy, n = 30) %>%
  # Reorder factor levels so the plot prints from highest to lowest score
  mutate(Variable = reorder(Variable, MeanDecreaseAccuracy))

print(top_importance)

# 3. Plot the Clean Feature Importance Chart
p_important_feature <- ggplot(top_importance, aes(x = MeanDecreaseAccuracy, y = Variable)) +
  geom_segment(aes(x = 0, xend = MeanDecreaseAccuracy, y = Variable, yend = Variable), 
               color = "grey70", linewidth = 0.6) +
  geom_point(color = "darkslateblue", size = 2.5) +
  theme_bw(base_size = 9) + 
  labs(
    x = NULL,
    y = NULL
  )

ggsave(plot = p_important_feature, 
       filename = "C:/Project/5_16s_thai_population/figure/revision_figures/16082026_p_important_feature_removeGenetic.png", 
       width = 2.1, height = 3.5, units = "in", dpi = 300)
