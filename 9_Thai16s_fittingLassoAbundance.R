library(tidyverse)
library(glmnet)
library(ggplot2)

# Load the data
alpha_diversity_2 <- readRDS("C:/Project/5_16s_thai_population/revise_aftermSystems/16082026_diversityMetadataSampleDF_pairwise.complete.rds")

asv_table <- readRDS("C:/Project/5_16s_thai_population/seqtab_final.rds")
asv_table <- t(asv_table)
taxa_table <- readRDS("C:/Project/5_16s_thai_population/tax_final.rds")

ASV <- otu_table(asv_table, taxa_are_rows = TRUE)
TAX <- tax_table(taxa_table)

physeq <- phyloseq(ASV, TAX)

taxonomic_df_all <- psmelt(physeq)

# CLR transformation
# 1. Create Sample × OTU count matrix

otu_mat <- taxonomic_df_all %>%
  select(Sample, OTU, Abundance) %>%
  pivot_wider(
    names_from = OTU,
    values_from = Abundance,
    values_fill = 0
  ) %>%
  column_to_rownames("Sample") %>%
  as.matrix()

# 2. Filter samples: library size >= 2500

lib_size <- rowSums(otu_mat)
otu_mat <- otu_mat[lib_size >= 2500, , drop = FALSE]

# 3. Calculate relative abundance

rel_abund <- otu_mat / rowSums(otu_mat)

# 4. Filter OTUs by prevalence >= 10%

otu_prevalence <- colMeans(otu_mat > 0)
keep_prevalence <- otu_prevalence >= 0.10

# 5. Filter OTUs reaching >= 1% relative abundance
#    in at least one sample

otu_max_rel_abund <- apply(rel_abund, 2, max)
keep_abundance <- otu_max_rel_abund >= 0.01

# 6. Apply both OTU filters

keep_otu <- keep_prevalence & keep_abundance
otu_mat_filt <- otu_mat[, keep_otu, drop = FALSE]

# 7. CLR transformation

otu_clr <- compositions::clr(otu_mat_filt + 0.5)
otu_clr <- as.data.frame(otu_clr) %>%
  rownames_to_column("Sample")

# Combine two tables together
taxa_info <- taxonomic_df_all %>% select(OTU, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% distinct()

abundance_df <- otu_clr %>%  
  pivot_longer(cols = -Sample, names_to = "OTU", values_to = "CLR_val") %>% 
  left_join(., alpha_diversity_2, by =c("Sample"))

write_rds(x = abundance_df, file = "C:/Project/5_16s_thai_population/revise_aftermSystems/24082026_CLR_abundace_table.rds")

abundance_df <- readRDS(file = "C:/Project/5_16s_thai_population/revise_aftermSystems/24082026_CLR_abundace_table.rds")

unq_otu <- abundance_df %>% pull(OTU) %>% unique()

loc <- unique(alpha_diversity_2$locality)

all_results <- list()

nn <- 1

for (i in loc) {
  for (j in unq_otu) {
    
    message("Processing ", nn, i, " - ", j)
    
    nn <- nn+1
    
    df <- abundance_df %>% 
      filter(locality == i, OTU == j) %>% 
      #Remove columns containing any NA
      select(where(~ !any(is.na(.))))
    
    y <- df %>% 
      select(-c("OTU", "Sample", "CLR_val", "Richness", "Shannon","Simpson", "Evenness", "locality", "PD", "PSVs")) %>% 
      # Keep a column if it is non-numeric OR if it is numeric and has at least 3 distinct values.
      select(where(~ !is.numeric(.) || n_distinct(.) >= 3)) %>% 
      # Remove variables with only 1 level
      select(where(~ n_distinct(., na.rm = TRUE) >= 2)) %>% 
      model.matrix(~ . - 1, data = .)
    
    features <- colnames(y)
    
    selection_matrix <- matrix(0, nrow = length(features), ncol = 400)
    rownames(selection_matrix) <- features
    
    discovered <- c()
    cumulative <- c()
    
    # Here the iteration is increased from 300 to 400 times since the number of taxa increase when analysis on ASV level.
    # With 300 interactions, the cumulative feature is still reach steady state nicely, but to be even more certain, 
    # I increase 100 iterations.
    for (k in 1:400) {
      fit_ok <- TRUE
      
      # When using ASV level, there are more chance for a ASV to have 0 count. Thus, it is highly likely for a random sampling in LASSO will 
      # pick all ZERO ASVs leading to error.
      # So in this version I use tryCatch to prevent the loop from breaking.
      cv_fit <- tryCatch({
        cv.glmnet(y, df$CLR_val, alpha = 1)
      }, error = function(e) {
        message("Iteration ", k, " failed: ", e$message)
        fit_ok <<- FALSE
        NULL
      })
      if (fit_ok) {
        final_model <- glmnet(y, df$CLR_val, lambda = cv_fit$lambda.min, alpha = 1)
        
        selected <- setdiff(
          rownames(coef(final_model))[coef(final_model)[,1] != 0],
          "(Intercept)")
        
        selection_matrix[, k] <- as.integer(features %in% selected)
        
        discovered <- union(discovered, selected)
        
      } else {
        
        # Fill failed iteration with NA
        selection_matrix[, k] <- NA
        
      }
      
      cumulative[k] <- length(discovered)
    }
    
    # Feature selection frequency
    freq <- rowSums(selection_matrix, na.rm = TRUE) / 400
    
    # Cutoff summary
    cutoff_df <- data.frame(
      Cutoff = seq(0.1, 1, 0.1) * 100,
      NumFeatures = sapply(seq(0.1, 1, 0.1), function(c) sum(freq >= c)),
      Locality = i,
      Variable = j
    )
    
    # Cumulative features found over time
    cumulative_df <- data.frame(
      Run = 1:400,
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

# Depending on location, some ASV with certain fail iteration (filled with NA) should be filtered out

# Combine results across all runs
cutoff_all <- bind_rows(lapply(all_results, `[[`, "cutoff"))
cumulative_all <- bind_rows(lapply(all_results, `[[`, "cumulative"))
features_all <- bind_rows(lapply(all_results, `[[`, "features"))
consensus_features_all <- bind_rows(lapply(all_results, `[[`, "consensus"))

# Optionally group features for downstream use
consensus_feature_list <- consensus_features_all %>%
  group_by(Locality, Variable) %>%
  summarise(Features = list(Feature), .groups = "drop")

write.csv(cutoff_all, file = "C:/Project/5_16s_thai_population/revise_aftermSystems/24082026_Lassocutoff_all_abundance_asv_pairwise.complete_CLR.csv")
write.csv(cumulative_all, file = "C:/Project/5_16s_thai_population/revise_aftermSystems/24082026_Lassocumulative_all_abundance_asv_pairwise.complete_CLR.csv")
write.csv(features_all, file = "C:/Project/5_16s_thai_population/revise_aftermSystems/24082026_Lassofeatures_all_abundance_asv_pairwise.complete_CLR.csv")
write.csv(consensus_features_all, file = "C:/Project/5_16s_thai_population/revise_aftermSystems/24082026_Lassoconsensus_features_all_abundance_asv_pairwise.complete_CLR.csv")

p_cutoff <- ggplot(cutoff_all, aes(x= Cutoff, y = NumFeatures, colour = Variable))+
  geom_line()+
  facet_wrap(~ Locality)+
  theme_bw()+
  theme(legend.position = "none")


p_cumulative <- ggplot(cumulative_all, aes(x = Run, y=CumulativeFeatures, colour = Variable))+
  geom_line()+
  facet_wrap(~ Locality)+
  theme_bw()+
  theme(legend.position = "none")

ggsave(plot = p_cutoff,
       filename = "C:/Project/5_16s_thai_population/figure/revision_figures/24082026_p_cutoff_featureLasso_abundance_asv_pairwise.complete_CLR.png", 
       width = 8, height = 5, units = "in")

ggsave(plot = p_cumulative,
       filename = "C:/Project/5_16s_thai_population/figure/revision_figures/24082026_p_cumulative_featureLasso_abundance_asv_pairwise.complete_CLR.png", 
       width = 8, height = 5, units = "in")
