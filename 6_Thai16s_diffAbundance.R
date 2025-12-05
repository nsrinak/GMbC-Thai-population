library(tidyverse)
library(vegan)
library(ggplot2)
library(cowplot)
library(ANCOMBC)


asv_table <- readRDS("C:/Project/5_16s_thai_population/seqtab_final.rds")
asv_table <- t(asv_table)

taxa_table <- readRDS("C:/Project/5_16s_thai_population/tax_final.rds")

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
ancombc_res <- ancombc2(
  data = physeq,              
  tax_level = "Genus",
  fix_formula = "locality + age + sex + bmi",
  rand_formula = NULL,
  prv_cut = 0.1, 
  lib_cut = 2500, 
  s0_perc = 0.05,
  group = "locality", 
  p_adj_method = "BH",
  pseudo_sens = TRUE,
  struc_zero = TRUE, 
  neg_lb = TRUE,
  alpha = 0.05, 
  n_cl = 8, 
  verbose = TRUE,
  global = TRUE, 
  pairwise = TRUE, 
  dunnet = TRUE, 
  trend = TRUE,
  iter_control = list(tol = 1e-2, max_iter = 20, 
                      verbose = TRUE),
  em_control = list(tol = 1e-5, max_iter = 100),
  lme_control = lme4::lmerControl(),
  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
  trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                              nrow = 2, 
                                              byrow = TRUE),
                                       matrix(c(-1, 0, 1, -1),
                                              nrow = 2, 
                                              byrow = TRUE),
                                       matrix(c(1, 0, 1, -1),
                                              nrow = 2, 
                                              byrow = TRUE)),
                       node = list(2, 2, 1),
                       solver = "ECOS",
                       B = 100))

saveRDS(ancombc_res, file = "C:/Project/5_16s_thai_population/25032025_ancombc_results_tryBMIAsFixEffect.rds")

ancombc_res <- readRDS("C:/Project/5_16s_thai_population/25032025_ancombc_results_tryBMIAsFixEffect.rds")

res_pairwise <- ancombc_res$res_pair %>% 
  mutate(qq_bkk_phat = if_else(q_localityphatthalung<0.05,"S", "NS")) %>% 
  mutate(qq_bkk_tak = if_else(q_localitytak<0.05,"S", "NS")) %>% 
  mutate(qq_tak_phat = if_else(q_localitytak_localityphatthalung<0.05,"S", "NS")) 
glimpse(res_pairwise)

p_bkk_phat <- res_pairwise %>% 
  filter(qq_bkk_phat == "S") %>% 
  ggplot(aes(y=reorder(taxon, lfc_localityphatthalung, FUN = "max"), x=lfc_localityphatthalung))+
  geom_bar(stat = "identity", fill = "#FFB7C5")+
  theme_classic()+
  xlim(-3,3)+
  labs(x= expression(Log[2]~"fold change"),
       y="",
       title="Phatthalung vs Bangkok")+
  theme(legend.position = "none")


p_bkk_tak <- res_pairwise %>% 
  filter(qq_bkk_tak == "S") %>% 
  ggplot(aes(y=reorder(taxon, lfc_localitytak, FUN = "max"), x=lfc_localitytak))+
  geom_bar(stat = "identity", fill = "#B2EBF2")+
  theme_classic()+
  xlim(-3,3)+
  labs(x=expression(Log[2]~"fold change"),
       y="",
       title="Tak vs Bangkok")+
  theme(legend.position = "none")

p_tak_phatthalung <- res_pairwise %>% 
  filter(qq_tak_phat == "S") %>% 
  ggplot(aes(y=reorder(taxon, lfc_localitytak_localityphatthalung, FUN = "max"), x=lfc_localitytak_localityphatthalung, fill = "#AEC6CF"))+
  geom_bar(stat = "identity", fill = "#AEC6CF")+
  theme_classic()+
  xlim(-3,3)+
  labs(x=expression(Log[2]~"fold change"),
       y="",
       title="Tak vs Phatthalung")+
  theme(legend.position = "none")

ggsave(plot = p_bkk_phat,
       filename = "C:/Project/5_16s_thai_population/figure/02092025_ANCOM_bkk_phat.png",
       height = 12, width = 10, units = "cm")

ggsave(plot = p_bkk_tak,
       filename = "C:/Project/5_16s_thai_population/figure/02092025_ANCOM_bkk_tak.png",
       height = 7.3, width = 10, units = "cm")

ggsave(plot = p_tak_phatthalung,
       filename = "C:/Project/5_16s_thai_population/figure/02092025_ANCOM_tak_phat.png",
       height = 4.2, width = 10, units = "cm")

view(ancombc_res$res_pair)

write.csv(ancombc_res$res_pair, file = "C:/Project/5_16s_thai_population/24072025_ANCOMBC2_result_pairwise_BMI.csv")
