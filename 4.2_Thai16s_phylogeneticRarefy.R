library(tidyverse)
library(vegan)
library(phyloseq)
library(cowplot)

asv_table <- readRDS("seqtab_final.rds")
asv_table <- t(asv_table)

taxa_table <- readRDS("tax_final.rds")


meta_table <- read_tsv("metadata.tsv")
meta_table <- meta_table %>% column_to_rownames(var = "donor_id") 

## Rarefying ----

# Rarefaction curve plot all ----
spnumber_all <- specnumber(t(asv_table))
Srare_all <- rarefy(t(asv_table), sample = min(rowSums(t(asv_table))))

png("all_species_rarefaction.png", width = 3.5, height = 3.5, res = 120, units = "in")

plot(spnumber_all, Srare_all,
     xlab = "Observed No. of Species",
     ylab = "Rarefied No. of Species",
     col = "black", pch = 19, 
     cex.lab = 0.7,     # axis label size
     cex.axis = 0.7,)
abline(0, 1, col = "black", lwd = 1.5)
dev.off()

rarecurve_all <- rarecurve(t(asv_table),
                            sample = min(rowSums(t(asv_table))),
                            step = 10,
                            tidy = TRUE)

p_all_rarecurv <- ggplot(rarecurve_all, aes(x = Sample, y = Species))+
  geom_point(colour = "black", size = 0.4, alpha = 0.7)+
  geom_vline(xintercept = min(rowSums(t(asv_table))))+
  theme_bw()+
  labs(x = "Sequencing read", y = "ASV")+
  ylim(0, 300)+xlim(0, 62000)

ggsave(plot = p_all_rarecurv,
       filename = "p_all_rarecurv.png", 
       width = 3.5, height = 3.5, units = "in")

# Rarefaction curve plots for each location ----
sampleID_bkk <- meta_table %>% filter(locality == "bangkok") %>% rownames()
sampleID_tak <- meta_table %>% filter(locality == "tak") %>% rownames()
sampleID_phat <- meta_table %>% filter(locality == "phatthalung") %>% rownames()

matrix_bkk <- t(asv_table[, sampleID_bkk])
spnumber_bkk <- specnumber(matrix_bkk)
Srare_bkk <- rarefy(matrix_bkk, sample = min(rowSums(matrix_bkk)))

matrix_tak <- t(asv_table[, sampleID_tak])
spnumber_tak <- specnumber(matrix_tak)
Srare_tak <- rarefy(matrix_tak, sample = min(rowSums(matrix_tak)))

matrix_phat <- t(asv_table[, sampleID_phat])
spnumber_phat <- specnumber(matrix_phat)
Srare_phat <- rarefy(matrix_phat, sample = min(rowSums(matrix_phat)))

png("bkk_species_rarefaction.png", width = 3.5, height = 3.5, res = 120, units = "in")
plot(spnumber_bkk, Srare_bkk,
     xlab = "Observed No. of Species",
     ylab = "Rarefied No. of Species",
     col = "#F8766D", pch = 19, 
     cex.lab = 0.7,     # axis label size
     cex.axis = 0.7,)
abline(0, 1, col = "black", lwd = 1.5)
dev.off()

png("tak_species_rarefaction.png", width = 3.5, height = 3.5, res = 120, units = "in")
plot(spnumber_tak, Srare_tak,
     xlab = "Observed No. of Species",
     ylab = "Rarefied No. of Species",
     col = "#619CFF", pch = 19, 
     cex.lab = 0.7,     # axis label size
     cex.axis = 0.7,)
abline(0, 1, col = "black", lwd = 1.5)
dev.off()

png("phat_species_rarefaction.png", width = 3.5, height = 3.5, res = 120, units = "in")
plot(spnumber_phat, Srare_phat,
     xlab = "Observed No. of Species",
     ylab = "Rarefied No. of Species",
     col = "#00BA38", pch = 19, 
     cex.lab = 0.7,     # axis label size
     cex.axis = 0.7,)
abline(0, 1, col = "black", lwd = 1.5)
dev.off()

rarecurve_bkk <- rarecurve(matrix_bkk,
                           sample = min(rowSums(matrix_bkk)),
                           step = 10,
                           tidy = TRUE)

rarecurve_tak <- rarecurve(matrix_tak,
                           sample = min(rowSums(matrix_tak)),
                           step = 10,
                           tidy = TRUE)

rarecurve_phat <- rarecurve(matrix_phat,
                           sample = min(rowSums(matrix_phat)),
                           step = 10,
                           tidy = TRUE)

p_bkk_rarecurv <- ggplot(rarecurve_bkk, aes(x = Sample, y = Species))+
  geom_point(colour = "#F8766D", size = 0.4, alpha = 0.7)+
  geom_vline(xintercept = min(rowSums(matrix_bkk)))+
  theme_bw()+
  labs(x = "Sequencing read", y = "ASV")+
  ylim(0, 300)+xlim(0, 62000)
  
p_tak_rarecurv <- ggplot(rarecurve_tak, aes(x = Sample, y = Species))+
  geom_point(colour = "#619CFF", size = 0.4, alpha = 0.7)+
  geom_vline(xintercept = min(rowSums(matrix_tak)))+
  theme_bw()+
  labs(x = "Sequencing read", y = "ASV")+
  ylim(0, 300)+xlim(0, 62000)

p_phat_rarecurv <- ggplot(rarecurve_phat, aes(x = Sample, y = Species))+
  geom_point(colour = "#00BA38", size = 0.4, alpha = 0.7)+
  geom_vline(xintercept = min(rowSums(matrix_phat)))+
  theme_bw()+
  labs(x = "Sequencing read", y = "ASV")+
  ylim(0, 300)+xlim(0, 62000)

ggsave(plot = p_bkk_rarecurv,
       filename = "p_bkk_rarecurv.png", 
       width = 3.5, height = 3.5, units = "in")

ggsave(plot = p_phat_rarecurv,
       filename = "p_phat_rarecurv.png", 
       width = 3.5, height = 3.5, units = "in")

ggsave(plot = p_tak_rarecurv,
       filename = "p_tak_rarecurv.png", 
       width = 3.5, height = 3.5, units = "in")