library(tidyverse)
library(phyloseq)
library(ggplot2)
library(cowplot)
library(ggvenn)
library(RColorBrewer)

alpha_diversity <- read.csv("C:/Project/5_16s_thai_population/22052005_diversityMetadataSampleDF.csv")

asv_table <- readRDS("C:/Project/5_16s_thai_population/seqtab_final.rds")
asv_table <- t(asv_table)
taxa_table <- readRDS("C:/Project/5_16s_thai_population/tax_final.rds")
meta_table <- read_tsv("C:/Project/5_16s_thai_population/thai_16s/metaData/metadata.tsv")
meta_table <- meta_table %>% column_to_rownames(var = "donor_id") 

# Convert to phyloseq components
ASV <- otu_table(asv_table, taxa_are_rows = TRUE)
TAX <- tax_table(taxa_table)
SAMP <- sample_data(meta_table)
physeq <- phyloseq(ASV, TAX, SAMP)

## Composition ----

## Plot relative abundance (filtered) ----
taxonomic_df_all <- psmelt(physeq)
#length(unique(taxonomic_df_all$Genus))

unique_genus_bkk <- taxonomic_df_all %>% filter(locality=="bangkok") %>% filter(Abundance!=0) %>% pull(Genus)
unique_genus_phat <- taxonomic_df_all %>% filter(locality=="phatthalung") %>% filter(Abundance!=0) %>% pull(Genus)
unique_genus_tak <- taxonomic_df_all %>% filter(locality=="tak") %>% filter(Abundance!=0) %>% pull(Genus)

unique_genus_df <- data.frame()
unique_genus_df_bkk <- taxonomic_df_all %>% filter(locality=="bangkok") %>% filter(Abundance!=0) %>% select(Genus,locality) %>% unique()
unique_genus_df_phat <- taxonomic_df_all %>% filter(locality=="phatthalung") %>% filter(Abundance!=0) %>% select(Genus,locality) %>% unique()
unique_genus_df_tak <- taxonomic_df_all %>% filter(locality=="tak") %>% filter(Abundance!=0) %>% select(Genus,locality) %>% unique()

unique_genus_df <- rbind(unique_genus_df, unique_genus_df_tak, unique_genus_df_phat, unique_genus_df_bkk)
rare_genus <- unique_genus_df %>% group_by(Genus) %>% summarise(count = n()) %>% filter( count == 1)

rare_genus_df <- unique_genus_df %>% filter(Genus %in% rare_genus$Genus) 

write.csv(x = rare_genus_df, file = "C:/Project/5_16s_thai_population/24072025_rare_genus.csv")

venn_list_all <- list("Bangkok"=unique_genus_bkk,
                      "Phatthalung"=unique_genus_phat,
                      "Tak"=unique_genus_tak)

p_venn_afterfiltered<-ggvenn(
  venn_list_all,
  fill_color = c("#F8766D", "#00BA38", "#619CFF"), # Pastel colors
  stroke_color = "black",        # Remove border lines
  #set_name_size = 2,        # Adjust set label size
  #text_size = 2,             # Adjust text inside the regions
  show_percentage = FALSE
)

ggsave(plot = p_venn_afterfiltered,
       filename = "C:/Project/5_16s_thai_population/figure/03092025_venn_afterfiltered_genus.png",
       width = 3.5, height = 3.5)

cute_colors <- colorRampPalette(brewer.pal(12, "Paired"))(82)
p_bar_family_phyloseq<-taxonomic_df_all %>% 
  group_by(locality) %>%
  mutate(avg_ra = Abundance/sum(Abundance)) %>%
  ggplot(aes(x=locality, y = avg_ra, fill = Family))+
  geom_bar(stat = "identity")+
  theme_bw() +
  theme(
    #legend.text = element_text(size = 7, , family = "Arial"),
    legend.position = "none",
    #legend.key.size = unit(0.4, "cm"),
    #text = element_text(size = 9, family = "Arial"),
    axis.text.x = element_text(angle = 315, hjust = 0)
  ) +
  scale_fill_manual(values = cute_colors) +  # Use custom colors
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung")) +
  labs(
    x = "",
    y = ""
  )
#guides(fill = guide_legend(nrow = 30))

p_bar_family_phyloseq_axis<-taxonomic_df_all %>% 
  group_by(locality) %>%
  mutate(avg_ra = Abundance/sum(Abundance)) %>%
  ggplot(aes(x=locality, y = avg_ra, fill = Family))+
  geom_bar(stat = "identity")+
  theme_bw() +
  theme(
    #legend.text = element_text(size = 7, , family = "Arial"),
    legend.position = "bottom",
    #legend.key.size = unit(0.4, "cm"),
    #text = element_text(size = 9, family = "Arial"),
    axis.text.x = element_text(angle = 315, hjust = 0)
  ) +
  scale_fill_manual(values = cute_colors) +  # Use custom colors
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung")) +
  labs(
    x = "",
    y = "Relative Abundance"
  )+ guides(fill = guide_legend(nrow = 10))

ggsave(plot = p_bar_family_phyloseq,
       filename = "C:/Project/5_16s_thai_population/figure/28032025_p_bar_family_phyloseq.png", 
       width = 17, height = 15, units = "cm")

ggsave(plot = p_bar_family_phyloseq_axis,
       filename = "C:/Project/5_16s_thai_population/figure/28072025_p_bar_family_phyloseq_axis.png", 
       width = 45, height = 15, units = "cm")

p_bar_order_phyloseq<-taxonomic_df_all %>% 
  group_by(locality) %>%
  mutate(avg_ra = Abundance/sum(Abundance)) %>%
  ggplot(aes(x=locality, y = avg_ra, fill = Order))+
  geom_bar(stat = "identity")+
  theme_classic() +
  theme_classic() +
  theme(
    legend.text = element_text(size = 7, , family = "Arial"),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    text = element_text(size = 9, family = "Arial"),
    axis.text.x = element_text(angle = 315, hjust = 0)
  ) +
  scale_fill_manual(values = cute_colors) +  # Use custom colors
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung")) +
  labs(
    x = "",
    y = "Relative Abundance"
  )+
  guides(fill = guide_legend(nrow = 30))

ggsave(plot = p_bar_order_phyloseq,
       filename = "C:/Project/5_16s_thai_population/figure/17022025_p_bar_order_phyloseq.png", 
       width = 10, height = 8)

cute_colors <- colorRampPalette(brewer.pal(12, "Paired"))(16)

p_bar_phylum_phyloseq<-taxonomic_df_all %>% 
  group_by(locality) %>%
  mutate(avg_ra = Abundance/sum(Abundance)) %>%
  ggplot(aes(x=locality, y = avg_ra, fill = Phylum))+
  geom_bar(stat = "identity")+
  theme_bw() +
  theme(
    #legend.text = element_text(size = 7, , family = "Arial"),
    legend.position = "none",
    #legend.key.size = unit(0.4, "cm"),
    #text = element_text(size = 9, family = "Arial"),
    axis.text.x = element_text(angle = 315, hjust = 0)
  ) +
  scale_fill_manual(values = cute_colors) +  # Use custom colors
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung")) +
  labs(
    x = "",
    y = "Relative Abundance"
  )
#guides(fill = guide_legend(nrow = 30))

p_bar_phylum_phyloseq_axis<-taxonomic_df_all %>% 
  group_by(locality) %>%
  mutate(avg_ra = Abundance/sum(Abundance)) %>%
  ggplot(aes(x=locality, y = avg_ra, fill = Phylum))+
  geom_bar(stat = "identity")+
  theme_bw() +
  theme(
    #legend.text = element_text(size = 7, , family = "Arial"),
    legend.position = "bottom",
    #legend.key.size = unit(0.4, "cm"),
    #text = element_text(size = 9, family = "Arial"),
    axis.text.x = element_text(angle = 315, hjust = 0)
  ) +
  scale_fill_manual(values = cute_colors) +  # Use custom colors
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung")) +
  labs(
    x = "",
    y = "Relative Abundance"
  )+ guides(fill = guide_legend(nrow = 2))

ggsave(plot = p_bar_phylum_phyloseq,
       filename = "C:/Project/5_16s_thai_population/figure/25032025_p_bar_phylum_phyloseq.png", 
       width = 7, height = 10, units = "cm")
ggsave(plot = p_bar_phylum_phyloseq_axis,
       filename = "C:/Project/5_16s_thai_population/figure/28072025_p_bar_phylum_phyloseq_axis.png", 
       width = 30, height = 15, units = "cm")
p_comb_bar_family_phylum_axis <- plot_grid(p_bar_phylum_phyloseq_axis, p_bar_family_phyloseq_axis, ncol = 2, align = "xy")
ggsave(plot = p_comb_bar_family_phylum_axis,
       filename = "C:/Project/5_16s_thai_population/figure/28072025_p_comb_bar_family_phylum_axis.png", 
       width =90, height = 22, units = "cm")
p_comb_bar_family_phylum <- plot_grid(p_bar_phylum_phyloseq, p_bar_family_phyloseq, ncol = 2, align = "y")
ggsave(plot = p_comb_bar_family_phylum,
       filename = "C:/Project/5_16s_thai_population/figure/28072025_p_comb_bar_family_phylum.png", 
       width = 9, height = 10, units = "cm")

## Alpha PD ----

#Higher Faith’s PD → More evolutionary history is preserved in the sample.
#Lower Faith’s PD → The species in the community are closely related (evolutionarily constrained).
library(rstatix)

stat_faith <- alpha_diversity %>% 
  pairwise_wilcox_test(PD ~ locality, p.adjust.method = "BH")

#.y.   group1      group2         n1    n2 statistic            p       p.adj p.adj.signif
#* <chr> <chr>       <chr>       <int> <int>     <dbl>        <dbl>       <dbl> <chr>       
#  1 PD    bangkok     phatthalung    47    29       200 0.0000000426 0.000000128 ****        
#  2 PD    bangkok     tak            47    30       446 0.006        0.006       **          
#  3 PD    phatthalung tak            29    30       676 0.000172     0.000258    *** 

p_faith <- alpha_diversity %>%  
  ggplot(aes(x = locality, y=PD, fill = locality))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 315, hjust = 0))+
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung")) +
  labs(
    x = "",
    y = "Faith's PD"
  )

print(p_faith)

#PSV ≈ 1 → Species are spread across deep evolutionary branches.
#PSV ≈ 0 → Species are closely related.

stat_psv <- alpha_diversity %>% 
  pairwise_wilcox_test(PSVs ~ locality, p.adjust.method = "BH")

#.y.   group1      group2         n1    n2 statistic         p    p.adj p.adj.signif
#* <chr> <chr>       <chr>       <int> <int>     <dbl>     <dbl>    <dbl> <chr>       
#  1 PSVs  bangkok     phatthalung    47    29       429 0.007     0.01     **          
#  2 PSVs  bangkok     tak            47    30       336 0.0000758 0.000227 ***         
#  3 PSVs  phatthalung tak            29    30       345 0.176     0.176    ns

p_psv <- alpha_diversity %>% 
  ggplot(aes(x = locality, y=PSVs, fill=locality))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 315, hjust = 0))+
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung")) +
  labs(
    x = "",
    y = "PSV"
  )


## Alpha diversity ----

richness_plot <- ggplot(alpha_diversity, aes(x = locality, y = Richness, fill = locality)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "", y = "Richness") +
  theme(legend.position = "none", axis.text.x = element_blank())+
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))

shannon_plot <- ggplot(alpha_diversity, aes(x = locality, y = Shannon, fill = locality)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "", y = "Shannon Index") +
  theme(legend.position = "none", axis.text.x = element_blank())+
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))

simpson_plot <- ggplot(alpha_diversity, aes(x = locality, y = Simpson, fill = locality)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "", y = "Simpson Index") +
  theme(legend.position = "none", axis.text.x = element_blank())+
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))

evenness_plot <- ggplot(alpha_diversity, aes(x = locality, y = Evenness, fill = locality)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "", y = "Evenness") +
  theme(legend.position = "none", axis.text.x = element_blank())+ #, axis.text.x = element_blank()
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))

p_alpha <- plot_grid(richness_plot, evenness_plot, shannon_plot, simpson_plot, p_psv,p_faith, 
                     ncol = 2, align = "v",rel_heights = c(3,3,4))

ggsave(plot = p_alpha,
       filename = "C:/Project/5_16s_thai_population/figure/28072025_alpha_diversity.png", 
       width = 2.7, height = 6, units = "in")


stat_richness <- alpha_diversity %>% 
  pairwise_wilcox_test(Richness ~ locality, p.adjust.method = "BH")

#.y.      group1      group2         n1    n2 statistic          p      p.adj p.adj.signif
#* <chr>    <chr>       <chr>       <int> <int>     <dbl>      <dbl>      <dbl> <chr>       
#  1 Richness bangkok     phatthalung    47    29      236  0.00000195 0.00000585 ****        
#  2 Richness bangkok     tak            47    30      538. 0.083      0.083      ns          
#  3 Richness phatthalung tak            29    30      708  0.000036   0.000054   ****   

stat_evenness <- alpha_diversity %>% 
  pairwise_wilcox_test(Evenness ~ locality, p.adjust.method = "BH")

#.y.      group1      group2         n1    n2 statistic         p     p.adj p.adj.signif
#* <chr>    <chr>       <chr>       <int> <int>     <dbl>     <dbl>     <dbl> <chr>       
#  1 Evenness bangkok     phatthalung    47    29       789 0.254     0.254     ns          
#  2 Evenness bangkok     tak            47    30      1096 0.0000251 0.0000753 ****        
#  3 Evenness phatthalung tak            29    30       632 0.002     0.004     **   

stat_shannon <- alpha_diversity %>% 
  pairwise_wilcox_test(Shannon ~ locality, p.adjust.method = "BH")

#.y.     group1      group2         n1    n2 statistic        p p.adj p.adj.signif
#* <chr>   <chr>       <chr>       <int> <int>     <dbl>    <dbl> <dbl> <chr>       
#  1 Shannon bangkok     phatthalung    47    29       474 0.026    0.037 *           
#  2 Shannon bangkok     tak            47    30       905 0.037    0.037 *           
#  3 Shannon phatthalung tak            29    30       663 0.000402 0.001 **

stat_simpson <- alpha_diversity %>% 
  pairwise_wilcox_test(Simpson ~ locality, p.adjust.method = "BH")

#.y.     group1      group2         n1    n2 statistic     p p.adj p.adj.signif
#* <chr>   <chr>       <chr>       <int> <int>     <dbl> <dbl> <dbl> <chr>       
#  1 Simpson bangkok     phatthalung    47    29       665 0.865 0.865 ns          
#  2 Simpson bangkok     tak            47    30       998 0.002 0.004 **          
#  3 Simpson phatthalung tak            29    30       629 0.003 0.004 **   
