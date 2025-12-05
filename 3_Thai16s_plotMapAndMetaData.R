library(tidyverse)
library(sf)
library(rnaturalearth)
library(ggplot2)

## Get Thailand map --
cities <- data.frame(
  name = c("Bangkok", "Phatthalung", "Tak"),
  lon  = c(100.523186, 99.937307, 98.691235),
  lat  = c(13.7563, 7.614167, 16.878056)
)

th_provinces <- ne_states(country = "Thailand", returnclass = "sf")

Yesp_map <- ggplot() +
  geom_sf(data = th_provinces, fill = "grey90", color = "white") +  # provinces
  geom_sf(data = st_union(th_provinces), fill = NA, color = "black", size = 0.5) + # country outline
  geom_point(data = cities, aes(x = lon, y = lat), color = "steelblue", size = 1.5) +
  geom_text(data = cities, aes(x = lon, y = lat, label = name), vjust = -1, size = 2.5) +
  theme_void()

ggsave(plot = Yesp_map,
       filename = "C:/Project/5_16s_thai_population/figure/19102025_Yesp_map.png", 
       width = 9, height = 12, units = "cm")

## Plot meta data (main paper) --

meta_table <- read_tsv("C:/Project/5_16s_thai_population/thai_16s/metaData/metadata.tsv")
meta_table <- meta_table %>% column_to_rownames(var = "donor_id") 

#meta_table %>% group_by(locality) %>% summarise(age_mean = mean(age),
#                                                age_sd = sd(age))

#    locality    age_mean age_sd
#    <chr>          <dbl>  <dbl>
#  1 bangkok         30.8   6.22
#  2 phatthalung     37.3   9.78
#  3 tak             30.9   9.01

p_sex = meta_table %>% 
  filter(sex != "na") %>%
  ggplot(aes(x=locality, fill = sex))+
  geom_bar(position=position_dodge())+
  theme_bw()+
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))+
  labs(x = "", y = "Number individual", fill = "Sex")+
  theme(legend.position = "none", axis.text.x = element_text(angle = 315, hjust = 0))

p_bmi = meta_table %>% 
  filter(bmi != "na") %>% 
  ggplot(aes(x=locality, y=bmi, fill=locality))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none", axis.text.x = element_blank()) +
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))+
  labs(x = "", y = "BMI")

p_age = meta_table %>% 
  filter(age != "na") %>% 
  ggplot(aes(x=locality, y=age, fill=locality))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none", axis.text.x = element_blank()) +
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))+
  labs(x = "", y = "Age (year)")

p_temp = meta_table %>% 
  filter(AnnualMeanTemp != "na") %>%
  mutate(AnnualMeanTemp = AnnualMeanTemp) %>% 
  select(locality, AnnualMeanTemp) %>% unique() %>% 
  ggplot(aes(x=locality, y = AnnualMeanTemp, fill = locality))+
  geom_bar(stat = "identity")+
  theme_bw()+
  theme(legend.position = "none", axis.text.x = element_blank()) +
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))+
  labs(x = "", y = "Annual temperature (ºC)")

p_precip = meta_table %>% 
  filter(AnnualPrecipitation != "na") %>%
  mutate(AnnualPrecipitation = AnnualPrecipitation) %>% 
  select(locality, AnnualPrecipitation) %>% unique() %>% 
  ggplot(aes(x=locality, y = AnnualPrecipitation, fill = locality))+
  geom_bar(stat = "identity")+
  theme_bw()+
  theme(legend.position = "none", axis.text.x = element_blank()) +
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))+
  labs(x = "", y = "Annual precipitation (mm)")

p_locdens = meta_table %>% 
  filter(locality_density != "na") %>%
  select(locality, locality_density) %>% unique() %>% 
  ggplot(aes(x=locality, y = locality_density, fill = locality))+
  geom_bar(stat = "identity")+
  theme_bw()+
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))+
  labs(x = "", y = "Local density")+
  theme(legend.position = "none", axis.text.x = element_text(angle = 315, hjust = 0))

p_metadata_basic <- plot_grid(p_temp, p_age, p_precip, p_bmi,p_locdens, p_sex, ncol = 2, align = "v", rel_heights = c(1,1,1.3))

ggsave(plot = p_metadata_basic,
       filename = "C:/Project/5_16s_thai_population/figure/01092025_p_metadata_basic.png", 
       width = 2.8, height = 6.2, units = "in")

## Plot meta data (supplementary) --

p_weight = meta_table %>% 
  filter(weight_kg != "na") %>% 
  ggplot(aes(x=locality, y=weight_kg, fill=locality))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none", axis.text.x = element_blank()) +
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))+
  labs(x = "", y = "Weight (kg)")

p_height = meta_table %>% 
  filter(height_cm != "na") %>% 
  ggplot(aes(x=locality, y=height_cm, fill=locality))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none", axis.text.x = element_blank()) +
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))+
  labs(x = "", y = "Height (cm)")

p_water = meta_table %>% 
  filter(water_source != "na") %>%
  select(locality, water_source) %>%
  ggplot(aes(x=locality, fill = water_source))+
  geom_bar(position = position_dodge(preserve = "single"), width = 0.6)+
  theme_bw()+
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))+
  theme( axis.text.x = element_blank(), legend.position = "none") +
  labs(x = "", y = "Number", fill = "Water source")

p_house = meta_table %>% 
  filter(household_size != "na") %>%
  ggplot(aes(x=locality, y=household_size, fill=locality))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none", axis.text.x = element_blank()) +
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))+
  labs(x = "", y = "Household size")

p_exercise = meta_table %>% 
  filter(exercice_per_week != "na") %>%
  ggplot(aes(x=locality, y=exercice_per_week, fill=locality))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none", axis.text.x = element_blank()) +
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))+
  labs(x = "", y = "Exercise per week")

p_electricity = meta_table %>% 
  filter(access_to_electricity != "na") %>%
  mutate(access_to_electricity = as.character(access_to_electricity)) %>% 
  select(locality, access_to_electricity) %>%
  ggplot(aes(x=locality, fill = access_to_electricity))+
  geom_bar(position = position_dodge(preserve = "single"), width = 0.6)+
  theme_bw()+
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))+
  theme(legend.position = "none", axis.text.x = element_blank()) +
  labs(x = "", y = "Number", fill = "Access to\nelectricity")

p_csection <- meta_table %>% 
  filter(c_section != "na") %>%
  mutate(c_section = as.character(c_section)) %>% 
  select(locality, c_section) %>%
  ggplot(aes(x=locality, fill = c_section))+
  geom_bar(position = position_dodge(preserve = "single"), width = 0.6)+
  theme_bw()+
  theme(legend.position = "none", axis.text.x = element_blank()) +
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))+
  labs(x = "", y = "Number", fill = "C-section")

p_breastfed <- meta_table %>% 
  filter(breast_fed != "na") %>%
  mutate(breast_fed = as.character(breast_fed)) %>% 
  select(locality, breast_fed) %>%
  ggplot(aes(x=locality, fill = breast_fed))+
  geom_bar(position = position_dodge(preserve = "single"), width = 0.6)+
  theme_bw()+
  theme(legend.position = "none", axis.text.x = element_blank()) +
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))+
  labs(x = "", y = "Number", fill = "Breast\nfeeding")

p_antibiotic <- meta_table %>% 
  filter(antibiotic != "na") %>%
  mutate(antibiotic = as.character(antibiotic)) %>% 
  select(locality, antibiotic) %>%
  ggplot(aes(x=locality, fill = antibiotic))+
  geom_bar(position = position_dodge(preserve = "single"), width = 0.6)+
  theme_bw()+
  theme(legend.position = "none", axis.text.x = element_blank()) +
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))+
  labs(x = "", y = "Number", fill = "Antibiotic")

p_supprement <- meta_table %>% 
  filter(supplement_traditional_medicine != "na") %>%
  mutate(supplement_traditional_medicine = as.character(supplement_traditional_medicine)) %>% 
  select(locality, supplement_traditional_medicine) %>%
  ggplot(aes(x=locality, fill = supplement_traditional_medicine))+
  geom_bar(position = position_dodge(preserve = "single"), width = 0.6)+
  theme_bw()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 315, hjust = 0)) +
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))+
  labs(x = "", y = "Number", fill = "Traditional\nmedicine")

p_oral <- meta_table %>% 
  filter(oral_topical_medication != "na") %>%
  mutate(oral_topical_medication = as.character(oral_topical_medication)) %>% 
  select(locality, oral_topical_medication) %>%
  ggplot(aes(x=locality, fill = oral_topical_medication))+
  geom_bar(position = position_dodge(preserve = "single"), width = 0.6)+
  theme_bw()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 315, hjust = 0)) +
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))+
  labs(x = "", y = "Number", fill = "Oral\ntreatment")

p_laxative <- meta_table %>% 
  filter(laxative != "na") %>%
  mutate(laxative = as.character(laxative)) %>% 
  select(locality, laxative) %>%
  ggplot(aes(x=locality, fill = laxative))+
  geom_bar(position = position_dodge(preserve = "single"), width = 0.6)+
  theme_bw()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 315, hjust = 0)) +
  scale_x_discrete(labels = c("bangkok" = "Bangkok", "tak" = "Tak", "phatthalung" = "Phatthalung"))+
  labs(x = "", y = "Number", fill = "Laxative")

p_com_supple <- plot_grid(p_weight, p_height, p_exercise,
                          p_house, p_csection, p_breastfed,
                          p_electricity, p_water, p_antibiotic,
                          p_supprement, p_oral, p_laxative, ncol = 3, align = "XY", rel_heights = c(1,1,1,1.2))

ggsave(plot = p_com_supple,
       filename = "C:/Project/5_16s_thai_population/figure/01092025_p_com_supple.png", 
       width = 7, height = 10, units = "in")
