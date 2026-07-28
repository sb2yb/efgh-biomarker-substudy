
#  boxplots characterizing the distribution of biomarker concentrations by diarrhea etiology
# this is a supplemental figure that has Shigella-attributable only, Shigella + co-bacterial attribution, Shigella + co-viral attribution, Shigella + Cryptosporidium attribution

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# to prevent scientific notation on the plot
options(scipen=999)

library(dplyr)
library(tidyr)
library(ggplot2)

# load in ELISA and TAC data 
load("R output/EFGH_substudy_rowperpid_ALLDATA_LOCKEDv3.Rda")


# keep necessary variables (pid, country, biomarkers, pathogen_attributable where pathogen_35 >5% of samples)
df <- EFGH_substudy_rowperpid_ALLDATA %>%
  select(pid, country, clean_conc_log10_CAL, clean_conc_log10_HAEM, clean_conc_log10_MPO, clean_conc_log10_NGAL, shigella_attributable,
         adenovirus_40_41_attributable, astrovirus_attributable, c_jejuni_coli_attributable, cryptosporidium_attributable, norovirus_gii_attributable, 
         rotavirus_attributable, sapovirus_attributable, ST.ETEC_attributable, tEPEC_attributable)

# Transform the dataframe so that there is only 1 column with the biomarker concentrations and 1 column indicating the name of the biomarker
long_df <- df %>%
  pivot_longer(
    cols = starts_with("clean_conc_log10_"),
    names_to = "biomarker",
    values_to = "log10_conc") %>% # Extract the biomarker names from the column names
  mutate(biomarker = gsub("clean_conc_log10_", "", biomarker))

# filter out where concentrations are NA
long_df <- long_df %>%
  filter(!is.na(log10_conc))



# Shigella only
shig_only <- long_df %>%
  filter(
    shigella_attributable == 1,
    c_jejuni_coli_attributable == 0,
    ST.ETEC_attributable == 0,
    tEPEC_attributable == 0,
    adenovirus_40_41_attributable == 0,
    astrovirus_attributable == 0,
    norovirus_gii_attributable == 0,
    rotavirus_attributable == 0,
    sapovirus_attributable == 0,
    cryptosporidium_attributable == 0) %>%
  mutate(pathogen = "Shigella only")

# Shigella + co-bacterial
shig_bacterial <- long_df %>%
  filter(
    shigella_attributable == 1,
    c_jejuni_coli_attributable == 1 |
      ST.ETEC_attributable == 1 |
      tEPEC_attributable == 1) %>%
  mutate(pathogen = "Shigella + co-bacterial")

# Shigella + co-viral
shig_viral <- long_df %>%
  filter(
    shigella_attributable == 1,
    adenovirus_40_41_attributable == 1 |
      astrovirus_attributable == 1 |
      norovirus_gii_attributable == 1 |
      rotavirus_attributable == 1 |
      sapovirus_attributable == 1) %>%
  mutate(pathogen = "Shigella + co-viral")

# Shigella + Cryptosporidium
shig_crypto <- long_df %>%
  filter(shigella_attributable == 1, cryptosporidium_attributable == 1) %>%
  mutate(pathogen = "Shigella + Cryptosporidium")

# Combine datasets
bind <- bind_rows(shig_only, shig_bacterial, shig_viral, shig_crypto)

# Set pathogen order
bind$pathogen <- factor(bind$pathogen,levels = c("Shigella only", "Shigella + co-bacterial", "Shigella + Cryptosporidium", "Shigella + co-viral"))

# set biomarker order
bind$biomarker <- factor(
  bind$biomarker,
  levels = c("HAEM", "NGAL", "MPO", "CAL"),
  labels = c("Hemoglobin", "Lipocalin-2", "Myeloperoxidase", "Calprotectin"))



# Figure 
plot = ggplot(bind, aes(x=biomarker, y=log10_conc, fill=pathogen)) 
plot = plot + theme_bw()
plot = plot + geom_boxplot() 
plot = plot + 
  labs(
    x = "",
    y = "Log10 (biomarker concentration)")
plot = plot + theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank())
plot = plot + scale_fill_manual(values=c("#AC1408", "#F8766D", "#FFC107", "#09B9FF"), 
                                name="",
                                breaks = c(
                                  "Shigella only",
                                  "Shigella + co-bacterial",
                                  "Shigella + Cryptosporidium",
                                  "Shigella + co-viral"))
plot = plot + theme(legend.position = c(.85, .28))
plot = plot +  theme(axis.title=element_text(size=16))
plot = plot +  theme(axis.text=element_text(size=14))
plot = plot +  theme(axis.line = element_line(color="black"))
plot

# Save Figure 
ggsave(file="Figures/Figure SX - Distribution of biomarkers by diarrhea etiology - shigella and co-attributuions - FINAL.pdf",plot,width=12,height=6,dpi=300)



