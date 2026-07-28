
# scatterplots characterizing the distribution of biomarker concentrations by diarrhea etiology

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# to prevent scientific notation on the plot
options(scipen=999)

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# load in ELISA and TAC data 
load("R output/EFGH_substudy_rowperpid_ALLDATA_LOCKEDv3.Rda")

# keep necessary variables (pid, country, biomarkers, pathogen_bin where pathogen_35 >5% of samples)
df <- EFGH_substudy_rowperpid_ALLDATA %>%
  select(pid, country, clean_conc_log10_CAL, clean_conc_log10_HAEM, clean_conc_log10_MPO, clean_conc_log10_NGAL, shigella_ct,
         adenovirus_40_41_ct, astrovirus_ct, c_jejuni_coli_ct, cryptosporidium_ct, norovirus_gii_ct, rotavirus_ct, sapovirus_ct, ST.ETEC_ct, tEPEC_ct)

# transform the original pathogen variable (continuous Ct value), to the log10 values (35-Ct)/3.322. 
# transforming the continuous Ct variable will make the coefficients more interpretable.
df1 <- df %>%
  mutate(shigella_log= (35-shigella_ct)/3.322) %>%
  mutate(adenovirus_40_41_log= (35-adenovirus_40_41_ct)/3.322) %>%
  mutate(astrovirus_log= (35-astrovirus_ct)/3.322) %>%
  mutate(c_jejuni_coli_log= (35-c_jejuni_coli_ct)/3.322) %>%
  mutate(cryptosporidium_log= (35-cryptosporidium_ct)/3.322) %>%
  mutate(norovirus_gii_log= (35-norovirus_gii_ct)/3.322) %>%
  mutate(rotavirus_log= (35-rotavirus_ct)/3.322) %>%
  mutate(sapovirus_log= (35-sapovirus_ct)/3.322) %>%
  mutate(ST_ETEC_log= (35-ST.ETEC_ct)/3.322) %>%
  mutate(tEPEC_log= (35-tEPEC_ct)/3.322)

shigella_df <- df1 %>%
  filter(shigella_log!=0)

adeno_df <- df1 %>%
  filter(adenovirus_40_41_log!=0)

astro_df <- df1 %>%
  filter(astrovirus_log!=0)

campy_df <- df1 %>%
  filter(c_jejuni_coli_log!=0)

crypto_df <- df1 %>%
  filter(cryptosporidium_log!=0)

norovirus_df <- df1 %>%
  filter(norovirus_gii_log!=0)

rotavirus_df <- df1 %>%
  filter(rotavirus_log!=0)

sapovirus_df <- df1 %>%
  filter(sapovirus_log!=0)

ST_ETEC_df <- df1 %>%
  filter(ST_ETEC_log!=0)

tEPEC_df <- df1 %>%
  filter(tEPEC_log!=0)

class(shigella_df$shigella_log)
class(shigella_df$clean_conc_log10_CAL)

# create function to loop through pathogens and biomarkers #

plot_function <- function(dat, pathogen, biomarker, biomarker_title, pathogen_title) {
  
  # Create the plot
  plot <- ggplot(dat, aes(x = !!rlang::sym(pathogen), y = !!rlang::sym(biomarker))) +
    theme_bw() +
    geom_point(size=0.5, alpha=0.5) + 
    geom_smooth(method = 'lm', formula = y ~ x) +
    labs(
      x = pathogen_title,
      y = biomarker_title) 
  plot = plot + theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank())
    
  plot = plot +  theme(axis.title=element_text(size=10))
  plot = plot +  theme(axis.text=element_text(size=10))
  plot = plot +  theme(axis.line = element_line(color="black"))
    
    return(plot)
}

# Call the function #

# shigella 
shigella_haem <- plot_function(shigella_df, "shigella_log", "clean_conc_log10_HAEM", "log10 concentration of hemoglobin", "log10 copies of Shigella")
shigella_haem

shigella_ngal <- plot_function(shigella_df, "shigella_log", "clean_conc_log10_NGAL", "log10 concentration of lipocalin-2", "log10 copies of Shigella")
shigella_ngal

shigella_mpo <- plot_function(shigella_df, "shigella_log", "clean_conc_log10_MPO", "log10 concentration of myeloperoxidase", "log10 copies of Shigella")
shigella_mpo

shigella_cal <- plot_function(shigella_df, "shigella_log", "clean_conc_log10_CAL", "log10 concentration of calprotectin", "log10 copies of Shigella")
shigella_cal

shigella <- ggarrange(shigella_haem, shigella_ngal, shigella_mpo, shigella_cal,   nrow=1, ncol=4)
shigella

# tEPEC 
tEPEC_haem<- plot_function(tEPEC_df, "tEPEC_log", "clean_conc_log10_HAEM", "log10 concentration of hemoglobin", "log10 copies of tEPEC")
tEPEC_haem

tEPEC_ngal <- plot_function(tEPEC_df, "tEPEC_log", "clean_conc_log10_NGAL", "log10 concentration of lipocalin-2", "log10 copies of tEPEC")
tEPEC_ngal

tEPEC_mpo <- plot_function(tEPEC_df, "tEPEC_log", "clean_conc_log10_MPO", "log10 concentration of myeloperoxidase", "log10 copies of tEPEC")
tEPEC_mpo

tEPEC_cal <- plot_function(tEPEC_df, "tEPEC_log", "clean_conc_log10_CAL", "log10 concentration of calprotectin", "log10 copies of tEPEC")
tEPEC_cal

tEPEC <- ggarrange(tEPEC_haem,  tEPEC_ngal,  tEPEC_mpo, tEPEC_cal,  nrow=1, ncol=4)
tEPEC

# campy 
campy_haem <- plot_function(campy_df, "c_jejuni_coli_log", "clean_conc_log10_HAEM", "log10 concentration of hemoglobin", "log10 copies of Campylobacter jejuni/coli")
campy_haem

campy_ngal <- plot_function(campy_df, "c_jejuni_coli_log", "clean_conc_log10_NGAL", "log10 concentration of lipocalin-2", "log10 copies of Campylobacter jejuni/coli")
campy_ngal

campy_mpo <- plot_function(campy_df, "c_jejuni_coli_log", "clean_conc_log10_MPO", "log10 concentration of myeloperoxidase", "log10 copies of Campylobacter jejuni/coli")
campy_mpo

campy_cal <- plot_function(campy_df, "c_jejuni_coli_log", "clean_conc_log10_CAL", "log10 concentration of calprotectin", "log10 copies of Campylobacter jejuni/coli")
campy_cal

campy <- ggarrange(campy_haem, campy_ngal, campy_mpo, campy_cal,   nrow=1, ncol=4)
campy

# ST-ETEC 
ST_ETEC_haem <- plot_function(ST_ETEC_df, "ST_ETEC_log", "clean_conc_log10_HAEM", "log10 concentration of hemoglobin", "log10 copies of ST-ETEC")
ST_ETEC_haem

ST_ETEC_ngal <- plot_function(ST_ETEC_df, "ST_ETEC_log", "clean_conc_log10_NGAL", "log10 concentration of lipocalin-2", "log10 copies of ST-ETEC")
ST_ETEC_ngal

ST_ETEC_mpo <- plot_function(ST_ETEC_df, "ST_ETEC_log", "clean_conc_log10_MPO", "log10 concentration of myeloperoxidase", "log10 copies of ST-ETEC")
ST_ETEC_mpo

ST_ETEC_cal <- plot_function(ST_ETEC_df, "ST_ETEC_log", "clean_conc_log10_CAL", "log10 concentration of calprotectin", "log10 copies of ST-ETEC")
ST_ETEC_cal

ST_ETEC <- ggarrange(ST_ETEC_haem,  ST_ETEC_ngal, ST_ETEC_mpo, ST_ETEC_cal,   nrow=1, ncol=4)
ST_ETEC

# crypto
crypto_haem <- plot_function(crypto_df, "cryptosporidium_log", "clean_conc_log10_HAEM", "log10 concentration of hemoglobin", "log10 copies of Cryptosporidium")
crypto_haem

crypto_ngal <- plot_function(crypto_df, "cryptosporidium_log", "clean_conc_log10_NGAL", "log10 concentration of lipocalin-2", "log10 copies of Cryptosporidium")
crypto_ngal

crypto_mpo <- plot_function(crypto_df, "cryptosporidium_log", "clean_conc_log10_MPO", "log10 concentration of myeloperoxidase", "log10 copies of Cryptosporidium")
crypto_mpo

crypto_cal <- plot_function(crypto_df, "cryptosporidium_log", "clean_conc_log10_CAL", "log10 concentration of calprotectin", "log10 copies of Cryptosporidium")
crypto_cal

crypto <- ggarrange( crypto_haem, crypto_ngal, crypto_mpo, crypto_cal, nrow=1, ncol=4)
crypto

#adenovirus 
adeno_haem <- plot_function(adeno_df, "adenovirus_40_41_log", "clean_conc_log10_HAEM", "log10 concentration of hemoglobin", "log10 copies of Adenovirus 40/41")
adeno_haem

adeno_ngal <- plot_function(adeno_df, "adenovirus_40_41_log", "clean_conc_log10_NGAL", "log10 concentration of lipocalin-2", "log10 copies of Adenovirus 40/41")
adeno_ngal

adeno_mpo <- plot_function(adeno_df, "adenovirus_40_41_log", "clean_conc_log10_MPO", "log10 concentration of myeloperoxidase", "log10 copies of Adenovirus 40/41")
adeno_mpo

adeno_cal <- plot_function(adeno_df, "adenovirus_40_41_log", "clean_conc_log10_CAL", "log10 concentration of calprotectin", "log10 copies of Adenovirus 40/41")
adeno_cal

adeno <- ggarrange(adeno_haem, adeno_ngal, adeno_mpo,  adeno_cal,  nrow=1, ncol=4)
adeno

# astrovirus
astro_haem <- plot_function(astro_df, "astrovirus_log", "clean_conc_log10_HAEM", "log10 concentration of hemoglobin", "log10 copies of Astrovirus")
astro_haem

astro_ngal <- plot_function(astro_df, "astrovirus_log", "clean_conc_log10_NGAL", "log10 concentration of lipocalin-2", "log10 copies of Astrovirus")
astro_ngal

astro_mpo <- plot_function(astro_df, "astrovirus_log", "clean_conc_log10_MPO", "log10 concentration of myeloperoxidase", "log10 copies of Astrovirus")
astro_mpo

astro_cal <- plot_function(astro_df, "astrovirus_log", "clean_conc_log10_CAL", "log10 concentration of calprotectin", "log10 copies of Astrovirus")
astro_cal

astro <- ggarrange(astro_haem, astro_ngal,  astro_mpo, astro_cal,  nrow=1, ncol=4)
astro

# norovirus 
norovirus_haem <- plot_function(norovirus_df, "norovirus_gii_log", "clean_conc_log10_HAEM", "log10 concentration of hemoglobin", "log10 copies of Norovirus GII")
norovirus_haem

norovirus_ngal <- plot_function(norovirus_df, "norovirus_gii_log", "clean_conc_log10_NGAL", "log10 concentration of lipocalin-2", "log10 copies of Norovirus GII")
norovirus_ngal

norovirus_mpo <- plot_function(norovirus_df, "norovirus_gii_log", "clean_conc_log10_MPO", "log10 concentration of myeloperoxidase", "log10 copies of Norovirus GII")
norovirus_mpo

norovirus_cal <- plot_function(norovirus_df, "norovirus_gii_log", "clean_conc_log10_CAL", "log10 concentration of calprotectin", "log10 copies of Norovirus GII")
norovirus_cal

norovirus <- ggarrange(norovirus_haem, norovirus_ngal, norovirus_mpo,  norovirus_cal,   nrow=1, ncol=4)
norovirus

#rotavirus
rotavirus_haem <- plot_function(rotavirus_df, "rotavirus_log", "clean_conc_log10_HAEM", "log10 concentration of hemoglobin", "log10 copies of Rotavirus")
rotavirus_haem

rotavirus_ngal <- plot_function(rotavirus_df, "rotavirus_log", "clean_conc_log10_NGAL", "log10 concentration of lipocalin-2", "log10 copies of Rotavirus")
rotavirus_ngal

rotavirus_mpo <- plot_function(rotavirus_df, "rotavirus_log", "clean_conc_log10_MPO", "log10 concentration of myeloperoxidase", "log10 copies of Rotavirus")
rotavirus_mpo

rotavirus_cal <- plot_function(rotavirus_df, "rotavirus_log", "clean_conc_log10_CAL", "log10 concentration of calprotectin", "log10 copies of Rotavirus")
rotavirus_cal

rotavirus <- ggarrange(rotavirus_haem, rotavirus_ngal, rotavirus_mpo, rotavirus_cal,    nrow=1, ncol=4)
rotavirus

#sapovirus
sapovirus_haem <- plot_function(sapovirus_df, "sapovirus_log", "clean_conc_log10_HAEM", "log10 concentration of hemoglobin", "log10 copies of Sapovirus")
sapovirus_haem

sapovirus_ngal <- plot_function(sapovirus_df, "sapovirus_log", "clean_conc_log10_NGAL", "log10 concentration of lipocalin-2", "log10 copies of Sapovirus")
sapovirus_ngal

sapovirus_mpo <- plot_function(sapovirus_df, "sapovirus_log", "clean_conc_log10_MPO", "log10 concentration of myeloperoxidase", "log10 copies of Sapovirus")
sapovirus_mpo

sapovirus_cal <- plot_function(sapovirus_df, "sapovirus_log", "clean_conc_log10_CAL", "log10 concentration of calprotectin", "log10 copies of Sapovirus")
sapovirus_cal

sapovirus <- ggarrange(sapovirus_haem, sapovirus_ngal, sapovirus_mpo, sapovirus_cal,  nrow=1, ncol=4)
sapovirus


all <- ggarrange(shigella, tEPEC, campy, ST_ETEC, crypto, adeno, astro, norovirus, rotavirus, sapovirus, nrow=10, ncol=1)
all


# Save Figure 
ggsave(file="Figures/Figure S6 - Biomarker distribution by etiology - Scatterplot - FINAL.pdf",all,width=13,height=30,dpi=300)

# split up figure 
firstfive <- ggarrange(shigella, tEPEC, campy, ST_ETEC, crypto, nrow=5, ncol=1)
lastfive <- ggarrange(adeno, astro, norovirus, rotavirus, sapovirus, nrow=5, ncol=1)

ggsave(file="Figures/Figure S6 - Biomarker distribution by etiology - Scatterplot - First5.pdf",firstfive,width=13,height=15,dpi=300)
ggsave(file="Figures/Figure S6 - Biomarker distribution by etiology - Scatterplot - Last5.pdf",lastfive,width=13,height=15,dpi=300)



