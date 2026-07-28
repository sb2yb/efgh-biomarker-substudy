
# linear regression analysis where the pathogen Ct value is the predictor, and the log10 biomarker concentration is the outcome. 
# Run the models with and without a quadratic term (I(x^2)). 

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# to prevent scientific notation on the plot
options(scipen=999)

library(dplyr)
library(tidyr)
library(gee)
library(geepack)

# load in elisa data
load("R output/EFGH_substudy_rowperpid_ALLDATA_LOCKEDv3.Rda")

# check NAs
table(EFGH_substudy_rowperpid_ALLDATA$clean_conc_log10_HAEM, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$clean_conc_log10_CAL, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$clean_conc_log10_MPO, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$clean_conc_log10_NGAL, useNA = "always")

table(EFGH_substudy_rowperpid_ALLDATA$shigella_ct, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$adenovirus_40_41_ct, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$astrovirus_ct, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$c_jejuni_coli_ct, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$cryptosporidium_ct, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$norovirus_gii_ct, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$rotavirus_ct, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$sapovirus_ct, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$ST.ETEC_ct, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$tEPEC_ct, useNA = "always")

table(EFGH_substudy_rowperpid_ALLDATA$country, useNA = "always")
class(EFGH_substudy_rowperpid_ALLDATA$child_id)

# keep necessary variables (pid, country, biomarkers, pathogen_bin where pathogen_35 >5% of samples)
df <- EFGH_substudy_rowperpid_ALLDATA %>%
  select(pid, child_id, country, clean_conc_log10_CAL, clean_conc_log10_HAEM, clean_conc_log10_MPO, clean_conc_log10_NGAL, shigella_ct,
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

#filter out NAs and make individual datasets
shigella_df <- df1 %>%
  filter(!is.na(shigella_log)) %>%
  filter(shigella_log!=0)

adeno_df <- df1 %>%
  filter(!is.na(adenovirus_40_41_log)) %>%
  filter(adenovirus_40_41_log!=0)

astrovirus_df <- df1 %>%
  filter(!is.na(astrovirus_log)) %>%
  filter(astrovirus_log!=0)

campy_df <- df1 %>%
  filter(!is.na(c_jejuni_coli_log)) %>%
  filter(c_jejuni_coli_log!=0)

crypto_df <- df1 %>%
  filter(!is.na(cryptosporidium_log)) %>%
  filter(cryptosporidium_log!=0)

norovirus_df <- df1 %>%
  filter(!is.na(norovirus_gii_log)) %>%
  filter(norovirus_gii_log!=0)

rotavirus_df <- df1 %>%
  filter(!is.na(rotavirus_log)) %>%
  filter(rotavirus_log!=0)

sapovirus_df <- df1 %>%
  filter(!is.na(sapovirus_log)) %>%
  filter(sapovirus_log!=0)

ST_ETEC_df <- df1 %>%
  filter(!is.na(ST_ETEC_log)) %>%
  filter(ST_ETEC_log!=0)

tEPEC_df <- df1 %>%
  filter(!is.na(tEPEC_log)) %>%
  filter(tEPEC_log!=0)

# function that will keep the log10 and antilog values that were run through the model 
risk_func <-function(biomarker, pathogen, dat) {
  model <- geeglm(biomarker ~ pathogen  + country, family = gaussian, data = dat, id = child_id, corstr = "exchangeable")
  beta <- summary(model)$coefficients[2,1]
  lb <- summary(model)$coefficients[2,1] - 1.96 * summary(model)$coefficients[2,2]
  ub <- summary(model)$coefficients[2,1] + 1.96 * summary(model)$coefficients[2,2]
  
  # Compute and format antilog values (base 10)
  beta_antilog <- 10^beta
  lb_antilog <- 10^lb
  ub_antilog <- 10^ub
  
  # Format both versions
  formatted_log10 <- paste0(
    sprintf("%.2f", beta), " (", sprintf("%.2f", lb), ", ", sprintf("%.2f", ub), ")")
  formatted_antilog <- paste0(
    sprintf("%.2f", beta_antilog), " (", sprintf("%.2f", lb_antilog), ", ", sprintf("%.2f", ub_antilog), ")")
  
  # Return everything in one row
  result <- data.frame(
    formatted_log10 = formatted_log10,
    formatted_antilog = formatted_antilog)
  
  return(result)
}

model <- geeglm(clean_conc_log10_HAEM ~ shigella_log  + country, family = gaussian, data = shigella_df, id = child_id, corstr = "exchangeable")
beta <- summary(model)$coefficients[2,1]
lb <- summary(model)$coefficients[2,1] - 1.96 * summary(model)$coefficients[2,2]
ub <- summary(model)$coefficients[2,1] + 1.96 * summary(model)$coefficients[2,2]

summary(model)


# call the function
# CAL
CAL_shigella <- risk_func(shigella_df$clean_conc_log10_CAL, shigella_df$shigella_log, dat=shigella_df)
CAL_tEPEC <- risk_func(tEPEC_df$clean_conc_log10_CAL, tEPEC_df$tEPEC_log, dat=tEPEC_df)
CAL_campy <- risk_func(campy_df$clean_conc_log10_CAL, campy_df$c_jejuni_coli_log, dat=campy_df)
CAL_ST_ETEC <- risk_func(ST_ETEC_df$clean_conc_log10_CAL, ST_ETEC_df$ST_ETEC_log, dat=ST_ETEC_df)
CAL_crypto <- risk_func(crypto_df$clean_conc_log10_CAL, crypto_df$cryptosporidium_log, dat=crypto_df)
CAL_adeno <- risk_func(adeno_df$clean_conc_log10_CAL, adeno_df$adenovirus_40_41_log, dat=adeno_df)
CAL_astrovirus <- risk_func(astrovirus_df$clean_conc_log10_CAL, astrovirus_df$astrovirus_log, dat=astrovirus_df)
CAL_norovirus <- risk_func(norovirus_df$clean_conc_log10_CAL, norovirus_df$norovirus_gii_log, dat=norovirus_df)
CAL_rotavirus <- risk_func(rotavirus_df$clean_conc_log10_CAL, rotavirus_df$rotavirus_log, dat=rotavirus_df)
CAL_sapovirus <- risk_func(sapovirus_df$clean_conc_log10_CAL, sapovirus_df$sapovirus_log, dat=sapovirus_df)

# HAEM
HAEM_shigella <- risk_func(shigella_df$clean_conc_log10_HAEM, shigella_df$shigella_log, dat=shigella_df)
HAEM_tEPEC <- risk_func(tEPEC_df$clean_conc_log10_HAEM, tEPEC_df$tEPEC_log, dat=tEPEC_df)
HAEM_campy <- risk_func(campy_df$clean_conc_log10_HAEM, campy_df$c_jejuni_coli_log, dat=campy_df)
HAEM_ST_ETEC <- risk_func(ST_ETEC_df$clean_conc_log10_HAEM, ST_ETEC_df$ST_ETEC_log, dat=ST_ETEC_df)
HAEM_crypto <- risk_func(crypto_df$clean_conc_log10_HAEM, crypto_df$cryptosporidium_log, dat=crypto_df)
HAEM_adeno <- risk_func(adeno_df$clean_conc_log10_HAEM, adeno_df$adenovirus_40_41_log, dat=adeno_df)
HAEM_astrovirus <- risk_func(astrovirus_df$clean_conc_log10_HAEM, astrovirus_df$astrovirus_log, dat=astrovirus_df)
HAEM_norovirus <- risk_func(norovirus_df$clean_conc_log10_HAEM, norovirus_df$norovirus_gii_log, dat=norovirus_df)
HAEM_rotavirus <- risk_func(rotavirus_df$clean_conc_log10_HAEM, rotavirus_df$rotavirus_log, dat=rotavirus_df)
HAEM_sapovirus <- risk_func(sapovirus_df$clean_conc_log10_HAEM, sapovirus_df$sapovirus_log, dat=sapovirus_df)

# MPO
MPO_shigella <- risk_func(shigella_df$clean_conc_log10_MPO, shigella_df$shigella_log, dat=shigella_df)
MPO_tEPEC <- risk_func(tEPEC_df$clean_conc_log10_MPO, tEPEC_df$tEPEC_log, dat=tEPEC_df)
MPO_campy <- risk_func(campy_df$clean_conc_log10_MPO, campy_df$c_jejuni_coli_log, dat=campy_df)
MPO_ST_ETEC <- risk_func(ST_ETEC_df$clean_conc_log10_MPO, ST_ETEC_df$ST_ETEC_log, dat=ST_ETEC_df)
MPO_crypto <- risk_func(crypto_df$clean_conc_log10_MPO, crypto_df$cryptosporidium_log, dat=crypto_df)
MPO_adeno <- risk_func(adeno_df$clean_conc_log10_MPO, adeno_df$adenovirus_40_41_log, dat=adeno_df)
MPO_astrovirus <- risk_func(astrovirus_df$clean_conc_log10_MPO, astrovirus_df$astrovirus_log, dat=astrovirus_df)
MPO_norovirus <- risk_func(norovirus_df$clean_conc_log10_MPO, norovirus_df$norovirus_gii_log, dat=norovirus_df)
MPO_rotavirus <- risk_func(rotavirus_df$clean_conc_log10_MPO, rotavirus_df$rotavirus_log, dat=rotavirus_df)
MPO_sapovirus <- risk_func(sapovirus_df$clean_conc_log10_MPO, sapovirus_df$sapovirus_log, dat=sapovirus_df)

# NGAL
NGAL_shigella <- risk_func(shigella_df$clean_conc_log10_NGAL, shigella_df$shigella_log, dat=shigella_df)
NGAL_tEPEC <- risk_func(tEPEC_df$clean_conc_log10_NGAL, tEPEC_df$tEPEC_log, dat=tEPEC_df)
NGAL_campy <- risk_func(campy_df$clean_conc_log10_NGAL, campy_df$c_jejuni_coli_log, dat=campy_df)
NGAL_ST_ETEC <- risk_func(ST_ETEC_df$clean_conc_log10_NGAL, ST_ETEC_df$ST_ETEC_log, dat=ST_ETEC_df)
NGAL_crypto <- risk_func(crypto_df$clean_conc_log10_NGAL, crypto_df$cryptosporidium_log, dat=crypto_df)
NGAL_adeno <- risk_func(adeno_df$clean_conc_log10_NGAL, adeno_df$adenovirus_40_41_log, dat=adeno_df)
NGAL_astrovirus <- risk_func(astrovirus_df$clean_conc_log10_NGAL, astrovirus_df$astrovirus_log, dat=astrovirus_df)
NGAL_norovirus <- risk_func(norovirus_df$clean_conc_log10_NGAL, norovirus_df$norovirus_gii_log, dat=norovirus_df)
NGAL_rotavirus <- risk_func(rotavirus_df$clean_conc_log10_NGAL, rotavirus_df$rotavirus_log, dat=rotavirus_df)
NGAL_sapovirus <- risk_func(sapovirus_df$clean_conc_log10_NGAL, sapovirus_df$sapovirus_log, dat=sapovirus_df)


# binding together 
CAL <- as.data.frame(rbind(CAL_shigella, CAL_tEPEC, CAL_campy, CAL_ST_ETEC, CAL_crypto, CAL_adeno, CAL_astrovirus, CAL_norovirus, CAL_rotavirus, CAL_sapovirus))
HAEM <- as.data.frame(rbind(HAEM_shigella, HAEM_tEPEC, HAEM_campy, HAEM_ST_ETEC, HAEM_crypto, HAEM_adeno, HAEM_astrovirus, HAEM_norovirus, HAEM_rotavirus, HAEM_sapovirus))
MPO <- as.data.frame(rbind(MPO_shigella, MPO_tEPEC, MPO_campy, MPO_ST_ETEC, MPO_crypto, MPO_adeno, MPO_astrovirus, MPO_norovirus, MPO_rotavirus, MPO_sapovirus))
NGAL <- as.data.frame(rbind(NGAL_shigella, NGAL_tEPEC, NGAL_campy, NGAL_ST_ETEC, NGAL_crypto, NGAL_adeno, NGAL_astrovirus, NGAL_norovirus, NGAL_rotavirus, NGAL_sapovirus))

# bind biomarkers together
biomarkers <- as.data.frame(cbind(HAEM, NGAL, MPO, CAL))

colnames(biomarkers) <- c("HAEM log10 Mean Difference (95% CI)", "HAEM Mean Difference (95% CI)", "NGAL log10 Mean Difference (95% CI)", "NGAL Mean Difference (95% CI)",
                          "MPO log10 Mean Difference (95% CI)", "MPO Mean Difference (95% CI)", "CAL log10 Mean Difference (95% CI)", "CAL Mean Difference (95% CI)")
rownames(biomarkers) <- c("Shigella", "tEPEC", "Campylobacter", "ST-ETEC", "Cryptosporidium", "Adenovirus 40/41", "Astrovirus", "Norovirus GII", "Rotavirus", "Sapovirus")

biomarkers <- biomarkers %>%
  select(`HAEM log10 Mean Difference (95% CI)`, `NGAL log10 Mean Difference (95% CI)`, `MPO log10 Mean Difference (95% CI)`, `CAL log10 Mean Difference (95% CI)`,
         `HAEM Mean Difference (95% CI)`, `NGAL Mean Difference (95% CI)`, `MPO Mean Difference (95% CI)`, `CAL Mean Difference (95% CI)`)

# save table
write.csv(biomarkers, "Figures/Table 1 - Pathogen and Biomarker Mean Differences - FINAL.csv")



