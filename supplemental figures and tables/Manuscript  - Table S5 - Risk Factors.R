
# risk factor regressions

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# to prevent scientific notation on the plot
options(scipen=999)

library(dplyr)
library(tidyr)
library(openxlsx)
library(gee)
library(geepack)

# load in elisa data
load("R output/EFGH_substudy_rowperpid_ALLDATA_LOCKEDv3.Rda")

# check NAs
table(EFGH_substudy_rowperpid_ALLDATA$clean_conc_log10_HAEM, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$clean_conc_log10_CAL, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$clean_conc_log10_MPO, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$clean_conc_log10_NGAL, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$dysentery, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$enr_age_months, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$enroll_diar_vom, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$enroll_diar_fever, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$enr_haz, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$enr_muac_before_cm, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$enr_waz, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$enr_whz, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$country, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$sex, useNA = "always")

table(EFGH_substudy_rowperpid_ALLDATA$tEPEC_attributable, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$c_jejuni_coli_attributable, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$ST.ETEC_attributable, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$cryptosporidium_attributable, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$adenovirus_40_41_attributable, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$astrovirus_attributable, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$norovirus_gii_attributable, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$rotavirus_attributable, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$sapovirus_attributable, useNA = "always")

#filter out NAs
small <-EFGH_substudy_rowperpid_ALLDATA %>%
  filter(!is.na(clean_conc_log10_HAEM)) %>%
  filter(!is.na(clean_conc_log10_CAL)) %>%
  filter(!is.na(clean_conc_log10_NGAL)) %>%
  filter(!is.na(clean_conc_log10_MPO)) %>%
  filter(!is.na(enr_age_months)) %>%
  filter(!is.na(sex)) %>%
  filter(!is.na(enroll_diar_fever)) %>%
  filter(!is.na(enroll_diar_vom)) %>%
  filter(!is.na(dysentery)) %>%
  filter(!is.na(country)) %>%
  filter(!is.na(enr_muac_before_cm)) %>%
  filter(!is.na(enr_haz)) %>%
  filter(!is.na(enr_waz)) %>%
  filter(!is.na(enr_whz)) %>%
  
  filter(!is.na(tEPEC_attributable)) %>%
  filter(!is.na(c_jejuni_coli_attributable)) %>%
  filter(!is.na(ST.ETEC_attributable)) %>%
  filter(!is.na(cryptosporidium_attributable)) %>%
  filter(!is.na(adenovirus_40_41_attributable)) %>%
  filter(!is.na(astrovirus_attributable)) %>%
  filter(!is.na(norovirus_gii_attributable)) %>%
  filter(!is.na(rotavirus_attributable)) %>%
  filter(!is.na(sapovirus_attributable)) 

# filter to where shigella is attributable
small <- small %>% filter(shigella_attributable==1)

###########################
# ONE MODEL PER BIOMARKER #
#        FINAL MODEL      #
###########################

# one model per biomarker with mutual adjustment
# do NOT include WAZ and WHZ in the same model as they're highly correlated - remove WAZ
risk_func <-function(biomarker, dat) {
cal_model <- geeglm(biomarker ~ sex + enr_age_months + enroll_diar_fever + enroll_diar_vom + dysentery + enr_muac_before_cm + enr_haz + enr_whz +
                      country + astrovirus_attributable + c_jejuni_coli_attributable + cryptosporidium_attributable + norovirus_gii_attributable + 
                      rotavirus_attributable + sapovirus_attributable + ST.ETEC_attributable + tEPEC_attributable, family = gaussian, data = dat, id = child_id, corstr = "exchangeable")

beta <- summary(cal_model)$coefficients[2:9,1]
lb <- summary(cal_model)$coefficients[2:9,1] - 1.96 * summary(cal_model)$coefficients[2:9,2]
ub <- summary(cal_model)$coefficients[2:9,1] + 1.96 * summary(cal_model)$coefficients[2:9,2]

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

CAL_model <- risk_func(small$clean_conc_log10_CAL, dat=small)
HAEM_model <- risk_func(small$clean_conc_log10_HAEM, dat=small)
MPO_model <- risk_func(small$clean_conc_log10_MPO, dat=small)
NGAL_model <- risk_func(small$clean_conc_log10_NGAL, dat=small)

biomarker_adj <- as.data.frame(cbind(HAEM_model, NGAL_model, MPO_model, CAL_model))

colnames(biomarker_adj) <- c("HAEM log10 Mean Difference (95% CI)", "HAEM Mean Difference (95% CI)", "NGAL log10 Mean Difference (95% CI)", "NGAL Mean Difference (95% CI)",
                          "MPO log10 Mean Difference (95% CI)", "MPO Mean Difference (95% CI)", "CAL log10 Mean Difference (95% CI)", "CAL Mean Difference (95% CI)")
rownames(biomarker_adj) <- c("Female", "Age (in months)", "Fever", "Vomiting", "Dysentery", "MUAC", "LAZ/HAZ", "WHZ")

biomarker_adj <- biomarker_adj %>%
  select(`HAEM log10 Mean Difference (95% CI)`, `NGAL log10 Mean Difference (95% CI)`, `MPO log10 Mean Difference (95% CI)`, `CAL log10 Mean Difference (95% CI)`,
         `HAEM Mean Difference (95% CI)`, `NGAL Mean Difference (95% CI)`, `MPO Mean Difference (95% CI)`, `CAL Mean Difference (95% CI)`)

# save table
write.csv(biomarker_adj, "Figures/Table S5 - Risk Factor Mean Differences - Shigella Only - Mutually Adjusted - FINAL.csv")





###########################################
# ONE MODEL PER BIOMARKER PER RISK FACTOR #
#     NOT MOVING FORWARD WITH THIS        #
###########################################

# create function to go through biomarkers and risk factors 
# dichomotous and continuous predictors ONLY (country cannot be used below)
risk_func <-function(biomarker, risk, dat) {
  model <- geeglm(biomarker ~ risk + country, family = gaussian, data = dat, id = child_id, corstr = "exchangeable")
  beta <- summary(model)$coefficients[2,1]
  lb <- summary(model)$coefficients[2,1] - 1.96 * summary(model)$coefficients[2,2]
  ub <- summary(model)$coefficients[2,1] + 1.96 * summary(model)$coefficients[2,2]
  est_CI<-paste0(round(beta,2)," (", round(lb,2),", ", round(ub,2),")")
  return(est_CI)
}


model <- geeglm(clean_conc_log10_HAEM ~ dysentery, family = gaussian, data = small, id = child_id, corstr = "exchangeable")
summary(model)

# call the function
# CAL
CAL_age <- risk_func(small$clean_conc_log10_CAL, small$enr_age_months, dat=small)
CAL_sex <- risk_func(small$clean_conc_log10_CAL, small$sex, dat=small)
CAL_fever <- risk_func(small$clean_conc_log10_CAL, small$enroll_diar_fever, dat=small)
CAL_vomit <- risk_func(small$clean_conc_log10_CAL, small$enroll_diar_vom, dat=small)
CAL_dysentery <- risk_func(small$clean_conc_log10_CAL, small$dysentery, dat=small)
CAL_muac <- risk_func(small$clean_conc_log10_CAL, small$enr_muac_before_cm, dat=small)
CAL_haz <- risk_func(small$clean_conc_log10_CAL, small$enr_haz, dat=small)
CAL_waz <- risk_func(small$clean_conc_log10_CAL, small$enr_waz, dat=small)
CAL_whz <- risk_func(small$clean_conc_log10_CAL, small$enr_whz, dat=small)

# HAEM
HAEM_age <- risk_func(small$clean_conc_log10_HAEM, small$enr_age_months, dat=small)
HAEM_sex <- risk_func(small$clean_conc_log10_HAEM, small$sex, dat=small)
HAEM_fever <- risk_func(small$clean_conc_log10_HAEM, small$enroll_diar_fever, dat=small)
HAEM_vomit <- risk_func(small$clean_conc_log10_HAEM, small$enroll_diar_vom, dat=small)
HAEM_dysentery <- risk_func(small$clean_conc_log10_HAEM, small$dysentery, dat=small)
HAEM_muac <- risk_func(small$clean_conc_log10_HAEM, small$enr_muac_before_cm, dat=small)
HAEM_haz <- risk_func(small$clean_conc_log10_HAEM, small$enr_haz, dat=small)
HAEM_waz <- risk_func(small$clean_conc_log10_HAEM, small$enr_waz, dat=small)
HAEM_whz <- risk_func(small$clean_conc_log10_HAEM, small$enr_whz, dat=small)

# MPO
MPO_age <- risk_func(small$clean_conc_log10_MPO, small$enr_age_months, dat=small)
MPO_sex <- risk_func(small$clean_conc_log10_MPO, small$sex, dat=small)
MPO_fever <- risk_func(small$clean_conc_log10_MPO, small$enroll_diar_fever, dat=small)
MPO_vomit <- risk_func(small$clean_conc_log10_MPO, small$enroll_diar_vom, dat=small)
MPO_dysentery <- risk_func(small$clean_conc_log10_MPO, small$dysentery, dat=small)
MPO_muac <- risk_func(small$clean_conc_log10_MPO, small$enr_muac_before_cm, dat=small)
MPO_haz <- risk_func(small$clean_conc_log10_MPO, small$enr_haz, dat=small)
MPO_waz <- risk_func(small$clean_conc_log10_MPO, small$enr_waz, dat=small)
MPO_whz <- risk_func(small$clean_conc_log10_MPO, small$enr_whz, dat=small)

# NGAL
NGAL_age <- risk_func(small$clean_conc_log10_NGAL, small$enr_age_months, dat=small)
NGAL_sex <- risk_func(small$clean_conc_log10_NGAL, small$sex, dat=small)
NGAL_fever <- risk_func(small$clean_conc_log10_NGAL, small$enroll_diar_fever, dat=small)
NGAL_vomit <- risk_func(small$clean_conc_log10_NGAL, small$enroll_diar_vom, dat=small)
NGAL_dysentery <- risk_func(small$clean_conc_log10_NGAL, small$dysentery, dat=small)
NGAL_muac <- risk_func(small$clean_conc_log10_NGAL, small$enr_muac_before_cm, dat=small)
NGAL_haz <- risk_func(small$clean_conc_log10_NGAL, small$enr_haz, dat=small)
NGAL_waz <- risk_func(small$clean_conc_log10_NGAL, small$enr_waz, dat=small)
NGAL_whz <- risk_func(small$clean_conc_log10_NGAL, small$enr_whz, dat=small)

CAL <- as.data.frame(rbind(CAL_sex, CAL_age, CAL_fever, CAL_vomit, CAL_dysentery, CAL_muac, CAL_haz, CAL_waz, CAL_whz))
HAEM <- as.data.frame(rbind(HAEM_sex, HAEM_age, HAEM_fever, HAEM_vomit, HAEM_dysentery, HAEM_muac, HAEM_haz, HAEM_waz, HAEM_whz))
MPO <- as.data.frame(rbind(MPO_sex, MPO_age, MPO_fever, MPO_vomit, MPO_dysentery, MPO_muac, MPO_haz, MPO_waz, MPO_whz))
NGAL <- as.data.frame(rbind(NGAL_sex, NGAL_age, NGAL_fever, NGAL_vomit, NGAL_dysentery, NGAL_muac, NGAL_haz, NGAL_waz, NGAL_whz))

# bind biomarkers together
biomarkers <- as.data.frame(cbind(HAEM, NGAL, MPO, CAL))

colnames(biomarkers) <- c("HAEM Mean Difference (95% CI)", "NGAL Mean Difference (95% CI)", "MPO Mean Difference (95% CI)", "CAL Mean Difference (95% CI)")
rownames(biomarkers) <- c("Female", "Age (in months)", "Fever", "Vomiting", "Dysentery", "MUAC", "LAZ/HAZ", "WAZ", "WHZ")

# save table
# write.csv(biomarkers, "Figures/Risk Factor Mean Differences - Shigella Only - FINAL.csv")




