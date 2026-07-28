
# sensitivity analysis requested by patty: 
# step 1: run geeglm with dysentery as my outcome and and scr_diar_days as my predictor
# step 2: run 3 ROC curves (group days 1-2, 3, and 4+ ) with log10 HAEM and dysentery

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

## ROC curves
library(pROC)
library(cutpointr) # The following objects are masked from ‘package:pROC’: auc, roc

# load in elisa data
load("R output/EFGH_substudy_rowperpid_ALLDATA_LOCKEDv3.Rda")

# keep necessary variables (pid, country, biomarkers, shigella_attributable)
small <- EFGH_substudy_rowperpid_ALLDATA %>%
  select(pid, country, child_id, clean_conc_log10_CAL, clean_conc_log10_HAEM, clean_conc_log10_MPO, clean_conc_log10_NGAL, dysentery, scr_diar_days)


# mean, SD, and counts for each biomarker 
# it is imperative that the mean and sd are calculate as they are below 
# the below ensures that 1) we are using the GEOMETIC mean when antilogging and 2) the entire logged values (and not rounded versions) are antilogged 

# List of biomarkers
biomarkers <- c("clean_conc_log10_CAL", 
                "clean_conc_log10_HAEM", 
                "clean_conc_log10_MPO", 
                "clean_conc_log10_NGAL")

# Function to summarize one biomarker (no rounding)
summarize_biomarker <- function(df, var) {
  mean_log10 <- mean(df[[var]], na.rm = TRUE)
  sd_log10   <- sd(df[[var]], na.rm = TRUE)
  geom_mean  <- 10^mean_log10
  geom_sd    <- 10^sd_log10
  n          <- sum(!is.na(df[[var]]))
  
  tibble(
    biomarker = gsub("clean_conc_log10_", "", var),
    n = n,
    mean_log10 = round(mean_log10, 2),
    sd_log10   = round(sd_log10, 2),
    geom_mean  = round(geom_mean, 2),
    geom_sd    = round(geom_sd, 2)
  )
}

# Apply to all biomarkers and combine
summary_table <- bind_rows(lapply(biomarkers, summarize_biomarker, df = small))



######################################################################################################################################
# STEP 1: are kids who are enrolled later in their diarrheal more likely to have dysentery than kids enrolled early in their illness #
######################################################################################################################################
# outcome = dysentery 
# predictor = number of days the child had diarrhea prior to enrollment 

model <- geeglm(dysentery ~ scr_diar_days + country, family = binomial (link="log"), data = small, id = child_id, corstr = "exchangeable")
risk <- exp(summary(model)$coefficients[2,1])
lb <- exp(summary(model)$coefficients[2,1] - 1.96 * summary(model)$coefficients[2,2])
ub <- exp(summary(model)$coefficients[2,1] + 1.96 * summary(model)$coefficients[2,2])
est_CI<-paste0(sprintf("%.2f", risk), " (", sprintf("%.2f", lb), ", ", sprintf("%.2f", ub), ")")



###########################################################################################################################################################################
# STEP 2: assuming that dysentery is more likley to be reported later, does detection of hemoglobin outperform maternal report of blood when kids present to care earlier #
###########################################################################################################################################################################
# using "any reported" because we want all of the data 

# keeping only the needed variables
haem <- EFGH_substudy_rowperpid_ALLDATA %>%
  select(pid, clean_conc_log10_HAEM, country, dysentery, scr_diar_days) %>%
  filter(!is.na(pid)) %>% # filter out any missing PIDs 
  filter(!is.na(clean_conc_log10_HAEM)) %>% # filter out the concentrations with no values
  filter(!is.na(dysentery)) # remove missing dysentery variables

# convert dysentery to a factor
haem$dysentery <- as.factor(haem$dysentery)

# we need 3 datasets (one with 1-2 days of diarrhea, 3 days of diarrhea, and 4+ days of diarrhea )
days1_2 <- haem %>%
  filter(scr_diar_days==1 | scr_diar_days==2)

day3 <- haem %>%
  filter(scr_diar_days==3)

days4_plus <- haem %>%
  filter(scr_diar_days==4 | scr_diar_days==5 | scr_diar_days==6 | scr_diar_days==7)

# ROC curve (one curve per biomarker)
# x = continuous biomarker
# y = dichotomous dysentery
# cutpointr(data, x, y)

# ROC curve plot
# The ROC curve plots the TPR against the FPR at various threshold settings. Each point on the curve corresponds to a different threshold used to classify log_10_biomarker_HAEM as predicting dysentery
all.cut_days1_2 <- cutpointr(days1_2$clean_conc_log10_HAEM, days1_2$dysentery, method = maximize_metric, metric = youden, na.rm=TRUE) %>% add_metric(metric=list(ppv, npv))

all.cut_day3 <- cutpointr(day3$clean_conc_log10_HAEM, day3$dysentery, method = maximize_metric, metric = youden, na.rm=TRUE) %>% add_metric(metric=list(ppv, npv))

all.cut_days4_plus <- cutpointr(days4_plus$clean_conc_log10_HAEM, days4_plus$dysentery, method = maximize_metric, metric = youden, na.rm=TRUE) %>% add_metric(metric=list(ppv, npv))


######################################################################################################################################
# STEP 3: i am curious now if you assocaite haem ~ day of presenatation #
######################################################################################################################################
# outcome = log10 HAEM 
# predictor = number of days the child had diarrhea prior to enrollment 

# LOOK AT ALL 4 BIOMARKERS
risk_func <-function(biomarker) {
model <- geeglm(biomarker ~ scr_diar_days + country, family = gaussian, data = small, id = child_id, corstr = "exchangeable")
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

CAL_cat <- risk_func(small$clean_conc_log10_CAL)
HAEM_cat <- risk_func(small$clean_conc_log10_HAEM)
MPO_cat <- risk_func( small$clean_conc_log10_MPO)
NGAL_cat <- risk_func(small$clean_conc_log10_NGAL)

bind <- data.frame(rbind(HAEM_cat, NGAL_cat, MPO_cat, CAL_cat))

