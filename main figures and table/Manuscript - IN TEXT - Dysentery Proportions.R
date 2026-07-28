
# This script answers the following questions:
#    1. the proportion of kids with dysentery that had blood in their stool prior to enrollment (and therefore would be experiencing dysentery at the time of or within a few days 
#       prior to stool collection) using scr_diar_blood, enroll_diar_blood, enroll_cond_comor___23, disch_fin_diagnoses___9
#    2. the proportion of kids with dysentery who did not have any blood seen prior to or at enrollment but only developed dysentery after stool collection using foll_stool_blood 
#       and dd_blood_0 :dd_blood_14


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

# dysentery dataset
dysentery <- read.csv("Data drops/dysentery_variables.csv")

# keep only pid and country from EFGH dataset
substudy <- EFGH_substudy_rowperpid_ALLDATA %>% select(pid, country)

# merge dysentery dataset into the EFGH sub study dataset 
substudy_dysentery <- left_join(substudy, dysentery, by="pid")

# the proportion of kids with dysentery that had blood in their stool prior to enrollment: scr_diar_blood, enroll_diar_blood, enroll_cond_comor___23, disch_fin_diagnoses___9 - "dysentery_at_or_before_enrollment"
table(substudy_dysentery$dysentery_at_or_before_enrollment, useNA = "always")


# the proportion of kids with dysentery who did not have any blood seen prior to or at enrollment but only developed dysentery after stool collection: foll_stool_blood and dd_blood_0 :dd_blood_14
substudy_small <- substudy_dysentery %>%
  filter(dysentery_at_or_before_enrollment==0)

table(substudy_small$dd_blood_any)
table(substudy_small$foll_stool_blood_mo3)
table(substudy_small$foll_stool_blood_wk4)


substudy_small <- substudy_small %>%
  mutate(dysentery_after_enroll = case_when(
    dd_blood_any==1 | foll_stool_blood_wk4==1 | foll_stool_blood_mo3==1 | (disch_dt != enroll_date & disch_fin_diagnoses___9 == 1) ~ 1,
    TRUE ~ 0))


table(substudy_small$dysentery_after_enroll, useNA = "always")
