
# this script is to further clean the original (i.e., not the repeats) ELISA dataset 
# 12 Jan 2026: This script is updated based on the new model for the ELISA data - this is a new list of "originals" that do not require further correction

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

library(dplyr)
library(tidyr)
library(xlsx)
library(stringr)
library(janitor)

# load in elisa data
load("R output/elisa_wholestool - dysentery - updatedmodel - 12Jan2026.Rda")

# remove standards and controls
ELISA_wholestool <- ELISA_wholestool %>%
  filter(sid_elisa!="CONTROL1" & sid_elisa!="CONTROL2" & sid_elisa!="STANDARD1" & sid_elisa!= "STANDARD2" & sid_elisa!="STANDARD3" & sid_elisa!= "STANDARD4" & sid_elisa!="STANDARD5" & sid_elisa!="STANDARD6") 

# remove samples with high concs - these should all be accounted for in the repeat script
ELISA_originals <- ELISA_wholestool %>%
  filter(flag=="")

# creating a seperate dataset with the flags - this should match the repeat lists sent to individual sites 
ELISA_originals_flags <- ELISA_wholestool %>%
  filter(flag!="")

# remove columns we don't need 
ELISA_originals <- ELISA_originals %>%
  select(-c(cv, st_whole_id, hmst_specid, st_whole_dt, st_whole_tm, spec_whole_dt, spec_whole_tm, st_whole_blood))

# we need to create a "clean_concentration" varaible to match with the repeated dataset - it will be the same as corrected_concentration for the originals
ELISA_originals <- ELISA_originals %>%
  mutate(clean_concentration = corrected_concentration)

# re-ordering ELISA_originals 
ELISA_originals <- ELISA_originals %>%
  select(country, sid_elisa, biomarker, well, od, dilution, raw_concentration, corrected_concentration, flag, plate, assay_date, pid, clean_concentration) %>%
  mutate(imputed_repeated= "no") %>%
  mutate(run_type= "original")

# save the original dataset - high all high concs removed
save(ELISA_originals, file="C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/R output/ELISA_originals_updatedmodel_12Jan2026.Rda")

# save the flags
save(ELISA_originals_flags, file="C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/R output/ELISA_originals_flags_updatedmodel - 12Jan2026.Rda")



# extra stuff 
bg <- ELISA_originals %>%
  filter(country=="BANGLADESH" & biomarker=="HAEM")


table(bg$corrected_concentration>1.66)

