# this script is to combine the cleaned original elisa data and the cleaned repeated elisa data 
# 16 Jan 2026: This script is updated based on the new model for the ELISA data 

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

library(dplyr)
library(tidyr)
library(xlsx)
library(stringr)
library(janitor)
library(lubridate)

# load in original data and elisa repeats 
load("R output/ELISA_originals_updatedmodel_12Jan2026.Rda")
load("R output/ELISA_repeats_updatedmodel - 27Jan2026.Rda")

tac <- readRDS("Data drops/tac_processed.Rds")
screening <- readRDS("Data drops/DCS_all_screening.Rds")
enrollment <- readRDS("Data drops/DCS_all_enrollment.Rds")
# don't need rectal_swab for now
rectal_swabs <- readRDS("Data drops/DCS_16_rectal_swab_results.Rds")
medical_management <- readRDS("Data drops/DCS_04a_medical_management.Rds")
# contains culture data at the isolate level and includes media type, processing times, an indicator of whether shigella was identified, and species/serotype. 
# we don't need the isolate data for now
isolates <- readRDS("Data drops/DCS_isolate_level.Rds")
# contains both tac and culture data at the participant level and includes variables for species/serotype, and whether shigella was identified via culture and/or attributable by TAC. 
shigella <- readRDS("Data drops/DCS_shigella.Rds")
# we have a variable "dysentery" in the DCS_severity_scores.Rds dataset that collates dysentery variables as reported by caregiver or clinician at screening, enrollment, discharge, and the diarrhea diary. 
severity <- readRDS("Data drops/DCS_severity_scores.Rds")
anthro <- readRDS("Data drops/DCS_anthro.Rds")
wealth <- readRDS("Data drops/DCS_wealth_index.Rds")
multiple_enroll <- readRDS("Data drops/DCS_multiple_enrollments.Rds")

# bind original and repeat dataset together
ELISA_originals_repeats <- rbind(ELISA_originals, ELISA_repeats_cleaned)

######################
# clean the datasets #
######################

# clean up elisa data
ELISA_clean1 <- ELISA_originals_repeats %>%
  mutate(clean_conc_log10 = log10(clean_concentration)) # create log10 cleaned concentration variable

# checking pid count
pid_counts_elisa <- ELISA_clean1 %>%
  count(pid, name = "pid_count")

# SOMETHING NOT RIGHT WITH PID COUNTS - everything should equal 4
ELISA_clean2 <- ELISA_clean1 %>%
  left_join(pid_counts_elisa, by = "pid") %>%
  arrange(sid_elisa)

# create an indicator for duplicated sid_elisa with multiple biomarkers
ELISA_clean3 <- ELISA_clean2 %>%
  group_by(sid_elisa, biomarker) %>%
  mutate(dup_biomarker_within_sid = row_number()) %>%
  ungroup()

# remove columns not needed
ELISA_clean4 <- ELISA_clean3 %>%
  select(-c(pid_count, dup_biomarker_within_sid))

# the value that remains as n=-3 in ELISA_clean4 (EF0006037 from peru for MPO) was accidentally duplicated during the run. however, the ods different substantially to the point that we could not confirm which od 
# truly corresponded to EF0006037 from peru. however, on re-testing, this sample was forgotten so we do not have data for this sid. 

# NEED TO ADD VARIABLE FOR KIDS WHO WERE ENROLLED MORE THAN ONCE 
# different than the multiple dataset cleaning below because each pid has 4 rows for the ELISA_clean1
# clean multiple dataset
multiple_enroll_clean_elisa <- multiple_enroll %>%
  mutate(child_id = row_number())  # Create a column with row numbers (so that we can group them analytically)

# Pivot data into longer format
multiple_enroll_clean_elisa1 <- multiple_enroll_clean_elisa %>%
  select(child_id, pid1, pid2, pid3, pid4) %>%
  pivot_longer(cols = starts_with("pid"),  # All pid columns
               names_to = "pid_column",   # New column name to identify pid columns
               values_to = "pid_value",   # The actual pid value
               values_drop_na = TRUE)     # Remove NAs (if needed)

# merged into ELISA_clean4 below
ELISA_multiple <- merge(ELISA_clean4, multiple_enroll_clean_elisa1, by.x = "pid", by.y = "pid_value")

########################################
# NEW version of transforming the data # 
#           ONE ROW PER PID            #
########################################

# Check for duplicates
duplicates <- ELISA_clean4 %>%
  group_by(pid, biomarker) %>%
  filter(n() > 1)

# Summarize the data for each pid and biomarker, taking the maximum of each concentration variable
elisa_summarized <- ELISA_clean4 %>%
  group_by(pid, country, biomarker) %>%
  summarise(
    clean_concentration = max(clean_concentration, na.rm = TRUE),
    clean_conc_log10 = max(clean_conc_log10, na.rm = TRUE),
    imputed_repeated = first(imputed_repeated),
    .groups = "drop")  # To drop the grouping after summarizing

# Now pivot the summarized data to have one row per pid with separate columns for each biomarker
elisa_1row <- elisa_summarized %>%
  pivot_wider(
    names_from = biomarker,
    values_from = c(clean_concentration, clean_conc_log10, imputed_repeated),
    names_glue = "{.value}_{biomarker}")  # Creating new column names with biomarker suffix

# Verify the number of rows per 'pid'
check <- elisa_1row %>%
  group_by(pid) %>%
  summarise(count = n(), .groups = 'drop') %>%
  filter(count != 1) 


##################
# clean tac data #
##################

# checking pid count
pid_counts_tac <- tac %>%
  count(pid, name = "pid_count")

# confirming only 1 pid per child
pid_counts_shigella <- shigella %>%
  count(pid, name = "pid_count")

# confirming only 1 pid per child
pid_counts_severity <- severity %>%
  count(pid, name = "pid_count")

# confirming only 1 pid per child
pid_counts_anthro <- anthro %>%
  count(pid, name = "pid_count")

# wealth index 
# confirming only 1 pid per child
pid_counts_wealth <- wealth %>%
  count(pid, name = "pid_count")

# only keep 
wealth_clean <- wealth %>%
  select(final_quintile, final_index, pid)

# clean mutliple dataset
multiple_enroll_clean <- multiple_enroll %>%
  mutate(child_id = row_number())  # Create a column with row numbers (so that we can group them analytically)

# Pivot data into longer format
multiple_enroll_clean1 <- multiple_enroll_clean %>%
  select(child_id, pid1, pid2, pid3, pid4) %>%
  pivot_longer(cols = starts_with("pid"),  # All pid columns
               names_to = "pid_column",   # New column name to identify pid columns
               values_to = "pid",   # The actual pid value
               values_drop_na = TRUE)     # Remove NAs (if needed)

multiple_enroll_clean1$child_id <- as.factor(multiple_enroll_clean1$child_id)

##################
# merge datasets #
##################

# ONE ROW PER PID
elisa_tac <- merge(elisa_1row, tac, by="pid", all.x = T)
elisa_tac_cul <- merge(elisa_tac, shigella, by="pid", all.x = T)
elisa_tac_cul_sev <- merge(elisa_tac_cul, severity, by="pid", all.x = T)
elisa_tac_cul_sev_scr <- merge(elisa_tac_cul_sev, screening, by="pid", all.x = T)
elisa_tac_cul_sev_scr_enr <- merge(elisa_tac_cul_sev_scr, enrollment, by="pid", all.x = T)
elisa_tac_cul_sev_scr_enr_anth <- merge(elisa_tac_cul_sev_scr_enr, anthro, by="pid", all.x = T)
elisa_tac_cul_sev_scr_enr_anth_wlth <- merge(elisa_tac_cul_sev_scr_enr_anth, wealth_clean, by="pid", all.x = T)
elisa_tac_cul_sev_scr_enr_anth_wlth_mult <- merge(elisa_tac_cul_sev_scr_enr_anth_wlth, multiple_enroll_clean1, by="pid", all.x = T)

# removing the antilog values last minute
elisa_tac_cul_sev_scr_enr_anth_wlth_mult <- elisa_tac_cul_sev_scr_enr_anth_wlth_mult %>%
  select(-c(clean_concentration_CAL, clean_concentration_HAEM, clean_concentration_MPO, clean_concentration_NGAL, groupid))

#################
# save datasets #
#################

# CREATE GATES DATASET (pid, sid, country, biomarker, clean_conc, assay_date) - in datalock folder 
EFGH_substudy_rowperpid_gates <- elisa_1row %>%
  select(pid, country, clean_conc_log10_CAL, clean_conc_log10_HAEM, clean_conc_log10_MPO, clean_conc_log10_NGAL)
save(EFGH_substudy_rowperpid_gates, file = "C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/R output/EFGH_substudy_rowperpid_GATES_LOCKEDv3.rda")

# one row per pid (all of the data)
EFGH_substudy_rowperpid_ALLDATA <- elisa_tac_cul_sev_scr_enr_anth_wlth_mult
save(EFGH_substudy_rowperpid_ALLDATA, file = "C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/R output/EFGH_substudy_rowperpid_ALLDATA_LOCKEDv3.rda")

# one row per biomarker (only ELISA data)
EFGH_substudy_rowperbiomarker_ELISA <- ELISA_multiple
save(EFGH_substudy_rowperbiomarker_ELISA, file = "C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/R output/EFGH_substudy_rowperbiomarker_ELISA_LOCKEDv3.rda")

# datsets by SITE 
EFGH_substudy_rowperbiomarker_ELISA_Bangladesh <- EFGH_substudy_rowperbiomarker_ELISA %>%
  filter(country=="BANGLADESH")
write.csv(EFGH_substudy_rowperbiomarker_ELISA_Bangladesh, "C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/R output/EFGH_substudy_rowperbiomarker_ELISA_Bangladesh.csv", row.names = FALSE)

EFGH_substudy_rowperbiomarker_ELISA_Kenya <- EFGH_substudy_rowperbiomarker_ELISA %>%
  filter(country=="KENYA")
write.csv(EFGH_substudy_rowperbiomarker_ELISA_Kenya, "C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/R output/EFGH_substudy_rowperbiomarker_ELISA_Kenyav2.csv", row.names = FALSE)

EFGH_substudy_rowperbiomarker_ELISA_Malawi <- EFGH_substudy_rowperbiomarker_ELISA %>%
  filter(country=="MALAWI")
write.csv(EFGH_substudy_rowperbiomarker_ELISA_Malawi, "C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/R output/EFGH_substudy_rowperbiomarker_ELISA_Malawi.csv", row.names = FALSE)

EFGH_substudy_rowperbiomarker_ELISA_Pakistan <- EFGH_substudy_rowperbiomarker_ELISA %>%
  filter(country=="PAKISTAN")
write.csv(EFGH_substudy_rowperbiomarker_ELISA_Pakistan, "C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/R output/EFGH_substudy_rowperbiomarker_ELISA_Pakistan.csv", row.names = FALSE)

EFGH_substudy_rowperbiomarker_ELISA_Peru <- EFGH_substudy_rowperbiomarker_ELISA %>%
  filter(country=="PERU")
write.csv(EFGH_substudy_rowperbiomarker_ELISA_Peru, "C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/R output/EFGH_substudy_rowperbiomarker_ELISA_Peru.csv", row.names = FALSE)

EFGH_substudy_rowperbiomarker_ELISA_TheGambia <- EFGH_substudy_rowperbiomarker_ELISA %>%
  filter(country=="THE GAMBIA")
write.csv(EFGH_substudy_rowperbiomarker_ELISA_TheGambia, "C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/R output/EFGH_substudy_rowperbiomarker_ELISA_TheGambia.csv", row.names = FALSE)



