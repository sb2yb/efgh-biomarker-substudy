# this script is to bring in the current ELISA data that has been repeated and check it against the repeat list sent to the sites
# 14 Jan 2026: This script is updated based on the new model for the ELISA data - we can no longer go off of the list of repeats as the flagged samples have changed 

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# avoid scientific notation
options(scipen = 999)

library(dplyr)
library(tidyr)
library(openxlsx)
library(stringr)
library(janitor)

# load in elisa repeats - these are the repeat results - where we choose 1 for each sample - the samples that the teams re-ran according to the old results
# 14 Jan 2026: not all of these will apply anymore - some of the samples we asked them to repeat didn't need to be repeated according to the NEW model 
#              we will only keep the repeated results where the OLD data and the NEW data indicated that the sample needed to be repeated 
load("R output/elisa_repeats - updatedmodel - 12Jan2026.Rda") # NEEDS TO BE NEW ONE SO THAT THE CONCENTRATIONS ARE CORRECT! 

# the samples that would require being repeated (based on the new model)
bg_list_NEW <- read.csv("C:/Users/sb2yb/Box/ELISA Inflammatory Markers/Forms and SOPs/Repeats/bangladesh_updatedmodel_newrepeats.csv")
kenya_list_NEW <- read.csv("C:/Users/sb2yb/Box/ELISA Inflammatory Markers/Forms and SOPs/Repeats/kenya_repeats_updatedmodel_newrepeats.csv")
malawi_list_NEW <- read.csv("C:/Users/sb2yb/Box/ELISA Inflammatory Markers/Forms and SOPs/Repeats/malawi_repeats_updatedmodel_newrepeats.csv")
gambia_list_NEW <- read.csv("C:/Users/sb2yb/Box/ELISA Inflammatory Markers/Forms and SOPs/Repeats/gambia_repeats_updatedmodel_newrepeats.csv")
pakistan_list_NEW <- read.csv("C:/Users/sb2yb/Box/ELISA Inflammatory Markers/Forms and SOPs/Repeats/pakistan_repeats_updatedmodel_newrepeats.csv")
peru_list_NEW <- read.csv("C:/Users/sb2yb/Box/ELISA Inflammatory Markers/Forms and SOPs/Repeats/peru_repeats_updatedmodel_newrepeats.csv")

# create column in repeat_df that removes D, D2, and D4 so that we can get back the original sid
# create column that stores D,D2,D4,D8
repeat_df <- repeat_df %>%
  separate(label, into = c("original_sid", "suffix"), sep = "(?=D\\d?$)", remove = FALSE, fill= "right") %>% # fill=right will leave the suffix column as NA rather than causing errors (for the 13 from BG without D,D2,D4,D8 on the end
  select(-c(cv))

# the above code did not catch D16 - manually correcting
repeat_df <- repeat_df %>%
  mutate(original_sid = case_when(label == "201353D16" ~ "201353",
      TRUE ~ original_sid),
         suffix = case_when(label == "201353D16" ~ "D16",
      TRUE ~ suffix ))

# first: format country lists to have the same columns 
bg_list_NEW <- bg_list_NEW %>%
  select(-c(X, cv, st_whole_id, st_whole_dt, st_whole_tm, hmst_specid, spec_whole_dt, spec_whole_tm, st_whole_blood))

gambia_list_NEW <- gambia_list_NEW %>%
  select(-c(X, cv, st_whole_id, st_whole_dt, st_whole_tm, hmst_specid, spec_whole_dt, spec_whole_tm, st_whole_blood))

kenya_list_NEW <- kenya_list_NEW %>%
  select(-c(X, cv, st_whole_id, st_whole_dt, st_whole_tm, hmst_specid, spec_whole_dt, spec_whole_tm, st_whole_blood))

malawi_list_NEW <- malawi_list_NEW %>%
  select(-c(X, cv, st_whole_id, st_whole_dt, st_whole_tm, hmst_specid, spec_whole_dt, spec_whole_tm, st_whole_blood))

pakistan_list_NEW <- pakistan_list_NEW %>%
  select(-c(X, cv, st_whole_id, st_whole_dt, st_whole_tm, hmst_specid, spec_whole_dt, spec_whole_tm, st_whole_blood))

peru_list_NEW <- peru_list_NEW %>%
  select(-c(X, cv, st_whole_id, st_whole_dt, st_whole_tm, hmst_specid, spec_whole_dt, spec_whole_tm, st_whole_blood))

# make sid_elisa character so that they'll bind
bg_list_NEW <- bg_list_NEW %>% mutate(sid_elisa = as.character(sid_elisa))
kenya_list_NEW <- kenya_list_NEW %>% mutate(sid_elisa = as.character(sid_elisa))
malawi_list_NEW <- malawi_list_NEW %>% mutate(sid_elisa = as.character(sid_elisa))
gambia_list_NEW <- gambia_list_NEW %>% mutate(sid_elisa = as.character(sid_elisa))
pakistan_list_NEW <- pakistan_list_NEW %>% mutate(sid_elisa = as.character(sid_elisa))
peru_list_NEW <- peru_list_NEW %>% mutate(sid_elisa = as.character(sid_elisa))


# add tracking id to peru list in place of sample id so that it links up
peru_crossref = read.xlsx("C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/stool_db_02072025.xlsx",1)
peru_crossref <- peru_crossref %>% select(SampleID, TrackingID)

# i think the best way forward is to remove peru from repeat_df - merge in the tracking id - and then merge peru back into repeat_df
repeat_df_peru <- repeat_df %>%
  filter(country=="PERU")

# getting the sample ID into the repeat dataset
repeat_df_peru1 <- merge(repeat_df_peru, peru_crossref, by.x = "original_sid", by.y = "TrackingID")
repeat_df_peru1 <- repeat_df_peru1 %>%
  select(-original_sid) %>%
  rename(original_sid = SampleID)

# re-order so that i can be properly bound
repeat_df_peru1 <- repeat_df_peru1[ , c(1,13,2:12)] 

# remove peru from repeat dataset
repeat_df_no_peru <- repeat_df %>%
  filter(country!="PERU")

# bringing the repeated data back together 
repeat_df1 <- rbind(repeat_df_no_peru, repeat_df_peru1)

# First let's see of those that were repeated, what still needs repeating? compare ELISA_repeats_cleaned (n=680) to the combined country "NEW" list (need to bind these)
all_sites_NEW <- bind_rows(bg_list_NEW, kenya_list_NEW,malawi_list_NEW,gambia_list_NEW,pakistan_list_NEW,peru_list_NEW)

# BEFORE THIS MERGE WE LOSE THE CALRPROTECTIN DATA FROM EF0006037 (PERU) BECAUSE WE ONLY HAVE DILUTION DATA (D, D2, AND D4) SO WHEN IT'S NOT NEEDED ON THE JOIN BELOW, IT'S LOST 
# we will separate it out here and add it at the end of this script - we will keep the first row as there is no flag associated with it 
PERU_special_case <- repeat_df1 %>%
  filter(original_sid=="EF0006037" & suffix=="D")

# merging the flagged data with the NEW model with the repeat data that we had done 
# we are looking to see of hte 2148 samples that need repeating how many of the 1830 previously repeated samples fit here 
overlap <- left_join(all_sites_NEW, repeat_df1,
                     by = c("country" = "country", "sid_elisa" = "original_sid", "biomarker" = "biomarker"),
                     suffix = c("_NEW", "_OLD"))


# Define the desired suffix order
# everyone has a suffix except the 1 instance of of the BG samples that we repeated off the bat, but it's fine because the D2 dilution is the one we want to keep and it gets selected
suffix_order <- c("D", "D2", "D4", "D8", "D16")

# identifying those with a value of 1 that are to be kept 
overlap1 <- overlap %>%
  group_by(sid_elisa, biomarker, plate_NEW) %>%
  arrange(factor(suffix, levels = suffix_order), .by_group = TRUE) %>%
  mutate(
  keep_flag = case_when(
    # Case 1: sid_elisa present (NEW), label missing (OLD) → keep
    !is.na(sid_elisa) & is.na(label) ~ "1",
    
    # Case 2a: both present → first blank/missing flag_OLD
    !is.na(sid_elisa) & !is.na(label) &
      row_number() == which(is.na(flag_OLD) | flag_OLD == "")[1] ~ "1",
    
    # Case 2b: both present → if no blanks/missings, flag last row
    !is.na(sid_elisa) & !is.na(label) &
      all(!is.na(flag_OLD) & flag_OLD != "") &
      row_number() == n() ~ "1",
    
    # Otherwise
    TRUE ~ "")) %>%
  ungroup()

# keep only the ones flagged with a "1"
overlap2 <- overlap1 %>%
  filter(keep_flag==1)

# now that we have one row per sample that needed repearing (n=2148) we need to clean the repeat data
# NEXT: match flag_OLD with flag_NEW - confirm that it's high concs for both 
table(overlap2$flag_NEW, overlap2$flag_OLD)

overlap3 <- overlap2 %>%
  mutate(clean_concentration= case_when(
    flag_OLD=="Concentration higher than highest standard" ~ (dilution_OLD * 50), # correcting the 18 samples that still had high concs after further dilutuion - they are all HAEM (highest standard is 50)

# ANY NEW sample that is flagged as too low or unable to generate does not need repeated data - this must come before the rest of the cleaning code     
# each biomarker will be assigned the value of it's lowest concentration/2
    (flag_NEW=="Concentration lower than lowest standard" | flag_NEW=="Unable to generate concentration despite OD being in range") & biomarker== "HAEM" ~ 0.00001831894 / 2,
    (flag_NEW=="Concentration lower than lowest standard" | flag_NEW=="Unable to generate concentration despite OD being in range") & biomarker== "CAL" ~ 64.10851 / 2,
    (flag_NEW=="Concentration lower than lowest standard" | flag_NEW=="Unable to generate concentration despite OD being in range") & biomarker== "MPO" ~ 1.520137 / 2,
    (flag_NEW=="Concentration lower than lowest standard" | flag_NEW=="Unable to generate concentration despite OD being in range") & biomarker== "NGAL" ~ 0.002186617 / 2,

    
# NEXT: clean up and finalize the repeat data that we can use     
    !is.na(label) ~ corrected_concentration_OLD,

# NEXT: clean up the data that does not have repeat data 
    is.na(label) & flag_NEW=="Concentration higher than highest standard" & biomarker== "HAEM" ~ dilution_NEW * 50,
    is.na(label) & flag_NEW=="Concentration higher than highest standard" & biomarker == "CAL" ~ dilution_NEW * 840,
    is.na(label) & flag_NEW=="Concentration higher than highest standard" & biomarker == "NGAL" ~ dilution_NEW * 118.8,
    is.na(label) & flag_NEW=="Concentration higher than highest standard" & biomarker == "MPO" ~ dilution_NEW * 100))

# need an indicator for repeats, imputed high values, imputed low values
overlap4 <- overlap3 %>%
  mutate(imputed_repeated= case_when(
    flag_OLD=="Concentration higher than highest standard" ~ "repeated_imputed_high",
    is.na(label) & flag_NEW=="Concentration higher than highest standard" ~ "imputed_high",
    flag_NEW=="Concentration lower than lowest standard" | flag_NEW=="Unable to generate concentration despite OD being in range" ~ "imputed_low",
    !is.na(label) ~ "repeated"))

# a series of checks to make sure that the mutate above works 
table(overlap4$imputed_repeated, useNA = "always")

check_unable <- overlap4 %>%
  filter(flag_NEW=="Unable to generate concentration despite OD being in range")
table(check_unable$clean_concentration, useNA = "always")
table(check_unable$imputed_repeated, useNA = "always")

check_lower <- overlap4 %>%
  filter(flag_NEW=="Concentration lower than lowest standard")
table(check_lower$clean_concentration, useNA = "always")
table(check_lower$imputed_repeated, useNA = "always")

check_higher <- overlap4 %>%
  filter(flag_NEW=="Concentration higher than highest standard" & is.na(label))
table(check_higher$clean_concentration, check_higher$biomarker, useNA = "always")
table(check_higher$imputed_repeated, useNA = "always")

check_higher_repeat <- overlap4 %>%
  filter(flag_NEW=="Concentration higher than highest standard" & !is.na(label))
table(check_higher_repeat$imputed_repeated, useNA = "always")

# FINALIZE THE DATA (if imputed use _NEW data if repeat use _OLD data)
# i think i want to break the dataset in half to do this? 
new_data <- overlap4 %>%
  select(country, sid_elisa, biomarker, well_NEW, od_NEW, dilution_NEW, raw_concentration_NEW, corrected_concentration_NEW, flag_NEW, plate_NEW, assay_date_NEW, pid, clean_concentration, imputed_repeated) %>%
  filter(imputed_repeated=='imputed_high' | imputed_repeated=='imputed_low') %>%
  rename_with(~ gsub("_NEW$", "", .x))

old_data <- overlap4 %>%
  select(country, sid_elisa, biomarker, well_OLD, od_OLD, dilution_OLD, raw_concentration_OLD, corrected_concentration_OLD, flag_OLD, plate_OLD, assay_date_OLD, pid, clean_concentration, imputed_repeated) %>%
  filter(imputed_repeated=='repeated' | imputed_repeated=='repeated_imputed_high') %>%
  rename_with(~ gsub("_OLD$", "", .x))

# bind new_data and old_data on top of eachother
ELISA_repeats_cleaned <- rbind(old_data, new_data)

# add a column called run_type
ELISA_repeats_cleaned <- ELISA_repeats_cleaned %>%
  mutate(run_type = case_when(
    imputed_repeated=='imputed_high' | imputed_repeated=='imputed_low' ~ "original",
    imputed_repeated=='repeated' | imputed_repeated=='repeated_imputed_high' ~ "repeated"))

# confirming same subset as before 
table(ELISA_repeats_cleaned$imputed_repeated)

# clean PERU_special_case and adding it to the cleaned dataset
PERU_special_case <- PERU_special_case %>%
  mutate(pid=6501770) %>%
  mutate(clean_concentration=corrected_concentration) %>%
  mutate(imputed_repeated="no") %>%
  mutate(run_type="original") %>%
  select(country, original_sid, biomarker, well, od, dilution, raw_concentration, corrected_concentration, flag, plate, assay_date, pid, clean_concentration, imputed_repeated, run_type) %>%
  rename(sid_elisa=original_sid)

# adding PERU special case to the cleaned dataset
ELISA_repeats_cleaned <- rbind(ELISA_repeats_cleaned, PERU_special_case)

# remove sample from KENYA - MPO - SID: 201588 as the flag generated does not make sense - OD between STD4 and 5 but no concentration generated
ELISA_repeats_cleaned <- ELISA_repeats_cleaned %>%
  filter(!(country == "KENYA" & biomarker == "MPO" & sid_elisa == 201588))

# save the cleaned repeated dataset
save(ELISA_repeats_cleaned, file="C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/R output/ELISA_repeats_updatedmodel - 27Jan2026.Rda")








