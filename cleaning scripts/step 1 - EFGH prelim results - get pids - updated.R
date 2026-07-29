# This script is for ongoing cleaning of the EFGH data. This is NOT for grant re-submissions. 

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

library(dplyr)
library(tidyr)
library(xlsx)
library(stringr)
library(janitor)

#########################
# Import the ELISA data #
#########################

metadata <- function(x) { # x = country name
  
  file.list <- list.files(path = paste0("C:/Users/sb2yb/Dropbox/EFGH_ELISA/",x,"/"),pattern='*.xlsx') ## get a list of all xlsx FILE NAMES (not the file path) in the Dropbox where the automatically uploaded files are 
  file.list <- paste0("C:/Users/sb2yb/Dropbox/EFGH_ELISA/",x,"/", file.list) # this is getting us the file paths for the file names in the line above

  df <- data.frame() # makes an empty dataframe
  df.list<-rep(list(df), length(file.list)) # this line creates XX (however many there are in the folder) empty dataframes - essentially replicating the blank df above however many times we need it 
  names(df.list) <- file.list # naming each of the empty dataframes with the file path (keeps from overwriting the dataframes)
               for (i in file.list){ # for loop - i is a file name, reads it in, takes tab6, and gets tab1, gets the biomarker data to have tab6 be useful
    
  tab6 <- read.xlsx(i,6) # 6 means the 6th tab in the excel file, i means do it however many times there's an xlsx file in that directory
  tab1 <- as.data.frame(t(read.xlsx(i,1))) # 1 means the 1st tab in the excel file, i means do it however many times there's an xlsx file in that directory and then transposing so that it's easier to call
  tab1 <- tab1 %>% row_to_names(row_number=1) # makes first row the column names
  
  # assigning what we pulled from tab1 and putting it into tab6
  tab6$biomarker<-tab1$Biomarker
  tab6$country<-tab1$Country
  tab6$plate<-tab1$`Plate ID` 
  tab6$assay_date<-tab1$`Assay date (DD MONTH YYYY)`
 
  df.list[[i]] <- tab6 # replaces a specific empty dataframe with the dataframe that we just made (so if there are 16 empty shells in df.list, this 
                       # replaces them iteratively until all 16 are filled) the for loop is repeated however many times there are excel files in the folder, but the function is run once

  } # end of the for loop
  
  country <- do.call("rbind", df.list) # binds all the dataframes within a country folder to one country dataframe 
  return(country) # returning each individual country dataframe
} # end of the function 


# Calling the country dataframes
Bangladesh <- metadata("Bangladesh")
Kenya <- metadata("Kenya")
Malawi <- metadata("Malawi")
Pakistan <- metadata("Pakistan")
Peru <- metadata("Peru")
Gambia <- metadata("The Gambia")

# fix the peru data as the IDs presented at tracking IDs and not PIDs or SIDs
peru_crossref = read.xlsx("C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/stool_db_02072025.xlsx",1)

# keeping only needed variables - WILL NEED TO USE LATER BUT CAN'T ADD NOW AS IT WILL PREVENT THE RBIND BELOW FROM HAPPENING
peru_crossref <- peru_crossref %>% dplyr::select(TrackingID, SampleID, PID)

# cleaning up peru data 
Peru <- Peru %>%
  mutate(label = gsub("^0+", "", label)) %>% # remove proceeding 0s - this is the TRACKING ID - not actually the SID
  mutate(plate = ifelse(plate %in% c("HEM10042024", "CAL10052024", "MPO10072024", "NGAL10082024"), 1, plate)) %>% # re-numbering plate #1
  mutate(plate = ifelse(plate %in% c("HEM10102024", "CAL10102024", "MPO10112024", "NGAL10122024"), 2, plate)) %>% # re-numbering plate #2
  mutate(plate = ifelse(plate %in% c("MPO10142024"), 3, plate)) %>% # re-numbering plate #3
  mutate(plate = ifelse(grepl("_3$", plate), 3, plate)) %>% # if the plate ends in _3 then make 3
  mutate(plate = ifelse(grepl("_4$", plate), 4, plate)) %>% # if the plate ends in _4 then make 4
  mutate(plate = ifelse(grepl("_5$", plate), 5, plate)) %>% # if the plate ends in _5 then make 5
  mutate(plate = ifelse(grepl("_6$", plate), 6, plate)) %>% # if the plate ends in _6 then make 6
  mutate(plate = ifelse(grepl("_7$", plate), 7, plate)) %>% # if the plate ends in _7 then make 7
  mutate(plate = ifelse(grepl("_8$", plate), 8, plate)) %>% # if the plate ends in _8 then make 8
  mutate(plate = ifelse(grepl("_9$", plate), 9, plate)) # if the plate ends in _9 then make 9

# merge peru_crossref into Peru

# SAMPLE EF0006037 GETS LOST HERE FOR CALPROTECTIN BECAUSE IT WAS A REPEAT ????
# I THINK WE HAVE TO CHANGE THE MERGE SO THAT IT RECOGNIZES THE ID EVEN WITH D,D2,D4, ETC ON THE END ????

Peru <- merge(Peru, peru_crossref, by.x = "label", by.y = "TrackingID", all.x = T)

# replace label (which is the tracking ID) with the sampleID when the sampleID is not missing 
Peru <- Peru %>%
  mutate(label = ifelse(!is.na(SampleID), SampleID, label))

# remove sampleID
Peru <- Peru %>%
  select(-c(SampleID, PID)) 

# PERU ACCIDENTALLY RAN 28 SAMPLES TWICE. SOME HAVE ODS THAT ARE NOT AT ALL SIMILIAR WHILE SOME ARE. THEY WERE ASKED TO REPEAT THOSE THAT WERE NOT SIMILIAR. 
# GOING THROUGH THE GUIDANCE BELOW ON WHICH ONES NEED TO BE REPEATED + SELECTING ONE OF THE TWO RUNS TO BE KEPT FOR THOSE THAT HAD SIMILIAR ODS ON THE DUPLICATE RUN.

# peru_duplicates: compares the ODs across the 28 accidential duplicates - those in green do not need repeating. those in red/brown do. 
# peru_duplicates_retest: the red/brown colored ones from peru_duplicates that need to be re-run because the ODs between the accidentical duplicates runs do not match.  
# peru_duplicates_DONOTretest: the green colored ones from peru_duplicates with similiar ODs - simply choose the OD from the first run and delete the second run 

# plan: bring in peru_duplicates_retest and merge the peru data - create a new column identifying them. 

# peru_duplicates = read.xlsx("C:/Users/sb2yb/Box/ELISA Inflammatory Markers/Forms and SOPs/Repeats/Peru/peru_duplicates.xlsx",1)
peru_duplicates_retest = read.xlsx("C:/Users/sb2yb/Box/ELISA Inflammatory Markers/Forms and SOPs/Repeats/Peru/peru_duplicates_retest.xlsx",1)
peru_duplicates_DONOTretest = read.xlsx("C:/Users/sb2yb/Box/ELISA Inflammatory Markers/Forms and SOPs/Repeats/Peru/peru_duplicates_DONOTretest.xlsx",1)

# first step, remove the Dec runs from peru_duplicates_DONOTretest
peru_duplicates_DONOTretest_small <- peru_duplicates_DONOTretest %>%
  filter(plate!=5) %>% # we want to select the ones we want to REMOVE -then match them up with the Peru dataset and delete the matches. 
  mutate(remove = "yes") %>%
  select(label, well, od, raw_concentration, biomarker, assay_date, remove)

# merge this up with peru data
peru_clean <- merge(Peru, peru_duplicates_DONOTretest_small, by=c("label", "well", "od"), all = T)

# remove the accidental duplicates that did not require a re-run
peru_clean <- peru_clean %>%
  filter(is.na(remove))

# clean up dataset
peru_clean <- peru_clean %>%
  select(label, well, od, raw_concentration.x, dilution, corrected_concentration, cv, flag, biomarker.x, country, plate, assay_date.x) %>%
  rename(biomarker = biomarker.x) %>%
  rename(raw_concentration = raw_concentration.x) %>%
  rename(assay_date = assay_date.x)

# now to deal with the accidental duplicates that had to be run a third time. 
# remove the first and second runs - keep the third run. 

peru_duplicates_retest <- peru_duplicates_retest %>% 
  mutate(remove = "yes") %>%
  select(label, well, od, raw_concentration, biomarker, assay_date, remove)

peru_clean1 <- merge(peru_clean, peru_duplicates_retest, by=c("label", "well", "od"), all = T)

# remove the accidental duplicates that needed a third run
peru_clean1 <- peru_clean1 %>%
  filter(is.na(remove))

# clean up dataset
peru_clean1 <- peru_clean1 %>%
  select(label, well, od, raw_concentration.x, dilution, corrected_concentration, cv, flag, biomarker.x, country, plate, assay_date.x) %>%
  rename(biomarker = biomarker.x) %>%
  rename(raw_concentration = raw_concentration.x) %>%
  rename(assay_date = assay_date.x)



##########################
# BELOW STARTS THE CHECK #
##########################

# checking in on peru data 
peru_check <- Peru %>%
  filter((plate==7 & biomarker=="Calprotectin") | (plate==7 & biomarker=="HAEM") | (plate==7 & biomarker=="MPO") | (plate==8 & biomarker=="Calprotectin") |
           (plate==8 & biomarker=="NGAL") |  (plate==9 & biomarker=="NGAL"))

# need to separate out D, D2, and D4 into separate column
peru_check1 <- peru_check %>%
  mutate(
    serial_dilution = case_when(
      grepl("D2$", label) ~ "D2",
      grepl("D4$", label) ~ "D4",
      grepl("D$", label) ~ "D",
      TRUE ~ NA_character_),
    tracking_id = gsub("D2$|D4$|D$", "", label)  # Remove D, D2, or D4 from label
  )

# remove standards and controls 
# rename columns to make it clear that these are repeats 
peru_check2 <- peru_check1 %>%
  filter(!grepl("STANDARD|CONTROL", label, ignore.case = TRUE)) %>%
  mutate(biomarker = ifelse(biomarker == "Calprotectin", "CAL", biomarker)) %>%
  rename(repeat_well = well) %>%
  rename(repeat_plate = plate) %>%
  rename(repeat_od = od) %>%
  rename(repeat_raw_concentration = raw_concentration) %>%
  rename(repeat_corrected_concentration = corrected_concentration) %>%
  rename(repeat_assay_date = assay_date) %>%
  rename(repeat_flag = flag) %>%
  select(-c(country, cv))

# check for duplicates in the serial dilution data:
problem <- peru_check2 %>%
  mutate(duplicate_flag = ifelse(duplicated(.) | duplicated(., fromLast = TRUE), "Duplicate", "Unique"))

# remerge peru_crossref in (bringing in sample IDs and PIDs linked to the tracking ID)
peru_check3 <- merge(peru_check2, peru_crossref, by.x = "tracking_id", by.y = "TrackingID", all.x = T)

# Ensure well format is consistent (e.g., A1, B1, ..., H12)
peru_check4 <- peru_check3 %>%
  mutate(
    well_row = substr(repeat_well, 1, 1),  # Extract row (A–H)
    well_col = as.integer(substr(repeat_well, 2, nchar(repeat_well)))  # Extract column (1–12)
  ) %>%
  arrange(repeat_plate, biomarker, well_col, well_row) %>%
  select(-well_row, -well_col)   # Remove helper columns if not needed


# MERGE IN A LIST WITH IDS THAT NEED 1) SERIAL DILUTIONS 
# bring in my repeat csv file (these are the samples that I said needed to be further diluted)
peru_repeat_list <- read.csv("C:/Users/sb2yb/Box/ELISA Inflammatory Markers/Forms and SOPs/Repeats/Peru/peru_repeats_FINAL.csv")

# only keep wanted columns
peru_repeat_list1 <- peru_repeat_list %>% 
  select(sid_elisa, well, biomarker, plate) %>%
  mutate(needs_serial_dilution = "Serial Dilution Required") %>%
  rename(SampleID = sid_elisa) %>%
  rename(original_well = well) %>%
  rename(original_plate = plate)

peru_repeat_list1$original_plate <- as.character(peru_repeat_list1$original_plate)

# merging the list of samples that needed to be repeated into the list of samples that was repeated
peru_check5 <- merge(peru_check4, peru_repeat_list1, by = c("SampleID", "biomarker"), all.x= T)

# Ensure well format is consistent (e.g., A1, B1, ..., H12)
peru_check6 <- peru_check5 %>%
  mutate(
    well_row = substr(repeat_well, 1, 1),  # Extract row (A–H)
    well_col = as.integer(substr(repeat_well, 2, nchar(repeat_well)))  # Extract column (1–12)
  ) %>%
  arrange(repeat_plate, biomarker, well_col, well_row) %>%
  select(-well_row, -well_col)

# MERGE IN A LIST WITH IDS THAT NEED 2) TO BE REPEATED BECAUSE THEY APPEARED AS DUPLICATES
peru_duplicates_list <- read.xlsx("C:/Users/sb2yb/Box/ELISA Inflammatory Markers/Forms and SOPs/Repeats/Peru/peru_duplicates_retest.xlsx",1)

peru_duplicates_list1 <- peru_duplicates_list %>%
  mutate(plate = ifelse(grepl("6$", plate), 6, plate)) %>% # if the plate ends in _6 then make 6
  mutate(plate = ifelse(grepl("7$", plate), 7, plate)) %>% # if the plate ends in _7 then make 7
  mutate(duplicate_run = "Duplicate Run") %>%
  rename(SampleID = label) %>%
  rename(original_duplicate_well = well) %>%
  rename(original_duplicate_plate = plate) %>%
  select(SampleID, original_duplicate_well, biomarker, original_duplicate_plate, duplicate_run) 

# assign unique grouping id to samples with the same sample ID
peru_duplicates_list1 <- peru_duplicates_list1 %>%
  mutate(group_id = as.integer(factor(SampleID)))

# AT THIS POINT EVERYTHING SEEMS TO HAVE MERGED SEAMLESSLLY - 19mar25

# bring in peru pending samples (samples that were left over and not run since there wasn't a full plate)
peru_pending_list <- read.xlsx("C:/Users/sb2yb/Box/ELISA Inflammatory Markers/Forms and SOPs/Repeats/Peru/peru_pending_samples.xlsx",1)
peru_pending_list <- peru_pending_list %>%
  select(SampleID) %>%
  mutate(pending_sample = "Pending Sample")

# try merging the pending sample list into peru_check6 to figure out the samples run that were not serial dilutions
peru_check7 <- merge(peru_check6, peru_pending_list, by= 'SampleID', all = T)

# create a column that identifies samples where neither a serial dilution or pending sample is flagged
peru_check8 <- peru_check7 %>%
  mutate(check_duplicate = ifelse(is.na(needs_serial_dilution) & is.na(pending_sample), "Duplicate?", NA_character_))

# create a dataset that will check the "check duplicate" column against the duplicate list provided by bri'anna
peru_check9 <- merge(peru_check8, peru_duplicates_list1, by = c("SampleID", "biomarker"), all= T)

peru_check10 <- peru_check9 %>%
  select(SampleID, biomarker, tracking_id, repeat_well, repeat_plate, PID, original_well, original_plate, needs_serial_dilution, pending_sample, check_duplicate, original_duplicate_well, original_duplicate_plate, duplicate_run, group_id)

# Ensure well format is consistent (e.g., A1, B1, ..., H12)
peru_check11 <- peru_check10 %>%
  mutate(
    well_row = substr(repeat_well, 1, 1),  # Extract row (A–H)
    well_col = as.integer(substr(repeat_well, 2, nchar(repeat_well)))  # Extract column (1–12)
  ) %>%
  arrange(repeat_plate, biomarker, well_col, well_row) %>%
  select(-well_row, -well_col)



# NOT SURE THAT THIS CODE BELOW IS FOR??

# the issue is merging in the repeat data 
# merge peru_duplicates_list1 (the list of samples that need 1 more run) into peru_check6 (the samples that were repeated)
peru_check7 <- merge(peru_check6, peru_duplicates_list1, by = c("SampleID"), all.x= T)


# single out the IDs that don't appear on the repeat or duplicate list, but were run anyway
peru_check8 <- peru_check7 %>%
  mutate(unknown = case_when(
    is.na(needs_serial_dilution) & is.na(duplicate_run) ~ "neither",
    TRUE ~ NA_character_  # Leave other rows as NA or you can specify a different response option here
  ))

# THIS LIST IS NOT RIGHT - SOMETHING WENT WRONG ON A MERGE ABOVE !
#paul_confirm <- peru_check8 %>% 
# filter(!is.na(unknown)) %>%
#  select(SampleID, tracking_id, PID, biomarker, diluted_well, diluted_plate, od) 


# find out where there are samples on both the repeat and duplicate lists - there's only 1 - EF0006037
# duplicate_repeat_lists <- merge(peru_repeats, peru_duplicates, by = c("SampleID", "biomarker", "well"))


#################
# END THE CHECK #
#################


# The Gambia needs ST added on to the end of their IDs
Gambia$label <- ifelse(nchar(Gambia$label) == 7, paste0(Gambia$label, "ST"), Gambia$label)

# bind the country datasets together
ELISA <- rbind(Bangladesh, Gambia, Kenya, Malawi, Pakistan, peru_clean1)

####################################
# clean the combined ELISA dataset #
####################################

# standardize the biomarker, country, and plate variables across sites
ELISA$biomarker <- toupper(ELISA$biomarker) # clean biomarker data - everything capitalized and in one acronym

ELISA <- ELISA %>%
  mutate(biomarker = case_when(
    biomarker == "MPO" | biomarker == "MYELOPEROXIDASE" ~ "MPO",
    biomarker == "CALPROTECTIN" | biomarker == "CAL" ~ "CAL",
    biomarker == "NGAL" | biomarker== "LIPOCALIN" | biomarker=="LIPOCALIN-2" ~ "NGAL",
    biomarker == "HAEM" | biomarker =="HEMOGLOBIN" | biomarker =="HAEMOGLOBIN" ~ "HAEM"))

ELISA$country <- toupper(ELISA$country) # clean country data - everything captialized 

ELISA$plate <- as.numeric(gsub("\\D", "", ELISA$plate)) # clean plate data - pull out only the plate number 

# cleaning a typo in duplicate data so that it goes to the right place
ELISA <- ELISA %>%
  mutate(label= case_when(
    label=="202198D1" & country== "KENYA" ~ "202198D",
    TRUE ~ label)) 

################################################
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#          SEPARATE OUT THE REPEATS            #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
################################################

# creating a repeats dataframe that only contains those that end in D, D2, or D4 (excludes any label that contains "STANDARD")
repeat_df <- ELISA %>%
  filter((grepl("D$", label) | grepl("D2$", label) | grepl("D4$", label) | grepl("D8$", label) | grepl("D16$", label)) & 
      !grepl("STANDARD", label) |
# keep all the 13 leftover samples from the BG CAL plate because we ran repeats on them right off the bat and need the original 1:1 dilution for repeat checks
  ((label == 1510626 | label == 1402348 | label == 1510555 | 
      label == 1510570 | label == 1510581 | label == 1510562 | 
      label == 1510596 | label == 1510615 | label == 1510603 | 
      label == 1510632 | label == 1510659 | label == 1510644 | 
      label == 1510667) & biomarker == "CAL"))


# removing SIDs that end in D, D2, or D4 (excludes any label that contains "STANDARD")
ELISA <- ELISA %>%
  filter(!((grepl("D$", label) | grepl("D2$", label) | grepl("D4$", label) | grepl("D8$", label) | grepl("D16$", label)) & 
        !grepl("^STANDARD", label)))


################################################
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
################################################

# double checking biomarker and country were appropriately cleaned
table(ELISA$biomarker)
table(ELISA$country)

# remove any accidental spaces in the datasets to be merged
ELISA$label <- str_trim(ELISA$label)

# The Gambia added the date to their plate 1, which caused an issue with the naming
ELISA <- ELISA %>%
  mutate(plate = case_when(
    plate==124042024 ~ 1,
    TRUE ~ plate))


############################################################
# Modifications to the ELISA dataset based on ELISA checks #
############################################################

# This is where all of the corrections to the ELISA dataset go
# Once an error is caught in the "check" below, it should be moved up here so that I don't keep "checking" the same issues

# removing the duplicates from kenya (determined above)
ELISA <- ELISA %>%
  filter(plate!=3 | country!="KENYA" | label!="200060")

ELISA <- ELISA %>%
  filter(plate!=3 | country!="KENYA" | label!="200688")

ELISA <- ELISA %>%
  filter(plate!=3 | country!="KENYA" | label!="200830")

ELISA <- ELISA %>%
  filter(plate!=3 | country!="KENYA" | label!="200947")

# BG sent an email on 28 Nov 2023 stating: 
# "We have carried out assays for MPO (20th November), Calprotectin (21st November) and NGAL (22nd November) and uploaded the data in Shinyapp. 
# Unintentionally, we made a mistake for one PID in the plate layout. The correct ID will be 1200566 instead of 1506066 at G2 position of the 96 well plate."
ELISA <- ELISA %>%
  mutate(label = case_when(
    label=="1506066" & well=="G2" & plate==9 & country=="BANGLADESH" ~ "1200566", # being as specific as possible in case this is an ID used later
    TRUE ~ label)) # keep everything else the same 

# The following 2 SIDs  could not be located for Malawi. Mutsai sent an email on 19 Dec indicating the PIDs that go along with the ELISA SIDs:
# 4.	CXN1W0F1 (3106628) – I will change this to: CXN1WOF1 as there was a mistype between O and 0.
# 8.	CXN1BVF1 (3106197) – I will use the spec_whole_id: CXN1VBF1 as there was a mistype between BV and VB.
# The other 6 errors will be fixed below as there is at least one of the three SIDs that match the PID given. 
ELISA <- ELISA %>%
  mutate(label = case_when(
    label=="CXN1W0F1" & well=="A8" & plate==2 & country=="MALAWI" ~ "CXN1WOF1", # being as specific as possible in case this is an ID used later
    TRUE ~ label)) # keep everything else the same 

ELISA <- ELISA %>%
  mutate(label = case_when(
    label=="CXN1BVF1" & well=="G6" & plate==2 & country=="MALAWI" ~ "CXN1VBF1", # being as specific as possible in case this is an ID used later
    TRUE ~ label)) # keep everything else the same 

# From Kenya on 1/9/23 - SID 211791 on plate 5, well E8 should be SID 201791
ELISA <- ELISA %>%
  mutate(label = case_when(
    label=="211791" & well=="E8" & plate==5 & country=="KENYA" ~ "201791", # being as specific as possible in case this is an ID used later
    TRUE ~ label)) # keep everything else the same 

# From Bangladesh on 1/9/23 - SID 1606432 on plate 10, well C2 should be SID 1506432 
ELISA <- ELISA %>%
  mutate(label = case_when(
    label=="1606432" & well=="C2" & plate==10 & country=="BANGLADESH" ~ "1506432", # being as specific as possible in case this is an ID used later
    TRUE ~ label)) # keep everything else the same 

# From Malawi on 13 Feb 24, for all Plate 3 assays, G6 should be corrected to CXN18HF1(3102413) and G9 should remain as CXN1AJF1(3103632) 
ELISA <- ELISA %>%
  mutate(label = case_when(
    label=="CXN1AJF1" & well=="G6" & plate==3 & country=="MALAWI" ~ "CXN18HF1", # being as specific as possible in case this is an ID used later
    TRUE ~ label)) # keep everything else the same 

ELISA <- ELISA %>%
  mutate(label = case_when(
    label=="CXNILNF1" & well=="G10" & plate==3 & country=="MALAWI" ~ "CXN1LN", # being as specific as possible in case this is an ID used later
    TRUE ~ label)) # keep everything else the same 

# From Malawi on 7 May 2024 - For plate04, well A8 should be CXN205F1 whilst well H9 should remain CXN20SF1.
ELISA <- ELISA %>%
  mutate(label = case_when(
    label=="CXN20SF1" & well=="A8" & plate==4 & country=="MALAWI" ~ "CXN205F1", # being as specific as possible in case this is an ID used later
    TRUE ~ label)) # keep everything else the same 

# From Malawi on 7 May 2024 - Indeed had accidentally run sample CXN1ZFF1 in plate1 and plate04. CXN1ZFF1 is to be removed from Plate 4 series (well A6)
ELISA <- ELISA %>%
  filter(label!="CXN1ZFF1" | well!="A6" | plate!=4 | country!="MALAWI") 

# From The Gambia on 7 May 2024 - The SIDs for  HAEM and NGAL within the plate 1 series in well B5 and H7 were typos. The correct SIDs should be 7100383 in well B5 and 7100762 in well H7.
# Within the MPO and CAL plate 1 series, the SIDs are correct in those respective wells.
ELISA <- ELISA %>%
mutate(label = case_when(
  label=="7100390ST" & well=="B5" & plate==1 & country=="THE GAMBIA" & biomarker=="HAEM" ~ "7100383ST", # being as specific as possible in case this is an ID used later
  TRUE ~ label)) # keep everything else the same 

ELISA <- ELISA %>%
  mutate(label = case_when(
    label=="7100390ST" & well=="B5" & plate==1 & country=="THE GAMBIA" & biomarker=="NGAL" ~ "7100383ST", # being as specific as possible in case this is an ID used later
    TRUE ~ label)) # keep everything else the same 

ELISA <- ELISA %>%
  mutate(label = case_when(
    label=="7100768ST" & well=="H7" & plate==1 & country=="THE GAMBIA" & biomarker=="HAEM" ~ "7100762ST", # being as specific as possible in case this is an ID used later
    TRUE ~ label)) # keep everything else the same 

ELISA <- ELISA %>%
  mutate(label = case_when(
    label=="7100768ST" & well=="H7" & plate==1 & country=="THE GAMBIA" & biomarker=="NGAL" ~ "7100762ST", # being as specific as possible in case this is an ID used later
    TRUE ~ label)) # keep everything else the same 

# From the Gambia on 10 May 2024 - The SIDs you shared for Plate 3 are typos. The correct ones are as follows: 7200848 instead of 7100848, 7102798 instead of 7102795, 7200830 instead of 7210830
ELISA <- ELISA %>%
  mutate(label = case_when(
    label=="7100848ST" & well=="C3" & plate==3 & country=="THE GAMBIA" ~ "7200848ST", # being as specific as possible in case this is an ID used later
    TRUE ~ label)) # keep everything else the same 

ELISA <- ELISA %>%
  mutate(label = case_when(
    label=="7102795ST" & well=="E6" & plate==3 & country=="THE GAMBIA" ~ "7102798ST", # being as specific as possible in case this is an ID used later
    TRUE ~ label)) # keep everything else the same 

ELISA <- ELISA %>%
  mutate(label = case_when(
    label=="7210830ST" & well=="A3" & plate==3 & country=="THE GAMBIA" ~ "7200830ST", # being as specific as possible in case this is an ID used later
    TRUE ~ label)) # keep everything else the same 

# Email from The Gambia on 14 June 24: SID 710584ST to be 7105846ST and 7106425TC to be 7106425ST
# follow up email on 25 June: I would like to inform you that our initial correction was incorrect. Please be advised that SID 7105846 in cell A2 is correct, whereas SID 710584ST in cell E3 contains a typo error. The correct SID for cell E3 should be 7105984ST
ELISA <- ELISA %>%
  mutate(label= case_when(
    label=="710584ST" & well=="E3" & country== "THE GAMBIA" & plate==7 ~ "7105984ST",
    TRUE ~ label)) %>%
  
  mutate(label= case_when(
    label=="7106425TC" & well=="C9" & country== "THE GAMBIA" & plate==7 ~ "7106425ST",
    TRUE ~ label)) 

# Malawi plate 6,7,8 upload has -F1 instead of F1 which is causing issues merging in PIDs 
ELISA <- ELISA %>%
  mutate(label = ifelse(country == "MALAWI" & plate == 6, gsub("-", "", label), label))

ELISA <- ELISA %>%
  mutate(label = ifelse(country == "MALAWI" & plate == 7, gsub("-", "", label), label))

ELISA <- ELISA %>%
  mutate(label = ifelse(country == "MALAWI" & plate == 8, gsub("-", "", label), label))


# Malawi has 3 SIDs that need revision - for each SID 2/4 uploads had a dash and the other 2 did not - keeping it consistent to remove them from being flagged
ELISA <- ELISA %>%
mutate(label= case_when(
  label=="CXN2AQ-F1" & country== "MALAWI" ~ "CXN2AQF1",
  TRUE ~ label)) 

ELISA <- ELISA %>%
mutate(label= case_when(
  label=="CXN2AT-F1" & country== "MALAWI" ~ "CXN2ATF1",
  TRUE ~ label)) 

ELISA <- ELISA %>%
mutate(label= case_when(
  label=="CXN2BU-F1" & country== "MALAWI" ~ "CXN2BUF1",
  TRUE ~ label)) 

# Email from Olivia on 5 Aug 24:
# We found that there were two cases of overlapping enrollments where the same child was enrolled before they had been closed out of a previous enrollment. 
# Since these children should not have been eligible for enrollment, we decided to exclude them from all analysis and removed the PIDs (6302484 and 7104928) from all analysis datasets.

# PID 7104928 is linked to SID 7104928ST from the Gambia - Removing the SID data 
ELISA <- ELISA %>%
  filter(!(label == "7104928ST" & country == "THE GAMBIA" & plate == 6))

# Email from Malawi - Throughout the plate 7 series, a duplicate SID was run: CXN2C8-F1 (well A4 and G3)
ELISA <- ELISA %>%
  filter(!(label == "CXN2C8F1" & country == "MALAWI" & plate == 7 & well=="A4"))

# BG has a blank label appearing for well D10 on plate 16 biomarker CAL - removing this so that it doesn't pop up in the checks
ELISA <- ELISA %>%
  filter(!(label == "" & country == "BANGLADESH" & plate == 16 & well=="D10"))


# for malawi plate 9 - remove the - in spec_whole_id
ELISA <- ELISA %>%
  mutate(label = ifelse(country == "MALAWI" & plate == 9, gsub("-", "", label), label))

##############################################
# Standardizing "checking" the ELISA dataset #
##############################################

# step 1: within a plate make sure there are no duplicate sample ids 
# it's possible that an sid could be listed 4 times across all 4 biomarkers in 1 plate series, but there be 2 in 1 plate and 0 in 1 plate, which would give a false assumption that everything is correct
step1_ELISA <- ELISA %>%
  filter(label!="CONTROL1" & label!="CONTROL2" & label!="STANDARD1" & label!= "STANDARD2" & label!="STANDARD3" & label!= "STANDARD4" & label!="STANDARD5" & label!="STANDARD6") %>%
  group_by(country, plate, biomarker, label)%>%
  summarise(n=length(label)) %>%
  ungroup()

step1_ELISA_fix <- step1_ELISA %>% # create a list of what i need to troubleshoot (make sure to move all "fixes" above this code so that I don't have to keep re-fixing)
  filter(n!=1)

# step 2 : making sure each SID is listed 4 times within a country 
# it's possible that what is found in step 1 could be listed here as well, but we need to keep step 1 based on the above scenario given
step2_ELISA <- ELISA %>%
  filter(label!="CONTROL1" & label!="CONTROL2" & label!="STANDARD1" & label!= "STANDARD2" & label!="STANDARD3" & label!= "STANDARD4" & label!="STANDARD5" & label!="STANDARD6") %>%
  group_by(country, label)  %>%
  summarise(n=length(label)) %>%
  ungroup()

# note on step2_ELISA_fix: sample 6037 is one of the samples that was accidentally duplicated. i have a 3rd run for NGAL and HAEM. CAL also needed a serial dilution - so it's clocked on the repeat dataset. 
#                          it was never run on the MPO plate and we decided to cut our losses instead of running a plate for 1 sample.
step2_ELISA_fix <- step2_ELISA %>% # create a list of what i need to troubleshoot (make sure to move all "fixes" above this code so that I don't have to keep re-fixing)
  filter(n!=4)


# ---------------------------------------------------------------------------- #

##############################
# Import Whole Stool dataset #
##############################

# merge in whole stool data so that we can get PIDs
wholestool <- read.csv("Data drops/WholeStool_Data_2025-04-10.csv") # FINAL LOCKED DATASET

#################################
# clean the Whole Stool dataset #
#################################
# Cleaning up the imported PID/SID list from Olivia

# st_whole_id, hmst_specid, spec_whole_id are the IDs we need to check 
wholestool_small <- wholestool %>%
  filter(enroll_site=="Kenya" | enroll_site=="Bangladesh" | enroll_site=="Malawi" | enroll_site=="Pakistan" | enroll_site=="Peru" | enroll_site=="The Gambia") %>% # bring in data from sites who have ELISA data 
  select(enroll_site, pid, st_whole_id, hmst_specid, spec_whole_id, st_whole_dt, st_whole_tm, spec_whole_dt, spec_whole_tm, st_whole_blood)

test <- wholestool_small %>%
  filter(pid== 3107893 | pid== 3106429 | pid== 3106441 | pid== 3106628 | pid==3108238 | pid==3106108 | pid==3106151 | pid==3106197)

# spec_whole_id is supposed to be the "final call" regarding discrepancies between the 3 TAC SIDs
# modify wholestool_small so that if spec_whole_id is missing, but something is in st_whole_id or hmst_special, it is populated

# spec_whole_id - specimen accession form (15)
# hmst_specid - home stool collection form (13b)
# st_whole_id - stool collection form (13)
wholestool_small1 <- wholestool_small %>%
  mutate(spec_whole_id= case_when(
    is.na(spec_whole_id) & is.na(hmst_specid) ~ st_whole_id,
    is.na(spec_whole_id) & is.na(st_whole_id) ~ hmst_specid,
    TRUE ~ spec_whole_id))

# removing all values in the () that were generated from the spec_whole_id mutate statement above as they will not merge on the sid below
wholestool_small1$spec_whole_id <- gsub("\\s*\\([^\\)]+\\)","",as.character(wholestool_small1$spec_whole_id))

# remove any accidental spaces in the datasets to be merged (does not remove interal space)
wholestool_small1$spec_whole_id <- str_trim(wholestool_small1$spec_whole_id)

# clean country data - everything captialized 
wholestool_small1$enroll_site <- toupper(wholestool_small1$enroll_site) 

# remove spaces and dashes within a response in R
wholestool_small1$spec_whole_id <- gsub(" ", "", wholestool_small1$spec_whole_id)
wholestool_small1$spec_whole_id <- gsub("-", "", wholestool_small1$spec_whole_id)

########################################################################
# Modifications to the Whole Stool dataset based on Whole Stool checks #
########################################################################

# This is where all of the corrections to the Whole Stool dataset go
# Once an error is caught in the "check" below, it should be moved up here so that I don't keep "checking" the same issues

# Peru on 23 Jan - I discovered that there are 2 different SIDs for PID 6500326. The one in spec_whole_id has an extra 0, which is preventing linkage: EF00001013.
# The SID for hmst_specid is the correct one: EF0001013. The hmst_specid also matches what we have in the ELISA dataset. It's hte wholestool dataset that is wrong. 
wholestool_small1 <- wholestool_small1 %>%
  mutate(spec_whole_id= case_when(
    spec_whole_id=="EF00001013" ~ "EF0001013",
    TRUE ~ spec_whole_id))

# for malawi on 9 sept 2024:
# PIDS 3111841 and 3112429 have CXN2D9F1 listed as their spec_whole_id (which is what I merge on). 
      # PID 3111841 = SID CXN2D9F1
#     # PID 3112429 = SID CXN2F9F1
# PIDS 3112670 and 3112356 have CXN2F2F1 listed as their spec_whole_id
      # PID 3112670 = SID CXN2FZF1
      # PID 3112356 = SID CXN2F2F1
# PIDS 3112681 and 3112922 have CXN2G0F1 listed as their spec_whole_id
      # PID 3112681 = SID CXN2G0F1
      # PID 3112922 = SID CXN2GQF1
wholestool_small1 <- wholestool_small1 %>%
mutate(spec_whole_id= case_when(
  pid==3112429 ~ "CXN2F9F1", 
  TRUE ~ spec_whole_id)) %>%
mutate(spec_whole_id= case_when(
  pid==3112670  ~ "CXN2FZF1", 
  TRUE ~ spec_whole_id)) %>%
mutate(spec_whole_id= case_when(
  pid==3112922 ~ "CXN2GQF1", 
  TRUE ~ spec_whole_id)) 

# Malawi: PID 3113553 is linked to the wrong SID (we are going to change the spec_whole_id to reflect the true SID)
wholestool_small1 <- wholestool_small1 %>%
  mutate(spec_whole_id= case_when(
    pid==3113553 ~ "CXN2IH", 
    TRUE ~ spec_whole_id)) 


####################################################
# Standardizing "checking" the Whole Stool dataset #
####################################################

# step 1: make sure there are no duplicates 
step1_wholestool_small1 <- wholestool_small1 %>%
  select(enroll_site, pid) %>%
  group_by(enroll_site, pid) %>%
  summarise(n=length(pid)) 

step1_wholestool_small1_fix <- step1_wholestool_small1 %>% # create a list of what i need to troubleshoot (make sure to move all "fixes" above this code so that I don't have to keep re-fixing)
  filter(n!=1)

# ---------------------------------------------------------------------------- #

################################################
# Merge Whole Stool dataset into ELISA dataset #
################################################

# merging on SIDs
# also merging on site in case there are duplicate IDs between sites
ELISA_wholestool <- merge(ELISA, wholestool_small1, by.x = c("country", "label"), by.y = c("enroll_site", "spec_whole_id"), all.x =TRUE) # keeping all ELISA data

ELISA_wholestool <- ELISA_wholestool %>% 
  dplyr::rename(sid_elisa=label) # renaming to make it clear this is the sid that was entered for the elisa

####################################################################################
# Modifications to the ELISA Whole Stool dataset based on ELISA Whole Stool checks #
####################################################################################

# This is where all of the corrections to the ELISA Whole Stool dataset go
# Once an error is caught in the "check" below, it should be moved up here so that I don't keep "checking" the same issues

# SID 201008, 201011, 1502231 DON'T HAVE PIDs

# Troubleshooting: go to wholestool_small1 and search for the above SIDs amongst st_whole_stool, hmst_specid, and spec_whole_id
#                 also check ELISA_wholestool


# On 26 Jan 2023, Kenya provided the following update:
#    SID 201008 = PID	2900379 ON PLATE 1
#    SID 202008	= PID 2901122	ON PLATE 6
#    SID 201011	= PID 2900393	ON PLATE 3
#    SID 202011	= PID 2901131	ON PLATE 6

# SID 201008 = PID 2900379 ON PLATE 1
ELISA_wholestool <- ELISA_wholestool %>%
  mutate(pid= case_when(
    sid_elisa=="201008" & country== "KENYA" & plate== 1 ~ 2900379,
    TRUE ~ pid)) 

# SID 202008 = PID 2901122	ON PLATE 6
# this SID was linked up with 2 different PIDs - removing the rows with the wrong PIDs
ELISA_wholestool <- ELISA_wholestool %>%
  filter(!sid_elisa=="202008" | !pid==2900379)

# SID 201011 = PID 2900393	ON PLATE 3
ELISA_wholestool <- ELISA_wholestool %>%
  mutate(pid= case_when(
    sid_elisa=="201011" & country== "KENYA" & plate== 3 ~ 2900393,
    TRUE ~ pid)) 

# SID 202011 = PID 2901131 ON PLATE 6
# this SID was linked up with 2 different PIDs - removing the rows with the wrong PIDs
ELISA_wholestool <- ELISA_wholestool %>%
  filter(!sid_elisa=="202011" | !pid==2900393)


# For sid 1502231 (can be found in ELISA_wholestool), but nothing can be found in wholestool_small
# Bangladesh uses the PID as the SID. 1502231 is supposed to be 1502238.
ELISA_wholestool <- ELISA_wholestool %>%
  mutate(pid= case_when(
    sid_elisa=="1502231" & plate==4 ~ 1502238,
    TRUE ~ pid)) %>%
  
  mutate(st_whole_id= case_when(
    sid_elisa=="1502231" & plate==4 ~ "1502238", # this is what it actually was in wholestool: 1502238
    TRUE ~ st_whole_id)) %>%

# the sid needs to change as well
mutate(sid_elisa= case_when(
  sid_elisa=="1502231" & plate==4 ~ "1502238", # this is what it actually was in wholestool: 1502238
  TRUE ~ sid_elisa)) 


# There were 8 SIDs that could not be located for Malawi. Mutsai sent an email on 18 Dec indicating the PIDs that go along with the ELISA SIDs:
# 4 and 8 were fixed above as there were typoes with SID naming - the below has at least one SID that matches the PID given by Mutsai
# 1.	CXN1YZF1 (3107893) – I will use the hmst_specid (CXN1YZF1) as this is what matches the ELISA SID
ELISA_wholestool <- ELISA_wholestool %>%
  mutate(pid= case_when(
    sid_elisa=="CXN1YZF1" & country== "MALAWI"  ~ 3107893,
    TRUE ~ pid)) %>%
  
  mutate(hmst_specid= case_when(
    sid_elisa=="CXN1YZF1" & country== "MALAWI" ~ "CXN1YZF1", 
    TRUE ~ hmst_specid)) %>%

# 2.	CXN1V0F1 (3106429) - I will use the hmst_ specid (CXN1V0F1) as this is what matches the ELISA SID
  mutate(pid= case_when(
    sid_elisa=="CXN1V0F1" & country== "MALAWI" ~ 3106429,
    TRUE ~ pid)) %>%
  
  mutate(hmst_specid= case_when(
    sid_elisa=="CXN1V0F1" & country== "MALAWI"~ "CXN1V0F1", 
    TRUE ~ hmst_specid)) %>%

# 3.	CXN1V1F1 (3106441) - hmst_specid and spec_whole_id both seem wrong as the IDs always start with CXN and end in F1. So I think I'll keep the ELISA ID just as is.
  mutate(pid= case_when(
    sid_elisa=="CXN1V1F1" & country== "MALAWI" ~ 3106441,
    TRUE ~ pid)) %>%

# 5.	CXN1ZRF1 (3108238) - I will use the hmst_ specid (CXN1ZRF1) as this is what matches the ELISA SID
  mutate(pid= case_when(
    sid_elisa=="CXN1ZRF1" & country== "MALAWI" ~ 3108238,
    TRUE ~ pid)) %>%
  
  mutate(hmst_specid= case_when(
    sid_elisa=="CXN1ZRF1" & country== "MALAWI" ~ "CXN1ZRF1", 
    TRUE ~ hmst_specid)) %>%

# 6.	CXN1S7F1 (3106108) – I will use the st_whole_id (CXN1S7F1) as this is what matches the ELISA SID
  mutate(pid= case_when(
    sid_elisa=="CXN1S7F1" & country== "MALAWI" ~ 3106108,
    TRUE ~ pid)) %>%
  
  mutate(hmst_specid= case_when(
    sid_elisa=="CXN1S7F1" & country== "MALAWI" ~ "CXN1S7F1", 
    TRUE ~ hmst_specid)) %>%

# 7.	CXN1SCF1 (3106151) - I will use the st_whole_id (CXN1SCF1) as this is what matches the ELISA SID
  mutate(pid= case_when(
    sid_elisa=="CXN1SCF1" & country== "MALAWI" ~ 3106151,
  TRUE ~ pid)) %>%
  
  mutate(hmst_specid= case_when(
    sid_elisa=="CXN1SCF1" & country== "MALAWI" ~ "CXN1SCF1", 
    TRUE ~ hmst_specid))

# As of 23 Jan - I discovered that there is no SID in the wholestool dataset for PID 6400178, so for now, we are adding it here. 
ELISA_wholestool <- ELISA_wholestool %>%
  mutate(pid= case_when(
    sid_elisa=="EF0000772" & country== "PERU" ~ 6400178,
    TRUE ~ pid)) 

# Email from Malawi on 13 Feb 24, SID CXN14VF1 = PID 3101265, SID CXN16JF1 = PID 3101937, SID CXN1A0F1  = PID 3103440, SID CXN1LN = PID 3103883
ELISA_wholestool <- ELISA_wholestool %>%
  mutate(pid= case_when(
    sid_elisa=="CXN14VF1" & country== "MALAWI" ~ 3101265,
    TRUE ~ pid)) %>%

  mutate(pid= case_when(
    sid_elisa=="CXN16JF1" & country== "MALAWI"~ 3101937,
    TRUE ~ pid)) %>%
  
  mutate(pid= case_when(
    sid_elisa=="CXN1A0F1" & country== "MALAWI" ~ 3103440,
    TRUE ~ pid)) %>%
  
  mutate(pid= case_when(
    sid_elisa=="CXN1LN" & country== "MALAWI"~ 3103883,
    TRUE ~ pid)) 

# Kenya plate 7 has SID 202060, which didn't link with a PID. In wholestool there is a typo in the spec_whole_id which caused the merge to fail: 2020260
# Adding in PID 2401601 for SID 202060
ELISA_wholestool <- ELISA_wholestool %>%
  mutate(pid= case_when(
    sid_elisa=="202060" & country== "KENYA" & plate==7 ~ 2401601,
    TRUE ~ pid))

# Email from Malawi on 13 May 24
# Plate 4 SID = PID corrections: CXN1M9 = 3104096, CXN1XZ = 3107512, CXN1YH = 3107599, CXN1YK = 3107634, CXN1Z2 = 3107724, CXN1Z7 = 3107751, CXN200 = 3108139
# Plate 5 SID = PID corrections: CXN21U = 3108778, CXN224 = 3108862

ELISA_wholestool <- ELISA_wholestool %>%
  mutate(pid= case_when(
    sid_elisa=="CXN1M9F1" & country== "MALAWI" & plate==4 ~ 3104096,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="CXN1XZF1" & country== "MALAWI" & plate==4 ~ 3107512,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="CXN1YHF1" & country== "MALAWI" & plate==4 ~ 3107599,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="CXN1YKF1" & country== "MALAWI" & plate==4 ~ 3107634,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="CXN1Z2S1" & country== "MALAWI" & plate==4 ~ 3107724,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="CXN1Z7S1" & country== "MALAWI" & plate==4 ~ 3107751,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="CXN200F1" & country== "MALAWI" & plate==4 ~ 3108139,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="CXN21UF1" & country== "MALAWI" & plate==5 ~ 3108778,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="CXN224F1" & country== "MALAWI" & plate==5 ~ 3108862,
    TRUE ~ pid))
  
# Email from Kenya on 16 June 24, The correct SID for PID 2901360 is 202316
ELISA_wholestool <- ELISA_wholestool %>%
mutate(pid= case_when(
  sid_elisa=="202316" & country== "KENYA" & plate==9 ~ 2901360,
  TRUE ~ pid)) 

# Email from Bangladesh on 14 June 24: The correct PID of H1 well is 1509759
ELISA_wholestool <- ELISA_wholestool %>%
  mutate(pid= case_when(
    sid_elisa=="1509749" & country== "BANGLADESH" & plate==14 ~ 1509759,
    TRUE ~ pid)) 

# Enter in the PID for Gambia plate 8 well A2
ELISA_wholestool <- ELISA_wholestool %>%
  mutate(pid= case_when(
    sid_elisa=="7106615ST" & country== "THE GAMBIA" & plate==8 ~ 7106615,
    TRUE ~ pid)) 

# Adding in PIDs for 5 samples from plate 6 from Malawi - there were typos/omissions in spec_whole_id that preventing linkage with PID
ELISA_wholestool <- ELISA_wholestool %>%
  mutate(pid= case_when(
    sid_elisa=="CXN29BF1" & country== "MALAWI" & plate==6 ~ 3110404,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="CXN2ADF1" & country== "MALAWI" & plate==6 ~ 3110772,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="CXN2AEF1" & country== "MALAWI" & plate==6 ~ 3110760,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="CXN2AJF1" & country== "MALAWI" & plate==6 ~ 3111013,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="CXN2E9F1" & country== "MALAWI" & plate==6 ~ 3110822,
    TRUE ~ pid)) 

# there are 4 gambia pids that need to be imported that didn't take on the merge
ELISA_wholestool <- ELISA_wholestool %>%
  mutate(pid= case_when(
    sid_elisa=="7109855ST" & country== "THE GAMBIA" & plate==12 ~ 7109855,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="7110471ST" & country== "THE GAMBIA" & plate==13 ~ 7110471,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="7110848ST" & country== "THE GAMBIA" & plate==13 ~ 7110848,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="7111172ST" & country== "THE GAMBIA" & plate==13 ~ 7111172,
    TRUE ~ pid)) 

# there are 10 SIDs from Malawi with  PIDs that did not populate in the check below 
# 7	B1	CXN2ATF1	3111217; 
# 7	D6	CXN2D5F1	3111807; 
# 7	E4	CXN2DVF1	3111537; 
# 7	C5	CXN2E0F1	3111615
# 8	C2	CXN2F5F1	3112388; 
# 8	G2	CXN2F9F1	 3112429; 
# 8	E3	CXN2FGF1	3112483; 
# 8	C7	CXN2GLF1	3112891; 
# 8	E7	CXN2GQF1	 3112922
# 7	C3	CXN2WBF1	 3111373

ELISA_wholestool <- ELISA_wholestool %>%
  mutate(pid= case_when(
    sid_elisa=="CXN2ATF1" & country== "MALAWI" & plate==7 & well=="B1" ~ 3111217,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="CXN2D5F1" & country== "MALAWI" & plate==7 & well=="D6" ~ 3111807,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="CXN2DVF1" & country== "MALAWI" & plate==7 & well=="E4" ~ 3111537,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="CXN2E0F1" & country== "MALAWI" & plate==7 & well=="C5" ~ 3111615,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="CXN2F5F1" & country== "MALAWI" & plate==8 & well=="C2" ~ 3112388,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="CXN2F9F1" & country== "MALAWI" & plate==8 & well=="G2" ~ 3112429,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="CXN2FGF1" & country== "MALAWI" & plate==8 & well=="E3" ~ 3112483,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="CXN2GLF1" & country== "MALAWI" & plate==8 & well=="C7" ~ 3112891,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="CXN2GQF1" & country== "MALAWI" & plate==8 & well=="E7" ~ 3112922,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="CXN2WBF1" & country== "MALAWI" & plate==7 & well=="C3" ~ 3111373,
    TRUE ~ pid))

# there are typos in spec_whole_id for the peru sites that kept the PID data from merging
ELISA_wholestool <- ELISA_wholestool %>%
  mutate(pid= case_when(
    sid_elisa=="EF0002357" & country== "PERU" & plate==2 & well=="H3" ~ 6301149,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="EF0002495" & country== "PERU" & plate==2 & well=="A5" ~ 6301247,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="EF0003185" & country== "PERU" & plate==2 & well=="D9" ~ 6500972,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="EF0004706" & country== "PERU" & plate==4 & well=="B4" ~ 6400824,
    TRUE ~ pid)) %>%
  mutate(pid= case_when(
    sid_elisa=="EF0005208" & country== "PERU" & plate==4 & well=="D9" ~ 6501568,
    TRUE ~ pid)) #%>%
 # mutate(pid= case_when(
 #   sid_elisa=="EF0005366" & country== "PERU" & plate==5 & well=="A1" ~ 6302484,
 #   TRUE ~ pid))

# there is a typo in spec_whole_id that kept a PID from merging from malawi 
ELISA_wholestool <- ELISA_wholestool %>%
  mutate(pid= case_when(
    sid_elisa=="CXN2BIF1" & country== "MALAWI" & plate==8 & well=="A1" ~ 3112264,
    TRUE ~ pid))

# we need to remove sid EF0005366 from the peru site- there is no corresponding pid in the wholestool dataset - we are cutting our losses 
ELISA_wholestool <- ELISA_wholestool %>%
  filter(sid_elisa!="EF0005366")

#################################################################
# Standardizing "checking" the merged ELISA Whole Stool dataset #
#################################################################

# step 1 : filter to missing PIDs
step1_ELISA_wholestool_fix <- ELISA_wholestool %>%
  filter(sid_elisa!="CONTROL1" & sid_elisa!="CONTROL2" & sid_elisa!="STANDARD1" & sid_elisa!= "STANDARD2" & sid_elisa!="STANDARD3" & sid_elisa!= "STANDARD4" & sid_elisa!="STANDARD5" & sid_elisa!="STANDARD6") %>%
  filter(is.na(pid))

# first select only variables of interest and only distinct rows
pid_sid <- ELISA_wholestool %>%
  select(country, pid, sid_elisa) %>%
  filter(sid_elisa!="CONTROL1" & sid_elisa!="CONTROL2" & sid_elisa!="STANDARD1" & sid_elisa!= "STANDARD2" & sid_elisa!="STANDARD3" & sid_elisa!= "STANDARD4" & sid_elisa!="STANDARD5" & sid_elisa!="STANDARD6") %>%
  distinct()

# step 2: filter to duplicate PIDs
step2_ELISA_wholestool_fix <- pid_sid %>%
  filter(duplicated(pid) | duplicated(pid, fromLast = TRUE)) 

# old code that didn't seem to work
# step2_ELISA_wholestool_fix <- pid_sid %>%
#   mutate(duplicate_sid_elisa = duplicated(sid_elisa) | duplicated(sid_elisa, fromLast = TRUE)) %>%
#   filter(duplicate_sid_elisa=="TRUE")

# ---------------------------------------------------------------------------- #

#############################################
# Merging dysentery variable to the dataset #
#############################################

# On 11 Sept - adding in dysentery variable
# blood <- read.csv("dysentery variable.csv") # old marker

# merge in blood data to bind_wholestool data - removes standards and controls
# ELISA_wholestool <- merge(ELISA_wholestool, blood, by = "pid", all.x = TRUE)

#########################
# Creating PID/SID list #
#########################

# grabbing PIDs to send to Erika 
pid_sid <- ELISA_wholestool %>%
  select(country, pid, sid_elisa) %>%
  filter(sid_elisa!="CONTROL1" & sid_elisa!="CONTROL2" & sid_elisa!="STANDARD1" & sid_elisa!= "STANDARD2" & sid_elisa!="STANDARD3" & sid_elisa!= "STANDARD4" & sid_elisa!="STANDARD5" & sid_elisa!="STANDARD6") %>%
  distinct()

check_kenya <- pid_sid %>% 
  filter(country=="KENYA")

####################################
# CHECK ON KENYA DATA - 11 July 24 #
####################################
Kenya_check <- ELISA_wholestool %>%
  filter(country=="KENYA") %>%
  filter(sid_elisa!="CONTROL1" & sid_elisa!="CONTROL2" & sid_elisa!="STANDARD1" & sid_elisa!= "STANDARD2" & sid_elisa!="STANDARD3" & sid_elisa!= "STANDARD4" & sid_elisa!="STANDARD5" & sid_elisa!="STANDARD6") %>%
  group_by(biomarker) %>%
  summarise(count = n())

kenya_check1 <- check_kenya %>%
mutate(duplicate_sid_elisa = duplicated(sid_elisa) | duplicated(sid_elisa, fromLast = TRUE))

##########################
# Mothly Meeting Metrics #
##########################

# BG has 3402 obs, to check this, 94 well plates * 3 (MPO, HAEM, CAL) + 96 = 378 * 9 plates = 3402
# But this would be 80*9 = 720 children

# distinct number of plates per biomarker by site (the # of plates for each biomarker uploaded by country
country_counts <- ELISA_wholestool %>%
  select(country, biomarker, plate) %>%
  distinct()
 
country_counts1 <- country_counts %>%
 group_by(country, biomarker) %>%
  summarise(num_plates=length(biomarker)) 

# number of children per site 
table(pid_sid$country)

# mean log10 biomarker concentration by site 
mean_conc <- ELISA_wholestool %>%
  filter(sid_elisa!="CONTROL1" & sid_elisa!="CONTROL2" & sid_elisa!="STANDARD1" & sid_elisa!= "STANDARD2" & sid_elisa!="STANDARD3" & sid_elisa!= "STANDARD4" & sid_elisa!="STANDARD5" & sid_elisa!="STANDARD6") %>%
  select(country, biomarker, corrected_concentration) 

# checking NAs for conc (if the conc is too high, NAs will be produced)
mean_conc_NA <- mean_conc %>%
  filter(is.na(corrected_concentration))

# back to mean log10 biomarker concentration by site 
mean_conc1 <- mean_conc %>%
  filter(!is.na(corrected_concentration)) %>% # filter out the NAs
  mutate(log10_conc = log10(corrected_concentration)) %>%
  group_by(country, biomarker) %>%
  summarise(mean(log10_conc))

mean_conc2 <- mean_conc %>%
  filter(!is.na(corrected_concentration)) %>% # filter out the NAs
  group_by(country, biomarker) %>%
  summarise(mean(corrected_concentration))

# filter to the flag where "concentration is higher than the highest standard"
high_conc <- ELISA_wholestool %>%
  filter(flag=="Concentration higher than highest standard") 

high_conc1 <- high_conc %>%
group_by(country, biomarker) %>%
  summarise(num_high_conc=length(flag)) 

# getting number of SIDs per country
pid_sid1 <- pid_sid %>%
  group_by(country) %>%
  summarise(sid_total=length(sid_elisa)) 

# merging high_conc1 and pid_sid1 to get total # of samples they've tested into a percentage
sid_high_conc <- merge(high_conc1, pid_sid1, by= "country")

# total # of samples they've tested into a percentage
sid_high_conc <- sid_high_conc %>%
  mutate(high_conc_perc = (num_high_conc/sid_total)*100) %>%
  mutate(high_conc_perc=round(high_conc_perc,2))

#################
# REPORT BY SID #
#################

# in the event that a site does not a complete plate series - we need a break down on number of children per plate per biomarker
sid_report <- ELISA_wholestool %>%
  filter(sid_elisa!="CONTROL1" & sid_elisa!="CONTROL2" & sid_elisa!="STANDARD1" & sid_elisa!= "STANDARD2" & sid_elisa!="STANDARD3" & sid_elisa!= "STANDARD4" & sid_elisa!="STANDARD5" & sid_elisa!="STANDARD6") %>%
  select(country, biomarker, plate, sid_elisa) %>%
#  filter(!is.na(pid)) %>% # we have some missing PIDs for BG and removing them here to get accurate numbers
  distinct()

sid_report1 <- sid_report %>%
  group_by(country, biomarker, plate) %>%
  summarise(num_samples=length(sid_elisa)) 

sid_report2 <- sid_report1 %>%
  group_by(country, biomarker) %>%
  summarise(total_samples=sum(num_samples)) 

# merge in country_counts1 with test2
sid_report3 <- merge(country_counts1, sid_report2, by = c("country", "biomarker"), all.x = T)

# merge in high_conc1
sid_report4 <- merge(sid_report3, high_conc1, by = c("country", "biomarker"), all.x = T)

# make NAs = 0
sid_report4 <- sid_report4 %>% mutate(num_high_conc = ifelse(is.na(num_high_conc),0,num_high_conc))

# total # of samples they've tested into a percentage
sid_report5 <- sid_report4 %>%
  mutate(high_conc_perc = (num_high_conc/total_samples)*100) %>%
  mutate(high_conc_perc=round(high_conc_perc,2))

# getting total number of samples per site
sid_report6 <- sid_report %>%
  group_by(country) %>%
  summarise(num_samples=length(sid_elisa)) 

# and then geting number of high conc by site
sid_report7 <- high_conc %>%
  group_by(country) %>%
  summarise(num_high_conc=length(flag)) 

# and then merging
sid_report8 <- merge(sid_report6, sid_report7, by="country")

# total # of samples they've tested into a percentage
sid_report8 <- sid_report8 %>%
  mutate(high_conc_perc = (num_high_conc/num_samples)*100) %>%
  mutate(high_conc_perc=round(high_conc_perc,2))

#################
# REPORT BY PID #
#################

# in the even that a site does not have a complete plate series - we need a break down on number of children per plate per biomarker
pid_report <- ELISA_wholestool %>%
  filter(sid_elisa!="CONTROL1" & sid_elisa!="CONTROL2" & sid_elisa!="STANDARD1" & sid_elisa!= "STANDARD2" & sid_elisa!="STANDARD3" & sid_elisa!= "STANDARD4" & sid_elisa!="STANDARD5" & sid_elisa!="STANDARD6") %>%
  select(country, biomarker, plate, pid) %>%
  filter(!is.na(pid)) %>% # we have some missing PIDs for BG and removing them here to get accurate numbers
  distinct()

pid_report1 <- pid_report %>%
  group_by(country, biomarker, plate) %>%
  summarise(num_kids=length(pid)) 

pid_report2 <- pid_report1 %>%
  group_by(country, biomarker) %>%
  summarise(total_kids=sum(num_kids)) 

# merge in country_counts1 with test2
pid_report3 <- merge(country_counts1, pid_report2, by = c("country", "biomarker"), all.x = T)

# merge in high_conc1
pid_report4 <- merge(pid_report3, high_conc1, by = c("country", "biomarker"), all.x = T)

# make NAs = 0
pid_report4 <- pid_report4 %>% mutate(num_high_conc = ifelse(is.na(num_high_conc),0,num_high_conc))

# total # of samples they've tested into a percentage
pid_report5 <- pid_report4 %>%
  mutate(high_conc_perc = (num_high_conc/total_kids)*100) %>%
  mutate(high_conc_perc=round(high_conc_perc,2))


###########################
# List of repeats by site #
###########################

bangladesh_repeats <- ELISA_wholestool %>%
  filter(flag=="Concentration higher than highest standard" | flag=="Concentration lower than lowest standard" | flag=="Unable to generate concentration despite OD being in range")%>%
  filter(country=="BANGLADESH") %>%
  arrange(plate)

kenya_repeats <- ELISA_wholestool %>%
  filter(flag=="Concentration higher than highest standard" | flag=="Concentration lower than lowest standard" | flag=="Unable to generate concentration despite OD being in range")%>%
# (well=="H8" & plate=="9" & biomarker=="NGAL" & country=="KENYA")) %>% # added in where the OD was unreadable for one sample - No longer needed after being re-run
  filter(country=="KENYA") %>%
  arrange(plate)

malawi_repeats <- ELISA_wholestool %>%
  filter(flag=="Concentration higher than highest standard" | flag=="Concentration lower than lowest standard" | flag=="Unable to generate concentration despite OD being in range") %>%
  filter(country=="MALAWI") %>%
  arrange(plate)

pakistan_repeats <- ELISA_wholestool %>%
  filter(flag=="Concentration higher than highest standard" | flag=="Concentration lower than lowest standard" | flag=="Unable to generate concentration despite OD being in range") %>%
  filter(country=="PAKISTAN") %>%
  arrange(plate)

gambia_repeats <- ELISA_wholestool %>%
  filter(flag=="Concentration higher than highest standard" | flag=="Concentration lower than lowest standard" | flag=="Unable to generate concentration despite OD being in range") %>%
  filter(country=="THE GAMBIA") %>%
  arrange(plate)

peru_repeats <- ELISA_wholestool %>%
  filter(flag=="Concentration higher than highest standard" | flag=="Concentration lower than lowest standard" | flag=="Unable to generate concentration despite OD being in range") %>%
  filter(country=="PERU") %>%
  arrange(plate)

 write.csv(bangladesh_repeats, file="C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/bangladesh_updatedmodel_newrepeats.csv")
 write.csv(kenya_repeats, file="C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/kenya_repeats_updatedmodel_newrepeats.csv")
 write.csv(malawi_repeats, file="C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/malawi_repeats_updatedmodel_newrepeats.csv")
 write.csv(pakistan_repeats, file="C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/pakistan_repeats_updatedmodel_newrepeats.csv")
 write.csv(gambia_repeats, file="C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/gambia_repeats_updatedmodel_newrepeats.csv")
 write.csv(peru_repeats, file="C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/peru_repeats_updatedmodel_newrepeats.csv")

#############################################################
# List of standards by plate, biomarker, standard, and site #
#############################################################

bangladesh_standards <- ELISA_wholestool %>%
  filter(sid_elisa=="STANDARD1" | sid_elisa=="STANDARD2" | sid_elisa=="STANDARD3" | sid_elisa=="STANDARD4" | sid_elisa=="STANDARD5" |sid_elisa=="STANDARD6") %>%
  filter(country=="BANGLADESH") %>%
  group_by(plate, biomarker, sid_elisa) %>%
  summarize(avg_od = mean(od, na.rm = TRUE), .groups = 'drop') %>%
  arrange(plate)

kenya_standards <- ELISA_wholestool %>%
  filter(sid_elisa=="STANDARD1" | sid_elisa=="STANDARD2" | sid_elisa=="STANDARD3" | sid_elisa=="STANDARD4" | sid_elisa=="STANDARD5" |sid_elisa=="STANDARD6") %>%
  filter(country=="KENYA") %>%
  group_by(plate, biomarker, sid_elisa) %>%
  summarize(avg_od = mean(od, na.rm = TRUE), .groups = 'drop') %>%
  arrange(plate)

malawi_standards <- ELISA_wholestool %>%
  filter(sid_elisa=="STANDARD1" | sid_elisa=="STANDARD2" | sid_elisa=="STANDARD3" | sid_elisa=="STANDARD4" | sid_elisa=="STANDARD5" |sid_elisa=="STANDARD6") %>%
  filter(country=="MALAWI") %>%
  group_by(plate, biomarker, sid_elisa) %>%
  summarize(avg_od = mean(od, na.rm = TRUE), .groups = 'drop') %>%
  arrange(plate)

pakistan_standards <- ELISA_wholestool %>%
  filter(sid_elisa=="STANDARD1" | sid_elisa=="STANDARD2" | sid_elisa=="STANDARD3" | sid_elisa=="STANDARD4" | sid_elisa=="STANDARD5" |sid_elisa=="STANDARD6") %>%
  filter(country=="PAKISTAN") %>%
  group_by(plate, biomarker, sid_elisa) %>%
  summarize(avg_od = mean(od, na.rm = TRUE), .groups = 'drop') %>%
  arrange(plate)

gambia_standards <- ELISA_wholestool %>%
  filter(sid_elisa=="STANDARD1" | sid_elisa=="STANDARD2" | sid_elisa=="STANDARD3" | sid_elisa=="STANDARD4" | sid_elisa=="STANDARD5" |sid_elisa=="STANDARD6") %>%
  filter(country=="THE GAMBIA") %>%
  group_by(plate, biomarker, sid_elisa) %>%
  summarize(avg_od = mean(od, na.rm = TRUE), .groups = 'drop') %>%
  arrange(plate)

peru_standards <- ELISA_wholestool %>%
  filter(sid_elisa=="STANDARD1" | sid_elisa=="STANDARD2" | sid_elisa=="STANDARD3" | sid_elisa=="STANDARD4" | sid_elisa=="STANDARD5" |sid_elisa=="STANDARD6") %>%
  filter(country=="PERU") %>%
  group_by(plate, biomarker, sid_elisa) %>%
  summarize(avg_od = mean(od, na.rm = TRUE), .groups = 'drop') %>%
  arrange(plate)

# write.csv(bangladesh_standards, file="C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/bangladesh_standards_FINAL.csv")
# write.csv(kenya_standards, file="C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/kenya_standards_FINAL.csv")
# write.csv(malawi_standards, file="C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/malawi_standards_FINAL.csv")
# write.csv(pakistan_standards, file="C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/pakistan_standards_FINAL.csv")
# write.csv(gambia_standards, file="C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/gambia_standards_FINAL.csv")
# write.csv(peru_standards, file="C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/peru_standards_FINAL.csv")

#####################
# Saving dataframes #
#####################

# save dataframes
write.csv(pid_sid, file="C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/PID and SID lists/PIDs_SIDs_FINAL.csv")
#write.csv(sid, file="C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/SIDs_11Oct2023.csv")
#write.csv(pid, file="C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/PIDs_11Oct2023.csv")

# data to be merged into TAC and culture data from Erika - THERE ARE NO REPEATS IN THIS DATASET
save(ELISA_wholestool, file = "C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/R output/elisa_wholestool - dysentery - updatedmodel - 12Jan2026.rda")

# save the repeat dataframe
save(repeat_df, file = "C:/Users/sb2yb/Box/ELISA Inflammatory Markers/EFGH Analyses/R output/elisa_repeats - updatedmodel - 12Jan2026.rda")
