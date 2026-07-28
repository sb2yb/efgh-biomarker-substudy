
# This is the code for Calculating Pathogen Prevalence
# Prevalence: total number of cases / total number of samples tested 
# Prevalence percentage: prevalence * 100

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# to prevent scientific notation on the plot
options(scipen=999)

library(dplyr)
library(tidyr)

# load in ELISA and TAC data 
load("R output/EFGH_substudy_rowperpid_ALLDATA_LOCKEDv3.Rda") 

# keep dichotomous variables
path_prev <- EFGH_substudy_rowperpid_ALLDATA %>%
  select(adenovirus_40_41_detect, aeromonas_detect, astrovirus_detect, c_infans_detect, c_jejuni_coli_detect, cryptosporidium_detect, c_upsalensis_detect,
         cyclospora_detect, e_bieneusi_detect, EAEC_detect, e_histolytica_detect, tEPEC_detect, aEPEC_detect, ST.ETEC_detect, LT.ETEC_detect,
         giardia_detect, h_pylori_detect, isospora_detect, norovirus_gi_detect, norovirus_gii_detect, rotavirus_detect, salmonella_detect, sapovirus_detect,
         v_cholerae_detect, shigella_detect, 
         
         adenovirus_40_41_attributable, aeromonas_attributable, astrovirus_attributable, c_jejuni_coli_attributable, cryptosporidium_attributable,
         cyclospora_attributable, e_histolytica_attributable, tEPEC_attributable, ST.ETEC_attributable, 
         isospora_attributable, norovirus_gii_attributable, rotavirus_attributable, salmonella_attributable, sapovirus_attributable, 
         v_cholerae_attributable, shigella_attributable)


# Function to calculate prevalence for dichotomous variables ending in '_35'
calculate_prevalence <- function(df, suffix) {
  # Select columns ending with _35
  variables <- grep(paste0(suffix, "$"), colnames(df), value = TRUE)

  # Calculate the mean (prevalence) for each specified dichotomous variable
  prevalence_df <- df %>%
    summarise(across(all_of(variables), 
                     ~ mean(., na.rm = TRUE),
                     .names = "prevalence_{.col}")) # renames the variable
  
  # Convert to long format
  prevalence_long_df <- prevalence_df %>%
    pivot_longer(cols = everything(),
                 names_to = "variable",
                 values_to = "prevalence")
  
  return(prevalence_long_df)
}


# Calculate detect and attributable prevalences
prev_detect <- calculate_prevalence(path_prev, "_detect") %>%
  mutate(variable = gsub("_detect", "", variable))

prev_attributable <- calculate_prevalence(path_prev, "_attributable") %>%
  mutate(variable = gsub("_attributable", "", variable))

# Join side-by-side
prevalence_wide <- full_join(prev_detect, prev_attributable, by = "variable", suffix = c("_detect", "_attributable"))

# check that we calculated the above correctly
table(path_prev$v_cholerae_detect, useNA = "always")

# Sort the prevalence results by the 'prevalence' column from highest to lowest
prevalence_results_sorted <- prevalence_wide %>%
  arrange(desc(prevalence_detect)) %>%
  mutate(
    prevalence_detect = prevalence_detect * 100,
    prevalence_attributable = prevalence_attributable * 100,
    prevalence_detect_rounded = round(prevalence_detect, 1),
    prevalence_attributable_rounded = round(prevalence_attributable, 1))


# save dataframe for later reference
write.csv(prevalence_results_sorted, file = "Figures/Table S2 - Pathogen Prevalence - Option 1 Coding - FINAL.csv")

