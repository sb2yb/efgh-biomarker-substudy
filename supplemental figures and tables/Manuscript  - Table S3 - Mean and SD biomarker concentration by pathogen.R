
# Calculate the mean (SD) log10 concentration for each biomarker 
#       (HAEM, NGAL, MPO, and CAL) for each pathogen (amongst positives only).

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

# keep necessary variables (pid, country, biomarkers, pathogen_attributable where pathogen_35 >5% of samples)
df <- EFGH_substudy_rowperpid_ALLDATA %>%
  select(pid, country, clean_conc_log10_CAL, clean_conc_log10_HAEM, clean_conc_log10_MPO, clean_conc_log10_NGAL, shigella_attributable,
         adenovirus_40_41_attributable, astrovirus_attributable, c_jejuni_coli_attributable, cryptosporidium_attributable, norovirus_gii_attributable, 
         rotavirus_attributable, ST.ETEC_attributable, tEPEC_attributable, sapovirus_attributable)


# Function to calculate mean and standard deviation for biomarker levels
calculate_biomarker_stats <- function(df, biomarker_cols, pathogen_cols) {
  # Convert wide data to long format
  long_df <- df %>%
    pivot_longer(
      cols = all_of(pathogen_cols),
      names_to = "pathogen",
      values_to = "presence")
  
  # Filter for presence of pathogen
  present_df <- long_df %>%
    filter(presence == 1)
  
  # Calculate mean and standard deviation for each biomarker
  result <- present_df %>%
    group_by(pathogen) %>%
    summarise(
      across(
        all_of(biomarker_cols), 
        list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE))
      )
    )
  
  return(result)
}

# Define the biomarker and pathogen column names
biomarker_cols <- c("clean_conc_log10_HAEM","clean_conc_log10_NGAL", "clean_conc_log10_MPO",  "clean_conc_log10_CAL")
pathogen_cols <- c("shigella_attributable", 'tEPEC_attributable', "c_jejuni_coli_attributable", 'ST.ETEC_attributable', "cryptosporidium_attributable",
                   "adenovirus_40_41_attributable", "astrovirus_attributable", "norovirus_gii_attributable", "rotavirus_attributable",  "sapovirus_attributable")

# Calculate the stats
biomarker_stats_per_pathogen <- calculate_biomarker_stats(df, biomarker_cols, pathogen_cols)

# round mean and sd to 2 decimal places
biomarker_stats_per_pathogen <- biomarker_stats_per_pathogen %>%
  mutate(
    # Log10 mean (SD)
    CAL_mean_sd = paste0(round(clean_conc_log10_CAL_mean, 2), " (", round(clean_conc_log10_CAL_sd, 2), ")"),
    HAEM_mean_sd = paste0(round(clean_conc_log10_HAEM_mean, 2), " (", round(clean_conc_log10_HAEM_sd, 2), ")"),
    MPO_mean_sd = paste0(round(clean_conc_log10_MPO_mean, 2), " (", round(clean_conc_log10_MPO_sd, 2), ")"),
    NGAL_mean_sd = paste0(round(clean_conc_log10_NGAL_mean, 2), " (", round(clean_conc_log10_NGAL_sd, 2), ")"),
    
    # Antilog mean (SD)
    CAL_mean_sd_antilog = paste0(round(10^clean_conc_log10_CAL_mean, 2), " (", round(10^clean_conc_log10_CAL_sd, 2), ")"),
    HAEM_mean_sd_antilog = paste0(round(10^clean_conc_log10_HAEM_mean, 2), " (", round(10^clean_conc_log10_HAEM_sd, 2), ")"),
    MPO_mean_sd_antilog = paste0(round(10^clean_conc_log10_MPO_mean, 2), " (", round(10^clean_conc_log10_MPO_sd, 2), ")"),
    NGAL_mean_sd_antilog = paste0(round(10^clean_conc_log10_NGAL_mean, 2), " (", round(10^clean_conc_log10_NGAL_sd, 2), ")"),
    pathogen = factor(pathogen, levels = c(
      "shigella_attributable", "tEPEC_attributable", "c_jejuni_coli_attributable", "ST.ETEC_attributable",
      "cryptosporidium_attributable", "adenovirus_40_41_attributable", "astrovirus_attributable",
      "norovirus_gii_attributable", "rotavirus_attributable", "sapovirus_attributable"))) %>%
    arrange(pathogen)

biomarker_stats_per_pathogen <- biomarker_stats_per_pathogen %>%
  select(pathogen, HAEM_mean_sd, NGAL_mean_sd, MPO_mean_sd, CAL_mean_sd, HAEM_mean_sd_antilog, NGAL_mean_sd_antilog, MPO_mean_sd_antilog, CAL_mean_sd_antilog)


# save dataframes
write.csv(biomarker_stats_per_pathogen, file="Figures/Table S3 - Mean and SD of biomarker concentrations - FINAL.csv")














  
  