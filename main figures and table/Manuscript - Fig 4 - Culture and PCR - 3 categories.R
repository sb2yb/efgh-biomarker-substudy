
# a boxplot characterizing the distribution of biomarker log10 concentration for Shigella depending on if the sample was culture/PCR +/-.
# IN THIS SCRIPT THE CULTURE POSITIVE GROUPS ARE COMBINED

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# to prevent scientific notation on the plot
options(scipen=999)

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(purrr)
library(gee)
library(geepack)

# load in ELISA and TAC data 
load("R output/EFGH_substudy_rowperpid_ALLDATA_LOCKEDv3.Rda")

# keep necessary variables (pid, country, biomarkers, shigella_attributable)
df <- EFGH_substudy_rowperpid_ALLDATA %>%
  select(pid, country, child_id, clean_conc_log10_CAL, clean_conc_log10_HAEM, clean_conc_log10_MPO, clean_conc_log10_NGAL, shigella_attributable, culture_shigella_positive)

# check out the data
table(df$shigella_attributable, useNA = "always")
table(df$culture_shigella_positive, useNA = "always")
table(df$shigella_attributable, df$culture_shigella_positive, useNA = "always")

# remove NAs for shigella_attributable and culture_shigella_positive
df <- df %>%
  filter(!is.na(shigella_attributable) & !is.na(culture_shigella_positive))

# create categories 
df <- df %>%
  mutate(culture_pcr_category = case_when(
    culture_shigella_positive == 1 ~ "culture+",
    culture_shigella_positive == 0 & shigella_attributable == 1 ~ "culture-/PCR+",
    culture_shigella_positive == 0 & shigella_attributable == 0 ~ "culture-/PCR-"))

# how we want the plot to look 
df$culture_pcr_category <- factor(df$culture_pcr_category, levels = c("culture-/PCR-", "culture-/PCR+", "culture+"))


# transform the dataframe so that there is only 1 column with the biomarker concentrations and 1 column indicating the name of the biomarker
long_df <- df %>%
  pivot_longer(
    cols = starts_with("clean_conc_log10_"),
    names_to = "biomarker",
    values_to = "log10_conc") %>% # Extract the biomarker names from the column names
  mutate(biomarker = gsub("clean_conc_log10_", "", biomarker))

# filter out where concentrations are NA
long_df <- long_df %>%
  filter(!is.na(log10_conc))

# set up biomarker levels
long_df$biomarker <- factor(long_df$biomarker, levels = c("HAEM", "NGAL", "MPO","CAL"), 
                                    labels = c("Hemoglobin", "Lipocalin-2", "Myeloperoxidase","Calprotectin"))


# df per biomarker
CAL_df <- long_df %>% filter(biomarker=="Calprotectin")
HAEM_df <- long_df %>% filter(biomarker=="Hemoglobin")
MPO_df <- long_df %>% filter(biomarker=="Myeloperoxidase")
NGAL_df <- long_df %>% filter(biomarker=="Lipocalin-2")


#### checking the model 
# set referent group in the model - we don't want to mess with the levels above because it changes our plot 
model <- geeglm(log10_conc ~ relevel(culture_pcr_category, ref = "culture-/PCR+"), family = gaussian, data = NGAL_df, id = child_id, corstr = "exchangeable")
beta <- summary(model)$coefficients[2:3,1]
lb <- summary(model)$coefficients[2:3,1] - 1.96 * summary(model)$coefficients[2:3,2]
ub <- summary(model)$coefficients[2:3,1] + 1.96 * summary(model)$coefficients[2:3,2]
pval_num <- summary(model)$coefficients[2:3, 4]
pval_str <- sprintf("%.5f", pval_num)

summary(model)


# --- Run GEE per biomarker ---

# function that will keep the log10 and antilog values that were run through the model 
risk_func <-function(dat, biomarker) {
  model <- geeglm(log10_conc ~ relevel(culture_pcr_category, ref = "culture-/PCR+"), family = gaussian, data = dat, id = child_id, corstr = "exchangeable")
  beta <- summary(model)$coefficients[2:3,1]
  lb <- summary(model)$coefficients[2:3,1] - 1.96 * summary(model)$coefficients[2:3,2]
  ub <- summary(model)$coefficients[2:3,1] + 1.96 * summary(model)$coefficients[2:3,2]
  pval_num <- summary(model)$coefficients[2:3, 4]
  pval_str <- sprintf("%.5f", pval_num)
  
  # Assign group names
  group <- c("culture-/PCR-", "culture+")
  
  result <- data.frame(
    biomarker = biomarker,
    group = group,
    beta = beta,
    lb = lb,
    ub = ub,
    p.value = pval_str,
    p.value_num = pval_num)
  
  return(result)
}


# call the function
CAL_cat <- risk_func(dat=CAL_df, "Calprotectin")
HAEM_cat <- risk_func(dat=HAEM_df, "Hemoglobin")
MPO_cat <- risk_func( dat=MPO_df, "Myeloperoxidase")
NGAL_cat <- risk_func(dat=NGAL_df, "Lipocalin-2")


# --- Combine results ---
gee_results <- data.frame(rbind(CAL_cat, HAEM_cat, MPO_cat, NGAL_cat))

# --- Add significance stars ---
gee_results <- gee_results %>%
  mutate(label = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE ~ ""))

# Prepare bar positions for significance using gee_results
pvals_manual <- gee_results %>%
  left_join(
    long_df %>% group_by(biomarker) %>% summarise(y_max = max(log10_conc, na.rm = TRUE)),
    by = "biomarker") %>%
  mutate(
    biomarker_num = as.numeric(factor(biomarker, levels = c("Hemoglobin", "Lipocalin-2", "Myeloperoxidase","Calprotectin"))),
    # Set line start/end for two comparisons vs reference (culture-/PCR+)
    xmin_num = case_when(
      group == "culture-/PCR-" ~ biomarker_num - 0.25,
      group == "culture+"      ~ biomarker_num + 0.05),
    xmax_num = case_when(
      group == "culture-/PCR-" ~ biomarker_num - 0.05,
      group == "culture+"      ~ biomarker_num + 0.25),
    y.position = y_max + 0.3 + row_number() * 0.15)  # stagger lines to avoid overlap
  

# Build figure in your preferred plot = plot + format
plot = ggplot(long_df, aes(x=biomarker, y=log10_conc, fill=culture_pcr_category)) 
plot = plot + theme_bw()
plot = plot + geom_boxplot(position = position_dodge(width = 0.8)) 
plot = plot + labs(x = "", y = "Log10(biomarker concentration)")
plot = plot + scale_fill_manual(values=c("#AC1408","#FFC107","#005D82"), 
                                name="",
                                breaks=c("culture-/PCR-", "culture-/PCR+", "culture+"))
plot = plot + theme(legend.position = c(.85, .25))
plot = plot + theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank())
plot = plot +  theme(axis.title=element_text(size=16))
plot = plot +  theme(axis.text=element_text(size=14))
plot = plot +  theme(axis.line = element_line(color="black"))
plot = plot + theme(legend.text = element_text(size = 14))   
plot = plot +  scale_y_continuous(limits = c(-5.5, 7.5), breaks = c(-5, -2.5, 0, 2.5, 5, 7.5))

# significance bars with offset
plot = plot + geom_segment(data = pvals_manual,
                           aes(x = xmin_num, xend = xmax_num, y = y.position, yend = y.position),
                           inherit.aes = FALSE, color="black")
plot = plot + geom_text(data = pvals_manual,
                        aes(x = (xmin_num + xmax_num)/2, y = y.position + 0.1, label = label),
                        inherit.aes = FALSE, size = 5)

# Display figure
plot


# Save Figure - Shigella Ct 31.1
ggsave(file="Figures/Figure 4 - Shigella culture and PCR positivity -FINAL.pdf",plot,width=8,height=6,dpi=300)
