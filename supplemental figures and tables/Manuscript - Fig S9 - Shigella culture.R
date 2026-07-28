
# a boxplot characterizing the distribution of biomarker concentration for Shigella based on culture+/-. 

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# to prevent scientific notation on the plot
options(scipen=999)

library(dplyr)
library(tidyr)
library(ggplot2)
library(gee)
library(geepack)

# load in ELISA and TAC data 
load("R output/EFGH_substudy_rowperpid_ALLDATA_LOCKEDv3.Rda")

# keep necessary variables (pid, country, biomarkers, sw_isolate)
df <- EFGH_substudy_rowperpid_ALLDATA %>%
  select(pid, child_id, country, clean_conc_log10_CAL, clean_conc_log10_HAEM, clean_conc_log10_MPO, clean_conc_log10_NGAL, culture_shigella_positive)

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

# filter out where shigella culture is NA
long_df <- long_df %>%
  filter(!is.na(culture_shigella_positive))

# make shigella culture variable a factor
long_df$culture_shigella_positive <- as.factor(long_df$culture_shigella_positive)

# create new yes/no variable so that the legend works
long_df <- long_df %>%
  mutate(shigella_culture = case_when(
    culture_shigella_positive == "0" ~ "No",
    culture_shigella_positive == "1" ~ "Yes"))

# make shigella culture variable a factor
long_df$shigella_culture <- as.factor(long_df$shigella_culture)

# set up biomarker levels
long_df$biomarker <- factor(long_df$biomarker, levels = c("HAEM", "NGAL", "MPO","CAL"), 
                              labels = c("Hemoglobin", "Lipocalin-2", "Myeloperoxidase","Calprotectin"))

# df per biomarker
CAL_df <- long_df %>% filter(biomarker=="Calprotectin")
HAEM_df <- long_df %>% filter(biomarker=="Hemoglobin")
MPO_df <- long_df %>% filter(biomarker=="Myeloperoxidase")
NGAL_df <- long_df %>% filter(biomarker=="Lipocalin-2")

# --- Run GEE per biomarker ---

  # function that will keep the log10 and antilog values that were run through the model 
  risk_func <-function(dat) {
    model <- geeglm(log10_conc ~ shigella_culture, family = gaussian, data = dat, id = child_id, corstr = "exchangeable")
    beta <- summary(model)$coefficients[2,1]
    lb <- summary(model)$coefficients[2,1] - 1.96 * summary(model)$coefficients[2,2]
    ub <- summary(model)$coefficients[2,1] + 1.96 * summary(model)$coefficients[2,2]
    pval_num <- summary(model)$coefficients[2, 4]
    pval_str <- sprintf("%.5f", pval_num)
    
    
    result <- data.frame(
      beta = beta,
      lb = lb,
      ub = ub,
      p.value = pval_str,
      p.value_num = pval_num)
    
    return(result)
  }
  

# call the function
# CAL
CAL_culture <- risk_func(dat=CAL_df)
HAEM_culture <- risk_func(dat=HAEM_df)
MPO_culture <- risk_func( dat=MPO_df)
NGAL_culture <- risk_func(dat=NGAL_df)


# --- Combine results ---
gee_results <- data.frame(
  biomarker = c("Calprotectin", "Hemoglobin", "Myeloperoxidase", "Lipocalin-2"),
  rbind(CAL_culture, HAEM_culture, MPO_culture, NGAL_culture))

# --- Add significance stars ---
gee_results <- gee_results %>%
  mutate(label = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE ~ ""))

# --- Prepare for plotting significance bars ---
pvals_manual <- long_df %>%
  group_by(biomarker) %>%
  summarise(y_max = max(log10_conc, na.rm = TRUE)) %>%
  left_join(gee_results, by = "biomarker") %>%
  mutate(
    biomarker_num = as.numeric(factor(biomarker, levels = c("Hemoglobin", "Lipocalin-2", "Myeloperoxidase","Calprotectin"))),
    xmin_num = biomarker_num - 0.15,
    xmax_num = biomarker_num + 0.15,
    y.position = y_max + 0.3)


# Figure - Shigella culture
plot = ggplot(long_df, aes(x=biomarker, y=log10_conc, fill=shigella_culture)) 
plot = plot + theme_bw()
plot = plot + geom_boxplot() 
plot = plot + 
  labs(
    x = "",
    y = "Log10(biomarker concentration)")
plot = plot + scale_fill_manual(values=c("#AC1408","#005D82"), 
                                   name="Shigella \nculture positive",
                                breaks=c("No", "Yes"))
plot = plot + theme(legend.position = c(.85, .25))
plot = plot + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank())
plot = plot +  theme(axis.title=element_text(size=16))
plot = plot +  theme(axis.text=element_text(size=12))
plot = plot + theme(legend.title = element_text(size = 14),  # Adjust title font size
                    legend.text = element_text(size = 12))   # Adjust labels font size
plot = plot +  theme(axis.line = element_line(color="black"))
plot = plot +  scale_y_continuous(limits = c(-5.5, 7), breaks = c(-5, -2.5, 0, 2.5, 5, 7.5))


# Add significance bars and labels
plot = plot +
  geom_segment(data = pvals_manual,
               aes(x = xmin_num, xend = xmax_num, y = y.position, yend = y.position),
               inherit.aes = FALSE, color = "black") +
  geom_text(data = pvals_manual,
            aes(x = biomarker_num, y = y.position + 0.05, label = label),
            inherit.aes = FALSE, size = 5)

plot

# Save Figure - Shigella Ct 31.1
ggsave(file="Figures/Figure S7 - Shigella Culture - FINAL.pdf",plot,width=8,height=6,dpi=300)

