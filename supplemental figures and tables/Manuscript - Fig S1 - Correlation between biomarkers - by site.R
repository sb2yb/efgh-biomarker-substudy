
# a scatterplot characterizing the relationship between the log10 concentrations of each of the 4 biomarkers - by site

# There will be 6 plots showing the interactions between each of the 4 biomarkers
#    - CAL x MPO
#    - CAL x HAEM
#    - CAL x NGAL
#    - MPO x HAEM
#    - MPO x NGAL
#    - HAEM x NGAL

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# to prevent scientific notation on the plot
options(scipen=999)

library(dplyr)
library(tidyr)
library(xlsx)
library(ggplot2)
library(ggpubr)
library(stringr)

# load in elisa data
load("R output/EFGH_substudy_rowperpid_ALLDATA_LOCKEDv3.Rda")

# keeping only the needed variables
biomarkers <- EFGH_substudy_rowperpid_ALLDATA %>%
  select(pid, country, clean_conc_log10_CAL, clean_conc_log10_MPO, clean_conc_log10_HAEM, clean_conc_log10_NGAL)%>%
  filter(!is.na(pid)) # filter out any missing PIDs

# we need to re-configure the dataframe where there is only 1 row per pid 
# ifelse(): Added to handle cases where all values are NA for a given group. It checks if all values are NA and returns NA instead of -Inf.
#       If all values in clean_conc_log10_CAL are NA for a group, the new column clean_conc_log10_CAL will be NA.
#       If not all values are NA, it calculates the maximum of the non-NA values in clean_conc_log10_CAL and assigns that as the value.
#       This approach ensures that you only get a meaningful maximum value if there are non-missing values, and it avoids returning -Inf when all values are NA.
df_consolidated <- biomarkers %>%
  group_by(pid, country) %>%
  summarize(
    clean_conc_log10_CAL = ifelse(all(is.na(clean_conc_log10_CAL)), NA, max(clean_conc_log10_CAL, na.rm = TRUE)),
    clean_conc_log10_MPO = ifelse(all(is.na(clean_conc_log10_MPO)), NA, max(clean_conc_log10_MPO, na.rm = TRUE)),
    clean_conc_log10_HAEM = ifelse(all(is.na(clean_conc_log10_HAEM)), NA, max(clean_conc_log10_HAEM, na.rm = TRUE)),
    clean_conc_log10_NGAL = ifelse(all(is.na(clean_conc_log10_NGAL)), NA, max(clean_conc_log10_NGAL, na.rm = TRUE)),
    .groups = 'drop') # Ensures that the final result is not grouped by pid or country, meaning the dataframe is treated as a standard dataframe with no grouping structure.
  

# create a function to loop through the plot code 6 times
scatter <- function(data,biomarker_x,biomarker_y,xtitle,ytitle,country_col, xmin, xmax, ymin, ymax) {

# Figure 
plot = ggplot(data, aes(x=biomarker_x, y=biomarker_y)) 
plot = plot + theme_bw()
plot = plot + geom_point() 
plot = plot + geom_smooth(method='lm', formula= y~x)
plot = plot + stat_cor(method = "spearman", aes(label = ..r.label..), size=4)
plot = plot + facet_wrap(as.formula(paste("~", country_col)))
plot = plot + 
  labs(
    x = xtitle,
    y = ytitle)
plot = plot + theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank())
plot = plot +  theme(axis.title=element_text(size=12))
plot = plot +  theme(axis.text=element_text(size=12))
plot = plot +  theme(axis.line = element_line(color="black"))
plot = plot +  scale_x_continuous(breaks = seq(xmin, xmax, by = 1)) 
plot = plot +  scale_y_continuous(breaks = seq(ymin, ymax, by = 1)) 


return(plot)
}

# call the function 
cal_mpo <- scatter(df_consolidated,df_consolidated$clean_conc_log10_CAL, df_consolidated$clean_conc_log10_MPO,"log10 concentration of calprotectin", "log10 concentration of myeloperoxidase", 'country', 1, 7, -1, 5)
cal_haem <- scatter(df_consolidated,df_consolidated$clean_conc_log10_CAL, df_consolidated$clean_conc_log10_HAEM,"log10 concentration of calprotectin", "log10 concentration of hemoglobin", 'country', 1, 7, -6, 3)
cal_ngal <- scatter(df_consolidated,df_consolidated$clean_conc_log10_CAL, df_consolidated$clean_conc_log10_NGAL,"log10 concentration of calprotectin", "log10 concentration of lipocalin-2", 'country', 1, 7, -3, 3)
mpo_haem <- scatter(df_consolidated,df_consolidated$clean_conc_log10_MPO, df_consolidated$clean_conc_log10_HAEM,"log10 concentration of myeloperoxidase", "log10 concentration of hemoglobin", 'country', -1, 5, -6, 3)
mpo_ngal <- scatter(df_consolidated,df_consolidated$clean_conc_log10_MPO, df_consolidated$clean_conc_log10_NGAL,"log10 concentration of myeloperoxidase", "log10 concentration of lipocalin-2", 'country', -1, 5, -3, 3)
haem_ngal <- scatter(df_consolidated,df_consolidated$clean_conc_log10_HAEM, df_consolidated$clean_conc_log10_NGAL,"log10 concentration of hemoglobin", "log10 concentration of lipocalin-2", 'country', -6, 3, -3, 3)

# visualize the plots
cal_mpo
cal_haem
cal_ngal
mpo_haem
mpo_ngal
haem_ngal

allplots <- ggarrange(cal_mpo, cal_haem, cal_ngal, mpo_haem, mpo_ngal, haem_ngal, nrow=2, ncol=3)
allplots

# Save Figure 
ggsave(file="Figures/Figure S1 - Correlation between biomarkers by site- FINAL.pdf",allplots,width=16,height=8,dpi=300)

