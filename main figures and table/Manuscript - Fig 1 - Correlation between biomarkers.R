
#  A scatterplot characterizing the relationship between the log10 concentrations of each of the 4 biomarkers

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

# load in elisa data
load("R output/EFGH_substudy_rowperpid_ALLDATA_LOCKEDv3.Rda")

# keeping only the needed variables
biomarkers <- EFGH_substudy_rowperpid_ALLDATA %>%
  select(pid, country, clean_conc_log10_CAL, clean_conc_log10_MPO, clean_conc_log10_HAEM, clean_conc_log10_NGAL)%>%
  filter(!is.na(pid)) # filter out any missing PIDs

# create a function to loop through the plot code 6 times
scatter <- function(data,xaxis,yaxis,xtitle,ytitle) {

# Remove rows where xaxis or yaxis values are NA- not needed but keeps warnings from appearing
data <- data %>%
  filter(!is.na(!!rlang::sym(xaxis)) & !is.na(!!rlang::sym(yaxis))) 
  
# Figure 
plot = ggplot(data, aes(x=!!rlang::sym(xaxis), y=!!rlang::sym(yaxis))) 
plot = plot + theme_bw()
plot = plot + geom_point() 
plot = plot + geom_smooth(method='lm', formula= y~x)
plot = plot + stat_cor(method = "spearman", aes(label = ..r.label..), size=6)
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

return(plot)
}

# call the function 
cal_mpo <- scatter(biomarkers,"clean_conc_log10_CAL", "clean_conc_log10_MPO","log10 concentration of calprotectin", "log10 concentration of myeloperoxidase")
cal_haem <- scatter(biomarkers,'clean_conc_log10_CAL', 'clean_conc_log10_HAEM',"log10 concentration of calprotectin", "log10 concentration of hemoglobin")
cal_ngal <- scatter(biomarkers,'clean_conc_log10_CAL', 'clean_conc_log10_NGAL',"log10 concentration of calprotectin", "log10 concentration of lipocalin-2")
mpo_haem <- scatter(biomarkers,'clean_conc_log10_MPO', 'clean_conc_log10_HAEM',"log10 concentration of myeloperoxidase", "log10 concentration of hemoglobin")
mpo_ngal <- scatter(biomarkers,'clean_conc_log10_MPO', 'clean_conc_log10_NGAL',"log10 concentration of myeloperoxidase", "log10 concentration of lipocalin-2")
haem_ngal <- scatter(biomarkers,'clean_conc_log10_HAEM', 'clean_conc_log10_NGAL',"log10 concentration of hemoglobin", "log10 concentration of lipocalin-2")

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
ggsave(file="Figures/Figure 1 - Correlation between biomarkers - FINAL.pdf",allplots,width=13,height=8,dpi=300)



