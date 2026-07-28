
# density plots by dysentery variable: ANY reported blood

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# to prevent scientific notation on the plot
options(scipen=999)

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

# load in ELISA and TAC data 
load("R output/EFGH_substudy_rowperpid_ALLDATA_LOCKEDv3.Rda")

# visualize my data
table(EFGH_substudy_rowperpid_ALLDATA$clean_conc_log10_HAEM, useNA = "always")  
table(EFGH_substudy_rowperpid_ALLDATA$dysentery, useNA = "always")  


###################
# set up for plot #
###################

plot1_df <- EFGH_substudy_rowperpid_ALLDATA %>%
  filter(!is.na(dysentery)) %>%
  filter(!is.na(clean_conc_log10_HAEM)) 

# check that NAs were removed
table(plot1_df$dysentery, useNA = "always")
table(plot1_df$clean_conc_log10_HAEM, useNA = "always")


########
# plot #
########

  # this is our plot
  blood_plot <- ggplot(plot1_df, aes(x = clean_conc_log10_HAEM, fill=factor(dysentery)))
  blood_plot = blood_plot + geom_density(alpha=0.2)
  blood_plot = blood_plot + theme_bw()
  blood_plot = blood_plot + 
    labs(x= "Log10 concentration of hemoglobin",
         y= "Density",
         title = "")
  blood_plot = blood_plot + theme(panel.grid.major = element_blank())
  blood_plot = blood_plot + theme(panel.grid.minor = element_blank())
  blood_plot = blood_plot + scale_fill_manual(values=c("#AC1408", "#09B9FF"),
                                              name=c(""),
                                              label=c("No visible blood in stool", "Visible blood in stool"),
                                              breaks = c(0,1))
  blood_plot = blood_plot + theme(legend.position = c(0.90, 0.80))
  blood_plot = blood_plot + theme(axis.title=element_text(size=16)) # x and y axes titles
  blood_plot = blood_plot + theme(axis.text=element_text(size=14)) # numbers on the axes
  blood_plot = blood_plot + theme(title=element_text(size=16)) # numbers on the axes
  blood_plot = blood_plot + xlim(-7, 4) # Set x-axis range to go from 0 to 5

  blood_plot

  ggsave(file="Figures/Figure 2 - HAEM density - Any reported - FINAL.pdf",blood_plot,width=10,height=6,dpi=300)
  