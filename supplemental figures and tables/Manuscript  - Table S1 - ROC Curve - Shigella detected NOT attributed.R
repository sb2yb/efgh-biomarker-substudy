
# ROC curve examining the sensitivity/specificity of HAEM biomarker and reported observed blood in stool 
# SHIGELLA DETECTED NOT ATTRIBUTED

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# to prevent scientific notation on the plot
options(scipen=999)

library(dplyr)
library(tidyr)
library(ggplot2)

## ROC curves
library(pROC)
library(cutpointr) # The following objects are masked from ‘package:pROC’: auc, roc

# load in ELISA and TAC data 
load("R output/EFGH_substudy_rowperpid_ALLDATA_LOCKEDv3.Rda")

# keeping only the needed variables
haem <- EFGH_substudy_rowperpid_ALLDATA %>%
  select(pid, clean_conc_log10_HAEM, country, dysentery, shigella_ct, shigella_detect, shigella_attributable) %>%
  filter(!is.na(pid)) %>% # filter out any missing PIDs 
  filter(!is.na(clean_conc_log10_HAEM)) %>% # filter out the concentrations with no values
  filter(!is.na(dysentery)) %>% # remove missing dysentery variables
  filter(shigella_detect==1 & shigella_attributable==0)

# check the class of dysentery
class(haem$dysentery)

# visualize the dysentery data
table(haem$dysentery, useNA = "always")

# convert dysentery to a factor
haem$dysentery <- as.factor(haem$dysentery)

# check the class of clean_conc_log10_HAEM
class(haem$clean_conc_log10_HAEM)

# visualize the dysentery data
table(haem$clean_conc_log10_HAEM, useNA = "always")


haem$dysentery <- as.factor(haem$dysentery)

# ROC curve (one curve per biomarker)
# x = continuous biomarker
# y = dichotomous dysentery
# cutpointr(data, x, y)

# ROC curve plot
# The ROC curve plots the TPR against the FPR at various threshold settings. Each point on the curve corresponds to a different threshold used to classify log_10_biomarker_HAEM as predicting dysentery
all.cut <- cutpointr(haem$clean_conc_log10_HAEM, haem$dysentery, method = maximize_metric, metric = youden, na.rm=TRUE)  %>%
  add_metric(metric=list(ppv, npv))

# Add PPV and NPV metrics
# all.cut <- add_metric(all.cut, metric = list(ppv, npv))

all.cut
plot(all.cut)
rocplot <-all.cut[[15]][[1]] # calls the 15th column and the 1st row


h <- ggplot(data=rocplot,aes(x=fpr,y=tpr)) +
  geom_line() +
  geom_point(aes(x=0.195423,y=0.470588),size=3, color="red") # x = 1- specificity y = sensitivity 
h = h + theme_bw()
h = h + scale_x_continuous(expand = c(0,0)) #gets rid of space on x-axis
h = h + scale_y_continuous(expand = c(0,0)) #gets rid of space on x-axis
h = h + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())
h = h + theme(panel.border = element_blank(), axis.line=element_line())
h = h +  theme(axis.title=element_text(size=16))
h = h + theme(axis.text=element_text(size=12))
h = h + labs(title = "Log 10 concentration of HAEM vs. Observed blood in stool", y="Sensitivity", x="1-Specificity")+ theme(axis.text = element_text(colour = "black"))
h


# save figure
# ggsave(file="Figures/Figure 2 - HAEM and dysentery ROC curve - Shigella detected NOT attributed - FINAL.pdf",h,width=8,height=6,dpi=300)
