
# risk factors forest plot - where co-attributions have been removed

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# to prevent scientific notation on the plot
options(scipen=999)

library(dplyr)
library(tidyr)
library(openxlsx)
library(gee)
library(geepack)
library(ggplot2)
library(forcats)

# load in elisa data
load("R output/EFGH_substudy_rowperpid_ALLDATA_LOCKEDv3.Rda")

# check NAs
table(EFGH_substudy_rowperpid_ALLDATA$clean_conc_log10_HAEM, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$clean_conc_log10_CAL, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$clean_conc_log10_MPO, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$clean_conc_log10_NGAL, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$dysentery, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$enr_age_months, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$enroll_diar_vom, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$enroll_diar_fever, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$enr_haz, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$enr_muac_before_cm, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$enr_waz, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$enr_whz, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$country, useNA = "always")
table(EFGH_substudy_rowperpid_ALLDATA$sex, useNA = "always")

#filter out NAs
small <-EFGH_substudy_rowperpid_ALLDATA %>%
  filter(!is.na(clean_conc_log10_HAEM)) %>%
  filter(!is.na(clean_conc_log10_CAL)) %>%
  filter(!is.na(clean_conc_log10_NGAL)) %>%
  filter(!is.na(clean_conc_log10_MPO)) %>%
  filter(!is.na(enr_age_months)) %>%
  filter(!is.na(sex)) %>%
  filter(!is.na(enroll_diar_fever)) %>%
  filter(!is.na(enroll_diar_vom)) %>%
  filter(!is.na(dysentery)) %>%
  filter(!is.na(country)) %>%
  filter(!is.na(enr_muac_before_cm)) %>%
  filter(!is.na(enr_haz)) %>%
  filter(!is.na(enr_waz)) %>%
  filter(!is.na(enr_whz))

# filter to where shigella is attributable
small <- small %>% filter(
  shigella_attributable == 1,
  tEPEC_attributable == 0,
  c_jejuni_coli_attributable == 0,
  ST.ETEC_attributable == 0,
  cryptosporidium_attributable == 0,
  adenovirus_40_41_attributable == 0,
  astrovirus_attributable == 0,
  norovirus_gii_attributable == 0,
  rotavirus_attributable == 0,
  sapovirus_attributable == 0)

table(small$tEPEC_attributable)
table(small$c_jejuni_coli_attributable)
table(small$ST.ETEC_attributable)
table(small$cryptosporidium_attributable)
table(small$adenovirus_40_41_attributable)
table(small$astrovirus_attributable)
table(small$norovirus_gii_attributable)
table(small$rotavirus_attributable)
table(small$sapovirus_attributable)


# one model per biomarker with mutual adjustment
risk_func <-function(biomarker, dat) {
cal_model <- geeglm(biomarker ~ sex + enr_age_months + enroll_diar_fever + enroll_diar_vom + dysentery + enr_muac_before_cm + enr_haz + enr_whz +
                    country, family = gaussian, data = dat, id = child_id, corstr = "exchangeable")

beta <- summary(cal_model)$coefficients[2:9,1]
lb <- summary(cal_model)$coefficients[2:9,1] - 1.96 * summary(cal_model)$coefficients[2:9,2]
ub <- summary(cal_model)$coefficients[2:9,1] + 1.96 * summary(cal_model)$coefficients[2:9,2]

# Save RR and CIs as a dataframe 
beta <- as.data.frame(round(beta,2))
lb <- as.data.frame(round(lb,2))
ub <- as.data.frame(round(ub,2))

# Bind RR and CIs  
return(abxRRplot <- as.data.frame(cbind(beta, lb, ub)))

}

HAEM_model <- risk_func(small$clean_conc_log10_HAEM, dat=small)
NGAL_model <- risk_func(small$clean_conc_log10_NGAL, dat=small)
MPO_model <- risk_func(small$clean_conc_log10_MPO, dat=small)
CAL_model <- risk_func(small$clean_conc_log10_CAL, dat=small)


biomarker_adj <- as.data.frame(rbind(HAEM_model, NGAL_model, MPO_model, CAL_model))

# Formatting the dataframe
colnames(biomarker_adj) <- c("MeanDifference", "LowerCI", "UpperCI")

# create biomarker column
biomarker_adj$biomarker <- c("HAEM","HAEM","HAEM","HAEM","HAEM","HAEM","HAEM","HAEM",
                             "NGAL","NGAL","NGAL","NGAL","NGAL","NGAL","NGAL","NGAL",
                             "MPO","MPO","MPO","MPO","MPO","MPO","MPO","MPO",
                             "CAL","CAL","CAL","CAL","CAL","CAL","CAL","CAL")

# create risk factor column
biomarker_adj$risk_factor <- c("Female", "Age (in months)", "Fever", "Vomiting", "Dysentery", "MUAC", "LAZ/HAZ", "WHZ")

# This will keep the order consistent with all of the other figures 
biomarker_adj$risk_factor <- factor(biomarker_adj$risk_factor, levels = c("Female", "Age (in months)", "Fever", "Vomiting", "Dysentery", "MUAC", "LAZ/HAZ", "WHZ"),
                          labels = c("Female", "Age (in months)", "Fever", "Vomiting", "Dysentery", "MUAC", "LAZ/HAZ",  "WHZ"))

biomarker_adj$biomarker <- factor(biomarker_adj$biomarker, levels = c("HAEM", "NGAL", "MPO", "CAL"),
                                    labels = c("Hemoglobin", "Lipocalin-2", "Myeloperoxidase","Calprotectin"))

#fct_review puts abx above no abx in the plot
q =  ggplot(biomarker_adj, aes(x = fct_rev(biomarker),y = MeanDifference, ymin = LowerCI, ymax = UpperCI, shape=biomarker))
q = q + geom_pointrange(aes(col=risk_factor))
#q = q + scale_colour_manual("dataset",values=c("#9590FF", "#00BF7D", "#00B0F6"))
q = q + theme_bw()
q = q + theme(strip.background = element_blank())
q = q + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())
q = q + theme(panel.border = element_blank(), axis.line.x = element_line(),
              axis.line.y = element_blank())
q = q +   xlab("")
q = q +   ylab("Mean Differences")
q = q +   geom_errorbar(aes(ymin=LowerCI, ymax=UpperCI,col=risk_factor),width=1,cex=0.5)
q = q +   facet_wrap(~risk_factor,strip.position="left",nrow=9) 
q = q +   theme(strip.text = element_text(size = 12))
q = q +   theme(axis.text.y=element_blank(),
                axis.ticks.y=element_blank())
q = q + theme(strip.text.y.left = element_text(angle = 0))
q = q + theme(axis.text.x = element_text(angle = 0,hjust=0.5,vjust=1,size=16, color="black"))
q = q + theme(axis.title=element_text(size=16))
q = q + theme(legend.text=element_text(size=16))
q = q + scale_y_continuous(limits=c(-1,2))
q = q + scale_color_manual(values=c("#AC1408", "#F42414", "#F8766D", "#FA9F98", "#FFC107","#005D82", "#0083B8", "#09B9FF", "#61D2FF"), 
                                name="",
                                breaks=c("Female", "Age (in months)", "Fever", "Vomiting", "Dysentery", "MUAC", "LAZ/HAZ", "WHZ"))
#q = q + theme(axis.ticks = element_blank())
q = q + guides(col = FALSE)
q = q + theme(
  legend.position = c(.85, .10),
  #  legend.justification = c("right", "bottom"),
  legend.title = element_blank())
q = q + geom_hline(yintercept =0, linetype=2)
q = q +   coord_flip()
q

# Save Figure 
ggsave(file="Figures/Figure SX - Risk Factors - no co-attribution FINAL.pdf",q,width=10,height=12,dpi=300)



