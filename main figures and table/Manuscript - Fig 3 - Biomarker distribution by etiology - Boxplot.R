
# boxplots characterizing the distribution of biomarker concentrations by diarrhea etiology

# Tidy up the workspace
rm(list=ls()) 

# Expand print options to avoid truncation
options(max.print = 1000000)

# to prevent scientific notation on the plot
options(scipen=999)

library(dplyr)
library(tidyr)
library(ggplot2)

# load in ELISA and TAC data 
load("R output/EFGH_substudy_rowperpid_ALLDATA_LOCKEDv3.Rda")


# keep necessary variables (pid, country, biomarkers, pathogen_attributable where pathogen_35 >5% of samples)
df <- EFGH_substudy_rowperpid_ALLDATA %>%
  select(pid, country, clean_conc_log10_CAL, clean_conc_log10_HAEM, clean_conc_log10_MPO, clean_conc_log10_NGAL, shigella_attributable,
         adenovirus_40_41_attributable, astrovirus_attributable, c_jejuni_coli_attributable, cryptosporidium_attributable, norovirus_gii_attributable, 
         rotavirus_attributable, sapovirus_attributable, ST.ETEC_attributable, tEPEC_attributable)

# Transform the dataframe so that there is only 1 column with the biomarker concentrations and 1 column indicating the name of the biomarker
long_df <- df %>%
  pivot_longer(
    cols = starts_with("clean_conc_log10_"),
    names_to = "biomarker",
    values_to = "log10_conc") %>% # Extract the biomarker names from the column names
  mutate(biomarker = gsub("clean_conc_log10_", "", biomarker))

# filter out where concentrations are NA
long_df <- long_df %>%
  filter(!is.na(log10_conc))

# create a function to create separate dataframe for each pathogen so that we can select where that pathogen is positive
# 
pathogen_function <- function(df, pathogen_attributable, pathogen_name) {
  result <- df %>%
    select(pid, !!sym(pathogen_attributable), biomarker, log10_conc) %>% # select the columns dynamically using sym() and !! within select()
    filter(!!sym(pathogen_attributable) == 1) %>% # filter rows where the pathogen column equals 1
    mutate(pathogen = pathogen_name) %>% # add a new column with the pathogen name
    rename(pathogen_present = !!sym(pathogen_attributable)) # rename the pathogen column to pathogen_present
  
  return(result)
}

# Example usage with long_df
adeno <- pathogen_function(long_df, "adenovirus_40_41_attributable", "Adenovirus")
astro <- pathogen_function(long_df, "astrovirus_attributable", "Astrovirus")
campy <- pathogen_function(long_df, "c_jejuni_coli_attributable", "Campylobacter")
crypto <- pathogen_function(long_df, "cryptosporidium_attributable", "Cryptosporidium")
noro <- pathogen_function(long_df, "norovirus_gii_attributable", "Norovirus")
rota <- pathogen_function(long_df, "rotavirus_attributable", "Rotavirus")
sapo <- pathogen_function(long_df, "sapovirus_attributable", "Sapovirus")
shigella <- pathogen_function(long_df, "shigella_attributable", "Shigella")
ST_ETEC <- pathogen_function(long_df, "ST.ETEC_attributable", "ST_ETEC")
tEPEC <- pathogen_function(long_df, "tEPEC_attributable", "t_EPEC")

# bind datasets together
bind <- rbind(adeno, astro, campy, crypto, noro, rota, sapo, shigella, ST_ETEC, tEPEC)

# set up pathogen levels
bind$pathogen <- factor(bind$pathogen, levels = c("Shigella","t_EPEC","Campylobacter","ST_ETEC","Cryptosporidium","Adenovirus","Astrovirus","Norovirus","Rotavirus","Sapovirus"),
                        labels = c("Shigella","tEPEC","Campylobacter jejuni/coli","ST-ETEC","Cryptosporidium","Adenovirus 40/41","Astrovirus","Norovirus GII","Rotavirus","Sapovirus"))

# set up biomarker levels
bind$biomarker <- factor(bind$biomarker, levels = c("HAEM", "NGAL", "MPO","CAL"), 
                         labels = c("Hemoglobin", "Lipocalin-2", "Myeloperoxidase","Calprotectin"))



# Figure 
plot = ggplot(bind, aes(x=biomarker, y=log10_conc, fill=pathogen)) 
plot = plot + theme_bw()
plot = plot + geom_boxplot() 
plot = plot + 
  labs(
    x = "",
    y = "Log10 (biomarker concentration)")
plot = plot + theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank())
plot = plot + scale_fill_manual(values=c("#AC1408", "#F42414", "#F8766D", "#FA9F98", "#FFC107","#005D82", "#0083B8", "#09B9FF", "#61D2FF", "#B3E9FF"), 
                                name="",
                                breaks=c("Shigella","tEPEC","Campylobacter jejuni/coli","ST-ETEC","Cryptosporidium","Adenovirus 40/41","Astrovirus","Norovirus GII","Rotavirus","Sapovirus"))
plot = plot + theme(legend.position = c(.85, .28))
plot = plot +  theme(axis.title=element_text(size=16))
plot = plot +  theme(axis.text=element_text(size=14))
plot = plot +  theme(axis.line = element_line(color="black"))
plot

# Save Figure 
ggsave(file="Figures/Figure 3 - Distribution of biomarkers by diarrhea etiology - FINAL.pdf",plot,width=12,height=6,dpi=300)


