in_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/'

library(tidyverse)
library(magrittr)

load(paste0(in_path,'outcomes.RData'))
load(paste0(in_path,'exposures_red_merged.RData'))

exposures_red_merged  <- merge(exposures_red_merged, outcomes[,c("patid", "start_date")], by='patid')

# first measurement of time-dependent variables after entry
exposures_first <- exposures_red_merged %>% filter(exp_date > start_date) %>% 
  group_by(patid) %>% filter(!duplicated(exposure)) %>% ungroup
# slice_min(order_by = exp_age)

#####################
# overall
outcomes %>% 
  summarise(N=n(),age_mean=mean((start_date - d_yob )/365.25), age_sd=sd((start_date - d_yob )/365.25))
((outcomes$end_date-outcomes$start_date)/365.25) %>%  as.numeric %>% summary()
outcomes %>% summarise(n=sum(cvd_ind, na.rm = T),mean=sum(cvd_ind, na.rm = T)/length(cvd_ind))

#####################
# Male

# age
outcomes %>% filter(gender=='Male')%>% 
  summarise(N=n(),age_mean=mean((start_date - d_yob )/365.25), age_sd=sd((start_date - d_yob )/365.25))
# CVD
outcomes %>% filter(gender=='Male')%>%   
  summarise(n=sum(cvd_ind, na.rm = T),mean=sum(cvd_ind, na.rm = T)/length(cvd_ind))
# first BMI after entry
exposures_first %>% filter(gender=='Male')%>% 
  filter(exposure=='bmi') %>% summarise(N=n(),mean=mean(original), sd=sd(original))
# first SBP after entry
exposures_first %>% filter(gender=='Male') %>% 
  filter(exposure=='sbp') %>% summarise(N=n(),mean=mean(original), sd=sd(original))
# first tchol after entry
exposures_first %>% filter(gender=='Male')  %>% 
  filter(exposure=='tchol') %>% summarise(N=n(),mean=mean(original), sd=sd(original))
# first HDL after entry
exposures_first %>% filter(gender=='Male')  %>% 
  filter(exposure=='hdl') %>% summarise(N=n(),mean=mean(original), sd=sd(original))
# first HDL after entry
exposures_first %>% filter(gender=='Male')%>% 
  filter(exposure=='smokbin') %>% summarise(N=n(),total=sum(original),mean=mean(original), sd=sd(original))
# Diabetes
outcomes %>% filter(gender=='Male') %>% summarise(N=sum(!is.na(diab_ind)), total=sum(diab_ind),mean=mean(diab_ind))
# renal 
outcomes %>% filter(gender=='Male') %>% summarise(N=sum(!is.na(renal_ind)), total=sum(renal_ind, na.rm = T),mean=mean(renal_ind, na.rm = T))
#migraine
outcomes %>% filter(gender=='Male') %>% summarise(N=sum(!is.na(Migraine_ind)), total=sum(Migraine_ind, na.rm = T),mean=mean(Migraine_ind, na.rm = T))
#rheumatoid arthritis
outcomes %>% filter(gender=='Male') %>% summarise(N=sum(!is.na(rheumatoid_arthritis_ind)), total=sum(rheumatoid_arthritis_ind, na.rm = T),mean=mean(rheumatoid_arthritis_ind, na.rm = T))
#depression
outcomes %>% filter(gender=='Male') %>% summarise(N=sum(!is.na(depression_ind)), total=sum(depression_ind, na.rm = T),mean=mean(depression_ind, na.rm = T))
#severe mental illness
outcomes %>% filter(gender=='Male') %>% summarise(N=sum(!is.na(Severe_mental_illness_ind)), total=sum(Severe_mental_illness_ind, na.rm = T),mean=mean(Severe_mental_illness_ind, na.rm = T))
#atrial fibrillation
outcomes %>% filter(gender=='Male') %>% summarise(N=sum(!is.na(Atrial_fibrillation_ind)), total=sum(Atrial_fibrillation_ind, na.rm = T),mean=mean(Atrial_fibrillation_ind, na.rm = T))
#Townsend 20 
outcomes %>% filter(gender=='Male') %>% summarise(N=sum(!is.na(townsend_20)), mean=mean(townsend_20, na.rm = T),sd=sd(townsend_20, na.rm = T))
#BP
outcomes %>% filter(gender=='Male') %>% summarise(N=sum(!is.na(townsend_20)), mean=mean(townsend_20, na.rm = T),sd=sd(townsend_20, na.rm = T))


# Female

# age
outcomes %>% filter(gender=='Female')%>% 
  summarise(N=n(),age_mean=mean((start_date - d_yob )/365.25), age_sd=sd((start_date - d_yob )/365.25))
# CVD
outcomes %>% filter(gender=='Female')%>%   
  summarise(n=sum(cvd_ind, na.rm = T),mean=sum(cvd_ind, na.rm = T)/length(cvd_ind))
# first BMI after entry
exposures_first %>% filter(gender=='Female')%>% 
  filter(exposure=='bmi') %>% summarise(N=n(),mean=mean(original), sd=sd(original))
# first SBP after entry
exposures_first %>% filter(gender=='Female') %>% 
  filter(exposure=='sbp') %>% summarise(N=n(),mean=mean(original), sd=sd(original))
# first tchol after entry
exposures_first %>% filter(gender=='Female')  %>% 
  filter(exposure=='tchol') %>% summarise(N=n(),mean=mean(original), sd=sd(original))
# first HDL after entry
exposures_first %>% filter(gender=='Female')  %>% 
  filter(exposure=='hdl') %>% summarise(N=n(),mean=mean(original), sd=sd(original))
# first HDL after entry
exposures_first %>% filter(gender=='Female')%>% 
  filter(exposure=='smokbin') %>% summarise(N=n(),total=sum(original),mean=mean(original), sd=sd(original))
# Diabetes
outcomes %>% filter(gender=='Female') %>% summarise(N=sum(!is.na(diab_ind)), total=sum(diab_ind),mean=mean(diab_ind))
# renal 
outcomes %>% filter(gender=='Female') %>% summarise(N=sum(!is.na(renal_ind)), total=sum(renal_ind, na.rm = T),mean=mean(renal_ind, na.rm = T))
#migraine
outcomes %>% filter(gender=='Female') %>% summarise(N=sum(!is.na(Migraine_ind)), total=sum(Migraine_ind, na.rm = T),mean=mean(Migraine_ind, na.rm = T))
#rheumatoid arthritis
outcomes %>% filter(gender=='Female') %>% summarise(N=sum(!is.na(rheumatoid_arthritis_ind)), total=sum(rheumatoid_arthritis_ind, na.rm = T),mean=mean(rheumatoid_arthritis_ind, na.rm = T))
#depression
outcomes %>% filter(gender=='Female') %>% summarise(N=sum(!is.na(depression_ind)), total=sum(depression_ind, na.rm = T),mean=mean(depression_ind, na.rm = T))
#severe mental illness
outcomes %>% filter(gender=='Female') %>% summarise(N=sum(!is.na(Severe_mental_illness_ind)), total=sum(Severe_mental_illness_ind, na.rm = T),mean=mean(Severe_mental_illness_ind, na.rm = T))
#atrial fibrillation
outcomes %>% filter(gender=='Female') %>% summarise(N=sum(!is.na(Atrial_fibrillation_ind)), total=sum(Atrial_fibrillation_ind, na.rm = T),mean=mean(Atrial_fibrillation_ind, na.rm = T))
#Townsend 20 
outcomes %>% filter(gender=='Female') %>% summarise(N=sum(!is.na(townsend_20)), mean=mean(townsend_20, na.rm = T),sd=sd(townsend_20, na.rm = T))

#####################
# Summary with derivation and validation
#####################
# Male derivation

# age
outcomes %>% filter(gender=='Male') %>% filter(derivation == 'derivation') %>% 
  summarise(N=n(),age_mean=mean((start_date - d_yob )/365.25), age_sd=sd((start_date - d_yob )/365.25))
# CVD
outcomes %>% filter(gender=='Male') %>% filter(derivation == 'derivation') %>%   
  summarise(n=sum(cvd_ind, na.rm = T),mean=sum(cvd_ind, na.rm = T)/length(cvd_ind))
# first BMI after entry
exposures_first %>% filter(gender=='Male') %>% filter(derivation == 'derivation') %>% 
  filter(exposure=='bmi') %>% summarise(N=n(),mean=mean(original), sd=sd(original))
# first SBP after entry
exposures_first %>% filter(gender=='Male') %>% filter(derivation == 'derivation') %>% 
  filter(exposure=='sbp') %>% summarise(N=n(),mean=mean(original), sd=sd(original))
# first tchol after entry
exposures_first %>% filter(gender=='Male') %>% filter(derivation == 'derivation') %>% 
  filter(exposure=='tchol') %>% summarise(N=n(),mean=mean(original), sd=sd(original))
# first HDL after entry
exposures_first %>% filter(gender=='Male') %>% filter(derivation == 'derivation') %>% 
  filter(exposure=='hdl') %>% summarise(N=n(),mean=mean(original), sd=sd(original))
# first HDL after entry
exposures_first %>% filter(gender=='Male') %>% filter(derivation == 'derivation') %>% 
  filter(exposure=='smokbin') %>% summarise(N=n(),mean=mean(original), sd=sd(original))






# Male validation

# age
outcomes %>% filter(gender=='Male') %>% filter(derivation == 'validation') %>% 
  summarise(N=n(),age_mean=mean((start_date - d_yob )/365.25), age_sd=sd((start_date - d_yob )/365.25))
# CVD
outcomes %>% filter(gender=='Male') %>% filter(derivation == 'validation') %>%   
  summarise(n=sum(cvd_ind, na.rm = T),mean=sum(cvd_ind, na.rm = T)/length(cvd_ind))
# first BMI after entry
exposures_first %>% filter(gender=='Male') %>% filter(derivation == 'validation') %>% 
  filter(exposure=='bmi') %>% summarise(N=n(),mean=mean(original), sd=sd(original))
# first SBP after entry
exposures_first %>% filter(gender=='Male') %>% filter(derivation == 'validation') %>% 
  filter(exposure=='sbp') %>% summarise(N=n(),mean=mean(original), sd=sd(original))
# first tchol after entry
exposures_first %>% filter(gender=='Male') %>% filter(derivation == 'validation') %>% 
  filter(exposure=='tchol') %>% summarise(N=n(),mean=mean(original), sd=sd(original))
# first HDL after entry
exposures_first %>% filter(gender=='Male') %>% filter(derivation == 'validation') %>% 
  filter(exposure=='hdl') %>% summarise(N=n(),mean=mean(original), sd=sd(original))
# first HDL after entry
exposures_first %>% filter(gender=='Male') %>% filter(derivation == 'validation') %>% 
  filter(exposure=='smokbin') %>% summarise(N=n(),mean=mean(original), sd=sd(original))


#####################
# Female derivation

# age
outcomes %>% filter(gender=='Female') %>% filter(derivation == 'derivation') %>% 
  summarise(N=n(),age_mean=mean((start_date - d_yob )/365.25), age_sd=sd((start_date - d_yob )/365.25))
# CVD
outcomes %>% filter(gender=='Female') %>% filter(derivation == 'derivation') %>%   
  summarise(n=sum(cvd_ind, na.rm = T),mean=sum(cvd_ind, na.rm = T)/length(cvd_ind))
# first BMI after entry
exposures_first %>% filter(gender=='Female') %>% filter(derivation == 'derivation') %>% 
  filter(exposure=='bmi') %>% summarise(N=n(),mean=mean(original), sd=sd(original))
# first SBP after entry
exposures_first %>% filter(gender=='Female') %>% filter(derivation == 'derivation') %>% 
  filter(exposure=='sbp') %>% summarise(N=n(),mean=mean(original), sd=sd(original))
# first tchol after entry
exposures_first %>% filter(gender=='Female') %>% filter(derivation == 'derivation') %>% 
  filter(exposure=='tchol') %>% summarise(N=n(),mean=mean(original), sd=sd(original))
# first HDL after entry
exposures_first %>% filter(gender=='Female') %>% filter(derivation == 'derivation') %>% 
  filter(exposure=='hdl') %>% summarise(N=n(),mean=mean(original), sd=sd(original))
# first HDL after entry
exposures_first %>% filter(gender=='Female') %>% filter(derivation == 'derivation') %>% 
  filter(exposure=='smokbin') %>% summarise(N=n(),mean=mean(original), sd=sd(original))

# Female validation

# age
outcomes %>% filter(gender=='Female') %>% filter(derivation == 'validation') %>% 
  summarise(N=n(),age_mean=mean((start_date - d_yob )/365.25), age_sd=sd((start_date - d_yob )/365.25))
# CVD
outcomes %>% filter(gender=='Female') %>% filter(derivation == 'validation') %>%   
  summarise(n=sum(cvd_ind, na.rm = T),mean=sum(cvd_ind, na.rm = T)/length(cvd_ind))
# first BMI after entry
exposures_first %>% filter(gender=='Female') %>% filter(derivation == 'validation') %>% 
  filter(exposure=='bmi') %>% summarise(N=n(),mean=mean(original), sd=sd(original))
# first SBP after entry
exposures_first %>% filter(gender=='Female') %>% filter(derivation == 'validation') %>% 
  filter(exposure=='sbp') %>% summarise(N=n(),mean=mean(original), sd=sd(original))
# first tchol after entry
exposures_first %>% filter(gender=='Female') %>% filter(derivation == 'validation') %>% 
  filter(exposure=='tchol') %>% summarise(N=n(),mean=mean(original), sd=sd(original))
# first HDL after entry
exposures_first %>% filter(gender=='Female') %>% filter(derivation == 'validation') %>% 
  filter(exposure=='hdl') %>% summarise(N=n(),mean=mean(original), sd=sd(original))
# first HDL after entry
exposures_first %>% filter(gender=='Female') %>% filter(derivation == 'validation') %>% 
  filter(exposure=='smokbin') %>% summarise(N=n(),mean=mean(original), sd=sd(original))