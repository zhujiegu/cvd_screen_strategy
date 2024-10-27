in_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/'

library(tidyverse)
library(magrittr)
select <- dplyr::select

load(paste0(in_path,'outcomes_new.RData'))
load(paste0(in_path,'exposures_merged_new_popmean.RData'))

# align exposure and outcome patid
patid_unique <- unique(exposures_merged$patid)
stopifnot(all(patid_unique %in% outcomes$patid))
outcomes %<>% filter(patid %in% patid_unique)

exposures_merged  <- merge(exposures_merged, outcomes[,c("patid", "start_date")], by='patid')

# first measurement of time-dependent variables after entry
exposures_first <- exposures_merged %>% filter(exp_date > start_date) %>% 
  group_by(patid) %>% filter(!duplicated(exposure)) %>% ungroup
# slice_min(order_by = exp_age)

#####################
# overall
outcomes %>% summarise(N=n(),age_mean=mean((start_date - d_yob )/365.25), age_sd=sd((start_date - d_yob )/365.25))
((outcomes$end_date-outcomes$start_date)/365.25) %>%  as.numeric %>% summary()
outcomes %>% summarise(n=sum(cvd_ind, na.rm = T),mean=sum(cvd_ind, na.rm = T)/length(cvd_ind))

outcomes$gender %>% table
is.na(outcomes$Townsend) %>% mean


#####################
# Nr of measurements
exposures_count <- exposures_merged %>% select(patid, exposure) %>% group_by(patid, exposure) %>%
  summarise(count = n(), .groups = "drop")

exposures_count %<>% group_by(patid) %>% mutate(type_exposure = n()) %>% ungroup

# Histogram for the number of patid with each type_exposure value
p1 <- ggplot(exposures_count %>% distinct(patid, type_exposure), aes(x = type_exposure)) +
  geom_histogram(binwidth = 1, color = "black", fill = "skyblue") +
  labs(title = "Count of individuals by Nr. measured exposure types",
       x = "Count of exposure types",
       y = "Number of individuals") +
  theme_minimal()

# Histogram for the number of patid with each count for each exposure
p2 <- ggplot(exposures_count, aes(x = count)) +
  geom_histogram(binwidth = 1, color = "black", fill = "lightgreen") +
  facet_wrap(~ factor(exposure, levels = c("bmi", "sbp", "tchol", "hdl", "smokbin")), 
             scales = "free", 
             labeller = as_labeller(c(bmi = "BMI", sbp = "SBP", 
                                      tchol = "Tot. Chol.", hdl = "HDL", 
                                      smokbin = "Smoke"))) +
  xlim(c(0,20)) +
  labs(title = "Count of individuals by number of each exposure (truncaated at 20)",
       x = "Count of measurements",
       y = "Number of individuals") +
  theme_minimal()

jpeg(file = '~/epi_paper/outp_figs/Nr_type.jpeg', units="in", width=5, height=2.5, res=600)
p1
dev.off()

jpeg(file = '~/epi_paper/outp_figs/Nr_each.jpeg', units="in", width=8, height=5, res=600)
p2
dev.off()
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
# first smk after entry
exposures_first %>% filter(gender=='Male')%>% 
  filter(exposure=='smokbin') %>% summarise(N=n(),total=sum(original),mean=mean(original), sd=sd(original))

#BP
outcomes %>% filter(gender=='Male') %>% summarise(total=sum(bp_med_user, na.rm = T),mean=total/n())
# Diabetes
outcomes %>% filter(gender=='Male') %>% summarise(total=sum(diab_ind, na.rm = T),mean=total/n())
# renal 
outcomes %>% filter(gender=='Male') %>% summarise(total=sum(renal_ind, na.rm = T),mean=total/n())
#migraine
outcomes %>% filter(gender=='Male') %>% summarise(total=sum(Migraine_ind, na.rm = T),mean=total/n())
#rheumatoid arthritis
outcomes %>% filter(gender=='Male') %>% summarise(total=sum(rheumatoid_arthritis_ind, na.rm = T),mean=total/n())
#depression
outcomes %>% filter(gender=='Male') %>% summarise(total=sum(depression_ind, na.rm = T),mean=total/n())
#severe mental illness
outcomes %>% filter(gender=='Male') %>% summarise(total=sum(Severe_mental_illness_ind, na.rm = T),mean=total/n())
#atrial fibrillation
outcomes %>% filter(gender=='Male') %>% summarise(total=sum(Atrial_fibrillation_ind, na.rm = T),mean=total/n())
#Townsend
nrow(outcomes %>% filter(gender=='Male', !is.na(Townsend)))
(outcomes %>% filter(gender=='Male') %>% select(Townsend) %>% table)/nrow(outcomes %>% filter(gender=='Male', !is.na(Townsend)))


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
# first smk after entry
exposures_first %>% filter(gender=='Female')%>% 
  filter(exposure=='smokbin') %>% summarise(N=n(),total=sum(original),mean=mean(original), sd=sd(original))

#BP
outcomes %>% filter(gender=='Female') %>% summarise(total=sum(bp_med_user, na.rm = T),mean=total/n())
# Diabetes
outcomes %>% filter(gender=='Female') %>% summarise(total=sum(diab_ind, na.rm = T),mean=total/n())
# renal 
outcomes %>% filter(gender=='Female') %>% summarise(total=sum(renal_ind, na.rm = T),mean=total/n())
#migraine
outcomes %>% filter(gender=='Female') %>% summarise(total=sum(Migraine_ind, na.rm = T),mean=total/n())
#rheumatoid arthritis
outcomes %>% filter(gender=='Female') %>% summarise(total=sum(rheumatoid_arthritis_ind, na.rm = T),mean=total/n())
#depression
outcomes %>% filter(gender=='Female') %>% summarise(total=sum(depression_ind, na.rm = T),mean=total/n())
#severe mental illness
outcomes %>% filter(gender=='Female') %>% summarise(total=sum(Severe_mental_illness_ind, na.rm = T),mean=total/n())
#atrial fibrillation
outcomes %>% filter(gender=='Female') %>% summarise(total=sum(Atrial_fibrillation_ind, na.rm = T),mean=total/n())
#Townsend
nrow(outcomes %>% filter(gender=='Female', !is.na(Townsend)))
(outcomes %>% filter(gender=='Female') %>% select(Townsend) %>% table)/nrow(outcomes %>% filter(gender=='Female', !is.na(Townsend)))

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


################################################
# in_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/data_created/'
# j=60
# gender='female'
# load(paste0(in_path, gender, "/data_lm_",j, "_", gender,".RData"))
# max(data_ext$end_age - j)
# which.max(data_ext$end_age - j)
