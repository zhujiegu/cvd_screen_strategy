library(haven)
library(dplyr)
library(magrittr)

out_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/'
##########################################
# Load Townsend20 covariate from the old data and merge with new Townsend5.
##########################################
# old
load('/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/patients_outcomes.RData')
townsend20 <- outcomes %>% select(patid, townsend_20)
# new
townsend <- read_dta('/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/dataprep/hes_2020/patient_townsend_new.dta')
townsend %<>% select(patid, townsend2001_5)
colnames(townsend)[2] <- 'Townsend'
# merge
townsend <- merge(townsend, townsend20, all = T)
townsend <- townsend %>% mutate(Townsend = ifelse(is.na(Townsend), ceiling(townsend_20/4), Townsend))
townsend %<>% select(patid, Townsend)
rm(outcomes, townsend20)
##########################################
# Load ethnicity
##########################################
ethnicity <- read_dta('/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/dataprep/cprd_2020/ethnic_cprd_new.dta')
ethnicity %<>% select(patid, ethnic_cprd)
colnames(ethnicity)[2] <- 'ethnicity'
# remove duplicated
ethnicity %<>% as.data.table()
ethnicity %<>% filter(ethnicity!=21)
dup <- duplicated(ethnicity$patid)
ethnicity <- ethnicity[!dup,]
##########################################
# Load data
##########################################
outcomes <- read_dta('/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/dataprep/cprd_2020/patientslist_outcomes_2020.dta')
outcomes %>% colnames
##########################################
# add death dates information
##########################################
death_date <- read_dta('/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/dataprep/hes_2020/death_dates_new.dta')
colnames(death_date) <- c('patid', 'death_date')
outcomes <- merge(outcomes,  death_date, by = 'patid', all.x = T)
# check
stopifnot(sum(!is.na(outcomes$death_date)) == sum(outcomes$death_ons_ind))
stopifnot(length(which( !is.na(outcomes$death_date) & outcomes$death_ons_ind!=1))==0)
##########################################
# Combine sources of information
##########################################
date_origin = '1960-01-01'

outcomes %<>% mutate(cvd_ind = ifelse(cvd_CPRD_ind==1|cvd_HES_primary_ind==1|cvd_HES_ind==1|cvd_ONS_primary_ind==1|cvd_ONS_ind==1,1,NA),
                     cvd_date = pmin(cvd_CPRD_date,cvd_HES_primary_date,cvd_HES_date,as.Date(cvd_ONS_primary_date,origin = date_origin),
                                     as.Date(cvd_ONS_date,origin = date_origin), na.rm=T),
                     death_ind = death_ons_ind,
                     diab_age_cprd = ifelse(diab_type_cprd %in% c(1,2), diab_age_cprd, NA),
                     diab_ind = ifelse(diab_type_cprd %in% c(1,2)|!is.na(diab_type2_date_hes)|!is.na(diab_type1_date_hes),1,NA),
                     diab_date = pmin(diab_type1_date_hes, diab_type2_date_hes, diab_age_cprd*365.25 + d_yob, na.rm=T), 
                     rheumatoid_arthritis_ind = rheumatoid_arthritis_CPRD_ind,
                     rheumatoid_arthritis_date = as.Date(rheumatoid_arthritis_CPRD_date,origin = date_origin),
                     Atrial_fibrillation_ind = Atrial_fibrillation_CPRD_ind,
                     Atrial_fibrillation_date = as.Date(Atrial_fibrillation_CPRD_date,origin = date_origin),
                     renal_ind = renal_CPRD_ind,
                     renal_date = as.Date(renal_CPRD_date,origin = date_origin),
                     depression_ind = depression_CPRD_ind,
                     depression_date = as.Date(depression_CPRD_date,origin = date_origin),
                     Migraine_ind = Migraine_CPRD_ind,
                     Migraine_date = as.Date(Migraine_CPRD_date,origin = date_origin),
                     Severe_mental_illness_ind = Severe_mental_illness_CPRD_ind,
                     Severe_mental_illness_date = as.Date(Severe_mental_illness_CPRD_date,origin = date_origin))
outcomes <- merge(outcomes, townsend, all.x = T)
outcomes <- merge(outcomes, ethnicity, all.x = T)
##########################################
# correct CVD and death post censoring
##########################################
#correct CVD post censoring
cvd_post_censoring = which( !is.na(outcomes$cvd_date) & outcomes$end_date < outcomes$cvd_date )
cvd_during_censoring = which( !is.na(outcomes$cvd_date) & outcomes$end_date == outcomes$cvd_date & is.na( outcomes$death_date )   )

length( which( !is.na(outcomes$cvd_date) & outcomes$end_date == outcomes$cvd_date))

length( cvd_post_censoring )
outcomes$cvd_date[ cvd_post_censoring ] = NA
outcomes$cvd_ind[ cvd_post_censoring ] = 0

outcomes$end_date[ cvd_during_censoring ] = outcomes$end_date[ cvd_during_censoring ] + 1
outcomes$cvd_ind[ cvd_during_censoring ] = 1 #it was already 1

cvd_post_censoring = which( !is.na(outcomes$cvd_date) & outcomes$end_date < outcomes$cvd_date )
stopifnot( length( cvd_post_censoring ) == 0 )
rm( cvd_post_censoring )

#uniform the notation for indices of events:
length( which(is.na( outcomes$cvd_ind)))
length( which(is.na( outcomes$cvd_ind ) & is.na(outcomes$cvd_date)))

outcomes$cvd_ind[ which(is.na( outcomes$cvd_ind)) ] = 0


#correct Death post censoring
death_post_censoring = which( !is.na(outcomes$death_date) & outcomes$end_date < outcomes$death_date )
length( death_post_censoring )
outcomes$death_date[ death_post_censoring ] = NA
outcomes$death_ind[ death_post_censoring ] = 0

death_post_censoring = which( !is.na(outcomes$death_date) & outcomes$end_date < outcomes$death_date )
stopifnot( length( death_post_censoring ) == 0 )
rm( death_post_censoring )

#correct CVD post Death
cvd_post_death = which( !is.na(outcomes$cvd_date) & outcomes$death_date < outcomes$cvd_date )
stopifnot( length( cvd_post_death ) == 0 )
rm( cvd_post_death )

##########################################
# CVD/Statin before start date
##########################################
cvd_before_start <- which( !is.na(outcomes$cvd_date) & outcomes$cvd_date <= outcomes$start_date )
statin_before_start <- which(!is.na(outcomes$statins_prscd) & outcomes$statins_prscd <= outcomes$start_date)
death_before_start <- which( !is.na(outcomes$death_date) & outcomes$death_date <= outcomes$start_date )
length(unique(c(cvd_before_start, statin_before_start))) # 371067

outcomes <- outcomes[-c(cvd_before_start, statin_before_start), ]

stopifnot( length(which( !is.na(outcomes$cvd_date) & outcomes$cvd_date < outcomes$start_date )) == 0 )
stopifnot( length(which( !is.na(outcomes$statins_prscd) & outcomes$statins_prscd < outcomes$start_date )) == 0 )
stopifnot( length(which( !is.na(outcomes$death_date) & outcomes$death_date < outcomes$start_date )) == 0 )

##########################################
# Adjust covariates diagnosis after end date (for summary statistics)
##########################################
outcomes %<>% select(patid, gender, d_yob, start_date, end_date, derivation, cvd_ind, cvd_date, death_ind, death_date,
                     diab_ind, diab_date, rheumatoid_arthritis_ind, rheumatoid_arthritis_date,
                     Atrial_fibrillation_ind, Atrial_fibrillation_date, renal_ind, renal_date,
                     depression_ind, depression_date, Migraine_ind, Migraine_date,
                     Severe_mental_illness_ind, Severe_mental_illness_date, statins_prscd, statinuser,
                     bp_med_prscd, bp_med_user, Townsend, ethnicity)

for(var in c('diab', 'rheumatoid_arthritis', 'Atrial_fibrillation', 'renal', 'depression', 'Migraine', 'Severe_mental_illness')){
  var_post_exit = which(!is.na(outcomes[,paste0(var,'_date')]) & outcomes[,paste0(var,'_date')] > outcomes$end_date)
  print(paste(var, length( var_post_exit )))
  outcomes[,paste0(var,'_date')][ var_post_exit ] = NA
  outcomes[,paste0(var,'_ind')][ var_post_exit ] = 0
}

var_post_exit = which(!is.na(outcomes$statins_prscd) & outcomes$statins_prscd > outcomes$end_date)
print(paste('statin', length( var_post_exit )))
outcomes$statins_prscd[ var_post_exit ] = NA
outcomes$statinuser[ var_post_exit ] = 0

var_post_exit = which(!is.na(outcomes$bp_med_prscd) & outcomes$bp_med_prscd > outcomes$end_date)
print(paste('bp_med', length( var_post_exit )))
outcomes$bp_med_prscd[ var_post_exit ] = NA
outcomes$bp_med_user[ var_post_exit ] = 0

# outcomes <- outcomes %>% mutate(diab_ind = ifelse(diab_date>end_date, NA, diab_ind),
#                                 diab_date = ifelse(is.na(diab_ind), NA, diab_date),
#                                 rheumatoid_arthritis_ind = ifelse(rheumatoid_arthritis_date>end_date, NA, rheumatoid_arthritis_ind),
#                                 rheumatoid_arthritis_date = ifelse(is.na(rheumatoid_arthritis_ind), NA, rheumatoid_arthritis_date),
#                                 Atrial_fibrillation_ind = ifelse(Atrial_fibrillation_date>end_date, NA, Atrial_fibrillation_ind),
#                                 Atrial_fibrillation_date = ifelse(is.na(Atrial_fibrillation_ind), NA, Atrial_fibrillation_date),
#                                 renal_ind = ifelse(renal_date>end_date, NA, renal_ind),
#                                 renal_date = ifelse(is.na(renal_ind), NA, renal_date),
#                                 depression_ind = ifelse(depression_date>end_date, NA, depression_ind),
#                                 depression_date = ifelse(is.na(depression_ind), NA, depression_date),
#                                 Migraine_ind = ifelse(Migraine_date>end_date, NA, Migraine_ind),
#                                 Migraine_date = ifelse(is.na(Migraine_ind), NA, Migraine_date),
#                                 Severe_mental_illness_ind = ifelse(Severe_mental_illness_date>end_date, NA, Severe_mental_illness_ind),
#                                 Severe_mental_illness_date = ifelse(is.na(Severe_mental_illness_ind), NA, Severe_mental_illness_date))
# check
stopifnot( length( which( outcomes$cvd_ind == 1 ) ) ==  length( which( !is.na( outcomes$cvd_date ) ) ) )
stopifnot( length( which( outcomes$death_ind == 1 ) ) ==  length( which( !is.na( outcomes$death_date ) ) ) )
stopifnot( length( which( outcomes$diab_ind == 1 ) ) ==  length( which( !is.na( outcomes$diab_date ) ) ) )
stopifnot( length( which( outcomes$rheumatoid_arthritis_ind == 1 ) ) ==  length( which( !is.na( outcomes$rheumatoid_arthritis_date ) ) ) )
stopifnot( length( which( outcomes$Atrial_fibrillation_ind == 1 ) ) ==  length( which( !is.na( outcomes$Atrial_fibrillation_date ) ) ) )
stopifnot( length( which( outcomes$renal_ind == 1 ) ) ==  length( which( !is.na( outcomes$renal_date ) ) ) )
stopifnot( length( which( outcomes$depression_ind == 1 ) ) ==  length( which( !is.na( outcomes$depression_date ) ) ) )
stopifnot( length( which( outcomes$Migraine_ind == 1 ) ) ==  length( which( !is.na( outcomes$Migraine_date ) ) ) )
stopifnot( length( which( outcomes$Severe_mental_illness_ind == 1 ) ) ==  length( which( !is.na( outcomes$Severe_mental_illness_date ) ) ) )
stopifnot( length( which( outcomes$statinuser == 1 ) ) ==  length( which( !is.na( outcomes$statins_prscd ) ) ) )
stopifnot( length( which( outcomes$bp_med_user == 1 ) ) ==  length( which( !is.na( outcomes$bp_med_prscd ) ) ) )

outcomes$gender %<>% as_factor(levels = "labels")
outcomes$ethnicity %<>% as_factor(levels = "labels")
outcomes$Townsend %<>% as_factor()

save(outcomes, file = (paste0(out_path,'outcomes_new.RData')))
