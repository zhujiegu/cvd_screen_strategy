#print(commandArgs(trailingOnly = TRUE))

library(data.table)
library(nlme)
library(survival)
library(tidyverse)
library(haven)
library(magrittr)

out_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/'
in_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/'

load(paste0(in_path,'outcomes_new.RData'))
load(paste0(in_path,'exposures_new.RData'))

# #keep variables of interest
# exp_oi = c("bmi","hdl", "tchol", "sbp", "smokbin")
# exposures_red = exposures[which( exposures$exposure %in% exp_oi ), ]

# #remove the exposures data
# rm(exposures); gc()

# Remove patients who do not have outcome data
exposures %<>% filter(patid %in% outcomes$patid)
stopifnot(all(unique(exposures$patid) %in% outcomes$patid))

sum(unique(exposures$patid) %in% outcomes$patid)
exposures_merged = merge( exposures, outcomes[,c("patid", "d_yob", "cvd_date", "end_date", "cvd_ind", "gender", "derivation") ], by = "patid")  
exposures_merged %<>% mutate(exp_age = as.numeric(exp_date-d_yob)/365.25)

#keep only the exposures before CVD date 
cvd_exp_cond = ( exposures_merged$cvd_date - exposures_merged$d_yob )/365.25 > exposures_merged$exp_age | is.na( exposures_merged$cvd_date )
exposures_merged = exposures_merged[which(cvd_exp_cond), ]
exposures_merged$patid %>% unique %>% length()

#keep only the exposures before end date 
exposures_merged %<>% filter(exp_age <= (end_date - d_yob)/365.25)
exposures_merged$patid %>% unique %>% length()

#creating the correct scaled variable
exposures_merged$scaled_corr = 0
exposures_merged$scaled_corr[exposures_merged$exposure == "smokbin"] = exposures_merged$original[exposures_merged$exposure == "smokbin"] 

female_id = which( exposures_merged$gender == "Female" )
male_id = which( exposures_merged$gender == "Male" )

exposures_merged %<>% as.data.frame
print("here1")

for(var in c("bmi","hdl", "tchol", "sbp")){
  # #creating scaled variable for the derivation set
  # ioi_female_deriv = which( exposures_merged$gender == "Female" & exposures_merged$exposure == var & exposures_merged$derivation == "derivation" )
  # ioi_male_deriv   = which( exposures_merged$gender == "Male"   & exposures_merged$exposure == var & exposures_merged$derivation == "derivation" )
  # 
  # female_mean_deriv = mean( exposures_merged[ioi_female_deriv, "original" ] ) 
  # male_mean_deriv   = mean( exposures_merged[ioi_male_deriv, "original" ] ) 
  # 
  # #overall_sd = sd( exposures_merged$original )
  # female_sd_deriv = sd( exposures_merged[ioi_female_deriv, "original" ] ) 
  # male_sd_deriv   = sd( exposures_merged[ioi_male_deriv, "original" ] ) 
  # 
  # exposures_merged$scaled_corr[ioi_female_deriv ] = (exposures_merged$original[ioi_female_deriv ] - female_mean_deriv)/female_sd_deriv
  # exposures_merged$scaled_corr[ioi_male_deriv ]   = (exposures_merged$original[ioi_male_deriv ] - male_mean_deriv)/male_sd_deriv
  # 
  # #creating scaled variable for the validation set
  # ioi_female_valid = which( exposures_merged$gender == "Female" & exposures_merged$exposure == var & exposures_merged$derivation == "validation" )
  # ioi_male_valid   = which( exposures_merged$gender == "Male"   & exposures_merged$exposure == var & exposures_merged$derivation == "validation" )
  # 
  # #standardize using mean and sd from derivation set
  # exposures_merged$scaled_corr[ioi_female_valid ] = (exposures_merged$original[ioi_female_valid ] - female_mean_deriv)/female_sd_deriv
  # exposures_merged$scaled_corr[ioi_male_valid ]   = (exposures_merged$original[ioi_male_valid ] - male_mean_deriv)/male_sd_deriv

  ioi_female = which( exposures_merged$gender == "Female" & exposures_merged$exposure == var)
  ioi_male   = which( exposures_merged$gender == "Male"   & exposures_merged$exposure == var)
  
  female_mean = mean( exposures_merged[ioi_female, "original" ] ) 
  male_mean   = mean( exposures_merged[ioi_male, "original" ] ) 
  
  #overall_sd = sd( exposures_merged$original )
  female_sd = sd( exposures_merged[ioi_female, "original"]) 
  male_sd   = sd( exposures_merged[ioi_male, "original"]) 
  
  exposures_merged$scaled_corr[ioi_female] = (exposures_merged$original[ioi_female] - female_mean)/female_sd
  exposures_merged$scaled_corr[ioi_male]   = (exposures_merged$original[ioi_male] - male_mean)/male_sd
  }

# align exposure and outcome patid
patid_unique <- unique(exposures_merged$patid)
stopifnot(all(patid_unique %in% outcomes$patid))
outcomes %<>% filter(patid %in% patid_unique)

print("People that died with CVD (during or after diagnosis)")
length( which( outcomes$cvd_ind == 1 & outcomes$death_ind == 1 ) )/length(outcomes$cvd_ind)
print("People that died without a diagnosis of CVD")
length( which( outcomes$cvd_ind == 0 & outcomes$death_ind == 1 ) )/length(outcomes$cvd_ind)
print("People censored with a diagnosis of CVD")
length( which( outcomes$cvd_ind == 1 & outcomes$death_ind == 0 ) )/length(outcomes$cvd_ind)
print("People censored without a diagnosis of CVD")
length( which( outcomes$cvd_ind == 0 & outcomes$death_ind == 0 ) )/length(outcomes$cvd_ind)
# table( exposures_merged$original[exposures_merged$exposure == "smokbin"] )
# table( exposures_merged$scaled[exposures_merged$exposure == "smokbin"] )
# table( exposures_merged$scaled_corr[exposures_merged$exposure == "smokbin"] )

# hist( exposures_merged$scaled[which( exposures_merged$exposure == "bmi" & exposures_merged$gender == "Female" )] )
# hist( exposures_merged$scaled_corr[which( exposures_merged$exposure == "bmi" & exposures_merged$gender == "Female" )] )
# 
# hist( exposures_merged$scaled[which( exposures_merged$exposure == "bmi" & exposures_merged$gender == "Male" )] )
# hist( exposures_merged$scaled_corr[which( exposures_merged$exposure == "bmi" & exposures_merged$gender == "Male" )] )
# 
# 
# hist( exposures_merged$scaled[which( exposures_merged$exposure == "hdl" & exposures_merged$gender == "Female" )] )
# hist( exposures_merged$scaled_corr[which( exposures_merged$exposure == "hdl" & exposures_merged$gender == "Female" )] )
# 
# hist( exposures_merged$scaled[which( exposures_merged$exposure == "hdl" & exposures_merged$gender == "Male" )] )
# hist( exposures_merged$scaled_corr[which( exposures_merged$exposure == "hdl" & exposures_merged$gender == "Male" )] )
# 
# 
# hist( exposures_merged$scaled[which( exposures_merged$exposure == "sbp" & exposures_merged$gender == "Female" )] )
# hist( exposures_merged$scaled_corr[which( exposures_merged$exposure == "sbp" & exposures_merged$gender == "Female" )] )
# 
# hist( exposures_merged$scaled[which( exposures_merged$exposure == "sbp" & exposures_merged$gender == "Male" )] )
# hist( exposures_merged$scaled_corr[which( exposures_merged$exposure == "sbp" & exposures_merged$gender == "Male" )] )
# 
# table( exposures_merged$scaled[which( exposures_merged$exposure == "smokbin" & exposures_merged$gender == "Female" )] )
# table( exposures_merged$scaled_corr[which( exposures_merged$exposure == "smokbin" & exposures_merged$gender == "Female" )] )
# 
# table( exposures_merged$scaled[which( exposures_merged$exposure == "smokbin" & exposures_merged$gender == "Male" )] )
# table( exposures_merged$scaled_corr[which( exposures_merged$exposure == "smokbin" & exposures_merged$gender == "Male" )] )
# 
# hist( exposures_merged$scaled[which( exposures_merged$exposure == "tchol" & exposures_merged$gender == "Female" )] )
# hist( exposures_merged$scaled_corr[which( exposures_merged$exposure == "tchol" & exposures_merged$gender == "Female" )] )
# 
# hist( exposures_merged$scaled[which( exposures_merged$exposure == "tchol" & exposures_merged$gender == "Male" )] )
# hist( exposures_merged$scaled_corr[which( exposures_merged$exposure == "tchol" & exposures_merged$gender == "Male" )] )

save(exposures_merged, file = (paste0(out_path,'exposures_merged_new_popmean.RData')))
# save(outcomes, file = (paste0(out_path,'outcomes.RData')))
