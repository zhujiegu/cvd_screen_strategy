#print(commandArgs(trailingOnly = TRUE))

library(data.table)
library(nlme)
library(survival)
library(tidyverse)
library(haven)

out_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/'
in_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/'

load(paste0(in_path,'patients_outcomes.RData'))
load(paste0(in_path,'/cprd_2M_exposures.RData'))

#keep variables of interest
exp_oi = c("bmi","hdl", "tchol", "sbp", "smokbin")
exposures_red = exposures[which( exposures$exposure %in% exp_oi ), ]

#remove the exposures data
rm(exposures); gc()

#remove unrealistic data
real_sbp     = which( exposures_red$exposure == "sbp"   & exposures_red$original >= 60   & exposures_red$original <= 250 )
real_tchol   = which( exposures_red$exposure == "tchol" & exposures_red$original >= 1.75 & exposures_red$original <= 20 )
real_hdl     = which( exposures_red$exposure == "hdl"   & exposures_red$original >= 0.3  & exposures_red$original <= 3.1 )
real_bmi     = which( exposures_red$exposure == "bmi"   & exposures_red$original <= 80 )
real_smokbin = which( exposures_red$exposure == "smokbin" )

real_values   = sort( c( real_sbp, real_tchol, real_hdl, real_bmi, real_smokbin ) )
exposures_red = exposures_red[real_values, ]

rm( real_bmi, real_sbp, real_smokbin, real_tchol, real_hdl, real_values )

stopifnot(all(unique(exposures_red$patid) %in% outcomes$patid))


#check on correctness CVD, DEATH and END times
# colnames(outcomes)
# outcomes[, c("cvd_date", "death_date", "end_date")]

#correct CVD post censoring
cvd_post_censoring = which( !is.na(outcomes$cvd_date) & outcomes$end_date < outcomes$cvd_date )
length( cvd_post_censoring )
outcomes$cvd_date[cvd_post_censoring ] = NA
outcomes$cvd_ind[cvd_post_censoring ] = 0

cvd_post_censoring = which( !is.na(outcomes$cvd_date) & outcomes$end_date < outcomes$cvd_date )
stopifnot( length( cvd_post_censoring ) == 0 )
rm( cvd_post_censoring )

#uniform the notation for indices of events:
length( which(is.na( outcomes$cvd_ind)))
length( which(is.na( outcomes$cvd_ind ) & is.na(outcomes$cvd_date)))

outcomes$cvd_ind[which(is.na( outcomes$cvd_ind)) ] = 0
length( which(is.na( outcomes$death_ind)))


#correct Death post censoring
death_post_censoring = which( !is.na(outcomes$death_date) & outcomes$end_date < outcomes$death_date )
length( death_post_censoring )
outcomes$death_date[death_post_censoring ] = NA
outcomes$death_ind[death_post_censoring ] = 0

death_post_censoring = which( !is.na(outcomes$death_date) & outcomes$end_date < outcomes$death_date )
stopifnot( length( death_post_censoring ) == 0 )
rm( death_post_censoring )

#correct CVD post Death
cvd_post_death = which( !is.na(outcomes$cvd_date) & outcomes$death_date < outcomes$cvd_date )
stopifnot( length( cvd_post_death ) == 0 )
rm(cvd_post_death)


#CVD/Statin/death before entering
cvd_before_start <- which( !is.na(outcomes$cvd_date) & outcomes$cvd_date < outcomes$start_date )
statin_before_start <- which(!is.na(outcomes$statins_prscd) & outcomes$statins_prscd < outcomes$start_date)
death_before_start <- which( !is.na(outcomes$death_date) & outcomes$death_date < outcomes$start_date )
length(unique(c(cvd_before_start, statin_before_start)))

outcomes <- outcomes[-c(cvd_before_start, statin_before_start), ]

stopifnot( length(which( !is.na(outcomes$cvd_date) & outcomes$cvd_date < outcomes$start_date )) == 0 )
stopifnot( length(which( !is.na(outcomes$statins_prscd) & outcomes$statins_prscd < outcomes$start_date )) == 0 )
stopifnot( length(which( !is.na(outcomes$death_date) & outcomes$death_date < outcomes$start_date )) == 0 )

#keep only the exposures before CVD date 
exposures_red = exposures_red[, c("patid", "exp_date", "exp_age", "exposure", "original", "scaled")]
exposures_red_merged = merge( exposures_red, outcomes[,c("patid", "d_yob", "cvd_date", "cvd_ind", "gender", "derivation") ], by = "patid")  
cvd_exp_cond = ( exposures_red_merged$cvd_date - exposures_red_merged$d_yob )/365.25 > exposures_red_merged$exp_age | is.na( exposures_red_merged$cvd_date )

exposures_red_merged = exposures_red_merged[which(cvd_exp_cond), ]

print("here0")

#creating the correct scaled variable
exposures_red_merged$scaled_corr = 0
exposures_red_merged$scaled_corr[exposures_red_merged$exposure == "smokbin"] = exposures_red_merged$original[exposures_red_merged$exposure == "smokbin"] 

female_id = which( exposures_red_merged$gender == "Female" )
male_id = which( exposures_red_merged$gender == "Male" )

print("here1")

for(var in c("bmi","hdl", "tchol", "sbp")){
  #creating scaled variable for the derivation set
  ioi_female_deriv = which( exposures_red_merged$gender == "Female" & exposures_red_merged$exposure == var & exposures_red_merged$derivation == "derivation" )
  ioi_male_deriv   = which( exposures_red_merged$gender == "Male"   & exposures_red_merged$exposure == var & exposures_red_merged$derivation == "derivation" )
  
  female_mean_deriv = mean( exposures_red_merged[ioi_female_deriv, "original" ] ) 
  male_mean_deriv   = mean( exposures_red_merged[ioi_male_deriv, "original" ] ) 
  
  #overall_sd = sd( exposures_red_merged$original )
  female_sd_deriv = sd( exposures_red_merged[ioi_female_deriv, "original" ] ) 
  male_sd_deriv   = sd( exposures_red_merged[ioi_male_deriv, "original" ] ) 
  
  exposures_red_merged$scaled_corr[ioi_female_deriv ] = (exposures_red_merged$original[ioi_female_deriv ] - female_mean_deriv)/female_sd_deriv
  exposures_red_merged$scaled_corr[ioi_male_deriv ]   = (exposures_red_merged$original[ioi_male_deriv ] - male_mean_deriv)/male_sd_deriv
  
  #creating scaled variable for the validation set
  ioi_female_valid = which( exposures_red_merged$gender == "Female" & exposures_red_merged$exposure == var & exposures_red_merged$derivation == "validation" )
  ioi_male_valid   = which( exposures_red_merged$gender == "Male"   & exposures_red_merged$exposure == var & exposures_red_merged$derivation == "validation" )
  
  #standardize using mean and sd from derivation set
  exposures_red_merged$scaled_corr[ioi_female_valid ] = (exposures_red_merged$original[ioi_female_valid ] - female_mean_deriv)/female_sd_deriv
  exposures_red_merged$scaled_corr[ioi_male_valid ]   = (exposures_red_merged$original[ioi_male_valid ] - male_mean_deriv)/male_sd_deriv
  }

# align exposure and outcome patid
patid_unique <- unique(exposures_red_merged$patid)
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
# table( exposures_red_merged$original[exposures_red_merged$exposure == "smokbin"] )
# table( exposures_red_merged$scaled[exposures_red_merged$exposure == "smokbin"] )
# table( exposures_red_merged$scaled_corr[exposures_red_merged$exposure == "smokbin"] )

# hist( exposures_red_merged$scaled[which( exposures_red_merged$exposure == "bmi" & exposures_red_merged$gender == "Female" )] )
# hist( exposures_red_merged$scaled_corr[which( exposures_red_merged$exposure == "bmi" & exposures_red_merged$gender == "Female" )] )
# 
# hist( exposures_red_merged$scaled[which( exposures_red_merged$exposure == "bmi" & exposures_red_merged$gender == "Male" )] )
# hist( exposures_red_merged$scaled_corr[which( exposures_red_merged$exposure == "bmi" & exposures_red_merged$gender == "Male" )] )
# 
# 
# hist( exposures_red_merged$scaled[which( exposures_red_merged$exposure == "hdl" & exposures_red_merged$gender == "Female" )] )
# hist( exposures_red_merged$scaled_corr[which( exposures_red_merged$exposure == "hdl" & exposures_red_merged$gender == "Female" )] )
# 
# hist( exposures_red_merged$scaled[which( exposures_red_merged$exposure == "hdl" & exposures_red_merged$gender == "Male" )] )
# hist( exposures_red_merged$scaled_corr[which( exposures_red_merged$exposure == "hdl" & exposures_red_merged$gender == "Male" )] )
# 
# 
# hist( exposures_red_merged$scaled[which( exposures_red_merged$exposure == "sbp" & exposures_red_merged$gender == "Female" )] )
# hist( exposures_red_merged$scaled_corr[which( exposures_red_merged$exposure == "sbp" & exposures_red_merged$gender == "Female" )] )
# 
# hist( exposures_red_merged$scaled[which( exposures_red_merged$exposure == "sbp" & exposures_red_merged$gender == "Male" )] )
# hist( exposures_red_merged$scaled_corr[which( exposures_red_merged$exposure == "sbp" & exposures_red_merged$gender == "Male" )] )
# 
# table( exposures_red_merged$scaled[which( exposures_red_merged$exposure == "smokbin" & exposures_red_merged$gender == "Female" )] )
# table( exposures_red_merged$scaled_corr[which( exposures_red_merged$exposure == "smokbin" & exposures_red_merged$gender == "Female" )] )
# 
# table( exposures_red_merged$scaled[which( exposures_red_merged$exposure == "smokbin" & exposures_red_merged$gender == "Male" )] )
# table( exposures_red_merged$scaled_corr[which( exposures_red_merged$exposure == "smokbin" & exposures_red_merged$gender == "Male" )] )
# 
# hist( exposures_red_merged$scaled[which( exposures_red_merged$exposure == "tchol" & exposures_red_merged$gender == "Female" )] )
# hist( exposures_red_merged$scaled_corr[which( exposures_red_merged$exposure == "tchol" & exposures_red_merged$gender == "Female" )] )
# 
# hist( exposures_red_merged$scaled[which( exposures_red_merged$exposure == "tchol" & exposures_red_merged$gender == "Male" )] )
# hist( exposures_red_merged$scaled_corr[which( exposures_red_merged$exposure == "tchol" & exposures_red_merged$gender == "Male" )] )

save( exposures_red_merged, file = (paste0(out_path,'exposures_red_merged.RData')))
save(outcomes, file = (paste0(out_path,'outcomes.RData')))
