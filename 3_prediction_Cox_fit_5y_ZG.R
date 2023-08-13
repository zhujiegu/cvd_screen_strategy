# #print(commandArgs(trailingOnly = TRUE))

args=commandArgs(trailingOnly = TRUE)
j      = as.integer(args[1])
gender = as.character(args[2]) #g can be male or female

in_path2 <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/2_lmfit_outp/'
in_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/data_created/'
out_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/3_cox_outp/'
# j=75
# gender='male'

library(data.table)
library(nlme)
library(survival)
library(tidyverse)
library(haven)
library(feather)
library(dplyr)
library(miceadds)
library(dplyr)
library(pec)
library(rms)
library(riskRegression)
library(risksetROC)
library(survAUC)
library(dynpred)

source("~/epi_paper/mixoutsamp_v3.R")
source("~/epi_paper/pred_error_check.R")

load(paste0(in_path, gender, "/data_lm_",j, "_", gender,".RData"))
mod_type = "RIRS_"
load(paste0(in_path2,j,"_", gender,".RData"))
#######################
#
# Type of analysis
#
#######################

#Number of investigated outcomes
outcomes_names = c('bmi','hdl','sbp', 'smokbin','tchol')
N_outcomes = length(outcomes_names)

#GLOBAL Horizon
GLOB_horizon = 8


#N subjects
N_patid_all = length(unique(data_ext$patid)) 

patid_all = unique(data_ext$patid) 

#setting covariates for the cox model
#no more SLE ind and family ind nor Dementia
cox_cov_less60 = c("bp_bin + diab_ind + hdl_blup_FE4 + sbp_blup_FE4 + smoke_blup_FE4 + 
                   tchol_blup_FE4 + bmi_blup_FE4 + Townsend  +  depression_ind + 
                   Severe_mental_illness_ind + Migraine_ind ")
cox_cov_more60 = c("bp_bin + diab_ind + hdl_blup_FE4 + sbp_blup_FE4 + smoke_blup_FE4 + 
                   tchol_blup_FE4 + bmi_blup_FE4 + Townsend + renal_disease + 
                   atrial_fibrillation + rheumatoid_arthritis +  depression_ind + 
                   Severe_mental_illness_ind + Migraine_ind")


if(j >= 60){
  cox_cov = cox_cov_more60
  cox_cov_list = c("bp_bin", "diab_ind", "hdl_blup_FE4", "sbp_blup_FE4", "smoke_blup_FE4", 
                   "tchol_blup_FE4", "bmi_blup_FE4", "Townsend", "renal_disease", 
                   "atrial_fibrillation", "rheumatoid_arthritis",   "depression_ind", 
                   "Severe_mental_illness_ind",  "Migraine_ind") #"Dementia_ind"
}else{
  cox_cov = cox_cov_less60
  cox_cov_list = c("bp_bin", "diab_ind", "hdl_blup_FE4", "sbp_blup_FE4", "smoke_blup_FE4", "tchol_blup_FE4", 
                   "bmi_blup_FE4", "Townsend",  "depression_ind", "Severe_mental_illness_ind", 
                   "Migraine_ind" ) # bp_bin removed
}

N_cox_cov = length(cox_cov_list)

#######################
#
# Risk prediction est
#
#######################

risk_prediction =  data_ext[ exp_count == 1, .(patid, exp_age, statin_age, cvd_ind, cvd_age, death_ind, death_age, end_age) ]

#at the baseline nobody was still prescribed statins

#######################
#
# Cumulative hazard est
#
#######################

cum_haz_est = matrix(0, nrow = length(unique(data_ext$exp_age)), ncol = 6)
cum_haz_est = as.data.table(cum_haz_est)

#######################
#
# Covariates matrix
#
#######################

covariates_matrix = as.data.table(unique(data_ext$patid))
#colnames(covariates_matrix) = "patid"

cov_rec_before_lm_name = c("bp_bin", "statin_bin", "renal_disease",
                           "atrial_fibrillation", "rheumatoid_arthritis",
                           "Severe_mental_illness_ind", "Migraine_ind", "depression_ind","diab_ind")

cov_rec_before_lm = c("bp_med_age", "statin_age", "renal_age",
                       "Atrial_fibrillation_age", "rheumatoid_arthritis_age",
                       "Severe_mental_illness_age", "Migraine_age", "depression_age","diab_age")

covariates_matrix = data_ext[exp_count == 1, lapply(.SD, function(x) ifelse(x <= j & !is.na(x), 1, 0)),
                             .SDcols = cov_rec_before_lm]
colnames(covariates_matrix) = cov_rec_before_lm_name
covariates_matrix = cbind(covariates_matrix, 
                          data_ext[exp_count == 1, 
                                   .(patid, Townsend, family_history, derivation) ])

setcolorder(covariates_matrix, c( "patid", "derivation", "diab_ind", "bp_bin",
                                    "statin_bin", "Townsend",
                                    "renal_disease", "atrial_fibrillation",
                                    "rheumatoid_arthritis", "Severe_mental_illness_ind",
                                    "Migraine_ind", 
                                    "depression_ind", "family_history"))


#adding space for predicted time-varying covariates
covariates_matrix = cbind(covariates_matrix, 
                          as.data.table(matrix(0, nrow = length(unique(data_ext$patid)), 
                                               ncol = N_outcomes*(GLOB_horizon + 1))))
colnames(covariates_matrix) = c("patid", "derivation", "diab_ind", "bp_med", "statin_bin",
                                   'Townsend', 'renal_disease','atrial_fibrillation','rheumatoid_arthritis',
                                   'Severe_mental_illness_ind','Migraine_ind', 'depression_ind', 'family_history',
                                   'bmi_blup_0', 'hdl_blup_0',  'sbp_blup_0', 'smoke_blup_0', 'tchol_blup_0',
                                   'bmi_blup_1', 'hdl_blup_1',   'sbp_blup_1',  'smoke_blup_1', 'tchol_blup_1',
                                   'bmi_blup_2', 'hdl_blup_2',   'sbp_blup_2',  'smoke_blup_2', 'tchol_blup_2',
                                   'bmi_blup_3', 'hdl_blup_3',   'sbp_blup_3',  'smoke_blup_3', 'tchol_blup_3',
                                   'bmi_blup_4', 'hdl_blup_4',   'sbp_blup_4',  'smoke_blup_4', 'tchol_blup_4',
                                   'bmi_blup_5', 'hdl_blup_5',   'sbp_blup_5',  'smoke_blup_5', 'tchol_blup_5',
                                   'bmi_blup_6', 'hdl_blup_6',   'sbp_blup_6',  'smoke_blup_6', 'tchol_blup_6',
                                   'bmi_blup_7', 'hdl_blup_7',   'sbp_blup_7',  'smoke_blup_7', 'tchol_blup_7',
                                   'bmi_blup_8', 'hdl_blup_8',   'sbp_blup_8',  'smoke_blup_8', 'tchol_blup_8'
)




pt_analysed = matrix(0, nrow = GLOB_horizon + 1, ncol = 11)


#computing fixed part of prediction
length(unique(data_ext$patid[ which(data_ext$exp_age_corr <= 0) ]))
#SOME PEOPLE HAVE OBSERVATION ONLY IN THE FUTURE!!!! (25% at 60 yo)
length(unique(data_ext$patid[ which(data_ext$exp_age_corr <= 0) ]))/ length(unique(data_ext$patid))

#All the predictions will be carried out from t_0 = Landamark age
dati_baseline = data.table(
  patid = unique(data_ext$patid),
  #bp_bin_base = as.numeric(ifelse(data_ext[ exp_count == 1, !is.na(bp_med_age) & bp_med_age  <= j ], 1, 0)),
  #statin_bin_base = as.numeric(ifelse(data_ext[ exp_count == 1, !is.na(statin_age) & statin_age  <= j ], 1, 0)),
  diab_ind = data_ext[ exp_count == 1, ifelse(diab_age <= j & !is.na(diab_age), 1, 0) ],
  bp_bin = data_ext[ exp_count == 1, ifelse(bp_med_age <= j & !is.na(bp_med_age), 1, 0) ],
  statin_ind = data_ext[ exp_count == 1, ifelse(statin_age <= j  & !is.na(statin_age), 1, 0) ],
  Townsend = data_ext[ exp_count == 1, Townsend ],
  renal_disease = data_ext[ exp_count == 1, ifelse(renal_age <= j  & !is.na(renal_age), 1, 0)  ],
  atrial_fibrillation = data_ext[ exp_count == 1, ifelse(Atrial_fibrillation_age <= j  & !is.na(Atrial_fibrillation_age), 1, 0)  ],
  rheumatoid_arthritis = data_ext[ exp_count == 1, ifelse(rheumatoid_arthritis_age <= j  & !is.na(rheumatoid_arthritis_age), 1, 0)  ],
  Severe_mental_illness_ind = data_ext[ exp_count == 1, ifelse(Severe_mental_illness_age <= j  & !is.na(Severe_mental_illness_age), 1, 0)  ],
  Migraine_ind = data_ext[ exp_count == 1, ifelse(Migraine_age <= j  & !is.na(Migraine_age), 1, 0)  ],
  depression_ind = data_ext[ exp_count == 1, ifelse(depression_age <= j  & !is.na(depression_age), 1, 0)  ],
  family_history = data_ext[ exp_count == 1, "family_history" ],
  status_cvd_death = -1000,
  status_death = -1000,
  status_cvd = -1000,
  status_statin = -1000,
  time_cvd_death = -1000,
  time_death = -1000,
  time_cvd = -1000,
  time_statin = -1000, 
  cvd_age_glob_hor = -1000,
  cvd_status_glob_hor = -1000,
  death_age_glob_hor = -1000,
  death_status_glob_hor = -1000,
  composite_age_glob_hor = -1000,
  composite_status_glob_hor = -1000,
  derivation = data_ext[ exp_count ==1, derivation ],
  lm_subcohort_index = 1)

#ref_all is composed of RE related only to those people that have measures 
# before landmark age j
ref_all = mixoutsamp(fit_rirs, as.data.frame(data_ext[exp_age <= j,c("patid","exp_age_corr","exp_age",
                                                                      "exposure", "sbp_ind","tchol_ind",
                                                                      "bp_bin", "statin_bin", "scaled_corr")]))$random

#setting as 0 al RE that are associated to those people that
# have no measurements before landmark age j
ref_all_complete = data_ext %>%
  distinct(patid, derivation) %>%
  left_join(ref_all, by = "patid") %>%
  mutate_if(is.numeric, ~replace_na(.,0))

for(l in 0:GLOB_horizon){
  blup_BMI_tmp = fit_rirs$coefficients$fixed["exposurebmi"] + ref_all_complete[,"exposurebmi"] +
    (fit_rirs$coefficients$fixed["exposurebmi:exp_age_corr"] + 
       ref_all_complete[,"exposurebmi:exp_age_corr"])  * l 

  blup_HDL_tmp =  fit_rirs$coefficients$fixed["exposurehdl"] + ref_all_complete[,"exposurehdl"] +
    (fit_rirs$coefficients$fixed["exposurehdl:exp_age_corr"] + ref_all_complete[,"exposurehdl:exp_age_corr"])  * l 

  blup_SBP_tmp = fit_rirs$coefficients$fixed["exposuresbp"] + 
    fit_rirs$coefficients$fixed["bp_bin:sbp_ind"] * dati_baseline[,bp_bin] + 
    ref_all_complete[,"exposuresbp"] +
    (fit_rirs$coefficients$fixed["exposuresbp:exp_age_corr"] + 
       ref_all_complete[,"exposuresbp:exp_age_corr"])  * l 

  blup_SMOKE_tmp = fit_rirs$coefficients$fixed["exposuresmokbin"] + ref_all_complete[,"exposuresmokbin"] +
    (fit_rirs$coefficients$fixed["exposuresmokbin:exp_age_corr"] + 
       ref_all_complete[,"exposuresmokbin:exp_age_corr"])  * l 

  blup_TCHOL_tmp = fit_rirs$coefficients$fixed["exposuretchol"] + 
    fit_rirs$coefficients$fixed["statin_bin:tchol_ind"] * dati_baseline[,statin_ind] + 
    ref_all_complete[,"exposuretchol"] +
    (fit_rirs$coefficients$fixed["exposuretchol:exp_age_corr"] + 
       ref_all_complete[,"exposuretchol:exp_age_corr"])  * l 

  blups_tmp = cbind(blup_BMI_tmp, blup_HDL_tmp, blup_SBP_tmp, blup_SMOKE_tmp, blup_TCHOL_tmp)
  head(blups_tmp)

  idx_cov = which(endsWith(colnames(covariates_matrix), paste0("_", l)))
  covariates_matrix[ , idx_cov ] = blups_tmp
}


for(l in 0:GLOB_horizon)
{
  time_origin = j + l
  time_horizon = j + l + 5
  
  #remove who died or who experienced cvd event before time_origin
  patid_out_time = data_ext[ (!is.na(death_date) & death_age <= time_origin) | 
                               (!is.na(cvd_date) & cvd_age  <= time_origin) | 
                               end_age <= time_origin , patid ]
  
  data_ext_in_time = data_ext[ !patid %in% patid_out_time, ]
  
  patid_in_time = unique(data_ext_in_time$patid)
  n_patid = length(unique(data_ext_in_time$patid))
  
  pt_analysed[ l+1, 1 ] = n_patid
  pt_analysed[ l+1, 2 ] = data_ext_in_time[ exp_count == 1, length(which(death_age <= time_origin + 1 & (cvd_age > time_origin + 1| is.na(cvd_age)))) ]
  pt_analysed[ l+1, 3 ] = data_ext_in_time[ exp_count == 1, length(which(cvd_age <= time_origin + 1 & (death_age > time_origin + 1| is.na(death_age))))]
  pt_analysed[ l+1, 4 ] = data_ext_in_time[ exp_count == 1, length(which(cvd_age <= time_origin + 1 & death_age <= time_origin + 1))]
  pt_analysed[ l+1, 5 ] = data_ext_in_time[ exp_count == 1, length(which(end_age <= time_origin + 1 & (death_age > time_origin + 1| is.na(death_age)) & (cvd_age > time_origin + 1| is.na(cvd_age)))) ]
  pt_analysed[ l+1, 6 ] = data_ext_in_time[ exp_count == 1, length(which(death_age <= time_origin + 1 |
                                                                             cvd_age <= time_origin + 1 |
                                                                             end_age <= time_origin + 1)) ]
  pt_analysed[ l+1, 7 ] = data_ext_in_time[ exp_count == 1, length(which(death_age <= time_horizon & (cvd_age > time_horizon | is.na(cvd_age))))]
  pt_analysed[ l+1, 8 ] = data_ext_in_time[ exp_count == 1, length(which(cvd_age <= time_horizon & (death_age > time_horizon | is.na(death_age)))) ]
  pt_analysed[ l+1, 9 ] = data_ext_in_time[ exp_count == 1, length(which(cvd_age <= time_horizon & death_age <= time_horizon)) ]
  pt_analysed[ l+1, 10 ] = data_ext_in_time[ exp_count == 1, length(which(end_age <= time_horizon & (death_age > time_horizon| is.na(death_age)) & (cvd_age > time_horizon| is.na(cvd_age)))) ]
  pt_analysed[ l+1, 11 ] = data_ext_in_time[ exp_count == 1, length(which(death_age <= time_horizon |
                                                                              cvd_age <= time_horizon |
                                                                              end_age <= time_horizon)) ]
  
  
   
  #CREATING DATA for COX model
  
  #setting the outcomes of interest
  outcome_time_death  = data_ext_in_time[exp_count == 1, lapply(.SD, function(x){
    ifelse(is.na(x), min(time_horizon, end_age) - time_origin, min(time_horizon, x) - time_origin)
  }), .SDcols = 'death_age',  by = patid ]
  outcome_time_cvd    = data_ext_in_time[exp_count == 1, lapply(.SD, function(x){
    ifelse(is.na(x), min(time_horizon, end_age) - time_origin, min(time_horizon, x) - time_origin)
  }), .SDcols = 'cvd_age',    by = patid ]
  outcome_time_statin = data_ext_in_time[ exp_count == 1, lapply(.SD, function(x) ifelse(is.na(x), min(time_horizon, end_age) - time_origin, min(time_horizon, x) - time_origin)), .SDcols = 'statin_age', by = patid ]
  
  outcome_status_death  = data_ext_in_time[exp_count == 1, lapply(.SD, function(x){
    ifelse(x <= time_horizon & x > time_origin & !is.na(x), 1, 0)
  }), .SDcols = 'death_age' ]
  outcome_status_cvd    = data_ext_in_time[ exp_count == 1, lapply(.SD, function(x) ifelse(x <= time_horizon & x > time_origin & !is.na(x), 1, 0)), .SDcols = 'cvd_age' ]
  outcome_status_statin = data_ext_in_time[ exp_count == 1, lapply(.SD, function(x) ifelse(x <= time_horizon & x > time_origin & !is.na(x), 1, 0)), .SDcols = 'statin_age' ]
  
  outcome_time_death_cvd   = apply(cbind(outcome_time_cvd$cvd_age, outcome_time_death$death_age), 1, min)
  outcome_status_death_cvd = rep(0, length(outcome_time_cvd$cvd_age))
  #Fatal cvd
  outcome_status_death_cvd[outcome_time_cvd$cvd_age == outcome_time_death$death_age & 
                             outcome_status_cvd$cvd_age == 1 & 
                             outcome_status_death$death_age == 1] = 2
  #Non fatal CVD
  outcome_status_death_cvd[outcome_time_cvd$cvd_age < outcome_time_death$death_age & 
                             outcome_status_cvd$cvd_age == 1] = 1
  ############################################
  # # Disagreement
  # sum(outcome_time_cvd$cvd_age < outcome_time_death$death_age & 
  #       outcome_status_cvd$cvd_age == 1)
  # 
  # # note here that when death is NA, the death_age is the end_age (see Line 290)
  # # There are people who's end_age is the same as cvd_age
  # sum(outcome_time_cvd$cvd_age == outcome_time_death$death_age & 
  #       outcome_status_cvd$cvd_age == 1 &  outcome_status_death$death_age == 0)
  # 
  # sum((outcome_time_cvd$cvd_age < outcome_time_death$death_age & outcome_status_cvd$cvd_age == 1)|
  #       (outcome_time_cvd$cvd_age == outcome_time_death$death_age & 
  #          outcome_status_cvd$cvd_age == 1 &  outcome_status_death$death_age == 0))
  
  outcome_status_death_cvd[(outcome_time_cvd$cvd_age < outcome_time_death$death_age & outcome_status_cvd$cvd_age == 1)|
                             (outcome_time_cvd$cvd_age == outcome_time_death$death_age & 
                                outcome_status_cvd$cvd_age == 1 &  outcome_status_death$death_age == 0)] = 1
  ############################################
  #Death for other causes
  outcome_status_death_cvd[  outcome_status_death$death_age == 1 & outcome_status_cvd$cvd_age == 0  ] = 3
  
  #SETTING AS 0 THOSE PEOPLE THAT DO NOT BELONG TO THE LANDMARK SUBCOHORT
  dati_baseline[ !patid %in% patid_in_time, 'lm_subcohort_index' ]      = 0
  dati_baseline[ patid %in% patid_in_time, 'status_cvd_death' ]  = outcome_status_death_cvd
  dati_baseline[ patid %in% patid_in_time, 'status_death'  ]     = outcome_status_death
  dati_baseline[ patid %in% patid_in_time, 'status_cvd' ]       = outcome_status_cvd
  dati_baseline[ patid %in% patid_in_time, 'status_statin' ]    = outcome_status_statin
  dati_baseline[ patid %in% patid_in_time, 'time_cvd_death' ]   = outcome_time_death_cvd
  dati_baseline[ patid %in% patid_in_time, 'time_death' ]       = outcome_time_death$death_age
  dati_baseline[ patid %in% patid_in_time, 'time_cvd' ]         = outcome_time_cvd$cvd_age
  dati_baseline[ patid %in% patid_in_time, 'time_statin' ]      = outcome_time_statin$statin_age
  
  idx_smoke = which(endsWith( colnames(covariates_matrix), paste0("smoke_blup_", l)))
  idx_hdl = which(endsWith( colnames(covariates_matrix), paste0("hdl_blup_", l)))
  idx_sbp = which(endsWith( colnames(covariates_matrix), paste0("sbp_blup_", l)))
  idx_tchol = which(endsWith( colnames(covariates_matrix), paste0("tchol_blup_", l)))
  idx_bmi = which(endsWith( colnames(covariates_matrix), paste0("bmi_blup_", l)))
  
  dati_baseline[, 'smoke_blup_FE4' ] = covariates_matrix[ , ..idx_smoke ] 
  dati_baseline[, 'hdl_blup_FE4' ]   = covariates_matrix[ , ..idx_hdl ]
  dati_baseline[, 'sbp_blup_FE4'   ] = covariates_matrix[ , ..idx_sbp ]
  dati_baseline[, 'tchol_blup_FE4' ] = covariates_matrix[ , ..idx_tchol ]
  dati_baseline[, 'bmi_blup_FE4'  ]  = covariates_matrix[ , ..idx_bmi ]
  
  ##################################
  ## explore the distribution of SBP within BP medication groups
  # library(ggplot2)
  # dati_baseline %>% group_by(bp_bin) %>% 
  #   summarise(mean=mean(sbp_blup_FE4), median=median(sbp_blup_FE4), sd=sd(sbp_blup_FE4))
  # pdf('temp.pdf')
  # ggplot(data = dati_baseline, aes(x=sbp_blup_FE4, color=as.factor(bp_bin)))+ geom_density()
  # dev.off()
  ##################################
  
  cox_formula_5CVD = as.formula(paste("Surv(time_cvd, status_cvd)~", cox_cov))
  model_cox_5CVD = coxph(cox_formula_5CVD, data = dati_baseline[lm_subcohort_index == 1,], x = T)
  #summary(model_cox_5CVD_new)
  
  # cox_formula_5CVD = as.formula(paste("Surv(time_cvd, status_cvd)~", cox_cov))
  # model_cox_5CVD = coxph(cox_formula_5CVD, data = data_cox[derivation == "derivation",], x = T)
  # summary(model_cox_5CVD)
  # 
  
  save(model_cox_5CVD, file = paste0(out_path,"cox_model_", mod_type, j, "_", l,"_", gender, ".RData"))
  
  print(paste(l, "Cox model fitted"))
  
  base_haz_5CVD  = basehaz(model_cox_5CVD, center = T)
  
  #Estimates for 5-year CVD risk prediction at landmark age (= on the whole landmark cohort) 
  lp_der_5CVD_all_pt_lm = predict(model_cox_5CVD, newdata = dati_baseline[, ], type = "lp")    
  estimated_surv = exp(- matrix(base_haz_5CVD$hazard, ncol = 1) %*% matrix(exp(lp_der_5CVD_all_pt_lm), nrow = 1))
  dim(estimated_surv)
  
  
  # if(l == 0){
  #   #CVD outcome
  #   cvd_age_glob_hor = data_ext_in_time[ exp_count == 1, lapply(.SD, function(x) ifelse(is.na(x), min(j + GLOB_horizon, end_age) - j, min(j + GLOB_horizon, x) -j)), .SDcols = 'cvd_age', by = patid ]
  #   cvd_status_glob_hor = data_ext_in_time[ exp_count == 1, lapply(.SD, function(x) ifelse(x <= (j + GLOB_horizon) & x > j & !is.na(x), 1, 0)), .SDcols = 'cvd_age' ]
  #   dati_baseline[ patid %in% patid_in_time, 'cvd_status_glob_hor'] = cvd_status_glob_hor# outcome_status_cvd_cumhaz
  #   dati_baseline[ patid %in% patid_in_time, 'cvd_age_glob_hor'] = cvd_age_glob_hor$cvd_age #outcome_time_cvd_cumhaz$cvd_age
  #   cox_formula_10CVD = as.formula(paste("Surv(cvd_age_glob_hor, cvd_status_glob_hor)~", cox_cov))
  #   model_cox_10CVD = coxph(cox_formula_10CVD, data = dati_baseline[lm_subcohort_index == 1,], x = T)
  #   base_haz_10CVD = basehaz(model_cox_10CVD, centered = T)
  #   cum_haz_est[ 1:dim(base_haz_10CVD)[ 1 ], 1:2 ] = base_haz_10CVD
  #   beta_l0_cvd = model_cox_10CVD$coef 
  #   
  #   save(beta_l0_cvd, file = paste0(out_path,"beta_l0_cvd_outcome_", mod_type, j,"_", gender, ".RData"))
  # 
  #   #Death outcome
  #   death_age_glob_hor = data_ext_in_time[ exp_count == 1, lapply(.SD, function(x) ifelse(is.na(x), min(j + GLOB_horizon, end_age) - j, min(j + GLOB_horizon, x) -j)), .SDcols = 'death_age', by = patid ]
  #   death_status_glob_hor = data_ext_in_time[ exp_count == 1, lapply(.SD, function(x) ifelse(x <= (j + GLOB_horizon) & x > j & !is.na(x), 1, 0)), .SDcols = 'death_age' ]
  #   dati_baseline[ patid %in% patid_in_time, 'death_status_glob_hor'] = death_status_glob_hor# outcome_status_death_cumhaz
  #   dati_baseline[ patid %in% patid_in_time, 'death_age_glob_hor'] = death_age_glob_hor$death_age #outcome_time_death_cumhaz$death_age
  #   cox_formula_10death = as.formula(paste("Surv(death_age_glob_hor, death_status_glob_hor)~", cox_cov))
  #   model_cox_10death = coxph(cox_formula_10death, data = dati_baseline[lm_subcohort_index == 1,], x = T)
  #   base_haz_10death = basehaz(model_cox_10death, centered = T)
  #   cum_haz_est[ 1:dim(base_haz_10death)[1], 3:4 ] = base_haz_10death
  #   beta_l0_death = model_cox_10death$coef 
  #   
  #   save(beta_l0_death, file = paste0(out_path,"beta_l0_death_outcome_", mod_type, j,"_", gender, ".RData"))
  # 
  #   #Composite outcome
  #   death_age_glob_hor = data_ext_in_time[ exp_count == 1, lapply(.SD, function(x) ifelse(is.na(x), min(j + GLOB_horizon, end_age) - j, min(j + GLOB_horizon, x) - j)), .SDcols = 'death_age', by = patid ]
  #   death_status_glob_hor = data_ext_in_time[ exp_count == 1, lapply(.SD, function(x) ifelse(x <= (j + GLOB_horizon) & x > j & !is.na(x), 1, 0)), .SDcols = 'death_age' ]
  #   
  #   composite_age_glob_hor = apply(cbind(cvd_age_glob_hor$cvd_age, death_age_glob_hor$death_age), 1, min)
  #   #composite_status_glob_hor = ifelse(cvd_status_glob_hor == 1, 1, ifelse(death_status_glob_hor == 1, 2, 0))
  #   composite_status_glob_hor = ifelse(cvd_status_glob_hor == 1 | death_status_glob_hor == 1, 1, 0)
  #   
  #   
  #   dati_baseline[ patid %in% patid_in_time, 'composite_status_glob_hor']  = composite_status_glob_hor 
  #   dati_baseline[ patid %in% patid_in_time, 'composite_age_glob_hor']  =  composite_age_glob_hor  
  #   cox_formula_10composite = as.formula(paste("Surv(composite_age_glob_hor, composite_status_glob_hor)~", cox_cov))
  #   model_cox_10composite = coxph(cox_formula_10composite, data = dati_baseline[lm_subcohort_index == 1,], x = T)
  #   base_haz_10composite = basehaz(model_cox_10composite, centered = T)
  #   cum_haz_est[ 1:dim(base_haz_10composite)[1], 5:6 ] = base_haz_10composite
  #   beta_l0_composite = model_cox_10composite$coef 
  #   
  #   save(beta_l0_composite, file = paste0(out_path,"beta_l0_composite_outcome_", mod_type, j,"_", gender, ".RData"))
  # }
  
  stopifnot(min(dati_baseline[ patid %in% patid_in_time, time_cvd ]) > 0)
 
  five_year_risk = estimated_surv[ dim(estimated_surv)[1],  ] #selecting the 5-year CVD risk
  risk_prediction[ patid %in% dati_baseline[, patid], paste0("5.y.risk.",l) := five_year_risk  ]
  risk_prediction[ patid %in% dati_baseline[, patid], paste0("5.y.risk.status.",l) := ifelse(five_year_risk <= 0.95, 1, 0) ]
  # hitting_time = apply(estimated_surv, 2, function(x) which(x  <= 0.95)[1] )
}


colnames(cum_haz_est)= c("hazard_cvd_outcome", "time_cvd_outcome", "hazard_death_outcome", "time_death_outcome", "hazard_composite_outcome", "time_composite_outcome")


warnings()

save(risk_prediction, file = paste0(out_path, "risk_pred_FE4_", mod_type, j,"_", gender, ".RData"))
save(cum_haz_est, file = paste0(out_path,"cum_haz_est_", mod_type, j,"_", gender, ".RData"))
save(covariates_matrix, file = paste0(out_path, "covariates_matrix_", mod_type, j,"_", gender, ".RData"))
