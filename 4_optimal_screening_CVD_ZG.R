# #print(commandArgs(trailingOnly = TRUE))

args=commandArgs(trailingOnly = TRUE)
j      = as.integer(args[1])
gender = as.character(args[2]) #g can be male or female

print(j)
print(gender)

in_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/3_cox_outp/'
out_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/4_screening_outp/'
# j=65
# gender='male'

library(tidyr)
library(data.table)
library(dplyr)
source("~/epi_paper/survival_integral_nb_functions.R" )
select <- dplyr::select

load( paste0(in_path, "risk_pred_FE4_RIRS_",j,"_", gender, ".RData" ) ) 
risk_pred_now = risk_prediction[ !is.na(`5.y.risk.0`) ] #keeping only the derivation set
rm( risk_prediction )
load( paste0(in_path, "covariates_matrix_RIRS_",j,"_", gender, ".RData") )
covariates_now = covariates_matrix#[ derivation == "derivation", ]
rm( covariates_matrix )
load( paste0(in_path, "cum_haz_est_RIRS_",j,"_", gender, ".RData") )
cum_haz_now = cum_haz_est
rm( cum_haz_est )

print("ora qui 0")

#optimal_matrix = matrix( 0, length(lm_age), 5)
#comparing_schemes_list = list(NULL)
#comparing_next_fup_list = list(NULL)
colnames(risk_pred_now)
GLOB_horizon = 10

risk_pred_now$status_cvd_visible = risk_pred_now$cvd_ind
risk_pred_now$status_cvd_visible[ risk_pred_now$cvd_age > (j + GLOB_horizon + 5) ] = 0 # 14 = 9+5 which is the wider observed time-window.


baseline_risk = risk_pred_now$risk_y10
# hist(baseline_risk)
baseline_class = ifelse( baseline_risk <= 0.025, "green",
                         ifelse( baseline_risk <= 0.05, "yellow",
                                 ifelse( baseline_risk <= 0.075, "orange",
                                         ifelse( baseline_risk <= 0.10, "red", "darkred"))))
risk_pred_now$risk_class = factor(baseline_class, levels =  c("darkred", "red", "orange", "yellow", "green"))
risk_pred_now$end_corr = ifelse(is.na( risk_pred_now$cvd_age ), pmin( risk_pred_now$end_age, j + GLOB_horizon + 5 ),  risk_pred_now$cvd_age )


# seq_index_status = grep("y.risk.status.", colnames(risk_pred_now))
# colnames(risk_pred_now)[seq_index_status]= paste0("5.y.status.", 0:GLOB_horizon)

# seq_index_status = grep("y.status.", colnames(risk_pred_now))
seq_index_risk = grep("y.risk.", colnames(risk_pred_now))

########################################################################
# 5-year threshold
load(file = paste0("~/epi_paper/risk_ratio.RData"))
risk_ratio_j_gender <- risk_ratio %>% filter(lm_age==j, gender==!!gender) %>% select(risk_ratio, risk_class)
risk_pred_now %<>% left_join(risk_ratio_j_gender, by='risk_class') %>% mutate(threshold_y5=0.1/risk_ratio)

#Detecting the first time 5y CVD risk crosses the threshold
risk_pred_now %<>% rowwise() %>%
  mutate( cross_year = which(c_across(starts_with("5.y.risk.")) > threshold_y5)[1] - 1,  # Subtract 1 to align with year 0-10
  ) %>% ungroup()

risk_pred_now %<>%  rowwise() %>%
  mutate(
    risk_crossed = ifelse(is.na(cross_year), NA, get(paste0("5.y.risk.", cross_year))),
    risk_before_crossed = ifelse(is.na(cross_year), NA, get(paste0("5.y.risk.", ifelse(cross_year == 0, 0, cross_year - 1))))
  ) %>% ungroup()

#people that never crosses the threshold
risk_pred_now$cross_year %>% is.na %>% sum

risk_pred_now %<>% mutate(
  pred_time = case_when(
    cross_year == 0 ~ 0,
    is.na(cross_year) ~ GLOB_horizon,
    TRUE ~ (threshold_y5 - risk_before_crossed) / (risk_crossed - risk_before_crossed) + cross_year -1)
) %>% mutate(pred_status=ifelse(is.na(cross_year), 0, 1))


comparing_schemes <- risk_pred_now %>% select (patid, risk_class, pred_time, pred_status, threshold_y5, statin_age, cvd_age, cvd_ind, 
                                               death_age, death_ind, end_age)

save(comparing_schemes, file = paste0(out_path,"comparing_schemes_RIRS_apply_", j, "_", gender, ".RData") )
