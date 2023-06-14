# #print(commandArgs(trailingOnly = TRUE))

args=commandArgs(trailingOnly = TRUE)
i      = as.integer(args[1])
gender = as.character(args[2]) #g can be male or female

in_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/3_cox_outp/'
out_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/4_screening_outp/'
# i=45
# gender='male'

library(tidyr)
library(data.table)
source("~/epi_paper/survival_integral_nb_functions.R" )


load( paste0(in_path, "risk_pred_FE4_RIRS_",i,"_", gender, ".RData" ) ) 
risk_pred_now = risk_prediction[ !is.na(`5.y.risk.0`) ] #keeping only the derivation set
rm( risk_prediction )
load( paste0(in_path, "covariates_matrix_RIRS_",i,"_", gender, ".RData") )
covariates_now = covariates_matrix#[ derivation == "derivation", ]
rm( covariates_matrix )
load( paste0(in_path, "cum_haz_est_RIRS_",i,"_", gender, ".RData") )
cum_haz_now = cum_haz_est
rm( cum_haz_est )


load(paste0(in_path,"beta_l0_cvd_outcome_RIRS_",i,"_", gender, ".RData") )

print("ora qui 0")


#setting parameters
# count = 1
# 
# bs = 60
# cv = 20
# cs = 20
# theta = 0.8
# 
# a_vec = seq(0.25, 10, length.out = 40)




#optimal_matrix = matrix( 0, length(lm_age), 5)
#comparing_schemes_list = list(NULL)
#comparing_next_fup_list = list(NULL)
colnames(risk_pred_now)
GLOB_horizon = 10

risk_pred_now$status_cvd_visible = risk_pred_now$cvd_ind
risk_pred_now$status_cvd_visible[ risk_pred_now$cvd_age > (i + GLOB_horizon + 5) ] = 0 # 14 = 9+5 which is the wider observed time-window.


baseline_risk = 1 - risk_pred_now$`5.y.risk.0` 
hist(baseline_risk)
baseline_class = ifelse( baseline_risk <= 0.025/2, "Low risk",
                         ifelse( baseline_risk <= 0.05/2, "Med-low risk",
                                 ifelse( baseline_risk <= 0.075/2, "Med-high risk",
                                         ifelse( baseline_risk <= 0.10/2, "High risk", "Very high risk"))))
risk_pred_now$risk_class = factor(baseline_class, levels =  c("Very high risk", "High risk", "Med-high risk", "Med-low risk", "Low risk"))
risk_pred_now$end_corr = ifelse(is.na( risk_pred_now$cvd_age ), pmin( risk_pred_now$end_age, i + GLOB_horizon + 5 ),  risk_pred_now$cvd_age )


seq_index_status = grep("y.risk.status.", colnames(risk_pred_now))
colnames(risk_pred_now)[seq_index_status]= paste0("5.y.status.", 0:10)

seq_index_status = grep("y.status.", colnames(risk_pred_now))
seq_index_risk = grep("y.risk.", colnames(risk_pred_now))

#Detecting the first status = 1, that is the first time 
# 5y CVD risk crosses the threshold
pred_time = apply( risk_pred_now[ , ..seq_index_status ], 1, function(x) which( x == 1 )[1] )
#people that crosses the theshold
length(which(!is.na(pred_time)))
#people that never crosses the threshold
length(which(is.na(pred_time)))

pred_risk = apply( risk_pred_now[ , ..seq_index_risk ], 1, function(x) x[ which( x < 0.95 )[1] ] )

stopifnot( length(which(!is.na(pred_risk))) == length(which(!is.na(pred_time))) )

risk_pred_now$pred_status = ifelse( is.na(pred_time), 0, 1)

#RISK is not monotonic!!!
pred_risk_before = apply(risk_pred_now[!is.na( pred_time ) & pred_time != 1 , 
                                        ..seq_index_risk ], 1, function(x) x[ which( x < 0.95 )[1] - 1 ] )


# #visualization of risk
# seq_to_melt = c(1,seq_index_risk)
# risk_pred_now_melt = reshape2::melt(risk_pred_now[sample(1:90000,1000, replace = F),..seq_to_melt], id = "patid")
# risk_pred_now_melt$value = 1 - risk_pred_now_melt$value
# ggplot( risk_pred_now_melt, aes(x = variable, y = as.factor(patid), fill = value))+
#   geom_tile() +
#   labs( x= "", y ="", title = paste0("Landmark ", i)) +
#   scale_fill_gradient(low="white", high="red") +
#   theme(axis.text.y = element_blank())
#   
# 
# length(pred_risk_before)
# pred_risk_before[1:10]
# names_cross =  paste0("5.y.risk.", pred_time[!is.na(pred_time)]-1)
# idx_cross = match( names_cross, colnames(risk_pred_now) )
# complete_matr =cbind(which(!is.na(pred_time)),idx_cross )
# risk_pred_now[ complete_matr[,1], complete_matr[,2], with = F]

# pred_risk_B = 3*( pred_time - 1 ) + 1
# pred_risk_A = 3*( pred_time - 1 ) - 2

xB = pred_time[ !is.na( pred_time ) ] - 1 + i
xA = pred_time[ !is.na( pred_time ) ] - 2 + i

yB = pred_risk[ !is.na( pred_time ) ]
#as.numeric( risk_pred_now[ cbind( which(!is.na( pred_time )), pred_risk_B[ which( !is.na(pred_risk_B) )] ) ] )
yA = rep( 1, length( xB ) )
yA[ which( xB != i ) ]  = pred_risk_before 
#as.numeric( risk_pred_now[ cbind( which(!is.na( pred_time ) & pred_risk_A > 0 ), pred_risk_A[ which( pred_risk_A > 0 ) ] ) ] )
m = (yB-yA)/(xB-xA)
real_pred_time = xB + (0.95-yB)/m

length(real_pred_time)
#risk_pred_delta_time_per_pt = i + pred_time[ which( !is.na( pred_time ) ) ] - 1 + as.numeric( risk_pred_now[ cbind( which( !is.na(pred_time) ) , 3*pred_time[!is.na(pred_time)] ) ] )
risk_pred_now$pred_time = pmin( risk_pred_now$end_corr, i + GLOB_horizon)
risk_pred_now$pred_time[ which( !is.na( pred_time ) ) ] = real_pred_time
pred_time[!is.na(pred_time)][1:10]
yA[1:10]
yB[1:10]
real_pred_time[1:10]

summary(risk_pred_now$`5.y.risk.1`)
# summary(risk_pred_now$`5.y.risk.1`[ which( risk_pred_now$`5.y.risk.1` != 0 ) ])
summary(risk_pred_now$`5.y.risk.2`[ which( risk_pred_now$`5.y.risk.2` != 0 ) ])
summary(risk_pred_now$`5.y.risk.3`[ which( risk_pred_now$`5.y.risk.3` != 0 ) ])
summary(risk_pred_now$`5.y.risk.4`[ which( risk_pred_now$`5.y.risk.4` != 0 ) ])
summary(risk_pred_now$`5.y.risk.5`[ which( risk_pred_now$`5.y.risk.5` != 0 ) ])
summary(risk_pred_now$`5.y.risk.6`[ which( risk_pred_now$`5.y.risk.6` != 0 ) ])
summary(risk_pred_now$`5.y.risk.7`[ which( risk_pred_now$`5.y.risk.7` != 0 ) ])
summary(risk_pred_now$`5.y.risk.8`[ which( risk_pred_now$`5.y.risk.8` != 0 ) ])
summary(risk_pred_now$`5.y.risk.9`[ which( risk_pred_now$`5.y.risk.9` != 0 ) ])
summary(risk_pred_now$`5.y.risk.10`[ which( risk_pred_now$`5.y.risk.10` != 0 ) ])


comparing_schemes <- data.frame(patid = risk_pred_now$patid,
                                risk_class = risk_pred_now$risk_class,
                                pred_time = risk_pred_now$pred_time,
                                pred_status = risk_pred_now$pred_status,
                                statin_age = risk_pred_now$statin_age,
                                cvd_age = risk_pred_now$cvd_age,
                                cvd_ind = risk_pred_now$cvd_ind,
                                death_age = risk_pred_now$death_age,
                                death_ind = risk_pred_now$death_ind,
                                end_age = risk_pred_now$end_age)




#comparing_schemes = cbind( patid = risk_pred_now$patid, risk_class = risk_pred_now$risk_class, 
#                           pred_time = risk_pred_now$pred_time, pred_status = risk_pred_now$pred_status,
#                           comparing_schemes )


save(comparing_schemes, file = paste0(out_path,"comparing_schemes_RIRS_apply_", i, "_", gender, ".RData") )
