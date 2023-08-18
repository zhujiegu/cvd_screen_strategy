library(data.table)
library(tidyverse)
library(gridExtra)
library(magrittr)
library(scales)
library(RColorBrewer)
library(nlme)
library(dplyr)
library(survival)

summarize <- dplyr::summarize
select <- dplyr::select

in_path2 <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/2_lmfit_outp/'
in_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/3_cox_outp/'

lm_age = seq(40,80,5)
# gender = "male" #female or male

##########################################################
# LMEM
##########################################################
beta_LMEM = NULL
#beta LMEM CVD
for( gender in c("female", "male")){
  for( i in lm_age ){
      load(paste0(in_path2, i,"_", gender, ".RData"))
      conf_int_temp = intervals(fit_rirs, which = "fixed")
      beta_LMEM_temp = data.frame( fitted_val = conf_int_temp$fixed[,2] )
      beta_LMEM_temp$lower  = conf_int_temp$fixed[,1]
      beta_LMEM_temp$upper  = conf_int_temp$fixed[,3]
      beta_LMEM_temp$lm_age = i
      beta_LMEM_temp$sex = gender
      beta_LMEM_temp$var = names( fit_rirs$coefficients$fixed )
      beta_LMEM_temp$type_exp = c("bmi", "hdl", "sbp", "smoke", "tchol", "bmi", "hdl", "sbp", "smoke", "tchol", "sbp", "tchol")
      beta_LMEM_temp$type_model = c( rep("intercept", 5), rep("slope", 5), "bp_bin", "statin_bin" )
      beta_LMEM  = rbind( beta_LMEM, beta_LMEM_temp )
      rm(beta_LMEM_temp)
  }
}

beta_LMEM %<>% mutate(sex=ifelse(sex=='male', 'Men', 'Women'))
beta_LMEM$sex %<>%  factor(levels = c('Women','Men'))

#Fixed Intercepts
int_BMI = ggplot( subset( beta_LMEM, var %in% c( "exposurebmi" ) ), aes( x = lm_age, y = fitted_val, col = sex)) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "", y = "BMI", title = "Intercepts", col = "") +
  lims( x = c(40,80), y = c(-0.85,0.6))+ #, y = c(-1.4,0.65)
  geom_hline( yintercept = 0, col = 1) +
  theme_bw() +
  theme( legend.position = "none" , plot.margin=margin(l = 0.1, b = -2.5, r = 0, t = 0, unit="cm") ) 

int_BMI_legend = ggplot( subset( beta_LMEM, var %in% c( "exposurebmi" ) ), aes( x = lm_age, y = fitted_val, col = sex)) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "", y = "BMI",  col = "") +
  lims( x = c(40,80), y = c(-0.85,0.6))+
  geom_hline( yintercept = 0, col = 1) +
  theme_bw() +
  theme(legend.text = element_text(size = 12)) 

int_HDL = ggplot( subset( beta_LMEM, var %in% c( "exposurehdl" ) ), aes( x = lm_age, y = fitted_val, col = sex)) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "", y = "HDL", col = "") +
  lims( x = c(40,80), y = c(-0.85,0.6))+
  geom_hline( yintercept = 0, col = 1) +
  theme_bw() +
  theme( legend.position = "none" , plot.margin=margin( l = 0.1, t = -2.5, b = -2.5, r = 0, unit="cm") ) 

int_smokebin = ggplot( subset( beta_LMEM, var %in% c( "exposuresmokbin" ) ), aes( x = lm_age, y = fitted_val, col = sex)) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "", y = "Smoke", col = "") +
  lims( x = c(40,80), y = c(-0.85,0.6))+
  geom_hline( yintercept = 0, col = 1) +
  theme_bw() +
  theme( legend.position = "none", plot.margin=margin( l = 0.1, t = -2.5, b = -2.5, r = 0, unit="cm") ) 

int_sbp = ggplot( subset( beta_LMEM, var %in% c( "exposuresbp" ) ), aes( x = lm_age, y = fitted_val, col = sex)) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "", y = "SBP",  col = "") +
  lims( x = c(40,80), y = c(-0.85,0.6))+
  geom_hline( yintercept = 0, col = 1) +
  theme_bw() +
  theme( legend.position = "none" , plot.margin=margin( l = 0.1, t = -2.5, b = -2.5, r = 0, unit="cm") ) 

int_tchol = ggplot( subset( beta_LMEM, var %in% c( "exposuretchol" ) ), aes( x = lm_age, y = fitted_val, col = sex)) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "Landmark age", y = "Tot. Chol.", col = "") +
  lims( x = c(40,80), y = c(-0.85,0.6))+
  geom_hline( yintercept = 0, col = 1) +
  theme_bw() +
  theme( legend.position = "none" , plot.margin=margin( l = 0.1, t = -2.5, b = -2.5, r = 0, unit="cm") ) 

#Fixed slopes
slope_BMI = ggplot( subset( beta_LMEM, var %in% c( "exposurebmi:exp_age_corr" ) ), aes( x = lm_age, y = fitted_val, col = sex)) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "", y = "", title = "Age effect", col = "") +
  lims( x = c(40,80), y = c(-0.06,0.05))+
  geom_hline( yintercept = 0, col = 1) +
  theme_bw() +
  theme( legend.position = "none", plot.margin=margin( l = -0.3,  b = -2.5, t= 0, r = 0, unit="cm")  ) 

slope_HDL = ggplot( subset( beta_LMEM, var %in% c( "exposurehdl:exp_age_corr" ) ), aes( x = lm_age, y = fitted_val, col = sex)) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "", y = "",  col = "") +
  lims( x = c(40,80), y = c(-0.06,0.05))+
  geom_hline( yintercept = 0, col = 1) +
  theme_bw() +
  theme( legend.position = "none" , plot.margin=margin(l=-0.3, t = -2.5, b = -2.5,r = 0, unit="cm")  ) 

slope_smokebin = ggplot( subset( beta_LMEM, var %in% c( "exposuresmokbin:exp_age_corr" ) ), aes( x = lm_age, y = fitted_val, col = sex)) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "", y = "", col = "") +
  lims( x = c(40,80), y = c(-0.06,0.05))+
  geom_hline( yintercept = 0, col = 1) +
  theme_bw() +
  theme( legend.position = "none", plot.margin=margin(l=-0.3, t = -2.5, b = -2.5, r = 0, unit="cm")  ) 

slope_sbp = ggplot( subset( beta_LMEM, var %in% c( "exposuresbp:exp_age_corr" ) ), aes( x = lm_age, y = fitted_val, col = sex)) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "", y = "",  col = "") +
  lims( x = c(40,80), y = c(-0.06,0.05))+
  geom_hline( yintercept = 0, col = 1) +
  theme_bw() +
  theme( legend.position = "none", plot.margin=margin(l=-0.3, t = -2.5, b = -2.5, r = 0, unit="cm") ) 

slope_tchol = ggplot( subset( beta_LMEM, var %in% c( "exposuretchol:exp_age_corr" ) ), aes( x = lm_age, y = fitted_val, col = sex)) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "Landmark age", y = "", col = "") +
  lims( x = c(40,80), y = c(-0.06,0.05))+
  geom_hline( yintercept = 0, col = 1) +
  theme_bw() +
  theme( legend.position = "none", plot.margin=margin(l=-0.3, t = -2.5, b = -2.5, r = 0, unit="cm") ) 

#Adjusting the intercepts
bp_bin_LMEM = ggplot( subset( beta_LMEM, var %in% c( "bp_bin:sbp_ind") ), aes( x = lm_age, y = fitted_val, col = sex)) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "", y = "", title = "BP med. effect") +
  lims( x = c(40,80))+
  geom_hline( yintercept = 0, col = 1) +
  theme_bw() +
  theme( legend.position = "none", plot.margin=margin(l=-0.3, t= 0, r = 0, b = 0, unit="cm") ) 

statin_LMEM = ggplot( subset( beta_LMEM, var %in% c( "statin_bin:tchol_ind") ), aes( x = lm_age, y = fitted_val, col = sex)) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "Landmark age", y = "", title = "Statin effect") +
  lims( x = c(40,80))+
  geom_hline( yintercept = 0, col = 1) +
  theme_bw() +
  theme( legend.position = "none", plot.margin=margin(l=-0.3, t= 0, r = 0, b = 0, unit="cm")) 


shared_legend = cowplot::get_legend( int_BMI_legend )



all_beta_LMEM = cowplot::plot_grid(int_BMI, slope_BMI, NULL,
                                   int_HDL, slope_HDL, NULL,
                   int_smokebin, slope_smokebin, NULL,
                   int_sbp, slope_sbp, bp_bin_LMEM,
                   int_tchol, slope_tchol, statin_LMEM,
                   align = 'hv', ncol = 3, nrow = 5, rel_heights = c(1,1,1,1,1))

final_plot_beta_LMEM = all_beta_LMEM + cowplot::draw_grob(shared_legend, 2/3.3, 0.2, 1.5/3.3, 1)


ggsave(final_plot_beta_LMEM, 
       file = paste0( "~/epi_paper/outp_figs/beta_LMEM_all.jpeg" ), width = 30, height = 40, unit = "cm",
       bg="white")


##########################################################
# beta 5 year CVD
##########################################################
beta_5_cox = NULL 
for( gender in c("female", "male")){
  for( i in lm_age ){
    for( j in 0:8){
      load(paste0(in_path, "cox_model_RIRS_", i,"_", j, "_", gender, ".RData"))
      beta_5_cox_temp = data.frame( exp_var = exp(model_cox_5CVD$coefficients) )
      beta_5_cox_temp$lower  = exp( confint(model_cox_5CVD)[,1] )
      beta_5_cox_temp$upper  = exp( confint(model_cox_5CVD)[,2] )
      beta_5_cox_temp$lm_age_j = i + j
      beta_5_cox_temp$lm_age = i
      beta_5_cox_temp$sex = gender
      beta_5_cox = rbind( beta_5_cox, beta_5_cox_temp )
      rm(beta_5_cox_temp)
    }
  }
}

beta_5_cox %<>% mutate(sex=ifelse(sex=='male', 'Men', 'Women'))
beta_5_cox$sex %<>%  factor(levels = c('Women','Men'))

#TIME DEPENDENT VARIABLES

#HDL
hdl = ggplot( beta_5_cox[ grep("hdl", rownames(beta_5_cox) ), ], aes(x = lm_age_j, y = exp_var, col = factor(lm_age)) ) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "Landmark age", y = "HR - HDL", col = "Risk class") +
  lims( x = c(40,88))+
  geom_hline( yintercept = 1, col = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(rows = vars(sex))
#ggsave(hdl, file = paste0( "img_epi_paper/beta_cox_hdl.jpeg" ), width = 15, height = 15, unit = "cm"  )

#sbp
sbp = ggplot( beta_5_cox[ grep("sbp", rownames(beta_5_cox) ), ], aes(x = lm_age_j, y = exp_var, col = factor(lm_age)) ) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "Landmark age", y = "HR - SBP", col = "Risk class") +
  lims( x = c(40,88))+
  geom_hline( yintercept = 1, col = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(rows = vars(sex))
#ggsave(sbp, file = paste0( "img_epi_paper/beta_cox_SBP.jpeg" ), width = 15, height = 15, unit = "cm"  )

#smoke
smoke = ggplot( beta_5_cox[ grep("smoke", rownames(beta_5_cox) ), ], aes(x = lm_age_j, y = exp_var, col = factor(lm_age)) ) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "Landmark age", y = "HR - Smoke", col = "Risk class") +
  lims( x = c(40,88))+
  geom_hline( yintercept = 1, col = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(rows = vars(sex))
#ggsave(smoke, file = paste0( "img_epi_paper/beta_cox_smoke.jpeg" ), width = 15, height = 15, unit = "cm"  )

#tchol
tchol = ggplot( beta_5_cox[ grep("tchol", rownames(beta_5_cox) ), ], aes(x = lm_age_j, y = exp_var, col = factor(lm_age)) ) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "Landmark age", y = "HR - Tot. Chol.", col = "Risk class") +
  lims( x = c(40,88))+
  geom_hline( yintercept = 1, col = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(rows = vars(sex))
#ggsave(tchol, file = paste0( "img_epi_paper/beta_cox_tchol.jpeg" ), width = 15, height = 15, unit = "cm"  )

#bmi
bmi = ggplot( beta_5_cox[ grep("bmi", rownames(beta_5_cox) ), ], aes(x = lm_age_j, y = exp_var, col = factor(lm_age)) ) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "Landmark age", y = "HR - BMI", col = "Risk class") +
  lims( x = c(40,88))+
  geom_hline( yintercept = 1, col = 1) +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank(), legend.text = element_text(size = 12)) +
  facet_grid(rows = vars(sex))
#ggsave(bmi, file = paste0( "img_epi_paper/beta_cox_bmi.jpeg" ), width = 15, height = 15, unit = "cm"  )

shared_legend_time_dep = cowplot::get_legend( bmi )

all_beta_cox_time_dep = cowplot::plot_grid(bmi + theme(legend.position = "none"),
                                           hdl, sbp, smoke, tchol, NULL, align = 'hv',  labels=c('A', 'B',
                                                                                               'C', 'D',
                                                                                               'E',  ''), ncol = 2, nrow = 3)

final_plot_beta_cox_time_dep = all_beta_cox_time_dep + cowplot::draw_grob(shared_legend_time_dep, 2/3.3, -0.3, 1.0/3.3, 1)


ggsave(final_plot_beta_cox_time_dep, file = paste0( "~/epi_paper/outp_figs/beta_cox_time_dep.jpeg" ), width = 30, height = 40, unit = "cm" ,
       bg="white")

#TIME FIXED covariates
#townsend
townsend2 = ggplot( beta_5_cox[ grep("Townsend2", rownames(beta_5_cox) ), ], aes(x = lm_age_j, y = exp_var, col = factor(lm_age)) ) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "Landmark age", y = "HR - Townsend2", col = "Risk class") +
  lims( x = c(40,90))+
  geom_hline( yintercept = 1, col = 1) +
  theme_bw() +
  theme(legend.position = "top", legend.box = "horizontal", legend.title = element_blank(), legend.text = element_text(size = 12)) +
  guides(color = guide_legend(nrow = 1))+
  facet_grid(rows = vars(sex))
townsend3 = ggplot( beta_5_cox[ grep("Townsend3", rownames(beta_5_cox) ), ], aes(x = lm_age_j, y = exp_var, col = factor(lm_age)) ) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "Landmark age", y = "HR - Townsend3", col = "Risk class") +
  lims( x = c(40,90))+
  geom_hline( yintercept = 1, col = 1) +
  theme_bw() +
  theme(legend.position = "top", legend.box = "horizontal", legend.title = element_blank(), legend.text = element_text(size = 12)) +
  guides(color = guide_legend(nrow = 1))+
  facet_grid(rows = vars(sex))
townsend4 = ggplot( beta_5_cox[ grep("Townsend4", rownames(beta_5_cox) ), ], aes(x = lm_age_j, y = exp_var, col = factor(lm_age)) ) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "Landmark age", y = "HR - Townsend4", col = "Risk class") +
  lims( x = c(40,90))+
  geom_hline( yintercept = 1, col = 1) +
  theme_bw() +
  theme(legend.position = "top", legend.box = "horizontal", legend.title = element_blank(), legend.text = element_text(size = 12)) +
  guides(color = guide_legend(nrow = 1))+
  facet_grid(rows = vars(sex))
townsend5 = ggplot( beta_5_cox[ grep("Townsend5", rownames(beta_5_cox) ), ], aes(x = lm_age_j, y = exp_var, col = factor(lm_age)) ) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "Landmark age", y = "HR - Townsend5", col = "Risk class") +
  lims( x = c(40,90))+
  geom_hline( yintercept = 1, col = 1) +
  theme_bw() +
  theme(legend.position = "top", legend.box = "horizontal", legend.title = element_blank(), legend.text = element_text(size = 12)) +
  guides(color = guide_legend(nrow = 1))+
  facet_grid(rows = vars(sex))
#ggsave(townsend, file = paste0( "img_epi_paper/beta_cox_townsend.jpeg" ), width = 15, height = 15, unit = "cm"  )

#blood pressure medication
bp_bin = ggplot( beta_5_cox[ grep("bp_bin", rownames(beta_5_cox) ), ], aes(x = lm_age_j, y = exp_var, col = factor(lm_age)) ) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "Landmark age", y = "HR - BP medication", col = "Risk class") +
  lims( x = c(40,90))+
  geom_hline( yintercept = 1, col = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(rows = vars(sex))
#ggsave(bp_bin, file = paste0( "img_epi_paper/beta_cox_bp_bin.jpeg" ), width = 15, height = 15, unit = "cm"  )

#diabetes
diab = ggplot( beta_5_cox[ grep("diab_ind", rownames(beta_5_cox) ), ], aes(x = lm_age_j, y = exp_var, col = factor(lm_age)) ) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "Landmark age", y = "HR - Diabetes", col = "Risk class") +
  lims( x = c(40,90))+
  geom_hline( yintercept = 1, col = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(rows = vars(sex))
#ggsave(diab, file = paste0( "img_epi_paper/beta_cox_diab.jpeg" ), width = 15, height = 15, unit = "cm"  )


#depression
depression = ggplot( beta_5_cox[ grep("depression", rownames(beta_5_cox) ), ], aes(x = lm_age_j, y = exp_var, col = factor(lm_age)) ) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "Landmark age", y = "HR - Depression", col = "Risk class") +
  lims( x = c(40,90))+
  geom_hline( yintercept = 1, col = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(rows = vars(sex))
#ggsave(depression, file = paste0( "img_epi_paper/beta_cox_depression.jpeg" ), width = 15, height = 15, unit = "cm"  )

#Severe_mental
Severe_mental = ggplot( beta_5_cox[ grep("Severe_mental", rownames(beta_5_cox) ), ], aes(x = lm_age_j, y = exp_var, col = factor(lm_age)) ) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "Landmark age", y = "HR - Sev. Mental ill", col = "Risk class") +
  lims( x = c(40,90))+
  geom_hline( yintercept = 1, col = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(rows = vars(sex))
#ggsave(Severe_mental, file = paste0( "img_epi_paper/beta_cox_Severe_mental.jpeg" ), width = 15, height = 15, unit = "cm"  )

#Migraine
Migraine = ggplot( beta_5_cox[ grep("Migraine", rownames(beta_5_cox) ), ], aes(x = lm_age_j, y = exp_var, col = factor(lm_age)) ) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "Landmark age", y = "HR - Migraine", col = "Risk class") +
  lims( x = c(40,90))+
  geom_hline( yintercept = 1, col = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(rows = vars(sex))
#ggsave(Migraine, file = paste0( "img_epi_paper/beta_cox_Migraine.jpeg" ), width = 15, height = 15, unit = "cm"  )

shared_legend_time_fixed1 = cowplot::get_legend( townsend2 )

all_beta_cox_time_fixed1 = cowplot::plot_grid(townsend2 + theme(legend.position = "none"),
                                              townsend3 + theme(legend.position = "none"),
                                              townsend4 + theme(legend.position = "none"),
                                              townsend5 + theme(legend.position = "none"),
                                              bp_bin, 
                                              diab, 
                                              depression, 
                                              Migraine, 
                                              Severe_mental, 
                                              align = 'h',  
                                              labels=c('A', 'B','C', 'D', 'E',  'F','G','H','I'), ncol = 2, nrow = 6, 
                                              rel_heights = c(1,1,1,1,1,0.2))

final_plot_beta_cox_time_fixed1 = all_beta_cox_time_fixed1 + cowplot::draw_grob(shared_legend_time_fixed1, 1.2/3.3, -0.48, 1.0/3.3, 1)


ggsave(final_plot_beta_cox_time_fixed1, 
       file = paste0("~/epi_paper/outp_figs/beta_cox_time_fixed_1.jpeg" ), width = 30, height = 40, unit = "cm" ,
       bg="white")



#Time fixed variables part 2

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

tot_col = gg_color_hue(9)

#atrial_f
atrial_f = ggplot( beta_5_cox[ grep("atrial_f", rownames(beta_5_cox) ), ], aes(x = lm_age_j, y = exp_var, col = factor(lm_age)) ) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "Landmark age", y = "HR - Atrial fibrillation", col = "Risk class") +
  lims( x = c(40,90))+
  geom_hline( yintercept = 1, col = 1) +
  scale_color_manual( values = tot_col[5:9]) +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank(), legend.text = element_text(size = 12)) +
  facet_grid(rows = vars(sex))
#ggsave(atrial_f, file = paste0( "img_epi_paper/beta_cox_atrial_f.jpeg" ), width = 15, height = 15, unit = "cm"  )

#renal
renal = ggplot( beta_5_cox[ grep("renal", rownames(beta_5_cox) ), ], aes(x = lm_age_j, y = exp_var, col = factor(lm_age)) ) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "Landmark age", y = "HR - Renal disease", col = "Risk class") +
  lims( x = c(40,90))+
  geom_hline( yintercept = 1, col = 1) +
  scale_color_manual( values = tot_col[5:9]) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(rows = vars(sex))
#ggsave(renal, file = paste0( "img_epi_paper/beta_cox_renal.jpeg" ), width = 15, height = 15, unit = "cm"  )

#rheumatoid
rheumatoid = ggplot( beta_5_cox[ grep("rheumatoid", rownames(beta_5_cox) ), ], aes(x = lm_age_j, y = exp_var, col = factor(lm_age)) ) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "Landmark age", y = "HR - Rheumatoid Arthritis", col = "Risk class") +
  lims( x = c(40,90))+
  geom_hline( yintercept = 1, col = 1) +
  scale_color_manual( values = tot_col[5:9]) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_grid(rows = vars(sex))
#ggsave(rheumatoid, file = paste0( "img_epi_paper/beta_cox_rheumatoid.jpeg" ), width = 15, height = 15, unit = "cm"  )


shared_legend_time_fixed2 = cowplot::get_legend( atrial_f )

all_beta_cox_time_fixed2 = cowplot::plot_grid( atrial_f + theme(legend.position = "none"), rheumatoid, renal, NULL, align = 'h',  labels=c('A', 'B', 'C', ''), ncol = 2, nrow = 2)

final_plot_beta_cox_time_fixed2 = all_beta_cox_time_fixed2 + cowplot::draw_grob(shared_legend_time_fixed2, 2/3.3, -0.3, 1.0/3.3, 1)


ggsave(final_plot_beta_cox_time_fixed2, file = paste0( "~/epi_paper/outp_figs/beta_cox_time_fixed_2.jpeg" ), width = 30, height = 40, unit = "cm" ,
       bg="white")



# 
# #plot risk
# risk_tot = NULL
# 
# for(j in c("female", "male"))
# {
#   for( i in lm_age )
#   {
#     load(paste0("results/results_", j, "/risk_pred_FE4_RIRS_",i,"_",j,".RData"))
#     risk_temp = risk_prediction[ !is.na( risk_prediction$`5.y.risk.0` ), ]
#     risk_temp$lm_age = i
#     risk_temp$sex = j
#     risk_tot = rbind( risk_tot, risk_temp )
#     rm(risk_temp)
#   }
# }
# 
# risk_tot[lm_age == 40 & `5.y.risk.0` <= 0.97, ]
# risk_tot[patid =="71091229", ]
# 
# 
# risk_to_plot = risk_tot %>% 
#   filter( patid %in% c("15187","71091229", "254205", "79435") & lm_age == 40 ) %>%
#   select( patid, sex, lm_age, cvd_age, death_age, end_age, paste0( "5.y.risk.", 0:8)) %>%
#   pivot_longer( cols = contains( "5.y.risk"), names_to = "risk_date", values_to = "value")
#  
# risk_to_plot$pred_age = risk_to_plot$lm_age + rep(0:8,4) 
# risk_to_plot$event_type = ifelse( !is.na(risk_to_plot$cvd_age), "CVD",
#                                   ifelse( !is.na( risk_to_plot$death_age ), "Death", "Censored" ) )
# 
# risk_to_plot$first_event_time = apply( risk_to_plot[, c("cvd_age", "death_age", "end_age")],1, function(x) min(x, na.rm = T)) 
#   
# risk_traj = ggplot( risk_to_plot, aes( x = pred_age, y = 100*(1-value), col = factor(patid))) +
#   geom_point() +
#   geom_hline(yintercept = 5) +
#   geom_vline(aes( xintercept = first_event_time, col = factor(patid), linetype = event_type)) +
#   labs(x = "Prediction time", y = "5-year CVD risk [%]") +
#   theme_bw()
# 
# ggsave(risk_traj, file = paste0( "img_epi_paper/risk_trajectory_40.jpeg" ), width = 15, height = 15, unit = "cm"  )
# 
# 
# #EPI paper: plotting scheme for males 60
# gender = "male"
# i = 60
# load(paste0("results/results_", gender, "/fit_LMEM_rirs_", i,"_", gender, ".RData"))
# load(paste0("data_created/", gender, "/data_lm_", i,"_", gender, ".RData"))
# 
# data_lines = data_ext[patid %in% c("1027", "1329", "1349", "2018"), c("patid", "exp_date","exp_age", "exposure", "bp_bin", "statin_bin")]
# data_lines$estim_LMEM = 0
# 
# #BMI 
# bmi_idx = which(data_lines$exposure == "bmi")
# data_lines$estim_LMEM[ bmi_idx ] = fit_rirs$coefficients$fixed[1] + fit_rirs$coefficients$fixed[6] * data_lines$exp_age[ bmi_idx ]
# #HDL 
# hdl_idx = which(data_lines$exposure == "hdl")
# data_lines$estim_LMEM[ hdl_idx ] = fit_rirs$coefficients$fixed[2] + fit_rirs$coefficients$fixed[7] * data_lines$exp_age[ hdl_idx ]
# #SBP 
# sbp_idx = which(data_lines$exposure == "sbp")
# data_lines$estim_LMEM[ sbp_idx ] = fit_rirs$coefficients$fixed[3] + fit_rirs$coefficients$fixed[8] * data_lines$exp_age[ sbp_idx ] + fit_rirs$coefficients$fixed[11] * data_lines$bp_bin[ sbp_idx ]
# #smoke 
# smoke_idx = which(data_lines$exposure == "smoke")
# data_lines$estim_LMEM[ smoke_idx ] = fit_rirs$coefficients$fixed[4] + fit_rirs$coefficients$fixed[9] * data_lines$exp_age[ smoke_idx ] 
# #TCHOL 
# tchol_idx = which(data_lines$exposure == "tchol")
# data_lines$estim_LMEM[ tchol_idx ] = fit_rirs$coefficients$fixed[5] + fit_rirs$coefficients$fixed[10] * data_lines$exp_age[ tchol_idx ] + fit_rirs$coefficients$fixed[12] * data_lines$statin_bin[ tchol_idx ] 
# 
# unique(data_ext$patid)
# ggplot( subset( data_ext, patid %in% c("1027", "1329", "1349", "2018")), aes(x = exp_date, y = scaled_corr, col = factor(patid) ) ) +
#   geom_point() +
#   labs(x = "Exposure date", y = "Scaled values") +
#   geom_line(data = data_lines, aes(x = exp_date, y = estim_LMEM, col = as.factor(patid))) +
#   theme_bw() +
#   facet_grid(rows = vars(exposure), scales = "free_y") 
#   
# 
# #Plot
# plot(0:8, fit_rirs$coefficients$fixed[1] + fit_rirs$coefficients$fixed[6] * 0:8, type= "l" )
# lines(0:8, fit_rirs$coefficients$fixed[1] + fit_rirs$coefficients$fixed[6] * 0:8  +
#         fit_rirs$coefficients$random$patid[1,1] + fit_rirs$coefficients$random$patid[1,6] * 0:8, col = 2 )
# 
# 
# 
# 
# data_female = NULL
# data_male = NULL
# 
# 
# 
# for( i in lm_age )
# {
#   load(paste0("data_created/female/data_lm_", i,"_female.RData"))
#   data_ext$lm_age = i
#   data_ext$sex = "female"
#   data_female = rbind( data_female, data_ext )
#   rm(data_ext)
#   
#   load(paste0("data_created/male/data_lm_", i,"_male.RData"))
#   data_ext$lm_age = i
#   data_ext$sex = "male"
#   data_male = rbind( data_male, data_ext )
#   rm(data_ext)
#   
#   print(i)
# }
# 
# #Table
# 
# max(data_female$end_age)
# summary(data_female$end_age[which(data_female$lm_age == 80)] - data_female$lm_age[ which(data_female$lm_age == 80)])
# mean(data_female$end_age[which(data_female$lm_age == 80)] - data_female$lm_age[ which(data_female$lm_age == 80)])
# sd(data_female$end_age[which(data_female$lm_age == 80)] - data_female$lm_age[ which(data_female$lm_age == 80)])
# summary(data_female$end_age[which(data_female$lm_age == 75)] - data_female$lm_age[ which(data_female$lm_age == 75)])
# mean(data_female$end_age[which(data_female$lm_age == 75)] - data_female$lm_age[ which(data_female$lm_age == 75)])
# sd(data_female$end_age[which(data_female$lm_age == 75)] - data_female$lm_age[ which(data_female$lm_age == 75)])
# 
# age = 80
# idx_old_old = which( data_female$end_age[which(data_female$lm_age == age)] - data_female$lm_age[ which(data_female$lm_age == age)] >13.3)
# data_female[which(data_female$lm_age == age)[idx_old_old[10]],]
# 
# max(data_male$end_age)
# summary(data_male$end_age[which(data_male$lm_age == 80)] - data_male$lm_age[ which(data_male$lm_age == 80)])
# mean(data_male$end_age[which(data_male$lm_age == 80)] - data_male$lm_age[ which(data_male$lm_age == 80)])
# sd(data_male$end_age[which(data_male$lm_age == 80)] - data_male$lm_age[ which(data_male$lm_age == 80)])
# 
# summary(data_male$end_age[which(data_male$lm_age == 75)] - data_male$lm_age[ which(data_male$lm_age == 75)])
# mean(data_male$end_age[which(data_male$lm_age == 75)] - data_male$lm_age[ which(data_male$lm_age == 75)])
# sd(data_male$end_age[which(data_male$lm_age == 75)] - data_male$lm_age[ which(data_male$lm_age == 75)])
# 
# age = 80
# idx_old_old = which( data_male$end_age[which(data_male$lm_age == age)] - data_male$lm_age[ which(data_male$lm_age == age)] >13.3)
# data_male[which(data_male$lm_age == age)[idx_old_old[10]],]
# 
# #load("../../patients_outcomes.RData")
# #load("exposures_red_merged.RData")
# #family_data = read_dta("../familyhistory.dta")
# 
# 
# #table( outcomes$pracid[outcomes$derivation == "derivation"] )
