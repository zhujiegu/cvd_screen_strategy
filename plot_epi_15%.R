library(data.table)
library(nlme)
library(survival)
library(tidyverse)
library(haven)
library(miceadds)
library(gridExtra)
library(magrittr)
library(RColorBrewer)
summarize <- dplyr::summarize
select <- dplyr::select

in_path2 <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/2_lmfit_outp/'
in_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/3_cox_outp/'
in_path3 <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/4_screening_outp/'
out_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/3_cox_outp/'

lm_age = seq(40,80,5)
gender = "male" #female or male
deriv_tmp = "derivation" #derivation or validation
cov = NULL
outcome = NULL
tstar = NULL
beta_5_cox = NULL 
beta_LMEM = NULL

#Covariates matrix
for( i in lm_age ){
  load(paste0(in_path, "/covariates_matrix_RIRS_", i,"_", gender, ".RData"))
  covariates_matrix = covariates_matrix[ derivation == deriv_tmp, 1:21]
  covariates_matrix$lm_age = i
  covariates_matrix$sex = gender
  cov = rbind( cov, covariates_matrix )
  rm(covariates_matrix)
}

binary_variables = c("diab_ind", "bp_med", "renal_disease", "atrial_fibrillation", "rheumatoid_arthritis",
                     "Severe_mental_illness_ind", "Migraine_ind", "Dementia_ind", "depression_ind", "SLE_ind", "family_history")

continuous_variables = c("townsend_20", "bmi_blup_0", "hdl_blup_0", "sbp_blup_0", "smoke_blup_0", "tchol_blup_0" )

cov = data.table(cov)

cont_var_summary = function(x){ paste0(  round(mean(x),2), " (", round(sd(x),2), ")" ) }

cont_desc = cov[, lapply(.SD, cont_var_summary), by=lm_age, .SDcols = continuous_variables ]

binary_var_summary = function(x){ paste0( length(which(x == 1)), " (", round( length(which(x == 1))/length(x)*100,2), "%)"  ) }

binary_desc = cov[,lapply(.SD, binary_var_summary), by=lm_age, .SDcols = binary_variables ]

final_descriptive = left_join( cont_desc, binary_desc, by = c( "lm_age" ) )

# write.table( t(final_descriptive), file = paste0( "img_epi_paper/final_descriptive_", gender, "_covariates_", deriv_tmp, ".txt"), quote = F, col.names = F, sep = ",")

#Outcome description
for( i in lm_age ){
  load(paste0(in_path, "risk_pred_FE4_RIRS_", i,"_", gender, ".RData"))
  outcome_matrix = risk_prediction[ !is.na(`5.y.risk.0`), 1:8 ]
  outcome_matrix$lm_age = i
  outcome_matrix$sex = gender
  outcome = rbind( outcome, outcome_matrix )
  rm(outcome_matrix)
}

outcome
outcome[ , fatal_cvd := ifelse( cvd_age == death_age & cvd_ind == 1 & death_ind == 1, 1, 0 ) ]
outcome[ , non_fatal_cvd := ifelse( (cvd_age < death_age | is.na(death_age)) & cvd_ind == 1, 1, 0 ) ]
outcome[ , death_other_causes := ifelse( (cvd_age < death_age | is.na(cvd_age)) & death_ind == 1, 1, 0 ) ]
outcome[ , followup := end_age - lm_age  ]

outcome_binary_desc = outcome[,lapply(.SD, binary_var_summary), by=lm_age, .SDcols = c("fatal_cvd", "non_fatal_cvd", "death_other_causes") ] #"cvd_ind", "death_ind", 
outcome_binary_desc$total = outcome[,lapply(.SD, length), by=lm_age, .SDcols = "patid"]$patid

outcome_cont_desc = outcome[, lapply(.SD, cont_var_summary), by=lm_age, .SDcols = "followup" ]


final_descriptive_outcome = left_join( outcome_cont_desc, outcome_binary_desc, by = "lm_age" )



#t*
for(j in c("female", "male")){
  for( i in lm_age ){
    load(paste0(in_path3, "comparing_schemes_RIRS_apply_", i,"_", j, "_15%.RData"))
    tstar_matrix = comparing_schemes[ , c('patid',"risk_class", "pred_time", "pred_status") ]
    tstar_matrix$lm_age = i
    tstar_matrix$sex = j
    tstar = rbind( tstar, tstar_matrix )
    rm(tstar_matrix)
  }
}
get_limits_mean_pred_time = function(x){ list(mean  = mean(x),
                                              lower = quantile( x, 0.25 ), 
                                              upper = quantile( x, 0.75 ),
                                              size = length(x)) }

tstar = data.table(tstar)
tstar$delta = tstar$pred_time-tstar$lm_age

# Proportion that did not cross the threshold in 8 years
tstar %<>% mutate(pred_status = replace(pred_status, delta > 8, 0))
levels(tstar$risk_class) <- c('Very_high','High','Medium_high',
                              'Medium_low','Low')
# levels(tstar$risk_class) <- c('Very_high(>5%)','High(3.75%-5%)','Medium_high(2.5%-3.75%)',
#                                        'Medium_low(1.25%-2.5%)','Low(<1.25%)')

# Total: n per (lm_age,sex,risk,crossed). 
# N: n per (lm_age, sex)
# riskgr_pct_total: % per (lm_age, sex, risk)/ (lm_age, sex)
# cross_pct: % per (lm_age, sex, risk,crossed)/ (lm_age, sex, risk)
N_cross_riskgr <- tstar %>% 
  group_by(lm_age, sex, risk_class, pred_status) %>% summarize(Total = n()) %>% 
  group_by(lm_age, sex) %>% mutate(N = sum(Total), riskgr_pct = Total/sum(Total)) %>% 
  group_by(lm_age, sex, risk_class) %>% mutate(riskgr_pct_total = sum(riskgr_pct), 
                                               cross_pct = Total/sum(Total))

N_cross_riskgr$pred_status <- N_cross_riskgr$pred_status==1
colnames(N_cross_riskgr)[4] <- 'crossed'

N_cross_riskgr %>% filter(risk_class=='Very_high') %>% filter(sex=='male')



# N_cross_riskgr <- tstar %>% group_by(lm_age, sex, risk_class) %>% summarize(Total = length(pred_status), 
#                                                          N_cross = sum(pred_status),
#                                                          Mean_cross = mean(pred_status),
#                                                          Mean_nocross = mean(1-pred_status),
#                                                          check = Mean_cross + Mean_nocross)
# N_cross_riskgr <- N_cross_riskgr %>% group_by(lm_age, sex) %>% 
#   mutate(riskgr_pct = Total/sum(Total)) %>% ungroup
##########################################################
## Major updated by ZG Dec 2022
##########################################################
# # Same subject appears only in one lm cohort (randomly selected)
# # This will result in biased results
# tstar_uniqueid <- tstar %>%  group_by(patid) %>% sample_n(1) %>% ungroup
# nocross_uniqueid <- tstar_uniqueid %>% group_by(lm_age, sex) %>% summarize(Total = length(pred_status), 
#                                                          N_nocross = sum(1-pred_status), Mean = 1- mean(pred_status))
# nocross_uniqueid
# 
# tstar_uniqueid_plot = tstar_uniqueid[ pred_status == 1 & delta <=8, 
#                     unlist( recursive = F, lapply( .SD, get_limits_mean_pred_time )), 
#                     .SDcols = "delta", by = c("lm_age", "risk_class", "sex") ]
# 
# tstar_uniqueid_distr = ggplot( tstar_uniqueid_plot, aes(x = lm_age, y = delta.mean, col = risk_class)) +
#   geom_pointrange(aes( ymin = delta.lower, ymax = delta.upper), position = position_dodge(width = 3)) +
#   labs( x = "Landmark age (unique ID)", y = "Time between landmark age and t*", col = "") +
#   geom_hline(yintercept = 5) +
#   facet_grid(~sex) +
#   theme_bw() +
#   theme(legend.position = "top", legend.justification = "center", legend.margin = margin(0,0,0,0), legend.box.margin = margin(-5,-5,-5,-5) ) 
# 
# ggsave(tstar_uniqueid_distr, file = paste0( "~/epi_paper/outp_figs/tstar_uniqueid_distr.jpeg" ), width = 25, height = 15, unit = "cm"  )

##########################################################
## Old tstar plot
# tstar_plot = tstar[ pred_status == 1 & delta <=8, 
#                     unlist( recursive = F, lapply( .SD, get_limits_mean_pred_time )), 
#                     .SDcols = "delta", by = c("lm_age", "risk_class", "sex") ]
# 
# tstar_distr = ggplot( tstar_plot, aes(x = lm_age, y = delta.mean, col = risk_class)) +
#   geom_pointrange(aes( ymin = delta.lower, ymax = delta.upper), position = position_dodge(width = 3)) +
#   labs( x = "Landmark age", y = "Time between landmark age and t*", col = "") +
#   geom_hline(yintercept = 5) +
#   facet_grid(~sex) +
#   theme_bw() +
#   theme(legend.position = "top", legend.justification = "center", legend.margin = margin(0,0,0,0), legend.box.margin = margin(-5,-5,-5,-5) ) 
# 
# ggsave(tstar_distr, file = paste0( "~/epi_paper/outp_figs/tstar_distr.jpeg" ), width = 25, height = 15, unit = "cm"  )

## New tstar plot
tstar_crossed = tstar %>% filter(pred_status==1) %>% filter(risk_class != 'Very_high') %>% 
  group_by(lm_age, sex, risk_class) %>% mutate(mean = mean(delta)) %>% ungroup
# ggplot(tstar_crossed, aes(x=delta, color=risk_class, fill = risk_class)) +
#   geom_density(alpha = 0.2)

tstar_all <- tstar %>% filter(risk_class != 'Very_high') %>% group_by(lm_age, sex, risk_class) %>% 
  mutate(delta = ifelse(pred_status==1, delta, 8)) %>% mutate(mean = mean(delta)) %>% ungroup

N_tstar_all <- tstar_all %>% 
  group_by(lm_age, sex, risk_class) %>% summarize(Total = n())

risk_gr_colors <- c("red4", "red","darkorange", "royalblue3", "springgreen3")
names(risk_gr_colors)<- N_cross_riskgr$risk_class %>% levels

# keep empty facet rows to be comparable
tstar_all$lm_age %<>% as.factor()
N_tstar_all$lm_age %<>% as.factor()
tstar_all$risk_class <- droplevels(tstar_all$risk_class)
N_tstar_all$risk_class <- droplevels(N_tstar_all$risk_class)

p_tstar_plot <- function(gender, risk_labels){
  ggplot(tstar_all %>% filter(sex==gender), aes(x=delta, color=risk_class, fill = risk_class)) + 
    geom_histogram(alpha = 0.2,  position = 'identity', breaks=seq(0,9,by=1), closed='left', aes(y = ..density..)) +
    geom_vline(aes(xintercept=mean, color='black'), linetype = "dashed") +
    geom_vline(aes(xintercept=8), linetype = "solid") +
    # scale_y_cut(breaks=c(5000, 10000), which=c(1, 3), scales=c(0.01, 1)) +
    # geom_text(aes(label=mean, y=Inf), color="black", angle=90, size=2, vjust = 0.1) +
    facet_grid(lm_age~risk_class) +
    scale_x_continuous(breaks=c(0:8, 9), labels=c(0:8, '>8')) + 
    scale_y_continuous(labels = percent) +
    xlab('Years') +
    ylab('Percentage') +
    ggtitle(gender)+
    # ggtitle(paste('Years to cross for people who crossed (males)')) +
    scale_color_manual(values=risk_gr_colors[-1], labels=risk_labels) +
    scale_fill_manual(values=risk_gr_colors[-1], labels=risk_labels) +
    geom_text(data = N_tstar_all %>% filter(sex==gender), 
              aes(x = 4, y = 0.9, label = paste0('n=', Total)), size = 3) +
    theme_bw()+
    theme(legend.position="bottom")
} 

risk_labels <- c('High(11.25%-15%)','Medium_high(7.5%-11.25%)',
                 'Medium_low(3.75%-7.5%)','Low(<3.75%)')
p_tstar_female <- p_tstar_plot('female', risk_labels)
p_tstar_male <- p_tstar_plot('male', risk_labels)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p_tstar_male)


jpeg(file = '~/epi_paper/outp_figs/t_detail_15.jpeg', units="in", width=9, height=12, res=300)
grid.arrange(arrangeGrob(p_tstar_female + theme(legend.position="none"),
                         p_tstar_male + theme(legend.position="none"), nrow=1, widths=c(1,1)),
             mylegend, nrow=2,heights=c(20, 1))
dev.off()

##########################################################
# Expected crossing time aggregated across age groups within each risk catogory

# INPUT AGE PROPORTION HERE!
# the input of age proportion does not necessarily add up to 1
# female and male weights can be on different scales
age_proportion_male <- tibble(lm_age=seq(40,80,by=5), sex='male',
                              w=c(436763,400353,442966,461041,400251,335745,316491,258938,174151))
age_proportion_female <- tibble(lm_age=seq(40,80,by=5), sex='female',
                                w=c(445589,407424,458714,477697,415958,353453,344326,291001,212739))
age_proportion <- rbind(age_proportion_female, age_proportion_male)
##############

t_summary_age_risk <- tstar %>% filter(risk_class!='Very_high') %>% 
  mutate(delta=ifelse(pred_status==0, 8, delta)) %>% group_by(lm_age, risk_class, sex) %>% 
  summarise(mean=mean(delta), median=median(delta), Q1=quantile(delta,0.25), Q3=quantile(delta,0.75))

weight <- merge(N_cross_riskgr, age_proportion, by=c('lm_age', 'sex'))
weight <- weight %>% mutate(w=w*riskgr_pct_total) %>% select(lm_age,sex,risk_class,w) %>% 
  mutate(risk_class= as.factor(risk_class)) %>% distinct %>% as_tibble()

t_summary_risk <- merge(t_summary_age_risk, weight, by=c('lm_age', 'sex', 'risk_class'))
t_summary_risk <- t_summary_risk %>% group_by(sex, risk_class) %>% mutate(w=w/sum(w)) %>% 
  summarise(mean=sum(mean*w),median=sum(median*w), Q1=sum(Q1*w), Q3=sum(Q3*w))

levels(t_summary_risk$risk_class) <- c('Very_high(>15%)','High(11.25%-15%)','Medium_high(7.5%-11.25%)',
                                       'Medium_low(3.75%-7.5%)','Low(<3.75%)')

jpeg(file = '~/epi_paper/outp_figs/t_aggregated_15.jpeg', units="in", width=4, height=4, res=300)
ggplot(t_summary_risk, aes(x=risk_class, y=median, color=sex, group=sex)) + 
  geom_point(position=position_dodge(0.2))+
  geom_errorbar(aes(ymin=Q1, ymax=Q3), width=.2,
                position=position_dodge(0.2))+
  xlab('Risk class') +
  ylab('Expected crossing time') +
  ylim(c(0,8)) +
  theme_bw()  +
  theme(axis.text.x = element_text(angle = 25, size =7))
dev.off()