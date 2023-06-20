library(data.table)
library(nlme)
library(survival)
library(tidyverse)
library(haven)
library(miceadds)
library(gridExtra)
library(magrittr)
library(scales)
library(RColorBrewer)
library(dplyr)
summarize <- dplyr::summarize
select <- dplyr::select

in_path2 <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/2_lmfit_outp/'
in_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/3_cox_outp/'
in_path3 <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/4_screening_outp/'
out_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/3_cox_outp/'

lm_age = seq(40,80,5)
gender = "male" #female or male
# deriv_tmp = "derivation" #derivation or validation
cov = NULL
outcome = NULL
tstar = NULL
beta_5_cox = NULL 
beta_LMEM = NULL

#Covariates matrix
for( i in lm_age ){
  load(paste0(in_path, "/covariates_matrix_RIRS_", i,"_", gender, ".RData"))
  covariates_matrix = covariates_matrix[, 1:21]
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
    load(paste0(in_path3, "comparing_schemes_RIRS_apply_", i,"_", j, ".RData"))
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

N_cross_riskgr
N_cross_riskgr %>% filter(risk_class=='Very_high') %>% filter(sex=='female')
# write.csv(N_cross_riskgr, file = '~/epi_paper/riskgr_pct.csv')
  
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
    facet_grid(lm_age~risk_class, drop=F) +
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

risk_labels <- c('High(3.75%-5%)','Medium_high(2.5%-3.75%)', 'Medium_low(1.25%-2.5%)','Low(<1.25%)')
p_tstar_female <- p_tstar_plot('female', risk_labels)
p_tstar_male <- p_tstar_plot('male', risk_labels)

  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  
  mylegend<-g_legend(p_tstar_male)
  

jpeg(file = '~/epi_paper/outp_figs/t_detail_5.jpeg', units="in", width=9, height=12, res=300)
grid.arrange(arrangeGrob(p_tstar_female + theme(legend.position="none"),
                         p_tstar_male + theme(legend.position="none"), nrow=1, widths=c(1,1)),
             mylegend, nrow=2,heights=c(20, 1))
dev.off()

# # proportion of risk class within lm_age, and proportion of crossing
# p1 <- ggplot(N_cross_riskgr %>% filter(sex=='male') %>% filter(lm_age <=ifelse(gender=='male',70,75)), 
#              aes(x=riskgr_pct, y=risk_class, color=risk_class, fill = risk_class, alpha = crossed))+ 
#   geom_bar(stat="identity")+
#   geom_text(aes(x= 1, y=-Inf, label = paste0('N=',N)), size=4, color='black',
#             hjust = 1, vjust =-3, check_overlap = TRUE) +
#   facet_grid(lm_age~.)+
#   scale_y_discrete(limits=rev)+
#   scale_color_manual(values=risk_gr_colors) +
#   scale_fill_manual(values=risk_gr_colors) +
#   scale_alpha_discrete(range = c(0.35, 0.9)) +
#   geom_text(data=N_cross_riskgr %>% filter(sex=='male') %>% filter(crossed==TRUE)%>% 
#               filter(lm_age <=ifelse(gender=='male',70,75)),
#             aes(x= riskgr_pct_total, label = paste0(round(100*cross_pct, digit=1), '%')), size=3,
#             hjust = -0.1, vjust =0.4) + 
#   theme_bw()+
#   ggtitle(paste('')) +
#   xlab('Percentage')+
#   theme(axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         legend.position="bottom")+ 
#   guides(alpha = 'none')

##########################################################

tstar_qt = tstar %>% filter(risk_class != 'Very_high') %>% 
  mutate(delta = ifelse(pred_status==1, delta, 8)) %>% 
  group_by(sex, lm_age, risk_class) %>% summarise(mean = mean(delta), qt_50 = quantile(delta, probs=0.5),
                                          qt_75 = quantile(delta, probs=0.25),qt_90 = quantile(delta, probs=0.1),
                                          qt_95 = quantile(delta, probs=0.05)) %>% ungroup
tstar_qt %>% filter(sex=='male', risk_class=='Medium_low')
  
# Proportion of t* < first assessment
prop_inneed <- function(stg){
  prop = tstar %>% filter(risk_class != 'Very_high') %>% 
    mutate(delta = ifelse(pred_status==1, delta, 8)) %>% 
    group_by(sex, lm_age, risk_class) %>% 
    summarise(mean = mean(delta <= case_when(risk_class == 'Low' & lm_age<65 ~ stg[1],
                                                     risk_class == 'Low' & lm_age>=65 ~ stg[2],
                                                     risk_class == 'Medium_low' & lm_age<65 ~ stg[3],
                                                     risk_class == 'Medium_low' & lm_age>=65 ~ stg[4],
                                                     risk_class == 'Medium_high' & lm_age<65 ~ stg[5],
                                                     risk_class == 'Medium_high' & lm_age>=65 ~ stg[6],
                                                     risk_class == 'High' & lm_age<65 ~ stg[7],
                                                     risk_class == 'High' & lm_age>=65 ~ stg[8]))) %>% 
    ungroup %>% pivot_wider(names_from = risk_class, values_from = mean)
  return(prop)
}


# input strategy, a vector of length length 8 for low(<65), low(>=65), ... high(>=65)

# Avg time of staying in high risk before screen
tstar_undetect <- function(stg){
  tstar_inneed <- tstar %>% filter(risk_class != 'Very_high') %>% 
    mutate(delta = ifelse(pred_status==1, delta, 8)) %>% 
    filter(delta <= case_when(risk_class == 'Low' & lm_age<65 ~ stg[1],
                              risk_class == 'Low' & lm_age>=65 ~ stg[2],
                              risk_class == 'Medium_low' & lm_age<65 ~ stg[3],
                              risk_class == 'Medium_low' & lm_age>=65 ~ stg[4],
                              risk_class == 'Medium_high' & lm_age<65 ~ stg[5],
                              risk_class == 'Medium_high' & lm_age>=65 ~ stg[6],
                              risk_class == 'High' & lm_age<65 ~ stg[7],
                              risk_class == 'High' & lm_age>=65 ~ stg[8]))
  
  tstar_undet = tstar_inneed %>% group_by(sex, lm_age, risk_class) %>% 
    mutate(undetect = case_when(risk_class == 'Low' & lm_age<65 ~ stg[1],
                                risk_class == 'Low' & lm_age>=65 ~ stg[2],
                                risk_class == 'Medium_low' & lm_age<65 ~ stg[3],
                                risk_class == 'Medium_low' & lm_age>=65 ~ stg[4],
                                risk_class == 'Medium_high' & lm_age<65 ~ stg[5],
                                risk_class == 'Medium_high' & lm_age>=65 ~ stg[6],
                                risk_class == 'High' & lm_age<65 ~ stg[7],
                                risk_class == 'High' & lm_age>=65 ~ stg[8])-delta) %>% 
    summarise(mean = mean(undetect)) %>% ungroup %>% 
    pivot_wider(names_from = risk_class, values_from = mean)
  return(tstar_undet)
}

# Proportion of T_wait > 1
tstar_uncovered <- function(stg){
  # Conditional on t* <= first assessment
  tstar_inneed <- tstar %>% filter(risk_class != 'Very_high') %>% 
    mutate(delta = ifelse(pred_status==1, delta, 8)) %>% 
    filter(delta <= case_when(risk_class == 'Low' & lm_age<65 ~ stg[1],
                              risk_class == 'Low' & lm_age>=65 ~ stg[2],
                              risk_class == 'Medium_low' & lm_age<65 ~ stg[3],
                              risk_class == 'Medium_low' & lm_age>=65 ~ stg[4],
                              risk_class == 'Medium_high' & lm_age<65 ~ stg[5],
                              risk_class == 'Medium_high' & lm_age>=65 ~ stg[6],
                              risk_class == 'High' & lm_age<65 ~ stg[7],
                              risk_class == 'High' & lm_age>=65 ~ stg[8]))
  # window of 1 year
  stg = stg -1
  tstar_uncover = tstar_inneed %>% group_by(sex, lm_age, risk_class) %>% 
    summarise(mean = 1-mean(ifelse(delta < case_when(risk_class == 'Low' & lm_age<65 ~ stg[1],
                                                   risk_class == 'Low' & lm_age>=65 ~ stg[2],
                                                   risk_class == 'Medium_low' & lm_age<65 ~ stg[3],
                                                   risk_class == 'Medium_low' & lm_age>=65 ~ stg[4],
                                                   risk_class == 'Medium_high' & lm_age<65 ~ stg[5],
                                                   risk_class == 'Medium_high' & lm_age>=65 ~ stg[6],
                                                   risk_class == 'High' & lm_age<65 ~ stg[7],
                                                   risk_class == 'High' & lm_age>=65 ~ stg[8]), 0, 1))) %>% 
    ungroup %>% pivot_wider(names_from = risk_class, values_from = mean) %>%
    select(sex, lm_age, Low, Medium_low, Medium_high, High)
  return(tstar_uncover)
}


prop_inneed(rep(5,8))
prop_inneed(c(8,8,8,8,6,6,2,2))
prop_inneed(c(8,8,8,8,5,5,2,2))
prop_inneed(c(8,8,8,8,4,4,1,1))
prop_inneed(c(8,8,8,8,5,4,2,1))
prop_inneed(c(8,8,8,8,4,3,1,1))
prop_inneed(c(8,8,8,7,8,5,3,2))
prop_inneed(c(8,8,8,6,5,3,1,1))



tstar_undetect(rep(5,8))
tstar_undetect(c(8,8,8,8,6,6,2,2))
tstar_undetect(c(8,8,8,8,5,5,2,2))
tstar_undetect(c(8,8,8,8,4,4,1,1))
tstar_undetect(c(8,8,8,8,5,4,2,1))
tstar_undetect(c(8,8,8,8,4,3,1,1))
tstar_undetect(c(8,8,8,7,8,5,3,2))
tstar_undetect(c(8,8,8,6,5,3,1,1))


# tstar_uncovered(rep(5,8))
# tstar_uncovered(c(8,8,8,8,6,6,2,2))
# tstar_uncovered(c(8,8,8,8,5,5,2,2))
# tstar_uncovered(c(8,8,8,8,4,4,1,1))
# tstar_uncovered(c(8,8,8,6,5,3,1,1))
# tstar_uncovered(c(8,8,8,8,4,3,1,1))
# tstar_uncovered(c(8,8,8,7,8,5,3,2))
# tstar_uncovered(c(8,8,8,8,5,4,2,1))





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

levels(t_summary_risk$risk_class) <- c('Very_high(>5%)','High(3.75%-5%)','Medium_high(2.5%-3.75%)',
                              'Medium_low(1.25%-2.5%)','Low(<1.25%)')

jpeg(file = '~/epi_paper/outp_figs/t_aggregated.jpeg', units="in", width=4, height=4, res=300)
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

##########################################################
# Same subjects in two cohorts (40, 45)
tstar_overlap_40 <- tstar %>% group_by(patid) %>% filter(any(lm_age==40) & any(lm_age==45)) %>% ungroup %>% 
  filter(lm_age==40)
tstar_overlap_45 <- tstar %>% group_by(patid) %>% filter(any(lm_age==40) & any(lm_age==45)) %>% ungroup %>% 
  filter(lm_age==45)

tstar_overlap_40 %<>% mutate(cross_gr = cut(delta, breaks=c(-Inf, 0, 5, 8, Inf), 
                                            labels=c("already","0-5","5-8",'>8y')))
tstar_overlap_40 %<>% mutate(cross_gr = replace(cross_gr, pred_status == 0, '>8y'))

tstar_overlap_45 %<>% mutate(cross_gr = cut(delta, breaks=c(-Inf, 0, 3, Inf), 
                                            labels=c("already","0-3",'>3y')))
tstar_overlap_45 %<>% mutate(cross_gr = replace(cross_gr, pred_status == 0, '>3y'))

compare_40_45 <- merge(tstar_overlap_40, tstar_overlap_45, by='patid') %>% select(starts_with('cross_gr'))
compare_40_45 %>% group_by(cross_gr.x,cross_gr.y) %>% summarise(N=n())

##########################################################
# LMEM
##########################################################

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



#beta 5 year CVD
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


rownames(beta_5_cox)

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
townsend = ggplot( beta_5_cox[ grep("townsend", rownames(beta_5_cox) ), ], aes(x = lm_age_j, y = exp_var, col = factor(lm_age)) ) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  labs( x = "Landmark age", y = "HR - Townsend20", col = "Risk class") +
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
  labs( x = "Landmark age", y = "HR - diabetes", col = "Risk class") +
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

shared_legend_time_fixed1 = cowplot::get_legend( townsend )

all_beta_cox_time_fixed1 = cowplot::plot_grid(townsend + theme(legend.position = "none"),
                                           bp_bin, diab, depression, Migraine, Severe_mental, align = 'h',  labels=c('A', 'B','C', 'D', 'E',  'F'), ncol = 2, nrow = 4, rel_heights = c(1,1,1,0.2))

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




#plot risk
risk_tot = NULL

for(j in c("female", "male"))
{
  for( i in lm_age )
  {
    load(paste0("results/results_", j, "/risk_pred_FE4_RIRS_",i,"_",j,".RData"))
    risk_temp = risk_prediction[ !is.na( risk_prediction$`5.y.risk.0` ), ]
    risk_temp$lm_age = i
    risk_temp$sex = j
    risk_tot = rbind( risk_tot, risk_temp )
    rm(risk_temp)
  }
}

risk_tot[lm_age == 40 & `5.y.risk.0` <= 0.97, ]
risk_tot[patid =="71091229", ]


risk_to_plot = risk_tot %>% 
  filter( patid %in% c("15187","71091229", "254205", "79435") & lm_age == 40 ) %>%
  select( patid, sex, lm_age, cvd_age, death_age, end_age, paste0( "5.y.risk.", 0:8)) %>%
  pivot_longer( cols = contains( "5.y.risk"), names_to = "risk_date", values_to = "value")
 
risk_to_plot$pred_age = risk_to_plot$lm_age + rep(0:8,4) 
risk_to_plot$event_type = ifelse( !is.na(risk_to_plot$cvd_age), "CVD",
                                  ifelse( !is.na( risk_to_plot$death_age ), "Death", "Censored" ) )

risk_to_plot$first_event_time = apply( risk_to_plot[, c("cvd_age", "death_age", "end_age")],1, function(x) min(x, na.rm = T)) 
  
risk_traj = ggplot( risk_to_plot, aes( x = pred_age, y = 100*(1-value), col = factor(patid))) +
  geom_point() +
  geom_hline(yintercept = 5) +
  geom_vline(aes( xintercept = first_event_time, col = factor(patid), linetype = event_type)) +
  labs(x = "Prediction time", y = "5-year CVD risk [%]") +
  theme_bw()

ggsave(risk_traj, file = paste0( "img_epi_paper/risk_trajectory_40.jpeg" ), width = 15, height = 15, unit = "cm"  )


#EPI paper: plotting scheme for males 60
gender = "male"
i = 60
load(paste0("results/results_", gender, "/fit_LMEM_rirs_", i,"_", gender, ".RData"))
load(paste0("data_created/", gender, "/data_lm_", i,"_", gender, ".RData"))

data_lines = data_ext[patid %in% c("1027", "1329", "1349", "2018"), c("patid", "exp_date","exp_age", "exposure", "bp_bin", "statin_bin")]
data_lines$estim_LMEM = 0

#BMI 
bmi_idx = which(data_lines$exposure == "bmi")
data_lines$estim_LMEM[ bmi_idx ] = fit_rirs$coefficients$fixed[1] + fit_rirs$coefficients$fixed[6] * data_lines$exp_age[ bmi_idx ]
#HDL 
hdl_idx = which(data_lines$exposure == "hdl")
data_lines$estim_LMEM[ hdl_idx ] = fit_rirs$coefficients$fixed[2] + fit_rirs$coefficients$fixed[7] * data_lines$exp_age[ hdl_idx ]
#SBP 
sbp_idx = which(data_lines$exposure == "sbp")
data_lines$estim_LMEM[ sbp_idx ] = fit_rirs$coefficients$fixed[3] + fit_rirs$coefficients$fixed[8] * data_lines$exp_age[ sbp_idx ] + fit_rirs$coefficients$fixed[11] * data_lines$bp_bin[ sbp_idx ]
#smoke 
smoke_idx = which(data_lines$exposure == "smoke")
data_lines$estim_LMEM[ smoke_idx ] = fit_rirs$coefficients$fixed[4] + fit_rirs$coefficients$fixed[9] * data_lines$exp_age[ smoke_idx ] 
#TCHOL 
tchol_idx = which(data_lines$exposure == "tchol")
data_lines$estim_LMEM[ tchol_idx ] = fit_rirs$coefficients$fixed[5] + fit_rirs$coefficients$fixed[10] * data_lines$exp_age[ tchol_idx ] + fit_rirs$coefficients$fixed[12] * data_lines$statin_bin[ tchol_idx ] 

unique(data_ext$patid)
ggplot( subset( data_ext, patid %in% c("1027", "1329", "1349", "2018")), aes(x = exp_date, y = scaled_corr, col = factor(patid) ) ) +
  geom_point() +
  labs(x = "Exposure date", y = "Scaled values") +
  geom_line(data = data_lines, aes(x = exp_date, y = estim_LMEM, col = as.factor(patid))) +
  theme_bw() +
  facet_grid(rows = vars(exposure), scales = "free_y") 
  

#Plot
plot(0:8, fit_rirs$coefficients$fixed[1] + fit_rirs$coefficients$fixed[6] * 0:8, type= "l" )
lines(0:8, fit_rirs$coefficients$fixed[1] + fit_rirs$coefficients$fixed[6] * 0:8  +
        fit_rirs$coefficients$random$patid[1,1] + fit_rirs$coefficients$random$patid[1,6] * 0:8, col = 2 )




data_female = NULL
data_male = NULL



for( i in lm_age )
{
  load(paste0("data_created/female/data_lm_", i,"_female.RData"))
  data_ext$lm_age = i
  data_ext$sex = "female"
  data_female = rbind( data_female, data_ext )
  rm(data_ext)
  
  load(paste0("data_created/male/data_lm_", i,"_male.RData"))
  data_ext$lm_age = i
  data_ext$sex = "male"
  data_male = rbind( data_male, data_ext )
  rm(data_ext)
  
  print(i)
}

#Table

max(data_female$end_age)
summary(data_female$end_age[which(data_female$lm_age == 80)] - data_female$lm_age[ which(data_female$lm_age == 80)])
mean(data_female$end_age[which(data_female$lm_age == 80)] - data_female$lm_age[ which(data_female$lm_age == 80)])
sd(data_female$end_age[which(data_female$lm_age == 80)] - data_female$lm_age[ which(data_female$lm_age == 80)])
summary(data_female$end_age[which(data_female$lm_age == 75)] - data_female$lm_age[ which(data_female$lm_age == 75)])
mean(data_female$end_age[which(data_female$lm_age == 75)] - data_female$lm_age[ which(data_female$lm_age == 75)])
sd(data_female$end_age[which(data_female$lm_age == 75)] - data_female$lm_age[ which(data_female$lm_age == 75)])

age = 80
idx_old_old = which( data_female$end_age[which(data_female$lm_age == age)] - data_female$lm_age[ which(data_female$lm_age == age)] >13.3)
data_female[which(data_female$lm_age == age)[idx_old_old[10]],]

max(data_male$end_age)
summary(data_male$end_age[which(data_male$lm_age == 80)] - data_male$lm_age[ which(data_male$lm_age == 80)])
mean(data_male$end_age[which(data_male$lm_age == 80)] - data_male$lm_age[ which(data_male$lm_age == 80)])
sd(data_male$end_age[which(data_male$lm_age == 80)] - data_male$lm_age[ which(data_male$lm_age == 80)])

summary(data_male$end_age[which(data_male$lm_age == 75)] - data_male$lm_age[ which(data_male$lm_age == 75)])
mean(data_male$end_age[which(data_male$lm_age == 75)] - data_male$lm_age[ which(data_male$lm_age == 75)])
sd(data_male$end_age[which(data_male$lm_age == 75)] - data_male$lm_age[ which(data_male$lm_age == 75)])

age = 80
idx_old_old = which( data_male$end_age[which(data_male$lm_age == age)] - data_male$lm_age[ which(data_male$lm_age == age)] >13.3)
data_male[which(data_male$lm_age == age)[idx_old_old[10]],]

#load("../../patients_outcomes.RData")
#load("exposures_red_merged.RData")
#family_data = read_dta("../familyhistory.dta")


#table( outcomes$pracid[outcomes$derivation == "derivation"] )
