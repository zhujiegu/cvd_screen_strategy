library(data.table)
library(tidyverse)
library(gridExtra)
library(magrittr)
library(scales)
library(dplyr)
library(openxlsx)

in_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/3_cox_outp/'
in_path3 <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/4_screening_outp/'
out_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/5_crosstime_summary/'

lm_age = seq(40,80,5)
gender = "male" #female or male
# deriv_tmp = "derivation" #derivation or validation
##########################################################
# INPUT AGE PROPORTION HERE!
# the input of age proportion does not necessarily add up to 1
# female and male weights can be on different scales
age_proportion_male <- tibble(lm_age=seq(40,80,by=5), sex='male',
                              n=c(436763,400353,442966,461041,400251,335745,316491,258938,174151))
age_proportion_female <- tibble(lm_age=seq(40,80,by=5), sex='female',
                                n=c(445589,407424,458714,477697,415958,353453,344326,291001,212739))
age_proportion <- rbind(age_proportion_female, age_proportion_male)
##############################################################################
# load crossing time
tstar = NULL

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
# get_limits_mean_pred_time = function(x){ list(mean  = mean(x),
#                                               lower = quantile( x, 0.25 ), 
#                                               upper = quantile( x, 0.75 ),
#                                               size = length(x)) }

tstar = data.table(tstar)
tstar$delta = tstar$pred_time-tstar$lm_age

# Proportion that did not cross the threshold in 8 years
tstar %<>% mutate(pred_status = replace(pred_status, delta > 8, 0))
levels(tstar$risk_class) <- c('Very_high','High','Medium_high',
                              'Medium_low','Low')
# levels(tstar$risk_class) <- c('Very_high(>5%)','High(3.75%-5%)','Medium_high(2.5%-3.75%)',
#                                        'Medium_low(1.25%-2.5%)','Low(<1.25%)')
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
##############################################################################
# Number/percentages of people crossed in each subgroup
##############################################################################
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
N_cross_riskgr %>% filter(risk_class=='Very_high') %>% filter(sex=='male')
# write.csv(N_cross_riskgr, file = '~/epi_paper/riskgr_pct.csv')
riskgr_pct <- N_cross_riskgr %>% select(lm_age,sex,risk_class,riskgr_pct_total) %>% unique

##############################################################################
# tstar plot
##############################################################################
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
    ggtitle(ifelse(gender=='male', 'Men', 'Women'))+
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

mylegend<-g_legend(p_tstar_male)


jpeg(file = '~/epi_paper/outp_figs/t_detail_5.jpeg', units="in", width=9, height=12, res=600)
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


##############################################################################
# remove male 65 medium low where only 1 sample available
##############################################################################
tstar %<>% filter(sex!='male'|lm_age!=65|risk_class!='Medium_low')
riskgr_pct <- tstar %>% 
  group_by(lm_age, sex, risk_class, pred_status) %>% summarize(Total = n()) %>% 
  group_by(lm_age, sex) %>% mutate(N = sum(Total), riskgr_pct = Total/sum(Total)) %>% 
  group_by(lm_age, sex, risk_class) %>% mutate(riskgr_pct_total = sum(riskgr_pct), 
                                               cross_pct = Total/sum(Total)) %>% 
  select(lm_age,sex,risk_class, riskgr_pct_total) %>% unique

##############################################################################
# Quantile plot
##############################################################################
tstar_qt = tstar %>% filter(risk_class != 'Very_high') %>% 
  mutate(delta = ifelse(pred_status==1, delta, 8)) %>% 
  group_by(sex, lm_age, risk_class) %>% summarise(mean = mean(delta), Q50 = quantile(delta, probs=0.5),
                                                  Q25 = quantile(delta, probs=0.25),Q10 = quantile(delta, probs=0.1),
                                                  Q05 = quantile(delta, probs=0.05)) %>% ungroup

tstar_qt <- merge(tstar_qt, riskgr_pct)
tstar_qt <- merge(tstar_qt, age_proportion)
tstar_qt %<>% mutate(w=riskgr_pct_total*n)
# age range
tstar_qt %<>% mutate(age_gr = ifelse(lm_age<65, '<65', '>=65'))
tstar_qt %>% head

# aggregation 2 age groups
tstar_qt_agg <- tstar_qt %>% group_by(sex, risk_class, age_gr) %>% mutate(Q50 = sum(Q50*w)/sum(w),
                                                          Q25 = sum(Q25*w)/sum(w),
                                                          Q10 = sum(Q10*w)/sum(w),
                                                          Q05 = sum(Q05*w)/sum(w)) %>% 
  select(-lm_age,-mean,-riskgr_pct_total,-n,-w) %>% unique
# aggregation all ages
tstar_qt_agg2 <- tstar_qt %>% group_by(sex, risk_class) %>% mutate(Q50 = sum(Q50*w)/sum(w),
                                                                  Q25 = sum(Q25*w)/sum(w),
                                                                  Q10 = sum(Q10*w)/sum(w),
                                                                  Q05 = sum(Q05*w)/sum(w)) %>% 
  mutate(age_gr='All') %>% select(-lm_age,-mean,-riskgr_pct_total,-n,-w) %>% unique

tstar_qt_agg <- rbind(tstar_qt_agg, tstar_qt_agg2)

tstar_qt_agg %<>% pivot_longer(cols = starts_with('Q'), names_to = 'Quantile', values_to = 'value')
tstar_qt_agg %<>% mutate_at(vars(risk_class, age_gr, Quantile), factor)
tstar_qt_agg %>% head

p_qt_female <- ggplot(tstar_qt_agg %>% filter(sex=='female'), aes(x=risk_class, y=value, color=Quantile, group=Quantile))+
  # geom_line(linetype = 2)+
  geom_point(aes(shape = Quantile), size=3)+
  facet_grid(~ age_gr)+
  ylab('Time')+
  xlab('Risk class')+
  ggtitle('Quantiles of t* (women)')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p_qt_male <- ggplot(tstar_qt_agg %>% filter(sex=='male'), aes(x=risk_class, y=value, color=Quantile, group=Quantile))+
  # geom_line(linetype = 2)+
  geom_point(aes(shape = Quantile), size=3)+
  facet_grid(~ age_gr)+
  ylab('Time')+
  xlab('Risk class')+
  ggtitle('Quantiles of t* (men)')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

mylegend<-g_legend(p_qt_male)

jpeg(file = '~/epi_paper/outp_figs/t_quantiles.jpeg', units="in", width=6, height=6, res=600)
grid.arrange(arrangeGrob(p_qt_female + theme(legend.position="none"),
                         p_qt_male + theme(legend.position="none"), nrow=1, widths=c(1,1)),
             mylegend, ncol=2, widths=c(8,1))
dev.off()
##############################################################################
# Strategies
##############################################################################
# based on Q50
tstar_qt_agg %>% filter(Quantile=='Q50', age_gr=='All') #

# based on Q10
tstar_qt_agg %>% filter(Quantile=='Q10', age_gr=='All') #4
# age stratified
tstar_qt_agg %>% filter(Quantile=='Q10', age_gr!='All') #

stg_l <- list('female'=list('5-5-5-5'=rep(5,8),
                            '2-6-8-8'=c(2,2,6,6,8,8,8,8),
                            '1-4-8-8'=c(1,1,4,4,8,8,8,8),
                            '1-5(3)-8(6)-8'=c(1,1,5,3,8,6,8,8)),
              'male'=list('5-5-5-5'=rep(5,8),
                            '2-5-8-8'=c(2,2,5,5,8,8,8,8),
                            '1-4-8-8'=c(1,1,4,4,8,8,8,8),
                            '1-4(2)-8-8'=c(1,1,4,2,8,8,8,8)))
#########################################
# Proportion of t* < first assessment
#########################################
# input strategy, a vector of length length 8 for low(<65), low(>=65), ... high(>=65)
prop_inneed <- function(stg){
  prop = tstar %>% filter(risk_class != 'Very_high') %>% 
    mutate(delta = ifelse(pred_status==1, delta, 8)) %>% 
    group_by(sex, lm_age, risk_class) %>% 
    summarise(mean = mean(delta <= case_when(risk_class == 'High' & lm_age<65 ~ stg[1],
                                             risk_class == 'High' & lm_age>=65 ~ stg[2],
                                             risk_class == 'Medium_high' & lm_age<65 ~ stg[3],
                                             risk_class == 'Medium_high' & lm_age>=65 ~ stg[4],
                                             risk_class == 'Medium_low' & lm_age<65 ~ stg[5],
                                             risk_class == 'Medium_low' & lm_age>=65 ~ stg[6],
                                             risk_class == 'Low' & lm_age<65 ~ stg[7],
                                             risk_class == 'Low' & lm_age>=65 ~ stg[8]))) %>% 
    ungroup %>% pivot_wider(names_from = risk_class, values_from = mean)
  return(prop)
}

eff_tables <- function(stg, gender){
  prop <- prop_inneed(stg) %>% filter(sex==gender)
  # write.csv(prop_, file = paste0(out_path, 'prop_inneed_',paste0(stg, collapse = ''),'.csv'))
  prop_ <- prop %>% pivot_longer(cols = High:Low, names_to = 'risk_class', values_to = 'Proportion') 
  
  prop_ <- merge(prop_, riskgr_pct)
  prop_ <- merge(prop_, age_proportion)
  prop_ %<>% mutate(w=riskgr_pct_total*n)
  # age range
  prop_ %<>% mutate(age_gr = ifelse(lm_age<65, '<65', '>=65'))
  
  # aggregation 2 age groups
  prop_agg <- prop_ %>% group_by(sex, risk_class, age_gr) %>% mutate(Proportion = sum(Proportion*w)/sum(w)) %>% 
    select(-lm_age,-riskgr_pct_total,-n,-w) %>% unique
  # aggregation all ages
  prop_agg2 <- prop_ %>% group_by(sex, risk_class) %>% mutate(Proportion = sum(Proportion*w)/sum(w)) %>% 
    mutate(age_gr='All') %>% select(-lm_age,-riskgr_pct_total,-n,-w) %>% unique
  
  prop_agg <- rbind(prop_agg, prop_agg2)
  prop_agg$risk_class %<>% factor(levels = c('High', 'Medium_high', 'Medium_low', 'Low'))
  prop_agg$age_gr %<>% factor
  return(list(prop = prop, prop_agg=prop_agg))
}

# Initiation
inneed_table_female <- vector(mode = 'list')
inneed_table_male <- vector(mode = 'list')
outp_exl_female <- createWorkbook()
outp_exl_male <- createWorkbook()

for (j in 1:length(stg_l$female)) {
  eff <- eff_tables(stg_l$female[[j]], 'female')
  addWorksheet(outp_exl_female, names(stg_l$female)[j])
  writeData(outp_exl_female, sheet = names(stg_l$female)[j], x = eff$prop)
  inneed_table_female[[j]] <- eff$prop_agg
}
for (j in 1:length(stg_l$male)) {
  eff <- eff_tables(stg_l$male[[j]], 'male')
  addWorksheet(outp_exl_male, names(stg_l$male)[j])
  writeData(outp_exl_male, sheet = names(stg_l$male)[j], x = eff$prop)
  inneed_table_male[[j]] <- eff$prop_agg
}

# save table for each stategy
saveWorkbook(outp_exl_female, paste0(out_path, 'prop_inneed_female.xlsx'))
saveWorkbook(outp_exl_male, paste0(out_path, 'prop_inneed_male.xlsx'))

# plot
names(inneed_table_female) <- names(stg_l$female)
names(inneed_table_male) <- names(stg_l$male)
inneed_table_female <- do.call(rbind, lapply(names(inneed_table_female), 
                                             function(e) inneed_table_female[[e]] %>% mutate(Strategy = e)))
inneed_table_male <- do.call(rbind, lapply(names(inneed_table_male), 
                                             function(e) inneed_table_male[[e]] %>% mutate(Strategy = e)))
inneed_table_female$Strategy %<>% factor(levels = names(stg_l$female))
inneed_table_male$Strategy %<>% factor(levels = names(stg_l$male))


p_inneed_male <- ggplot(inneed_table_male, aes(x=risk_class, y=Proportion, color=Strategy, fill=Strategy))+
  geom_point(aes(shape = Strategy), size=3, alpha=0.8, position=position_jitter(h=0,w=0.2))+
  scale_shape_manual(values=c(21:24))+
  # geom_segment(data=filter(inneed_table_male, Strategy=='5-5-5-5'),
  #              aes(x=as.numeric(risk_class)-0.25, xend=as.numeric(risk_class)+0.25, y=Proportion, yend=Proportion))+
  facet_grid(~ age_gr)+
  ylab('Proportion')+
  xlab('Risk class')+
  ggtitle('Proportion of t*<t1 (men)')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p_inneed_female <- ggplot(inneed_table_female, aes(x=risk_class, y=Proportion, color=Strategy, fill=Strategy))+
  geom_point(aes(shape = Strategy), size=3, alpha=0.8, position=position_jitter(h=0,w=0.2))+
  scale_shape_manual(values=c(21:24))+
  # geom_segment(data=filter(inneed_table_female, Strategy=='5-5-5-5'),
  #              aes(x=as.numeric(risk_class)-0.25, xend=as.numeric(risk_class)+0.25, y=Proportion, yend=Proportion))+
  facet_grid(~ age_gr)+
  ylab('Proportion')+
  xlab('Risk class')+
  ggtitle('Proportion of t*<t1 (women)')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

jpeg(file = '~/epi_paper/outp_figs/inneed.jpeg', units="in", width=10, height=5, res=600)
grid.arrange(p_inneed_female, p_inneed_male, ncol=2)
dev.off()

#########################################
# Avg time of staying in high risk: mean (t1-t* | t*<t1)
#########################################
tstar_undetect <- function(stg){
  tstar_inneed <- tstar %>% filter(risk_class != 'Very_high') %>% 
    mutate(delta = ifelse(pred_status==1, delta, 8)) %>% 
    filter(delta <= case_when(risk_class == 'High' & lm_age<65 ~ stg[1],
                              risk_class == 'High' & lm_age>=65 ~ stg[2],
                              risk_class == 'Medium_high' & lm_age<65 ~ stg[3],
                              risk_class == 'Medium_high' & lm_age>=65 ~ stg[4],
                              risk_class == 'Medium_low' & lm_age<65 ~ stg[5],
                              risk_class == 'Medium_low' & lm_age>=65 ~ stg[6],
                              risk_class == 'Low' & lm_age<65 ~ stg[7],
                              risk_class == 'Low' & lm_age>=65 ~ stg[8]))
  
  tstar_undet = tstar_inneed %>% group_by(sex, lm_age, risk_class) %>% 
    mutate(undetect = case_when(risk_class == 'High' & lm_age<65 ~ stg[1],
                                risk_class == 'High' & lm_age>=65 ~ stg[2],
                                risk_class == 'Medium_high' & lm_age<65 ~ stg[3],
                                risk_class == 'Medium_high' & lm_age>=65 ~ stg[4],
                                risk_class == 'Medium_low' & lm_age<65 ~ stg[5],
                                risk_class == 'Medium_low' & lm_age>=65 ~ stg[6],
                                risk_class == 'Low' & lm_age<65 ~ stg[7],
                                risk_class == 'Low' & lm_age>=65 ~ stg[8])-delta) %>% 
    summarise(mean = mean(undetect)) %>% ungroup %>% 
    pivot_wider(names_from = risk_class, values_from = mean)
  return(tstar_undet)
}

wtime_tables <- function(stg_name, gender){
  stg <- stg_l[[gender]][[stg_name]]
  wtime <- tstar_undetect(stg) %>% filter(sex==gender)
  # write.csv(wtime_, file = paste0(out_path, 'wtime_inneed_',paste0(stg, collapse = ''),'.csv'))
  wtime_ <- wtime %>% pivot_longer(cols = High:Low, names_to = 'risk_class', values_to = 'wait_time') 
  
  cross_pct <- read.xlsx(paste0(out_path, 'prop_inneed_', gender, '.xlsx'), sheet = stg_name) %>% 
    pivot_longer(cols = High:Low, names_to = 'risk_class', values_to = 'crossed')
  
  wtime_ <- merge(wtime_, riskgr_pct)
  wtime_ <- merge(wtime_, age_proportion)
  wtime_ <- merge(wtime_, cross_pct)
  wtime_ %<>% mutate(w=n*riskgr_pct_total*crossed)
  # age range
  wtime_ %<>% mutate(age_gr = ifelse(lm_age<65, '<65', '>=65'))
  wtime_[is.na(wtime_)] <- 0
  
  # aggregation 2 age groups
  wtime_agg <- wtime_ %>% group_by(sex, risk_class, age_gr) %>% mutate(wait_time = sum(wait_time*w)/sum(w)) %>% 
    select(-lm_age,-riskgr_pct_total,-n,-w, -crossed) %>% unique
  # aggregation all ages
  wtime_agg2 <- wtime_ %>% group_by(sex, risk_class) %>% mutate(wait_time = sum(wait_time*w)/sum(w)) %>% 
    mutate(age_gr='All') %>% select(-lm_age,-riskgr_pct_total,-n,-w, -crossed) %>% unique
  
  wtime_agg <- rbind(wtime_agg, wtime_agg2)
  wtime_agg$risk_class %<>% factor(levels = c('High', 'Medium_high', 'Medium_low', 'Low'))
  wtime_agg$age_gr %<>% factor
  return(list(wtime = wtime, wtime_agg=wtime_agg))
}

# Initiation
wtime_table_female <- vector(mode = 'list')
wtime_table_male <- vector(mode = 'list')
outp_wt_exl_female <- createWorkbook()
outp_wt_exl_male <- createWorkbook()

for (j in 1:length(stg_l$female)) {
  eff <- wtime_tables(names(stg_l$female)[j], 'female')
  addWorksheet(outp_wt_exl_female, names(stg_l$female)[j])
  writeData(outp_wt_exl_female, sheet = names(stg_l$female)[j], x = eff$wtime)
  wtime_table_female[[j]] <- eff$wtime_agg
}
for (j in 1:length(stg_l$male)) {
  eff <- wtime_tables(names(stg_l$male)[j], 'male')
  addWorksheet(outp_wt_exl_male, names(stg_l$male)[j])
  writeData(outp_wt_exl_male, sheet = names(stg_l$male)[j], x = eff$wtime)
  wtime_table_male[[j]] <- eff$wtime_agg
}

# save table for each stategy
saveWorkbook(outp_wt_exl_female, paste0(out_path, 'wtime_inneed_female.xlsx'))
saveWorkbook(outp_wt_exl_male, paste0(out_path, 'wtime_inneed_male.xlsx'))

# plot
names(wtime_table_female) <- names(stg_l$female)
names(wtime_table_male) <- names(stg_l$male)
wtime_table_female <- do.call(rbind, lapply(names(wtime_table_female), 
                                             function(e) wtime_table_female[[e]] %>% mutate(Strategy = e)))
wtime_table_male <- do.call(rbind, lapply(names(wtime_table_male), 
                                           function(e) wtime_table_male[[e]] %>% mutate(Strategy = e)))
wtime_table_female$Strategy %<>% factor(levels = names(stg_l$female))
wtime_table_male$Strategy %<>% factor(levels = names(stg_l$male))


p_wtime_male <- ggplot(wtime_table_male, aes(x=risk_class, y=wait_time, color=Strategy, fill=Strategy))+
  geom_point(aes(shape = Strategy), size=3, alpha=0.8, position=position_jitter(h=0,w=0.2))+
  scale_shape_manual(values=c(21:24))+
  # geom_segment(data=filter(wtime_table_male, Strategy=='5-5-5-5'),
  #              aes(x=as.numeric(risk_class)-0.25, xend=as.numeric(risk_class)+0.25, y=wtimeortion, yend=wtimeortion))+
  facet_grid(~ age_gr)+
  ylab('Waiting time')+
  xlab('Risk class')+
  ggtitle('Average waiting time given t*<t1 (men)')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p_wtime_female <- ggplot(wtime_table_female, aes(x=risk_class, y=wait_time, color=Strategy, fill=Strategy))+
  geom_point(aes(shape = Strategy), size=3, alpha=0.8, position=position_jitter(h=0,w=0.2))+
  scale_shape_manual(values=c(21:24))+
  # geom_segment(data=filter(wtime_table_female, Strategy=='5-5-5-5'),
  #              aes(x=as.numeric(risk_class)-0.25, xend=as.numeric(risk_class)+0.25, y=wtimeortion, yend=wtimeortion))+
  facet_grid(~ age_gr)+
  ylab('Waiting time')+
  xlab('Risk class')+
  ggtitle('Average waiting time given t*<t1 (women)')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

jpeg(file = '~/epi_paper/outp_figs/wtime.jpeg', units="in", width=10, height=5, res=600)
grid.arrange(p_wtime_female, p_wtime_male, ncol=2)
dev.off()

#########################################
# Proportion of T_wait > 1
#########################################
# tstar_uncovered <- function(stg){
#   # Conditional on t* <= first assessment
#   tstar_inneed <- tstar %>% filter(risk_class != 'Very_high') %>% 
#     mutate(delta = ifelse(pred_status==1, delta, 8)) %>% 
#     filter(delta <= case_when(risk_class == 'Low' & lm_age<65 ~ stg[1],
#                               risk_class == 'Low' & lm_age>=65 ~ stg[2],
#                               risk_class == 'Medium_low' & lm_age<65 ~ stg[3],
#                               risk_class == 'Medium_low' & lm_age>=65 ~ stg[4],
#                               risk_class == 'Medium_high' & lm_age<65 ~ stg[5],
#                               risk_class == 'Medium_high' & lm_age>=65 ~ stg[6],
#                               risk_class == 'High' & lm_age<65 ~ stg[7],
#                               risk_class == 'High' & lm_age>=65 ~ stg[8]))
#   # window of 1 year
#   stg = stg -1
#   tstar_uncover = tstar_inneed %>% group_by(sex, lm_age, risk_class) %>% 
#     summarise(mean = 1-mean(ifelse(delta < case_when(risk_class == 'Low' & lm_age<65 ~ stg[1],
#                                                      risk_class == 'Low' & lm_age>=65 ~ stg[2],
#                                                      risk_class == 'Medium_low' & lm_age<65 ~ stg[3],
#                                                      risk_class == 'Medium_low' & lm_age>=65 ~ stg[4],
#                                                      risk_class == 'Medium_high' & lm_age<65 ~ stg[5],
#                                                      risk_class == 'Medium_high' & lm_age>=65 ~ stg[6],
#                                                      risk_class == 'High' & lm_age<65 ~ stg[7],
#                                                      risk_class == 'High' & lm_age>=65 ~ stg[8]), 0, 1))) %>% 
#     ungroup %>% pivot_wider(names_from = risk_class, values_from = mean) %>%
#     select(sex, lm_age, Low, Medium_low, Medium_high, High)
#   return(tstar_uncover)
# }
# 






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











##########################################################################
# #Covariates matrix
# cov = NULL
# for( i in lm_age ){
#   load(paste0(in_path, "/covariates_matrix_RIRS_", i,"_", gender, ".RData"))
#   covariates_matrix = covariates_matrix[, 1:21]
#   covariates_matrix$lm_age = i
#   covariates_matrix$sex = gender
#   cov = rbind( cov, covariates_matrix )
#   rm(covariates_matrix)
# }
# 
# binary_variables = c("diab_ind", "bp_med", "renal_disease", "atrial_fibrillation", "rheumatoid_arthritis",
#                      "Severe_mental_illness_ind", "Migraine_ind", "depression_ind", "family_history")
# continuous_variables = c("Townsend", "bmi_blup_0", "hdl_blup_0", "sbp_blup_0", "smoke_blup_0", "tchol_blup_0" )
# cov = data.table(cov)
# cont_var_summary = function(x){ paste0(  round(mean(x),2), " (", round(sd(x),2), ")" ) }
# cont_desc = cov[, lapply(.SD, cont_var_summary), by=lm_age, .SDcols = continuous_variables ]
# binary_var_summary = function(x){ paste0( length(which(x == 1)), " (", round( length(which(x == 1))/length(x)*100,2), "%)"  ) }
# binary_desc = cov[,lapply(.SD, binary_var_summary), by=lm_age, .SDcols = binary_variables ]
# final_descriptive = left_join( cont_desc, binary_desc, by = c( "lm_age" ) )
# # write.table( t(final_descriptive), file = paste0( "img_epi_paper/final_descriptive_", gender, "_covariates_", deriv_tmp, ".txt"), quote = F, col.names = F, sep = ",")

##########################################################################
# #Outcome description
# outcome = NULL
# for( i in lm_age ){
#   load(paste0(in_path, "risk_pred_FE4_RIRS_", i,"_", gender, ".RData"))
#   outcome_matrix = risk_prediction[ !is.na(`5.y.risk.0`), 1:8 ]
#   outcome_matrix$lm_age = i
#   outcome_matrix$sex = gender
#   outcome = rbind( outcome, outcome_matrix )
#   rm(outcome_matrix)
# }
# 
# outcome
# outcome[ , fatal_cvd := ifelse( cvd_age == death_age & cvd_ind == 1 & death_ind == 1, 1, 0 ) ]
# outcome[ , non_fatal_cvd := ifelse( (cvd_age < death_age | is.na(death_age)) & cvd_ind == 1, 1, 0 ) ]
# outcome[ , death_other_causes := ifelse( (cvd_age < death_age | is.na(cvd_age)) & death_ind == 1, 1, 0 ) ]
# outcome[ , followup := end_age - lm_age  ]
# 
# outcome_binary_desc = outcome[,lapply(.SD, binary_var_summary), by=lm_age, .SDcols = c("fatal_cvd", "non_fatal_cvd", "death_other_causes") ] #"cvd_ind", "death_ind", 
# outcome_binary_desc$total = outcome[,lapply(.SD, length), by=lm_age, .SDcols = "patid"]$patid
# 
# outcome_cont_desc = outcome[, lapply(.SD, cont_var_summary), by=lm_age, .SDcols = "followup" ]
# final_descriptive_outcome = left_join( outcome_cont_desc, outcome_binary_desc, by = "lm_age" )
