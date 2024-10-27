library(data.table)
library(tidyverse)
library(gridExtra)
library(magrittr)
library(scales)
library(dplyr)
library(openxlsx)
library(grid)
library(ggnewscale)
library(wesanderson)
summarise <- dplyr::summarise
select <- dplyr::select

in_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/3_cox_outp/'
in_path3 <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/4_screening_outp/'
out_path <- '/home/zg304/epi_paper/5_crosstime_summary/'

lm_age = seq(40,80,5)
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
  group_by(lm_age, sex, risk_class, pred_status) %>% summarise(Total = n()) %>% 
  group_by(lm_age, sex) %>% mutate(N = sum(Total), riskgr_pct = Total/sum(Total)) %>% 
  group_by(lm_age, sex, risk_class) %>% mutate(riskgr_pct_total = sum(riskgr_pct), 
                                               cross_pct = Total/sum(Total))

N_cross_riskgr$pred_status <- N_cross_riskgr$pred_status==1
colnames(N_cross_riskgr)[4] <- 'crossed'

N_cross_riskgr
N_cross_riskgr %>% filter(risk_class=='darkred') %>% filter(sex=='male')
write.csv(N_cross_riskgr, file = '~/epi_paper/riskgr_pct.csv')
riskgr_pct <- N_cross_riskgr %>% select(lm_age,sex,risk_class,riskgr_pct_total) %>% unique

##############################################################################
# tstar plot
##############################################################################
tstar_crossed = tstar %>% filter(pred_status==1) %>% filter(risk_class != 'darkred') %>% 
  group_by(lm_age, sex, risk_class) %>% mutate(mean = mean(pred_time)) %>% ungroup
# ggplot(tstar_crossed, aes(x=pred_time, color=risk_class, fill = risk_class)) +
#   geom_density(alpha = 0.2)

tstar_all <- tstar %>% filter(risk_class != 'darkred') %>% group_by(lm_age, sex, risk_class) %>% 
  mutate(mean = mean(pred_time)) %>% ungroup

N_tstar_all <- tstar_all %>% 
  group_by(lm_age, sex, risk_class) %>% summarise(Total = n())

cross_riskgr_pct <- N_cross_riskgr %>% select(lm_age, sex, risk_class, crossed, cross_pct) %>% filter(risk_class!='darkred')
cross_riskgr_pct <- cross_riskgr_pct %>%  # Group by lm_age, sex, and risk_class to ensure we consider each combination
  group_by(lm_age, sex, risk_class) %>%
  # Complete the dataset for crossed within each group
  complete(crossed = c(TRUE, FALSE), fill = list(cross_pct = NA)) %>%
  # Ensure that cross_pct for TRUE and FALSE adds up to 1
  mutate(cross_pct = ifelse(is.na(cross_pct), 1 - sum(cross_pct, na.rm = TRUE), cross_pct)) %>%
  # Ungroup to avoid grouping-related issues later on
  ungroup()
cross_riskgr_pct$risk_class <- droplevels(cross_riskgr_pct$risk_class)


risk_gr_colors <- c("red4", "red","darkorange", "#ffd200", "springgreen3")
names(risk_gr_colors)<- N_cross_riskgr$risk_class %>% levels

# keep empty facet rows to be comparable
tstar_all$lm_age %<>% as.factor()
N_tstar_all$lm_age %<>% as.factor()
tstar_all$risk_class <- droplevels(tstar_all$risk_class)
N_tstar_all$risk_class <- droplevels(N_tstar_all$risk_class)

p_tstar_plot <- function(gender, risk_labels){
  ggplot(tstar_all %>% filter(sex==gender), aes(x=pred_time, color=risk_class, fill = risk_class)) + 
    geom_histogram(alpha = 0.4,  position = 'identity', breaks=seq(0,11,by=1), closed='left', aes(y = ..density..)) +
    geom_vline(aes(xintercept=mean, color=risk_class), linetype = "dashed") +
    # scale_y_cut(breaks=c(5000, 10000), which=c(1, 3), scales=c(0.01, 1)) +
    # geom_text(aes(label=mean, y=Inf), color="black", angle=90, size=2, vjust = 0.1) +
    facet_grid(lm_age ~ risk_class, drop = FALSE, 
      labeller = labeller(risk_class = setNames(risk_labels, levels(tstar_all$risk_class)))
    )+
    scale_x_continuous(breaks=1:10, labels = c('',2,'',4,'',6,'',8,'',10)) + 
    scale_y_continuous(labels = percent) +
    xlab('Years') +
    ylab('Percentage') +
    labs(fill="Risk group", color='Risk group')+
    ggtitle(ifelse(gender=='male', 'Men', 'Women'))+
    # ggtitle(paste('Years to cross for people who crossed (males)')) +
    scale_color_manual(values=risk_gr_colors[-1], labels=risk_labels) +
    scale_fill_manual(values=risk_gr_colors[-1], labels=risk_labels) +
    geom_vline(aes(xintercept=10), linetype = "solid") +
    geom_text(data = N_tstar_all %>% filter(sex==gender), 
              aes(x = 4, y = 0.9, label = paste0('n=', Total)), color = "black", size = 3) +
    geom_text(data = cross_riskgr_pct %>% filter(sex==gender, crossed==FALSE),
              aes(x = 4, y = 0.75, label = paste0(sprintf("%.1f",round(cross_pct *100, digits = 1)), '% >10')),
              color = "black", size = 3) +
    theme_bw()+
    theme(legend.position="bottom")
} 

risk_labels <- c('7.5%-10%','5%-7.5%', '2.5%-5%','<2.5%')
p_tstar_female <- p_tstar_plot('female', risk_labels)
p_tstar_male <- p_tstar_plot('male', risk_labels)

mylegend<-g_legend(p_tstar_male)


jpeg(file = '~/epi_paper/outp_figs/t_detail_10.jpeg', units="in", width=9, height=12, res=600)
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
# remove groups that contain less than 100 samples
##############################################################################
tstar %<>% filter(sex!='male'|lm_age!=55|risk_class!='green')
tstar %<>% filter(sex!='male'|lm_age!=65|risk_class!='orange')
tstar %<>% filter(sex!='female'|lm_age!=75|risk_class!='red')
##############################################################################

riskgr_pct <- tstar %>% 
  group_by(lm_age, sex, risk_class, pred_status) %>% summarise(Total = n()) %>% 
  group_by(lm_age, sex) %>% mutate(N = sum(Total), riskgr_pct = Total/sum(Total)) %>% 
  group_by(lm_age, sex, risk_class) %>% mutate(riskgr_pct_total = sum(riskgr_pct), 
                                               cross_pct = Total/sum(Total)) %>% 
  select(lm_age,sex,risk_class, riskgr_pct_total) %>% unique

weight <- merge(riskgr_pct, age_proportion)
weight_all <- weight %>% mutate(w=riskgr_pct_total*n) %>% group_by(sex, risk_class) %>% summarise(N_riskgr=sum(w)) %>% 
  mutate(age_gr='All')
weight_str <- weight %>% mutate(w=riskgr_pct_total*n) %>% mutate(age_gr = ifelse(lm_age<65, '<65', '>=65')) %>% 
  group_by(sex, risk_class, age_gr) %>% summarise(N_riskgr=sum(w))
weight <- rbind(weight_all, weight_str)
##############################################################################
# Percentile plot
##############################################################################
tstar_qt = tstar %>% filter(risk_class != 'darkred') %>% 
  group_by(sex, lm_age, risk_class) %>% summarise(mean = mean(pred_time), P50 = quantile(pred_time, probs=0.5),
                                                  P25 = quantile(pred_time, probs=0.25),P10 = quantile(pred_time, probs=0.1),
                                                  P5 = quantile(pred_time, probs=0.05)) %>% ungroup
tstar_qt %<>% mutate(risk_class = factor(risk_class, levels = c("red", "orange", "yellow", "green"))) %>%
  arrange(sex, lm_age, risk_class)
# write.xlsx(tstar_qt, file = '~/epi_paper/5_crosstime_summary/tstar_qt.xlsx')

tstar_qt <- merge(tstar_qt, riskgr_pct)
tstar_qt <- merge(tstar_qt, age_proportion)
tstar_qt %<>% mutate(w=riskgr_pct_total*n)
# age range
tstar_qt %<>% mutate(age_gr = ifelse(lm_age<65, '<65', '>=65'))
tstar_qt %>% head

# aggregation 2 age groups
tstar_qt_agg <- tstar_qt %>% group_by(sex, risk_class, age_gr) %>% mutate(P50 = sum(P50*w)/sum(w),
                                                          P25 = sum(P25*w)/sum(w),
                                                          P10 = sum(P10*w)/sum(w),
                                                          P5 = sum(P5*w)/sum(w)) %>% 
  select(-lm_age,-mean,-riskgr_pct_total,-n,-w) %>% unique
# aggregation all ages
tstar_qt_agg2 <- tstar_qt %>% group_by(sex, risk_class) %>% mutate(P50 = sum(P50*w)/sum(w),
                                                                  P25 = sum(P25*w)/sum(w),
                                                                  P10 = sum(P10*w)/sum(w),
                                                                  P5 = sum(P5*w)/sum(w)) %>% 
  mutate(age_gr='All') %>% select(-lm_age,-mean,-riskgr_pct_total,-n,-w) %>% unique

tstar_qt_agg <- rbind(tstar_qt_agg, tstar_qt_agg2)

tstar_qt_agg %<>% pivot_longer(cols = starts_with('P'), names_to = 'Percentile', values_to = 'value')
tstar_qt_agg %<>% mutate_at(vars(risk_class, age_gr, Percentile), factor)
tstar_qt_agg %>% head
tstar_qt_agg$Percentile %<>%  factor(levels = c('P50','P25','P10','P5'))

p_qt_female <- ggplot(tstar_qt_agg %>% filter(sex=='female',Percentile!='P5'), aes(x=risk_class, y=value, color=age_gr, fill=age_gr))+
  geom_bar(stat='identity', position = 'dodge', width=.6, alpha=0.3)+
  facet_grid(~ Percentile)+
  ylab('Time (years)')+
  xlab('Risk group')+
  ggtitle('Women')+
  theme_bw()+
  scale_x_discrete(labels=c("red" = "7.5%-10%", "orange" = "5%-7.5%",
                            "yellow" = "2.5%-5%", "green" = "<2.5%"))+
  scale_color_manual(values = wes_palette("Moonrise3", n = 3)) +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 3)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  labs(fill="Age", color='Age') +
  # Change the y-axis breaks so the top break is "8+"
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10), labels = c("0", "2", "4", "6", "8", "10+")) 

arrows_df <- tstar_qt_agg %>% 
  filter(sex == 'female', value == 10)

# Dodge width for arrows
dodge_width <- 0.6 / length(unique(tstar_qt_agg$age_gr))

p_qt_female <- p_qt_female + geom_segment(data = arrows_df,
                           aes(x = as.numeric(risk_class) - 0.2 + (as.numeric(age_gr) - 1) * dodge_width, 
                               xend = as.numeric(risk_class) - 0.2 + (as.numeric(age_gr) - 1) * dodge_width,
                               y = 10, yend = 10.5, group = age_gr, color=age_gr),
                           arrow = arrow(length = unit(2, "mm")), linewidth = 0.5, inherit.aes = FALSE, 
                           show.legend = FALSE)


p_qt_male <- ggplot(tstar_qt_agg %>% filter(sex=='male',Percentile!='P5'), aes(x=risk_class, y=value, color=age_gr, fill=age_gr))+
  geom_bar(stat='identity', position = 'dodge', width=.6, alpha=0.3)+
  # geom_line(linetype = 2)+
  # geom_point(aes(shape = Percentile), size=3)+
  facet_grid(~ Percentile)+
  ylab('Time (years)')+
  xlab('Risk group')+
  ggtitle('Men')+
  theme_bw()+
  scale_x_discrete(labels=c("red" = "7.5%-10%", "orange" = "5%-7.5%",
                            "yellow" = "2.5%-5%", "green" = "<2.5%"))+
  scale_color_manual(values = wes_palette("Moonrise3", n = 3)) +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 3)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  labs(fill="Age", color='Age') +
  # Change the y-axis breaks so the top break is "8+"
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10), labels = c("0", "2", "4", "6", "8", "10+")) 

arrows_df <- tstar_qt_agg %>% 
  filter(sex == 'male', value == 10)

# Dodge width for arrows
dodge_width <- 0.6 / length(unique(tstar_qt_agg$age_gr))

p_qt_male <- p_qt_male + geom_segment(data = arrows_df,
                                      aes(x = as.numeric(risk_class) - 0.2 + (as.numeric(age_gr) - 1) * dodge_width, 
                                          xend = as.numeric(risk_class) - 0.2 + (as.numeric(age_gr) - 1) * dodge_width,
                                          y = 10, yend = 10.5, group = age_gr, color=age_gr),
                                      arrow = arrow(length = unit(2, "mm")), linewidth = 0.5, inherit.aes = FALSE, 
                                      show.legend = FALSE)

# mylegend<-g_legend(p_qt_male+ labs(fill="Age", color='Age'))

jpeg(file = '~/epi_paper/outp_figs/t_quantiles.jpeg', units="in", width=8.5, height=6.5, res=600)
grid.arrange(p_qt_female, p_qt_male, nrow=2, heights=c(1,1)
             # top=textGrob('Average percentiles of expected crossing time', gp=gpar(fontsize=18,font=8))
             )
dev.off()
##############################################################################
# Strategies
##############################################################################
# based on P50
tstar_qt_agg %>% filter(Percentile=='P50', age_gr=='All') 
tstar_qt_agg %>% filter(Percentile=='P50', age_gr!='All') #age stratified

# based on P25
tstar_qt_agg %>% filter(Percentile=='P25', age_gr=='All') 
tstar_qt_agg %>% filter(Percentile=='P25', age_gr!='All') #age stratified

# based on P10
tstar_qt_agg %>% filter(Percentile=='P10', age_gr=='All') 
tstar_qt_agg %>% filter(Percentile=='P10', age_gr!='All') 

# based on P5
tstar_qt_agg %>% filter(Percentile=='P5', age_gr=='All') 
tstar_qt_agg %>% filter(Percentile=='P5', age_gr!='All') 

# # rounding up
# stg_l <- list('female'=list('5-5-5-5'=rep(5,8), # Benchmark Non-stratified
#                             '2-6-10-10'=c(2,2,6,6,10,10,10,10), #P50 Non-stratified
#                             '3(2)-7(4)-10(6)-10'=c(3,2,7,4,10,6,10,10), #P50 age-stratified
#                             '1-5-9-10'=c(1,1,5,5,9,9,10,10), #P25 Non-stratified
#                             '2(1)-5(3)-10(6)-10'=c(2,1,5,3,10,6,10,10), #P25 age-stratified
#                             '1-4-8-10'=c(1,1,4,4,8,8,10,10), #P10 Non-stratified #P5 Non-stratified
#                             '1-4(3)-8(5)-10'=c(1,1,4,3,8,5,10,10), #P10 age-stratified
#                             '1-4(2)-8(5)-10'=c(1,1,4,2,8,5,10,10)), #P5 age-stratified
#               'male'=list('5-5-5-5'=rep(5,8), # Benchmark Non-stratified
#                           '2-5-10-10'=c(2,2,5,5,10,10,10,10), #P50 Non-stratified
#                           '2(1)-5-10-10'=c(2,1,5,5,10,10,10,10), #P50 age-stratified
#                           '1-4-8-10'=c(1,1,4,4,8,8,10,10), #P25
#                           '1-3-6-10'=c(1,1,3,3,6,6,10,10))) #P10 #P5


# round down if first decimal is 0
stg_l <- list('female'=list('5-5-5-5'=rep(5,8), # Benchmark Non-stratified
                            '2-6-10-10'=c(2,2,6,6,10,10,10,10), #P50 Non-stratified
                            '3(1)-7(4)-10(6)-10'=c(3,1,7,4,10,6,10,10), #P50 age-stratified
                            '1-5-9-10'=c(1,1,5,5,9,9,10,10), #P25 Non-stratified
                            '1-5(3)-9(6)-10'=c(1,1,5,3,9,6,10,10), #P25 age-stratified
                            '1-4-8-10'=c(1,1,4,4,8,8,10,10), #P10 Non-stratified
                            '1-3-8-10'=c(1,1,3,3,8,8,10,10), #P5 Non-stratified
                            '1-4(2)-8(5)-10'=c(1,1,4,2,8,5,10,10)), #P10 age-stratified #P5 age-stratified
              'male'=list('5-5-5-5'=rep(5,8), # Benchmark Non-stratified
                          '2-5-9-10'=c(2,2,5,5,9,9,10,10), #P50 Non-stratified
                          '2(1)-5-9-10'=c(2,1,5,5,9,9,10,10), #P50 age-stratified
                          '1-4-8-10'=c(1,1,4,4,8,8,10,10), #P25
                          '1-3-6-10'=c(1,1,3,3,6,6,10,10))) #P10 #P5

str_type <- expand.grid(Percentile = c('P50','P25','P10','P5'), Stratified= c('Non-stratified', 'Age-stratified'), 
                        Sex=c('female','male'))


# Create a data frame based on stg_l (modified manually from output of ChatGPT, check!)
# strategies_df <- data.frame(
#   sex = rep(c("female", "male"), each = 9),
#   quantile = rep(c("Benchmark", "P50", "P50", "P25", "P25", "P10", "P10", "P5", "P5"),2),
#   stratified = rep(c("Non-stratified", rep(c("Non-stratified", "Age-stratified"),4)), 2),
#   Strategy = c("5-5-5-5", "2-6-10-10", "3(2)-7(4)-10(6)-10", "1-5-9-10", "2(1)-5(3)-10(6)-10", 
#                "1-4-8-10", "1-4(3)-8(5)-10","1-4-8-10", "1-4(2)-8(5)-10",
#                "5-5-5-5", "2-5-10-10", "2(1)-5-10-10", "1-4-8-10", "1-4-8-10", "1-3-6-10","1-3-6-10",
#                "1-3-6-10", "1-3-6-10"),
#   stringsAsFactors = FALSE
# )

strategies_df <- data.frame(
  sex = rep(c("female", "male"), each = 9),
  quantile = rep(c("Benchmark", "P50", "P50", "P25", "P25", "P10", "P10", "P5", "P5"),2),
  stratified = rep(c("Non-stratified", rep(c("Non-stratified", "Age-stratified"),4)), 2),
  Strategy = c("5-5-5-5", "2-6-10-10", "3(1)-7(4)-10(6)-10", "1-5-9-10", "1-5(3)-9(6)-10", 
               "1-4-8-10", "1-4(2)-8(5)-10","1-3-8-10", "1-4(2)-8(5)-10",
               "5-5-5-5", "2-5-9-10", "2(1)-5-9-10", "1-4-8-10", "1-4-8-10", "1-3-6-10","1-3-6-10",
               "1-3-6-10", "1-3-6-10"),
  stringsAsFactors = FALSE
)

strategies_df


#########################################
# Number of assessments needed per i
#########################################
nb_ass <- function(stg, gender){
  weight <- merge(riskgr_pct, age_proportion)
  weight %<>% mutate(w=riskgr_pct_total*n)
  weight %<>% mutate(age_gr = ifelse(lm_age<65, '<65', '>=65'))
  
  if(gender!='both'){
    stg_ <- tibble(age_gr=rep(c('<65', '>=65'),4), 
                   risk_class=rep(c('red', 'orange', 'yellow', 'green'), each=2),
                   int=stg)
    stg_ <- merge(stg_, weight)
    # aggregation 2 age groups
    nb_agg <- stg_ %>% filter(sex%in%gender) %>% summarise(nb = sum((10/int)*w)/sum(w))
    return(nb_agg)
  }else{
    stg_ <- tibble(age_gr=rep(c('<65', '>=65'),8), 
                   risk_class=rep(rep(c('red', 'orange', 'yellow', 'green'), each=2),2),
                   sex=rep(c('female','male'), each=8),
                   int=c(stg$female, stg$male))
    stg_ <- merge(stg_, weight)
    # aggregation 2 age groups
    nb_agg <- stg_ %>% summarise(nb = sum((10/int)*w)/sum(w))
    return(nb_agg)
  }
}
# Sex specific
lapply(stg_l$female, function(e) nb_ass(e, 'female'))
lapply(stg_l$male, function(e) nb_ass(e, 'male'))

# whole population
split_data <- split(strategies_df, strategies_df$sex)
unique_combinations <- unique(strategies_df[, c("quantile", "stratified")])
# Loop over combinations and create the list
nb_population <- lapply(1:nrow(unique_combinations), function(i) {
  combo <- unique_combinations[i, ]
  quantile <- combo$quantile
  stratified <- combo$stratified
    female_stg_name = split_data$female[split_data$female$quantile == quantile & 
                                 split_data$female$stratified == stratified, "Strategy"]
    male_stg_name = split_data$male[split_data$male$quantile == quantile & 
                             split_data$male$stratified == stratified, "Strategy"]
    
    stg_subl <- list(
      female = stg_l$female[[female_stg_name]],
      male = stg_l$male[[male_stg_name]]
    )
    nb_ass(stg_subl, gender = 'both')
})
names(nb_population) <- apply(unique_combinations, 1, function(x) paste(x, collapse = "_"))
nb_population
#########################################
# Proportion of t* < first assessment
#########################################
eff_tables <- function(stg_subl){
  prop <- tibble()
  for(k in 1:length(stg_subl)){
    gender=names(stg_subl)[k]
    stg=stg_subl[[k]]
    prop_ = tstar %>% filter(risk_class != 'darkred', sex==gender) %>% 
      group_by(lm_age, risk_class) %>% 
      summarise(mean = mean(pred_time < case_when(risk_class == 'red' & lm_age<65 ~ stg[1],
                                                  risk_class == 'red' & lm_age>=65 ~ stg[2],
                                                  risk_class == 'orange' & lm_age<65 ~ stg[3],
                                                  risk_class == 'orange' & lm_age>=65 ~ stg[4],
                                                  risk_class == 'yellow' & lm_age<65 ~ stg[5],
                                                  risk_class == 'yellow' & lm_age>=65 ~ stg[6],
                                                  risk_class == 'green' & lm_age<65 ~ stg[7],
                                                  risk_class == 'green' & lm_age>=65 ~ stg[8]))) %>% 
      ungroup %>% pivot_wider(names_from = risk_class, values_from = mean) %>% mutate(sex=gender)
    prop %<>% bind_rows(prop_)
  }
  
    # write.csv(prop_, file = paste0(out_path, 'prop_inneed_',paste0(stg, collapse = ''),'.csv'))
    prop_ <- prop %>% pivot_longer(cols = red:green, names_to = 'risk_class', values_to = 'Proportion') 
    prop_ <- merge(prop_, riskgr_pct)
    prop_ <- merge(prop_, age_proportion)
    prop_ %<>% mutate(w=riskgr_pct_total*n)
    # age range
    prop_ %<>% mutate(age_gr = ifelse(lm_age<65, '<65', '>=65'))
    
    # aggregation 2 age groups
    prop_agg <- prop_ %>% group_by(risk_class, age_gr) %>% 
      summarise(Proportion = sum(Proportion*w)/sum(w))
    # aggregation all ages
    prop_agg2 <- prop_ %>% group_by(risk_class) %>% 
      summarise(Proportion = sum(Proportion*w)/sum(w)) %>% mutate(age_gr='All')
    
    prop_agg <- rbind(prop_agg, prop_agg2)
    prop_agg$risk_class %<>% factor(levels = c('red', 'orange', 'yellow', 'green'))
    prop_agg$age_gr %<>% factor
  return(list(prop = prop, prop_agg=prop_agg))
}

# Initiation
inneed_table_female <- vector(mode = 'list')
inneed_table_male <- vector(mode = 'list')
outp_exl_female <- createWorkbook()
outp_exl_male <- createWorkbook()

for (j in 1:length(stg_l$female)) {
  eff <- eff_tables(list('female'=stg_l$female[[j]]))
  addWorksheet(outp_exl_female, names(stg_l$female)[j])
  writeData(outp_exl_female, sheet = names(stg_l$female)[j], x = eff$prop)
  inneed_table_female[[j]] <- eff$prop_agg
}
for (j in 1:length(stg_l$male)) {
  eff <- eff_tables(list('male'=stg_l$male[[j]]))
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

inneed_table_female <- strategies_df %>% filter(sex=='female') %>% select(-sex) %>% merge(inneed_table_female)
inneed_table_male <- strategies_df %>% filter(sex=='male') %>% select(-sex) %>% merge(inneed_table_male)

inneed_table_female$quantile %<>% factor(levels = c('Benchmark','P50','P25','P10','P5'))
inneed_table_male$quantile %<>% factor(levels = c('Benchmark','P50','P25','P10','P5'))
inneed_table_female$stratified %<>% factor(levels = c('Non-stratified','Age-stratified'))
inneed_table_male$stratified %<>% factor(levels = c('Non-stratified','Age-stratified'))
inneed_table_female$Strategy %<>% factor(levels = names(stg_l$female))
inneed_table_male$Strategy %<>% factor(levels = names(stg_l$male))


p_inneed_male <- ggplot(inneed_table_male %>% filter(quantile!='P5'), aes(x=risk_class, y=Proportion, color=age_gr, fill=age_gr))+
  geom_bar(stat='identity', position = 'dodge', width=.6, alpha=0.4)+
  # geom_segment(data=filter(inneed_table_male, Strategy=='5-5-5-5'),
  #              aes(x=as.numeric(risk_class)-0.25, xend=as.numeric(risk_class)+0.25, y=Proportion, yend=Proportion))+
  facet_grid(stratified~quantile)+
  ylab('Proportion')+
  xlab('Risk group')+
  ggtitle('Men')+
  theme_bw()+
  labs(fill="Age", color='Age')+
  scale_x_discrete(labels=c("red" = "7.5%-10%", "orange" = "5%-7.5%",
                            "yellow" = "2.5%-5%", "green" = "<2.5%"))+
  scale_color_manual(values = wes_palette("Moonrise3", n = 3))+
  scale_fill_manual(values = wes_palette("Moonrise3", n = 3))+
  geom_text(aes(x = 3, y = 1, label = paste0(Strategy)), size = 3, color='black') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p_inneed_female <- ggplot(inneed_table_female %>% filter(quantile!='P5'), aes(x=risk_class, y=Proportion, color=age_gr, fill=age_gr))+
  geom_bar(stat='identity', position = 'dodge', width=.6, alpha=0.4)+
  # geom_segment(data=filter(inneed_table_female, Strategy=='5-5-5-5'),
  #              aes(x=as.numeric(risk_class)-0.25, xend=as.numeric(risk_class)+0.25, y=Proportion, yend=Proportion))+
  facet_grid(stratified~quantile)+
  ylab('Proportion')+
  xlab('Risk group')+
  ggtitle('Women')+
  theme_bw()+
  scale_x_discrete(labels=c("red" = "7.5%-10%", "orange" = "5%-7.5%",
                            "yellow" = "2.5%-5%", "green" = "<2.5%"))+
  labs(fill="Age", color='Age')+
  scale_color_manual(values = wes_palette("Moonrise3", n = 3))+
  scale_fill_manual(values = wes_palette("Moonrise3", n = 3))+
  geom_text(aes(x = 3, y = 1, label = paste0(Strategy)), size = 3, color='black') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# mylegend<-g_legend(p_inneed_female+ labs(fill="Age", color='Age'))

jpeg(file = '~/epi_paper/outp_figs/inneed.jpeg', units="in", width=8.5, height=10, res=600)
grid.arrange(p_inneed_female, p_inneed_male, nrow=2, heights=c(1,1)
             # top=textGrob('Proportion of individuals crossed at the first assessment', gp=gpar(fontsize=18,font=8))
             )
# grid.arrange(arrangeGrob(p_inneed_female + theme(legend.position="none"),
#                          p_inneed_male + theme(legend.position="none"), nrow=2, heights=c(1,1)),
#              mylegend, ncol=1, heights=c(10,1), top=textGrob('Proportion of t*<t1', gp=gpar(fontsize=20,font=8)))
dev.off()

####################
# # P50 only
# inneed_table_P50 <- bind_rows(inneed_table_female, inneed_table_male)
# 
# inneed_table_P50 <- inneed_table_P50 %>% filter(quantile %in% c('Benchmark','P50')) %>% 
#   mutate(quantile=ifelse(quantile=='P50', 
#                          ifelse(stratified=='Non-stratified', 'P50 (non-stratified)', 'P50 (age-stratified)'),
#                          'Benchmark'))
# inneed_table_P50$quantile %<>% factor(levels = c('Benchmark', 'P50 (non-stratified)', 'P50 (age-stratified)'))
# inneed_table_P50$sex %<>% factor(levels = c('female', 'male'), labels = c('Women','Men'))
# 
# p_inneed_P50 <- ggplot(inneed_table_P50, aes(x=risk_class, y=Proportion, color=age_gr, fill=age_gr))+
#   geom_bar(stat='identity', position = 'dodge', width=.6, alpha=0.4)+
#   # geom_segment(data=filter(inneed_table_male, Strategy=='5-5-5-5'),
#   #              aes(x=as.numeric(risk_class)-0.25, xend=as.numeric(risk_class)+0.25, y=Proportion, yend=Proportion))+
#   facet_grid(sex~quantile)+
#   ylab('Proportion')+
#   xlab('Risk group')+
#   theme_bw()+
#   labs(fill="Age", color='Age')+
#   scale_x_discrete(labels=c("red" = "7.5%-10%", "orange" = "5%-7.5%",
#                             "yellow" = "2.5%-5%", "green" = "<2.5%"))+
#   scale_color_manual(values = wes_palette("Moonrise3", n = 3))+
#   scale_fill_manual(values = wes_palette("Moonrise3", n = 3))+
#   geom_text(aes(x = 3, y = 1, label = paste0(Strategy)), size = 3, color='black') +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# jpeg(file = '~/epi_paper/outp_figs/inneed_P50.jpeg', units="in", width=6, height=5, res=600)
# p_inneed_P50
# dev.off()

#aggregate across risk gr
inneed_table_female %<>% mutate(sex='female')
inneed_table_male %<>% mutate(sex='male')
inneed_table_female_agg <- merge(inneed_table_female, weight)
inneed_table_male_agg <- merge(inneed_table_male, weight)
inneed_table_agg <- rbind(inneed_table_female_agg, inneed_table_male_agg)

prop_outp <- lapply(1:nrow(unique_combinations), function(i) {
  combo <- unique_combinations[i, ]
  each_sex = inneed_table_agg %>% filter(age_gr=='All', stratified==combo$stratified, quantile==combo$quantile) %>% group_by(sex)%>% 
    summarise(prop_agg = sum(Proportion*N_riskgr)/sum(N_riskgr))
  population = inneed_table_agg %>% filter(age_gr=='All', stratified==combo$stratified, quantile==combo$quantile) %>% 
    summarise(prop_agg = sum(Proportion*N_riskgr)/sum(N_riskgr))
  bind_rows(each_sex, tibble(sex='both', prop_agg=as.numeric(population)))
})
names(prop_outp) <- apply(unique_combinations, 1, function(x) paste(x, collapse = "_"))
prop_outp

#########################################
# Avg time of staying in high risk: mean (t1-t* | t*<t1)
#########################################
tstar_undetect <- function(stg, gender){
  tstar_inneed <- tstar %>% filter(risk_class != 'darkred', sex==gender) %>% 
    filter(pred_time < case_when(risk_class == 'red' & lm_age<65 ~ stg[1],
                              risk_class == 'red' & lm_age>=65 ~ stg[2],
                              risk_class == 'orange' & lm_age<65 ~ stg[3],
                              risk_class == 'orange' & lm_age>=65 ~ stg[4],
                              risk_class == 'yellow' & lm_age<65 ~ stg[5],
                              risk_class == 'yellow' & lm_age>=65 ~ stg[6],
                              risk_class == 'green' & lm_age<65 ~ stg[7],
                              risk_class == 'green' & lm_age>=65 ~ stg[8]))
  
  tstar_undet = tstar_inneed %>% group_by(sex, lm_age, risk_class) %>% 
    mutate(undetect = case_when(risk_class == 'red' & lm_age<65 ~ stg[1],
                                risk_class == 'red' & lm_age>=65 ~ stg[2],
                                risk_class == 'orange' & lm_age<65 ~ stg[3],
                                risk_class == 'orange' & lm_age>=65 ~ stg[4],
                                risk_class == 'yellow' & lm_age<65 ~ stg[5],
                                risk_class == 'yellow' & lm_age>=65 ~ stg[6],
                                risk_class == 'green' & lm_age<65 ~ stg[7],
                                risk_class == 'green' & lm_age>=65 ~ stg[8])-pred_time) %>% 
    summarise(mean = mean(undetect)) %>% ungroup %>% 
    pivot_wider(names_from = risk_class, values_from = mean) 
  
  # Find missing columns
  required_cols <- c("red", "orange", "yellow", "green")
  missing_cols <- setdiff(required_cols, names(tstar_undet))
  tstar_undet %<>% mutate(!!!set_names(rep(list(NA_real_), length(missing_cols)), missing_cols))
  return(tstar_undet)
}

wtime_tables <- function(stg_name_l){
  wtime <- tibble()
  cross_pct <- tibble()
  for(k in 1:length(stg_name_l)){
    gender=names(stg_name_l)[k]
    stg_name=stg_name_l[[k]]
    stg <- stg_l[[gender]][[stg_name]]
    wtime_ <- tstar_undetect(stg, gender)
    wtime %<>% bind_rows(wtime_)
    cross_pct_ <- read.xlsx(paste0(out_path, 'prop_inneed_', gender, '.xlsx'), sheet = stg_name) %>% 
      pivot_longer(cols = red:green, names_to = 'risk_class', values_to = 'crossed')
    cross_pct %<>% bind_rows(cross_pct_)
  }
  # write.csv(wtime_, file = paste0(out_path, 'wtime_inneed_',paste0(stg, collapse = ''),'.csv'))
  wtime_ <- wtime %>% pivot_longer(cols = red:green, names_to = 'risk_class', values_to = 'wait_time') 
  wtime_ <- merge(wtime_, riskgr_pct)
  wtime_ <- merge(wtime_, age_proportion)
  wtime_ <- merge(wtime_, cross_pct)
  wtime_ %<>% mutate(w=n*riskgr_pct_total*crossed)
  # age range
  wtime_ %<>% mutate(age_gr = ifelse(lm_age<65, '<65', '>=65'))
  wtime_[is.na(wtime_)] <- 0
  
  # aggregation 2 age groups
  wtime_agg <- wtime_ %>% group_by(risk_class, age_gr) %>% summarise(wait_time = sum(wait_time*w)/sum(w))
  # aggregation all ages
  wtime_agg2 <- wtime_ %>% group_by(risk_class) %>% summarise(wait_time = sum(wait_time*w)/sum(w)) %>% 
    mutate(age_gr='All')
  
  wtime_agg <- rbind(wtime_agg, wtime_agg2)
  wtime_agg$risk_class %<>% factor(levels = c('red', 'orange', 'yellow', 'green'))
  wtime_agg$age_gr %<>% factor
  return(list(wtime = wtime, wtime_agg=wtime_agg))
}

# Initiation
wtime_table_female <- vector(mode = 'list')
wtime_table_male <- vector(mode = 'list')
outp_wt_exl_female <- createWorkbook()
outp_wt_exl_male <- createWorkbook()

for (j in 1:length(stg_l$female)) {
  eff <- wtime_tables(list('female'=names(stg_l$female)[j]))
  addWorksheet(outp_wt_exl_female, names(stg_l$female)[j])
  writeData(outp_wt_exl_female, sheet = names(stg_l$female)[j], x = eff$wtime)
  wtime_table_female[[j]] <- eff$wtime_agg
}
for (j in 1:length(stg_l$male)){
  eff <- wtime_tables(list('male'=names(stg_l$male)[j]))
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

wtime_table_female <- strategies_df %>% filter(sex=='female') %>% select(-sex) %>% merge(wtime_table_female)
wtime_table_male <- strategies_df %>% filter(sex=='male') %>% select(-sex) %>% merge(wtime_table_male)

wtime_table_female$quantile %<>% factor(levels = c('Benchmark','P50','P25','P10','P5'))
wtime_table_male$quantile %<>% factor(levels = c('Benchmark','P50','P25','P10','P5'))
wtime_table_female$stratified %<>% factor(levels = c('Non-stratified','Age-stratified'))
wtime_table_male$stratified %<>% factor(levels = c('Non-stratified','Age-stratified'))
wtime_table_female$Strategy %<>% factor(levels = names(stg_l$female))
wtime_table_male$Strategy %<>% factor(levels = names(stg_l$male))

p_wtime_male <- ggplot(wtime_table_male %>% filter(quantile!='P5'), aes(x=risk_class, y=wait_time, color=age_gr, fill=age_gr))+
  geom_bar(stat='identity', position = 'dodge', width=.6, alpha=0.4)+
  # geom_point(aes(shape = Strategy), size=3, alpha=0.8, position=position_jitter(h=0,w=0.2))+
  # scale_shape_manual(values=c(21:24))+
  # geom_segment(data=filter(wtime_table_male, Strategy=='5-5-5-5'),
  #              aes(x=as.numeric(risk_class)-0.25, xend=as.numeric(risk_class)+0.25, y=wtimeortion, yend=wtimeortion))+
  facet_grid(stratified~quantile)+
  ylab('Waiting time (years)')+
  xlab('Risk group')+
  ggtitle('Men')+
  theme_bw()+
  labs(fill="Age", color='Age')+
  scale_x_discrete(labels=c("red" = "7.5%-10%", "orange" = "5%-7.5%",
                            "yellow" = "2.5%-5%", "green" = "<2.5%"))+
  scale_color_manual(values = wes_palette("Moonrise3", n = 3))+
  scale_fill_manual(values = wes_palette("Moonrise3", n = 3))+
  geom_text(aes(x = 3, y = Inf, label = paste0(Strategy)), size = 3, vjust =1.5, color='black') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p_wtime_female <- ggplot(wtime_table_female %>% filter(quantile!='P5'), aes(x=risk_class, y=wait_time, color=age_gr, fill=age_gr))+
  geom_bar(stat='identity', position = 'dodge', width=.6, alpha=0.4)+
  # scale_shape_manual(values=c(21:24))+
  # geom_segment(data=filter(wtime_table_female, Strategy=='5-5-5-5'),
  #              aes(x=as.numeric(risk_class)-0.25, xend=as.numeric(risk_class)+0.25, y=wtimeortion, yend=wtimeortion))+
  facet_grid(stratified~quantile)+
  ylab('Waiting time (years)')+
  xlab('Risk group')+
  ggtitle('Women')+
  theme_bw()+
  labs(fill="Age", color='Age')+
  scale_x_discrete(labels=c("red" = "7.5%-10%", "orange" = "5%-7.5%",
                            "yellow" = "2.5%-5%", "green" = "<2.5%"))+
  scale_color_manual(values = wes_palette("Moonrise3", n = 3))+
  scale_fill_manual(values = wes_palette("Moonrise3", n = 3))+
  geom_text(aes(x = 3, y = Inf, label = paste0(Strategy)), size = 3, vjust =1.5, color='black') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# mylegend<-g_legend(p_wtime_female+ labs(fill="Age", color='Age'))

jpeg(file = '~/epi_paper/outp_figs/wtime.jpeg', units="in", width=8.5, height=10, res=600)
grid.arrange(p_wtime_female, p_wtime_male, nrow=2, heights=c(1,1)
             # top=textGrob('Average waiting time after crossing', gp=gpar(fontsize=18,font=8))
             )
# grid.arrange(arrangeGrob(p_wtime_female + theme(legend.position="none"),
#                          p_wtime_male + theme(legend.position="none"), nrow=1, widths=c(1,1)),
#              mylegend, ncol=2, widths=c(8,1))
dev.off()

####################
# # P50 only
# wtime_table_P50 <- bind_rows(wtime_table_female, wtime_table_male)
# 
# wtime_table_P50 <- wtime_table_P50 %>% filter(quantile %in% c('Benchmark','P50')) %>% 
#   mutate(quantile=ifelse(quantile=='P50', 
#                          ifelse(stratified=='Non-stratified', 'P50 (non-stratified)', 'P50 (age-stratified)'),
#                          'Benchmark'))
# wtime_table_P50$quantile %<>% factor(levels = c('Benchmark', 'P50 (non-stratified)', 'P50 (age-stratified)'))
# wtime_table_P50$sex %<>% factor(levels = c('female', 'male'), labels = c('Women','Men'))
# 
# p_wtime_P50 <- ggplot(wtime_table_P50, aes(x=risk_class, y=wait_time, color=age_gr, fill=age_gr))+
#   geom_bar(stat='identity', position = 'dodge', width=.6, alpha=0.4)+
#   # scale_shape_manual(values=c(21:24))+
#   # geom_segment(data=filter(wtime_table_female, Strategy=='5-5-5-5'),
#   #              aes(x=as.numeric(risk_class)-0.25, xend=as.numeric(risk_class)+0.25, y=wtimeortion, yend=wtimeortion))+
#   facet_grid(sex~quantile)+
#   ylab('Waiting time (years)')+
#   xlab('Risk group')+
#   theme_bw()+
#   labs(fill="Age", color='Age')+
#   scale_x_discrete(labels=c("red" = "7.5%-10%", "orange" = "5%-7.5%",
#                             "yellow" = "2.5%-5%", "green" = "<2.5%"))+
#   scale_color_manual(values = wes_palette("Moonrise3", n = 3))+
#   scale_fill_manual(values = wes_palette("Moonrise3", n = 3))+
#   new_scale_color()+new_scale_fill()+
#   geom_text(aes(x = 3, y = Inf, label = paste0(Strategy)), size = 3, vjust =1.5) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# jpeg(file = '~/epi_paper/outp_figs/wtime_P50.jpeg', units="in", width=6, height=5, res=600)
# p_wtime_P50
# dev.off()

#aggregate across risk gr
n_inneed <- inneed_table_agg %>% mutate(N_inneed=Proportion*N_riskgr) %>% select(-Proportion, -N_riskgr)
wtime_table_female %<>% mutate(sex='female')
wtime_table_male %<>% mutate(sex='male')
wtime_table <- rbind(wtime_table_female, wtime_table_male)
wtime_table_agg <- merge(wtime_table, n_inneed)

wtime_outp <- lapply(1:nrow(unique_combinations), function(i) {
  combo <- unique_combinations[i, ]
  each_sex = wtime_table_agg %>% filter(age_gr=='All', stratified==combo$stratified, quantile==combo$quantile) %>% 
    group_by(sex) %>% summarise(wtime_agg = sum(wait_time*N_inneed, na.rm = T)/sum(N_inneed))
  population = wtime_table_agg %>% filter(age_gr=='All', stratified==combo$stratified, quantile==combo$quantile) %>% 
    summarise(wtime_agg = sum(wait_time*N_inneed, na.rm = T)/sum(N_inneed))
  bind_rows(each_sex, tibble(sex='both', wtime_agg=as.numeric(population)))
})
names(wtime_outp) <- apply(unique_combinations, 1, function(x) paste(x, collapse = "_"))
wtime_outp

#########################################
# Proportion of T_wait > 1
#########################################
# tstar_uncovered <- function(stg){
#   # Conditional on t* <= first assessment
#   tstar_inneed <- tstar %>% filter(risk_class != 'darkred') %>% 
#     mutate(pred_time = ifelse(pred_status==1, pred_time, 8)) %>% 
#     filter(pred_time <= case_when(risk_class == 'green' & lm_age<65 ~ stg[1],
#                               risk_class == 'green' & lm_age>=65 ~ stg[2],
#                               risk_class == 'yellow' & lm_age<65 ~ stg[3],
#                               risk_class == 'yellow' & lm_age>=65 ~ stg[4],
#                               risk_class == 'orange' & lm_age<65 ~ stg[5],
#                               risk_class == 'orange' & lm_age>=65 ~ stg[6],
#                               risk_class == 'red' & lm_age<65 ~ stg[7],
#                               risk_class == 'red' & lm_age>=65 ~ stg[8]))
#   # window of 1 year
#   stg = stg -1
#   tstar_uncover = tstar_inneed %>% group_by(sex, lm_age, risk_class) %>% 
#     summarise(mean = 1-mean(ifelse(pred_time < case_when(risk_class == 'green' & lm_age<65 ~ stg[1],
#                                                      risk_class == 'green' & lm_age>=65 ~ stg[2],
#                                                      risk_class == 'yellow' & lm_age<65 ~ stg[3],
#                                                      risk_class == 'yellow' & lm_age>=65 ~ stg[4],
#                                                      risk_class == 'orange' & lm_age<65 ~ stg[5],
#                                                      risk_class == 'orange' & lm_age>=65 ~ stg[6],
#                                                      risk_class == 'red' & lm_age<65 ~ stg[7],
#                                                      risk_class == 'red' & lm_age>=65 ~ stg[8]), 0, 1))) %>% 
#     ungroup %>% pivot_wider(names_from = risk_class, values_from = mean) %>%
#     select(sex, lm_age, Low, Medium_low, Medium_high, High)
#   return(tstar_uncover)
# }
# 






####################################################################################################
## Explorative below

t_summary_age_risk <- tstar %>% filter(risk_class!='darkred') %>% 
  mutate(pred_time=ifelse(pred_status==0, 8, pred_time)) %>% group_by(lm_age, risk_class, sex) %>% 
  summarise(mean=mean(pred_time), median=median(pred_time), Q1=quantile(pred_time,0.25), Q3=quantile(pred_time,0.75))

weight <- merge(N_cross_riskgr, age_proportion, by=c('lm_age', 'sex'))
weight <- weight %>% mutate(w=w*riskgr_pct_total) %>% select(lm_age,sex,risk_class,w) %>% 
  mutate(risk_class= as.factor(risk_class)) %>% distinct %>% as_tibble()

t_summary_risk <- merge(t_summary_age_risk, weight, by=c('lm_age', 'sex', 'risk_class'))
t_summary_risk <- t_summary_risk %>% group_by(sex, risk_class) %>% mutate(w=w/sum(w)) %>% 
  summarise(mean=sum(mean*w),median=sum(median*w), Q1=sum(Q1*w), Q3=sum(Q3*w))

levels(t_summary_risk$risk_class) <- c('darkred(>5%)','High(3.75%-5%)','Medium_high(2.5%-3.75%)',
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

tstar_overlap_40 %<>% mutate(cross_gr = cut(pred_time, breaks=c(-Inf, 0, 5, 8, Inf), 
                                            labels=c("already","0-5","5-8",'>8y')))
tstar_overlap_40 %<>% mutate(cross_gr = replace(cross_gr, pred_status == 0, '>8y'))

tstar_overlap_45 %<>% mutate(cross_gr = cut(pred_time, breaks=c(-Inf, 0, 3, Inf), 
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
