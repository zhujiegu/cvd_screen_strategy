in_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/'

library(data.table)
library(tidyverse)
library(magrittr)
library(parallel)

load(paste0(in_path,'patients_outcomes.RData'))

outcomes <- outcomes %>% rowwise %>% mutate(end_date_cor = min(end_date, cvd_date, death_date, na.rm = T))

# ages distribution of population at a specific date which does not consider statin
age_dist <- function(date){
  age_vec <- ((date - outcomes$d_yob)/365.25) %>% as.numeric()
  include = outcomes$end_date >= date & age_vec>=40 & age_vec<85 & outcomes$start_date <= date
  return(age_vec[include])
}

# ages of person in the subcohort which considers statin
age_nostatin_dist <- function(date){
  age_vec <- ((date - outcomes$d_yob)/365.25) %>% as.numeric()
  include = outcomes$end_date >= date & age_vec>=40 & age_vec<85 & 
    (is.na(outcomes$statins_prscd) | outcomes$statins_prscd>date) & outcomes$start_date <= date
  return(age_vec[include])
}

dates <- paste0(seq(2004,2017, by=4), '-4-1') %>% as.Date

ages_list <- lapply(dates, age_dist)
ages_nostatin_list <- lapply(dates, age_nostatin_dist)
names(ages_list) <- as.character(seq(2004,2017, by=4))
names(ages_nostatin_list) <- as.character(seq(2004,2017, by=4))

jpeg(file = '~/epi_paper/outp_figs/age_snapshot_marginal.jpeg', units="in", width=12, height=4, res=300)
par(mfrow=c(1,4))
hist(ages_list$`2004`, main = '2004', ylim = c(0, 90000))
hist(ages_list$`2008`, main = '2008', ylim = c(0, 90000))
hist(ages_list$`2012`, main = '2012', ylim = c(0, 90000))
hist(ages_list$`2016`, main = '2016', ylim = c(0, 90000))
par(mfrow=c(1,1))
dev.off()


jpeg(file = '~/epi_paper/outp_figs/age_snapshot.jpeg', units="in", width=12, height=4, res=300)
par(mfrow=c(1,4))
hist(ages_nostatin_list$`2004`, main = '2004', ylim = c(0, 90000))
hist(ages_nostatin_list$`2008`, main = '2008', ylim = c(0, 90000))
hist(ages_nostatin_list$`2012`, main = '2012', ylim = c(0, 90000))
hist(ages_nostatin_list$`2016`, main = '2016', ylim = c(0, 90000))
par(mfrow=c(1,1))
dev.off()

ages_nostatin_long <- stack(ages_nostatin_list)
colnames(ages_nostatin_long) <- c('age', 'year')

jpeg(file = '~/epi_paper/outp_figs/age_snapshot_density.jpeg', units="in", width=6, height=4, res=300)
ggplot(ages_nostatin_long, aes(x=age, color=year)) +
  geom_density()
dev.off()

######################################################
# lm dataset
N_female <- c(155497,159975,144907,	123135,	109173,	86403,	63257,	46787,	37178)
N_male <- c(131548,	141304,	128547,	106788,	92200,	67790,	45512,	30978,	23189)
N_lm <- N_female + N_male
names(N_lm) <- seq(40,80, by=5)
barplot(N_lm)

# # check relationship with snapshots
# dates_check <- paste0(seq(2004,2017, by=1), '-4-1') %>% as.Date
# 
# ages_nostatin_list_check <- lapply(dates_check, age_nostatin_dist)
# N_40 <- sapply(ages_nostatin_list_check, function(e) sum(e>=50 & e<51))
# sum(N_40)
