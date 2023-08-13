library(haven)
library(dplyr)
library(magrittr)
library(data.table)

in_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/dataprep/cprd_2020/'
out_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/'

exp_sbp <- read_dta(paste0(in_path,'exposures_sbp_new.dta'))
exp_chol <- read_dta(paste0(in_path,'exposures_cholesterol_new.dta'))
exp_smk <- read_dta(paste0(in_path,'exposures_smoking_new.dta'))
exp_bmi <- read_dta(paste0(in_path,'exposures_clinical_raw.dta')) # this needs to select BMI and preprocessing for outliers

##########################################
# Compute BMI and remove outliers
##########################################
# select weight and height for BMI
exp_bmi <- exp_bmi %>% group_by(patid) %>%
  filter(all(c("height", "weight") %in% exposure)) %>% filter(exposure %in% c("height", "weight"))
exp_bmi <- exp_bmi %>% group_by(patid) %>% mutate(bmi = value/(first(value[exposure == "height"]))^2) %>% 
  filter(exposure == 'weight') %>% select(patid, exp_date, bmi)
colnames(exp_bmi)[3] <- 'value'
exp_bmi %<>% mutate(exposure='bmi')

# outliers
exp_sbp %<>% filter(value>=60, value<=250)
exp_tchol <- exp_chol %>% filter(exposure == "tchol", value>= 1.75, value <= 20)
exp_hdl <- exp_chol %>% filter(exposure == "hdl", value>= 0.3, value <= 3.1)
exp_bmi %<>% filter(value<=80)
rm(exp_chol)
##########################################
# Merge, remove duplicated and average conflicts
##########################################
exposures <- rbind(exp_bmi, exp_sbp, exp_hdl, exp_tchol, exp_smk)
exposures <- as.data.table(exposures)
colnames(exposures)[3] <- 'original'
# remove NAs
exposures <- exposures[!is.na(exposures$original),]
rm(exp_sbp, exp_hdl, exp_tchol, exp_smk, exp_bmi); gc()
# Duplicated and conflicted measurements
# remove duplicated
dup <- duplicated(exposures)
sum(dup)
exposures$exposure[dup] %>% table
exposures <- exposures[!dup,]
# take average conflicted
exposures <- exposures[, original := mean(original), by = .(patid, exp_date, exposure)][, .SD[1], by = .(patid, exp_date, exposure)]
exposures %<>% dplyr::filter(exposure!='smokbin' | original!=0.5)
save(exposures, file = (paste0(out_path,'exposures_new.RData')))
