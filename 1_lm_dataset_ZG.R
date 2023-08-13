#creating derivation landmark datasets 

# #print(commandArgs(trailingOnly = TRUE))

args=commandArgs(trailingOnly = TRUE)
j = as.integer(args[1])
gender = as.character(args[2])

out_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/'
in_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/'
path_to_save = paste0(out_path, "data_created/", gender,"/")


library(data.table)
library(nlme)
library(survival)
library(tidyverse)
library(magrittr)
library(haven)

load(paste0(in_path,'outcomes_new.RData'))
load(paste0(in_path,'exposures_merged_new_popmean.RData'))
# align exposure and outcome patid
patid_unique <- unique(exposures_merged$patid)
stopifnot(all(patid_unique %in% outcomes$patid))
outcomes %<>% filter(patid %in% patid_unique)

family_data <- read_dta(paste0(in_path,'familyhistory.dta'))

#selecting correct patients 
cat(j,"\n")
#j = 80

## Descriptives
# patid_unique <- unique(exposures_merged$patid)
# length(patid_unique)
# length(unique(exposures_merged$patid[exposures_merged$gender == "Female"]))
# length(unique(exposures_merged$patid[exposures_merged$gender == "Male"]))
# length(unique(outcomes$pracid))
# length(unique(outcomes$region))
# length(unique(outcomes$pracid[outcomes$derivation == 'derivation']))
# length(unique(outcomes$pracid[outcomes$derivation == 'validation']))
# length(unique(outcomes$region[outcomes$derivation == 0]))

# length(which(outcomes$derivation == 'derivation'))
# length(which(outcomes$derivation == 'validation'))
# dist_fup = as.numeric((outcomes$end_date - outcomes$start_date)/365.25)
# summary(dist_fup)

#Filter on patients


#Focus on end, start, cvd, statins
end_age_cond     = ( outcomes$end_date - outcomes$d_yob )/365.25 > j
statin_age_cond  = ( outcomes$statins_prscd - outcomes$d_yob )/365.25 > j | is.na( outcomes$statins_prscd )
start_age_cond   = ( outcomes$start_date - outcomes$d_yob )/365.25 <= j
cvd_age_cond     = ( outcomes$cvd_date - outcomes$d_yob )/365.25 > j | is.na( outcomes$cvd_date )

if(gender == "male"){
  gender_cond      = outcomes$gender == 1 | outcomes$gender == "Male"
}else{
  gender_cond      = outcomes$gender == 0 | outcomes$gender == "Female"
}

selected_id = end_age_cond & statin_age_cond & start_age_cond & cvd_age_cond & gender_cond
selected_pt = outcomes$patid[ selected_id ]


#check on correctness CVD, DEATH and END times
colnames(outcomes)
# outcomes[, c("cvd_date", "death_date", "end_date")]

print("People that died with CVD (during or after diagnosis)")
length( which( outcomes$cvd_ind == 1 & outcomes$death_ind == 1 ) )/length(outcomes$cvd_ind)
print("People that died without a diagnosis of CVD")
length( which( outcomes$cvd_ind == 0 & outcomes$death_ind == 1 ) )/length(outcomes$cvd_ind)
print("People censored with a diagnosis of CVD")
length( which( outcomes$cvd_ind == 1 & outcomes$death_ind == 0 ) )/length(outcomes$cvd_ind)
print("People censored without a diagnosis of CVD")
length( which( outcomes$cvd_ind == 0 & outcomes$death_ind == 0 ) )/length(outcomes$cvd_ind)


outcomes_selected = outcomes[ selected_id, ]

#Keep only the exposures of those patients that satisfy all the requirements
exposures_selected = exposures_merged[ exposures_merged$patid %in% selected_pt, ]

rm( exposures_merged, outcomes, selected_id, selected_pt, end_age_cond,
    cvd_age_cond, gender_cond, start_age_cond,
    statin_age_cond)

#adding index for exposures in the future (after L_a)
exposures_selected$future = ifelse( exposures_selected$exp_age >=j, 1, 0 )
#adding age corrected (centralised)
exposures_selected$exp_age_corr = exposures_selected$exp_age - j

#Creating the final data: we can have people with NO exposures
data_ext = merge( exposures_selected, outcomes_selected, all = F )

dim( data_ext )

#statin_index
data_ext$statin_bin = ifelse( data_ext$statins_prscd <= data_ext$exp_date & !is.na(data_ext$statins_prscd), 1, 0 ) 

#bp index
data_ext$bp_bin = ifelse( data_ext$bp_med_prscd <= data_ext$exp_date & !is.na(data_ext$bp_med_prscd), 1, 0 ) 

#sbp measure ind
data_ext$sbp_ind = ifelse( data_ext$exposure == "sbp", 1, 0 ) 

#tchol measure ind
data_ext$tchol_ind = ifelse( data_ext$exposure == "tchol", 1, 0 ) 


print("almost creating data_table")
data_ext = data.table( data_ext )
print("created data_table")

setkey( data_ext, patid )

data_ext[ ,exp_count := 1:.N, by = patid ]
print("created exp_count")

#data_ext[ death_ind == 1, ]
#data_ext[ death_ind == 1 & death_date != end_date, lapply(.SD, function(x) length(unique(x))), .SDcols = "patid" ]
#2902

#warning decison!!!!!!!
#data_ext[ death_ind == 1 & death_date != end_date, death_ind := 0 ] #censoring death
# table( data_ext$SLE_ind )

stopifnot( length( which( data_ext$SLE_ind == 1 ) ) ==  length( which( !is.na( data_ext$SLE_date ) ) ) )
stopifnot( length( which( data_ext$Migraine_ind == 1 ) ) ==  length( which( !is.na( data_ext$Migraine_date ) ) ) )
stopifnot( length( which( data_ext$depression_ind == 1 ) ) ==  length( which( !is.na( data_ext$depression_date ) ) ) )
stopifnot( length( which( data_ext$Dementia_ind == 1 ) ) ==  length( which( !is.na( data_ext$Dementia_date ) ) ) )
stopifnot( length( which( data_ext$Severe_mental_illness_ind == 1 ) ) ==  length( which( !is.na( data_ext$Severe_mental_illness_date ) ) ) )
stopifnot( length( which( data_ext$diab_ind == 1 ) ) ==  length( which( !is.na( data_ext$diab_date ) ) ) )
stopifnot( length( which( data_ext$Atrial_fibrillation_ind == 1 ) ) ==  length( which( !is.na( data_ext$Atrial_fibrillation_date ) ) ) )
stopifnot( length( which( data_ext$rheumatoid_arthritis_ind == 1 ) ) ==  length( which( !is.na( data_ext$rheumatoid_arthritis_date ) ) ) )
stopifnot( length( which( data_ext$renal_ind == 1 ) ) ==  length( which( !is.na( data_ext$renal_date ) ) ) )
stopifnot( length( which( data_ext$death_ind == 1 ) ) ==  length( which( !is.na( data_ext$death_date ) ) ) )
stopifnot( length( which( data_ext$cvd_ind == 1 ) ) ==  length( which( !is.na( data_ext$cvd_date ) ) ) )


#creating age variables
data_ext[ ,statin_age := ifelse( is.na( statins_prscd ), NA, ( statins_prscd - d_yob )/365.25 ) ]
data_ext[ ,end_age := ifelse( is.na( end_date ), NA, ( end_date - d_yob )/365.25 ) ]
data_ext[ ,cvd_age := ifelse( is.na( cvd_date ), NA, ( cvd_date - d_yob )/365.25 ) ]
data_ext[ ,death_age := ifelse( is.na( death_date ), NA, ( death_date - d_yob )/365.25 ) ]
data_ext[ ,diab_age := ifelse( is.na( diab_date ), NA, ( diab_date - d_yob )/365.25 ) ]
data_ext[ ,bp_med_age := ifelse( is.na( bp_med_prscd ), NA, ( bp_med_prscd - d_yob )/365.25 ) ]
data_ext[ ,renal_age := ifelse( is.na( renal_date ), NA, ( renal_date - d_yob )/365.25 ) ]
data_ext[ ,rheumatoid_arthritis_age := ifelse( is.na( rheumatoid_arthritis_date ), NA, ( rheumatoid_arthritis_date - d_yob )/365.25 ) ]
data_ext[ ,Atrial_fibrillation_age := ifelse( is.na( Atrial_fibrillation_date ), NA, ( Atrial_fibrillation_date - d_yob )/365.25 ) ]
data_ext[ ,Severe_mental_illness_age := ifelse( is.na( Severe_mental_illness_date ), NA, ( Severe_mental_illness_date - d_yob )/365.25 ) ]
data_ext[ ,Migraine_age := ifelse( is.na( Migraine_date ), NA, ( Migraine_date - d_yob )/365.25 ) ]
data_ext[ ,depression_age := ifelse( is.na( depression_date ), NA, ( depression_date - d_yob )/365.25 ) ]

#imputing !!!!!!!!!!!!!!!!!!!!!!!!
data_ext$Townsend[ which( is.na( data_ext$Townsend ) ) ] = median( data_ext$Townsend, na.rm = T )

#family hisotry 
data_ext[ ,family_history := 0 ]
CHD_patid = family_data$patid[grep( "FH: premature coronary heart disease", family_data$desc ) ]
data_ext[ patid %in% CHD_patid, family_history := 1 ]

stopifnot( length( which( data_ext$exp_age > data_ext$end_age ) ) == 0 )

save( data_ext, file =  paste0( path_to_save, "data_lm_",j, "_", gender,".RData" ) )
print('1_lm_dataset.R completed')
