in_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/3_cox_outp/'

mod_type="RIRS_"

sample_size <- tibble()
for(gender in c('female','male')){
  for(j in seq(40,80,by=5)){
    for(l in 0:10){
      load(paste0(out_path,"cox_model_", mod_type, j, "_", l,"_", gender, ".RData"))
      sample_size %<>% bind_rows(tibble(gender=gender, lm_age=j, t_interest=l, n=model_cox_5CVD$n))
    }
  }
}

sample_size %>% print(n=Inf)

