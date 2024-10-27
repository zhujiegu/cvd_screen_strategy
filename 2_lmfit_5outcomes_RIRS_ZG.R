# #print(commandArgs(trailingOnly = TRUE))

#LOCAL = FALSE

# gender='male'
# j=40

args=commandArgs(trailingOnly = TRUE)
j      = as.integer(args[1])
gender = as.character(args[2]) #g can be male or female

out_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/2_lmfit_outp/'
in_path <- '/rds/project/jmmh2/rds-jmmh2-hes_data/electronic_health_records/cprd/DataFiles/analysis/zhujie/data_created/'

library(data.table)
library(nlme)
library(survival)
library(tidyverse)
library(haven)

load(paste0(in_path, gender, "/data_lm_",j, "_", gender,".RData"))

# data_ext <- data_ext[derivation == 'derivation', ]

cat(j,"\n")

# test multi-cores
print(Sys.time())
fit_rirs = lme(scaled_corr ~ 0 + exposure + exp_age_corr:exposure + bp_bin:sbp_ind + statin_bin:tchol_ind,  
                random=~ 0 + exposure + exp_age_corr:exposure|patid,
                weights=varIdent(form=~1|exposure),
                data = data_ext,
                control = lmeControl(maxIter=1000,msMaxIter=1000,msVerbose=TRUE,rel.tol=1e-4,
                                     msMaxEval=1000, niterEM = 50))
print(Sys.time())
save(fit_rirs, file = paste0(out_path,j,"_", gender,".RData"))