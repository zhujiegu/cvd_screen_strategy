library(dplyr)

diff <- readRDS('scale_corr_diff.rds')
diff$exposure %>% unique

# Smoking all the same
diff_smk <- diff %>% filter(exposure=='smokbin')
all.equal(diff_smk$fr, diff_smk$zg)


diff_sbp <- diff %>% filter(exposure=='sbp')
diff_sbp
diff_sbp %>% filter(fr==0)

miss_sbp <- diff_sbp %>% filter(fr==0)
miss_sbp$patid %>% unique() %>% length

diff_sbp$patid %>% unique %>% length
