

#setwd("C:/Users/ftian1/PhD program/PhD program/Research/Dr. Ying Yuan/Dose_ranging_trial/Code/Ord_simulation")

#!/usr/bin/env Rscript
setwd("/rsrch3/scratch/biostatistics/ftian1/dose_ranging_trial/simulation22")


#------------------------------------

source("ANOVA_func_new.R")
source("simulation_setting.R")

#------------------------------------


res_ANOVA = list()
for (i in 1:nrow(Eff_curve)) {
  res_ANOVA[[i]] = ANOVA_sim(N=200, Eff_mean=Eff_ord_mean[i,], Eff_std=4, 
                             P_Tox=Tox_curve[i,], rho=0.3,
                             mu_low=1.6, Cf_dr=0.05,
                             nsim=5000, seed_number=1234567,
                             gamma_index, score)
  
}


save(res_ANOVA, file="res_ANOVA.RData")













