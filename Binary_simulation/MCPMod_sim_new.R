

#setwd("C:/Users/ftian1/PhD program/PhD program/Research/Dr. Ying Yuan/Dose_ranging_trial/Code/Ord_simulation")

#!/usr/bin/env Rscript
setwd("/rsrch3/scratch/biostatistics/ftian1/dose_ranging_trial/simulation23")

#---------------------------------------

source("MCPMod_func_new.R")
source("simulation_setting.R")

#------------------------------------------


res_MM=list()
for (i in 1:11) {
    res_MM[[i]]=MCPMod_sim(N=160,d=c(0,1,2,3,4), Eff_mean=Eff_bin_mean[i,],Eff_std=2,
                                    P_Tox=Tox_curve[i,], rho=0.3, w=1, mu_low=0.2, pi_up=0.3,
                                    nsim=5000,seed_number=1234567)
}


save(res_MM,file="res_MM.RData")






