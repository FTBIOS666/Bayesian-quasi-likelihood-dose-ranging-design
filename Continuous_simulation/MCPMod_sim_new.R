

#setwd("C:/Users/ftian1/PhD program/PhD program/Research/Dr. Ying Yuan/Dose_ranging_trial/Code/Main_simulation")

#!/usr/bin/env Rscript
setwd("/rsrch3/scratch/biostatistics/ftian1/dose_ranging_trial/simulation21")

#---------------------------------------

source("MCPMod_func_new.R")
source("simulation_setting.R")

#------------------------------------------


res_MM=list()
for (i in 1:11) {
    res_MM[[i]]=MCPMod_sim(N=200,d=c(0,1,2,3,4), Eff_mean=Eff_curve[i,],Eff_std=1,
                                    P_Tox=Tox_curve[i,], rho=0.3, w=2, mu_low=0.4, pi_up=0.3,
                                    nsim=5000,seed_number=1234567)
}


save(res_MM,file="res_MM.RData")






