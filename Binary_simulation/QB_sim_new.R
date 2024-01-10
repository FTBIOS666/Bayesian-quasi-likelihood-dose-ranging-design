

#setwd("C:/Users/ftian1/PhD program/PhD program/Research/Dr. Ying Yuan/Dose_ranging_trial/Code/Bin_simulation")

#!/usr/bin/env Rscript
setwd("/rsrch3/scratch/biostatistics/ftian1/dose_ranging_trial/simulation23")

#-----------------------------
source("simulation_setting.R")
source("QB_func_new.R")
source("Functions_for_sim.R")
#------------------------------


res_QB=list()
for (i in 1:11) {
  res_QB[[i]]=QB_sim(N1=80,ncohort=20,
                     Eff_mean=Eff_bin_mean[i,], Eff_std=2,
                     P_Tox=Tox_curve[i,], rho=0.3,
                     mu_low=0.2,pi_up=0.3,K=5,add_w=0,tox_cut=0.25,
                     Cf_tox=0.9,Cf_fu=0.9,Cf_dr=0.991, Cf_med=0,
                     w=1,tau=0.5, prob_cap=0.25,
                     detect_DR=TRUE,assurance=FALSE,
                     nsample=1000,nsim=5000,seed_number=1234567)
}


save(res_QB,file="res_QB.RData")


#Association between std and dr_cut
#0.1: 0.953 
#0.5: 0.954 
#1:0.958
#2: 0.953
#3:0.955
#4:0.952
#5: 0.957
    

