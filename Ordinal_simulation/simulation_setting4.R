#tempdir()
#dir.create(tempdir())
library(truncnorm) #for generating random truncated normal values
library(Iso)  #for isotonic regression with weight
library(dplyr) #for the between function
library(mvtnorm)
library(PMCMRplus)
library(DescTools)
library("DoseFinding")

#---------------------------------
#Exponential curve: s10
d = c(0, 1, 2, 3, 4)
e0=0.2
e1=1.010e-06 
delta=3.103e-01

eff_exp=e0+e1*(exp(d/delta)-1)
tox_exp=c(0.05, 0.06,0.08, 0.1, 0.24)


#Logistic curve: 



#s2
d = c(0, 1, 2, 3, 4)
e0=0.2
emax=0.57
ed50=1.84
delta=0.23

eff_logis=e0+emax/(1+exp((ed50-d)/delta))
tox_logis=c(0.05, 0.1,0.18,0.2,0.45)




#s7
d = c(0, 1, 2, 3, 4)
e0=0.195
emax=0.6
ed50=2.3
delta=0.5

eff_logis7=e0+emax/(1+exp((ed50-d)/delta))
tox_logis7=c(0.05, 0.06,0.1, 0.12, 0.32)



#Emax curve: s1
d = c(0, 1, 2, 3, 4)
e0=0.2
emax=0.75
ed50=1

eff_emax=e0+(emax*d)/(ed50+d)
tox_emax=c(0.05, 0.1,0.11, 0.3, 0.34)


#Linear log-dose curve: s5
d = c(0, 1, 2, 3, 4)
e0=0.2
delta=0.34

eff_linlog=e0+delta*log(d+1)
tox_linlog=c(0.05, 0.07, 0.22, 0.34, 0.45)


###########################################
Eff_curve = rbind(
  c(0.2,0.2,0.2,0.2,0.2),
  #1
  eff_emax,
  #2
  eff_linlog,
  #3
  eff_logis,
  #4
  eff_exp,
  
  #5
  c(0.20,0.34,0.68,0.76,0.78),
  #6
  c(0.2, 0.21, 0.72, 0.75, 0.80),
  #7
  eff_logis7,
  
  #8
  c(0.2, 0.23, 0.32, 0.65, 0.79),
  #9
  c(0.20,0.23,0.25,0.72,0.80),
  #10
  c(0.20,0.20,0.22,0.54,0.80)
)


Tox_curve = rbind(
  c(0.05,0.05,0.05,0.05,0.05),
  #1
  tox_emax,
  #2
  tox_linlog,
  #3
  tox_logis,
  #4
  tox_exp,
  
  #5
  c(0.05, 0.12,0.14,0.35,0.45),
  #6
  c(0.05, 0.06, 0.15, 0.24,0.28),
  #7
  tox_logis7,
  
  #8
  c(0.05, 0.08, 0.1, 0.32, 0.45),
  #9
  c(0.05, 0.06, 0.08, 0.15, 0.34),
  #10
  c(0.05,0.06, 0.08, 0.18, 0.2)
  
)

###########

Eff_ord_mean = round(2*Eff_curve+0.4,2)


##########

Ord_mean = sort(unique(as.vector(Eff_ord_mean))) 
score=c(0,1,2,3,4)



setwd("C:/Users/ftian1/PhD program/PhD program/Research/Dr. Ying Yuan/Dose_ranging_trial/Code/Ord_simulation")
#setwd("/rsrch3/scratch/biostatistics/ftian1/dose_ranging_trial/simulation22")

#setting4=data.frame(read_excel("setting_decreasing2.xlsx"))

#save(setting4, file="setting4.RData")

load("setting4.RData")

gamma_cut=setting4

gamma_index=cbind(Ord_mean,gamma_cut)


