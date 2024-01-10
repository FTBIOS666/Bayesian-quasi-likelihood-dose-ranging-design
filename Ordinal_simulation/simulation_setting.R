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

Eff_ord_mean = round((Eff_curve+0.8)*(10/8),2)


##########

Ord_mean = sort(unique(as.vector(Eff_ord_mean))) 
score=c(0,1,2,3,4)

Ord_mean[1]
p1=c(0.31,0.3,0.27,0.07,0.05)   
p1%*%score


Ord_mean[2]
p2=c(0.31,0.29,0.28,0.07,0.05)
p2%*%score

Ord_mean[3]
p3=c(0.3,0.3,0.28,0.07,0.05)
p3%*%score

Ord_mean[4]
p4=c(0.28,0.32,0.28,0.07,0.05)
p4%*%score

Ord_mean[5]
p5=c(0.27,0.33,0.28,0.07,0.05)
p5%*%score

Ord_mean[6]
p6=c(0.26,0.34,0.28,0.07,0.05)
p6%*%score

Ord_mean[7]
p7=c(0.26,0.28,0.31,0.1,0.05)
p7%*%score

Ord_mean[8]
p8=c(0.26,0.28,0.28,0.13,0.05)
p8%*%score

Ord_mean[9]
p9=c(0.26,0.25,0.26,0.18,0.05)
p9%*%score

Ord_mean[10]
p10=c(0.25,0.25,0.26,0.19,0.05)
p10%*%score


Ord_mean[11]
p11=c(0.2,0.2,0.37,0.18,0.05)
p11%*%score

Ord_mean[12]
p12=c(0.2,0.2,0.33,0.22,0.05)
p12%*%score

Ord_mean[13]
p13=c(0.2,0.2,0.32,0.23,0.05)
p13%*%score


Ord_mean[14]
p14=c(0.2,0.2,0.3,0.25,0.05)
p14%*%score

Ord_mean[15]
p15=c(0.2,0.2,0.24,0.31,0.05)
p15%*%score

Ord_mean[16]
p16=c(0.2,0.2,0.21,0.34,0.05)
p16%*%score

Ord_mean[17]
p17=c(0.2,0.2,0.2,0.35,0.05)
p17%*%score


Ord_mean[18]
p18=c(0.2,0.2,0.2,0.32,0.08)
p18%*%score

Ord_mean[19]
p19=c(0.2,0.2,0.2,0.3,0.1)
p19%*%score


Ord_mean[20]
p20=c(0.2,0.2,0.2,0.27,0.13)
p20%*%score


Ord_mean[21]
p21=p20=c(0.2,0.2,0.2,0.26,0.14)
p21%*%score


Ord_mean[22]
p22=p20=c(0.2,0.2,0.2,0.25,0.15)
p22%*%score


Ord_mean[23]
p23=p20=c(0.2,0.2,0.2,0.24,0.16)
p23%*%score


Ord_mean[24]
p24=p20=c(0.2,0.2,0.2,0.23,0.17)
p24%*%score


Ord_mean[25]
p25=p20=c(0.2,0.2,0.2,0.22,0.18)
p25%*%score

Ord_mean[26]
p26=c(0.2,0.2,0.2,0.21,0.19)
p26%*%score

Ord_mean[27]
p27=c(0.2,0.2,0.2,0.2,0.2)
p27%*%score

gamma_cut=rbind(p1,p2,p3,p4,p5,
                p6,p7,p8,p9,p10,
                p11, p12, p13, p13, p15,
                p16, p17, p18, p19, p20,
                p21, p22, p23, p24, p25, p26, p27)

gamma_index=cbind(Ord_mean,gamma_cut)


