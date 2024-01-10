
r=0.3
sdev=1
Eff=Eff_curve[2,]
n=rep(10000,5)
Tox=c(0.05,0.1,0.15,0.2,0.25)

#----------------------------------

R<-list()
mu<-list()
ob_data<-list()

for (j in 1:length(n)) {
  if(n[j]==0){
    ob_data[[j]]<-data.frame(X=integer(0),S=integer(0),Z=integer(0))
  }else if(n[j]>0){
    R[[j]] <- matrix(c(1,(r*sdev),(r*sdev), (sdev)^2), nrow = 2, ncol = 2)
    mu[[j]] <- c(X = 0, S = Eff[j])
    ob_data[[j]] <- as.data.frame(mvtnorm::rmvnorm(n[j], mean = mu[[j]], sigma = R[[j]]))
    ob_data[[j]]$Z <- ifelse(ob_data[[j]]$X>qnorm(1-Tox[j]),1,0)
    ob_data[[j]]$Y=ifelse(ob_data[[j]]$S < qnorm(p=Eff[j], mean=Eff[j], sd=sdev), 1,0)
  }
}
dat=list(sim_data = ob_data)

mean(ob_data[[1]]$Y)
mean(ob_data[[2]]$Y)
mean(ob_data[[3]]$Y)
mean(ob_data[[4]]$Y)
mean(ob_data[[5]]$Y)

  
