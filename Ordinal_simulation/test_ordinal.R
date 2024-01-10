
r=0.3
sdev=1
Eff=Eff_ord_mean[2,]
n=c(100000,100000,100000,100000,100000)
Tox=c(0.05,0.1,0.15,0.2,0.25)

#----------------------------------

match_ord=function(val,bnd, score){
  if(val<bnd[1]){
    out=score[1]
  }else if(val>=bnd[1] & val<bnd[2]){
    out=score[2]
  }else if(val>=bnd[2] & val<bnd[3]){
    out=score[3]
  }else if(val>=bnd[3] & val<bnd[4]){
    out=score[4]
  }else if(val>=bnd[4]){
    out=score[5]
  }
  return(out)
}

match_ord_vec = Vectorize(match_ord, vectorize.args = "val")


get_ord=function(vec, Eff, gamma_index, sdev, score){
  lab = which(gamma_index[,1]==Eff[j])
  interv = gamma_index[lab,c(2:6)]
  qintv = cumsum(interv)
  pcent =sapply(qintv, FUN=function(x) qnorm(x, mean=Eff[j], sd=sdev))
  res=match_ord_vec(ob_data[[j]]$S, bnd=pcent, score)
}


#------------------------------------


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
    ob_data[[j]]$Y=get_ord(vec=ob_data[[j]]$Y, Eff, gamma_index, sdev, score)
  }
}

Eff

mean(ob_data[[1]]$Y)
mean(ob_data[[2]]$Y)
mean(ob_data[[3]]$Y)
mean(ob_data[[4]]$Y)
mean(ob_data[[5]]$Y)

  
