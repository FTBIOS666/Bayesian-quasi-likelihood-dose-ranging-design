
ET_data <- function(r, sdev, Eff, n, Tox){
    R<-list()
    mu<-list()
    ob_data<-list()
    
    

    for (j in 1:length(n)) {
        if(n[j]==0){
            ob_data[[j]]<-data.frame(X=integer(0),Y=integer(0),Z=integer(0))
        }else if(n[j]>0){
            R[[j]] <- matrix(c(1,(r*sdev),(r*sdev), (sdev)^2), nrow = 2, ncol = 2)
            mu[[j]] <- c(X = 0, Y = Eff[j])
            ob_data[[j]] <- as.data.frame(mvtnorm::rmvnorm(n[j], mean = mu[[j]], sigma = R[[j]]))
            ob_data[[j]]$Z <- ifelse(ob_data[[j]]$X>qnorm(1-Tox[j]),1,0)
        }
    }
    
    dat=list(sim_data = ob_data)
    
    return(dat)
    
}



#------------------------------------------------------------------------------

Tox_monitor = function(dat, n, Tox_cap, Cf, n_post){
    
    Tox_outcome <- list()
    a_Tox<-c()
    b_Tox<-c()
    unrestricted_post_pi<-matrix(data=NA, nrow = n_post, ncol = length(n))
    for(j in 1:length(n)){
        Tox_outcome[[j]] <- c(dat[[j]]$Z)
        a_Tox[j]<-1 + length(which(Tox_outcome[[j]]==1))
        b_Tox[j]<-1+ n[j] - length(which(Tox_outcome[[j]]==1))
        unrestricted_post_pi[,j]<-rbeta(n=n_post, shape1=a_Tox[j], shape2=b_Tox[j])
    }
    
    var_pi<-(a_Tox*b_Tox)/(((a_Tox+b_Tox)^2)*(a_Tox+b_Tox+1))
    
    #Then do isotonic regression under Bayesian framework:
    restricted_post_pi<-matrix(data=NA, nrow = n_post, ncol = length(n))
    for (i in 1:n_post) {
        restricted_post_pi[i,]<-pava(y = unrestricted_post_pi[i,], w=1/var_pi)
    }
    
    
    marginal_post_pi<-restricted_post_pi
    
    ADM_set_tox<-rep(0,length(n))
    for(j in 1:length(n)){
        if((mean((marginal_post_pi[,j]-marginal_post_pi[,1])> Tox_cap)) > Cf){
            ADM_set_tox[j]<-0
        }else{ADM_set_tox[j]<-1}
    }
    ADM_set_tox[1]<-1
    
    Tox_result = list(Tox_data=Tox_outcome,
                      post_pi=marginal_post_pi,
                      retain_arm=ADM_set_tox,
                      aT = a_Tox,
                      bT=b_Tox)
    
    return(Tox_result )
    
}


#-------------------------------------------------------------------------------
Eff_monitor=function(dat, n,post_tox, eff_cap, Cf, n_post,Tox_data){
        
    OG_Eff0_outcome<-list()
    OG_Eff1_outcome<-list()
    OG_Eff_outcome<-list()
    num_T0<-c()
    num_T1<-c()
    Eff0_outcome<-list()
    Eff1_outcome<-list()
    a_Eff0<-c()
    b_Eff0<-c()
    a_Eff1<-c()
    b_Eff1<-c()
    unrestricted_post_mu0<-matrix(data=NA, nrow = n_post, ncol = length(n))
    unrestricted_post_mu1<-matrix(data=NA, nrow = n_post, ncol = length(n))
    var_mu0<-c()
    var_mu1<-c()
    for (j in 1:length(n)) {
        
        num_T0[j]<-length(which(Tox_data[[j]]==0))
        if(num_T0[j]==0){
            OG_Eff0_outcome[[j]]<-integer(0)
        }else{
            OG_Eff0_outcome[[j]]<-dat[[j]]$Y[dat[[j]]$Z==0]
        }
        
        num_T1[j]<-length(which(Tox_data[[j]]==1))
        if(num_T1[j]==0){
            OG_Eff1_outcome[[j]]<-integer(0)
        }else{
            OG_Eff1_outcome[[j]]<-dat[[j]]$Y[dat[[j]]$Z==1]
        }
        
        OG_Eff_outcome[[j]]<-c(OG_Eff0_outcome[[j]],OG_Eff1_outcome[[j]])
        
        #####
        
        if(num_T0[j]==0){
            Eff0_outcome[[j]]<-integer(0)
        }else{
            Eff0_outcome[[j]]<-((OG_Eff0_outcome[[j]] - min(OG_Eff_outcome[[j]]))/
                                    (max(OG_Eff_outcome[[j]]) - min(OG_Eff_outcome[[j]])))
        }
        
        if(num_T1[j]==0){
            Eff1_outcome[[j]]<-integer(0)
        }else{
            Eff1_outcome[[j]]<-((OG_Eff1_outcome[[j]] - min(OG_Eff_outcome[[j]]))/
                                    (max(OG_Eff_outcome[[j]]) - min(OG_Eff_outcome[[j]])))
        }
        
        
        #####
        a_Eff0[j]<-1 + sum(Eff0_outcome[[j]])
        b_Eff0[j]<-1 + num_T0[j] -sum(Eff0_outcome[[j]])
        a_Eff1[j]<-1 + sum(Eff1_outcome[[j]])
        b_Eff1[j]<-1 + num_T1[j] -sum(Eff1_outcome[[j]])
        
        #####
        unrestricted_post_mu0[,j]<-rbeta(n=n_post, shape1=a_Eff0[j], shape2=b_Eff0[j])
        unrestricted_post_mu1[,j]<-rbeta(n=n_post, shape1=a_Eff1[j], shape2=b_Eff1[j])
        unrestricted_post_mu0[,j]<-(unrestricted_post_mu0[,j]* (max(OG_Eff_outcome[[j]]) - min(OG_Eff_outcome[[j]])) +
                                        min(OG_Eff_outcome[[j]]) )
        unrestricted_post_mu1[,j]<-(unrestricted_post_mu1[,j]* (max(OG_Eff_outcome[[j]]) - min(OG_Eff_outcome[[j]])) +
                                        min(OG_Eff_outcome[[j]]) )
        var_mu0[j]<-(((a_Eff0[j]*b_Eff0[j])/(((a_Eff0[j]+b_Eff0[j] )^2)*(a_Eff0[j]+b_Eff0[j]+1)))*  
                         (max(OG_Eff_outcome[[j]]) - min(OG_Eff_outcome[[j]]))^2 )
        var_mu1[j]<-(((a_Eff1[j]*b_Eff1[j])/(((a_Eff1[j]+b_Eff1[j] )^2)*(a_Eff1[j]+b_Eff1[j]+1)))*  
                         (max(OG_Eff_outcome[[j]]) - min(OG_Eff_outcome[[j]]))^2 )
        
        
    }
    
    #Then do isotonic regresison under Bayesian framework:
    unrestricted_post_mu<-matrix(data=NA, nrow = n_post, ncol = length(n))
    restricted_post_mu<-matrix(data=NA, nrow = n_post, ncol = length(n))
    for (i in 1:n_post){
        
        unrestricted_post_mu[i,]<-(post_tox[i,]*unrestricted_post_mu1[i,]+
                                       (1-post_tox[i,])*unrestricted_post_mu0[i,])
        restricted_post_mu[i,]<-pava(y = unrestricted_post_mu[i,], 
                                     w= 1/(((post_tox[i,])^2)*var_mu1 + 
                                               ((1-post_tox[i,])^2)*var_mu0) )
    }
    
    marginal_post_mu<-restricted_post_mu
    
    ADM_set_two<-rep(1,length(n))
    for(j in 1:length(n)){
        if((mean((marginal_post_mu[,j]<=marginal_post_mu[,1]))) > Cf){
            ADM_set_two[j]<-0
        }
    }
    ADM_set_two[1]<-1
    
    
    Eff_result = list(retain_arm = ADM_set_two,
                      post_mu = marginal_post_mu,
                      og_eff_dat = OG_Eff_outcome,
                      aE0=a_Eff0,
                      bE0=b_Eff0,
                      aE1=a_Eff1,
                      bE1=b_Eff1)
    
    return(Eff_result)
    
}

#-------------------------------------------------------------------------------

adap_allo = function(adm_temp, post_pi, post_mu, pi_accp, add_pnl, pnl, t, n,ncht, n_post, pb_cap,eff_cap){
    
    post_u_matrix = matrix(data=NA, nrow=n_post, ncol=length(n))
    exp_post_u_vec=c()
    for(j in 1:length(n)){
        wegt=ifelse(post_pi[,j]>(post_pi[,1]+pi_accp),(add_pnl*pnl+pnl),pnl)
        post_u_matrix[,j] = post_mu[,j]-wegt*post_pi[,j]
        exp_post_u_vec[j]=mean(post_u_matrix[,j])
        
    }
    
    #Only arms in the admissible set will be counted for adaptive allocation
    #post_u_matrix[,which(adm_temp==0)] = -999
    #marginal_post_mu[,which(adm_temp==0)] = -999
    
    max_u_dbn <- (as.vector(tabulate(apply(post_u_matrix,1,which.max),nbins=length(n))))/n_post
    
    close_med_dbn <- (as.vector(tabulate(apply(post_mu,1,function(x) 
        which.min(abs(x- (x[1]+eff_cap)))),nbins=length(n))))/n_post
    
    P_max <- max(max(max_u_dbn[-1]), max(close_med_dbn[-1]))
    
    P_adap = c()
    if(t==1){
        close_med_dbn[1] = max(close_med_dbn[-1])
    }else if(t==0){
        max_u_dbn[1] = max(max_u_dbn[-1])
    }else if(t==0.5){
        temp_u<-max_u_dbn
        temp_e<-close_med_dbn
        max_u_dbn[c((which.max(temp_e[-1])+1), (which.max(temp_u[-1])+1))] <- P_max
        close_med_dbn[c((which.max(temp_e[-1])+1), (which.max(temp_u[-1])+1))] <- P_max
    }
    P_adap <- t*close_med_dbn + (1-t)*max_u_dbn
    
    #the standardization is only for dose arms in the admissible set.
    P_adap[adm_temp==0]<-0
    
    if(sum(P_adap[-1])==0){
        P_adap<-rep(1/(length(n)), (length(n)))
    }
    
    denom<-sum(P_adap)-P_adap[1]+max(P_adap[-1])
    P_adap[1]<-max(P_adap[-1])
    P_adap<-P_adap/denom
    
    if(P_adap[1]>pb_cap){
        P_adap[1]=pb_cap
    }
    
    npts<-c(rmultinom(1,ncht,P_adap))
    
    return(list(npt = npts))
    
}

#-------------------------------------------------------------------------------

dose_sel=function(a_Tox, b_Tox, a_Eff0, b_Eff0, a_Eff1, b_Eff1,og_dat, Eff_cap, J, adm_d,pi_accp, add_pnl, pnl){
    
    #First get the closed-form posterior mean of toxicity rate:
    #--------------------------
    mean_pi<-c()
    var_pi<-c()
    mean_pi<-a_Tox/(a_Tox+b_Tox)
    var_pi<-(a_Tox*b_Tox)/(((a_Tox+b_Tox)^2)*(a_Tox+b_Tox+1))
    restricted_mean_pi<-pava(mean_pi, w=1/var_pi)
    #---------------------------
    
    #Then get the closed-form posterior mean of mean efficacy
    #and select MED
    #-------------------------------------
    mean_mu0<-c()
    mean_mu1<-c()
    mean_mu<-c()
    var_mu<-c()
    var_mu0<-c()
    var_mu1<-c()
    for (j in 1:J) {
        
        mean_mu0[j]<-(a_Eff0[j]/(a_Eff0[j]+b_Eff0[j])) * (max(og_dat[[j]]) - min(og_dat[[j]])) +
            min(og_dat[[j]])
        mean_mu1[j]<-(a_Eff1[j]/(a_Eff1[j]+b_Eff1[j])) * (max(og_dat[[j]]) - min(og_dat[[j]])) +
            min(og_dat[[j]])
        mean_mu[j]<-mean_mu0[j]*(1-restricted_mean_pi[j]) + mean_mu1[j]*restricted_mean_pi[j]
        
        
        var_mu0[j]<-(((a_Eff0[j]*b_Eff0[j])/(((a_Eff0[j]+b_Eff0[j] )^2)*(a_Eff0[j]+b_Eff0[j]+1)))*  
                         (max(og_dat[[j]]) - min(og_dat[[j]]))^2 )
        var_mu1[j]<-(((a_Eff1[j]*b_Eff1[j])/(((a_Eff1[j]+b_Eff1[j] )^2)*(a_Eff1[j]+b_Eff1[j]+1)))*  
                         (max(og_dat[[j]]) - min(og_dat[[j]]))^2 )
        
        
        
    }
    restricted_mean_mu<-pava(y = mean_mu, 
                             w= 1/(((restricted_mean_pi)^2)*var_mu1 + 
                                       ((1-restricted_mean_pi)^2)*var_mu0))
    
    min_dist = min(abs(restricted_mean_mu[adm_d==1& c(1:J)!=1] - restricted_mean_mu[1]  - Eff_cap))
    #---------------------------------------------
    
    
    #Select mud
    #------------------------------------------------
    mean_u_vec=c()
    for(j in 1:J){
        wegt=ifelse(restricted_mean_pi[j] > (restricted_mean_pi[1]+pi_accp), (add_pnl*pnl+pnl),pnl)
        mean_u_vec[j] = restricted_mean_mu[j] - wegt*restricted_mean_pi[j]
    }
    
    med=min(which(abs(restricted_mean_mu[-1] - restricted_mean_mu[1] - Eff_cap)==min_dist)+1)
    mud=min(which(mean_u_vec[-1]  == max(mean_u_vec[adm_d==1& c(1:J)!=1]))+1)
    
    sel_res = list(MED=med, MUD=mud)
    
    return(sel_res)
    
}


