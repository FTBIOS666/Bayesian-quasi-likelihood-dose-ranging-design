
#Project name:
#---------------------------------------------------------
#A Bayesian quasi-likelihood design for dose-ranging studies 
#with bivariate safety and efficacy endpoints
#---------------------------------------------------------

#Author information:
#---------------------------------------------------------
#Author: Feng Tian
#Contact: ftian1@mdanderson.org
#---------------------------------------------------------


#Instruction:
#----------------------------------------------------------
#The code corresponds to the ANOVE approach 
#as a dose-ranging trial design

#The function here can produce PoC power/Correct MED selection percentages/so on

#Below are the description for each of the argument input 
#used in the simulation functions:

#N: the total sample size

#Eff_mean: the true mean efficacy for each dose level (vector)

#Eff_std: the true standard deviation for efficacy response (scalar)

#P_Tox: the true toxicity rate for each dose level (vector)

#rho: the correlation coefficient between the efficacy outcome and the latent continuous toxicity outcome

#mu_low: the clinically relevant effect of dose

#Cf_dr: the cutoff to control the type 1 error of PoC

#nsim: number of simulations for each scenarios

#seed_number: set the seed number for reproduction of results

#---------------------------------------------------------

ANOVA_sim<-function(N, Eff_mean, Eff_std, P_Tox, rho,
                    mu_low, Cf_dr=0.05,
                    nsim=10000, seed_number=1234567,
                    gamma_index, score){
    
    
    #Count the number of dose levels as J
    J=length(Eff_mean)
    
    # Get the true mean marginal efficacy:
    true_eff<-Eff_mean
    
    #Find the true MED location
    j_med_true<-which.min(abs(true_eff[-1]-mu_low-true_eff[1]))+1
    
    #Vector to store the test results for PoC
    PoC_vec<-c()
    
    #Vector to store the estimated MED
    j_med_vec<-c()
    
    set.seed(seed_number)
    for (ii in 1:nsim){
        
        
        #Decide number of patients at each arm
        #--------------------------
        if(N%%J==0){
            #When N1/J as integer: Permutated block randomization
            N_dose <- rep(N/J,J)
        }else{
            #When N1/J is not integer: simple randomization/equal randomization
            N_dose <- rep(N%/%J, J)
            N_dose<-N_dose + c(rmultinom(1,N-sum(N_dose),rep(1/J, J)))
        }
        #--------------------------
        
        
        #Simulate the efficacy and toxicity data
        #------------------------------------------
        #Obtain the mean vector and covariance matrix for each dose level
      
      
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
        res=match_ord_vec(val=ob_data[[j]]$S, bnd=pcent, score)
        return(res)
      }
      
      r=rho
      sdev=Eff_std
      Eff=Eff_mean
      Tox=P_Tox
      n=N_dose
      
      R<-list()
      mu<-list()
      ob_data<-list()
      for (j in 1:J) {
        R[[j]] <- matrix(c(1,(r*sdev),(r*sdev), (sdev)^2), nrow = 2, ncol = 2)
        mu[[j]] <- c(X = 0, S = Eff[j])
        ob_data[[j]] <- as.data.frame(mvtnorm::rmvnorm(n[j], mean = mu[[j]], sigma = R[[j]]))
        ob_data[[j]]$Z <- ifelse(ob_data[[j]]$X>qnorm(1-Tox[j]),1,0)
        ob_data[[j]]$Y=get_ord(vec=ob_data[[j]]$S, Eff, gamma_index, sdev, score)
        
      }
      
      lg = length(ob_data)
      er_index=0
      for (i in 1:lg) {
        if(length(unique(ob_data[[i]]$Y))==1){
          er_index=er_index+1
        }
      }
      
      if(er_index>=1){
        next
      }
       
        #Simulate marginal efficacy outcome 
        OG_Eff_outcome<-list()
        for (j in 1:J) {
            OG_Eff_outcome[[j]]<-ob_data[[j]]$Y
        }
        
        #--------------------------------
        
        
        #ANOVA TEST to detect dose-response
        #-------------------------------------------
        marg_ob_eff<-data.frame(Y=unlist(OG_Eff_outcome), Treatment = rep(seq(1:J),N_dose))
        test_result<-dunnettTest(x=marg_ob_eff$Y,g=marg_ob_eff$Treatment,alternative="greater")
        PoC_vec[ii]<-ifelse(sum(test_result$p.value<=Cf_dr)==0,0,1)
        #-------------------------------------------
        
        
        #If PoC is established, select the MED:
        #-------------------------------------
        if(PoC_vec[ii]==1){
            ob_mu<-sapply(OG_Eff_outcome, mean)
            min_dist<-min(abs(ob_mu[-1] - mu_low - ob_mu[1]))
            j_med_vec[ii]<-which(abs(ob_mu[-1] - mu_low - ob_mu[1]) == min_dist)+1
        }
        #-------------------------------------
    }
    
    
    #Summary simulation results
    #----------------------------
    POC_power <- mean(PoC_vec[!is.na(PoC_vec)])
    correct_med_p <- mean(j_med_vec[!is.na(j_med_vec)]==j_med_true)
    med_pc=tabulate(j_med_vec, nbins = J)/sum(tabulate(j_med_vec,nbins = J))
    p_err <- 100*(mean(abs( j_med_vec[!is.na(j_med_vec)]-j_med_true)))/j_med_true
    #----------------------------
    
    result = list(true_eff,
                  PoC_vec,
                  POC_power,
                  j_med_true,
                  j_med_vec,
                  correct_med_p,
                  p_err,
                  med_pc)
    
    return(result)
    
}

