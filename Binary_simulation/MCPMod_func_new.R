
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
#---------------------------------------------------------

#The code corresponds to the MCPMod approach 
#as a dose-ranging trial design

#The function here can produce PoC power/Correct MED& MUD selection percentages/so on


#Below are the description for each of the argument input 
#used in the simulation functions:

#N: the total sample size

#d: the true dose levels (vector)

#Eff_mean: the true mean efficacy for each dose level (vector)

#Eff_std: the true standard deviation for efficacy response (scalar)

#P_Tox: the true toxicity rate for each dose level (vector)

#rho: the correlation coefficient between the efficacy outcome and the latent continuous toxicity outcome

#w: the weight value for the utility function

#mu_low: the clinically relevant effect of dose

#pi_up: the maximum tolerated dose change compared to the placebo

#nsim: number of simulations for each scenarios

#seed_number: set the seed number for reproduction of results

#---------------------------------------------------------


#function
##############################
MCPMod_sim= function(N,d, Eff_mean,Eff_std,P_Tox,rho, 
                    w,mu_low,pi_up,
                    nsim=10000,seed_number=1234567){
    
    
    #Vector to store simulation results
    #-------------------------------------------
    #Count the number of dose levels as J
    J=length(d)
    
    #Vector to store the test results for PoC
    PoC_vec=c()
    
    #Vector to store the test results for estimated MED (continuous)
    est_dose = c()
    
    #suppress warnings for unstable computation of
    #mvtnorm::pmvt: Completion with error > abseps.
    options(warn = -1)
    
    
    #To store the selected MED with no final med admissible rule
    j_med_vec<-c() 
    
    #To store the selected MUD for each simulation
    j_mud_vec<-c()
    
    #-------------------------------------------
    
    #Find the true MED& MUD location
    #------------------------------
    true_eff<-Eff_mean
    true_tox <- P_Tox
    true_U<-true_eff-w*true_tox 
    
    j_med_true<-which.min(abs(true_eff[-1]-mu_low-true_eff[1]))+1
    j_mud_true<-which.max(true_U[-1])+1
    #------------------------------------
    
    
    #start the simulation
    #--------------------------
    set.seed(seed_number)
    for (ii in 1:nsim) {
        
        
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
        
        
        
        #Obtain the mean vector and covariance matrix for each dose level
        #get the observed data for each dose
        #--------------------------------
        r=rho
        sdev=Eff_std
        Eff=Eff_mean
        Tox=P_Tox
        n=N_dose
        
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
            ob_data[[j]]$Y=ifelse(ob_data[[j]]$S > qnorm(p=(1-Eff[j]), mean=Eff[j], sd=sdev), 1,0)
          }
        }
      
      
       
        OG_Tox_outcome <-list()
        for(j in 1:J){
            OG_Tox_outcome[[j]]<-c(ob_data[[j]]$Z)
        }
        
        Tox_outcome <- unlist(OG_Tox_outcome)
        
        #Simulate marginal efficacy outcome 
        OG_Eff_outcome<-list()
        for (j in 1:J) {
            OG_Eff_outcome[[j]]<-ob_data[[j]]$Y
        }
        
        Eff_outcome <- unlist(OG_Eff_outcome)
        Eff_data<-data.frame(dose=rep(d, N_dose), resp=Eff_outcome)
        Tox_data<-data.frame(dose=rep(d, N_dose), resp=Tox_outcome)
        #-----------------------------------------------------------------------
        
        
        
        #run MCP-MOD to get the estimated MED:
        #-----------------------------------------------------------------------
        
        logit <- function(p) log(p / (1 - p))
        logit<-Vectorize(logit)
        inv_logit <- function(y) 1 / (1 + exp(-y))
        inv_logit <-Vectorize( inv_logit )
        
        Eff_models <- Mods(doses=d, placEff=logit(0.2), maxEff=logit(0.5)-logit(0.2),
                           linear = NULL, linlog = NULL, emax = 0.5,
                           exponential = c(1), logistic=c(1, 0.4),
                           fullMod = FALSE,
                           addArgs=list(off=1), direction = "increasing")
        
        
        fit_nocov <- glm(resp~factor(dose) + 0, data = Eff_data, family = binomial)
        mu_hat <- coef(fit_nocov)
        S_hat <- vcov(fit_nocov)
        
        Eff_dfe <- DoseFinding::MCPMod(d,mu_hat, S=S_hat ,  models=Eff_models, 
                                       alpha = 0.05, type = "general", selModel = "maxT",
                                       pVal = TRUE, doseType = "TD", Delta = mu_low, 
                                       alternative="one.sided")
        
        
        
        if(is.null(Eff_dfe$selMod)){
            #Null results means dose-response is not detected
            PoC_vec[ii]<-0
        }else{
            PoC_vec[ii]<-1
            #est_dose[ii] <- Eff_dfe$doseEst[Eff_dfe$selMod]
            
            fitted_logit_mu = predict(Eff_dfe, predType = "ls-means", doseSeq = d)[Eff_dfe$selMod]
            
            #Transform back from logit-sclae to range of [0,1] as toxicity rates
            fitted_mu <- as.vector(inv_logit(fitted_logit_mu))
            
            j_med_vec[ii]<- which.min(abs(fitted_mu[-1]-mu_low - fitted_mu[1]))+1
            
            #-----------------------------------------------
            
            
            
            #Now we perform MCPMod to binary toxicity data
            #The model parameters here should be based on logit scale
            #-----------------------------------------------
            logit <- function(p) log(p / (1 - p))
            logit<-Vectorize(logit)
            inv_logit <- function(y) 1 / (1 + exp(-y))
            inv_logit <-Vectorize( inv_logit )
            
            Tox_models <- Mods(doses=d, placEff=logit(0.05), maxEff=logit(0.5)-logit(0.05),
                               linear = NULL, linlog = NULL, emax = 0.25,
                               exponential = c(1), logistic=c(1, 0.4),
                               fullMod = FALSE,
                               addArgs=list(off=1), direction = "increasing")
            
            
            fit_nocov <- glm(resp~factor(dose) + 0, data = Tox_data, family = binomial)
            mu_hat <- coef(fit_nocov)
            S_hat <- vcov(fit_nocov)
            
            Tox_dfe <- DoseFinding::MCPMod(d,mu_hat, S=S_hat ,  models=Tox_models, 
                                           alpha = 1, type = "general", selModel = "maxT",
                                           pVal = TRUE, doseType = "TD", Delta = pi_up, 
                                           alternative="one.sided")
            
            
            fitted_logit_pi = predict(Tox_dfe, predType = "ls-means", doseSeq = d)[Tox_dfe$selMod]
            
            #Transform back from logit-sclae to range of [0,1] as toxicity rates
            fitted_pi <- as.vector(inv_logit(fitted_logit_pi))
            #-----------------------------------------------
            
            #Get the estimated utility values
            #-----------------------------------------------
            fitted_u = fitted_mu-w*fitted_pi
            j_mud_vec[ii]<-which.max(fitted_u[-1])+1
            #-----------------------------------------------
            
        }
        
    }
    
    #end of simulation
    #--------------------------
    
    #Summary simulation results
    #----------------------------
    POC_power <- mean(PoC_vec[!is.na(PoC_vec)])
    
    correct_med_p = mean(j_med_vec[!is.na(j_med_vec)]==j_med_true)
    med_pc=tabulate(j_med_vec, nbins = J)/sum(tabulate(j_med_vec,nbins = J))
    
    correct_mud_p = mean(j_mud_vec[!is.na(j_mud_vec)]==j_mud_true)
    mud_pc=tabulate(j_mud_vec, nbins = J)/sum(tabulate(j_mud_vec,nbins = J))
    
    p_err <- 100*(mean(abs( j_med_vec[!is.na(j_med_vec)]-j_med_true)))/j_med_true
    #----------------------------
    
    
    result = list(true_eff,
                  true_tox,
                  true_U,
                  ###
                  PoC_vec,
                  POC_power,
                  ###
                  j_med_true,
                  j_med_vec,
                  correct_med_p,
                  p_err,
                  med_pc,
                  ###
                  j_mud_true,
                  j_mud_vec,
                  correct_mud_p,
                  mud_pc)
    
    return(result)
}









