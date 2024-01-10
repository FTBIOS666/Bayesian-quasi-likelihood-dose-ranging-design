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
#Below are the description for each of the argument input
#used in the simulation functions:

#N1: total number of patients at the first stage
#ncohort: number of patients for each sub-sequent stage
#Eff_mean: the true mean efficacy for each dose level (vector)
#Eff_std: the true standard deviation for efficacy response (scalar)
#P_Tox: the true toxicity rate for each dose level (vector)
#rho: the correlation coefficient between the efficacy outcome and the latent continuous toxicity outcome


#mu_low: the clinically relevant effect of dose
#pi_up: the maximum tolerated dose change compared to the placebo
#K: number of stages
#add_w: the additional penalty if toxicity rates excess a cut
#tox_cut: the cutoff value beyond which the penalty for toxicity will increase

#Cf_tox: the probability cutoff for safety monitoring step
#Cf_fu: the probability cutoff for futility monitoring step
#Cf_dr: the probability cutoff for the detection of dose-response activity
#Cf_med: the probability cutoff to retain dose that are likely to be MED at final decision stage

#w: the weight value for the utility function
#tau: the value reflecting the primary target dose as MED, OBD or both.

#prob_cap: the allocation probability cap for the control arm
#detect_DR: whether to detect dose-response in the trial
#assurance: whether to add a final admissible rule to exclude possible futile dose levels

#nsim: number of simulations for each scenarios
#nsample: number of posterior sample to be drew for parameters of interests
#seed_number: set the seed number for reproduction of results
#---------------------------------------------------------


#function
##############################
QB_sim<-function(N1,ncohort,
                 Eff_mean, Eff_std,P_Tox,rho,
                 mu_low,pi_up,K,add_w,tox_cut,
                 Cf_tox,Cf_fu,Cf_dr, Cf_med=0,
                 w,tau=0.5, prob_cap=0.25,
                 detect_DR=TRUE,assurance=FALSE,
                 nsample=1000,nsim=10000,seed_number=1234567){
    


    #Set up empty vector and matrix to store simulation results
    #-----------------------------------------------
    #Count the number of dose levels as J
    J=length(Eff_mean)

    #To store the selected MED with no final med admissible rule
    j_med_vec<-c()

    #To store the selected MUD for each simulation
    j_mud_vec<-c()

    #Retained dose at the final analysis
    retain_dose_matrix<-matrix(data=NA, nrow = nsim, ncol=J)

    #record the test significance for PoC
    stat_pwr_vec<-c()

    #number of patients at each arm in the final analysis stage
    patient_mat<-matrix(data=NA, nrow = nsim, ncol=J)

    #record the cutoff value for PoC
    DR_cutoff_vec<-c()

    #record early stopped trial
    nstage_vec<-c()

    ADM_one_mat<-matrix(data=NA, nrow = nsim, ncol=J)  #dose-response
    ADM_two_mat<-matrix(data=NA, nrow = nsim, ncol=J)  #futility monitoring
    ADM_tox_mat<-matrix(data=NA, nrow = nsim, ncol=J)  #toxicity monitoring
    ADM_med_mat<-matrix(data=NA, nrow = nsim, ncol=J)  #MED admissible set
    #----------------------------------------------



    #Get the true MED/MUD location
    #---------------------------------------------
    true_eff<-Eff_mean
    true_tox<-P_Tox

    #See if it is necessary to increase the penalty of increasing toxicity
    updated_W<-ifelse(true_tox>(true_tox[1]+tox_cut),(add_w*w+w),w)
    true_U<-true_eff-updated_W*true_tox

    j_med_true<-which.min(abs(true_eff[-1]-mu_low-true_eff[1]))+1
    j_mud_true<-which.max(true_U[-1])+1
    #--------------------------------------------


    #Start simulation:
    #------------------------------------------------
    set.seed(seed_number)
    for (ii in 1:nsim){


        ###############
        #Initial stage#
        ###############

        #Assign patients to each arm
        #--------------------------------------------------------
        if(N1%%J==0){
            #To deal with N1/J as integer
            #Permutated block randomization
            N_dose <- rep(N1/J,J)
        }else{
            #To deal with N1/J as not integer:
            #simple randomization/equal randomization
            N_dose <- rep(N1%/%J, J)
            N_dose<-N_dose + c(rmultinom(1,N1-sum(N_dose),rep(1/J, J)))
        }
        #-----------------------------------------------------

        #Obtain the mean vector and covariance matrix for each dose level
        #get the observed data for each dose
        #-----------------------------------------------------

        ob_data = ET_data(r=rho, sdev=Eff_std, Eff=Eff_mean, n=N_dose, Tox=P_Tox)$sim_data
        
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

        #--------------------------------------------


        #Do the initial data analysis: first step as the safety monitoring
        #----------------------------------------------

        Tox_info = Tox_monitor(dat=ob_data, n=N_dose, Tox_cap=pi_up, Cf=Cf_tox, n_post=nsample)
        #Tox_info
        #---------------------------------------------------------------------------


        #now we do the second part of initial analysis as futility monitoring
        #-----------------------------------------------------------------------

         Eff_info= Eff_monitor(dat=ob_data, n=N_dose,
                               post_tox=Tox_info$post_pi, eff_cap=mu_low, Cf=Cf_fu,
                               n_post=nsample,Tox_data=Tox_info$Tox_data)

        #---------------------------------------------------------------------------


        #############################
        #End of initial stage#
        #############################



        ########################
        #Interim analysis stage#
        ########################
        #-------------------------------------
        k_flag<-1
        nstage<-1
        while(k_flag<K){

            #Count the number of stage
            nstage<-nstage+1

            pre_retain_dose<-ifelse((Tox_info$retain_arm==1 & Eff_info$retain_arm==1),1,0)

            if(sum(pre_retain_dose[-1])==0 ){

                k_flag<-K

            }else{


                #Conduct adaptive allocation
                #----------------------------------------------------

                allo_res<-adap_allo(adm_temp=pre_retain_dose, post_pi=Tox_info$post_pi, post_mu=Eff_info$post_mu,
                                    pi_accp=tox_cut, add_pnl=add_w, pnl=w, t=tau, n=N_dose,ncht=ncohort,
                                    n_post=nsample, pb_cap=prob_cap,eff_cap=mu_low)

                #allo_res
                #--------------------------------------
                N_dose <- N_dose + allo_res$npt
                #-------------------------------------


                #generate new data
                #--------------------------------
                new_ob_data = ET_data(r=rho, sdev=Eff_std, Eff=Eff_mean, n=allo_res$npt, Tox=P_Tox)$sim_data
                for(i in 1:J){
                    ob_data[[i]] = rbind(ob_data[[i]], new_ob_data[[i]])
                }
                #------------------------------------

                #Update tox info and perform safety monitoring
                #------------------------------------
                Tox_info = Tox_monitor(dat=ob_data, n=N_dose, Tox_cap=pi_up, Cf=Cf_tox, n_post=nsample)
                #-----------------------------------

                #Now we update efficacy info and perform futility monitoring
                #-----------------------------------
                Eff_info = Eff_monitor(dat=ob_data, n=N_dose,
                                       post_tox=Tox_info$post_pi, eff_cap=mu_low, Cf=Cf_fu,
                                       n_post=nsample,Tox_data=Tox_info$Tox_data)
                #---------------------------------------------

                #Update the flag to check whether to stop while loop
                #---------------
                k_flag<-k_flag+1
                #---------------
            }

        }

        #-------------------------------------------


        #Record number of stages for each simulation
        #Record patient allocation result
        #-----------------------------------
        nstage_vec[i]<-nstage
        patient_mat[ii,]<-N_dose
        #-----------------------------------

        ###############################
        #End of interim analysis stage#
        ###############################


        ######################
        #Final decision stage#
        ######################

        #Now we perform the detection of dose-response:
        #------------------------------------------------
        ADM_set_one<-rep(0,J)
        for(j in 1:J){
            if(mean((Eff_info$post_mu[,j] >Eff_info$post_mu[,1])) > Cf_dr ){
                ADM_set_one[j]<-1
            }
        }
        ADM_set_one[1]<-1

        #if we do not have this step of PoC, then all dose are retained
        if(!detect_DR){
            ADM_set_one<-rep(1,J)
        }
        #----------------------------------------



        #Now we get the MDE admissible set:
        #------------------------------------------------
        ADM_set_med<-rep(0,J)
        for(j in 1:J){
            if(mean((Eff_info$post_mu[,j]>(Eff_info$post_mu[,1]+mu_low))) > Cf_med ){
                ADM_set_med[j]<-1
            }
        }
        ADM_set_med[1]<-1
        #----------------------


        #record the all the results of retained dose
        #----------------------------------
        ADM_one_mat[ii,]<-ADM_set_one
        ADM_tox_mat[ii,]<-Tox_info$retain_arm
        ADM_two_mat[ii,]<-Eff_info$retain_arm
        ADM_med_mat[ii,]<-ADM_set_med


        if(!assurance){
            retain_dose<-ifelse((Tox_info$retain_arm==1 & Eff_info$retain_arm==1),1,0)
        }else{retain_dose<-ifelse((Tox_info$retain_arm==1 & Eff_info$retain_arm==1 & ADM_set_med==1),1,0)}

        #----------------------------------------


        #record the DR cutoff percentile:
        #--------------------------------------
        DR_cutoff_vec[ii]<-mean(Eff_info$post_mu[,J]>Eff_info$post_mu[,1])
        #-------------------------------------


        #select med and mud
        #------------------------------

        #At the final look we check the existence of med and select med

        if( sum(ADM_set_one[-1])==0 ){

            stat_power = 0
            retain_dose=rep(0,J)

        }else if(sum(ADM_set_one[-1])>0 & sum(retain_dose[-1])==0 ){

            stat_power = 1
            retain_dose=rep(0,J)

        }else{

            stat_power=1

            #---------------------------------------------------
            #Please the source code for function "dose_sel"
            sel_d = dose_sel(a_Tox=Tox_info$aT, b_Tox=Tox_info$bT, a_Eff0=Eff_info$aE0, b_Eff0=Eff_info$bE0,
                             a_Eff1=Eff_info$aE1, b_Eff1=Eff_info$bE1, og_dat=Eff_info$og_eff_dat,
                             Eff_cap=mu_low, J=J, adm_d=retain_dose,pi_accp=tox_cut, add_pnl=add_w, pnl=w)

            #-----------------------------------------------------

            j_med_vec[ii]=sel_d$MED
            j_mud_vec[ii]=sel_d$MUD

            #------------------------------------------------------

        }

        retain_dose_matrix[ii,] = retain_dose
        stat_pwr_vec[ii] = stat_power

        #############################
        #End of final decision stage#
        #############################


    }
    #end of for loop for simulation
    #------------------------------------------------


    #get the correct selection percentages of med/mud
    #--------------------------------------------------------
    #***: ignore the NA value
    correct_med_p = mean(j_med_vec[!is.na(j_med_vec)]==j_med_true)
    med_pc=tabulate(j_med_vec, nbins = J)/sum(tabulate(j_med_vec,nbins = J))
    correct_mud_p = mean(j_mud_vec[!is.na(j_mud_vec)]==j_mud_true)
    mud_pc=tabulate(j_mud_vec, nbins = J)/sum(tabulate(j_mud_vec,nbins = J))
    #--------------------------------------------------------



    #get the pErr for med
    #------------------------------------------------
    #percentage_bias <- 100*(mean(j_med_vec[!is.na(j_med_vec)]-j_med_true))/j_med_true
    percentage_ab_error <- 100*(mean(abs(j_med_vec[!is.na(j_med_vec)]-j_med_true)))/j_med_true
    #----------------------------------------------

    #get the PoC power
    #-------------------------------------
    PoC_power = sum(stat_pwr_vec,na.rm=TRUE)/(nsim-sum(is.na(stat_pwr_vec)))
    #-------------------------------------

    #probability cutoff for PoC
    #--------------------------------------
    dr_cut = quantile(DR_cutoff_vec[!is.na(DR_cutoff_vec)], 0.95)
    #--------------------------------------


    #percentages of patients at each arm
    #--------------------------------------
    patient_mat=na.omit(patient_mat)
    for (i in 1: nrow(patient_mat)) {
        patient_mat[i,]<-patient_mat[i,]/sum(patient_mat[i,])
    }
    avg_allo = apply(patient_mat,2,mean)
    #--------------------------------------

    final_result = list(Efficacy = true_eff,
                        Toxicity=true_tox,
                        Utility = true_U,

                        MED=j_med_true,
                        MED_vector=j_med_vec,
                        MED_percentage=med_pc,
                        PCS_MED=correct_med_p,
                        pError=percentage_ab_error,

                        MUD=j_mud_true,
                        MUD_vector=j_mud_vec,
                        MUD_percentage=mud_pc,
                        PCS_MUD=correct_mud_p,

                        #ADM_one_mat,
                        #ADM_tox_mat,
                        #ADM_two_mat,
                        #ADM_med_mat,
                        retain_dose_matrix,

                        Patients=avg_allo,
                        
                        Cf_vector=DR_cutoff_vec,
                        Cf_dr=dr_cut,
                        PoC = PoC_power)

}
