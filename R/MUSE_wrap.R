#' IPW and AIPW Methods for Multi-cohort Selection Bias in Non-probability Samples
#' 
#' @param K Necessary Input. Number of cohorts. Should be a numeric positive integer.
#' @param Z_names Necessary Input. A character vector containing the names of the Z 
#' variables or disease model covariates.
#' @param intdata_list Necessary Input. A list of size K where each element of 
#' list corresponds to the data for each of the K multiple cohorts including the disease indicator D, the Z variables and the selection variables in the disease model. Each data should be of the form of a data frame. Please include a column named "id" to indicate unique identifiers to the units.
#' @param N Target Population Size.
#' @param UW_CS An indicator variable (TRUE or FALSE) for using the unweighted logistic 
#' regression model with cohort-specific intercepts.
#' @param IPW An indicator variable (TRUE or FALSE) for using the Inverse Probability Weighted Methods (IPW) methods.
#' @param weights_user User specified weights. A numeric vector of weights for the 
#' combined data (not duplicated).
#' @param AIPW An indicator variable (TRUE or FALSE) for using the Joint Augmented 
#' IPW method. If AIPW is TRUE, please also input IPW to be TRUE.
#' @param ipw_method If IPW is TRUE, specify the IPW method to be used. A character 
#' variable. Default is PL (Pseudolikelihood). Other options are Simplex Regression 
#' (SR), Calibration (CL) or User Specified (US).
#' @param extdata If IPW method is set to PL or SR or AIPW is TRUE, please provide 
#' a data frame containing individual level external data which is a probability 
#' sample like NHANES. The external data should contain all the selection variables 
#' and if AIPW is TRUE, then also all the auxiliary score model variables.  Please 
#' include a column named "id" to indicate unique identifiers to the units.
#' @param marginals_list If IPW method is set to CL, please provide a list of size 
#' K in which each element is a numeric vector containing the marginal sums of the selection variables of each of the K cohorts. Please ensure the first element for each of the K numeric vector should be the population size, N. 
#' @param select_var_list If IPW is set to be TRUE, please provide a list of size 
#' K in which each element is a character vector corresponding to the selection 
#' variables' names for each of the K cohorts.
#' @param aux_var_list If AIPW is set to be TRUE, please provide a list of size K 
#' in which each element is a character vector corresponding to the auxiliary score 
#' model variables' names for each of the K cohorts.
#' @param Weights_e If IPW method is set to PL or SR or AIPW is TRUE, please provide 
#' the known selection weights for the external probability sample. The input should 
#' be a numeric vector.
#' @param aux_model If AIPW is true, please provide the auxiliary score model. 
#' Default is XGboost.
#' @param variance An indicator variable (TRUE or FALSE) whether variance should 
#' be computed or not, along with the point estimate.
#' @param type_var If variance is true, indicate the method type to be used for 
#' computing the variance. For the unweighted method, do not provide any type. 
#' For IPW methods, PL and CL, we have two options, asy (asymptotic variance 
#' incorporating the variance from nuisance parameters) and "approx" ignoring the 
#' variance from nuisance parameters. For SR and AIPW, we have only the approx method. 
#' The default for IPW and AIPW methods is approx.
#'
#' @return If variance=TRUE, it will return a list of estimate vector and variance vector.
#' If variance=FALSE, it will return an estimate vector.
#'
#' @examples
#' #library(MASS)
#' 
#' K=3 ## Number of Cohorts
#' set.seed(100)
#' mean_w_p=0
#' mean_z_1=0
#' mean_z_2=0
#' mean_z_3=0
#' corr=0.5
#' var_z_w_p=matrix(c(1,corr,corr,corr,
#'                    corr,1,corr,corr,
#'                    corr,corr,1,corr,
#'                    corr,corr,corr,1),
#'                  nrow=4,ncol=4)
#' 
#' theta=c(-2,0.35,0.45,0.25) ## Theta_Z vector
#' N=5e4 ## Population size
#' 
#' ### selection models
#' dw=1
#' dwz1=c(1,0.8,0.6)
#' dwz2=c(0.6,0.8,1)
#' dwz3=rep(1,3)
#' 
#' gamma_ext=c(-0.6,1.2,0.4,-0.2,0.5)
#' gamma_int_1=c(-1,1.5,0.2,0.8,-0.3)
#' gamma_int_2=c(-1,1.25,0.4,0.6)
#' gamma_int_3=c(-3,0.8,0.5)
#' 
#' ## Generation of population level data
#' simu_popu<-function(N,mean_w_p,mean_z_1,mean_z_2,mean_z_3,
#'                     var_z_w_p,theta,dw){
#'     cov<- MASS::mvrnorm(n = N, mu = c(mean_w_p,mean_z_1,mean_z_2,mean_z_3), Sigma = var_z_w_p)
#'     data <- data.frame(Z1 = cov[, 2], Z2 = cov[, 3], Z3=cov[,4])
#'     W_p=cov[,1]
#'     # Generate random uniforms
#'     #set.seed(5678)
#'     U1 <- runif(N)
#'     #set.seed(4321)
#'     # Generate Disease Status
#'     DISEASE <- expit(theta[1] + theta[2] * data$Z1 + theta[3]*data$Z2 +theta[4]*data$Z3)
#'     data$D   <- ifelse(DISEASE > U1, 1, 0)
#'     # Relate W_p and D
#'     data$W_1 <- W_p + dw* data$D + dwz1[1]*data$Z1 +
#'         dwz2[1]*data$Z2 + dwz3[1]*data$Z3 +
#'         rnorm(n=N,0,1)
#' 
#'     data$W_2 <- W_p + dw* data$D + dwz1[2]*data$Z1 +
#'         dwz2[2]*data$Z2 + dwz3[2]*data$Z3 +
#'         rnorm(n=N,0,1)
#' 
#'     data$W_3 <- W_p + dw* data$D + dwz1[3]*data$Z1 +
#'         dwz2[3]*data$Z2 + dwz3[3]*data$Z3 +
#'         rnorm(n=N,0,1)
#' 
#'     data$id=c(1:N)
#'     return(data)
#' }
#' ## Generation of external individual level data
#' simu_ext<-function(data,gamma_ext){
#'     U2e <- runif(N)
#'     # Generate Sampling Status
#'     SELECT <-0.75*expit(gamma_ext[1] +
#'                             gamma_ext[2]* data$D +
#'                             gamma_ext[3] * data$Z1 +
#'                             gamma_ext[4]* data$Z2 +
#'                             gamma_ext[5] * data$Z3)
#'     S_e  <- ifelse(SELECT > U2e, TRUE, FALSE)
#'     # Observed Data
#'     data_e <- data[which(S_e==1),]
#'     data_e$Select_Weights = 0.75*expit(gamma_ext[1] +
#'                                            gamma_ext[2]* data_e$D +
#'                                            gamma_ext[3] * data_e$Z1 +
#'                                            gamma_ext[4]* data_e$Z2 +
#'                                            gamma_ext[5] * data_e$Z3)
#'     return(data_e)
#' }
#' ## Generation of internal data 1
#' simu_int_1<-function(data,gamma_int_1){
#'     U2i <- runif(N)
#'     # Generate Sampling Status
#'     SELECT <- expit(cbind(1,data$D,data$W_1,data$Z2,data$Z3)
#'                     %*% gamma_int_1)
#'     S_i  <- ifelse(SELECT > U2i, TRUE, FALSE)
#'     # Observed Data
#'     data_i <- data[which(S_i==1),]
#'     return(data_i)
#' }
#' 
#' ## Generation of internal data 2
#' simu_int_2<-function(data,gamma_int_2){
#'     U2i <- runif(N)
#'     # Generate Sampling Status
#'     SELECT <- expit(cbind(1,data$D,data$W_2,data$Z3)
#'                     %*% gamma_int_2)
#'     S_i  <- ifelse(SELECT > U2i, TRUE, FALSE)
#'     # Observed Data
#'     data_i <- data[which(S_i==1),]
#'     return(data_i)
#' }
#' 
#' 
#' ## Generation of internal data 3
#' simu_int_3<-function(data,gamma_int_3){
#'     U2i <- runif(N)
#'     # Generate Sampling Status
#'     SELECT <- expit(cbind(1,data$W_3,data$Z2)
#'                     %*% gamma_int_3)
#'     S_i  <- ifelse(SELECT > U2i, TRUE, FALSE)
#'     # Observed Data
#'     data_i <- data[which(S_i==1),]
#'     return(data_i)
#' }
#' 
#' data=simu_popu(N,mean_w_p,mean_z_1,mean_z_2,mean_z_3,
#'                var_z_w_p,theta,dw)
#' 
#' extdata=simu_ext(data,gamma_ext)
#' 
#' 
#' intdata1=simu_int_1(data,gamma_int_1)
#' intdata2=simu_int_2(data,gamma_int_2)
#' intdata3=simu_int_3(data,gamma_int_3)
#' 
#' ## names of selection variables in each cohort
#' select_var_list=list(c("D","W_1","Z2","Z3"),c("D","W_2","Z3"),c("W_3","Z2"))
#' 
#' ## names of auxiliary variables in each cohort
#' aux_var_list=list(c("D","W_1","Z2","Z3"),c("D","W_2","Z3"),c("W_3","Z2"))
#' 
#' ## list of internal data
#' intdata_list=list(intdata1,intdata2,intdata3)
#' ## names of Z variables
#' Z_names=c("Z1","Z2","Z3")
#' 
#' theta ## actual theta_z
#' res_uw=EHRmuse(K=K,N=N,Z_names=Z_names,
#'                intdata_list=intdata_list,variance = TRUE)
#'
#' @export
EHRmuse<-function(K,Z_names,intdata_list,
                    N=NULL,
                    UW_CS=FALSE,
                    IPW=FALSE,
                    weights_user=NULL,
                    AIPW=FALSE,
                    ipw_method="PL",
                    extdata=NULL,
                    marginals_list=NULL,
                    select_var_list=NULL,
                    aux_var_list=NULL,
                    Weights_e=NULL,
                    aux_model="XGBoost",
                    variance=FALSE,
                    type_var="approx"){
  # require(nleqslv)
  # require(xgboost)
  # require(simplexreg)
  # require(survey)
  # require(dplyr)
  # require(nnet)
  
  message("Please rename the outcome variable to be D.")
  
  if(!is.numeric(K)||!(K>0))
    stop("K, the number of cohorts should be a positive integer.")
  
  if (!is.list(intdata_list)||length(intdata_list)!=K) 
    stop("The input for data for different cohorts
         should be in a list of length K")
  
  if(!is.vector(Z_names)||!is.character(Z_names))
    stop("Z_names should be a character vector")
  
  if(IPW==FALSE){
    if(AIPW==TRUE){
      stop("Please set both IPW and AIPW to be true to use AIPW")
    }
    if(UW_CS==FALSE){
      # source("Unweighted.R")
      est_uw=unweighted(K,intdata_list,Z_names)
      if(variance==TRUE){
        return(list(est=est_uw$est,var=est_uw$var))
      }else{
        return(est_uw$est)
      }
    }else{
      # source("Unweighted_cohort_specific.R")
      est_uw=unweighted_cs(K,intdata_list,Z_names)
      if(variance==TRUE){
        return(list(est=est_uw$est,var=est_uw$var))
      }else{
        return(est_uw$est)
      }
    }
    
  }
  
  if(IPW==TRUE){
    if(AIPW==FALSE){
      if(ipw_method=="PL"){
        if (!is.list(select_var_list)||length((select_var_list))!=K) 
          stop("select_var_list should be a list of length K.")
        
        if (!is.data.frame(extdata)) stop("individual level 
                                        extdata must be a data frame.")
        
        if (!is.numeric(Weights_e) || !is.vector(Weights_e))
          stop("Weights_e must be a numeric vector.")
        
        if(!is.numeric(N)||!(N>0))
          stop("N, population size should be a positive integer.")
        
        # source("PL_estimation.R")
        
        if(variance==TRUE){
          if(!is.character(type_var))
            stop("if variance is true for PL method, please choose asy or approx.")
          
          if(type_var=="asy"){
            # source("PL_var.R")
            
            est_PL_var=variance_pl_actual(K,N,intdata_list,extdata,select_var_list,Weights_e,Z_names)
            
            return(list(est=est_PL_var$final_est,var=est_PL_var$var_est))
            
          }else if(type_var=="approx"){
            
            est_PL_var=PL_est_multi_cohort(K,intdata_list,extdata,select_var_list,Weights_e,Z_names)
            
            return(list(est=est_PL_var$final_est,var=est_PL_var$variance_est))
            
          }
          
        }else{
          
          est_PL=PL_est_multi_cohort(K,intdata_list,extdata,select_var_list,Weights_e,Z_names)
          
          return(est_PL$final_est)
        }
        
      }else if(ipw_method=="SR"){
        
        if (!is.list(select_var_list)||length((select_var_list))!=K) 
          stop("select_var_list should be a list of length K.")
        
        if (!is.data.frame(extdata)) stop("individual level 
                                        extdata must be a data frame.")
        
        if (!is.numeric(Weights_e) || !is.vector(Weights_e))
          stop("Weights_e must be a numeric vector.")
        
        if(!is.numeric(N)||!(N>0))
          stop("N, population size should be a positive integer.")
        
        # source("SR_estimation.R")
        
        if(variance==TRUE){
          message("Only approx method of variance is available for SR")
          est_SR_var=SR_est_var_multi_cohort(K,intdata_list,extdata,select_var_list,Weights_e,Z_names)
          
          return(list(est=est_SR_var$final_est,var=est_SR_var$variance_est))
        }else{
          
          est_SR=SR_est_var_multi_cohort(K,intdata_list,extdata,select_var_list,Weights_e,Z_names)
          
          return(est_SR$final_est)
        }
        
      }else if(ipw_method=="CL"){
        if (!is.list(marginals_list)||length(marginals_list)!=K)
          stop("The input for marginal statistics for 
             selection variables should be in a list of length K")
        
        if (!is.list(select_var_list)||length((select_var_list))!=K) 
          stop("select_var_list should be a list of length K.")
        
        
        if(!is.numeric(N)||!(N>0))
          stop("N, population size should be a positive integer.")
        
        # source("CL_estimation.R")
        
        if(variance==TRUE){
          if(!is.character(type_var))
            stop("if variance is true for CL method, please choose asy or approx.")
          
          if(type_var=="asy"){
            # source("CL_var.R")
            
            est_CL_var=variance_cl_actual(K,N,intdata_list,
                                          select_var_list,marginals_list,Z_names)
            
            return(list(est=est_CL_var$final_est,var=est_CL_var$var_est))
            
          }else if(type_var=="approx"){
            
            est_CL_var=CL_est_multi_cohort(K,intdata_list,select_var_list,marginals_list,Z_names)
            
            return(list(est=est_CL_var$final_est,var=est_CL_var$variance_est))
            
          }
          
        }else{
          
          est_CL=CL_est_multi_cohort(K,intdata_list,select_var_list,marginals_list,Z_names)
          
          return(est_CL$final_est)
        }
        
        
      }else if(ipw_method=="US"){
        if(!is.vector(weights_user)||is.numeric(weights_user))
          stop("If ipw method is US, please provided the user 
               specified weights for the combined data 
               without any repition which is a numeric vector.")
        
        intdata_comb=NULL
        for(i in 1:K){
          intdata=intdata_list[[i]]
          if (!is.data.frame(intdata)) stop("All Internal data should be a data frame.")
          intdata_comb=rbind.data.frame(intdata_comb,intdata)
        }
        intdata_comb1=intdata_comb[!duplicated(intdata_comb$id),]
        
        if(variance==TRUE){
          est=weighted(intdata,weights_user,Z_names)
          
          return(list(est=est_us$final,var=est_us$var))
        }else{
          est=weighted(intdata,weights_user,Z_names)
          
          return(est_us$final)
        }
        
      }else{
        stop("If IPW is true, then input ipw_method as a 
           character from the options, PL, SR, CL 
           or US which stands for user specified weights 
           in the argument weights_user.")
      }
    }else{
      if (!is.list(select_var_list)||length((select_var_list))!=K) 
        stop("select_var_list should be a list of length K.")
      
      if (!is.list(aux_var_list)||length((aux_var_list))!=K) 
        stop("aux_var_list should be a list of length K.")
      
      if (!is.data.frame(extdata)) stop("individual level 
                                        extdata must be a data frame.")
      
      if (!is.numeric(Weights_e) || !is.vector(Weights_e))
        stop("Weights_e must be a numeric vector.")
      
      if(!is.numeric(N)||!(N>0))
        stop("N, population size should be a positive integer.")
      
      if(!(ipw_method %in% c("PL","SR","CL","US")))
        stop("Please choose ipw_method from the four options.")
      
      if(ipw_method=="CL"){
        if (!is.list(marginals_list)||length(marginals_list)!=K)
          stop("The input for marginal statistics for 
             selection variables should be in a list of length K")
        
        if (!is.list(select_var_list)||length((select_var_list))!=K) 
          stop("select_var_list should be a list of length K.")
        
        
        if(!is.numeric(N)||!(N>0))
          stop("N, population size should be a positive integer.")
      }
      
     if(variance==TRUE){
        # source("AIPW_var_approx.R")
        var_AIPW=AIPW_var_approx(K,N,intdata_list,extdata,
                                           select_var_list,
                                           aux_var_list,
                                           Weights_e,Z_names,
                                           aux_model,ipw_method,
                                           marginals_list,
                                           weights_user)
        return(list(est=var_AIPW$est,var=var_AIPW$var))
      }else{
        # source("AIPW_estimation.R")
        AIPW_est=AIPW_joint(K,N,intdata_list,extdata,
                            select_var_list,
                            aux_var_list,
                            Weights_e,Z_names,
                            aux_model,ipw_method,
                            marginals_list,
                            weights_user)
        return(est=AIPW_est$est)
      }
      
      
    }
   
      
  }
}

