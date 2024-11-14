SR_est_multi_cohort <- function(K,Select_list,
                                Select_e_list,Weights_e,
                                intdata_id_list,
                                extdata_id){
  
  if (!is.list(Select_list) || !is.list(Select_e_list) 
      || !is.list(intdata_id_list)) 
    stop("All selection input data should be lists for 
         multiple cohorts and the external data.")
  
  if (length((Select_list))!=K || length((Select_e_list))!=K 
      || length((intdata_id_list))!=K) 
    stop("All input lists should have length K.")
  
  ## Sample weights for the external data
  if (!is.numeric(Weights_e) || !is.vector(Weights_e))
    stop("'Weights_e' must be a numeric vector.")
  
  if (!is.numeric(extdata_id) || !is.vector(extdata_id))
    stop("Id of extdata must be a numeric vector.")
  
  # Initialize lists to store results for each cohort
  estweights_list <- list()
  modsimp <-list()
  modmult <-list()

  # Loop over each cohort
  for (k in 1:K){
    select_k=Select_list[[k]]
    select_e_k=Select_e_list[[k]]
    intdata_k_id=intdata_id_list[[k]]
    
    if (!is.numeric(intdata_k_id) || !is.vector(intdata_k_id))
      stop("Id of intdata must be a numeric vector.")
    
    if (!is.list(select_k)) 
      stop("Selection_variable_input_for_each_cohort_should_be_a_list")
    
    if (length((select_k))!=K) 
      stop("Selection_variable_input_for_each_cohort_should_be_a_list of length K")
    
    
    
    ## Selection variables exclusive to cohort K
    
    select_k_var=select_k[[k]]
    select_k_e=Select_e_list[[k]]
    
    # Check select_k_var and select_k_e to be data frames or not
    if (!is.data.frame(select_k_var)) stop("select_k_var must be a data frame.")
    
    if (!is.data.frame(select_k_e)) stop("select_k_e must be a data frame.")
    
    formula_string <- paste("1/Weights_e ~", 
                            paste(colnames(select_k_e), collapse = " + "))
    
    # Convert the string to a formula
    formula <- as.formula(formula_string)
    
    modsimp_k=simplexreg(data=data.frame(select_k_e),
                         formula= formula,
                         link="logit")
    
    # Find common and exclusive IDs
    inter <- intersect(intdata_k_id, extdata_id)
    intju <- setdiff(intdata_k_id, inter)
    extju <- setdiff(extdata_id, inter)
    
    # Subset the data and add the 'co' column in one step
    bothie <- transform(select_k_e[extdata_id %in% inter, ], co = 1)
    colnames(bothie)=c(colnames(select_k_e),"co")
    justint <- transform(select_k_var[intdata_k_id %in% intju, ], co = 2)
    colnames(justint)=c(colnames(select_k_e),"co")
    justext <- transform(select_k_e[extdata_id %in% extju, ], co = 3)
    colnames(justext)=c(colnames(select_k_e),"co")
    
    combdata=rbind(justint,bothie,justext)
    
    combdata$group=relevel(as.factor(combdata$co),ref=1)
    formula_string <- paste("group ~", 
                            paste(colnames(select_k_e), collapse = " + "))
    formula <- as.formula(formula_string)
    model1=multinom(formula, data=combdata)
    
    modsimp[[k]]=modsimp_k
    
    modmult[[k]]=model1
  } 
  
  comb_est_weights_list=list()
  for(i in 1:K){
    select_i_i=Select_list[[i]][[i]]
    estweights_list[[i]]=matrix(0,nrow(select_i_i),K)
    for(j in 1:K){
      select_i_j=Select_list[[i]][[j]]
      wtintsimp=predict(modsimp[[j]],newdata=select_i_j,type="response")
      Pmult1=predict(modmult[[j]],select_i_j,type="prob")
      
      prob=rep(0,times=nrow(select_i_j))
      for(h in 1:nrow(select_i_j)){
        prob[h]= wtintsimp[h] * (Pmult1[h,1] + Pmult1[h,2])/
          (Pmult1[h,1] + Pmult1[h,3])
      }
      
      prob[prob>1]=1
      estweights_list[[i]][,j]=1/prob
    }
    comb_est_weights=rep(1,nrow(select_i_i))
    for(l in 1:K){
      comb_est_weights=comb_est_weights*(1-1/estweights_list[[i]][,l])
    }
    comb_est_weights_list[[i]]=1/(1-comb_est_weights)
  }
  return(list(combined_weights = comb_est_weights_list))
}



###
library(nleqslv)
library(nloptr)
library(MASS)
library(xgboost)
expit<-function(x){
  return(exp(x)/(1+exp(x)))
}
K=3
set.seed(100)
mean_w_p=0
mean_z_1=0
mean_z_2=0
mean_z_3=0
corr=0.5
var_z_w_p=matrix(c(1,corr,corr,corr,
                   corr,1,corr,corr,
                   corr,corr,1,corr,
                   corr,corr,corr,1),
                 nrow=4,ncol=4)

theta=c(-2,0.35,0.45,0.25)
N=5e4
dw=1
dwz1=c(1,0.8,0.6)
dwz2=c(0.6,0.8,1)
dwz3=rep(1,3)

gamma_ext=c(-0.6,1.2,0.4,-0.2,0.5)
gamma_int_1=c(-1,1.5,0.2,0.8,-0.3)
gamma_int_2=c(-1,1.25,0.4,0.6)
gamma_int_3=c(-3,0.8,0.5)

simu_popu<-function(N,mean_w_p,mean_z_1,mean_z_2,mean_z_3,
                    var_z_w_p,theta,dw){
  cov<- mvrnorm(n = N, mu = c(mean_w_p,mean_z_1,mean_z_2,mean_z_3), Sigma = var_z_w_p)
  data <- data.frame(Z1 = cov[, 2], Z2 = cov[, 3], Z3=cov[,4])
  W_p=cov[,1]
  # Generate random uniforms
  #set.seed(5678)
  U1 <- runif(N)
  #set.seed(4321)
  # Generate Disease Status
  DISEASE <- expit(theta[1] + theta[2] * data$Z1 + theta[3]*data$Z2 +theta[4]*data$Z3)
  data$D   <- ifelse(DISEASE > U1, 1, 0)
  # Relate W_p and D
  data$W_1 <- W_p + dw* data$D + dwz1[1]*data$Z1 + 
    dwz2[1]*data$Z2 + dwz3[1]*data$Z3 +
    rnorm(n=N,0,1)
  
  data$W_2 <- W_p + dw* data$D + dwz1[2]*data$Z1 + 
    dwz2[2]*data$Z2 + dwz3[2]*data$Z3 +
    rnorm(n=N,0,1)
  
  data$W_3 <- W_p + dw* data$D + dwz1[3]*data$Z1 + 
    dwz2[3]*data$Z2 + dwz3[3]*data$Z3 +
    rnorm(n=N,0,1)
  
  data$id=c(1:N)
  return(data)
}

simu_ext<-function(data,gamma_ext){
  U2e <- runif(N)
  # Generate Sampling Status
  SELECT <-0.75*expit(gamma_ext[1] +  
                        gamma_ext[2]* data$D + 
                        gamma_ext[3] * data$Z1 +
                        gamma_ext[4]* data$Z2 +
                        gamma_ext[5] * data$Z3)
  S_e  <- ifelse(SELECT > U2e, T, F)
  # Observed Data
  data_e <- data[which(S_e==1),]
  data_e$Select_Weights = 0.75*expit(gamma_ext[1] +  
                                       gamma_ext[2]* data_e$D + 
                                       gamma_ext[3] * data_e$Z1 +
                                       gamma_ext[4]* data_e$Z2 +
                                       gamma_ext[5] * data_e$Z3)
  return(data_e)
}

simu_int_1<-function(data,gamma_int_1){
  U2i <- runif(N)
  # Generate Sampling Status
  SELECT <- expit(cbind(1,data$D,data$W_1,data$Z2,data$Z3)
                  %*% gamma_int_1)
  S_i  <- ifelse(SELECT > U2i, T, F)
  # Observed Data
  data_i <- data[which(S_i==1),]
  return(data_i)
}

simu_int_2<-function(data,gamma_int_2){
  U2i <- runif(N)
  # Generate Sampling Status
  SELECT <- expit(cbind(1,data$D,data$W_2,data$Z3)
                  %*% gamma_int_2)
  S_i  <- ifelse(SELECT > U2i, T, F)
  # Observed Data
  data_i <- data[which(S_i==1),]
  return(data_i)
}



simu_int_3<-function(data,gamma_int_3){
  U2i <- runif(N)
  # Generate Sampling Status
  SELECT <- expit(cbind(1,data$W_3,data$Z2)
                  %*% gamma_int_3)
  S_i  <- ifelse(SELECT > U2i, T, F)
  # Observed Data
  data_i <- data[which(S_i==1),]
  return(data_i)
}
data=simu_popu(N,mean_w_p,mean_z_1,mean_z_2,mean_z_3,
               var_z_w_p,theta,dw)
extdata=simu_ext(data,gamma_ext)
intdata_comb=NULL
intdata1=simu_int_1(data,gamma_int_1)
intdata2=simu_int_2(data,gamma_int_2)
intdata3=simu_int_3(data,gamma_int_3)
intdata_comb=rbind(intdata1,intdata2,intdata3)

intdata_comb1=intdata_comb[!duplicated(intdata_comb$id),]

## 
library(dplyr)
select_1=list(intdata1 %>%
                dplyr::select("D","W_1","Z2","Z3"),
              intdata1 %>%
                dplyr::select("D","W_2","Z3"),
              intdata1 %>%
                dplyr::select("W_3","Z2"))


select_2=list(intdata2 %>%
                dplyr::select("D","W_1","Z2","Z3"),
              intdata2 %>%
                dplyr::select("D","W_2","Z3"),
              intdata2 %>%
                dplyr::select("W_3","Z2"))



select_3=list(intdata3 %>%
                dplyr::select("D","W_1","Z2","Z3"),
              intdata3 %>%
                dplyr::select("D","W_2","Z3"),
              intdata3 %>%
                dplyr::select("W_3","Z2"))
select_e=list(extdata %>%
                dplyr::select("D","W_1","Z2","Z3"),
              extdata %>%
                dplyr::select("D","W_2","Z3"),
              extdata %>%
                dplyr::select("W_3","Z2"))
select=list(select_1,select_2,select_3)

K=3

res=SR_est_multi_cohort(K=3,Select_list=select,
                                    Select_e_list=select_e,
                                    Weights_e=1/extdata$Select_Weights,
                                    intdata_id_list=
                                    list(intdata1$id,intdata2$id,intdata3$id),
                                    extdata_id=extdata$id)

comb_est_weights=NULL
for(i in 1:K){
  comb_est_weights=c(comb_est_weights,res$combined_weights[[i]])
}
comb_est_weights1=comb_est_weights[!duplicated(intdata_comb$id)]


weighted<-function(intdata,estweights){
  modelinit=glm(D ~ Z1 +Z2 +Z3, family = stats::quasibinomial(),data=intdata)
  start <- coef(modelinit)
  design <- survey::svydesign(data = intdata,ids = 1:length(intdata$D), strata = NULL,
                              weights = estweights)
  mod=survey::svyglm(data=intdata,stats::formula(D ~ Z1 + Z2 +Z3), design = design,
                     family = quasibinomial(),
                     start = as.numeric(start))
  final<-coef(mod)
  return(final)
}
weighted(intdata_comb1,comb_est_weights1)


