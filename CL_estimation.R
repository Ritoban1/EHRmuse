CL_est_multi_cohort <- function(K,Select_list,marginals_list){
  
  if (!is.list(Select_list) || !is.list(marginals_list)) 
    stop("All selection input data should be lists for 
         multiple cohorts and the external data.")
  
  if (length((Select_list))!=K || length((marginals_list))!=K) 
    stop("All input lists should have length K.")
  

  # Initialize lists to store results for each cohort
  estweights_list <- list()
  gamma_estimates <- list()
  
  # Define expit function if not already defined
  expit <- function(x) 1 / (1 + exp(-x))
  
  for (k in 1:K){
    select_k=Select_list[[k]]
    
    if (!is.list(select_k)) 
      stop("Selection_variable_input_for_each_cohort_should_be_a_list")
    
    if (length((select_k))!=K) 
      stop("Selection_variable_input_for_each_cohort_should_be_a_list of length K")
    
    ## Selection variables exclusive to cohort K
    
    select_k_var=select_k[[k]]
    
    ##marginals for cohort K
    marg_k=marginals_list[[k]]
    
    # Convert select_k_var to matrices if necessary
    if (!is.data.frame(select_k_var)) stop("select_k_var must be a data frame.")
    
    ## marg_k should be a numeric vector
    if (!is.numeric(marg_k)) stop("marg_k must be a numeric vector.")
    
    
    prop<-function(gamma){
      y <- c(rep(0,(ncol(select_k_var)+1)))
      for(i in 1:nrow(select_k_var)){
        vec=c(1, as.numeric(select_k_var[i,]))
        y = y + 1/(as.vector(expit(gamma %*% vec))) * vec 
      }
      y= y - marg_k
      y
    }
    
    # Set starting values for gamma based on dimensions
    start <- c(rep(0,(ncol(select_k_var)+1)))
    
    # Optimization
    z <- nleqslv(x = start, fn = prop, method = "Newton", global = "dbldog",
                 control = list(trace = 1, allowSingular = TRUE))
    
    # gamma_values_for_cohort_k
    gamma_estimates[[k]] <- z$x
  }
  
  comb_est_weights_list=list()
  
  for(i in 1:K){
    select_i_i=Select_list[[i]][[i]]
    estweights_list[[i]]=matrix(0,nrow(select_i_i),K)
    for(j in 1:K){
      select_i_j=Select_list[[i]][[j]]
      estweights_list[[i]][,j]= 1/expit(as.matrix(cbind(1,select_i_j)) %*% 
                                          as.numeric(gamma_estimates[[j]]))
    }
    comb_est_weights=rep(1,nrow(select_i_i))
    for(l in 1:K){
      comb_est_weights=comb_est_weights*(1-1/estweights_list[[i]][,l])
    }
    comb_est_weights_list[[i]]=1/(1-comb_est_weights)
  }  
  
  return(list(gamma_estimates = gamma_estimates, 
              combined_weights = comb_est_weights_list))
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

totals<-function(data){
  total_D=sum(data$D)
  total_W_1=sum(data$W_1)
  total_W_2=sum(data$W_2)
  total_W_3=sum(data$W_3)
  total_Z2=sum(data$Z2)
  total_Z3=sum(data$Z3)
  return(c(nrow(data),total_D,total_W_1,
           total_W_2,total_W_3,total_Z2,total_Z3))
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


select=list(select_1,select_2,select_3)

marginal=totals(data)
margs=list(marginal[c(1,2,3,6,7)],
           marginal[c(1,2,4,7)],
           marginal[c(1,5,6)])

K=3



res=CL_est_multi_cohort(K,Select_list=select,
                        marginals_list=margs)


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
