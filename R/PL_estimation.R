PL_wt_multi_cohort <- function(K,intdata_list,extdata,select_var_list,Weights_e){
  # expit<-function(x){
  #   return(exp(x)/(1+exp(x)))
  # }
  # Initialize lists to store results for each cohort
  estweights_list <- list()
  gamma_estimates <- list()
  
  # Define expit function if not already defined
  # expit <- function(x) 1 / (1 + exp(-x))
  
  # Loop over each cohort
  for (k in 1:K){
    intdata=intdata_list[[k]]
    
    # Check select_k_var and select_k_e to be data frames or not
    if (!is.data.frame(intdata)) stop("intdata must be a data frame.")
    
    intdata_select=intdata%>%
      dplyr::select(select_var_list[[k]])
    
    extdata_select=extdata%>%
      dplyr::select(select_var_list[[k]])
    
    prop <- function(gamma) {
      y <- c(rep(0,(ncol(intdata_select)+1)))
      vec1=c(nrow(intdata_select),colSums(intdata_select))
      
      for(i in 1:nrow(extdata_select)){
        vec=c(1,as.numeric(extdata_select[i,]))
        y = y + as.vector(expit(gamma %*% vec)) * vec *
          Weights_e[i]
      }
      y <- vec1 - y
      y
    }
    
    # Set starting values for gamma based on dimensions
    start <- c(rep(0,(ncol(intdata_select)+1)))
    
    # Optimization
    z <- nleqslv(x = start, fn = prop, method = "Newton", global = "dbldog",
                 control = list(trace = 1, allowSingular = TRUE))
    
    # gamma_values_for_cohort_k
    gamma_estimates[[k]] <- z$x
  }
  
  comb_est_weights_list=NULL
  
  for(i in 1:K){
    select_i_i=intdata_list[[i]]%>%
      dplyr::select(select_var_list[[i]])
    
    estweights_list[[i]]=matrix(0,nrow(select_i_i),K)
    for(j in 1:K){
      select_i_j=intdata_list[[i]]%>%
        dplyr::select(select_var_list[[j]])
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

# source("weighted.R")

PL_est_multi_cohort <- function(K,intdata_list,extdata,select_var_list,
                                Weights_e,Z_names){
  wts_fun=PL_wt_multi_cohort(K,intdata_list,extdata,select_var_list,Weights_e)
  wts=wts_fun$combined_weights
  intdata_comb=NULL
  comb_est_weights=NULL
  for(i in 1:K){
    intdata=intdata_list[[i]]
    intdata_comb=rbind.data.frame(intdata_comb,intdata)
    comb_est_weights=c(comb_est_weights,wts[[i]])
  }
  intdata_comb1=intdata_comb[!duplicated(intdata_comb$id),]
  comb_est_weights1=comb_est_weights[!duplicated(intdata_comb$id)]
  
  fin_est=weighted(intdata=intdata_comb1,estweights=comb_est_weights1,Z_names=Z_names)
  
  return(list(gamma_estimates = wts_fun$gamma_estimates, 
              combined_weights = wts,
              final_est=fin_est$final,
              variance_est=fin_est$var))
}



