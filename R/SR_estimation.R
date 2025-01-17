SR_wt_multi_cohort <- function(K,intdata_list,extdata,select_var_list,Weights_e){
  
  # expit<-function(x){
  #   return(exp(x)/(1+exp(x)))
  # }

  
  # Initialize lists to store results for each cohort
  estweights_list <- list()
  modsimp <-list()
  modmult <-list()
  
  # Loop over each cohort
  for (k in 1:K){
    intdata=intdata_list[[k]]
    
    # Check select_k_var and select_k_e to be data frames or not
    if (!is.data.frame(intdata)) stop("intdata must be a data frame.")
    
    intdata_select=intdata%>%
      dplyr::select(select_var_list[[k]])
    
    extdata_select=extdata%>%
      dplyr::select(select_var_list[[k]])
    
    formula_string <- paste("1/Weights_e ~", 
                            paste(colnames(extdata_select), collapse = " + "))
    
    # Convert the string to a formula
    formula <- as.formula(formula_string)
    
    modsimp_k=simplexreg(data=data.frame(extdata),
                         formula= formula,
                         link="logit")
    
    # Find common and exclusive IDs
    inter <- intersect(intdata$id, extdata$id)
    intju <- setdiff(intdata$id, inter)
    extju <- setdiff(extdata$id, inter)
    
    # Subset the data and add the 'co' column in one step
    bothie <- transform(extdata_select[extdata$id %in% inter, ], co = 1)
    colnames(bothie)=c(colnames(extdata_select),"co")
    justint <- transform(intdata_select[intdata$id %in% intju, ], co = 2)
    colnames(justint)=c(colnames(extdata_select),"co")
    justext <- transform(extdata_select[extdata$id %in% extju, ], co = 3)
    colnames(justext)=c(colnames(extdata_select),"co")
    
    combdata=rbind(justint,bothie,justext)
    
    combdata$group=relevel(as.factor(combdata$co),ref=1)
    formula_string <- paste("group ~", 
                            paste(colnames(extdata_select), collapse = " + "))
    formula <- as.formula(formula_string)
    model1=multinom(formula, data=combdata)
    
    modsimp[[k]]=modsimp_k
    
    modmult[[k]]=model1
  } 
  
  comb_est_weights_list=list()
  for(i in 1:K){
    select_i_i=intdata_list[[i]]%>%
      dplyr::select(select_var_list[[i]])
    
    estweights_list[[i]]=matrix(0,nrow(select_i_i),K)
    for(j in 1:K){
      select_i_j=intdata_list[[i]]%>%
        dplyr::select(select_var_list[[j]])
      
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

# source("weighted.R")

SR_est_var_multi_cohort <- function(K,intdata_list,extdata,
                                    select_var_list,Weights_e,Z_names,
                                    intdata_comb_boot=NULL,
                                    AIPW=NULL,type_var=NULL,variance=NULL){
  wts_fun=SR_wt_multi_cohort(K,intdata_list,extdata,select_var_list,Weights_e)
  wts=wts_fun$combined_weights
  intdata_comb=NULL
  comb_est_weights=NULL
  for(i in 1:K){
    intdata=intdata_list[[i]]
    if (!is.data.frame(intdata)) stop("All Internal data should be a data frame.")
    intdata_comb=rbind.data.frame(intdata_comb,intdata)
    comb_est_weights=c(comb_est_weights,wts[[i]])
  }
  intdata_comb1=intdata_comb[!duplicated(intdata_comb$id),]
  comb_est_weights1=comb_est_weights[!duplicated(intdata_comb$id)]
  
  fin_est=weighted(intdata=intdata_comb1,estweights=comb_est_weights1,Z_names=Z_names)
  
  return(list(combined_weights = wts,
              final_est=fin_est$final,
              variance_est=fin_est$var))
}


