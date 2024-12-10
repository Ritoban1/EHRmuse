AIPW_joint<-function(K,N,intdata_list,extdata,
                     select_var_list,
                     aux_var_list,
                     Weights_e,Z_names,
                     aux_model="XGBoost",ipw_method="PL",
                     marginals_list=NULL,
                     weights_user=NULL){
 
  if(ipw_method=="PL"){
    # source("PL_estimation.R")
    
    IPW_est=PL_est_multi_cohort(K,intdata_list,extdata,
                                select_var_list,Weights_e,Z_names)

    wts=IPW_est$combined_weights
    intdata_comb=NULL
    comb_est_weights=NULL
    start=IPW_est$final_est
    select_var_comb_name=NULL
    aux_var_comb_name=NULL
    for(i in 1:K){
      intdata=intdata_list[[i]]
      if (!is.data.frame(intdata)) stop("All Internal data should be a data frame.")
      intdata_comb=rbind.data.frame(intdata_comb,intdata)
      comb_est_weights=c(comb_est_weights,wts[[i]])
      select_var_comb_name=union(select_var_comb_name,select_var_list[[i]])
      aux_var_comb_name=union(aux_var_comb_name,aux_var_list[[i]])
    }
    intdata_comb1=intdata_comb[!duplicated(intdata_comb$id),]
    estweights_est=comb_est_weights[!duplicated(intdata_comb$id)]
    
    
  }else if(ipw_method=="SR"){
    # source("SR_estimation.R")
    
    IPW_est=SR_est_var_multi_cohort(K,intdata_list,extdata,
                                    select_var_list,Weights_e,Z_names)

    wts=IPW_est$combined_weights
    intdata_comb=NULL
    comb_est_weights=NULL
    start=IPW_est$final_est
    select_var_comb_name=NULL
    aux_var_comb_name=NULL
    for(i in 1:K){
      intdata=intdata_list[[i]]
      if (!is.data.frame(intdata)) stop("All Internal data should be a data frame.")
      intdata_comb=rbind.data.frame(intdata_comb,intdata)
      comb_est_weights=c(comb_est_weights,wts[[i]])
      select_var_comb_name=union(select_var_comb_name,select_var_list[[i]])
      aux_var_comb_name=union(aux_var_comb_name,aux_var_list[[i]])
    }
    intdata_comb1=intdata_comb[!duplicated(intdata_comb$id),]
    estweights_est=comb_est_weights[!duplicated(intdata_comb$id)]
    

  }else if(ipw_method=="CL"){
    # source("CL_estimation.R")
    IPW_est=CL_est_multi_cohort(K,intdata_list,select_var_list,
                                marginals_list,Z_names)
    wts=IPW_est$combined_weights
    intdata_comb=NULL
    comb_est_weights=NULL
    start=IPW_est$final_est
    select_var_comb_name=NULL
    
    aux_var_comb_name=NULL
    for(i in 1:K){
      intdata=intdata_list[[i]]
      if (!is.data.frame(intdata)) stop("All Internal data should be a data frame.")
      intdata_comb=rbind.data.frame(intdata_comb,intdata)
      comb_est_weights=c(comb_est_weights,wts[[i]])
      select_var_comb_name=union(select_var_comb_name,select_var_list[[i]])
      aux_var_comb_name=union(aux_var_comb_name,aux_var_list[[i]])
    }
    intdata_comb1=intdata_comb[!duplicated(intdata_comb$id),]
    estweights_est=comb_est_weights[!duplicated(intdata_comb$id)]
  }else if(ipw_method=="US"){
    estweights_est=weights_users
    intdata_comb=NULL
    aux_var_comb_name=NULL
    for(i in 1:K){
      intdata=intdata_list[[i]]
      if (!is.data.frame(intdata)) stop("All Internal data should be a data frame.")
      intdata_comb=rbind.data.frame(intdata_comb,intdata)
      aux_var_comb_name=union(aux_var_comb_name,aux_var_list[[i]])
    }
    intdata_comb1=intdata_comb[!duplicated(intdata_comb$id),]
    weight_us=weighted(intdata=intdata_comb1,
                       estweights=estweights_est,Z_names=Z_names)
    start=weight_us$final
  }
  
  
  D_int=intdata_comb1$D
  D_ext=extdata$D
  Z_1_int_name=setdiff(Z_names,aux_var_comb_name)
  Z_n_1_name=setdiff(Z_names,Z_1_int_name)
  W_name=setdiff(aux_var_comb_name,Z_names)
  
  Z_1_int=intdata_comb1%>%
    dplyr::select(Z_1_int_name)
  Z_n_1_int=intdata_comb1%>%
    dplyr::select(Z_n_1_name)
  Z_n_1_ext=extdata%>%
    dplyr::select(Z_n_1_name)
  
  W_int=intdata_comb1%>%
    dplyr::select(W_name)
  W_ext=extdata%>%
    dplyr::select(W_name)
  
 if(aux_model=="XGBoost"){
    
    prop_xg<-function(theta_vec){
      extpred1=matrix(0,length(D_ext),ncol(Z_1_int))
      intpred1=matrix(0,length(D_int),ncol(Z_1_int))
      N_Z=1+ncol(Z_1_int)+ncol(Z_n_1_int)
      y1=c(rep(0,times=N_Z))
      for(d in 1:ncol(Z_1_int)){
        response1=Z_1_int[,d]*(D_int-as.vector(expit(as.matrix(cbind(1,
                                                    Z_1_int,
                                      Z_n_1_int)) %*% theta_vec)))
        
        xgb_int1=xgb.DMatrix(data=data.matrix(cbind(Z_n_1_int,
                                                    W_int)),
                             label=response1)
        xgbc1=xgboost(data=xgb_int1,max.depth=3,nrounds=100,verbose = 0)
        
        xgb_ext=xgb.DMatrix(data=data.matrix(cbind(Z_n_1_ext,
                                                   W_ext)))
  
        extpred1[,d]=predict(xgbc1,xgb_ext)
        intpred1[,d]=predict(xgbc1,xgb_int1)
        
      }
      
      response2=(D_int-as.vector(expit(as.matrix(cbind(1,
                                                       Z_1_int,
                                    Z_n_1_int)) %*% theta_vec)))
      
      
      xgb_int2=xgb.DMatrix(data=data.matrix(cbind(Z_n_1_int,
                                                  W_int)),
                           label=response2)
      xgbc2=xgboost(data=xgb_int2,max.depth=3,nrounds=100,verbose = 0)
      xgb_ext=xgb.DMatrix(data=data.matrix(cbind(Z_n_1_ext,
                                                 W_ext)))
      extpred2=predict(xgbc2,xgb_ext)
      intpred2=predict(xgbc2,xgb_int2)
      

      
      #####
      
      for(i in 1:length(D_ext)){
        vec=c(as.numeric(extpred2[i]),as.numeric(extpred1[i,]),
                as.vector(extpred2[i]*as.numeric(Z_n_1_ext[i,])))
        y1 = y1 + vec*Weights_e[i]
      }
      
      y2=c(rep(0,times=N_Z))
      
      for(i in 1:length(D_int)){
        vec1=as.numeric(c(1, Z_1_int[i,],Z_n_1_int[i,]))
        vec2=as.numeric(c(intpred2[i],intpred1[i,],
                          intpred2[i]*Z_n_1_int[i,]))
        y2 = y2 + ((D_int[i]*vec1-
                    as.vector(expit(theta_vec %*% vec1))*vec1)-
                     vec2)*estweights_est[i]
        
      }
      y=y1+y2
      y/N
     }
     start=as.numeric(start)
     z = nleqslv(x=start, fn=prop_xg,method="Newton",
                 global = "dbldog",control=list(trace=1,allowSingular=TRUE),
                 jacobian = TRUE)
     return(list(est=z$x,jac=z$jac,estweights_est=wts,start=start,
                 init_var=IPW_est$variance_est))
      
  }
  
  
}
