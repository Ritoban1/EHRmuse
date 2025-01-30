AIPW_var_approx<-function(K,N,intdata_list,extdata,
                          select_var_list,
                          aux_var_list,
                          Weights_e,Z_names,
                          aux_model="XGBoost",
                          ipw_method="PL",
                          marginals_list=NULL,
                          weights_user=NULL){
  # source("AIPW_estimation.R")
  
  est_AIPW=AIPW_joint(K=K,N=N,intdata_list=intdata_list,
                      extdata=extdata,
                      select_var_list=select_var_list,
                      aux_var_list=aux_var_list,
                      Weights_e=Weights_e,Z_names=Z_names,
                      aux_model=aux_model,ipw_method=ipw_method,
                      marginals_list=marginals_list,
                      weights_user=weights_user)
  
  
  
  wts=est_AIPW$estweights_est
  intdata_comb=NULL
  comb_est_weights=NULL
  start=est_AIPW$start
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
    
    ### part 1
    sum_1=matrix(0,N_Z,N_Z)
    
    for(i in 1:length(D_int)){
      vec1=as.numeric(c(1, Z_1_int[i,],Z_n_1_int[i,]))
      vec2=as.numeric(c(intpred2[i],intpred1[i,],
                        intpred2[i]*Z_n_1_int[i,]))
      vec3=((D_int[i]*vec1-
               as.vector(expit(theta_vec %*% vec1))*vec1)-
              vec2)
      sum_1 = sum_1 + vec3 %*% t(vec3)*(estweights_est[i])^2
    }
    
    #### external part
    sum_4=matrix(0,N_Z,N_Z)
    
    for(i in 1:length(D_ext)){
      vec=c(as.numeric(extpred2[i]),as.numeric(extpred1[i,]),
            as.vector(extpred2[i]*as.numeric(Z_n_1_ext[i,])))
      sum_4 = sum_4 + (Weights_e[i])^2*vec %*% t(vec)
    }
    
    ### common part
    
    common_id=intersect(intdata_comb1$id,extdata$id)
    
    if(sum(is.na(common_id))==0){
      sum_2=matrix(0,N_Z,N_Z)
      ind_common=which(intdata_comb1$id %in% common_id)
      D_com=D_int[ind_common]
      intpred1_com=intpred1[ind_common]
      intpred2_com=intpred2[ind_common]
      Z_1_com=data.frame(Z_1_int[ind_common,])
      Z_n_1_com=data.frame(Z_n_1_int[ind_common,])
      estweights_est_com=estweights_est[ind_common]
      
      for(i in 1:length(D_com)){
        vec1=as.numeric(c(1,Z_1_com[i,],Z_n_1_com[i,]))
        vec2=as.numeric(c(intpred2_com[i],intpred1_com[i],
                          intpred2_com[i]*Z_n_1_com[i,]))
        vec3=((D_com[i]*vec1-
                 as.vector(expit(theta_vec %*% vec1))*vec1)-
                vec2)
        sum_2=sum_2 + vec3%*%t(vec2)*estweights_est_com[i]*Weights_e[i]
      }
      
      sum_3=t(sum_2)
      
      E=sum_1+sum_4+sum_2+sum_3
      
    }else{
      E=sum_1+sum_4
    }
    return(E/N)
  }
  
    
  E=prop_xg(est_AIPW$est)
  
  G_theta=(est_AIPW$jac)
  
  var_est=1/N *solve(G_theta) %*% E %*% t(solve(G_theta))
  
  for(i in 1:nrow(var_est)){
    if(var_est[i,i]<0){
      var_est[i,i]=(est_AIPW$init_var)[i]
    }
  }
  return(list(est=est_AIPW$est,var=diag(var_est)))
}
