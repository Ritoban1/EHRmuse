AIPW_joint<-function(D_int,Z_1_int,
                     Z_n_1_int,W_int,
                     D_ext,
                     Z_n_1_ext,W_ext,
                     Weights_e,
                     estweights_est,aux_model,start){
  
  if (!is.numeric(D_int) || !is.vector(D_int))
    stop("'D_int must be a numeric vector.")
  
  if (!is.numeric(D_ext) || !is.vector(D_ext))
    stop("'D_ext must be a numeric vector.")
  
  # Check select_k_var and select_k_e to be data frames or not
  if (!is.data.frame(Z_1_int)) stop("Z_1_int must be a data frame.")
  
  if (!is.data.frame(Z_n_1_int)) stop("Z_n_1_int must be a data frame.")
  
  if (!is.data.frame(W_int)) stop("W_int must be a data frame.")
  
  if (!is.data.frame(Z_n_1_ext)) stop("Z_n_1_ext must be a data frame.")
  
  if (!is.data.frame(W_ext)) stop("W_ext must be a data frame.")
  
  ## Estimated weights for the combined internal data
  if (!is.numeric(estweights_est) || !is.vector(estweights_est))
    stop("estweights_est must be a numeric vector.")
  
  ## Sample weights for the external data
  if (!is.numeric(Weights_e) || !is.vector(Weights_e))
    stop("estweights_est must be a numeric vector.")
  
  if(aux_model=="XGBoost"){
    
     extpred1=matrix(0,length(D_ext),ncol(Z_1_int))
     intpred1=matrix(0,length(D_int),ncol(Z_1_int))
     
     prop_xg<-function(theta_vec){
      N_Z=1+ncol(Z_1_int)+ncol(Z_n_1_int)
      y1=c(rep(0,times=N_Z))
      for(d in 1:ncol(Z_1_int)){
        response1=Z_1_int[,d]*(D_int-as.vector(expit(as.matrix(cbind(1,
                                                    Z_1_int,
                                      Z_n_1_int)) %*% theta_vec)))
        
        xgb_int1=xgb.DMatrix(data=data.matrix(cbind(D_int,Z_n_1_int,
                                                    W_int)),
                             label=response1)
        colnames(xgb_int1)[1]="D"
        xgbc1=xgboost(data=xgb_int1,max.depth=3,nrounds=100,verbose = 0)
        
        xgb_ext=xgb.DMatrix(data=data.matrix(cbind(D_ext,Z_n_1_ext,
                                                   W_ext)))
        colnames(xgb_ext)[1]="D"
        extpred1[,d]=predict(xgbc1,xgb_ext)
        intpred1[,d]=predict(xgbc1,xgb_int1)
        
      }
      
      response2=(D_int-as.vector(expit(as.matrix(cbind(1,
                                                       Z_1_int,
                                    Z_n_1_int)) %*% theta_vec)))
      
      
      xgb_int2=xgb.DMatrix(data=data.matrix(cbind(D_int,Z_n_1_int,
                                                  W_int)),
                           label=response2)
      colnames(xgb_int2)[1]="D"
      xgbc2=xgboost(data=xgb_int2,max.depth=3,nrounds=100,verbose = 0)
      xgb_ext=xgb.DMatrix(data=data.matrix(cbind(D_ext,Z_n_1_ext,
                                                 W_ext)))
      colnames(xgb_ext)[1]="D"
      extpred2=predict(xgbc2,xgb_ext)
      intpred2=predict(xgbc2,xgb_int2)
      
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
                 global = "dbldog",control=list(trace=1,allowSingular=TRUE))
     return(z$x)
      
    }
}
 
estweights_est=comb_est_weights1

D_int=intdata_comb1$D
Z_1_int=data.frame(intdata_comb1$Z1)
Z_n_1_int=intdata_comb1%>%
  dplyr::select("Z2","Z3")

W_int=intdata_comb1%>%
  dplyr::select("W_1","W_2","W_3")

D_ext=extdata$D

Z_n_1_ext=extdata%>%
  dplyr::select("Z2","Z3")

W_ext=extdata%>%
  dplyr::select("W_1","W_2","W_3")
Weights_e=1/extdata$Select_Weights
aux_model="XGBoost"
start=weighted(intdata_comb1,comb_est_weights1)
AIPW_joint(D_int,Z_1_int,
                     Z_n_1_int,W_int,
                     D_ext,
                     Z_n_1_ext,W_ext,
                     Weights_e,
                     estweights_est,aux_model,start)


