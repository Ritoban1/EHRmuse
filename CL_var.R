variance_cl_actual<-function(K,N,intdata_list,
                             select_var_list,
                             marginals_list,Z_names){
  
  wts_fun=CL_est_multi_cohort(K,intdata_list,select_var_list,
                              marginals_list,Z_names)
  wts=wts_fun$combined_weights
  N_est=N
  theta_hat=as.numeric(wts_fun$final_est)
  intdata_comb=NULL
  comb_est_weights=NULL
  ind_list=NULL
  
  for(i in 1:K){
    intdata=intdata_list[[i]]
    intdata_comb=rbind.data.frame(intdata_comb,intdata)
    comb_est_weights=c(comb_est_weights,wts[[i]])
  }
  
  intdata_comb1=intdata_comb[!duplicated(intdata_comb$id),]
  comb_est_weights1=comb_est_weights[!duplicated(intdata_comb$id)]
  estweights=comb_est_weights1
  
  for(i in 1:K){
    intdata=intdata_list[[i]]
    ind_list[[i]]=as.numeric(intdata_comb1$id %in% intdata$id)
  }
  
  estweights_indi=NULL
  for(j in 1:K){
    intdata_comb_select=intdata_comb1 %>%
      dplyr::select(select_var_list[[j]])
    
    intdata_comb_select=cbind.data.frame(1,intdata_comb_select)
    
    prob_est=expit(as.matrix(intdata_comb_select) %*%  
                     wts_fun$gamma_estimates[[j]])
    estweights_indi[[j]]=1/prob_est
  }
  
  
  ## g_theta
  z_len=length(Z_names)+1
  
  g_theta = matrix(0,nrow=z_len,ncol=z_len)
  intprob = 1/comb_est_weights1
  intdata_Z=intdata_comb1%>%
    dplyr::select(all_of(Z_names))
  
  
  for(i in 1:nrow(intdata_Z)){
    z=as.numeric(c(1,intdata_Z[i,]))
    g_theta = g_theta + as.numeric(1/intprob[i] *
                                     exp(theta_hat%*%z)/(1+ exp(theta_hat%*%z))^2) * 
      z%*%t(z)
  }
  g_theta=g_theta/N_est
  
  ## g_alpha
  
  x_dim=K
  
  for(i in 1:K){
    x_dim=x_dim+length(select_var_list[[i]])
  }
  
  g_alpha = matrix(0,nrow=z_len,ncol=x_dim)
  for(i in 1:nrow(intdata_comb1)){
    z=as.numeric(c(1,intdata_Z[i,]))
    
    temp_alpha=NULL
    for(j in 1:K){
      intdata_comb_select=intdata_comb1 %>%
        dplyr::select(select_var_list[[j]])
      select_x=as.numeric(c(1,intdata_comb_select[i,]))
      temp_alpha=cbind(temp_alpha,1/estweights_indi[[j]][i]*t(select_x))
    }
    
    temp=1
    
    for(j in 1:K){
      temp=temp * (1-1/(estweights_indi[[j]][i]))
    }
    
    g_alpha=g_alpha + estweights[i]^2*temp*
      as.numeric(intdata_comb1$D[i] - expit(theta_hat%*%z))* z %*% temp_alpha
    
  }
  g_alpha = -g_alpha/N_est
  
  ## M
  H = matrix(0,nrow=(x_dim),ncol=(x_dim))
  
  cum_length=rep(0,K)
  len_select=rep(0,K)
  for(j in 1:K){
    len_select[j]=length(select_var_list[[j]])+1
  }
  cum_length=cumsum(len_select)
  
  for(i in 1:nrow(intdata_comb1)){
    temp=matrix(0,nrow=x_dim,ncol=(x_dim))
    x_1=intdata_comb1%>%
      dplyr::select(select_var_list[[1]])
    x_1=as.numeric(c(1,x_1[i,]))
    temp[c(1:length(x_1)),c(1:length(x_1))]=ind_list[[1]][i]*as.numeric((estweights_indi[[1]][i])
                                                       *(1-1/(estweights_indi[[1]][i])))*
                                           x_1%*%t(x_1)
    if(K>1){
      for(j in 2:K){
        x=intdata_comb1 %>%
          dplyr::select(select_var_list[[j]])
        x=as.numeric(c(1,x[i,]))
        
        temp[c((cum_length[j-1]+1):(cum_length[j])),
             c((cum_length[j-1]+1):(cum_length[j]))]=
          ind_list[[j]][i]*as.numeric((estweights_indi[[j]][i])
                                      *(1-1/(estweights_indi[[j]][i])))*
          x%*%t(x)
      }
    }
    
    H=H+temp
  }
  H= -H/N_est
  
  ## E1
  
  E1 = matrix(0,nrow=z_len,ncol=z_len)
  
  for(i in 1:nrow(intdata_comb1)){
    z=as.numeric(c(1,intdata_Z[i,]))
    g = as.numeric(1/intprob[i] * (intdata_comb1$D[i] - 
                                     exp(theta_hat%*%z)/(1+ exp(theta_hat%*%z)))) * z
    E1 = E1 + g %*% t(g)
  }
  E1 = E1/N_est
  
  ## E2
  E2 = matrix(0,nrow=z_len,ncol=z_len)
  sum1 = matrix(0,nrow=x_dim,ncol=z_len)
  
  
  for(i in 1:nrow(intdata_comb1)){
    x=NULL
    for(j in 1:K){
      intdata_comb_select=intdata_comb1 %>%
        dplyr::select(select_var_list[[j]])
      select_x=as.numeric(c(1,intdata_comb_select[i,]))
      x=c(x,ind_list[[j]][i]*select_x*estweights_indi[[j]][i])
    }
    z=as.numeric(c(1,intdata_Z[i,]))
    g = as.numeric(1/intprob[i] * (intdata_comb1$D[i] - 
                                     exp(theta_hat%*%z)/(1+ exp(theta_hat%*%z)))) * z
    sum1 = sum1 + x %*% t(g)
  }
  sum2 = matrix(0,nrow=x_dim,ncol=z_len)
  for(i in 1:nrow(intdata_comb1)){
    x=NULL
    for(j in 1:K){
      intdata_comb_select=intdata_comb1 %>%
        dplyr::select(select_var_list[[j]])
      select_x=as.numeric(c(1,intdata_comb_select[i,]))
      x=c(x,select_x)
    }
    z=as.numeric(c(1,intdata_Z[i,]))
    g = as.numeric(1/intprob[i] * (intdata_comb1$D[i] - 
                                     exp(theta_hat%*%z)/(1+ exp(theta_hat%*%z)))) * z
    sum2 = sum2 + x %*% t(g)
  }
  E2= (g_alpha) %*%  (solve(H)) %*% (sum1-sum2)
  E2= E2/N_est
  
  ## E3
  E3 = t(E2)
  
  # E4
  
  E4 = matrix(0,nrow=z_len,ncol=z_len)
  sum1 = matrix(0,nrow=x_dim,ncol=x_dim)
  sum2 = matrix(0,nrow=x_dim,ncol=x_dim)
  sum3= matrix(0,nrow=x_dim,ncol=x_dim)
  
  for(i in 1:nrow(intdata_comb1)){
    x=NULL
    for(j in 1:K){
      intdata_comb_select=intdata_comb1 %>%
        dplyr::select(select_var_list[[j]])
      select_x=as.numeric(c(1,intdata_comb_select[i,]))
      x=c(x,ind_list[[j]][i]*select_x*estweights_indi[[j]][i])
    }
    
    sum1 = sum1 + x %*% t(x)
  }
  
  for(i in 1:nrow(intdata_comb1)){
    x=NULL
    for(j in 1:K){
      intdata_comb_select=intdata_comb1 %>%
        dplyr::select(select_var_list[[j]])
      select_x=as.numeric(c(1,intdata_comb_select[i,]))
      x=c(x,select_x)
    }
    
    sum2 = sum2 + x %*% t(x)*comb_est_weights1[i]
  }
  
  for(i in 1:nrow(intdata_comb1)){
    x_a=NULL
    x_b=NULL
    for(j in 1:K){
      intdata_comb_select=intdata_comb1 %>%
        dplyr::select(select_var_list[[j]])
      select_x=as.numeric(c(1,intdata_comb_select[i,]))
      x_b=c(x_b,ind_list[[j]][i]*select_x*estweights_indi[[j]][i])
      x_a=c(x_a,select_x)
    }
    sum3 = sum3 + x_a %*% t(x_b)
  }
  
  E4 =  (g_alpha) %*%  (solve(H)) %*% (sum1+sum2-2*sum3) %*% t((solve(H))) %*% t(g_alpha) 
  
  E4= E4/N_est
  
  E= E1-E2-E3+E4
  
  V = solve(g_theta) %*% E %*% t(solve(g_theta))
  V = V/N_est
  return(list(final_est=wts_fun$final_est,var_est=diag(V)))
}



