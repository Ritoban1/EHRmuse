unweighted<-function(K,intdata_list,Z_names){
  intdata_comb=NULL
  for(i in 1:K){
    intdata=intdata_list[[i]]
    if (!is.data.frame(intdata)) stop("All Internal data should be a data frame.")
    intdata_comb=rbind.data.frame(intdata_comb,intdata)
  }
  intdata_comb1=intdata_comb[!duplicated(intdata_comb$id),]
  formula <- as.formula(paste("D ~", paste(Z_names, collapse = " + ")))
  modelinit=glm(formula, family = stats::quasibinomial(),data=intdata_comb1)
  
  start <- coef(modelinit)
  mod<-summary(modelinit)
  var<-(mod$coefficients[,2])^2
  return(list(est=start,var=var))
}
