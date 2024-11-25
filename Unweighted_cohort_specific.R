unweighted_cs<-function(K,intdata_list,Z_names){
  intdata_comb=NULL
  cohort=NULL
  for(i in 1:K){
    intdata=intdata_list[[i]]
    if (!is.data.frame(intdata)) stop("All Internal data should be a data frame.")
    intdata_comb=rbind.data.frame(intdata_comb,intdata)
    cohort=c(cohort,rep(i,nrow(intdata)))
  }
  intdata_comb$cohort=as.factor(cohort)
  intdata_comb1=intdata_comb[!duplicated(intdata_comb$id),]
  
  Z_names=c(Z_names,"cohort")
  formula <- as.formula(paste("D ~", paste(Z_names, collapse = " + ")))
  modelinit=glm(formula, family = stats::quasibinomial(),data=intdata_comb1)
  
  start <- coef(modelinit)[c(1:length(Z_names))]
  mod<-summary(modelinit)
  var<-(mod$coefficients[c(1:length(Z_names)),2])^2
  return(list(est=start,var=var))
}
