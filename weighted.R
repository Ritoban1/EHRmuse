weighted<-function(intdata,estweights,Z_names){
  formula <- as.formula(paste("D ~", paste(Z_names, collapse = " + ")))
  modelinit=glm(formula, family = stats::quasibinomial(),data=intdata)
  start <- coef(modelinit)
  design <- survey::svydesign(data = intdata,ids = 1:length(intdata$D), strata = NULL,
                              weights = estweights)
  mod=survey::svyglm(data=intdata,formula=formula, design = design,
                     family = quasibinomial(),
                     start = as.numeric(start))
  final<-coef(mod)
  var=diag(summary(mod)$cov.scaled)
  return(list(final=final,var=var))
}
