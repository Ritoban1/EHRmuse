
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EHRmuse

## Installation

``` r
# From CRAN:
install.packages("EHRmuse")
```

### Development version

``` r
# install.packages("pak")
pak::pak("Ritoban1/EHRmuse")
```

This R package, EHRmuse implements Joint Inverse Probability Weighted
(JIPW) and Joint Augmented Inverse Probability Weighted (JAIPW) methods
to address selection bias in non-probability samples, such as Electronic
Health Records (EHRs). These methods leverage data integration
techniques, incorporating either individual-level data or summary-level
statistics from external sources, to improve the estimation of
association parameters in binary disease risk models. The link to the
methods paper is: <https://arxiv.org/abs/2412.00228>.

Selection bias poses significant challenges when working with EHR data,
particularly when participants are recruited from multiple clinics or
centers, each with distinct and potentially outcome-dependent selection
mechanisms. Standard inverse-probability-weighted (IPW) methods, which
rely on parametric selection models, often suffer from inaccuracies due
to model misspecification in such settings.

The JAIPW method enhances robustness by integrating individual-level
data from multiple cohorts with diverse selection mechanisms and
supplementing it with data from an external probability sample. By
incorporating a flexible auxiliary score model, JAIPW achieves double
robustness, effectively mitigating biases caused by misspecified
selection models and improving the reliability of analysis in complex,
non-probability sampling scenarios.

Let $D$ be a binary disease indicator, and let $\mathbf{Z}$ denote a set
of covariates. The primary disease model of interest in the target
population is:

$$
\text{logit}(P(D=1|\mathbf{Z})) = \theta_0 + \mathbf{\theta}_{\mathbf{Z}}' \mathbf{Z}.
$$

We analyze data from $K$ internal non-probability samples (cohorts)
drawn from the same target population. Each cohort $k$
($k = 1, 2, \ldots, K$) is characterized by a binary selection indicator
$S_k$, which represents inclusion in the cohort. For each cohort $k$,
$Z_{1k}$ denotes the subset of covariates that appear exclusively in the
disease model and do not directly influence the corresponding selection
indicator $S_k$. $Z_{2k}$ includes covariates shared by both the disease
and selection models. Across all cohorts, the set of covariates
$\mathbf{Z}$ is the union of $Z_{1k}$ and $Z_{2k}$. Additionally, $W_k$
represents variables specific to the selection model for cohort $k$,
which can vary across cohorts. The probability of selection into cohort
$k$, given covariates, is modeled as $P(S_k=1|X_k) = \pi_k(X_k)$, where
$X_k$ includes $D$, $Z_{2k}$, and $W_k$.

## Data Generation

The data is generated using the R script **“Data_gen_indi.R”**. For this
example, we set the number of cohorts to $K = 3$, the population size to
$N = 50{,}000$, and the dimension of the covariate vector to
$\text{dim}(Z) = 3$. Different selection mechanisms are employed for the
three cohorts to reflect cohort-specific variability. Additionally, we
simulated individual-level external data to evaluate the methods JPL
(Joint Pseudolikelihood), JSR (Joint Simplex Regression), and JAIPW
(Joint Augmented Inverse Probability Weighted), and generated marginal
totals for JCL.

``` r
expit<-function(x){
  return(exp(x)/(1+exp(x)))
}

K=3 ## Number of Cohorts
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

theta=c(-2,0.35,0.45,0.25) ## Theta_Z vector
N=5e4 ## Population size

### selection models
dw=1
dwz1=c(1,0.8,0.6)
dwz2=c(0.6,0.8,1)
dwz3=rep(1,3)

gamma_ext=c(-0.6,1.2,0.4,-0.2,0.5)
gamma_int_1=c(-1,1.5,0.2,0.8,-0.3)
gamma_int_2=c(-1,1.25,0.4,0.6)
gamma_int_3=c(-3,0.8,0.5)

## Generation of population level data
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
## Generation of external individual level data
simu_ext<-function(data,gamma_ext){
  U2e <- runif(N)
  # Generate Sampling Status
  SELECT <-0.75*expit(gamma_ext[1] +  
                        gamma_ext[2]* data$D + 
                        gamma_ext[3] * data$Z1 +
                        gamma_ext[4]* data$Z2 +
                        gamma_ext[5] * data$Z3)
  S_e  <- ifelse(SELECT > U2e, TRUE, FALSE)
  # Observed Data
  data_e <- data[which(S_e==1),]
  data_e$Select_Weights = 0.75*expit(gamma_ext[1] +  
                                       gamma_ext[2]* data_e$D + 
                                       gamma_ext[3] * data_e$Z1 +
                                       gamma_ext[4]* data_e$Z2 +
                                       gamma_ext[5] * data_e$Z3)
  return(data_e)
}
## Generation of internal data 1
simu_int_1<-function(data,gamma_int_1){
  U2i <- runif(N)
  # Generate Sampling Status
  SELECT <- expit(cbind(1,data$D,data$W_1,data$Z2,data$Z3)
                  %*% gamma_int_1)
  S_i  <- ifelse(SELECT > U2i, TRUE, FALSE)
  # Observed Data
  data_i <- data[which(S_i==1),]
  return(data_i)
}

## Generation of internal data 2
simu_int_2<-function(data,gamma_int_2){
  U2i <- runif(N)
  # Generate Sampling Status
  SELECT <- expit(cbind(1,data$D,data$W_2,data$Z3)
                  %*% gamma_int_2)
  S_i  <- ifelse(SELECT > U2i, TRUE, FALSE)
  # Observed Data
  data_i <- data[which(S_i==1),]
  return(data_i)
}


## Generation of internal data 3
simu_int_3<-function(data,gamma_int_3){
  U2i <- runif(N)
  # Generate Sampling Status
  SELECT <- expit(cbind(1,data$W_3,data$Z2)
                  %*% gamma_int_3)
  S_i  <- ifelse(SELECT > U2i, TRUE, FALSE)
  # Observed Data
  data_i <- data[which(S_i==1),]
  return(data_i)
}

data=simu_popu(N,mean_w_p,mean_z_1,mean_z_2,mean_z_3,
               var_z_w_p,theta,dw)

extdata=simu_ext(data,gamma_ext)


intdata1=simu_int_1(data,gamma_int_1)
intdata2=simu_int_2(data,gamma_int_2)
intdata3=simu_int_3(data,gamma_int_3)

## names of selection variables in each cohort
select_var_list=list(c("D","W_1","Z2","Z3"),c("D","W_2","Z3"),c("W_3","Z2"))

## names of auxiliary variables in each cohort
aux_var_list=list(c("D","W_1","Z2","Z3"),c("D","W_2","Z3"),c("W_3","Z2"))

## list of internal data
intdata_list=list(intdata1,intdata2,intdata3)
## names of Z variables
Z_names=c("Z1","Z2","Z3")

theta ## actual theta_z
```

## Unweighted Logistic Method

At first, we implement an unweighted logistic method.

``` r
## Unweighted Logistic Method with variance estimate

res_uw=EHRmuse(K=K,N=N,Z_names=Z_names,
                 intdata_list=intdata_list,variance = TRUE)
res_uw

## set UW_CS=TRUE for cohort specific intercepts
```

## Inverse Probability Weighted Regression Using Individual level External data

## Pseudolikelihood (PL)

``` r
## Approximate variance

res_pl= EHRmuse(K=K,N=N,Z_names=Z_names,intdata_list=intdata_list,IPW=TRUE,
              select_var_list=select_var_list,
              extdata=extdata,Weights_e = 1/extdata$Select_Weights,
              variance = TRUE,ipw_method ="PL")
res_pl

## Asymptotic variance

# res_pl= EHRmuse(K=K,N=N,Z_names=Z_names,intdata_list=intdata_list,IPW=TRUE,
#               select_var_list=select_var_list,
#               extdata=extdata,Weights_e = 1/extdata$Select_Weights,
#               variance = TRUE,ipw_method ="PL", type_var = "asy")

## Just change to ipw_method="SR" for Simplex Regression
```

## Inverse Probability Weighted Regression Using Marginal level external data

## Calibration (CL)

The data is generated using the R script **“Data_gen_marg.R”**

``` r
expit<-function(x){
  return(exp(x)/(1+exp(x)))
}
K=3 ## Number of cohorts
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

theta=c(-2,0.35,0.45,0.25) ## Theta_Z vector
N=5e4 ## Population size
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
  S_i  <- ifelse(SELECT > U2i, TRUE, FALSE)
  # Observed Data
  data_i <- data[which(S_i==1),]
  return(data_i)
}

simu_int_2<-function(data,gamma_int_2){
  U2i <- runif(N)
  # Generate Sampling Status
  SELECT <- expit(cbind(1,data$D,data$W_2,data$Z3)
                  %*% gamma_int_2)
  S_i  <- ifelse(SELECT > U2i, TRUE, FALSE)
  # Observed Data
  data_i <- data[which(S_i==1),]
  return(data_i)
}



simu_int_3<-function(data,gamma_int_3){
  U2i <- runif(N)
  # Generate Sampling Status
  SELECT <- expit(cbind(1,data$W_3,data$Z2)
                  %*% gamma_int_3)
  S_i  <- ifelse(SELECT > U2i, TRUE, FALSE)
  # Observed Data
  data_i <- data[which(S_i==1),]
  return(data_i)
}
data=simu_popu(N,mean_w_p,mean_z_1,mean_z_2,mean_z_3,
               var_z_w_p,theta,dw)
intdata1=simu_int_1(data,gamma_int_1)
intdata2=simu_int_2(data,gamma_int_2)
intdata3=simu_int_3(data,gamma_int_3)

## 
marginal=totals(data) ## generation of marginal totals

##list of marginals
margs=list(marginal[c(1,2,3,6,7)],
           marginal[c(1,2,4,7)],
           marginal[c(1,5,6)])

## names of selection variables in each cohort
select_var_list=list(c("D","W_1","Z2","Z3"),c("D","W_2","Z3"),c("W_3","Z2"))

## names of auxiliary variables in each cohort
aux_var_list=list(c("D","W_1","Z2","Z3"),c("D","W_2","Z3"),c("W_3","Z2"))

## list of internal data
intdata_list=list(intdata1,intdata2,intdata3)

## names of Z variables
Z_names=c("Z1","Z2","Z3")
```

``` r

res_cl=EHRmuse(K=K,N=N,Z_names=Z_names,intdata_list=intdata_list,IPW=TRUE,
            select_var_list=select_var_list,
            marginals_list = margs,
            variance =TRUE,ipw_method = "CL")
```

## Joint Augmented Inverse Probability Weighted Method (JAIPW)

``` r

res_aipw=EHRmuse(K=K,N=N,Z_names=Z_names,intdata_list=intdata_list,IPW=TRUE,
              select_var_list=select_var_list,aux_var_list=aux_var_list,
              extdata=extdata,Weights_e = 1/extdata$Select_Weights,
              variance = TRUE,ipw_method ="PL",AIPW=TRUE)
```

## Getting help

Contact the author at <kundur@umich.edu>.
