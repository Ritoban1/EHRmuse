#' Expit
#'
#' @param x numeric vector
#' @return exp(x)/(1+exp(x))
expit<-function(x){
    return(exp(x)/(1+exp(x)))
}