#' Expit
#'
#' @param x numeric vector
#' @return exp(x)/(1+exp(x))
#' @examples expit(1)
expit<-function(x){
    return(exp(x)/(1+exp(x)))
}