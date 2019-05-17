#' @export
approx_zono <- function(P, fit_ratio = NULL, Parameters = NULL){
  
  ret_list = zono_approx(P, fit_ratio, Parameters)
  
  Mat = ret_list$Mat
  
  # first column is the vector b
  b = Mat[,1]
  
  # remove first column
  A = Mat[,-c(1)]
  PP = list("P" = Hpolytope$new(A,b), "fit_ratio" = ret_list$fit_ratio)
  
  return(PP)
  
}