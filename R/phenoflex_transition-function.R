#just for plotting, though
calc_py <- function(yc, s1, y){
  sy <- exp(s1 * yc * ((y-yc)/y))
  return((sy)/(sy+1))
}

get_transition_fun <- function(yc, s1, y_vector = NULL, p = 0.5){
  if(is.null(y_vector)){
    y_vector <- seq(yc * (1-p), yc *(1+p), 1)
  }
  
  purrr::map_dbl(y_vector, function(y) calc_py(yc, s1, y))
  
}