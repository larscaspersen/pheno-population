#sampling population
get_skewed_dist <- function(mean, sd, skew, n=100){
  
  delta <- skew / sqrt(1 + skew^2)
  
  omega <- sd /
    sqrt(1 - 2 * delta^2 / pi)
  
  xi <- mean -
    omega * delta * sqrt(2 / pi)
  
  set.seed(123)
  
  return(sn::rsn(
    n = n,
    xi = xi,
    omega = omega,
    alpha = skew
  ))
}

helper_sample_dist <- function(n, par, yc_sd, zc_sd, dist_chill, dist_heat, add_par, seed = 12345){
  if(n > 1){
    set.seed(seed)
    if(yc_sd > 0){
      if(dist_chill == 'normal'){
        yc_pop <- rnorm(n = n, mean = par[1], sd =  yc_sd)
      } else if(dist_chill == 'normal_skewed'){
        yc_pop <- get_skewed_dist(mean = par[1], sd = yc_sd, skew = add_par[1],  n = n)
      }
    }
    if(zc_sd > 0){
      if(dist_heat == 'normal'){
        zc_pop <- rnorm(n = n, mean = par[2], sd = zc_sd)
      } else if(dist_heat == 'normal_skewed'){
        zc_pop <- get_skewed_dist(mean = par[2], sd = zc_sd, skew = add_par[2],  n = n)
      }
    }
  }
  
  return(list('yc_pop' = yc_pop |> as.vector(),
              'zc_pop' = zc_pop |> as.vector()))
}
