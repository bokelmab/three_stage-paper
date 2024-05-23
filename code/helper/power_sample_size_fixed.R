get_N_fixed = function(p_delta, p_beta = 0.1){
  
  return(2*((qnorm(p=0.975)-qnorm(p_beta))/p_delta)^2)
}
get_Pow_fixed = function(p_N, p_delta){
  return(1-pnorm(qnorm(p=0.975)-p_delta*sqrt(p_N/2)))
}



