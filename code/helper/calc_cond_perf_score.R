calc_cond_perf_score <- function(p_ocp_values, p_EN, p_VN, p_delta, p_design_param){
  
  ## conditional performance score values
  if(get_N_fixed(p_delta = p_delta, p_beta = p_design_param$beta) > p_design_param$n_max){
    cp_targ <- p_design_param$alpha
    n_target <- p_design_param$n1
  }else{
    cp_targ <- 1-p_design_param$beta
    n_target <- get_N_fixed(p_delta = delta, p_beta = p_design_param$beta)
  }
  
  ## calculate observed conditional power components
  #ocp_values <- mapply(p_rec_obj$calc_cp, p_z1_values, p_n2_values, p_n3_values, p_effect_est)
  E_ocp <- mean(p_ocp_values)
  Var_ocp <- var(p_ocp_values)
  l_cp <- 1-abs((E_ocp-cp_targ))/(1-p_design_param$alpha) 
  v_cp <- 1-sqrt(Var_ocp/(1/4))
  
  ## calculate sample size components
  l_n <- 1-abs(p_EN-n_target)/(p_design_param$n_max-p_design_param$n1)
  v_n <- 1-sqrt(p_VN/((p_design_param$n_max-p_design_param$n1)/2)^2)
  cond_score_results_ocp <- c(E_ocp = E_ocp, Var_ocp = Var_ocp, E_n = p_EN, Var_n = p_VN, l_cp = l_cp,
                                      v_cp = v_cp, l_n = l_n, v_n = v_n, cond_score = (l_cp+v_cp+l_n+v_n)/4)
  
  ## output score results
  return(cond_score_results_ocp)
  
}