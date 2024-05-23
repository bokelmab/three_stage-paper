create_recalculation_function <- function(p_design_param, p_delta_values){
  
  calc_env <- environment()
  ## take relevant parameters from design
  c1 <- p_design_param$c1
  c2 <- p_design_param$c2
  c3 <- p_design_param$c3
  f1 <- p_design_param$f1
  f2 <- p_design_param$f2
  w1 <- p_design_param$w1
  w2 <- p_design_param$w2
  w3 <- p_design_param$w3
  n1 <- p_design_param$n1
  gamma <- p_design_param$gamma
  delta_values <- p_delta_values
  n_approx <- 100
  n_max <- p_design_param$n_max
  cp_target <- p_design_param$cp_target
  
  output <- list(
    
    env = calc_env,
    
    
    calc_crp2 = function(p_z, p_n, p_delta){
      return(1-pnorm((sqrt(w1^2+w2^2)*c2-w1*p_z)/w2-p_delta*sqrt(p_n/2)))
    },
    calc_crp3 = function(p_z, p_n, p_delta){
      return(1-pnorm((c3-sqrt(w1^2+w2^2)*p_z)/w3-p_delta*sqrt(p_n/2)))
    },
    calc_cp = function(p_z1, p_n2, p_n3, p_delta){
      
      crp2 <- output$calc_crp2(p_z1, p_n2, p_delta) ## probability to reject at stage 2
      
      func_crp3 <- function(p_z2){
        return(output$calc_crp3(p_z = p_z2, p_n = p_n3, p_delta = p_delta)*output$cond_dens_z2(p_z2 = p_z2, p_n2 = p_n2, p_delta = p_delta, p_z1 = p_z1))
      }
      crp3 <- integrate(func_crp3, lower = f2, upper = c2)$value ## probability to reject at stage 3
      
      return(crp2+crp3) ## return conditional power
    },
    recalc_ocp = function(p_z1){
      
      delta_est <- sqrt(2/n1)*p_z1
      n_values <- 0:floor((n_max-n1)/2) ## try these sample sizes for stage 2 and 3
      
      ## calculate for each possible sample size the ocp value
      ocp_values <- mapply(output$calc_cp, p_z1, n_values, n_values, delta_est)
      
      ## choose sample size, such that ocp meets target
      idx_cp <- which(ocp_values >= cp_target)
      if(length(idx_cp) == 0){
        return(max(n_values))
      }else{
        return(n_values[min(idx_cp)])
      }
      
    },
    calc_En = function(p_z1, p_n2, p_n3, p_delta){
      func_dens_z2 <- function(p_z2){
        return(output$cond_dens_z2(p_z2 = p_z2, p_n2 = p_n2, p_delta = p_delta, p_z1 = p_z1))
      }
      return(p_n2+p_n3*integrate(func_dens_z2, lower = f2, upper = c2)$value)
    },
    # calc_condVar_n = function(p_z1, p_n2, p_n3, p_delta){ ## variance of n, given z1
    #   
    #   ## first case: trial stops
    #   if(p_z1 <= f1 || p_z1 >= c1){
    #     return(0)
    #   }
    #   
    #   ## second case: trial continues in second stage
    #   ## E[N^2]-E[N]^2
    #   func_dens_z2 <- function(p_z2){
    #     return(output$cond_dens_z2(p_z2 = p_z2, p_n2 = p_n2, p_delta = p_delta, p_z1 = p_z1))
    #   }
    #   prob_third_stage <- integrate(func_dens_z2, lower = f2, upper = c2)$value
    #   E_N2 <- (1-prob_third_stage)*(n1+p_n2)^2 + prob_third_stage*(n1+p_n2+p_n3)^2
    #   E_N <- (1-prob_third_stage)*(n1+p_n2) + prob_third_stage*(n1+p_n2+p_n3)
    #   return(E_N2-E_N^2)
    # },
    calc_prob_third_stage = function(p_z1, p_n2, p_n3, p_delta){
        func_dens_z2 <- function(p_z2){
          return(output$cond_dens_z2(p_z2 = p_z2, p_n2 = p_n2, p_delta = p_delta, p_z1 = p_z1))
        }
        prob_third_stage <- integrate(func_dens_z2, lower = f2, upper = c2)$value
        return(prob_third_stage)
    },
    recalc_opt_n = function(p_z1){
      
      n_values <- 0:floor((n_max-n1)/2)
      score_per_n <- NULL
      
      ## optimal 2nd stage sample size at z1
      for(i_delta in 1:length(delta_values)){
        delta <- delta_values[i_delta]
        cp_per_n <- mapply(output$calc_cp, p_z1, n_values, n_values, delta)
        En_per_n <- mapply(output$calc_En, p_z1, n_values, n_values, delta)
        post_dens_delta <- output$cond_dens_z1(p_z1 = p_z1, p_n1 = n1, p_delta = delta)
        score_per_n %<>% cbind((cp_per_n - gamma[i_delta]*En_per_n)*post_dens_delta)
      }
      
      ## obtain optimal n
      idx_max <- which.max(apply(score_per_n,1,mean))
      return(n_values[idx_max])
    },
    norm_dens <- function(p_x, p_mu, p_sigma){
      
      return(1/sqrt(2*pi)/p_sigma*exp(-1/2*((p_x-p_mu)/p_sigma)^2))
    },
    cond_dens_z2 = function(p_z2, p_n2, p_delta, p_z1){
      
      mu <- (w1*p_z1+w2*sqrt(p_n2/2)*p_delta)/sqrt(w1^2+w2^2)
      sigma <- (w2)/sqrt(w1^2+w2^2)
      return(norm_dens(p_x = p_z2, p_mu = mu, p_sigma = sigma))
    },
    cond_dens_z1 = function(p_z1, p_n1, p_delta){
      return(1/sqrt(2*pi)*exp(-((p_z1-p_delta*sqrt(p_n1/2))^2)/2))
    },
    obtain_rec_function = function(p_algorithm){
      
      ## get a grid of z values in recalculation area
      z_values <- f2+(0:n_approx)/n_approx*(c2-f2)
      
      ## get optimal n for a grid of z values
      data_recalc <- data.frame(z_value = z_values, n = sapply(z_values, p_algorithm))
      
      ## build recalculation function
      rec_n <- function(p_z){
        
        ## pick nearest neighbor to z value
        b <- min(which(data_recalc$z_value >= p_z))
        if(b == 1){
          return(data_recalc$n[b])
        }
        a <- b-1
        return(data_recalc$n[a])
        # if(data_recalc$n[a] == data_recalc$n[b]){
        #   return(data_recalc$n[a])
        # }else{
        #   ## calculate optimal n for p_Zint
        #   return(output$calculate_opt_n(p_z))
        # }
      }
      
      ## return recalculation function
      return(rec_n)
      
    },
    fut_sec_stage = function(p_z1, p_delta, p_n2){
      
      func <- function(p_z2){
        z1 <- p_z1
        delta <- p_delta
        n2 <- p_n2
        return(output$cond_dens_z2(p_z2, p_n2 = n2, p_delta = delta, p_z1 = z1))
      }
      
      result <- integrate(func, lower = -Inf, upper = f2)
      return(result$value)
    },
    pow_third_stage = function(p_z1, p_delta, p_n2, p_n3){
      
      func_z2 <- function(p_z2){
        z1 <- p_z1
        delta <- p_delta
        n2 <- p_n2
        return(output$cond_dens_z2(p_z2 = p_z2, p_n2 = n2, p_delta = delta, p_z1 = z1))
      }
      
      func_cp3 = function(p_z2){
        z1 <- p_z1
        n3 <- p_n3
        delta <- p_delta
        return(1-pnorm((c3-w1*z1-w2*p_z2-w3*sqrt(n3/2)*delta)/w3))
      }
      
      func <- function(p_z2){
        return(func_z2(p_z2)*func_cp3(p_z2))
      }
      result <- integrate(func, lower = f2, upper = c2)
      return(result$value)
      
      
    }
    
  )
  
}


