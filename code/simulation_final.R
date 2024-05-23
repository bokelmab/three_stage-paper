## required libraries
library(rpact)
library(dplyr)
library(magrittr)
library(data.table)

## required source files
source('code/helper/power_sample_size_fixed.R')
source('code/helper/create_recalculation_function_final.R')
source('code/helper/calc_cond_perf_score.R')

## simmulation parameters
inf_rates <- c(1/3, 2/3, 1)
beta = 0.2
alpha = 0.025
assumed_delta <- 0.3
interval_delta_eval <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6) ## delta values on which to evaluate
interval_delta_rel <- setdiff(interval_delta_eval, c(0,0.05,0.1,0.15)) ## clinically relevant delta
n_iter <- 10000
sel_seed <- 11234

## define design 
designIN <- getDesignInverseNormal(typeOfDesign = "P", futilityBounds = c(0, 0), informationRates = inf_rates, alpha = alpha, beta = beta)

## calculate sample size for three-stage design based on assumed delta
sample_size_fixed <- getSampleSizeMeans(design = designIN, groups = 2,
                                        alternative = assumed_delta*10, stDev = 10)
sample_size_stages_per_group <- sample_size_fixed%$% c(numberOfSubjects1[1], numberOfSubjects1[2], numberOfSubjects1[3])
sample_size_stages_per_group %<>% round
n1 <- sample_size_stages_per_group[1]
n2_gs <- sample_size_stages_per_group[2]-sample_size_stages_per_group[1]
n3_gs <- sample_size_stages_per_group[3]-sample_size_stages_per_group[2]

## calculate sample size for two-stage design based on assumed delta
design_two_stage <- getDesignInverseNormal(typeOfDesign = "P", futilityBounds = c(0), informationRates = c(1/3,1), alpha = alpha, beta = beta)
sample_size_two_stage_stages_per_group <- c(n1, n1+n2_gs+n3_gs)
n2_two_stage <- n2_gs+n3_gs
design_param_gs2 <- list(f1 = design_two_stage$futilityBounds[1], c1 = design_two_stage$criticalValues[1], c2 = design_two_stage$criticalValues[2],
                     n1 = n1, n2_gs = n2_gs+n3_gs, 
                     w1 = sqrt(design_two_stage$informationRates[1]),
                     w2 = sqrt(design_two_stage$informationRates[2]-design_two_stage$informationRates[1]),
                     gamma = gamma, n_max = 393, alpha = alpha, beta = beta, cp_target = 1-beta)
spec_gs2 <- create_recalculation_function(p_design_param = design_param_gs2, p_delta_values = interval_delta_eval)

gamma <- (get_Pow_fixed(get_N_fixed(interval_delta_rel, p_beta = beta)+0.01, interval_delta_rel)-get_Pow_fixed(get_N_fixed(interval_delta_rel, p_beta = beta),interval_delta_rel))/0.01


## design parameter
design_param <- list(f1 = designIN$futilityBounds[1], c1 = designIN$criticalValues[1], f2 = designIN$futilityBounds[2], c2 = designIN$criticalValues[2], c3 = designIN$criticalValues[3],
                     n1 = sample_size_stages_per_group[1], n2_gs = sample_size_stages_per_group[2]-sample_size_stages_per_group[1], n3_gs = sample_size_stages_per_group[3]-sample_size_stages_per_group[2], 
                     w1 = sqrt(designIN$informationRates[1]),
                     w2 = sqrt(designIN$informationRates[2]-designIN$informationRates[1]), w3 = sqrt(designIN$informationRates[3]-designIN$informationRates[2]),
                     gamma = gamma, n_max = 393, alpha = alpha, beta = beta, cp_target = 1-beta)

## create recalculation function
spec_rec <- create_recalculation_function(p_design_param = design_param, p_delta_values = interval_delta_rel)
rec_func_ocp <- spec_rec$obtain_rec_function(spec_rec$recalc_ocp) 
rec_func_opt <- spec_rec$obtain_rec_function(spec_rec$recalc_opt_n)
design_param$gamma <- rep(0.00295, length(interval_delta_rel))
spec_rec_jt <- create_recalculation_function(p_design_param = design_param, p_delta_values = interval_delta_rel)
rec_func_jt <- spec_rec_jt$obtain_rec_function(spec_rec_jt$recalc_opt_n)

## results to obtain in simulation
Pow_gs <- c()
EN_gs <- c()
Pow_gs2 <- c()
EN_gs2 <- c()
Pow_opt <- c()
EN_opt <- c()
Pow_ocp <- c()
EN_ocp <- c()
Pow_jt <- c()
EN_jt <- c()
cond_score_results_ocp <- NULL
cond_score_results_opt <- NULL
cond_score_results_gs <- NULL
cond_score_results_gs2 <- NULL
cond_score_results_jt <- NULL

for(i_delta in 1:length(interval_delta_eval)){
  
  delta <- interval_delta_eval[i_delta]
  
  ## simulate fixed three stage design
  three_stage <- getSimulationMeans(designIN, alternative = delta*10, stDev = 10,
                                    plannedSubjects = 2*sample_size_stages_per_group,
                                    maxNumberOfIterations = n_iter)
  prob_second_stage <- 1-three_stage$futilityPerStage[1,]-three_stage$rejectPerStage[1,]
  prob_third_stage <- prob_second_stage-three_stage$futilityPerStage[2,]-three_stage$rejectPerStage[2,]
  
  data_sim <- getData(three_stage) %>% as.data.table()
  #stage_three_results[[as.character(delta)]] <- data_sim
  stage_two_gs <- (data_sim$stageNumber == 1 & data_sim$trialStop == F)
  
  ## interim results for recalculation
  z1_values <- data_sim[stageNumber == 1 & trialStop == F,]$testStatistic
  effect_est <- sqrt(2/n1)*z1_values
  
  ## sample size and power three stage
  ocp_values_gs <- mapply(spec_rec$calc_cp, z1_values, n2_gs, n3_gs, effect_est)
  Pow_gs <- c(Pow_gs, three_stage$overallReject)
  EN_gs_cond <- mean(data_sim[stageNumber>1 & trialStop,]$numberOfCumulatedSubjects/2)
  EN_gs <- c(EN_gs, (1-prob_second_stage)*n1 + prob_second_stage*EN_gs_cond)
  VarN_gs_cond <- var(data_sim[stageNumber>1 & trialStop,]$numberOfCumulatedSubjects/2)
  
  ## recalculation OCP
  n2_n3_rec_ocp <- sapply(z1_values, rec_func_ocp)
  ocp_values_ocp <- mapply(spec_rec$calc_cp, z1_values, n2_n3_rec_ocp, n2_n3_rec_ocp, effect_est)
  prob_third_stage_ocp <- mapply(spec_rec$calc_prob_third_stage, z1_values, n2_n3_rec_ocp, n2_n3_rec_ocp, delta)
  EN_ocp_cond <- mean(n1+n2_n3_rec_ocp+prob_third_stage_ocp*n2_n3_rec_ocp)
  EN2_ocp_cond <- mean((1-prob_third_stage_ocp)*(n1+n2_n3_rec_ocp)^2+prob_third_stage_ocp*(n1+2*n2_n3_rec_ocp)^2)
  VarN_ocp_cond <- EN2_ocp_cond-EN_ocp_cond^2
  Pow_ocp <- c(Pow_ocp, three_stage$rejectPerStage[1] + prob_second_stage*mean(mapply(spec_rec$calc_cp, z1_values, n2_n3_rec_ocp, n2_n3_rec_ocp, delta)))
  EN_ocp <- c(EN_ocp, (1-prob_second_stage)*n1 + prob_second_stage*EN_ocp_cond)
  
  ## recalculation score-optimized
  n2_n3_rec_opt <- sapply(z1_values, rec_func_opt)
  ocp_values_opt <- mapply(spec_rec$calc_cp, z1_values, n2_n3_rec_opt, n2_n3_rec_opt, effect_est)
  prob_third_stage_opt <- mapply(spec_rec$calc_prob_third_stage, z1_values, n2_n3_rec_opt, n2_n3_rec_opt, delta)
  EN_opt_cond <- mean(n1+n2_n3_rec_opt+prob_third_stage_opt*n2_n3_rec_opt)
  EN2_opt_cond <- mean((1-prob_third_stage_opt)*(n1+n2_n3_rec_opt)^2+prob_third_stage_opt*(n1+2*n2_n3_rec_opt)^2)
  VarN_opt_cond <- EN2_opt_cond-EN_opt_cond^2
  Pow_opt <- c(Pow_opt, three_stage$rejectPerStage[1] + prob_second_stage*mean(mapply(spec_rec$calc_cp, z1_values, n2_n3_rec_opt, n2_n3_rec_opt, delta)))
  EN_opt <- c(EN_opt, (1-prob_second_stage)*n1 + prob_second_stage*EN_opt_cond)
  
  ## recalculation Jennison & Turnbull
  n2_n3_rec_jt <- sapply(z1_values, rec_func_jt)
  ocp_values_jt <- mapply(spec_rec$calc_cp, z1_values, n2_n3_rec_jt, n2_n3_rec_jt, effect_est)
  prob_third_stage_jt <- mapply(spec_rec$calc_prob_third_stage, z1_values, n2_n3_rec_jt, n2_n3_rec_jt, delta)
  EN_jt_cond <- mean(n1+n2_n3_rec_jt+prob_third_stage_jt*n2_n3_rec_jt)
  EN2_jt_cond <- mean((1-prob_third_stage_jt)*(n1+n2_n3_rec_jt)^2+prob_third_stage_jt*(n1+2*n2_n3_rec_jt)^2)
  VarN_jt_cond <- EN2_jt_cond-EN_jt_cond^2
  Pow_jt <- c(Pow_jt, three_stage$rejectPerStage[1] + prob_second_stage*mean(mapply(spec_rec$calc_cp, z1_values, n2_n3_rec_jt, n2_n3_rec_jt, delta)))
  EN_jt <- c(EN_jt, (1-prob_second_stage)*n1 + prob_second_stage*EN_jt_cond)
  
  ## two-stage group sequential
  two_stage <- getSimulationMeans(design_two_stage, alternative = delta*10, stDev = 10,
                                    plannedSubjects = 2*sample_size_two_stage_stages_per_group,
                                    maxNumberOfIterations = n_iter)
  data_sim_two_stage <- getData(two_stage) %>% as.data.table()
  z1_values_two_stage <- data_sim_two_stage[stageNumber == 1 & trialStop == F,]$testStatistic
  effect_est_two_stage <- sqrt(2/n1)*z1_values_two_stage
  ocp_values_gs2 <- mapply(spec_gs2$calc_crp2, z1_values_two_stage, n2_two_stage, effect_est_two_stage)
  Pow_gs2 <- c(Pow_gs2, two_stage$overallReject)
  EN_gs2_cond <- mean(data_sim_two_stage[stageNumber>1 & trialStop,]$numberOfCumulatedSubjects/2)
  EN_gs2 <- c(EN_gs2, (1-prob_second_stage)*n1 + prob_second_stage*EN_gs2_cond)
  VarN_gs2_cond <- var(data_sim_two_stage[stageNumber>1 & trialStop,]$numberOfCumulatedSubjects/2)
  
  ## conditional performance score
  cond_score_ocp <- calc_cond_perf_score(p_ocp_values = ocp_values_ocp, p_EN = EN_ocp_cond, p_VN = VarN_ocp_cond, p_delta = delta, p_design_param = design_param)
  cond_score_results_ocp %<>% rbind(cond_score_ocp)
  cond_score_opt <- calc_cond_perf_score(p_ocp_values = ocp_values_opt, p_EN = EN_opt_cond, p_VN = VarN_opt_cond, p_delta = delta, p_design_param = design_param)
  cond_score_results_opt %<>% rbind(cond_score_opt)
  cond_score_gs <- calc_cond_perf_score(p_ocp_values = ocp_values_gs, p_EN = EN_gs_cond, p_VN = VarN_gs_cond, p_delta = delta, p_design_param = design_param)
  cond_score_results_gs %<>% rbind(cond_score_gs)
  cond_score_gs2 <- calc_cond_perf_score(p_ocp_values = ocp_values_gs2, p_EN = EN_gs2_cond, p_VN = VarN_gs2_cond, p_delta = delta, p_design_param = design_param_gs2)
  cond_score_results_gs2 %<>% rbind(cond_score_gs2)
  cond_score_jt <- calc_cond_perf_score(p_ocp_values = ocp_values_jt, p_EN = EN_jt_cond, p_VN = VarN_jt_cond, p_delta = delta, p_design_param = design_param)
  cond_score_results_jt %<>% rbind(cond_score_jt)
  
  ## print progress
  print(paste0('Delta ', delta))
}
cond_score_results_gs %<>% as.data.table()
cond_score_results_gs2 %<>% as.data.table()
cond_score_results_ocp %<>% as.data.table()
cond_score_results_opt %<>% as.data.table()
cond_score_results_jt %<>% as.data.table()

## get ingredients of score
Nadap <- get_N_fixed(p_delta = interval_delta_rel, p_beta = beta)

## calculate global score values 
gscore_fix <- 0.8-gamma*Nadap
gscore_gs <- Pow_gs[-(1:4)] - gamma*EN_gs[-(1:4)]
gscore_gs2 <- Pow_gs2[-(1:4)] - gamma*EN_gs2[-(1:4)]
gscore_ocp <- Pow_ocp[-(1:4)] - gamma*EN_ocp[-(1:4)]
gscore_opt <- Pow_opt[-(1:4)] - gamma*EN_opt[-(1:4)]
gscore_jt <- Pow_jt[-(1:4)] - gamma*EN_jt[-(1:4)]

### Plot score values, power and sample size ###################################
layout(matrix(c(1,2,3,4,4,4), ncol=3, byrow=TRUE), heights = c(0.9, 0.1))
par(mar=c(3, 3, 3, 1), mgp=c(1.5, 0.4, 0))

## plot global score
plot(interval_delta_eval, c(NA,NA,NA,NA,gscore_gs), lwd = 3, ylim = c(0,0.7),
     type = 'l', main = 'Global Performance Score', ylab = expression(S^G), xlab = expression(Delta), col = 'blue')
lines(interval_delta_eval, c(NA,NA,NA,NA,gscore_ocp), col = 'blueviolet', lwd = 3, lty = 3)
lines(interval_delta_eval, c(NA,NA,NA,NA,gscore_opt), col = 'orange', lwd = 3, lty = 3)
lines(interval_delta_eval, c(NA,NA,NA,NA,gscore_gs2), col = 'green', lwd = 3, lty = 3)
lines(interval_delta_eval, c(NA,NA,NA,NA,gscore_jt), col = 'cadetblue', lwd = 3, lty = 3)

## plot power
plot(interval_delta_eval, Pow_gs, lwd = 3, ylim = c(0, 1), 
     type = 'l', main = 'Power', ylab = 'Power', xlab = expression(Delta), col = 'blue')
polygon(c(interval_delta_rel, rev(c(interval_delta_rel))),
        c(rep(1-beta, length(interval_delta_rel)), rep(0, length(interval_delta_rel))), col=rgb(1, 0, 0,0.3), border=NA)
lines(interval_delta_eval, Pow_ocp, col = 'blueviolet', lwd = 3, lty = 3)
lines(interval_delta_eval, Pow_opt, col = 'orange', lwd = 3, lty = 3)
lines(interval_delta_eval, Pow_gs2, col = 'green', lwd = 3, lty = 3)
lines(interval_delta_eval, Pow_jt, col = 'cadetblue', lwd = 3, lty = 3)

## plot sample size
plot(interval_delta_eval, EN_gs, lwd = 3, 
     type = 'n', main = 'Sample size', ylab = 'E[N]', xlab = expression(Delta), col = 'black',
     ylim = c(min(Nadap), max(Nadap)))
polygon(c(interval_delta_rel, rev(interval_delta_rel)), c(Nadap,rep(max(Nadap), length(interval_delta_rel))),
        col=rgb(1, 0, 0,0.3), border=NA)
lines(interval_delta_eval, EN_gs, col = 'blue', lwd = 3)
lines(interval_delta_eval, EN_ocp, col = 'blueviolet', lwd = 3, lty = 3)
lines(interval_delta_eval, EN_opt, col = 'orange', lwd = 3, lty = 3)
lines(interval_delta_eval, EN_gs2, col = 'green', lwd = 3, lty = 3)
lines(interval_delta_eval, EN_jt, col = 'cadetblue', lwd = 3, lty = 3)

## legend
par(mar=c(1, 1, 0.5, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c('green','blue', 'blueviolet','orange', 'cadetblue')
legend(x = "bottom",inset = 0,
       legend = c('gs2','gs3', 'ocp', expression(optS^G),expression(optN)), 
       col=plot_colors, lty = c(3,1,3,3,3,3), lwd=3, cex=1, horiz = TRUE, seg.len = 3, text.width = 0.1, bty = 'n')
plot.new()

### save results table
results_power <- data.frame(Pow_gs2,Pow_gs,Pow_ocp,Pow_opt,Pow_jt)
row.names(results_power) <- interval_delta_eval
#saveRDS(results_power, 'results/results_power.RDS')

results_N <- data.frame(EN_gs2,EN_gs,EN_ocp,EN_opt,EN_jt)
row.names(results_N) <- interval_delta_eval
#saveRDS(results_N, 'results/results_N.RDS')

results_Sg <- data.frame(gscore_gs2,gscore_gs,gscore_ocp,gscore_opt,gscore_jt)
row.names(results_Sg) <- interval_delta_rel
#saveRDS(results_Sg, 'results/results_Sg.RDS')

######################## Plot Conditional evaluation ###########################
layout(matrix(c(1,2,3,1,4,5,6,6,6), ncol=3, byrow=TRUE), heights = c(0.45, 0.45,0.1))
par(mar=c(3, 3, 3, 1), mgp=c(1.5, 0.4, 0))

## CP Score
plot(interval_delta_eval, cond_score_results_gs$cond_score, main = 'Cond. Perf. -Score', 
     type = 'l', lwd = 3, col = 'blue', ylab = expression(S^C), xlab = expression(Delta), ylim = c(0,1))
lines(interval_delta_eval, cond_score_results_ocp$cond_score, col = 'blueviolet', lwd = 3, lty = 3)
lines(interval_delta_eval, cond_score_results_opt$cond_score, col = 'orange', lwd = 3, lty = 3)
lines(interval_delta_eval, cond_score_results_gs2$cond_score, col = 'green', lwd = 3, lty = 3)
lines(interval_delta_eval, cond_score_results_jt$cond_score, col = 'cadetblue', lwd = 3, lty = 3)

## l_cp
plot(interval_delta_eval, cond_score_results_gs$l_cp, main = 'Location Cond. Power', 
     type = 'l', lwd = 3, col = 'blue', ylab = expression(l[cp]), xlab = expression(Delta), ylim = c(0,1))
lines(interval_delta_eval, cond_score_results_ocp$l_cp, col = 'blueviolet', lwd = 3, lty = 3)
lines(interval_delta_eval, cond_score_results_opt$l_cp, col = 'orange', lwd = 3, lty = 3)
lines(interval_delta_eval, cond_score_results_gs2$l_cp, col = 'green', lwd = 3, lty = 3)
lines(interval_delta_eval, cond_score_results_jt$l_cp, col = 'cadetblue', lwd = 3, lty = 3)

## l_n
plot(interval_delta_eval, cond_score_results_gs$l_n, main = 'Location N', 
     type = 'l', lwd = 3, col = 'blue', ylab = expression(l[n]), xlab = expression(Delta), ylim = c(0,1))
lines(interval_delta_eval, cond_score_results_ocp$l_n, col = 'blueviolet', lwd = 3, lty = 3)
lines(interval_delta_eval, cond_score_results_opt$l_n, col = 'orange', lwd = 3, lty = 3)
lines(interval_delta_eval, cond_score_results_gs2$l_n, col = 'green', lwd = 3, lty = 3)
lines(interval_delta_eval, cond_score_results_jt$l_n, col = 'cadetblue', lwd = 3, lty = 3)


## v_cp
plot(interval_delta_eval, cond_score_results_gs$v_cp, main = 'Variation Cond. Power', 
     type = 'l', lwd = 3, col = 'blue', ylab = expression(v[cp]), xlab = expression(Delta), ylim = c(0,1))
lines(interval_delta_eval, cond_score_results_ocp$v_cp, col = 'blueviolet', lwd = 3, lty = 3)
lines(interval_delta_eval, cond_score_results_opt$v_cp, col = 'orange', lwd = 3, lty = 3)
lines(interval_delta_eval, cond_score_results_gs2$v_cp, col = 'green', lwd = 3, lty = 3)
lines(interval_delta_eval, cond_score_results_jt$v_cp, col = 'cadetblue', lwd = 3, lty = 3)

## v_n
plot(interval_delta_eval, cond_score_results_gs$v_n, main = 'Variation N', 
     type = 'l', lwd = 3, col = 'blue', ylab = expression(v[n]), xlab = expression(Delta), ylim = c(0,1))
lines(interval_delta_eval, cond_score_results_ocp$v_n, col = 'blueviolet', lwd = 3, lty = 3)
lines(interval_delta_eval, cond_score_results_opt$v_n, col = 'orange', lwd = 3, lty = 3)
lines(interval_delta_eval, cond_score_results_gs2$v_n, col = 'green', lwd = 3, lty = 3)
lines(interval_delta_eval, cond_score_results_jt$v_n, col = 'cadetblue', lwd = 3, lty = 3)

## legend
par(mar=c(1, 1, 0.5, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c('green','blue', 'blueviolet','orange', 'cadetblue')
legend(x = "bottom",inset = 0,
       legend = c('gs2','gs3', 'ocp', expression(optS^G),expression(optN)), 
       col=plot_colors, lty = c(3,1,3,3,3,3), lwd=3, cex=1, horiz = TRUE, seg.len = 3, text.width = 0.1, bty = 'n')
plot.new()

## results to save ####
results_Sc <- data.frame(gs2 = cond_score_results_gs2$cond_score, gs3 = cond_score_results_gs$cond_score,
                         ocp = cond_score_results_ocp$cond_score, opt = cond_score_results_opt$cond_score,
                         jt = cond_score_results_jt$cond_score)
row.names(results_Sc) <- interval_delta_eval
#saveRDS(results_Sc, 'results/results_Sc.RDS')

results_lcp <- data.frame(gs2 = cond_score_results_gs2$l_cp, gs3 = cond_score_results_gs$l_cp,
                         ocp = cond_score_results_ocp$l_cp, opt = cond_score_results_opt$l_cp,
                         jt = cond_score_results_jt$l_cp)
row.names(results_lcp) <- interval_delta_eval
#saveRDS(results_lcp, 'results/results_lcp.RDS')

results_vcp <- data.frame(gs2 = cond_score_results_gs2$v_cp, gs3 = cond_score_results_gs$v_cp,
                          ocp = cond_score_results_ocp$v_cp, opt = cond_score_results_opt$v_cp,
                          jt = cond_score_results_jt$v_cp)
row.names(results_vcp) <- interval_delta_eval
#saveRDS(results_vcp, 'results/results_vcp.RDS')

results_ln <- data.frame(gs2 = cond_score_results_gs2$l_n, gs3 = cond_score_results_gs$l_n,
                         ocp = cond_score_results_ocp$l_n, opt = cond_score_results_opt$l_n,
                         jt = cond_score_results_jt$l_n)
row.names(results_ln) <- interval_delta_eval
#saveRDS(results_ln, 'results/results_ln.RDS')

results_vn <- data.frame(gs2 = cond_score_results_gs2$v_n, gs3 = cond_score_results_gs$v_n,
                          ocp = cond_score_results_ocp$v_n, opt = cond_score_results_opt$v_n,
                          jt = cond_score_results_jt$v_n)
row.names(results_vn) <- interval_delta_eval
#saveRDS(results_vn, 'results/results_vn.RDS')

####### sample size and power plot #############################################
layout(matrix(c(1,2,3,3), ncol=2, byrow=TRUE), heights = c(0.9, 0.1))
par(mar=c(3, 3, 3, 1), mgp=c(1.5, 0.4, 0))

## plot recalculation function n2
z1_values <- design_param$f1 + (0:1000)/1000*(design_param$c1-design_param$f1)
n2_rec_values <- sapply(z1_values, rec_func_opt)
plot(z1_values, n2_rec_values, type = 'l', ylim = c(0, 5*n3_gs), lwd = 3, 
     main = expression(n[2](.)), ylab = expression(n[2]), xlab = expression(z[1]))

## plot recalculation function n3
z2_values <- design_param$f2 + (0:1000)/1000*(design_param$c2-design_param$f2)
n3_rec_values <- sapply(z2_values, rec_func3)
plot(z2_values, n3_rec_values, type = 'l', lwd = 3, main = expression(n[3](.)),
     ylab = expression(n[3]), xlab = expression(z[2]*"*"), ylim = c(0, 5*n3_gs))

## plot observed condition power recalculation function
# data_plot_rec_cp <- data_plot_rec_cp[order(z2),]
# data_plot_rec_cp %$% plot(z2, n3/2, type = 'l', lwd = 3, main = 'OCP stage 3',
#                           ylab = expression(n[3]), xlab = expression(z[2]*"*"), ylim = c(0, 5*n3_gs))


