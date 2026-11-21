################################################################################
####       Confounding Heterogeneity in Matched Observational Studies       ####
####                 Reproducible Code for Main Calculations                ####
####                Periodontal Disease & Smoking Application               ####
################################################################################

#### README --------------------------------------------------------------------

# running this file produces the following outputs:

# 1) "sensitivity_values_main.csv": Reproduces table 3 in the main text
# (Sensitivity Values)

# 2) "sensitivity_values_supplement.csv": Reproduces table S6 in the Supplementary
# Materials.

# all files are placed in the folder "Confounding_Heterogeneity_Application".

#### Libraries and functions ---------------------------------------------------

# installs any missing packages, then loads them into the library
required_packages <- c("quadprog", "ggplot2", "latex2exp", "lqmm", "mosaic", 
                       "stats", "tidyverse", "mvtnorm", "MASS", "VGAM", 
                       "Rsolnp", "nloptr", "sensitivitymult", "doSNOW", 
                       "foreach", "doParallel")
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages) > 0) {
  install.packages(missing_packages)
}
invisible(lapply(required_packages, library, character.only = TRUE))

# creates folder holding all simulation output
dir.create("Confounding_Heterogeneity_Application")
setwd("Confounding_Heterogeneity_Application")

# Confounder acts first
chibar_noconfounding_classic <- function(data, kappa, alpha = 0.05, qfunc = "perm"){
  
  if(qfunc == "perm"){
    qf <- function(Y){
      return(abs(Y))
    }
  }
  I <- nrow(data)/2
  Z <- data[, 1]
  outcomes <- data[, 2:ncol(data)]
  K <- ncol(outcomes)
  Yc <- outcomes[seq(1, nrow(data), by = 2), ] - outcomes[seq(2, nrow(data), by = 2), ] # treated - control
  Yckappa <- base::t(base::t(Yc) * kappa)
  
  # calculate the test statistic
  observed_test_stat <- observed_test_statistic(data, qfunc)
  observed_cov <- Covariance_Rosenbaum(data, rep(0, I), qfunc)
  scaled_test_stat <- abs(kappa %*% base::t(as.matrix(observed_test_stat))) / sqrt((kappa %*% observed_cov %*% kappa))
  
  # calculate the critical value
  cor <- Yckappa[, 1] %*% Yckappa[, 2] / sqrt(sum(Yckappa[, 1]^2 * sum(Yckappa[, 2]^2)))
  M <- matrix(c(1, cor, cor, 1), nrow = 2, byrow = TRUE)
  Minv <- solve(M)
  rho <- Minv[1,2] / sqrt(Minv[1,1]*Minv[2,2])
  w2 <- 1/4 + asin(rho)/(2*pi)
  w0 <- 1/4 - asin(rho)/(2*pi)
  w1 <- 1 - w0 - w2
  
  Fmix <- function(x){
    w0 +
      w1 * pchisq(x,1) +
      w2 * pchisq(x,2)
  }
  
  crit_val <- uniroot(function(x) Fmix(x)-(1-alpha), interval=c(0,100))$root
  pval <- 1 - Fmix(scaled_test_stat^2)
  
  output <- list("Test_stat" = scaled_test_stat,
                 "p_value" = pval, 
                 "critical_val" = crit_val)
  return(output)
  
  
}
chibar_noconfounding_amp <- function(data, kappa, alpha = 0.05, qfunc = "perm"){
  
  if(qfunc == "perm"){
    qf <- function(Y){
      return(abs(Y))
    }
  }
  I <- nrow(data)/2
  Z <- data[, 1]
  outcomes <- data[, 2:ncol(data)]
  K <- ncol(outcomes)
  Yc <- outcomes[seq(1, nrow(data), by = 2), ] - outcomes[seq(2, nrow(data), by = 2), ] # treated - control
  Yckappa <- base::t(base::t(Yc) * kappa)
  
  # calculate the test statistic
  observed_test_stat <- observed_test_statistic(data, qfunc)
  observed_cov <- Covariance_CH2(data, rep(1, I), 1, c(1, 1), qfunc)
  scaled_test_stat <- abs(kappa %*% base::t(as.matrix(observed_test_stat))) / sqrt((kappa %*% observed_cov %*% kappa))
  
  # calculate the critical value
  cor <- Yckappa[, 1] %*% Yckappa[, 2] / sqrt(sum(Yckappa[, 1]^2 * sum(Yckappa[, 2]^2)))
  M <- matrix(c(1, cor, cor, 1), nrow = 2, byrow = TRUE)
  Minv <- solve(M)
  rho <- Minv[1,2] / sqrt(Minv[1,1]*Minv[2,2])
  w2 <- 1/4 + asin(rho)/(2*pi)
  w0 <- 1/4 - asin(rho)/(2*pi)
  w1 <- 1 - w0 - w2
  
  Fmix <- function(x){
    w0 +
      w1 * pchisq(x,1) +
      w2 * pchisq(x,2)
  }
  
  crit_val <- uniroot(function(x) Fmix(x)-(1-alpha), interval=c(0,100))$root
  pval <- 1 - Fmix(scaled_test_stat^2)
  
  output <- list("Test_stat" = scaled_test_stat,
                 "p_value" = pval, 
                 "critical_val" = crit_val)
  return(output)
  
  
}
chibar_classic_ttest <- function(data, Gamma, kappa, alpha = 0.05, qfunc = "perm"){
  
  if(qfunc == "perm"){
    qf <- function(Y){
      return(abs(Y))
    }
  }
  
  I <- nrow(data)/2
  Z <- data[, 1]
  outcomes <- data[, 2:ncol(data)]
  K <- ncol(outcomes)
  test_stat <- observed_test_statistic(data, qfunc = "perm")
  Yc <- outcomes[seq(1, nrow(data), by = 2), ] - outcomes[seq(2, nrow(data), by = 2), ] # treated - control
  Yckappa <- base::t(base::t(Yc) * kappa)
  
  # calculating worst-case u
  
  Dmat <- matrix(0, nrow = I, ncol = I)
  for(i in 1:(I-1)){
    for(j in (i+1):I){
      Dmat[i, j] <- sum(base::t((qf(as.matrix(Yckappa[i, ])) * sign(as.matrix(Yckappa[i, ])))) %*% (qf(as.matrix(Yckappa[j, ])) * sign(as.matrix(Yckappa[j, ]))))
      Dmat[j, i] <- Dmat[i, j]
    }
  }
  
  for(i in 1:I){
    Dmat[i, i] <- sum(base::t((qf(as.matrix(Yckappa[i, ])) * sign(as.matrix(Yckappa[i, ])))) %*% (qf(as.matrix(Yckappa[i, ])) * sign(as.matrix(Yckappa[i, ]))))
  }
  
  Dmat_pd <- make.positive.definite(Dmat)
  dvec <- as.numeric(rowSums(qf(Yckappa)*sign(Yckappa)) * sum(test_stat * kappa))
  Amat <- cbind(diag(I), -diag(I))
  bvec <- rep(-(Gamma-1)/(Gamma+1), 2*I)
  rough_solution_u <- solve.QP(Dmat_pd, dvec, Amat, bvec, meq = 0)$solution
  
  worst_ev <- Expected_value_Rosenbaum(data, rough_solution_u, qfunc = "perm")
  worst_cov <- Covariance_Rosenbaum(data, rough_solution_u, qfunc = "perm")
  worst_test_stat <- abs(kappa %*% base::t(as.matrix(test_stat - worst_ev))) / sqrt((kappa %*% worst_cov %*% kappa))
  
  # kappa_func <- function(kappa){
  #   abs(kappa %*% base::t(as.matrix(test_stat - worst_ev))) / sqrt((kappa %*% worst_cov %*% kappa))
  # }
  # find_worst_test <- optim(rep(1, K), kappa_func, method = "L-BFGS-B", lower = c(0, 0), control = list(fnscale = -1))
  # worst_test_stat <- find_worst_test$value
  
  # finding the lower bound on the correlation between outcomes
  
  correlation_objective <- function(p){
    return((Yckappa[, 1] %*% Yckappa[, 2] - (Yckappa[, 1] * p) %*% (Yckappa[, 2] * p))/(sqrt(sum(Yckappa[, 1]^2 * (1-p^2))) * sqrt(sum(Yckappa[, 2]^2 * (1-p^2)))))
  }
  
  lb_cor <- optim(rep((Gamma-1)/(1+Gamma), I), correlation_objective, method = "L-BFGS-B", lower = (1-Gamma)/(1+Gamma), upper = (Gamma-1)/(1+Gamma))$value
  
  # calculating distribution function with worse-case crit val
  
  M <- matrix(c(1, lb_cor, lb_cor, 1), nrow = 2, byrow = TRUE)
  Minv <- solve(M)
  rho <- Minv[1,2] / sqrt(Minv[1,1]*Minv[2,2])
  w2 <- 1/4 + asin(rho)/(2*pi)
  w0 <- 1/4 - asin(rho)/(2*pi)
  w1 <- 1 - w0 - w2
  
  Fmix <- function(x){
    w0 +
      w1 * pchisq(x,1) +
      w2 * pchisq(x,2)
  }
  
  crit_val <- uniroot(function(x) Fmix(x)-(1-alpha), interval=c(0,100))$root
  
  pval <- 1 - Fmix(worst_test_stat^2)
  
  output <- list("sensitivity_analysis_stat" = worst_test_stat, 
                 "worst_case_pvalue" = pval,
                 "worst_case_critval" = crit_val)
  return(output)
}
chibar_T2CH_ttest <- function(data, Lambda, Delta_vec, kappa, alpha = 0.05, qfunc = "perm"){
  
  if(qfunc == "perm"){
    qf <- function(Y){
      return(abs(Y))
    }
  }
  
  lambda <- log(Lambda)
  delta_vec <- log(Delta_vec)
  I <- nrow(data)/2
  Z <- data[, 1]
  outcomes <- data[, 2:ncol(data)]
  K <- ncol(outcomes)
  test_stat <- observed_test_statistic(data, qfunc)
  Yc <- outcomes[seq(1, nrow(data), by = 2), ] - outcomes[seq(2, nrow(data), by = 2), ]
  qYc <- qf(base::t(base::t(Yc) * kappa))
  
  program <- function(u){
    # u is an I-length vector of exp\{u_{i1}-u_{i2}\}, i=1,...,I  
    # u <- runif(I, 1/exp(1), exp(1))
    
    worst_ev <- Expected_value_CH2(data, u, Lambda, Delta_vec, qfunc)
    worst_cov <- Covariance_CH2(data, u, Lambda, Delta_vec, qfunc)
    worst_test_stat <- abs(kappa %*% base::t(as.matrix(test_stat - worst_ev))) / sqrt((kappa %*% worst_cov %*% kappa))
    return(worst_test_stat)
  }
  
  uinit <- rep(exp(1), I)
  #uinit <- c(rep(exp(0.5), 8), rep(exp(-0.5), 7))
  lower_bound <- rep(1/exp(1), I)
  upper_bound <- rep(exp(1), I)
  
  local_sol <- nlminb(uinit, program, lower = lower_bound, upper = upper_bound)
  
  worst_ev <- Expected_value_CH2(data, local_sol$par, Lambda, Delta_vec, qfunc)
  worst_cov <- Covariance_CH2(data, local_sol$par, Lambda, Delta_vec, qfunc)
  worst_test_stat <- abs(kappa %*% base::t(as.matrix(test_stat - worst_ev))) / sqrt((kappa %*% worst_cov %*% kappa))
  
  # kappa_func <- function(kappa){
  #   abs(kappa %*% base::t(as.matrix(test_stat - worst_ev))) / sqrt((kappa %*% worst_cov %*% kappa))
  # }
  # find_worst_test <- optim(rep(1, K), kappa_func, method = "L-BFGS-B", lower = c(0, 0), control = list(fnscale = -1))
  # worst_test_stat <- find_worst_test$value
  
  # finding the lower bound on the correlation between outcomes
  
  correlation_objective <- function(p){
    return((Yc[, 1] %*% Yc[, 2] - (Yc[, 1] * p[1:I]) %*% (Yc[, 2] * p[(I+1):(2*I)]))/(sqrt(sum(Yc[, 1]^2 * (1-p[1:I]^2))) * sqrt(sum(Yc[, 2]^2 * (1-p[(I+1):(2*I)]^2)))))
  }
  
  lb_cor <- optim(c(rep(((Lambda-1)/(1+Lambda))*((Delta_vec[1]-1)/(1+Delta_vec[1])), I), rep(((Lambda-1)/(1+Lambda))*((Delta_vec[2]-1)/(1+Delta_vec[2])), I)), 
                  correlation_objective, method = "L-BFGS-B", 
                  lower = rep(0, 2*I), 
                  upper = c(rep(((Lambda-1)/(1+Lambda))*((Delta_vec[1]-1)/(1+Delta_vec[1])), I), rep(((Lambda-1)/(1+Lambda))*((Delta_vec[2]-1)/(1+Delta_vec[2])), I)))$value
  
  # calculating distribution function with worse-case crit val
  
  M <- matrix(c(1, lb_cor, lb_cor, 1), nrow = 2, byrow = TRUE)
  Minv <- solve(M)
  rho <- Minv[1,2] / sqrt(Minv[1,1]*Minv[2,2])
  w2 <- 1/4 + asin(rho)/(2*pi)
  w0 <- 1/4 - asin(rho)/(2*pi)
  w1 <- 1 - w0 - w2
  
  Fmix <- function(x){
    w0 +
      w1 * pchisq(x,1) +
      w2 * pchisq(x,2)
  }
  
  crit_val <- uniroot(function(x) Fmix(x)-(1-alpha), interval=c(0,100))$root
  
  pval <- 1 - Fmix(worst_test_stat^2)
  
  output <- list("sensitivity_analysis_stat" = worst_test_stat, 
                 "worst_case_pvalue" = pval,
                 "worst_case_critval" = crit_val)
  return(output)
}
chibar_infCH_ttest <- function(data, Gamma, Lambda, kappa, alpha = 0.05, qfunc = "perm"){
  
  if(qfunc == "perm"){
    qf <- function(Y){
      return(abs(Y))
    }
  }
  
  I <- nrow(data)/2
  Z <- data[, 1]
  outcomes <- data[, 2:ncol(data)]
  K <- ncol(outcomes)
  test_stat <- observed_test_statistic(data, qfunc)
  Yc <- outcomes[seq(1, nrow(data), by = 2), ] - outcomes[seq(2, nrow(data), by = 2), ]
  qYc <- qf(base::t(base::t(Yc) * kappa))
  
  xL <- (Lambda - 1)/(Lambda + 1)
  lambda <- log(Lambda)
  obj_func <- function(x1){
    
    x2 <- ((K * (Gamma - 1)/(Gamma + 1)) - xL*x1)/xL
    
    Delta_vec <- c((x1 + 1)/(1 - x1), (x2 + 1)/(1 - x2))
    Delta_vec[which(Delta_vec == Inf)] <- 1000000
    delta_vec <- log(Delta_vec)
    
    worst_ev <- Expected_value_CH2(data, rep(exp(1), I), Lambda, Delta_vec, qfunc)
    worst_cov <- Covariance_CH2(data, rep(exp(1), I), Lambda, Delta_vec, qfunc)
    worst_test_stat <- abs(kappa %*% base::t(as.matrix(test_stat - worst_ev))) / sqrt((kappa %*% worst_cov %*% kappa))
    return(worst_test_stat)
  }
  
  lb <- max(0, ((K * (Gamma - 1)/(Gamma + 1)) - xL)/xL)
  ub <- min(1, (K * (Gamma - 1)/(Gamma + 1))/xL)
  opts <- list("algorithm" = "NLOPT_LN_BOBYQA", "xtol_rel" = 1.0e-5)
  
  res_list <- list()
  obj_evals <- c()
  search <- seq(lb, ub, length.out = 3)
  for(i in 1:3){
    res_list[[i]] <- nloptr(x0 = search[i], eval_f = obj_func, opts = opts, lb = lb, ub = ub)
    obj_evals[i] <- res_list[[i]]$objective
  }
  best <- which(obj_evals == min(obj_evals))[1]
  solution <- res_list[[best]]$solution
  
  edel <- c(solution, ((K * (Gamma - 1)/(Gamma + 1)) - xL*solution)/xL)
  Delta_vec <- (1+edel)/(1-edel)
  Delta_vec[which(Delta_vec == Inf)] <- 1000000
  lambda <- log(Lambda)
  delta_vec <- log(Delta_vec)
  
  program <- function(u){
    # u is an I-length vector of exp\{u_{i1}-u_{i2}\}, i=1,...,I  
    # u <- runif(I, 1/exp(1), exp(1))
    
    worst_ev <- Expected_value_CH2(data, u, Lambda, Delta_vec, qfunc)
    worst_cov <- Covariance_CH2(data, u, Lambda, Delta_vec, qfunc)
    worst_test_stat <- abs(kappa %*% base::t(as.matrix(test_stat - worst_ev))) / sqrt((kappa %*% worst_cov %*% kappa))
    return(worst_test_stat)
  }
  
  uinit <- rep(exp(1), I)
  #uinit <- c(rep(exp(0.5), 8), rep(exp(-0.5), 7))
  lower_bound <- rep(1/exp(1), I)
  upper_bound <- rep(exp(1), I)
  
  local_sol <- nlminb(uinit, program, lower = lower_bound, upper = upper_bound)
  worst_test_stat <- local_sol$objective
  
  correlation_objective <- function(p){
    return((Yc[, 1] %*% Yc[, 2] - (Yc[, 1] * p[1:I]) %*% (Yc[, 2] * p[(I+1):(2*I)]))/(sqrt(sum(Yc[, 1]^2 * (1-p[1:I]^2))) * sqrt(sum(Yc[, 2]^2 * (1-p[(I+1):(2*I)]^2)))))
  }
  
  upper_lim <- c(((Lambda-1)/(1+Lambda))*((Delta_vec[1]-1)/(1+Delta_vec[1])), ((Lambda-1)/(1+Lambda))*((Delta_vec[2]-1)/(1+Delta_vec[2])))
  upper_lim[which(upper_lim == 0)] <- 0.0001
  
  lb_cor <- optim(c(rep(((Lambda-1)/(1+Lambda))*((Delta_vec[1]-1)/(1+Delta_vec[1])), I), rep(((Lambda-1)/(1+Lambda))*((Delta_vec[2]-1)/(1+Delta_vec[2])), I)), 
                  correlation_objective, method = "L-BFGS-B", 
                  lower = rep(0, 2*I), 
                  upper = c(rep(upper_lim[1], I), rep(upper_lim[2], I)))$value
  
  # calculating distribution function with worse-case crit val
  
  M <- matrix(c(1, lb_cor, lb_cor, 1), nrow = 2, byrow = TRUE)
  Minv <- solve(M)
  rho <- Minv[1,2] / sqrt(Minv[1,1]*Minv[2,2])
  w2 <- 1/4 + asin(rho)/(2*pi)
  w0 <- 1/4 - asin(rho)/(2*pi)
  w1 <- 1 - w0 - w2
  
  Fmix <- function(x){
    w0 +
      w1 * pchisq(x,1) +
      w2 * pchisq(x,2)
  }
  
  crit_val <- uniroot(function(x) Fmix(x)-(1-alpha), interval=c(0,100))$root
  
  pval <- 1 - Fmix(worst_test_stat^2)
  
  output <- list("sensitivity_analysis_stat" = worst_test_stat, 
                 "worst_case_pvalue" = pval,
                 "worst_case_critval" = crit_val,
                 "Delta_vec" = Delta_vec)
  return(output)
}
chibar_supCH_ttest <- function(data, Gamma, Lambda, kappa, alpha = 0.05, qfunc = "perm"){
  
  if(qfunc == "perm"){
    qf <- function(Y){
      return(abs(Y))
    }
  }
  
  I <- nrow(data)/2
  Z <- data[, 1]
  outcomes <- data[, 2:ncol(data)]
  K <- ncol(outcomes)
  test_stat <- observed_test_statistic(data, qfunc)
  Yc <- outcomes[seq(1, nrow(data), by = 2), ] - outcomes[seq(2, nrow(data), by = 2), ]
  qYc <- qf(base::t(base::t(Yc) * kappa))
  
  xL <- (Lambda - 1)/(Lambda + 1)
  lambda <- log(Lambda)
  obj_func <- function(x1){
    
    x2 <- ((K * (Gamma - 1)/(Gamma + 1)) - xL*x1)/xL
    
    Delta_vec <- c((x1 + 1)/(1 - x1), (x2 + 1)/(1 - x2))
    Delta_vec[which(Delta_vec == Inf)] <- 1000000
    delta_vec <- log(Delta_vec)
    
    program <- function(u){
      worst_ev <- Expected_value_CH2(data, u, Lambda, Delta_vec, qfunc)
      worst_cov <- Covariance_CH2(data, u, Lambda, Delta_vec, qfunc)
      worst_test_stat <- abs(kappa %*% base::t(as.matrix(test_stat - worst_ev))) / sqrt((kappa %*% worst_cov %*% kappa))
      return(worst_test_stat)
    }
    
    uinit <- rep(exp(1), I)
    lower_bound <- rep(1/exp(1), I)
    upper_bound <- rep(exp(1), I)
    
    local_sol <- nlminb(uinit, program, lower = lower_bound, upper = upper_bound)
    
    worst_ev <- Expected_value_CH2(data, local_sol$par, Lambda, Delta_vec, qfunc)
    worst_cov <- Covariance_CH2(data, local_sol$par, Lambda, Delta_vec, qfunc)
    worst_test_stat <- abs(kappa %*% base::t(as.matrix(test_stat - worst_ev))) / sqrt((kappa %*% worst_cov %*% kappa))
    return(-worst_test_stat)
  }
  
  lb <- max(0, ((K * (Gamma - 1)/(Gamma + 1)) - xL)/xL)
  ub <- min(1, (K * (Gamma - 1)/(Gamma + 1))/xL)
  opts <- list("algorithm" = "NLOPT_LN_BOBYQA", "xtol_rel" = 1.0e-5)
  
  res_list <- list()
  obj_evals <- c()
  search <- seq(lb, ub, length.out = 3)
  for(i in 1:3){
    res_list[[i]] <- nloptr(x0 = search[i], eval_f = obj_func, opts = opts, lb = lb, ub = ub)
    obj_evals[i] <- res_list[[i]]$objective
  }
  best <- which(obj_evals == min(obj_evals))[1]
  solution <- res_list[[best]]$solution
  
  edel <- c(solution, ((K * (Gamma - 1)/(Gamma + 1)) - xL*solution)/xL)
  Delta_vec <- (1+edel)/(1-edel)
  Delta_vec[which(Delta_vec == Inf)] <- 1000000
  lambda <- log(Lambda)
  delta_vec <- log(Delta_vec)
  
  program <- function(u){
    worst_ev <- Expected_value_CH2(data, u, Lambda, Delta_vec, qfunc)
    worst_cov <- Covariance_CH2(data, u, Lambda, Delta_vec, qfunc)
    worst_test_stat <- abs(kappa %*% base::t(as.matrix(test_stat - worst_ev))) / sqrt((kappa %*% worst_cov %*% kappa))
    return(worst_test_stat)
  }
  
  uinit <- rep(exp(1), I)
  lower_bound <- rep(1/exp(1), I)
  upper_bound <- rep(exp(1), I)
  
  local_sol <- nlminb(uinit, program, lower = lower_bound, upper = upper_bound)
  worst_test_stat <- local_sol$objective
  
  correlation_objective <- function(p){
    return((Yc[, 1] %*% Yc[, 2] - (Yc[, 1] * p[1:I]) %*% (Yc[, 2] * p[(I+1):(2*I)]))/(sqrt(sum(Yc[, 1]^2 * (1-p[1:I]^2))) * sqrt(sum(Yc[, 2]^2 * (1-p[(I+1):(2*I)]^2)))))
  }
  
  upper_lim <- c(((Lambda-1)/(1+Lambda))*((Delta_vec[1]-1)/(1+Delta_vec[1])), ((Lambda-1)/(1+Lambda))*((Delta_vec[2]-1)/(1+Delta_vec[2])))
  upper_lim[which(upper_lim == 0)] <- 0.0001
  
  lb_cor <- optim(c(rep(((Lambda-1)/(1+Lambda))*((Delta_vec[1]-1)/(1+Delta_vec[1])), I), rep(((Lambda-1)/(1+Lambda))*((Delta_vec[2]-1)/(1+Delta_vec[2])), I)), 
                  correlation_objective, method = "L-BFGS-B", 
                  lower = rep(0, 2*I), 
                  upper = c(rep(upper_lim[1], I), rep(upper_lim[2], I)))$value
  
  # calculating distribution function with worse-case crit val
  
  M <- matrix(c(1, lb_cor, lb_cor, 1), nrow = 2, byrow = TRUE)
  Minv <- solve(M)
  rho <- Minv[1,2] / sqrt(Minv[1,1]*Minv[2,2])
  w2 <- 1/4 + asin(rho)/(2*pi)
  w0 <- 1/4 - asin(rho)/(2*pi)
  w1 <- 1 - w0 - w2
  
  Fmix <- function(x){
    w0 +
      w1 * pchisq(x,1) +
      w2 * pchisq(x,2)
  }
  
  crit_val <- uniroot(function(x) Fmix(x)-(1-alpha), interval=c(0,100))$root
  
  pval <- 1 - Fmix(worst_test_stat^2)
  
  output <- list("sensitivity_analysis_stat" = worst_test_stat, 
                 "worst_case_pvalue" = pval,
                 "worst_case_critval" = crit_val,
                 "Delta_vec" = Delta_vec)
  return(output)
}

# Practitioner acts first
chibar_kappa_finder <- function(data){
  
  I <- nrow(data)/2
  K <- ncol(data) - 1
  
  test_stat <- observed_test_statistic(data, qfunc = "perm")
  worst_ev <- Expected_value_Rosenbaum(data, rep(0, I), qfunc = "perm")
  worst_cov <- Covariance_Rosenbaum(data, rep(0, I), qfunc = "perm")
  
  kappa_func <- function(kappa){
    abs(kappa %*% base::t(as.matrix(test_stat - worst_ev))) / sqrt((kappa %*% worst_cov %*% kappa))
  }
  find_worst_test <- optim(rep(0.1, K), kappa_func, method = "L-BFGS-B", lower = c(0, 0), control = list(fnscale = -1))
  kappa_unscaled <- find_worst_test$par
  kappa <- kappa_unscaled / sum(abs(kappa_unscaled))
  return(kappa)
  
}
fixed_kappa_noconfounding_classic <- function(data, kappa, alpha = 0.05, qfunc = "perm"){
  
  if(qfunc == "perm"){
    qf <- function(Y){
      return(abs(Y))
    }
  }
  I <- nrow(data)/2
  # Z <- data[, 1]
  # outcomes <- data[, 2:ncol(data)]
  # K <- ncol(outcomes)
  # Yc <- outcomes[seq(1, nrow(data), by = 2), ] - outcomes[seq(2, nrow(data), by = 2), ] # treated - control
  # Yckappa <- base::t(base::t(Yc) * kappa)
  
  # calculate the test statistic
  observed_test_stat <- observed_test_statistic(data, qfunc)
  observed_cov <- Covariance_Rosenbaum(data, rep(0, I), qfunc)
  scaled_test_stat <- abs(kappa %*% base::t(as.matrix(observed_test_stat))) / sqrt((kappa %*% observed_cov %*% kappa))
  
  # calculate the critical value
  crit_val <- sqrt(qchisq(1-alpha, 1))
  pval <- 1-pchisq(scaled_test_stat^2, 1)
  
  output <- list("Test_stat" = scaled_test_stat,
                 "p_value" = pval, 
                 "critical_val" = crit_val)
  return(output)
  
}
fixed_kappa_noconfounding_amp <- function(data, kappa, alpha = 0.05, qfunc = "perm"){
  
  if(qfunc == "perm"){
    qf <- function(Y){
      return(abs(Y))
    }
  }
  
  I <- nrow(data)/2
  #Z <- data[, 1]
  #outcomes <- data[, 2:ncol(data)]
  #K <- ncol(outcomes)
  #Yc <- outcomes[seq(1, nrow(data), by = 2), ] - outcomes[seq(2, nrow(data), by = 2), ] # treated - control
  #Yckappa <- base::t(base::t(Yc) * kappa)
  
  # calculate the test statistic
  observed_test_stat <- observed_test_statistic(data, qfunc)
  observed_cov <- Covariance_CH2(data, rep(1, I), 1, c(1, 1), qfunc)
  scaled_test_stat <- abs(kappa %*% base::t(as.matrix(observed_test_stat))) / sqrt((kappa %*% observed_cov %*% kappa))
  
  # calcualte the crit val
  crit_val <- sqrt(qchisq(1-alpha, 1))
  pval <- 1-pchisq(scaled_test_stat^2, 1)
  
  output <- list("Test_stat" = scaled_test_stat,
                 "p_value" = pval, 
                 "critical_val" = crit_val)
  return(output)
  
}
fixed_kappa_T2CH_ttest <- function(data, Lambda, Delta_vec, kappa, alpha = 0.05, qfunc = "perm"){
  
  if(qfunc == "perm"){
    qf <- function(Y){
      return(abs(Y))
    }
  }
  
  lambda <- log(Lambda)
  delta_vec <- log(Delta_vec)
  I <- nrow(data)/2
  # Z <- data[, 1]
  # outcomes <- data[, 2:ncol(data)]
  # K <- ncol(outcomes)
  test_stat <- observed_test_statistic(data, qfunc)
  # Yc <- outcomes[seq(1, nrow(data), by = 2), ] - outcomes[seq(2, nrow(data), by = 2), ]
  # qYc <- qf(base::t(base::t(Yc) * kappa))
  
  program <- function(u){
    # u is an I-length vector of exp\{u_{i1}-u_{i2}\}, i=1,...,I  
    # u <- runif(I, 1/exp(1), exp(1))
    
    worst_ev <- Expected_value_CH2(data, u, Lambda, Delta_vec, qfunc)
    worst_cov <- Covariance_CH2(data, u, Lambda, Delta_vec, qfunc)
    worst_test_stat <- abs(kappa %*% base::t(as.matrix(test_stat - worst_ev))) / sqrt((kappa %*% worst_cov %*% kappa))
    return(worst_test_stat)
  }
  
  uinit <- rep(exp(1), I)
  #uinit <- c(rep(exp(0.5), 8), rep(exp(-0.5), 7))
  lower_bound <- rep(1/exp(1), I)
  upper_bound <- rep(exp(1), I)
  
  local_sol <- nlminb(uinit, program, lower = lower_bound, upper = upper_bound)
  
  worst_ev <- Expected_value_CH2(data, local_sol$par, Lambda, Delta_vec, qfunc)
  worst_cov <- Covariance_CH2(data, local_sol$par, Lambda, Delta_vec, qfunc)
  worst_test_stat <- abs(kappa %*% base::t(as.matrix(test_stat - worst_ev))) / sqrt((kappa %*% worst_cov %*% kappa))
  
  # kappa_func <- function(kappa){
  #   abs(kappa %*% base::t(as.matrix(test_stat - worst_ev))) / sqrt((kappa %*% worst_cov %*% kappa))
  # }
  # find_worst_test <- optim(rep(1, K), kappa_func, method = "L-BFGS-B", lower = c(0, 0), control = list(fnscale = -1))
  # worst_test_stat <- find_worst_test$value
  
  # finding the lower bound on the correlation between outcomes
  
  # calculate the critical value
  crit_val <- sqrt(qchisq(1-alpha, 1))
  pval <- 1-pchisq(worst_test_stat^2, 1)
  
  output <- list("sensitivity_analysis_stat" = worst_test_stat, 
                 "worst_case_pvalue" = pval,
                 "worst_case_critval" = crit_val)
  return(output)
}
fixed_kappa_infCH_ttest <- function(data, Gamma, Lambda, kappa, alpha = 0.05, qfunc = "perm"){
  
  if(qfunc == "perm"){
    qf <- function(Y){
      return(abs(Y))
    }
  }
  
  I <- nrow(data)/2
  # Z <- data[, 1]
  # outcomes <- data[, 2:ncol(data)]
  K <- ncol(outcomes)
  test_stat <- observed_test_statistic(data, qfunc)
  # Yc <- outcomes[seq(1, nrow(data), by = 2), ] - outcomes[seq(2, nrow(data), by = 2), ]
  # qYc <- qf(base::t(base::t(Yc) * kappa))
  
  xL <- (Lambda - 1)/(Lambda + 1)
  lambda <- log(Lambda)
  obj_func <- function(x1){
    
    x2 <- ((K * (Gamma - 1)/(Gamma + 1)) - xL*x1)/xL
    
    Delta_vec <- c((x1 + 1)/(1 - x1), (x2 + 1)/(1 - x2))
    Delta_vec[which(Delta_vec == Inf)] <- 1000000
    delta_vec <- log(Delta_vec)
    
    worst_ev <- Expected_value_CH2(data, rep(exp(1), I), Lambda, Delta_vec, qfunc)
    worst_cov <- Covariance_CH2(data, rep(exp(1), I), Lambda, Delta_vec, qfunc)
    worst_test_stat <- abs(kappa %*% base::t(as.matrix(test_stat - worst_ev))) / sqrt((kappa %*% worst_cov %*% kappa))
    return(worst_test_stat)
  }
  
  lb <- max(0, ((K * (Gamma - 1)/(Gamma + 1)) - xL)/xL)
  ub <- min(1, (K * (Gamma - 1)/(Gamma + 1))/xL)
  opts <- list("algorithm" = "NLOPT_LN_BOBYQA", "xtol_rel" = 1.0e-5)
  
  res_list <- list()
  obj_evals <- c()
  search <- seq(lb, ub, length.out = 3)
  for(i in 1:3){
    res_list[[i]] <- nloptr(x0 = search[i], eval_f = obj_func, opts = opts, lb = lb, ub = ub)
    obj_evals[i] <- res_list[[i]]$objective
  }
  best <- which(obj_evals == min(obj_evals))[1]
  solution <- res_list[[best]]$solution
  
  edel <- c(solution, ((K * (Gamma - 1)/(Gamma + 1)) - xL*solution)/xL)
  Delta_vec <- (1+edel)/(1-edel)
  Delta_vec[which(Delta_vec == Inf)] <- 1000000
  lambda <- log(Lambda)
  delta_vec <- log(Delta_vec)
  
  program <- function(u){
    # u is an I-length vector of exp\{u_{i1}-u_{i2}\}, i=1,...,I  
    # u <- runif(I, 1/exp(1), exp(1))
    
    worst_ev <- Expected_value_CH2(data, u, Lambda, Delta_vec, qfunc)
    worst_cov <- Covariance_CH2(data, u, Lambda, Delta_vec, qfunc)
    worst_test_stat <- abs(kappa %*% base::t(as.matrix(test_stat - worst_ev))) / sqrt((kappa %*% worst_cov %*% kappa))
    return(worst_test_stat)
  }
  
  uinit <- rep(exp(1), I)
  #uinit <- c(rep(exp(0.5), 8), rep(exp(-0.5), 7))
  lower_bound <- rep(1/exp(1), I)
  upper_bound <- rep(exp(1), I)
  
  local_sol <- nlminb(uinit, program, lower = lower_bound, upper = upper_bound)
  worst_test_stat <- local_sol$objective
  
  # calculate the critical value
  crit_val <- sqrt(qchisq(1-alpha, 1))
  pval <- 1-pchisq(worst_test_stat^2, 1)
  
  output <- list("sensitivity_analysis_stat" = worst_test_stat, 
                 "worst_case_pvalue" = pval,
                 "worst_case_critval" = crit_val,
                 "Delta_vec" = Delta_vec)
  return(output)
}
fixed_kappa_supCH_ttest <- function(data, Gamma, Lambda, kappa, alpha = 0.05, qfunc = "perm"){
  
  if(qfunc == "perm"){
    qf <- function(Y){
      return(abs(Y))
    }
  }
  
  I <- nrow(data)/2
  # Z <- data[, 1]
  # outcomes <- data[, 2:ncol(data)]
  K <- ncol(outcomes)
  test_stat <- observed_test_statistic(data, qfunc)
  # Yc <- outcomes[seq(1, nrow(data), by = 2), ] - outcomes[seq(2, nrow(data), by = 2), ]
  # qYc <- qf(base::t(base::t(Yc) * kappa))
  
  xL <- (Lambda - 1)/(Lambda + 1)
  lambda <- log(Lambda)
  obj_func <- function(x1){
    
    x2 <- ((K * (Gamma - 1)/(Gamma + 1)) - xL*x1)/xL
    
    Delta_vec <- c((x1 + 1)/(1 - x1), (x2 + 1)/(1 - x2))
    Delta_vec[which(Delta_vec == Inf)] <- 1000000
    delta_vec <- log(Delta_vec)
    
    program <- function(u){
      worst_ev <- Expected_value_CH2(data, u, Lambda, Delta_vec, qfunc)
      worst_cov <- Covariance_CH2(data, u, Lambda, Delta_vec, qfunc)
      worst_test_stat <- abs(kappa %*% base::t(as.matrix(test_stat - worst_ev))) / sqrt((kappa %*% worst_cov %*% kappa))
      return(worst_test_stat)
    }
    
    uinit <- rep(exp(1), I)
    lower_bound <- rep(1/exp(1), I)
    upper_bound <- rep(exp(1), I)
    
    local_sol <- nlminb(uinit, program, lower = lower_bound, upper = upper_bound)
    
    # worst_ev <- Expected_value_CH2(data, local_sol$par, Lambda, Delta_vec, qfunc)
    # worst_cov <- Covariance_CH2(data, local_sol$par, Lambda, Delta_vec, qfunc)
    # worst_test_stat <- abs(kappa %*% base::t(as.matrix(test_stat - worst_ev))) / sqrt((kappa %*% worst_cov %*% kappa))
    return(-local_sol$objective)
  }
  
  lb <- max(0, ((K * (Gamma - 1)/(Gamma + 1)) - xL)/xL)
  ub <- min(1, (K * (Gamma - 1)/(Gamma + 1))/xL)
  
  low_sol <- obj_func(lb)
  up_sol <- obj_func(ub)
  solution <- lb
  if(up_sol < low_sol){
    solution <- ub
  }
  
  # opts <- list("algorithm" = "NLOPT_LN_BOBYQA", "xtol_rel" = 1.0e-5)
  # 
  # res_list <- list()
  # obj_evals <- c()
  # search <- seq(lb, ub, length.out = 3)
  # for(i in 1:3){
  #   res_list[[i]] <- nloptr(x0 = search[i], eval_f = obj_func, opts = opts, lb = lb, ub = ub)
  #   obj_evals[i] <- res_list[[i]]$objective
  # }
  # best <- which(obj_evals == min(obj_evals))[1]
  # solution <- res_list[[best]]$solution
  
  edel <- c(solution, ((K * (Gamma - 1)/(Gamma + 1)) - xL*solution)/xL)
  Delta_vec <- (1+edel)/(1-edel)
  Delta_vec[which(Delta_vec == Inf)] <- 1000000
  lambda <- log(Lambda)
  delta_vec <- log(Delta_vec)
  
  program <- function(u){
    worst_ev <- Expected_value_CH2(data, u, Lambda, Delta_vec, qfunc)
    worst_cov <- Covariance_CH2(data, u, Lambda, Delta_vec, qfunc)
    worst_test_stat <- abs(kappa %*% base::t(as.matrix(test_stat - worst_ev))) / sqrt((kappa %*% worst_cov %*% kappa))
    return(worst_test_stat)
  }
  
  uinit <- rep(exp(1), I)
  lower_bound <- rep(1/exp(1), I)
  upper_bound <- rep(exp(1), I)
  
  local_sol <- nlminb(uinit, program, lower = lower_bound, upper = upper_bound)
  worst_test_stat <- local_sol$objective
  
  # calculate the critical value
  crit_val <- sqrt(qchisq(1-alpha, 1))
  pval <- 1-pchisq(worst_test_stat^2, 1)
  
  output <- list("sensitivity_analysis_stat" = worst_test_stat, 
                 "worst_case_pvalue" = pval,
                 "worst_case_critval" = crit_val,
                 "Delta_vec" = Delta_vec)
  return(output)
}
fixed_kappa_supCH_unitkappa_ttest <- function(data, Gamma, Lambda, kappa, alpha = 0.05, qfunc = "perm"){
  
  if(qfunc == "perm"){
    qf <- function(Y){
      return(abs(Y))
    }
  }
  
  I <- nrow(data)/2
  # Z <- data[, 1]
  # outcomes <- data[, 2:ncol(data)]
  K <- ncol(outcomes)
  test_stat <- observed_test_statistic(data, qfunc)
  # Yc <- outcomes[seq(1, nrow(data), by = 2), ] - outcomes[seq(2, nrow(data), by = 2), ]
  # qYc <- qf(base::t(base::t(Yc) * kappa))
  
  xL <- (Lambda - 1)/(Lambda + 1)
  lambda <- log(Lambda)
  
  lb <- max(0, ((K * (Gamma - 1)/(Gamma + 1)) - xL)/xL)
  ub <- min(1, (K * (Gamma - 1)/(Gamma + 1))/xL)
  if(kappa[1] == 0){
    solution <- ub
  }
  if(kappa[1] == 1){
    solution <- lb
  }
  
  edel <- c(solution, ((K * (Gamma - 1)/(Gamma + 1)) - xL*solution)/xL)
  Delta_vec <- (1+edel)/(1-edel)
  Delta_vec[which(Delta_vec == Inf)] <- 1000000
  delta_vec <- log(Delta_vec)
  
  program <- function(u){
    worst_ev <- Expected_value_CH2(data, u, Lambda, Delta_vec, qfunc)
    worst_cov <- Covariance_CH2(data, u, Lambda, Delta_vec, qfunc)
    worst_test_stat <- abs(kappa %*% base::t(as.matrix(test_stat - worst_ev))) / sqrt((kappa %*% worst_cov %*% kappa))
    return(worst_test_stat)
  }
  
  uinit <- rep(exp(1), I)
  lower_bound <- rep(1/exp(1), I)
  upper_bound <- rep(exp(1), I)
  
  local_sol <- nlminb(uinit, program, lower = lower_bound, upper = upper_bound)
  worst_test_stat <- local_sol$objective
  
  # calculate the critical value
  crit_val <- sqrt(qchisq(1-alpha, 1))
  pval <- 1-pchisq(worst_test_stat^2, 1)
  
  output <- list("sensitivity_analysis_stat" = worst_test_stat, 
                 "worst_case_pvalue" = pval,
                 "worst_case_critval" = crit_val,
                 "Delta_vec" = Delta_vec)
  return(output)
}
fixed_kappa_classic_ttest <- function(data, Gamma, kappa, alpha = 0.05, qfunc = "perm"){
  
  if(qfunc == "perm"){
    qf <- function(Y){
      return(abs(Y))
    }
  }
  
  I <- nrow(data)/2
  Z <- data[, 1]
  outcomes <- data[, 2:ncol(data)]
  K <- ncol(outcomes)
  test_stat <- observed_test_statistic(data, qfunc = "perm")
  Yc <- outcomes[seq(1, nrow(data), by = 2), ] - outcomes[seq(2, nrow(data), by = 2), ] # treated - control
  Yckappa <- base::t(base::t(Yc) * kappa)
  
  # calculating worst-case u
  
  program <- function(pu){
    # u is an I-length vector of exp\{u_{i1}-u_{i2}\}, i=1,...,I  
    # u <- runif(I, 1/exp(1), exp(1))
    
    worst_ev <- Expected_value_Rosenbaum(data, pu, qfunc)
    worst_cov <- Covariance_Rosenbaum(data, pu, qfunc)
    worst_test_stat <- abs(kappa %*% base::t(as.matrix(test_stat - worst_ev))) / sqrt((kappa %*% worst_cov %*% kappa))
    return(worst_test_stat)
  }
  
  #uinit <- c(rep(exp(0.5), 8), rep(exp(-0.5), 7))
  lower_bound <- rep((1 - Gamma)/(Gamma + 1), I)
  upper_bound <- rep((Gamma - 1)/(Gamma + 1), I)
  uinit <- upper_bound*sign(rowSums(Yckappa))
  
  local_sol <- nlminb(uinit, program, lower = lower_bound, upper = upper_bound)
  worst_test_stat <- local_sol$objective
  
  #worst_ev <- Expected_value_Rosenbaum(data, pu, qfunc)
  #worst_cov <- Covariance_Rosenbaum(data, pu, qfunc)
  #worst_test_stat <- abs(kappa %*% base::t(as.matrix(test_stat - worst_ev))) / sqrt((kappa %*% worst_cov %*% kappa))
  
  # kappa_func <- function(kappa){
  #   abs(kappa %*% base::t(as.matrix(test_stat - worst_ev))) / sqrt((kappa %*% worst_cov %*% kappa))
  # }
  # find_worst_test <- optim(rep(1, K), kappa_func, method = "L-BFGS-B", lower = c(0, 0), control = list(fnscale = -1))
  # worst_test_stat <- find_worst_test$value
  
  # calculate the critical value
  crit_val <- sqrt(qchisq(1-alpha, 1))
  pval <- 1-pchisq(worst_test_stat^2, 1)
  
  output <- list("sensitivity_analysis_stat" = worst_test_stat, 
                 "worst_case_pvalue" = pval,
                 "worst_case_critval" = crit_val)
  return(output)
}

#### Sensitivity Value, Unobserved Confounder Acts First (main) ----------------

summary_df <- data.frame(matrix(nrow = 0, ncol = 9))

# classic: gamma = 2.76, kappa = c(0, 1)
classic <- chibar_classic_ttest_adv(data, 2.76)
summary_df <- rbind(summary_df, c("Classic", "NA", classic$kappa, 2.76, "NA", "NA", "NA", round(classic$worst_case_pvalue, 3)))


# amp (a=1/3): Gamma = 2.52, Lambda = 9.77, Delta1 = Delta2 = 3.26, kappa = c(0, 1)
Gamma <- 2.52
a <- 1/3
Pars <- Scaled_amp_parameters(Gamma, a)
amp_athird <- chibar_T2CH_ttest_adv(data, Pars[1], rep(Pars[2], 2))
summary_df <- rbind(summary_df, c("Amplification", 0.33, amp_athird$kappa, 2.52, 9.77, 3.26, 3.26, round(amp_athird$worst_case_pvalue, 3)))

# amp (a=1): Gamma = 2.52, Lambda = 4.83, Delta1 = Delta2 = 4.83, kappa = c(0, 1)
Gamma <- 2.52
a <- 1
Pars <- Scaled_amp_parameters(Gamma, a)
amp_a1 <- chibar_T2CH_ttest_adv(data, Pars[1], rep(Pars[2], 2))
summary_df <- rbind(summary_df, c("Amplification", 1, amp_a1$kappa, 2.52, 4.83, 4.83, 4.83, round(amp_a1$worst_case_pvalue, 3)))

# amp (a=3): Gamma = 2.52, Lambda = 3.26, Delta1 = Delta2 = 9.77, kappa = c(0, 1)
Gamma <- 2.52
a <- 3
Pars <- Scaled_amp_parameters(Gamma, a)
amp_a3 <- chibar_T2CH_ttest_adv(data, Pars[1], rep(Pars[2], 2))
summary_df <- rbind(summary_df, c("Amplification", 3, amp_a1$kappa, 2.52, 3.26, 9.77, 9.77, round(amp_a3$worst_case_pvalue, 3)))

# inf (a=1/3): Gamma = 2.38, Lambda = 9.19, Delta1 = 2.50, Delta2 = 3.85, kappa = c(0.55, 0.45)
Gamma <- 2.38
a <- 1/3
Pars <- Scaled_amp_parameters(Gamma, a)
inf_athird <- chibar_infCH_ttest_adv(data, Gamma, Pars[1])
summary_df <- rbind(summary_df, c("Inf_CH", 0.33, inf_athird$kappa, 2.38, 9.19, 2.50, 3.85, round(inf_athird$worst_case_pvalue, 3)))

# inf (a=1): Gamma = 2.35, Lambda = 4.48, Delta1 = 3.31, Delta2 = 6.51, kappa = c(0.54, 0.46)
Gamma <- 2.35
a <- 1
Pars <- Scaled_amp_parameters(Gamma, a)
inf_a1 <- chibar_infCH_ttest_adv(data, Gamma, Pars[1])
summary_df <- rbind(summary_df, c("Inf_CH", 1, inf_a1$kappa, 2.35, 4.48, 3.31, 6.51, round(inf_a1$worst_case_pvalue, 3)))

# inf (a=3): Gamma = 2.30, Lambda = 2.95, Delta1 = 5.11, Delta2 = 24.49, kappa = c(0.54, 0.46)
Gamma <- 2.30
a <- 3
Pars <- Scaled_amp_parameters(Gamma, a)
inf_a3 <- chibar_infCH_ttest_adv(data, Gamma, Pars[1])
summary_df <- rbind(summary_df, c("Inf_CH", 3, inf_a3$kappa, 2.30, 2.95, 5.11, 24.49, round(inf_a3$worst_case_pvalue, 3)))

# sup (a=1/3): Gamma = 8.89, Lambda = 35.48, Delta1 = infty, Delta2 = 5.41, kappa = c(0, 1)
Gamma <- 8.89
a <- 1/3
Pars <- Scaled_amp_parameters(Gamma, a)
sup_athird <- chibar_supCH_ttest_adv_biggamma(data, Gamma, Pars[1])
summary_df <- rbind(summary_df, c("Sup_CH", 0.33, sup_athird$kappa, 8.89, 35.48, "infinity", 5.41, round(sup_athird$worst_case_pvalue, 3)))

# sup (a=1): Gamma = 5.92, Lambda = 11.75, Delta1 = infty, Delta2 = 5.38, kappa = c(0, 1)
Gamma <- 5.92
a <- 1
Pars <- Scaled_amp_parameters(Gamma, a)
sup_a1 <- hibar_supCH_ttest_adv_biggamma(data, Gamma, Pars[1])
summary_df <- rbind(summary_df, c("Sup_CH", 1, sup_a1$kappa, 5.92, 11.75, "infinity", 5.38, round(sup_a1$worst_case_pvalue, 3)))

# sup (a=3): Gamma = 3.29, Lambda = 4.31, Delta1 = infty, Delta2 = 5.96, kappa = c(0, 1)
Gamma <- 3.29
a <- 3
Pars <- Scaled_amp_parameters(Gamma, a)
sup_a3 <- chibar_supCH_ttest_adv_biggamma(data, Gamma, Pars[1])
summary_df <- rbind(summary_df, c("Sup_CH", 3, sup_a3$kappa, 3.29, 4.31, "infinity", 5.96, round(sup_a3$worst_case_pvalue, 3)))

names(summary_df) <- c("Framework", "a", "kappa1", "kappa2", "Gamma", "Lambda", "Delta1", "Delta2", "pvalue")
write.csv(summary_df, "sensitivity_values_main.csv")

#### Sensitivity Value, Practitioner Acts First (supplement) -------------------

summary_df <- data.frame(matrix(nrow = 0, ncol = 9))
kappa <- chibar_kappa_finder(data)

# classic: gamma = 2.60
classic <- chibar_classic_ttest(data, 2.60, kappa)
summary_df <- rbind(summary_df, c("Classic", "NA", kappa, 2.60, "NA", "NA", "NA", round(classic$worst_case_pvalue, 3)))

# amp (a=1/3): Gamma = 2.51, Lambda = 9.73, Delta1 = Delta2 = 3.24
a <- 1/3
Gamma <- 2.51
pars <- Scaled_amp_parameters(Gamma, a)
amp_athird <- chibar_T2CH_ttest(data, pars[1], rep(pars[2], 2), kappa)
summary_df <- rbind(summary_df, c("Amplification", 0.33, kappa, 2.51, 9.73, 3.24, 3.24, round(amp_athird$worst_case_pvalue, 3)))

# amp (a=1): Gamma = 2.49, Lambda = 4.79, Delta1 = Delta2 = 4.79
a <- 1
Gamma <- 2.49
pars <- Scaled_amp_parameters(Gamma, a)
amp_a1 <- chibar_T2CH_ttest(data, pars[1], rep(pars[2], 2), kappa)
summary_df <- rbind(summary_df, c("Amplification", 1, kappa, 2.49, 4.79, 4.79, 4.79, round(amp_a1$worst_case_pvalue, 3)))

# amp (a=3): Gamma = 2.46, Lambda = 3.18, Delta1 = Delta2 = 9.53
a <- 3
Gamma <- 2.46
pars <- Scaled_amp_parameters(Gamma, a)
amp_a3 <- chibar_T2CH_ttest(data, pars[1], rep(pars[2], 2), kappa)
summary_df <- rbind(summary_df, c("Amplification", 3, kappa, 2.46, 3.18, 9.53, 9.53, round(amp_a3$worst_case_pvalue, 3)))

# inf (a=1/3): Gamma = 1.72, Lambda = 6.41, Delta1 = 1, Delta2 = 6.27
a <- 1/3
Gamma <- 1.72
pars <- Scaled_amp_parameters(Gamma, a)
inf_athird <- chibar_infCH_ttest(data, Gamma, pars[1], kappa)
summary_df <- rbind(summary_df, c("Inf_CH", 0.33, kappa, 1.72, 6.41, 1, 6.27, round(inf_athird$worst_case_pvalue, 3)))

# inf (a=1): Gamma = 1.76, Lambda = 3.208, Delta1 = 1.10, Delta2 = infty
a <- 1
Gamma <- 1.76
pars <- Scaled_amp_parameters(Gamma, a)
inf_a1 <- chibar_infCH_ttest(data, Gamma, pars[1], kappa)
summary_df <- rbind(summary_df, c("Inf_CH", 1, kappa, 1.76, 3.208, 1.10, "Infinity", round(inf_a1$worst_case_pvalue, 3)))

# inf (a=3): Gamma = 2.12, Lambda = 2.70, Delta1 = 3.56, Delta2 = infty
a <- 3
Gamma <- 2.12
pars <- Scaled_amp_parameters(Gamma, a)
inf_a3 <- chibar_infCH_ttest(data, Gamma, pars[1], kappa)
summary_df <- rbind(summary_df, c("Inf_CH", 3, kappa, 2.12, 2.70, 3.56, "Infinity", round(inf_a3$worst_case_pvalue, 3)))

# sup (a=1/3): Gamma = 4.10, Lambda = 16.21, Delta1 = infty, Delta2 = 2.20
a <- 1/3
Gamma <- 4.1
pars <- Scaled_amp_parameters(Gamma, a)
sup_athird <- chibar_supCH_ttest(data, Gamma, pars[1], kappa)
summary_df <- rbind(summary_df, c("Sup_CH", 0.33, kappa, 4.10, 16.21, "Infinity", 2.20, round(sup_athird$worst_case_pvalue, 3)))

# sup (a=1): Gamma = 3.45, Lambda = 6.75, Delta1 = infty, Delta2 = 2.88
a <- 1
Gamma <- 3.45
pars <- Scaled_amp_parameters(Gamma, a)
sup_a1 <- chibar_supCH_ttest(data, Gamma, pars[1], kappa)
summary_df <- rbind(summary_df, c("Sup_CH", 1, kappa, 3.45, 6.75, "Infinity", 2.88, round(sup_a1$worst_case_pvalue, 3)))

# sup (a=3): Gamma = 2.87, Lambda = 3.73, Delta1 = infty, Delta2 = 5.11
a <- 3
Gamma <- 2.87
pars <- Scaled_amp_parameters(Gamma, a)
sup_a3 <- chibar_supCH_ttest(data, Gamma, pars[1], kappa)
summary_df <- rbind(summary_df, c("Sup_CH", 3, kappa, 2.87, 3.73, "Infinity", 5.1, round(sup_a3$worst_case_pvalue, 3)))

names(summary_df) <- c("Framework", "a", "kappa1", "kappa2", "Gamma", "Lambda", "Delta1", "Delta2", "pvalue")
write.csv(summary_df, "sensitivity_values_supplement.csv")









