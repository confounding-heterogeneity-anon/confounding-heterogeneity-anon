################################################################################
####       Confounding Heterogeneity in Matched Observational Studies       ####
####                 Reproducible Code for Main Calculations                ####
####                          Finite-Sample Power                           ####
################################################################################

#### README --------------------------------------------------------------------

# running this file produces the following outputs:

# 1) "power_I300_kappa(.5.5 / .25.75 / 01)_tau.1.5.csv": Recreates Table S4
# in the Supplementary Materials.

# 2) "power_(athird / a1 / a3)_kappa(.5.5 / .25.75 / 01).png": Recreates the 
# power figures in the main text and Supplementary Materials.

# all files are placed in the folder "Confounding_Heterogeneity_Finite_Sample_Power".

#### Libraries and functions ---------------------------------------------------

# installs any missing packages, then loads them into the library
required_packages <- c("quadprog", "ggplot2", "latex2exp", "lqmm", "MASS",
                       "stats", "tidyverse", "mvtnorm", "MASS", "VGAM", 
                       "Rsolnp", "nloptr", "sensitivitymult", "doSNOW", 
                       "foreach", "doParallel")
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages) > 0) {
  install.packages(missing_packages)
}
invisible(lapply(required_packages, library, character.only = TRUE))

# creates folder holding all simulation output
dir.create("Confounding_Heterogeneity_Finite_Sample_Power")
setwd("Confounding_Heterogeneity_Finite_Sample_Power")

gen_data_normal <- function(tau, I, rho = 0){
  # tau: vector of effect sizes
  # rho: correlation between potential outcomes under control
  # I: number of matched pairs
  
  Sigma <- matrix(c(0.5, rho, rho, 0.5), nrow = 2)
  Rc <- mvrnorm(n = I, mu = c(0, 0), Sigma = Sigma)
  Rt <- mvrnorm(n = I, mu = tau, Sigma = Sigma)
  Z <- rep(c(1, 0), I)
  
  Y1 <- c(rbind(Rt[, 1], Rc[, 1]))
  Y2 <- c(rbind(Rt[, 2], Rc[, 2]))
  
  output <- data.frame(do.call("cbind", list(Z, Y1, Y2)))
  names(output) <- c("Z", "Y1", "Y2")
  return(output)
  
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
  crit_val <- qchisq(1-2*alpha, 1)
  reject <- ifelse(scaled_test_stat^2 > crit_val, 1, 0)
  
  output <- list("Test_stat" = scaled_test_stat,
                 "reject" = reject, 
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
  K <- ncol(outcomes)
  #Yc <- outcomes[seq(1, nrow(data), by = 2), ] - outcomes[seq(2, nrow(data), by = 2), ] # treated - control
  #Yckappa <- base::t(base::t(Yc) * kappa)
  
  # calculate the test statistic
  observed_test_stat <- observed_test_statistic(data, qfunc)
  observed_cov <- Covariance_CH2(data, rep(1, I), 1, c(1, 1), qfunc)
  scaled_test_stat <- abs(kappa %*% base::t(as.matrix(observed_test_stat))) / sqrt((kappa %*% observed_cov %*% kappa))
  
  # calcualte the crit val
  crit_val <- qchisq(1-2*alpha, 1)
  reject <- ifelse(scaled_test_stat^2 > crit_val, 1, 0)
  
  output <- list("Test_stat" = scaled_test_stat,
                 "reject" = reject, 
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
  K <- ncol(outcomes)
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
  crit_val <- qchisq(1-2*alpha, 1)
  reject <- ifelse(worst_test_stat^2 > crit_val, 1, 0)
  
  output <- list("sensitivity_analysis_stat" = worst_test_stat, 
                 "reject" = reject, 
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
  crit_val <- qchisq(1-2*alpha, 1)
  reject <- ifelse(worst_test_stat^2 > crit_val, 1, 0)
  
  output <- list("sensitivity_analysis_stat" = worst_test_stat, 
                 "reject" = reject, 
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
  crit_val <- qchisq(1-2*alpha, 1)
  reject <- ifelse(worst_test_stat^2 > crit_val, 1, 0)
  
  output <- list("sensitivity_analysis_stat" = worst_test_stat, 
                 "reject" = reject, 
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
  crit_val <- qchisq(1-2*alpha, 1)
  reject <- ifelse(worst_test_stat^2 > crit_val, 1, 0)
  
  output <- list("sensitivity_analysis_stat" = worst_test_stat, 
                 "reject" = reject, 
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
  crit_val <- qchisq(1-2*alpha, 1)
  reject <- ifelse(worst_test_stat^2 > crit_val, 1, 0)
  
  output <- list("sensitivity_analysis_stat" = worst_test_stat, 
                 "reject" = reject, 
                 "worst_case_critval" = crit_val)
  return(output)
}

#### Finite-Sample Power Simulation --------------------------------------------

I <- 300
iter <- 1000
tau <- c(0.1, 0.5)
rho <- 0
set.seed(12345)
data_list <- list()
for(i in 1:iter){
  data_list[[i]] <- gen_data_normal(tau, I, rho)
}
a_grid <- c(1/3, 1, 3)
alpha <- 0.05
qfunc <- "perm"

spacing <- 0.2
Gamma_grid <- seq(1.2, 4, by = spacing)
kappa <- c(0.5, 0.5)

# tau = c(0.1, 0.5), kappa = c(0.5, 0.5)

# setting up progress bar
cl <- makeCluster(10)
registerDoSNOW(cl)
pb <- txtProgressBar(max = iter, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

#registerDoParallel(10)
output <- foreach(i = 1:iter, .combine = rbind,
                  .packages = c("nloptr", "Rsolnp", "mvtnorm", "MASS", "lqmm", "quadprog"),
                  .options.snow = opts) %dopar%
  {
    
    out <- data.frame(matrix(ncol = 5))
    data <- data_list[[i]]
    noconf_classic <- fixed_kappa_noconfounding_classic(data, kappa, alpha, qfunc)
    noconf_amp <- fixed_kappa_noconfounding_amp(data, kappa, alpha, qfunc)
    
    out <- do.call("rbind", list(c(1, 0, "classic", noconf_classic$Test_stat, noconf_classic$reject),
                                 c(1, 0, "amp", noconf_amp$Test_stat, noconf_amp$reject),
                                 c(1, 0, "inf_CH", noconf_amp$Test_stat, noconf_amp$reject),
                                 c(1, 0, "sup_CH", noconf_amp$Test_stat, noconf_amp$reject)))
    
    for(g in Gamma_grid){
      
      if(as.numeric(out[which(out[, 1] == as.character(g-spacing) & out[, 3] == "classic"), 5]) < 0.5){
        out <- rbind(out, c(g, 0, "classic", 0, 0))
      }else{
        classic <- fixed_kappa_classic_ttest(data, g, kappa, alpha, qfunc)
        out <- rbind(out, c(g, 0, "classic", classic$sensitivity_analysis_stat, classic$reject))
      }
      
      
      for(a in a_grid){
        
        pars <- Scaled_amp_parameters(g, a)
        Lambda <- pars[1]
        
        
        if(g < Gamma_grid[3]){
          
          amp <- fixed_kappa_T2CH_ttest(data, Lambda, rep(pars[2], 2), kappa, alpha, qfunc)
          inf_CH <- fixed_kappa_infCH_ttest(data, g, Lambda, kappa, alpha, qfunc)
          sup_CH <- fixed_kappa_supCH_ttest(data, g, Lambda, kappa, alpha, qfunc)
          
          out <- do.call("rbind", list(out, c(g, a, "amp", amp$sensitivity_analysis_stat, amp$reject),
                                       c(g, a, "inf_CH", inf_CH$sensitivity_analysis_stat, inf_CH$reject),
                                       c(g, a, "sup_CH", sup_CH$sensitivity_analysis_stat, sup_CH$reject)))
          
        }
        if(g >= Gamma_grid[3]){
          
          
          if(as.numeric(out[which(out[, 1] == as.character(g-spacing) & out[, 2] == a & out[, 3] == "amp"), 5]) < 0.5){
            out <- rbind(out, c(g, a, "amp", 0, 0))
          }else{
            amp <- fixed_kappa_T2CH_ttest(data, Lambda, rep(pars[2], 2), kappa, alpha, qfunc)
            out <- rbind(out, c(g, a, "amp", amp$sensitivity_analysis_stat, amp$reject))
          }
          
          if(as.numeric(out[which(out[, 1] == as.character(g-spacing) & out[, 2] == a & out[, 3] == "inf_CH"), 5]) < 0.5){
            out <- rbind(out, c(g, a, "inf_CH", 0, 0))
          }else{
            inf_CH <- fixed_kappa_infCH_ttest(data, g, Lambda, kappa, alpha, qfunc)
            out <- rbind(out, c(g, a, "inf_CH", inf_CH$sensitivity_analysis_stat, inf_CH$reject))
          }
          
          if(as.numeric(out[which(out[, 1] == as.character(g-spacing) & out[, 2] == a & out[, 3] == "sup_CH"), 5]) < 0.5){
            out <- rbind(out, c(g, a, "sup_CH", 0, 0))
          }else{
            sup_CH <- fixed_kappa_supCH_ttest(data, g, Lambda, kappa, alpha, qfunc)
            out <- rbind(out, c(g, a, "sup_CH", sup_CH$sensitivity_analysis_stat, sup_CH$reject))
          }
          
          
        }
        
        
        
      }
      
      print(g)
    }
    
    return(out)
  }
close(pb)
stopCluster(cl)
output <- data.frame(output)
output[, c(1:2, 4:5)] <- lapply(output[, c(1:2, 4:5)], as.numeric)
output[, c(1:2, 4:5)] <- round(output[, c(1:2, 4:5)], 3)

summary_df <- data.frame(matrix(nrow = 0, ncol = 4))
Gamma_grid_full <- round(c(1, Gamma_grid), 1)
framework <- c("classic", "amp", "inf_CH", "sup_CH")
for(gidx in 1:length(Gamma_grid_full)){
  for(f in framework){
    g <- Gamma_grid_full[gidx]
    subdata <- output[which(output$X1 == g & output$X3 == f), ]
    unique_a <- unique(subdata$X2)
    
    for(a in unique_a){
      
      subdata_a <- subdata[which(subdata$X2 == a), ]
      num_reject <- mean(subdata_a[, 5] == 1)
      summary_df <- rbind(summary_df, c(g, a, f, num_reject))
      
    }
  }
}
names(summary_df) <- c("Gamma", "a", "Framework", "Num_Reject")
summary_df[, c(1, 2, 4)] <- lapply(summary_df[, c(1, 2, 4)], as.numeric)
summary_df
write.csv(summary_df, "power_I300_kappa.5.5_tau.1.5.csv", row.names = FALSE)


summary_df <- read.csv("power_I300_kappa.5.5_tau.1.5.csv")
summary_df$Framework[which(summary_df$Framework == "classic")] <- "Classic"
summary_df$Framework[which(summary_df$Framework == "amp")] <- "Amplification"
summary_df$Framework[which(summary_df$Framework == "inf_CH")] <- "Inf_CH"
summary_df$Framework[which(summary_df$Framework == "sup_CH")] <- "Sup_CH"

ggplot(summary_df[which(summary_df$a %in% c(0, 0.333)), ]) + 
  geom_point(aes(x = Gamma, y = Num_Reject, shape = Framework, color = Framework), size = 3) +
  geom_line(aes(x = Gamma, y = Num_Reject, color = Framework), linewidth = 0.75) + 
  theme_bw() +
  labs(y = "Power") +
  geom_vline(xintercept = 2.91, color = custom_colors["classic"]) +
  geom_vline(xintercept = 2.09, color = custom_colors["amp"]) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 13),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.background = element_rect(colour = "black"),
        legend.direction = "horizontal",
        legend.position = "inside",
        legend.position.inside = c(1, 1),
        legend.justification = c("right", "top"),
        panel.grid.minor = element_blank()) + 
  guides(shape = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_shape_manual(values = c("Classic" = 1, "Amplification" = 2,
                                "Inf_CH" = 0, "Sup_CH" = 5)) + 
  scale_color_manual(values = custom_colors_slope)
ggsave("power_athird_kappa.5.5.png", width = 6, height = 3)


ggplot(summary_df[which(summary_df$a %in% c(0, 1)), ]) + 
  geom_point(aes(x = Gamma, y = Num_Reject, shape = Framework, color = Framework), size = 3) +
  geom_line(aes(x = Gamma, y = Num_Reject, color = Framework), linewidth = 0.75) + 
  theme_bw() +
  labs(y = "Power") +
  geom_vline(xintercept = 2.91, color = custom_colors["classic"]) +
  geom_vline(xintercept = 2.09, color = custom_colors["amp"]) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 13),
        legend.title = element_blank(),
        legend.key.size = unit(1, "cm"),
        legend.background = element_rect(colour = "black"),
        legend.direction = "horizontal",
        legend.position = "inside",
        legend.position.inside = c(1, 1),
        legend.justification = c("right", "top"),
        panel.grid.minor = element_blank()) + 
  guides(shape = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_shape_manual(values = c("Classic" = 1, "Amplification" = 2,
                                "Inf_CH" = 0, "Sup_CH" = 5)) + 
  scale_color_manual(values = custom_colors_slope)
ggsave("power_a1_kappa.5.5.png", width = 6, height = 3)


ggplot(summary_df[which(summary_df$a %in% c(0, 3)), ]) + 
  geom_point(aes(x = Gamma, y = Num_Reject, shape = Framework, color = Framework), size = 3) +
  geom_line(aes(x = Gamma, y = Num_Reject, color = Framework), linewidth = 0.75) + 
  theme_bw() +
  labs(y = "Power") +
  geom_vline(xintercept = 2.91, color = custom_colors["classic"]) +
  geom_vline(xintercept = 2.09, color = custom_colors["amp"]) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.position = "none",
        panel.grid.minor = element_blank()) + 
  guides(linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_shape_manual(values = c("Classic" = 1, "Amplification" = 2,
                                "Inf_CH" = 0, "Sup_CH" = 5)) + 
  scale_color_manual(values = custom_colors_slope)
ggsave("power_a3_kappa.5.5.png", width = 6, height = 3)



# tau = c(0.1, 0.5), kappa = c(0.25, 0.75)
spacing <- 0.20
Gamma_grid <- seq(1.20, 4, by = spacing)
kappa = c(0.25, 0.75)

# setting up progress bar
cl <- makeCluster(10)
registerDoSNOW(cl)
pb <- txtProgressBar(max = iter, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

#registerDoParallel(10)
output <- foreach(i = 1:iter, .combine = rbind,
                  .packages = c("nloptr", "Rsolnp", "mvtnorm", "MASS", "lqmm", "quadprog"),
                  .options.snow = opts) %dopar%
  {
    
    out <- data.frame(matrix(ncol = 5))
    data <- data_list[[i]]
    noconf_classic <- fixed_kappa_noconfounding_classic(data, kappa, alpha, qfunc)
    noconf_amp <- fixed_kappa_noconfounding_amp(data, kappa, alpha, qfunc)
    
    out <- do.call("rbind", list(c(1, 0, "classic", noconf_classic$Test_stat, noconf_classic$reject),
                                 c(1, 0, "amp", noconf_amp$Test_stat, noconf_amp$reject),
                                 c(1, 0, "inf_CH", noconf_amp$Test_stat, noconf_amp$reject),
                                 c(1, 0, "sup_CH", noconf_amp$Test_stat, noconf_amp$reject)))
    
    for(g in Gamma_grid){
      
      if(as.numeric(out[which(out[, 1] == as.character(g-spacing) & out[, 3] == "classic"), 5]) < 0.5){
        out <- rbind(out, c(g, 0, "classic", 0, 0))
      }else{
        classic <- fixed_kappa_classic_ttest(data, g, kappa, alpha, qfunc)
        out <- rbind(out, c(g, 0, "classic", classic$sensitivity_analysis_stat, classic$reject))
      }
      
      
      for(a in a_grid){
        
        pars <- Scaled_amp_parameters(g, a)
        Lambda <- pars[1]
        
        
        if(g < Gamma_grid[3]){
          
          amp <- fixed_kappa_T2CH_ttest(data, Lambda, rep(pars[2], 2), kappa, alpha, qfunc)
          inf_CH <- fixed_kappa_infCH_ttest(data, g, Lambda, kappa, alpha, qfunc)
          sup_CH <- fixed_kappa_supCH_ttest(data, g, Lambda, kappa, alpha, qfunc)
          
          out <- do.call("rbind", list(out, c(g, a, "amp", amp$sensitivity_analysis_stat, amp$reject),
                                       c(g, a, "inf_CH", inf_CH$sensitivity_analysis_stat, inf_CH$reject),
                                       c(g, a, "sup_CH", sup_CH$sensitivity_analysis_stat, sup_CH$reject)))
          
        }
        if(g >= Gamma_grid[3]){
          
          
          if(as.numeric(out[which(out[, 1] == as.character(g-spacing) & out[, 2] == a & out[, 3] == "amp"), 5]) < 0.5){
            out <- rbind(out, c(g, a, "amp", 0, 0))
          }else{
            amp <- fixed_kappa_T2CH_ttest(data, Lambda, rep(pars[2], 2), kappa, alpha, qfunc)
            out <- rbind(out, c(g, a, "amp", amp$sensitivity_analysis_stat, amp$reject))
          }
          
          if(as.numeric(out[which(out[, 1] == as.character(g-spacing) & out[, 2] == a & out[, 3] == "inf_CH"), 5]) < 0.5){
            out <- rbind(out, c(g, a, "inf_CH", 0, 0))
          }else{
            inf_CH <- fixed_kappa_infCH_ttest(data, g, Lambda, kappa, alpha, qfunc)
            out <- rbind(out, c(g, a, "inf_CH", inf_CH$sensitivity_analysis_stat, inf_CH$reject))
          }
          
          if(as.numeric(out[which(out[, 1] == as.character(g-spacing) & out[, 2] == a & out[, 3] == "sup_CH"), 5]) < 0.5){
            out <- rbind(out, c(g, a, "sup_CH", 0, 0))
          }else{
            sup_CH <- fixed_kappa_supCH_ttest(data, g, Lambda, kappa, alpha, qfunc)
            out <- rbind(out, c(g, a, "sup_CH", sup_CH$sensitivity_analysis_stat, sup_CH$reject))
          }
          
          
        }
        
        
        
      }
      
      print(g)
    }
    
    return(out)
  }
close(pb)
stopCluster(cl)
output <- data.frame(output)
output[, c(1:2, 4:5)] <- lapply(output[, c(1:2, 4:5)], as.numeric)
output[, c(1:2, 4:5)] <- round(output[, c(1:2, 4:5)], 3)

summary_df <- data.frame(matrix(nrow = 0, ncol = 4))
Gamma_grid_full <- round(c(1, Gamma_grid), 1)
framework <- c("classic", "amp", "inf_CH", "sup_CH")
for(gidx in 1:length(Gamma_grid_full)){
  for(f in framework){
    g <- Gamma_grid_full[gidx]
    subdata <- output[which(output$X1 == g & output$X3 == f), ]
    unique_a <- unique(subdata$X2)
    
    for(a in unique_a){
      
      subdata_a <- subdata[which(subdata$X2 == a), ]
      num_reject <- mean(subdata_a[, 5] == 1)
      summary_df <- rbind(summary_df, c(g, a, f, num_reject))
      
    }
  }
}
names(summary_df) <- c("Gamma", "a", "Framework", "Num_Reject")
summary_df
summary_df[, c(1, 2, 4)] <- lapply(summary_df[, c(1, 2, 4)], as.numeric)
write.csv(summary_df, "power_I300_kappa.25.75_tau.1.5.csv", row.names = FALSE)


summary_df <- read.csv("power_I300_kappa.25.75_tau.1.5.csv")
summary_df$Framework[which(summary_df$Framework == "classic")] <- "Classic"
summary_df$Framework[which(summary_df$Framework == "amp")] <- "Amplification"
summary_df$Framework[which(summary_df$Framework == "inf_CH")] <- "Inf_CH"
summary_df$Framework[which(summary_df$Framework == "sup_CH")] <- "Sup_CH"

ggplot(summary_df[which(summary_df$a %in% c(0, 0.333)), ]) + 
  geom_point(aes(x = Gamma, y = Num_Reject, shape = Framework, color = Framework), size = 3) +
  geom_line(aes(x = Gamma, y = Num_Reject, color = Framework), linewidth = 0.75) + 
  theme_bw() +
  labs(y = "Power") +
  geom_vline(xintercept = 3.58, color = custom_colors["classic"]) +
  geom_vline(xintercept = 2.69, color = custom_colors["amp"]) +
  geom_vline(xintercept = 2.11, color = custom_colors["inf_CH"]) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.position = "none",
        panel.grid.minor = element_blank()) + 
  guides(linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_shape_manual(values = c("Classic" = 1, "Amplification" = 2,
                                "Inf_CH" = 0, "Sup_CH" = 5)) + 
  scale_color_manual(values = custom_colors_slope)
ggsave("power_athird_kappa.25.75.png", width = 6, height = 3)


ggplot(summary_df[which(summary_df$a %in% c(0, 1)), ]) + 
  geom_point(aes(x = Gamma, y = Num_Reject, shape = Framework, color = Framework), size = 3) +
  geom_line(aes(x = Gamma, y = Num_Reject, color = Framework), linewidth = 0.75) + 
  theme_bw() +
  labs(y = "Power") +
  geom_vline(xintercept = 3.58, color = custom_colors["classic"]) +
  geom_vline(xintercept = 2.69, color = custom_colors["amp"]) +
  geom_vline(xintercept = 2.11, color = custom_colors["inf_CH"]) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.position = "none",
        panel.grid.minor = element_blank()) + 
  guides(linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_shape_manual(values = c("Classic" = 1, "Amplification" = 2,
                                "Inf_CH" = 0, "Sup_CH" = 5)) + 
  scale_color_manual(values = custom_colors_slope)
ggsave("power_a1_kappa.25.75.png", width = 6, height = 3)


ggplot(summary_df[which(summary_df$a %in% c(0, 3)), ]) + 
  geom_point(aes(x = Gamma, y = Num_Reject, shape = Framework, color = Framework), size = 3) +
  geom_line(aes(x = Gamma, y = Num_Reject, color = Framework), linewidth = 0.75) + 
  theme_bw() +
  labs(y = "Power") +
  geom_vline(xintercept = 3.58, color = custom_colors["classic"]) +
  geom_vline(xintercept = 2.69, color = custom_colors["amp"]) +
  geom_vline(xintercept = 2.11, color = custom_colors["inf_CH"]) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.position = "none",
        panel.grid.minor = element_blank()) + 
  guides(linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_shape_manual(values = c("Classic" = 1, "Amplification" = 2,
                                "Inf_CH" = 0, "Sup_CH" = 5)) + 
  scale_color_manual(values = custom_colors_slope)
ggsave("power_a3_kappa.25.75.png", width = 6, height = 3)

# tau = c(0.1, 0.5), kappa = c(0, 1)
spacing <- 0.20
Gamma_grid <- seq(1.20, 4, by = spacing)
kappa = c(0, 1)

# setting up progress bar
cl <- makeCluster(10)
registerDoSNOW(cl)
pb <- txtProgressBar(max = iter, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

#registerDoParallel(10)
output <- foreach(i = 1:iter, .combine = rbind,
                  .packages = c("nloptr", "Rsolnp", "mvtnorm", "MASS", "lqmm", "quadprog"),
                  .options.snow = opts) %dopar%
  {
    
    out <- data.frame(matrix(ncol = 5))
    data <- data_list[[i]]
    noconf_classic <- fixed_kappa_noconfounding_classic(data, kappa, alpha, qfunc)
    noconf_amp <- fixed_kappa_noconfounding_amp(data, kappa, alpha, qfunc)
    
    out <- do.call("rbind", list(c(1, 0, "classic", noconf_classic$Test_stat, noconf_classic$reject),
                                 c(1, 0, "amp", noconf_amp$Test_stat, noconf_amp$reject),
                                 c(1, 0, "inf_CH", noconf_amp$Test_stat, noconf_amp$reject),
                                 c(1, 0, "sup_CH", noconf_amp$Test_stat, noconf_amp$reject)))
    
    for(g in Gamma_grid){
      
      if(as.numeric(out[which(out[, 1] == as.character(g-spacing) & out[, 3] == "classic"), 5]) < 0.5){
        out <- rbind(out, c(g, 0, "classic", 0, 0))
      }else{
        classic <- fixed_kappa_classic_ttest(data, g, kappa, alpha, qfunc)
        out <- rbind(out, c(g, 0, "classic", classic$sensitivity_analysis_stat, classic$reject))
      }
      
      
      for(a in a_grid){
        
        pars <- Scaled_amp_parameters(g, a)
        Lambda <- pars[1]
        
        
        if(g < Gamma_grid[3]){
          
          amp <- fixed_kappa_T2CH_ttest(data, Lambda, rep(pars[2], 2), kappa, alpha, qfunc)
          inf_CH <- fixed_kappa_infCH_ttest(data, g, Lambda, kappa, alpha, qfunc)
          sup_CH <- fixed_kappa_supCH_unitkappa_ttest(data, g, Lambda, kappa, alpha, qfunc)
          
          out <- do.call("rbind", list(out, c(g, a, "amp", amp$sensitivity_analysis_stat, amp$reject),
                                       c(g, a, "inf_CH", inf_CH$sensitivity_analysis_stat, inf_CH$reject),
                                       c(g, a, "sup_CH", sup_CH$sensitivity_analysis_stat, sup_CH$reject)))
          
        }
        if(g >= Gamma_grid[3]){
          
          
          if(as.numeric(out[which(out[, 1] == as.character(g-spacing) & out[, 2] == a & out[, 3] == "amp"), 5]) < 0.5){
            out <- rbind(out, c(g, a, "amp", 0, 0))
          }else{
            amp <- fixed_kappa_T2CH_ttest(data, Lambda, rep(pars[2], 2), kappa, alpha, qfunc)
            out <- rbind(out, c(g, a, "amp", amp$sensitivity_analysis_stat, amp$reject))
          }
          
          if(as.numeric(out[which(out[, 1] == as.character(g-spacing) & out[, 2] == a & out[, 3] == "inf_CH"), 5]) < 0.5){
            out <- rbind(out, c(g, a, "inf_CH", 0, 0))
          }else{
            inf_CH <- fixed_kappa_infCH_ttest(data, g, Lambda, kappa, alpha, qfunc)
            out <- rbind(out, c(g, a, "inf_CH", inf_CH$sensitivity_analysis_stat, inf_CH$reject))
          }
          
          if(as.numeric(out[which(out[, 1] == as.character(g-spacing) & out[, 2] == a & out[, 3] == "sup_CH"), 5]) < 0.5){
            out <- rbind(out, c(g, a, "sup_CH", 0, 0))
          }else{
            sup_CH <- fixed_kappa_supCH_unitkappa_ttest(data, g, Lambda, kappa, alpha, qfunc)
            out <- rbind(out, c(g, a, "sup_CH", sup_CH$sensitivity_analysis_stat, sup_CH$reject))
          }
          
          
        }
        
        
        
      }
      
      print(g)
    }
    
    return(out)
  }
close(pb)
stopCluster(cl)
output <- data.frame(output)
output[, c(1:2, 4:5)] <- lapply(output[, c(1:2, 4:5)], as.numeric)
output[, c(1:2, 4:5)] <- round(output[, c(1:2, 4:5)], 3)

summary_df <- data.frame(matrix(nrow = 0, ncol = 4))
Gamma_grid_full <- round(c(1, Gamma_grid), 1)
framework <- c("classic", "amp", "inf_CH", "sup_CH")
for(gidx in 1:length(Gamma_grid_full)){
  for(f in framework){
    g <- Gamma_grid_full[gidx]
    subdata <- output[which(output$X1 == g & output$X3 == f), ]
    unique_a <- unique(subdata$X2)
    
    for(a in unique_a){
      
      subdata_a <- subdata[which(subdata$X2 == a), ]
      num_reject <- mean(subdata_a[, 5] == 1)
      summary_df <- rbind(summary_df, c(g, a, f, num_reject))
      
    }
  }
}
names(summary_df) <- c("Gamma", "a", "Framework", "Num_Reject")
summary_df
summary_df[, c(1, 2, 4)] <- lapply(summary_df[, c(1, 2, 4)], as.numeric)
write.csv(summary_df, "power_I300_kappa01_tau.1.5.csv", row.names = FALSE)


summary_df <- read.csv("power_I300_kappa01_tau.1.5.csv")
summary_df$Framework[which(summary_df$Framework == "classic")] <- "Classic"
summary_df$Framework[which(summary_df$Framework == "amp")] <- "Amplification"
summary_df$Framework[which(summary_df$Framework == "inf_CH")] <- "Inf_CH"
summary_df$Framework[which(summary_df$Framework == "sup_CH")] <- "Sup_CH"

ggplot(summary_df[which(summary_df$a %in% c(0, 0.333)), ]) + 
  geom_point(aes(x = Gamma, y = Num_Reject, shape = Framework, color = Framework), size = 3) +
  geom_line(aes(x = Gamma, y = Num_Reject, color = Framework), linewidth = 0.75) + 
  theme_bw() +
  labs(y = "Power") +
  geom_vline(xintercept = 3.53, color = custom_colors["classic"]) +
  geom_vline(xintercept = 1.77, color = custom_colors["inf_CH"]) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.position = "none",
        panel.grid.minor = element_blank()) + 
  guides(linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_shape_manual(values = c("Classic" = 1, "Amplification" = 2,
                                "Inf_CH" = 0, "Sup_CH" = 5)) + 
  scale_color_manual(values = custom_colors_slope)
ggsave("power_athird_kappa01.png", width = 6, height = 3)


ggplot(summary_df[which(summary_df$a %in% c(0, 1)), ]) + 
  geom_point(aes(x = Gamma, y = Num_Reject, shape = Framework, color = Framework), size = 3) +
  geom_line(aes(x = Gamma, y = Num_Reject, color = Framework), linewidth = 0.75) + 
  theme_bw() +
  labs(y = "Power") +
  geom_vline(xintercept = 3.53, color = custom_colors["classic"]) +
  geom_vline(xintercept = 1.77, color = custom_colors["inf_CH"]) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.position = "none",
        panel.grid.minor = element_blank()) + 
  guides(linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_shape_manual(values = c("Classic" = 1, "Amplification" = 2,
                                "Inf_CH" = 0, "Sup_CH" = 5)) + 
  scale_color_manual(values = custom_colors_slope)
ggsave("power_a1_kappa01.png", width = 6, height = 3)


ggplot(summary_df[which(summary_df$a %in% c(0, 3)), ]) + 
  geom_point(aes(x = Gamma, y = Num_Reject, shape = Framework, color = Framework), size = 3) +
  geom_line(aes(x = Gamma, y = Num_Reject, color = Framework), linewidth = 0.75) + 
  theme_bw() +
  labs(y = "Power") +
  geom_vline(xintercept = 3.53, color = custom_colors["classic"]) +
  geom_vline(xintercept = 1.77, color = custom_colors["inf_CH"]) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.position = "none",
        panel.grid.minor = element_blank()) + 
  guides(linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_shape_manual(values = c("Classic" = 1, "Amplification" = 2,
                                "Inf_CH" = 0, "Sup_CH" = 5)) + 
  scale_color_manual(values = custom_colors_slope)
ggsave("power_a3_kappa01.png", width = 6, height = 3)












