################################################################################
####       Confounding Heterogeneity in Matched Observational Studies       ####
####                 Reproducible Code for Main Calculations                ####
####                           Bahadur Efficiency                           ####
################################################################################

#### README --------------------------------------------------------------------

# running this file produces the following outputs:

# 1) "relative_efficiency_(.5.5 / .25.75 / 01).csv": Reproduces 
# the relative Bahadur efficiency results reported in Table S5 of the 
# Supplementary Materials.

# 2) "relative_n_required_(.5.5 / .25.75 / 01).csv": Reproduces 
# the simulated number of matched pairs required to achieve 95% power at 
# significance level 0.01, reported as the parenthesized values of Table S5 in 
# the Supplementary Materials.

# 3) "Bahadur_Slope_Comparison_(athird / a1 / a3)_(.5.5 / .25.75 / 01).png": Reproduces the 
# Bahadur efficiency figures in both the main text and the supplement

# all files are placed in the folder "Confounding_Heterogeneity_Bahadur_Efficiency".

#### Libraries and functions ---------------------------------------------------

# installs any missing packages, then loads them into the library
required_packages <- c("quadprog", "ggplot2", "latex2exp", "lqmm", # "mosaic", 
                       "stats", "tidyverse", "mvtnorm", "MASS", "VGAM", 
                       "Rsolnp", "nloptr", "sensitivitymult", "doSNOW", 
                       "foreach", "doParallel")
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages) > 0) {
  install.packages(missing_packages)
}
invisible(lapply(required_packages, library, character.only = TRUE))

# creates folder holding all simulation output
dir.create("Confounding_Heterogeneity_Bahadur_Efficiency")
setwd("Confounding_Heterogeneity_Bahadur_Efficiency")

# Quality of life
round_to_reference <- function(main, ref){
  #? rounds each entry in main vector to closest entry in reference vector
  # main: vector of doubles to round
  # ref: reference vector
  
  # calculates how different each vector in main is from the entries in ref
  diffmat <- data.frame(matrix(0, nrow = length(ref), ncol = length(main)))
  for(r in 1:length(ref)){
    diffmat[r, ] <- abs(main - ref[r])
  }
  
  # converts each entry in main to its closest match in ref
  output <- rep(0, length(main))
  for(m in 1:length(main)){
    output[m] <- ref[which(diffmat[, m] == min(diffmat[, m]))]
  }
  return(output)
}
erf <- function(x) {
  return(2 * pnorm(x * sqrt(2)) - 1)
}
Scaled_amp_parameters <- function(Gamma, a){
  # a > 0, we are finding \Lambda, \Delta that satisifes \Delta = a\Lambda 
  # and (\Lambda \Delta + 1)/(\Lambda + \Delta) = \Gamma / (1 + \Gamma)
  
  Lambda <- (Gamma*(a+1) + sqrt((Gamma*(a+1))^2 - (4*a)))/(2*a)
  Delta <- a*Lambda
  return(c(Lambda, Delta))
  
}

# phi and eta calculations for qik = |YikVi|, normally distributed
expected_abs_normal <- function(kappa, mu, sigma) {
  # calculates expected value of |X| if X is gaussian
  expected_value <- sigma * sqrt(2/pi) * exp(-mu^2 / (2 * sigma^2)) + mu * (1 - 2 * pnorm(-mu/sigma))
  return(kappa * expected_value)
}
expected_abs_correlated_normals <- function(kappa, muv, sigv, rho){
  # calculates expected value of |kappa_1 x + kappa_2 y| for correlated gaussians X and Y
  mu <- sum(kappa * muv)
  sigma <- sqrt((kappa[1]*sigv[1])^2 + (kappa[2]*sigv[2])^2 + (2*kappa[1]*kappa[2]*rho*sigv[1]*sigv[2]))
  expected_value <- sigma * sqrt(2/pi) * exp(-mu^2 / (2 * sigma^2)) + mu * (1 - 2 * pnorm(-mu/sigma))
  return(expected_value)  
}
mu_kappa <- function(kappa, muv){
  # calculates expected value of T_kappa
  # kappa is a K-length real-valued vector
  # muv is a K-length vector means from K outcomes. 
  return(sum(kappa * muv))
}
phi_kappa_vector_normal <- function(kappa, muv, sigv){
  # calculates expected value of \sum_{k=1}^{K}\kappa_k q_{k}, where q_{k} in our
  # case will be the expected value of |X| for a gaussian X. 
  # kappa is a K-length real-valued vector
  # muv is a K-length vector of the means of the K normally distributed pairwise
  # treatment effects
  # sigv is also a K-length vector, holding the standard deviations for the K normally
  # distributed pairwise treatment effects
  # Deltav is the K-length vector of sensitivity parameters for the K different 
  # outcomes. This is for type 2 confounding heterogeneity. 
  
  exp_absv <- c()
  for(k in 1:length(kappa)){
    exp_absv[k] <- expected_abs_normal(kappa[k], muv[k], sigv[k])
  }
  return(exp_absv)
}
eta_kappa_normal <- function(kappa, muv, sigv, rho = 0){
  eta <- expected_abs_correlated_normals(kappa, muv, sigv, rho)
  return(eta)
}
infimum_sensitivity_CH_normal <- function(kappa, muv, sigv){
  
  # stores \phi_{\kappa, k}
  mu <- mu_kappa(kappa, muv)
  phiv <- phi_kappa_vector_normal(kappa, muv, sigv)
  phi <- sum(phiv)
  
  # asymptote for Lambda
  Lambda_asymptotev <- rep(((phi + mu)/(phi - mu))^(1/length(kappa)), length(kappa))
  
  # worst-case delta asymptote
  phi_ordered <- sort(phiv, decreasing = TRUE)
  nu <- 0
  for(m in 1:(length(phi_ordered)-1)){
    active_set <- seq(1, m)
    nu_m <- 2*(sum(sqrt(phi_ordered[1:m]))/(sum(phi_ordered[1:m]) - mu))^2
    passed <- ifelse(nu_m * phi_ordered[m] > 2 && nu_m * phi_ordered[m+1] <= 2, TRUE, FALSE)
    if(passed){
      nu <- nu_m
    }
  }
  if(nu == 0){
    nu <- 2*(sum(sqrt(phiv))/(sum(phiv) - mu))^2
  }
  Delta_asymptotev <- pmax(1, sqrt(2*phiv*nu) - 1)
  return(c(Lambda_asymptotev, Delta_asymptotev))
}

# bahadur efficiency for qik = |Yik|, normally distributed
bahadur_slope_rosenbaum_normal <- function(kappa, muv, sigv, rho = 0, Gamma, alpha = 0.01, power = 0.95, n = 1e6){
  
  # generate correlated normals
  Sigma <- matrix(
    c(sigv[1]^2,
      rho*sigv[1]*sigv[2],
      rho*sigv[1]*sigv[2],
      sigv[2]^2),
    2,2)
  
  Z <- mvrnorm(
    n = n,
    mu = muv,
    Sigma = Sigma)
  
  # define c_1 function
  omegav <- abs(kappa[1] * Z[, 1] + kappa[2] * Z[, 2])
  c1 <- function(t){
    mean(omegav * (Gamma * exp(omegav * t) - exp(-omegav * t))/(Gamma * exp(omegav * t) + exp(-omegav * t)))
  }
  
  # find the root that solves c1(t) = mu
  mu <- sum(kappa * muv)
  c1_root <- function(t){
    c1(t) - mu
  }
  
  xi_tilde <- uniroot(c1_root, interval = c(0, 100))$root
  
  # calculate bahadur slope
  c0 <- function(t){
    log(1/(1 + Gamma)) + mean(log(Gamma * exp(omegav * t) + exp(-omegav * t)))
  }
  
  slope <- 2 * (xi_tilde * mu - c0(xi_tilde))
  
  # simulated sample size power analysis
  
  F_dist <- kappa[1] * Z[, 1] + kappa[2] * Z[, 2]
  mean_F <- mean(F_dist)
  var_F <- var(F_dist)
  worst_case_mean <- mean((Gamma - 1)/(Gamma + 1)*omegav)
  worst_case_var <- (4*Gamma)/(Gamma + 1)^2*var(omegav)
  
  power_func <- function(n){
    power - pnorm(((mean_F - worst_case_mean) - qnorm(1 - alpha)*sqrt(worst_case_var)/sqrt(n))/(sqrt(var_F)/sqrt(n)))
  }
  min_sample_size <- ceiling(uniroot(power_func, interval = c(1, 100000000000000))$root)
  
  output <- list("slope" = slope, "min_n" = min_sample_size)
  return(output)
}
bahadur_slope_CH_normal <- function(kappa, muv, sigv, Lambdav, Deltav, alpha = 0.01, power = 0.95, n = 1e6){
  # bahadur slope calculation for amplification (prod(Lambdav) = Lambda, Deltav = rep(Delta)),
  # CH2 (prod(Lambav) = Lambda), and CH1.
  
  Z <- matrix(rnorm(2*n, mean = rep(muv, each = n), sd = rep(sigv, each = n)),
              nrow = n, ncol = 2)
  
  omegav <- matrix(cbind(abs(kappa[1] * Z[, 1]), abs(kappa[2] * Z[, 2])), ncol = 2)
  prod_Lambda <- prod(Lambdav)
  pv <- ((prod_Lambda * Deltav + 1))/((prod_Lambda + 1)*(Deltav + 1))
  
  c1 <- function(t){
    mean(omegav[, 1] * ((pv[1] * exp(omegav[, 1] * t) - (1 - pv[1])*exp(-omegav[, 1] * t))/(pv[1] * exp(omegav[, 1] * t) + (1 - pv[1])*exp(-omegav[, 1] * t)))) + mean(omegav[, 2] * ((pv[2] * exp(omegav[, 2] * t) - (1 - pv[2])*exp(-omegav[, 2] * t))/(pv[2] * exp(omegav[, 2] * t) + (1 - pv[2])*exp(-omegav[, 2] * t))))   
  }
  
  # find the root that solves c1(t) = mu
  mu <- sum(kappa * muv)
  c1_root <- function(t){
    c1(t) - mu
  }
  
  xi_tilde <- uniroot(c1_root, interval = c(0, 100))$root
  
  # calculate bahadur slope
  c0 <- function(t){
    mean(log(pv[1]*exp(omegav[, 1] * t) + (1-pv[1])*exp(-omegav[, 1] * t))) + mean(log(pv[2]*exp(omegav[, 2] * t) + (1-pv[2])*exp(-omegav[, 2] * t)))
  }
  
  slope <- 2 * (xi_tilde * mu - c0(xi_tilde))
  
  F_dist <- kappa[1] * Z[, 1] + kappa[2] * Z[, 2]
  mean_F <- mean(F_dist)
  var_F <- var(F_dist)
  worst_case_mean <- mean((2*pv[1] - 1)*omegav[, 1] + (2*pv[2] - 1)*omegav[, 2])
  worst_case_var <- 4*pv[1]*(1-pv[1])*var(omegav[, 1]) + 4*pv[2]*(1-pv[2])*var(omegav[, 2])
  power_func <- function(n){
    power - pnorm(((mean_F - worst_case_mean) - qnorm(1 - alpha)*sqrt(worst_case_var)/sqrt(n))/(sqrt(var_F)/sqrt(n)))
  }
  min_sample_size <- ceiling(uniroot(power_func, interval = c(1, 100000000000000))$root)
  
  output <- list("slope" = slope, "min_n" = min_sample_size)
  return(output)
  
}
bahadur_slope_CH_proj <- function(kappa, muv, sigv, Gamma, a, alpha = 0.01, power = 0.95, n = 1e6){
  
  Z <- matrix(rnorm(2*n, mean = rep(muv, each = n), sd = rep(sigv, each = n)),
              nrow = n, ncol = 2)
  
  omegav <- matrix(cbind(abs(kappa[1] * Z[, 1]), abs(kappa[2] * Z[, 2])), ncol = 2)
  F_dist <- kappa[1] * Z[, 1] + kappa[2] * Z[, 2]
  mean_F <- mean(F_dist)
  var_F <- var(F_dist)
  
  # lambda and delta parameters
  lamdel <- Scaled_amp_parameters(Gamma, a)
  
  #### first, the bahadur slope for amplification framework
  Lambda <- lamdel[1]
  Delta <- lamdel[2]
  p_amp <- (Lambda * Delta + 1)/((Lambda + 1)*(Delta + 1))
  K <- length(kappa)
  pv <- rep(p_amp, K)
  mean_vector <- kappa * muv
  mu <- sum(mean_vector)
  
  
  m0v <- function(mat, t){
    rowSums(log((Delta*exp(mat * t) + exp(-mat * t))/(Delta + 1)))
  }
  
  Rv <- function(mat, t){
    log((Lambda + exp(m0v(mat, -t) - m0v(mat, t)))/(Lambda + 1))
  }
  
  m1v <- function(mat, t){
    rowSums(mat*((Delta * exp(mat * t) - exp(-mat * t))/(Delta * exp(mat * t) + exp(-mat * t))))
  }
  
  Qv <- function(mat, t){
    exp(m0v(mat, -t) - m0v(mat, t))/(Lambda + exp(m0v(mat, -t) - m0v(mat, t)))
  }
  
  c1 <- function(t){
    mean(m1v(omegav, t) - Qv(omegav, t)*(m1v(omegav, -t) + m1v(omegav, t)))
  }
  
  # find the root that solves c1(t) = mu
  c1_root <- function(t){
    c1(t) - mu
  }
  
  xi_tilde <- uniroot(c1_root, interval = c(0, 100))$root
  
  # calculate bahadur slope
  c0 <- function(t){
    mean(m0v(omegav, t) + Rv(omegav, t))
  }
  
  amp_slope <- 2 * (xi_tilde * mu - c0(xi_tilde))
  worst_case_mean <- mean((2*pv[1] - 1)*omegav[, 1] + (2*pv[2] - 1)*omegav[, 2])
  off_diag <- ((Delta - 1)/(Delta + 1))^2 * (1 - ((Lambda - 1)/(Lambda + 1))^2)
  worst_case_var <- 4*pv[1]*(1-pv[1])*var(omegav[, 1]) + 4*pv[2]*(1-pv[2])*var(omegav[, 2]) + 2*mean(omegav[, 1] * omegav[, 2])*off_diag
  power_func <- function(n){
    power - pnorm(((mean_F - worst_case_mean) - qnorm(1 - alpha)*sqrt(worst_case_var)/sqrt(n))/(sqrt(var_F)/sqrt(n)))
  }
  min_sample_size_amp <- ceiling(uniroot(power_func, interval = c(1, 100000000000000))$root)
  
  
  xL <- (Lambda - 1)/(Lambda + 1)
  obj_func <- function(x1){
    
    clip_omega <- omegav[1:1e4, ]
    
    if(xL != 0){
      x2 <- ((K * (Gamma - 1)/(Gamma + 1)) - xL*x1)/xL
    }
    if(xL == 0){
      x2 <- ((K * (Gamma - 1)/(Gamma + 1)) - 0.0001*x1)/0.0001
    }
    
    Delta_vec <- c((x1 + 1)/(1 - x1), (x2 + 1)/(1 - x2))
    Delta_vec[which(Delta_vec == Inf)] <- 1000000
    
    m0v <- function(mat, t){
      num <- (t(t(exp(mat * t)) * Delta_vec) + exp(-mat * t))
      t(t(num) * (1/(Delta_vec + 1)))
      return(rowSums(log(t(t(num) * (1/(Delta_vec + 1))))))
    }
    
    Rv <- function(mat, t){
      log((Lambda + exp(m0v(mat, -t) - m0v(mat, t)))/(Lambda + 1))
    }
    
    m1v <- function(mat, t){
      num <- t(t(exp(mat * t)) * Delta_vec) - exp(-mat * t)
      denom <- t(t(exp(mat * t)) * Delta_vec) + exp(-mat * t)
      return(rowSums((mat * num)/denom))
    }
    
    Qv <- function(mat, t){
      exp(m0v(mat, -t) - m0v(mat, t))/(Lambda + exp(m0v(mat, -t) - m0v(mat, t)))
    }
    
    c1 <- function(t){
      mean(m1v(clip_omega, t) - Qv(clip_omega, t)*(m1v(clip_omega, -t) + m1v(clip_omega, t)))
    }
    
    # find the root that solves c1(t) = mu
    c1_root <- function(t){
      c1(t) - mu
    }
    
    xi_tilde <- uniroot(c1_root, interval = c(0, 100))$root
    
    # calculate bahadur slope
    c0 <- function(t){
      mean(m0v(clip_omega, t) + Rv(clip_omega, t))
    }
    
    slope <- 2 * (xi_tilde * mu - c0(xi_tilde))
    return(slope)
  }
  
  if(Lambda != 1){
    lb <- max(0, ((K * (Gamma - 1)/(Gamma + 1)) - xL)/xL)
    ub <- min(1, (K * (Gamma - 1)/(Gamma + 1))/xL)
  }
  if(Lambda == 1){
    lb <- 0
    ub <- 1
  }
  opts <- list("algorithm" = "NLOPT_LN_BOBYQA", "xtol_rel" = 1.0e-5)
  inf_sol <- nloptr(x0 = (lb+ub)/2, eval_f = obj_func, opts = opts, lb = lb, ub = ub)
  edel <- c(inf_sol$solution, ((K * (Gamma - 1)/(Gamma + 1)) - xL*inf_sol$solution)/xL)
  if(inf_sol$solution == 0 & Gamma == 1){
    edel <- c(0, 1)
  }
  Delta_vec <- (1+edel)/(1-edel)
  Delta_vec[which(Delta_vec > 1e6)] <- 1e6
  
  m0v <- function(mat, t){
    num <- (t(t(exp(mat * t)) * Delta_vec) + exp(-mat * t))
    t(t(num) * (1/(Delta_vec + 1)))
    return(rowSums(log(t(t(num) * (1/(Delta_vec + 1))))))
  }
  
  Rv <- function(mat, t){
    log((Lambda + exp(m0v(mat, -t) - m0v(mat, t)))/(Lambda + 1))
  }
  
  m1v <- function(mat, t){
    num <- t(t(exp(mat * t)) * Delta_vec) - exp(-mat * t)
    denom <- t(t(exp(mat * t)) * Delta_vec) + exp(-mat * t)
    return(rowSums((mat * num)/denom))
  }
  
  Qv <- function(mat, t){
    exp(m0v(mat, -t) - m0v(mat, t))/(Lambda + exp(m0v(mat, -t) - m0v(mat, t)))
  }
  
  c1 <- function(t){
    mean(m1v(omegav, t) - Qv(omegav, t)*(m1v(omegav, -t) + m1v(omegav, t)))
  }
  
  # find the root that solves c1(t) = mu
  c1_root <- function(t){
    c1(t) - mu
  }
  
  xi_tilde <- uniroot(c1_root, interval = c(0, 100))$root
  
  # calculate bahadur slope
  c0 <- function(t){
    mean(m0v(omegav, t) + Rv(omegav, t))
  }
  
  inf_slope <- 2 * (xi_tilde * mu - c0(xi_tilde))
  pv <- c((Lambda * Delta_vec[1] + 1)/((Lambda + 1)*(Delta_vec[1] + 1)), (Lambda * Delta_vec[2] + 1)/((Lambda + 1)*(Delta_vec[2] + 1)))
  worst_case_mean <- mean((2*pv[1] - 1)*omegav[, 1] + (2*pv[2] - 1)*omegav[, 2])
  off_diag <- ((Delta_vec[1] - 1)/(Delta_vec[1] + 1))*((Delta_vec[2] - 1)/(Delta_vec[2] + 1)) * (1 - ((Lambda - 1)/(Lambda + 1))^2)
  worst_case_var <- 4*pv[1]*(1-pv[1])*var(omegav[, 1]) + 4*pv[2]*(1-pv[2])*var(omegav[, 2]) + 2*mean(omegav[, 1] * omegav[, 2])*off_diag
  power_func <- function(n){
    power - pnorm(((mean_F - worst_case_mean) - qnorm(1 - alpha)*sqrt(worst_case_var)/sqrt(n))/(sqrt(var_F)/sqrt(n)))
  }
  min_sample_size_inf <- ceiling(uniroot(power_func, interval = c(1, 100000000000000))$root)
  
  lower <- obj_func(lb)
  upper <- obj_func(ub)
  solution <- lb
  if(upper > lower){
    solution <- ub
  }
  
  edel <- c(solution, ((K * (Gamma - 1)/(Gamma + 1)) - xL*solution)/xL)
  Delta_vec <- (1+edel)/(1-edel)
  Delta_vec[which(Delta_vec == Inf)] <- 1e6
  
  m0v <- function(mat, t){
    num <- (t(t(exp(mat * t)) * Delta_vec) + exp(-mat * t))
    t(t(num) * (1/(Delta_vec + 1)))
    return(rowSums(log(t(t(num) * (1/(Delta_vec + 1))))))
  }
  
  Rv <- function(mat, t){
    log((Lambda + exp(m0v(mat, -t) - m0v(mat, t)))/(Lambda + 1))
  }
  
  m1v <- function(mat, t){
    num <- t(t(exp(mat * t)) * Delta_vec) - exp(-mat * t)
    denom <- t(t(exp(mat * t)) * Delta_vec) + exp(-mat * t)
    return(rowSums((mat * num)/denom))
  }
  
  Qv <- function(mat, t){
    exp(m0v(mat, -t) - m0v(mat, t))/(Lambda + exp(m0v(mat, -t) - m0v(mat, t)))
  }
  
  c1 <- function(t){
    mean(m1v(omegav, t) - Qv(omegav, t)*(m1v(omegav, -t) + m1v(omegav, t)))
  }
  
  # find the root that solves c1(t) = mu
  c1_root <- function(t){
    c1(t) - mu
  }
  
  xi_tilde <- uniroot(c1_root, interval = c(0, 100))$root
  
  # calculate bahadur slope
  c0 <- function(t){
    mean(m0v(omegav, t) + Rv(omegav, t))
  }
  
  sup_slope <- 2 * (xi_tilde * mu - c0(xi_tilde))
  pv <- c((Lambda * Delta_vec[1] + 1)/((Lambda + 1)*(Delta_vec[1] + 1)), (Lambda * Delta_vec[2] + 1)/((Lambda + 1)*(Delta_vec[2] + 1)))
  worst_case_mean <- mean((2*pv[1] - 1)*omegav[, 1] + (2*pv[2] - 1)*omegav[, 2])
  off_diag <- ((Delta_vec[1] - 1)/(Delta_vec[1] + 1))*((Delta_vec[2] - 1)/(Delta_vec[2] + 1)) * (1 - ((Lambda - 1)/(Lambda + 1))^2)
  worst_case_var <- 4*pv[1]*(1-pv[1])*var(omegav[, 1]) + 4*pv[2]*(1-pv[2])*var(omegav[, 2]) + 2*mean(omegav[, 1] * omegav[, 2])*off_diag
  power_func <- function(n){
    power - pnorm(((mean_F - worst_case_mean) - qnorm(1 - alpha)*sqrt(worst_case_var)/sqrt(n))/(sqrt(var_F)/sqrt(n)))
  }
  min_sample_size_sup <- ceiling(uniroot(power_func, interval = c(1, 100000000000000))$root)
  
  output <- list("slopes" = c(inf_slope, amp_slope, sup_slope),
                 "min_n" = c(min_sample_size_inf, min_sample_size_amp, min_sample_size_sup))
  
  return(output)
  
}
bahadur_slope_CH_initial <- function(kappa, muv, sigv, a, alpha = 0.01, power = 0.95, n = 1e6){
  # bahadur slope calculation for amplification (prod(Lambdav) = Lambda, Deltav = rep(Delta)),
Z <- matrix(rnorm(2*n, mean = rep(muv, each = n), sd = rep(sigv, each = n)),
              nrow = n, ncol = 2)
  
  omegav <- matrix(cbind(abs(kappa[1] * Z[, 1]), abs(kappa[2] * Z[, 2])), ncol = 2)
  F_dist <- kappa[1] * Z[, 1] + kappa[2] * Z[, 2]
  mean_F <- mean(F_dist)
  var_F <- var(F_dist)
  
  #### first, the bahadur slope for amplification framework
  
  lamdel <- Scaled_amp_parameters(Gamma, a)
  Lambda <- lamdel[1]
  Delta <- lamdel[2]
  p_amp <- (Lambda * Delta + 1)/((Lambda + 1)*(Delta + 1))
  K <- length(kappa)
  pv <- rep(p_amp, K)
  mean_vector <- kappa * muv
  mu <- sum(mean_vector)
  
  
  m0v <- function(mat, t){
    rowSums(log((Delta*exp(mat * t) + exp(-mat * t))/(Delta + 1)))
  }
  
  Rv <- function(mat, t){
    log((Lambda + exp(m0v(mat, -t) - m0v(mat, t)))/(Lambda + 1))
  }
  
  m1v <- function(mat, t){
    rowSums(mat*((Delta * exp(mat * t) - exp(-mat * t))/(Delta * exp(mat * t) + exp(-mat * t))))
  }
  
  Qv <- function(mat, t){
    exp(m0v(mat, -t) - m0v(mat, t))/(Lambda + exp(m0v(mat, -t) - m0v(mat, t)))
  }
  
  c1 <- function(t){
    mean(m1v(omegav, t) - Qv(omegav, t)*(m1v(omegav, -t) + m1v(omegav, t)))
  }
  
  # find the root that solves c1(t) = mu
  c1_root <- function(t){
    c1(t) - mu
  }
  
  xi_tilde <- uniroot(c1_root, interval = c(0, 100))$root
  
  # calculate bahadur slope
  c0 <- function(t){
    mean(m0v(omegav, t) + Rv(omegav, t))
  }
  
  amp_slope <- 2 * (xi_tilde * mu - c0(xi_tilde))
  worst_case_mean <- mean((2*pv[1] - 1)*omegav[, 1] + (2*pv[2] - 1)*omegav[, 2])
  off_diag <- ((Delta - 1)/(Delta + 1))^2 * (1 - ((Lambda - 1)/(Lambda + 1))^2)
  worst_case_var <- 4*pv[1]*(1-pv[1])*var(omegav[, 1]) + 4*pv[2]*(1-pv[2])*var(omegav[, 2]) + 2*mean(omegav[, 1] * omegav[, 2])*off_diag
  power_func <- function(n){
    power - pnorm(((mean_F - worst_case_mean) - qnorm(1 - alpha)*sqrt(worst_case_var)/sqrt(n))/(sqrt(var_F)/sqrt(n)))
  }
  min_sample_size <- ceiling(uniroot(power_func, interval = c(1, 100000000000000))$root)
  
  output <- list("slope" = amp_slope, "min_n" = min_sample_size)
  return(output)
  
}
slope_grid_amp_normal <- function(kappa, muv, sigv, LambdaMax, refine){
  
  # LambdaMax = largest value of lambda I want to investigate
  # refine = grid spacing
  
  mu <- mu_kappa(kappa, muv)
  phi <- sum(phi_kappa_vector_normal(kappa, muv, sigv))
  asymptote <- (phi + mu)/(phi - mu)
  
  amp_curve <- function(Lambda){
    rhs <- ((Lambda + 1)/(Lambda - 1))*(mu / phi)
    return((rhs + 1)/(1 - rhs))
  }
  
  Lambda_grid <- seq(1.1, LambdaMax, by = refine)
  Lam_Del_grid <- data.frame(matrix(cbind(Lambda_grid[1], Lambda_grid), ncol = 2))
  for(Lambda in Lambda_grid[-1]){
    
    if(Lambda < asymptote){
      Lam_Del_grid <- rbind(Lam_Del_grid, data.frame(matrix(cbind(Lambda, Lambda_grid), ncol = 2)))
    }
    
    if(Lambda > asymptote){
      DeltaMax <- min(LambdaMax, amp_curve(Lambda))
      Delta_grid <- seq(1.1, DeltaMax, by = refine)
      Lam_Del_grid <- rbind(Lam_Del_grid, data.frame(matrix(cbind(Lambda, Delta_grid), ncol = 2)))
    }
    
  }
  names(Lam_Del_grid) <- c("Lambda", "Delta")
  return(Lam_Del_grid)
  
}
CH2_delta_grid <- function(Gamma, Lambda, length_out){
  # K = 2
  
  min_LHS <- ((Gamma / (Gamma + 1))*2 - (Lambda / (Lambda + 1))) * (Lambda + 1)
  min_Delta <- max((1 - min_LHS)/(min_LHS - Lambda), 1) + 1
  
  LHS <- ((Gamma / (Gamma + 1))*2 - ((Lambda * min_Delta + 1)/((Lambda + 1)*(min_Delta + 1)))) * (Lambda + 1)
  max_Delta <- (1 - LHS) / (LHS - Lambda)
  
  grid <- seq(min_Delta, max_Delta, length.out = length_out)
  Delta2 <- c()
  for(Delta in grid){
    LHS <- ((Gamma / (Gamma + 1))*2 - ((Lambda * Delta + 1)/((Lambda + 1)*(Delta + 1)))) * (Lambda + 1)
    Delta2 <- c(Delta2, (1 - LHS) / (LHS - Lambda))
  }
  
  return(matrix(cbind(grid, Delta2), ncol = 2))
  
}

#### Bahadur Slope, kappa = 0.5, 0.5 --------------------

# kappa = c(0.5, 0.5)
set.seed(12345)
kappa <- c(0.5, 0.5)
muv <- c(0.1, 0.5)
sigv <- c(1, 1)
#inf_gamma <- round(mean(infimum_sensitivity_CH_normal(kappa, muv, sigv)[(length(kappa) + 1):(2*length(kappa))]), 2)
gamma_grid <- c(1.1, 1.3, 1.5, 1.7)
a_grid <- c(1/3, 1, 3)

# initializer

slope_mat <- matrix(ncol = 6, nrow = 0)
n_mat <- matrix(ncol = 6, nrow = 0)

init_slope_rosenbaum <- bahadur_slope_rosenbaum_normal(kappa, muv, sigv, rho = 0, Gamma = 1)
for(a in a_grid){
  
  if(a <= 1){
    slope_CH <- bahadur_slope_CH_initial(kappa, muv, sigv, a)
    slope_mat <- rbind(slope_mat, c(1, a, init_slope_rosenbaum$slope, rep(slope_CH$slope, 3)))
    n_mat <- rbind(n_mat, c(1, a, init_slope_rosenbaum$min_n, rep(slope_CH$min_n, 3))) 
  }
  if(a > 1){
    slope_CH <- bahadur_slope_CH_initial(kappa, muv, sigv, a)
    slope_mat <- rbind(slope_mat, c(1, a, init_slope_rosenbaum$slope, rep(slope_CH$slope, 2), init_slope_rosenbaum$slope))
    n_mat <- rbind(n_mat, c(1, a, init_slope_rosenbaum$min_n, rep(slope_CH$min_n, 2), init_slope_rosenbaum$min_n)) 
  }
  
}


for(g in gamma_grid){
  
  slope_rosenbaum <- bahadur_slope_rosenbaum_normal(kappa, muv, sigv, rho = 0, g)
  
  for(a in a_grid){
    
    slope_CH <- bahadur_slope_CH_proj(kappa, muv, sigv, g, a)
    slope_mat <- rbind(slope_mat, c(g, a, slope_rosenbaum$slope, slope_CH$slope))
    n_mat <- rbind(n_mat, c(g, a, slope_rosenbaum$min_n, slope_CH$min_n))
    print(a)
  }
  print(g)
}

slope_mat <- data.frame(slope_mat)
names(slope_mat) <- c("Gamma", "a", "classic", "infCH", "amp", "supCH")
n_mat <- data.frame(n_mat)
names(n_mat) <- c("Gamma", "a", "classic", "infCH", "amp", "supCH")

rel_eff_mat <- round(do.call("cbind", list(slope_mat[, 1:2], slope_mat[, 3:6] / slope_mat[, 3])), 2)
names(rel_eff_mat) <- c("Gamma", "a", "Classic", "Inf_Slope_CH", "Amplification", "Sup_Slope_CH")
rel_n_mat <- round(do.call("cbind", list(n_mat[, 1:2], n_mat[, 3] / n_mat[, 3:6])), 2)
names(rel_n_mat) <- c("Gamma", "a", "Classic", "Inf_Slope_CH", "Amplification", "Sup_Slope_CH")
rel_eff_mat

write.csv(rel_eff_mat, "relative_efficiency_.5.5.csv")
write.csv(rel_n_mat, "relative_n_required_.5.5.csv")
write.csv(slope_mat, "slope_mat_.5.5.csv")

# visual of the slopes
slope_mat <- read.csv("slope_mat_.5.5.csv")[, -1]
slopes_df <- round(data.frame(slope_mat), 3)
names(slopes_df) <- c("Gamma", "a", "Classic", "Inf_CH", "Amplification", "Sup_CH")
slopes_df_long <- slopes_df %>%
  pivot_longer(cols = Classic:Sup_CH, 
               names_to = "Framework",
               values_to = "Slope")

# a = 1/3
slopes_df_long_third <- slopes_df_long[which(slopes_df_long$a < 0.5), ]

ggplot(slopes_df_long_third) + 
  geom_point(aes(x = Gamma, y = Slope, shape = Framework, color = Framework), size = 3) +
  geom_line(aes(x = Gamma, y = Slope, color = Framework), linewidth = 0.75) + 
  theme_bw() + 
  labs(y = "Bahadur Slope") + 
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
  coord_cartesian(xlim = c(1, 1.7), ylim = c(0, 0.25)) +
  scale_x_continuous(breaks = c(1, 1.1, 1.3, 1.5, 1.7)) + 
  scale_y_continuous(breaks = seq(0, 0.25, by = 0.05)) + 
  scale_shape_manual(values = c("Classic" = 1, "Amplification" = 2,
                                "Inf_CH" = 0, "Sup_CH" = 5)) + 
  scale_color_manual(values = custom_colors_slope)
ggsave("Bahadur_slope_comparison_athird_.5.5.png", width = 6, height = 3)

# a = 1
slopes_df_long_1 <- slopes_df_long[which(slopes_df_long$a == 0 | slopes_df_long$a == 1), ]

ggplot(slopes_df_long_1) + 
  geom_point(aes(x = Gamma, y = Slope, shape = Framework, color = Framework), size = 3) +
  geom_line(aes(x = Gamma, y = Slope, color = Framework), linewidth = 0.75) + 
  theme_bw() + 
  labs(y = "Bahadur Slope") + 
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
  coord_cartesian(xlim = c(1, 1.7), ylim = c(0, 0.25)) +
  scale_x_continuous(breaks = c(1, 1.1, 1.3, 1.5, 1.7)) + 
  scale_y_continuous(breaks = seq(0, 0.25, by = 0.05)) + 
  scale_shape_manual(values = c("Classic" = 1, "Amplification" = 2,
                                "Inf_CH" = 0, "Sup_CH" = 5)) + 
  scale_color_manual(values = custom_colors_slope)
ggsave("Bahadur_slope_comparison_a1_.5.5.png", width = 6, height = 3)

# a = 3
slopes_df_long_3 <- slopes_df_long[which(slopes_df_long$a == 0 | slopes_df_long$a == 3), ]

ggplot(slopes_df_long_3) + 
  geom_point(aes(x = Gamma, y = Slope, shape = Framework, color = Framework), size = 3) +
  geom_line(aes(x = Gamma, y = Slope, color = Framework), linewidth = 0.75) + 
  theme_bw() + 
  labs(y = "Bahadur Slope") + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.position = "blank",
        panel.grid.minor = element_blank()) + 
  guides(shape = guide_legend(nrow = 2, byrow = TRUE)) + 
  coord_cartesian(xlim = c(1, 1.7), ylim = c(0, 0.25)) +
  scale_x_continuous(breaks = c(1, 1.1, 1.3, 1.5, 1.7)) + 
  scale_y_continuous(breaks = seq(0, 0.25, by = 0.05)) + 
  scale_shape_manual(values = c("Classic" = 1, "Amplification" = 2,
                                "Inf_CH" = 0, "Sup_CH" = 5)) + 
  scale_color_manual(values = custom_colors_slope)
ggsave("Bahadur_slope_comparison_a3_.5.5.png", width = 6, height = 3)

#### Bahadur Slope, kappa = 0.25, 0.75 --------------------

# kappa = c(0.25, 0.75)
set.seed(12345)
kappa <- c(0.25, 0.75)
muv <- c(0.1, 0.5)
sigv <- c(1, 1)
#inf_gamma <- round(mean(infimum_sensitivity_CH_normal(kappa, muv, sigv)[(length(kappa) + 1):(2*length(kappa))]), 2)
gamma_grid <- c(1.1, 1.3, 1.5, 1.7)
a_grid <- c(1/3, 1, 3)

# initializer

slope_mat <- matrix(ncol = 6, nrow = 0)
n_mat <- matrix(ncol = 6, nrow = 0)

init_slope_rosenbaum <- bahadur_slope_rosenbaum_normal(kappa, muv, sigv, rho = 0, Gamma = 1)
for(a in a_grid){
  
  if(a <= 1){
    slope_CH <- bahadur_slope_CH_initial(kappa, muv, sigv, a)
    slope_mat <- rbind(slope_mat, c(1, a, init_slope_rosenbaum$slope, rep(slope_CH$slope, 3)))
    n_mat <- rbind(n_mat, c(1, a, init_slope_rosenbaum$min_n, rep(slope_CH$min_n, 3))) 
  }
  if(a > 1){
    slope_CH <- bahadur_slope_CH_initial(kappa, muv, sigv, a)
    slope_mat <- rbind(slope_mat, c(1, a, init_slope_rosenbaum$slope, rep(slope_CH$slope, 2), init_slope_rosenbaum$slope))
    n_mat <- rbind(n_mat, c(1, a, init_slope_rosenbaum$min_n, rep(slope_CH$min_n, 2), init_slope_rosenbaum$min_n)) 
  }
  
}

for(g in gamma_grid){
  
  for(a in a_grid){
    
    slope_rosenbaum <- bahadur_slope_rosenbaum_normal(kappa, muv, sigv, rho = 0, g)
    slope_CH <- bahadur_slope_CH_proj(kappa, muv, sigv, g, a)
    slope_mat <- rbind(slope_mat, c(g, a, slope_rosenbaum$slope, slope_CH$slope))
    n_mat <- rbind(n_mat, c(g, a, slope_rosenbaum$min_n, slope_CH$min_n))
    print(a)
  }
  print(g)
}

slope_mat <- data.frame(slope_mat)
names(slope_mat) <- c("Gamma", "a", "classic", "infCH", "amp", "supCH")
n_mat <- data.frame(n_mat)
names(n_mat) <- c("Gamma", "a", "classic", "infCH", "amp", "supCH")

rel_eff_mat <- round(do.call("cbind", list(slope_mat[, 1:2], slope_mat[, 3:6] / slope_mat[, 3])), 2)
names(rel_eff_mat) <- c("Gamma", "a", "Classic", "Inf_Slope_CH", "Amplification", "Sup_Slope_CH")
rel_n_mat <- round(do.call("cbind", list(n_mat[, 1:2], n_mat[, 3] / n_mat[, 3:6])), 2)
names(rel_n_mat) <- c("Gamma", "a", "Classic", "Inf_Slope_CH", "Amplification", "Sup_Slope_CH")
rel_eff_mat

write.csv(rel_eff_mat, "relative_efficiency_.25.75.csv")
write.csv(rel_n_mat, "relative_n_required_.25.75.csv")
write.csv(slope_mat, "slope_mat_.25.75.csv")

# visual of the slopes
slope_mat <- read.csv("slope_mat_.25.75.csv")[, -1]
slopes_df <- round(data.frame(slope_mat), 3)
names(slopes_df) <- c("Gamma", "a", "Classic", "Inf_CH", "Amplification", "Sup_CH")
slopes_df_long <- slopes_df %>%
  pivot_longer(cols = Classic:Sup_CH, 
               names_to = "Framework",
               values_to = "Slope")

# a = 1/3
slopes_df_long_third <- slopes_df_long[which(slopes_df_long$a < 0.5), ]

ggplot(slopes_df_long_third) + 
  geom_point(aes(x = Gamma, y = Slope, shape = Framework, color = Framework), size = 3) +
  geom_line(aes(x = Gamma, y = Slope, color = Framework), linewidth = 0.75) + 
  theme_bw() + 
  labs(y = "Bahadur Slope") + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.position = "blank",
        panel.grid.minor = element_blank()) + 
  guides(shape = guide_legend(nrow = 2, byrow = TRUE)) + 
  coord_cartesian(xlim = c(1, 1.7), ylim = c(0, 0.25)) +
  scale_x_continuous(breaks = c(1, 1.1, 1.3, 1.5, 1.7)) + 
  scale_y_continuous(breaks = seq(0, 0.25, by = 0.05)) + 
  scale_shape_manual(values = c("Classic" = 1, "Amplification" = 2,
                                "Inf_CH" = 0, "Sup_CH" = 5)) + 
  scale_color_manual(values = custom_colors_slope)
ggsave("Bahadur_slope_comparison_athird_.25.75.png", width = 6, height = 3)

# a = 1
slopes_df_long_1 <- slopes_df_long[which(slopes_df_long$a == 0 | slopes_df_long$a == 1), ]

ggplot(slopes_df_long_1) + 
  geom_point(aes(x = Gamma, y = Slope, shape = Framework, color = Framework), size = 3) +
  geom_line(aes(x = Gamma, y = Slope, color = Framework), linewidth = 0.75) + 
  theme_bw() + 
  labs(y = "Bahadur Slope") + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.position = "blank",
        panel.grid.minor = element_blank()) + 
  guides(shape = guide_legend(nrow = 2, byrow = TRUE)) + 
  coord_cartesian(xlim = c(1, 1.7), ylim = c(0, 0.25)) +
  scale_x_continuous(breaks = c(1, 1.1, 1.3, 1.5, 1.7)) + 
  scale_y_continuous(breaks = seq(0, 0.25, by = 0.05)) + 
  scale_shape_manual(values = c("Classic" = 1, "Amplification" = 2,
                                "Inf_CH" = 0, "Sup_CH" = 5)) + 
  scale_color_manual(values = custom_colors_slope)
ggsave("Bahadur_slope_comparison_a1_.25.75.png", width = 6, height = 3)

# a = 3
slopes_df_long_3 <- slopes_df_long[which(slopes_df_long$a == 0 | slopes_df_long$a == 3), ]

ggplot(slopes_df_long_3) + 
  geom_point(aes(x = Gamma, y = Slope, shape = Framework, color = Framework), size = 3) +
  geom_line(aes(x = Gamma, y = Slope, color = Framework), linewidth = 0.75) + 
  theme_bw() + 
  labs(y = "Bahadur Slope") + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.position = "blank",
        panel.grid.minor = element_blank()) + 
  guides(shape = guide_legend(nrow = 2, byrow = TRUE)) + 
  coord_cartesian(xlim = c(1, 1.7), ylim = c(0, 0.25)) +
  scale_x_continuous(breaks = c(1, 1.1, 1.3, 1.5, 1.7)) + 
  scale_y_continuous(breaks = seq(0, 0.25, by = 0.05)) + 
  scale_shape_manual(values = c("Classic" = 1, "Amplification" = 2,
                                "Inf_CH" = 0, "Sup_CH" = 5)) + 
  scale_color_manual(values = custom_colors_slope)
ggsave("Bahadur_slope_comparison_a3_.25.75.png", width = 6, height = 3)

#### Bahadur Slope, kappa = 0, 1 --------------------

# kappa = c(0, 1)
set.seed(12345)
kappa <- c(0.0, 0.1)
muv <- c(0.1, 0.5)
sigv <- c(1, 1)
#inf_gamma <- round(mean(infimum_sensitivity_CH_normal(kappa, muv, sigv)[(length(kappa) + 1):(2*length(kappa))]), 2)
gamma_grid <- c(1.1, 1.3, 1.5, 1.7)
a_grid <- c(1/3, 1, 3)

# initializer

slope_mat <- matrix(ncol = 6, nrow = 0)
n_mat <- matrix(ncol = 6, nrow = 0)

init_slope_rosenbaum <- bahadur_slope_rosenbaum_normal(kappa, muv, sigv, rho = 0, Gamma = 1)
for(a in a_grid){
  
  if(a <= 1){
    slope_CH <- bahadur_slope_CH_initial(kappa, muv, sigv, a)
    slope_mat <- rbind(slope_mat, c(1, a, init_slope_rosenbaum$slope, rep(slope_CH$slope, 3)))
    n_mat <- rbind(n_mat, c(1, a, init_slope_rosenbaum$min_n, rep(slope_CH$min_n, 3))) 
  }
  if(a > 1){
    slope_CH <- bahadur_slope_CH_initial(kappa, muv, sigv, a)
    slope_mat <- rbind(slope_mat, c(1, a, init_slope_rosenbaum$slope, rep(slope_CH$slope, 2), init_slope_rosenbaum$slope))
    n_mat <- rbind(n_mat, c(1, a, init_slope_rosenbaum$min_n, rep(slope_CH$min_n, 2), init_slope_rosenbaum$min_n)) 
  }
  
}
for(g in gamma_grid){
  
  for(a in a_grid){
    
    slope_rosenbaum <- bahadur_slope_rosenbaum_normal(kappa, muv, sigv, rho = 0, g)
    slope_CH <- bahadur_slope_CH_proj(kappa, muv, sigv, g, a)
    slope_mat <- rbind(slope_mat, c(g, a, slope_rosenbaum$slope, slope_CH$slope))
    n_mat <- rbind(n_mat, c(g, a, slope_rosenbaum$min_n, slope_CH$min_n))
    print(a)
  }
  print(g)
}

slope_mat <- data.frame(slope_mat)
names(slope_mat) <- c("Gamma", "a", "classic", "infCH", "amp", "supCH")
n_mat <- data.frame(n_mat)
names(n_mat) <- c("Gamma", "a", "classic", "infCH", "amp", "supCH")

rel_eff_mat <- round(do.call("cbind", list(slope_mat[, 1:2], slope_mat[, 3:6] / slope_mat[, 3])), 2)
names(rel_eff_mat) <- c("Gamma", "a", "Classic", "Inf_Slope_CH", "Amplification", "Sup_Slope_CH")
rel_n_mat <- round(do.call("cbind", list(n_mat[, 1:2], n_mat[, 3] / n_mat[, 3:6])), 2)
names(rel_n_mat) <- c("Gamma", "a", "Classic", "Inf_Slope_CH", "Amplification", "Sup_Slope_CH")
rel_eff_mat

write.csv(rel_eff_mat, "relative_efficiency_01.csv")
write.csv(rel_n_mat, "relative_n_required_01.csv")
write.csv(slope_mat, "slope_mat_01.csv")

# visual of the slopes
slope_mat <- read.csv("slope_mat_01.csv")[, -1]
slopes_df <- round(data.frame(slope_mat), 3)
names(slopes_df) <- c("Gamma", "a", "Classic", "Inf_CH", "Amplification", "Sup_CH")
slopes_df_long <- slopes_df %>%
  pivot_longer(cols = Classic:Sup_CH, 
               names_to = "Framework",
               values_to = "Slope")

# a = 1/3
slopes_df_long_third <- slopes_df_long[which(slopes_df_long$a < 0.5), ]

ggplot(slopes_df_long_third) + 
  geom_point(aes(x = Gamma, y = Slope, shape = Framework, color = Framework), size = 3) +
  geom_line(aes(x = Gamma, y = Slope, color = Framework), linewidth = 0.75) + 
  theme_bw() + 
  labs(y = "Bahadur Slope") + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.position = "blank",
        panel.grid.minor = element_blank()) + 
  guides(shape = guide_legend(nrow = 2, byrow = TRUE)) + 
  coord_cartesian(xlim = c(1, 1.7), ylim = c(0, 0.25)) +
  scale_x_continuous(breaks = c(1, 1.1, 1.3, 1.5, 1.7)) + 
  scale_y_continuous(breaks = seq(0, 0.25, by = 0.05)) + 
  scale_shape_manual(values = c("Classic" = 1, "Amplification" = 2,
                                "Inf_CH" = 0, "Sup_CH" = 5)) + 
  scale_color_manual(values = custom_colors_slope)
ggsave("Bahadur_slope_comparison_athird_01.png", width = 6, height = 3)

# a = 1
slopes_df_long_1 <- slopes_df_long[which(slopes_df_long$a == 0 | slopes_df_long$a == 1), ]

ggplot(slopes_df_long_1) + 
  geom_point(aes(x = Gamma, y = Slope, shape = Framework, color = Framework), size = 3) +
  geom_line(aes(x = Gamma, y = Slope, color = Framework), linewidth = 0.75) + 
  theme_bw() + 
  labs(y = "Bahadur Slope") + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.position = "blank",
        panel.grid.minor = element_blank()) + 
  guides(shape = guide_legend(nrow = 2, byrow = TRUE)) + 
  coord_cartesian(xlim = c(1, 1.7), ylim = c(0, 0.25)) +
  scale_x_continuous(breaks = c(1, 1.1, 1.3, 1.5, 1.7)) + 
  scale_y_continuous(breaks = seq(0, 0.25, by = 0.05)) + 
  scale_shape_manual(values = c("Classic" = 1, "Amplification" = 2,
                                "Inf_CH" = 0, "Sup_CH" = 5)) + 
  scale_color_manual(values = custom_colors_slope)
ggsave("Bahadur_slope_comparison_a1_01.png", width = 6, height = 3)

# a = 3
slopes_df_long_3 <- slopes_df_long[which(slopes_df_long$a == 0 | slopes_df_long$a == 3), ]

ggplot(slopes_df_long_3) + 
  geom_point(aes(x = Gamma, y = Slope, shape = Framework, color = Framework), size = 3) +
  geom_line(aes(x = Gamma, y = Slope, color = Framework), linewidth = 0.75) + 
  theme_bw() + 
  labs(y = "Bahadur Slope") + 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.position = "blank",
        panel.grid.minor = element_blank()) + 
  guides(shape = guide_legend(nrow = 2, byrow = TRUE)) + 
  coord_cartesian(xlim = c(1, 1.7), ylim = c(0, 0.25)) +
  scale_x_continuous(breaks = c(1, 1.1, 1.3, 1.5, 1.7)) + 
  scale_y_continuous(breaks = seq(0, 0.25, by = 0.05)) + 
  scale_shape_manual(values = c("Classic" = 1, "Amplification" = 2,
                                "Inf_CH" = 0, "Sup_CH" = 5)) + 
  scale_color_manual(values = custom_colors_slope)
ggsave("Bahadur_slope_comparison_a3_01.png", width = 6, height = 3)










