################################################################################
####       Confounding Heterogeneity in Matched Observational Studies       ####
####                 Reproducible Code for Main Calculations                ####
####                           Design Sensitivity                           ####
################################################################################

#### README --------------------------------------------------------------------

# running this file produces the following outputs:

# 1) "Design_Sensitivity_(Normal / t)_(Perm / Rank)_Setting(1 / 2).csv":
# Reproduces the calculated / simulated (via Monte Carlo) design sensitivity
# for all simulation settings in the main and supplementary files.

# 2) "Design_Sensitivity_amp.png", "Design_Sensitivity_crosssection_CH2.png",
# and "Design_Sensitivity_crosssection_CH1.png": Reproduces the three subfigures
# shown in Figure 2 of the main text.

# all files are placed in the folder "Confounding_Heterogeneity_Design_Sensitivity".

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
dir.create("Confounding_Heterogeneity_Design_Sensitivity")
setwd("Confounding_Heterogeneity_Design_Sensitivity")

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
expected_abs_sum_tdist <- function(kappa, muv, degree = 5, n = 1e7) {
  # monte-carlo approximation to absolute value of sum of two t-distributions
  # kappa = vector of kappas
  # muv = vector of means of the t-distributions
  # degree = degrees of freedom for both t-distributions
  # n = number of samples for monte carlo approximation
  
  x1 <- muv[1] + rt(n, df = degree)
  x2 <- muv[2] + rt(n, df = degree)
  z <- kappa[1] * x1 + kappa[2] * x2
  return(mean(abs(z)))
}
expected_abs_tdist <- function(kappa, mu, degree = 5, n = 1e7) {
  # monte-carlo approximation of absolute value of a single t-distribution
  # kappa = a singular kappa
  # mu = a singular mean 
  # degree = degrees of freedom for the t-distribution
  # n = number of samples for monte carlo approximation
  
  x <- mu + rt(n, df = degree)
  z <- kappa * x
  return(mean(abs(z)))
}
expected_abs_sum_correlated_tdist <- function(kappa, muv, rho, degree = 5, n = 1e7){
  
  Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  L <- chol(Sigma)
  Z <- matrix(rnorm(2 * n), ncol = 2) %*% L
  U <- rchisq(n, df = degree)
  scale <- sqrt(degree / U)
  
  x1 <- muv[1] + Z[, 1] * scale
  x2 <- muv[2] + Z[, 2] * scale
  
  return(mean(abs(kappa[1] * x1 + kappa[2] * x2)))
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
phi_kappa_vector_tdist <- function(kappa, muv, degree){
  exp_absv <- c()
  for(k in 1:length(kappa)){
    exp_absv[k] <- expected_abs_tdist(kappa[k], muv[k], degree)
  }
  return(exp_absv)
}
eta_kappa_tdist <- function(kappa, muv, degree){
  eta <- expected_abs_sum_tdist(kappa, muv, degree)
  return(eta)
}
eta_kappa_correlated_tdist <- function(kappa, muv, rho, degree){
  eta <- expected_abs_sum_correlated_tdist(kappa, muv, rho, degree)
  return(eta)
}

design_sensitivity_rosenbaum_normal <- function(kappa, muv, sigv, rho = 0){
  mu <- mu_kappa(kappa, muv)
  eta <- eta_kappa_normal(kappa, muv, sigv, rho)
  Gamma <- (eta + mu)/(eta - mu)
  return(Gamma)
}
infimum_sensitivity_amp_normal <- function(kappa, muv, sigv){
  mu <- mu_kappa(kappa, muv)
  phiv <- phi_kappa_vector_normal(kappa, muv, sigv)
  phi <- sum(phiv)
  return(rep((phi + mu)/(phi - mu), length(kappa)))
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

design_sensitivity_rosenbaum_tdist <- function(kappa, muv, degree){
  mu <- mu_kappa(kappa, muv)
  eta <- eta_kappa_tdist(kappa, muv, degree)
  Gamma <- (eta + mu)/(eta - mu)
  return(Gamma)
}
design_sensitivity_rosenbaum_correlated_tdist <- function(kappa, muv, rho, degree){
  mu <- mu_kappa(kappa, muv)
  eta <- eta_kappa_correlated_tdist(kappa, muv, rho, degree)
  Gamma <- (eta + mu)/(eta - mu)
  return(Gamma)
}
infimum_sensitivity_amp_tdist <- function(kappa, muv, degree){
  mu <- mu_kappa(kappa, muv)
  phiv <- phi_kappa_vector_tdist(kappa, muv, degree)
  phi <- sum(phiv)
  return(rep((phi + mu)/(phi - mu), length(kappa)))
}
infimum_sensitivity_CH_tdist <- function(kappa, muv, degree){
  
  # stores \phi_{\kappa, k}
  mu <- mu_kappa(kappa, muv)
  phiv <- phi_kappa_vector_tdist(kappa, muv, degree)
  phi <- sum(phiv)
  
  # asymptote for Lambda
  Lambda_asymptotev <- rep(((phi + mu)/(phi - mu))^(1/2), length(kappa))
  
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

design_curve_amp_normal <- function(kappa, muv, sigv, Lam_min = -1){
  
  mu <- mu_kappa(kappa, muv)
  phi <- sum(phi_kappa_vector_normal(kappa, muv, sigv))
  square <- sqrt(mu/phi)
  square_corner <- (square + 1)/(1 - square)
  
  if(Lam_min == -1){
    Lam_min <- ceiling(asymptote)
  }
  
  asymptote <- (phi + mu)/(phi - mu)
  Lam_grid <- seq(Lam_min, square_corner, length.out = 25)
  rhs <- ((Lam_grid + 1)/(Lam_grid - 1))*(mu / phi)
  Delta_grid <- (rhs + 1)/(1 - rhs)
  first_half <- cbind(Lam_grid, Delta_grid)
  second_half <- cbind(rev(Delta_grid)[-1], rev(Lam_grid)[-1])
  curve <- rbind(first_half, second_half)
  return(list(curve = curve, asymptote = asymptote))
}
design_curve_crossection_CH2 <- function(kappa, muv, sigv, Lambda, maxv){
  
  # maxv: the largest value of delta1 and delta2 you want to see
  
  mu <- mu_kappa(kappa, muv)
  phiv <- phi_kappa_vector_normal(kappa, muv, sigv)
  Lambda_const <- (Lambda - 1)/(Lambda + 1)
  rh_const <- mu/Lambda_const
  
  # calculate local asymptotes
  Del1_const <- (rh_const - phiv[2])/phiv[1]
  Del1_asymptote <- (Del1_const + 1)/(1-Del1_const)
  Del2_const <- (rh_const - phiv[1])/phiv[2]
  Del2_asymptote <- (Del2_const + 1)/(1-Del2_const)
  
  # construct points along the design sensitivity
  Del1_grid <- seq(ceiling(Del1_asymptote), ceiling(Del1_asymptote)+20, length.out = 50)
  Del2_grid <- seq(ceiling(Del2_asymptote), ceiling(Del2_asymptote)+20, length.out = 50)
  
  first_half_rhs <- (rh_const - (Del1_grid - 1)/(Del1_grid + 1)*phiv[1])/phiv[2]
  first_half <- (first_half_rhs + 1)/(1 - first_half_rhs)
  first_half_df <- matrix(cbind(Del1_grid, first_half), ncol = 2)
  
  second_half_rhs <- (rh_const - (Del2_grid - 1)/(Del2_grid + 1)*phiv[2])/phiv[1]
  second_half <- (second_half_rhs + 1)/(1 - second_half_rhs)
  second_half_df <- matrix(cbind(second_half, Del2_grid), ncol = 2)
  
  Design_df <- rbind(first_half_df, second_half_df)
  Design_df <- Design_df[which(Design_df[, 1] < maxv[1] & Design_df[, 2] < maxv[2]), ]
  Design_df <- data.frame(Design_df[order(Design_df[, 1]), ])
  names(Design_df) <- c("Delta1", "Delta2")
  
  output <- list(curve = Design_df, local_asymptotes = c(Del1_asymptote, Del2_asymptote))
  return(output)
  
}
design_curve_crossection_CH1 <- function(kappa, muv, sigv, Deltav){
  
  mu <- mu_kappa(kappa, muv)
  phiv <- phi_kappa_vector_normal(kappa, muv, sigv)
  phi_const <- ((Deltav[1] - 1)/(Deltav[1] + 1))*phiv[1] + ((Deltav[2] - 1)/(Deltav[2] + 1))*phiv[2]
  rh_const <- mu/phi_const
  
  asymptote <- (rh_const + 1)/(1 - rh_const)
  
  Lambda1_first_half <- seq(1, sqrt(asymptote), length.out = 25)
  Lambda2_first_half <- (rh_const + 1)/(Lambda1_first_half  * (1 - rh_const))
  first_half <- cbind(Lambda1_first_half, Lambda2_first_half)
  second_half <- cbind(rev(Lambda2_first_half)[-1], rev(Lambda1_first_half)[-1])
  
  Design_df <- data.frame(rbind(first_half, second_half))
  names(Design_df) <- c("Lambda1", "Lambda2")
  
  output <- list(curve = Design_df, upper_boundary = asymptote)
  return(output)
  
} 

mu_kappa_rank_normal <- function(kappa, muv, sigv, n = 1e7){
  
  Z1 <- rnorm(n, muv[1], sigv[1])
  Z2 <- rnorm(n, muv[2], sigv[2])
  mu1 <- kappa[1] * mean(pfoldnorm(abs(Z1), muv[1], sigv[1]) * sign(Z1))
  mu2 <- kappa[2] * mean(pfoldnorm(abs(Z2), muv[2], sigv[2]) * sign(Z2))
  return(mu1 + mu2)
  
}
eta_kappa_rank_normal <- function(kappa, muv, sigv, rho = 0, n = 1e7){
  
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
  
  UZ1 <- kappa[1] * pfoldnorm(abs(Z[, 1]), muv[1], sigv[1]) * sign(Z[, 1])
  UZ2 <- kappa[2] * pfoldnorm(abs(Z[, 2]), muv[2], sigv[2]) * sign(Z[, 2])
  
  return(mean(abs(UZ1 + UZ2)))
  
}
mu_kappa_rank_tdist <- function(kappa, muv, degree, n = 1e7){
  
  t1 <- muv[1] + rt(n, degree)
  t2 <- muv[2] + rt(n, degree)
  mu1 <- kappa[1]*mean((pt(abs(t1) - muv[1], degree) - pt(-abs(t1) - muv[1], degree)) * sign(t1))
  mu2 <- kappa[2]*mean((pt(abs(t2) - muv[2], degree) - pt(-abs(t2) - muv[2], degree)) * sign(t2))
  return(mu1 + mu2)
  
} 
eta_kappa_rank_tdist <- function(kappa, muv, degree, n = 1e7){
  
  t1 <- muv[1] + rt(n, degree)
  t2 <- muv[2] + rt(n, degree)
  UT1 <- kappa[1]*(pt(abs(t1) - muv[1], degree) - pt(-abs(t1) - muv[1], degree)) * sign(t1)
  UT2 <- kappa[2]*(pt(abs(t2) - muv[2], degree) - pt(-abs(t2) - muv[2], degree)) * sign(t2)
  
  return(mean(abs(UT1 + UT2)))
  
  
}
eta_kappa_rank_correlated_tdist <- function(kappa, muv, degree, rho = 0, n = 1e7){
  
  # generating correlated t dists
  Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  L <- chol(Sigma)
  Z <- matrix(rnorm(2 * n), ncol = 2) %*% L
  U <- rchisq(n, df = degree)
  scale <- sqrt(degree / U)
  
  t1 <- muv[1] + Z[, 1] * scale
  t2 <- muv[2] + Z[, 2] * scale
  
  UT1 <- kappa[1]*(pt(abs(t1) - muv[1], degree) - pt(-abs(t1) - muv[1], degree)) * sign(t1)
  UT2 <- kappa[2]*(pt(abs(t2) - muv[2], degree) - pt(-abs(t2) - muv[2], degree)) * sign(t2)
  
  return(mean(abs(UT1 + UT2)))
}
phi_kappa_vector_rank <- function(kappa){
  return(kappa*0.5)
}

design_sensitivity_rosenbaum_rank_normal <- function(kappa, muv, sigv, rho = 0, n = 1e7){
  mu <- mu_kappa_rank_normal(kappa, muv, sigv, n)
  eta <- eta_kappa_rank_normal(kappa, muv, sigv, rho, n)
  Gamma <- (eta + mu)/(eta - mu)
  return(Gamma)
}
infimum_sensitivity_amp_rank_normal <- function(kappa, muv, sigv, n = 1e7){
  mu <- mu_kappa_rank_normal(kappa, muv, sigv, n)
  phiv <- phi_kappa_vector_rank(kappa)
  phi <- sum(phiv)
  return(rep((phi + mu)/(phi - mu), length(kappa)))
}
infimum_sensitivity_CH_rank_normal <- function(kappa, muv, sigv, n = 1e7){
  
  # stores \phi_{\kappa, k}
  mu <- mu_kappa_rank_normal(kappa, muv, sigv, n)
  phiv <- phi_kappa_vector_rank(kappa)
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

design_sensitivity_rosenbaum_rank_tdist <- function(kappa, muv, degree, n = 1e7){
  mu <- mu_kappa_rank_tdist(kappa, muv, degree, n)
  eta <- eta_kappa_rank_tdist(kappa, muv, degree, n)
  Gamma <- (eta + mu)/(eta - mu)
  return(Gamma)
}
design_sensitivity_rosenbaum_rank_correlated_tdist <- function(kappa, muv, rho, degree, n = 1e7){
  mu <- mu_kappa_rank_tdist(kappa, muv, degree, n)
  eta <- eta_kappa_rank_correlated_tdist(kappa, muv, degree, rho, n)
  Gamma <- (eta + mu)/(eta - mu)
  return(Gamma)
}
infimum_sensitivity_amp_rank_tdist <- function(kappa, muv, sigv, n = 1e7){
  mu <- mu_kappa_rank_tdist(kappa, muv, sigv, n)
  phiv <- phi_kappa_vector_rank(kappa)
  phi <- sum(phiv)
  return(rep((phi + mu)/(phi - mu), length(kappa)))
}
infimum_sensitivity_CH_rank_tdist <- function(kappa, muv, sigv, n = 1e7){
  
  # stores \phi_{\kappa, k}
  mu <- mu_kappa_rank_tdist(kappa, muv, sigv, n)
  phiv <- phi_kappa_vector_rank(kappa)
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


#### Design Sensitivity, Normally Distributed Outcome --------------------------

# normal dist, permutation test stat, setting 1
muv <- c(0.25, 0.25)
sigv <- c(1, 1)
kappa_mat <- matrix(c(1, 0, 0.75, 0.25, 0.5, 0.5, 0.25, 0.75, 0, 1), nrow = 5, byrow = TRUE)

output <- data.frame(matrix(nrow = 0, ncol = 11))
for(k in 1:nrow(kappa_mat)){
  kappa <- kappa_mat[k, ]
  rosenbaum_negcorr <- design_sensitivity_rosenbaum_normal(kappa, muv, sigv, rho = -0.5)
  rosenbaum <- design_sensitivity_rosenbaum_normal(kappa, muv, sigv)
  rosenbaum_poscorr <- design_sensitivity_rosenbaum_normal(kappa, muv, sigv, rho = 0.5)
  amp <- infimum_sensitivity_amp_normal(kappa, muv, sigv)
  CH <- infimum_sensitivity_CH_normal(kappa, muv, sigv)
  output[k, ] <- c(kappa, rosenbaum_negcorr, rosenbaum, rosenbaum_poscorr, amp, CH)
}

names(output) <- c("kappa1", "kappa2", "Gamma_rho-0.5", "Gamma_rho0", "Gamma_rho0.5", "Lambda", "Delta", "Lambda1", "Lambda2", "Delta1", "Delta2")
round(output, 2)
write.csv(output, "Design_Sensitivity_Normal_Perm_setting1.csv")

# normal dist, permutation test stat, setting 2
muv <- c(0.1, 0.5)
sigv <- c(1, 1)
kappa_mat <- matrix(c(1, 0, 0.75, 0.25, 0.5, 0.5, 0.25, 0.75, 0, 1), nrow = 5, byrow = TRUE)

output <- data.frame(matrix(nrow = 0, ncol = 11))
for(k in 1:nrow(kappa_mat)){
  kappa <- kappa_mat[k, ]
  rosenbaum_negcorr <- design_sensitivity_rosenbaum_normal(kappa, muv, sigv, rho = -0.5)
  rosenbaum <- design_sensitivity_rosenbaum_normal(kappa, muv, sigv)
  rosenbaum_poscorr <- design_sensitivity_rosenbaum_normal(kappa, muv, sigv, rho = 0.5)
  amp <- infimum_sensitivity_amp_normal(kappa, muv, sigv)
  CH <- infimum_sensitivity_CH_normal(kappa, muv, sigv)
  output[k, ] <- c(kappa, rosenbaum_negcorr, rosenbaum, rosenbaum_poscorr, amp, CH)
}

names(output) <- c("kappa1", "kappa2", "Gamma_rho-0.5", "Gamma_rho0", "Gamma_rho0.5", "Lambda", "Delta", "Lambda1", "Lambda2", "Delta1", "Delta2")
round(output, 2)
write.csv(output, "Design_Sensitivity_Normal_Perm_setting2.csv")

set.seed(12345)

# normal dist, rank stat, setting 1
muv <- c(0.25, 0.25)
sigv <- c(1, 1)
kappa_mat <- matrix(c(1, 0, 0.75, 0.25, 0.5, 0.5, 0.25, 0.75, 0, 1), nrow = 5, byrow = TRUE)

output <- data.frame(matrix(nrow = 0, ncol = 11))
for(k in 1:nrow(kappa_mat)){
  kappa <- kappa_mat[k, ]
  rosenbaum_negcorr <- design_sensitivity_rosenbaum_rank_normal(kappa, muv, sigv, rho = -0.5)
  rosenbaum <- design_sensitivity_rosenbaum_rank_normal(kappa, muv, sigv)
  rosenbaum_poscorr <- design_sensitivity_rosenbaum_rank_normal(kappa, muv, sigv, rho = 0.5)
  amp <- infimum_sensitivity_amp_rank_normal(kappa, muv, sigv)
  CH <- infimum_sensitivity_CH_rank_normal(kappa, muv, sigv)
  output[k, ] <- c(kappa, rosenbaum_negcorr, rosenbaum, rosenbaum_poscorr, amp, CH)
}

names(output) <- c("kappa1", "kappa2", "Gamma_rho-0.5", "Gamma_rho0", "Gamma_rho0.5", "Lambda", "Delta", "Lambda1", "Lambda2", "Delta1", "Delta2")
round(output, 2)
write.csv(output, "Design_Sensitivity_Normal_Rank_setting1.csv")

# normal dist, rank stat, setting 2
muv <- c(0.1, 0.5)
sigv <- c(1, 1)
kappa_mat <- matrix(c(1, 0, 0.75, 0.25, 0.5, 0.5, 0.25, 0.75, 0, 1), nrow = 5, byrow = TRUE)

output <- data.frame(matrix(nrow = 0, ncol = 11))
for(k in 1:nrow(kappa_mat)){
  kappa <- kappa_mat[k, ]
  rosenbaum_negcorr <- design_sensitivity_rosenbaum_rank_normal(kappa, muv, sigv, rho = -0.5)
  rosenbaum <- design_sensitivity_rosenbaum_rank_normal(kappa, muv, sigv)
  rosenbaum_poscorr <- design_sensitivity_rosenbaum_rank_normal(kappa, muv, sigv, rho = 0.5)
  amp <- infimum_sensitivity_amp_rank_normal(kappa, muv, sigv)
  CH <- infimum_sensitivity_CH_rank_normal(kappa, muv, sigv)
  output[k, ] <- c(kappa, rosenbaum_negcorr, rosenbaum, rosenbaum_poscorr, amp, CH)
}

names(output) <- c("kappa1", "kappa2", "Gamma_rho-0.5", "Gamma_rho0", "Gamma_rho0.5", "Lambda", "Delta", "Lambda1", "Lambda2", "Delta1", "Delta2")
round(output, 2)
write.csv(output, "Design_Sensitivity_Normal_Rank_setting2.csv")

#### Design Sensitivity Figures ------------------------------------------------

# figure a) amplification case

amp_design <- design_curve_amp_normal(kappa, muv, sigv)
ampdf <- amp_design$curve

ggplot(ampdf) + 
  geom_line(aes(x = Lam_grid, y = Delta_grid), linewidth = 0.75) + 
  geom_vline(xintercept = amp_design$asymptote, linetype = "dashed") +
  geom_hline(yintercept = amp_design$asymptote, linetype = "dashed") +
  geom_ribbon(aes(x = Lam_grid, y = Delta_grid, xmin = 1, xmax = Lam_grid), fill = "blue3", alpha = 0.1) +
  geom_ribbon(aes(x = Lam_grid, y = Delta_grid, ymin = 1, ymax = 3), fill = "blue3", alpha = 0.1) +
  annotate("rect", xmin = 1, xmax = 3, ymin = 1, ymax = 3, fill = "blue3", alpha = 0.1) +
  theme_bw() + 
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 12))
#scale_x_continuous(expand = c(0, 0)) +
#scale_y_continuous(expand = c(0, 0))
ggsave("Design_Sensitivity_amp.png", width = 4, height = 4)

# figure b) type II confounding heterogeneity, cross section at \Lambda = 3.5

Lambda <- 3.5
CH2_design <- design_curve_crossection_CH2(kappa, muv, sigv, Lambda, c(30, 30))
ch2df <- CH2_design$curve
infimum_sensitivity <- round(infimum_sensitivity_CH_normal(kappa, muv, sigv)[3:4], 2)

ggplot(ch2df) +
  geom_line(aes(x = Delta1, y = Delta2), linewidth = 0.75) +
  geom_vline(xintercept = infimum_sensitivity[1], linetype = "dashed") +
  geom_hline(yintercept = infimum_sensitivity[2], linetype = "dashed") +
  geom_ribbon(aes(x = Delta1, y = Delta2, xmin = 1, xmax = Delta1), fill = "lightblue3", alpha = 0.1) +
  geom_ribbon(aes(x = Delta1, y = Delta2, ymin = 1, ymax = min(ch2df[, 2])), fill = "lightblue3", alpha = 0.1) +
  annotate("rect", xmin = 1, xmax = min(ch2df[, 1]), ymin = 1, ymax = min(ch2df[,2]), fill = "lightblue3", alpha = 0.1) +
  theme_bw() + 
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 12)) # +
# scale_x_continuous(expand = c(0, 0)) +
# scale_y_continuous(expand = c(0, 0))
ggsave("Design_Sensitivity_crosssection_CH2.png", width = 4, height = 4)


# figure c) Type I confounding heterogeneity, cross section at \Delta1 = 2, \Delta2 = 5
Deltav <- c(2, 5)
CH1_design <- design_curve_crossection_CH1(kappa, muv, sigv, Deltav)
ch1df <- CH1_design$curve
infimum_sensitivity_lambda <- infimum_sensitivity_CH_normal(kappa, muv, sigv)[1:2]

ggplot(ch1df) + 
  geom_line(aes(x = Lambda1, y = Lambda2), linewidth = 0.75) + 
  annotate("point", x = infimum_sensitivity_lambda[1], infimum_sensitivity_lambda[2], shape = 1, size = 3) +
  geom_ribbon(aes(x = Lambda1, y = Lambda2, xmin = 1, xmax = Lambda1), fill = "lightblue3", alpha = 0.1) + 
  geom_vline(xintercept = CH1_design$upper_boundary, linetype = "dotted") + 
  geom_hline(yintercept = CH1_design$upper_boundary, linetype = "dotted") +
  theme_bw() + 
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 12))
# scale_x_continuous(expand = c(0, 0)) +
# scale_y_continuous(expand = c(0, 0))
ggsave("Design_Sensitivity_crosssection_CH1.png", width = 4, height = 4)












#### Design Sensitivity, t-Distributed Outcome ---------------------------------

set.seed(12345)
# t dist, setting 1
muv <- c(0.25, 0.25) * sqrt(5/3)
degree <- 5
kappa_mat <- matrix(c(1, 0, 0.75, 0.25, 0.5, 0.5, 0.25, 0.75, 0, 1), nrow = 5, byrow = TRUE)

output <- data.frame(matrix(nrow = 0, ncol = 11))
for(k in 1:nrow(kappa_mat)){
  kappa <- kappa_mat[k, ]
  rosenbaum_negcorr <- design_sensitivity_rosenbaum_correlated_tdist(kappa, muv, -0.5, degree)  
  rosenbaum <- design_sensitivity_rosenbaum_tdist(kappa, muv, degree)
  rosenbaum_poscorr <- design_sensitivity_rosenbaum_correlated_tdist(kappa, muv, 0.5, degree)
  amp <- infimum_sensitivity_amp_tdist(kappa, muv, degree)
  CH <- infimum_sensitivity_CH_tdist(kappa, muv, degree)
  output[k, ] <- c(kappa, rosenbaum_negcorr, rosenbaum, rosenbaum_poscorr, amp, CH)
}

names(output) <- c("kappa1", "kappa2", "Gamma_rho-0.5", "Gamma_rho0", "Gamma_rho0.5", "Lambda", "Delta", "Lambda1", "Lambda2", "Delta1", "Delta2")
round(output, 2)
write.csv(output, "Design_Sensitivity_t_Perm_setting1.csv")

# t dist, setting 2
muv <- c(0.1, 0.5) * sqrt(5/3)
degree <- 5
kappa_mat <- matrix(c(1, 0, 0.75, 0.25, 0.5, 0.5, 0.25, 0.75, 0, 1), nrow = 5, byrow = TRUE)

output <- data.frame(matrix(nrow = 0, ncol = 11))
for(k in 1:nrow(kappa_mat)){
  kappa <- kappa_mat[k, ]
  rosenbaum_negcorr <- design_sensitivity_rosenbaum_correlated_tdist(kappa, muv, -0.5, degree)  
  rosenbaum <- design_sensitivity_rosenbaum_tdist(kappa, muv, degree)
  rosenbaum_poscorr <- design_sensitivity_rosenbaum_correlated_tdist(kappa, muv, 0.5, degree)
  amp <- infimum_sensitivity_amp_tdist(kappa, muv, degree)
  CH <- infimum_sensitivity_CH_tdist(kappa, muv, degree)
  output[k, ] <- c(kappa, rosenbaum_negcorr, rosenbaum, rosenbaum_poscorr, amp, CH)
}

names(output) <- c("kappa1", "kappa2", "Gamma_rho-0.5", "Gamma_rho0", "Gamma_rho0.5", "Lambda", "Delta", "Lambda1", "Lambda2", "Delta1", "Delta2")
round(output, 2)
write.csv(output, "Design_Sensitivity_t_Perm_setting2.csv")

#### rank

set.seed(12345)

# t dist, setting 1
muv <- c(0.25, 0.25) * sqrt(5/3)
degree <- 5
kappa_mat <- matrix(c(1, 0, 0.75, 0.25, 0.5, 0.5, 0.25, 0.75, 0, 1), nrow = 5, byrow = TRUE)

output <- data.frame(matrix(nrow = 0, ncol = 11))
for(k in 1:nrow(kappa_mat)){
  kappa <- kappa_mat[k, ]
  rosenbaum_negcorr <- design_sensitivity_rosenbaum_rank_correlated_tdist(kappa, muv, -0.5, degree)  
  rosenbaum <- design_sensitivity_rosenbaum_rank_tdist(kappa, muv, degree)
  rosenbaum_poscorr <- design_sensitivity_rosenbaum_rank_correlated_tdist(kappa, muv, 0.5, degree)
  amp <- infimum_sensitivity_amp_rank_tdist(kappa, muv, degree)
  CH <- infimum_sensitivity_CH_rank_tdist(kappa, muv, degree)
  output[k, ] <- c(kappa, rosenbaum_negcorr, rosenbaum, rosenbaum_poscorr, amp, CH)
}

names(output) <- c("kappa1", "kappa2", "Gamma_rho-0.5", "Gamma_rho0", "Gamma_rho0.5", "Lambda", "Delta", "Lambda1", "Lambda2", "Delta1", "Delta2")
round(output, 2)
write.csv(output, "Design_Sensitivity_t_Rank_setting1.csv")

# t dist, setting 2
muv <- c(0.1, 0.5) * sqrt(5/3)
degree <- 5
kappa_mat <- matrix(c(1, 0, 0.75, 0.25, 0.5, 0.5, 0.25, 0.75, 0, 1), nrow = 5, byrow = TRUE)

output <- data.frame(matrix(nrow = 0, ncol = 11))
for(k in 1:nrow(kappa_mat)){
  kappa <- kappa_mat[k, ]
  rosenbaum_negcorr <- design_sensitivity_rosenbaum_rank_correlated_tdist(kappa, muv, -0.5, degree)  
  rosenbaum <- design_sensitivity_rosenbaum_rank_tdist(kappa, muv, degree)
  rosenbaum_poscorr <- design_sensitivity_rosenbaum_rank_correlated_tdist(kappa, muv, 0.5, degree)
  amp <- infimum_sensitivity_amp_rank_tdist(kappa, muv, degree)
  CH <- infimum_sensitivity_CH_rank_tdist(kappa, muv, degree)
  output[k, ] <- c(kappa, rosenbaum_negcorr, rosenbaum, rosenbaum_poscorr, amp, CH)
}

names(output) <- c("kappa1", "kappa2", "Gamma_rho-0.5", "Gamma_rho0", "Gamma_rho0.5", "Lambda", "Delta", "Lambda1", "Lambda2", "Delta1", "Delta2")
round(output, 2)
write.csv(output, "Design_Sensitivity_t_Rank_setting2.csv")






















