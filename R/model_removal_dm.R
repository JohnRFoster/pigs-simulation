modelCode <- nimbleCode({

  # priors
  for(i in 1:n_method){
    log_rho[i] ~ dnorm(0, tau = 0.1)
  }

  for(i in 1:3){
    p_mu[i] ~ dnorm(0, tau = 1)
    logit(p_unique[i]) <- p_mu[i]
  }

  for(i in 1:2){
    log_gamma[i] ~ dnorm(0, tau = 0.1)
  }

  # non time varying coefficients - observation model
  for(m in 1:n_method){
    beta1[m] ~ dnorm(0, tau = 1)
    for(i in 1:m_p){
      beta_p[m, i] ~ dnorm(0, tau = 1)
    }
  }

  # estimate apparent survival
  phi_mu ~ dbeta(phi_mu_a, phi_mu_b)
  psi_phi ~ dgamma(1, 0.1)
  a_phi <- phi_mu * psi_phi
  b_phi <- (1 - phi_mu) * psi_phi

  log_nu ~ dnorm(2, tau = 1)  # mean litter size
  log(nu) <- log_nu

  ## convert to expected number of pigs per primary period
  log_zeta <- log(pp_len) + log_nu - log(365)
  log(zeta) <- log_zeta
  for(i in 1:n_ls){
    J[i] ~ dpois(nu)
  }

  for(i in 1:n_survey){

    log_potential_area[i] <- calc_log_potential_area(
      log_rho = log_rho[1:n_method],
      log_gamma = log_gamma[1:2],
      p_unique = p_unique[1:3],
      log_effort_per = log_effort_per[i],
      effort_per = effort_per[i],
      n_trap_m1 = n_trap_m1[i],
      log_pi = log_pi,
      method = method[i]
    )

    # probability of capture, given that an individual is in the surveyed area
    log_theta[i] <- log(
      ilogit(beta1[method[i]] + inprod(X_p[county[i], 1:m_p], beta_p[method[i], 1:m_p]))
    ) +
      min(0, log_potential_area[i] - log_survey_area_km2[i])

    # likelihood
    y[i] ~ dpois(p[i] * (N[nH_p[i]] - y_sum[i]))

  }

  # the probability an individual is captured on the first survey
  for(i in 1:n_first_survey){
    log(p[first_survey[i]]) <- log_theta[first_survey[i]]
  }

  # the probability an individual is captured after the first survey
  for(i in 1:n_not_first_survey){
    log(p[not_first_survey[i]]) <- log_theta[start[not_first_survey[i]]] +
      sum(log(1 - exp(log_theta[start[not_first_survey[i]]:end[not_first_survey[i]]])))
  }

  for(i in 1:n_property){

    log_lambda_1[i] ~ dunif(0, 10)
    # log(N[i, all_pp[i, 1]]) <- log_lambda_1[i]
    log(N[nH[i, 1]]) <- log_lambda_1[i]

    # population growth across time steps
    for(j in 2:n_time_prop[i]){ # loop through every PP, including missing ones

      lambda[nH[i, j-1]] <- (N[nH[i, j-1]] - rem[i, j-1]) * zeta / 2 +
        (N[nH[i, j-1]] - rem[i, j-1]) * phi[nH[i, j-1]]

      N[nH[i, j]] ~ dpois(lambda[nH[i, j-1]])
      phi[nH[i, j-1]] ~ dbeta(a_phi, b_phi)

    }

  }

  # # for easier monitoring of abundance - long format
  # for(i in 1:n_units){
  #   xn[i] <- N[property_x[i], pp_x[i]]
  # }

})
