
reference_pars_generic_vax <- function() {
  dem_pars <- parameters_demographic()
  n_group <- dem_pars$n_group
  n_vax <- dem_pars$n_vax
  Ea0 <- matrix(0, n_group, n_vax)
  S0 <- Ea0
  S0[,1] <- dem_pars$N
  vaccination_campaign_length <- 10
  n_vaccination <- array(1000,
                        dim=c(n_group,
                              n_vax,
                              (vaccination_campaign_length)
                        ))
  n_vaccination[,,vaccination_campaign_length] <- 0
  n_vaccination[,n_vax,] <- 0

  out <- list(dt = 1,
       N = dem_pars$N,
       S0 = S0,
       Ea0 = Ea0,
       Eb0 = matrix(0, n_group, n_vax),
       Ir0 = matrix(0, n_group, n_vax),
       Id0 = matrix(0, n_group, n_vax),
       R0 =  matrix(0, n_group, n_vax),
       D0 =  matrix(0, n_group, n_vax),
       beta_h = 0.2 / 12.11,
       beta_z = rep(0.4 / 12.11, n_group),
       gamma_E = 0.05,
       #gamma_I = 0.1,
       gamma_Ir = 0.1,
       gamma_Id = 0.05,
       CFR = matrix(c(1 / c(seq(2.5, 77.5, 5), 25, 25),
                      (0.5*1) / c(seq(2.5, 77.5, 5), 25, 25)),
                    n_group,n_vax),
       ve_T = c(0,0.9),
       ve_I = c(0,0.8),
       m = dem_pars$m,
       n_group = n_group,
       n_vax = n_vax,
       #n_vaccination_allocation_SER = c(0.85,0.05,0.05,0.05),
       vaccination_campaign_length = vaccination_campaign_length,
       n_vaccination = n_vaccination
       )
#
#   if(sum(out$n_vaccination_allocation_SER)!=1){
#     stop("Entries in n_vaccination_allocation_SER must sum to 1.")
#   }

  return(out)
}

reference_pars_targeted_vax <- function() {
  dem_pars <- parameters_demographic()
  n_group <- dem_pars$n_group
  n_vax <- dem_pars$n_vax
  Ea0 <- matrix(0, n_group, n_vax)
  S0 <- Ea0
  S0[,1] <- dem_pars$N

  vaccination_campaign_length <- 10
  daily_doses <- matrix(c(rep(1000,vaccination_campaign_length),
                          rep(0,vaccination_campaign_length)),
                        ncol=n_vax,
                        nrow=vaccination_campaign_length)
  daily_doses[vaccination_campaign_length,] <- 0

  N_prioritisation_steps <- 3
  prioritisation_strategy <- cbind(c(1,1,1,rep(0,n_group-3)),# kids
                                   c(1,1,1,rep(0,n_group-5),1,1), # kids + CSW
                                   c(rep(1,n_group))) # all
  vaccination_coverage_target_prop <- 0.8
  vaccination_coverage_target <- round(prioritisation_strategy * vaccination_coverage_target_prop * dem_pars$N0)


  out <- list(dt = 1,
              N = dem_pars$N,
              S0 = S0,
              Ea0 = Ea0,
              Eb0 = matrix(0, n_group, n_vax),
              Ir0 = matrix(0, n_group, n_vax),
              Id0 = matrix(0, n_group, n_vax),
              R0 =  matrix(0, n_group, n_vax),
              D0 =  matrix(0, n_group, n_vax),
              beta_h = 0.2 / 12.11,
              beta_s = 0.2 / 12.11,
              beta_z = rep(0.4 / 12.11, n_group),
              gamma_E = 0.05,
              gamma_I = 0.1,
              gamma_Ir = 0.1,
              gamma_Id = 0.05,
              CFR = matrix(c(1 / c(seq(2.5, 77.5, 5), 25, 25),
                             (0.5*1) / c(seq(2.5, 77.5, 5), 25, 25)),
                           n_group,n_vax),
              ve_T = c(0,0.9),
              ve_I = c(0,0.8),
              m_gen_pop = dem_pars$m_gen_pop,
              m_sex = dem_pars$m_sex,
              n_group = n_group,
              n_vax = n_vax,
              #n_vaccination_allocation_SER = c(0.85,0.05,0.05,0.05),
              vaccination_campaign_length = vaccination_campaign_length,
              daily_doses = daily_doses,
              N_prioritisation_steps = N_prioritisation_steps,
              prioritisation_strategy = prioritisation_strategy,
              vaccination_coverage_target = vaccination_coverage_target,
              vaccine_uptake = rep(1,n_group)
  )

  return(out)
}



reference_names <- function() {
  n_group <- parameters_demographic()$n_group
  n_vax <- parameters_demographic()$n_vax
  states <- c("S", "Ea", "Eb", "Ir", "Id", "R", "D")
  list(
    states = paste0(rep(states, each = n_group*n_vax), seq_len(n_group*n_vax))
  )
}
