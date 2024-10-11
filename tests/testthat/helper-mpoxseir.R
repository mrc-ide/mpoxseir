
reference_pars_targeted_vax <- function() {
  dem_pars <- parameters_demographic()
  n_group <- dem_pars$n_group
  n_vax <- dem_pars$n_vax
  Ea0 <- matrix(0, n_group, n_vax)
  S0 <- Ea0
  #S0[,2] <- dem_pars$N
  S0[,2] <- rbinom(size=dem_pars$N,
                   prob=dem_pars$sus_prop,n=length(dem_pars$N))
  S0[,1] <- dem_pars$N-S0[,2]
  # rowSums(S0)==dem_pars$N0
  dem_pars$m_sex["SW","PBS"] <- dem_pars$m_sex["PBS","SW"] <- max(dem_pars$m_gen_pop)

  daily_doses_time <- c(1, 11)
  daily_doses_value <- matrix(0, ncol = n_vax, nrow = 2)
  # nothing happens for j=1 and 2 doses are the maximum
  daily_doses_value[1, 2] <- daily_doses_value[1, 3] <- 1000

  N_prioritisation_steps <- 3
  prioritisation_strategy <- cbind(c(1,1,1,rep(0,n_group-3)),# kids
                                   c(1,1,1,rep(0,n_group-5),1,1), # kids + CSW
                                   c(rep(1,n_group))) # all
  vaccination_coverage_target_1st_dose_prop <- 0.8
  vaccination_coverage_target_2nd_dose_prop <- 0.5
  #n_group,N_prioritisation_steps,n_vax
  # vaccination_coverage_target <- array(0,dim=c(n_group,N_prioritisation_steps,
  #                                              n_vax))
  #
  # vaccination_coverage_target[,,3] <- round(prioritisation_strategy * vaccination_coverage_target_1st_dose_prop * dem_pars$N0)
  # vaccination_coverage_target[,,4] <- round(prioritisation_strategy * vaccination_coverage_target_2nd_dose_prop * dem_pars$N0)



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
              ve_T = c(0.8,0,0.8,0.9),
              ve_I = c(0.8,0,0.8,0.9),
              m_gen_pop = dem_pars$m_gen_pop,
              m_sex = dem_pars$m_sex,
              n_group = n_group,
              n_vax = n_vax,
              #n_vaccination_allocation_SER = c(0.85,0.05,0.05,0.05),
              daily_doses_time = daily_doses_time,
              daily_doses_value = daily_doses_value,
              N_prioritisation_steps = N_prioritisation_steps,
              prioritisation_strategy = prioritisation_strategy,
              #vaccination_coverage_target = vaccination_coverage_target,
              vaccine_uptake = rep(1,n_group),
              vaccination_coverage_target_1st_dose_prop = vaccination_coverage_target_1st_dose_prop,
              vaccination_coverage_target_2nd_dose_prop = vaccination_coverage_target_2nd_dose_prop
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
