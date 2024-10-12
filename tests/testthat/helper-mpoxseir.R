
reference_pars_targeted_vax <- function(region = "equateur", uptake = 0) {
  pars <- parameters_fixed(region, initial_infections = 10)
  n_group <- pars$n_group
  n_vax <- pars$n_vax
  idx <- get_compartment_indices()
  
  pars$beta_h <- pars$beta_s <- 0.2 / 12.11
  pars$beta_z <- rep(0.4 / 12.11, n_group)

  daily_doses_time <- c(1, 11)
  daily_doses_value <- matrix(0, nrow = n_vax, ncol = 2)
  daily_doses_value[idx$vax$unvaccinated, 1] <- 
    daily_doses_value[idx$vax$one_dose, 1] <- 1000
  
  daily_doses <- t(vapply(seq_len(pars$n_vax),
                          function (i) {
                            approx(daily_doses_time, 
                                   daily_doses_value[i, ],
                                   xout = seq(daily_doses_time[1],
                                              daily_doses_time[2]), 
                                   method = "constant")$y},
                          numeric(max(daily_doses_time))))

  pars$daily_doses_time <- daily_doses_time
  pars$daily_doses_value <- daily_doses_value
  pars$daily_doses <- daily_doses
  pars$N_prioritisation_steps <- 3
  pars$prioritisation_strategy <- cbind(c(1,1,1,rep(0,n_group-3)),# kids
                                        c(1,1,1,rep(0,n_group-5),1,1), # kids + CSW
                                        c(rep(1,n_group))) # all
  pars$vaccination_coverage_target_1st_dose_prop <- 0.8
  pars$vaccination_coverage_target_2nd_dose_prop <- 0.5
  pars$vaccine_uptake[] <- uptake

  pars
}



reference_names <- function() {
  n_group <- parameters_demographic()$n_group
  n_vax <- parameters_demographic()$n_vax
  states <- c("S", "Ea", "Eb", "Ir", "Id", "R", "D")
  list(
    states = paste0(rep(states, each = n_group*n_vax), seq_len(n_group*n_vax))
  )
}
