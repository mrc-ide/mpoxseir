
reference_pars <- function() {
  dem_pars <- parameters_demographic()
  N_age <- dem_pars$N_age
  E10 <- rep(0, N_age)
  list(dt = 1,
       N = dem_pars$population,
       S0 = dem_pars$population - E10,
       E10 = E10,
       E20 = rep(0, N_age),
       Ir0 = rep(0, N_age),
       Id0 = rep(0, N_age),
       R0 = rep(0, N_age),
       D0 = rep(0, N_age),
       beta_h = 0.2 / 12.11,
       beta_z = 0.4 / 12.11,
       gamma_E = 0.05,
       gamma_I = 0.1,
       gamma_Ir = 0.1,
       gamma_Id = 0.05,
       CFR = 1 / seq(2.5, 77.5, 5),
       m = dem_pars$m,
       N_age = dem_pars$N_age)
}

reference_names <- function() {
  N_age <- parameters_demographic()$N_age
  states <- c("S", "E1", "E2", "Ir", "Id", "R", "D")
  list(
    states = paste0(rep(states, each = N_age), seq_len(N_age))
  )
}
