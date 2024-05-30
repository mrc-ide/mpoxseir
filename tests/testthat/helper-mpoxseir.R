
reference_pars <- function() {
  dem_pars <- parameters_demographic()
  n_group <- dem_pars$n_group
  Ea0 <- rep(0, n_group)
  list(dt = 1,
       N = dem_pars$N,
       S0 = dem_pars$N - Ea0,
       Ea0 = Ea0,
       Eb0 = rep(0, n_group),
       Ir0 = rep(0, n_group),
       Id0 = rep(0, n_group),
       R0 =  rep(0, n_group),
       D0 =  rep(0, n_group),
       beta_h = 0.2 / 12.11,
       beta_z = rep(0.4 / 12.11,n_group),
       gamma_E = 0.05,
       gamma_I = 0.1,
       gamma_Ir = 0.1,
       gamma_Id = 0.05,
       CFR = 1 / c(seq(2.5, 77.5, 5), 25, 25),
       m = dem_pars$m,
       n_group = n_group)
}

reference_names <- function() {
  n_group <- parameters_demographic()$n_group
  states <- c("S", "Ea", "Eb", "Ir", "Id", "R", "D")
  list(
    states = paste0(rep(states, each = n_group), seq_len(n_group))
  )
}
