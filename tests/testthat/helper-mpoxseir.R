
reference_pars <- function() {
  dem_pars <- parameters_demographic()
  n_group <- dem_pars$n_group
  n_vax <- dem_pars$n_vax
  Ea0 <- matrix(0, n_group, n_vax)
  list(dt = 1,
       N = dem_pars$N,
       S0 = dem_pars$N - Ea0,
       Ea0 = Ea0,
       Eb0 = matrix(0, n_group, n_vax),
       Ir0 = matrix(0, n_group, n_vax),
       Id0 = matrix(0, n_group, n_vax),
       R0 =  matrix(0, n_group, n_vax),
       D0 =  matrix(0, n_group, n_vax),
       beta_h = 0.2 / 12.11,
       beta_z = rep(0.4 / 12.11, n_group),
       gamma_E = 0.05,
       gamma_I = 0.1,
       gamma_Ir = 0.1,
       gamma_Id = 0.05,
       CFR = matrix(c(1 / c(seq(2.5, 77.5, 5), 25, 25),
                      (0.5*1) / c(seq(2.5, 77.5, 5), 25, 25),
                      (0.25*1) / c(seq(2.5, 77.5, 5), 25, 25)),
                    n_group,n_vax),
       m = dem_pars$m,
       n_group = n_group,
       n_vax = n_vax)
}

reference_names <- function() {
  n_group <- parameters_demographic()$n_group
  n_vax <- parameters_demographic()$n_vax
  states <- c("S", "Ea", "Eb", "Ir", "Id", "R", "D")
  list(
    states = paste0(rep(states, each = n_group*n_vax), seq_len(n_group*n_vax))
  )
}
