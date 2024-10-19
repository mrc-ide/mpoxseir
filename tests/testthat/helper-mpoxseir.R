

reference_pars_targeted_vax <- function(region = "equateur", uptake = 0) {
  pars <- parameters_fixed(region, initial_infections = 10)
  n_group <- pars$n_group
  n_vax <- pars$n_vax
  idx <- get_compartment_indices()
  age_bins <- get_age_bins()
  
  pars$beta_h <- pars$beta_s <- 0.2 / 12.11
  pars$beta_z <- rep(0.4 / 12.11, n_group)
  
  # vaccination
  
  ### these are the params we need to consider 
  
  # prioritisation_strategy_children = matrix(1, nrow = n_group, ncol = N_prioritisation_steps_children),
  # prioritisation_strategy_adults = matrix(1, nrow = n_group, ncol = N_prioritisation_steps_adults),
  # vaccine_uptake = rep(0, n_group),
  # ve_I = ve_I,
  # ve_T = ve_T,
  # vaccination_coverage_target_1st_dose_children_prop = 0,
  # vaccination_coverage_target_1st_dose_adults_prop = 0,
  # vaccination_coverage_target_1st_dose_adults_prop = 0,
  # vaccination_campaign_length_children = vaccination_campaign_length_children,
  # vaccination_campaign_length_adults = vaccination_campaign_length_adults,
  # daily_doses_children = matrix(0, nrow = vaccination_campaign_length_children,
  #                               ncol = n_vax),
  # daily_doses_adults = matrix(0, nrow = vaccination_campaign_length_adults,
  #                             ncol = n_vax))
  
  # children
  
  pars$m_sex["CSW", "PBS"] <- pars$m_sex["PBS", "CSW"] <-
    pars$m_sex["ASW", "PBS"] <- pars$m_sex["PBS", "ASW"] <- max(pars$m_gen_pop)

  pars$vaccination_campaign_length <- 10
  daily_doses <- matrix(0, nrow = pars$vaccination_campaign_length, ncol = n_vax)
  # nothing happens for j=1 and 2 doses are the maximum
  
  daily_doses[, idx$vax$unvaccinated] <- daily_doses[, idx$vax$one_dose] <- 1000
  daily_doses[pars$vaccination_campaign_length, ] <- 0


  pars$daily_doses <- daily_doses
  pars$N_prioritisation_steps <- 3
  
  vax_kids <- setNames(rep(0, n_group), names(idx$group))
  vax_kids[which(age_bins$end < 15)] <- 1
  vax_kids_SW <- vax_kids
  vax_kids_SW[c("CSW", "ASW")] <- 1
  
  pars$prioritisation_strategy <- cbind(
    vax_kids, # kids
    vax_kids_SW, # kids + SW
    rep(1, n_group)) # all
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
