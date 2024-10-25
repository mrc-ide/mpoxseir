

reference_pars_targeted_vax <- function(region = "equateur") {
  
  pars <- parameters_fixed(region, initial_infections = 10)
  n_group <- pars$n_group
  n_vax <- pars$n_vax
  idx <- get_compartment_indices()
  age_bins <- get_age_bins()
  
  # these parameters are fitted so give dummy values for testing purposes
  pars$beta_h <- pars$beta_s <- 0.2 / 12.11
  pars$beta_z <- rep(0.4 / 12.11, n_group)
  pars$m_sex["CSW", "PBS"] <- pars$m_sex["PBS", "CSW"] <- 
    pars$m_sex["ASW", "PBS"] <- pars$m_sex["PBS", "ASW"] <- max(pars$m_gen_pop)
  
  # vaccination
  # some basic parameters so vaccination runs for all tests but then can develop further for vax specific ones in the tests themselves
  
  # daily doses
  pars$vaccination_campaign_length_children <- 10
  pars$daily_doses_children <- matrix(0,ncol=n_vax,
                                 nrow=pars$vaccination_campaign_length_children)
  pars$daily_doses_children[1:(pars$vaccination_campaign_length_children-1),
                       2] <- 1000
  
  pars$vaccination_campaign_length_adults <- 15
  pars$daily_doses_adults <- matrix(0,ncol=n_vax,
                               nrow=pars$vaccination_campaign_length_adults)
  pars$daily_doses_adults[1:(pars$vaccination_campaign_length_adults-5),2] <- 1000
  pars$daily_doses_adults[5:(pars$vaccination_campaign_length_adults-1),3] <- 1000

  pars$N_prioritisation_steps_children <- 3
  pars$N_prioritisation_steps_adults <- 2
  
  priority_children <- matrix(rep(pars$prioritisation_strategy_children,
                                  times=pars$N_prioritisation_steps_children),
                              ncol=pars$N_prioritisation_steps_children,
                              byrow=FALSE)
  priority_children[c(3,4,17),1] <- 0
  priority_children[c(4,17),2] <- 0
  
  pars$prioritisation_strategy_children <- priority_children * 0.5
  
  priority_adults <- matrix(rep(pars$prioritisation_strategy_adults,
                                times=pars$N_prioritisation_steps_adults),
                            ncol=pars$N_prioritisation_steps_adults,
                            byrow=FALSE)
  priority_adults[c(4:16,20),1] <- 0
  
  pars$prioritisation_strategy_adults <- priority_adults * 0.5

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
