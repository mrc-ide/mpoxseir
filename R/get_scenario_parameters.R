## vaccine scenario parameters
### The parameters in here need to be updated/pulled from the other functions so use with caution! but the scenario information is included in here

get_scenario_parameters <- function(location = c("equateur",
                                        "sudkivu"),
                           vaccine_dose_scenario = c("1st_doses_only",
                                                     "2nd_doses_only"),
                           total_doses = 2000000,
                           vaccine_uptake_scenario = c("all_uptake",
                                                       "uptake_by_group"),
                           vaccine_prioritisation_scenario = c("kids_1st",
                                                               "kids_CSW_1st",
                                                               "all_equal"),
                           doses_per_day,
                           vaccine_coverage_target_prop = 0.8){

  if(!is.integer(doses_per_day)|!is.integer(total_doses)){
    stop("doses_per_day and total_doses must both be integers")
  }

  if(vaccine_coverage_target_prop<=0|vaccine_coverage_target_prop>1){
    stop("vaccine_coverage_target_prop must be greater than 0 and less than or equal to 1")
  }

  dem_pars <- parameters_demographic()
  n_group <- dem_pars$n_group
  n_vax <- dem_pars$n_vax

  ## location

  if(location=="equateur"){

    N <- round(dem_pars$age_prop*dem_pars$province_pop$equateur)
    beta_z <- rep(1/100,times=dem_pars$n_group) # complete nonsense for now
    doses_to_give <- rbinom(n=1,size=total_doses,
                            prob=dem_pars$province_pop$equateur/sum(unlist(dem_pars$province_pop)))

  } else if(location=="sudkivu"){

    N <- round(dem_pars$age_prop*dem_pars$province_pop$sudkivu)
    beta_z <- rep(1/1000,times=dem_pars$province_pop$sudkivu) # complete nonsense for now
    doses_to_give <- rbinom(n=1,size=total_doses,
                            prob=dem_pars$province_pop$sudkivu/sum(unlist(dem_pars$province_pop)))

  } else{
    stop("location not valid")
  }

  ## vaccine prioritisation scenario
  ## must be a matrix
  if(vaccine_prioritisation_scenario=="kids_1st"){
    N_prioritisation_steps <- 3
    prioritisation_strategy <- cbind(c(1,1,1,rep(0,n_group-3)),# kids
                                     c(1,1,1,rep(0,n_group-5),1,1), # kids + CSW
                                     c(rep(1,n_group))) # all
  } else if(vaccine_prioritisation_scenario=="kids_CSW_1st"){
    N_prioritisation_steps <- 2
    prioritisation_strategy <- cbind(c(1,1,1,rep(0,n_group-5),1,1), # kids + CSW
                                     c(rep(1,n_group))) # all
  } else if(vaccine_prioritisation_scenario=="all_equal"){
    N_prioritisation_steps <- 1
    prioritisation_strategy <- matrix((rep(1,n_group)),ncol=1) # all
  } else{
    stop("vaccine_prioritisation_scenario not valid")
  }

  if(!is.matrix(prioritisation_strategy)){
    stop("prioritisation_strategy must be a matrix")
  }

  ## vaccine targets (which determine progression through prioritisation scenario)
  vaccine_coverage_target <- round(prioritisation_strategy * vaccine_coverage_target_prop * dem_pars$N0)

  if(!is.matrix(vaccine_coverage_target)){
    stop("vaccine_coverage_target must be a matrix")
  }

  ## vaccine uptake
  ## uncertainty on this if we need (see spreadsheet)
  if(vaccine_uptake_scenario=="all_uptake"){
    vaccine_uptake <- rep(1,n_group)
  } else if(vaccine_uptake_scenario=="uptake_by_group"&region=="equateur"){
    vaccine_uptake <- c(rep(1,3),rep(0.478,n_group-3))
  } else if(vaccine_uptake_scenario=="uptake_by_group"&region=="sudkivu"){
    vaccine_uptake <- c(rep(1,3),rep(0.558,n_group-3))
  } else{
    stop("vaccine_uptake_scenario or location not valid")
  }


  ## set up the vaccination strategy
  vaccination_campaign_length_proposal <- floor(doses_to_give/doses_per_day)

  daily_doses <- matrix(0,
                        nrow=vaccination_campaign_length_proposal+2,
                        ncol=n_vax)
  daily_doses[1:vaccination_campaign_length_proposal,1] <- doses_per_day
  # allocate whatever is leftover
  daily_doses[vaccination_campaign_length_proposal+1,1] <- (doses_to_give - sum(daily_doses[1:vaccination_campaign_length_proposal,1]))
  ## the last row must be 0

  if(sum(daily_doses)!=doses_to_give){
    warning("total vaccines allocated does not equal total vaccines planned")
  }

  ## vaccine dosage scenarios

  if(vaccine_dose_scenario=="1st_doses_only"){
    ve_I <- c(1,0.736)
  } else if(vaccine_dose_scenario=="2nd_doses_only"){
    ve_I <- c(1,0.818)
    daily_doses <- round(daily_doses * 0.5)
  } else{
    stop("vaccine_dose_scenario not valid")
  }


  # number of unvaccinated susceptibles given historic smallpox vaccination
  S0 <- matrix(0,nrow=dem_pars$n_group,ncol=dem_pars$n_vax)
  S0[,1] <- rbinom(size=N,prob=dem_pars$sus_prop,n=length(N))
  S0[,n_vax] <- N-S0[,1]


  scenario <- list(
    location = location,
    vaccine_dose_scenario = vaccine_dose_scenario,
    vaccine_uptake_scenario = vaccine_uptake_scenario,
    vaccine_prioritisation_scenario = vaccine_prioritisation_scenario,
    doses_per_day = doses_per_day,
    total_doses = total_doses,
    total_doses_for_location = doses_to_give,
    vaccine_coverage_target_prop = vaccine_coverage_target_prop)

  params <- list(
    dt = 1,
    n_group = n_group,
    n_vax = n_vax,
    N = dem_pars$N,
    CFR = matrix(c(0.119,0.082,0.06,0.047,0.039,0.033,0.028,0.025,
                   rep(0.025,10),
                   0.048,0.033,0.024,0.019,0.015,0.013,0.011,0.01,
                   rep(0.01,10)),
                 nrow=dem_pars$n_group,ncol=dem_pars$n_vax),
    ve_I = ve_I,
    ve_T = c(1,1),#dummy for now
    m = dem_pars$m,
    gamma_E = 0.5*1/7,
    gamma_I = 0.5*1/7,
    gamma_Ir = 1/19,
    gamma_Id = 1/10,
    beta_h = 1/1000,
    beta_z = beta_z,
    S0 = S0,
    Ea0 = matrix(0, dem_pars$n_group, dem_pars$n_vax),
    Eb0 = matrix(0, dem_pars$n_group, dem_pars$n_vax),
    Ir0 = matrix(0, dem_pars$n_group, dem_pars$n_vax),
    Id0 = matrix(0, dem_pars$n_group, dem_pars$n_vax),
    R0 = matrix(0, dem_pars$n_group, dem_pars$n_vax),
    D0 = matrix(0, dem_pars$n_group, dem_pars$n_vax),
    vaccination_campaign_length = nrow(daily_doses),
    daily_doses = daily_doses,
    N_prioritisation_steps = N_prioritisation_steps,
    prioritisation_strategy = prioritisation_strategy,
    vaccine_coverage_target = vaccine_coverage_target,
    vaccine_uptake = vaccine_uptake
  )

  return(list(scenario = scenario,
              params = params))

}


