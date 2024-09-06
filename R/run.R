
#' @export
parameters_fixed <- function(N, overrides = list()) {

  ## Initialising variable that other parameters depend on
  demographic_params <- parameters_demographic()
  n_group <- demographic_params$n_group
  n_vax <- demographic_params$n_vax
  N0 <- round(N * demographic_params$N0 / sum(demographic_params$N0)) # total number in each age-group
  N_prioritisation_steps <- 1

  # seed infections
  Ea0 <- matrix(5, nrow = n_group, ncol = n_vax)
  X0 <- matrix(0, nrow = n_group, ncol = n_vax)

  # CFR from Whittles 2024, 5-year bands to 40
  age_bins <- get_age_bins()
  CFR <- rep(0, n_group)
  names(CFR) <- names(demographic_params$N0)
  CFR[which(age_bins$end < 40)] <- c(0.102, 0.054, 0.035, 0.026, 0.02, 0.016, 0.013, 0.012)
  CFR[which(age_bins$start >= 40)] <- 0.01
  CFR["SW"] <- CFR["20-24"]
  CFR["PBS"] <- CFR["35-39"]

  vaccination_campaign_length <- 1

  ## note this needs to be updated with proper assignment of people into unvax vs vax
  ## NOTE THIS REALLY NEEDS TO BE UPDATED WITH PROPER ASSIGNMENT OF PEOPLE INTO UNVAX VS VAX
  params_list = list(
    n_group = n_group,
    n_vax = n_vax,
    N_prioritisation_steps = N_prioritisation_steps,
    S0 = round(N0 / n_vax) - Ea0,
    Ea0 = Ea0,
    Eb0 = X0,
    Ir0 = X0,
    Id0 = X0,
    R0 = X0,
    D0 = X0,
    N0 = N0,
    R0_hh = 0.67, # Jezek 1988 SAR paper - will be fitted
    R0_sw_st = 1.3, # Will be fitted
    beta_z_max = 0.01, # Will be fitted
    RR_z = c(0.977, 1, 0.444, rep(0.078, n_group - 3)), # Jezek 1988 zoonotic + Jezek 1987
    gamma_E = 1 / 8,  # WHO Shiny - Kröger et al (clade II) (RM comment: propose 1/7 based on Besombes et al. on 29 clade I patients)
    gamma_Ir = 1 / 18, # Jezek 1988 "clinical features of 282.."
    gamma_Id = 1 / 10, # Jezek 1988
    CFR = matrix(CFR, nrow = n_group, ncol = n_vax, byrow = FALSE),
    m_sex = demographic_params$m_sex,
    m_gen_pop = demographic_params$m_gen_pop,
    prioritisation_strategy = matrix(1, nrow = n_group, ncol = N_prioritisation_steps),
    vaccination_coverage_target = matrix(0.01, nrow = n_group, ncol = N_prioritisation_steps),
    vaccine_uptake = rep(0.8, n_group),
    ve_T = rep(0, n_vax),
    ve_I = rep(0, n_vax),
    vaccination_campaign_length = vaccination_campaign_length,
    daily_doses = matrix(1, nrow = vaccination_campaign_length, ncol = n_vax))

  # Ensure overridden parameters are passed as a list
  if (!is.list(overrides)) {
    stop('overrides must be a list')
  }

  # Override parameter values in the overrides input
  for (name in names(overrides)) {
    if (!(name %in% names(params_list))) {
      stop(paste('unknown parameter', name, sep=' '))
    }
    params_list[[name]] <- overrides[[name]]
  }

  return(params_list)
}

#' @export
create_transform_params <- function(N, overrides = list()) {
  
  pars <- parameters_fixed(N = N, overrides = overrides)
  
  function(beta_z_max,  ## zoonotic beta for the age-group with the highest risk (used to calculate RR_z)
           R0_hh,       ## R0 for household transmission (used to calculate beta_h)
           R0_sw_st) {  ## R0 for sex workers to people who buy sex)
    
    nms_group <- names(pars$N0)
    
    ## Updating values in pars with parameters passed into transform_params 
    pars$beta_z_max <- beta_z_max
    pars$R0_hh <- R0_hh
    pars$R0_sw_st <- R0_sw_st
    
    # Converting R0_hh to the beta_hh parameter given mixing matrix (excluding SW & PBS)
    
    ## Calc. duration of infectiousness by age, weighted by disease severity
    ## We assume here that the R0 calculation that duration_infectious_by_age is used in
    ## is based on duration_infectious_by_age in unvaccinated individuals (i.e. with the unvaxxed CFR)
    if (pars$n_vax > 1) {
      duration_infectious_by_age <-
        pars$CFR[, 1] * (1 / pars$gamma_Id) + (1 - pars$CFR[, 1]) * (1 / pars$gamma_Ir)
    } else{
      duration_infectious_by_age <-
        pars$CFR * (1 / pars$gamma_Id) + (1 - pars$CFR) * (1 / pars$gamma_Ir)
    }
    
    ## Get indices of mixing matrix for general pop (those we assume hh transmission predominates)
    age_groups <- get_age_bins()
    idx_gen_pop <- seq_len(nrow(age_groups))
    
    ## Calculate beta_household given R0, mixing matrix and duration of infectiousness
    ## Note this isn't quite correct and we need to sort this
    ## We should be include household contacts of PBS and SWs
    pars$beta_h <- pars$R0_hh / Re(eigen(pars$m_gen_pop * duration_infectious_by_age)$values[1])
    
    # Converting the relative risk age-spline to the age-specific beta_z
    pars$beta_z <- pars$RR_z * pars$beta_z_max
    
    # Converting the inputted R0_sw_st into the rates required for the mixing matrix
    idx_SW <- which(nms_group == "SW")
    idx_PBS <- which(nms_group == "PBS")
    N_SW <- pars$N0[idx_SW]
    N_PBS <- pars$N0[idx_PBS]
    
    m_sw_pbs <- (pars$R0_sw_st / duration_infectious_by_age[idx_SW])
    m_pbs_sw <- m_sw_pbs  * (N_SW / N_PBS)
    pars$m_sex[idx_SW, idx_PBS] <- m_sw_pbs
    pars$m_sex[idx_PBS, idx_SW] <- m_pbs_sw
    pars$beta_s <- 1
    
    pars
    
  }   
  
}

# Run the model with single set of parameter values
#' @export
run_mpoxSEIR_targetedVax_single <- function(
  
    N,       # population size to run the model with
    n_weeks, # number of weeks to run for
    
    ## Transmission related parameters
    R0_hh,      ## R0 for the household
    R0_sw_st,   ## R0 for sex workers to people who buy sex
    beta_z_max, ## beta for the age-group with highest zoonotic transmission (a number)

    ## Vaccination related parameters
    n_vax = 0,                       ## number of vaccination compartments (integer, basis for an additional dimension in odin states)
    daily_doses = 0,                 ## the daily number of doses administered (matrix of vaccination_campaign_length * number of vaccination compartments)
    N_prioritisation_steps = 0,      ## the number of different vaccination prioritisation categories we're considering
    prioritisation_strategy = 0,     ## what each step corresponds to in terms of strategy (matrix of n_group * N_prioritisation_steps, with 1s and 0s indicating whether a group is included in prioritisation step)
    vaccination_coverage_target = 0, ## vacccination coverage target for each group and prioritisation step (matrix of n_group * N_prioritisation_steps)
    vaccine_uptake = 0,              ## max achievable coverage for each group (vector of length n_group) CHECK WITH RUTH THIS IS RIGHT
    ve_T = 0,                        ## vaccine efficacy against onwards transmissibility for each vaccinated compartment (vector of length n_vax)
    ve_I = 0,                        ## vaccine efficacy against infection for each vaccinated compartment (vector of length n_vax)
    ve_D = 0,
    vaccination_campaign_length = 0, ## length of the vaccination campaign (in timesteps NOT days - CHECK WITH RUTH THIS IS RIGHT)

    ## Other model parameters
    overrides = list(),              ## list of other model parameters which if specified will overwrite the defaults - IN PROGRESS
    
    ## Model simulation related parameters
    n_particles = 1,
    n_threads = 1,
    seed = 42,
    deterministic = TRUE,
    outputs_retained = NULL) {

  ########################
  ### checks to go here
  ########################

  ## Transforming inputted parameters into those required for model running
  transform_params <- create_transform_params(N = N, overrides = list())
  pars <- transform_params(beta_z_max = beta_z_max,
                           R0_hh = R0_hh,
                           R0_sw_st = R0_sw_st)

  ## Setting up the model with the inputted parameters
  mod <- mpoxseir:::model_targeted_vax$new(pars = pars,
                                           time = 0,
                                           n_particles = n_particles,
                                           n_threads = n_threads,
                                           seed = seed,
                                           deterministic = deterministic)

  ## Running the model
  days_per_week <- 7
  output_times <- c(0, seq_len(n_weeks) * days_per_week)
  output <- mod$simulate(output_times)

  ## Postprocessing of model output - subsetting relevant outputs to retain
  indices <- mod$info()$index
  if (any(!outputs_retained %in% names(indices))) {
    stop("outputs_retained must be in the named outputs of transformed model output")
  }
  if (is.null(outputs_retained)) {
    outputs_retained <- names(indices)
  }
  indices_keep <- unlist(indices[which(names(indices) %in% outputs_retained)])

  ## Postprocessing of model output - creating tidy dataframe of outputs
  grid <- list(state = outputs_retained,
               replicate = seq_len(n_particles),
               time = output_times)

  ret <- do.call(expand.grid, grid)
  ret$value <- c(output[indices_keep, , ])

  ret
}

# Run multiple iterations of the model with different parameter values
#' @export
run_mpoxSEIR_targetedVax_multiple <- function(
    N,
    n_weeks,
    beta_z_max,
    R0_hh,
    R0_sw_st,
    ## Vaccination related parameters
    n_vax = 0,                       ## number of vaccination compartments (integer, basis for an additional dimension in odin states)
    daily_doses = 0,                 ## the daily number of doses administered (matrix of vaccination_campaign_length * number of vaccination compartments)
    N_prioritisation_steps = 0,      ## the number of different vaccination prioritisation categories we're considering
    prioritisation_strategy = 0,     ## what each step corresponds to in terms of strategy (matrix of n_group * N_prioritisation_steps, with 1s and 0s indicating whether a group is included in prioritisation step)
    vaccination_coverage_target = 0, ## vacccination coverage target for each group and prioritisation step (matrix of n_group * N_prioritisation_steps)
    vaccine_uptake = 0,              ## max achievable coverage for each group (vector of length n_group) CHECK WITH RUTH THIS IS RIGHT
    ve_T = 0,                        ## vaccine efficacy against onwards transmissibility for each vaccinated compartment (vector of length n_vax)
    ve_I = 0,                        ## vaccine efficacy against infection for each vaccinated compartment (vector of length n_vax)
    ve_D = 0,
    vaccination_campaign_length = 0, ## length of the vaccination campaign (in timesteps NOT days - CHECK WITH RUTH THIS IS RIGHT)
    overrides = list(),              ## list of other model parameters which if specified will overwrite the defaults - IN PROGRESS
    
    ## Model simulation related parameters
    n_particles = 1,
    n_threads = 1,
    seed = 42,
    deterministic = TRUE,
    outputs_retained = NULL) {

      # Check that the input vectors are the same length
      if (!(all.equal(length(beta_z_max), length(R0_hh), length(R0_sw_st), length(ve_I)))) {
        stop("beta_z_max, R0_hh, R0_sw_st and ve_I must be the same length")
      }

      # Use lapply to iterate over the indices of the input vectors
      ret <- Map(run_mpoxSEIR_targetedVax_single,
                 beta_z_max = beta_z_max,
                 R0_hh = R0_hh,
                 R0_sw_st = R0_sw_st,
                 MoreArgs = list(
                   N = N,                               ## size of the population to run the model with
                   n_weeks = n_weeks,
                   n_vax = n_vax,                       ## number of vaccination compartments (integer, basis for an additional dimension in odin states)
                   daily_doses = daily_doses,                 ## the daily number of doses administered (matrix of vaccination_campaign_length * number of vaccination compartments)
                   N_prioritisation_steps = N_prioritisation_steps,      ## the number of different vaccination prioritisation categories we're considering
                   prioritisation_strategy = prioritisation_strategy,     ## what each step corresponds to in terms of strategy (matrix of n_group * N_prioritisation_steps, with 1s and 0s indicating whether a group is included in prioritisation step)
                   vaccination_coverage_target = vaccination_coverage_target, ## vacccination coverage target for each group and prioritisation step (matrix of n_group * N_prioritisation_steps)
                   vaccine_uptake = vaccine_uptake,              ## max achievable coverage for each group (vector of length n_group) CHECK WITH RUTH THIS IS RIGHT
                   ve_T = ve_T,                        ## vaccine efficacy against onwards transmissibility for each vaccinated compartment (vector of length n_vax)
                   ve_I = ve_I,                        ## vaccine efficacy against infection for each vaccinated compartment (vector of length n_vax)
                   ve_D = ve_D,
                   vaccination_campaign_length = vaccination_campaign_length, ## length of the vaccination campaign (in timesteps NOT days - CHECK WITH RUTH THIS IS RIGHT)
                   overrides = overrides,              ## list of other model parameters which if specified will overwrite the defaults - IN PROGRESS
                   
                   ## Model simulation related parameters
                   n_particles = n_particles,
                   n_threads = n_threads,
                   seed = seed,
                   deterministic = deterministic,
                   outputs_retained = outputs_retained))

      dplyr::bind_rows(ret, .id = "parameter_set")
}
