##' A function that creates the transform function for use in the fitting
##' 
##' @title Create transform function
##' @inheritParams parameters_fixed
##' 
##' @return A transform function
##'   
##' @export
##'
create_transform_params <- function(region, initial_infections,
                                    use_ve_D = FALSE, overrides = list()) {

  pars <- parameters_fixed(region = region,
                           initial_infections = initial_infections,
                           use_ve_D = use_ve_D,
                           overrides = overrides)

  function(transformed_pars) {
    ## fitted params are:
    ## beta_z_max = zoonotic beta for the age-group with the highest risk (used to calculate RR_z)
    ## R0_hh = R0 for household transmission (used to calculate beta_h)
    ## R0_sw_st = R0 for sex workers to people who buy sex)
    nms_group <- names(pars$N0)

    ## Updating values in pars with parameters passed into transform_params
    pars$beta_z_max <- transformed_pars["beta_z_max"]
    pars$R0_hh <- transformed_pars["R0_hh"]
    pars$R0_sw_st <- transformed_pars["R0_sw_st"]
    
    pars$alpha_cases <- transformed_pars["alpha_cases"]
    pars$alpha_cases_00_04 <- pars$alpha_cases
    pars$alpha_cases_05_14 <- pars$alpha_cases
    pars$alpha_cases_15_plus <- pars$alpha_cases
    
    pars$alpha_deaths <- transformed_pars["alpha_deaths"]
    pars$alpha_deaths_00_04 <- pars$alpha_deaths
    pars$alpha_deaths_05_14 <- pars$alpha_deaths
    pars$alpha_deaths_15_plus <- pars$alpha_deaths

    # Converting R0_hh to the beta_hh parameter given mixing matrix (excluding SW & PBS)

    ## Calc. duration of infectiousness by age, weighted by disease severity
    ## We assume here that the R0 calculation that duration_infectious_by_age is used in
    ## is based on duration_infectious_by_age in unvaccinated individuals (i.e. with the unvaxxed CFR)
    k <- 1 # number of compartments per infectious disease state
    dt <- 1 # timestep
    idx <- get_compartment_indices()
    
    if (pars$n_vax > 1) {
      CFR <- pars$CFR[, idx$vax$unvaccinated]
    } else {
     CFR <- pars$CFR
    }

    duration_infectious_by_age <-
      CFR * ((k * dt) / (1 - exp(-pars$gamma_Id * dt))) +
      (1 - CFR) * ((k * dt) / (1 - exp(-pars$gamma_Ir * dt)))
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
    # we never resolve the product m[i, j] = beta_s * c[i, j] here, so beta_s is
    # set to 1. I.e. we are working with transmission rates, rather than transmission
    # probabilities + contact rates.

    N_CSW <- pars$N0["CSW"]
    N_ASW <- pars$N0["ASW"]
    N_PBS <- pars$N0["PBS"]

    Delta_PBS <- duration_infectious_by_age[idx$group$PBS]
    
    # transmission rate from CSW / ASW to PBS is the same
    m_csw_pbs <- pars$R0_sw_st / Delta_PBS * N_PBS / (N_CSW + N_ASW)
    m_asw_pbs <- m_csw_pbs
    # transmission from PBS -> CSW / ASW is proportional to their group sizes
    m_pbs_csw <- m_csw_pbs * N_CSW / N_PBS
    m_pbs_asw <- m_csw_pbs * N_ASW / N_PBS

    pars$m_sex[idx$group$CSW, idx$group$PBS] <- m_csw_pbs
    pars$m_sex[idx$group$ASW, idx$group$PBS] <- m_asw_pbs
    pars$m_sex[idx$group$PBS, idx$group$CSW] <- m_pbs_csw
    pars$m_sex[idx$group$PBS, idx$group$ASW] <- m_pbs_asw

    pars$beta_s <- 1

    pars

  }

}

##' Run the model with single set of parameter values
##'
##' @title Run the model with single set of parameter values
##' 
##' @param region The region to run the model for, must be either `"equateur"` 
##'   or `"sudkivu"`
##'   
##' @param initial_infections The initial number of infections to seed with
##' 
##' @param n_weeks number of weeks to run for
##' 
##' @param R0_hh R0 for the household
##' 
##' @param R0_sw_st R0 for sex workers to people who buy sex
##' 
##' @param beta_z_max beta for the age-group with highest zoonotic 
##'   transmission (a number)
##' 
##' @param n_vax number of vaccination compartments (integer, basis for an 
##'   additional dimension in odin states)
##' 
##' @param daily_doses the daily number of doses administered (matrix of 
##'   vaccination_campaign_length * number of vaccination compartments)
##' 
##' @param N_prioritisation_steps the number of different vaccination 
##'   prioritisation categories we're considering
##' 
##' @param prioritisation_strategy what each step corresponds to in terms of 
##'   strategy (matrix of n_group * N_prioritisation_steps, with 1s and 0s
##'   indicating whether a group is included in prioritisation step)
##' 
##' @param vaccination_coverage_target vaccination coverage target for each 
##'   group and prioritisation step (matrix of n_group * N_prioritisation_steps)
##' 
##' @param vaccine_uptake max achievable coverage for each group (vector of 
##'   length n_group)
##' 
##' @param ve_T vaccine efficacy against onwards transmissibility for each 
##'   vaccinated compartment (vector of length n_vax)
##' 
##' @param ve_I vaccine efficacy against infection for each vaccinated 
##'   compartment (vector of length n_vax)
##' 
##' @param ve_D vaccine efficacy against death for each vaccinated 
##'   compartment (vector of length n_vax)
##' 
##' @param vaccination_campaign_length length of the vaccination campaign
##'   (in timesteps NOT days - CHECK WITH RUTH THIS IS RIGHT)
##' 
##' @param overrides list of other model parameters which if specified will
##'   overwrite the defaults
##' 
##' @param n_particles Number of particles
##' 
##' @param n_threads Number of threads
##' 
##' @param seed The random seed to use
##' 
##' @param deterministic Logical, whether to run the model deterministically or
##'   not
##' 
##' @param outputs_retained The outputs to retain
##' 
##' @return The output from running the model
##'   
##' @export
##'
run_mpoxSEIR_targetedVax_single <- function(

    region,  # region to run the model for (either Sud Kivu or Equateur)
    initial_infections, # number of initial infections to seed with
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
  transform_params <- create_transform_params(region = region, initial_infections = initial_infections, overrides = list())
  params_for_transform <- c("beta_z_max" = beta_z_max, "R0_hh" = R0_hh, "R0_sw_st" = R0_sw_st)
  pars <- transform_params(params_for_transform)

  ## Setting up the model with the inputted parameters
  mod <- model_targeted_vax$new(pars = pars,
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

##' Run multiple iterations of the model with different parameter values
##'
##' @title Run multiple iterations of the model with different parameter values
##' 
##' @param region The region to run the model for, must be either `"equateur"` 
##'   or `"sudkivu"`
##'   
##' @param initial_infections The initial number of infections to seed with
##' 
##' @param n_weeks number of weeks to run for
##' 
##' @param R0_hh R0 for the household
##' 
##' @param R0_sw_st R0 for sex workers to people who buy sex
##' 
##' @param beta_z_max beta for the age-group with highest zoonotic 
##'   transmission (a number)
##' 
##' @param n_vax number of vaccination compartments (integer, basis for an 
##'   additional dimension in odin states)
##' 
##' @param daily_doses the daily number of doses administered (matrix of 
##'   vaccination_campaign_length * number of vaccination compartments)
##' 
##' @param N_prioritisation_steps the number of different vaccination 
##'   prioritisation categories we're considering
##' 
##' @param prioritisation_strategy what each step corresponds to in terms of 
##'   strategy (matrix of n_group * N_prioritisation_steps, with 1s and 0s
##'   indicating whether a group is included in prioritisation step)
##' 
##' @param vaccination_coverage_target vaccination coverage target for each 
##'   group and prioritisation step (matrix of n_group * N_prioritisation_steps)
##' 
##' @param vaccine_uptake max achievable coverage for each group (vector of 
##'   length n_group)
##' 
##' @param ve_T vaccine efficacy against onwards transmissibility for each 
##'   vaccinated compartment (vector of length n_vax)
##' 
##' @param ve_I vaccine efficacy against infection for each vaccinated 
##'   compartment (vector of length n_vax)
##' 
##' @param ve_D vaccine efficacy against death for each vaccinated 
##'   compartment (vector of length n_vax)
##' 
##' @param vaccination_campaign_length length of the vaccination campaign
##'   (in timesteps NOT days - CHECK WITH RUTH THIS IS RIGHT)
##' 
##' @param overrides list of other model parameters which if specified will
##'   overwrite the defaults
##' 
##' @param n_particles Number of particles
##' 
##' @param n_threads Number of threads
##' 
##' @param seed The random seed to use
##' 
##' @param deterministic Logical, whether to run the model deterministically or
##'   not
##' 
##' @param outputs_retained The outputs to retain
##' 
##' @return The output from running the model
##'   
##' @export
##'
run_mpoxSEIR_targetedVax_multiple <- function(
    region,
    initial_infections, 
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
                   region = region,                     ## region to run the model for
                   initial_infections = initial_infections, 
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
