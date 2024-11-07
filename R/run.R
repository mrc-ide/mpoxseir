
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
