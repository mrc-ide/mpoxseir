
#' @export
parameters_fixed <- function(overrides = list()) {

  ## Initialising variable that other parameters depend on
  demographic_params <- parameters_demographic()
  n_group <- demographic_params$n_group
  n_vax <- demographic_params$n_vax
  N0 <- demographic_params$N0
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
  params_list = list(
    n_group = n_group,
    n_vax = n_vax,
    N_prioritisation_steps = N_prioritisation_steps,
    Ea0 = Ea0,
    S0 = N0 - Ea0,
    Eb0 = X0,
    Ir0 = X0,
    Id0 = X0,
    R0 = X0,
    D0 = X0,
    R0_hh = 0.67, # Jezek 1988 SAR paper - will be fitted
    R0_sw_st = 1.3, # Will be fitted
    beta_z_max = 0.01, # Will be fitted
    RR_z = c(0.977, 1, 0.444, rep(0.078, n_group - 3)), # Jezek 1988 zoonotic + Jezek 1987
    gamma_E = 1 / 8,  # WHO Shiny - Kröger et al (clade II)
    gamma_Ir = 1 / 18, # Jezek 1988 "clinical features of 282.."
    gamma_Id = 1 / 10, # Jezek 1988
    CFR = matrix(CFR, nrow = n_group, ncol = n_vax, byrow = TRUE),
    m = demographic_params$m,
    # dt = 0.5,
    # runtime = 150,
    # particles = 1,
    # threads = 1,
    # seed = 42,
    # deterministic = TRUE,
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
transform_params <- function(
    beta_z_max,  ## zoonotic beta for the age-group with the highest risk (used to calculate RR_z)
    R0_hh,       ## R0 for household transmission (used to calculate beta_h)
    R0_sw_st    ## R0 for sex workers to people who buy sex
) {

  pars <- parameters_fixed()

  # Converting R0_hh to the beta_hh parameter given mixing matrix (excluding SW & PBS)

  ## Calc. duration of infectiousness by age, weighted by disease severity
  duration_infectious_by_age <-
    pars$CFR * (1 / pars$gamma_Id) + (1 - pars$CFR) * (1 / pars$gamma_Ir)

  ## Get indices of mixing matrix for general pop (those we assume hh transmission predominates)
  age_groups <- get_age_bins()
  idx_gen_pop <- seq_len(nrow(age_groups))

  ## Calculate beta_household given R0, mixing matrix and duration of infectiousness
  pars$beta_h <- R0_hh / Re(eigen(pars$m[idx_gen_pop, idx_gen_pop] * duration_infectious_by_age[idx_gen_pop])$values[1])

  # Converting the relative risk age-spline to the age-specific beta_z
  pars$beta_z <- pars$RR_z * pars$beta_z_max

  # Converting the inputted R0_sw_st into the rates required for the mixing matrix
  idx_SW <- which(colnames(m) == "SW")
  idx_PBS <- which(colnames(m) == "PBS")
  N_SW <- pars$N0[idx_SW]
  N_PBS <- pars$N0[idx_PBS]

  m_sw_pbs <- (R0_sw_st / duration_infectious_by_age[idx_SW]) / beta_h # multiplied by beta_h within in the model
  m_pbs_sw <- m_sw_pbs  * (N_SW / N_PBS) / pars$beta_h  # multiplied by beta_h within in the model
  pars$m[idx_SW, idx_PBS] <- m_sw_pbs
  pars$m[idx_PBS, idx_SW] <- m_pbs_sw

  pars
}

# Run the model with single set of parameter values
#' @export
run_mpoxSEIR_targetedVax_single <- function(
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

    ## Model simulation related parameters
    runtime = 150,
    n_particles = 1,
    n_threads = 1,
    seed = 42,
    deterministic = TRUE,
    outputs_retained = NULL,
) {

  ########################
  ### checks to go here
  ########################

  ## Transforming inputted parameters into those required for model running
  pars <- transform_params(beta_z_max = beta_z_max,
                                         R0_hh = R0_hh,
                                         R0_sw_st = R0_sw_st)

  ## Setting up the model with the inputted parameters
  mod <- mpoxseir:::model_targeted_vax$new(pars = pars,
                             time = 1,
                             n_particles = n_particles,
                             n_threads = n_threads,
                             seed = seed,
                             deterministic = deterministic)

  ## Running the model
  output <- mod$simulate(1:runtime)

  ## Postprocessing of model output - subsetting relevant outputs to retain
  indices <- mod$info()$idx
  if (any(!outputs_retained %in% names(indices))) {
    stop("outputs_retained must be in the named outputs of transformed model output")
  }
  indices_keep <- unlist(indices[which(names(indices) %in% outputs_retained)])
  output2 <- output[indices_keep, , ,drop = FALSE]

  ## Postprocessing of model output - creating tidy dataframe of outputs
  dimnames(output2) <- list(
    Output = names(indices_keep),
    Replicate = 1:particles,
    Timestep = dt * (1:runtime))
  grid <- list(state = dimnames(output2)[[1]] %||% seq_len(dim(output2)[1]),
               replicate = dimnames(output2)[[2]] %||% seq_len(dim(output2)[2]),
               time = dimnames(output2)[[3]] %||% seq_len(dim(output2)[3]))
  state <- output2
  ret <- do.call(expand.grid, grid)
  ret$time <- as.numeric(ret$time)
  ret$value <- c(state)

  ## Returning the output
  return(ret)

}

# Run multiple iterations of the model with different parameter values
#' @export
run_mpoxSEIR_targetedVax_multiple <- function(beta_z_max,
                                              R0_hh,
                                              R0_sw_st,
                                              ve_I,
                                              fixed_params) {

  # Check that the input vectors are the same length
  if (!(all.equal(length(beta_z_max), length(R0_hh), length(R0_sw_st), length(ve_I)))) {
    stop("beta_z_max, R0_hh, R0_sw_st and ve_I must be the same length")
  }

  # Use lapply to iterate over the indices of the input vectors
  output_list <- lapply(seq_along(beta_z_max), function(i) {

    # Extract the current values for this iteration
    current_beta_z_max <- beta_z_max[i]
    current_R0_hh <- R0_hh[i]
    current_R0_sw_st <- R0_sw_st[i]
    current_ve_I <- c(0, ve_I[i]) # 0 efficacy for unvaccinated

    # Run the single simulation function with the current parameters
    temp <- run_mpoxSEIR_targetedVax_single(n_group = fixed_params$n_group,
                                            S0 = fixed_params$S0,
                                            Ea0 = fixed_params$Ea0,
                                            Eb0 = fixed_params$Eb0,
                                            Ir0 = fixed_params$Ir0,
                                            Id0 = fixed_params$Id0,
                                            R0 = fixed_params$R0,
                                            D0 = fixed_params$D0,
                                            R0_hh = current_R0_hh,
                                            R0_sw_st = current_R0_sw_st,
                                            beta_z_max = current_beta_z_max,
                                            RR_z = fixed_params$RR_z,
                                            gamma_E = fixed_params$gamma_E,
                                            gamma_Ir = fixed_params$gamma_Ir,
                                            gamma_Id = fixed_params$gamma_Id,
                                            CFR = fixed_params$CFR,
                                            m = fixed_params$m,
                                            n_vax = fixed_params$n_vax,
                                            daily_doses = fixed_params$daily_doses,
                                            N_prioritisation_steps = fixed_params$N_prioritisation_steps,
                                            prioritisation_strategy = fixed_params$prioritisation_strategy,
                                            vaccination_coverage_target = fixed_params$vaccination_coverage_target,
                                            vaccine_uptake = fixed_params$vaccine_uptake,
                                            ve_T = fixed_params$ve_T,
                                            ve_I = current_ve_I,
                                            vaccination_campaign_length = fixed_params$vaccination_campaign_length,
                                            dt = fixed_params$dt,
                                            runtime = fixed_params$runtime,
                                            particles = fixed_params$particles,
                                            threads = fixed_params$threads,
                                            seed = fixed_params$seed,
                                            deterministic = fixed_params$deterministic,
                                            outputs_retained = fixed_params$outputs_retained)

    # Appending the parameters used in this particular simulation to the output
    temp$parameter_set_idx <- i
    temp$beta_z_max <- current_beta_z_max
    temp$R0_sw_st <- current_R0_sw_st
    temp$R0_hh <- current_R0_hh
    temp$ve_Is <- current_ve_I
    temp
  })

  return(data.table::rbindlist(output_list))
}
