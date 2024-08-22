## NOTE THAT THESE NEED UPDATING IN A MAJOR WAY
#' @export
default_params <- function(overrides = list()) {
  
  ## Initialising variable that other parameters depend on
  demographic_params <- parameters_demographic()
  n_group <- demographic_params$n_group
  n_vax <- demographic_params$n_vax
  N0 <- demographic_params$N0
  N_prioritisation_steps <- 1
  Ea0 <- matrix(rep(rep(5, n_group), n_vax), ncol = n_vax)
  vaccination_campaign_length <- 1
    
  ## note this needs to be updated with proper assignment of people into unvax vs vax
  params_list = list(n_group = n_group,
                     n_vax = n_vax,
                     N_prioritisation_steps = N_prioritisation_steps,
                     Ea0 = Ea0,
                     S0 = matrix(N0 - Ea0, nrow = n_group, ncol = n_vax),   
                     Eb0 = matrix(rep(rep(0, n_group), n_vax), ncol = n_vax),
                     Ir0 = matrix(rep(rep(round(255 * 0.9/ n_group), n_group), n_vax), ncol = n_vax),      
                     Id0 = matrix(rep(rep(round(255 * 0.1/ n_group), n_group), n_vax), ncol = n_vax),   
                     R0 = matrix(rep(rep(0, n_group), n_vax), ncol = n_vax),       
                     D0 = matrix(rep(rep(0, n_group), n_vax), ncol = n_vax),       
                     R0_hh = 0.75,
                     R0_sw_st = 2,
                     beta_z_max = 0.01,
                     RR_z = c(0.875, 1, 0.514, rep(0.101, n_group - 3)),      
                     gamma_E = 1 / 8,  
                     gamma_Ir = 1 / 4, 
                     gamma_Id = 1 / 4,
                     CFR = matrix(rep(c(0.105, rep(0.048, 2), rep(0.03, n_group - 3)), n_vax), ncol = n_vax),      
                     m = demographic_params$m,         
                     dt = 0.5, 
                     runtime = 150, 
                     particles = 1,
                     threads = 1,
                     seed = 42,
                     deterministic = TRUE,
                     prioritisation_strategy = matrix(rep(1, n_group), ncol = N_prioritisation_steps),  
                     vaccination_coverage_target = matrix(rep(0.01, n_group), ncol = N_prioritisation_steps),  
                     vaccine_uptake = rep(0.8, n_group),
                     ve_T = rep(0, n_vax),
                     ve_I = rep(0, n_vax),
                     vaccination_campaign_length = vaccination_campaign_length,
                     daily_doses = matrix(rep(1, vaccination_campaign_length * n_vax),
                                          nrow = vaccination_campaign_length, ncol = n_vax))
  
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
  S0,          ## number of initial Susceptible individuals (vector of length n_group)
  Ea0,         ## number of initial Exposed individuals in first E compartment (vector of length n_group)
  Eb0,         ## number of initial Exposed individuals in second E compartment  (vector of length n_group)
  Ir0,         ## number of initial Infectious individuals who will go on to recover (vector of length n_group) 
  Id0,         ## number of initial Infectious individuals who will go on to die (vector of length n_group) 
  R0,          ## number of initial Recovered individuals (vector of length n_group) 
  D0,          ## number of initial Dead individuals (vector of length n_group) 
  RR_z,        ## relative rate scaling factor for zoonotic transmission in other groups (vector of length n_group)
  beta_z_max,  ## zoonotic beta for the age-group with the highest risk (used to calculate RR_z)
  R0_hh,       ## R0 for household transmission (used to calculate beta_h)
  R0_sw_st,    ## R0 for sex workers to people who buy sex
  m,           ## mixing matrix
  CFR,         ## CFR
  gamma_Ir,    ## gamma_Id
  gamma_Id     ## gamma_Ir
) {

  
  # Converting R0_hh to the beta_hh parameter given mixing matrix (excluding SW & PBS)
  if (dim(CFR)[2] > 1) {
    CFR <- CFR[, 1]
  }
  duration_infectious_by_age <- (CFR * (1 / gamma_Id)) + ((1 - CFR) * (1 / gamma_Ir)) ## Calc. duration of infectiousness by age, weighted by disease severity
  m_dim <- dim(m)[1]                                                                              
  index_gen_pop <- 1:(m_dim-2) ## Get indices of mixing matrix for general pop (those we assume hh transmission predominates)
  beta_h <- R0_hh / Re(eigen(m[index_gen_pop, index_gen_pop] * duration_infectious_by_age[index_gen_pop])$values[1]) ## Calculate beta_household given R0, mixing matrix and duration of infectiousness
  
  # Converting the relative risk age-spline to the age-specific beta_z
  beta_z <- RR_z * beta_z_max
  
  # Converting the inputted R0_sw_st into the rates required for the mixing matrix
  index_SW <- which(colnames(m) == "SW")
  index_PBS <- which(colnames(m) == "PBS") 
  N_SW <- S0[index_SW] + Ea0[index_SW] + Eb0[index_SW] + Ir0[index_SW] + Id0[index_SW] + R0[index_SW] + D0[index_SW]
  N_PBS <- S0[index_PBS] + Ea0[index_PBS] + Eb0[index_PBS] + Ir0[index_PBS] + Id0[index_PBS] + R0[index_PBS] + D0[index_PBS]
  
  beta_sw_pbs <- (R0_sw_st / duration_infectious_by_age[index_SW]) / beta_h # multiplied by beta_h within in the model
  beta_pbs_sw <- ((R0_sw_st * N_SW) / (duration_infectious_by_age[index_SW] * N_PBS)) / beta_h  # multiplied by beta_h within in the model
  m[index_SW, index_PBS] <- beta_sw_pbs
  m[index_PBS, index_SW] <- beta_pbs_sw
  
  return(list(beta_h = beta_h,
              beta_z = beta_z,
              m = m))
}

# Run the model with single set of parameter values
#' @export
run_mpoxSEIR_targetedVax_single <- function(
    
    ## Number of groups and initial states
    n_group,    ## number of groups (integer, includes age-groups, SWs and PBS)
    S0,         ## number of initial Susceptible individuals (vector of length n_group)
    Ea0,        ## number of initial Exposed individuals in first E compartment (vector of length n_group)
    Eb0,        ## number of initial Exposed individuals in second E compartment  (vector of length n_group)
    Ir0,        ## number of initial Infectious individuals who will go on to recover (vector of length n_group) 
    Id0,        ## number of initial Infectious individuals who will go on to die (vector of length n_group) 
    R0,         ## number of initial Recovered individuals (vector of length n_group) 
    D0,         ## number of initial Dead individuals (vector of length n_group) 
  
    ## Transmission related parameters
    R0_hh,      ## R0 for the household
    R0_sw_st,   ## R0 for sex workers to people who buy sex
    beta_z_max, ## beta for the age-group with highest zoonotic transmission (a number)
    RR_z,       ## relative rate scaling factor for zoonotic transmission in other groups (vector of length n_group)
    gamma_E,    ## rate of transition from E->I (a number, 1/gamma_E = duration of incubation period) 
    gamma_Ir,   ## rate of transition from I->R (a number, 1/gamma_Ir = duration of infectiousness in those who will go on to recover) 
    gamma_Id,   ## rate of transition from I->D (a number, 1/gamma_Id = duration of infectiousness in those who will go on to die) 
    CFR,        ## case fatality ratio (matrix of dim n_group * n_vax)
    m,          ## mixing matrix of contact rates per day (matrix of dim n_group * n_group)
    
    ## Vaccination related parameters
    n_vax,                       ## number of vaccination compartments (integer, basis for an additional dimension in odin states)
    daily_doses,                 ## the daily number of doses administered (matrix of vaccination_campaign_length * number of vaccination compartments)
    N_prioritisation_steps,      ## the number of different vaccination prioritisation categories we're considering
    prioritisation_strategy,     ## what each step corresponds to in terms of strategy (matrix of n_group * N_prioritisation_steps, with 1s and 0s indicating whether a group is included in prioritisation step)
    vaccination_coverage_target, ## vacccination coverage target for each group and prioritisation step (matrix of n_group * N_prioritisation_steps)
    vaccine_uptake,              ## max achievable coverage for each group (vector of length n_group) CHECK WITH RUTH THIS IS RIGHT
    ve_T,                        ## vaccine efficacy against onwards transmissibility for each vaccinated compartment (vector of length n_vax)
    ve_I,                        ## vaccine efficacy against infection for each vaccinated compartment (vector of length n_vax)
    vaccination_campaign_length, ## length of the vaccination campaign (in timesteps NOT days - CHECK WITH RUTH THIS IS RIGHT)
    
    ## Model simulation related parameters
    dt = 0.5, 
    runtime = 150, 
    particles = 1,
    threads = 1,
    seed = 42,
    deterministic = TRUE,
    outputs_retained
) {
  
  ########################
  ### checks to go here
  ########################
  
  ## Transforming inputted parameters into those required for model running
  transformed_params <- transform_params(S0 = S0, Ea0 = Ea0, Eb0 = Eb0, Ir0 = Ir0, Id0 = Id0, R0 = R0, D0 = D0,
                                         RR_z = RR_z, beta_z_max = beta_z_max, 
                                         R0_hh = R0_hh,  
                                         R0_sw_st = R0_sw_st,
                                         m = m, 
                                         CFR = CFR, 
                                         gamma_Ir = gamma_Ir, gamma_Id = gamma_Id)
  beta_h <- transformed_params$beta_h
  beta_z <- transformed_params$beta_z
  m <- transformed_params$m
  
  ## Setting up the model with the inputted parameters
  mod <- mpoxseir:::model_targeted_vax$new(pars = list(n_group = n_group,    
                                                       S0 = S0, Ea0 = Ea0, Eb0 = Eb0, Ir0 = Ir0, Id0 = Id0, R0 = R0, D0 = D0,         
                                                       beta_h = beta_h, beta_z = beta_z,       
                                                       gamma_E = gamma_E, gamma_Ir = gamma_Ir, gamma_Id = gamma_Id,
                                                       CFR = CFR, m = m,
                                                       n_vax = n_vax, 
                                                       N_prioritisation_steps = N_prioritisation_steps, prioritisation_strategy = prioritisation_strategy,
                                                       vaccination_coverage_target = vaccination_coverage_target, vaccine_uptake = vaccine_uptake,
                                                       vaccination_campaign_length = vaccination_campaign_length, daily_doses = daily_doses,
                                                       ve_T = ve_T, ve_I = ve_I,
                                                       dt = dt), 
                             time = 1,
                             n_particles = particles,
                             n_threads = threads,
                             seed = seed,
                             deterministic = deterministic)  
  
  ## Running the model
  output <- mod$simulate(1:runtime)
  
  ## Postprocessing of model output - subsetting relevant outputs to retain 
  indices <- mod$info()$index
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
    temp$parameter_set_index <- i
    temp$beta_z_max <- current_beta_z_max
    temp$R0_sw_st <- current_R0_sw_st
    temp$R0_hh <- current_R0_hh
    temp$ve_Is <- current_ve_I
    temp
  })
  
  return(data.table::rbindlist(output_list))
}