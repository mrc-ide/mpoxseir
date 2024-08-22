
transform_params <- function(
    
  
) {
  
  # Generating parameters list for model running (ret) and updating this with parameter values derived from inputted pars (e.g. R0)
  pars <- as.list(pars)
  ret <- fixed_pars()
  
  # 1) Converting R0_hh to the beta_hh parameter given mixing matrix (excluding SW & PBS)
  duration_infectious_by_age <- (ret$CFR * (1 / ret$gamma_Id)) + ((1 - ret$CFR) * (1 / ret$gamma_Ir)) ## Calc. duration of infectiousness by age, weighted by disease severity
  m_dim <- dim(ret$m)[1]                                                                              
  index_gen_pop <- 1:(m_dim-2) ## Get indices of mixing matrix for general pop (those we assume hh transmission predominates)
  ret$beta_h <- pars$R0_hh / Re(eigen(ret$m[index_gen_pop, index_gen_pop] * duration_infectious_by_age[index_gen_pop])$values[1]) ## Calculate beta_household given R0, mixing matrix and duration of infectiousness
  
  # 2) Converting the relative risk age-spline to the age-specific beta_z
  ret$beta_z <- ret$RR_z * pars$beta_z_max
  
  # 3) Converting the inputted R0_sw_st into the rates required for the mixing matrix
  model_compartments <- c("S0", "Ea0", "Eb0", "Ir0", "Id0", "R0", "D0")
  index_SW <- which(colnames(ret$m) == "SW")
  index_PBS <- which(colnames(ret$m) == "PBS") 
  
  N_SW <- sum(sapply(ret[model_compartments], function(x) x[index_SW]))
  N_PBS <- sum(sapply(ret[model_compartments], function(x) x[index_PBS]))
  
  beta_sw_pbs <- (pars$R0_sw_st / duration_infectious_by_age[index_SW]) / ret$beta_h
  beta_pbs_sw <- ((pars$R0_sw_st * N_SW) / (duration_infectious_by_age[index_SW] * N_PBS)) / ret$beta_h
  ret$m[index_SW, index_PBS] <- beta_sw_pbs
  ret$m[index_PBS, index_SW] <- beta_pbs_sw
  
}

#' @importFrom stats dmultinom
#' @importFrom squire get_mixing_matrix
#' 
#' @export model-targeted-vax
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
    beta_h,     ## beta for household transmission (a number)
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
    dt = 1, 
    runtime = 150, 
    particles = 1,
    threads = 1,
    seed = 42,
    deterministic = TRUE
) {
  
  ########################
  ### checks to go here
  ########################
  
  ## Setting up the model with the inputted parameters
  beta_z <- beta_z_max * RR_z
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
  output2 <- mod$transform_variables(output)
  
  ## Returning the output
  return(output2)
  
}



run_mpoxSEIR_targetedVax_multi <- function()

params <- parameters_demographic()
n_group <- params$n_group
n_vax <- 2
Ea0 <- matrix(rep(rep(5, n_group), n_vax), ncol = n_vax)
S0 <- matrix(params$N - Ea0, ncol = n_vax)   
Eb0 <- matrix(rep(rep(0, n_group), n_vax), ncol = n_vax)
Ir0 <- matrix(rep(rep(round(255 * 0.9/ n_group), n_group), n_vax), ncol = n_vax)      
Id0 <- matrix(rep(rep(round(255 * 0.1/ n_group), n_group), n_vax), ncol = n_vax)   
R0 <- matrix(rep(rep(0, n_group), n_vax), ncol = n_vax)       
D0 <- matrix(rep(rep(0, n_group), n_vax), ncol = n_vax)       
beta_h <- 0.2 / 12.11   
beta_z_max <- 0.01
RR_z <- c(0.875, 1, 0.514, rep(0.101, n_group - 3))      
gamma_E <- 1 / 8  
gamma_Ir <- 1 / 4 
gamma_Id <- 1 / 4
CFR <- matrix(rep(c(0.105, rep(0.048, 2), rep(0.03, n_group - 3)), n_vax), ncol = n_vax)      
m <- params$m         
dt <- 0.5 
runtime <- 150 
particles <- 1
threads <- 1
seed <- 42
deterministic <- TRUE
N_prioritisation_steps <- 1
prioritisation_strategy <- matrix(rep(1, n_group), ncol = N_prioritisation_steps)  
vaccination_coverage_target <- matrix(rep(0.01, n_group), ncol = N_prioritisation_steps)  
vaccine_uptake <- rep(0.8, n_group)
ve_T <- rep(0, n_vax)
ve_I <- rep(0, n_vax)
vaccination_campaign_length <- 1
daily_doses <- matrix(rep(1, vaccination_campaign_length * n_vax),
                      nrow = vaccination_campaign_length, ncol = n_vax)
