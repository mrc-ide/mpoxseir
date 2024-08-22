library(mpoxseir); library(ggplot2)

par <- default_params()
transformed_params <- transform_params(S0 = par$S0, Ea0 = par$Ea0, Eb0 = par$Eb0, 
                                       Ir0 = par$Ir0, Id0 = par$Id0, R0 = par$R0, D0 = par$D0,
                                       RR_z = par$RR_z, beta_z_max = par$beta_z_max, 
                                       R0_hh = par$R0_hh, R0_sw_st = par$R0_sw_st,
                                       m = par$m,  
                                       CFR = par$CFR, 
                                       gamma_Ir = par$gamma_Ir, gamma_Id = par$gamma_Id)
single_model_output <- run_mpoxSEIR_targetedVax_single(n_group = par$n_group,    
                                                       S0 = par$S0,         
                                                       Ea0 = par$Ea0,        
                                                       Eb0 = par$Eb0,        
                                                       Ir0 = par$Ir0,        
                                                       Id0 = par$Id0,        
                                                       R0 = par$R0,         
                                                       D0 = par$D0,         
                                                       R0_hh = par$R0_hh,      
                                                       R0_sw_st = par$R0_sw_st,   
                                                       beta_z_max = par$beta_z_max, 
                                                       RR_z = par$RR_z,       
                                                       gamma_E = par$gamma_E,    
                                                       gamma_Ir = par$gamma_Ir,   
                                                       gamma_Id = par$gamma_Id,  
                                                       CFR = par$CFR,        
                                                       m = par$m,          
                                                       n_vax = par$n_vax,                       
                                                       daily_doses = par$daily_doses,                 
                                                       N_prioritisation_steps = par$N_prioritisation_steps,      
                                                       prioritisation_strategy = par$prioritisation_strategy,      
                                                       vaccination_coverage_target = par$vaccination_coverage_target, 
                                                       vaccine_uptake = par$vaccine_uptake,              
                                                       ve_T = par$ve_T,                        
                                                       ve_I = par$ve_I,                        
                                                       vaccination_campaign_length = par$vaccination_campaign_length, 
                                                       dt = par$dt, 
                                                       runtime = par$runtime, 
                                                       particles = 10, #par$particles,
                                                       threads = par$threads,
                                                       seed = par$seed,
                                                       deterministic = FALSE, # par$deterministic,
                                                       outputs_retained = c("cases", "deaths"))

ggplot(single_model_output, aes(x = time, y = value, col = factor(replicate))) +
  geom_line() +
  facet_wrap(.~state, scales = "free_y")

iterations <- 10
beta_z_max_vector <- rep(0.1, iterations)
R0_hh_vector <- rep(0.75, iterations)
R0_sw_st_vector <- rep(1.75, iterations)
ve_I_vector <- rep(0.8, iterations)
fixed_params <- par
fixed_params$outputs_retained <- c("cases", "deaths")
fixed_params$particles <- 5
multiple_model_outputs <- run_mpoxSEIR_targetedVax_multiple(beta_z_max = beta_z_max_vector,
                                                            R0_hh = R0_hh_vector,
                                                            R0_sw_st = R0_sw_st_vector,
                                                            ve_I = ve_I_vector,
                                                            fixed_params = fixed_params)

ggplot(multiple_model_outputs, aes(x = time, y = value, col = factor(replicate))) +
  geom_line() +
  facet_grid(parameter_set_index~state, scales = "free_y")


# par <- default_params()
# n_group = par$n_group    
# S0 = par$S0         
# Ea0 = par$Ea0        
# Eb0 = par$Eb0        
# Ir0 = par$Ir0        
# Id0 = par$Id0        
# R0 = par$R0         
# D0 = par$D0         
# R0_hh = par$R0_hh      
# R0_sw_st = par$R0_sw_st   
# beta_z_max = par$beta_z_max 
# RR_z = par$RR_z       
# gamma_E = par$gamma_E    
# gamma_Ir = par$gamma_Ir   
# gamma_Id = par$gamma_Id  
# CFR = par$CFR        
# m = par$m          
# n_vax = par$n_vax                       
# daily_doses = par$daily_doses                 
# N_prioritisation_steps = par$N_prioritisation_steps      
# prioritisation_strategy = par$prioritisation_strategy      
# vaccination_coverage_target = par$vaccination_coverage_target 
# vaccine_uptake = par$vaccine_uptake              
# ve_T = par$ve_T                        
# ve_I = par$ve_I                        
# vaccination_campaign_length = par$vaccination_campaign_length 
# dt = par$dt 
# runtime = par$runtime 
# particles = par$particles
# threads = par$threads
# seed = par$seed
# deterministic = par$deterministic
# outputs_retained <- c("cases", "deaths")
