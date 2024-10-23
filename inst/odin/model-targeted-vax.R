# State indices: [i: age; j: vaccination]

# Time Steps
dt <- user(1)
steps_per_week <- 7 / dt
initial(time) <- step
update(time) <- (step + 1) * dt
#output(time) <- TRUE

#### Vaccination

## j = 1: previous smallpox vaccine (efficacy roughly at one dose)
## j = 2: unvaccinated
## j = 3: 1 dose
## j = 4: 2 doses

## Emergency use of the MVA-BN has not been granted for children (u18s) in DRC and so we split vaccination into two components depending on age
## below are indicators for children vs adult groups
## use of _raw indicates that we are in prop form for the boundary category of 15 - 19
children_ind_raw[] <- user()
dim(children_ind_raw) <- c(n_group)
adults_ind_raw[] <- user()
dim(adults_ind_raw) <- c(n_group)

## specify the daily number of vaccinations that can happen at each time point
## first and last columns must be 0; last row must be 0 (pre-processing job)
## column j corresponds to the people in j receiving a vaccine (so moving
## to j + 1)
daily_doses_children[, ] <- user()
dim(daily_doses_children) <- c(vaccination_campaign_length_children, n_vax)

daily_doses_adults[, ] <- user()
dim(daily_doses_adults) <- c(vaccination_campaign_length_adults, n_vax)

## with this we need to ensure daily_doses for children and adults final time step is all 0 - done outside of the model in pre-processing
daily_doses_children_t[ ] <- if (as.integer(time) >= (vaccination_campaign_length_children))
  daily_doses_children[vaccination_campaign_length_children, i] else daily_doses_children[time, i]
dim(daily_doses_children_t) <- n_vax

daily_doses_adults_t[, ] <- if (as.integer(time) >= (vaccination_campaign_length_adults))
  daily_doses_adults[vaccination_campaign_length_adults, j] else daily_doses_adults[time, j]
dim(daily_doses_adults_t) <- c(1, n_vax)

## allocate the daily doses according to prioritisation strategy

## the number of different prioritisation strategies that we are considering
N_prioritisation_steps_children <- user()
N_prioritisation_steps_adults <- user()
## what this corresponds to in terms of groups is encoded here
## this is independent of dose / applies to both doses
prioritisation_strategy_children[, ] <- user()
dim(prioritisation_strategy_children) <- c(n_group, N_prioritisation_steps_children)

prioritisation_strategy_adults[, ] <- user()
dim(prioritisation_strategy_adults) <- c(n_group, N_prioritisation_steps_adults)

## vaccination targets: for each group, what is the level of coverage we want to
## see before we move onto the next set of groups to target?
## when we input this do the calc outside by the population size as part of
## pre-processing to save faffing around in here.
## here, j=1,...,n_vax corresponds to the coverage wanted in j, so for j=0
## and j=1 this stays at 0

## what is the current state of the vaccination targets? e.g. what
## prioritisation step are we on, this depends on whether we have met our
## targets

## has the target been met
## for cols 1 (smallpox vax) and 2 (unvaccinated) this doesn't matter for neither children or adults, for children col 4 (2nd dose) also doesn't matter
target_met_children_t[,] <- 0
dim(target_met_children_t) <- c(n_group, n_vax)

target_met_adults_t[, ] <- 0
dim(target_met_adults_t) <- c(n_group, n_vax)


## children
## 1st doses
## if you have a 2nd dose this implies you also have had a 1st dose so account
## for this in the 1st dose target - this now isn't strictly relevant for children but leaving it in in case we do expand this to 2 doses in future
target_met_children_t[, 3] <-
  if ((sum(N[i, 3:4])*children_ind_raw[i]) > # number of vaxxed kids in relevant compartments 
      round(prioritisation_strategy_children[i, prioritisation_step_1st_dose_children] * #children_ind_raw[i] *
            #vaccination_coverage_target_1st_dose_children_prop * 
            sum(N[i, 1:4])))
    1 else 0


## adults 
## 1st doses
## if you have a 2nd dose this implies you also have had a 1st dose so account
## for this in the 1st dose target
target_met_adults_t[, 3] <-
  if ((sum(N[i, 3:4])*adults_ind_raw[i])>
      round(prioritisation_strategy_adults[i, prioritisation_step_1st_dose_adults] * #adults_ind_raw[i] * 
            #vaccination_coverage_target_1st_dose_adults_prop * 
            sum(N[i, 1:4])))
    1 else 0


## 2nd doses
target_met_adults_t[, 4] <-
  if ((sum(N[i, 4])*adults_ind_raw[i]) >
      round(prioritisation_strategy_adults[i, prioritisation_step_2nd_dose_adults] * #adults_ind_raw[i] *
            #vaccination_coverage_target_2nd_dose_adults_prop * 
            sum(N[i, 1:4])))
    1 else 0

## prioritisation step proposal to account for the fact that this would update
## every single time step if we vaccinate quickly enough (unlikely but would
## like it to be done properly)

## the target goal for the current prioritisation strategy
## with ceiling basically saying if there is any target in that group to account for it
coverage_achieved_1st_dose_children[] <- ceiling(prioritisation_strategy_children[i, prioritisation_step_1st_dose_children])
dim(coverage_achieved_1st_dose_children) <- c(n_group)

coverage_achieved_1st_dose_adults[] <- ceiling(prioritisation_strategy_adults[i, prioritisation_step_1st_dose_adults])
dim(coverage_achieved_1st_dose_adults) <- c(n_group)

coverage_achieved_2nd_dose_adults[] <- ceiling(prioritisation_strategy_adults[i, prioritisation_step_2nd_dose_adults])
dim(coverage_achieved_2nd_dose_adults) <- c(n_group)


## children
prioritisation_step_1st_dose_children_proposal <-
  if (sum(target_met_children_t[, 3]) == sum(coverage_achieved_1st_dose_children[])) 
    prioritisation_step_1st_dose_children + 1 else
    prioritisation_step_1st_dose_children

update(prioritisation_step_1st_dose_children) <-
  if (prioritisation_step_1st_dose_children_proposal > N_prioritisation_steps_children)
    N_prioritisation_steps_children else prioritisation_step_1st_dose_children_proposal


## adults
prioritisation_step_1st_dose_adults_proposal <-
  if (sum(target_met_adults_t[, 3]) == sum(coverage_achieved_1st_dose_adults[]))
    prioritisation_step_1st_dose_adults + 1 else
    prioritisation_step_1st_dose_adults

update(prioritisation_step_1st_dose_adults) <-
  if (prioritisation_step_1st_dose_adults_proposal > N_prioritisation_steps_adults)
    N_prioritisation_steps_adults else prioritisation_step_1st_dose_adults_proposal

prioritisation_step_2nd_dose_adults_proposal <-
  if (sum(target_met_adults_t[, 4]) == sum(coverage_achieved_2nd_dose_adults[]))
    prioritisation_step_2nd_dose_adults + 1 else
    prioritisation_step_2nd_dose_adults

update(prioritisation_step_2nd_dose_adults) <-
  if (prioritisation_step_2nd_dose_adults_proposal > N_prioritisation_steps_adults)
    N_prioritisation_steps_adults else prioritisation_step_2nd_dose_adults_proposal


#prioritisation_step <- user()
#initial(prioritisation_step) <- 1
initial(prioritisation_step_1st_dose_children) <- 1
initial(prioritisation_step_1st_dose_adults) <- 1
initial(prioritisation_step_2nd_dose_adults) <- 1

## now that we know the step we are on, we know who is eligible per age group to
## be vaccinated


## who can get a children's dose
n_eligible_for_dose1_children[] <- (S[i, 2] + Ea[i, 2] + Eb[i, 2] + R[i, 2]) *
  prioritisation_strategy_children[i, prioritisation_step_1st_dose_children] 
dim(n_eligible_for_dose1_children) <- c(n_group)
## who can get a 1st dose
## now we have to look in the preceding j class (e.g. who in unvaccinated j = 1
## can get a 1st dose and move to j = 2) as these are the people who will be
## eligible
n_eligible_for_dose1_adults[] <- (S[i, 2] + Ea[i, 2] + Eb[i, 2] + R[i, 2]) *
  prioritisation_strategy_adults[i, prioritisation_step_1st_dose_adults] 
dim(n_eligible_for_dose1_adults) <- c(n_group)
## who can get a 2nd dose
n_eligible_for_dose2_adults[] <- (S[i, 3] + Ea[i, 3] + Eb[i, 3] + R[i, 3]) *
  prioritisation_strategy_adults[i, prioritisation_step_2nd_dose_adults] 
dim(n_eligible_for_dose2_adults) <- c(n_group)

## allocate the doses to the unvaccinated by age group, prioritisation strategy
## and across S, E, R
## hacky fix for now
n_vaccination_t_S[, ] <- 0
n_vaccination_t_Ea[, ] <- 0
n_vaccination_t_Eb[, ] <- 0
n_vaccination_t_R[, ] <- 0

n_vaccination_t_S_children[] <- 0
n_vaccination_t_Ea_children[] <- 0
n_vaccination_t_Eb_children[ ] <- 0
n_vaccination_t_R_children[ ] <- 0

n_vaccination_t_S_adults[] <- 0
n_vaccination_t_Ea_adults[ ] <- 0
n_vaccination_t_Eb_adults[ ] <- 0
n_vaccination_t_R_adults[ ] <- 0


### children 1st doses 
## allocate 1st doses
n_vaccination_t_S_children[] <- if (sum(n_eligible_for_dose1_children[]) == 0) 0 else
  min(floor((daily_doses_children_t[1, 2] * S[i, 2] *
               prioritisation_strategy_children[i, prioritisation_step_1st_dose_children]) / sum(n_eligible_for_dose1_children[])),
      S[i, 2])


n_vaccination_t_Ea_children[] <- if (sum(n_eligible_for_dose1_children[]) == 0) 0 else
  min(floor((daily_doses_children_t[1, 2] * Ea[i, 2] *
               prioritisation_strategy_children[i, prioritisation_step_1st_dose_children]) / sum(n_eligible_for_dose1_children[])),
      Ea[i, 2])


n_vaccination_t_Eb_children[] <- if (sum(n_eligible_for_dose1_children[]) == 0) 0 else
  min(floor((daily_doses_children_t[1, 2] * Eb[i, 2] *
               prioritisation_strategy_children[i, prioritisation_step_1st_dose_children]) / sum(n_eligible_for_dose1_children[])),
      Eb[i, 2])


n_vaccination_t_R_children[] <- if (sum(n_eligible_for_dose1_children[]) == 0) 0 else
  min(floor((daily_doses_children_t[1, 2] * R[i, 2] *
               prioritisation_strategy_children[i, prioritisation_step_1st_dose_children]) / sum(n_eligible_for_dose1_children[])),
      R[i, 2])


### adults

## allocate 1st doses
n_vaccination_t_S_adults[] <- if (sum(n_eligible_for_dose1_adults[]) == 0) 0 else
  min(floor((daily_doses_adults_t[1, 2] * S[i, 2] *
               prioritisation_strategy_adults[i, prioritisation_step_1st_dose_adults]) / sum(n_eligible_for_dose1_adults[])),
      S[i, 2])
  

n_vaccination_t_Ea_adults[] <- if (sum(n_eligible_for_dose1_adults[]) == 0) 0 else
  min(floor((daily_doses_adults_t[1, 2] * Ea[i, 2] *
               prioritisation_strategy_adults[i, prioritisation_step_1st_dose_adults]) / sum(n_eligible_for_dose1_adults[])),
      Ea[i, 2])
  

n_vaccination_t_Eb_adults[] <- if (sum(n_eligible_for_dose1_adults[]) == 0) 0 else
  min(floor((daily_doses_adults_t[1, 2] * Eb[i, 2] *
               prioritisation_strategy_adults[i, prioritisation_step_1st_dose_adults]) / sum(n_eligible_for_dose1_adults[])),
      Eb[i, 2])
  

n_vaccination_t_R_adults[] <- if (sum(n_eligible_for_dose1_adults[]) == 0) 0 else
  min(floor((daily_doses_adults_t[1, 2] * R[i, 2] *
               prioritisation_strategy_adults[i, prioritisation_step_1st_dose_adults] ) / sum(n_eligible_for_dose1_adults[])),
      R[i, 2] )


### combine total first doses

n_vaccination_t_S[,2] <- n_vaccination_t_S_children[i] + n_vaccination_t_S_adults[i]

n_vaccination_t_Ea[,2] <- n_vaccination_t_Ea_children[i] + n_vaccination_t_Ea_adults[i]

n_vaccination_t_Eb[,2] <- n_vaccination_t_Eb_children[i] + n_vaccination_t_Eb_adults[i]

n_vaccination_t_R[,2] <- n_vaccination_t_R_children[i] + n_vaccination_t_R_adults[i]

## for the boundary case do an extra check that we haven't gone over the number of people in each compartment 

n_vaccination_t_S[4,2] <- min(n_vaccination_t_S[4, 2], S[4, 2])

n_vaccination_t_Ea[4,2] <- if(n_vaccination_t_Ea[4,2]>Ea[4,2]) 
  Ea[4,2] else
    n_vaccination_t_Ea[4,2]

n_vaccination_t_Eb[4,2] <- if(n_vaccination_t_Eb[4,2]>Eb[4,2]) 
  Eb[4,2] else
    n_vaccination_t_Eb[4,2]

n_vaccination_t_R[4,2] <- if(n_vaccination_t_R[4,2]>R[4,2]) 
  R[4,2] else
    n_vaccination_t_R[4,2]


## allocate 2nd doses (adults only for now)
n_vaccination_t_S[, 3] <- if (sum(n_eligible_for_dose2_adults[]) == 0) 0 else
  min(floor((daily_doses_adults_t[1, 3] * S[i, 3] *
               prioritisation_strategy_adults[i, prioritisation_step_2nd_dose_adults]) / sum(n_eligible_for_dose2_adults[])),
      S[i, 3])


n_vaccination_t_Ea[, 3] <- if (sum(n_eligible_for_dose2_adults[]) == 0) 0 else
  min(floor((daily_doses_adults_t[1, 3] * Ea[i, 3] *
               prioritisation_strategy_adults[i, prioritisation_step_2nd_dose_adults]) / sum(n_eligible_for_dose2_adults[])),
      Ea[i, 3])


n_vaccination_t_Eb[, 3] <- if (sum(n_eligible_for_dose2_adults[]) == 0) 0 else
  min(floor((daily_doses_adults_t[1, 3] * Eb[i, 3] *
               prioritisation_strategy_adults[i, prioritisation_step_2nd_dose_adults]) / sum(n_eligible_for_dose2_adults[])),
      Eb[i, 3])


n_vaccination_t_R[, 3] <- if (sum(n_eligible_for_dose2_adults[]) == 0) 0 else
  min(floor((daily_doses_adults_t[1, 3] * R[i, 3] *
               prioritisation_strategy_adults[i, prioritisation_step_2nd_dose_adults]) / sum(n_eligible_for_dose2_adults[])),
      R[i,3])


## net vaccination change for relevant classes (S, Ea, Eb, R)
## logic here depends on vaccine class you are in (e.g. can only increase
## through j)

delta_S_n_vaccination[, ] <- if (j == 1) 0 else
  if (j == 2) (-n_vaccination_t_S[i, j]) else
    if (j == n_vax) (n_vaccination_t_S[i, j - 1]) else
      (-n_vaccination_t_S[i, j] + n_vaccination_t_S[i, j - 1])

delta_Ea_n_vaccination[, ] <- if (j == 1) 0 else
  if (j == 2) (-n_vaccination_t_Ea[i, j]) else
    if (j == n_vax) (n_vaccination_t_Ea[i, j - 1]) else
      (-n_vaccination_t_Ea[i, j] + n_vaccination_t_Ea[i, j - 1])

delta_Eb_n_vaccination[, ] <- if (j == 1) 0 else
  if (j == 2) (-n_vaccination_t_Eb[i, j]) else
    if (j == n_vax) (n_vaccination_t_Eb[i, j - 1]) else
      (-n_vaccination_t_Eb[i, j] + n_vaccination_t_Eb[i, j - 1])

delta_R_n_vaccination[, ] <- if (j == 1) 0 else
  if (j == 2) (-n_vaccination_t_R[i, j]) else
    if (j == n_vax) (n_vaccination_t_R[i, j - 1]) else
      (-n_vaccination_t_R[i, j] + n_vaccination_t_R[i, j - 1])

# use to test to see if vaccination working
update(vax_given_S) <- sum(n_vaccination_t_S[, ])
update(vax_given_Ea) <- sum(n_vaccination_t_Ea[, ])
update(vax_given_Eb) <- sum(n_vaccination_t_Eb[, ])
update(vax_given_R) <- sum(n_vaccination_t_R[, ])

# split by 1st and 2nd doses
update(vax_1stdose_given_S) <- sum(n_vaccination_t_S[, 2])
update(vax_1stdose_given_Ea) <- sum(n_vaccination_t_Ea[, 2])
update(vax_1stdose_given_Eb) <- sum(n_vaccination_t_Eb[, 2])
update(vax_1stdose_given_R) <- sum(n_vaccination_t_R[, 2])

update(vax_2nddose_given_S) <- sum(n_vaccination_t_S[, 3])
update(vax_2nddose_given_Ea) <- sum(n_vaccination_t_Ea[, 3])
update(vax_2nddose_given_Eb) <- sum(n_vaccination_t_Eb[, 3])
update(vax_2nddose_given_R) <- sum(n_vaccination_t_R[, 3])

## end of vaccination section


## Core equations for transitions between compartments:
# by age groups and vaccination class
# after vaccination has taken place
update(S[, ]) <- S[i, j] + delta_S_n_vaccination[i, j] - n_SEa[i, j]
update(Ea[, ]) <- Ea[i, j] + delta_Ea_n_vaccination[i, j] + delta_Ea[i, j]
update(Eb[, ]) <- Eb[i, j] + delta_Eb_n_vaccination[i, j] + delta_Eb[i, j]
update(Ir[, ]) <- Ir[i, j] + delta_Ir[i, j]
update(Id[, ]) <- Id[i, j] + delta_Id[i, j]
update(R[, ]) <- R[i, j] + delta_R_n_vaccination[i, j] + delta_R[i, j]
update(D[, ]) <- D[i, j] + delta_D[i, j]

## Additional outputs
update(E[, ]) <- Ea[i, j] + Eb[i, j]
update(I[, ]) <- Ir[i, j] + Id[i, j]
update(N[, ]) <- S[i, j] + Ea[i, j] + Eb[i, j] + Ir[i, j] + Id[i, j] + 
  R[i, j] +  D[i, j]
#update(N[, ]) <- S[i, j] + E[i, j] + I[i, j] + R[i, j] + D[i, j]

# weekly cases
# group indices 
# [1: 0-4,    2: 5-9,    3: 10-14,  4: 15-19,  5: 20-24,  6: 25-29,
#  7: 30-34,  8: 35-39,  9: 40-44, 10: 45-49, 11: 50-54, 12: 55-59,
# 13: 60-64, 14: 65-69, 15: 70-74, 16: 75+,   17: CSW,   18: ASW,
# 19: PBS,   20: HCW]
# Note that:
# X_05_14 includes 50% CSW (ages 12-14)
# X_15_plus includes all 50% CSW (ages 15-17), and all ASW, PBS, HCW

new_cases_00_04 <- sum(n_SEa[1, ])
new_cases_SW_12_14 <- rbinom(sum(n_SEa[17, ]), 0.5)
new_cases_SW_15_17 <- sum(n_SEa[17, ]) - new_cases_SW_12_14
new_cases_05_14 <- sum(n_SEa[2:3, ]) + new_cases_SW_12_14
new_cases_15_plus <-
  sum(n_SEa[4:16, ]) + new_cases_SW_15_17 + sum(n_SEa[18:20, ])
new_cases_SW <- sum(n_SEa[17:18, ])
new_cases_PBS <- sum(n_SEa[19, ])
new_cases_HCW <- sum(n_SEa[20, ])

is_same_week <- step %% steps_per_week > 0
update(cases_inc) <- cases_inc * is_same_week + sum(n_SEa[, ])
update(cases_inc_00_04) <- cases_inc_00_04 * is_same_week + new_cases_00_04
update(cases_inc_05_14) <- cases_inc_05_14 * is_same_week + new_cases_05_14
update(cases_inc_15_plus) <- cases_inc_15_plus * is_same_week +
  new_cases_15_plus
  
update(cases_inc_SW) <- cases_inc_SW * is_same_week + new_cases_SW
update(cases_inc_PBS) <- cases_inc_PBS * is_same_week + new_cases_PBS
update(cases_inc_HCW) <- cases_inc_SW * is_same_week + new_cases_HCW

# cumulative cases
update(cases_cumulative) <- cases_cumulative + sum(n_SEa[, ])
update(cases_cumulative_00_04) <- cases_cumulative_00_04 + new_cases_00_04
update(cases_cumulative_05_14) <- cases_cumulative_05_14 + new_cases_05_14
update(cases_cumulative_15_plus) <- cases_cumulative_15_plus + new_cases_15_plus
update(cases_cumulative_SW) <- cases_cumulative_SW + new_cases_SW
update(cases_cumulative_PBS) <- cases_cumulative_PBS + new_cases_PBS
update(cases_cumulative_HCW) <- cases_cumulative_HCW + new_cases_HCW


new_deaths_00_04 <- sum(n_IdD[1, ])
new_deaths_SW_12_14 <- rbinom(sum(n_IdD[17, ]), 0.5)
new_deaths_SW_15_17 <- sum(n_IdD[17, ]) - new_deaths_SW_12_14
new_deaths_05_14 <- sum(n_IdD[2:3, ]) + new_deaths_SW_12_14
new_deaths_15_plus <-
  sum(n_IdD[4:16, ]) + new_deaths_SW_15_17 + sum(n_IdD[18:20, ])
new_deaths_SW <- sum(n_IdD[17:18, ])
new_deaths_PBS <- sum(n_IdD[19, ])
new_deaths_HCW <- sum(n_IdD[20, ])
# weekly deaths
update(deaths_inc) <- deaths_inc * is_same_week + sum(n_IdD[, ])
update(deaths_inc_00_04) <- deaths_inc_00_04 * is_same_week + new_deaths_00_04
update(deaths_inc_05_14) <- deaths_inc_05_14 * is_same_week + new_deaths_05_14
update(deaths_inc_15_plus) <-
  deaths_inc_15_plus * is_same_week + new_deaths_15_plus

update(deaths_inc_SW) <- deaths_inc_SW * is_same_week + new_deaths_SW
update(deaths_inc_PBS) <- deaths_inc_PBS * is_same_week + new_deaths_PBS
update(deaths_inc_HCW) <- deaths_inc_HCW * is_same_week + new_deaths_HCW

# cumulative deaths
update(deaths_cumulative) <- deaths_cumulative + sum(n_IdD[, ])
update(deaths_cumulative_00_04) <- deaths_cumulative_00_04 + new_deaths_00_04
update(deaths_cumulative_05_14) <- deaths_cumulative_05_14 + new_deaths_05_14
update(deaths_cumulative_15_plus) <-
  deaths_cumulative_15_plus + new_deaths_15_plus

update(deaths_cumulative_SW) <- deaths_cumulative_SW + new_deaths_SW
update(deaths_cumulative_PBS) <- deaths_cumulative_PBS + new_deaths_PBS
update(deaths_cumulative_HCW) <- deaths_cumulative_HCW + new_deaths_HCW

update(S_tot) <- sum(S[, ])
update(E_tot) <- sum(E[, ])
update(I_tot) <- sum(I[, ])
update(R_tot) <- sum(R[, ])
update(D_tot) <- sum(D[, ])
update(N_tot) <- sum(N[, ])
update(total_vax) <- total_vax + vax_given_S + vax_given_Ea + vax_given_Eb +
  vax_given_R
update(total_vax_1stdose) <- total_vax_1stdose + vax_1stdose_given_S +
  vax_1stdose_given_Ea + vax_1stdose_given_Eb + vax_1stdose_given_R
update(total_vax_2nddose) <- total_vax_2nddose + vax_2nddose_given_S +
  vax_2nddose_given_Ea + vax_2nddose_given_Eb + vax_2nddose_given_R

## Individual probabilities of transition:
# S to E - age dependent
p_SE[, ] <- 1 - exp(-lambda[i, j] * dt)
# progression through latent period (2 subcompartments)
p_EE <- 1 - exp(-gamma_E * 2 * dt)
# progression to infection
p_EI <- 1 - exp(-gamma_E * 2 * dt)
# progression through infectious period to recovery
p_IrR <- 1 - exp(-gamma_Ir * dt)
# progression through infectious period to death
p_IdD <- 1 - exp(-gamma_Id * dt)

# Compute the force of infection

#  Mixing Matrix
m_gen_pop[, ] <- user()
m_sex[, ] <- user()
# I adjusted for reduced transmissibility
I_infectious[, ] <- I[i, j] * (1 - ve_T[j])

prop_infectious[] <-
  if (sum(N[i, ]) == 0) 0 else sum(I_infectious[i, ]) / sum(N[i, ])
# Generating Force of Infection
# for susceptible age i, % contacts infectious age j
s_ij_gen_pop[, ] <- m_gen_pop[i, j] * prop_infectious[j]
# as above but for the sexual contacts only
s_ij_sex[, ] <- m_sex[i, j] * prop_infectious[j]
lambda[, ] <- ((beta_h * sum(s_ij_gen_pop[i, ])) +
                 (beta_s * sum(s_ij_sex[i, ])) + beta_z[i]) * (1 - ve_I[i,j])

## Draws from binomial distributions for numbers changing between compartments
# accounting for vaccination:
n_SEa[, ] <- rbinom(S[i, j] + delta_S_n_vaccination[i, j], p_SE[i, j])
n_EaEb[, ] <- rbinom(Ea[i, j] + delta_Ea_n_vaccination[i, j], p_EE)
n_EbI[, ] <- rbinom(Eb[i, j] + delta_Eb_n_vaccination[i, j], p_EI)
# Proportion of the infections that will die, impact of vaccination included
# already as part of inputs
n_EbId[, ] <- rbinom(n_EbI[i, j], CFR[i, j])
# whatever infections don't die go to R
n_EbIr[, ] <- n_EbI[i, j] - n_EbId[i, j]

n_IrR[, ] <- rbinom(Ir[i, j], p_IrR)
n_IdD[, ] <- rbinom(Id[i, j], p_IdD)

## Calculate net change in each model state
delta_Ea[, ] <- n_SEa[i, j] - n_EaEb[i, j]
delta_Eb[, ] <- n_EaEb[i, j] - n_EbI[i, j]
delta_Ir[, ] <- n_EbIr[i, j] - n_IrR[i, j]
delta_Id[, ] <- n_EbId[i, j] - n_IdD[i, j]
delta_R[, ] <- n_IrR[i, j]
delta_D[, ] <- n_IdD[i, j]

## Initial states:
initial(S[, ]) <- S0[i, j]
initial(Ea[, ]) <- Ea0[i, j]
initial(Eb[, ]) <- Eb0[i, j]
initial(Ir[, ]) <- Ir0[i, j]
initial(Id[, ]) <- Id0[i, j]
initial(R[, ]) <- R0[i, j]
initial(D[, ]) <- D0[i, j]

initial(E[, ]) <- Ea0[i, j] + Eb0[i, j]
initial(I[, ]) <- Ir0[i, j] + Id0[i, j]
initial(N[, ]) <- S0[i, j] + Ea0[i, j] + Eb0[i, j] + Ir0[i, j] + Id0[i, j] +
  R0[i, j] + D0[i, j]
initial(cases_inc) <- 0
initial(deaths_inc) <- 0
initial(cases_cumulative) <- 0
initial(deaths_cumulative) <- 0

initial(cases_inc_00_04) <- 0
initial(cases_inc_05_14) <- 0
initial(cases_inc_15_plus) <- 0
initial(cases_inc_PBS) <- 0
initial(cases_inc_SW) <- 0
initial(cases_inc_HCW) <- 0
initial(deaths_inc_00_04) <- 0
initial(deaths_inc_05_14) <- 0
initial(deaths_inc_15_plus) <- 0
initial(deaths_inc_PBS) <- 0
initial(deaths_inc_SW) <- 0
initial(deaths_inc_HCW) <- 0

initial(cases_cumulative_00_04) <- 0
initial(cases_cumulative_05_14) <- 0
initial(cases_cumulative_15_plus) <- 0
initial(cases_cumulative_PBS) <- 0
initial(cases_cumulative_SW) <- 0
initial(cases_cumulative_HCW) <- 0
initial(deaths_cumulative_00_04) <- 0
initial(deaths_cumulative_05_14) <- 0
initial(deaths_cumulative_15_plus) <- 0
initial(deaths_cumulative_PBS) <- 0
initial(deaths_cumulative_SW) <- 0
initial(deaths_cumulative_HCW) <- 0

initial(vax_given_S) <- 0
initial(vax_given_Ea) <- 0
initial(vax_given_Eb) <- 0
initial(vax_given_R) <- 0

initial(vax_1stdose_given_S) <- 0
initial(vax_1stdose_given_Ea) <- 0
initial(vax_1stdose_given_Eb) <- 0
initial(vax_1stdose_given_R) <- 0

initial(vax_2nddose_given_S) <- 0
initial(vax_2nddose_given_Ea) <- 0
initial(vax_2nddose_given_Eb) <- 0
initial(vax_2nddose_given_R) <- 0

initial(S_tot) <- sum(S0[, ])
initial(E_tot) <- sum(Ea0[, ]) + sum(Eb0[, ])
initial(I_tot) <- sum(Ir0[, ]) + sum(Id0[, ])
initial(R_tot) <- sum(R0[, ])
initial(D_tot) <- sum(D0[, ])
initial(N_tot) <- sum(S0[, ]) + sum(Ea0[, ]) + sum(Eb0[, ]) + sum(Ir0[, ]) +
  sum(Id0[, ]) + sum(R0[, ]) + sum(D0[, ])
initial(total_vax) <- 0
initial(total_vax_1stdose) <- 0
initial(total_vax_2nddose) <- 0

##Initial vectors
S0[, ] <- user()
Ea0[, ] <- user()
Eb0[, ] <- user()
Ir0[, ] <- user()
Id0[, ] <- user()
R0[, ] <- user()
D0[, ] <- user()

##Parameters
beta_h <- user()
beta_s <- user()
beta_z[] <- user()
gamma_E <- user()
gamma_Ir <- user()
gamma_Id <- user()
CFR[, ] <- user()

#vaccine efficacy parameters
ve_T[] <- user()
ve_I[,] <- user()
#ve_D[,] <- user() # this is included within the CFR

#Number of age classes & number of transmissibility classes
n_vax <- user()
n_group <- user()

## Dimensions of the different "vectors" here vectors stand for
## multi-dimensional arrays
dim(N) <- c(n_group, n_vax)
dim(S) <- c(n_group, n_vax)
dim(S0) <- c(n_group, n_vax)
dim(p_SE) <- c(n_group, n_vax)
dim(n_SEa) <- c(n_group, n_vax)

dim(Ea) <- c(n_group, n_vax)
dim(Ea0) <- c(n_group, n_vax)
dim(Eb0) <- c(n_group, n_vax)
dim(delta_Ea) <- c(n_group, n_vax)
dim(n_EaEb) <- c(n_group, n_vax)

dim(Eb) <- c(n_group, n_vax)
dim(delta_Eb) <- c(n_group, n_vax)
dim(n_EbI) <- c(n_group, n_vax)

dim(n_EbId) <- c(n_group, n_vax)
dim(n_EbIr) <- c(n_group, n_vax)
dim(E) <- c(n_group, n_vax)

dim(Ir0) <- c(n_group, n_vax)
dim(Ir) <- c(n_group, n_vax)
dim(delta_Ir) <- c(n_group, n_vax)
dim(n_IrR) <- c(n_group, n_vax)

dim(Id0) <- c(n_group, n_vax)
dim(Id) <- c(n_group, n_vax)
dim(delta_Id) <- c(n_group, n_vax)
dim(n_IdD) <- c(n_group, n_vax)
dim(I) <- c(n_group, n_vax)

dim(R) <- c(n_group, n_vax)
dim(R0) <- c(n_group, n_vax)
dim(delta_R) <- c(n_group, n_vax)

dim(D) <- c(n_group, n_vax)
dim(D0) <- c(n_group, n_vax)
dim(delta_D) <- c(n_group, n_vax)

dim(lambda) <- c(n_group, n_vax)
dim(m_gen_pop) <- c(n_group, n_group)
dim(m_sex) <- c(n_group, n_group)
dim(I_infectious) <- c(n_group, n_vax)
dim(prop_infectious) <- c(n_group)
dim(s_ij_gen_pop) <- c(n_group, n_group)
dim(s_ij_sex) <- c(n_group, n_group)
dim(beta_z) <- c(n_group)

dim(CFR) <- c(n_group, n_vax)

dim(ve_T) <- c(n_vax)
dim(ve_I) <- c(n_group,n_vax)

vaccination_campaign_length_children <- user()
vaccination_campaign_length_adults <- user()

dim(n_vaccination_t_S) <- c(n_group, n_vax)
dim(n_vaccination_t_Ea) <- c(n_group, n_vax)
dim(n_vaccination_t_Eb) <- c(n_group, n_vax)
dim(n_vaccination_t_R) <- c(n_group, n_vax)

dim(n_vaccination_t_S_children) <- c(n_group)
dim(n_vaccination_t_Ea_children) <- c(n_group)
dim(n_vaccination_t_Eb_children) <- c(n_group)
dim(n_vaccination_t_R_children) <- c(n_group)

dim(n_vaccination_t_S_adults) <- c(n_group)
dim(n_vaccination_t_Ea_adults) <- c(n_group)
dim(n_vaccination_t_Eb_adults) <- c(n_group)
dim(n_vaccination_t_R_adults) <- c(n_group)

dim(delta_S_n_vaccination) <- c(n_group, n_vax)
dim(delta_Ea_n_vaccination) <- c(n_group, n_vax)
dim(delta_Eb_n_vaccination) <- c(n_group, n_vax)
dim(delta_R_n_vaccination) <- c(n_group, n_vax)


#### compare

exp_noise <- user(1e6)

# cases
cases <- data()
model_cases <- cases_inc + rexp(exp_noise)
compare(cases) ~ poisson(model_cases)

cases_00_04 <- data()
model_cases_00_04 <- cases_inc_00_04 + rexp(exp_noise)
compare(cases_00_04) ~ poisson(model_cases_00_04)

cases_05_14 <- data()
model_cases_05_14 <- cases_inc_05_14 + rexp(exp_noise)
compare(cases_05_14) ~ poisson(model_cases_05_14)

cases_15_plus <- data()
model_cases_15_plus <- cases_inc_15_plus + rexp(exp_noise)
compare(cases_15_plus) ~ poisson(model_cases_15_plus)

# deaths
deaths <- data()
model_deaths <- deaths_inc + rexp(exp_noise)
compare(deaths) ~ poisson(model_deaths)

deaths_00_04 <- data()
model_deaths_00_04 <- deaths_inc_00_04 + rexp(exp_noise)
compare(deaths_00_04) ~ poisson(model_deaths_00_04)

deaths_05_14 <- data()
model_deaths_05_14 <- deaths_inc_05_14 + rexp(exp_noise)
compare(deaths_05_14) ~ poisson(model_deaths_05_14)

deaths_15_plus <- data()
model_deaths_15_plus <- deaths_inc_15_plus + rexp(exp_noise)
compare(deaths_15_plus) ~ poisson(model_deaths_15_plus)
