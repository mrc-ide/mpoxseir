# State indices: [i: age; j: vaccination]

#### Vaccination

## j = 1: previous smallpox vaccine (efficacy roughly at one dose)
## j = 2: unvaccinated
## j = 3: 1 dose
## j = 4: 2 doses

## Emergency use of the MVA-BN has not been granted for children (u18s) in DRC
## and so we split vaccination into two components depending on age
## below are indicators for children vs adult groups
## use of _raw indicates that we are in prop form for the boundary category
## of 15 - 19
children_ind_raw <- parameter()
dim(children_ind_raw) <- c(n_group)
adults_ind_raw <- parameter()
dim(adults_ind_raw) <- c(n_group)

### setup the number of daily doses
##    daily_doses_value has dimensions n_vax x n_time
##    daily_doses_time has dimensions n_time
##    daily_doses_t then has dimensions n_vax
## we use interpolate such that 
##    daily_doses_t[i] = daily_doses_value[i, j] when 
##      daily_doses_time[j] <= time < daily_doses_time[j + 1]
daily_doses_children_value <- parameter()
dim(daily_doses_children_value) <- parameter(rank = 2)
daily_doses_children_time <- parameter()
dim(daily_doses_children_time) <- parameter(rank = 1)

daily_doses_adults_value <- parameter()
dim(daily_doses_adults_value) <- parameter(rank = 2)
daily_doses_adults_time <- parameter()
dim(daily_doses_adults_time) <- parameter(rank = 1)

daily_doses_children_t <- 
  interpolate(daily_doses_children_time, daily_doses_children_value, "constant")
dim(daily_doses_children_t) <- n_vax

daily_doses_adults_t <- 
  interpolate(daily_doses_adults_time, daily_doses_adults_value, "constant")
dim(daily_doses_adults_t) <- n_vax

## allocate the daily doses according to prioritisation strategy

## the number of different prioritisation strategies that we are considering
N_prioritisation_steps_children <- parameter()
N_prioritisation_steps_adults <- parameter()
## what this corresponds to in terms of groups is encoded here
## this is independent of dose / applies to both doses
prioritisation_strategy_children <- parameter()
dim(prioritisation_strategy_children) <-
  c(n_group, N_prioritisation_steps_children)

prioritisation_strategy_adults <- parameter()
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
## for cols 1 (smallpox vax) and 2 (unvaccinated) this doesn't matter for
## neither children or adults, for children col 4 (2nd dose) also doesn't matter

## children
dim(target_met_children_t) <- c(n_group, n_vax)

## 1st doses
## if you have a 2nd dose this implies you also have had a 1st dose so account
## for this in the 1st dose target - this now isn't strictly relevant for
## children but leaving it in in case we do expand this to 2 doses in future
target_met_children_t[, ] <- 0
target_met_children_t[, 3] <-
  ((sum(N[i, 3:4]) * children_ind_raw[i]) >
     prioritisation_strategy_children[
       i, prioritisation_step_1st_dose_children] * sum(N[i, ]))

## children_ind_raw -> is_child 
## is_adult <- 1 - is_child 

## adults
dim(target_met_adults_t) <- c(n_group, n_vax)

## 1st doses
## if you have a 2nd dose this implies you also have had a 1st dose so account
## for this in the 1st dose target
target_met_adults_t[, ] <- 0
target_met_adults_t[, 3] <-
  ((sum(N[i, 3:4]) * adults_ind_raw[i]) >
     prioritisation_strategy_adults[i, prioritisation_step_1st_dose_adults] *
     sum(N[i, ]))


## 2nd doses
target_met_adults_t[, 4] <-
  ((sum(N[i, 4]) * adults_ind_raw[i]) >
     prioritisation_strategy_adults[i, prioritisation_step_2nd_dose_adults] *
     sum(N[i, ]))


## prioritisation step proposal to account for the fact that this would update
## every single time step if we vaccinate quickly enough (unlikely but would
## like it to be done properly)

## the target goal for the current prioritisation strategy
## with ceiling basically saying if there is any target in that group to account
## for it
coverage_achieved_1st_dose_children[] <- ceiling(
  prioritisation_strategy_children[i, prioritisation_step_1st_dose_children])
dim(coverage_achieved_1st_dose_children) <- c(n_group)

coverage_achieved_1st_dose_adults[] <- ceiling(
  prioritisation_strategy_adults[i, prioritisation_step_1st_dose_adults])
dim(coverage_achieved_1st_dose_adults) <- c(n_group)

coverage_achieved_2nd_dose_adults[] <- ceiling(
  prioritisation_strategy_adults[i, prioritisation_step_2nd_dose_adults])
dim(coverage_achieved_2nd_dose_adults) <- c(n_group)

## children
prioritisation_step_1st_dose_children_proposal <-
  if (sum(target_met_children_t[, 3]) ==
      sum(coverage_achieved_1st_dose_children[]))
    prioritisation_step_1st_dose_children + 1 else
    prioritisation_step_1st_dose_children

update(prioritisation_step_1st_dose_children) <-
  if (prioritisation_step_1st_dose_children_proposal >
      N_prioritisation_steps_children)
    N_prioritisation_steps_children else
      prioritisation_step_1st_dose_children_proposal


## adults
prioritisation_step_1st_dose_adults_proposal <-
  if (sum(target_met_adults_t[, 3]) == sum(coverage_achieved_1st_dose_adults[]))
    prioritisation_step_1st_dose_adults + 1 else
    prioritisation_step_1st_dose_adults

update(prioritisation_step_1st_dose_adults) <-
  if (prioritisation_step_1st_dose_adults_proposal >
      N_prioritisation_steps_adults)
    N_prioritisation_steps_adults else
      prioritisation_step_1st_dose_adults_proposal

prioritisation_step_2nd_dose_adults_proposal <-
  if (sum(target_met_adults_t[, 4]) == sum(coverage_achieved_2nd_dose_adults[]))
    prioritisation_step_2nd_dose_adults + 1 else
    prioritisation_step_2nd_dose_adults

update(prioritisation_step_2nd_dose_adults) <-
  if (prioritisation_step_2nd_dose_adults_proposal >
      N_prioritisation_steps_adults)
    N_prioritisation_steps_adults else
      prioritisation_step_2nd_dose_adults_proposal


#prioritisation_step <- parameter()
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

### allocate to S

### if target_met == 1 then vax is 0; 
## if n_eligb == 0 then 0;
## else (daily_doses * S / n_elig) * ceiling(prioiritisation) (remove prioritisation)

## give_dose_X X = {dose_1_children,dose_1_adult,dose_2_adult} length n_group 

## is target met or elgiible or prioiritsation step - single indicator 

## children 1st doses
n_vaccination_t_S_children[] <- 0
n_vaccination_t_S_children[] <-
  if (sum(n_eligible_for_dose1_children[]) == 0) 0 else
    min(floor((daily_doses_children_t[2] * S[i, 2] *
                 prioritisation_strategy_children[
                   i, prioritisation_step_1st_dose_children]) /
                sum(n_eligible_for_dose1_children[])),
        S[i, 2])

n_vaccination_t_S_adults[] <- 0

## adults 1st doses
n_vaccination_t_S_adults[] <-
  if (sum(n_eligible_for_dose1_adults[]) == 0) 0 else
    min(floor((daily_doses_adults_t[2] * S[i, 2] *
                 prioritisation_strategy_adults[
                   i, prioritisation_step_1st_dose_adults]) /
                sum(n_eligible_for_dose1_adults[])),
        S[i, 2])

## combine total first doses
n_vaccination_t_S[, ] <- 0
n_vaccination_t_S[, 2] <-
  n_vaccination_t_S_children[i] + n_vaccination_t_S_adults[i]

## for the boundary case do an extra check that we haven't gone over the number
## of people in each compartment
n_vaccination_t_S[3, 2] <- min(n_vaccination_t_S[3, 2], S[3, 2])

## allocate 2nd doses (adults only for now)
n_vaccination_t_S[, 3] <- if (sum(n_eligible_for_dose2_adults[]) == 0) 0 else
  min(floor((daily_doses_adults_t[3] * S[i, 3] *
               prioritisation_strategy_adults[
                 i, prioritisation_step_2nd_dose_adults]) /
              sum(n_eligible_for_dose2_adults[])),
      S[i, 3])


### allocate to Ea

## children 1st doses
n_vaccination_t_Ea_children[] <- 0
n_vaccination_t_Ea_children[] <-
  if (sum(n_eligible_for_dose1_children[]) == 0) 0 else
    min(floor((daily_doses_children_t[2] * Ea[i, 2] *
                 prioritisation_strategy_children[
                   i, prioritisation_step_1st_dose_children]) /
                sum(n_eligible_for_dose1_children[])),
        Ea[i, 2])

## adults 1st doses
n_vaccination_t_Ea_adults[] <- 0
n_vaccination_t_Ea_adults[] <-
  if (sum(n_eligible_for_dose1_adults[]) == 0) 0 else
    min(floor((daily_doses_adults_t[2] * Ea[i, 2] *
                 prioritisation_strategy_adults[
                   i, prioritisation_step_1st_dose_adults]) /
                sum(n_eligible_for_dose1_adults[])),
        Ea[i, 2])

## combine total first doses
n_vaccination_t_Ea[, ] <- 0
n_vaccination_t_Ea[, 2] <-
  n_vaccination_t_Ea_children[i] + n_vaccination_t_Ea_adults[i]

## for the boundary case do an extra check that we haven't gone over the number
## of people in each compartment
n_vaccination_t_Ea[3, 2] <- min(n_vaccination_t_Ea[3, 2], Ea[3, 2])

## adults 2nd doses
n_vaccination_t_Ea[, 3] <- if (sum(n_eligible_for_dose2_adults[]) == 0) 0 else
  min(floor((daily_doses_adults_t[3] * Ea[i, 3] *
               prioritisation_strategy_adults[
                 i, prioritisation_step_2nd_dose_adults]) /
              sum(n_eligible_for_dose2_adults[])),
      Ea[i, 3])


### allocate to Eb

## children 1st doses
n_vaccination_t_Eb_children[] <- 0
n_vaccination_t_Eb_children[] <-
  if (sum(n_eligible_for_dose1_children[]) == 0) 0 else
    min(floor((daily_doses_children_t[2] * Eb[i, 2] *
                 prioritisation_strategy_children[
                   i, prioritisation_step_1st_dose_children]) /
                sum(n_eligible_for_dose1_children[])),
        Eb[i, 2])

## adults 1st doses
n_vaccination_t_Eb_adults[] <- 0
n_vaccination_t_Eb_adults[] <-
  if (sum(n_eligible_for_dose1_adults[]) == 0) 0 else
    min(floor((daily_doses_adults_t[2] * Eb[i, 2] *
                 prioritisation_strategy_adults[
                   i, prioritisation_step_1st_dose_adults]) /
                sum(n_eligible_for_dose1_adults[])),
        Eb[i, 2])

## combine total first doses
n_vaccination_t_Eb[, ] <- 0
n_vaccination_t_Eb[, 2] <-
  n_vaccination_t_Eb_children[i] + n_vaccination_t_Eb_adults[i]

## for the boundary case do an extra check that we haven't gone over the number
## of people in each compartment
n_vaccination_t_Eb[3, 2] <- min(n_vaccination_t_Eb[3, 2], Eb[3, 2])

## adults 2nd doses
n_vaccination_t_Eb[, 3] <- if (sum(n_eligible_for_dose2_adults[]) == 0) 0 else
  min(floor((daily_doses_adults_t[3] * Eb[i, 3] *
               prioritisation_strategy_adults[
                 i, prioritisation_step_2nd_dose_adults]) /
              sum(n_eligible_for_dose2_adults[])),
      Eb[i, 3])


### allocate to R

## children 1st doses
n_vaccination_t_R_children[] <- 0
n_vaccination_t_R_children[] <-
  if (sum(n_eligible_for_dose1_children[]) == 0) 0 else
    min(floor((daily_doses_children_t[2] * R[i, 2] *
                 prioritisation_strategy_children[
                   i, prioritisation_step_1st_dose_children]) /
                sum(n_eligible_for_dose1_children[])),
        R[i, 2])

## adults 1st doses
n_vaccination_t_R_adults[] <- 0
n_vaccination_t_R_adults[] <-
  if (sum(n_eligible_for_dose1_adults[]) == 0) 0 else
    min(floor((daily_doses_adults_t[2] * R[i, 2] *
                 prioritisation_strategy_adults[
                   i, prioritisation_step_1st_dose_adults]) /
                sum(n_eligible_for_dose1_adults[])),
        R[i, 2])

## combine total first doses
n_vaccination_t_R[, ] <- 0
n_vaccination_t_R[, 2] <-
  n_vaccination_t_R_children[i] + n_vaccination_t_R_adults[i]

## for the boundary case do an extra check that we haven't gone over the number
## of people in each compartment
n_vaccination_t_R[3, 2] <- min(n_vaccination_t_R[3, 2], R[3, 2])

## adults 2nd doses
n_vaccination_t_R[, 3] <- if (sum(n_eligible_for_dose2_adults[]) == 0) 0 else
  min(floor((daily_doses_adults_t[3] * R[i, 3] *
               prioritisation_strategy_adults[
                 i, prioritisation_step_2nd_dose_adults]) /
              sum(n_eligible_for_dose2_adults[])),
      R[i, 3])


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
new_S[, ] <- S[i, j] + delta_S_n_vaccination[i, j] - n_SEa[i, j]
new_Ea[, ] <- Ea[i, j] + delta_Ea_n_vaccination[i, j] + delta_Ea[i, j]
new_Eb[, ] <- Eb[i, j] + delta_Eb_n_vaccination[i, j] + delta_Eb[i, j]
new_Ir[, ] <- Ir[i, j] + delta_Ir[i, j]
new_Id[, ] <- Id[i, j] + delta_Id[i, j]
new_R[, ] <- R[i, j] + delta_R_n_vaccination[i, j] + delta_R[i, j]
new_D[, ] <- D[i, j] + delta_D[i, j]

update(S[, ]) <- new_S[i, j]
update(Ea[, ]) <- new_Ea[i, j]
update(Eb[, ]) <- new_Eb[i, j]
update(Ir[, ]) <- new_Ir[i, j]
update(Id[, ]) <- new_Id[i, j]
update(R[, ]) <- new_R[i, j]
update(D[, ]) <- new_D[i, j]

## Additional outputs
new_E[, ] <- new_Ea[i, j] + new_Eb[i, j]
new_I[, ] <- new_Ir[i, j] + new_Id[i, j]
new_N[, ] <- new_S[i, j] + new_Ea[i, j] + new_Eb[i, j] + new_Ir[i, j] +
  new_Id[i, j] + new_R[i, j] + new_D[i, j]

update(E[, ]) <- new_E[i, j] 
update(I[, ]) <- new_I[i, j]
update(N[, ]) <- new_N[i, j]

# cumulative cases by transmission route
update(cases_cumulative_hh)  <- cases_cumulative_hh + sum(n_SEa_hh[, ])
update(cases_cumulative_s)   <- cases_cumulative_s + sum(n_SEa_s[, ])
update(cases_cumulative_z)   <- cases_cumulative_z + sum(n_SEa_z[, ])
update(cases_cumulative_hc) <- cases_cumulative_hc + sum(n_SEa_hc[, ])

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
new_cases_SW_12_14 <- Binomial(sum(n_SEa[17, ]), 0.5)
new_cases_SW_15_17 <- sum(n_SEa[17, ]) - new_cases_SW_12_14
new_cases_05_14 <- sum(n_SEa[2:3, ]) + new_cases_SW_12_14
new_cases_15_plus <-
  sum(n_SEa[4:16, ]) + new_cases_SW_15_17 + sum(n_SEa[18:20, ])
new_cases_CSW <- sum(n_SEa[17, ])
new_cases_ASW <- sum(n_SEa[18, ])
new_cases_SW <- new_cases_CSW + new_cases_ASW
new_cases_PBS <- sum(n_SEa[19, ])
new_cases_HCW <- sum(n_SEa[20, ])

update(cases_inc) <- cases_inc + sum(n_SEa[, ])
update(cases_inc_00_04) <- cases_inc_00_04 + new_cases_00_04
update(cases_inc_05_14) <- cases_inc_05_14 + new_cases_05_14
update(cases_inc_15_plus) <- cases_inc_15_plus + new_cases_15_plus

update(cases_inc_CSW) <- cases_inc_CSW + new_cases_CSW
update(cases_inc_ASW) <- cases_inc_ASW + new_cases_ASW
update(cases_inc_SW) <- cases_inc_SW + new_cases_SW
update(cases_inc_PBS) <- cases_inc_PBS + new_cases_PBS
update(cases_inc_HCW) <- cases_inc_HCW + new_cases_HCW

# cumulative cases
update(cases_cumulative_by_age[]) <- cases_cumulative_by_age[i] + sum(n_SEa[i, ])
update(cases_cumulative) <- cases_cumulative + sum(n_SEa[, ])
update(cases_cumulative_00_04) <- cases_cumulative_00_04 + new_cases_00_04
update(cases_cumulative_05_14) <- cases_cumulative_05_14 + new_cases_05_14
update(cases_cumulative_15_plus) <- cases_cumulative_15_plus + new_cases_15_plus
update(cases_cumulative_CSW) <- cases_cumulative_CSW + new_cases_CSW
update(cases_cumulative_ASW) <- cases_cumulative_ASW + new_cases_ASW
update(cases_cumulative_SW) <- cases_cumulative_SW + new_cases_SW
update(cases_cumulative_PBS) <- cases_cumulative_PBS + new_cases_PBS
update(cases_cumulative_HCW) <- cases_cumulative_HCW + new_cases_HCW


new_deaths_00_04 <- sum(n_IdD[1, ])
new_deaths_SW_12_14 <- Binomial(sum(n_IdD[17, ]), 0.5)
new_deaths_SW_15_17 <- sum(n_IdD[17, ]) - new_deaths_SW_12_14
new_deaths_05_14 <- sum(n_IdD[2:3, ]) + new_deaths_SW_12_14
new_deaths_15_plus <-
  sum(n_IdD[4:16, ]) + new_deaths_SW_15_17 + sum(n_IdD[18:20, ])

new_deaths_CSW <- sum(n_IdD[17, ])
new_deaths_ASW <- sum(n_IdD[18, ])
new_deaths_SW <- new_deaths_CSW + new_deaths_ASW
new_deaths_PBS <- sum(n_IdD[19, ])
new_deaths_HCW <- sum(n_IdD[20, ])
# weekly deaths
update(deaths_inc) <- deaths_inc + sum(n_IdD[, ])
update(deaths_inc_00_04) <- deaths_inc_00_04 + new_deaths_00_04
update(deaths_inc_05_14) <- deaths_inc_05_14 + new_deaths_05_14
update(deaths_inc_15_plus) <- deaths_inc_15_plus + new_deaths_15_plus

update(deaths_inc_CSW) <- deaths_inc_CSW + new_deaths_CSW
update(deaths_inc_ASW) <- deaths_inc_ASW + new_deaths_ASW
update(deaths_inc_SW) <- deaths_inc_SW + new_deaths_SW
update(deaths_inc_PBS) <- deaths_inc_PBS + new_deaths_PBS
update(deaths_inc_HCW) <- deaths_inc_HCW + new_deaths_HCW

# cumulative deaths
update(deaths_cumulative) <- deaths_cumulative + sum(n_IdD[, ])
update(deaths_cumulative_00_04) <- deaths_cumulative_00_04 + new_deaths_00_04
update(deaths_cumulative_05_14) <- deaths_cumulative_05_14 + new_deaths_05_14
update(deaths_cumulative_15_plus) <-
  deaths_cumulative_15_plus + new_deaths_15_plus
update(deaths_cumulative_CSW) <- deaths_cumulative_CSW + new_deaths_CSW
update(deaths_cumulative_ASW) <- deaths_cumulative_ASW + new_deaths_ASW
update(deaths_cumulative_SW) <- deaths_cumulative_SW + new_deaths_SW
update(deaths_cumulative_PBS) <- deaths_cumulative_PBS + new_deaths_PBS
update(deaths_cumulative_HCW) <- deaths_cumulative_HCW + new_deaths_HCW

update(S_tot) <- sum(new_S[, ])
update(E_tot) <- sum(new_E[, ])
update(I_tot) <- sum(new_I[, ])
update(R_tot) <- sum(new_R[, ])
update(D_tot) <- sum(new_D[, ])
update(N_tot) <- sum(new_N[, ])

update(total_vax) <- total_vax + sum(n_vaccination_t_S[, ]) +
  sum(n_vaccination_t_Ea[, ]) + sum(n_vaccination_t_Eb[, ]) +
  sum(n_vaccination_t_R[, ])
update(total_vax_1stdose) <- total_vax_1stdose + sum(n_vaccination_t_S[, 2]) +
  sum(n_vaccination_t_Ea[, 2]) + sum(n_vaccination_t_Eb[, 2]) +
  sum(n_vaccination_t_R[, 2])
update(total_vax_2nddose) <- total_vax_2nddose + sum(n_vaccination_t_S[, 3]) +
  sum(n_vaccination_t_Ea[, 3]) + sum(n_vaccination_t_Eb[, 3]) +
  sum(n_vaccination_t_R[, 3])

## Vaccine doses by age / KP
# group indices
# [1: 0-4,    2: 5-9,    3: 10-14,  4: 15-19,  5: 20-24,  6: 25-29,
#  7: 30-34,  8: 35-39,  9: 40-44, 10: 45-49, 11: 50-54, 12: 55-59,
# 13: 60-64, 14: 65-69, 15: 70-74, 16: 75+,   17: CSW,   18: ASW,
# 19: PBS,   20: HCW]
# Note that:
# X_05_14 includes 50% CSW (ages 12-14)
# X_15_plus includes all 50% CSW (ages 15-17), and all ASW, PBS, HCW
# vaccination indices:
## j = 2: unvaccinated -> 1st dose
## j = 3: 1st dose -> second dose

n_vaccination_t[, ] <- n_vaccination_t_S[i, j] + n_vaccination_t_Ea[i, j] +
  n_vaccination_t_Eb[i, j] + n_vaccination_t_R[i, j]

new_dose1 <- sum(n_vaccination_t[, 2])
new_dose1_00_04 <- n_vaccination_t[1, 2]
new_dose1_SW_12_14 <- round(n_vaccination_t[17, 2] * 0.5)
new_dose1_SW_15_17 <- n_vaccination_t[17, 2] - new_dose1_SW_12_14
new_dose1_05_14 <- sum(n_vaccination_t[2:3, 2]) + new_dose1_SW_12_14
new_dose1_15_plus <-
  sum(n_vaccination_t[4:16, 2]) + new_dose1_SW_15_17 + sum(n_vaccination_t[18:20, 2])
new_dose1_CSW <- n_vaccination_t[17, 2]
new_dose1_ASW <- n_vaccination_t[18, 2]
new_dose1_SW <- new_dose1_CSW + new_dose1_ASW
new_dose1_PBS <- n_vaccination_t[19, 2]
new_dose1_HCW <- n_vaccination_t[20, 2]

update(dose1_inc) <- dose1_inc + new_dose1
update(dose1_inc_00_04) <- dose1_inc_00_04 + new_dose1_00_04
update(dose1_inc_05_14) <- dose1_inc_05_14 + new_dose1_05_14
update(dose1_inc_15_plus) <- dose1_inc_15_plus + new_dose1_15_plus
update(dose1_inc_SW) <- dose1_inc_SW + new_dose1_SW
update(dose1_inc_CSW) <- dose1_inc_CSW + new_dose1_CSW
update(dose1_inc_ASW) <- dose1_inc_ASW + new_dose1_ASW
update(dose1_inc_PBS) <- dose1_inc_PBS + new_dose1_PBS
update(dose1_inc_HCW) <- dose1_inc_HCW + new_dose1_HCW

new_dose2 <- sum(n_vaccination_t[, 3])
new_dose2_00_04 <- n_vaccination_t[1, 3]
new_dose2_SW_12_14 <- round(n_vaccination_t[17, 3] * 0.5)
new_dose2_SW_15_17 <- n_vaccination_t[17, 3] - new_dose2_SW_12_14
new_dose2_05_14 <- sum(n_vaccination_t[2:3, 3]) + new_dose2_SW_12_14
new_dose2_15_plus <-
  sum(n_vaccination_t[4:16, 3]) + new_dose2_SW_15_17 + sum(n_vaccination_t[18:20, 3])
new_dose2_CSW <- n_vaccination_t[17, 3]
new_dose2_ASW <- n_vaccination_t[18, 3]
new_dose2_SW <- new_dose2_CSW + new_dose2_ASW
new_dose2_PBS <- n_vaccination_t[19, 3]
new_dose2_HCW <- n_vaccination_t[20, 3]

update(dose2_inc) <- dose2_inc + new_dose2
update(dose2_inc_00_04) <- dose2_inc_00_04 + new_dose2_00_04
update(dose2_inc_05_14) <- dose2_inc_05_14 + new_dose2_05_14
update(dose2_inc_15_plus) <- dose2_inc_15_plus + new_dose2_15_plus
update(dose2_inc_SW) <- dose2_inc_SW + new_dose2_SW
update(dose2_inc_CSW) <- dose2_inc_CSW + new_dose2_CSW
update(dose2_inc_ASW) <- dose2_inc_ASW + new_dose2_ASW
update(dose2_inc_PBS) <- dose2_inc_PBS + new_dose2_PBS
update(dose2_inc_HCW) <- dose2_inc_HCW + new_dose2_HCW

update(dose1_cumulative) <- dose1_cumulative + new_dose1
update(dose1_cumulative_00_04) <- dose1_cumulative_00_04 + new_dose1_00_04
update(dose1_cumulative_05_14) <- dose1_cumulative_05_14 + new_dose1_05_14
update(dose1_cumulative_15_plus) <- dose1_cumulative_15_plus + new_dose1_15_plus
update(dose1_cumulative_CSW) <- dose1_cumulative_CSW + new_dose1_CSW
update(dose1_cumulative_ASW) <- dose1_cumulative_ASW + new_dose1_ASW
update(dose1_cumulative_SW) <- dose1_cumulative_SW + new_dose1_SW
update(dose1_cumulative_PBS) <- dose1_cumulative_PBS + new_dose1_PBS
update(dose1_cumulative_HCW) <- dose1_cumulative_HCW + new_dose1_HCW

update(dose2_cumulative) <- dose2_cumulative + new_dose2
update(dose2_cumulative_00_04) <- dose2_cumulative_00_04 + new_dose2_00_04
update(dose2_cumulative_05_14) <- dose2_cumulative_05_14 + new_dose2_05_14
update(dose2_cumulative_15_plus) <- dose2_cumulative_15_plus + new_dose2_15_plus
update(dose2_cumulative_CSW) <- dose2_cumulative_CSW + new_dose2_CSW
update(dose2_cumulative_ASW) <- dose2_cumulative_ASW + new_dose2_ASW
update(dose2_cumulative_SW) <- dose2_cumulative_SW + new_dose2_SW
update(dose2_cumulative_PBS) <- dose2_cumulative_PBS + new_dose2_PBS
update(dose2_cumulative_HCW) <- dose2_cumulative_HCW + new_dose2_HCW
  

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
m_gen_pop <- parameter()
m_sex <- parameter()
# I adjusted for reduced transmissibility
I_infectious[, ] <- I[i, j] * (1 - ve_T[j])

prop_infectious[] <-
  if (sum(N[i, ]) == 0) 0 else sum(I_infectious[i, ]) / sum(N[i, ])
# Generating Force of Infection
# for susceptible age i, % contacts infectious age j
s_ij_gen_pop[, ] <- m_gen_pop[i, j] * prop_infectious[j]
# as above but for the sexual contacts only
s_ij_sex[, ] <- m_sex[i, j] * prop_infectious[j]

lambda_hh[, ] <- beta_h * sum(s_ij_gen_pop[i, ]) * (1 - ve_I[i, j])
lambda_s[, ] <- beta_s * sum(s_ij_sex[i, ]) * (1 - ve_I[i, j])
# additional foi in HCW only (i = 20) homogeneous from infected as assumed equally
# likely to attend hospital
lambda_hc[, ] <-
  if (i == 20) beta_hcw * sum(I_infectious[, ]) * (1 - ve_I[i, j]) else 0
lambda_z[, ] <- beta_z[i] * (1 - ve_I[i, j])

lambda[, ] <- lambda_hh[i, j] + lambda_s[i, j] + lambda_hc[i, j] + lambda_z[i, j] 

## Draws from binomial distributions for numbers changing between compartments
# accounting for vaccination:
n_SEa[, ] <- Binomial(S[i, j] + delta_S_n_vaccination[i, j], p_SE[i, j])

p_hh[, ]  <- if (lambda[i, j] > 0) lambda_hh[i, j] / lambda[i, j] else 0
p_s[, ]   <- if (lambda[i, j] > 0) lambda_s[i, j] / lambda[i, j] else 0
p_hc[, ] <- if (lambda[i, j] > 0) lambda_hc[i, j] / lambda[i, j] else 0

## Split n_SEa by transmission route
n_SEa_hh[, ] <- Binomial(n_SEa[i, j], p_hh[i, j])
n_SEa_s[, ] <- Binomial(n_SEa[i, j] - n_SEa_hh[i, j], p_s[i, j])
n_SEa_hc[, ] <- Binomial(n_SEa[i, j] - n_SEa_hh[i, j] - n_SEa_s[i, j],
                          p_hc[i, j])
n_SEa_z[, ] <- n_SEa[i, j] - n_SEa_hh[i, j] - n_SEa_s[i, j] - n_SEa_hc[i, j]


n_EaEb[, ] <- Binomial(Ea[i, j] + delta_Ea_n_vaccination[i, j], p_EE)
n_EbI[, ] <- Binomial(Eb[i, j] + delta_Eb_n_vaccination[i, j], p_EI)
# Proportion of the infections that will die, impact of vaccination included
# already as part of inputs
n_EbId[, ] <- Binomial(n_EbI[i, j], CFR[i, j])
# whatever infections don't die go to R
n_EbIr[, ] <- n_EbI[i, j] - n_EbId[i, j]

n_IrR[, ] <- Binomial(Ir[i, j], p_IrR)
n_IdD[, ] <- Binomial(Id[i, j], p_IdD)

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
initial(cases_inc, zero_every = 7) <- 0
initial(deaths_inc, zero_every = 7) <- 0
initial(cases_cumulative) <- 0
initial(deaths_cumulative) <- 0

initial(cases_cumulative_hh)  <- 0
initial(cases_cumulative_s)   <- 0
initial(cases_cumulative_z)   <- 0
initial(cases_cumulative_hc) <- 0

initial(cases_inc_00_04, zero_every = 7) <- 0
initial(cases_inc_05_14, zero_every = 7) <- 0
initial(cases_inc_15_plus, zero_every = 7) <- 0
initial(cases_inc_PBS, zero_every = 7) <- 0
initial(cases_inc_CSW, zero_every = 7) <- 0
initial(cases_inc_ASW, zero_every = 7) <- 0
initial(cases_inc_SW, zero_every = 7) <- 0
initial(cases_inc_HCW, zero_every = 7) <- 0
initial(deaths_inc_00_04, zero_every = 7) <- 0
initial(deaths_inc_05_14, zero_every = 7) <- 0
initial(deaths_inc_15_plus, zero_every = 7) <- 0
initial(deaths_inc_PBS, zero_every = 7) <- 0
initial(deaths_inc_CSW, zero_every = 7) <- 0
initial(deaths_inc_ASW, zero_every = 7) <- 0
initial(deaths_inc_SW, zero_every = 7) <- 0
initial(deaths_inc_HCW, zero_every = 7) <- 0

initial(dose1_inc, zero_every = 7) <- 0
initial(dose1_inc_00_04, zero_every = 7) <- 0
initial(dose1_inc_05_14, zero_every = 7) <- 0
initial(dose1_inc_15_plus, zero_every = 7) <- 0
initial(dose1_inc_PBS, zero_every = 7) <- 0
initial(dose1_inc_CSW, zero_every = 7) <- 0
initial(dose1_inc_ASW, zero_every = 7) <- 0
initial(dose1_inc_SW, zero_every = 7) <- 0
initial(dose1_inc_HCW, zero_every = 7) <- 0

initial(dose2_inc, zero_every = 7) <- 0
initial(dose2_inc_00_04, zero_every = 7) <- 0
initial(dose2_inc_05_14, zero_every = 7) <- 0
initial(dose2_inc_15_plus, zero_every = 7) <- 0
initial(dose2_inc_PBS, zero_every = 7) <- 0
initial(dose2_inc_CSW, zero_every = 7) <- 0
initial(dose2_inc_ASW, zero_every = 7) <- 0
initial(dose2_inc_SW, zero_every = 7) <- 0
initial(dose2_inc_HCW, zero_every = 7) <- 0

initial(cases_cumulative_by_age[]) <- 0
initial(cases_cumulative_00_04) <- 0
initial(cases_cumulative_05_14) <- 0
initial(cases_cumulative_15_plus) <- 0
initial(cases_cumulative_PBS) <- 0
initial(cases_cumulative_CSW) <- 0
initial(cases_cumulative_ASW) <- 0
initial(cases_cumulative_SW) <- 0
initial(cases_cumulative_HCW) <- 0
initial(deaths_cumulative_00_04) <- 0
initial(deaths_cumulative_05_14) <- 0
initial(deaths_cumulative_15_plus) <- 0
initial(deaths_cumulative_PBS) <- 0
initial(deaths_cumulative_CSW) <- 0
initial(deaths_cumulative_ASW) <- 0
initial(deaths_cumulative_SW) <- 0
initial(deaths_cumulative_HCW) <- 0

initial(dose1_cumulative) <- 0
initial(dose1_cumulative_00_04) <- 0
initial(dose1_cumulative_05_14) <- 0
initial(dose1_cumulative_15_plus) <- 0
initial(dose1_cumulative_PBS) <- 0
initial(dose1_cumulative_CSW) <- 0
initial(dose1_cumulative_ASW) <- 0
initial(dose1_cumulative_SW) <- 0
initial(dose1_cumulative_HCW) <- 0

initial(dose2_cumulative) <- 0
initial(dose2_cumulative_00_04) <- 0
initial(dose2_cumulative_05_14) <- 0
initial(dose2_cumulative_15_plus) <- 0
initial(dose2_cumulative_PBS) <- 0
initial(dose2_cumulative_CSW) <- 0
initial(dose2_cumulative_ASW) <- 0
initial(dose2_cumulative_SW) <- 0
initial(dose2_cumulative_HCW) <- 0

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
S0 <- parameter()
Ea0 <- parameter()
Eb0 <- parameter()
Ir0 <- parameter()
Id0 <- parameter()
R0 <- parameter()
D0 <- parameter()

##Parameters
beta_h <- parameter()
beta_s <- parameter()
beta_z <- parameter()
beta_hcw <- parameter()
gamma_E <- parameter()
gamma_Ir <- parameter()
gamma_Id <- parameter()
CFR <- parameter()

#vaccine efficacy parameters
ve_T <- parameter()
ve_I <- parameter()
#ve_D[,] <- parameter() # this is included within the CFR

#Number of age classes & number of transmissibility classes
n_vax <- parameter()
n_group <- parameter()

## Dimensions of the different "vectors" here vectors stand for
## multi-dimensional arrays
dim(N, new_N) <- c(n_group, n_vax)
dim(S, new_S) <- c(n_group, n_vax)
dim(S0) <- c(n_group, n_vax)
dim(p_SE) <- c(n_group, n_vax)
dim(n_SEa) <- c(n_group, n_vax)

dim(Ea, new_Ea) <- c(n_group, n_vax)
dim(Ea0) <- c(n_group, n_vax)
dim(Eb0) <- c(n_group, n_vax)
dim(delta_Ea) <- c(n_group, n_vax)
dim(n_EaEb) <- c(n_group, n_vax)

dim(Eb, new_Eb) <- c(n_group, n_vax)
dim(delta_Eb) <- c(n_group, n_vax)
dim(n_EbI) <- c(n_group, n_vax)

dim(n_EbId) <- c(n_group, n_vax)
dim(n_EbIr) <- c(n_group, n_vax)
dim(E, new_E) <- c(n_group, n_vax)

dim(Ir0) <- c(n_group, n_vax)
dim(Ir, new_Ir) <- c(n_group, n_vax)
dim(delta_Ir) <- c(n_group, n_vax)
dim(n_IrR) <- c(n_group, n_vax)

dim(Id0) <- c(n_group, n_vax)
dim(Id, new_Id) <- c(n_group, n_vax)
dim(delta_Id) <- c(n_group, n_vax)
dim(n_IdD) <- c(n_group, n_vax)
dim(I, new_I) <- c(n_group, n_vax)

dim(R, new_R) <- c(n_group, n_vax)
dim(R0) <- c(n_group, n_vax)
dim(delta_R) <- c(n_group, n_vax)

dim(D, new_D) <- c(n_group, n_vax)
dim(D0) <- c(n_group, n_vax)
dim(delta_D) <- c(n_group, n_vax)

dim(lambda) <- c(n_group, n_vax)
dim(lambda_hh) <- c(n_group, n_vax)
dim(lambda_s) <- c(n_group, n_vax)
dim(lambda_hc) <- c(n_group, n_vax)
dim(lambda_z) <- c(n_group, n_vax)
dim(p_hh) <- c(n_group, n_vax)
dim(p_s) <- c(n_group, n_vax)
dim(p_hc) <- c(n_group, n_vax)
dim(n_SEa_hh) <- c(n_group, n_vax)
dim(n_SEa_s) <- c(n_group, n_vax)
dim(n_SEa_hc) <- c(n_group, n_vax)
dim(n_SEa_z) <- c(n_group, n_vax)
dim(m_gen_pop) <- c(n_group, n_group)
dim(m_sex) <- c(n_group, n_group)
dim(I_infectious) <- c(n_group, n_vax)
dim(prop_infectious) <- c(n_group)
dim(s_ij_gen_pop) <- c(n_group, n_group)
dim(s_ij_sex) <- c(n_group, n_group)
dim(beta_z) <- c(n_group)

dim(CFR) <- c(n_group, n_vax)

dim(ve_T) <- c(n_vax)
dim(ve_I) <- c(n_group, n_vax)

dim(cases_cumulative_by_age) <- n_group

dim(n_vaccination_t_S) <- c(n_group, n_vax)
dim(n_vaccination_t_Ea) <- c(n_group, n_vax)
dim(n_vaccination_t_Eb) <- c(n_group, n_vax)
dim(n_vaccination_t_R) <- c(n_group, n_vax)
dim(n_vaccination_t) <-  c(n_group, n_vax)

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


#### Compare functions
# Options are Negative Binomial [Aggregate | By-age] for cases + deaths
# Plus optional Binomial on % cases in HCW and/or SW

exp_noise <- parameter(1e+06)

## cases
# Aggregate
alpha_cases <- parameter()
cases <- data()
model_cases <- cases_inc + Exponential(exp_noise)
cases ~ NegativeBinomial(size = 1 / alpha_cases, mu = model_cases)

# By-age
alpha_cases_00_04 <- parameter()
cases_00_04 <- data()
model_cases_00_04 <- cases_inc_00_04 + Exponential(exp_noise)
cases_00_04 ~
  NegativeBinomial(size = 1 / alpha_cases_00_04, mu = model_cases_00_04)

alpha_cases_05_14 <- parameter()
cases_05_14 <- data()
model_cases_05_14 <- cases_inc_05_14 + Exponential(exp_noise)
cases_05_14 ~ 
  NegativeBinomial(size = 1 / alpha_cases_05_14, mu = model_cases_05_14)

alpha_cases_15_plus <- parameter()
cases_15_plus <- data()
model_cases_15_plus <- cases_inc_15_plus + Exponential(exp_noise)
cases_15_plus ~ 
  NegativeBinomial(size = 1 / alpha_cases_15_plus, mu = model_cases_15_plus)

## deaths
# Aggregate
alpha_deaths <- parameter()
deaths <- data()
model_deaths <- deaths_inc + Exponential(exp_noise)
deaths ~ NegativeBinomial(size = 1 / alpha_deaths, mu = model_deaths)

# By-age
alpha_deaths_00_04 <- parameter()
deaths_00_04 <- data()
model_deaths_00_04 <- deaths_inc_00_04 + Exponential(exp_noise)
deaths_00_04 ~ 
  NegativeBinomial(size = 1 / alpha_deaths_00_04, mu = model_deaths_00_04)

alpha_deaths_05_14 <- parameter()
deaths_05_14 <- data()
model_deaths_05_14 <- deaths_inc_05_14 + Exponential(exp_noise)
deaths_05_14 ~ 
  NegativeBinomial(size = 1 / alpha_deaths_05_14, mu = model_deaths_05_14)

alpha_deaths_15_plus <- parameter()
deaths_15_plus <- data()
model_deaths_15_plus <- deaths_inc_15_plus + Exponential(exp_noise)
deaths_15_plus ~ 
  NegativeBinomial(size = 1 / alpha_deaths_15_plus, mu = model_deaths_15_plus)

# Cumulative CFR
cfr_00_04 <- data()
cfr_00_04 ~ Beta(deaths_cumulative_00_04,
                 cases_cumulative_00_04 - deaths_cumulative_00_04)
cfr_05_14 <- data()
cfr_05_14 ~ Beta(deaths_cumulative_05_14,
                 cases_cumulative_05_14 - deaths_cumulative_05_14)
cfr_15_plus <- data()
cfr_15_plus ~ Beta(deaths_cumulative_15_plus,
                   cases_cumulative_15_plus - deaths_cumulative_15_plus)

# Proportion of cases in key pops
# create a data stream of aggregated cases that will work regardless of whether
# fitting is by age or in aggregate
cases_total <- data()

cases_HCW <- data()
model_cases_HCW <- cases_inc_HCW +  Exponential(exp_noise)
model_cases_non_HCW <- cases_inc - cases_inc_HCW +  Exponential(exp_noise)
model_prop_HCW <- model_cases_HCW / (model_cases_HCW + model_cases_non_HCW)
cases_HCW ~ Binomial(cases_total, model_prop_HCW)

cases_SW <- data()
model_cases_SW <- cases_inc_SW +  Exponential(exp_noise)
model_cases_non_SW <- cases_inc - cases_inc_SW +  Exponential(exp_noise)
model_prop_SW <- model_cases_SW / (model_cases_SW + model_cases_non_SW)
cases_SW ~ Binomial(cases_total, model_prop_SW)
