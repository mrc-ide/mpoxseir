# State indices: [i: age; j: vaccination]

#### Vaccination

## j = 1: previous smallpox vaccine (efficacy roughly at one dose)
## j = 2: unvaccinated
## j = 3: 1 dose
## j = 4: 2 doses

## specify the daily number of vaccinations that can happen at each time point
## first and last columns must be 0; last row must be 0 (pre-processing job)
## column j corresponds to the people in j receiving a vaccine (so moving
## to j + 1)
daily_doses_value <- parameter(constant = TRUE)
dim(daily_doses_value) <- parameter(rank = 2)
daily_doses_time <- parameter(constant = TRUE)
dim(daily_doses_time) <- parameter(rank = 1)

## with this we need to ensure daily_doses final time step is all 0 - done
## outside of the model in pre-processing
daily_doses_t[] <- interpolate(daily_doses_time, daily_doses_value, "constant")
dim(daily_doses_t) <- c(n_vax)

## allocate the daily doses according to prioritisation strategy

## the number of different prioritisation strategies that we are considering
N_prioritisation_steps <- parameter()
## what this corresponds to in terms of groups is encoded here
## this is independent of dose / applies to both doeses
prioritisation_strategy <- parameter()
dim(prioritisation_strategy) <- c(n_group, N_prioritisation_steps)

## vaccination targets: for each group, what is the level of coverage we want to
## see before we move onto the next set of groups to target?
## when we input this do the calc outside by the population size as part of
## pre-processing to save faffing around in here.
## here, j=1,...,n_vax corresponds to the coverage wanted in j, so for j=0
## and j=1 this stays at 0
vaccination_coverage_target_1st_dose_prop <- parameter()
vaccination_coverage_target_2nd_dose_prop <- parameter()

## what is the current state of the vaccination targets? e.g. what
## prioritisation step are we on, this depends on whether we have met our
## targets

## has the target been met
## for cols 1 (smallpox vax) and 2 (unvaccinated) this doesn't matter
dim(target_met_t) <- c(n_group, n_vax)
target_met_t[, ] <- 0

## 1st doses
## if you have a 2nd dose this implies you also have had a 1st dose so account
## for this in the 1st dose target
target_met_t[, 3] <-
  if (sum(N[, 3:4]) >
      round(prioritisation_strategy[i, prioritisation_step_1st_dose] *
            vaccination_coverage_target_1st_dose_prop * sum(N[i, 1:4])))
    1 else 0


## 2nd doses
target_met_t[, 4] <-
  if (sum(N[, 4]) >
      round(prioritisation_strategy[i, prioritisation_step_2nd_dose] *
            vaccination_coverage_target_2nd_dose_prop * sum(N[i, 1:4])))
    1 else 0

## prioritisation step proposal to account for the fact that this would update
## every single time step if we vaccinate quickly enough (unlikely but would
## like it to be done properly)
prioritisation_step_1st_dose_proposal <-
  if (sum(target_met_t[, 3]) == n_group) prioritisation_step_1st_dose + 1 else
    prioritisation_step_1st_dose

update(prioritisation_step_1st_dose) <-
  if (prioritisation_step_1st_dose_proposal > N_prioritisation_steps)
    N_prioritisation_steps else prioritisation_step_1st_dose_proposal

prioritisation_step_2nd_dose_proposal <-
  if (sum(target_met_t[, 4]) == n_group) prioritisation_step_2nd_dose + 1 else
    prioritisation_step_2nd_dose

update(prioritisation_step_2nd_dose) <-
  if (prioritisation_step_2nd_dose_proposal > N_prioritisation_steps)
    N_prioritisation_steps else prioritisation_step_2nd_dose_proposal


#prioritisation_step <- user()
#initial(prioritisation_step) <- 1
initial(prioritisation_step_1st_dose) <- 1
initial(prioritisation_step_2nd_dose) <- 1

## now that we know the step we are on, we know who is eligible per age group to
## be vaccinated

## who can get a 1st dose
## now we have to look in the preceding j class (e.g. who in unvaccinated j = 1
## can get a 1st dose and move to j = 2) as these are the people who will be
## eligible
n_eligible_for_dose1[] <- (S[i, 2] + Ea[i, 2] + Eb[i, 2] + R[i, 2]) *
  prioritisation_strategy[i, prioritisation_step_1st_dose]
dim(n_eligible_for_dose1) <- c(n_group)
## who can get a 2nd dose
n_eligible_for_dose2[] <- (S[i, 3] + Ea[i, 3] + Eb[i, 3] + R[i, 3]) *
  prioritisation_strategy[i, prioritisation_step_2nd_dose]
dim(n_eligible_for_dose2) <- c(n_group)

## vaccine uptake
## account for the fact that some groups (in n_group) may be less inclined to
## take the vaccine
vaccine_uptake <- parameter()
dim(vaccine_uptake) <- c(n_group)

## allocate the doses to the unvaccinated by age group, prioritisation strategy
## and across S, E, R
## hacky fix for now

## allocate doses to S
n_vaccination_t_S[, ] <- 0

## allocate 1st doses
n_vaccination_t_S[, 2] <- min(
  floor((daily_doses_t[2] * S[i, 2] *
           prioritisation_strategy[i, prioritisation_step_1st_dose] *
           vaccine_uptake[i]) / sum(n_eligible_for_dose1[])),
  S[i, 2])

## allocate 2nd doses
n_vaccination_t_S[, 3] <- min(
  floor((daily_doses_t[3] * S[i, 3] *
           prioritisation_strategy[i, prioritisation_step_2nd_dose] *
           vaccine_uptake[i]) / sum(n_eligible_for_dose2[])),
  S[i, 3])

## allocate doses to Ea
n_vaccination_t_Ea[, ] <- 0

## allocate 1st doses
n_vaccination_t_Ea[, 2] <- min(
  floor((daily_doses_t[2] * Ea[i, 2] *
           prioritisation_strategy[i, prioritisation_step_1st_dose] *
           vaccine_uptake[i]) / sum(n_eligible_for_dose1[])),
  Ea[i, 2])

## allocate 2nd doses
n_vaccination_t_Ea[, 3] <- min(
  floor((daily_doses_t[3] * Ea[i, 3] *
           prioritisation_strategy[i, prioritisation_step_2nd_dose] *
           vaccine_uptake[i]) / sum(n_eligible_for_dose2[])),
  Ea[i, 3])

## allocate doses to Eb
n_vaccination_t_Eb[, ] <- 0

## allocate 1st doses
n_vaccination_t_Eb[, 2] <- min(
  floor((daily_doses_t[2] * Eb[i, 2] *
           prioritisation_strategy[i, prioritisation_step_1st_dose] *
           vaccine_uptake[i]) / sum(n_eligible_for_dose1[])),
  Eb[i, 2])

## allocate 2nd doses
n_vaccination_t_Eb[, 3] <- min(
  floor((daily_doses_t[3] * Eb[i, 3] *
           prioritisation_strategy[i, prioritisation_step_2nd_dose] *
           vaccine_uptake[i]) / sum(n_eligible_for_dose2[])),
  Eb[i, 3])

## allocate doses to R
n_vaccination_t_R[, ] <- 0

## allocate 1st doses
n_vaccination_t_R[, 2] <- min(
  floor((daily_doses_t[2] * R[i, 2] *
           prioritisation_strategy[i, prioritisation_step_1st_dose] *
           vaccine_uptake[i]) / sum(n_eligible_for_dose1[])),
  R[i, 2])

## allocate 2nd doses
n_vaccination_t_R[, 3] <- min(
  floor((daily_doses_t[3] * R[i, 3] *
           prioritisation_strategy[i, prioritisation_step_2nd_dose] *
           vaccine_uptake[i]) / sum(n_eligible_for_dose2[])),
  R[i, 3])


## net vaccination change for relevant classes (S,Ea,Eb,R)
## logic here depends on vaccine class you are in (e.g. can only increase
## through j)
## j=1 corresponds to previous smallpox vaccination
## this is generic so accounts for n_vax whatever it is
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
update(N[, ]) <- S[i, j] + Ea[i, j] + Eb[i, j] + Ir[i, j] + Id[i, j] + R[i, j] +
  D[i, j]
#update(N[, ]) <- S[i, j] + E[i, j] + I[i, j] + R[i, j] + D[i, j]

# weekly cases
update(cases_inc) <- cases_inc + sum(n_SEa[, ])
update(cases_inc_00_04) <- cases_inc_00_04 + sum(n_SEa[1, ])
update(cases_inc_05_14) <- cases_inc_05_14 + sum(n_SEa[2:3, ])
update(cases_inc_15_plus) <- cases_inc_15_plus + sum(n_SEa[4:16, ])
update(cases_inc_PBS) <- cases_inc_PBS + sum(n_SEa[17, ])
update(cases_inc_SW) <- cases_inc_SW + sum(n_SEa[18, ])

# cumulative cases
update(cases_cumulative) <- cases_cumulative + sum(n_SEa[, ])
update(cases_cumulative_00_04) <- cases_cumulative_00_04 + sum(n_SEa[1, ])
update(cases_cumulative_05_14) <- cases_cumulative_05_14 + sum(n_SEa[2:3, ])
update(cases_cumulative_15_plus) <- cases_cumulative_15_plus +
  sum(n_SEa[4:16, ])
update(cases_cumulative_PBS) <- cases_cumulative_PBS + sum(n_SEa[17, ])
update(cases_cumulative_SW) <- cases_cumulative_SW + sum(n_SEa[18, ])

# weekly deaths
update(deaths_inc) <- deaths_inc + sum(n_IdD[, ])
update(deaths_inc_00_04) <- deaths_inc_00_04 + sum(n_IdD[1, ])
update(deaths_inc_05_14) <- deaths_inc_05_14 + sum(n_IdD[2:3, ])
update(deaths_inc_15_plus) <- deaths_inc_15_plus + sum(n_IdD[4:16, ])
update(deaths_inc_PBS) <- deaths_inc_PBS + sum(n_IdD[17, ])
update(deaths_inc_SW) <- deaths_inc_SW + sum(n_IdD[18, ])

# cumulative deaths
update(deaths_cumulative) <- deaths_cumulative + sum(n_IdD[, ])
update(deaths_cumulative_00_04) <- deaths_cumulative_00_04 + sum(n_IdD[1, ])
update(deaths_cumulative_05_14) <- deaths_cumulative_05_14 + sum(n_IdD[2:3, ])
update(deaths_cumulative_15_plus) <- deaths_cumulative_15_plus +
  sum(n_IdD[4:16, ])
update(deaths_cumulative_PBS) <- deaths_cumulative_PBS + sum(n_IdD[17, ])
update(deaths_cumulative_SW) <- deaths_cumulative_SW + sum(n_IdD[18, ])

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
m_gen_pop <- parameter()
m_sex <- parameter()
# I adjusted for reduced transmissibility
I_infectious[, ] <- I[i, j] * (1 - ve_T[j])

prop_infectious[] <- sum(I_infectious[i, ]) / sum(N[i, ])
# Generating Force of Infection
# for susceptible age i, % contacts infectious age j
s_ij_gen_pop[, ] <- m_gen_pop[i, j] * prop_infectious[j]
# as above but for the sexual contacts only
s_ij_sex[, ] <- m_sex[i, j] * prop_infectious[j]
lambda[, ] <- ((beta_h * sum(s_ij_gen_pop[i, ])) +
                 (beta_s * sum(s_ij_sex[i, ])) + beta_z[i]) * (1 - ve_I[j])

## Draws from binomial distributions for numbers changing between compartments
# accounting for vaccination:
n_SEa[, ] <- Binomial(S[i, j] + delta_S_n_vaccination[i, j], p_SE[i, j])
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
initial(cases_cumulative, zero_every = 7) <- 0
initial(deaths_cumulative, zero_every = 7) <- 0

initial(cases_inc_00_04, zero_every = 7) <- 0
initial(cases_inc_05_14, zero_every = 7) <- 0
initial(cases_inc_15_plus, zero_every = 7) <- 0
initial(cases_inc_PBS, zero_every = 7) <- 0
initial(cases_inc_SW, zero_every = 7) <- 0
initial(deaths_inc_00_04, zero_every = 7) <- 0
initial(deaths_inc_05_14, zero_every = 7) <- 0
initial(deaths_inc_15_plus, zero_every = 7) <- 0
initial(deaths_inc_PBS, zero_every = 7) <- 0
initial(deaths_inc_SW, zero_every = 7) <- 0

initial(cases_cumulative_00_04) <- 0
initial(cases_cumulative_05_14) <- 0
initial(cases_cumulative_15_plus) <- 0
initial(cases_cumulative_PBS) <- 0
initial(cases_cumulative_SW) <- 0
initial(deaths_cumulative_00_04) <- 0
initial(deaths_cumulative_05_14) <- 0
initial(deaths_cumulative_15_plus) <- 0
initial(deaths_cumulative_PBS) <- 0
initial(deaths_cumulative_SW) <- 0

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
gamma_E <- parameter()
gamma_Ir <- parameter()
gamma_Id <- parameter()
CFR <- parameter()

#vaccine efficacy parameters
ve_T <- parameter()
ve_I <- parameter()
#ve_D[,] <- user() # this is included within the CFR

#Number of age classes & number of transmissibility classes
n_vax <- parameter()
n_group <- parameter()

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
dim(ve_I) <- c(n_vax)

vaccination_campaign_length <- parameter()

dim(n_vaccination_t_S) <- c(n_group, n_vax)
dim(n_vaccination_t_Ea) <- c(n_group, n_vax)
dim(n_vaccination_t_Eb) <- c(n_group, n_vax)
dim(n_vaccination_t_R) <- c(n_group, n_vax)

dim(delta_S_n_vaccination) <- c(n_group, n_vax)
dim(delta_Ea_n_vaccination) <- c(n_group, n_vax)
dim(delta_Eb_n_vaccination) <- c(n_group, n_vax)
dim(delta_R_n_vaccination) <- c(n_group, n_vax)


#### compare

exp_noise <- parameter(1e+06)

# cases
cases <- data()
model_cases <- cases_inc + Exponential(exp_noise)
cases ~ Poisson(model_cases)

cases_00_04 <- data()
model_cases_00_04 <- cases_inc_00_04 + Exponential(exp_noise)
cases_00_04 ~ Poisson(model_cases_00_04)

cases_05_14 <- data()
model_cases_05_14 <- cases_inc_05_14 + Exponential(exp_noise)
cases_05_14 ~ Poisson(model_cases_05_14)

cases_15_plus <- data()
model_cases_15_plus <- cases_inc_15_plus + cases_inc_PBS + cases_inc_SW + Exponential(exp_noise)
cases_15_plus ~ Poisson(model_cases_15_plus)

# deaths
deaths <- data()
model_deaths <- deaths_inc + Exponential(exp_noise)
deaths ~ Poisson(model_deaths)

deaths_00_04 <- data()
model_deaths_00_04 <- deaths_inc_00_04 + Exponential(exp_noise)
deaths_00_04 ~ Poisson(model_deaths_00_04)

deaths_05_14 <- data()
model_deaths_05_14 <- deaths_inc_05_14 + Exponential(exp_noise)
deaths_05_14 ~ Poisson(model_deaths_05_14)

deaths_15_plus <- data()
model_deaths_15_plus <- deaths_inc_15_plus + deaths_inc_PBS + deaths_inc_SW + Exponential(exp_noise)
deaths_15_plus ~ Poisson(model_deaths_15_plus)
