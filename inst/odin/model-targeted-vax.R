# State indices: [i: age; j: vaccination]

# Time Steps
dt <- user(1)
steps_per_week <- 7 / dt
initial(time) <- step
update(time) <- (step + 1) * dt
#output(time) <- TRUE

#### Vaccination

## specify the daily number of vaccinations that can happen at each time point
daily_doses[,] <- user()
dim(daily_doses) <- c(vaccination_campaign_length,n_vax)

## with this we need to ensure daily_doses final time step is all 0 - done outside of the model in pre-processing
daily_doses_t[,] <- if (as.integer(time) >= (vaccination_campaign_length))
  daily_doses[vaccination_campaign_length,j] else daily_doses[time,j]
dim(daily_doses_t) <- c(1,n_vax)

## where are we allocating the daily doses for that time step?
## we have to look at the prioritisation strategy

## the number of different prioritisation strategies that we are considering
N_prioritisation_steps <- user()
#priotisation_step <- user() ### must intialise at 1
## what this corresponds to in terms of groups is encoded here
prioritisation_strategy[,] <- user()
dim(prioritisation_strategy) <- c(n_group,N_prioritisation_steps)

## vaccination targets: for each group, what is the level of coverage we want to see before we move onto the next set of groups to target?
## when we input this do the calc outside by the population size as part of pre-processing to save faffing around in here
vaccination_coverage_target[,] <- user()
dim(vaccination_coverage_target) <- c(n_group,N_prioritisation_steps)

## what is the current state of the vaccination targets? e.g. what prioritisation step are we on
# this depends on whether we have met our targets
## total vaccinated in the model comes from the columns of N for which j>1

## total number of people vaccinated per group
total_vaccinated_t[,] <- sum(N[i,j]) - N[i,1]
dim(total_vaccinated_t) <- c(n_group,1)

## has the target been met
target_met_t[,] <- if(total_vaccinated_t[i,1]>vaccination_coverage_target[i,prioritisation_step]) 1 else 0
dim(target_met_t) <- c(n_group,1)

update(prioritisation_step) <- if(sum(target_met_t[,])==n_group) prioritisation_step + 1 else prioritisation_step

#prioritisation_step <- user()
initial(prioritisation_step) <- 1

## now that we know the step we are on, we know who is eligible per age group to be vaccinated
## this targets only the unvaccinated but we have to change this if we do end up modelling 1 -> 2 doses
n_eligible[] <- (S[i,1] + Ea[i,1] + Eb[i,1] + R[i,1])*prioritisation_strategy[i,prioritisation_step]
dim(n_eligible) <- c(n_group)


## vaccine uptake
## account for the fact that some groups (in n_group) may be less inclined to take the vaccine
vaccine_uptake[] <- user()
dim(vaccine_uptake) <- c(n_group)


## allocate the doses to the unvaccinated by age group, prioritisation strategy and across S, E, R
## daily_doses_t should be a vector and not summed over when accounting for more than just vaccinated and unvaccinated (and then the columns of S etc. would also need to be expanded )
## hacky fix for now to only allocate a dose
n_vaccination_t_S[,] <- 0
n_vaccination_t_Ea[,] <- 0
n_vaccination_t_Eb[,] <- 0
n_vaccination_t_R[,] <- 0

n_vaccination_t_S[,1] <- min(
  floor((daily_doses_t[1,j] * S[i,1] * prioritisation_strategy[i,prioritisation_step]*vaccine_uptake[i])/sum(n_eligible[])),
  S[i,1])

n_vaccination_t_Ea[,1] <- min(
  floor((daily_doses_t[1,j] * Ea[i,1] * prioritisation_strategy[i,prioritisation_step]*vaccine_uptake[i])/sum(n_eligible[])),
  Ea[i,1])

n_vaccination_t_Eb[,1] <- min(
  floor((daily_doses_t[1,j] * Eb[i,1] * prioritisation_strategy[i,prioritisation_step]*vaccine_uptake[i])/sum(n_eligible[])),
  Eb[i,1])

n_vaccination_t_R[,1] <- min(
  floor((daily_doses_t[1,j] * R[i,1] * prioritisation_strategy[i,prioritisation_step]*vaccine_uptake[i])/sum(n_eligible[])),
  R[i,1])

## net vaccination change for relevant classes (S,Ea,Eb,R)
## logic here depends on vaccine class you are in (e.g. can only increase through j)
## this is generic so accounts for n_vax whatever it is
delta_S_n_vaccination[,] <- if(j==1) (-n_vaccination_t_S[i,j]) else if(j==n_vax) (n_vaccination_t_S[i,j-1]) else (-n_vaccination_t_S[i,j] + n_vaccination_t_S[i,j-1])

delta_Ea_n_vaccination[,] <- if(j==1) (-n_vaccination_t_Ea[i,j]) else if(j==n_vax) (n_vaccination_t_Ea[i,j-1]) else (-n_vaccination_t_Ea[i,j] + n_vaccination_t_Ea[i,j-1])

delta_Eb_n_vaccination[,] <- if(j==1) (-n_vaccination_t_Eb[i,j]) else if(j==n_vax) (n_vaccination_t_Eb[i,j-1]) else (-n_vaccination_t_Eb[i,j] + n_vaccination_t_Eb[i,j-1])

delta_R_n_vaccination[,] <- if(j==1) (-n_vaccination_t_R[i,j]) else if(j==n_vax) (n_vaccination_t_R[i,j-1]) else (-n_vaccination_t_R[i,j] + n_vaccination_t_R[i,j-1])

# use to test to see if vaccination working
update(vax_given_S) <- sum(n_vaccination_t_S[,])
update(vax_given_Ea) <- sum(n_vaccination_t_Ea[,])
update(vax_given_Eb) <- sum(n_vaccination_t_Eb[,])
update(vax_given_R) <- sum(n_vaccination_t_R[,])

## end of vaccination section


## Core equations for transitions between compartments:
# by age groups and vaccination class
# after vaccination has taken place
update(S[,]) <- S[i,j] + delta_S_n_vaccination[i,j] - n_SEa[i,j]
update(Ea[,]) <- Ea[i,j] + delta_Ea_n_vaccination[i,j] + delta_Ea[i,j]
update(Eb[,]) <- Eb[i,j] + delta_Eb_n_vaccination[i,j] + delta_Eb[i,j]
update(Ir[,]) <- Ir[i,j] + delta_Ir[i,j]
update(Id[,]) <- Id[i,j] + delta_Id[i,j]
update(R[,]) <- R[i,j] + delta_R_n_vaccination[i,j] + delta_R[i,j]
update(D[,]) <- D[i,j] + delta_D[i,j]

## Additional outputs
update(E[,]) <- Ea[i,j] + Eb[i,j]
update(I[,]) <- Ir[i,j] + Id[i,j]
update(N[,]) <- S[i,j] + Ea[i,j] + Eb[i,j] + Ir[i,j] + Id[i,j] + R[i,j] + D[i,j]
#update(N[,]) <- S[i,j] + E[i,j] + I[i,j] + R[i,j] + D[i,j]

is_same_week <- step %% steps_per_week > 0
update(cases) <- cases * is_same_week + sum(n_SEa[,])
update(cases_0_5) <- cases * is_same_week + sum(n_SEa[1,])
update(cases_5_15) <- cases * is_same_week + sum(n_SEa[2:3,])
update(cases_15_plus) <- cases * is_same_week + sum(n_SEa[4:16,])
update(cases_PBS) <- cases * is_same_week + sum(n_SEa[17,])
update(cases_SW) <- cases * is_same_week + sum(n_SEa[18,])

update(deaths) <- deaths * is_same_week + sum(n_IdD[,])
update(deaths_0_5) <- cases * is_same_week + sum(n_IdD[1,])
update(deaths_5_15) <- cases * is_same_week + sum(n_IdD[2:3,])
update(deaths_15_plus) <- cases * is_same_week + sum(n_IdD[4:16,])
update(deaths_PBS) <- cases * is_same_week + sum(n_IdD[17,])
update(deaths_SW) <- cases * is_same_week + sum(n_IdD[18,])

update(S_tot) <- sum(S[,])
update(E_tot) <- sum(E[,])
update(I_tot) <- sum(I[,])
update(R_tot) <- sum(R[,])
update(D_tot) <- sum(D[,])
update(N_tot) <- sum(N[,])
update(total_vax) <- total_vax + vax_given_S + vax_given_Ea + vax_given_Eb + vax_given_R


## Individual probabilities of transition:
p_SE[,] <- 1 - exp(-lambda[i,j] * dt) # S to E - age dependent
p_EE <- 1 - exp(-gamma_E * dt) # progression through latent period
p_EI <- 1 - exp(-gamma_E * dt) # progression to infection
p_IrR <- 1 - exp(-gamma_Ir * dt) # progression through infectious period to recovery
p_IdD <- 1 - exp(-gamma_Id * dt) # progression through infectious period to death

#Compute the force of infection

#  Mixing Matrix
m[,] <- user()
I_infectious[, ] <- I[i, j] * (1 - ve_T[j]) # I adjusted for reduced transmissibility

prop_infectious[] <- sum(I_infectious[i, ]) / sum(N[i, ])
# Generating Force of Infection
s_ij[, ] <- m[i, j] * prop_infectious[j] # for susceptible age i, % contacts infectious age j
lambda[,] <- ((beta_h * sum(s_ij[i, ])) + beta_z[i]) * (1-ve_I[j])

## Draws from binomial distributions for numbers changing between compartments accounting for vaccination:
n_SEa[,] <- rbinom(S[i,j] + delta_S_n_vaccination[i,j], p_SE[i,j])
n_EaEb[,] <- rbinom(Ea[i,j] +delta_Ea_n_vaccination[i,j], p_EE)
n_EbI[,] <- rbinom(Eb[i,j] +delta_Eb_n_vaccination[i,j], p_EI)
n_EbId[,] <- rbinom(n_EbI[i,j], CFR[i,j]) # Proportion of the infections that will die, impact of vaccination included already as part of inputs
n_EbIr[,] <- n_EbI[i,j] - n_EbId[i,j] # whatever infections don't die go to R

n_IrR[,] <- rbinom(Ir[i,j], p_IrR)
n_IdD[,] <- rbinom(Id[i,j], p_IdD)

## Calculate net change in each model state
delta_Ea[,] <- n_SEa[i,j] - n_EaEb[i,j]
delta_Eb[,] <- n_EaEb[i,j] - n_EbI[i,j]
delta_Ir[,] <- n_EbIr[i,j] - n_IrR[i,j]
delta_Id[,] <- n_EbId[i,j] - n_IdD[i,j]
delta_R[,] <- n_IrR[i,j]
delta_D[,] <- n_IdD[i,j]

## Initial states:
initial(S[,]) <- S0[i,j]
initial(Ea[,]) <- Ea0[i,j]
initial(Eb[,]) <- Eb0[i,j]
initial(Ir[,]) <- Ir0[i,j]
initial(Id[,]) <- Id0[i,j]
initial(R[,]) <- R0[i,j]
initial(D[,]) <- D0[i,j]

initial(E[,]) <- Ea0[i,j] + Eb0[i,j]
initial(I[,]) <- Ir0[i,j] + Id0[i,j]
initial(N[,]) <- S0[i,j] + Ea0[i,j] + Eb0[i,j] + Ir0[i,j] + Id0[i,j] + R0[i,j] + D0[i,j]
initial(cases) <- 0
initial(deaths) <- 0

initial(cases_0_5) <- 0
initial(cases_5_15) <- 0
initial(cases_15_plus) <- 0
initial(cases_PBS) <- 0
initial(cases_SW) <- 0
initial(deaths_0_5) <- 0
initial(deaths_5_15) <- 0
initial(deaths_15_plus) <- 0
initial(deaths_PBS) <- 0
initial(deaths_SW) <- 0

initial(vax_given_S) <- 0
initial(vax_given_Ea) <- 0
initial(vax_given_Eb) <- 0
initial(vax_given_R) <- 0

initial(S_tot) <- sum(S0[,])
initial(E_tot) <- sum(Ea0[,]) + sum(Eb0[,])
initial(I_tot) <- sum(Ir0[,]) + sum(Id0[,])
initial(R_tot) <- sum(R0[,])
initial(D_tot) <- sum(D0[,])
initial(N_tot) <- sum(S0[,]) + sum(Ea0[,]) + sum(Eb0[,]) + sum(Ir0[,]) +
  sum(Id0[,]) + sum(R0[,]) + sum(D0[,])
initial(total_vax) <- 0

##Initial vectors
S0[,] <- user()
Ea0[,] <- user()
Eb0[,] <- user()
Ir0[,] <- user()
Id0[,] <- user()
R0[,] <- user()
D0[,] <- user()

##Parameters
beta_h <- user()
beta_z[] <- user()
gamma_E <- user()
# gamma_I <- user()
gamma_Ir <- user()
gamma_Id <- user()
CFR[,] <- user()

#vaccine efficacy parameters
ve_T[] <- user()
ve_I[] <- user()
#ve_D[,] <- user() # this is included within the CFR

#Number of age classes & number of transmissibility classes
n_vax <- user()
n_group <- user()

##Dimensions of the different "vectors" here vectors stand for multi-dimensional arrays
dim(N) <- c(n_group,n_vax)
dim(S) <- c(n_group,n_vax)
dim(S0) <- c(n_group,n_vax)
dim(p_SE) <- c(n_group,n_vax)
dim(n_SEa) <- c(n_group,n_vax)

dim(Ea) <- c(n_group,n_vax)
dim(Ea0) <- c(n_group,n_vax)
dim(Eb0) <- c(n_group,n_vax)
dim(delta_Ea) <- c(n_group,n_vax)
dim(n_EaEb) <- c(n_group,n_vax)

dim(Eb) <- c(n_group,n_vax)
dim(delta_Eb) <- c(n_group,n_vax)
dim(n_EbI) <- c(n_group,n_vax)

dim(n_EbId) <- c(n_group,n_vax)
dim(n_EbIr) <- c(n_group,n_vax)
dim(E) <- c(n_group,n_vax)

dim(Ir0) <- c(n_group,n_vax)
dim(Ir) <- c(n_group,n_vax)
dim(delta_Ir) <- c(n_group,n_vax)
dim(n_IrR) <- c(n_group,n_vax)

dim(Id0) <- c(n_group,n_vax)
dim(Id) <- c(n_group,n_vax)
dim(delta_Id) <- c(n_group,n_vax)
dim(n_IdD) <- c(n_group,n_vax)
dim(I) <- c(n_group,n_vax)

dim(R) <- c(n_group,n_vax)
dim(R0) <- c(n_group,n_vax)
dim(delta_R) <- c(n_group,n_vax)

dim(D) <- c(n_group,n_vax)
dim(D0) <- c(n_group,n_vax)
dim(delta_D) <- c(n_group,n_vax)

dim(lambda) <- c(n_group,n_vax)
dim(m) <- c(n_group, n_group)
dim(I_infectious) <- c(n_group, n_vax)
dim(prop_infectious) <- c(n_group)
dim(s_ij) <- c(n_group,n_group)
dim(beta_z) <- c(n_group)

dim(CFR) <- c(n_group,n_vax)

dim(ve_T) <- c(n_vax)
dim(ve_I) <- c(n_vax)

vaccination_campaign_length <- user()
#dim(vaccination_campaign_length) <- 1L

#dim(n_vaccination) <- c(n_group,n_vax,vaccination_campaign_length) ##Lilith, the idea is that column j corresponds to the number of people leaving that class (e.g. going from unvaccinated to 1st dose). This makes the last column a dummy variable because it doesn't seem to like a dimension of (n_vax-1). And then dimension 3 corresponds to time
#dim(n_vaccination_t) <- c(n_group,n_vax)

dim(n_vaccination_t_S) <- c(n_group,n_vax)
dim(n_vaccination_t_Ea) <- c(n_group,n_vax)
dim(n_vaccination_t_Eb) <- c(n_group,n_vax)
dim(n_vaccination_t_R) <- c(n_group,n_vax)

dim(delta_S_n_vaccination) <- c(n_group,n_vax)
dim(delta_Ea_n_vaccination) <- c(n_group,n_vax)
dim(delta_Eb_n_vaccination) <- c(n_group,n_vax)
dim(delta_R_n_vaccination) <- c(n_group,n_vax)

#dim(n_eligible) <- c(n_group)

