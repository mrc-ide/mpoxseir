# State indices: [i: age; j: vaccination]
#restructure to follow help page for debugging
#https://mrc-ide.github.io/odin.dust/articles/sir_models.html

# Time Steps
dt <- user(1)
steps_per_week <- 7 / dt
initial(time) <- step
update(time) <- (step + 1) * dt
#output(time) <- TRUE

#### Vaccination

## specify the number of vaccinations that happen daily per age group and vaccine strata (n_vax: vaccine strata, n_vaccination: new vaccinations given where column j represents individuals leaving that strata)
n_vaccination[,,] <- user()

n_vaccination_allocation_SER[] <- user() #det calc of where the vaccines go to in terms of S,Ea,Eb,R, but could also be used as vector of probs for a multinomial draw - will need to be updated when do Erlangs

## with this we need to ensure n_vaccination final time step is all 0 - done outside of the model in pre-processing
n_vaccination_t[,] <- if (as.integer(time) >= (vaccination_campaign_length))
  n_vaccination[i,j,vaccination_campaign_length] else n_vaccination[i,j,time]
n_eligible[] <- S[i, 1] + Ea[i, 1] + Eb[i, 1] + R[i, 1]
## first allocate across the classes 
n_vaccination_t_S[,]   <- min(round(n_vaccination_t[i,j] * S[i, 1] / n_eligible[i]), S[i, 1])
n_vaccination_t_Ea[,] <- min(round(n_vaccination_t[i,j] * Ea[i, 1] / n_eligible[i]), Ea[i, 1])
n_vaccination_t_Eb[,] <- min(round(n_vaccination_t[i,j] * Eb[i, 1] / n_eligible[i]), Eb[i, 1])
n_vaccination_t_R[,]   <- min(round(n_vaccination_t[i,j] * R[i, 1] / n_eligible[i]), R[i, 1])

## need to check that we have the capacity to vaccinate according to the schedule in n_vaccination_t
## set equal to number of person in class if the proposed vaccination number is higher
n_vaccination_t_S[,] <- if(n_vaccination_t_S[i,j]>=S[i,j]) S[i,j] else n_vaccination_t_S[i,j]
n_vaccination_t_Ea[,] <- if(n_vaccination_t_Ea[i,j]>=Ea[i,j]) Ea[i,j] else n_vaccination_t_Ea[i,j]
n_vaccination_t_Eb[,] <- if(n_vaccination_t_Eb[i,j]>=Eb[i,j]) Eb[i,j] else n_vaccination_t_Eb[i,j]
n_vaccination_t_R[,] <- if(n_vaccination_t_R[i,j]>=R[i,j]) R[i,j] else n_vaccination_t_R[i,j]

## calculate net vaccination change for relevant classes (S,Ea,Eb,R)
## logic here depends on vaccine class you are in (e.g. can only increase through j)
delta_S_n_vaccination[,] <- if(j==1) (-n_vaccination_t_S[i,j]) else if(j==n_vax) (n_vaccination_t_S[i,j-1]) else (-n_vaccination_t_S[i,j] + n_vaccination_t_S[i,j-1]) 

delta_Ea_n_vaccination[,] <- if(j==1) (-n_vaccination_t_Ea[i,j]) else if(j==n_vax) (n_vaccination_t_Ea[i,j-1]) else (-n_vaccination_t_Ea[i,j] + n_vaccination_t_Ea[i,j-1]) 

delta_Eb_n_vaccination[,] <- if(j==1) (-n_vaccination_t_Eb[i,j]) else if(j==n_vax) (n_vaccination_t_Eb[i,j-1]) else (-n_vaccination_t_Eb[i,j] + n_vaccination_t_Eb[i,j-1]) 

delta_R_n_vaccination[,] <- if(j==1) (-n_vaccination_t_R[i,j]) else if(j==n_vax) (n_vaccination_t_R[i,j-1]) else (-n_vaccination_t_R[i,j] + n_vaccination_t_R[i,j-1]) 


## update states to reflect vaccination
## this is what was causing issues so don't have as a separate variable now
# update(S_after_vax[,]) <- S[i,j] + delta_S_n_vaccination[i,j]
# update(Ea_after_vax[,]) <- Ea[i,j] + delta_Ea_n_vaccination[i,j]
# update(Eb_after_vax[,]) <- Eb[i,j] + delta_Eb_n_vaccination[i,j]
# update(R_after_vax[,]) <- R[i,j] + delta_R_n_vaccination[i,j]

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
update(deaths) <- deaths * is_same_week + sum(n_IdD[,])

update(S_tot) <- sum(S[,])
update(E_tot) <- sum(E[,])
update(I_tot) <- sum(I[,])
update(R_tot) <- sum(R[,])
update(D_tot) <- sum(D[,])
update(N_tot) <- sum(N[,])

# update(S_after_vax_tot) <- sum(S_after_vax[,])
# update(Ea_after_vax_tot) <- sum(Ea_after_vax[,])
# update(Eb_after_vax_tot) <- sum(Eb_after_vax[,])
# update(R_after_vax_tot) <- sum(R_after_vax[,])

# ## troubleshooting
# update(S_tot_test) <- sum(S[,]) + sum(delta_S_n_vaccination[,]) - sum(n_SEa[,])
# update(Ea_tot_test) <- sum(Ea[,]) + sum(delta_Ea_n_vaccination[,]) + sum(delta_Ea[,])
# update(Eb_tot_test) <- sum(Eb[,]) + sum(delta_Eb_n_vaccination[,]) + sum(delta_Eb[,])
# update(Ir_tot_test) <- sum(Ir[,]) + sum(delta_Ir[,])
# update(Id_tot_test) <- sum(Id[,]) + sum(delta_Id[,])
# update(R_tot_test) <- sum(R[,]) + sum(delta_R_n_vaccination[,]) + sum(delta_R[,])
# update(D_tot_test) <- sum(D[,])  + sum(delta_D[,])
# 
# update(N_tot_test) <- (S_tot_test + Ea_tot_test + Eb_tot_test + Ir_tot_test + Id_tot_test + R_tot_test + D_tot_test) 



## Individual probabilities of transition:
p_SE[,] <- 1 - exp(-lambda[i,j] * dt) # S to E - age dependent
p_EE <- 1 - exp(-gamma_E * dt) # progression through latent period
p_EI <- 1 - exp(-gamma_I * dt) # progression to infection
p_IrR <- 1 - exp(-gamma_Ir * dt) # progression through infectious period to recovery
p_IdD <- 1 - exp(-gamma_Id * dt) # progression through infectious period to death

#Compute the force of infection

#  Mixing Matrix
m[, ,] <- user()

# Generating Force of Infection
s_ijk[, ,] <- m[i,j,k] * I[j,k] # for susceptible age i, % contacts infectious age k
# does zoonotic spillover need to be stratified by vaccination as well?
lambda[,] <- (beta_h * sum(s_ijk[i, j, ]) * (1 - ve_T[i,j]))+ beta_z[i]

## Draws from binomial distributions for numbers changing between compartments accounting for vaccination:
n_SEa[,] <- rbinom(S[i,j] + delta_S_n_vaccination[i,j], (1-ve_I[i,j])*p_SE[i,j])
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

# initial(S_after_vax[,]) <- 0
# initial(Ea_after_vax[,]) <- 0
# initial(Eb_after_vax[,]) <- 0
# initial(R_after_vax[,]) <- 0

initial(E[,]) <- Ea0[i,j] + Eb0[i,j]
initial(I[,]) <- Ir0[i,j] + Id0[i,j]
initial(N[,]) <- S0[i,j] + Ea0[i,j] + Eb0[i,j] + Ir0[i,j] + Id0[i,j] + R0[i,j] + D0[i,j]
initial(cases) <- 0
initial(deaths) <- 0

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

# initial(S_after_vax_tot) <- 0
# initial(Ea_after_vax_tot) <- 0
# initial(Eb_after_vax_tot) <- 0
# initial(R_after_vax_tot) <- 0

# initial(S_tot_test) <- 0
# initial(Ea_tot_test) <- 0
# initial(Eb_tot_test) <- 0
# initial(Ir_tot_test) <- 0
# initial(Id_tot_test) <- 0
# initial(R_tot_test) <- 0
# initial(D_tot_test) <- 0
# initial(N_tot_test) <- 0


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
gamma_I <- user()
gamma_Ir <- user()
gamma_Id <- user()
CFR[,] <- user() 
#vaccine efficacy parameters
ve_T[,] <- user()
ve_I[,] <- user()
#ve_D[,] <- user() # this is included within the CFR

#Number of age classes & number of transmissibility classes
n_vax <- user()
n_group <- user() 

##Dimensions of the different "vectors" here vectors stand for multi-dimensional arrays
dim(N) <- c(n_group,n_vax) 
dim(S) <- c(n_group,n_vax)
#dim(S_after_vax) <- c(n_group,n_vax)
dim(S0) <- c(n_group,n_vax)
dim(p_SE) <- c(n_group,n_vax)
dim(n_SEa) <- c(n_group,n_vax)

dim(Ea) <- c(n_group,n_vax)
#dim(Ea_after_vax) <- c(n_group,n_vax)
dim(Ea0) <- c(n_group,n_vax)
dim(Eb0) <- c(n_group,n_vax)
dim(delta_Ea) <- c(n_group,n_vax)
dim(n_EaEb) <- c(n_group,n_vax)

dim(Eb) <- c(n_group,n_vax)
#dim(Eb_after_vax) <- c(n_group,n_vax)
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
#dim(R_after_vax) <- c(n_group,n_vax)
dim(R0) <- c(n_group,n_vax)
dim(delta_R) <- c(n_group,n_vax)

dim(D) <- c(n_group,n_vax)
dim(D0) <- c(n_group,n_vax)
dim(delta_D) <- c(n_group,n_vax)

dim(lambda) <- c(n_group,n_vax)
dim(m) <- c(n_group, n_group,n_vax) 
dim(s_ijk) <- c(n_group,n_group,n_vax)
dim(beta_z) <- c(n_group)

dim(CFR) <- c(n_group,n_vax)

dim(ve_T) <- c(n_group,n_vax)
dim(ve_I) <- c(n_group,n_vax)

vaccination_campaign_length <- user()
#dim(vaccination_campaign_length) <- 1L

dim(n_vaccination) <- c(n_group,n_vax,vaccination_campaign_length) ##Lilith, the idea is that column j corresponds to the number of people leaving that class (e.g. going from unvaccinated to 1st dose). This makes the last column a dummy variable because it doesn't seem to like a dimension of (n_vax-1). And then dimension 3 corresponds to time
dim(n_vaccination_t) <- c(n_group,n_vax)

dim(n_vaccination_t_S) <- c(n_group,n_vax)
dim(n_vaccination_t_Ea) <- c(n_group,n_vax)
dim(n_vaccination_t_Eb) <- c(n_group,n_vax)
dim(n_vaccination_t_R) <- c(n_group,n_vax)

dim(delta_S_n_vaccination) <- c(n_group,n_vax) 
dim(delta_Ea_n_vaccination) <- c(n_group,n_vax) 
dim(delta_Eb_n_vaccination) <- c(n_group,n_vax) 
dim(delta_R_n_vaccination) <- c(n_group,n_vax) 

dim(n_vaccination_allocation_SER) <- 4L # need to be updated when implement Erlangs
