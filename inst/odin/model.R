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
# at a given prespecified time we give N number of doses - Lilith see note on dimensionality please towards end of script 
# time working solution: we set up a big array for every time point that is defaulted to 0 and then we just change the entries that we want to, so this means that this step takes place every time. But this could be very inefficient so can replace in future if needed?

n_vaccination_allocation_SER[] <- user() #det calc of where the vaccines go to in terms of S,Ea,Eb,R, but could also be used as vector of probs for a multinomial draw - will need to be updated when do Erlangs

## code adapted from https://github.com/mrc-ide/gasworks/blob/ce197d380e91bde536835cb9da3b09df784bb64f/inst/odin/model.R#L11 re vaccination
## with this we need to ensure n_vaccination final time step is all 0 - done outside of the model in pre-processing
n_vaccination_t[,] <- if (as.integer(step) >= (vaccination_campaign_length))
  n_vaccination[i,j,vaccination_campaign_length] else n_vaccination[i,j,step]

## need to check that we have the capacity to vaccinate according to the schedule in n_vaccination_t
## first see what the allocation is across the classes 
n_vaccination_t_S[,] <- round(n_vaccination_t[i,j]*n_vaccination_allocation_SER[1])
n_vaccination_t_Ea[,] <- round(n_vaccination_t[i,j]*n_vaccination_allocation_SER[2])
n_vaccination_t_Eb[,] <- round(n_vaccination_t[i,j]*n_vaccination_allocation_SER[3])
n_vaccination_t_R[,] <- round(n_vaccination_t[i,j]*n_vaccination_allocation_SER[4])

## set equal to number of person in class if the proposed vaccination number is higher
n_vaccination_t_S[,] <- if(abs(n_vaccination_t_S[i,j])>=S[i,j]) S[i,j] else n_vaccination_t_S[i,j]
n_vaccination_t_Ea[,] <- if(abs(n_vaccination_t_Ea[i,j])>=Ea[i,j]) Ea[i,j] else n_vaccination_t_Ea[i,j]
n_vaccination_t_Eb[,] <- if(abs(n_vaccination_t_Eb[i,j])>=Eb[i,j]) Eb[i,j] else n_vaccination_t_Eb[i,j]
n_vaccination_t_R[,] <- if(abs(n_vaccination_t_R[i,j])>=R[i,j]) R[i,j] else n_vaccination_t_R[i,j]

## calculate net vaccination change for relevant classes (S,Ea,Eb,R)
## logic here depends on vaccine class you are in (e.g. can only increase through j)
delta_S_n_vaccination[,] <- if(j==1) (n_vaccination_t_S[i,j]) else if(j==n_vax) (n_vaccination_t_S[i,j-1]) else (n_vaccination_t_S[i,j] + n_vaccination_t_S[i,j-1]) 

delta_Ea_n_vaccination[,] <- if(j==1) (n_vaccination_t_Ea[i,j]) else if(j==n_vax) (n_vaccination_t_Ea[i,j-1]) else (n_vaccination_t_Ea[i,j] + n_vaccination_t_Ea[i,j-1]) 

delta_Eb_n_vaccination[,] <- if(j==1) (n_vaccination_t_Eb[i,j]) else if(j==n_vax) (n_vaccination_t_Eb[i,j-1]) else (n_vaccination_t_Eb[i,j] + n_vaccination_t_Eb[i,j-1]) 

delta_R_n_vaccination[,] <- if(j==1) (n_vaccination_t_R[i,j]) else if(j==n_vax) (n_vaccination_t_R[i,j-1]) else (n_vaccination_t_R[i,j] + n_vaccination_t_R[i,j-1]) 


# ## calculate net vaccination change for relevant classes (S,Ea,Eb,R)
# ## logic here depends on vaccine class you are in (e.g. can only increase through j)
# delta_S_n_vaccination[,] <- if(j==1) (-round(n_vaccination_t[i,j]*n_vaccination_allocation_SER[1])) else if(j==n_vax) (round(n_vaccination_t[i,j-1]*n_vaccination_allocation_SER[1])) else (-round(n_vaccination_t[i,j]*n_vaccination_allocation_SER[1]) + round(n_vaccination_t[i,j-1]*n_vaccination_allocation_SER[1])) 
# 
# delta_Ea_n_vaccination[,] <- if(j==1) (-round(n_vaccination_t[i,j]*n_vaccination_allocation_SER[2])) else if(j==n_vax) (round(n_vaccination_t[i,j-1]*n_vaccination_allocation_SER[2])) else (-round(n_vaccination_t[i,j]*n_vaccination_allocation_SER[2]) + round(n_vaccination_t[i,j-1]*n_vaccination_allocation_SER[2]))
# 
# delta_Eb_n_vaccination[,] <- if(j==1) (-round(n_vaccination_t[i,j]*n_vaccination_allocation_SER[3])) else if(j==n_vax) (round(n_vaccination_t[i,j-1]*n_vaccination_allocation_SER[3])) else (-round(n_vaccination_t[i,j]*n_vaccination_allocation_SER[3]) + round(n_vaccination_t[i,j-1]*n_vaccination_allocation_SER[3])) 
# 
# delta_R_n_vaccination[,] <- if(j==1) (-round(n_vaccination_t[i,j]*n_vaccination_allocation_SER[4])) else if(j==n_vax) (round(n_vaccination_t[i,j-1]*n_vaccination_allocation_SER[4])) else (-round(n_vaccination_t[i,j]*n_vaccination_allocation_SER[4]) + round(n_vaccination_t[i,j-1]*n_vaccination_allocation_SER[4]))

# ## Logic in case number to be vaccinated exceeds number of people in compartment: set to the number of people in the compartment if planned vaccination exceeds this number 
# delta_S_n_vaccination[,] <- if(delta_S_n_vaccination[i,j]>=S[i,j]) S[i,j] else delta_S_n_vaccination[i,j]
# 
# delta_Ea_n_vaccination[,] <- if(delta_Ea_n_vaccination[i,j]>=Ea[i,j]) Ea[i,j] else delta_Ea_n_vaccination[i,j]
# 
# delta_Eb_n_vaccination[,] <- if(delta_Eb_n_vaccination[i,j]>=Eb[i,j]) Eb[i,j] else delta_Eb_n_vaccination[i,j]
# 
# delta_R_n_vaccination[,] <- if(delta_R_n_vaccination[i,j]>=R[i,j]) R[i,j] else delta_R_n_vaccination[i,j]


## update states to reflect vaccination
update(S_after_vax[,]) <- S[i,j] + delta_S_n_vaccination[i,j]
update(Ea_after_vax[,]) <- Ea[i,j] + delta_Ea_n_vaccination[i,j]
update(Eb_after_vax[,]) <- Eb[i,j] + delta_Eb_n_vaccination[i,j]
update(R_after_vax[,]) <- R[i,j] + delta_R_n_vaccination[i,j]

## end of vaccination section

## Core equations for transitions between compartments:
# by age groups and vaccination class 
# after vaccination has taken place
update(S[,]) <- S_after_vax[i,j] - n_SEa[i,j]
update(Ea[,]) <- Ea_after_vax[i,j] + delta_Ea[i,j]
update(Eb[,]) <- Eb_after_vax[i,j] + delta_Eb[i,j]
update(Ir[,]) <- Ir[i,j] + delta_Ir[i,j]
update(Id[,]) <- Id[i,j] + delta_Id[i,j]
update(R[,]) <- R_after_vax[i,j] + delta_R[i,j]
update(D[,]) <- D[i,j] + delta_D[i,j]

## Additional outputs
update(E[,]) <- Ea[i,j] + Eb[i,j]
update(I[,]) <- Ir[i,j] + Id[i,j]
update(N[,]) <- S[i,j] + Ea[i,j] + Eb[i,j] + Ir[i,j] + Id[i,j] + R[i,j] + D[i,j]

is_same_week <- step %% steps_per_week > 0
update(cases) <- cases * is_same_week + sum(n_SEa[,])
update(deaths) <- deaths * is_same_week + sum(n_IdD[,])

update(S_tot) <- sum(S[,])
update(E_tot) <- sum(E[,])
update(I_tot) <- sum(I[,])
update(R_tot) <- sum(R[,])
update(D_tot) <- sum(D[,])
update(N_tot) <- sum(N[,])

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

## Draws from binomial distributions for numbers changing between compartments:
n_SEa[,] <- rbinom(S_after_vax[i,j], (1-ve_I[i,j])*p_SE[i,j])
n_EaEb[,] <- rbinom(Ea_after_vax[i,j], p_EE)
n_EbI[,] <- rbinom(Eb_after_vax[i,j], p_EI)
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

initial(S_after_vax[,]) <- S0[i,j]
initial(Ea_after_vax[,]) <- Ea0[i,j]
initial(Eb_after_vax[,]) <- Eb0[i,j]
initial(R_after_vax[,]) <- R0[i,j]

initial(E[,]) <- Ea0[i,j] + Eb0[i,j]
initial(I[,]) <- Ir0[i,j] + Id0[i,j]
initial(N[,]) <- S0[i,j] + Ea0[i,j] + Eb0[i,j] + Ir0[i,j] + Id0[i,j] + R0[i,j] + D0[i,j]
initial(cases) <- 0
initial(deaths) <- 0

initial(S_tot) <- sum(S0[,])
initial(E_tot) <- sum(Ea0[,]) + sum(Eb0[,])
initial(I_tot) <- sum(Ir0[,]) + sum(Id0[,])
initial(R_tot) <- sum(R0[,])
initial(D_tot) <- sum(D0[,])
initial(N_tot) <- sum(S0[,]) + sum(Ea0[,]) + sum(Eb0[,]) + sum(Ir0[,]) +
  sum(Id0[,]) + sum(R0[,]) + sum(D0[,])

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
dim(S_after_vax) <- c(n_group,n_vax)
dim(S0) <- c(n_group,n_vax)
dim(p_SE) <- c(n_group,n_vax)
dim(n_SEa) <- c(n_group,n_vax)

dim(Ea) <- c(n_group,n_vax)
dim(Ea_after_vax) <- c(n_group,n_vax)
dim(Ea0) <- c(n_group,n_vax)
dim(Eb0) <- c(n_group,n_vax)
dim(delta_Ea) <- c(n_group,n_vax)
dim(n_EaEb) <- c(n_group,n_vax)

dim(Eb) <- c(n_group,n_vax)
dim(Eb_after_vax) <- c(n_group,n_vax)
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
dim(R_after_vax) <- c(n_group,n_vax)
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
