# State indices: [i: age]
#restructure to follow help page for debugging
#https://mrc-ide.github.io/odin.dust/articles/sir_models.html

# Time Steps
dt <- user(1)
steps_per_week <- 7 / dt
initial(time) <- step
update(time) <- (step + 1) * dt
#output(time) <- TRUE

## help page has a "total" set of eqs here eg. not split by age group

## Core equations for transitions between compartments:
# by age groups
update(S[]) <- S[i] - n_SEa[i]
update(Ea[]) <- Ea[i] + delta_Ea[i]
update(Eb[]) <- Eb[i] + delta_Eb[i]
update(Ir[]) <- Ir[i] + delta_Ir[i]
update(Id[]) <- Id[i] + delta_Id[i]
update(R[]) <- R[i] + delta_R[i]
update(D[]) <- D[i] + delta_D[i]

## Additional outputs
update(E[]) <- Ea[i] + Eb[i]
update(I[]) <- Ir[i] + Id[i]
update(N[]) <- S[i] + Ea[i] + Eb[i] + Ir[i] + Id[i] + R[i] + D[i]

is_same_week <- step %% steps_per_week > 0
update(cases) <- cases * is_same_week + sum(n_SEa[])
update(deaths) <- deaths * is_same_week + sum(n_IdD[])

update(S_tot) <- sum(S[])
update(E_tot) <- sum(E[])
update(I_tot) <- sum(I[])
update(R_tot) <- sum(R[])
update(D_tot) <- sum(D[])
update(N_tot) <- sum(N[])

## Individual probabilities of transition:
p_SE[] <- 1 - exp(-lambda[i] * dt) # S to E - age dependent
p_EE <- 1 - exp(-gamma_E * dt) # progression through latent period
p_EI <- 1 - exp(-gamma_I * dt) # progression to infection
p_IrR <- 1 - exp(-gamma_Ir * dt) # progression through infectious period to recovery
p_IdD <- 1 - exp(-gamma_Id * dt) # progression through infectious period to death

#Compute the force of infection

#  Mixing Matrix
m[, ] <- user()

# Generating Force of Infection
s_ij[, ] <- m[i, j] * I[j] # for susceptible age i, % contacts infectious age j
lambda[] <- beta_h * sum(s_ij[i,]) + beta_z[i]

## Draws from binomial distributions for numbers changing between compartments:
n_SEa[] <- rbinom(S[i], p_SE[i])
n_EaEb[] <- rbinom(Ea[i], p_EE)
n_EbI[] <- rbinom(Eb[i], p_EI)
n_EbId[] <- rbinom(n_EbI[i], CFR[i]) # Proportion of the infections that will die
n_EbIr[] <- n_EbI[i] - n_EbId[i] # whatever infections don't die go to R

n_IrR[] <- rbinom(Ir[i], p_IrR)
n_IdD[] <- rbinom(Id[i], p_IdD)


## Calculate net change in each model state
delta_Ea[] <- n_SEa[i] - n_EaEb[i]
delta_Eb[] <- n_EaEb[i] - n_EbI[i]
delta_Ir[] <- n_EbIr[i] - n_IrR[i]
delta_Id[] <- n_EbId[i] - n_IdD[i]
delta_R[] <- n_IrR[i]
delta_D[] <- n_IdD[i]

## Initial states:
initial(S[]) <- S0[i]
initial(Ea[]) <- Ea0[i]
initial(Eb[]) <- Eb0[i]
initial(Ir[]) <- Ir0[i]
initial(Id[]) <- Id0[i]
initial(R[]) <- R0[i]
initial(D[]) <- D0[i]

initial(E[]) <- Ea0[i] + Eb0[i]
initial(I[]) <- Ir0[i] + Id0[i]
initial(N[]) <- S0[i] + Ea0[i] + Eb0[i] + Ir0[i] + Id0[i] + R0[i] + D0[i]
initial(cases) <- 0
initial(deaths) <- 0

initial(S_tot) <- sum(S0[])
initial(E_tot) <- sum(Ea0[]) + sum(Eb0[])
initial(I_tot) <- sum(Ir0[]) + sum(Id0[])
initial(R_tot) <- sum(R0[])
initial(D_tot) <- sum(D0[])
initial(N_tot) <- sum(S0[]) + sum(Ea0[]) + sum(Eb0[]) + sum(Ir0[]) +
  sum(Id0[]) + sum(R0[]) + sum(D0[])

##Initial vectors
S0[] <- user()
Ea0[] <- user()
Eb0[] <- user()
Ir0[] <- user()
Id0[] <- user()
R0[] <- user()
D0[] <- user()

##Parameters
beta_h <- user()
beta_z[] <- user()
gamma_E <- user()
gamma_I <- user()
gamma_Ir <- user()
gamma_Id <- user()
CFR[] <- user() #CFR

#Number of age classes & number of transmissibility classes
n_group <- user()

##Dimensions of the different "vectors" here vectors stand for multi-dimensional arrays
dim(N) <- n_group
dim(S) <- n_group
dim(S0) <- n_group
dim(p_SE) <- n_group
dim(n_SEa) <- n_group

dim(Ea) <- c(n_group)
dim(Ea0) <- c(n_group)
dim(Eb0) <- c(n_group)
dim(delta_Ea) <- c(n_group)
dim(n_EaEb) <- c(n_group)

dim(Eb) <- c(n_group)
dim(delta_Eb) <- c(n_group)
dim(n_EbI) <- c(n_group)

dim(n_EbId) <- c(n_group)
dim(n_EbIr) <- c(n_group)
dim(E) <- n_group

dim(Ir0) <- c(n_group)
dim(Ir) <- c(n_group)
dim(delta_Ir) <- c(n_group)
dim(n_IrR) <- c(n_group)

dim(Id0) <- c(n_group)
dim(Id) <- c(n_group)
dim(delta_Id) <- c(n_group)
dim(n_IdD) <- c(n_group)
dim(I) <- n_group

dim(R) <- c(n_group)
dim(R0) <- c(n_group)
dim(delta_R) <- c(n_group)

dim(D) <- c(n_group)
dim(D0) <- c(n_group)
dim(delta_D) <- c(n_group)

dim(lambda) <- n_group
dim(m) <- c(n_group, n_group)
dim(s_ij) <- c(n_group,n_group)
dim(beta_z) <- n_group

dim(CFR) <- c(n_group)


