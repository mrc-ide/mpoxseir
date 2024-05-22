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
update(S[]) <- S[i] - n_SE1[i]
update(E1[]) <- E1[i] + delta_E1[i]
update(E2[]) <- E2[i] + delta_E2[i]
update(Ir[]) <- Ir[i] + delta_Ir[i]
update(Id[]) <- Id[i] + delta_Id[i]
update(R[]) <- R[i] + delta_R[i]
update(D[]) <- D[i] + delta_D[i]

## Additional outputs
update(E[]) <- E1[i] + E2[i]
update(I[]) <- Ir[i] + Id[i]
update(N[]) <- S[i] + E1[i] + E2[i] + Ir[i] + Id[i] + R[i] + D[i]

is_same_week <- step %% steps_per_week > 0
update(cases) <- cases * is_same_week + sum(n_SE1[])
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
lambda[] <- beta_h * sum(s_ij[i,]) + beta_z

## Draws from binomial distributions for numbers changing between compartments:
n_SE1[] <- rbinom(S[i], p_SE[i])
n_E1E2[] <- rbinom(E1[i], p_EE)
n_E2I[] <- rbinom(E2[i], p_EI)
n_E2Id[] <- rbinom(n_E2I[i], CFR[i]) # Proportion of the infections that will die
n_E2Ir[] <- n_E2I[i] - n_E2Id[i] # whatever infections don't die go to R

n_IrR[] <- rbinom(Ir[i], p_IrR)
n_IdD[] <- rbinom(Id[i], p_IdD)


## Calculate net change in each model state
delta_E1[] <- n_SE1[i] - n_E1E2[i]
delta_E2[] <- n_E1E2[i] - n_E2I[i]
delta_Ir[] <- n_E2Ir[i] - n_IrR[i]
delta_Id[] <- n_E2Id[i] - n_IdD[i]
delta_R[] <- n_IrR[i]
delta_D[] <- n_IdD[i]

## Initial states:
initial(S[]) <- S0[i]
initial(E1[]) <- E10[i]
initial(E2[]) <- E20[i]
initial(Ir[]) <- Ir0[i]
initial(Id[]) <- Id0[i]
initial(R[]) <- R0[i]
initial(D[]) <- D0[i]

initial(E[]) <- E10[i] + E20[i]
initial(I[]) <- Ir0[i] + Id0[i]
initial(N[]) <- S0[i] + E10[i] + E20[i] + Ir0[i] + Id0[i] + R0[i] + D0[i]
initial(cases) <- 0
initial(deaths) <- 0

initial(S_tot) <- sum(S0[])
initial(E_tot) <- sum(E10[]) + sum(E20[])
initial(I_tot) <- sum(Ir0[]) + sum(Id0[])
initial(R_tot) <- sum(R0[])
initial(D_tot) <- sum(D0[])
initial(N_tot) <- sum(S0[]) + sum(E10[]) + sum(E20[]) + sum(Ir0[]) +
  sum(Id0[]) + sum(R0[]) + sum(D0[])

##Initial vectors
S0[] <- user()
E10[] <- user()
E20[] <- user()
Ir0[] <- user()
Id0[] <- user()
R0[] <- user()
D0[] <- user()

##Parameters
beta_h <- user()
beta_z <- user()
gamma_E <- user()
gamma_I <- user()
gamma_Ir <- user()
gamma_Id <- user()
CFR[] <- user() #CFR

#Number of age classes & number of transmissibility classes
N_age <- user()

##Dimensions of the different "vectors" here vectors stand for multi-dimensional arrays
dim(N) <- N_age
dim(S) <- N_age
dim(S0) <- N_age
dim(p_SE) <- N_age
dim(n_SE1) <- N_age

dim(E1) <- c(N_age)
dim(E10) <- c(N_age)
dim(E20) <- c(N_age)
dim(delta_E1) <- c(N_age)
dim(n_E1E2) <- c(N_age)

dim(E2) <- c(N_age)
dim(delta_E2) <- c(N_age)
dim(n_E2I) <- c(N_age)

dim(n_E2Id) <- c(N_age)
dim(n_E2Ir) <- c(N_age)
dim(E) <- N_age

dim(Ir0) <- c(N_age)
dim(Ir) <- c(N_age)
dim(delta_Ir) <- c(N_age)
dim(n_IrR) <- c(N_age)

dim(Id0) <- c(N_age)
dim(Id) <- c(N_age)
dim(delta_Id) <- c(N_age)
dim(n_IdD) <- c(N_age)
dim(I) <- N_age

dim(R) <- c(N_age)
dim(R0) <- c(N_age)
dim(delta_R) <- c(N_age)

dim(D) <- c(N_age)
dim(D0) <- c(N_age)
dim(delta_D) <- c(N_age)

dim(lambda) <- N_age
dim(m) <- c(N_age, N_age)
dim(s_ij) <- c(N_age,N_age)

dim(CFR) <- c(N_age)


