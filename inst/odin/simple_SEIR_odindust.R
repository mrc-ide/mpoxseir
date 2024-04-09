#E and R stage indexed by i andj with
#i for the age group
#j for the progression (not exponential latent and infectious period)
#write in odin.dust
#restructure to follow help page for debugging
#https://mrc-ide.github.io/odin.dust/articles/sir_models.html

# Time Steps
dt <- user()
initial(time) <- 0
update(time) <- (step + 1) * dt
#output(time) <- TRUE

## help page has a "total" set of eqs here eg. not split by age group

## Core equations for transitions between compartments:
# by age groups
update(S[]) <- S[i] - n_SE[i]
update(E1[]) <- E1[i] + delta_E1[i]
update(E2[]) <- E2[i] + delta_E2[i]
update(I[]) <- I[i] + delta_I[i]
update(R[]) <- R[i] + delta_R[i]

## help page has cumulative incidence here 

## Individual probabilities of transition:
p_SE[] <- 1 - exp(-lambda[i] * dt) # S to E - age dependent
p_EE <- 1 - exp(-gamma_E * dt) # progression through latent period
p_EI <- 1 - exp(-gamma_E * dt) # progression through latent period
p_IR <- 1 - exp(-gamma_I * dt) # progression through infectious period

#Compute the force of infection

# Interpolation for Mixing Matrix
# get rid of interpolation as limited support for it
m[, ] <- user()
#m[, ] <- interpolate(tt_matrix, mix_mat_set, "constant")
#dim(m) <- c(N_age, N_age)
#tt_matrix[] <- user()
#mix_mat_set[, ,] <- user()
#dim(tt_matrix) <- user()
#dim(mix_mat_set) <- c(length(tt_matrix), N_age, N_age)

# Interpolation for beta
# get rid of interpolation as limited support for it
# beta <- interpolate(tt_beta, beta_set, "constant")
# tt_beta[] <- user()
# beta_set[] <- user()
# dim(tt_beta) <- user()
# dim(beta_set) <- length(tt_beta)

# Generating Force of Infection
#temp[] <- I[i]
s_ij[, ] <- m[i, j] * I[i] # why temp[j] not temp[i] RM corrected as with help page 9/4
lambda[] <- beta * sum(s_ij[i,]) ## in help page column i not row i


## Draws from binomial distributions for numbers changing between compartments:
n_SE[] <- rbinom(S[i], p_SE[i])
n_EE[] <- rbinom(E1[i], p_EE)
n_EI[] <- rbinom(E2[i], p_EI)
#output(n_EI[]) <- TRUE
n_IR[] <- rbinom(I[i], p_IR)


## these are just changes to the number in each state?
delta_E1[] <- n_SE[i] - n_EE[i]
delta_E2[] <- n_EE[i] - n_EI[i]
delta_I[] <- n_EI[i] - n_IR[i]
delta_R[] <- n_IR[i]

## Initial states:
initial(S[]) <- S0[i]
initial(E1[]) <- E0[i]
initial(E2[]) <- E02[i]
initial(I[]) <- I0[i]
initial(R[]) <- R0[i]

##Initial vectors
S0[] <- user()
E0[] <- user()
E02[] <- user()
I0[] <- user()
R0[] <- user()

##Parameters
beta <- user()
gamma_E <- user()
gamma_I <- user()


#Number of age classes & number of transmissibility classes
N_age <- user()


##Dimensions of the different "vectors" here vectors stand for multi-dimensional arrays
dim(S) <- N_age
dim(S0) <- N_age
dim(p_SE) <- N_age
dim(n_SE) <- N_age

dim(E1) <- c(N_age)
dim(E0) <- c(N_age)
dim(E02) <- c(N_age)
dim(delta_E1) <- c(N_age)
dim(n_EE) <- c(N_age)

dim(E2) <- c(N_age)
dim(delta_E2) <- c(N_age)
dim(n_EI) <- c(N_age)

dim(I0) <- c(N_age)
dim(I) <- c(N_age)
dim(delta_I) <- c(N_age)
dim(n_IR) <- c(N_age)

dim(R) <- c(N_age)
dim(R0) <- c(N_age)
dim(delta_R) <- c(N_age)

dim(lambda) <- N_age
dim(m) <- c(N_age, N_age)
dim(s_ij) <- c(N_age,N_age)
#dim(temp) <- N_age
