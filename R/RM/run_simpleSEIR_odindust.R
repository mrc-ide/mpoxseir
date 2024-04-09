## running odin.dust

gen_sir <- odin.dust::odin_dust("C:/Users/rom116/OneDrive - Imperial College London/mpox/Modelling/sir-demo.R")

## to run in package is a different thing 

sir_model <- gen_sir$new(pars = list(dt = 1,
                                     S_ini = 1000,
                                     I_ini = 10,
                                     beta = 0.2,
                                     gamma = 0.1),
                         time = 1,
                         n_particles = 10L,
                         n_threads = 4L,
                         seed = 1L)

sir_model$state()
sir_model$run(10)
sir_model$run(20)
sir_model$run(100)

dt <- 0.25
n_particles <- 10L
p_new <- list(dt = dt, S_ini = 2000, I_ini = 10, beta = 0.4, gamma = 0.1)
sir_model$update_state(pars = p_new, time = 0)
#sir_model$run(100)
sir_model$state()


n_times <- 200
x <- array(NA, dim = c(sir_model$info()$len, n_particles, n_times))

for (t in seq_len(n_times)) {
  x[ , , t] <- sir_model$run(t)
}
time <- x[1, 1, ]
x <- x[-1, , ]

par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
cols <- c(S = "#8c8cd9", I = "#cc0044", R = "#999966")
matplot(time, t(x[1, , ]), type = "l",
        xlab = "Time", ylab = "Number of individuals",
        col = cols[["S"]], lty = 1, ylim = range(x))
matlines(time, t(x[2, , ]), col = cols[["I"]], lty = 1)
matlines(time, t(x[3, , ]), col = cols[["R"]], lty = 1)
legend("left", lwd = 1, col = cols, legend = names(cols), bty = "n")


# run_simple_SEEIR_model in odin.dust

# test <- odin_dust("inst/odin/simple_SEIR.R") # doesn't like output
# test <- odin_dust("C:/Users/rom116/OneDrive - Imperial College London/mpox/Modelling/sir-age-demo.R") # compiled

first_attempt <- odin.dust::odin_dust("inst/odin/simple_SEIR_odindust.R") ##it worked!!! issue was defining dimension of the contact mixing matrix 

data(polymod, package = "socialmixr")
age.limits = seq(0, 70, 10)
contact <- socialmixr::contact_matrix(
  survey = polymod,
  countries = "United Kingdom",
  age.limits = age.limits,
  symmetric = TRUE)

transmission <- contact$matrix /
  rep(contact$demography$population, each = ncol(contact$matrix))
transmission

N_age <- length(age.limits)
n_particles <- 5L
dt <- 0.25

# ##Initial vectors
# S0[] <- user()
# E0[] <- user()
# E02[] <- user()
# I0[] <- user()
# R0[] <- user()
# 
# ##Parameters
# beta <- user()
# gamma_E <- user()
# gamma_I <- user()
# 
# 
# #Number of age classes & number of transmissibility classes
# N_age <- user()
# m[, ] <- user()
# 
# dt <- user()

model <- first_attempt$new(pars = list(dt = 0.25,
                                       S0 = contact$demography$population,
                                       E0 = c(0, 0, 0, 0, 0, 0, 0, 0),
                                       E02 = c(0, 0, 0, 0, 0, 0, 0, 0),
                                       I0 = c(0, 1000, 0, 0, 0, 0, 0, 0),
                                       R0 = c(0, 0, 0, 0, 0, 0, 0, 0),
                                       beta = 0.2 / 12.11,
                                       gamma_E = 0.05,
                                       gamma_I = 0.1,
                                       m = transmission,
                                       N_age = N_age),
                           time = 1,
                           n_particles = 5L,
                           n_threads = 1L,
                           seed = 1L)


#model$run(10)

n_times <- 10000

x <- array(NA, dim = c(model$info()$len, n_particles, n_times))

# For loop to run the model iteratively
for (t in seq_len(n_times)) {
  x[ , , t] <- model$run(t)
}
time <- x[1, 1, ]
x <- x[-1, , ]
# Plotting the trajectories
par(mfrow = c(2,4), oma=c(2,3,0,0))
for (i in 1:N_age) {
  par(mar = c(3, 4, 2, 0.5))
  cols <- c(S = "#8c8cd9", I = "#cc0044", R = "#999966")
  matplot(time, t(x[i + 3,, ]), type = "l", # Offset to access numbers in age compartment
          xlab = "", ylab = "", yaxt="none", main = paste0("Age ", contact$demography$age.group[i]),
          col = cols[["S"]], lty = 1, ylim=range(x[-1:-3,,]))
  matlines(time, t(x[i + 3 + N_age, , ]), col = cols[["I"]], lty = 1)
  matlines(time, t(x[i + 3 + 2*N_age, , ]), col = cols[["R"]], lty = 1)
  legend("right", lwd = 1, col = cols, legend = names(cols), bty = "n")
  axis(2, las =2)
}
mtext("Number of individuals", side=2,line=1, outer=T)
mtext("Time", side = 1, line = 0, outer =T)


