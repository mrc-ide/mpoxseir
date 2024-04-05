## can we run compartmental-mpox like we would squire?

library(odin)
odin::can_compile()

library(devtools)

Sys.getenv("GITHUB_PAT")
Sys.unsetenv("GITHUB_PAT")

devtools::install_github("mrc-ide/compartmental-mpox") ### because private?

## instead have to run manually to check is working

# from squire help page: https://github.com/mrc-ide/squire
#r <- run_explicit_SEEIR_model(country = "Afghanistan")

## want to run simple_SEIR.R (inst/odin/)

source("R/run.R")
source("R/parameters.R")
source("R/check.R")

#r <- run_simple_SEEIR_model(population=100000,contact_matrix_set = TRUE)
## look at get_mixing_matrix

## try to run the bigger model as a check 
source("R/population.R")
source("R/assertions.R")
library(tidyverse)
source("R/contact_matrices.R")
source("R/beta.R")
source("R/odin.R")


r2 <- run_explicit_SEEIR_model(country="Afghanistan")
plot(r2)
### THIS RAN YAY


## now back to simple_SEIR
## last error with contact_matrix_set
## from ?squire::run_simple_SEEIR_model

pop <- get_population("Afghanistan", simple_SEIR = TRUE)
# m1 <- run_simple_SEEIR_model(population = pop$n, dt = 1,
#                              R0 = 2,
#                              contact_matrix_set=contact_matrices[[1]])

contact_matrices <- get_mixing_matrix(country="Afghanistan")

r <- run_simple_SEEIR_model(population=pop$n,
                            dt=1,
                            R0=2,
                            contact_matrix_set = contact_matrices)
plot(r)


## THIS ALSO RAN!