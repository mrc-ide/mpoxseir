## can we run compartmental-mpox like we would squire?

"squire" %in% loadedNamespaces() ### want to be FALSE

library(odin)
odin::can_compile()

# library(devtools)
# 
# Sys.getenv("GITHUB_PAT")
# Sys.unsetenv("GITHUB_PAT")
# 
# devtools::install_github("mrc-ide/compartmental-mpox") ### didn't work, maybe because private?

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

"squire" %in% loadedNamespaces() ## TRUE (something must load the package when running the code?)
## would have a dependency if so
## happens during run_explicit_SEEIR_model


## now back to simple_SEIR
## last error with contact_matrix_set
## from ?squire::run_simple_SEEIR_model

pop <- get_population("Afghanistan", simple_SEIR = TRUE)
## this loads squire 
# m1 <- run_simple_SEEIR_model(population = pop$n, dt = 1,
#                              R0 = 2,
#                              contact_matrix_set=contact_matrices[[1]])
unloadNamespace("squire")

contact_matrices <- get_mixing_matrix(country="Afghanistan")
## this loads squire 
unloadNamespace("squire")

r <- run_simple_SEEIR_model(population=pop$n,
                            dt=1,
                            R0=2,
                            contact_matrix_set = contact_matrices)
## but this doesn't load squire
plot(r)

## THIS ALSO RAN!




