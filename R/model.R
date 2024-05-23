##' @name model
##' @title Stochastic compartmental transmission model of Mpox virus
##' @description The basic model; we may add more or adapt this one
##'   over time.
##' @export model
NULL

##' Index of model outputs needed for fitting and validation. This function
##' conforms to the mcstate interface.
##' @title Index of model variables
##' @param info The result of running the `$info()` method on an
##'   initialised [model]
##' @return A list with element `run`, indicating the locations of the
##'   compartments used to compare
##' @export
model_index <- function(info) {
  run <- c("cases", "deaths")
  state <- c(run, "S", "E", "I", "R")
  index <- unlist(info$index)
  list(run = index[run], state = index[state])
}
