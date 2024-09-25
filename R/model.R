##' @name model_targeted_vax
##' @title Stochastic compartmental transmission model of Mpox virus
##' @description The basic model; we may add more or adapt this one
##'   over time.
##' @export model_targeted_vax
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
  run <- c(
    "cases",  "deaths",
    "cases_0_5", "cases_05_15", "cases_15_plus", "cases_PBS", "cases_SW",
    "deaths_0_5", "deaths_05_15", "deaths_15_plus", "deaths_PBS", "deaths_SW")
  save <- c(run, "S_tot", "E_tot", "I_tot", "R_tot", "D_tot", "N_tot")
  index <- unlist(info$index)
  list(run = index[run], state = index[save])
}
