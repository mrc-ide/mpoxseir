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
    "cases_inc", "deaths_inc",
    "cases_inc_00_04", "cases_inc_05_14", "cases_inc_15_plus", 
    "cases_inc_PBS", "cases_inc_SW",
    "deaths_inc_00_04", "deaths_inc_05_14", "deaths_inc_15_plus", 
    "deaths_inc_PBS", "deaths_inc_SW")
  save <- c(run, "S_tot", "E_tot", "I_tot", "R_tot", "D_tot", "N_tot")
  index <- unlist(info$index)
  list(run = index[run], state = index[save])
}
