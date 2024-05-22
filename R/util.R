#' @noRd
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

#' @noRd
is_ptr_null <- function(pointer){
  a <- attributes(pointer)
  attributes(pointer) <- NULL
  out <- identical(pointer, methods::new("externalptr"))
  attributes(pointer) <- a
  return(out)
}

## Index locations of outputs in odin model
#' @noRd
odin_index <- function(model) {
  len <- length(model$.__enclos_env__$private$ynames)
  model$transform_variables(seq_len(len))
}


## Indices for cumulative cases total
#' @noRd
cases_total_index <- function(model) {

  index <- odin_index(model)
  indices <- c("IMild", "ICase1", "ICase2", "IOxGetLive1", "IOxGetLive2",
               "IOxGetDie1", "IOxGetDie2", "IOxNotGetLive1", "IOxNotGetLive2",
               "IOxNotGetDie1", "IOxNotGetDie2", "IMVGetLive1", "IMVGetLive2",
               "IMVGetDie1", "IMVGetDie2", "IMVNotGetLive1", "IMVNotGetLive2",
               "IMVNotGetDie1", "IMVNotGetDie2", "IRec1", "IRec2", "R1", "R2","D")
  return(unlist(index[indices]))
}

