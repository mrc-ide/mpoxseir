#!/usr/bin/env Rscript
if (packageVersion("odin2") < "0.3.34") {
  stop("Please upgrade odin2 to at least 0.3.34")
}
if (packageVersion("dust2") < "0.3.24") {
  stop("Please upgrade dust2 to at least 0.3.24")
}

odin2::odin_package(here::here(), check_bounds = "disabled")
devtools::load_all(here::here())
