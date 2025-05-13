#!/usr/bin/env Rscript
if (packageVersion("odin2") < "0.3.24") {
  stop("Please upgrade odin2 to at least 0.3.24")
}
if (packageVersion("dust2") < "0.3.22") {
  stop("Please upgrade dust2 to at least 0.3.22")
}

odin2::odin_package(here::here())
devtools::load_all(here::here())
