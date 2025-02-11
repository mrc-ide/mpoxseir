#!/usr/bin/env Rscript
if (packageVersion("odin2") < "0.3.17") {
  stop("Please upgrade odin2 to at least 0.3.17")
}

odin2::odin_package(here::here())
devtools::load_all(here::here())

