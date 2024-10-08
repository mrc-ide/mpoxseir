#!/usr/bin/env Rscript
if (packageVersion("odin.dust") < "0.2.29") {
  stop("Please upgrade odin.dust to at least 0.2.29")
}

withr::with_options(
  list(odin.no_check_naked_index = TRUE),
  withr::with_collate(
    "C",
    odin.dust::odin_dust_package(
      here::here()
      , options = odin.dust::odin_dust_options(debug_enable = FALSE)
    )
  )
)
devtools::load_all()
