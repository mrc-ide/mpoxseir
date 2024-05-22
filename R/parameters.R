
parameters_demographic <- function(country = "Democratic Republic of Congo") {

  data <- squire::get_population(country)
  population <- c(data$n[1:15], sum(data$n[16:17])) # combine 75+
  nms <- c(levels(data$age_group)[1:15], "75+")
  m <- squire::get_mixing_matrix(country)
  dimnames(m) <- list(nms, nms)

  list(
    N_age = length(population),
    population = setNames(population, nms),
    m = m
  )
}
