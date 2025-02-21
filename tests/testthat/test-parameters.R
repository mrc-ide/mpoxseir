test_that("assign seeds works", {
  tmp <- assign_seeds(10, rep(1 / 3, 3))
  expect_equal(sum(tmp), 10)
  expect_equal(assign_seeds(6, c(3, 2, 1)), c(3, 2, 1))
})

test_that("Clade Ib seeding in (Adult)SW actually happens", {
  initial_infections <- 150
  
  pars_sudkivu <- parameters_fixed("sudkivu", initial_infections = initial_infections)
  expect_equal(sum(pars_sudkivu$Ea0) , initial_infections)
  expect_equal(sum(pars_sudkivu$Ea0[18,2]) , initial_infections) # check all are in unvaccinated ASW
  
  pars_burundi <- parameters_fixed("burundi", initial_infections = initial_infections)
  expect_equal(sum(pars_burundi$Ea0) , initial_infections)
  expect_equal(sum(pars_burundi$Ea0[18,2]) , initial_infections) # check all are in unvaccinated ASW
  
  
  pars_bujumbura <- parameters_fixed("bujumbura", initial_infections = initial_infections)
  expect_equal(sum(pars_bujumbura$Ea0) , initial_infections)
  expect_equal(sum(pars_bujumbura$Ea0[18,2]) , initial_infections) # check all are in unvaccinated ASW
})

test_that("Mixing matrices are correct in Equateur", {
  pars <- parameters_demographic(region = "equateur")
  
  expect_equal(sum(pars$N0), sum(pars$N_age))
  
  M <- pars$total_contacts_nonsexual
  M_age <- pars$total_contacts_age
  
  m <- M / pars$N0
  m[is.na(m)] <- 0

  expect_true(isSymmetric(M))
  expect_equal(sum(M_age - t(M_age)), 0)
  
  expect_equal(sum(upper.tri(M, diag = TRUE) * M), 
               sum(upper.tri(M_age, diag = TRUE) * M_age))
  expect_equal(m, pars$m_gen_pop)
})

test_that("Mixing matrices are correct in Sud Kivu", {
  pars <- parameters_demographic(region = "sudkivu")
  
  expect_equal(sum(pars$N0), sum(pars$N_age))
  
  M <- pars$total_contacts_nonsexual
  M_age <- pars$total_contacts_age
  
  m <- M / pars$N0
  m[is.na(m)] <- 0
  
  expect_true(isSymmetric(M))
  expect_equal(sum(M_age - t(M_age)), 0)
  
  expect_equal(sum(upper.tri(M, diag = TRUE) * M), 
               sum(upper.tri(M_age, diag = TRUE) * M_age))
  expect_equal(m, pars$m_gen_pop)
})

test_that("can load other mixing matrices", {
  pars <- parameters_demographic(region = "equateur",
                                 mixing_matrix = "synthetic_home")
  old_pars <- parameters_demographic(region = "equateur")
  pars$total_contacts_age
  old_pars$total_contacts_age
  expect_equal(pars$N_age, old_pars$N_age)
  expect_false(any(pars$total_contacts_age == old_pars$total_contacts_age))

  pars_all <- parameters_demographic(region = "equateur",
                                      mixing_matrix = "synthetic_all")
  expect_true(all(pars$total_contacts_age < pars_all$total_contacts_age))

  expect_error(parameters_demographic(region = "equateur",
                                      mixing_matrix = "hello"))
})


test_that("unknown region throws error", {
  
  expect_error(parameters_fixed("burundi_not_bujumbura",
                                initial_infections = 10))
  
})


test_that("total sex workers in Bujumbura as we expect", {
  
  expect_equal(
    sum(parameters_fixed("bujumbura_mairie",
                         initial_infections = 10)$N0[c("CSW","ASW")]),
    3852
  )

})


test_that("Overrides works", {
  
  expect_equal(
    parameters_fixed("equateur", initial_infections = 10,
                     overrides = list(gamma_E = 1 / 8))$gamma_E,
    1 / 8
  )
  
  expect_error(
    parameters_fixed("equateur", initial_infections = 10,
                     overrides = 1 / 8),
    'overrides must be a list'
  )
  
  expect_error(
    parameters_fixed("equateur", initial_infections = 10,
                     overrides = list(gamma_X = 1 / 8)),
    'unknown parameter gamma_X'
  )
  
})

test_that("use_ve_D works", {
  
  p1 <- parameters_fixed("equateur", initial_infections = 10, use_ve_D = FALSE)
  p2 <- parameters_fixed("equateur", initial_infections = 10, use_ve_D = TRUE)
  
  idx <- get_compartment_indices()
  
  expect_equal(p1$CFR[, idx$vax$unvaccinated], p2$CFR[, idx$vax$unvaccinated])
  expect_vector_lt(p2$CFR[, -idx$vax$unvaccinated],
                   p1$CFR[, -idx$vax$unvaccinated])
  
})

test_that("create_age_bins works", {
  
  start <- seq(0, 75, 5)
  
  age_bins <- create_age_bins(start = seq(0, 75, 5))
  
  expect_equal(age_bins$start, start)
  expect_equal(age_bins$end, c(start[-1L] - 1, 100))
  
  expected_labels <- c(paste0(start[-length(start)], "-", start[-1L] - 1),
                       paste0(start[length(start)], "+"))
  expect_equal(age_bins$label, expected_labels)
  
  expect_error(create_age_bins(start = seq(0, 75, 5), max = 70),
               "max_age is too small")
  
})

test_that("Can not seed too high", {
  
  expect_error(parameters_fixed("equateur", initial_infections = 1e10),
               "population size and seeding infections is incompatible")
  
})
