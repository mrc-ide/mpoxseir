test_that("assign seeds works", {
  tmp <- assign_seeds(10, rep(1 / 3, 3))
  expect_equal(sum(tmp), 10)
  expect_equal(assign_seeds(6, c(3, 2, 1)), c(3, 2, 1))
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
