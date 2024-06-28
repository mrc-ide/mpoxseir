test_that("run is equal to reference", {
  pars <- reference_pars()
  nms <- reference_names()

  m <- model$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))

  expect_true(any(res["weekly_cases", , ] > 0))
  expect_true(any(res["weekly_deaths", , ] > 0))

  expect_equal(sum(res["N_tot", , ] - sum(pars$N)), 0)
})


test_that("when beta_h = beta_z = 0 there are no new infections", {
  pars <- reference_pars()
  pars$beta_h <- 0
  pars$beta_z<- rep(0,pars$n_group)

  m <- model$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))

  expect_true(all(res["weekly_cases", , ] == 0))
  expect_equal(sum(res["N_tot", , ] - sum(pars$N)), 0)
})

test_that("when CFR = 0 nobody dies", {
  pars <- reference_pars()
  pars$CFR[] <- 0

  m <- model$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))

  expect_true(any(res["weekly_cases", , ] > 0))
  expect_true(all(res["weekly_deaths", , ] == 0))
  expect_true(any(res["R_tot", , ] > 0))
  expect_true(all(res["D_tot", , ] == 0))
  expect_equal(sum(res["N_tot", , ] - sum(pars$N)), 0)
})

test_that("when CFR = 1 everybody dies", {
  pars <- reference_pars()
  pars$CFR[] <- 1

  m <- model$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))

  expect_true(any(res["weekly_cases", , ] > 0))
  expect_true(any(res["weekly_deaths", , ] > 0))
  expect_true(all(res["R_tot", , ] == 0))
  expect_true(any(res["D_tot", , ] > 0))
  expect_equal(sum(res["N_tot", , ] - sum(pars$N)), 0)
})


test_that("when beta_h = 0 there are only zoonotic infections", {
  pars <- reference_pars()
  pars$beta_h <- 0
  pars$beta_z <- c(rep(0, pars$n_group - 1), 0.4 / 12.11) # last group only for test purpose

  m <- model$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))

  expect_true(all(res[paste0("Ea", seq_len(pars$n_group - 1)), , ] == 0))
  expect_true(all(res[paste0("Eb", seq_len(pars$n_group - 1)), , ] == 0))
  expect_true(any(res["Ea18", , ] > 0))
  expect_true(any(res["Eb18", , ] > 0))
})

test_that("weekly outputs work as expected", {
  pars <- reference_pars()
  pars$beta_z[] <- 1e-4
  pars$beta_h <- 0
  m <- model$new(pars, 1, 3, seed = 1)
  t <- seq(1, 22)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))

  week_ends <- which(t %% 7 == 0)

  # weekly cases reset every 7 days
  expect_true(all(res["weekly_cases", , week_ends] >
                    res["weekly_cases", , week_ends - 1]))
  # cumulative cases match (N - S) = sum(weekly_cases)
  expect_equal(res["N_tot", , max(t)] - res["S_tot", , max(t)],
              rowSums(res["weekly_cases", , week_ends]))
  # weekly cases by age sum to total
  expect_equal(colSums(res[grep("weekly_cases", rownames(res)), , max(t)]),
               res["weekly_cases", , max(t)] * 2)

  # weekly deaths reset every 7 days
  expect_true(all(res["weekly_deaths", , week_ends] >
                    res["weekly_deaths", , week_ends - 1]))
  # cumulative deaths match (N - S) = sum(weekly_deaths)
  expect_equal(res["D_tot", , max(t)],
               rowSums(res["weekly_deaths", , week_ends]))
  # weekly deaths by age sum to total
  expect_equal(colSums(res[grep("weekly_deaths", rownames(res)), , max(t)]),
               res["weekly_deaths", , max(t)] * 2)
})
