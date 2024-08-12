test_that("run is equal to reference", {
  pars <- reference_pars()
  nms <- reference_names()

  m <- model$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))

  expect_true(any(res["cases", , ] > 0))
  expect_true(any(res["deaths", , ] > 0))

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

  expect_true(all(res["cases", , ] == 0))
  expect_equal(sum(res["N_tot", , ] - sum(pars$N)), 0)
})

test_that("when CFR = 0 nobody dies", {
  pars <- reference_pars()
  pars$CFR[] <- 0

  m <- model$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))

  expect_true(any(res["cases", , ] > 0))
  expect_true(all(res["deaths", , ] == 0))
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

  expect_true(any(res["cases", , ] > 0))
  expect_true(any(res["deaths", , ] > 0))
  expect_true(all(res["R_tot", , ] == 0))
  expect_true(any(res["D_tot", , ] > 0))
  expect_equal(sum(res["N_tot", , ] - sum(pars$N)), 0)
})


test_that("when beta_h = 0 there are only zoonotic infections", {
  pars <- reference_pars()
  pars$beta_h <- 0
  pars$beta_z<- c(rep(0,pars$n_group-1),0.4 / 12.11) # last group only for test purpose
  
  m <- model$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))
  
  ## isolate the compartments which can have infections (e.g. the last group)
  comps_of_interest_Ea <- paste0("Ea",
                                seq(from=pars$n_group,
                                    length.out=pars$n_vax,by=pars$n_group))
  comps_of_interest_Eb <- paste0("Eb",seq(from=pars$n_group,
                                         length.out=pars$n_vax,by=pars$n_group))
  
  expect_true(all(res[grep("^Ea",rownames(res)),,][!(rownames(res)[grep("^Ea",rownames(res))] %in% comps_of_interest_Ea)]==0))
  expect_true(all(res[grep("^Eb",rownames(res)),,][!(rownames(res)[grep("^Eb",rownames(res))] %in% comps_of_interest_Eb)]==0))
  
  expect_equal(sum(res["N_tot", , ] - sum(pars$N)), 0)
  
})

test_that("when n_vaccination>0, vaccinations are given", {
  pars <- reference_pars()

  m <- model$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))
  
  comps_of_interest_S <- paste0("S",
                                seq(from=pars$n_group+1,
                                    to=pars$n_group*pars$n_vax,by=1))
  comps_of_interest_Ea <- paste0("Ea",
                                 seq(from=pars$n_group+1,
                                     to=pars$n_group*pars$n_vax,by=1))
  comps_of_interest_Eb <- paste0("Eb",
                                 seq(from=pars$n_group+1,
                                     to=pars$n_group*pars$n_vax,by=1))
  comps_of_interest_R <- paste0("R",
                                seq(from=pars$n_group+1,
                                    to=pars$n_group*pars$n_vax,by=1))

  
  ## in first time point there should be no one vaccinated 
  expect_true(all(res[comps_of_interest_S, ,1] == 0))
  expect_true(all(res[comps_of_interest_Ea, ,1] == 0))
  expect_true(all(res[comps_of_interest_Eb, ,1] == 0))
  expect_true(all(res[comps_of_interest_R, ,1] == 0))
  
  ## for subsequent time points we should have people in these classes
  expect_true(any(res[comps_of_interest_S,,]>0))
  expect_true(any(res[comps_of_interest_Ea,,]>0))
  expect_true(any(res[comps_of_interest_Eb,,]>0))
  expect_true(any(res[comps_of_interest_R,,]>0))
  
  ## vaccines given should be positive 
  expect_true(any(res[grep("vax_given",rownames(res)),,]>0))

  ## population check 
  expect_equal(sum(res["N_tot", , ] - sum(pars$N)), 0)
  
})


test_that("when time is after vaccination_campaign_length, no vaccines are given", {
  pars <- reference_pars()
  
  m <- model$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))
  
  expect_true(all(res[grep("vax_given",rownames(res)),,seq(pars$vaccination_campaign_length+1,max(t),1)]==0))
  
  
})


