test_that("run is equal to reference", {
  pars <- reference_pars_targeted_vax()
  nms <- reference_names()

  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))

  expect_true(any(res["cases_inc", , ] > 0))
  expect_true(any(res["deaths_inc", , ] > 0))

  expect_equal(sum(res["N_tot", , ] - sum(pars$N0)), 0)
})


test_that("when beta_h = beta_z = beta_s = 0 there are no new infections", {
  pars <- reference_pars_targeted_vax()
  pars$beta_h <- 0
  pars$beta_s <- 0
  pars$beta_z<- rep(0,pars$n_group)

  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))

  expect_true(all(res["cases_inc", , ] == 0))
  expect_equal(sum(res["N_tot", , ] - sum(pars$N0)), 0)
})

test_that("when CFR = 0 nobody dies", {
  pars <- reference_pars_targeted_vax()
  pars$CFR[] <- 0

  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))

  expect_true(any(res["cases_inc", , ] > 0))
  expect_true(all(res["deaths_inc", , ] == 0))
  expect_true(any(res["R_tot", , ] > 0))
  expect_true(all(res["D_tot", , ] == 0))
  expect_equal(sum(res["N_tot", , ] - sum(pars$N0)), 0)
})

test_that("when CFR = 1 everybody dies", {
  pars <- reference_pars_targeted_vax()
  pars$CFR[] <- 1

  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))

  expect_true(any(res["cases_inc", , ] > 0))
  expect_true(any(res["deaths_inc", , ] > 0))
  expect_true(all(res["R_tot", , ] == 0))
  expect_true(any(res["D_tot", , ] > 0))
  expect_equal(sum(res["N_tot", , ] - sum(pars$N0)), 0)
})


test_that("when beta_h = 0 and beta_s = 0 there are only zoonotic infections", {
  pars <- reference_pars_targeted_vax()
  pars$beta_h <- 0
  pars$beta_s <- 0
  pars$beta_z[-pars$n_group] <- 0 # last group only for test purpose (i.e. SW)
  n_init <- sum(pars$Ea0)
  pars$Ea0[] <- 0

  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  y <- dust2::dust_unpack_state(sys, res)
  
  # check only infections are in SW
  expect_equal(y$cases_cumulative_SW, y$cases_cumulative)
  expect_equal(sum(y$I[-pars$n_group, , , ]), 0)
  expect_equal(sum(y$N_tot + n_init - sum(pars$N0)), 0)
})


test_that("when beta_h = 0 and beta_z = 0 infections only from sexual contact", {
  pars <- reference_pars_targeted_vax()
  pars$beta_h <- 0
  pars$beta_z[] <- 0
  pars$beta_s <- 0.2
  pars$m_sex["SW", "PBS"] <- pars$m_sex["PBS", "SW"] <- 0.5
  
  # need to seed infections as so low otherwise
  
  idx_comp <- get_compartment_indices()
  idx_kp <- unlist(idx_comp$group[c("SW", "PBS")])
  idx_unvax <- idx_comp$vax$unvaccinated
  
  pars$Ir0[idx_kp, idx_unvax] <- pars$Ir0[idx_kp, idx_unvax] + 10
  pars$Id0[idx_kp, idx_unvax] <- pars$Id0[idx_kp, idx_unvax] + 10
  pars$S0[idx_kp, idx_unvax] <- pars$S0[idx_kp, idx_unvax] -
    pars$Id0[idx_kp, idx_unvax] - pars$Ir0[idx_kp, idx_unvax]

  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))


  ## should have cases in CSW and PBS
  expect_true(all(res["cases_cumulative_PBS", , max(t)]>0))
  expect_true(all(res["cases_cumulative_SW",,max(t)]>0))

  ## shouldn't have cases in the age groups
  expect_true(all(res["cases_cumulative_00_04",,max(t)]==0))
  expect_true(all(res["cases_cumulative_05_14",,max(t)]==0))
  expect_true(all(res["cases_cumulative_15_plus",,max(t)]==0))

  ## make sure population size continues behaving
  expect_equal(sum(res["N_tot", , ] - sum(pars$N0)), 0)

})




test_that("when n_vaccination>0, vaccinations are given", {
  pars <- reference_pars_targeted_vax(uptake = 0.8)
  pars$prioritisation_strategy[] <- 1
  

  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  y <- dust2::dust_unpack_state(sys, res)
  
  idx_comp <- get_compartment_indices()
  idx_vax <- c(idx_comp$vax$one_dose, idx_comp$vax$two_dose)


  ## in first time point there should be no one vaccinated
  expect_true(all(y$N[, idx_vax, , 1] == 0))
  # at last time point we should have lots vaccinated
  expect_true(all(y$N[, idx_vax, , 21] > 0))

  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))
  ## vaccines given should be positive
  expect_true(any(res[grep("vax_given",rownames(res)),,]>0))

  ## population check
  expect_equal(sum(res["N_tot", , ] - sum(pars$N0)), 0)
  
  ## when time is after vaccination_campaign_length, no vaccines are given
  expect_equal(sum(y$vax_given_S[,  t >= max(pars$daily_doses_time)]), 0)
  expect_true(all(y$vax_given_S[,  t >= max(pars$daily_doses_time)] > 0))
  
  ## the vaccines given do not exceed the total set out in the strategy
  expect_true(max(res["total_vax", , ]) <= sum(pars$daily_doses))

})



test_that("if prioritisation_step==1, vaccines are only given in groups 1 - 3", {
  pars <- reference_pars_targeted_vax()

  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))

  if(all(res["prioritisation_step_1st_dose",,]==1)){
    expect_true(all(res[paste0("S",(2*pars$n_group+4):(3*pars$n_group)),,]==0))
    expect_true(all(res[paste0("Ea",
                               (2*pars$n_group+4):(3*pars$n_group)),,]==0))
    expect_true(all(res[paste0("Eb",
                               (2*pars$n_group+4):(3*pars$n_group)),,]==0))
    expect_true(all(res[paste0("R",
                               (2*pars$n_group+4):(3*pars$n_group)),,]==0))
  }

  if(all(res["prioritisation_step_2nd_dose",,]==1)){
    expect_true(all(res[paste0("S",(3*pars$n_group+4):(4*pars$n_group)),,]==0))
    expect_true(all(res[paste0("Ea",
                               (3*pars$n_group+4):(4*pars$n_group)),,]==0))
    expect_true(all(res[paste0("Eb",
                               (3*pars$n_group+4):(4*pars$n_group)),,]==0))
    expect_true(all(res[paste0("R",
                               (3*pars$n_group+4):(4*pars$n_group)),,]==0))
  }

})

test_that("if vaccine_uptake = 0.5, half the expected vaccines are given out", {
  pars <- reference_pars_targeted_vax(uptake = 0.5)

  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))

  # doses given out should be roughly half of those allocated (roughly comes from rounding due to lack of multinomial)
  expect_true(all(res["total_vax", , max(t)]<=0.5*sum(pars$daily_doses)))


})


test_that("no one moves in or out of j=1 (previous smallpox vaccine)", {
  pars <- reference_pars_targeted_vax(uptake = 1)

  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))

  # the total number of people in j=1 shouldn't change
  expect_true(all(res[paste0("N",seq(1:pars$n_group)),,]==pars$S0[,1]))

  ## population check
  expect_equal(sum(res["N_tot", , ] - sum(pars$N0)), 0)

})


test_that("1st and 2nd doses are given", {
  pars <- reference_pars_targeted_vax(uptake = 1)

  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))

  # 1st doses are given
  expect_true(all(res["total_vax_1stdose",,max(t)]>0))
  # 2nd doses are given
  expect_true(all(res["total_vax_2nddose",,max(t)]>0))

})


test_that("1st and 2nd prioritisation steps can be different depending on the targets", {

  pars <- reference_pars_targeted_vax(uptake = 1)
  # low vaccination target coverage
  pars$vaccination_coverage_target_1st_dose_prop <- 0.1
  pars$vaccination_campaign_length <- 15
  pars$daily_doses <- matrix(0,ncol=pars$n_vax,
                        nrow=pars$vaccination_campaign_length)
  # give loads of vaccines to push through prioiritsation quickly
  pars$daily_doses[1:10,2] <- pars$daily_doses[11:15,3] <- 1000000
  pars$daily_doses[pars$vaccination_campaign_length,] <- 0

  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))

  expect_false(all(res["prioritisation_step_1st_dose",,]==res["prioritisation_step_2nd_dose",,]))

  # extra check here for good measure
  expect_true(all((sum(pars$daily_doses) - res["total_vax",,max(t)])>0))

})

test_that("no 2nd doses are given if no 1st doses are given", {

  pars <- reference_pars_targeted_vax(uptake = 1)
  # no 1st doses allocated
  pars$daily_doses_value[2, ] <- 0

  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))

  ## check no 1st doses given
  expect_true(all(res["total_vax_1stdose",,]==0))

  ## check no 2nd doses given
  expect_true(all(res["total_vax_2nddose",,]==0))

})


test_that("Test compiled compare components", {
  pars <- reference_pars_targeted_vax(uptake = 1)
  pars$exp_noise <- Inf
  nms <- reference_names()
  
  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  time <- 350
  y <- dust2::dust_system_run_to_time(sys, time)
  
  d <- data.frame(time = 350,
                  cases = 150,
                  cases_00_04 = 30,
                  cases_05_14 = 40,
                  cases_15_plus = 80,
                  deaths = 50,
                  deaths_00_04 = 10,
                  deaths_05_14 = 15,
                  deaths_15_plus = 25)
  
  parts <- list(
    c("cases"),
    c("cases_00_04"),
    c("cases_05_14"),
    c("cases_15_plus"),
    c("deaths"),
    c("deaths_00_04"),
    c("deaths_05_14"),
    c("deaths_15_plus"))
  
  
  compare_part <- function(nms) {
    d_test <- d
    d_test[, setdiff(names(d), c(nms, "time"))] <- NA_real_
    dust2::dust_system_compare_data(sys, d_test)
  }
  
  ll_parts <- lapply(parts, compare_part)
  
  ll_all <- compare_part(do.call(cbind, parts))
  
  ## check that using each datastream individually sums to using them all
  expect_equal(ll_all, rowSums(do.call(cbind, ll_parts)))
    
  ## check that all datastreams are supplying non-zero likelihood contributions
  expect_true(all(sapply(ll_parts, function(x) any(x != 0))))
    
})
