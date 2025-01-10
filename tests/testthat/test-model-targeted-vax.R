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
  
  expect_equal(res["cases_cumulative", ,  ], 
               res["cases_cumulative_z", ,  ] +
                 res["cases_cumulative_hh", ,  ] +
                 res["cases_cumulative_hc", ,  ] +
                 res["cases_cumulative_s", ,  ])

  expect_equal(sum(res["N_tot", , ] - sum(pars$N0)), 0)
})

test_that("check cases and deaths are counted correctly", {
  pars <- reference_pars_targeted_vax()
  nms <- reference_names()
  
  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(7, 21, by = 7)
  res <- dust2::dust_system_simulate(sys, t)
  y <- dust2::dust_unpack_state(sys, res)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))
  
  expect_equal(y$cases_inc,
               y$cases_inc_00_04 + y$cases_inc_05_14 + y$cases_inc_15_plus)
  expect_equal(y$deaths_inc,
               y$deaths_inc_00_04 + y$deaths_inc_05_14 + y$deaths_inc_15_plus)
  expect_equal(y$cases_cumulative,
               y$cases_cumulative_00_04 + y$cases_cumulative_05_14 +
                 y$cases_cumulative_15_plus)
  expect_equal(y$deaths_cumulative,
               y$deaths_cumulative_00_04 + y$deaths_cumulative_05_14 + 
                 y$deaths_cumulative_15_plus)
  
  for (nm in grep("cases_inc", rownames(res), value = TRUE)) {
    expect_equal(t(apply(res[nm, , ], 1, cumsum)),
                 res[gsub("inc", "cumulative", nm), , ])
  }

  for (nm in grep("deaths_inc", rownames(res), value = TRUE)) {
    expect_equal(t(apply(res[nm, , ], 1, cumsum)),
                 res[gsub("inc", "cumulative", nm), , ])
  }

  expect_equal(apply(y$cases_cumulative_by_age, c(2, 3), sum),
               res["cases_cumulative", , ])
})


test_that("relevant states sum correctly", {
  pars <- reference_pars_targeted_vax()
  nms <- reference_names()
  
  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  y <- dust2::dust_unpack_state(sys, res)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))
  
  ## E
  expect_equal(y$E, y$Ea + y$Eb)
  
  ## I
  expect_equal(y$I, y$Ir + y$Id)
  
  ## N
  expect_equal(y$N, y$S + y$E + y$I + y$R + y$D)
  
  ## X_tot
  expect_equal(y$S_tot, apply(y$S, c(3, 4), sum))
  expect_equal(y$E_tot, apply(y$E, c(3, 4), sum))
  expect_equal(y$I_tot, apply(y$I, c(3, 4), sum))
  expect_equal(y$R_tot, apply(y$R, c(3, 4), sum))
  expect_equal(y$D_tot, apply(y$D, c(3, 4), sum))
  expect_equal(y$N_tot, apply(y$N, c(3, 4), sum))
  
  ## total_vax
  expect_equal(y$total_vax, y$total_vax_1stdose + y$total_vax_2nddose)
  cumulative_vax_given <- 
    t(apply(y$vax_given_S + y$vax_given_Ea + y$vax_given_Eb + y$vax_given_R,
            1, cumsum))
  expect_equal(y$total_vax, cumulative_vax_given)
  cumulative_vax_1stdose_given <- 
    t(apply(y$vax_1stdose_given_S + y$vax_1stdose_given_Ea +
              y$vax_1stdose_given_Eb + y$vax_1stdose_given_R, 1, cumsum))
  expect_equal(y$total_vax_1stdose, cumulative_vax_1stdose_given)
  cumulative_vax_2nddose_given <- 
    t(apply(y$vax_2nddose_given_S + y$vax_2nddose_given_Ea +
              y$vax_2nddose_given_Eb + y$vax_2nddose_given_R, 1, cumsum))
  expect_equal(y$total_vax_2nddose, cumulative_vax_2nddose_given)
  
})


test_that("when beta_h = beta_z = beta_s = 0 there are no new infections", {
  pars <- reference_pars_targeted_vax()
  pars$beta_h <- 0
  pars$beta_s <- 0
  pars$beta_z<- rep(0,pars$n_group)
  pars$beta_hcw <- 0

  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))

  expect_true(all(res["cases_inc", , ] == 0))
  expect_equal(sum(res["N_tot", , ] - sum(pars$N0)), 0)
})


test_that("when beta_hcw is > 0 there are infections in HCW only", {
  pars <- reference_pars_targeted_vax()
  pars$beta_h <- 0
  pars$beta_s <- 0
  pars$beta_z <- rep(0, pars$n_group)
  
  idx <- get_compartment_indices()
  n_init <- sum(pars$Ea0)
  pars$Ea0[] <- 0
  pars$Eb0[idx$group$HCW, 2] <- n_init
  
  
  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))
  
  expect_true(all(res["cases_inc", , ] == res["cases_inc_HCW", , ]))
  expect_equal(sum(res["N_tot", , ] - sum(pars$N0)), 0)
  expect_equal(res["cases_cumulative_hc", , ], res["cases_cumulative", , ])
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
  pars$beta_hcw <- 0
  pars$beta_z[-pars$n_group] <- 0 # last group only for test purpose (i.e. HCW)
  n_init <- sum(pars$Ea0)
  pars$Ea0[] <- 0

  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)

  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  y <- dust2::dust_unpack_state(sys, res)

  
  # check only infections are in final group
  expect_equal(y$cases_cumulative_z, y$cases_cumulative)
  expect_equal(y$cases_cumulative_HCW, y$cases_cumulative)
  expect_equal(sum(y$I[-pars$n_group, , , ]), 0)
  expect_equal(sum(y$N_tot + n_init - sum(pars$N0)), 0)
})


test_that("when beta_h = 0 and beta_z = 0 infections only from sexual contact", {
  pars <- reference_pars_targeted_vax()
  pars$beta_h <- 0
  pars$beta_z[] <- 0
  pars$beta_s <- 0.2
  pars$beta_hcw <- 0
  pars$Ea0[] <- 0
  pars$m_sex["CSW", "PBS"] <- pars$m_sex["PBS", "CSW"] <- 0.5
  pars$m_sex["ASW", "PBS"] <- pars$m_sex["PBS", "ASW"] <- 0.5
  
  # need to seed infections as so low otherwise
  
  idx_comp <- get_compartment_indices()
  idx_kp <- unlist(idx_comp$group[c("CSW", "ASW", "PBS")])
  idx_unvax <- idx_comp$vax$unvaccinated
  
  pars$Ea0[idx_kp, idx_unvax] <- pars$Ea0[idx_kp, idx_unvax] + 10
  pars$S0[idx_kp, idx_unvax] <- pars$S0[idx_kp, idx_unvax] - pars$Ea0[idx_kp, idx_unvax]

  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  y <- dust2::dust_unpack_state(sys, res)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))


  ## should have cases in CSW and PBS
  expect_true(all(res["cases_cumulative_PBS", , max(t)] > 0))
  expect_true(all(res["cases_cumulative_SW", , max(t)] > 0))


  ## shouldn't have cases in the age groups
  expect_true(all(res["cases_cumulative_00_04",,max(t)] == 0))
  expect_equal(res["cases_inc_05_14",,] + res["cases_inc_15_plus", , ],
               res["cases_inc_SW", , ] + res["cases_inc_PBS", , ])
  
  n_age <- nrow(get_age_bins())
  expect_equal(sum(y$Ea[seq_len(n_age), , , ]), 0)

  ## make sure population size continues behaving
  expect_equal(sum(res["N_tot", , ] - sum(pars$N0) + 10), 0)

})

test_that("when daily_doses==0, no vaccinations are given", {
  pars <- reference_pars_targeted_vax()
  
  # inputs
  pars$daily_doses_adults_value[, ] <- 0
  pars$daily_doses_children_value[, ] <- 0

  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  y <- dust2::dust_unpack_state(sys, res)
  idx_comp <- get_compartment_indices()
  idx_vax <- c(idx_comp$vax$one_dose, idx_comp$vax$two_dose)
  
  ## no one should be vaccinated at any time point
  expect_true(all(y$N[, idx_vax, , ] == 0))
  
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))
  
  ## vaccines given should be 0
  expect_true(any(res[grep("vax_given",rownames(res)),,]==0))
  
  ## population check
  expect_equal(sum(res["N_tot", , ] - sum(pars$N0)), 0)
  
})



test_that("when daily_doses>0, vaccinations are given", {
  pars <- reference_pars_targeted_vax()
  
  # inputs
  # check inputs
  expect_true(sum(pars$daily_doses_children_value)>0)
  expect_true(sum(pars$daily_doses_adults_value)>0)
  
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
  expect_true(any(y$N[, idx_vax, , 21] > 0))

  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))
  
  ## vaccines given should be positive
  expect_true(any(res[grep("vax_given",rownames(res)),,]>0))

  ## population check
  expect_equal(sum(res["N_tot", , ] - sum(pars$N0)), 0)
  
  ## when time is after final vaccine time, no vaccines are given
  vaccine_end_time <- 
    max(pars$daily_doses_adults_time, pars$daily_doses_children_time)
  expect_equal(sum(y$vax_given_S[,  t > vaccine_end_time]), 0)
  expect_true(all(y$vax_given_S[,  t == vaccine_end_time] > 0))
  
  ## the vaccines given do not exceed the total set out in the strategy
  daily_doses_children <- 
    interpolate_daily_doses(pars$daily_doses_children_time,
                            pars$daily_doses_children_value)
  daily_doses_adults <- 
    interpolate_daily_doses(pars$daily_doses_adults_time,
                            pars$daily_doses_adults_value)
  expect_true(max(res["total_vax", , ]) <= 
                (sum(daily_doses_children) + sum(daily_doses_adults)))

})



test_that("vaccines are only given in the prioritised groups", {
  pars <- reference_pars_targeted_vax()

  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  y <- dust2::dust_unpack_state(sys, res)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))
  
  idx <- get_compartment_indices()

  ## identify which child groups aren't prioritised for vaccination in the first step 
  idx_novax_children <- (pars$prioritisation_strategy_children[, 1] == 0) &
    (pars$is_child > 0)
  idx_vax <- c(idx$vax$one_dose, idx$vax$two_dose)
  
  if(all(res["prioritisation_step_1st_dose_children", , ] == 1)){
    expect_true(all(y$N[idx_novax_children, idx_vax, , ] == 0))
  }
  

  ## repeat above for adults, including first and second doses 
  
  idx_novax_adults <- (pars$prioritisation_strategy_adults[, 1] == 0) &
    ((1 - pars$is_child) > 0)
  
  if(all(res["prioritisation_step_1st_dose_adults",,] == 1)){
    expect_true(all(y$N[idx_novax_adults, idx_vax, , ] == 0))
  }

})

test_that("no one moves in or out of j=1 (previous smallpox vaccine)", {
  pars <- reference_pars_targeted_vax()

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
  pars <- reference_pars_targeted_vax()
  
  # check inputs
  expect_true(sum(pars$daily_doses_children_value)>0)
  expect_true(sum(pars$daily_doses_adults_value)>0)

  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))

  # check outputs
  # 1st doses are given
  expect_true(all(res["total_vax_1stdose",,max(t)]>0))
  # 2nd doses are given
  expect_true(all(res["total_vax_2nddose",,max(t)]>0))

})

test_that("2nd doses are not given to children", {
  pars <- reference_pars_targeted_vax()
  
  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))

  # 2nd doses are given
  expect_true(all(res["total_vax_2nddose",,max(t)]>0))
  
  # children_idx
  idx_children <- (rep(ceiling(pars$is_child),pars$n_vax)) *seq(1:(pars$n_group*pars$n_vax))
  # 2nd dose compartments
  idx_children <- idx_children[which(idx_children!=0&idx_children>=3*pars$n_group)]

  idx_children_2nd_dose <- paste0(rep("N", each = length(idx_children)),
                                  idx_children)
  
  expect_true(all(res[idx_children_2nd_dose,,]==0))
  
})

test_that("prioritisation steps can progress", {
  
  pars <- reference_pars_targeted_vax()
  
  # set up parameters to push things through quickly 
  pars$daily_doses_children_value <- pars$daily_doses_children_value * 100
  pars$daily_doses_adults_value <- pars$daily_doses_adults_value * 10
  
  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))
  
  expect_true(any(res["prioritisation_step_1st_dose_children",,max(t)]!=1))
  expect_true(any(res["prioritisation_step_1st_dose_adults",,max(t)]!=1))
  expect_true(any(res["prioritisation_step_2nd_dose_adults",,max(t)]!=1))
  
})


test_that("1st/2nd dose and adult/child prioritisation steps can be different depending on the targets", {

  pars <- reference_pars_targeted_vax()
  
  # low vaccination target coverage
  pars$prioritisation_strategy_children[,1] <- pars$prioritisation_strategy_children[,1]/4
  pars$prioritisation_strategy_adults[,1] <- pars$prioritisation_strategy_adults[,1]/4
  
  # give lots of vaccines to push through quickly
  pars$daily_doses_children_value <- pars$daily_doses_children_value * 10
  pars$daily_doses_adults_value <- pars$daily_doses_adults_value * 100
  
  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))
  
  # 1st doses children vs adults
  expect_false(all(res["prioritisation_step_1st_dose_children",,]==res["prioritisation_step_1st_dose_adults",,]))
  # adults 1st vs 2 doses
  expect_false(all(res["prioritisation_step_1st_dose_adults",,]==res["prioritisation_step_2nd_dose_adults",,]))
  # adults 2nd dose vs children 1st dose
  expect_false(all(res["prioritisation_step_1st_dose_children",,]==res["prioritisation_step_2nd_dose_adults",,]))

  daily_doses_children <- 
    interpolate_daily_doses(pars$daily_doses_children_time,
                            pars$daily_doses_children_value)
  daily_doses_adults <- 
    interpolate_daily_doses(pars$daily_doses_adults_time,
                            pars$daily_doses_adults_value)
  
  # extra check here for good measure
  expect_true(all(((sum(daily_doses_children)+sum(daily_doses_adults)) - res["total_vax",,max(t)])>=0))

})


test_that("no 2nd doses are given if no 1st doses are given", {

  pars <- reference_pars_targeted_vax()
  # no 1st doses allocated, but 2nd doses still positive
  pars$daily_doses_adults_value[2, ] <- 0

  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))

  idx_adults <- ceiling(1-pars$is_child)*seq(1:(pars$n_group*pars$n_vax))
  idx_adults <- idx_adults[which(idx_adults!=0&idx_adults>2*pars$n_group)]
  
  idx_adults_vax <- paste0(rep("N", each = length(idx_adults)),
                                idx_adults)
  
  ## no 1st doses given to adults
  expect_true(all(res[idx_adults_vax,,]==0))

  ## check no 2nd doses given at all
  expect_true(all(res["total_vax_2nddose",,]==0))

})

test_that("vaccination of children still occurs if adult vaccination is turned off", {
  
  pars <- reference_pars_targeted_vax()
  # no adult doses
  pars$daily_doses_adults_value[, ] <- 0
  # children still have doses
  expect_true(any(pars$daily_doses_children_value>0))
  
  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))
  
  # expect 1st doses given
  expect_true(all(res["total_vax_1stdose",,max(t)]>0))
  
  # expect no 2nd doses given
  expect_true(all(res["total_vax_2nddose",,max(t)]==0))
  
  # check that the doses aren't in the adult groups
  idx_adults <- ceiling(1-pars$is_child)*seq(1:(pars$n_group*pars$n_vax))
  idx_adults <- idx_adults[which(idx_adults!=0&idx_adults>2*pars$n_group)]
  
  idx_adults_vax <- paste0(rep("N", each = length(idx_adults)),
                           idx_adults)
  
  ## none of the 1st doses given are given to adults
  expect_true(all(res[idx_adults_vax,,]==0))
  
})


test_that("vaccination of adults still occurs if child vaccination is turned off", {
  
  pars <- reference_pars_targeted_vax()
  # adult doses
  expect_true(any(pars$daily_doses_adults_value > 0))
  # no children doses
  pars$daily_doses_children_value[, ] <- 0
  
  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))
  
  # expect 1st doses given
  expect_true(all(res["total_vax_1stdose",,max(t)]>0))
  
  # expect 2nd doses given
  expect_true(all(res["total_vax_2nddose",,max(t)]>0))
  
  # check that the doses aren't in the child groups
  idx_children <- ceiling(pars$is_child)*seq(1:(pars$n_group*pars$n_vax))
  idx_children <- idx_children[which(idx_children!=0&idx_children>2*pars$n_group)]
  
  idx_children_vax <- paste0(rep("N", each = length(idx_children)),
                           idx_children)
  
  ## none of the 1st doses given are given to children
  expect_true(all(res[idx_children_vax,,]==0))
  
})

test_that("children vax targets are being reached as we expect before prioritisation step moves", {
  
  pars <- reference_pars_targeted_vax()
  
  # give lots of vaccines to push through quickly
  pars$daily_doses_children_value <- pars$daily_doses_children_value * 100
  
  # confirm 3 child prioritisation steps
  expect_true(ncol(pars$prioritisation_strategy_children)==3)

  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))
  
  ## for the test to work we need to have moved to at least prioritisation step 3
  expect_true(any(res["prioritisation_step_1st_dose_children",,max(t)]==3))
  
  # when do children move to the next step
  time_2nd_step_children <- which(res["prioritisation_step_1st_dose_children",1,]==2)[1]
  time_3rd_step_children <- which(res["prioritisation_step_1st_dose_children",1,]==3)[1]
  
  ## moving from step 1 to step 2
  for(g in (1:pars$n_group)){
    
    group_idx <- seq(from = g, to = pars$n_group*pars$n_vax, by = pars$n_group)
    group_idx_vax <- group_idx[which(group_idx>(2*pars$n_group))]
    
    if(all(colSums(res[paste0("N",group_idx),,(time_2nd_step_children-1)])!=0)){
      # preceding time step expect target to be met
      expect_true(all(colSums(res[paste0("N",group_idx_vax),,(time_2nd_step_children-1)])/colSums(res[paste0("N",group_idx),,(time_2nd_step_children-1)])>=pars$prioritisation_strategy_children[g,1]))
      # different targets aren't all met at once so difficult to make this generic, so instead just check that the target wasn't met in the first time step
      expect_true(all(colSums(res[paste0("N",group_idx_vax),,1])/colSums(res[paste0("N",group_idx),,1])<=pars$prioritisation_strategy_children[g,1]))
      
    }
  }
  
  ## moving from step 2 to step 3 
  for(g in 1:pars$n_group){
    
    group_idx <- seq(from = g, to = pars$n_group*pars$n_vax, by = pars$n_group)
    group_idx_vax <- group_idx[which(group_idx>(2*pars$n_group))]
    
    if(all(colSums(res[paste0("N",group_idx),,(time_3rd_step_children-1)])!=0)){
    # preceding time step expect target to be met
    expect_true(all(colSums(res[paste0("N",group_idx_vax),,(time_3rd_step_children-1)])/colSums(res[paste0("N",group_idx),,(time_3rd_step_children-1)])>=pars$prioritisation_strategy_children[g,2]))
      # if the first part of the test passes then there isn't a need to reassess the first time step here
    
  }
  
  }
  
  })


test_that("adult vax 1st dose targets are being reached as we expect before prioritisation step moves", {
  
  pars <- reference_pars_targeted_vax()
  
  # give lots of vaccines to push through quickly
  pars$daily_doses_adults_value <- pars$daily_doses_adults_value * 10
  #pars$prioritisation_strategy_adults[20,] <- 0

  # confirm 2 adult prioritisation steps
  expect_true(ncol(pars$prioritisation_strategy_adult)==2)
  
  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))
  
  ## for the test to work we need to have moved to at least prioritisation step 2
  expect_true(any(res["prioritisation_step_1st_dose_adults",,max(t)]==2))
  
  # when do adult 1st doses move to the next step
  time_2nd_step_adults <- which(res["prioritisation_step_1st_dose_adults",1,]==2)[1]

  
  ## moving from step 1 to step 2
  for(g in 1:pars$n_group){
    
    group_idx <- seq(from = g, to = pars$n_group*pars$n_vax, by = pars$n_group)
    group_idx_vax <- group_idx[which(group_idx>(2*pars$n_group))]
    
    if(all(colSums(res[paste0("N",group_idx),,(time_2nd_step_adults-1)])!=0)){
      expect_true(all(colSums(res[paste0("N",group_idx_vax),,(time_2nd_step_adults-1)])/colSums(res[paste0("N",group_idx),,(time_2nd_step_adults-1)])>=pars$prioritisation_strategy_adults[g,1]))
      # different targets aren't all met at once so difficult to make this generic, so instead just check that the target wasn't met in the first time step
      expect_true(all(colSums(res[paste0("N",group_idx_vax),,1])/colSums(res[paste0("N",group_idx),,1])<=pars$prioritisation_strategy_adults[g,1]))
      
      
    }
    
  }
  
})


test_that("adult vax 2nd dose targets are being reached as we expect before prioritisation step moves", {
  
  pars <- reference_pars_targeted_vax()
  
  # give lots of vaccines to push through quickly
  pars$daily_doses_adults_value <- pars$daily_doses_adults_value * 10
  # and lower the target to make more achievable 
  #pars$prioritisation_strategy_adults <- pars$prioritisation_strategy_adults/2
  pars$prioritisation_strategy_adults[20,] <- 0
  
  # confirm 2 adult prioritisation steps
  expect_true(ncol(pars$prioritisation_strategy_adult)==2)
  
  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))
  
  ## for the test to work we need to have moved to at least prioritisation step 2
  expect_true(any(res["prioritisation_step_2nd_dose_adults",,max(t)]==2))
  
  # when do adult 2nd doses go to the next step
  time_2nd_step_adults <- which(res["prioritisation_step_2nd_dose_adults",1,]==2)[1]
  
  
  ## moving from step 1 to step 2
  for(g in 1:pars$n_group){
    
    group_idx <- seq(from = g, to = pars$n_group*pars$n_vax, by = pars$n_group)
    group_idx_vax <- group_idx[which(group_idx>(3*pars$n_group))]
    
    if(all(colSums(res[paste0("N",group_idx),,(time_2nd_step_adults-1)])!=0)){
      expect_true(all(sum(res[paste0("N",group_idx_vax),,(time_2nd_step_adults-1)])/colSums(res[paste0("N",group_idx),,(time_2nd_step_adults-1)])>=pars$prioritisation_strategy_adults[g,1]))
      # different targets aren't all met at once so difficult to make this generic, so instead just check that the target wasn't met in the first time step
      expect_true(all(sum(res[paste0("N",group_idx_vax),,1])/colSums(res[paste0("N",group_idx),,1])<=pars$prioritisation_strategy_adults[g,1]))
    }
    
  }
  
})

test_that("Test vaccine outputs sum correctly", {
  pars <- reference_pars_targeted_vax()
  
  # set up parameters to push things through quickly 
  pars$daily_doses_children_value <- pars$daily_doses_children_value * 100
  pars$daily_doses_adults_value <- pars$daily_doses_adults_value * 10
  
  sys <- dust2::dust_system_create(model_targeted_vax(), pars, time = 1,
                                   n_particles = 3, seed = 1, dt = 1)
  dust2::dust_system_set_state_initial(sys)
  
  t <- seq(1, 21)
  res <- dust2::dust_system_simulate(sys, t)
  y <-   dust2::dust_unpack_state(sys, res)
  rownames(res) <- names(unlist(dust2::dust_unpack_index(sys)))

  idx <- get_compartment_indices()
  idx_age <- seq_len(nrow(get_age_bins()))
  group_bins <- get_group_bins()
  idx_15_plus <- which(group_bins$start >=15)
  
  ## function that does a diff but appends 0 at the beginning
  diff0 <- function(x) c(0, diff(x))
  
  ## age/group based tests
  ## second doses
  expect_equal(y$N[idx$group$ASW, idx$vax$two_dose, , ], y$dose2_cumulative_ASW)
  expect_equal(y$N[idx$group$CSW, idx$vax$two_dose, , ], y$dose2_cumulative_CSW)
  expect_equal(y$N[idx$group$HCW, idx$vax$two_dose, , ], y$dose2_cumulative_HCW)
  expect_equal(y$N[idx$group$PBS, idx$vax$two_dose, , ], y$dose2_cumulative_PBS)
  
  ## This is a bit fiddly but we calculate daily incidence from cumulative,
  ## split in half (taking the floor) and then calculate cumulative of that
  dose2_cumulative_CSW_15_plus <- 
    t(apply(floor(0.5 * apply(y$dose2_cumulative_CSW, 1, diff0)), 2, cumsum))
  expect_equal(apply(y$N[idx_15_plus, idx$vax$two_dose, , ], c(2, 3), sum) + 
                 dose2_cumulative_CSW_15_plus,
               res["dose2_cumulative_15_plus", , ])
  expect_equal(apply(y$N[c(idx$group$`5-9`, idx$group$`10-14`), idx$vax$two_dose, , ], c(2, 3), sum) +
                 y$dose2_cumulative_CSW - dose2_cumulative_CSW_15_plus,
               res["dose2_cumulative_05_14", , ])
  
   
  ## first doses - note we have to account also for people who have had two doses!
  expect_equal(apply(y$N[idx$group$ASW, c(idx$vax$one_dose, idx$vax$two_dose), , ], c(2, 3), sum), y$dose1_cumulative_ASW)
  expect_equal(apply(y$N[idx$group$CSW, c(idx$vax$one_dose, idx$vax$two_dose), , ], c(2, 3), sum), y$dose1_cumulative_CSW)
  expect_equal(apply(y$N[idx$group$HCW, c(idx$vax$one_dose, idx$vax$two_dose), , ], c(2, 3), sum), y$dose1_cumulative_HCW)
  expect_equal(apply(y$N[idx$group$PBS, c(idx$vax$one_dose, idx$vax$two_dose), , ], c(2, 3), sum), y$dose1_cumulative_PBS)
  
  dose1_cumulative_CSW_15_plus <- 
    t(apply(floor(0.5 * apply(y$dose1_cumulative_CSW, 1, diff0)), 2, cumsum)) 
  expect_equal(apply(y$N[idx_15_plus, c(idx$vax$one_dose,idx$vax$two_dose), , ], c(3, 4), sum) +
                 dose1_cumulative_CSW_15_plus,
               res["dose1_cumulative_15_plus", , ])
  expect_equal(apply(y$N[c(idx$group$`5-9`, idx$group$`10-14`), c(idx$vax$one_dose, idx$vax$two_dose), , ], c(3, 4), sum) +
                 y$dose1_cumulative_CSW - dose1_cumulative_CSW_15_plus,
               res["dose1_cumulative_05_14", , ])

  expect_equal(y$N[idx$group$`0-4`, idx$vax$one_dose, , ],
               res["dose1_cumulative_00_04", , ])
  
  
  # check age outputs sum to total doses given
  expect_equal(res["total_vax_1stdose", , ],
               res["dose1_cumulative_00_04", , ] +
                 res["dose1_cumulative_05_14", , ] +
                 res["dose1_cumulative_15_plus", , ])
  
  expect_equal(res["total_vax_2nddose", , ],
               res["dose2_cumulative_00_04", , ] +
                 res["dose2_cumulative_05_14", , ] +
                 res["dose2_cumulative_15_plus", , ])
  
  # check key pop outputs sum to total doses given
  expect_equal(res["dose1_cumulative_SW", , ],
               res["dose1_cumulative_CSW", , ] +
                 res["dose1_cumulative_ASW", , ])
  expect_equal(res["dose2_cumulative_SW", , ],
               res["dose2_cumulative_CSW", , ] +
                 res["dose2_cumulative_ASW", , ])
  
  expect_equal(res["total_vax_1stdose", , ],
               res["dose1_cumulative_SW", , ] +
                 res["dose1_cumulative_PBS", , ] +
                 res["dose1_cumulative_HCW", , ] +
                 apply(y$N[idx_age, c(idx$vax$one_dose, idx$vax$two_dose), , ],
                       c(3, 4), sum))
  expect_equal(res["total_vax_2nddose", , ],
               res["dose2_cumulative_SW", , ] +
                 res["dose2_cumulative_PBS", , ] +
                 res["dose2_cumulative_HCW", , ] +
                 apply(y$N[idx_age, idx$vax$two_dose, , ], c(2, 3), sum))
  
  # check all incident and cumulative dose outputs match
  t_eow <- seq(7, max(t), 7)
  nms <- grep("dose[0-9]_inc", rownames(res), value = TRUE)
  
  for (nm in nms) {
    expect_equal(t(apply(res[nm, , t_eow], 1, cumsum)),
                 res[gsub("inc", "cumulative", nm), , t_eow])
  }

})

test_that("Test compiled compare components", {
  pars <- reference_pars_targeted_vax()
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
                  deaths_15_plus = 25,
                  cases_total = 150,
                  cases_HCW = 5,
                  cases_SW = 10,
                  cfr_00_04 = 0.01,
                  cfr_05_14 = 0.03,
                  cfr_15_plus = 0.05)
  
  parts <- list(
    c("cases"),
    c("cases_00_04"),
    c("cases_05_14"),
    c("cases_15_plus"),
    c("deaths"),
    c("deaths_00_04"),
    c("deaths_05_14"),
    c("deaths_15_plus"),
    c("cases_HCW", "cases_total"),
    c("cases_SW", "cases_total"),
    c("cfr_00_04"),
    c("cfr_05_14"),
    c("cfr_15_plus"))

  
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
