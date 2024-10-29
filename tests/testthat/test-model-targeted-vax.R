test_that("run is equal to reference", {
  pars <- reference_pars_targeted_vax()
  nms <- reference_names()

  m <- model_targeted_vax$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))

  expect_true(any(res["cases_inc", , ] > 0))
  expect_true(any(res["deaths_inc", , ] > 0))

  expect_equal(sum(res["N_tot", , ] - sum(pars$N0)), 0)
})

test_that("check cases and deaths are counted correctly", {
  pars <- reference_pars_targeted_vax()
  nms <- reference_names()
  
  m <- model_targeted_vax$new(pars, 1, 3, seed = 1)
  t <- seq(7, 21, by = 7)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))
  y <- m$transform_variables(res)
  
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
  expect_equal(apply(y$cases_inc, 2, sum),
               y$cases_cumulative[, , length(t)])
  expect_equal(apply(y$cases_inc_00_04, 2, sum),
               y$cases_cumulative_00_04[, , length(t)])
  expect_equal(apply(y$cases_inc_05_14, 2, sum),
               y$cases_cumulative_05_14[, , length(t)])
  expect_equal(apply(y$cases_inc_15_plus, 2, sum),
               y$cases_cumulative_15_plus[, , length(t)])
  expect_equal(apply(y$cases_inc_PBS, 2, sum),
               y$cases_cumulative_PBS[, , length(t)])
  expect_equal(apply(y$cases_inc_SW, 2, sum),
               y$cases_cumulative_SW[, , length(t)])
  expect_equal(apply(y$cases_inc_HCW, 2, sum),
               y$cases_cumulative_HCW[, , length(t)])
  
  expect_equal(apply(y$deaths_inc, 2, sum),
               y$deaths_cumulative[, , length(t)])
  expect_equal(apply(y$deaths_inc_00_04, 2, sum),
               y$deaths_cumulative_00_04[, , length(t)])
  expect_equal(apply(y$deaths_inc_05_14, 2, sum),
               y$deaths_cumulative_05_14[, , length(t)])
  expect_equal(apply(y$deaths_inc_15_plus, 2, sum),
               y$deaths_cumulative_15_plus[, , length(t)])
  expect_equal(apply(y$deaths_inc_PBS, 2, sum),
               y$deaths_cumulative_PBS[, , length(t)])
  expect_equal(apply(y$deaths_inc_SW, 2, sum),
               y$deaths_cumulative_SW[, , length(t)])
  expect_equal(apply(y$deaths_inc_HCW, 2, sum),
               y$deaths_cumulative_HCW[, , length(t)])

})


test_that("when beta_h = beta_z = beta_s = 0 there are no new infections", {
  pars <- reference_pars_targeted_vax()
  pars$beta_h <- 0
  pars$beta_s <- 0
  pars$beta_z<- rep(0,pars$n_group)

  m <- model_targeted_vax$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))

  expect_true(all(res["cases_inc", , ] == 0))
  expect_equal(sum(res["N_tot", , ] - sum(pars$N0)), 0)
})

test_that("when CFR = 0 nobody dies", {
  pars <- reference_pars_targeted_vax()
  pars$CFR[] <- 0

  m <- model_targeted_vax$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))

  expect_true(any(res["cases_inc", , ] > 0))
  expect_true(all(res["deaths_inc", , ] == 0))
  expect_true(any(res["R_tot", , ] > 0))
  expect_true(all(res["D_tot", , ] == 0))
  expect_equal(sum(res["N_tot", , ] - sum(pars$N0)), 0)
})

test_that("when CFR = 1 everybody dies", {
  pars <- reference_pars_targeted_vax()
  pars$CFR[] <- 1

  m <- model_targeted_vax$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))

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
  pars$beta_z[-pars$n_group] <- 0 # last group only for test purpose (i.e. HCW)
  n_init <- sum(pars$Ea0)
  pars$Ea0[] <- 0

  m <- model_targeted_vax$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  idx <- m$info()$index
  rownames(res) <- names(unlist(idx))
  y <- m$transform_variables(res)
  names(y) <- names(idx)
  
  # check only infections are in SW
  expect_equal(y$cases_cumulative_HCW, y$cases_cumulative)
  expect_equal(sum(y$I[-pars$n_group, , , ]), 0)
  expect_equal(sum(y$N_tot + n_init - sum(pars$N0)), 0)
})


test_that("when beta_h = 0 and beta_z = 0 infections only from sexual contact", {
  pars <- reference_pars_targeted_vax()
  pars$beta_h <- 0
  pars$beta_z[] <- 0
  pars$beta_s <- 0.2
  pars$Ea0[] <- 0
  pars$m_sex["CSW", "PBS"] <- pars$m_sex["PBS", "CSW"] <- 0.5
  pars$m_sex["ASW", "PBS"] <- pars$m_sex["PBS", "ASW"] <- 0.5
  
  # need to seed infections as so low otherwise
  
  idx_comp <- get_compartment_indices()
  idx_kp <- unlist(idx_comp$group[c("CSW", "ASW", "PBS")])
  idx_unvax <- idx_comp$vax$unvaccinated
  
  pars$Ea0[idx_kp, idx_unvax] <- pars$Ea0[idx_kp, idx_unvax] + 10
  pars$S0[idx_kp, idx_unvax] <- pars$S0[idx_kp, idx_unvax] - pars$Ea0[idx_kp, idx_unvax]

  m <- model_targeted_vax$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21) 
  res <- m$simulate(t)
  idx <- m$info()$index
  rownames(res) <- names(unlist(idx))


  ## should have cases in CSW and PBS
  expect_true(all(res["cases_cumulative_PBS", , max(t)] > 0))
  expect_true(all(res["cases_cumulative_SW", , max(t)] > 0))


  ## shouldn't have cases in the age groups
  expect_true(all(res["cases_cumulative_00_04",,max(t)] == 0))
  expect_equal(res["cases_inc_05_14",,] + res["cases_inc_15_plus", , ],
               res["cases_inc_SW", , ] + res["cases_inc_PBS", , ])
  y <- m$transform_variables(res)
  n_age <- nrow(get_age_bins())
  expect_equal(sum(y$Ea[seq_len(n_age), , , ]), 0)

  ## make sure population size continues behaving
  expect_equal(sum(res["N_tot", , ] - sum(pars$N0) + 10), 0)

})

test_that("when daily_doses==0, no vaccinations are given", {
  pars <- reference_pars_targeted_vax()
  
  # inputs
  pars$daily_doses_adults[,] <- 0
  pars$daily_doses_children[,] <- 0

  m <- model_targeted_vax$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  idx <- m$info()$index
  rownames(res) <- names(unlist(idx))
  y <- m$transform_variables(res)
  names(y) <- names(idx)
  idx_comp <- get_compartment_indices()
  idx_vax <- c(idx_comp$vax$one_dose, idx_comp$vax$two_dose)
  
  ## no one should be vaccinated at any time point
  expect_true(all(y$N[, idx_vax, , ] == 0))
  
  ## vaccines given should be 0
  expect_true(any(res[grep("vax_given",rownames(res)),,]==0))
  
  ## population check
  expect_equal(sum(res["N_tot", , ] - sum(pars$N0)), 0)
  
})



test_that("when daily_doses>0, vaccinations are given", {
  pars <- reference_pars_targeted_vax()
  
  # inputs
  # check inputs
  expect_true(sum(pars$daily_doses_children)>0)
  expect_true(sum(pars$daily_doses_adults)>0)
  
  m <- model_targeted_vax$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  idx <- m$info()$index
  rownames(res) <- names(unlist(idx))
  y <- m$transform_variables(res)
  names(y) <- names(idx)
  idx_comp <- get_compartment_indices()
  idx_vax <- c(idx_comp$vax$one_dose, idx_comp$vax$two_dose)

  ## in first time point there should be no one vaccinated
  expect_true(all(y$N[, idx_vax, , 1] == 0))
  # at last time point we should have lots vaccinated
  expect_true(any(y$N[, idx_vax, , 21] > 0))

  ## vaccines given should be positive
  expect_true(any(res[grep("vax_given",rownames(res)),,]>0))

  ## population check
  expect_equal(sum(res["N_tot", , ] - sum(pars$N0)), 0)
  
  ## when time is after vaccination_campaign_length, no vaccines are given
  expect_equal(sum(y$vax_given_S[, ,  t > pars$vaccination_campaign_length]), 0)
  expect_true(all(y$vax_given_S[, ,  t == pars$vaccination_campaign_length] > 0))
  
  ## the vaccines given do not exceed the total set out in the strategy
  expect_true(max(res["total_vax", , ]) <= (sum(pars$daily_doses_children)+sum(pars$daily_doses_adults)))

})



test_that("vaccines are only given in the prioritised groups", {
  pars <- reference_pars_targeted_vax()

  m <- model_targeted_vax$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))

  y <-  m$transform_variables(res)
  idx <- get_compartment_indices()

  pars$prioritisation_strategy_adults
  ## identify which child groups aren't prioritised for vaccination in the first step 
  idx_novax_children <- (pars$prioritisation_strategy_children[, 1] == 0) &
    (pars$children_ind_raw > 0)
  idx_vax <- c(idx$vax$one_dose, idx$vax$two_dose)
  
  if(all(res["prioritisation_step_1st_dose_children", , ] == 1)){
    expect_true(all(y$N[idx_novax_children, idx_vax, , ] == 0))
  }
  

  ## repeat above for adults, including first and second doses 
  
  idx_novax_adults <- (pars$prioritisation_strategy_adults[, 1] == 0) &
    (pars$adults_ind_raw > 0)
  
  if(all(res["prioritisation_step_1st_dose_adults",,] == 1)){
    expect_true(all(y$N[idx_novax_adults, idx_vax, , ] == 0))
  }

})

test_that("no one moves in or out of j=1 (previous smallpox vaccine)", {
  pars <- reference_pars_targeted_vax()

  m <- model_targeted_vax$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))

  # the total number of people in j=1 shouldn't change
  expect_true(all(res[paste0("N",seq(1:pars$n_group)),,]==pars$S0[,1]))

  ## population check
  expect_equal(sum(res["N_tot", , ] - sum(pars$N0)), 0)

})


test_that("1st and 2nd doses are given", {
  pars <- reference_pars_targeted_vax()
  
  # check inputs
  expect_true(sum(pars$daily_doses_children)>0)
  expect_true(sum(pars$daily_doses_adults)>0)

  m <- model_targeted_vax$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))

  # check outputs
  # 1st doses are given
  expect_true(all(res["total_vax_1stdose",,max(t)]>0))
  # 2nd doses are given
  expect_true(all(res["total_vax_2nddose",,max(t)]>0))

})

test_that("2nd doses are not given to children", {
  pars <- reference_pars_targeted_vax()
  
  m <- model_targeted_vax$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))

  # 2nd doses are given
  expect_true(all(res["total_vax_2nddose",,max(t)]>0))
  
  # children_idx
  idx_children <- (rep(ceiling(pars$children_ind_raw),pars$n_vax)) *seq(1:(pars$n_group*pars$n_vax))
  # 2nd dose compartments
  idx_children <- idx_children[which(idx_children!=0&idx_children>=3*pars$n_group)]

  idx_children_2nd_dose <- paste0(rep("N", each = length(idx_children)),
                                  idx_children)
  
  expect_true(all(res[idx_children_2nd_dose,,]==0))
  
})

test_that("prioritisation steps can progress", {
  
  pars <- reference_pars_targeted_vax()
  
  # set up parameters to push things through quickly 
  pars$daily_doses_children <- pars$daily_doses_children * 100
  pars$daily_doses_adults <- pars$daily_doses_adults * 10
  
  m <- model_targeted_vax$new(pars, 1, 3, seed = 1)
  t <- seq(1, pars$vaccination_campaign_length_adults)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))
  
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
  pars$daily_doses_children <- pars$daily_doses_children*100
  pars$daily_doses_adults <- pars$daily_doses_adults*100
  
  m <- model_targeted_vax$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))
  
  # 1st doses children vs adults
  expect_false(all(res["prioritisation_step_1st_dose_children",,]==res["prioritisation_step_1st_dose_adults",,]))
  # adults 1st vs 2 doses
  expect_false(all(res["prioritisation_step_1st_dose_adults",,]==res["prioritisation_step_2nd_dose_adults",,]))
  # adults 2nd dose vs children 1st dose
  expect_false(all(res["prioritisation_step_1st_dose_children",,]==res["prioritisation_step_2nd_dose_adults",,]))

  # extra check here for good measure
  expect_true(all(((sum(pars$daily_doses_children)+sum(pars$daily_doses_adults)) - res["total_vax",,max(t)])>=0))

})


test_that("no 2nd doses are given if no 1st doses are given", {

  pars <- reference_pars_targeted_vax()
  # no 1st doses allocated, but 2nd doses still positive
  pars$daily_doses_adults[,2] <- 0

  m <- model_targeted_vax$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))

  idx_adults <- ceiling(pars$adults_ind_raw)*seq(1:(pars$n_group*pars$n_vax))
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
  pars$daily_doses_adults[,] <- 0
  # children still have doses
  expect_true(any(pars$daily_doses_children>0))
  
  m <- model_targeted_vax$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))
  
  # expect 1st doses given
  expect_true(all(res["total_vax_1stdose",,max(t)]>0))
  
  # expect no 2nd doses given
  expect_true(all(res["total_vax_2nddose",,max(t)]==0))
  
  # check that the doses aren't in the adult groups
  idx_adults <- ceiling(pars$adults_ind_raw)*seq(1:(pars$n_group*pars$n_vax))
  idx_adults <- idx_adults[which(idx_adults!=0&idx_adults>2*pars$n_group)]
  
  idx_adults_vax <- paste0(rep("N", each = length(idx_adults)),
                           idx_adults)
  
  ## none of the 1st doses given are given to adults
  expect_true(all(res[idx_adults_vax,,]==0))
  
})


test_that("vaccination of adults still occurs if child vaccination is turned off", {
  
  pars <- reference_pars_targeted_vax()
  # adult doses
  expect_true(any(pars$daily_doses_adults>0))
  # no children doses
  pars$daily_doses_children[,] <- 0
  
  m <- model_targeted_vax$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))
  
  # expect 1st doses given
  expect_true(all(res["total_vax_1stdose",,max(t)]>0))
  
  # expect 2nd doses given
  expect_true(all(res["total_vax_2nddose",,max(t)]>0))
  
  # check that the doses aren't in the child groups
  idx_children <- ceiling(pars$children_ind_raw)*seq(1:(pars$n_group*pars$n_vax))
  idx_children <- idx_children[which(idx_children!=0&idx_children>2*pars$n_group)]
  
  idx_children_vax <- paste0(rep("N", each = length(idx_children)),
                           idx_children)
  
  ## none of the 1st doses given are given to children
  expect_true(all(res[idx_children_vax,,]==0))
  
})

test_that("children vax targets are being reached as we expect before prioritisation step moves", {
  
  pars <- reference_pars_targeted_vax()
  
  # give lots of vaccines to push through quickly
  pars$daily_doses_children <- pars$daily_doses_children*100
  
  # confirm 3 child prioritisation steps
  expect_true(ncol(pars$prioritisation_strategy_children)==3)

  m <- model_targeted_vax$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))
  
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
  pars$daily_doses_adults <- pars$daily_doses_adults *10
  #pars$prioritisation_strategy_adults[20,] <- 0

  # confirm 2 adult prioritisation steps
  expect_true(ncol(pars$prioritisation_strategy_adult)==2)
  
  m <- model_targeted_vax$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))
  
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
  pars$daily_doses_adults <- pars$daily_doses_adults *10
  # and lower the target to make more achievable 
  #pars$prioritisation_strategy_adults <- pars$prioritisation_strategy_adults/2
  pars$prioritisation_strategy_adults[20,] <- 0
  
  # confirm 2 adult prioritisation steps
  expect_true(ncol(pars$prioritisation_strategy_adult)==2)
  
  m <- model_targeted_vax$new(pars, 1, 3, seed = 1)
  t <- seq(1, 21)
  res <- m$simulate(t)
  rownames(res) <- names(unlist(m$info()$index))
  
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


test_that("Test compiled compare components", {
  pars <- reference_pars_targeted_vax()
  pars$exp_noise <- Inf
  nms <- reference_names()
  
  m <- model_targeted_vax$new(pars, 1, 3, seed = 1L)
  
  time <- 350
  y <- m$run(time)
  
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
    m$set_data(dust::dust_data(d_test, "time"))
    m$compare_data()
  }
  
  ll_parts <- lapply(parts, compare_part)
  
  ll_all <- compare_part(do.call(cbind, parts))
  
  ## check that using each datastream individually sums to using them all
  expect_equal(ll_all, rowSums(do.call(cbind, ll_parts)))
  
  ## check that all datastreams are supplying non-zero likelihood contributions
  expect_true(all(sapply(ll_parts, function(x) any(x != 0))))
  
})

