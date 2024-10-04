
##' A function that returns the demographic parameters for use in the model
##' 
##' @title Get demographic parameters
##' 
##' @return A list containing all the demographic parameters
##'   
##' @export
##' @importFrom stats dmultinom
##' @importFrom stats setNames
##' @importFrom squire get_population
##' @importFrom squire get_mixing_matrix
##' 
##' @export
parameters_demographic <- function() {
  
  age_bins <- get_age_bins()

  ## Set up population denominators
  country <- "Democratic Republic of Congo"
  data <- squire::get_population(country)
  N_age <- c(data$n[1:15], sum(data$n[16:17])) # combine 75+

  ## Add in Sex workers (SW) and People who Buy Sex (PBS)
  
  idx_12_49 <- age_bins$start>=10 & age_bins$end <= 49
  
  ## CSW
  ## splitting this for child (12 - 17) and adult SWs (18 - 49)
  idx_12_17 <- age_bins$start >=10 & age_bins$end <= 19
  idx_18_49  <- age_bins$start >= 15 & age_bins$end <= 49
  
  N_12_17 <- N_age * idx_12_17
  ## adjust for different age bands
  N_12_17[which(age_bins=="10-14")] <- 0.6*N_12_17[which(age_bins=="10-14")]
  N_12_17[which(age_bins=="15-19")] <- 0.6*N_12_17[which(age_bins=="15-19")]
  
  N_18_49 <- N_age * idx_18_49
  ## adjust for different age bands
  N_18_49[which(age_bins=="15-19")] <- 0.4*N_18_49[which(age_bins=="15-19")]

  p_SW <- 0.007 * 0.5 # 0.7% women (50%) 15-49 Laga et al
  N_SW_children <- round(p_SW * N_12_17)
  N_SW_adults <- round(p_SW * N_18_49)

  ## PBS
  
  p_PBS <- 0.11 * 0.5 # 11% men (50%) 15-49 DHS https://www.statcompiler.com/en/
  idx_15_49  <- age_bins$start >= 15 & age_bins$end <= 49
  N_15_49 <- N_age * idx_15_49
  
  N_PBS <- round(p_PBS * N_15_49)

  N <- c(N_age - N_SW_adults - N_PBS - N_SW_children, 
         sum(N_SW_adults), sum(N_PBS), sum(N_SW_children))
  
  # check allocated properly
  sum(N) - sum(N_age)
  # all fine

  # Set up mixing matrices: 1) general population and 2) sexual contact
  # Squire gives unbalanced per-capita daily rates, we need to
  # 1. scale these by the population size to get total daily contacts between
  #    age group i and age group j
  # 2. balance the matrix so total contacts i->j == j->i
  # 3. add in non-age groups (SW + PBS) and attribute contacts
  # 4. convert back to rates

  m_age <- squire::get_mixing_matrix(country)
  M_raw <- t(t(m_age) * N_age) # total daily contacts in population

  ## balance the matrix so M[i,j] == M[j,i]
  M_age <- (M_raw + t(M_raw)) / 2

  # Need to split total contacts i->j by gen pop / SW / PBS
  # assume homogenous mixing in day-to-day contacts (will add sex in fitting)

  # decompose total contact matrix by number in pair that could be KP (0/1/2)
  
  M0 <- M_age * outer(!idx_12_49, !idx_12_49) # neither could be KP
  ## border case of 10 - 14 (reduce contacts to 40% of the total contacts in this band e.g. 10 and 11 year olds)
  M0[which(age_bins=="10-14"),!idx_12_49] <- 0.4 * M_age[which(age_bins=="10-14"),!idx_12_49] 
  M0[!idx_12_49,which(age_bins=="10-14")] <- 0.4 * M_age[!idx_12_49,which(age_bins=="10-14")]
  M0[which(age_bins=="10-14"),which(age_bins=="10-14")] <- 0.4 * M_age[which(age_bins=="10-14"),which(age_bins=="10-14")]
  
  M1 <- M_age * outer(idx_12_49, idx_12_49, FUN = "xor") # only 1 could be KP
  ## border case of 10 - 14 (reduce contacts to 60% of this band, e.g. 12, 13, 14)
  M1[which(age_bins=="10-14"),] <- 0.6 * M1[which(age_bins=="10-14"),] 
  M1[,which(age_bins=="10-14")] <- 0.6 * M1[,which(age_bins=="10-14")]
  
  M2 <- M_age * outer(idx_12_49, idx_12_49) # both could be KP
  ## border case of 10 - 14
  M2[which(age_bins=="10-14"),] <- 0.6 * M2[which(age_bins=="10-14"),]
  M2[,which(age_bins=="10-14")] <- 0.6 * M2[,which(age_bins=="10-14")]

  p_KP <- p_SW + p_PBS
  p <- c(1 - p_KP, p_SW, p_PBS)

  # Split contact matrix by all 6 combinations of contact between groups
  M_gen_pop <- M0 + (1 - p_KP) * M1 + (1 - p_KP)^2 * M2 # gen pop x gen pop
  M_gen_SW  <- M1 * p_SW  + M2 * dmultinom(c(1, 1, 0), 2, p) # gen pop x SW 
  M_gen_PBS <- M1 * p_PBS + M2 * dmultinom(c(1, 0, 1), 2, p) # gen pop x PBS
  M_SW_SW <- M2 * p_SW ^ 2 # SW x SW
  M_PBS_PBS <- M2 * p_PBS ^ 2 # PBS x PBS
  M_SW_PBS <- M2 * dmultinom(c(0, 1, 1), 2, p) # SW x PBS

  ## check matrices are decomposed properly
  ## LILITH right now they are not sorry
  sum(M_gen_pop + M_gen_PBS + M_gen_SW + M_SW_SW + M_PBS_PBS + M_SW_PBS - M_age)

  marginalise <- function(M) (rowSums(M) + diag(M)) / 2

  n_age <- length(N_age)
  idx_age <- seq_len(n_age)
  n_group <- n_age + 3
  nms_group <- c(age_bins$label, "SW_adults", "PBS", "SW_children")


  ### General population
  # Construct new contact matrix including groups
  M_all <- matrix(0, n_group, n_group, dimnames = list(nms_group, nms_group))
  M_all[idx_age, idx_age] <- M_gen_pop
  M_all["SW_adults", idx_age] <- M_all[idx_age, "SW_adults"] <- marginalise(M_gen_SW)
  M_all["PBS", idx_age] <- M_all[idx_age, "PBS"] <- marginalise(M_gen_PBS)

  M_all["SW_adults", "SW_adults"] <- sum(marginalise(M_SW_SW))
  M_all["SW_children", "SW_children"] <- sum(marginalise(M_SW_SW))
  M_all["PBS", "PBS"] <- sum(marginalise(M_PBS_PBS))
  M_all["SW_adults", "PBS"] <- M_all["PBS", "SW_adults"] <- sum(marginalise(M_SW_PBS))
  M_all["SW_children", "PBS"] <- M_all["PBS", "SW_children"] <- sum(marginalise(M_SW_PBS))

  # check that total contacts are the same as original
  sum(marginalise(M_all)) - sum(marginalise(M_age)) ## check total number of contacts

  # Convert to per-capita rates by dividing by population
  # Resulting matrix is Asymmetric c_ij != c_ji
  # BUT total number of contacts i->j and j->i is balanced
  m_all <- M_all / N


  ## set up sexual contact matrix for parameterisation in transform function
  m_sex <- matrix(0, n_group, n_group, dimnames = list(nms_group, nms_group))

  # proportion of susceptibles estimated to have smallpox vaccine
  sus_prop <- c(rep(1,8),0.54,0.29,0.29,0.23,0.21,0.21,0.21,0.21,1,1,1)

  # province populations
  province_pop = list("equateur" = 1712000,
                      "sudkivu" = 6565000)


  list(
    n_group = n_group,
    N0 = setNames(N, nms_group),
    m_gen_pop = m_all,
    m_sex = m_sex,
    total_contacts_gen_pop = M_all,
    #total_contacts_sex = M_sex,
    n_vax = 4,
    sus_prop = sus_prop,
    province_pop = province_pop
  )
}


##' A function that gets the age bins used in the model.
##' 
##' @title Get age bins for use in the model
##' 
##' @return A data frame containing the labels for the age bins, and their
##'   start and end values
##'   
##' @export
get_age_bins <- function() {
  ## We always use these age bands, so rather than detect them, we will
  ## check that things conform to them.
  end <- c(seq(4, 75, by = 5), 100)
  start <- c(0, end[-length(end)] + 1L)
  label <- paste(start, end, sep = "-")
  label[length(label)] <- paste0(max(start), "+")

  data.frame(label = label, start = start, end = end)
}

# Function to allocate N individuals into m groups based on weights
assign_seeds <- function(N, w) {

    w_norm <- w / sum(w)
    raw_alloc <- N * w_norm
    int_alloc <- floor(raw_alloc)

    remainder <- N - sum(int_alloc)
    fractional_parts <- raw_alloc - int_alloc

    # Distribute the remainder to the groups with the largest fractional parts
    if (remainder > 0) {
      extra_alloc <- order(fractional_parts * w_norm, decreasing = TRUE)[1:remainder]
      int_alloc[extra_alloc] <- int_alloc[extra_alloc] + 1
    }

    int_alloc
  }

##' A function that gets the fixed parameters for use in the model
##' 
##' @title Get fixed parameters for use in the model
##' 
##' @param region The region for the parameters, must be either `"equateur"` or
##'   `"sudkivu"`
##'   
##' @param initial_infections The initial number of infections
##' 
##' @param overrides A list, containing any parameters for which you want to
##'   override the default values
##' 
##' @return A list of the fixed parameters
##'   
##' @export
##' 
#' @export
parameters_fixed <- function(region, initial_infections, overrides = list()) {

  ## Checking region
  if (!(region %in% c("equateur", "sudkivu"))) {
    stop("region must be equateur or sudkivu")
  }

  ## Initialising variable that other parameters depend on
  demographic_params <- parameters_demographic()
  age_bins <- get_age_bins()

  n_group <- demographic_params$n_group
  n_vax <- demographic_params$n_vax
  N <- demographic_params$province_pop[[region]]
  N0 <- round(N * demographic_params$N0 / sum(demographic_params$N0)) # total number in each age-group
  N_prioritisation_steps <- 1

  ## Seed infections in the unvaccinated group in a region-specific manner
  X0 <- matrix(0, nrow = n_group, ncol = n_vax)
  Ea0 <- matrix(0, nrow = n_group, ncol = n_vax)
  RR_z <- c(0.977, 1, 0.444, rep(0.078, n_group - 3)) # Jezek 1988 zoonotic + Jezek 1987
  if (region == "sudkivu") { # seeding in sex workers in Sud Kivu

    ## Extract sex-worker index and put initial infections in this group (unvaccinated strata)
    sw_index <- which(colnames(demographic_params$m_gen_pop) == "SW")
    Ea0[sw_index, 2] <- initial_infections

  } else if (region == "equateur") { # seeding in general pop in proportion to zoonotic risk in equateur

    ## Extract gen-pop index and put initial infections in this group (unvaccinated strata) in proportion to zoonotic risk
    gen_pop_index <- seq_len(nrow(age_bins))
    seeding_infections <- assign_seeds(initial_infections, w = RR_z[gen_pop_index])

    Ea0[gen_pop_index, 2] <- seeding_infections

  } else {
    stop("something is wrong with the name of the region - change to sudkivu or equateur")
  }

  # CFR from Whittles 2024, 5-year bands to 40

  CFR <- rep(0, n_group)
  names(CFR) <- names(demographic_params$N0)
  CFR[which(age_bins$end < 40)] <- c(0.102, 0.054, 0.035, 0.026, 0.02, 0.016, 0.013, 0.012)
  CFR[which(age_bins$start >= 40)] <- 0.01
  CFR["SW"] <- CFR["20-24"]
  CFR["PBS"] <- CFR["35-39"]

  vaccination_campaign_length <- 1

  ## update with assignment into smallpox vaccination compartment (j=1)
  ## if we have a small population then seeding with 5 infections per compt and strata could be a problem (e.g. with N=10000 too few SW)
  if(all(rowSums(Ea0)<N0)){
      S0 <- matrix(0,nrow=demographic_params$n_group,
                   ncol=demographic_params$n_vax)
      S0[,2] <- round((N0-rowSums(Ea0))*demographic_params$sus_prop)
      #rbinom(size=N0,prob=demographic_params$sus_prop,n=length(N0)) - Ea0[,2]
      S0[,1] <- N0-rowSums(Ea0)-S0[,2]
      ##rowSums(S0)+rowSums(Ea0)==N0
  }else{
    stop("combination of population size and seeding infections is incompatible, please review.")
  }

  params_list = list(
    region = region,
    n_group = n_group,
    n_vax = n_vax,
    N_prioritisation_steps = N_prioritisation_steps,
    S0 = S0,
    Ea0 = Ea0,
    Eb0 = X0,
    Ir0 = X0,
    Id0 = X0,
    R0 = X0,
    D0 = X0,
    N0 = N0,
    R0_hh = 0.67, # Jezek 1988 SAR paper - will be fitted
    R0_sw_st = 1.3, # Will be fitted
    beta_z_max = 0.01, # Will be fitted
    RR_z = RR_z,
    gamma_E = 1 / 7,  #  1/7 based on Besombes et al. on 29 clade I patients
    gamma_Ir = 1 / 18, # Jezek 1988 "clinical features of 282.."
    gamma_Id = 1 / 10, # Jezek 1988
    CFR = matrix(CFR, nrow = n_group, ncol = n_vax, byrow = FALSE),
    m_sex = demographic_params$m_sex,
    m_gen_pop = demographic_params$m_gen_pop,
    prioritisation_strategy = matrix(1, nrow = n_group, ncol = N_prioritisation_steps),
    # vaccination_coverage_target = matrix(0.01, nrow = n_group, ncol = N_prioritisation_steps),
    vaccine_uptake = rep(0.8, n_group),
    ve_T = rep(0, n_vax),
    ve_I = rep(0, n_vax),
    vaccination_coverage_target_1st_dose_prop = 0.8,
    vaccination_coverage_target_2nd_dose_prop = 0.5,
    vaccination_campaign_length = vaccination_campaign_length,
    daily_doses = matrix(1, nrow = vaccination_campaign_length, ncol = n_vax))

  # Ensure overridden parameters are passed as a list
  if (!is.list(overrides)) {
    stop('overrides must be a list')
  }

  # Override parameter values in the overrides input
  for (name in names(overrides)) {
    if (!(name %in% names(params_list))) {
      stop(paste('unknown parameter', name, sep=' '))
    }
    params_list[[name]] <- overrides[[name]]
  }

  return(params_list)
}

