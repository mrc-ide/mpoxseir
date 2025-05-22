
##' A function that returns the demographic parameters for use in the model
##' 
##' @title Get demographic parameters
##' 
##' @param region The region for the parameters, must be either `"equateur"`, 
##'   `"sudkivu"`, `"burundi"` or `"bujumbura"`
##' @param mixing_matrix The mixing matrix must be either `"Zimbabwe"`,
##'  `"synthetic_home"`, or `"synthetic_all"`
##' @param p_SW The proportion of SW-age groups that are sex workers. Note that
##'   e.g. a value of 0.01 means 1% of all SW-age groups are sex workers and not
##'   just 1% of women in those groups. Default is NULL, in which case we use
##'   the default value for the region given in the package
##' @param p_HCW The proportion of HCW-age groups that are healthcare workers. Default is NULL, in which the default values for DRC and Burundi are used.
##' @return A list containing all the demographic parameters
##'   
##' @export
##' @importFrom stats dmultinom
##' @importFrom stats setNames
##' @importFrom squire get_population
##' @importFrom squire get_mixing_matrix
##' 
##' @export
parameters_demographic <- function(region, mixing_matrix = "Zimbabwe",
                                   p_SW = NULL,
                                   p_HCW = NULL) {
  age_bins <- get_age_bins()
  squire_age_bins <- create_age_bins(start = seq(0, 75, 5))
  group_bins <- get_group_bins()
  row.names(group_bins) <- group_bins$label
  
  ## Set up population denominators
  if(region %in% c("equateur","sudkivu")){
  country <- "Democratic Republic of Congo"
  } else if(region %in% c("burundi","bujumbura","bujumbura_mairie")){
    country <- "Burundi"
  }

  data <- squire::get_population(country)
  
  ## combine 75+
  N_age_squire <- c(data$n[1:15], sum(data$n[16:17]))
  names(N_age_squire) <- squire_age_bins$label

  ## Calc new N_age
  N_age <- N_age_squire
  names(N_age) <- age_bins$label
  p_10to14_in_05to11 <- 0.4
  p_10to14_in_12to14 <- 1 - p_10to14_in_05to11
  N_age["5-11"] <- N_age_squire["5-9"] + p_10to14_in_05to11 * N_age_squire["10-14"]
  N_age["12-14"] <- p_10to14_in_12to14 * N_age_squire["10-14"]
  
  ## Add in (child/adult) sex workers, (C/ASW), people who buy sex (PBS),
  ## and healthcare workers (HCW)
  
  ## CSW: age 12-17
  
  w_CSW <- proportion_in_age_bins(group_bins["CSW", "start"],
                                  group_bins["CSW", "end"])
  # 60% 10-14 + 60% 15-19
  N_CSW <- N_age * w_CSW
  
  ## ASW: age 18-49
  w_ASW  <- proportion_in_age_bins(group_bins["ASW", "start"],
                                   group_bins["ASW", "end"])
  N_ASW <- N_age * w_ASW
  
  if (region == "equateur") {
    p_SW_default <- 0.007 * 0.5 
    # 0.7% women (50%) 15-49 Laga et al - assume this holds down to age 12
  } else if (region == "sudkivu"){
    p_SW_default <- 0.03 * 0.5 # WHO press release
  } else if (region == "burundi"){
    p_SW_default <- 0.028 * 0.5 # Laga et al
  } else if (region %in% c("bujumbura","bujumbura_mairie")){
    p_poss_SW  <-  sum(N_ASW + N_CSW)/ sum(N_age)
    p_SW_default <- 3852 / (792503 * p_poss_SW)
    ## key pops report, taken as median of Bujumbura Mairie in Figure 7 
  }
  
  p_SW <- p_SW %||% p_SW_default
  
  N_CSW <- round(p_SW * N_CSW)
  N_ASW <- round(p_SW * N_ASW)
  
  ## HCW / PBS: age 20-49
  w_PBS <- proportion_in_age_bins(group_bins["PBS", "start"],
                                  group_bins["PBS", "end"])
  N_PBS <- N_age * w_PBS
  
  ## PBS
  # not changing this based on age as we only heard about young SWs not young PBS
  p_PBS <- 0.11 * 0.5 # 11% men (50%) 15-49 DHS https://www.statcompiler.com/en/
  # we assume all greater than 20 to avoid vaccine issues
  ## for Burundi same source saying 0.6% which seems far too low
  N_PBS <- round(p_PBS * N_PBS)
  
  ## HCWs: age 20-69 (from https://apps.who.int/nhwaportal/)
  w_HCW <- proportion_in_age_bins(group_bins["HCW", "start"],
                                  group_bins["HCW", "end"])
  N_HCW <- N_age * w_HCW
  
  if(region %in% c("equateur","sudkivu")){
    p_HCW_default <- 136606 / sum(N_age)
    } else if(region %in% c("burundi","bujumbura","bujumbura_mairie")){
    p_HCW_default <- 11911 / sum(N_age)
    }
  
  p_HCW <- p_HCW %||% p_HCW_default

   # possibly want to reduce this further to account for fact that not every HCW will have contact with mpox patients? 
  N_HCW <- round(p_HCW * N_HCW)
  
  N <- c(N_age - N_ASW - N_PBS - N_CSW - N_HCW,
         CSW = sum(N_CSW),
         ASW = sum(N_ASW),
         PBS = sum(N_PBS),
         HCW = sum(N_HCW))

  # Set up mixing matrices: 1) general population and 2) sexual contact
  # Squire gives unbalanced per-capita daily rates, we need to
  # 1. scale these by the population size to get total daily contacts between
  #    age group i and age group j
  # 2. balance the matrix so total contacts i->j == j->i
  # 3. add in non-age groups (CSW, ASW, PBS, HCW) and attribute non-sexual contacts
  # 4. convert back to rates

  if (mixing_matrix == "Zimbabwe") {
    m_age <- squire::get_mixing_matrix(country)
  } else if (mixing_matrix == "synthetic_all") {
    path <- system.file("extdata", "prem2021_synthetic_contact_all.rds",
                        package = "mpoxseir")
    m_age <- readRDS(path)
  } else if (mixing_matrix == "synthetic_home") {
    path <- system.file("extdata", "prem2021_synthetic_contact_home.rds",
                        package = "mpoxseir")
    m_age <- readRDS(path)
  } else {
    stop(sprintf("mixing_matrix %s not recognised", mixing_matrix))
  }
  
  ## need to do this with the squire age groups, rather than the new ones
  M_raw <- t(t(m_age) * N_age_squire) # total daily contacts in population
  
  
  ## balance the matrix so M[i,j] == M[j,i]
  M_age_squire <- (M_raw + t(M_raw)) / 2
  rownames(M_age_squire) <- colnames(M_age_squire) <- squire_age_bins$label

  # Adjust age groups to allow for new partition
  # [0-4], [5-9],  [10-14], [15-19], ... ->
  # [0-4], [5-11], [12-14], [15-19], ... ->
  
  # partition contacts within 10-14 category between 10-11 = 0.4 and 12-14 = 0.6
  M_12to14_12to14 <- M_age_squire["10-14", "10-14"] * 0.6^2
  M_10to11_12to14 <- M_age_squire["10-14", "10-14"] * 2 * 0.4 * 0.6
  M_10to11_10to11 <- M_age_squire["10-14", "10-14"] * 0.4^2
  
  # check equal
  stopifnot(M_10to11_10to11 + M_10to11_12to14 + M_12to14_12to14 ==
              M_age_squire["10-14", "10-14"])
  
  # partition contacts between 5-9 and 10-14
  M_05to09_10to11 <- M_age_squire["5-9", "10-14"] * (2 / 5) 
  M_05to09_12to14 <- M_age_squire["5-9", "10-14"] - M_05to09_10to11
  
  # compile contacts between 5-11 and 12-14 
  M_05to11_12to14 <- M_05to09_12to14 + M_10to11_12to14
  
  # compile contacts between 5-11
  M_05to11_05to11 <- M_age_squire["5-9", "5-9"] + M_10to11_10to11 +
    M_05to09_10to11

  # check that all have been accounted for
  stopifnot(abs(M_05to11_05to11 + M_05to11_12to14 + M_12to14_12to14 -
              (M_age_squire["5-9", "5-9"] + M_age_squire["5-9", "10-14"] +
              M_age_squire["10-14", "10-14"])) < 1e-6)
  
  M_10to11 <- M_age_squire["10-14", ] * 2 / 5
  M_12to14 <- M_age_squire["10-14", ] - M_10to11
  M_05to11 <- M_age_squire["5-9", ] + M_10to11
  
  M_age <- M_age_squire
  rownames(M_age) <- colnames(M_age) <- age_bins$label
  M_age["5-11", ] <- M_age[, "5-11"] <- M_05to11
  M_age["12-14", ] <- M_age[, "12-14"] <- M_12to14
  M_age[c("5-11", "12-14"), c("5-11", "12-14")] <-
    matrix(c(M_05to11_05to11, M_05to11_12to14,
             M_05to11_12to14, M_12to14_12to14), nrow = 2)
  
  M_age[lower.tri(M_age)] <- t(M_age)[lower.tri(M_age)] # populate lower triangle
  # check the totals match
  stopifnot(abs(sum(M_age[upper.tri(M_age, diag = TRUE)]) -
              sum(M_age_squire[upper.tri(M_age_squire, diag = TRUE)])) < 1e-6)

  
  # Need to split total contacts i->j by gen pop / CSW / ASW / PBS
  # assume homogenous mixing in day-to-day contacts (will add sex in fitting)
  # create a matrix with the probability of being in each key pop by age
  
  p_kp <- cbind(CSW = p_SW * w_CSW,
                ASW = p_SW * w_ASW,
                PBS = p_PBS * w_PBS,
                HCW = p_HCW * w_HCW)
  nms_kp <- colnames(p_kp)
  p_kp <- cbind(gen = 1 - rowSums(p_kp), p_kp)
  rownames(p_kp) <- age_bins$label

  ## For each i,j pair in M_age, split contacts up

  idx_compartment <- get_compartment_indices()
  nms_group <- names(idx_compartment$group)
  n_group <- idx_compartment$dim$group
  n_age <- nrow(age_bins)
  
  M <- matrix(0, n_group, n_group,  dimnames = list(nms_group, nms_group))
  
  ## Take the unique entries of M_age (upper triangular including diagonal) and
  ## Consider one by one, splitting into gen <-> gen / gen <-> key / key <-> key
  ## contacts for each i,j combo
  for (i in seq_len(n_age)) {
    for (j in seq(i, n_age)) {
      # Calculate proportion of age: i -> j contacts that are between gen / key pops
      # p sums to 1 and can be used to allocate all of M_age[i,j]
      # can consider p as the proportion of total contacts that involve
      # person 1 being in k and person 2 being in l, where k,l = gen, kp1, ...
      p <- outer(p_kp[i, ], p_kp[j, ])
      n <- M_age[i, j] * p
      # Separate out gen pop -> gen pop contacts and assign to age group
      M[i, j] <- n["gen", "gen"]
      # Count contacts between gen pop age group i and key pops aged j
      # Add to total contacts between age group i and each key pop. 
      # Repeat for age j as we only consider upper triangular so must allocate
      # all of M_age[i, j]. This holds for i=j as p separates this out into
      # 1.gen aged i -> 2.key aged i; and 1.key aged i -> 2.gen aged i
      M[i, nms_kp] <- M[i, nms_kp] + n["gen", nms_kp] 
      M[j, nms_kp] <- M[j, nms_kp] + n[nms_kp, "gen"]

      # Count contacts between key pops aged i and key pops aged j
      n_kp <- n[nms_kp, nms_kp]
      # add in reverse direction contacts (as described above), omitting diagonal
      # to avoid double counting
      M[nms_kp, nms_kp] <- M[nms_kp, nms_kp] + n_kp + t(n_kp * lower.tri(n_kp)) 

    }
  }

  M[lower.tri(M)] <- t(M)[lower.tri(M)] # populate lower triangle

  # check the totals match
  # stopifnot(abs(sum(M_age[upper.tri(M_age, diag = TRUE)]) -
  #                 sum(M[upper.tri(M, diag = TRUE)])) < 1e-6)
  
  
  
  
  # Convert to per-capita rates by dividing by population
  # Resulting matrix is Asymmetric c_ij != c_ji
  # BUT total number of contacts i->j and j->i is balanced
  m <- M / N
  # correct for any zero population denominators (e.g. HCW)
  m[N == 0, ] <- m[, N == 0] <- 0
  


  ## set up sexual contact matrix for parameterisation in transform function
  m_sex <- matrix(0, n_group, n_group, dimnames = list(nms_group, nms_group))

  # province populations
  province_pop = list("equateur" = 1712000,
                      "sudkivu" = 6565000,
                      "burundi" = 11890781, ## taken from squire (above)
                      "bujumbura" = 1095302,## Annuaire statisitique 2022 (from Olivier, in folder in Teams) - Bujumbura Mairie + Isare
                      "bujumbura_mairie" = 792503) ## purely for testing purposes

  # proportion of susceptibles estimated to be unvaccinated (historically)
  # In Burundi, no-one born after 1970 thought to be historically (smallpox) vaccinated (source: Ruth's email from Jean-Claude)
  # At time of writing (2025), this corresponds to over 55s
  p_unvaccinated <- setNames(rep(0, n_group), nms_group)
  if (region %in% c("equateur", "sudkivu")) {
    p_unvaccinated[which(age_bins$end < 40)] <- 1
    p_unvaccinated[which(age_bins$start >= 40)] <-
      c(0.54, 0.29, 0.29, 0.23, 0.21, 0.21, 0.21, 0.21)
  } else if (region %in% c("burundi","bujumbura","bujumbura_mairie")) {
    p_unvaccinated[which(age_bins$end < 55)] <- 1
    p_unvaccinated[which(age_bins$start >= 55)] <-
      c(0.23, 0.21, 0.21, 0.21, 0.21)
  }
  
  p_unvaccinated[nms_kp] <- 1 # assume no prior vaccination in KPs

  list(
    n_group = n_group,
    N0 = N,
    N_age = N_age,
    m_gen_pop = m,
    m_sex = m_sex,
    total_contacts_nonsexual = M,
    total_contacts_age = M_age,
    #total_contacts_sex = M_sex,
    n_vax = idx_compartment$dim$vax,
    p_unvaccinated = p_unvaccinated,
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
  start <- c(0, 5, 12, seq(15, 75, by = 5))
  create_age_bins(start, max_age = 100)
}

create_age_bins <- function(start, max_age = 100) {
  if (max_age <= max(start)) stop("max_age is too small")
  end <- c(start[-1] - 1L, max_age)
  label <- paste(start, end, sep = "-")
  label[length(label)] <- paste0(max(start), "+")
  data.frame(label = label, start = start, end = end)
}

##' A function that calculates the proportion of each age group than lies inside
##' a given interval
##' 
##' @title calculate the proportion of each age group that lies between `min_age`
##' and `max_age` (inclusive) based on uniform distribution within each age band
##' @param min_age a scalar giving the bottom of the age range
##' @param max_age a scalar giving the top of the age range
##' @return a vector of length n_age = 16 giving the proportion for each age band
##'   
proportion_in_age_bins <- function(min_age, max_age) {
  bins <- get_age_bins()
  
  # Calculate the overlap range for each bin
  overlap_start <- pmax(bins$start, min_age)
  overlap_end <- pmin(bins$end, max_age)
  
  # Calculate the size of the bin and the overlap
  bin_size <- bins$end - bins$start + 1
  overlap_size <- pmax(0, overlap_end - overlap_start + 1)
  
  # Calculate the proportion of overlap for each bin
  proportions <- overlap_size / bin_size
  
  proportions
}

##' A function that gets the compartment indices used in the model
##' 
##' @title Get compartment indices used throughout the model
##' 
##' @return A list containing entries: `dim` giving the dimensions of each
##' compartment, currently: group = 18, vax = 4;
##' `group`, a named list of the array indices corresponding to the
##' first dimension of model compartments: 16 5-year age bands + SW + PBS;
##' and `vax`, a named list of the vaccine strata used in
##' the second dimension of model compartments:
##' 1. historic smallpox; 2. unvaccinated; 3. one-dose; 4. two-dose.
##' 
##' Use this function wherever you need to refer to these standards in your code
##' to avoid duplication, and make modifying universally easier.
##'   
##' @export
get_compartment_indices <- function() {
  group_bins <- get_group_bins()
  n_group <- nrow(group_bins)

  groups <- seq_len(n_group)
  names(groups) <- group_bins$label
  
  n_vax <- 4
  vax_strata <- seq_len(n_vax)
  names(vax_strata) <- c("historic", "unvaccinated", "one_dose", "two_dose")
  
  list(dim = list(group = n_group, vax = n_vax),
       group = as.list(groups),
       vax = as.list(vax_strata))
}


## A function that indicates whether the group is classed as an adult or child, and for the border case what the prop is
get_group_bins <- function() {
  
  age_bins <- get_age_bins()
  age_bins$children <- proportion_in_age_bins(0, 11)
  age_bins$fifteen_plus <- 1 - proportion_in_age_bins(0,14)
  
  groups <- data.frame(label = c("CSW", "ASW", "PBS", "HCW"),
                       start = c(12, 18, 20, 20),
                       end = c(17, 49, 49, 69),
                       children = c(0, 0, 0, 0),
                       fifteen_plus = c(0.5,1,1,1))
  ret <- rbind(age_bins, groups)
  ret$adults <- 1 - ret$children
  
  ret
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
##' @inheritParams parameters_demographic
##'   
##' @param initial_infections The initial number of infections
##' @param use_ve_D logical, indicating whether model should allow for vaccine
##' efficacy against death (above and beyond protection against infection)
##' @param overrides A list, containing any parameters for which you want to
##'   override the default values
##' 
##' @return A list of the fixed parameters
##'   
##' @export
##' 
#' @export
parameters_fixed <- function(region, initial_infections, use_ve_D = FALSE,
                             mixing_matrix = "Zimbabwe", p_SW = NULL,
                             p_HCW = NULL,
                             overrides = list()) {

  ## Checking region
  if (!(region %in% c("equateur", "sudkivu",
                      "burundi","bujumbura","bujumbura_mairie"))) {
    stop("region must be equateur, sudkivu, burundi, bujumbura or bujumbura_mairie")
  }

  ## Initialising variable that other parameters depend on
  demographic_params <- parameters_demographic(region = region,
                                               mixing_matrix = mixing_matrix,
                                               p_SW = p_SW,
                                               p_HCW = p_HCW)
  age_bins <- get_age_bins()
  idx_compartment <- get_compartment_indices()

  n_group <- demographic_params$n_group
  n_vax <- demographic_params$n_vax
  idx_unvax <- idx_compartment$vax$unvaccinated
  idx_historic_vax <- idx_compartment$vax$historic
  
  ## standard vaccination parameters, we may want to move this entire section to
  ## the transform in future if we want to allow for vaccine uncertainty
  ## VE against infection
  ve_I <- matrix(c(0.736, 0, 0.736, 0.818),
                 nrow = n_group, ncol = n_vax, byrow = TRUE) ##VALUES TO BE UPDATED
  ## VE against onward transmission
  ve_T <- rep(0, n_vax)
  
  
  N <- demographic_params$province_pop[[region]]
  N0 <- round(N * demographic_params$N0 / sum(demographic_params$N0)) # total number in each age-group
  ## ages: 0-4, 5-11, 12-14, 15-19, 20+
  RR_z <- c(1, 0.857, 0.438, 0.086, rep(0.073, n_group - 4)) # Jezek 1988 zoonotic + Jezek 1987
  

  ## Seed infections in the unvaccinated group in a region-specific manner
  X0 <- matrix(0, nrow = n_group, ncol = n_vax)
  seed_rate <- X0

  
  if (region %in% c("sudkivu","burundi","bujumbura","bujumbura_mairie")) { # seeding in sex workers in Clade Ib affected areas

    ## Extract sex-worker index and put initial infections in this group (unvaccinated strata)
    index_asw <- get_compartment_indices()$group$ASW
    seed_rate[index_asw, idx_unvax] <- initial_infections

  } else if (region == "equateur") { # seeding in general pop in proportion to zoonotic risk in equateur

    ## Extract gen-pop index and put initial infections in this group (unvaccinated strata) in proportion to zoonotic risk
    index_gen_pop <- seq_len(nrow(age_bins))
    seeding_infections <- initial_infections * 
      RR_z[index_gen_pop] / sum(RR_z[index_gen_pop])

    seed_rate[index_gen_pop, idx_unvax] <- seeding_infections

  }
  
  p_unvaccinated <- demographic_params$p_unvaccinated
  ## update with assignment into smallpox vaccination compartment (j = 1)
  ## if we have a small population then seeding with 5 infections per compt and
  ## strata could be a problem (e.g. with N = 10,000 too few SW)
  ## split population into historically- / un-vaccinated first before seeding
  ## need to think about how this will interact with equilibrium position for
  ## seeding in endemic areas.
  S0 <- X0
  S0[, idx_unvax] <- round(N0 * p_unvaccinated)
  S0[, idx_historic_vax] <- N0 - S0[, idx_unvax]
  if (any (S0 < seed_rate)) {
    stop("population size and seeding infections is incompatible")
  }

  CFR_unvax <- rep(0, n_group)
  names(CFR_unvax) <- names(demographic_params$N0)
  # CFR from Whittles 2024, 5-year bands to 40
  CFR_unvax[which(age_bins$end < 40)] <-
    c(0.102, 0.054, 0.035, 0.026, 0.02, 0.016, 0.013, 0.012)
  CFR_unvax[which(age_bins$start >= 40)] <- 0.01

  # Create matrix [group | vax strata]
  # Vax strata: 1. historic smallpox; 2. unvaccinated; 3. one-dose; 4. two-dose
  # initially all strata contain CFR by age for unvaccinated, which we then
  # modify to incorporate any VE against death 
  CFR <- matrix(CFR_unvax, nrow = n_group, ncol = n_vax)
  rownames(CFR) <- names(CFR_unvax)
  
  if (use_ve_D) {
  # If allowing for vaccine protection against death
    CFR_historic_vax <- CFR_unvax
    CFR_historic_vax[which(age_bins$end < 40)] <-
      c(0.048, 0.025, 0.016, 0.012, 0.009, 0.007, 0.006, 0.005)
    CFR_historic_vax[which(age_bins$start >= 40)] <- 0.004

    # use same CFR for all vaccinated regardless of efficacy
    CFR[, -idx_unvax] <- CFR_historic_vax
  }
  
  group_bins <- get_group_bins()
  
  CFR["ASW", ] <- colMeans(CFR[which(age_bins$start>=group_bins$start[which(group_bins$label=="CSW")]&age_bins$end<=group_bins$end[which(group_bins$label=="ASW")]),])
  CFR["CSW", ] <- CFR["ASW", ]
  CFR["PBS", ] <- colMeans(CFR[which(age_bins$start>=group_bins$start[which(group_bins$label=="PBS")]&age_bins$end<=group_bins$end[which(group_bins$label=="PBS")]),])
  CFR["HCW", ] <- colMeans(CFR[which(age_bins$start>=group_bins$start[which(group_bins$label=="HCW")]&age_bins$end<=group_bins$end[which(group_bins$label=="HCW")]),])
  
  ## vaccination default 
  
  ## Initially final coverage within children / adults is set to be the proportion
  ## in each age group (i.e. 100% coverage, given age). For example, 15-19yo 
  ## would have 60% target coverage of the child vaccines and 40% of the adult
  ## This may cause issues in future if vaccines that can be given to adults / 
  ## children have different efficacies.
  
  N_prioritisation_steps_children <- 1
  N_prioritisation_steps_adults <- 1
  
  group_bins <- get_group_bins()
  prioritisation_strategy_children <- matrix(group_bins$children,
                                             nrow = n_group,
                                             ncol = N_prioritisation_steps_children)
  prioritisation_strategy_adults <- matrix(group_bins$adults,
                                             nrow = n_group,
                                             ncol = N_prioritisation_steps_adults)
  
  params_list = list(
    region = region,
    n_group = n_group,
    n_vax = n_vax,
    S0 = S0,
    Ea0 = X0,
    Eb0 = X0,
    Ir0 = X0,
    Id0 = X0,
    R0 = X0,
    D0 = X0,
    N0 = N0,
    seed_rate = seed_rate,
    R0_hh = 0.67, # Jezek 1988 SAR paper - will be fitted
    R0_sw_st = 1.3, # Will be fitted
    beta_z_max = 0.01, # Will be fitted
    beta_hcw = 0.001, # Will be fitted
    alpha_cases = 0.5,
    alpha_cases_00_04 = 0.5,
    alpha_cases_05_14 = 0.5,
    alpha_cases_15_plus = 0.5,
    alpha_deaths = 0.5,
    alpha_deaths_00_04 = 0.5,
    alpha_deaths_05_14 = 0.5,
    alpha_deaths_15_plus = 0.5,
    phi_00_04 = 1,
    phi_05_14 = 1,
    phi_15_plus = 1,
    phi_CSW_12_14 = 1,
    phi_CSW_15_17 = 1,
    phi_ASW = 1,
    phi_PBS = 1,
    phi_HCW = 1,
    rho_00_04 = 0.1,
    rho_00_14 = 0.1,
    RR_z = RR_z,
    gamma_E = 1 / 7,  #  1/7 based on Besombes et al. on 29 clade I patients
    gamma_Ir = 1 / 18, # Jezek 1988 "clinical features of 282.."
    gamma_Id = 1 / 10, # Jezek 1988
    CFR = as.matrix(CFR),
    m_sex = demographic_params$m_sex,
    m_gen_pop = demographic_params$m_gen_pop,
    N_prioritisation_steps_children = N_prioritisation_steps_children,
    N_prioritisation_steps_adults = N_prioritisation_steps_adults,
    prioritisation_strategy_children = prioritisation_strategy_children,
    prioritisation_strategy_adults = prioritisation_strategy_adults,
    ve_I = ve_I,
    ve_T = ve_T,
    daily_doses_children_value = matrix(0, nrow = n_vax, ncol = 1),
    daily_doses_children_time = 1,
    daily_doses_adults_value = matrix(0, nrow = n_vax, ncol = 1),
    daily_doses_adults_time = 1,
    is_child = group_bins$children)

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
