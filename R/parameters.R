
#' @importFrom stats dmultinom
#' @importFrom stats setNames
#' @importFrom squire get_population
#' @importFrom squire get_mixing_matrix
#' @export
parameters_demographic <- function() {
  age_bins <- get_age_bins()

  ## Set up population denominators
  country <- "Democratic Republic of Congo"
  data <- squire::get_population(country)
  N_age <- c(data$n[1:15], sum(data$n[16:17])) # combine 75+

  ## Add in Sex workers (SW) and People who Buy Sex (PBS)
  idx_15_49  <- age_bins$start >= 15 & age_bins$end < 49
  N_15_49 <- N_age * idx_15_49

  p_SW <- 0.007 * 0.5 # 0.7% women (50%) 15-49 Laga et al
  N_SW <- round(p_SW * N_15_49)

  p_PBS <- 0.11 * 0.5 # 11% men (50%) 15-49 DHS https://www.statcompiler.com/en/
  N_PBS <- round(p_PBS * N_15_49)

  N <- c(N_age - N_SW - N_PBS, sum(N_SW), sum(N_PBS))

  # Set up mixing matrix
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
  M0 <- M_age * outer(!idx_15_49, !idx_15_49) # neither could be KP
  M1 <- M_age * outer(idx_15_49, idx_15_49, FUN = "xor") # only 1 could be KP
  M2 <- M_age * outer(idx_15_49, idx_15_49) # both could be KP


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
  sum(M_gen_pop + M_gen_PBS + M_gen_SW + M_SW_SW + M_PBS_PBS + M_SW_PBS - M_age)

  marginalise <- function(M) (rowSums(M) + diag(M)) / 2

  n_age <- length(N_age)
  idx_age <- seq_len(n_age)
  n_group <- n_age + 2
  nms_group <- c(age_bins$label, "SW", "PBS")

  # Construct new contact matrix including groups
  M <- matrix(0, n_group, n_group, dimnames = list(nms_group, nms_group))
  M[idx_age, idx_age] <- M_gen_pop
  M["SW", idx_age] <- M[idx_age, "SW"] <- marginalise(M_gen_SW)
  M["PBS", idx_age] <- M[idx_age, "PBS"] <- marginalise(M_gen_PBS)

  M["SW", "SW"] <- sum(marginalise(M_SW_SW))
  M["PBS", "PBS"] <- sum(marginalise(M_PBS_PBS))
  M["SW", "PBS"] <- M["PBS", "SW"] <- sum(marginalise(M_SW_PBS))

  # check that total contacts are the same as original
  sum(marginalise(M)) - sum(marginalise(M_age)) ## check total number of contacts

  # Convert to per-capita rates by dividing by population
  # Resulting matrix is Asymmetric c_ij != c_ji
  # BUT total number of contacts i->j and j->i is balanced
  m <- M / N

  # proportion of susceptibles estimated to have smallpox vaccine
  sus_prop <- c(rep(1,8),0.54,0.29,0.29,0.23,0.21,0.21,0.21,0.21,1,1)

  # province populations
  province_pop = list("equateur" = 1712000,
                      "sudkivu" = 6565000)


  list(
    n_group = n_group,
    N0 = setNames(N, nms_group),
    m = m,
    total_contacts = M,
    n_vax = 2,
    sus_prop = sus_prop,
    province_pop = province_pop
  )
}

## We always use these age bands, so rather than detect them, we will
## check that things conform to them.
get_age_bins <- function() {
  end <- c(seq(4, 75, by = 5), 100)
  start <- c(0, end[-length(end)] + 1L)
  label <- paste(start, end, sep = "-")
  label[length(label)] <- paste0(max(start), "+")

  data.frame(label = label, start = start, end = end)
}

