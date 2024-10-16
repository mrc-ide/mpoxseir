test_that("mpoxseir_date can convert to dates into 2023", {
  expect_equal(mpoxseir_date("2023-01-01"), 1)
  expect_equal(mpoxseir_date("2024-12-31"), 365 + 366) # 2024 is a leap year!
  r <- seq(as_date("2023-01-01"), as_date("2024-12-31"), by = 1)
  expect_equal(mpoxseir_date(r), 1:(365 + 366))
  
  expect_equal(mpoxseir_date_as_date(mpoxseir_date(r)), r)
  
  expect_error(mpoxseir_date("2019-01-01"),
               "Negative dates, mpoxseir_date likely applied twice")
  expect_error(mpoxseir_date(c("2023-01-01", "2019-01-01")),
               "Negative dates, mpoxseir_date likely applied twice")
})


test_that("assert mpoxseir date throws on non mpoxseir dates", {
  expect_silent(assert_mpoxseir_date(1))
  expect_error(assert_mpoxseir_date(as_date("2023-01-01")),
               "'date' must be numeric - did you forget mpoxseir_date()?",
               fixed = TRUE)
})


test_that("helper function can convert somewhat helpfully", {
  expect_equal(as_mpoxseir_date("2023-02-01"), 32)
  expect_equal(as_mpoxseir_date(as_date("2023-02-01")), 32)
  expect_equal(as_mpoxseir_date(32), 32)
})


test_that("negative numbers are not allowed", {
  expect_error(assert_mpoxseir_date(-1),
               "Negative dates, mpoxseir_date likely applied twice")
  expect_error(assert_mpoxseir_date(c(10, -1, 29)),
               "Negative dates, mpoxseir_date likely applied twice")
})

test_that("ISO dates only", {
  expect_error(
    as_date("02-05-2024"),
    "Expected ISO dates or R dates - please convert")
})
