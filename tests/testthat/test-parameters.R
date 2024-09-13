test_that("assign seeds works", {
  tmp <- assign_seeds(10, rep(1 / 3, 3))
  expect_equal(sum(tmp), 10)
  expect_equal(assign_seeds(6, c(3, 2, 1)), c(3, 2, 1))
})
