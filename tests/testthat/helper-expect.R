expect_vector_lt <- function(x, y, digits = 100, tol = 0) {
  expect_true(all(round(x, digits) - round(y, digits) < tol),
              sprintf("\nNot all '%s' less than '%s' (tol %s)",
                      deparse(substitute(x)), deparse(substitute(y)), tol))
}
