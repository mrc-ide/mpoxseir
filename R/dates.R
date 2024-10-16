##' We need to map "dates" onto [`dust::dust`]'s concept of model
##' "step" and we do this by mapping a date such as `2024-03-02` into
##' the number of days into (and beyond) 2023 (62 here, with the 1st of January
##' being day 1). We call this integer number an "mpoxseir date".
##'
##' There are several related functions here
##'
##' * `mpoxseir_date` converts its argument into an R `Date` object,
##'   then applies this tranformation. If the argument is not a `Date`
##'   object or a string representing one, an error will be thrown.
##'
##' * `mpoxseir_date_to_date` does the reverse conversion to
##'   `mpoxseir_date`, converting an integer mpoxseir date into an R
##'   `Date`
##'
##' * `as_mpoxseir_date` does the same conversion as `mpoxseir_date`
##'   but will assume that an integer *already* represents an mpoxseir
##'   date and will return it unmodified rather than erroring.
##'
##' * `as_date` does a string to date conversion, using [as.Date()]
##'   but requiring the dates are in ISO 8601 (YYYY-MM-DD) format (it
##'   is a helper that avoids conversion to `NA`, instead throwing an
##'   error)
##'
##' @title Date handling for mpoxseir
##'
##' @param date A Date object, or something that can be converted to
##'   one, or an "mpoxseir date"; see Details
##'
##' @return An integer, being the number of days into 2023
##' @export
##' @examples
##' # Convert dates into mpoxseir dates:
##' mpoxseir::mpoxseir_date("2023-01-01")
##' mpoxseir::mpoxseir_date(c("2024-03-01", "2024-10-01"))
##'
##' # Reverse the conversion:
##' mpoxseir::mpoxseir_date_as_date(1)
##' mpoxseir::mpoxseir_date_as_date(c(61, 275))
##'
##' # Double conversion not possible with mpoxseir_date...
##' try(mpoxseir::mpoxseir_date(61))
##' # ...but allowed with as_mpoxseir_date
##' mpoxseir::as_mpoxseir_date(61)
##'
##' # Strict date conversion with as_date
##' mpoxseir::as_date("2024-03-01")
##' try(mpoxseir::as_date("03-01-2024"))
mpoxseir_date <- function(date) {
  days_into_2023 <- as.numeric(as_date(date) - as_date("2022-12-31"))
  if (any(days_into_2023 < 0)) {
    stop("Negative dates, mpoxseir_date likely applied twice")
  }
  days_into_2023
}


##' @export
##' @rdname mpoxseir_date
mpoxseir_date_as_date <- function(date) {
  assert_mpoxseir_date(date)
  as_date("2022-12-31") + date
}


##' @export
##' @rdname mpoxseir_date
as_mpoxseir_date <- function(date) {
  if (is.character(date) || is_date(date)) {
    mpoxseir_date(as_date(date))
  } else {
    assert_mpoxseir_date(date)
  }
}

##' @export
##' @rdname mpoxseir_date
as_date <- function(date) {
  if (is_date(date)) {
    return(date)
  }
  if (!all(grepl("^[0-9]{4}-[0-9]{1,2}-[0-9]{1,2}$", date))) {
    stop("Expected ISO dates or R dates - please convert")
  }
  as.Date(date)
}


assert_mpoxseir_date <- function(date) {
  if (!is.numeric(date)) {
    stop("'date' must be numeric - did you forget mpoxseir_date()?")
  }
  if (any(date < 0)) {
    stop("Negative dates, mpoxseir_date likely applied twice")
  }
  date
}

is_date <- function(x) {
  inherits(x, "Date")
}
