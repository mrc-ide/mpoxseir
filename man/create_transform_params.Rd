% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run.R
\name{create_transform_params}
\alias{create_transform_params}
\title{Create transform function}
\usage{
create_transform_params(region, initial_infections, overrides = list())
}
\arguments{
\item{region}{The region for fitting, must be either \code{"equateur"} or
\code{"sudkivu"}}

\item{initial_infections}{The initial number of infections}

\item{overrides}{A list, containing any parameters for which you want to
override the default fixed values}
}
\value{
A transform function
}
\description{
A function that creates the transform function for use in the fitting
}
