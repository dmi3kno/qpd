% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chebyshev.R
\name{check_roots}
\alias{check_roots}
\alias{is_qdf_valid}
\title{Validate quantile function by checking the corresponding quantile density function}
\usage{
check_roots(fun, ..., n = 13, interval = c(0, 1), parts = 10)

is_qdf_valid(fun, ..., n = 13L, parts = 10)
}
\arguments{
\item{fun}{quantile density function (QDF) to check}

\item{...}{other parameters passed to QDF function}

\item{n}{integer polynomial degree to fit to the function. Default is 13}

\item{interval}{interval to check the roots on}

\item{parts}{integer number of parts to split the function range. Default is 10}
}
\value{
logical value indicating if the QDF is valid (i.e. returns only non-negative values)
}
\description{
Validate quantile function by checking the corresponding quantile density function
}
