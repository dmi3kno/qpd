% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gnh.R
\name{qgnh}
\alias{qgnh}
\alias{fgnh}
\alias{dqgnh}
\alias{rgnh}
\alias{pgnh}
\alias{is_gnh_valid}
\title{Functions for g-and-h distribution}
\usage{
qgnh(p, A, B, C = 0.8, g, h, zscale = FALSE)

fgnh(p, A, B, C = 0.8, g, h, log = FALSE, zscale = FALSE)

dqgnh(p, A, B, C = 0.8, g, h, log = FALSE, zscale = FALSE)

rgnh(n, A, B, C = 0.8, g, h)

pgnh(q, ..., lower = 0, upper = 1, tol = 1e-06, silent = TRUE, trace = 0)

is_gnh_valid(A, B, C = 0.8, g, h, n_grid = 100L, s_grid = 2L)
}
\arguments{
\item{p}{numeric vector of probabilities}

\item{A}{numeric location parameter of g-and-h distribution}

\item{B}{positive numeric scale parameter of g-and-h distribution}

\item{C}{numeric parameter of g-and-h distribution. Default is 0.8}

\item{g}{numeric skeweness parameter of g-and-h distribution}

\item{h}{non-negative numeric kurtosis parameter of g-and-h distribution}

\item{zscale}{is the probability input provided on scale of \code{z} values, i.e. as \code{z=qnorm(p)}. Default is FALSE}

\item{log}{logical should the result be returned as log(). Default is FALSE}

\item{n}{numeric number of samples to draw}

\item{q}{vector of quantiles}

\item{...}{used by method}

\item{lower, upper}{the \code{stats::uniroot} lower and upper end points of the interval to be searched. Defaults are 0 and 1, respectively}

\item{tol}{the \code{stats::uniroot} desired accuracy (convergence tolerance). Default value 1e-06}

\item{silent}{the \code{base::try} argument. Default is TRUE}

\item{trace}{integer number passed to \code{stats::uniroot}; if positive, tracing information is produced. Higher values giving more details.}

\item{n_grid}{integer size of helper grid to be passed to \code{make_pgrid}. Default is 50}

\item{s_grid}{integer beta shape of helper grid to be passed to \code{make_pgrid}. Default is 5}
}
\value{
vector of quantiles for g-and-h distribution

vector of quantiles for g-and-h distribution
}
\description{
Functions for g-and-h distribution
}
\examples{
qgnh(0.1, 3,1,0.8, 2,0.5)
is_gnh_valid(A=5, B=5, C=0.8, g=0.5, h=0.5)
}
