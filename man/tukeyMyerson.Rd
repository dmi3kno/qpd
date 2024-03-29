% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Myerson.R
\name{tukeyMyerson}
\alias{tukeyMyerson}
\alias{qtukeyMyerson}
\alias{ftukeyMyerson}
\alias{dqtukeyMyerson}
\alias{rtukeyMyerson}
\alias{ptukeyMyerson}
\title{Density function for Tukey-Myerson distribution (Myerson distribution with Tukey lambda base function).}
\usage{
qtukeyMyerson(
  p,
  q1,
  q2,
  q3,
  alpha = 0.25,
  tlambda = 0,
  lower.tail = TRUE,
  log.p = FALSE
)

ftukeyMyerson(
  p,
  q1,
  q2,
  q3,
  alpha = 0.25,
  tlambda = 0,
  log = FALSE,
  lower.tail = TRUE,
  log.p = FALSE
)

dqtukeyMyerson(
  p,
  q1,
  q2,
  q3,
  alpha = 0.25,
  tlambda = 0,
  log = FALSE,
  lower.tail = TRUE,
  log.p = FALSE
)

rtukeyMyerson(n, q1, q2, q3, alpha = 0.25, tlambda = 0)

ptukeyMyerson(
  q,
  ...,
  lower = 0,
  upper = 1,
  tol = 1e-06,
  silent = TRUE,
  trace = 0
)
}
\arguments{
\item{p}{vector of probabilities}

\item{q1, q2, q3}{numeric values representing bottom quantile, middle quantile (50th percentile)
and top quantile. Quantiles are assumed to be symmetrical.}

\item{alpha}{numerical fixed proportion of distribution below the bottom quantile.
Default value is 0.25}

\item{tlambda}{numerical Tukey lambda \code{lambda} parameter}

\item{lower.tail}{logical; if TRUE (default), probabilities are \eqn{P[X\leq x]}, otherwise, \eqn{P[X > x]}}

\item{log, log.p}{logical; if TRUE, probabilities p are given as log(p)}

\item{n}{number of observations. If `length(n)>1``, the length is taken to be the number required.}

\item{q}{vector of quantiles.}

\item{...}{used by method}

\item{lower, upper}{the \code{stats::uniroot} lower and upper end points of the interval to be searched. Defaults are 0 and 1, respectively}

\item{tol}{the \code{stats::uniroot} desired accuracy (convergence tolerance). Default value 1e-06}

\item{silent}{the \code{base::try} argument. Default is TRUE}

\item{trace}{integer number passed to \code{stats::uniroot}; if positive, tracing information is produced. Higher values giving more details.}
}
\value{
a vector of probabilities of length equal to \code{length(x)}

a vector of exceedance probabilities of length equal to \code{length(x)}.
}
\description{
Density function for Tukey-Myerson distribution (Myerson distribution with Tukey lambda base function).
}
\examples{
qtukeyMyerson(0.25, 10, 20, 40)
flogitMyerson(0.25, 10, 20, 40)
dqlogitMyerson(0.25, 10, 20, 40)
rlogitMyerson(1, 10, 20, 40)
ptukeyMyerson(0.25, 10, 20, 40)
}
