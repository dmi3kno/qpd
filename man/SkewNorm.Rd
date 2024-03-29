% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SkewNorm.R
\name{dSkewNorm}
\alias{dSkewNorm}
\alias{pSkewNorm}
\alias{qSkewNorm}
\alias{rSkewNorm}
\title{Generalized Normal distribution with skew parameter.}
\usage{
dSkewNorm(x, q1, q3, s = 0)

pSkewNorm(x, q1, q3, s = 0)

qSkewNorm(p, q1, q3, s = 0)

rSkewNorm(n, q1, q3, s = 0)
}
\arguments{
\item{x}{vector of quantiles. If `length(n)>1``, the length is taken to be the number required.}

\item{q1}{minimum value}

\item{q3}{maximum value}

\item{s}{skew parameter betweem -1 and 1. When skew parameter equal 0, distribution is normal.}

\item{p}{numeric vector of probabilities}

\item{n}{number of observations. If `length(n)>1``, the length is taken to be the number required.}
}
\value{
a vector of probabilities of length equal to \code{length(x)}.

a vector of exceedance probabilities of length equal to \code{length(x)}.

a vector of exceedance probabilities of length equal to \code{length(x)}.

a length \code{n} vector of random values.
}
\description{
Generalized Normal distribution with skew parameter.
}
\examples{
dSkewNorm(12, 10, 20, 0.5)
pSkewNorm(12, 10, 20, 0.5)
qSkewNorm(12, 10, 20, 0.5)
rSkewNorm(1, 10, 20, 0.5)
}
