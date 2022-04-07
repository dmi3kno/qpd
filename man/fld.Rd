% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fld.R
\name{fld}
\alias{fld}
\alias{qfld}
\alias{ffld}
\alias{dqfld}
\alias{pfld}
\title{Flattened Generalized Logistic Distribution (FLD)}
\usage{
qfld(p, bt, k, dlt = 0.5)

ffld(p, bt, k, dlt = 0.5)

dqfld(p, bt, k, dlt = 0.5)

pfld(q, bt, k, dlt = 0.5, tol = 1e-06)
}
\arguments{
\item{p}{vector of probabilities}

\item{bt}{parameter beta of FGLD distribution (scale). Should be non-negative.}

\item{k}{parameter k of FGLD distribution (shape). Should be non-negative.}

\item{dlt}{parameter delta of FGLD distribution(mixing parameter). Should be within the interval (0,1), default is 0.5}

\item{q}{vector of quantiles}

\item{tol}{tolerance value for optimization. Default value 1e-06}
}
\value{
vector
}
\description{
Defines quantile function (ICDF), probability density and quantile density function as well as probability function (CDF) for flattened generalized logistic distribution
}
\examples{
qfld(0.1, 0.5, 0.3, 0.5)
# centered fld
p <- runif(1e4)
x <- qfld(p, 0.25, 1)-qfld(0.5, 0.25, 1)
hist(x,30)
}