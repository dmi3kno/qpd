# qpd 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
* Added several distributions 
 - Metalog distribution (general metalog)
 - Keelin-Powley Simple Q-Normal distribution
 - Johnson Quantile Parameterized Distribution (JQPDB, JQPDS, JQPDS2)
 - Generalized Normal distribution with skew parameter (QPD)
 - Govindarajulu distribution
 - Generalized Lambda Distribution (RS parameterization)
 - Generalized Lambda Distribution (FKML parameterization)
 - Generalized Lambda Distribution (CSW parameterization)
 - Skew-Logistic Distribution
 - Flattened (Skew-) Logistic Distribution
 - Wakeby distribution
 - Exponential distribution (quantile density functions)
 - g-and-h distribution, g-and-k distribution
 - Singe and two-parameter Rayleigh distribution
 - Normal distribution (quantile density functions)
 - Generalized Exponential distribution
 - Myerson distribution and its variations (logit-Myerson, sech-Myerson and cauchy-Myerson)
 - Quantile-parameterized Skew-Logistic xPower distribution (x="exteded")
* Function for fitting beta distribution to 3 elicited quantiles
* Service functions for accumulating simplices
* HDR pseudo random numbers algorithm by Hubbard Decision Research
* Service functions for eliciting QDirichlet distributions
* Service functions for checking the validity of quantile functions by proxy root finding with Chebyshev polynomials
* Added generic inverse quantile function factory `iqf` for creating approximated CDFs for quantile distributions.
* Change default metalog bounds to (-Inf,Inf)
* Add `make_ecdf_df()` for creating ECDF from sample.
