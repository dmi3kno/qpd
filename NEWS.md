# qpd 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
* Added several distributions 
 - Metalog distribution (general metalog)
 - Johnson Quantile Parameterized Distribution (JQPDB, JQPDS, JQPDS2)
 - Generalized Normal distribution with skew parameter (QPD)
 - Govindarajulu distribution
 - Generalized Lambda Distribution (RS parameterization)
 - (Generalized) Flattened Logistic Distribution
 - Wakeby distribution
 - Exponential distribution (quantile density functions)
 - g-and-h distribution, g-and-k distribution
 - two-parameter Rayleigh distribution
 - Normal distribution (quantile density functions)
 - Generalized Exponential distribution
 - Myerson distribution
* Function for fitting beta distribution to 3 elicited quantiles
* Service functions for accumulating simplices
* HDR pseudo random numbers algorithm by Hubbard Decision Research
* Service functions for eliciting QDirichlet distributions
* Service functions for checking the validity of quantile functions by proxy root finding with Chebyshev polynomials
* Added generic inverse quantile function factory `iqf` for creating approximated CDFs for quantile distributions.
