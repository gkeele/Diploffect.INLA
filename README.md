Diploffect.INLA
===============

This package contains an importance sampling-based integrated nested Laplace approximation (INLA) implementation of the Diploffect model to estimate Bayesian posterior distributions of genetic effects while taking into account haplotype uncertainty.

# Required (or useful) non-CRAN R packages

This package requires the [INLA](http://www.r-inla.org/) R package, which is not available through CRAN. Visit the site, or more simply enter the following code at the prompt in R:
```r
install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/stable")
```  
If the haplotype data are stored in a genome cache, Diploffect.INLA can make use of convenience utilities contained in bagpipe.backend, the R backend to the [Bagpipe](http://valdarlab.unc.edu/software/bagpipe/_build/html/index.html) software that can be used to map QTL in multiparent populations.

# Example 

Example data is included in the package. The following code services as a simple vignette for using the package for now.

```r
library(devtools)
install_github(“gkeele/Diploffect.INLA”)
library(Diploffect.INLA)
data(exampleCC)
data(locusmatrix)
inla.diploffect <- run.diploffect.inla(formula=y~1+(1|strain)+(1|dose.date), add.on=FALSE, 				       
                                       data=exampleCC, K=NULL, prob.matrix=locusmatrix,
                                       num.draws=10, use.dip.lincomb=TRUE, seed=1, 
				       gamma.rate=1, impute.on="CCline")
inla.diploffect.summary <- run.diploffect.inla.summary.stats(inla.diploffect)
plot_straineff.ci(inla.diploffect.summary, flip=FALSE)
plot_deviation.ci(inla.diploffect.summary, flip=FALSE)
plot_diplotype.ci(inla.diploffect.summary, flip=FALSE)
plot_varexp.ci(inla.diploffect.summary, add.numbers=TRUE)
```
