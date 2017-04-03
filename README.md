Diploffect.INLA
===============

This package contains an importance sampling integrated nested Laplace approximation (INLA) implementation of the Diploffect model to estimate Bayesian posterior distributions of genetic effects while taking into account haplotype uncertainty. 

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
