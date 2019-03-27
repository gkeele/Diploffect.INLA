Diploffect.INLA
===============

This package contains an importance sampling-based integrated nested Laplace approximation (INLA) implementation of the Diploffect model to estimate Bayesian posterior distributions of genetic effects while taking into account haplotype uncertainty.

## Required (or useful) non-CRAN R packages

This package requires the [INLA](http://www.r-inla.org/) R package, which is not available through CRAN. Visit the site, or more simply enter the following code at the prompt in R:
```r
install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE) 
```  
In addition to phenotype values, the Diploffect package requires that each individual has a set of diplotype (founder haplotype pair identity) probabilities for the locus. These diplotype probabilities can be provided directly as a matrix, or, for backwards compatibility, they can be provided as a data directory as would be output from the haplotype reconstruction method HAPPY (Mott *et al.* 2000), henceforth described as a genome cache. If the diplotype probabilities are stored in a genome cache, Diploffect.INLA can make use of convenience utilities to run the Diploffect model.

## Example 

Example data is included in the package. The following code services as a simple vignette for using the package for now.

```r
library(devtools)
install_github("gkeele/Diploffect.INLA")
library(Diploffect.INLA)
library(INLA) # This is obnoxious - for some reason import within Diploffect.INLA does not work
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
