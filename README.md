
# Diploffect.INLA

This package contains an importance sampling-based integrated nested
Laplace approximation (INLA) implementation of the Diploffect model to
estimate Bayesian posterior distributions of genetic effects while
taking into account haplotype uncertainty.

## Required (or useful) non-CRAN R packages

This package requires the [INLA](http://www.r-inla.org/) R package,
which is not available through CRAN. Visit the site, or more simply
enter the following code at the prompt in R: In addition to phenotype
values, the Diploffect package requires that each individual has a set
of diplotype (founder haplotype pair identity) probabilities for the
locus. These diplotype probabilities can be provided directly as a
matrix, or, for backwards compatibility, they can be provided as a data
directory as would be output from the haplotype reconstruction method
HAPPY (Mott *et al.* 2000), henceforth described as a genome cache. If
the diplotype probabilities are stored in a genome cache,
Diploffect.INLA can make use of convenience utilities to run the
Diploffect model.

## Example

Example data is included in the package. The following code services as
a simple vignette for using the package for now.

``` r
library(Diploffect.INLA)
data(exampleCC)
data(locusmatrix)
```

``` r
inla.diploffect <- run.diploffect.inla(formula=y~1+(1|strain)+(1|dose.date), add.on=FALSE,                     
                                       data=exampleCC, K=NULL, prob.matrix=locusmatrix,
                                       num.draws=10, use.dip.lincomb=TRUE, seed=1, 
                       gamma.rate=1, impute.on="CCline")
```

    ## 1: -236.837262206168 
    ## 2: -238.057618245916 
    ## 3: -237.342529220994 
    ## 4: -236.848450287613 
    ## 5: -236.848450287613 
    ## 6: -236.400242074085 
    ## 7: -236.837262206168 
    ## 8: -236.400242074085 
    ## 9: -237.342529220994 
    ## 10: -237.302667466497

``` r
inla.diploffect.summary <- run.diploffect.inla.summary.stats(inla.diploffect)
```

    ## Loading nonqtl proportion of variance explained (var components and sums of squares)... 
    ## Loading strain effects...
    ## Loading additive heritability...
    ## Loading additive sums of squares...
    ## Loading deviation effects...
    ## Loading deviation means...
    ## Loading diplotype means...
    ## Loading dominant heritability...
    ## Loading total QTL heritability...
    ## Loading dominant sums of squares...
    ## Loading total QTL sums of squares...

``` r
plot.straineff.ci(inla.diploffect.summary)
```


