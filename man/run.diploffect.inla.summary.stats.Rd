% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run.diploffect.inla.R
\name{run.diploffect.inla.summary.stats}
\alias{run.diploffect.inla.summary.stats}
\title{Estimate approximate posterior credible intervals for model parameters from Diploffect.INLA posterior distributions}
\usage{
run.diploffect.inla.summary.stats(diploffect.object)
}
\arguments{
\item{diploffect.object}{A Diploffect.INLA summary object. Output from run.diploffect.inla.summary.stats(). Contains credible information
for all the model parameters.}
}
\value{
Diploffect.INLA summary object.
}
\description{
This function takes a Diploffect.INLA object and estimates the posterior credible intervals of model parameters
}
\examples{
library(Diploffect.INLA)
data(exampleCC)
data(locusmatrix)
inla.diploffect <- run.diploffect.inla(formula=y~1+(1|strain)+(1|dose.date), add.on=FALSE, data=exampleCC, K=NULL,
                                       prob.matrix=locusmatrix,
                                       num.draws=10, 
                                       use.dip.lincomb=TRUE, seed=1, gamma.rate=1, impute.on="CCline")
inla.diploffect.summary <- run.diploffect.inla.summary.stats(inla.diploffect)
}
