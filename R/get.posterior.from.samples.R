
#' @export
external.get.ci.from.inla.marginal <- function(marginals){
  effect.num = length(marginals)
  med = rep(0, effect.num)
  mu = rep(0, effect.num)
  quant.wide.left = rep(0, effect.num)
  quant.wide.right = rep(0, effect.num)
  quant.narrow.left = rep(0, effect.num)
  quant.narrow.right = rep(0, effect.num)
  identity <- function(x) {x}
  for (i in 1:effect.num) {
    med[i] = INLA::inla.qmarginal(0.5, marginals[[i]])
    mu[i] = INLA::inla.emarginal(identity, marginals[[i]])
    quant.wide.left[i] = INLA::inla.qmarginal(0.025, marginals[[i]])
    quant.wide.right[i] = INLA::inla.qmarginal(0.975, marginals[[i]])
    quant.narrow.left[i] = INLA::inla.qmarginal(0.25, marginals[[i]])
    quant.narrow.right[i] = INLA::inla.qmarginal(0.75, marginals[[i]])
  }
  return(list(med=med, mu=mu,
              quant.narrow=cbind(quant.narrow.left, quant.narrow.right),
              quant.wide=cbind(quant.wide.left, quant.wide.right)))
}

#' @export
get.posterior.from.samples <- function(samples.matrix.list, logmlik, variable){
  logmliks = as.vector(scale(result$results$logmliks, scale=F))
  mlik = exp(logmlik)
  
  make.density.matrix <- function(density.object) {
    density.matrix <- cbind(density.object$x, density.object$y)
    colnames(density.matrix) <- c("x", "y")
    return(density.matrix)
  }
  
  density.list <- lapply(1:length(mlik), function(i) make.density.matrix(density(samples.matrix.list[[i]][,variable], from=0, to=1)))
  marginal <- list(marginal=get.combined.marginal.emp(density.list, mliks))
  return (get.ci.from.inla.marginal(marginal))
}