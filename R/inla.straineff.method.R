################## General straineff method
StrainEffMethod <- setRefClass("StrainEffMethod",
                               fields = list(X_="matrix", Y_="matrix", scale_="matrix", model_="character",
                                             data_="matrix", N_="numeric",
                                             Z_="matrix",
                                             Z2_="list",
                                             K_="matrix", cholK_="matrix", rawZ_="factor", rawZ2_="data.frame"),
                               methods = list(
                                 init = function(model, data, X, Y, scale, Z, Z2, K){
                                   'TODO: call by data.additive[1:N,1:8]'
                                   Y_ <<- Y
                                   if(!is.null(X)) X_ <<- X
                                   data_ <<- data
                                   model_ <<- model

                                   N_ <<- dim(Y)[2]

                                   if(!is.null(Z)){
                                     Z_ <<- incidence.matrix(as.factor(Z))
                                     rawZ_ <<- Z
                                   }
                                   if(!is.null(Z2)){
                                     Z2_ <<- list()
                                     for(i in 1:ncol(Z2)){
                                       Z2_[[i]] <<- incidence.matrix(as.factor(Z2[,i]))
                                     }
                                     names(Z2_) <<- names(Z2)
                                     rawZ2_ <<- Z2
                                   }
                                   if(!is.null(K)){
                                     K_ <<- K
                                   }
                                 },
                                 can.estimate = function(model){
                                   'check whether the class can be used for estimation'
                                   stop("Not Implemented")
                                 },
                                 estimate = function(){
                                   'Estimate strain effects'
                                   stop("Not Implemented")
                                 },
                                 get.straineff = function(M, result){
                                   'Estimate strain effects'
                                   stop("Not Implemented")
                                 },
                                 is.straineff = function(){
                                   'Is strain effects model?'
                                   return(FALSE)
                                 },
                                 is.straineff.mcmc = function(){
                                   'Is strain effects model with MCMC sampling?'
                                   return(FALSE)
                                 },
                                 has.deviation.effects = function(){
                                   'Does the method supply deviation estimation?'
                                   return (FALSE)
                                 },
                                 get.deviation.effects = function(M, result){
                                   'Estimate strain effects'
                                   stop("Not Implemented")
                                 },
                                 has.deviation.effects = function(){
                                   return(FALSE)
                                 },
                                 get.diplotype.effects = function(M, result){
                                   beta = get.straineff(M, result)
                                   if (!has.deviation.effects()){
                                     deviation.effects <- mat.or.vec(1, M * (M + 1) / 2)
                                   }
                                   else{
                                     deviation.effects <- get.deviation.effects(M, result)
                                   }
                                   return(calculate.diplotype.effects(beta, deviation.effects))
                                 },
                                 get.haplotype.score = function(result){
                                   return(0)
                                 },
                                 has.haplotype.score = function(result){
                                   return(F)
                                 },
                                 get.display.name = function(result){
                                   return(model_)
                                 },
                                 get.extra = function(result){
                                   return(NULL)
                                 }
                               )
)

get.straineff.method <- function(model){
  candidate.methods = c(INLAMethod)
  for (m in candidate.methods){
    ret = m$new()
    if(ret$can.estimate(model)){
      return(ret)
    }
  }
  return(NULL)
}

straineff.method.factory <- function(model, data, data.additive, data.deviated,
                                     X, Y, scale, R, R2, K=NULL){
  method = get.straineff.method(model)
  if (!is.null(method)){
    method$init(model, data, data.additive, data.deviated, X, Y, scale, R, R2, K)
    return(method)
  }
  else{
    stop("Did you forget to put your method file in straineff.method.R?")
  }
}

if(exists(".inlaEnv") && is.environment(.inlaEnv)){
  ## then reuse it
} else{
  .inlaEnv = new.env()
}
requireNamespace("INLA")

################## INLA method
if(!isTRUE(requireNamespace("INLA", quietly = TRUE))){
  stop("You need to install the packages 'INLA'. Please run in your R terminal:\n install.packages('INLA', repos='https://www.math.ntnu.no/inla/R/stable')")
}
# If INLA is installed, then attach the Namespace (so that all the relevant functions are available)
if(isTRUE(requireNamespace("INLA", quietly = TRUE))){
  if(!is.element("INLA", (.packages()))){
    attachNamespace("INLA")
  }

  ## Annoying issue with INLA and its Namespace
INLAMethod <- setRefClass("INLAMethod",
                          contains = "StrainEffMethod",
                          methods = list(
                            can.estimate = function(model){
                              model_ <<- model
                              return(grepl('inla', model))
                            },
                            use.weight = function(){
                              return(!grepl('noweight', model_))
                            },
                            use.kinship = function(){
                              return(grepl('kinship', model_))
                            },
                            estimate = function(num.draws=200,
                                                family="gaussian",
                                                num.threads=4,
                                                use.dip.lincomb=FALSE,
                                                this.seed=NULL,
                                                founder.names,
                                                gamma.rate=1,
                                                latent.sample.num=1000,
                                                return.joint.posterior.samples=FALSE,
                                                add.only=FALSE,
                                                var.components,
                                                imputation.map){
                              S = num.draws

                              M <- length(founder.names)
                              ## The number of phenotypes
                              results = list(logmliks=rep(0, S), betas = matrix(0, S, M),
                                             summary.random=list(), inla.objects=list(), hyper.samples=list())
                              T = M * (M + 1) / 2
                              diplotypes = list()
                              case = 1
                              strain.names <- colnames(data_)

                              if(!is.null(this.seed)){
                                set.seed(this.seed)
                              }
                              if(use.kinship()){
                                eigen.K <- eigen(K_)
                              }
                              while(TRUE){
                                s = dim(data_)[2]
                                data = matrix(0, N_, s)

                                diplotype = rep(0, N_)
                                mapping = straineff.mapping.matrix.happy(M)
                                ## Impute the diplotype (haplotype pair) status
                                x = matrix(0, N_, M)

                                data <- run.imputation(diplotype.probs=data_, imputation.map=imputation.map)

                                diplotype <- apply(data, 1, function(x) which(x == 1))
                                x <- t(mapping[,diplotype])
                                X.add <- x
                                X.dom <- data[,-(1:M)]
                                X.full <- data

                                diplotypes[[case]] = diplotype
                                ## If the diplotype does not exist in the imputed data,
                                ## remove the corresponding diplotype effects.
                                selected = c()
                                for (i in (M + 1):T){
                                  if (sum(data[, i]) != 0){
                                    selected = c(selected, i)
                                  }
                                }
                                dom.M <- length(selected)
                                ## the following script is based on the example at the bottom
                                ## of http://www.r-inla.org/models/tools
                                if(dim(X_)[1]==0){
                                  x = cbind(rep(1, N_), x)
                                }
                                else{
                                  x = cbind(rep(1, N_), X_, x)
                                }
                                if(!add.only){
                                  x = cbind(x, data[,(M+1):T])
                                }
                                ## Kinship effects
                                poly.M <- 0
                                if(dim(Z_)[1] != 0){
                                  if(use.kinship()){
                                    ck = K_
                                    x = cbind(x, diag(ncol(K_)))
                                    poly.M = dim(K_)[1]
                                  }
                                  else{
                                    x = cbind(x, Z_)
                                    poly.M = dim(Z_)[2]
                                  }
                                }
                                ## Sparse random effects
                                random.M <- 0
                                if(length(Z2_)[1] > 0){
                                  for(i in 1:length(Z2_)){
                                    x = cbind(x, Z2_[[i]])
                                    random.M = random.M + dim(Z2_[[i]])[2]
                                  }
                                }
                                ## Fixed effects
                                fixed.M <- dim(X_)[2]
                                fixed <- NA

                                dev.M <- T - M
                                smart.index.return <- function(x){
                                  if(x == 0) { return(NULL) }
                                  if(x != 0) { return(1:x) }
                                }
                                if(add.only){
                                  y = Y_[1, ]
                                  data <- list(y=y,
                                               intercept=c(1, rep(NA, fixed.M + M + poly.M + random.M)),
                                               #fixed=c(NA, smart.index.return(fixed.M), rep(NA,  M + poly.M + random.M)),
                                               idx=c(rep(NA, 1 + fixed.M), smart.index.return(M), rep(NA, poly.M + random.M)),
                                               poly.idx=c(rep(NA, 1 + fixed.M + M), smart.index.return(poly.M), rep(NA, random.M)))
                                               #random.idx=c(rep(NA, 1 + fixed.M + M + poly.M), smart.index.return(random.M))
                                  if(dim(X_)[2] != 0){
                                    for(i in 1:ncol(X_)){
                                      this.fix.effect <- rep(NA, ncol(X_))
                                      this.fix.effect[i] <- 1
                                      data[[paste0("fixed.", gsub(pattern=" ", replacement=".", x=colnames(X_), fixed=TRUE)[i])]] <- c(NA, this.fix.effect, rep(NA, M + poly.M + random.M))
                                    }
                                  }
                                  Z2 <- Z2_
                                  if(length(Z2_) > 0){
                                    start.index <- 1
                                    for(i in 1:length(Z2_)){
                                      this.random.effect <- rep(NA, random.M)
                                      this.random.effect[start.index:(start.index + ncol(Z2_[[i]]) - 1)] <- 1:ncol(Z2_[[i]])
                                      data[[paste0("random.", gsub(pattern=" ", replacement=".", x=names(rawZ2_), fixed=TRUE)[i])]] <- c(rep(NA, 1 + fixed.M + M + poly.M), this.random.effect)
                                      start.index <- ncol(Z2_[[i]]) + 1
                                    }
                                    #Z2 <- Z2_
                                    names(Z2) <- paste0("random.", names(Z2_))
                                  }
                                }
                                if(!add.only){
                                  y = Y_[1, ]
                                  data <- list(y=y,
                                               intercept=c(1, rep(NA, fixed.M + M + dev.M + poly.M + random.M)),
                                               #fixed=c(NA, smart.index.return(fixed.M), rep(NA,  M + dev.M + poly.M + random.M)),
                                               idx=c(rep(NA, 1 + fixed.M), smart.index.return(M), rep(NA, dev.M + poly.M + random.M)),
                                               dom.idx=c(rep(NA, 1 + fixed.M + M), smart.index.return(dev.M), rep(NA, poly.M + random.M)),
                                               poly.idx=c(rep(NA, 1 + fixed.M + M + dev.M), smart.index.return(poly.M), rep(NA, random.M)))
                                               #random.idx=c(rep(NA, 1 + fixed.M + M + dev.M + poly.M), smart.index.return(random.M)))
                                  if(dim(X_)[2] != 0){
                                    for(i in 1:ncol(X_)){
                                      this.fix.effect <- rep(NA, ncol(X_))
                                      this.fix.effect[i] <- 1
                                      data[[paste0("fixed.", gsub(pattern=" ", replacement=".", x=colnames(X_), fixed=TRUE)[i])]] <- c(NA, this.fix.effect, rep(NA, M + dev.M + poly.M + random.M))
                                    }
                                  }
                                  Z2 <- Z2_
                                  if(length(Z2_) > 0){
                                    start.index <- 1
                                    for(i in 1:length(Z2_)){
                                      this.random.effect <- rep(NA, random.M)
                                      this.random.effect[start.index:(start.index + ncol(Z2_[[i]]) - 1)] <- 1:ncol(Z2_[[i]])
                                      data[[paste0("random.", gsub(pattern=" ", replacement=".", x=names(rawZ2_), fixed=TRUE)[i])]] <- c(rep(NA, 1 + fixed.M + M + dev.M + poly.M), this.random.effect)
                                      start.index <- ncol(Z2_[[i]]) + 1
                                    }
                                    #Z2 <- Z2_
                                    names(Z2) <- paste0("random.", names(Z2_))
                                  }
                                }
                                if(use.kinship()){
                                  data$inv.K = eigen.K$vectors %*% ((1/eigen.K$values) * t(eigen.K$vectors))
                                }
                                formula = NA
                                ## http://stackoverflow.com/questions/4951442/formula-with-dynamic-number-of-variables
                                pre = "y ~ -1"
                                intercept = "intercept"
                                if (fixed.M != 0) {
                                  ##fixed = "f(fixed, hyper=list(theta=list(prior=\"loggamma\", param=c(1,0.01))))"
                                  fixed = paste0("fixed.", gsub(pattern=" ", replacement=".", x=colnames(X_), fixed=TRUE))
                                }
                                dom.effect <- NULL
                                poly.effect <- NULL
                                random.effect <- NULL

                                effect = paste("f(idx, constr=TRUE, hyper=list(theta=list(prior=\"loggamma\", param=c(1,", gamma.rate, "))))", sep="")
                                if(!add.only){
                                  dom.effect = paste("f(dom.idx, constr=TRUE, hyper=list(theta=list(prior=\"loggamma\", param=c(1,", gamma.rate, "))))", sep="")
                                }
                                if(poly.M != 0){
                                  poly.effect = paste("f(poly.idx, constr=TRUE, model='generic0', Cmatrix=inv.K, hyper=list(theta=list(prior=\"loggamma\", param=c(1,", gamma.rate, "))))", sep="")
                                }
                                if(random.M != 0){
                                  random.effect <- paste(paste0("f(",
                                                                paste0("random.", gsub(pattern=" ", replacement=".", x=names(rawZ2_), fixed=TRUE)),
                                                                ", constr=TRUE, hyper=list(theta=list(prior=\"loggamma\", param=c(1,", gamma.rate, "))))"),
                                                         collapse=" + ")
                                }
                                all = c(pre, intercept, fixed, effect, dom.effect, poly.effect, random.effect)
                                all = all[!is.na(all)]
                                formula <- as.formula(paste(all, collapse = "+"))
                                this.lincomb <- NULL
                                if (!add.only & use.dip.lincomb){
                                  this.lincomb <- rep(NA, T)
                                  for (i in 1:M){
                                    this.homozygote <- rep(NA, M)
                                    this.homozygote[i] <- 2
                                    lc = INLA::inla.make.lincomb(idx=this.homozygote)
                                    this.lincomb[i] <- lc
                                    names(this.lincomb)[i] <- paste("lc", i, sep="")
                                  }
                                  for (i in (M+1):T){
                                    this.heterozygote.add <- rep(NA, M)
                                    this.heterozygote.add[which(mapping[,i] == 1)] <- 1
                                    this.heterozygote.dom <- rep(NA, dev.M)
                                    this.heterozygote.dom[i-8] <- 1
                                    lc = INLA::inla.make.lincomb(idx=this.heterozygote.add, dom.idx=this.heterozygote.dom)
                                    this.lincomb[i] <- lc
                                    names(this.lincomb)[i] <- paste("lc", i, sep="")
                                  }
                                }
                                if(family=="gaussian"){
                                  result <- INLA::inla(formula=formula,
                                                       data=data,
                                                       family=family,
                                                       control.family=list(hyper = list(prec = list(prior = "loggamma", param = c(1, gamma.rate)))),
                                                       scale=scale,
                                                       lincomb=this.lincomb,
                                                       control.predictor=list(A=x, compute=TRUE),
                                                       control.compute=list(config=TRUE),
                                                       quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975),
                                                       num.threads=num.threads,
                                                       verbose=FALSE)
                                }
                                if(family=="T"){
                                  result <- INLA::inla(formula=formula,
                                                       data=data,
                                                       family=family,
                                                       control.family=list(hyper = list(prec = list(prior="loggamma", param = c(1, gamma.rate)),
                                                                                                    dof = list(prior = "loggamma", param=c(30, 5)))),
                                                       lincomb=this.lincomb,
                                                       control.predictor=list(A=x, compute=TRUE),
                                                       quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975),
                                                       num.threads=num.threads, verbose=FALSE)
                                }
                                if(length(result$mlik) == 0){
                                  print("failed, reattempt..")
                                  next;
                                }
                                results$logmliks[case] = result$mlik[1]
                                cat(paste(paste0(case, ":"), result$mlik[1]), "\n")

                                results$betas[case, 1:M] = result$summary.random[["idx"]]$mean[1:M]
                                result$summary.random$lincomb <- result$summary.lincomb.derived
                                results$summary.random[[case]] = result$summary.random
                                results$inla.objects[[case]] = list(marginals.random = result$marginals.random, marginals.lincomb = result$marginals.lincomb.derived, summary.hyperpar=result$summary.hyperpar)
                                hyper.samples <- INLA::inla.hyperpar.sample(n=1000, result=result, intern=FALSE)
                                latent.samples <- INLA::inla.posterior.sample(n=latent.sample.num, result=result, use.improved.mean=TRUE)

                                full.par <- sapply(latent.samples, function(x) x$latent)
                                rownames(full.par) <- rownames(latent.samples[[1]]$latent)
                                results$joint.posterior.samples <- NULL
                                if(return.joint.posterior.samples){
                                  results$joint.posterior.samples[[case]] <- full.par
                                }

                                ## Calculating sums of squares
                                SS.mat <- t(apply(full.par, 2, function(x) get.partial.ss(x, y=y, X.fix=X_, X.add=X.add, X.dom=X.dom, Z2=Z2, var.components=var.components)))

                                hyper.samples.exp <- matrix(NA, nrow=nrow(hyper.samples), ncol=length(var.components))
                                colnames(hyper.samples.exp) <- var.components
                                for(i in 1:length(var.components)){
                                  hyper.samples.exp[,i] <- apply(hyper.samples, 1, function(x) (1/exp(x[grepl(colnames(hyper.samples), pattern=paste0(" ", var.components[i]), fixed=TRUE)]))/sum(1/exp(x)))
                                }
                                colnames(hyper.samples.exp)[colnames(hyper.samples.exp) == "idx"] <- "QTL.add.exp"
                                colnames(hyper.samples.exp)[colnames(hyper.samples.exp) == "dom.idx"] <- "QTL.dom.exp"
                                colnames(hyper.samples.exp)[colnames(hyper.samples.exp) == "poly.idx"] <- "QTL.kinship.exp"

                                results$hyper.samples[[case]] = hyper.samples
                                results$hyper.varexp.samples[[case]] = hyper.samples.exp
                                results$SS.samples[[case]] = SS.mat
                                #diplotypes
                                if(case == S){ break; }
                                case <- case + 1;
                              }
                              list(results=results, diplotypes=diplotypes, strain.names=strain.names)
                            },
                            is.straineff = function(){
                              "Is strain effects model?"
                              return(T)
                            },
                            get.sim.straineff = function(M, result){
                              ret = get.straineff(M, result$strain.effects.result)
                              return(ret)
                            },
                            get.samples.mean = function(samples, fieldname, rowidx, colname, mliks){
                              S = length(samples)
                              stopifnot(S == length(mliks))
                              s = rep(0, S)

                              num.effect = length(samples[[i]][[fieldname]][, 'mean'])
                              for (i in 1:S) {
                                shift = sum(samples[[i]][[fieldname]][, 'mean']) / num.effect
                                s[i] = samples[[i]][[fieldname]][rowidx, colname] - shift
                              }
                              if(use.weight()){
                                return(weighted.mean(s, mliks))
                              }
                              else{
                                return(mean(result$betas[, i]))
                              }
                            },
                            get.ci.from.inla.marginal = function(marginals){
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
                            },
                            get.straineff.ci = function(M, result) {
                              logmliks = as.vector(scale(result$results$logmliks, scale=F))
                              mliks = exp(logmliks)
                              marginal = list()
                              for (i in 1:M) {
                                marginal[[i]] =
                                  get.combined.marginal(result$results$inla.objects, "idx",
                                                        i, mliks)
                              }
                              return(get.ci.from.inla.marginal(marginal))
                            },
                            get.deviated.effects = function(M, result) {
                              get.deviated.effects.mean(M, result)
                            },
                            get.deviated.effects.mean = function(M, result) {
                              logmliks = as.vector(scale(result$results$logmliks, scale=F))
                              mliks = exp(logmliks)
                              rowSums(sapply(1:length(mliks), function(i) (result$results$summary.random[[i]][["dom.idx"]]$mean * mliks[i]))) / sum(mliks)
                            },
                            get.deviated.effects.ci = function(M, result){
                              logmliks = as.vector(scale(result$results$logmliks, scale=F))
                              mliks = exp(logmliks)
                              marginal = list()
                              for (i in 1:choose(M, 2)) {
                                marginal[[i]] =
                                  get.combined.marginal(result$results$inla.objects, "dom.idx",
                                                        i, mliks)
                              }
                              return (get.ci.from.inla.marginal(marginal))
                            },
                            get.diplotypes.ci = function(M, result) {
                              logmliks = as.vector(scale(result$results$logmliks, scale=F))
                              mliks = exp(logmliks)
                              marginal = list()
                              for (i in 1:(M + choose(M, 2))) {
                                marginal[[i]] =
                                  get.combined.marginal(result$results$inla.objects, paste('lc', i, sep=""), i,
                                                        mliks, is.lincomb=TRUE)
                              }
                              return (get.ci.from.inla.marginal(marginal))
                            },
                            get.straineff = function(M, result) {
                              logmliks = as.vector(scale(result$results$logmliks, scale=F))
                              if (use.weight()) {
                                mliks = exp(logmliks)
                                rowSums(sapply(1:length(mliks), function(i) (result$results$summary.random[[i]]$idx$mean * mliks[i]))) / sum(mliks)
                              }
                              else {
                                mliks = rep(1, length(logmliks))
                                rowSums(sapply(1:length(mliks), function(i) (result$results$summary.random[[i]]$idx$mean * mliks[i]))) / sum(mliks)
                              }
                            },
                            get.qtl.h2 = function(result, effect="add") {
                              logmliks = as.vector(scale(result$results$logmliks, scale=F))
                              if (use.weight()) {
                                mliks = exp(logmliks)
                              }
                              else {
                                mliks = rep(1, length(logmliks))
                              }
                              make.density.matrix <- function(density.object) {
                                density.matrix <- cbind(density.object$x, density.object$y)
                                colnames(density.matrix) <- c("x", "y")
                                return(density.matrix)
                              }
                              density.list <- lapply(1:length(mliks), function(i) make.density.matrix(density(result$results$hyper.varexp.samples[[i]][,paste("QTL.", effect, ".exp", sep="")], from=0, to=1)))
                              marginal = list(marginal=get.combined.marginal.emp(density.list, mliks))
                              return (get.ci.from.inla.marginal(marginal))
                            },
                            get.total.qtl.h2 = function(result) {
                              logmliks = as.vector(scale(result$results$logmliks, scale=F))
                              if (use.weight()) {
                                mliks = exp(logmliks)
                              }
                              else {
                                mliks = rep(1, length(logmliks))
                              }
                              make.density.matrix <- function(density.object) {
                                density.matrix <- cbind(density.object$x, density.object$y)
                                colnames(density.matrix) <- c("x", "y")
                                return(density.matrix)
                              }
                              h2.total.list <- lapply(1:length(mliks), function(i) rowSums(result$results$hyper.varexp.samples[[i]][, c("QTL.add.exp", "QTL.dom.exp")]))
                              density.list <- lapply(1:length(mliks), function(i) make.density.matrix(density(h2.total.list[[i]], from=0, to=1)))
                              marginal = list(marginal=get.combined.marginal.emp(density.list, mliks))
                              return (get.ci.from.inla.marginal(marginal))
                            },
                            get.nonqtl.h2 = function(result, variable) {
                              logmliks = as.vector(scale(result$results$logmliks, scale=F))
                              if (use.weight()) {
                                mliks = exp(logmliks)
                              }
                              else {
                                mliks = rep(1, length(logmliks))
                              }
                              make.density.matrix <- function(density.object) {
                                density.matrix <- cbind(density.object$x, density.object$y)
                                colnames(density.matrix) <- c("x", "y")
                                return(density.matrix)
                              }
                              density.list <- lapply(1:length(mliks), function(i) make.density.matrix(density(result$results$hyper.varexp.samples[[i]][,variable], from=0, to=1)))
                              marginal = list(marginal=get.combined.marginal.emp(density.list, mliks))
                              return (get.ci.from.inla.marginal(marginal))
                            },
                            get.nonqtl.ss.h2 = function(result, variable) {
                              logmliks = as.vector(scale(result$results$logmliks, scale=F))
                              if (use.weight()) {
                                mliks = exp(logmliks)
                              }
                              else {
                                mliks = rep(1, length(logmliks))
                              }
                              make.density.matrix <- function(density.object) {
                                density.matrix <- cbind(density.object$x, density.object$y)
                                colnames(density.matrix) <- c("x", "y")
                                return(density.matrix)
                              }
                              density.list <- lapply(1:length(mliks), function(i) make.density.matrix(density(result$results$SS.samples[[i]][,variable], from=0, to=1)))
                              marginal = list(marginal=get.combined.marginal.emp(density.list, mliks))
                              return (get.ci.from.inla.marginal(marginal))
                            },
                            get.SS.h2 = function(result, effect="add") {
                             logmliks = as.vector(scale(result$results$logmliks, scale=F))
                             if (use.weight()) {
                               mliks = exp(logmliks)
                             }
                             else {
                               mliks = rep(1, length(logmliks))
                             }
                             make.density.matrix <- function(density.object) {
                               density.matrix <- cbind(density.object$x, density.object$y)
                               colnames(density.matrix) <- c("x", "y")
                               return(density.matrix)
                             }
                             density.list <- lapply(1:length(mliks), function(i) make.density.matrix(density(result$results$SS.samples[[i]][,paste("S.", effect, sep="")], from=0, to=1)))
                             marginal = list(marginal=get.combined.marginal.emp(density.list, mliks))
                             return (get.ci.from.inla.marginal(marginal))
                            }
                          )
)

run.imputation <- function(diplotype.probs, imputation.map){
  diplotype.probs <- data.frame(original.order=1:nrow(diplotype.probs), SUBJECT.NAME=rownames(diplotype.probs), diplotype.probs, stringsAsFactors=FALSE)
  diplotype.probs <- merge(x=diplotype.probs, y=imputation.map, by="SUBJECT.NAME")
  diplotype.probs <- diplotype.probs[order(diplotype.probs$original.order),]
  diplotype.probs <- diplotype.probs[, names(diplotype.probs) != "original.order"]

  imputable.diplotype.probs <- diplotype.probs[as.integer(rownames(unique(data.frame(diplotype.probs[,"impute.on"])))),]
  rownames(imputable.diplotype.probs) <- imputable.diplotype.probs[, "impute.on"]
  imputable.diplotype.probs <- imputable.diplotype.probs[,!(names(imputable.diplotype.probs) %in% names(imputation.map))]
  imputation <- t(apply(imputable.diplotype.probs, 1, function(x) rmultinom(1, 1, x)))
  full.imputation <- imputation[as.character(imputation.map[, "impute.on"]),]
  rownames(full.imputation) <- imputation.map[, "SUBJECT.NAME"]
  return(full.imputation)
}

get.combined.marginal = function(samples, fieldname, rowidx, mliks, is.lincomb=FALSE) {
  S = length(samples)
  stopifnot(S == length(mliks))
  s = rep(0, S)
  get.max.from.marginal <- function(s, is.lincomb){
    if(is.lincomb){
      return(max(s$marginals.lincomb[[fieldname]][, 1]))
    }
    else if(is.null(fieldname) & is.null(rowidx)){
      return(max(s[, 1]))
    }
    else{
      return(max(s$marginals.random[[fieldname]][[rowidx]][, 1]))
    }
  }
  get.min.from.marginal <- function(s, is.lincomb){
    if(is.lincomb) {
      return(min(s$marginals.lincomb[[fieldname]][, 1]))
    }
    else if(is.null(fieldname) & is.null(rowidx)){
      return(min(s[, 1]))
    }
    else{
      return(min(s$marginals.random[[fieldname]][[rowidx]][, 1]))
    }
  }
  left = min(unlist(lapply(samples, function(x) get.max.from.marginal(x, is.lincomb=is.lincomb))))
  right = max(unlist(lapply(samples, function(x) get.min.from.marginal(x, is.lincomb=is.lincomb))))

  x = seq(left, right, length.out=200)
  get.weighted.density <- function(x, is.lincomb) {
    s = rep(0, S)
    for (i in 1:S) {
      if (is.lincomb) {
        s[i] = INLA::inla.dmarginal(x, samples[[i]]$marginals.lincomb[[fieldname]])
      }
      else {
        s[i] = INLA::inla.dmarginal(x, samples[[i]]$marginals.random[[fieldname]][[rowidx]])
      }
    }
    return(weighted.mean(s, mliks))
  }
  y = unlist(lapply(x, function(x) get.weighted.density(x, is.lincomb=is.lincomb)))
  return (data.frame(x, y))
}

get.combined.marginal.emp = function(samples, mliks) {
  S = length(samples)
  stopifnot(S == length(mliks))
  s = rep(0, S)
  get.max.from.marginal <- function(s) {
    return(max(s[, 1]))
  }
  get.min.from.marginal <- function(s) {
    return(min(s[, 1]))
  }
  left = min(unlist(lapply(samples, function(x) get.max.from.marginal(x))))
  right = max(unlist(lapply(samples, function(x) get.min.from.marginal(x))))

  x = seq(left, right, length.out=200)
  get.weighted.density <- function(x) {
    s = rep(0, S)
    for (i in 1:S) {
      s[i] = inla.dmarginal(x, samples[[i]])
    }
    return (weighted.mean(s, mliks))
  }
  y = unlist(lapply(x, function(x) get.weighted.density(x)))
  return (data.frame(x, y))
}

get.partial.ss <- function(par, y, X.fix=NULL, X.add=NULL, X.dom=NULL, Z2=NULL, var.components) {
  intercept.idx <- grepl(pattern="intercept", names(par))
  fix.par.idx <- grepl(pattern="fixed", names(par))
  add.par.idx <- grepl(pattern="^idx", names(par))
  dom.par.idx <- grepl(pattern="^dom", names(par))
  poly.par.idx <- grepl(pattern="poly", names(par))
  # Additional variance components
  extra.idx.var <- var.components[!(var.components %in% c("idx", "dom.idx", "poly.idx"))]

  SS.results <- SS.names <- NULL
  SS.0 <- sum((y - rep(1, length(y)) * par[intercept.idx])^2)

  #### Setting up fixed effects matrix
  fix.dim.from.inla <- dim(X.fix)
  if(fix.dim.from.inla[2] == 0){
    X.fix <- matrix(rep(1, length(y)), ncol=1)
  }
  if(fix.dim.from.inla[2] > 0){
    X.fix <-  cbind(rep(1, length(y)), X.fix)
  }
  fix.par <- c(par[intercept.idx], par[fix.par.idx])
  this.fit <- X.fix %*% fix.par
  SS.fix <- sum((y - this.fit)^2)
  recent.SS <- SS.fix
  if(any(add.par.idx)){
    this.fit <- this.fit + X.add %*% par[add.par.idx]
    SS.add <- sum((y - this.fit)^2)
    SS.names <- c(SS.names, "S.add")
    SS.results <- c(SS.results, (recent.SS - SS.add)/SS.fix)
    recent.SS <- SS.add
  }
  if(any(dom.par.idx)){
    this.fit <- this.fit + X.dom %*% par[dom.par.idx]
    SS.dom <- sum((y - this.fit)^2)
    SS.names <- c(SS.names, "S.dom")
    SS.results <- c(SS.results, (recent.SS - SS.dom)/SS.fix)
    recent.SS <- SS.dom
  }
  if(any(dom.par.idx) & any(add.par.idx)){
    SS.names <- c(SS.names, "S.qtl")
    SS.results <- c(SS.results, (SS.fix - SS.dom)/SS.fix)
  }
  if(any(poly.par.idx)){
    this.fit <- this.fit + diag(length(y)) %*% par[poly.par.idx]
    SS.poly <- sum((y - this.fit)^2)
    SS.results <- c(SS.results, (recent.SS - SS.poly)/SS.fix)
    SS.names <- c(SS.names, "S.poly")
    recent.SS <- SS.poly
  }
  if(length(extra.idx.var) > 0){
    for(i in 1:length(extra.idx.var)){
      extra.idx <- grepl(pattern=extra.idx.var[i], names(par))
      this.fit <- this.fit + Z2[[extra.idx.var[i]]] %*% par[extra.idx]
      SS.random <- sum((y - this.fit)^2)
      SS.names <- c(SS.names, extra.idx.var[i])
      SS.results <- c(SS.results, (recent.SS - SS.random)/SS.fix)
      recent.SS <- SS.random
    }
  }

  if(is.vector(SS.results)){
    SS.results <- matrix(SS.results, ncol=1)
  }
  names(SS.results) <- SS.names
  return(SS.results)
}

load.deviated.mean.from.inla <- function(M, diploffect.inla, scale=1) {
  method = get.straineff.method("inla")
  ci <- method$get.deviated.effects.mean(M, diploffect.inla)
  ci <- sapply(ci, function(t) {t * scale})
  return(ci)
}

load.ci.from.inla <- function(M, diploffect.inla, scale=1) {
  method = get.straineff.method("inla")
  ci <- method$get.straineff.ci(M, diploffect.inla)
  ci <- lapply(ci, function(t) {t * scale})
  return(ci)
}

load.deviated.ci.from.inla <- function(M, diploffect.inla, scale=1) {
  method = get.straineff.method("inla")
  ci <- method$get.deviated.effects.ci(M, diploffect.inla)
  ci <- lapply(ci, function(t) {t * scale})
  return(ci)
}

load.diplotypes.ci.from.inla <- function(M, diploffect.inla, scale=1) {
  method = get.straineff.method("inla")
  ci <- method$get.diplotypes.ci(M, diploffect.inla)
  ci <- lapply(ci, function(t) {t * scale})
  return(ci)
}

load.nonqtl.h2.from.inla <- function(diploffect.inla, variable, scale=1){
  method = get.straineff.method("inla")
  h2.stats <- method$get.nonqtl.h2(diploffect.inla, variable=variable)
  h2.stats <- lapply(h2.stats, function(t){ t * scale })
  return(h2.stats)
}

load.nonqtl.ss.h2.from.inla <- function(diploffect.inla, variable, scale=1){
  method = get.straineff.method("inla")
  h2.stats <- method$get.nonqtl.ss.h2(diploffect.inla, variable=variable)
  h2.stats <- lapply(h2.stats, function(t){ t * scale })
  return(h2.stats)
}

load.qtl.h2.from.inla <- function(diploffect.inla, effect="add", scale=1) {
  method = get.straineff.method("inla")
  h2.stats <- method$get.qtl.h2(diploffect.inla, effect=effect)
  h2.stats <- lapply(h2.stats, function(t) {t * scale})
  return(h2.stats)
}

load.total.qtl.h2.from.inla <- function(diploffect.inla, scale=1) {
  method = get.straineff.method("inla")
  h2.stats <- method$get.total.qtl.h2(diploffect.inla)
  h2.stats <- lapply(h2.stats, function(t) {t * scale})
  return(h2.stats)
}

load.SS.h2.from.inla <- function(diploffect.inla, effect="add", scale=1) {
  method = get.straineff.method("inla")
  SS.h2.stats <- method$get.SS.h2(diploffect.inla, effect=effect)
  SS.h2.stats <- lapply(SS.h2.stats, function(t) {t * scale})
  return(SS.h2.stats)
}


}

