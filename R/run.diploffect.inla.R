#' Run Diploffect model on diplotype (haplotype pair) probabilities
#'
#' This function primarily takes a formula, data frame, and the diplotype probability matrix for a locus
#' and runs the corresponding Diploffect model for a given number of samples.
#'
#' @param formula An lm style formula with functions of outcome and covariates in data frame.
#' @param data A data frame with outcome and potential covariates. Should also have individual IDs
#' that link to IDs in the diplotype probability matrix, commonly with this column named "SUBJECT.NAME".
#' @param K A positive semi-definite relationship matrix, usually a realized genetic relationship matrix (GRM)
#' based on SNP genotypes are the haplotype probabilistic reconstructions. Colnames and rownames should match
#' the SUBJECT.NAME column in the data frame.
#' @param num.founders DEFAULT: 8. Number of haplotype alleles. For Collaborative Cross (CC) derived populations, this is 8.
#' @param prob.matrix An \eqn{n \times num.founders + {num.founders \choose 2}} matrix of diplotype probabilities. Rownames should match
#' IDs in K and data.
#' @param add.only DEFAULT: FALSE. Allows only additive haplotype effects to be included in the model.
#' @param num.draws The number of times the diplotype (haplotype pair) probabilities are sampled, and
#' INLA is run.
#' @param use.dip.lincomb DEFAULT: TRUE. Diploffect is parameterized based on additive haplotype effects
#' and heterozygous diplotype dominance effects. If posteriors on diplotypes are wanted, this option is required.
#' It does add to the computation time.
#' @param seed DEFAULT: 1. Diploffect involves a sampling process of the diplotypes, thus a seed is necessary
#' to produce the same results over multiple runs and different machines.
#' @param gamma.rate DEFAULT: 1. Hyper parameter for the variance component priors.
#' @param impute.on DEFAULT: "SUBJECT.NAME". Samples of diplotypes are drawn for each level of this variable.
#' Must be a column in the data frame. Allows the same diplotypes to be sampled if there are repeat observations
#' of a genome, such as in CC or CC-RIX.
#' @param weights.on DEFAULT: NULL. If numeric column in data frame is specified, weights are applied to the
#' diagonal of the random noise (diagonal) covariance matrix. This is an efficient approach to modeling
#' repeat observations, such as in CC or CC-RIX.
#' @param scale.response DEFAULT: TRUE. It is good to scale the response with INLA, thus the posteriors of effects
#' have a consistent scale (in standard deviations of the outcome).
#' @param founders DEFAULT: NULL. Provide alternative labels for the founder alleles. Defaults to the colnames of prob.matrix.
#' @param locus.name DEFAULT: NULL. Provide a locus label.
#' @return Diploffect.INLA object. Approximate posterior distributions for model parameters.
#' @export
#' @examples
#' library(Diploffect.INLA)
#' data(exampleCC)
#' data(locusmatrix)
#' inla.diploffect <- run.diploffect.inla(formula=y~1+(1|strain)+(1|dose.date), add.on=FALSE, data=exampleCC, K=NULL,
#'                                        prob.matrix=locusmatrix,
#'                                        num.draws=10, 
#'                                        use.dip.lincomb=TRUE, seed=1, gamma.rate=1, impute.on="CCline")
run.diploffect.inla <- function(formula, data, K=NULL,
                                num.founders=8, prob.matrix, add.only=FALSE,
                                num.draws,
                                use.dip.lincomb=TRUE, seed=1, gamma.rate=1,
                                impute.on="SUBJECT.NAME", weights.on=NULL, #"NUM.OBS",
                                scale.response=TRUE,
                                founders=NULL, locus.name=NULL, return.joint.posterior.samples=FALSE){

  weights.on <- weights.on[1]

  formula.string.new.data <- process.formula(formula, action="make.new.data.formula", impute.on=c(impute.on, weights.on))
  fixed.formula.string <- process.formula(formula, action="make.fixed.formula")
  final.formula.string <- process.formula(formula, action="make.final.formula")
  new.data <- model.frame(formula(formula.string.new.data), data=data)
  rownames(new.data) <- new.data$SUBJECT.NAME

  # Setting up K structured random effect
  if(is.null(K)) { Z <- NULL }
  if(!is.null(K)){
    # taking overlap of data with kinship
    overlap <- intersect(new.data$SUBJECT.NAME, rownames(K))
    new.data <- new.data[overlap,]
    K <- K[overlap, overlap]

    Z <- factor(1:ncol(K))
  }

  Y <- new.data[,1]
  Y <- matrix(Y, 1, length(Y))
  rownames(Y) <- names(new.data)[1]
  colnames(Y) <- new.data$SUBJECT.NAME
  if(scale.response){
    Y[1,] <- scale(as.vector(Y), center=TRUE, scale=TRUE)
  }
  else{
    cat("Warning: Not centering and scaling the phenotype is not suggested.\n Doing so helps standardize hyperparameters across phenotypes\n")
  }
  X <- model.matrix(formula(fixed.formula.string), data=new.data)[,-1, drop=FALSE]
  rownames(X) <- new.data$SUBJECT.NAME
  if(dim(X)[2] == 0){
    X <- NULL
  }

  # Setting up heterscedasticity
  if(is.null(weights.on)){
    scale <- matrix(rep(1, nrow(Y)), 1, ncol(Y))
  }
  if(!is.null(weights.on)){
    scale <- matrix(new.data[, weights.on], 1, ncol(Y))
  }
  # Setting up sparse random effects
  Z2 <- NULL
  if(grepl(pattern="(1", x=final.formula.string, fixed=TRUE)){ # make random effect matrix to set up R2
    Z2.frame <- model.frame(formula(formula.string.new.data), data=data)
    random.var <- process.formula(formula, action="return.random.variable")
    Z2 <- Z2.frame[,random.var, drop=FALSE]
    rownames(Z2) <- Z2.frame$SUBJECT.NAME
  }

  # Founders
  if(is.null(founders)){
    founders <- LETTERS[1:num.founders]
  }

  # Prepping imputation in Diploffect
  imputation.map <- new.data[, c("SUBJECT.NAME", impute.on), drop=FALSE]
  imputation.map <- data.frame(SUBJECT.NAME=new.data$SUBJECT.NAME, impute.on=new.data[,impute.on])

  model <- ifelse(is.null(K), "inla", "inla.kinship")

  # Variance components
  random.var <- process.formula(formula, action="return.random.variable")
  if(length(random.var) > 0){ random.var <- paste0("random.", random.var) }
  var.components <- c(random.var, c("idx", "dom.idx", "poly.idx")[c(TRUE, !add.only, !is.null(K))])
  inla = INLAMethod$new()
  inla$init(model=model, data=prob.matrix, X=X, Y=Y, scale=scale, Z=Z, Z2=Z2, K=K)
  result = inla$estimate(num.draws=num.draws,
                         family="gaussian",
                         num.threads=1,
                         use.dip.lincomb=use.dip.lincomb,
                         gamma.rate=gamma.rate,
                         this.seed=seed,
                         founder.names=founders,
                         add.only=add.only,
                         var.components=var.components,
                         return.joint.posterior.samples=return.joint.posterior.samples,
                         imputation.map=imputation.map)
  # Adding analysis information to Diploffect object
  genetic.effects <- c("additive", "dominant", "polygene")[c(TRUE, !add.only, !is.null(K))]
  result$analysis.id <- list(formula=final.formula.string, locus=locus.name,
                             num.draws=num.draws, M=length(founders), founders=founders,
                             genetic.effects=genetic.effects,
                             var.components=var.components)
  return(result)
}

process.formula <- function(formula,
                            action=c("add.subjects", "return.effect", "make.final", "make.random.data.formula", "make.fixed.formula"),
                            impute.on="SUBJECT.NAME"){
  formula.string <- paste0(Reduce(paste, deparse(formula)))

  lh.formula <- trimws(unlist(strsplit(x=formula.string, split="~", fixed=TRUE)))[1]
  rh.formula <- trimws(unlist(strsplit(x=formula.string, split="~", fixed=TRUE)))[2]
  rh.formula <- trimws(unlist(strsplit(x=rh.formula, split="+", fixed=TRUE)))

  locus.index <- grepl(pattern="locus.", x=rh.formula, fixed=TRUE)
  random.effect.index <- grepl(pattern="\\(1\\s*\\|", x=rh.formula, perl=TRUE)

  if(action == "return.QTL.effect"){
    QTL.effect.type <- gsub(pattern="locus.", replacement="", x=rh.formula[locus.index])
    return(QTL.effect.type)
  }
  if(action == "make.fixed.formula"){
    fixed.formula.string <- paste(lh.formula,
                                    paste(rh.formula[!(locus.index | random.effect.index)], collapse=" + "),
                                    sep=" ~ ")
    return(fixed.formula.string)
  }
  if(action == "make.final.formula"){
    final.formula.string <- paste(lh.formula,
                                  paste(rh.formula[!locus.index], collapse=" + "),
                                  sep=" ~ ")
    return(final.formula.string)
  }
  if(action == "make.new.data.formula"){
    random.var <- trimws(gsub(pattern="\\(1\\s*\\|", replacement="", x=rh.formula[random.effect.index]))
    random.var <- gsub(pattern=")", replacement="", x=random.var)
    formula.string.new.data <- paste(lh.formula,
                                  paste(paste(rh.formula[!(locus.index | random.effect.index)], collapse=" + "),
                                        paste(random.var, collapse=" + "),
                                        paste(unique(c("SUBJECT.NAME", impute.on)), collapse=" + "),
                                        sep=" + "),
                                  sep=" ~ ")
    return(formula.string.new.data)
  }
  if(action == "return.random.variable"){
    random.var <- trimws(gsub(pattern="\\(1\\s*\\|", replacement="", x=rh.formula[random.effect.index]))
    random.var <- gsub(pattern=")", replacement="", x=random.var)
    return(random.var)
  }
}

#' Estimate approximate posterior credible intervals for model parameters from Diploffect.INLA posterior distributions
#'
#' This function takes a Diploffect.INLA object and estimates the posterior credible intervals of model parameters
#'
#' @param diploffect.object A Diploffect.INLA summary object. Output from run.diploffect.inla.summary.stats(). Contains credible information
#' for all the model parameters.
#' @return Diploffect.INLA summary object.
#' @export
#' @examples 
#' library(Diploffect.INLA)
#' data(exampleCC)
#' data(locusmatrix)
#' inla.diploffect <- run.diploffect.inla(formula=y~1+(1|strain)+(1|dose.date), add.on=FALSE, data=exampleCC, K=NULL,
#'                                        prob.matrix=locusmatrix,
#'                                        num.draws=10, 
#'                                        use.dip.lincomb=TRUE, seed=1, gamma.rate=1, impute.on="CCline")
#' inla.diploffect.summary <- run.diploffect.inla.summary.stats(inla.diploffect)
run.diploffect.inla.summary.stats <- function(diploffect.object){
  deviation.ci <- deviation.mean <- diplotype.ci <- qtl.dom.ci <- SS.dom.ci <- SS.total.ci <- qtl.kinship.ci <- SS.poly.ci <- qtl.total.ci <- nonqtl.ci.list <- nonqtl.ss.ci.list <- NULL
  M <- diploffect.object$analysis.id$M

  ##### Handling additional non-qtl variance components
  cat("Loading nonqtl proportion of variance explained (var components and sums of squares)... \n")
  var.components <- diploffect.object$analysis.id$var.components
  nonqtl.var.components <- var.components[grepl(var.components, pattern="random.")]
  if(length(nonqtl.var.components) > 0){
    nonqtl.ci.list <- nonqtl.ss.ci.list <- list()
    for(i in 1:length(nonqtl.var.components)){
      nonqtl.ci.list[[i]] <- load.nonqtl.h2.from.inla(diploffect.inla=diploffect.object, variable=nonqtl.var.components[i], scale=1)
      nonqtl.ss.ci.list[[i]] <- load.nonqtl.ss.h2.from.inla(diploffect.inla=diploffect.object, variable=nonqtl.var.components[i], scale=1)
    }
    names(nonqtl.ci.list) <- names(nonqtl.ss.ci.list) <- nonqtl.var.components
  }

  #### Additive
  if("additive" %in% diploffect.object$analysis.id$genetic.effects){
    strain.ci <- load.ci.from.inla(M=M, diploffect.inla=diploffect.object)
    cat("Loading strain effects...\n")
    qtl.add.ci <- load.qtl.h2.from.inla(diploffect.inla=diploffect.object, effect="add", scale=1)
    cat("Loading additive heritability...\n")
    SS.add.ci <- load.SS.h2.from.inla(diploffect.inla=diploffect.object, effect="add", scale=1)
    cat("Loading additive sums of squares...\n")
  }

  ##### Dominant
  if("dominant" %in% diploffect.object$analysis.id$genetic.effects){
    deviation.ci <- load.deviated.ci.from.inla(M, diploffect.inla=diploffect.object, scale=1)
    cat("Loading deviation effects...\n")
    deviation.mean <- load.deviated.mean.from.inla(M, diploffect.inla=diploffect.object, scale=1)
    cat("Loading deviation means...\n")
    diplotype.ci <- load.diplotypes.ci.from.inla(M, diploffect.inla=diploffect.object, scale=1)
    cat("Loading diplotype means...\n")
    qtl.dom.ci <- load.qtl.h2.from.inla(diploffect.inla=diploffect.object, effect="dom", scale=1)
    cat("Loading dominant heritability...\n")
    qtl.total.ci <- load.total.qtl.h2.from.inla(diploffect.inla=diploffect.object, scale=1)
    cat("Loading total QTL heritability...\n")
    SS.dom.ci <- load.SS.h2.from.inla(diploffect.inla=diploffect.object, effect="dom", scale=1)
    cat("Loading dominant sums of squares...\n")
    SS.total.ci <- load.SS.h2.from.inla(diploffect.inla=diploffect.object, effect="qtl", scale=1)
    cat("Loading total QTL sums of squares...\n")
  }
  if("polygene" %in% diploffect.object$analysis.id$genetic.effects){
    qtl.kinship.ci <- load.qtl.h2.from.inla(diploffect.inla=diploffect.object, effect="kinship", scale=1)
    cat("Loading polygene heritability...\n")
    SS.poly.ci <- load.SS.h2.from.inla(diploffect.inla=diploffect.object, effect="poly", scale=1)
    cat("Loading polygene sums of squares...\n")
  }

  ci.results <- list(strain.ci=strain.ci, deviation.ci=deviation.ci, deviation.mean=deviation.mean, diplotype.ci=diplotype.ci,
                     qtl.total.ci=qtl.total.ci, qtl.add.ci=qtl.add.ci, qtl.dom.ci=qtl.dom.ci, kinship.ci=qtl.kinship.ci,
                     SS.total.ci=SS.total.ci, SS.add.ci=SS.add.ci, SS.dom.ci=SS.dom.ci, SS.poly.ci=SS.poly.ci,
                     nonqtl.ci.list=nonqtl.ci.list, nonqtl.ss.ci.list=nonqtl.ss.ci.list,
                     analysis.id=diploffect.object$analysis.id)
  return(ci.results)
}

make.M <- function(M){
  j <- M
  k <- (-1 + sqrt(j))*(j - 1)^(-3/2)
  m <- 1/sqrt(j - 1)
  M <- diag(j - 1)
  M[M == 0] <- -k
  M <- rbind(M, rep(-m, j - 1))
  return(M)
}
