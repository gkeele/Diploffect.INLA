#' Run Diploffect model on diplotype (haplotype pair) probabilities stored in a genomecache directory
#'
#' This function primarily takes a formula, locus, data frame, and genomecache and runs the
#' corresponding Diploffect model for a given number of samples.
#'
#' @param formula An lm style formula with functions of outcome and covariates in data frame.
#' @param data A data frame with outcome and potential covariates. Should also have individual IDs
#' that link to IDs in the genomecache, commonly with this column named "SUBJECT.NAME".
#' @param K A positive semi-definite relationship matrix, usually a realized genetic relationship matrix (GRM)
#' based on SNP genotypes are the haplotype probabilistic reconstructions. Colnames and rownames should match
#' the SUBJECT.NAME column in the data frame.
#' @param genomecache The path to the genomecache directory. The genomecache is a particularly structured
#' directory that stores the haplotype probabilities/dosages at each locus. It has an additive model
#' subdirectory and a full model subdirectory. Each contains subdirectories for each chromosome, which then
#' store .RData files for the probabilities/dosages.
#' @param locus The locus for which to run Diploffect. It must be contained in the genomecache.
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
#' of a genome, such as in Collaborative Cross (CC) or CC-RIX.
#' @param weights.on DEFAULT: NULL. If numeric column in data frame is specified, weights are applied to the
#' diagonal of the random noise (diagonal) covariance matrix. This is an efficient approach to modeling
#' repeat observations, such as in CC or CC-RIX.
#' @param scale.response DEFAULT: TRUE. It is good to scale the response with INLA, thus the posteriors of effects
#' have a consistent scale (in standard deviations of the outcome).
#' @param do.augment.of.cache DEFAULT: FALSE. Not suggested. Could be used for CC genomecaches with only additive
#' dosages, though full probabilities are still greatly preferred.
#' @return Diploffect.INLA object. Approximate posterior distributions for model parameters.
#' @export
#' @examples run.diploffect.inla.through.genomecache()
run.diploffect.inla.through.genomecache <- function(formula, data, K=NULL,
                                                    genomecache, locus,
                                                    num.draws, use.dip.lincomb=TRUE,
                                                    seed=1, gamma.rate=1, impute.on="SUBJECT.NAME", weights.on=NULL,  #"NUM.OBS",
                                                    scale.response=TRUE,
                                                    do.augment.of.cache=FALSE){
  weights.on <- weights.on[1]

  QTL.effect <- process.formula(formula, action="return.QTL.effect")
  if(QTL.effect == "additive") { add.only <- TRUE}
  if(QTL.effect == "full") { add.only <- FALSE }

  formula.string.new.data <- process.formula(formula, action="make.new.data.formula", impute.on=c(impute.on, weights.on))
  new.data <- model.frame(formula(formula.string.new.data), data=data)
  h <- bagpipe.backend::DiploprobReader$new(genomecache)
  founders <- h$getFounders()
  num.founders <- length(founders)
  # cache has rows per individual, not line
  if(!do.augment.of.cache){
    diplotype.matrix <- h$getLocusMatrix(locus=locus, model="full", subjects=new.data$SUBJECT.NAME)
    diplotype.matrix[diplotype.matrix < 0] <- 0
    # pick overlap between data and genomecache
    data <- data[as.character(data$SUBJECT.NAME) %in% h$getSubjects(),]
    new.data <- new.data[as.character(new.data$SUBJECT.NAME) %in% h$getSubjects(),]
  }
  # cache has rows per line
  if(do.augment.of.cache){
    diplotype.matrix <- h$getLocusMatrix(locus=locus, model="full", subjects=new.data[, impute.on])
    diplotype.matrix[diplotype.matrix < 0] <- 0

    rownames(diplotype.matrix) <- new.data$SUBJECT.NAME
    # pick overlap between data and genomecache
    data <- data[as.character(data[,impute.on]) %in% h$getSubjects(),]
    new.data <- new.data[as.character(new.data[,impute.on]) %in% h$getSubjects(),]
  }
  # pick overlap between data and genomecache
  diplotype.matrix <- diplotype.matrix[as.character(new.data$SUBJECT.NAME),]
  diploffect.inla.object <- run.diploffect.inla(formula=formula, data=data, K=K,
                                                num.founders=num.founders, prob.matrix=diplotype.matrix, add.only=add.only,
                                                num.draws=num.draws, use.dip.lincomb=use.dip.lincomb,
                                                seed=seed, gamma.rate=gamma.rate,
                                                impute.on=impute.on, weights.on=weights.on, scale.response=scale.response, founders=founders, locus.name=locus)
  return(diploffect.inla.object)
}
