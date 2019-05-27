get.ci.from.mcmc <- function(dat, wanted) {
  dat[,wanted] = dat[,wanted]- rowMeans(dat[,wanted])
  which.wanted=ifow(is.integer(wanted), wanted, match(wanted, varnames(dat)))
  num.wanted=length(which.wanted)
  chain <- mcmc.stack(dat)
  mu    <- colMeans(chain[,which.wanted])
  med   <- apply(coda::HPDinterval(chain, prob=0.01)[which.wanted,],
                 1, mean)
  prob.wide=0.95
  prob.narrow=0.50
  hpd.wide    <- coda::HPDinterval(chain, prob=prob.wide)[which.wanted,]
  hpd.narrow  <- coda::HPDinterval(chain, prob=prob.narrow)[which.wanted,]
  return (list(med=med, mu=mu, hpd.narrow=hpd.narrow, hpd.wide=hpd.wide))
}

load.ci.from.mcmc <- function(M, file, H=1, h=1) {
  load(file)
  mcmc = result[[1]][[1]]
  extra <- straineff.extra(H, h)
  ci <- get.ci.from.mcmc(mcmc, wanted=paste('beta[', extra, 1:M, ']', sep=""))
  ci
}

plot.ci <- function(midvals, narrow.intervals, wide.intervals,
                    names=1:length(midvals),
                    add=FALSE,
                    xlab="Estimate", xlab.line=2.5, xlab.cex=1,
                    ylab="",
                    yaxis=TRUE,
                    name.margin=6,
                    name.line=4,
                    pch.midvals=19,
                    col=rep("black", length(midvals)),
                    col.midvals=col,
                    include.top.axis=TRUE,
                    shift=0,
                    type="p",
                    use.this.lim=NULL,
                    main="", 
                    main.cex=1, 
                    main.line=1,
                    ...)
{
  nvals <- length(midvals)
  col.midvals <- rep(col.midvals, length.out=nvals)
  y.pos <- (1:nvals)-0.5
  if(!add){
    if(is.null(use.this.lim)){
      lim <- range(c(wide.intervals,narrow.intervals,midvals), na.rm=TRUE)
      lim <- c(-1,1) * diff(lim)*0.1 + lim
    }
    if(!is.null(use.this.lim)){
      lim <- use.this.lim
    }

    mar <- c(5, name.margin, 4, 2)+0.1
    oldmar <- par(mar=mar); on.exit(par(mar=oldmar))
    plot(lim, c(0,nvals+0.5), type="n", axes=FALSE, ylab=ylab, xlab="", main=NA, ...)
    title(xlab=xlab, line=xlab.line, cex.main=xlab.cex)
    title(main=main, line=main.line, cex.main=main.cex)
    axis(1)
    if(include.top.axis){ axis(3, line=-1) }
    if(yaxis){
      axis(2, at=y.pos, labels=rev(names), las=1, lty=0, hadj=0, line=name.line)
    }
  }
  if("p"==type){
    for(i in 1:nvals){
      pos <- nvals-i + 0.5 + shift
      lines(wide.intervals[i,], rep(pos,2), col=col[i])
      lines(narrow.intervals[i,], rep(pos,2), lwd=3, col=col[i])
      points(midvals[i], pos, pch=pch.midvals, col=col.midvals[i])
    }
  }
  invisible(rev(y.pos))
}

prepare.additive.dominant.ratio.posterior <- function(files, num.draws=1000) {
  samples = c()
  weights = c()
  if (!is.vector(files)) {
    files = c(files)
  }
  for (file in files) {
    load(file)
    print(file)
    logmliks = as.vector(scale(result$results$logmliks, scale=F))
    mliks = exp(logmliks)
    ratio = lapply(1:length(mliks), function (i) {
      hyper.samples = result$results$hyper.samples[[i]]
      ## bug in inla, this is not log scale
      precision.add = hyper.samples[, 'Log precision for idx1 in user-scale']
      precision.dom = hyper.samples[, 'Log precision for dom.idx in user-scale']

      (1 / (precision.add ^ 2)) / (1 / (precision.add ^ 2) + 1 / (precision.dom ^2))
    })

    samples = as.vector(unlist(ratio))
    weights = rep(mliks, each=num.draws)
  }
  weights = weights / sum(weights)
  density(samples, weights=weights, from=0, to=1, kernel='cosine')
}

#' Plot approximate posterior credible intervals of additive haplotype effects as a caterpillar plot
#'
#' This function takes a Diploffect.INLA summary object and plots the posterior credible intervals of
#' the additive haplotype effects.
#'
#' @param inla.diploffect.ci A Diploffect.INLA summary object. Output from run.diploffect.inla.summary.stats(). Contains credible information
#' for various model parameters.
#' @param sn DEFAULT: NULL. Strain names to be used as labels in the plot. Will default to the names used in the Diploffect.INLA model probabilities.
#' @param xlab DEFAULT: "Haplotype Effects". Label for the x axis of the plot.
#' @param main DEFAULT: "". Title for the plot.
#' @param flip DEFAULT: TRUE. Flips the order of the haplotypes. Allows a little flexibility in how the effects are arranged.
#' @return Nothing. Produces plot.
#' @export plot.straineff.ci
#' @examples
#' library(Diploffect.INLA)
#' data(exampleCC)
#' data(locusmatrix)
#' inla.diploffect <- run.diploffect.inla(formula=y~1+(1|strain)+(1|dose.date), add.on=FALSE, data=exampleCC, K=NULL,
#'                                        prob.matrix=locusmatrix,
#'                                        num.draws=10, 
#'                                        use.dip.lincomb=TRUE, seed=1, gamma.rate=1, impute.on="CCline")
#' inla.diploffect.summary <- run.diploffect.inla.summary.stats(inla.diploffect)
#' plot.straineff.ci(inla.diploffect.summary, flip=FALSE)
plot.straineff.ci <- function(inla.diploffect.ci, 
                              sn=NULL, 
                              xlab="Haplotype Effects", 
                              main=NULL, 
                              main.cex=1, 
                              main.line=2,
                              include.top.axis=TRUE,
                              flip=TRUE, ...) {
  ci <- inla.diploffect.ci$strain.ci
  if(is.null(sn)){
    sn <- inla.diploffect.ci$analysis.id$founders
  }

  if(flip){ order <- rev(1:nrow(ci$quant.narrow)) }
  if(!flip){ order <- 1:nrow(ci$quant.narrow) }
  
  if(is.null(main)){
    main <- c(paste(inla.diploffect.ci$analysis.id$formula, paste0("locus(", inla.diploffect.ci$analysis.id$locus, ")"), sep=" + "),
              paste("INLA samples:", inla.diploffect.ci$analysis.id$num.draws))
  }
  
  ypos <- plot.ci(ci$med[order], ci$quant.narrow[order,], ci$quant.wide[order,], names=sn[order],
                  xlab=xlab, col.midvals="white",
                  pch.midvals="|", type="p", 
                  main=main, main.cex=main.cex, main.line=main.line,
                  include.top.axis=include.top.axis, ...)
  points(ci$mu[order], ypos, pch="|")
  abline(v=0, lty=2)
}

#' Plot approximate posterior credible intervals of heterozygous diplotype dominance effects as a caterpillar plot
#'
#' This function takes a Diploffect.INLA summary object and plots the posterior credible intervals of
#' the dominance effects, corresponding to deviations for heterozygous diplotypes.
#'
#' @param inla.diploffect.ci A Diploffect.INLA summary object. Output from run.diploffect.inla.summary.stats(). Contains credible information
#' for various model parameters.
#' @param sn DEFAULT: NULL. Strain names to be used as labels in the plot. Will default to the names used in the Diploffect.INLA model probabilities.
#' @param xlab DEFAULT: "Dominant Deviation Effects". Label for the x axis of the plot.
#' @param flip DEFAULT: TRUE. Flips the order of the haplotypes. Allows a little flexibility in how the effects are arranged.
#' @return Nothing. Produces plot.
#' @export plot.deviation.ci
#' @examples
#' library(Diploffect.INLA)
#' data(exampleCC)
#' data(locusmatrix)
#' inla.diploffect <- run.diploffect.inla(formula=y~1+(1|strain)+(1|dose.date), add.on=FALSE, data=exampleCC, K=NULL,
#'                                        prob.matrix=locusmatrix,
#'                                        num.draws=10, 
#'                                        use.dip.lincomb=TRUE, seed=1, gamma.rate=1, impute.on="CCline")
#' inla.diploffect.summary <- run.diploffect.inla.summary.stats(inla.diploffect)
#' plot.deviation.ci(inla.diploffect.summary, flip=FALSE)
plot.deviation.ci <- function(inla.diploffect.ci, sn=NULL, xlab="Dominant Deviation Effects", flip=TRUE, ...) {
  ci <- inla.diploffect.ci$deviation.ci
  if(is.null(sn)){
    founders <- inla.diploffect.ci$analysis.id$founders
    full.to.dosages <- t(straineff.mapping.matrix.happy(M=length(founders)))

    sn <- apply(full.to.dosages[-(1:8),], 1, function(x) paste(founders[sort(which(x==1), decreasing=FALSE)], collapse=" x "))
  }
  if(flip){ order <- rev(1:nrow(ci$quant.narrow)) }
  if(!flip){ order <- 1:nrow(ci$quant.narrow) }
  main <- c(paste(inla.diploffect.ci$analysis.id$formula, paste0("locus(", inla.diploffect.ci$analysis.id$locus, ")"), sep=" + "),
            paste("INLA samples:", inla.diploffect.ci$analysis.id$num.draws))
  ypos <- plot.ci(ci$med[order], ci$quant.narrow[order,], ci$quant.wide[order,], names=sn[order],
                  xlab=xlab, col.midvals="white",
                  pch.midvals="|", type="p", main=main, ...)
  points(ci$mu[order], ypos, pch="|")
  abline(v=0, lty=2)
}

#' Plot approximate posterior credible intervals of diplotype effects as a caterpillar plot
#'
#' This function takes a Diploffect.INLA summary object and plots the posterior credible intervals of
#' the diplotype effects, which are a linear combination of the additive and dominant effects.
#'
#' @param inla.diploffect.ci A Diploffect.INLA summary object. Output from run.diploffect.inla.summary.stats(). Contains credible information
#' for various model parameters.
#' @param sn DEFAULT: NULL. Strain names to be used as labels in the plot. Will default to the names used in the Diploffect.INLA model probabilities.
#' @param xlab DEFAULT: "Diplotype Effects". Label for the x axis of the plot.
#' @param flip DEFAULT: TRUE. Flips the order of the haplotypes. Allows a little flexibility in how the effects are arranged.
#' @return Nothing. Produces plot.
#' @export plot.diplotype.ci
#' @examples
#' library(Diploffect.INLA)
#' data(exampleCC)
#' data(locusmatrix)
#' inla.diploffect <- run.diploffect.inla(formula=y~1+(1|strain)+(1|dose.date), add.on=FALSE, data=exampleCC, K=NULL,
#'                                        prob.matrix=locusmatrix,
#'                                        num.draws=10, 
#'                                        use.dip.lincomb=TRUE, seed=1, gamma.rate=1, impute.on="CCline")
#' inla.diploffect.summary <- run.diploffect.inla.summary.stats(inla.diploffect)
#' plot_diplotype.ci(inla.diploffect.summary, flip=FALSE)
plot.diplotype.ci <- function(inla.diploffect.ci, sn=NULL, xlab="Diplotype Effects", flip=TRUE, ...) {
  ci <- inla.diploffect.ci$diplotype.ci
  if(is.null(sn)){
    founders <- inla.diploffect.ci$analysis.id$founders
    full.to.dosages <- t(straineff.mapping.matrix.happy(M=length(founders)))

    sn <- c(paste(founders, founders, sep=" x "),
            apply(full.to.dosages[-(1:8),], 1, function(x) paste(founders[sort(which(x==1), decreasing=FALSE)], collapse=" x ")))
  }
  if(flip){ order <- rev(1:nrow(ci$quant.narrow)) }
  if(!flip){ order <- 1:nrow(ci$quant.narrow) }
  main <- c(paste(inla.diploffect.ci$analysis.id$formula, paste0("locus(", inla.diploffect.ci$analysis.id$locus, ")"), sep=" + "),
            paste("INLA samples:", inla.diploffect.ci$analysis.id$num.draws))
  par(mar=c(3.2,3.2,2,1), mgp=c(2.2,.7,0), tck=-.01, las=1)
  ypos <- plot.ci(ci$med[order], ci$quant.narrow[order,], ci$quant.wide[order,], names=sn[order],
                  xlab=xlab, col.midvals="white",
                  pch.midvals="|", type="p", main=main, ...)
  points(ci$mu[order], ypos, pch="|")
  abline(v=0, lty=2)
}

plot.comparison.cis <- function(ci.list, 
                                analysis.id, 
                                labels=NULL, 
                                sn, 
                                main="", 
                                main.cex=1, 
                                main.line=1,
                                xlab="Haplotype Effects", 
                                add.numbers=FALSE,
                                use.this.lim=NULL, 
                                comp.col=c("black", "red"), ...)
  {
  par(mar=c(3.2,3.2,2,1), mgp=c(2.2,.7,0), tck=-.01, las=1)
  # main <- c(paste(analysis.id$formula, paste0("locus(", analysis.id$locus, ")"), sep=" + "),
  #           paste("INLA samples:", analysis.id$num.draws))
  add = FALSE
  step = 0.4 / (length(ci.list) - 1)
  shift = 0.2 ## start
  if(is.null(use.this.lim)){
    lim <- range(ci.list, na.rm=TRUE)
    lim <- c(-1,1) * diff(lim)*0.1 + lim
  }
  else{
    lim <- use.this.lim
  }
  for(i in 1:length(ci.list)){
    ci <- ci.list[[i]]
    if(add){
      ypos <- plot.ci(ci$med, ci[[3]], ci[[4]], names=sn, add=add,
                      xlab=xlab, col.midvals="white", col=rep(comp.col[i], length(ci$med)),
                      pch.midvals="|", type="p", shift=shift, main=main, use.this.lim=lim, ...)
    }
    if(!add){
      ypos <- plot.ci(ci$med, ci[[3]], ci[[4]], names=sn, add=add,
                      xlab=xlab, col.midvals="white", col=rep(comp.col[i], length(ci$med)),
                      pch.midvals="|", type="p", shift=shift, use.this.lim=lim, main=main, ...)
    }
    points(ci$mu, ypos + shift, pch="|")
    if(add.numbers){
      mtext(text=paste(round(ci$mu, 2), paste0("(", round(ci[[4]][, 1], 2), ", ", round(ci[[4]][, 2], 2), ")")), side=4, at=ypos + shift, adj=1, col=comp.col[i])
    }
    if(!add){ add <- TRUE }
    shift <- shift - step
  }
  abline(v=0, lty=2)

  if(!is.null(labels)){
    legend('bottomleft', labels,
           col=1:length(labels), lty=1, bty='n', lwd=3)
  }
}

#' Plot approximate posterior credible intervals of proportion variance explained (PVE) as a caterpillar plot
#'
#' This function takes a Diploffect.INLA summary object and plots the posterior credible intervals of
#' the proportion of the variance explained.
#'
#' @param inla.diploffect.ci A Diploffect.INLA summary object. Output from run.diploffect.inla.summary.stats(). Contains credible information
#' for various model parameters.
#' @param xlab DEFAULT: "Variance Explained". Label for the x axis of the plot.
#' @param add.numbers DEFAULT: FALSE. Adds numeric summaries next to intervals.
#' @return Nothing. Produces plot.
#' @export plot.varexp.ci
#' @examples
#' library(Diploffect.INLA)
#' data(exampleCC)
#' data(locusmatrix)
#' inla.diploffect <- run.diploffect.inla(formula=y~1+(1|strain)+(1|dose.date), add.on=FALSE, data=exampleCC, K=NULL,
#'                                        prob.matrix=locusmatrix,
#'                                        num.draws=10, 
#'                                        use.dip.lincomb=TRUE, seed=1, gamma.rate=1, impute.on="CCline")
#' inla.diploffect.summary <- run.diploffect.inla.summary.stats(inla.diploffect)
#' plot.varexp.ci(inla.diploffect.summary, add.numbers=TRUE)
plot.varexp.ci <- function(inla.diploffect.ci, 
                           xlab="Variance Explained",
                           main=NULL, 
                           main.cex=1, 
                           main.line=2,
                           add.numbers=FALSE, ...){
  join.ci <- function(ci.list) {
    combine.ci <- list()
    med <- NULL
    mu <- NULL
    quant.narrow <- NULL
    quant.wide <- NULL
    for (i in 1:length(ci.list)) {
      med <- c(med, ci.list[[i]]$med)
      mu <- c(mu, ci.list[[i]]$mu)
      quant.narrow <- rbind(quant.narrow, ci.list[[i]]$quant.narrow)
      quant.wide <- rbind(quant.wide, ci.list[[i]]$quant.wide)
    }
    combine.ci <- list(med=med, mu=mu, quant.narrow=quant.narrow, quant.wide=quant.wide)
    return(combine.ci)
  }

  SS.ci.list <- list()
  h2.ci.list <- list()
  
  if(is.null(main)){
    main <- c(paste(inla.diploffect.ci$analysis.id$formula, paste0("locus(", inla.diploffect.ci$analysis.id$locus, ")"), sep=" + "),
              paste("INLA samples:", inla.diploffect.ci$analysis.id$num.draws))
  }

  effect.labels <- NULL
  if(!is.null(inla.diploffect.ci$nonqtl.ci.list)){
    for(i in 1:length(inla.diploffect.ci$nonqtl.ci.list)){
      h2.ci.list[[i]] <- inla.diploffect.ci$nonqtl.ci.list[[i]]
      SS.ci.list[[i]] <- inla.diploffect.ci$nonqtl.ss.ci.list[[i]]
      effect.labels <- c(effect.labels, names(inla.diploffect.ci$nonqtl.ci.list)[i])
    }
  }

  SS.ci.list[[length(SS.ci.list)+1]] <- inla.diploffect.ci$SS.add.ci
  h2.ci.list[[length(h2.ci.list)+1]] <- inla.diploffect.ci$qtl.add.ci
  effect.labels <- c(effect.labels, "QTL Additive")
  if("dominant" %in% inla.diploffect.ci$analysis.id$genetic.effects){
    SS.ci.list[[length(SS.ci.list)+1]] <- inla.diploffect.ci$SS.dom.ci
    SS.ci.list[[length(SS.ci.list)+1]] <- inla.diploffect.ci$SS.total.ci

    h2.ci.list[[length(h2.ci.list)+1]] <- inla.diploffect.ci$qtl.dom.ci
    h2.ci.list[[length(h2.ci.list)+1]] <- inla.diploffect.ci$qtl.total.ci

    effect.labels <- c(effect.labels, "QTL Dominant", "QTL Total")
  }
  if("polygene" %in% inla.diploffect.ci$analysis.id$genetic.effects){
    SS.ci.list[[length(SS.ci.list)+1]] <- inla.diploffect.ci$SS.poly.ci
    h2.ci.list[[length(h2.ci.list)+1]] <- inla.diploffect.ci$kinship.ci

    effect.labels <- c(effect.labels, "Polygene")
  }

  SS.ci <- join.ci(SS.ci.list)
  h2.ci <- join.ci(h2.ci.list)

  plot.comparison.cis(ci.list=list(SS.ci, h2.ci),
                      analysis.id=inla.diploffect.ci$analysis.id,
                      labels=c("SS", "VC"), 
                      sn=effect.labels, 
                      xlab=xlab, 
                      add.numbers=add.numbers, 
                      use.this.lim=c(0, 1),
                      main=main, main.cex=main.cex, main.line=main.line,
                      ...)
  abline(v=1, lty=2)
}

plot.diallel <- function(inla.diploffect.ci, sn=NULL) {
  if(is.null(sn)){
    sn <- inla.diploffect.ci$analysis.id$founders
  }
  diplotype <- inla.diploffect.ci$diplotype.ci$mu
  M <- inla.diploffect.ci$analysis.id$M
  main <- c(paste(inla.diploffect.ci$analysis.id$formula, paste0("locus(", inla.diploffect.ci$analysis.id$locus, ")"), sep=" + "),
            paste("INLA samples:", inla.diploffect.ci$analysis.id$num.draws))
  par(mar=c(3.5, 3.2, 4.5, 1.0))
  mapping = straineff.mapping.matrix(M, 'happy')
  data.matrix = matrix(0, M, M)
  for (i in 1:ncol(mapping)) {
    matrix.idx <- which(mapping[,i] != 0)
    if (length(matrix.idx) == 1) {
      data.matrix[matrix.idx, matrix.idx] <- diplotype[i]
    }
    else if (length(matrix.idx) == 2) {
      data.matrix[matrix.idx[1], matrix.idx[2]] <- data.matrix[matrix.idx[2], matrix.idx[1]] <- diplotype[i]
    }
  }
  draw.diallel(data.matrix, sn=sn, main=main)
}
