## Return the nearest marker name to the position at chromosome chr.
## position is in bp
get.nearest.marker <- function(h, chr, position) {
  chromosomes = h$genotype$genome$chromosome
  bps = h$genotype$genome$bp
  markers = h$genotype$genome$marker
  chr.idx = which(chromosomes == chr)
  idx = which.min(abs(bps[chr.idx] - position))
  if (abs(bps[chr.idx[idx]] - position) > 1000000) {
    print("The nearest marker is still very far from the specified position.")
  }
  return(list(idx=idx, name=markers[chr.idx[idx]]))
}

incidence.matrix <- function(fact)
{
  m=diag(nlevels(fact))[fact,]
  colnames(m)=levels(fact)
  return(m)
}

shannon.entropy <- function(p){
  if (min(p) < 0 || sum(p) <= 0) {
    p[p<=0]=0
  }
  p.norm <- p[p>0]/sum(p)
  return(-sum(log2(p.norm)*p.norm))
}

straineff.mapping.matrix <- function(M=8, matrixtype){
  T=M*(M+1)/2
  TT=M*(M+1)/2
  mapping<-matrix(rep(0,T*M),M,T)
  if( matrixtype %in% c("hmm")){
    ## GAIN Matrix
    idx<-1;
    for (i in 1:M){
      for (j in 1:(i)){
        mapping[i,idx]<-mapping[i,idx]+1;
        mapping[j,idx]<-mapping[j,idx]+1;
        idx<-idx+1;
      }
    }
  }
  else {
    ## HAPPY matrix
    idx<-1;
    for (i in 1:M){
      mapping[i,idx]<- mapping[i,idx]+2
      idx<-idx+1;
    }
    for (i in 2:M){
      for (j in 1:(i-1)){
        mapping[i,idx]<-mapping[i,idx]+1;
        mapping[j,idx]<-mapping[j,idx]+1;
        idx<-idx+1;
      }
    }
  }
  return(mapping)
}

straineff.mapping.matrix.happy <- function(M=8){
  T=M*(M+1)/2
  TT=M*(M+1)/2
  mapping<-matrix(rep(0,T*M),M,T)
  ## HAPPY matrix
  idx<-1;
  for (i in 1:M){
    mapping[i,idx]<- mapping[i,idx]+2
    idx<-idx+1;
  }
  for (i in 2:M){
    for (j in 1:(i-1)){
      mapping[i,idx]<-mapping[i,idx]+1;
      mapping[j,idx]<-mapping[j,idx]+1;
      idx<-idx+1;
    }
  }
  return(mapping)
}

straineff.mapping.pair <- function(M, matrixtype){
  T = M * (M + 1) / 2
  mapping <- matrix(0, 2, T)
  if( matrixtype %in% c("happy")){
    ## HAPPY matrix
    idx<-1;
    for (i in 1:M){
      mapping[, idx] <- c(i, i)
      idx<-idx+1;
    }
    for (i in 2:M){
      for (j in 1:(i-1)){
        mapping[1, idx] <- i
        mapping[2, idx] <- j
        idx <- idx + 1;
      }
    }
  }
  else {
    stop("straineff.mapping.pair does not support other matrix type yet.")
  }
  mapping
}

straineff.smooth.probability.matrix <- function(N, data){
  p <- data[1:N,]
  for (i in 1:N){
    total = sum(p[i,]) + 0.0000036
    for (j in 1:36){
      p[i,j] = (p[i,j] + 0.0000001) / total
    }
  }
  p <- t(matrix(unlist(p), ncol=N, byrow=TRUE))
  p
}

straineff.get.posterior.matrix <- function(M, N, mcmc.matrix){
  TT=M*(M+1)/2
  data <- mat.or.vec(N,M)
  x <- mat.or.vec(N,TT)
  lmapping<-mat.or.vec(TT,1)
  rmapping<-mat.or.vec(TT,1)
  idx<-1;
  for (i in 1:M){
    lmapping[idx]<- i
    rmapping[idx]<- i
    idx<-idx+1;
  }
  for (i in 2:M){
    for (j in 1:(i-1)){
      lmapping[idx]<- i
      rmapping[idx]<- j
      idx<-idx+1;
    }
  }
  ss=mcmc.matrix[[1]]
  for (i in 1:N){
    r=table(ss[[1]][,paste('idx[',i,']',sep="")])
    x[i,as.numeric(names(r))]=r/(sum(r))
  }
  x
}

straineff.prior.entropy <- function(N,data){
  ret <- apply(data,1,shannon.entropy)
  ret <- sum(unlist(ret))
  ret
}

sic <- function(p) {
  if (min(p) < 0 || sum(p) <= 0) {
    p[p <= 0] = 0
  }
  p.norm <- p[p > 0] / sum(p)
  N <- length(p.norm)
  p.norm <- p.norm[which(p.norm >= 0.00000001)]
  sum(p.norm*log(p.norm/rep(1/N, length(p.norm))))
}

straineff.prior.sic <- function(N, data) {
  ret <- apply(data, 1, sic)
  sum(unlist(ret)) / N
}

straineff.prior.sic.vector <- function(N, data) {
  ret <- apply(data, 1, sic)
  unlist(ret)
}

straineff.posterior.entropy <- function(N,mcmc.matrix){
  straineff.prior.entropy(N,straineff.get.posterior.matrix(8,mcmc.matrix))
}

straineff.true.founder.prior <- function(Y,N,M,phenotype.name,Y1,data){
  Z1=Y1
  idx<-1;
  happymap <- mat.or.vec(M,M)
  for (i in 1:M){
    happymap[i,i]<-idx
    idx <- idx+1
  }
  for (i in 2:M){
    for (j in 1:(i-1)){
      happymap[i,j]<-idx
      idx<-idx+1;
    }
  }
  mapping=as.numeric(phenotype.name)
  p <- mat.or.vec(N,M)
  for ( i in 1:N){
    Z1[i,]=Y1[mapping[i],]
  }
  ret=0
  for (i in 1:N){
    x=as.numeric(Z1[i,3])
    y=as.numeric(Z1[i,5])

    if (x>y)
      ret=ret+data[i,happymap[x,y]]
    else
      ret=ret+data[i,happymap[y,x]]
  }
  ret
}
straineff.true.founder.posterior <- function(Y,N,M,phenotype.name,Y1,mcmc.matrix){
  straineff.true.founder.prior(Y,N,M,phenotype.name,Y1,straineff.get.posterior.matrix(8,mcmc.matrix))
}

straineff.extra <- function(H, h) {
  extra <- NULL
  if (H > 1) extra <- paste(h, ',', sep="")
  extra
}

get.extra <- function(H=1, h=1) {
  extra <- NULL
  if (H > 1) extra <- paste(h, ',', sep="")
  return(extra)
}

summarize.beta <- function(M, dat, H=1, h=1){
  beta = mat.or.vec(M,niter(dat))
  mbeta = mat.or.vec(M,1)
  extra <- get.extra(H, h)
  for (i in 1:M){
    beta[i,] = dat[,paste('beta[',extra, i,']',sep="")]
  }
  for (i in 1:niter(dat)){
    beta[,i] = beta[,i] - mean(beta[,i])
  }
  for (i in 1:M){
    mbeta[i] = mean(beta[i,])
  }
  return (mbeta)
}

summarize.weighted.beta <- function(M, N, mcmc.mat, haplotype.prior, H=1, h=1){
  beta = mat.or.vec(M, niter(mcmc.mat))
  w = mat.or.vec(niter(mcmc.mat), 1)
  mbeta = mat.or.vec(M, 1)
  extra <- get.extra(H, h)
  for (i in 1:M){
    beta[i,] = mcmc.mat[, paste('beta[', extra, i, ']', sep="")]
  }
  ## mcmc.mat[, 'deviance'] is the deviance defined in JAGS (not proper
  ## deviance though)
  ## defined as -2*logDensity(model)
  ## convert it back to log likelihood
  w = mcmc.mat[, 'deviance'] * (-0.5)
  lprior = log.haplotype.prior(N, haplotype.prior, mcmc.mat)
  w = exp(w)
  for (i in 1:niter(mcmc.mat)) {
    beta[,i] = beta[,i] - mean(beta[,i])
  }
  for (i in 1:M) {
    mbeta[i] = weighted.mean(beta[i,], w)
  }
  print(mbeta)
  print(summarize.beta(M, mcmc.mat))
  return (mbeta)
}

traces.deviated.effects <- function(M, dat, H=1, h=1) {
  total = M * (M + 1) / 2
  gamma = mat.or.vec(total, niter(dat))
  extra <- get.extra(H, h)
  for (j in (M + 1):total) {
    gamma[j, ] = dat[, paste('gamma[', extra, j ,']', sep="")]
  }
  return (gamma)
}

traces.diplotype.effects <- function(M, dat, deviated, H=1, h=1) {
  total = M * (M + 1) / 2
  beta = mat.or.vec(M, niter(dat))
  gamma = mat.or.vec(total, niter(dat))
  diplotype.effects = mat.or.vec(total, niter(dat))
  extra <- get.extra(H, h)
  for (j in 1:M){
    beta[j, ] = dat[, paste('beta[', extra, j, ']', sep="")]
  }
  if (deviated) {
    for (j in (M + 1):total){
      gamma[j, ] = dat[, paste('gamma[', extra, j ,']', sep="")]
    }
  }
  mapping.matrix <- straineff.mapping.matrix(M, "happy")
  for (i in 1:niter(dat)){
    for (j in 1:total) {
      diplotype.effects[j, i] = t(mapping.matrix[, j]) %*% beta[, i]
      if (deviated) {
        diplotype.effects[j, i] = diplotype.effects[j, i] + gamma[j, i]
      }
    }
    diplotype.effects[, i] = diplotype.effects[, i] - mean(diplotype.effects[, i])
  }
  return(diplotype.effects)
}

summarize.diplotype.effects <- function(M, dat, deviated, H=1, h=1) {
  total = M * (M + 1) / 2
  diplotype.effects <- traces.diplotype.effects(M, dat, deviated, H, h)
  mean.diplotype.effects = mat.or.vec(1, total)
  for (j in 1:total){
    mean.diplotype.effects[j] = mean(diplotype.effects[j, ])
  }
  return(mean.diplotype.effects)
}

summarize.haplotype <- function(N, data, true.haplotype, mcmc.matrix){
  n.iter = niter(mcmc.matrix)
  map2 = data
  delta = 0
  for ( i in 1:N){
    ret = sort.int(data[i,], index.return=T)
    map2[i,] = ret$ix
    s = mcmc.matrix[,paste('idx[', i, ']', sep="")]
    estimate = (sum(as.numeric(map2[i, s]) == true.haplotype[1])/n.iter)
    real = data[i, true.haplotype[1]]
    delta = delta + (estimate - real)
  }
  return (delta/N)
}

log.haplotype.prior <- function(N, data, mcmc.matrix){
  n.iter = niter(mcmc.matrix)
  map2 = data
  r = rep(0, n.iter)
  for ( i in 1:N) {
    ret = sort.int(data[i,], index.return=T)
    map2[i,] = ret$ix
    s = mcmc.matrix[,paste('idx[', i, ']', sep="")]
    r = r + log(data[i, map2[i, s]])
  }
  return (r)
}

calculate.diplotype.effects <- function(beta, deviation.effects) {
  M <- length(beta)
  mapping.matrix <- straineff.mapping.matrix(M, "happy")
  total = M * (M + 1) / 2
  diplotype.effects <- mat.or.vec(total, 1)
  for (j in 1:total) {
    diplotype.effects[j] = t(mapping.matrix[, j]) %*% beta + deviation.effects[j]
  }
  return(diplotype.effects)
}

remove.small.probability <- function(data, numprop=36){
  MM = dim(data)[2]
  for (i in 1:dim(data)[1]) {
    x = sort(data[i, ], index.return=T)
    data[i,  which(x$ix <= MM - numprop)]=0
    data[ i, ]=data[i, ]/sum(data[i,  ])
  }
  return(data)
}
