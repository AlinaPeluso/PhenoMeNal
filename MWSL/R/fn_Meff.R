#' Estimation of closed-form-expression metabolome-wide significance level (MWSL) and corresponding effective number of tests (ENT)
#'
#' Approximation methods based on the spectral decomposition of the correlation matrix of the metabolomics variates are employed to estimate a closed-form-expression significance level to be used for multiple tests adjustement in the case of correlated metabolomic variates.
#' @param \code{features} a data frame of \code{n} samples (rows) and \code{M} features (columns)
#' @param \code{method} a string with possible values \code{ecorr} (empirical correlation), or \code{scorr} (shrinkage correlation), default=\code{ecorr}
#' @param \code{alpha} an optional probability value, default=0.05
#' @param \code{big.mat} an optional logic value to be set to \code{TRUE} when dealing with a very large set of feature (M>>5000), default=\code{FALSE}
#'
#' @return \code{Meff_Nyholt} closed-form expression of the effective number of tests based on Nyholt (2004)
#' @return \code{Meff_Liji} closed-form expression of the effective number of tests based on Li & Ji (2005)
#' @return \code{Meff_Gao} closed-form expression of the effective number of tests based on Gao (2008)
#' @return \code{Meff_Galwey} closed-form expression of the effective number of tests based on Galwey (2009)
#' @return \code{Meff_Bonferroni} closed-form expression of the effective number of tests based on Bonferroni (1963)
#' @return \code{Meff_Sidak} closed-form expression of the effective number of tests based on Sidak (1967)
#' @return \code{Meff_MWSL} closed-form expression of the effective number of tests based on Peluso (2020)
#' @return \code{res.Meff_MWSL} is the vector of result estimates which includes:
#' \code{MWSL} = metabolome-wide significance level (MWSL),
#' \code{MWSL_CI.up} = upper value confidence interval MWSL,
#' \code{MWSL_CI.low} = lower value confidence interval MWSL,
#' \code{ENT} = effective number of tests (ENT),
#' \code{ENT_CI.up} = upper value confidence interval ENT,
#' \code{ENT_CI.low} = lower value confidence interval ENT,
#' \code{R.percent} = \code{ENT}/\code{M}
#' @export
#'

Meff <- function(features,n.permutation=NULL,method='ecorr',big.mat=FALSE,alpha=NULL,verbose=TRUE){

  n <- nrow(features)
  M <- ncol(features)
  K <- ifelse(is.null(n.permutation)==T,10000,n.permutation)
  method <- ifelse(is.null(method)==T,'ecorr',method)
  alpha <- ifelse(is.null(alpha)==T,0.05,alpha)

  ### compute the **shrinkage** correlation matrix
  if (method=='scorr'){
    if (big.mat==TRUE){
      big.scorr_Xmat = bigcor.shrink(x=features, size= M/10, fun = "cor")
      big.scorr_Xmat <- big.scorr_Xmat[1:nrow(big.scorr_Xmat), 1:ncol(big.scorr_Xmat)]
      #evs = eigen(big.scorr_Xmat)$values
      evs <- RSpectra::eigs(big.scorr_Xmat , k=ncol(big.scorr_Xmat), opts = list(retvec = FALSE))
      evs <- as.numeric(evs$values)
      abs_evs <- abs(evs)
      levs <- length(evs)
    }else{
    scorr_Xmat= corpcor::cor.shrink(features)
    evs = eigen(scorr_Xmat)$values
    abs_evs <- abs(evs)
    levs <- length(evs)
    }
  }

  ### compute the **empirical** correlation matrix
  if (method=='ecorr'){
    if (big.mat==TRUE){
      big.ecorr_Xmat = propagate::bigcor(x=features, size= M/10, fun = "cor")
      big.ecorr_Xmat <- big.ecorr_Xmat[1:nrow(big.ecorr_Xmat), 1:ncol(big.ecorr_Xmat)]
      #evs = eigen(big.ecorr_Xmat)$values
      evs <- RSpectra::eigs(big.ecorr_Xmat , k=ncol(big.ecorr_Xmat), opts = list(retvec = FALSE))
      evs <- as.numeric(evs$values)
      abs_evs <- abs(evs)
      levs <- length(evs)
    }else{
    ecorr_Xmat = cor(features, use = "pairwise.complete.obs")
    ecorr_Xmat = as.matrix( Matrix::nearPD(ecorr_Xmat)$mat)
    evs = eigen(ecorr_Xmat)$values
    abs_evs <- abs(evs)
    levs <- length(evs)
    }
  }

  # nyholt: effective number of tests (based on Nyholt, 2004)
  m <- 1 + (levs - 1) * (1 - var(evs) / levs)
  Meff_Nyholt <- floor(m)

  # liji: effective number of tests (based on Li & Ji, 2005)
  abs_evs <- abs_evs + sqrt(.Machine$double.eps)
  m <- sum(ifelse(abs_evs >= 1, 1, 0) + (abs_evs - floor(abs_evs)))
  Meff_Liji <- floor(m)

  # gao: effective number of tests (based on Gao, 2008)
  m <- which(cumsum(sort(abs_evs, decreasing = TRUE)) / sum(abs_evs) > 0.995)[1]
  Meff_Gao <- floor(m)

  # galwey: effective number of tests (based on Galwey, 2009)
  evs[evs < 0] <- 0
  m <- sum(sqrt(evs))^2 / sum(evs)
  Meff_Galwey <- floor(m)

  # Sidak
  Sidak <- 1-(1-alpha)^(1/M)
  ENT_Sidak <- alpha/Sidak
  R.percent_Sidak <- (ENT_Sidak/M)*100
  res.ENT_Sidak <- cbind(Sidak, ENT_Sidak, R.percent_Sidak)

  # Bonferroni
  Bonferroni <- alpha/M
  ENT_Bonferroni <- alpha/Bonferroni
  R.percent_Bonferroni <- (ENT_Bonferroni/M)*100
  res.ENT_Bonferroni <- cbind(Bonferroni, ENT_Bonferroni, R.percent_Bonferroni)

  # MWSL
  evs[evs < 0] <- 0
  shape2.star <- (sum(sqrt(evs))/log(evs[1]))^2 / (sum(evs)/evs[1] + sqrt(evs[1]))
  q <- sort(rbeta(K,shape1=1,shape2=floor(shape2.star)))
  MWSL <- q[ceiling(alpha*K)]
  MWSL_CI.up <- q[ceiling(alpha*K) + sqrt(ceiling(alpha*K)*(1-alpha))]
  MWSL_CI.low <- q[ceiling(alpha*K) - sqrt(ceiling(alpha*K)*(1-alpha))]
  ENT_MWSL <- alpha/MWSL
  ENT_MWSL_CI.up = alpha/MWSL_CI.low
  ENT_MWSL_CI.low = alpha/MWSL_CI.up
  R.percent = (ENT_MWSL/M)*100
  res.ENT_MWSL <- cbind(MWSL, MWSL_CI.up, MWSL_CI.low, ENT_MWSL, ENT_MWSL_CI.up, ENT_MWSL_CI.low, R.percent)


  res <- list()
  class(res) = "Meff"
  res$Meff_Nyholt <- Meff_Nyholt
  res$Meff_Liji <- Meff_Liji
  res$Meff_Gao <- Meff_Gao
  res$Meff_Galwey <- Meff_Galwey
  res$Meff_Bonferroni <- ENT_Bonferroni
  res$Meff_Sidak <- ENT_Sidak
  res$Meff_MWSL <- ENT_MWSL
  res$res.Meff_MWSL <- res.ENT_MWSL

  return(res)

}
