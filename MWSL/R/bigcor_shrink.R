bigcor.shrink = function (x, y = NULL, fun = c("cor", "cov"), size = 2000, verbose = TRUE, ...) 
{
  fun <- match.arg(fun)
  if (fun == "cor") 
    FUN <- cor.shrink
  else FUN <- cov.shrink
  if (fun == "cor") 
    STR <- "Correlation"
  else STR <- "Covariance"
  if (!is.null(y) & NROW(x) != NROW(y)) 
    stop("'x' and 'y' must have compatible dimensions!")
  NCOL <- ncol(x)
  if (!is.null(y)) 
    YCOL <- NCOL(y)
  REST <- NCOL%%size
  LARGE <- NCOL - REST
  NBLOCKS <- NCOL%/%size
  if (is.null(y)) 
    resMAT <- ff::ff(vmode = "double", dim = c(NCOL, NCOL))
  else resMAT <- ff::ff(vmode = "double", dim = c(NCOL, YCOL))
  GROUP <- rep(1:NBLOCKS, each = size)
  if (REST > 0) 
    GROUP <- c(GROUP, rep(NBLOCKS + 1, REST))
  SPLIT <- split(1:NCOL, GROUP)
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)
  if (!is.null(y)) 
    COMBS <- cbind(1:length(SPLIT), rep(1, length(SPLIT)))
  timeINIT <- proc.time()
  for (i in 1:nrow(COMBS)) {
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    if (is.null(y)) {
      if (verbose) 
        message("bigcor: ", sprintf("#%d: %s of Block %s and Block %s (%s x %s) ... ", 
                                    i, STR, COMB[1], COMB[2], length(G1), length(G2)))
      if (all.equal(G1,G2)==T) { 
        RES <- FUN(x[, c(G1)])
        resMAT[G1, G2] <- RES
      } else {
        RES <- FUN(x[, c(G1,G2)])
        dfRES <- matrix(RES[1:size,(size+1):(2*size)],ncol=size)
        resMAT[G1, G2] <- dfRES
        resMAT[G2, G1] <- t(dfRES)
      }
    }
    else {
      if (verbose) 
        message("bigcor: ", sprintf("#%d: %s of Block %s and 'y' (%s x %s) ... ", 
                                    i, STR, COMB[1], length(G1), YCOL))
      RES <- FUN(x[, G1], y, ...)
      resMAT[G1, ] <- RES
    }
    if (verbose) {
      timeNOW <- proc.time() - timeINIT
      message("bigcor: ", round(timeNOW[3], 2), " sec\n")
    }
    gc()
  }
  return(resMAT)
}
