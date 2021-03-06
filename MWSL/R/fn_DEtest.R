#' Identification of differentially regulated metabolomics variates
#'
#'Identification of differentially regulated metabolomics variates directly linked to a specified outcome (e.g. disease risk). This includes the case of continuous outcomes both symmetric and skewed, as well as the case of a binary (0/1) outcome, and a countable outcome taking discrete values 0,1,2,3,..., as well as a survival time-to-event outcome. Once the appropriate analysis has been performed the metabolites of interest are identified comparing the respective p-value to the the adjusted thresold (MWSL or Meff). When the raw p-value corresponding to a certain feature is smaller than the adjusted thresold we identify that metabolomic variate as significant. For the other metabolites we conclude that there is no association between the changes in the outcome variable and the shifts in these features.
#'
#' @param \code{outcome} a vector or a data frame of \code{n} samples. Possible outcome types are: continuous (both symmetric or skewed), discrete binary (0/1), discrete countable, or survival time-to-event
#' @param \code{features} a data frame of \code{n} samples (rows) and \code{M} features (columns)
#' @param \code{confounders} an optional data frame of \code{n} samples (rows) and \code{P} fixed effects confounders (columns), default=\code{NULL}
#' @param \code{MWSL} metabolome-wide significance level (MWSL) estimated via the permutation-based approach (see \code{MWSL::FWERperm}) or via the via the closed-form-expression approach (see \code{MWSL::Meff})
#' @param \code{alpha} an optional probability value, default=0.05
#' @param \code{vennPlot} an optional logic value to be set to \code{TRUE} for the visualisation of the Venn-plot for the total count of differrentally regulated metabolomics variates when employing the Benjamini-Hochberg correction (FDR), the Bonferroni and the MWSL (FWER) correction, respectively
#' @param \code{verbose} an optional logic value to suppress some output status messages, default=\code{TRUE}
#'
#' @return \code{res.DE_count} total count of differrentally regulated metabolomics variates when employing the Benjamini-Hochberg correction (FDR), the Bonferroni and the MWSL (FWER) correction, respectively
#' @return \code{res.DE_names} names of the of differrentally regulated metabolomics variates when employing the Benjamini-Hochberg correction (FDR), the Bonferroni and the MWSL (FWER) correction, respectively
#' @return \code{Venn_plot} Venn-plot visualisation tool for the total count of differrentally regulated metabolomics variates when employing the Benjamini-Hochberg correction (FDR), the Bonferroni and the MWSL (FWER) correction, respectively
#'
#' @export
#'
DEtest <- function(outcome,features,confounders=NULL,MWSL=NULL,alpha=NULL,vennPlot=TRUE,verbose=TRUE){

  Y <- outcome
  n <- nrow(features)
  M <- ncol(features)
  alpha <- ifelse(is.null(alpha)==T,0.05,alpha)

  if (verbose==T){
    print(paste("n.samples =", n))
    print(paste("n.features =", M))}

  # Y ~ survival
  if (length(Y)<n) {
    time <- Y[,1]
    status <- Y[,2]
    model_type <- 'Cox survival regression'
    if (verbose==T){print(paste("model_type =", model_type))}
  } else {
    outcome_type <- ifelse(length(Y) != length(which(Y==as.integer(Y))) | length(table(Y))>180,'continuous','discrete')
    if (verbose==T){print(paste("outcome_type =", outcome_type))}
    # Y ~ continuous
    if (outcome_type=='continuous') {
      at <-  moments::agostino.test(Y)
      model_type <- ifelse(at$p.value<0.001,'median regression','OLS regression')
      if (verbose==T){print(paste("model_type =", model_type))}
    }
    # Y ~ discrete
    if (outcome_type=='discrete') {
      outcome_type.discrete <- ifelse(length(table(Y))>2,'count','binary')
      if (verbose==T){print(paste("outcome_type.discrete =", outcome_type.discrete))}
      model_type <- 'GLM regression'
      if (verbose==T){print(paste("model_type =", model_type))}
      if (outcome_type.discrete=='binary') {family_type <- "binomial"}
      if (outcome_type.discrete=='count') {family_type <- ifelse(var(Y)/mean(Y)>1.5,"negative.binomial","poisson")}
      if (verbose==T){print(paste("family_type =", family_type))}
    }
  }


  if (is.null(confounders)==F){

    cl <- parallel::makeCluster(parallel::detectCores()-1)
    doParallel::registerDoParallel(cl)

    rawp <-
      foreach(m=1:ncol(features), .combine='c') %dopar% {
        new_data <- cbind(Y,features[,m],confounders)
        names(new_data)[1] <- "Y"
        names(new_data)[2] <- "X"

        if (model_type == 'OLS regression') {
         fmla <- as.formula(paste("Y ~ X +", paste(colnames(confounders),collapse="+")))
         reg.out <- lm(fmla,data=as.data.frame(new_data))
         pval <- coef(summary(reg.out))["X","Pr(>|t|)"]}

        if (model_type == 'median regression') {
          fmla <- as.formula(paste("Y ~ X +", paste(colnames(confounders),collapse="+")))
          reg.out <- quantreg::rq(formula=fmla,data=as.data.frame(new_data),tau = .5)
          pval <- coef(summary(reg.out, se = "boot"))["X","Pr(>|t|)"]}

        if (model_type == 'GLM regression') {
          fmla <- as.formula(paste("Y ~ X +", paste(colnames(confounders),collapse="+")))
          if (family_type %in% c('binomial')) {
            reg.out <- glm(fmla,data=as.data.frame(new_data),family="binomial")
            pval <- coef(summary(reg.out))["X","Pr(>|z|)"]}
          if (family_type %in% c('poisson')) {
            reg.out <- speedglm::speedglm(fmla,data=as.data.frame(new_data),family=stats::poisson(link="log"))
            pval <- as.numeric(as.character(coef(summary(reg.out))["X","Pr(>|z|)"]))}
          if (family_type =='negative.binomial'){
            reg.out <- speedglm::speedglm(fmla,data=as.data.frame(new_data),family=MASS::negative.binomial(theta=mean(Y)))
            pval <- as.numeric(as.character(coef(summary(reg.out))["X","Pr(>|t|)"]))}}

        if (model_type == 'Cox survival regression') {
          if (ncol(Y)>2){
            fmla <- as.formula(paste("survival::Surv(time1,time2,status) ~ X +", paste(colnames(confounders),collapse="+")))
          } else {
            fmla <- as.formula(paste("survival::Surv(time,status) ~ X +", paste(colnames(confounders),collapse="+")))}
          reg.out <- survival::coxph(formula=fmla,data=as.data.frame(new_data))
          pval <- summary(reg.out)$coefficients["X","Pr(>|z|)"]}

        pval
      }

    parallel::stopCluster(cl)
  }

  if (is.null(confounders)==T){

    cl <- parallel::makeCluster(parallel::detectCores()-1)
    doParallel::registerDoParallel(cl)

    rawp <-
      foreach(m=1:ncol(features), .combine='c') %dopar% {
        new_data <- cbind(Y,features[,m],confounders)
        names(new_data)[1] <- "Y"
        names(new_data)[2] <- "X"

        if (model_type == 'median reg.') {
          colnames(new_data) <- c('Y','X')
          fmla <- as.formula(paste("Y ~ X"))
          reg.out <- quantreg::rq(formula=fmla,data=as.data.frame(new_data),tau = .5)
          pval <- coef(summary(reg.out, se = "boot"))["X","Pr(>|t|)"]}

        if (model_type == 'OLS') {
          colnames(new_data) <- c('Y','X')
          fmla <- as.formula(paste("Y ~ X"))
          reg.out <- lm(fmla, data=as.data.frame(new_data))
          pval <- coef(summary(reg.out))["X","Pr(>|t|)"]}

        if (model_type == 'GLM') {
          colnames(new_data) <- c('Y','X')
          fmla <- as.formula(paste("Y ~ X"))
          if (family_type %in% c('binomial')) {
            reg.out <- glm(fmla,data=as.data.frame(new_data),family="binomial")
            pval <- coef(summary(reg.out))["X","Pr(>|z|)"]}
          if (family_type %in% c('poisson')) {
            reg.out <- speedglm::speedglm(fmla,data=as.data.frame(new_data),family=stats::poisson(link="log"))
            pval <- as.numeric(as.character(coef(summary(reg.out))["X","Pr(>|z|)"]))}
          if (family_type =='negative.binomial'){
            reg.out <- speedglm::speedglm(fmla,data=as.data.frame(new_data),family=MASS::negative.binomial(theta=mean(Y)))
            pval <- as.numeric(as.character(coef(summary(reg.out))["X","Pr(>|t|)"]))}}

        if (model_type == 'Cox') {
          if (ncol(Y)>2){
            colnames(new_data) <- c('time1','time2','status','X')
            fmla <- as.formula(paste("survival::Surv('time1','time2','status') ~ X"))
          } else {
            colnames(new_data) <- c('time','status','X')
            fmla <- as.formula(paste("survival::Surv(time, status) ~ X"))
          }
          reg.out <- survival::coxph(formula=fmla,data=as.data.frame(new_data))
          pval <- summary(reg.out)$coefficients["X","Pr(>|z|)"]}

        pval
      }

    parallel::stopCluster(cl)
  }

  MWSL_Sidak <- 1-(1-alpha)^(1/M)
  MWSL_Bonferroni <- alpha/M

  res <- as.data.frame(rawp); row.names(res) <- names(features)
  res$FWER.Bonf_de <- ifelse(rawp<MWSL_Bonferroni,1,0)
  res$FWER.Sidak_de <- ifelse(rawp<MWSL_Sidak,1,0)
  res$FWER.MWSL_de <- ifelse(rawp<MWSL,1,0)
  res$FDR.BH_pval <- p.adjust(rawp,method='BH',length(rawp))
  res$FDR.BH_de <- ifelse(res$FDR.BH_pval<alpha,1,0)

  res.DE_count <- cbind(
    FWER.Bonf=sum(res$FWER.Bonf_de),
    #FWER.Sidak=sum(res$FWER.Sidak_de),
    FWER.MWSL=sum(res$FWER.MWSL_de),
    FDR.BH=sum(res$FDR.BH_de))

  res.DE_names <- list(
    FWER.Bonf=row.names(res[res$FWER.Bonf_de>0,]),
    #FWER.Sidak=row.names(res[res$FWER.Sidak_de>0,]),
    FWER.MWSL=row.names(res[res$FWER.MWSL_de>0,]),
    FDR.BH=row.names(res[res$FDR.BH_de>0,]))

  if (vennPlot==TRUE){
    universe <- sort(unique(c(res.DE_names[[1]], res.DE_names[[2]], res.DE_names[[3]])))
    if (length(universe)>0){
      Counts <- matrix(0, nrow=length(universe), ncol=3)
      for (i in 1:length(universe)) {
        Counts[i,1] <- universe[i] %in% res.DE_names[[1]]
        Counts[i,2] <- universe[i] %in% res.DE_names[[2]]
        Counts[i,3] <- universe[i] %in% res.DE_names[[3]]
      }
      colnames(Counts) <- c('FWER.Bonf','FWER.MWSL','FDR.BH')
      cols<-c("Red", "Green", "Blue", "Gray")
      Venn_plot <- limma::vennDiagram(limma::vennCounts(Counts), circle.col=cols)
    }
  }


  res <- list()
  class(res) = "MWSL"
  res$res.DE_count <- res.DE_count
  res$res.DE_names <- res.DE_names
  res$Venn_plot
  return(res)
}
