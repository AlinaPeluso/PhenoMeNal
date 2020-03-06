
options("scipen"=100, "digits"=4)

devtools::install_github("AlinaPeluso/PhenoMeNal", subdir="MWSL")
library(MWSL)


#==================================================================
#------------------------------------------------------------------
### *** MESA Binned data ***
#------------------------------------------------------------------
#==================================================================

#--------------------------
### Import MESA Binned data
#--------------------------
data("MESA_binned")

#### 4 outcomes
glucose <- MESA_binned[,1]; log_glucose <- MESA_binned[,2]; bmi <- MESA_binned[,7]; log_bmi <- MESA_binned[,8];
outcomes <- MESA_binned[,c(1:4)]

#### 13 confounders: age,sex,height,ethnicityC,ethnicityH,ethnicityAA,smokingF,smokingC,ldl_chol,hdl_chol,sbp,bp_treatment,diabetes,lipids_treatment
MESA_binned$male <- ifelse(MESA_binned$sex<2,1,0)
confounders <- MESA_binned[,c("age","male","height","ethnicityH","ethnicityAA","ethnicityCA","smokingF","smokingC","ldl_chol","hdl_chol","sbp","bp_treatment","diabetes","lipids_treatment")]

#### 655 features
features <- MESA_binned[,23:(ncol(MESA_binned)-1)]

#--------------------------
### Descriptive statistics
#--------------------------
#### features
head(t(round(sapply(features, function(x) c(mean=mean(x),sd=sd(x),median=median(x),min=min(x),max=max(x))),2)),5)
tail(t(round(sapply(features, function(x) c(mean=mean(x),sd=sd(x),median=median(x),min=min(x),max=max(x))),2)),5)

#### confounders
t(round(sapply(MESA_binned[,c(7,9:22,ncol(MESA_binned))], function(x) c(mean=mean(x),sd=sd(x),median=median(x),min=min(x),max=max(x))),3))

#### clinical outcomes
t(round(sapply(outcomes, function(x) c(mean=mean(x),sd=sd(x),median=median(x),min=min(x),max=max(x))),2))

par(mfrow=c(2,4))
for (i in 1:ncol(outcomes)){hist(outcomes[,i],main=names(outcomes)[i],xlab=NULL)}
for (i in 1:ncol(outcomes)){boxplot(outcomes[,i],main=names(outcomes)[i],xlab=NULL)}



#--------------------------------------------
### Permutation-based MWSL and ENT estimation
#--------------------------------------------
methods <- c('identity','mN','mlogN')
mat <- matrix(NA,3,8)
colnames(mat) <- c('MWSL','MWSL_CI.up','MWSL_CI.low','ENT','ENT_CI.up','ENT_CI.low','R.percent','t1err.percent')
rownames(mat) <- methods
rmesa_FWERperm <- list(glucose=mat,log_glucose=mat,bmi=mat,log_bmi=mat)

rmesa_pval <- list(glucose=mat,log_glucose=mat,bmi=mat,log_bmi=mat)
rmesa_FWERperm <- list(glucose=mat,log_glucose=mat,bmi=mat,log_bmi=mat)
allres_mesa <- list()
for (j in 1:length(methods)){
  for (i in 1:ncol(outcomes)){
    rmesa <- FWERperm(outcome=outcomes[,i],
                      features=features,
                      confounders=confounders,
                      n.permutation=10000,
                      method=methods[j],
                      verbose=F)
    allres_mesa[[3*(i-1)+j]] <- rmesa
    rmesa_FWERperm[[i]][j,1:7] <- rmesa$res
    rmesa_FWERperm[[i]][j,8] <- rmesa$t1err.percent
  }
}

df.names <- expand.grid(methods, names(outcomes))
names(allres_mesa) <- paste(df.names$Var1, df.names$Var2,sep='.')
rmesa_FWERperm


# plot of p-values when the features are simulated via multivariate Normal distribution
hist(allres_mesa[["mN.glucose"]][["matPvals"]],main="Plot p-values under the null",breaks=50,xlab=NULL)


# plot of minimum p-values when the features are simulated via multivariate Normal distribution
hist(allres_mesa[["mN.glucose"]][["q"]],main="Plot minimum p-values",breaks=150,xlab=NULL) # ,xlim = c(0,0.035)
op <- par(cex = 1.5); alpha=0.05
abline(v=rmesa_FWERperm$glucose[1,1],col="red",lwd=5)
abline(v=rmesa_FWERperm$glucose[2,1],col="blue",lwd=5)
abline(v=rmesa_FWERperm$glucose[3,1],col="brown",lwd=5)
abline(v=1-(1-alpha)^(1/ncol(features)),col="green",lwd=5)
abline(v=alpha/ncol(features),col="orange",lwd=5)
legend("topright",c('perm_id','perm_mN','perm_mlogN','Sidak','Bonferroni'),fill=c("red","blue","brown","green","orange"))


df_rmesa_FWERperm <- do.call(rbind,rmesa_FWERperm)
df1_rmesa_FWERperm<- data.frame(
  outcome = c('glucose','glucose','glucose',
              'logGlucose','logGlucose','logGlucose',
              'BMI','BMI','BMI',
              'logBMI','logBMI','logBMI'),
  type = c('identity','multivariate Normal','multivariate log-Normal',
           'identity','multivariate Normal','multivariate log-Normal',
           'identity','multivariate Normal','multivariate log-Normal',
           'identity','multivariate Normal','multivariate log-Normal'),
  ENT = c(df_rmesa_FWERperm[,4]),
  ENT.ciUP = c(df_rmesa_FWERperm[,5]),
  ENT.ciLOW = c(df_rmesa_FWERperm[,6]))
(plot_res.MESA_co <- ggplot(data=df1_rmesa_FWERperm,aes(x=outcome,y=ENT)) +
    facet_grid(~ type) +
    geom_hline(yintercept=ncol(features)) +
    annotate("text",x='BMI',y=(ncol(features)+15),label='ANT=655') +
    geom_text(mapping=aes(label=round(ENT,0)),hjust=-.5)+
    geom_point(size=3) +
    geom_errorbar(aes(ymin=ENT.ciLOW,ymax=ENT.ciUP),size=1) +
    theme(text = element_text(size=20)) +
    theme(legend.position="bottom") +
    theme(axis.text.x = element_text(angle=30, hjust=1)) +
    theme(axis.text.y = element_blank()) +
    ggtitle("MESA_binned data - all clinical outcomes")
)



#-------------------------------------------------------------------
### Closed form expression eigenvalues-based MWSL and ENT estimation
#-------------------------------------------------------------------

### Empirical correlation
rmesa_Meff_ecorr <- Meff(features=features,
                         n.permutation=100000,
                         method='scorr',
                         #big.mat=TRUE,
                         alpha=0.05)
mat.rmesa_Meff_ecorr <- rbind(Meff_Nyholt=rmesa_Meff_ecorr$Meff_Nyholt,
                              Meff_Liji=rmesa_Meff_ecorr$Meff_Liji,
                              Meff_Gao=rmesa_Meff_ecorr$Meff_Gao,
                              Meff_Galwey=rmesa_Meff_ecorr$Meff_Galwey,
                              Meff_Bonferroni=rmesa_Meff_ecorr$Meff_Bonferroni,
                              Meff_Sidak=rmesa_Meff_ecorr$Meff_Sidak,
                              Meff_MWSL=rmesa_Meff_ecorr$Meff_MWSL);
colnames(mat.rmesa_Meff_ecorr) <- "Estimate"
mat.rmesa_Meff_ecorr
rmesa_Meff_ecorr$res.Meff_MWSL


### Shrinkage correlation
rmesa_Meff_scorr <- Meff(features=features,
                         n.permutation=100000,
                         method='scorr',
                         #big.mat=TRUE,
                         alpha=0.05)
mat.rmesa_Meff_scorr <- rbind(Meff_Nyholt=rmesa_Meff_scorr$Meff_Nyholt,
                              Meff_Liji=rmesa_Meff_scorr$Meff_Liji,
                              Meff_Gao=rmesa_Meff_scorr$Meff_Gao,
                              Meff_Galwey=rmesa_Meff_scorr$Meff_Galwey,
                              Meff_Bonferroni=rmesa_Meff_scorr$Meff_Bonferroni,
                              Meff_Sidak=rmesa_Meff_scorr$Meff_Sidak,
                              Meff_MWSL=rmesa_Meff_scorr$Meff_MWSL);
colnames(mat.rmesa_Meff_scorr) <- "Estimate"
mat.rmesa_Meff_scorr
rmesa_Meff_scorr$res.Meff_MWSL


#---------------------------------------------------------------------------
### Identification of outcome-specific differentially regulated  metabolites
#---------------------------------------------------------------------------

rmesa_DEtest <- rmesa_DEtest_count  <- list()
names(rmesa_DEtest) <- names(rmesa_DEtest_count) <- c('glucose','logGlucose','BMI','logBMI')
for (i in 1:ncol(outcomes)){
  res_DEtest <- DEtest(outcome=outcomes[,i],
                      features=features,
                      confounders=confounders,
                      MWSL=0.000145102,
                      alpha=0.05,
                      vennPlot=TRUE)
  rmesa_DEtest[[i]]<- res_DEtest
  rmesa_DEtest_count[[i]]<- res_DEtest$res.DE_count
}

res.mesa_DEtest_count <- do.call(rbind, rmesa_DEtest_count)
rownames(res.mesa_DEtest_count) <- c('glucose','logGlucose','BMI','logBMI')
res.mesa_DEtest_count


rmesa_DEtest[["glucose"]][["res.DE_names"]][["FWER.MWSL"]]

