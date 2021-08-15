#' @title MI-aware Modifed eBayes Function
#' 
#' @description Modified eBayes function to be used instead of the one in the limma package
#'
#' @param fit an MArrayLM fitted model object produced by lmFit or contrasts.fit. For ebayes only, fit can alternatively be an unclassed list produced by lm.series, gls.series or mrlm containing components coefficients, stdev.unscaled, sigma and df.residual.
#' @param VarRubin a variance-covariance matrix.
#' @param proportion numeric value between 0 and 1, assumed proportion of genes which are differentially expressed
#' @param stdev.coef.lim numeric vector of length 2, assumed lower and upper limits for the standard deviation of log2-fold-changes for differentially expressed genes
#' @param trend logical, should an intensity-trend be allowed for the prior variance? Default is that the prior variance is constant.
#' @param robust logical, should the estimation of df.prior and var.prior be robustified against outlier sample variances?
#' @param winsor.tail.p numeric vector of length 1 or 2, giving left and right tail proportions of x to Winsorize. Used only when robust=TRUE.
#'
#' @return eBayes produces an object of class MArrayLM (see MArrayLM-class) containing everything found in \code{fit} plus the following added components:
#' \describe{ 
#' \item{t}{numeric matrix of moderated t-statistics.}
#' \item{p.value}{numeric matrix of two-sided p-values corresponding to the t-statistics.}
#' \item{lods}{numeric matrix giving the log-odds of differential expression (on the natural log scale).}
#' \item{s2.prior}{estimated prior value for sigma^2. A row-wise vector if covariate is non-NULL, otherwise a single value.}
#' \item{df.prior}{degrees of freedom associated with s2.prior. A row-wise vector if robust=TRUE, otherwise a single value.}
#' \item{df.total}{row-wise numeric vector giving the total degrees of freedom associated with the t-statistics for each gene. Equal to df.prior+df.residual or sum(df.residual), whichever is smaller.}
#' \item{s2.post}{row-wise numeric vector giving the posterior values for sigma^2.}
#' \item{var.prior}{column-wise numeric vector giving estimated prior values for the variance of the log2-fold-changes for differentially expressed gene for each constrast. Used for evaluating lods.}
#' \item{F}{row-wise numeric vector of moderated F-statistics for testing all contrasts defined by the columns of fit simultaneously equal to zero.}
#' \item{F.p.value}{row-wise numeric vector giving p-values corresponding to F.}
#' }
#'
#'The matrices t, p.value and lods have the same dimensions as the input object fit, with rows corresponding to genes and columns to coefficients or contrasts. The vectors s2.prior, df.prior, df.total, F and F.p.value correspond to rows, with length equal to the number of genes. The vector var.prior corresponds to columns, with length equal to the number of contrasts. If s2.prior or df.prior have length 1, then the same value applies to all genes.
#'
#'s2.prior, df.prior and var.prior contain empirical Bayes hyperparameters used to obtain df.total, s2.post and lods.
#' 
#' @author Modified by M. Chion and F. Bertrand. Original by Gordon Smyth and Davis McCarthy
#' @export
#'
#' @examples
#' library(mi4p)
#' data(datasim)
#' datasim_imp <- multi.impute(data = datasim[,-1], conditions = 
#' attr(datasim,"metadata")$Condition, method = "MLE")
#' VarRubin.matrix <- rubin2.all(datasim_imp[1:5,,],
#' attr(datasim,"metadata")$Condition)
#' set.seed(2016)
#' sigma2 <- 0.05 / rchisq(100, df=10) * 10
#' y <- datasim_imp[,,1]
#' design <- cbind(Intercept=1,Group=as.numeric(
#' attr(datasim,"metadata")$Condition)-1)
#' fit.model <- limma::lmFit(y,design)
#' eBayes.mod(fit=fit.model,VarRubin.matrix[[1]])
eBayes.mod<-function (fit, VarRubin, proportion = 0.01, stdev.coef.lim = c(0.1, 4), 
          trend = FALSE, robust = FALSE, winsor.tail.p = c(0.05, 0.1)) 
{
  if (trend) 
    if (is.null(fit$Amean)) 
      stop("Need Amean component in fit to estimate trend")
    eb <- hid.ebayes(fit = fit, VarRubin = VarRubin, proportion = proportion, stdev.coef.lim = stdev.coef.lim, 
                     trend = trend, robust = robust, winsor.tail.p = winsor.tail.p)
    fit$df.prior <- eb$df.prior
    fit$s2.prior <- eb$s2.prior
    fit$var.prior <- eb$var.prior
    fit$proportion <- proportion
    fit$s2.post <- eb$s2.post
    fit$t <- eb$t
    fit$df.total <- eb$df.total
    fit$p.value <- eb$p.value
    fit$lods <- eb$lods
  if (!is.null(fit$design) && limma::is.fullrank(fit$design)) {
    F.stat <- limma::classifyTestsF(fit, fstat.only = TRUE)
    fit$F <- as.vector(F.stat)
    df1 <- attr(F.stat, "df1")
    df2 <- attr(F.stat, "df2")
    if (df2[1] > 1e+06) 
      fit$F.p.value <- pchisq(df1 * fit$F, df1, lower.tail = FALSE)
    else fit$F.p.value <- pf(fit$F, df1, df2, lower.tail = FALSE)
  }
  fit
}
