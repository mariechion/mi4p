# After the eBayes function of the limma package

#' Title
#'
#' @param fit 
#' @param VarRubin 
#' @param proportion 
#' @param stdev.coef.lim 
#' @param trend 
#' @param robust 
#' @param winsor.tail.p 
#'
#' @return
#' @export
#'
#' @examples
#' 1+1
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
