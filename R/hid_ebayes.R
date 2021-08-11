# Adapted from the .ebayes function of the limma package

#' Title
#'
#' @param fit 
#' @param VarRubin 
#' @param mod 
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
hid.ebayes<-function (fit, VarRubin, mod=T, proportion = 0.01, stdev.coef.lim = c(0.1, 4), 
          trend = FALSE, robust = FALSE, winsor.tail.p = c(0.05, 0.1)) 
{
  coefficients <- fit$coefficients
  stdev.unscaled <- fit$stdev.unscaled
  sigma <- VarRubin #fit$sigma
  df.residual <- fit$df.residual
  if (is.null(coefficients) || is.null(stdev.unscaled) || is.null(sigma) || 
      is.null(df.residual)) 
    stop("No data, or argument is not a valid lmFit object")
  if (all(df.residual == 0)) 
    stop("No residual degrees of freedom in linear model fits")
  if (all(!is.finite(sigma))) 
    stop("No finite residual standard deviations")
  if (trend) {
    covariate <- fit$Amean
    if (is.null(covariate)) 
      stop("Need Amean component in fit to estimate trend")
  }
  else {
    covariate <- NULL
  }
  if (mod) {
    out <- limma::squeezeVar(sigma^2, df.residual, covariate = covariate, 
                    robust = robust, winsor.tail.p = winsor.tail.p)
    out$s2.prior <- out$var.prior
    out$s2.post <- out$var.post
    out$var.prior <- out$var.post <- NULL
    out$t <- coefficients/stdev.unscaled/sqrt(out$s2.post)
    df.total <- df.residual + out$df.prior
    df.pooled <- sum(df.residual, na.rm = TRUE)
    df.total <- pmin(df.total, df.pooled)
    out$df.total <- df.total
    out$p.value <- 2 * pt(-abs(out$t), df = df.total)
    var.prior.lim <- stdev.coef.lim^2/median(out$s2.prior)
    out$var.prior <- limma::tmixture.matrix(out$t, stdev.unscaled, df.total, 
                                   proportion, var.prior.lim)
  }
  else {
    #A completer
  }
  
  if (any(is.na(out$var.prior))) {
    out$var.prior[is.na(out$var.prior)] <- 1/out$s2.prior
    warning("Estimation of var.prior failed - set to default value")
  }
  r <- rep(1, NROW(out$t)) %o% out$var.prior
  r <- (stdev.unscaled^2 + r)/stdev.unscaled^2
  t2 <- out$t^2
  Infdf <- out$df.prior > 10^6
  if (any(Infdf)) {
    kernel <- t2 * (1 - 1/r)/2
    if (any(!Infdf)) {
      t2.f <- t2[!Infdf]
      r.f <- r[!Infdf]
      df.total.f <- df.total[!Infdf]
      kernel[!Infdf] <- (1 + df.total.f)/2 * log((t2.f + 
                                                    df.total.f)/(t2.f/r.f + df.total.f))
    }
  }
  else kernel <- (1 + df.total)/2 * log((t2 + df.total)/(t2/r + 
                                                           df.total))
  out$lods <- log(proportion/(1 - proportion)) - log(r)/2 + 
    kernel
  out
}


