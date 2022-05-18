#' @title Multiple imputation of quantitative proteomics datasets
#' 
#' @description \code{\link{multi.impute}} performs multiple imputation on a
#' given quantitative proteomics dataset.
#'
#' @param data A quantitative matrix to be imputed, with proteins/peptides in
#' rows and samples in columns.
#' @param conditions A vector of length the number of samples where each element
#' corresponds to the condition the sample belongs to.
#' @param nb.imp The number of imputation to perform.
#' @param method A single character string describing the imputation method to 
#' be used. See details.
#' @param parallel Logical, whether or not use parallel computing
#' (with \code{\link[foreach]{foreach}}).
#'
#' @return A numeric array of dimension c(dim(data),nb.imp).
#' 
#' @details Multiple imputation consists in imputing several times a given
#' dataset using a given method. Here, imputation methods can be chosen either 
#' from \code{\link[mice]{mice}}, \code{\link[imp4p]{imp4p-package}} or 
#' \code{\link[impute]{impute.knn}}:
#' \itemize{
#'     \item "pmm", "midastouch", "sample", "cart", "rf","mean", "norm", 
#'     "norm.nob", "norm.boot", "norm.predict": imputation methods as described 
#'     in \code{\link[mice]{mice}}.
#'     \item "RF" imputes missing values using random forests algorithm as 
#'     described in \code{\link[imp4p]{impute.RF}}.
#'     \item "MLE" imputes missing values using maximum likelihood estimation
#'     as described in \code{\link[imp4p]{impute.mle}}.
#'     \item "PCA" imputes missing values using principal component analysis as 
#'     described in \code{\link[imp4p]{impute.PCA}}.
#'     \item "SLSA" imputes missing values using structured least squares 
#'     algorithm as described in \code{\link[imp4p]{impute.slsa}}.
#'     \item "kNN" imputes missing values using k nearest neighbors as 
#'     described in \code{\link[impute]{impute.knn}}.
#' }
#' 
#' 
#' @references 
#' M. Chion, Ch. Carapito and F. Bertrand (2021). \emph{Accounting for multiple 
#' imputation-induced variability for differential analysis in mass 
#' spectrometry-based label-free quantitative proteomics}. arxiv:2108.07086. 
#' \url{https://arxiv.org/abs/2108.07086}.
#' 
#' @export
#'
#' @examples
#' library(mi4p)
#' data(datasim)
#' multi.impute(data = datasim[,-1], conditions = attr(datasim,"metadata")$Condition, method = "MLE")
multi.impute <- function(data, conditions, nb.imp = NULL, method, parallel = FALSE){
  conditions <- as.factor(conditions)
  # Determine the number of imputations to be done
  # Rule of thumb (White et al. 2011) : I = % MV
  if (is.null(nb.imp)) {
    nb.imp <- ceiling(sum(is.na(data))/prod(dim(data))*100)
  }
  if (nb.imp<2) {nb.imp <- 2}
  # Imputation method from mice package
  if (method %in% c("pmm", "midastouch", "sample", "cart", "rf",
                        "mean", "norm", "norm.nob", "norm.boot", "norm.predict")) {
    data.imp <- mice::mice(data, m = nb.imp, method = method, printFlag = FALSE)
    data.imp.list <- mice::complete(data.imp, action = "all")
    data.array <- array(unlist(data.imp.list), 
                        dim = c(nrow(data.imp.list[[1]]), ncol(data.imp.list[[1]]), length(data.imp.list)))
  }
  # Imputation method from imp4p package
  if (method == "RF") {
    if (parallel) {
      iforeach<-NA
      requireNamespace("foreach",quietly = TRUE)
      data.array <- simplify2array(foreach::foreach(iforeach =  1:nb.imp,
                                           .errorhandling = 'remove', .verbose = TRUE) %dopar% 
                                     imp4p::impute.RF(tab = data, conditions = conditions))
    }
    else {
      data.array <- simplify2array(lapply(1:nb.imp, function(aaa) {
        imp4p::impute.RF(tab = data, conditions = conditions)}
      ))
    }
  }
  if (method == "MLE") {
    if (parallel) {
      iforeach<-NA
      requireNamespace("foreach",quietly = TRUE)
      data.array <- simplify2array(foreach::foreach(iforeach =  1:nb.imp,
                                           .errorhandling = 'remove', .verbose = TRUE) %dopar% 
                                     imp4p::impute.mle(tab = data, conditions = conditions))
    }
    else {
      data.array <- simplify2array(lapply(1:nb.imp, function(aaa) {
        imp4p::impute.mle(tab = data, conditions = conditions)}
      ))
    }
  }
  if (method == "PCA") {
    if (parallel) {
      iforeach<-NA
      requireNamespace("foreach",quietly = TRUE)
      data.array <- simplify2array(foreach::foreach(iforeach =  1:nb.imp,
                                           .errorhandling = 'remove', .verbose = TRUE) %dopar% 
                                     imp4p::impute.PCA(tab = data, conditions = conditions))
    }
    else {
      data.array <- simplify2array(lapply(1:nb.imp, function(aaa) {
        imp4p::impute.PCA(tab = data, conditions = conditions)}
      ))
    }
  }
  if (method == "SLSA") {
    if (parallel) {
      iforeach<-NA
      requireNamespace("foreach",quietly = TRUE)
      data.array <- simplify2array(foreach::foreach(iforeach =  1:nb.imp,
                                           .errorhandling = 'remove', .verbose = TRUE) %dopar% 
                                     imp4p::impute.slsa(tab = data, conditions = conditions))
    }
    else {
      data.array <- simplify2array(lapply(1:nb.imp, function(aaa) {
        imp4p::impute.slsa(tab = data, conditions = conditions)}
      ))
    }
  }
  # Imputation using kNN
  if (method == "kNN") {
    if (parallel) {
      iforeach<-NA
      requireNamespace("foreach",quietly = TRUE)
      data.array <- simplify2array(foreach::foreach(iforeach =  1:nb.imp,
                                           .errorhandling = 'remove', .verbose = TRUE) %dopar% 
                                     impute::impute.knn(data = as.matrix(data))$data)
    }
    else {
      data.array <- simplify2array(lapply(1:nb.imp, function(aaa) {
        impute::impute.knn(data = as.matrix(data))$data}
      ))
    }
  }
  # Returrn arrray of imputed matrices
  return(data.array)
}
