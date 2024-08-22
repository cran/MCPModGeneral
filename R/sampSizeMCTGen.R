#' Sample Size Calculations (General Case)
#'
#' This function build on the \code{sampSizeMCT} function in the
#' \code{DoseFinding} package, allowing the procedure to work with the
#' \code{powMCTGen} function for the general case.
#'
#' @param family A character string containing the error distribution to be used
#'   in the model.
#' @param link A character string for the model link function.
#' @param modelPar A numeric vector containing the additional parameters for the
#'   family argument. If the family is negative binomial, the dispersion
#'   parameter should be supplied. If the family is binomial, no model parameter
#'   should be supplied.
#' @param theoResp A numerical vector of theoretical response values, on the
#'   transformed scale (e.g. on the log-scale for the negative binomial family).
#'   This should be the same length as the doses argument.
#' @param doses A numerical vector of doses, corresponding to the theoretical
#'   response values provided.
#' @param upperN,lowerN Upper and lower bound for the power sample size.
#'   \code{lowerN} defaults to \code{floor(upperN/2)}.
#' @param Ntype One of "arm", "total", or 'actual". See documentation for
#'   \code{Ntype} in \code{\link[DoseFinding:powMCT]{powMCT}} for descriptions
#'   of "arm" and 'total". For "actual", the nSample should be a numerical
#'   vector containing the actual patient allocation for each dose provided.
#' @param alRatio A numeric vector specifying the ratios between the patient
#'   allocation for the specified doses.
#' @param altModels An object of class \code{Mods}, defining the mean vectors
#'   under which the power should be calculated.
#' @param alpha A numeric value specifying the significance level
#' @param power A numeric value specifying the power power of \code{sumFct}
#' @param sumFct Either an included character vector or a function that combines
#'   the power values under the different alternative into one value.
#' @param verbose A logical specifying whether the patient allocation should be
#'   printed, in addition to the results.
#' @param tol A positive numeric value specifying the tolerance level for the
#'   bisection search algorithm. Bisection is stopped if the \code{targFunc}
#'   value is within \code{tol} of power.
#' @return Numeric containing the calculated power values
#' @examples
#' \donttest{
#' dose.vec = c(0, 5, 10, 20, 30, 40)
#' models.full = Mods(doses = dose.vec, linear = NULL,
#'       sigEmax = rbind(c(9, 2), c(6, 3)),
#'       emax = 0.8,
#'       quadratic = -0.02,
#'       placEff = 0, maxEff = 2)
#' ## Now we can calculate the sample sizes needed in order to achieve a certain power
#' sampSizeMCTGen("negative binomial", "log", modelPar = 0.1, upperN = 50, Ntype = "arm",
#'       altModels = models.full, alpha = 0.05, sumFct = "min", power = 0.8)
#' }
#' @export
sampSizeMCTGen <- function(family = c("negative binomial", "binomial", "poisson"),
                           link = c("log", "logit", "probit", "cauchit", "cloglog", "log risk ratio", "risk ratio"),
                           modelPar = NULL, theoResp = NULL, doses = NULL, upperN, lowerN = floor(upperN/2),
                           Ntype = c("arm", "total"), alRatio = NULL, altModels, alpha = 0.025, power = 0.8,
                           sumFct = c("min", "mean", "max"), verbose = FALSE, tol = 0.001) {

  ## This function matches the implementation of sampSizeMCT from the DoseFinding
  ## package, and simple adjusted the input arguments for the general case and the
  ## powMCTGen function

  link <- match.arg(link)
  if (is.null(doses)) {
    n.doses <- length(attr(altModels, "doses"))
  } else {
    n.doses <- length(doses)
  }
  if (is.null(alRatio)) {
    alRatio = rep(1, n.doses)
  }

  sumFct <- match.arg(sumFct)
  if (sumFct == "min") {
    targFunc <- function(n) {
      min(powMCTGen(n, Ntype = Ntype, alRatio = alRatio, altModels = altModels,
                    modelPar = modelPar, theoResp = theoResp, doses = doses, alpha = alpha,
                    family = family, link = link))
    }
  } else if (sumFct == "mean") {
    targFunc <- function(n) {
      mean(powMCTGen(n, Ntype = Ntype, alRatio = alRatio, altModels = altModels,
                     modelPar = modelPar, theoResp = theoResp, doses = doses, alpha = alpha,
                     family = family, link = link))
    }
  } else if (sumFct == "max") {
    targFunc <- function(n) {
      max(powMCTGen(n, Ntype = Ntype, alRatio = alRatio, altModels = altModels,
                    modelPar = modelPar, theoResp = theoResp, doses = doses, alpha = alpha,
                    family = family, link = link))
    }
  } else if (inherits(sumFct, "function")) {
    targFunc <- sumFct
  } else {
    stop("invalid sumFct input")
  }

  upper <- targFunc(upperN) - power
  while (upper < 0) {
    upperN <- upperN * 2
    upper <- targFunc(upperN) - power
  }
  lower <- targFunc(lowerN) - power
  while (lower > 0) {
    lowerN <- lowerN/2
    lower <- targFunc(lowerN) - power
  }

  current <- tol + 1
  niter <- 0
  while (abs(current) > tol & (upperN > lowerN + 1)) {
    currN <- round((upperN + lowerN)/2)
    current <- targFunc(round(currN)) - power
    if (current > 0) {
      upperN <- currN
    } else {
      lowerN <- currN
    }
    niter <- niter + 1
    if (verbose) {
      cat("Iter: ", niter, ", N = ", currN, ", current value = ",
          round(current + power, 4), "\n", sep = "")
    }
  }
  while (current < 0) {
    currN <- currN + 1
    current <- targFunc(round(currN)) - power
  }
  if (Ntype == "total") {
    arm.sizes <- round(currN * alRatio/sum(alRatio))
  } else {
    arm.sizes <- round(currN * alRatio/min(alRatio))
  }
  res <- list(samp.size = arm.sizes, target = round(current + power, 4))
  attr(res, "alRatio") <- round(alRatio/min(alRatio), 4)
  attr(res, "power") <- power
  attr(res, "Ntype") <- Ntype
  class(res) <- "sampSize"
  # message("\nUsing \"", sumFct, "\" of power\n", sep = "")
  res

}
