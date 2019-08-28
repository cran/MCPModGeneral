#' Covariance Matrix Simulation
#'
#' For non-canonical links, simulating the covariance matrix is sometimes the
#' easiest way to get an estimate of the covariance matrix. Even for the
#' canonical links, simulating the covariance matrix may be desirable, as
#' theoretical covariance matrices are based off of asymptotic properties which
#' may not hold for small sample sizes.
#'
#' @param doses A numerical vector of doses, corresponding to the theoretical
#'   response values provided.
#' @param resp A numeric vector of response values corresponding to the doses.
#'   This should be on the link scale (e.g. \code{resp} should be on the
#'   log-scale if using the log-link).
#' @param nSample An integer if \code{Ntype} is "arm" or "total", or a numerical
#'   vector of patient allocations for each arm if \code{Ntype} is "actual".
#' @param Ntype One of "arm", "total", or 'actual". See documentation for
#'   \code{Ntype} in \code{\link[DoseFinding:powMCT]{powMCT}} for descriptions
#'   of "arm" and 'total". For "actual", the nSample should be a numerical
#'   vector containing the actual patient allocation for each dose provided.
#' @param nSim An integer specifying the number of simulations used to estimate
#'   the covariance matrix.
#' @param alRatio Vector describing the relative patient allocations to the dose
#'   groups up to proportionality, e.g. \code{rep(1, length(doses))} corresponds
#'   to balanced allocations.
#' @param family A character string containing the error distribution to be used
#'   in the model.
#' @param link A character string specifying the link to be using when modeling
#'   the glm.
#' @param modelPar A numeric vector containing the additional parameters for the
#'   family argument. If the family is negative binomial, the dispersion
#'   parameter should be supplied. If the family is binomial, no model parameter
#'   should be supplied.
#' @param placEff A numeric value specifying the mean response at the placebo
#'   This is required if \code{link} = "risk ratio" or "log risk ratio" and
#'   ignored otherwise.
#' @param verbose A logical specifying whether the patient allocation should be
#'   printed, in addition to the results.
#' @param ... For all other arguments, see the documentation in
#'   \code{\link[DoseFinding:powMCT]{powMCT}}.
#' @return Numeric containing the estimated covariance matrix.
#' @examples
#' data(migraine)
#' models = Mods(linear = NULL, emax = 1, quadratic = c(-0.004), doses = migraine$dose)
#' simS(migraine$dose, getResp(models)[,1], 30, "arm", 10, family = "binomial",
#'      link = "logit", verbose = TRUE)
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
simS <- function(doses, resp, nSample, Ntype = c("arm", "total", "actual"),
                 nSim = 1000, alRatio = NULL, family, link, modelPar = NULL,
                 placEff = NULL, verbose = FALSE) {


  n.doses <- length(doses)

  if (length(doses) != length(resp)) {
    stop("doses must be the same length as resp")
  }
  if (!is.null(alRatio)) {
    if (length(alRatio) != length(doses)) {
      stop("alRatio must be the same length as doses and resp")
    }
  }

  Ntype <- match.arg(Ntype)
  if (Ntype == "actual") {
    dose.samples <- nSample
  } else if (Ntype == "arm" & !is.null(alRatio)) {
    dose.samples <- nSample * alRatio/min(alRatio)
  } else if (Ntype == "total" & !is.null(alRatio)) {
    dose.samples <- nSample * alRatio/sum(alRatio)
  } else if (Ntype == "total" & is.null(alRatio)) {
    dose.samples <- nSample/n.doses * rep(1, n.doses)
  } else {
    dose.samples <- rep(nSample, n.doses)
  }
  if (!all(dose.samples == floor(dose.samples))) {
    dose.samples = round(dose.samples)
  }




  if (family == "negative binomial") {
    if (!(link %in% c("log", "identity", "risk ratio", "sqrt", "log risk ratio"))) {
      stop("invalid link function")
    }
    if (link == "log") {
      resp <- exp(resp)
    } else if (link == "identity") {
      resp <- resp
    } else if (link == "sqrt") {
      resp <- resp^2
    } else if (link == "risk ratio") {
      resp <- resp * placEff
    } else if (link == "log risk ratio") {
      resp <- exp(resp) * placEff
    }

    S.hat.sum <- matrix(0, nrow = n.doses, ncol = n.doses)
    for (j in 1:nSim) {
      dose.vec <- c()
      resp.vec <- c()
      for (i in 1:n.doses) {
        dose.vec <- c(dose.vec, rep(doses[i], dose.samples[i]))
        resp.vec <- c(resp.vec, rnbinom(dose.samples[i], mu = resp[i], size = modelPar))
      }
      dat <- data.frame(dose = dose.vec, resp = resp.vec)

      muS <- prepareGen(family, link, NULL, "dose", "resp", dat)
      S.hat.sum <- S.hat.sum + muS$S

      if (verbose) {
        pb <- txtProgressBar(min = 0, max = 1, initial = 0, style = 3)
        setTxtProgressBar(pb, j/nSim)
      }
    }

  } else if (family == "poisson") {
    if (!(link %in% c("log", "identity", "risk ratio", "sqrt", "log risk ratio"))) {
      stop("invalid link function")
    }
    if (link == "log") {
      resp <- exp(resp)
    } else if (link == "identity") {
      resp <- resp
    } else if (link == "sqrt") {
      resp <- resp^2
    } else if (link == "risk ratio") {
      resp <- resp * placEff
    } else if (link == "log risk ratio") {
      resp <- exp(resp) * placEff
    }

    S.hat.sum <- matrix(0, nrow = n.doses, ncol = n.doses)
    for (j in 1:nSim) {
      dose.vec <- c()
      resp.vec <- c()
      for (i in 1:n.doses) {
        dose.vec <- c(dose.vec, rep(doses[i], dose.samples[i]))
        resp.vec <- c(resp.vec, rpois(dose.samples[i], lambda = resp[i]))
      }

      dat <- data.frame(dose = dose.vec, resp = resp.vec)
      muS <- prepareGen(family, link, NULL, "dose", "resp", dat)
      S.hat.sum <- S.hat.sum + muS$S

      if (verbose) {
        pb <- txtProgressBar(min = 0, max = 1, initial = 0, style = 3)
        setTxtProgressBar(pb, j/nSim)
      }
    }



  } else if (family == "binomial") {
    if (link == "logit") {
      resp.mean <- plogis(resp)
      start <- resp
    } else if (link == "probit") {
      resp.mean <- pnorm(resp)
      start <- resp
    } else if (link == "cauchit") {
      resp.mean <- pcauchy(resp)
      start <- resp
    } else if (link == "cloglog") {
      resp.mean <- 1 - exp(-exp(resp))
      start <- resp
    } else if (link == "log") {
      resp.mean <- exp(resp)
      start <- resp
    } else if (link == "identity") {
      resp.mean <- resp
      start <- resp
    } else if (link == "risk ratio") {
      resp.mean <- resp * placEff
      start <- NA
    } else if (link == "log risk ratio") {
      resp.mean <- exp(resp) * placEff
      start <- NA
    } else {
      stop("invalid link function")
    }

    S.hat.sum <- matrix(0, nrow = n.doses, ncol = n.doses)
    for (j in 1:nSim) {
      dose.vec <- c()
      resp.vec <- c()
      for (i in 1:n.doses) {
        dose.vec <- c(dose.vec, rep(doses[i], dose.samples[i]))
        resp.vec <- c(resp.vec, rbinom(dose.samples[i], size = 1, prob = resp.mean[i]))
      }

      dat <- data.frame(dose = dose.vec, resp = resp.vec)
      muS <- prepareGen(family, link, NULL, "dose", "resp", dat, start = start)
      S.hat.sum <- S.hat.sum + muS$S

      if (verbose) {
        pb <- txtProgressBar(min = 0, max = 1, initial = 0, style = 3)
        setTxtProgressBar(pb, j/nSim)
      }


    }



  }
  return(S.hat.sum/nSim)
}
