## Helper function for powMCT.gen SCalc calculates the theoretical covariance
## matrix for a given family and link Note, the theoretical covariance matrices
## for probit and cauchit link rely on transforming the identity link using the
## delta method All derivations are provided in TheoCov.pdf in the package folder

SCalc <- function(doses, dose.samples, n.doses, family, link, mean.vec, modelPar,
                  placEff = NULL, offset = 1) {

  if (family == "negative binomial") {

    if (length(offset) != 1 | class(offset) != "numeric") {
      stop("invalid offset")
    }
    if (!(link %in% c("log", "identity", "risk ratio", "sqrt", "log risk ratio"))) {
      stop("invalid link function")
    }
    if (link == "log") {
      dose.vec <- c()
      resp.vec <- c()
      for (i in 1:n.doses) {
        dose.vec <- c(dose.vec, rep(doses[i], dose.samples[i]))
        resp.vec <- c(resp.vec, rep(mean.vec[i], dose.samples[i]))
      }
      resp.vec <- exp(resp.vec)

      dat <- data.frame(dose = dose.vec, resp = resp.vec)
      dat$dose <- factor(dat$dose)
      x.mat <- model.matrix(resp ~ dose - 1, data = dat)
      S.hat.inv <- matrix(0, nrow = n.doses, ncol = n.doses)
      for (i in 1:nrow(dat)) {
        S.hat.inv <- S.hat.inv + dat$resp[i]/(1 + dat$resp[i]/modelPar) *
          t(t(x.mat[i, ])) %*% t(x.mat[i, ])
      }
      S.hat <- solve(S.hat.inv)/offset

    } else if (link == "identity") {
      dose.vec <- c()
      resp.vec <- c()
      for (i in 1:n.doses) {
        dose.vec <- c(dose.vec, rep(doses[i], dose.samples[i]))
        resp.vec <- c(resp.vec, rep(mean.vec[i], dose.samples[i]))
      }

      dat <- data.frame(dose = dose.vec, resp = resp.vec)
      dat$dose <- factor(dat$dose)
      x.mat <- model.matrix(resp ~ dose - 1, data = dat)
      S.hat.inv <- matrix(0, nrow = n.doses, ncol = n.doses)
      for (i in 1:nrow(dat)) {
        S.hat.inv <- S.hat.inv + (-1/modelPar^2 * (resp.vec[i] + modelPar)/(1 +
                      resp.vec[i]/modelPar)^2 + 1/resp.vec[i]) * t(t(x.mat[i, ])) %*%
          t(x.mat[i, ])
      }
      S.hat <- solve(S.hat.inv)/offset

    } else if (link == "risk ratio") {
      if (is.null(placEff)) {
        stop("must provide placEff for risk ratio link")
      }

      x.mat = diag(n.doses)
      S.hat.inv = matrix(0, nrow = n.doses, ncol = n.doses)

      dose.vec <- c()
      resp.vec <- c()
      mean.vec[1] <- placEff
      rr.tmp <- mean.vec
      mean.vec[-1] <- mean.vec[-1] * placEff
      mean.vec <- log(mean.vec)  # Put it on the log-scale
      for (i in 1:n.doses) {
        dose.vec <- c(dose.vec, rep(doses[i], dose.samples[i]))
        resp.vec <- c(resp.vec, rep(mean.vec[i], dose.samples[i]))
      }
      resp.vec <- exp(resp.vec)


      dat <- data.frame(dose = dose.vec, resp = resp.vec)
      dat$dose <- factor(dat$dose)
      x.mat <- model.matrix(resp ~ dose, data = dat)
      S.hat.inv <- matrix(0, nrow = n.doses, ncol = n.doses)
      for (i in 1:nrow(dat)) {
        S.hat.inv <- S.hat.inv + dat$resp[i]/(1 + dat$resp[i]/modelPar) *
          t(t(x.mat[i, ])) %*% t(x.mat[i, ])
      }
      S.hat <- solve(S.hat.inv)

      ## Now apply the delta method
      S.hat <- diag(rr.tmp) %*% S.hat %*% t(diag(rr.tmp))
      S.hat <- S.hat/offset




    } else if (link == "log risk ratio") {
      if (is.null(placEff)) {
        stop("must provide placEff for log risk ratio link")
      }

      dose.vec <- c()
      resp.vec <- c()
      mean.vec[1] <- placEff
      mean.vec[-1] <- mean.vec[-1] + placEff
      for (i in 1:n.doses) {
        dose.vec <- c(dose.vec, rep(doses[i], dose.samples[i]))
        resp.vec <- c(resp.vec, rep(mean.vec[i], dose.samples[i]))
      }
      resp.vec <- exp(resp.vec)


      dat <- data.frame(dose = dose.vec, resp = resp.vec)
      dat$dose <- factor(dat$dose)
      x.mat <- model.matrix(resp ~ dose, data = dat)
      S.hat.inv <- matrix(0, nrow = n.doses, ncol = n.doses)
      for (i in 1:nrow(dat)) {
        S.hat.inv <- S.hat.inv + dat$resp[i]/(1 + dat$resp[i]/modelPar) *
          t(t(x.mat[i, ])) %*% t(x.mat[i, ])
      }
      S.hat = solve(S.hat.inv)/offset

    } else if (link == "sqrt") {
      dose.vec <- c()
      resp.vec <- c()
      for (i in 1:n.doses) {
        dose.vec <- c(dose.vec, rep(doses[i], dose.samples[i]))
        resp.vec <- c(resp.vec, rep(mean.vec[i], dose.samples[i]))
      }
      resp.vec <- resp.vec^2
      dat <- data.frame(dose = dose.vec, resp = resp.vec)
      dat$dose <- factor(dat$dose)
      x.mat <- model.matrix(resp ~ dose - 1, data = dat)

      S.hat.inv <- matrix(0, nrow = n.doses, ncol = n.doses)
      for (i in 1:nrow(dat)) {
        S.hat.inv <- S.hat.inv + 2 * modelPar * (2 * modelPar * dat$resp[i] +
                            2 * dat$resp[i]^2)/(dat$resp[i] * (modelPar + dat$resp[i])^2) *
          t(t(x.mat[i, ])) %*% t(x.mat[i, ])

      }
      S.hat <- solve(S.hat.inv)/offset

    }

  } else if (family == "binomial") {
    if (link == "logit") {
      resp.vec = plogis(mean.vec)

      dat <- data.frame(dose = doses, resp = resp.vec)
      dat$dose <- factor(dat$dose)
      x.mat <- model.matrix(resp ~ dose - 1, data = dat)
      w.mat <- diag(dose.samples * dat$resp * (1 - dat$resp))
      S.hat <- solve(t(x.mat) %*% w.mat %*% x.mat)/offset


    } else if (link == "probit") {
      mean.vec <- pnorm(mean.vec)
      dat <- data.frame(dose = doses, resp = mean.vec, w = dose.samples)
      dat$dose <- factor(dat$dose)
      x.mat <- model.matrix(resp ~ dose - 1, data = dat)
      S.hat.inv <- matrix(0, nrow = n.doses, ncol = n.doses)
      for (i in 1:n.doses) {
        S.hat.inv <- S.hat.inv + (dat$w[i]/dat$resp[i] + (dat$w[i] - dat$w[i] *
                       dat$resp[i])/(1 - dat$resp[i])^2) * t(t(x.mat[i, ])) %*% t(x.mat[i,])
      }
      S.hat <- solve(S.hat.inv)/offset

      ## Now transform S.hat using Delta method
      S.hat <- diag((1/dnorm(qnorm(mean.vec)))^2) %*% S.hat

    } else if (link == "cauchit") {
      mean.vec <- pcauchy(mean.vec)
      dat <- data.frame(dose = doses, resp = mean.vec, w = dose.samples)
      dat$dose <- factor(dat$dose)
      x.mat <- model.matrix(resp ~ dose - 1, data = dat)
      S.hat.inv <- matrix(0, nrow = n.doses, ncol = n.doses)
      for (i in 1:n.doses) {
        S.hat.inv <- S.hat.inv + (dat$w[i]/dat$resp[i] + (dat$w[i] - dat$w[i] *
                             dat$resp[i])/(1 - dat$resp[i])^2) * t(t(x.mat[i, ])) %*% t(x.mat[i,])
      }
      S.hat <- solve(S.hat.inv)

      ## Now transform S.hat using Delta method
      S.hat <- diag((1/dcauchy(qcauchy(mean.vec)))^2) %*% S.hat/offset

    } else if (link == "cloglog") {
      dat <- data.frame(dose = doses, resp = 1 - exp(-exp(mean.vec)), w = dose.samples)
      dat$dose <- factor(dat$dose)
      x.mat <- model.matrix(resp ~ dose - 1, data = dat)
      S.hat.inv <- matrix(0, nrow = n.doses, ncol = n.doses)
      for (i in 1:n.doses) {
        S.hat.inv <- S.hat.inv + exp(mean.vec[i]) * (-2 * dat$w[i] * exp(exp(mean.vec[i])) +
                       dat$w[i] * exp(2 * exp(mean.vec[i])) + dat$w[i] * dat$resp[i] *
                       exp(exp(mean.vec[i])) - dat$w[i] * dat$resp[i] * exp(2 * exp(mean.vec[i])) +
                       dat$w[i] * dat$resp[i] * exp(exp(mean.vec[i]) + mean.vec[i]) +
                       dat$w[i])/(exp(exp(mean.vec[i])) - 1)^2 * t(t(x.mat[i, ])) %*%
          t(x.mat[i, ])

      }
      S.hat <- solve(S.hat.inv)/offset


    } else if (link == "log") {
      dat <- data.frame(dose = doses, resp = exp(mean.vec), w = dose.samples)
      dat$dose <- factor(dat$dose)
      x.mat <- model.matrix(resp ~ dose - 1, data = dat)
      S.hat.inv <- matrix(0, nrow = n.doses, ncol = n.doses)
      for (i in 1:n.doses) {
        S.hat.inv <- S.hat.inv + ((1 - dat$resp[i]) * (dat$w[i] - dat$w[i] *
                          dat$resp[i]) * dat$resp[i] + (dat$w[i] - dat$w[i] * dat$resp[i]) *
                          dat$resp[i]^2)/(1 - dat$resp[i])^2 * t(t(x.mat[i, ])) %*% t(x.mat[i,])
      }
      S.hat <- solve(S.hat.inv)/offset

    } else if (link == "identity") {

      dat <- data.frame(dose = doses, resp = mean.vec, w = dose.samples)
      dat$dose <- factor(dat$dose)
      x.mat <- model.matrix(resp ~ dose - 1, data = dat)
      S.hat.inv <- matrix(0, nrow = n.doses, ncol = n.doses)
      for (i in 1:n.doses) {
        S.hat.inv <- S.hat.inv + (dat$w[i]/dat$resp[i] + (dat$w[i] - dat$w[i] *
                        dat$resp[i])/(1 - dat$resp[i])^2) * t(t(x.mat[i, ])) %*% t(x.mat[i,])
      }
      S.hat <- solve(S.hat.inv)/offset

    } else if (link == "risk ratio") {

      if (is.null(placEff)) {
        stop("must provide placEff for risk ratio link")
      }
      mean.vec <- log(mean.vec)
      mean.vec[1] <- log(placEff)
      rr.tmp <- exp(mean.vec)
      mean.vec[-1] <- mean.vec[-1] + log(placEff)
      dat <- data.frame(dose = doses, resp = exp(mean.vec), w = dose.samples)
      dat$dose <- factor(dat$dose)
      x.mat <- model.matrix(resp ~ dose, data = dat)
      S.hat.inv <- matrix(0, nrow = n.doses, ncol = n.doses)
      for (i in 1:n.doses) {
        S.hat.inv <- S.hat.inv + ((1 - dat$resp[i]) * (dat$w[i] - dat$w[i] *
                        dat$resp[i]) * dat$resp[i] + (dat$w[i] - dat$w[i] * dat$resp[i]) *
                        dat$resp[i]^2)/(1 - dat$resp[i])^2 * t(t(x.mat[i, ])) %*% t(x.mat[i,])
      }
      S.hat <- solve(S.hat.inv)/offset

      ## Now apply the Delta method
      S.hat <- diag(rr.tmp) %*% S.hat %*% t(diag(rr.tmp))
      S.hat <- S.hat/offset


    } else if (link == "log risk ratio") {
      if (is.null(placEff)) {
        stop("must provide placEff for log risk ratio link")
      }
      mean.vec[1] <- placEff
      mean.vec[-1] <- mean.vec[-1] + placEff
      dat <- data.frame(dose = doses, resp = exp(mean.vec), w = dose.samples)
      dat$dose <- factor(dat$dose)
      x.mat <- model.matrix(resp ~ dose, data = dat)
      S.hat.inv <- matrix(0, nrow = n.doses, ncol = n.doses)
      for (i in 1:n.doses) {
        S.hat.inv <- S.hat.inv + ((1 - dat$resp[i]) * (dat$w[i] - dat$w[i] *
                      dat$resp[i]) * dat$resp[i] + (dat$w[i] - dat$w[i] * dat$resp[i]) *
                      dat$resp[i]^2)/(1 - dat$resp[i])^2 * t(t(x.mat[i, ])) %*% t(x.mat[i,])
      }
      S.hat <- solve(S.hat.inv)/offset


    } else {
      stop("invalid link function")
    }

  } else if (family == "poisson") {
    if (length(offset) != 1 | class(offset) != "numeric") {
      stop("invalid offset")
    }
    if (!(link %in% c("log", "identity", "risk ratio", "sqrt", "log risk ratio"))) {
      stop("invalid link function")
    }

    if (link == "log") {
      dose.vec <- doses
      resp.vec <- exp(mean.vec)

      dat <- data.frame(dose = dose.vec, resp = resp.vec)
      dat$dose <- factor(dat$dose)
      x.mat <- model.matrix(resp ~ dose - 1, data = dat)
      S.hat.inv <- matrix(0, nrow = n.doses, ncol = n.doses)
      for (i in 1:n.doses) {
        S.hat.inv <- S.hat.inv + dose.samples[i] * dat$resp[i] *
          t(t(x.mat[i,])) %*% t(x.mat[i, ])
      }
      S.hat <- solve(S.hat.inv)/offset




    } else if (link == "identity") {
      dose.vec <- doses
      resp.vec <- mean.vec

      dat <- data.frame(dose = dose.vec, resp = resp.vec)
      dat$dose <- factor(dat$dose)
      x.mat <- model.matrix(resp ~ dose - 1, data = dat)
      S.hat.inv <- matrix(0, nrow = n.doses, ncol = n.doses)
      for (i in 1:n.doses) {
        S.hat.inv <- S.hat.inv + dose.samples[i] * t(t(x.mat[i, ])) %*%
          t(x.mat[i,])/dat$resp[i]
      }
      S.hat <- solve(S.hat.inv)/offset


    } else if (link == "sqrt") {
      dose.vec <- doses
      resp.vec <- mean.vec^2

      dat <- data.frame(dose = dose.vec, resp = resp.vec)
      dat$dose <- factor(dat$dose)
      x.mat <- model.matrix(resp ~ dose - 1, data = dat)
      S.hat.inv <- matrix(0, nrow = n.doses, ncol = n.doses)
      for (i in 1:n.doses) {
        S.hat.inv <- S.hat.inv + dose.samples[i] * (4) * t(t(x.mat[i, ])) %*%
          t(x.mat[i, ])
      }
      S.hat <- solve(S.hat.inv)/offset


    } else if (link == "log risk ratio") {

      dose.vec <- doses
      mean.vec[1] <- placEff
      mean.vec[-1] <- mean.vec[-1] + placEff
      resp.vec <- exp(mean.vec)

      dat <- data.frame(dose = dose.vec, resp = resp.vec)
      dat$dose <- factor(dat$dose)
      x.mat <- model.matrix(resp ~ dose, data = dat)
      S.hat.inv <- matrix(0, nrow = n.doses, ncol = n.doses)
      for (i in 1:n.doses) {
        S.hat.inv <- S.hat.inv + dose.samples[i] * dat$resp[i] *
          t(t(x.mat[i,])) %*% t(x.mat[i, ])
      }
      S.hat <- solve(S.hat.inv)/offset



    } else if (link == "risk ratio") {

      dose.vec <- doses
      mean.vec[1] <- placEff
      rr.tmp <- mean.vec
      mean.vec[-1] <- mean.vec[-1] * placEff
      resp.vec <- mean.vec

      dat <- data.frame(dose = dose.vec, resp = resp.vec)
      dat$dose <- factor(dat$dose)
      x.mat <- model.matrix(resp ~ dose, data = dat)
      S.hat.inv <- matrix(0, nrow = n.doses, ncol = n.doses)
      for (i in 1:n.doses) {
        S.hat.inv <- S.hat.inv + dose.samples[i] * dat$resp[i] *
          t(t(x.mat[i,])) %*% t(x.mat[i, ])
      }
      S.hat <- solve(S.hat.inv)/offset
      S.hat <- diag(rr.tmp) %*% S.hat %*% diag(rr.tmp)

    }


  } else {
    stop("invalid family argument")
  }
  return(S.hat)
}




## Slightly modified powMCT function from the DoseFinding package.  I've changed
## the way that muMat get values, by supplying a doses argument This allows me to
## change the doses
powMCTDose <- function(contMat, alpha = 0.025, altModels, n, sigma, S, placAdj = FALSE,
                       alternative = c("one.sided", "two.sided"), df, critV = TRUE,
                       control = mvtnorm.control(), doses = NULL) {
  alternative <- match.arg(alternative)
  if (inherits(contMat, "optContr")) {
    if (attr(contMat, "placAdj") != placAdj) {
      message("using 'placAdj' specification from contMat object")
      placAdj <- attr(contMat, "placAdj")
    }
    contMat <- contMat$contMat
  }
  if (!is.matrix(contMat))
    stop("contMat needs to be a matrix")
  nD <- nrow(contMat)  # nr of doses
  nC <- ncol(contMat)  # nr of contrasts
  ## extract covariance matrix
  if (missing(S)) {
    if (missing(n) | missing(sigma))
      stop("either S or n and sigma need to be specified")
    if (length(n) == 1)
      n <- rep(n, nD)
    if (length(n) != nD)
      stop("n needs to be of length nrow(contMat)")
    S <- sigma^2 * diag(1/n)
    df <- sum(n) - nD
  } else {
    if (!missing(n) | !missing(sigma))
      stop("need to specify exactly one of \"S\" or \"n\" and \"sigma\"")
    if (nrow(S) != ncol(S))
      stop("S needs to be a square matrix")
    if (nrow(S) != nD)
      stop("S needs to have as many rows&cols as there are doses")
    if (missing(df))
      stop("need to specify degrees of freedom in \"df\", when specifying \"S\"")
  }
  ## extract means under the alternative
  if (missing(altModels))
    stop("altModels argument needs to be specified")
  if (!is.null(doses)) {
    muMat <- getResp(altModels, doses)
  } else {
    muMat <- getResp(altModels)
  }
  if (placAdj) {
    muMat <- sweep(muMat, 2, muMat[1, ], "-")  # remove placebo column
    muMat <- muMat[-1, , drop = FALSE]
  }
  if (nrow(muMat) != nD)
    stop("incompatible contMat and muMat")
  ## extract df
  if (missing(S)) {
    if (missing(df))
      stop("degrees of freedom need to be specified in df")
    df <- sum(n) - nD
  }
  ## calculate non-centrality parameter
  deltaMat <- t(contMat) %*% muMat
  covMat <- t(contMat) %*% S %*% contMat
  den <- sqrt(diag(covMat))
  deltaMat <- deltaMat/den
  if (alternative == "two.sided") {
    deltaMat <- abs(deltaMat)
  }
  corMat <- cov2cor(covMat)

  if (!is.finite(df)) {
    df <- 0
  }
  ## calculate critical value
  if (is.logical(critV) & critV == TRUE)
  {
    critV <- critVal(corMat, alpha, df, alternative, control)
  }  # else assume critV already contains critical value
  res <- powCalc(alternative, critV, df, corMat, deltaMat, control)
  ## class(res) <- 'powMCT' attr(res, 'type') <- ifelse(missing(n), 'S', 'n&sigma')
  ## attr(res, 'contMat') <- contMat attr(res, 'muMat') <- muMat
  res
}


## A copy of the powCalc function from DoseFinding found in powMCT.R
powCalc <- function(alternative, critV, df, corMat, deltaMat, control) {
  nC <- nrow(corMat)  # number of contrasts
  if (alternative[1] == "two.sided") {
    lower <- rep(-critV, nC)
  } else {
    lower <- rep(-Inf, nC)
  }
  upper <- rep(critV, nC)
  if (!missing(control)) {
    if (!is.list(control)) {
      stop("when specified, 'control' must be a list")
    }
    ctrl <- do.call("mvtnorm.control", control)
  } else {
    ctrl <- control
  }
  ctrl$interval <- NULL  # not used with pmvt
  nScen <- ncol(deltaMat)
  res <- numeric(nScen)
  for (i in 1:nScen) {
    pmvtCall <- c(list(lower, upper, df = df, corr = corMat,
                       delta = deltaMat[,i], algorithm = ctrl))
    res[i] <- as.vector(1 - do.call("pmvt", pmvtCall))
  }
  names(res) <- colnames(deltaMat)
  res
}

## A copy of the critVal function from DoseFinding found in MCTtest.R
critVal <- function(corMat, alpha = 0.025, df = NULL,
                    alternative = c("one.sided","two.sided"),
                    control = mvtnorm.control()) {
  ## calculate critical value
  alternative <- match.arg(alternative)
  if (missing(corMat))
    stop("corMat needs to be specified")
  if (is.null(df))
    stop("degrees of freedom need to be specified")
  tail <- ifelse(alternative[1] == "two.sided", "both.tails", "lower.tail")
  if (!missing(control)) {
    if (!is.list(control)) {
      stop("when specified, 'control' must be a list")
    }
    ctrl <- do.call("mvtnorm.control", control)
  } else {
    ctrl <- control
  }
  if (!is.finite(df)) {
    # normal case
    df <- 0
  }
  qmvtCall <- c(list(1 - alpha, tail = tail, df = df, corr = corMat, algorithm = ctrl,
                     interval = ctrl$interval))
  do.call("qmvt", qmvtCall)$quantile
}





## Power Calculation at Specific Sample Sizes

#' Calculate Power for Multiple Contrast Test (General Case)
#'
#' Like the powMCT function, this function allows the user to calculate power
#' for a multiple contrast test for a set of specified alternatives for the
#' general case, but without specifying and S matrix. The user supplies the
#' patient allocation, the alternative models, and any parameters needed for the
#' distribution family (e.g. the dispersion parameter for the negative binomial
#' distribution). The function works by calculating the \eqn{\mu} and \eqn{S} for each
#' model in the alternative models and supplying the calculated values to the
#' powMCT function from the DoseFinding package, forwarding relevant arguments.
#' \cr This function also allows a new Ntype, namely 'actual'. If nSample is a
#' vector and Ntype is 'actual', the function interprets nSample to be the exact
#' patient allocation. This is useful for slightly modifying patient allocation
#' and avoiding messy ratios. \cr Furthermore, the function also accepts a
#' theoResp and doses argument, which together describe the theoretical
#' dose-response relationship. The returned power is the probability of
#' accepting at least one of the models specified in altModels.
#'
#'
#' @param nSample An integer if \code{Ntype} is 'arm' or 'total', or a numerical
#'   vector of patient allocations for each arm if \code{Ntype} is 'actual'.
#' @param family A character string containing the error distribution to be used
#'   in the model.
#' @param modelPar A numeric vector containing the additional parameters for the
#'   family argument. If the family is negative binomial, the dispersion
#'   parameter should be supplied. If the family is binomial, no model parameter
#'   should be supplied.
#' @param placEff A numeric value specifying the mean response at the placebo
#'   This is required if \code{link} = 'risk ratio' and ignored otherwise.
#' @param theoResp A numerical vector of theoretical response values, on the
#'   transformed scale (e.g. on the log-scale for the negative binomial family).
#'   This should be the same length as the doses argument.
#' @param doses A numerical vector of doses, corresponding to the theoretical
#'   response values provided.
#' @param Ntype One of 'arm', 'total', or 'actual'. See documentation for
#'   \code{Ntype} in \code{\link[DoseFinding:powMCT]{powMCT}} for descriptions
#'   of 'arm' and 'total'. For 'actual', the nSample should be a numerical
#'   vector containing the actual patient allocation for each dose provided.
#' @param verbose A logical specifying whether the patient allocation should be
#'   printed, in addition to the results.
#' @param offset A positive numeric value specifying the offset term for the
#'   negative binomial distribution. If offset = 1 (the default), then the
#'   offset has no effect. Theoretically, the offset should be a numeric vector
#'   the same length as the number of observations, but for planning purposes,
#'   it is unlikely to know the individual offsets in advance.
#' @param link A character string for the model link function.
#' @param alRatio Vector describing the relative patient allocations to the dose
#'   groups up to proportionality, e.g. \code{rep(1, length(doses))} corresponds
#'   to balanced allocations.
#' @param altModels An object of class \code{Mods}, defining the mean vectors
#'   under which the power should be calculated.
#' @param alpha Significance level to use.
#' @param df Degrees of freedom to assume.
#' @param critV Critical value, if equal to \code{TRUE} the critical value will
#'   be calculated. Otherwise one can directly specify the critical value here.
#' @param alternative Character determining the alternative for the multiple
#'   contrast trend test.
#' @return Numeric containing the calculated power values
#' @examples
#' dose.vec = c(0, 5, 10, 20, 30, 40)
#' models.full = Mods(doses = dose.vec, linear = NULL,
#'       sigEmax = rbind(c(9, 2), c(6, 3)),
#'       emax = 0.8,
#'       quadratic = -0.02,
#'       placEff = 0, maxEff = 2)
#' ## Calculate the power using the responses and doses specified in Mods
#' powMCTGen(30, 'negative binomial', 'log', modelPar = 0.1, Ntype = 'arm',
#'       alpha = 0.05, altModels = models.full)
#' ## Calculate the power at theoretical dose-response values
#' powMCTGen(30, 'negative binomial', 'log', modelPar = 0.1,
#'       theoResp = c(0, 0.01, 0.02, 1, 1.6, 1.8), doses = c(0, 10, 20, 30, 40, 50),
#'       alpha = 0.05, altModels = models.full)
#' @importFrom DoseFinding optContr getResp mvtnorm.control
#' @export
powMCTGen <- function(nSample, family = c("negative binomial", "binomial", "poisson"),
                      link = c("log", "logit", "sqrt", "probit", "cauchit", "cloglog", "identity",
                               "risk ratio", "log risk ratio"), modelPar = NULL, placEff = NULL, theoResp = NULL,
                      doses = NULL, Ntype = c("arm", "total", "actual"), alRatio = NULL, altModels,
                      alpha = 0.025, df = NULL, critV = TRUE, alternative = c("one.sided", "two.sided"),
                      verbose = FALSE, offset = NULL) {

  if (is.null(offset)) {
    offset <- 1
  }

  # Check user-supplied data, format nSample, and handle patient allocation
  if (!is.null(theoResp)) {
    if (length(theoResp) != length(doses)) {
      stop("theoResp and doses must have the same length")
    }
    mean.mat <- matrix(theoResp, ncol = 1)
    n.doses <- length(doses)
  } else {
    if (is.null(doses)) {
      doses <- attr(altModels, "doses")
    }
    n.doses <- length(doses)
    mean.mat <- getResp(altModels, doses)
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
    dose.samples <- round(dose.samples)
  }
  if (verbose) {
    print("the patient allocation is given by:")
    pat.df = data.frame(Dose = doses, nPatients = dose.samples)
    print(pat.df)
  }

  family <- match.arg(family)
  link <- match.arg(link)
  if (link == "risk ratio" & is.null(placEff)) {
    stop("must supply \"placEff\" for link \"risk ratio\"")
  }
  alternative <- match.arg(alternative)

  if (family == "negative binomial" & is.null(modelPar)) {
    stop("must supply modelPar for family \"negative binomial \"")
  }

  ## Need to do this for all alternative models
  n.models <- ncol(mean.mat)
  powVals <- rep(NA, n.models)
  doses.or <- doses

  for (j in 1:n.models) {
    doses <- doses.or
    n.doses <- length(doses)
    if (link == "risk ratio") {
      S.hat <- SCalc(doses, dose.samples, n.doses, family, link, mean.mat[,j],
                     modelPar, placEff, offset)
      S.hat <- S.hat[-1, -1]
      contMat <- optContr(altModels, doses = doses[-1], S = S.hat, placAdj = TRUE)
      n.doses <- n.doses - 1
      if (!is.null(theoResp)) {
        theoResp <- theoResp[-1]
      }
    } else if (link == "log risk ratio") {
      S.hat <- SCalc(doses, dose.samples, n.doses, family, link, mean.mat[,j],
                     modelPar, placEff, offset)
      S.hat <- S.hat[-1, -1]
      contMat <- optContr(altModels, doses = doses[-1], S = S.hat, placAdj = TRUE)
      n.doses <- n.doses - 1
      if (!is.null(theoResp)) {
        theoResp <- theoResp[-1]
      }
    } else {
      S.hat <- SCalc(doses, dose.samples, n.doses, family, link, mean.mat[,j],
                     modelPar, placEff, offset)
      contMat <- optContr(altModels, doses = doses, S = S.hat)
    }

    if (!is.null(theoResp)) {
      if (link == "risk ratio" | link == "log risk ratio") {
        doses <- doses[-1]
      }
      contMat.use = contMat$contMat
      nD <- nrow(contMat.use)
      nC <- ncol(contMat.use)
      muMat <- matrix(theoResp, nrow = n.doses)
      if (link == "risk ratio") {
        muMat <- muMat - 1
      }
      deltaMat <- t(contMat.use) %*% muMat
      covMat <- t(contMat.use) %*% S.hat %*% contMat.use
      den <- sqrt(diag(covMat))
      deltaMat <- deltaMat/den
      if (alternative == "two.sided") {
        deltaMat <- abs(deltaMat)
      }
      corMat <- cov2cor(covMat)
      if (is.null(df)) {
        df <- 0
      } else if (!is.finite(df)) {
        df <- 0
      }
      if (is.logical(critV) & critV == TRUE) {
        critV <- critVal(corMat, alpha, df, alternative = alternative)
      } else if (is.numeric(critV)) {
        critV <- critVal(corMat, alpha, df, alternative = alternative)
      } else {
        stop("critV must either be numeric or \"TRUE\"")
      }
      res <- powCalc(alternative = alternative, critV, df, corMat, deltaMat,
                     control = list())
      powVals <- res
    } else {
      if (is.null(df)) {
        df <- Inf
      }

      if (link == "log risk ratio") {
        powVals[j] <- powMCTDose(contMat, altModels = altModels, df = df,
                                 S = S.hat, placAdj = TRUE, alpha = alpha,
                                 doses = doses)[j]
      } else {
        powVals[j] <- powMCTDose(contMat, altModels = altModels, df = df,
                                 S = S.hat, alpha = alpha, doses = doses)[j]
      }
    }

  }


  if (is.null(theoResp)) {
    names(powVals) <- colnames(mean.mat)
  }

  return(powVals)
}



#' Return the S Matrix for a Theoretical DR-Curve
#'
#' This function is useful for several \code{DoseFinding} functions, but
#' particular for \code{planMod}. Given the true dose-response curve at
#' specified doses, this function will calculate and return the S matrix
#' associated with the specified distribution. If an object of class \code{Mods}
#' is provided in the \code{models} argument, then a list of S matrices will be
#' returned.
#'
#' @param nSample An integer if \code{Ntype} is 'arm' or 'total', or a numerical
#'   vector of patient allocations for each arm if \code{Ntype} is 'actual'.
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
#' @param Ntype One of 'arm', 'total', or 'actual'. See documentation for
#'   \code{Ntype} in \code{\link[DoseFinding:powMCT]{powMCT}} for descriptions
#'   of 'arm' and 'total'. For 'actual', the nSample should be a numerical
#'   vector containing the actual patient allocation for each dose provided.
#' @param alRatio A numeric vector specifying the ratios between the patient
#'   allocation for the specified doses.
#' @param placEff A numeric value of the placebo effect. This is needed only
#'   when the link is risk ratio.
#' @param models Instead of supplying a theoretical dose-response curve and
#'   doses, an object of class \code{Mods} may be provided. The doses will be
#'   pulled from this object, along with the responses at the doses.
#' @param verbose A logical specifying whether the patient allocation should be
#'   printed, in addition to the results.
#' @param offset A positive numeric value specifying the offset term for the
#'   negative binomial distribution. If offset = NULL (the default), then the
#'   offset has no effect. Theoretically, the offset should be a numeric vector
#'   the same length as the number of observations, but for planning purposes,
#'   it is unlikely to know the individual offsets in advance.
#' @return A numeric S matrix, or a list S matrices, for each model is
#'   \code{models}.
#' @examples
#' dose.vec = c(0, 5, 10, 20, 30, 40)
#' models.full = Mods(doses = dose.vec, linear = NULL,
#'       sigEmax = rbind(c(9, 2), c(6, 3)),
#'       emax = 0.8,
#'       quadratic = -0.02,
#'       placEff = 0, maxEff = 2)
#' planModPrepare(30, 'negative binomial', 'log', 0.3, getResp(models.full)[,3],
#'       dose.vec, 'arm')
#' @export
planModPrepare <- function(nSample, family = c("negative binomial", "binomial", "poisson"),
                           link = c("log", "logit", "sqrt", "probit", "cauchit", "cloglog", "identity",
                                    "risk ratio", "log risk ratio"), modelPar = NULL, theoResp = NULL, doses = NULL,
                           Ntype = c("arm", "total", "actual"), alRatio = NULL, placEff = NULL, models = NULL,
                           verbose = FALSE, offset = NULL) {

  if (is.null(offset)) {
    offset <- 1
  }

  if (is.null(models)) {
    if (length(theoResp) != length(doses)) {
      stop("theoResp and doses must have the same length")
    }
    mean.mat <- matrix(theoResp, ncol = 1)
  } else {
    if (class(models) != "Mods") {
      stop("models must be of class Mods if provided")
    }
    mean.mat <- getResp(models)
    doses <- attr(models, "doses")
  }

  n.doses <- length(doses)

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
    message("non-integer patient allocations. Rounding patient allocation")
    dose.samples <- round(dose.samples)
  }

  if (verbose) {
    message("the patient allocation is given by:\n")
    pat.df <- data.frame(Dose = doses, nPatients = dose.samples)
    message(pat.df)
  }

  family <- match.arg(family)
  link <- match.arg(link)
  S.hat <- list()
  for (i in 1:ncol(mean.mat)) {
    S.hat[[i]] <- SCalc(doses, dose.samples, n.doses, family, link,
                        c(mean.mat[,i]), modelPar, placEff, offset)
  }
  if (length(S.hat) > 1) {
    names(S.hat) <- attr(models, "names")
  } else {
    S.hat = S.hat[[1]]
  }

  return(S.hat)
}
