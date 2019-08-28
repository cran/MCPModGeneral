#' Prepare General Data for the MCPMod Function
#'
#' This function serves as an alternative for using the MCPModGen function
#' directly for general data. The function returns the estimates for \eqn{\mu} and
#' \eqn{S}, which are needed for MCPMod.
#'
#' @param family A character string containing the error distribution to be used
#'   in the model.
#' @param link A character string for the model link function.
#' @param w Either a numeric vector of the same length as dose and resp, or a
#'   character vector denoting the column name in the data.
#' @param offset Either a numeric vector of the same length as dose and resp, or a
#'   character vector denoting the column name in the data.
#' @param ... Additional arguments to be passed to \code{glm} or \code{glm.nb}.
#'   This is especially useful when a fitting error is returned. In these cases,
#'   it may be useful to supply a \code{start} vector for the parameters.
#' @inheritParams MCPModGen
#' @return A list containing the \eqn{\mu} vector and \eqn{S} matrix.
#' @examples
#' # Analyze the binary migraine data from the DoseFinding package.
#' data(migraine)
#' models = Mods(linear = NULL, emax = 1, quadratic = c(-0.004), doses = migraine$dose)
#'
#' # Now analyze using binomial weights
#' PFrate <- migraine$painfree/migraine$ntrt
#' migraine$pfrat = migraine$painfree / migraine$ntrt
#' muS = prepareGen("binomial", "logit", w = "ntrt", dose = "dose",
#'                  resp = "pfrat", data = migraine)
#' ## Look at the elements of muS
#' muS
#' MCPMod(muS$data$dose, muS$data$resp, models = models, S = muS$S,
#'        type = "general", selModel = "aveAIC",Delta = 0.2)
#' @export
prepareGen <- function(family = c("negative binomial", "binomial", "poisson"),
                       link = c("log", "logit", "probit", "cauchit", "cloglog",
                                "identity", "log risk ratio", "risk ratio",
                                "sqrt"), w = NULL, dose, resp, data = NULL,
                       addCovars = ~1, placAdj = FALSE,
                       offset = NULL, ...) {

  ## Format user data
  if (!is.null(data)) {
    if (class(data) != "data.frame") {
      stop("data must be of class \"data.frame\"")
    }
    if (class(dose) != "character" | class(resp) != "character") {
      stop("dose and resp must be character vectors")
    }
    dose.vec <- data[[dose]]
    resp.vec <- data[[resp]]
    if (!is.null(w)) {
      if (class(w) != "character") {
        stop("w must be a character vector")
      }
      w.vec <- data[[w]]
    }
  } else {
    if (length(dose) != length(resp)) {
      stop("dose and resp must be of equal length")
    }
    dose.vec <- dose
    resp.vec <- resp
    if (!is.null(w)) {
      if (length(w) != length(resp)) {
        stop("w must be of equal length as dose and resp")
      }
      w.vec <- w
    }
    if (addCovars != ~1) {
      stop("must supply data when addCovars is given")
    }
  }

  if (addCovars == ~1) {
    if (!is.null(w)) {
      dat.glm <- data.frame(dose = dose.vec, resp = resp.vec, w = w.vec)
    } else {
      dat.glm <- data.frame(dose = dose.vec, resp = resp.vec)
    }
  } else {
    dat.glm <- data
    ## Need to rename dose and resp columns
    which.dose <- which(names(dat.glm) == dose)
    which.resp <- which(names(dat.glm) == resp)
    names(dat.glm)[which.dose] = "dose"
    names(dat.glm)[which.resp] = "resp"
  }


  if (is.null(offset)) {
    ## Do nothing
  } else if (!is.null(data) & class(offset) == "character") {
    offset <- data[[offset]]
  } else if (length(offset) != nrow(dat.glm)) {
    stop("offset must be the same length as data")
  }


  doses <- sort(unique(dat.glm$dose))
  n.doses <- length(doses)

  family <- match.arg(family)
  link.in <- match.arg(link)

  ## Negative binomial, binomial, and Poisson are similar, so comments are only
  ## included for the negative binomial chunk
  if (family == "negative binomial") {
    if (!is.null(w)) {
      message("ignoring w argument for negative binomial")
    }

    ## For log risk ratio and risk ratio, we use a log link We create the formula then
    ## call MASS::glm.nb. Log risk ratio is identical to placebo-adjusted
    ## on log-scale, so fit that addCovars[2] gets the arguments part of the addCovars
    if (link.in == "log risk ratio" | link.in == "risk ratio") {
      if (is.null(offset)) {
        ## If there is no offset
        form <- as.formula(paste("resp~factor(dose)+", addCovars[2], sep = ""))
      } else {
        form <- as.formula(paste("resp~factor(dose)+", addCovars[2],
                                 "+ offset(log(offset))", sep = ""))
      }
      # call.glm <- list(form, dat.glm, "log")
      # names(call.glm) <- c("formula", "data", "link")
      # glm.fit <- do.call(MASS::glm.nb, call.glm)
      glm.fit = MASS::glm.nb(form, data = dat.glm, link = "log", ...)

      mu.hat <- coef(glm.fit)
      S.hat <- vcov(glm.fit)
      if (link == "risk ratio") {
        ## Convert to risk-ratio using Delta method
        mu.hat <- exp(mu.hat)
        S.hat <- diag(mu.hat) %*% S.hat %*% diag(mu.hat)
      }
      mu.hat <- mu.hat[1:n.doses]
      S.hat <- S.hat[1:n.doses, 1:n.doses]
    } else {
      if (placAdj == FALSE) {
        if (is.null(offset)) {
          form <- as.formula(paste("resp~factor(dose)+", addCovars[2], "-1",
                                   sep = ""))
        } else {
          form <- as.formula(paste("resp~factor(dose)+", addCovars[2],
                                   "-1 + offset(log(offset))", sep = ""))
        }
        call.glm <- list(form, dat.glm, link.in)
        names(call.glm) <- c("formula", "data", "link")
        glm.fit <- do.call(MASS::glm.nb, call.glm)
        ## I would like to use the MASS::glm.nb version so I can pass the ...
        ## argument, or find a way to pass ... through do.call. However,
        ## MASS::glm.nb fails to pass link.in inside the Poisson call
        # glm.fit <- MASS::glm.nb(form, data = dat.glm, link = link.in, ...)

        mu.hat <- coef(glm.fit)[1:n.doses]
        S.hat <- vcov(glm.fit)[1:n.doses, 1:n.doses]
      } else if (placAdj) {
        if (is.null(offset)) {
          form <- as.formula(paste("resp~factor(dose)+", addCovars[2], sep = ""))
        } else {
          form <- as.formula(paste("resp~factor(dose)+", addCovars[2],
                                   "+ offset(log(offset))", sep = ""))
        }
        call.glm <- list(form, dat.glm, link.in)
        names(call.glm) <- c("formula", "data", "link")
        glm.fit <- do.call(MASS::glm.nb, call.glm)
        # glm.fit <- MASS::glm.nb(form, data = dat.glm, link = link.in, ...)
        mu.hat <- coef(glm.fit)[2:n.doses]
        S.hat <- vcov(glm.fit)[2:n.doses, 2:n.doses]
        doses <- doses[2:n.doses]
      } else {
        stop("placAdj must be logical")
      }
    }

  } else if (family == "binomial") {


    if (link.in == "log risk ratio" | link.in == "risk ratio") {
      if (is.null(w)) {
        risks.vec = rep(NA, n.doses)
        k <- 1
        for (i in doses) {
          risks.vec[k] <- mean(dat.glm$resp[dat.glm$dose == i])
          k <- k + 1
        }
        emp.rr <- risks.vec/risks.vec[1]
        if (length(labels(terms(addCovars))) > 0) {
          emp.rr <- c(emp.rr, rep(1, length(labels(terms(addCovars)))))
        }

        ## Calculate log-risk ratio glm.fit <- glm(resp ~ factor(dose), data = dat.glm,
        ## family = binomial(link = 'log'), start = log(c(mean(dat.glm$resp[dat.glm$dose
        ## == 0]), emp.rr[-1])))
        form <- as.formula(paste("resp~factor(dose)+", addCovars[2], sep = ""))
        # call.glm <- list(form, dat.glm, binomial(link = "log"),
        #                  log(c(mean(dat.glm$resp[dat.glm$dose == 0]), emp.rr[-1])))
        # names(call.glm) <- c("formula", "data", "family", "start")
        # glm.fit <- do.call(glm, call.glm)
        glm.fit <- glm(form, data = dat.glm, family = binomial(link = "log"),
                       start = log(c(mean(dat.glm$resp[dat.glm$dose == 0]),
                                     emp.rr[-1])))

        mu.hat <- coef(glm.fit)
        S.hat <- vcov(glm.fit)
        if (link.in == "risk ratio") {
          ## Convert to risk-ratio using Delta method
          mu.hat <- exp(mu.hat)
          S.hat <- diag(mu.hat) %*% S.hat %*% t(diag(mu.hat))
        }
        mu.hat <- mu.hat[1:n.doses]
        S.hat <- S.hat[1:n.doses, 1:n.doses]
      } else {
        risks.vec <- rep(NA, n.doses)
        k <- 1
        for (i in doses) {
          risks.vec[k] <- dat.glm$resp[dat.glm$dose == i]
          k <- k + 1
          ## there is a faster way to do this, but this ensures the order otherwise use
          ## emp.rr = dat.glm$resp/dat.glm$resp[dat.glm$dose == i]
        }
        emp.rr <- risks.vec/risks.vec[1]
        if (length(labels(terms(addCovars))) > 0) {
          emp.rr <- c(emp.rr, rep(1, length(labels(terms(addCovars)))))
        }

        # glm.fit <- glm(resp ~ factor(dose), data = dat.glm, family = binomial(link =
        # 'log'), weights = w, start = log(c(dat.glm$resp[dat.glm$dose == 0],
        # emp.rr[-1])))
        form <- as.formula(paste("resp~factor(dose)+", addCovars[2], sep = ""))
        # call.glm <- list(form, dat.glm, binomial(link = "log"),
        #                  log(c(dat.glm$resp[dat.glm$dose == 0], emp.rr[-1])), dat.glm$w)
        # names(call.glm) <- c("formula", "data", "family", "start", "weights")
        # glm.fit <- do.call(glm, call.glm)
        glm.fit <- glm(form, data = dat.glm, family = binomial(link = "log"),
                       weights = dat.glm$w,
                       start = log(c(dat.glm$resp[dat.glm$dose == 0], emp.rr[-1])))

        mu.hat <- coef(glm.fit)
        S.hat <- vcov(glm.fit)
        if (link.in == "risk ratio") {
          ## Convert to risk-ratio using Delta method
          mu.hat <- exp(mu.hat)
          S.hat <- diag(mu.hat) %*% S.hat %*% t(diag(mu.hat))
        }
        mu.hat <- mu.hat[1:n.doses]
        S.hat <- S.hat[1:n.doses, 1:n.doses]
      }
    } else {
      if (placAdj == FALSE) {
        form <- as.formula(paste("resp~factor(dose)+", addCovars[2], "-1",
                                 sep = ""))
        glm.fit <- glm(form, data = dat.glm, weights = w,
                       family = binomial(link = link.in), ...)
        mu.hat <- coef(glm.fit)[1:n.doses]
        S.hat <- vcov(glm.fit)[1:n.doses, 1:n.doses]

      } else if (placAdj) {
        form <- as.formula(paste("resp~factor(dose)+", addCovars[2], sep = ""))
        glm.fit <- glm(form, data = dat.glm, weights = w,
                       family = binomial(link = link.in), ...)
        mu.hat <- coef(glm.fit)[2:n.doses]
        S.hat <- vcov(glm.fit)[2:n.doses, 2:n.doses]
        doses <- doses[2:n.doses]
      } else {
        stop("placAdj must be logical")
      }
    }

  } else if (family == "poisson") {
    if (!is.null(w)) {
      message("ignoring w argument for Poisson")
    }

    if (link.in == "log risk ratio" | link.in == "risk ratio") {
      if (is.null(offset)) {
        form <- as.formula(paste("resp~factor(dose)+", addCovars[2], sep = ""))
      } else {
        form <- as.formula(paste("resp~factor(dose)+", addCovars[2],
                                 " + offset(log(offset))", sep = ""))
      }
      glm.fit <- glm(form, data = dat.glm, family = poisson(link = "log"),
                     ...)

      mu.hat <- coef(glm.fit)
      S.hat <- vcov(glm.fit)
      if (link.in == "risk ratio") {
        ## Convert to risk-ratio using Delta method
        mu.hat <- exp(mu.hat)
        S.hat <- diag(mu.hat) %*% S.hat %*% diag(mu.hat)
      }
      mu.hat <- mu.hat[1:n.doses]
      S.hat <- S.hat[1:n.doses, 1:n.doses]
    } else {
      if (placAdj == FALSE) {
        if (is.null(offset)) {
          form <- as.formula(paste("resp~factor(dose)+", addCovars[2], "-1",
                                   sep = ""))
        } else {
          form <- as.formula(paste("resp~factor(dose)+", addCovars[2],
                                   "-1 + offset(log(offset))", sep = ""))
        }
        glm.fit <- glm(form, data = dat.glm, family = poisson(link = link.in),
                       ...)
        mu.hat <- coef(glm.fit)[1:n.doses]
        S.hat <- vcov(glm.fit)[1:n.doses, 1:n.doses]
      } else if (placAdj) {
        form <- as.formula(paste("resp~factor(dose)+", addCovars[2], sep = ""))
        glm.fit <- glm(form, data = dat.glm, family = poisson(link = link.in),
                       ...)
        mu.hat <- coef(glm.fit)[2:n.doses]
        S.hat <- vcov(glm.fit)[2:n.doses, 2:n.doses]
        doses <- doses[2:n.doses]
      } else {
        stop("placAdj must be logical")
      }
    }



  } else {
    stop("invalid family argument")
  }

  data.df <- data.frame(dose = doses, resp = mu.hat)
  return.list <- list(data = data.df, S = S.hat)

  return(return.list)
}
