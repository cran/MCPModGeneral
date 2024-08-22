#' MCPModSurv - Multiple Comparison and Modeling for Coxph and Parametric
#' Survival Models
#'
#' This function allows the user to implement the MCPMod function on a Cox
#' proportional hazards regression model and a parametric survival model. The
#' function works very similarly to
#' \code{\link[MCPModGeneral:MCPModGen]{MCPModGen}}, but is unique enough in
#' terms of the data and the parameters to warrant its own function.
#'
#' `MCPModSurv` works by making calls to `coxph`, `survreg`, and `Surv` from the
#' `survival` package. After retrieving coefficient estimates and the estimated
#' covariance matrix, values are passed into the `MCPMod` function from the
#' `DoseFinding` package.
#'
#' @param model A character string containing the survival regression model.
#' @param dist A character string for the distribution, in the case when
#'   \code{model} is "parametric". Must be one of "\code{weibull}",
#'   "\code{exponential}", "\code{gaussian}", "\code{logistic}",
#'   "\code{lognormal}", or "\code{loglogistic}".
#' @param returnS Logical determining whether muHat and SHat should be returned,
#'   in additional to the MCPMod output.
#' @param dose,resp,status Either character strings specifying the names of the
#'   respective columns in the \code{data} data frame, or numeric vectors of
#'   equal length containing their respective values. \code{status} refers to
#'   whether an observation was censored or not. If no observations were
#'   censored, \code{status} should be a vector of 1s.
#' @inheritParams MCPModGen
#' @param ... Additional arguments to be passed to \code{coxph} or
#'   \code{survreg}. This is especially useful when a fitting error is returned.
#' @return An object of class MCPMod if returnS = FALSE. Otherwise, a list
#'   containing an object of class MCPMod, the numeric vector \eqn{\mu}, and the
#'   numeric matrix \eqn{S}.
#' @export
MCPModSurv <- function(model = c("coxph", "parametric"), dist = NULL,
                       returnS = FALSE, dose, resp, status, data = NULL, models,
                       placAdj = FALSE, selModel = c("AIC", "maxT", "aveAIC"),
                       alpha = 0.025, df = NULL, critV = NULL,
                       doseType = c("TD", "ED"), Delta, p, pVal = TRUE,
                       alternative = c("one.sided", "two.sided"),
                       na.action = na.fail, mvtcontrol = mvtnorm.control(),
                       bnds, control = NULL, ...) {

  ## Check and format the incoming data.
  if (!is.null(data)) {
    if (!inherits(data, "data.frame")) {
      stop("data must be of class \"data.frame\"")
    }
    if (!inherits(dose, "character") | !inherits(resp, "character") | 
        !inherits(status, "character")) {
      stop("dose, resp, and status must be of class \"character\" when supplying data")
    }
    dose.vec <- data[, dose]
    resp.vec <- data[, resp]
    status.vec <- data[, status]
  } else {
    if (inherits(dose, "character") | inherits(resp, "character") | 
        inherits(status, "character")) {
      stop("Must supply data when dose and resp are of class \"character\"")
    }
    if ((length(dose) != length(resp)) | (length(dose) != length(status))) {
      stop("dose, resp, and status must be of equal length")
    }
    dose.vec <- dose
    resp.vec <- resp
    status.vec <- status
  }

  dat <- data.frame(dose = dose.vec, resp = resp.vec, status = status.vec)
  dat$dose <- as.factor(dat$dose)
  n.doses <- length(unique(dat$dose))
  doses <- sort(unique(dose.vec))


  model <- match.arg(model)
  if (model == "coxph") {
    # Fit a Cox proportional hazards regression model By assumptions, this is a
    # placebo-adjusted estimate
    cox.mod <- survival::coxph(survival::Surv(resp, status) ~ dose, data = dat,
                               ...)
    mu.hat <- coef(cox.mod)  # Retrive the estimates coefficients
    # The cox model doesn't return an intercept term from coef()
    S.hat <- vcov(cox.mod)

    doses <- doses[-1]  # For placebo adjusted, remove the placebo dose
    # Rely on the DoseFinding package for fitting
    mod.out <- MCPMod(doses, mu.hat, models = models, S = S.hat, type = "general",
                      placAdj = TRUE, selModel = selModel, alpha = alpha, df = df, critV = critV,
                      doseType = doseType, Delta = Delta, p = p, pVal = pVal, alternative = alternative,
                      na.action = na.action, mvtcontrol = mvtcontrol, bnds = bnds, control = control,
                      ...)
  } else {
    # Otherwise fit a parametric survival model Check potential distribution
    # assumptions for survreg This could potentially include hand-written links, but
    # not currently available
    if (!(dist %in% c("weibull", "exponential", "gaussian", "logistic", "lognormal",
                      "loglogistic"))) {
      stop("dist must be one of \"weibull\", \"exponential\", \"gaussian\",
           \"logistic\", \"lognormal\", \"loglogistic\".")
    }
    # Unlike the CoxPH model, we can now have an intercept term Fit the non-placebo
    # adjusted case
    if (placAdj == FALSE) {
      surv.mod <- survival::survreg(survival::Surv(resp, status) ~ dose - 1,
                                    data = dat, dist = dist)
      mu.hat <- coef(surv.mod)  # Retrieve all coefficient estimates
      S.hat <- vcov(surv.mod)

      mod.out <- MCPMod(doses, mu.hat, models = models,
                        S = S.hat[1:n.doses, 1:n.doses], type = "general",
                        placAdj = FALSE, selModel = selModel, alpha = alpha,
                        df = df, critV = critV, doseType = doseType, Delta = Delta,
                        p = p, pVal = pVal, alternative = alternative,
                        na.action = na.action, mvtcontrol = mvtcontrol,
                        bnds = bnds, ...)

    } else {
      surv.mod <- survival::survreg(survival::Surv(resp, status) ~ dose, data = dat,
                                    dist = dist)

      mu.hat <- coef(surv.mod)  # I want to return even the intercept, so
      # still retrieve all of them
      S.hat <- vcov(surv.mod)

      # Only pass forward the non-intercept terms to MCPMod.  S is passed as
      # S.hat[2:n.doses, 2:n.doses] instead of S.hat[-1,-1], as depending on the
      # distribution, there will be additional terms in the covariance matrix
      mod.out <- MCPMod(doses[-1], mu.hat[-1], models = models,
                        S = S.hat[2:n.doses, 2:n.doses], type = "general",
                        placAdj = FALSE, selModel = selModel, alpha = alpha,
                        df = df, critV = critV, doseType = doseType,
                        Delta = Delta, p = p,
                        pVal = pVal, alternative = alternative,
                        na.action = na.action, mvtcontrol = mvtcontrol,
                        bnds = bnds, ...)
    }
  }



  if (returnS == FALSE) {
    return(mod.out)
  } else {
    data.df <- data.frame(dose = doses, resp = mu.hat)
    return.list <- list(MCPMod = mod.out, data = data.df, S = S.hat)
    return(return.list)
  }


}
