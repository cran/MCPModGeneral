#' MCPModGen - Multiple Comparison and Modeling (General Case)
#'
#' This function allows the user to implement the MCPMod function on negative
#' binomial, Poisson, and binary data, without having to write any additional
#' code. If analyzing survival data, see the
#' \code{\link[MCPModGeneral:MCPModSurv]{MCPModSurv}} function.
#'
#' This function works by first fitting a glm with the chosen family and link.
#' The \eqn{\mu} vector and \eqn{S} matrix are extracted from the glm, and these
#' values are supplied to the MCPMod function, along with all user-defined
#' arguments. \cr Currently, the function can take the negative binomial and
#' Poisson family with a log, sqrt, identity, risk ratio, and log risk ratio
#' links, or a bernoulli family with a log, logit, probit, cauchit,
#' cloglog-link, identity, risk ratio, and log risk ratio links.
#'
#' @param family A character string containing the error distribution to be used
#'   in the model.
#' @param link A character string for the model link function.
#' @param w Either a numeric vector of the same length as dose and resp, or a
#'   character vector denoting the column name in the data.
#' @param returnS Logical determining whether muHat and SHat should be returned,
#'   in additional to the MCPMod output.
#' @param offset Either a numeric vector of the same length as dose and resp, or
#'   a character vector denoting the column name in the data.
#' @param ... Additional arguments to be passed to \code{glm} or \code{glm.nb}.
#'   This is especially useful when a fitting error is returned. In these cases,
#'   it may be useful to supply a \code{start} vector for the parameters.
#' @param dose,resp Either vectors of equal length specifying dose and response
#'   values, or character vectors specifying the names of variables in the data
#'   frame specified in \code{data}.
#' @param data Data frame with names specified in `dose`, `resp`, and optionally
#'   `w`. If data is not specified, it is assumed that `dose` and `resp` are
#'   numerical vectors
#' @param addCovars Formula specifying additive linear covariates (e.g. `~
#'   factor(gender)`).
#' @param placAdj Logical specifying whether the provided by `resp` are to be
#'   treated as placebo-adjusted estimates.
#' @inheritParams DoseFinding::MCPMod
#' @param df An optional numeric value specifying the degrees of freedom.
#'   Infinite degrees of freedom (`df=Inf`, the default), correspond to the
#'   multivariate normal distribution.
#' @return An object of class MCPMod if `returnS = FALSE`. Otherwise, a list
#'   containing an object of class MCPMod, the numeric vector \eqn{\mu}, and the
#'   numeric matrix \eqn{S}.
#' @references Buckland, S. T., Burnham, K. P. and Augustin, N. H. (1997). Model
#'   selection an integral part of inference, \emph{Biometrics}, \bold{53},
#'   603--618
#' @examples
#' # Analyze the binary migraine data from the DoseFinding package.
#' data(migraine)
#' models = Mods(linear = NULL, emax = 1, quadratic = c(-0.004), doses = migraine$dose)
#'
#' # Now analyze using binomial weights
#' PFrate <- migraine$painfree/migraine$ntrt
#' migraine$pfrat = migraine$painfree / migraine$ntrt
#' MCPModGen("binomial","logit",returnS = TRUE, w = "ntrt", dose = "dose",
#'    resp = "pfrat", data = migraine, models = models, selModel = "aveAIC",
#'    Delta = 0.2)
#' @import stats
#' @importFrom DoseFinding MCPMod Mods
#' @export
MCPModGen <- function(family = c("negative binomial", "binomial"),
                      link = c("log", "logit", "probit", "cauchit", "cloglog", "identity",
                               "log risk ratio", "risk ratio"), returnS = FALSE,
                      w = NULL, dose, resp, data = NULL, models, addCovars = ~1,
                      placAdj = FALSE, selModel = c("AIC", "maxT", "aveAIC"),
                      alpha = 0.025, df = NULL, critV = NULL,
                      doseType = c("TD", "ED"), Delta, p, pVal = TRUE,
                      alternative = c("one.sided", "two.sided"),
                      na.action = na.fail, mvtcontrol = mvtnorm.control(), bnds,
                      control = NULL, offset = NULL, ...) {

  link = match.arg(link)
  if (link == "risk ratio" | link == "log risk ratio") {
    if (placAdj != TRUE) {
      placAdj <- TRUE
      message("Forcing 'placAdj = TRUE' for risk ratio links")
    }
  }

  if (family == "negative binomial" | family == "poisson" | family == "binomial") {

    if (link == "risk ratio" | link == "log risk ratio") {
      ## Call prepareGen to retrieve the mu.hat and S.hat vector
      muS.hat <- prepareGen(family, link = link, w = w, dose = dose, resp = resp,
                            data = data, placAdj = placAdj, addCovars = addCovars,
                            offset = offset)
      doses <- muS.hat$data$dose
      mu.hat <- muS.hat$data$resp
      S.hat <- muS.hat$S

      if (link == "risk ratio") {
        ## Call the MCPMod function with the placebo-adjusted values For the risk ratio, we subtract mu.hat by 1 for the fitting
        mod.out <- MCPMod(doses[-1], mu.hat[-1] - 1, models = models,
                          S = S.hat[-1, -1], type = "general", placAdj = TRUE,
                          selModel = selModel, alpha = alpha, df = df,
                          critV = critV, doseType = doseType, Delta = Delta,
                          p = p, pVal = pVal, alternative = alternative,
                          na.action = na.action, mvtcontrol = mvtcontrol,
                          bnds = bnds, control = control, ...)
      } else {
        mod.out <- MCPMod(doses[-1], mu.hat[-1], models = models,
                          S = S.hat[-1, -1], type = "general", placAdj = TRUE,
                          selModel = selModel, alpha = alpha, df = df,
                          critV = critV, doseType = doseType, Delta = Delta,
                          p = p, pVal = pVal, alternative = alternative,
                          na.action = na.action, mvtcontrol = mvtcontrol,
                          bnds = bnds, control = control, ...)
      }

    } else {
      ## prepareGen fits the glm and returns muHat and SHat
      muS.hat <- prepareGen(family, link = link, w = w, dose = dose, resp = resp,
                            data = data, placAdj = placAdj, addCovars = addCovars,
                            offset = offset)
      doses <- muS.hat$data$dose
      mu.hat <- muS.hat$data$resp
      S.hat <- muS.hat$S

      ## Rely on the DoseFinding package to implement the procedure with the provided arguments
      mod.out <- MCPMod(doses, mu.hat, models = models, S = S.hat,
                        type = "general", placAdj = placAdj, selModel = selModel,
                        alpha = alpha, df = df, critV = critV,
                        doseType = doseType, Delta = Delta, p = p, pVal = pVal,
                        alternative = alternative, na.action = na.action,
                        mvtcontrol = mvtcontrol, bnds = bnds, control = control,
                        ...)
    }

    ## Everything is equivalent for the binomial case
  } else {
    stop("Invalid 'family' argument")
  }

  ## If user wants to returnS, return a list with the MCPMod and other vals
  if (returnS == FALSE) {
    return(mod.out)
  } else {
    data.df <- data.frame(dose = doses, resp = mu.hat)
    return.list <- list(MCPMod = mod.out, data = data.df, S = S.hat)
    return(return.list)
  }

}
