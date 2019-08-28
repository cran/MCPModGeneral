\name{MCPModGeneral-package}
\alias{MCPModGeneral-package}
\alias{MCPModGeneral}
\docType{package}
\title{
  A Supplement to the \code{DoseFinding} Package for the General Case
}
\description{
  Analyzes non-normal data via the Multiple Comparison Procedures and Modeling approach ('MCP-Mod'). Many functions rely on the `DoseFinding` package. This package makes it so the user does not need to prespecify or calculate the \eqn{\mu} vector and \eqn{S} matrix. Instead, the user typically supplies the data in its raw form, and this package will calculate the needed objects and pass them into the `DoseFinding` functions. If the user wishes to primarily use the functions provided in the `DoseFinding` package, a singular function (`prepareGen`) will provide mu and S. The package currently handles power analysis and the 'MCP-Mod' procedure for negative binomial, Poisson, and binomial data. The 'MCP-Mod' procedure can also be applied to survival data, but power analysis is not available.

}
\details{
  \tabular{ll}{
    Package: \tab MCPModGeneral\cr
    Type: \tab Package\cr
    Version: \tab 0.1-0 \cr
    Date: \tab 2019-8-27\cr
    License: \tab GPL-3
  }
  The main functions are:\cr
  \bold{\code{\link[MCPModGeneral:prepareGen]{prepareGen}}}: Calculates the \eqn{\mu} vector and \eqn{S} matrix to be supplied to regular \code{MCPMod} functions (e.g. \code{MCPMod}, \code{MCTtest}, \code{planMod})\cr
  \bold{\code{\link[MCPModGeneral:MCPModGen]{MCPModGen}}}: Implements the full \code{MCPMod} procedure for raw negative binomial and binary data.\cr
  \bold{\code{\link[MCPModGeneral:planModPrepare]{planModPrepare}}}: Calculate the theoretical covariance matrix \eqn{S}.\cr
  \bold{\code{\link[MCPModGeneral:powMCTGen]{powMCTGen}}}: Calculates the power of the multiple contrast test.\cr
  \bold{\code{\link[MCPModGeneral:sampSizeMCTGen]{sampSizeMCTGen}}}: Calculates the sample size needed to reach the target power.\cr

  The secondary functions are:\cr
  \bold{\code{\link[MCPModGeneral:MCPModSurv]{MCPModSurv}}}: Implements the full \code{MCPMod} procedure for basic survival data. This includes a Cox-PH model and parametric survival regression models. Power analysis is not available for the survival data.\cr
  \bold{\code{\link[MCPModGeneral:simS]{simS}}}: A simulation based method for estimating the theoretical covariance matrices.

}
\author{
  Ian Laga

  Maintainer: Ian Laga <ilaga25@gmail.com>
}
\references{
  Pinheiro, J. C., Bornkamp, B., Glimm, E. and Bretz, F. (2014)
  Model-based dose finding under model uncertainty using general
  parametric models, \emph{Statistics in Medicine}, \bold{33},
  1646--1661

  Buckland, S. T., Burnham, K. P. and Augustin, N. H. (1997). Model selection an integral part
  of inference, \emph{Biometrics}, \bold{53}, 603--618
}

\keyword{ package }
\examples{
  # Analyze the binary migraine data from the DoseFinding package.
  data(migraine)
  models = Mods(linear = NULL, emax = 1, quadratic = c(-0.004),
                doses = migraine$dose)

  powMCTGen(migraine$ntrt, "binomial", "logit",
            Ntype = "actual", altModels = models)
  sampSizeMCTGen("binomial", "logit", altModels = models, power = 0.8,
                 Ntype = "arm", upperN = 30, verbose = TRUE)

  # Now analyze using binomial weights
  PFrate <- migraine$painfree/migraine$ntrt
  migraine$pfrat = migraine$painfree / migraine$ntrt
  MCPModGen("binomial","logit",returnS = TRUE, w = "ntrt", dose = "dose",
     resp = "pfrat", data = migraine, models = models, selModel = "aveAIC",
     Delta = 0.2)
}
