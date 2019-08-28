
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MCPModGeneral

<!-- badges: start -->

<!-- badges: end -->

`MCPModGeneral` is an extension of the `DoseFinding` package,
streamlining the analysis of non-normal endpoints. Almost all of the
functions in the `DoseFinding` package rely on the user supplying
\(\mu\), the estimated dose-response coefficients, and \(S\), the
variance-covariance matrix of \(\mu\). However, \(S\) is difficult to
know before-hand, and for functions like `powMCT`, differ for each
alternative model. The `MCPModGeneral` package does not require the user
to supply a matrix for \(S\) and instead calculates the theoretical
variance-covariance matrix for the negative binomial, binomial, and
Poisson distributions. Alternatively, the emperical covariance matrices
can be estimated via simulation.

Users can also use the `MCPModGeneral` package to fit the full MCPMod
procedure on negative binomial, Poisson, and binomial data, as well as
basic survival data. The relevant functions for fitting the models are
`prepareGen` and `MCPModGen`. The full capabilities of these functions
will be explored later in the vignette.

As with the `DoseFinding` package, the `MCPModGeneral` package still
requires the user to create and specify the potential dose-response
curves. The `DoseFinding` package provides the `guesst` and `Mods`
functions to create these models, but ultimately, the dose-response
curves come from prior knowledge or discussion with clinicians. For the
entirety of this vignette, it is assumed that the user is able to
construct the potential models. Recall that dose-response curves must be
constructed on the same scale as the ANOVA output, meaning if negative
binomial data is to be analyzed via a GLM with a log-link, the
dose-response curve should represent the log of the means at each dose.
The `MCPModGeneral` package allows usage with the most commonly used
links. See the `family` page from the `stats` package for a list of the
common links.

## Installation

You can install the released version of MCPModGeneral from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("MCPModGeneral")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(MCPModGeneral)
#> Loading required package: DoseFinding
#> Loading required package: lattice
#> Loading required package: mvtnorm
# Analyze the binary migraine data from the DoseFinding package.
data(migraine)
models = Mods(linear = NULL, emax = 1, quadratic = c(-0.004),
              doses = migraine$dose)

powMCTGen(migraine$ntrt, "binomial", "logit",
          Ntype = "actual", altModels = models)
#>    linear      emax quadratic 
#> 0.8637783 0.9893745 0.9148810
sampSizeMCTGen("binomial", "logit", altModels = models, power = 0.8,
               Ntype = "arm", upperN = 30, verbose = TRUE)
#> Iter: 1, N = 45, current value = 0.7357
#> Iter: 2, N = 52, current value = 0.799
#> Iter: 3, N = 56, current value = 0.829
#> Iter: 4, N = 54, current value = 0.8154
#> Iter: 5, N = 53, current value = 0.8067
#> 
#> Using "min" of power
#> Sample size calculation
#> 
#> alRatio: 1 1 1 1 1 1 1 1 
#> Total sample size: 424 
#> Sample size per arm: 53 53 53 53 53 53 53 53 
#> targFunc:

# Now analyze using binomial weights
PFrate <- migraine$painfree/migraine$ntrt
migraine$pfrat = migraine$painfree / migraine$ntrt
MCPModGen("binomial","logit",returnS = TRUE, w = "ntrt", dose = "dose",
   resp = "pfrat", data = migraine, models = models, selModel = "aveAIC",
   Delta = 0.2)
#> $MCPMod
#> MCPMod
#> 
#> Multiple Contrast Test:
#>           t-Stat   adj-p
#> linear     3.703 < 0.001
#> emax       3.636 < 0.001
#> quadratic  3.079 0.00278
#> 
#> Estimated Dose Response Models:
#> linear model
#>     e0  delta 
#> -1.710  0.006 
#> 
#> emax model
#>     e0   eMax   ed50 
#> -2.219  1.387  8.473 
#> 
#> quadratic model
#>     e0     b1     b2 
#> -1.776  0.010  0.000 
#> 
#> Model weights (AIC):
#>    linear      emax quadratic 
#>    0.3388    0.5071    0.1541 
#> 
#> Estimated TD, Delta=0.2
#>    linear      emax quadratic 
#>   33.8758    1.4274   20.9810 
#> 
#> $data
#>    dose       resp
#> 1   0.0 -2.2225424
#> 2   2.5 -1.9459101
#> 3   5.0 -2.0541237
#> 4  10.0 -1.0775589
#> 5  20.0 -1.4469190
#> 6  50.0 -1.2927683
#> 7 100.0 -1.1676052
#> 8 200.0 -0.5663955
#> 
#> $S
#>                 factor(dose)0 factor(dose)2.5 factor(dose)5 factor(dose)10
#> factor(dose)0       0.0852564       0.0000000     0.0000000      0.0000000
#> factor(dose)2.5     0.0000000       0.2857133     0.0000000      0.0000000
#> factor(dose)5       0.0000000       0.0000000     0.2256406      0.0000000
#> factor(dose)10      0.0000000       0.0000000     0.0000000      0.0837766
#> factor(dose)20      0.0000000       0.0000000     0.0000000      0.0000000
#> factor(dose)50      0.0000000       0.0000000     0.0000000      0.0000000
#> factor(dose)100     0.0000000       0.0000000     0.0000000      0.0000000
#> factor(dose)200     0.0000000       0.0000000     0.0000000      0.0000000
#>                 factor(dose)20 factor(dose)50 factor(dose)100
#> factor(dose)0        0.0000000     0.00000000      0.00000000
#> factor(dose)2.5      0.0000000     0.00000000      0.00000000
#> factor(dose)5        0.0000000     0.00000000      0.00000000
#> factor(dose)10       0.0000000     0.00000000      0.00000000
#> factor(dose)20       0.1029412     0.00000000      0.00000000
#> factor(dose)50       0.0000000     0.09103641      0.00000000
#> factor(dose)100      0.0000000     0.00000000      0.09365079
#> factor(dose)200      0.0000000     0.00000000      0.00000000
#>                 factor(dose)200
#> factor(dose)0        0.00000000
#> factor(dose)2.5      0.00000000
#> factor(dose)5        0.00000000
#> factor(dose)10       0.00000000
#> factor(dose)20       0.00000000
#> factor(dose)50       0.00000000
#> factor(dose)100      0.00000000
#> factor(dose)200      0.07464607
```
