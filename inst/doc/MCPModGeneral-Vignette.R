## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include = FALSE, fig.width = 3, fig.height = 3--------------------
library(DoseFinding)
library(MCPModGeneral)

## ----powMCTGen, fig.height=3, fig.width=5-------------------------------------
dose.vec = c(0, 5, 10, 20, 30, 40)
models.full = Mods(doses = dose.vec, linear = NULL,
                   sigEmax = rbind(c(9, 4), c(20, 3)), 
                   emax = 1.25, quadratic = -0.044/2.667,
                   placEff = 0, maxEff = 2)
plot(models.full)

## ----powMCTGen2, include = TRUE-----------------------------------------------
## Look at the power for each possible DR-curve
powMCTGen(30, "negative binomial", "log", modelPar = 0.1,
          Ntype = "arm", alpha = 0.05, altModels = models.full, verbose = T)

## ----eval=FALSE---------------------------------------------------------------
#  powMCTGen(180, "negative binomial", "log", modelPar = 0.1,
#            Ntype = "total", alpha = 0.05, altModels = models.full, verbose = T)

## ----eval=FALSE---------------------------------------------------------------
#  powMCTGen(c(30,30,30,30,30,30), "negative binomial", "log", modelPar = 0.1,
#            Ntype = "actual", alpha = 0.05, altModels = models.full, verbose = T)

## ----eval=FALSE---------------------------------------------------------------
#  powMCTGen(30, "binomial", "probit",
#            Ntype = "arm", alpha = 0.05, altModels = models.full)

## -----------------------------------------------------------------------------
powMCTGen(30, "binomial", "probit",
          Ntype = "arm", alpha = 0.05, doses = c(0, 1, 2, 36, 38, 40),
          altModels = models.full, verbose = TRUE)

## -----------------------------------------------------------------------------
## Now consider power at some theoretical DR-values
powMCTGen(30, "negative binomial", "log", modelPar = 0.1,
          theoResp = c(0, 0.2, 1.8), doses = c(0, 20, 40),
          alpha = 0.05, altModels = models.full)

## -----------------------------------------------------------------------------
## Can also check type-1 error
powMCTGen(30, "negative binomial", "log", modelPar = 0.01, theoResp = rep(0, 5),
          doses = c(0, 50, 10, 20, 30),
          alpha = 0.05, altModels = models.full)



## -----------------------------------------------------------------------------
sampSizeMCTGen("binomial", "logit", upperN = 50, Ntype = "arm",
               altModels = models.full, alpha = 0.05, alRatio = c(3/2, 1/2, 1, 1, 1, 1),
               sumFct = "min", power = 0.8)

## -----------------------------------------------------------------------------
sampSizeMCTGen("negative binomial", "log", modelPar = 0.1, upperN = 50, Ntype = "arm",
               altModels = models.full, alpha = 0.05,
               sumFct = "max", power = 0.8, verbose = T)

## -----------------------------------------------------------------------------
sampSizeMCTGen("negative binomial", "log", modelPar = 0.1, upperN = 100, Ntype = "total",
               alRatio = c(3/2, 1/2, 1),
               theoResp = c(0, 0.2, 1.8), doses = c(0, 20, 40),
               altModels = models.full, alpha = 0.05)

## ----fig.height=3, fig.width=5------------------------------------------------
data(migraine)
migraine$pfrat = migraine$painfree / migraine$ntrt
migraine

models = Mods(linear = NULL, emax = 10, quadratic = c(-0.004), doses = migraine$dose)
plot(models)

## ----fig.height=3, fig.width=5------------------------------------------------
mu.S = prepareGen(family = "binomial", link = "logit", w = "ntrt", dose = "dose",
                  resp = "pfrat", data = migraine)
mcp.hand = MCPMod(dose = mu.S$data$dose, resp = mu.S$data$resp, models = models,
                  S = mu.S$S, Delta = 0.2, type = "general")
plot(mcp.hand)
mcp.hand

## -----------------------------------------------------------------------------
mcp.gen = MCPModGen("binomial", "logit", returnS = F, w = "ntrt", dose = "dose",
            resp = "pfrat", data = migraine, models = models, Delta = 0.2)
mcp.gen


## -----------------------------------------------------------------------------
## Simulate some negative binomial data according to one of the models
set.seed(188)
mean.vec = getResp(models)[,2]
dose.dat = c()
resp.dat = c()
gender.dat = c()
for(i in 1:length(migraine$dose)){
  gender.tmp = rbinom(300, 1, prob = 0.3)
  gender.dat = c(gender.dat, gender.tmp)
  dose.dat = c(dose.dat, rep(migraine$dose[i], 300))
  resp.dat = c(resp.dat, rnbinom(300, mu = exp(mean.vec[i] + 5*gender.tmp), size = 1))
}
nb.dat = data.frame(dose = dose.dat, resp = resp.dat, gender = gender.dat)
nb.dat[sample(1:nrow(nb.dat), 5),]

## ----fig.height=3, fig.width=5------------------------------------------------
mcp.nb1 = MCPModGen("negative binomial", link = "log", returnS = T,
          dose = "dose", resp = "resp", data = nb.dat, models = models, Delta = 0.6)

mcp.nb2 = MCPModGen("negative binomial", link = "log", returnS = T,
          dose = dose.dat, resp = resp.dat, models = models, Delta = 0.6)

mcp.nb3 = MCPModGen("negative binomial", link = "log", returnS = T, placAdj = T,
          dose = "dose", resp = "resp", data = nb.dat, models = models, Delta = 0.6)

mcp.nb4 = MCPModGen("negative binomial", link = "log", returnS = T, placAdj = T,
          dose = dose.dat, resp = resp.dat, models = models, Delta = 0.6)

names(mcp.nb1)
mcp.nb1$data


mcp.nb1$MCPMod$doseEst
mcp.nb4$MCPMod$doseEst

plot(mcp.nb1$MCPMod)
plot(mcp.nb4$MCPMod)

## ----fig.height=3, fig.width=5------------------------------------------------
mcp.covars = MCPModGen("negative binomial", link = "log", returnS = F, addCovars = ~ factor(gender),
                    dose = "dose", resp = "resp", data = nb.dat, models = models, Delta = 0.6)

mcp.covars
plot(mcp.covars)



## -----------------------------------------------------------------------------
TD(models, Delta = 0.6)[2]

mcp.nb1$MCPMod$doseEst[mcp.nb1$MCPMod$selMod]
mcp.covars$doseEst[mcp.covars$selMod]



## ----fig.height=3, fig.width=5------------------------------------------------
set.seed(1786)
doses = c(0, 0.1, 0.5, 0.75, 1)
n.vec = c(30, 20, 23, 19, 32)
n.doses = length(doses)
models = Mods(doses = doses, linear = NULL, emax = 0.1, exponential = 0.2,
              quadratic = -0.75, placEff = 1, maxEff = -0.3)

## Perform power-analysis
powMCTGen(n.vec, "binomial", "risk ratio", altModels = models, placEff = 0.9,
          Ntype = "actual")




## Simulate the data according to the exponential curve
means = getResp(models)[,3]*0.9
cbind(Dose = doses, Means = means)

resp.dat = c()
for(i in 1:n.doses){
  resp.dat = c(resp.dat, rbinom(1, size = n.vec[i], prob = means[i]))
}

bin.dat = data.frame(dose = doses, resp = resp.dat/n.vec, w = n.vec)




## Fit using the package

mod.pack = MCPModGen("binomial", "risk ratio", returnS = F, w = "w", dose = "dose", resp = "resp",
                     data = bin.dat, models = models, Delta = 0.1)

plot(mod.pack)
## Look at the MED estimate
mod.pack$doseEst[mod.pack$selMod]
TD(models, Delta = 0.1, direction = "decreasing")[3]



