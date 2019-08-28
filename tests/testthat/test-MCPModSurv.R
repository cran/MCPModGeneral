
test_that("Parametric Model", {

  doses = c(0, 0.05, 0.2, 0.5, 0.8, 1)
  n.doses = length(doses)
  models.shape = Mods(doses = doses, exponential = 0.2, emax = 0.05,
                      quadratic = -1.1797/1.8876, direction = "increasing")
  models = Mods(doses = doses, exponential = 0.2, emax = 0.05, quadratic = -1.1797/1.8876,
                placEff = 0.2, maxEff = -0.8)


  ## For time-to-event data, the models specify log-means

  ## Simulate from the quadratic model
  stand.resp = getResp(models)[,3]

  n = 30
  dose.dat = c()
  resp.dat = c()
  for(i in 1:n.doses){
    dose.dat = c(dose.dat, rep(doses[i], n))
    resp.dat = c(resp.dat, rexp(n, rate = 1/exp(stand.resp[i])))
  }
  dat = data.frame(dose = dose.dat, resp = resp.dat)
  dat$status = 1
  dat$status[dat$resp > 10] = 0 ## censor values above 10


  weib.mod = survival::survreg(survival::Surv(resp, status) ~ factor(dose) - 1, data = dat, dist = "weibull")
  expo.mod = survival::survreg(survival::Surv(resp, status) ~ factor(dose) - 1, data = dat, dist = "exponential")
  gaus.mod = survival::survreg(survival::Surv(resp, status) ~ factor(dose) - 1, data = dat, dist = "gaussian")
  logi.mod = survival::survreg(survival::Surv(resp, status) ~ factor(dose) - 1, data = dat, dist = "logistic")
  logn.mod = survival::survreg(survival::Surv(resp, status) ~ factor(dose) - 1, data = dat, dist = "lognormal")
  logl.mod = survival::survreg(survival::Surv(resp, status) ~ factor(dose) - 1, data = dat, dist = "loglogistic")

  weib.out = MCPMod(doses, coef(weib.mod), models = models.shape,
                   S = vcov(weib.mod)[1:n.doses,1:n.doses], type = "general",
                   Delta = 0.4, placAdj = F)
  weib.out.pack = MCPModSurv("parametric", dist = "weibull", returnS = T, dose = "dose", resp = "resp", status = "status",
                            data = dat, models = models.shape, Delta = 0.4)

  expo.out = MCPMod(doses, coef(expo.mod), models = models.shape,
                    S = vcov(expo.mod)[1:n.doses,1:n.doses], type = "general",
                    Delta = 0.4, placAdj = F)
  expo.out.pack = MCPModSurv("parametric", dist = "exponential", returnS = T, dose = "dose", resp = "resp", status = "status",
                             data = dat, models = models.shape, Delta = 0.4)

  gaus.out = MCPMod(doses, coef(gaus.mod), models = models.shape,
                    S = vcov(gaus.mod)[1:n.doses,1:n.doses], type = "general",
                    Delta = 0.4, placAdj = F)
  gaus.out.pack = MCPModSurv("parametric", dist = "gaussian", returnS = T, dose = "dose", resp = "resp", status = "status",
                             data = dat, models = models.shape, Delta = 0.4)

  logi.out = MCPMod(doses, coef(logi.mod), models = models.shape,
                    S = vcov(logi.mod)[1:n.doses,1:n.doses], type = "general",
                    Delta = 0.4, placAdj = F)
  logi.out.pack = MCPModSurv("parametric", dist = "logistic", returnS = T, dose = "dose", resp = "resp", status = "status",
                             data = dat, models = models.shape, Delta = 0.4)

  logn.out = MCPMod(doses, coef(logn.mod), models = models.shape,
                    S = vcov(logn.mod)[1:n.doses,1:n.doses], type = "general",
                    Delta = 0.4, placAdj = F)
  logn.out.pack = MCPModSurv("parametric", dist = "lognormal", returnS = T, dose = "dose", resp = "resp", status = "status",
                             data = dat, models = models.shape, Delta = 0.4)

  logl.out = MCPMod(doses, coef(logl.mod), models = models.shape,
                    S = vcov(logl.mod)[1:n.doses,1:n.doses], type = "general",
                    Delta = 0.4, placAdj = F)
  logl.out.pack = MCPModSurv("parametric", dist = "loglogistic", returnS = T, dose = "dose", resp = "resp", status = "status",
                             data = dat, models = models.shape, Delta = 0.4)


  weib.mod.plac = survival::survreg(survival::Surv(resp, status) ~ factor(dose), data = dat, dist = "weibull")
  weib.out.plac = MCPMod(doses[-1], coef(weib.mod.plac)[-1], models = models.shape,
                    S = vcov(weib.mod.plac)[2:n.doses,2:n.doses], type = "general",
                    Delta = 0.4, placAdj = F)
  weib.out.plac.pack = MCPModSurv("parametric", dist = "weibull", returnS = T, dose = "dose", resp = "resp", status = "status",
                             data = dat, models = models.shape, placAdj = T, Delta = 0.4)

  expect_equivalent(weib.out, weib.out.pack$MCPMod, tolerance = 0.001, scale = 1)
  expect_equivalent(expo.out, expo.out.pack$MCPMod, tolerance = 0.001, scale = 1)
  expect_equivalent(gaus.out, gaus.out.pack$MCPMod, tolerance = 0.001, scale = 1)
  expect_equivalent(logi.out, logi.out.pack$MCPMod, tolerance = 0.001, scale = 1)
  expect_equivalent(logn.out, logn.out.pack$MCPMod, tolerance = 0.001, scale = 1)
  expect_equivalent(logl.out, logl.out.pack$MCPMod, tolerance = 0.001, scale = 1)
  expect_equivalent(weib.out.plac, weib.out.plac.pack$MCPMod, tolerance = 0.001, scale = 1)
})




test_that("Coxph Model", {

  doses = c(0, 0.05, 0.2, 0.5, 0.8, 1)
  n.doses = length(doses)
  models.shape = Mods(doses = doses, exponential = 0.2, emax = 0.05, quadratic = -1.1797/1.8876, direction = "increasing")
  models = Mods(doses = doses, exponential = 0.2, emax = 0.05, quadratic = -1.1797/1.8876, placEff = 0.2, maxEff = -0.8)


  ## For time-to-event data, the models specify log-means

  ## Simulate from the quadratic model
  stand.resp = getResp(models)[,3]

  n = 30
  dose.dat = c()
  resp.dat = c()
  for(i in 1:n.doses){
    dose.dat = c(dose.dat, rep(doses[i], n))
    resp.dat = c(resp.dat, rexp(n, rate = 1/exp(stand.resp[i])))
  }
  dat = data.frame(dose = dose.dat, resp = resp.dat)
  dat$status = 1
  dat$status[dat$resp > 10] = 0 ## censor values above 10


  cox.mod = survival::coxph(survival::Surv(resp, status) ~ factor(dose), data = dat, model = T)

  mu.hat = coef(cox.mod)
  S.hat = vcov(cox.mod)

  mod.out = MCPMod(doses[-1], mu.hat, models = models.shape, S = S.hat, type = "general", Delta = 0.4, placAdj = T)

  mod.out.pack = MCPModSurv("coxph", returnS = T, dose = "dose", resp = "resp", status = "status",
                            data = dat, models = models.shape, Delta = 0.4)

  expect_equivalent(mod.out, mod.out.pack$MCPMod, tolerance = 0.001, scale = 1)
})
