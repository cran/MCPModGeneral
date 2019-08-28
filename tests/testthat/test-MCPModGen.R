## These tests also cover prepareGen, since they rely on them
test_that("MCPModGen Binomial 1", {
  ## First do the Pinheiro example

  data(migraine)
  models <- DoseFinding::Mods(linear = NULL, emax = 1, quadratic = c(-0.004),
                doses = migraine$dose)
  PFrate <- migraine$painfree/migraine$ntrt

  doseVec <- migraine$dose
  doseVecFac <- as.factor(migraine$dose)
  ## fit logistic regression with dose as factor
  fitBin <- glm(PFrate~doseVecFac-1, family = binomial,
                weights = migraine$ntrt)
  drEst <- coef(fitBin)
  vCov <- vcov(fitBin)

  pin.mod <- DoseFinding::MCPMod(doseVec, drEst, models = models, S = vCov,
                                type = "general",
                   Delta = 0.2, selModel = "aveAIC")

  ## Now using MCPModGen with data frame
  migraine$pfrat <- migraine$painfree / migraine$ntrt
  gen.1 <- MCPModGen("binomial","logit",returnS = T, w = "ntrt", dose = "dose",
                    resp = "pfrat", data = migraine,
                    models = models, selModel = "aveAIC", Delta = 0.2)

  ## Now with raw vectors
  gen.2 <- MCPModGen("binomial", "logit", returnS = T, w = migraine$ntrt,
                    dose = migraine$dose, resp = migraine$pfrat,
                    models = models, selModel = "aveAIC", Delta = 0.2)

  expect_equal(gen.1$MCPMod$doseEst, pin.mod$doseEst, tolerance = 1e-2)
  expect_equal(gen.2$MCPMod$doseEst, pin.mod$doseEst, tolerance = 1e-2)
  expect_equivalent(gen.1$S, vCov)
  expect_equivalent(gen.2$S, vCov)
})


test_that("MCPModGen Binomial 1 (addCovars)", {
  ## First do the Pinheiro example

  set.seed(48)
  data(migraine)
  n.doses <- nrow(migraine)
  models <- DoseFinding::Mods(linear = NULL, emax = 1, quadratic = c(-0.004),
                doses = migraine$dose)

  ## Convert to long form
  dose.vec <- c()
  resp.vec <- c()
  for(i in 1:nrow(migraine)){
    dose.vec <- c(dose.vec, rep(migraine$dose[i], migraine$ntrt[i]))
    resp.vec <- c(resp.vec, rep(1, migraine$painfree[i]))
    resp.vec <- c(resp.vec, rep(0, migraine$ntrt[i] - migraine$painfree[i]))
  }
  dat <- data.frame(dose = dose.vec, resp = resp.vec)
  ## Add random gender
  dat$gender <- rbinom(nrow(dat), 1, prob = 0.3)


  ## fit logistic regression with dose as factor
  fitBin <- glm(resp~factor(dose) - 1 + gender, family = binomial, data = dat)
  drEst <- coef(fitBin)[1:n.doses]
  vCov <- vcov(fitBin)[1:n.doses, 1:n.doses]

  DF.reg <- DoseFinding::MCPMod(sort(unique(migraine$dose)), drEst,
                               models = models,
                   S = vCov, type = "general", Delta = 0.2, selModel = "aveAIC")

  fitBin <- glm(resp~factor(dose) + gender, family = binomial, data = dat)
  drEst <- coef(fitBin)[2:n.doses]
  vCov <- vcov(fitBin)[2:n.doses, 2:n.doses]

  DF.pj <- DoseFinding::MCPMod(sort(unique(migraine$dose))[-1], drEst,
                              models = models,
           placAdj = T, S = vCov, type = "general", Delta = 0.2, selModel = "aveAIC")

  ## Now using MCPModGen with data frame
  gen.reg <- MCPModGen("binomial","logit", returnS = T, dose = "dose",
                    resp = "resp", data = dat, addCovars = ~ gender,
                    models = models, selModel = "aveAIC", Delta = 0.2)

  ## Now with raw vectors
  gen.pj <- MCPModGen("binomial","logit", returnS = T, dose = "dose", placAdj = T,
                    resp = "resp", data = dat, addCovars = ~ gender,
                    models = models, selModel = "aveAIC", Delta = 0.2)

  expect_equal(gen.reg$MCPMod$MCTtest, DF.reg$MCTtest, tolerance = 1e-2)
  expect_equal(gen.pj$MCPMod$MCTtest, DF.pj$MCTtest, tolerance = 1e-2)
})




test_that("MCPModGen NB 1", {
  set.seed(31)
  ## Simulate data and use DoseFinding
  dose.vec <- c(0, 5, 10, 20, 30, 40)
  models <- DoseFinding::Mods(doses = dose.vec, linear = NULL,
                              sigEmax = rbind(c(9,1), c(15,3)),
                              emax = 0.7, quadratic = -0.05,
                              placEff = 0, maxEff = 2)
  lin.means <- DoseFinding::getResp(models)[,1]

  n <- 30
  dat.0 <- rnbinom(n*3/2, mu = exp(lin.means[1]), size = .1)
  dat.5 <- rnbinom(n/2, mu = exp(lin.means[2]), size = .1)
  dat.10 <- rnbinom(n, mu = exp(lin.means[3]), size = .1)
  dat.20 <- rnbinom(n, mu = exp(lin.means[4]), size = .1)
  dat.30 <- rnbinom(n, mu = exp(lin.means[5]), size = .1)
  dat.40 <- rnbinom(n, mu = exp(lin.means[6]), size = .1)

  dat <- data.frame(dose = c(rep(0, n*3/2), rep(5, n/2), rep(c(10,20,30,40), each = n)),
                   resp = c(dat.0, dat.5, dat.10, dat.20, dat.30, dat.40))

  ## Now get the initial glm estimates and covariance matrix
  glm.fit <- MASS::glm.nb(resp ~ factor(dose) - 1, data = dat)
  mu.hat <- coef(glm.fit)
  S.hat <- vcov(glm.fit)
  mcp.dat <- data.frame(dose = dose.vec, resp = mu.hat)

  ## Now that we have models, perform the test and find the optimal dose
  set.seed(51)
  pin.mod <- DoseFinding::MCPMod(dose, resp, mcp.dat, models = models,
                                S = S.hat, type = "general",
                   selModel = "aveAIC", alpha = 0.06, Delta = 0.2)

  ## Now run my models

  set.seed(51)
  gen.1 <- MCPModGen(family = "negative binomial", link = "log",
              returnS = T, dose = "dose", resp = "resp", data = dat,
              models = models, alpha = 0.06, Delta = 0.2, selModel = 'aveAIC')

  set.seed(51)
  gen.2 <- MCPModGen(family = "negative binomial", link = "log",
                         returnS = T, dose = dat$dose, resp = dat$resp,
            models = models, alpha = 0.06, Delta = 0.2, selModel = 'aveAIC')

  expect_equal(gen.1$MCPMod, pin.mod, tolerance = 1e-1)
  expect_equal(gen.2$MCPMod, pin.mod, tolerance = 1e-1)
  expect_equal(gen.1$S, S.hat)
  expect_equal(gen.2$S, S.hat)
})


test_that("MCPModGen NB 1 (addCovars)", {
  ## Simulate data and use DoseFinding
  set.seed(57)
  dose.vec <- c(0, 5, 10, 20, 30, 40)
  n.doses <- length(dose.vec)
  models <- DoseFinding::Mods(doses = dose.vec, linear = NULL,
                sigEmax = rbind(c(9,1), c(15,3)),
                emax = 0.7, quadratic = -0.05,
                placEff = 0, maxEff = 2)
  lin.means <- DoseFinding::getResp(models)[,1]

  n <- 30
  dat.0 <- rnbinom(n*3/2, mu = exp(lin.means[1]), size = .1)
  dat.5 <- rnbinom(n/2, mu = exp(lin.means[2]), size = .1)
  dat.10 <- rnbinom(n, mu = exp(lin.means[3]), size = .1)
  dat.20 <- rnbinom(n, mu = exp(lin.means[4]), size = .1)
  dat.30 <- rnbinom(n, mu = exp(lin.means[5]), size = .1)
  dat.40 <- rnbinom(n, mu = exp(lin.means[6]), size = .1)

  dat <- data.frame(dose = c(rep(0, n*3/2), rep(5, n/2), rep(c(10,20,30,40), each = n)),
                   resp = c(dat.0, dat.5, dat.10, dat.20, dat.30, dat.40))
  ## Randomly add a covariate
  dat$gender <- rbinom(nrow(dat), 1, prob = 0.3)

  ## Now get the initial glm estimates and covariance matrix
  glm.fit <- MASS::glm.nb(resp ~ factor(dose) - 1 + gender, data = dat)
  mu.hat <- coef(glm.fit)[1:n.doses]
  S.hat <- vcov(glm.fit)[1:n.doses, 1:n.doses]
  mcp.dat <- data.frame(dose = dose.vec, resp = mu.hat)

  ## Now that we have models, perform the test and find the optimal dose
  set.seed(51)
  DF.reg <- DoseFinding::MCPMod(dose, resp, mcp.dat, models = models, S = S.hat,
                   type = "general",
                   selModel = "aveAIC", alpha = 0.06, Delta = 0.2)

  glm.fit <- MASS::glm.nb(resp ~ factor(dose) + gender, data = dat)
  mu.hat <- coef(glm.fit)[2:n.doses]
  S.hat <- vcov(glm.fit)[2:n.doses, 2:n.doses]
  mcp.dat <- data.frame(dose = dose.vec[-1], resp = mu.hat)
  DF.pj <- DoseFinding::MCPMod(dose, resp, mcp.dat, models = models, S = S.hat,
                  type = "general", placAdj = T,
                  selModel = "aveAIC", alpha = 0.06, Delta = 0.2)

  ## Now run my models

  set.seed(51)
  gen.reg <- MCPModGen(family = "negative binomial", link = "log", returnS = T,
              dose = "dose", resp = "resp", data = dat, models = models,
              alpha = 0.06, Delta = 0.2, selModel = 'aveAIC', addCovars = ~gender)

  set.seed(51)
  gen.pj <- MCPModGen(family = "negative binomial", link = "log", returnS = T,
                    dose = "dose", resp = "resp", data = dat,
                    models = models, alpha = 0.06, Delta = 0.2, selModel = 'aveAIC',
                    addCovars = ~gender, placAdj = T)

  expect_equal(gen.reg$MCPMod$MCTtest, DF.reg$MCTtest, tolerance = 1e-1)
  expect_equal(gen.pj$MCPMod$MCTtest, DF.pj$MCTtest, tolerance = 1e-1)
})





test_that("MCPModGen NB Log Risk Ratio", {
  ## First generate a negative binomial data set
  doses <- c(0, 1, 2, 3, 4)
  n.doses <- length(doses)
  n.vec <- c(60, 40, 30, 40, 53)
  models <- DoseFinding::Mods(doses = doses, linear = NULL, emax = 1.5,
                             sigEmax = c(2, 8), exponential = 0.7,
                placEff = 0, maxEff = -0.7)
  means <- exp(getResp(models)[,2]) * 8

  dose.vec <- c()
  resp.vec <- c()
  for(i in 1:n.doses){
    dose.vec <- c(dose.vec, rep(doses[i], n.vec[i]))
    resp.vec <- c(resp.vec, rnbinom(n.vec[i], size = 0.7, mu = means[i]))
  }
  dat <- data.frame(dose = dose.vec, resp = resp.vec)

  ## Now fit the log risk ratio by hand

  glm.fit <- MASS::glm.nb(resp ~ factor(dose), data = dat)

  ## Do log-risk ratio first
  mu.hat <- coef(glm.fit)[-1]
  S.hat <- vcov(glm.fit)[-1,-1]
  dose.in <- doses[-1]

  set.seed(51)
  mod.DF.log <- DoseFinding::MCPMod(dose.in, mu.hat, models = models,
                  S = S.hat, type = "general", placAdj = T, Delta = 0.3)

  ## And now using the package
  set.seed(51)
  mod.pack.log <- MCPModGen(family = "negative binomial", link = "log risk ratio",
          dose = "dose", resp = "resp",
          data = dat, models = models, Delta = 0.3, returnS = T)


  expect_equal(mod.DF.log$MCTtest, mod.pack.log$MCPMod$MCTtest)
})




test_that("MCPModGen NB Risk Ratio", {
  doses <- c(0, 1, 2, 3, 4)
  n.doses <- length(doses)
  n.vec <- c(60, 40, 30, 40, 53)
  models <- DoseFinding::Mods(doses = doses, linear = NULL, emax = 1.5,
                             sigEmax = c(2, 8), exponential = 0.7,
                placEff = 1, maxEff = -0.7)
  means <- DoseFinding::getResp(models)[,2] * 8

  dose.vec <- c()
  resp.vec <- c()
  for(i in 1:n.doses){
    dose.vec <- c(dose.vec, rep(doses[i], n.vec[i]))
    resp.vec <- c(resp.vec, rnbinom(n.vec[i], size = 0.7, mu = means[i]))
  }
  dat <- data.frame(dose = dose.vec, resp = resp.vec)

  ## Now fit the risk ratio by hand

  glm.fit <- MASS::glm.nb(resp ~ factor(dose), data = dat)

  mu.hat <- coef(glm.fit)
  S.hat <- vcov(glm.fit)
  S.hat <- diag(exp(mu.hat))%*%S.hat%*%diag(exp(mu.hat))
  mu.hat <- exp(mu.hat)

  S.hat <- S.hat[-1,-1]
  mu.hat <- mu.hat[-1] - 1
  dose.in <- doses[-1]

  mod.DF <- DoseFinding::MCPMod(dose.in, mu.hat, models = models,
                               S = S.hat, type = "general",
                  placAdj = T, Delta = 0.3)

  mod.pack <- MCPModGen(family = "negative binomial", link = "risk ratio",
                 dose = "dose", resp = "resp",
                 data = dat, models = models, Delta = 0.3, returnS = T)


  expect_equal(mod.DF$MCTtest, mod.pack$MCPMod$MCTtest, tolerance = 1e-2, scale = 1)

})





test_that("MCPModGen Binomial Log Risk Ratio", {
  ## First generate a negative binomial data set
  doses <- c(0, 1, 2, 3, 4)
  n.doses <- length(doses)
  n.vec <- c(60, 40, 30, 40, 53)
  models <- DoseFinding::Mods(doses = doses, linear = NULL, emax = 1.5,
                             sigEmax = c(2, 8), exponential = 0.7,
                placEff = 0, maxEff = -0.7)
  means <- exp(DoseFinding::getResp(models)[,2]) * 0.8

  dose.vec <- doses
  resp.vec <- c()
  for(i in 1:n.doses){
    resp.vec <- c(resp.vec, rbinom(1, size = n.vec[i], prob = means[i]))
  }
  dat <- data.frame(dose = dose.vec, resp = resp.vec/n.vec, w = n.vec)

  ## Now fit the log risk ratio by hand

  risks.vec <- rep(NA, n.doses)
  k <- 1
  for(i in doses){
    risks.vec[k] <- mean(dat$resp[dat$dose == i])
    k <- k + 1
  }
  emp.rr <- risks.vec/risks.vec[1]

  glm.fit <- glm(resp ~ factor(dose), data = dat, weights = w,
                family = binomial(link = "log"),
                start = log(c(mean(dat$resp[dat$dose == 0]), emp.rr[-1])))

  ## Do log-risk ratio first
  mu.hat <- coef(glm.fit)[-1]
  S.hat <- vcov(glm.fit)[-1,-1]
  dose.in <- doses[-1]

  set.seed(51)
  mod.DF.log <- DoseFinding::MCPMod(dose.in, mu.hat, models = models,
                                   S = S.hat, type = "general",
                      placAdj = T, Delta = 0.3)

  ## And now using the package
  set.seed(51)
  mod.pack.log <- MCPModGen(family = "binomial", link = "log risk ratio",
                           dose = "dose", resp = "resp", w = "w",
                           data = dat, models = models, Delta = 0.3, returnS = T)


  expect_equal(mod.DF.log$MCTtest, mod.pack.log$MCPMod$MCTtest)
})

test_that("MCPModGen Binomial Log Risk Ratio - long form", {
  ## First generate a negative binomial data set
  doses <- c(0, 1, 2, 3, 4)
  n.doses <- length(doses)
  n.vec <- c(60, 40, 30, 40, 53)
  models <- DoseFinding::Mods(doses = doses, linear = NULL, emax = 1.5,
                             sigEmax = c(2, 8), exponential = 0.7,
                placEff = 0, maxEff = -0.7)
  means <- exp(DoseFinding::getResp(models)[,2]) * 0.8

  dose.vec <- c()
  resp.vec <- c()
  for(i in 1:n.doses){
    dose.vec <- c(dose.vec, rep(doses[i], n.vec[i]))
    resp.vec <- c(resp.vec, rbinom(n.vec[i], size = 1, prob = means[i]))
  }
  dat <- data.frame(dose = dose.vec, resp = resp.vec)

  ## Now fit the log risk ratio by hand

  risks.vec <- rep(NA, n.doses)
  k <- 1
  for(i in doses){
    risks.vec[k] <- mean(dat$resp[dat$dose == i])
    k <- k + 1
  }
  emp.rr <- risks.vec/risks.vec[1]

  glm.fit <- glm(resp ~ factor(dose), data = dat,
                family = binomial(link = "log"),
                start = log(c(mean(dat$resp[dat$dose == 0]), emp.rr[-1])))

  ## Do log-risk ratio first
  mu.hat <- coef(glm.fit)[-1]
  S.hat <- vcov(glm.fit)[-1,-1]
  dose.in <- doses[-1]

  set.seed(51)
  mod.DF.log <- DoseFinding::MCPMod(dose.in, mu.hat, models = models,
                                   S = S.hat, type = "general",
                      placAdj = T, Delta = 0.3)

  ## And now using the package
  set.seed(51)
  mod.pack.log <- MCPModGen(family = "binomial", link = "log risk ratio",
                           dose = "dose", resp = "resp",
                           data = dat, models = models, Delta = 0.3, returnS = T)


  expect_equal(mod.DF.log$MCTtest, mod.pack.log$MCPMod$MCTtest)
})




test_that("MCPModGen NB Risk Ratio", {
  doses <- c(0, 1, 2, 3, 4)
  n.doses <- length(doses)
  n.vec <- c(60, 40, 30, 40, 53)
  models <- DoseFinding::Mods(doses = doses, linear = NULL, emax = 1.5, sigEmax = c(2, 8), exponential = 0.7,
                placEff = 1, maxEff = -0.7)
  means <- DoseFinding::getResp(models)[,2] * 8

  dose.vec <- c()
  resp.vec <- c()
  for(i in 1:n.doses){
    dose.vec <- c(dose.vec, rep(doses[i], n.vec[i]))
    resp.vec <- c(resp.vec, rnbinom(n.vec[i], size = 0.7, mu = means[i]))
  }
  dat <- data.frame(dose = dose.vec, resp = resp.vec)

  ## Now fit the risk ratio by hand

  glm.fit <- MASS::glm.nb(resp ~ factor(dose), data = dat)

  mu.hat <- coef(glm.fit)
  S.hat <- vcov(glm.fit)
  S.hat <- diag(exp(mu.hat))%*%S.hat%*%diag(exp(mu.hat))
  mu.hat <- exp(mu.hat)

  S.hat <- S.hat[-1,-1]
  mu.hat <- mu.hat[-1] - 1
  dose.in <- doses[-1]

  set.seed(51)
  mod.DF <- DoseFinding::MCPMod(dose.in, mu.hat, models = models,
                               S = S.hat, type = "general",
                  placAdj = T, Delta = 0.3)

  set.seed(51)
  mod.pack <- MCPModGen(family = "negative binomial", link = "risk ratio",
                       dose = "dose", resp = "resp",
                       data = dat, models = models, Delta = 0.3, returnS = T)


  expect_equivalent(mod.DF$MCTtest, mod.pack$MCPMod$MCTtest)

})








test_that("MCPModGen Binomial Risk Ratio", {
  ## First generate a negative binomial data set
  doses <- c(0, 1, 2, 3, 4)
  n.doses <- length(doses)
  n.vec <- c(60, 40, 30, 40, 53)
  models <- DoseFinding::Mods(doses = doses, linear = NULL, emax = 1.5,
                             sigEmax = c(2, 8), exponential = 0.7,
                placEff = 0, maxEff = -0.7)
  means <- exp(DoseFinding::getResp(models)[,2]) * 0.8

  dose.vec <- doses
  resp.vec <- c()
  for(i in 1:n.doses){
    resp.vec <- c(resp.vec, rbinom(1, size = n.vec[i], prob = means[i]))
  }
  dat <- data.frame(dose = dose.vec, resp = resp.vec/n.vec, w = n.vec)

  ## Now fit the log risk ratio by hand

  risks.vec <- rep(NA, n.doses)
  k <- 1
  for(i in doses){
    risks.vec[k] <- mean(dat$resp[dat$dose == i])
    k <- k + 1
  }
  emp.rr <- risks.vec/risks.vec[1]

  glm.fit <- glm(resp ~ factor(dose), data = dat, weights = w,
                family = binomial(link = "log"),
                start = log(c(mean(dat$resp[dat$dose == 0]), emp.rr[-1])))

  mu.hat <- coef(glm.fit)
  S.hat <- vcov(glm.fit)
  S.hat <- diag(exp(mu.hat))%*%S.hat%*%diag(exp(mu.hat))
  mu.hat <- exp(mu.hat)

  S.hat <- S.hat[-1,-1]
  mu.hat <- mu.hat[-1] - 1
  dose.in <- doses[-1]

  set.seed(51)
  mod.DF.log <- DoseFinding::MCPMod(dose.in, mu.hat, models = models,
                                   S = S.hat, type = "general",
                      placAdj = T, Delta = 0.3)

  ## And now using the package
  set.seed(51)
  mod.pack.log <- MCPModGen(family = "binomial", link = "risk ratio",
                           dose = "dose", resp = "resp", w = "w",
                           data = dat, models = models, Delta = 0.3, returnS = T)


  expect_equal(mod.DF.log$MCTtest, mod.pack.log$MCPMod$MCTtest)
})

test_that("MCPModGen Binomial Risk Ratio - long form", {
  ## First generate a negative binomial data set
  doses <- c(0, 1, 2, 3, 4)
  n.doses <- length(doses)
  n.vec <- c(60, 40, 30, 40, 53)
  models <- DoseFinding::Mods(doses = doses, linear = NULL,
                             emax = 1.5, sigEmax = c(2, 8), exponential = 0.7,
                placEff = 0, maxEff = -0.7)
  means <- exp(DoseFinding::getResp(models)[,2]) * 0.8

  dose.vec <- c()
  resp.vec <- c()
  for(i in 1:n.doses){
    dose.vec <- c(dose.vec, rep(doses[i], n.vec[i]))
    resp.vec <- c(resp.vec, rbinom(n.vec[i], size = 1, prob = means[i]))
  }
  dat <- data.frame(dose = dose.vec, resp = resp.vec)

  ## Now fit the log risk ratio by hand

  risks.vec <- rep(NA, n.doses)
  k <- 1
  for(i in doses){
    risks.vec[k] <- mean(dat$resp[dat$dose == i])
    k <- k + 1
  }
  emp.rr <- risks.vec/risks.vec[1]

  glm.fit <- glm(resp ~ factor(dose), data = dat,
                family = binomial(link = "log"),
                start = log(c(mean(dat$resp[dat$dose == 0]), emp.rr[-1])))

  mu.hat <- coef(glm.fit)
  S.hat <- vcov(glm.fit)
  S.hat <- diag(exp(mu.hat))%*%S.hat%*%diag(exp(mu.hat))
  mu.hat <- exp(mu.hat)

  S.hat <- S.hat[-1,-1]
  mu.hat <- mu.hat[-1] - 1
  dose.in <- doses[-1]

  set.seed(51)
  mod.DF.log <- DoseFinding::MCPMod(dose.in, mu.hat, models = models,
                                   S = S.hat, type = "general",
                      placAdj = T, Delta = 0.3)

  ## And now using the package
  set.seed(51)
  mod.pack.log <- MCPModGen(family = "binomial", link = "risk ratio",
                           dose = "dose", resp = "resp",
                           data = dat, models = models, Delta = 0.3, returnS = T)


  expect_equal(mod.DF.log$MCTtest, mod.pack.log$MCPMod$MCTtest)
})










test_that("MCPModGen Poisson Log Risk Ratio", {
  ## First generate a negative binomial data set
  doses <- c(0, 1, 2, 3, 4)
  n.doses <- length(doses)
  n.vec <- c(60, 40, 30, 40, 53)
  models <- DoseFinding::Mods(doses = doses, linear = NULL, emax = 1.5,
                             sigEmax = c(2, 8), exponential = 0.7,
                placEff = 0, maxEff = -0.7)
  means <- exp(DoseFinding::getResp(models)[,2]) * 5

  dose.vec <- c()
  resp.vec <- c()
  for(i in 1:n.doses){
    dose.vec <- c(dose.vec, rep(doses[i], n.vec[i]))
    resp.vec <- c(resp.vec, rpois(n.vec[i], lambda = means[i]))
  }
  dat <- data.frame(dose = dose.vec, resp = resp.vec)



  glm.fit <- glm(resp ~ factor(dose), data = dat,
                family = poisson(link = "log"))

  ## Do log-risk ratio first
  mu.hat <- coef(glm.fit)[-1]
  S.hat <- vcov(glm.fit)[-1,-1]
  dose.in <- doses[-1]

  set.seed(51)
  mod.DF.log <- DoseFinding::MCPMod(dose.in, mu.hat, models = models,
                                   S = S.hat, type = "general",
                      placAdj = T, Delta = 0.3)

  ## And now using the package
  set.seed(51)
  mod.pack.log <- MCPModGen(family = "poisson", link = "log risk ratio",
                           dose = "dose", resp = "resp",
                           data = dat, models = models, Delta = 0.3, returnS = T)


  expect_equal(mod.DF.log$MCTtest, mod.pack.log$MCPMod$MCTtest)
})




test_that("MCPModGen Poisson Risk Ratio", {
  doses <- c(0, 1, 2, 3, 4)
  n.doses <- length(doses)
  n.vec <- c(60, 40, 30, 40, 53)
  models <- DoseFinding::Mods(doses = doses, linear = NULL, emax = 1.5,
                             sigEmax = c(2, 8), exponential = 0.7,
                placEff = 1, maxEff = -0.7)
  means <- exp(DoseFinding::getResp(models)[,2]) * 5

  dose.vec <- c()
  resp.vec <- c()
  for(i in 1:n.doses){
    dose.vec <- c(dose.vec, rep(doses[i], n.vec[i]))
    resp.vec <- c(resp.vec, rpois(n.vec[i], lambda = means[i]))
  }
  dat <- data.frame(dose = dose.vec, resp = resp.vec)

  ## Now fit the risk ratio by hand

  glm.fit <- glm(resp ~ factor(dose), data = dat,
                family = poisson(link = "log"))

  mu.hat <- coef(glm.fit)
  S.hat <- vcov(glm.fit)
  S.hat <- diag(exp(mu.hat))%*%S.hat%*%diag(exp(mu.hat))
  mu.hat <- exp(mu.hat)

  S.hat <- S.hat[-1,-1]
  mu.hat <- mu.hat[-1] - 1
  dose.in <- doses[-1]

  set.seed(51)
  mod.DF <- DoseFinding::MCPMod(dose.in, mu.hat, models = models,
                               S = S.hat, type = "general",
                  placAdj = T, Delta = 0.3)

  set.seed(51)
  mod.pack <- MCPModGen(family = "poisson", link = "risk ratio",
                       dose = "dose", resp = "resp",
                       data = dat, models = models, Delta = 0.3, returnS = T)


  expect_equal(mod.DF$MCTtest, mod.pack$MCPMod$MCTtest)

})







test_that("MCPMod Binomial addCovars", {

  dose.vec <- c(0, 2.5, 5, 10, 20, 50)
  dat <- data.frame(dose = c(rep(0, 133), rep(2.5, 32), rep(5, 44),
                            rep(10, 63), rep(20, 63), rep(50, 65)))
  dat$gender <- rbinom(nrow(dat), size = 1, prob = 1/3)
  dat$painfree <- rbinom(nrow(dat), size = 1,
              prob = plogis(-2 + 0.2*dat$dose - 0.5*dat$gender - 6*(dat$dose == 50)))

  ## Now define the model
  models <- DoseFinding::Mods(doses = dose.vec, linear = NULL,
                             sigEmax = rbind(c(30,4), c(100,3)),
                emax = 3, quadratic = -0.009/2.667)

  ## Covariate and placAdj
  fitBin.cp <- glm(painfree ~ factor(dose) + factor(gender), data = dat, family = "binomial")
  mu.hat.hand.cp <- coef(fitBin.cp)[2:6]
  S.hat.hand.cp <- vcov(fitBin.cp)[2:6,2:6]
  pin.mod.cp <- DoseFinding::MCPMod(dose.vec[-1], mu.hat.hand.cp, models = models,
                      S = S.hat.hand.cp, type = "general", placAdj = T,
                      selModel = "aveAIC", alpha = 0.074, Delta = 0.2)

  ## No covariate and placAdj
  fitBin.ncp <- glm(painfree ~ factor(dose), data = dat, family = "binomial")
  mu.hat.hand.ncp <- coef(fitBin.ncp)[2:6]
  S.hat.hand.ncp <- vcov(fitBin.ncp)[2:6,2:6]
  pin.mod.ncp <- DoseFinding::MCPMod(dose.vec[-1], mu.hat.hand.ncp,
            models = models, S = S.hat.hand.ncp, type = "general", placAdj = T,
            selModel = "aveAIC", alpha = 0.074, Delta = 0.2)

  ## covariate and no placAdj
  fitBin.cnp <- glm(painfree ~ factor(dose) + factor(gender) - 1, data = dat, family = "binomial")
  mu.hat.hand.cnp <- coef(fitBin.cnp)[1:6]
  S.hat.hand.cnp <- vcov(fitBin.cnp)[1:6,1:6]
  pin.mod.cnp <- DoseFinding::MCPMod(dose.vec, mu.hat.hand.cnp, models = models,
                        S = S.hat.hand.cnp, type = "general", placAdj = F,
                       selModel = "aveAIC", alpha = 0.074, Delta = 0.2)

  ## no covariate and no placAdj
  fitBin.ncnp <- glm(painfree ~ factor(dose) - 1, data = dat, family = "binomial")
  mu.hat.hand.ncnp <- coef(fitBin.ncnp)[1:6]
  S.hat.hand.ncnp <- vcov(fitBin.ncnp)[1:6,1:6]
  pin.mod.ncnp <- DoseFinding::MCPMod(dose.vec, mu.hat.hand.ncnp, models = models,
                        S = S.hat.hand.ncnp, type = "general", placAdj = F,
                        selModel = "aveAIC", alpha = 0.074, Delta = 0.2)


  gen.cp <- MCPModGen("binomial", "logit", returnS = T, dose = "dose", resp = "painfree", data = dat, models = models,
                     addCovars = ~ factor(gender), placAdj = T, alpha = 0.074, selModel = "aveAIC", Delta = 0.2)

  gen.ncp <- MCPModGen("binomial", "logit", returnS = T, dose = "dose", resp = "painfree", data = dat, models = models,
                      placAdj = T, alpha = 0.074, selModel = "aveAIC", Delta = 0.2)

  gen.cnp <- MCPModGen("binomial", "logit", returnS = T, dose = "dose", resp = "painfree", data = dat, models = models,
                      addCovars = ~ factor(gender), alpha = 0.074, selModel = "aveAIC", Delta = 0.2)

  gen.ncnp <- MCPModGen("binomial", "logit", returnS = T, dose = "dose", resp = "painfree", data = dat, models = models,
                       alpha = 0.074, selModel = "aveAIC", Delta = 0.2)

  ## Also make sure it's okay if the data.frame is formatted correctly anyways.
  dat.factor <- dat
  dat.factor$dose <- factor(dat.factor$dose)
  dat.factor$gender <- factor(dat.factor$gender)
  gen.cp2 <- MCPModGen("binomial", "logit", returnS = T, dose = "dose", resp = "painfree", data = dat, models = models,
                      addCovars = ~ gender, placAdj = T, alpha = 0.074, selModel = "aveAIC", Delta = 0.2)

  ## This one should fail
  expect_error(gen.wrong.cov = MCPModGen("binomial", "logit", returnS = T, dose = "dose", resp = "painfree", data = dat, models = models,
                                         addCovars = ~ hello, alpha = 0.074, selModel = "aveAIC", Delta = 0.2))



  expect_equivalent(gen.cp$S, S.hat.hand.cp)
  expect_equivalent(gen.ncp$S, S.hat.hand.ncp)
  expect_equivalent(gen.cnp$S, S.hat.hand.cnp)
  expect_equivalent(gen.ncnp$S, S.hat.hand.ncnp)
  expect_equivalent(gen.cp2$S, S.hat.hand.cp)

  expect_equivalent(gen.cp$MCPMod, pin.mod.cp)
  expect_equivalent(gen.ncp$MCPMod, pin.mod.ncp)
  expect_equivalent(gen.cnp$MCPMod, pin.mod.cnp)
  expect_equivalent(gen.ncnp$MCPMod, pin.mod.ncnp)
  expect_equivalent(gen.cp2$MCPMod, pin.mod.cp)


})



# Do Three tests for Risk Ratio link with covariates ----------------------



test_that("MCPModGen NB Risk Ratio (addCovars)", {
  doses <- c(0, 1, 2, 3, 4)
  n.doses <- length(doses)
  n.vec <- c(60, 40, 30, 40, 53)
  models <- DoseFinding::Mods(doses = doses, linear = NULL, emax = 1.5,
                             sigEmax = c(2, 8), exponential = 0.7,
                placEff = 1, maxEff = -0.7)
  means <- DoseFinding::getResp(models)[,2] * 8

  dose.vec <- c()
  resp.vec <- c()
  for(i in 1:n.doses){
    dose.vec <- c(dose.vec, rep(doses[i], n.vec[i]))
    resp.vec <- c(resp.vec, rnbinom(n.vec[i], size = 0.7, mu = means[i]))
  }
  dat <- data.frame(dose = dose.vec, resp = resp.vec)
  dat$gender <- rbinom(nrow(dat), 1, prob = 0.3)

  ## Now fit the risk ratio by hand

  glm.fit <- MASS::glm.nb(resp ~ factor(dose) + gender, data = dat)

  mu.hat <- coef(glm.fit)
  S.hat <- vcov(glm.fit)
  S.hat <- diag(exp(mu.hat))%*%S.hat%*%diag(exp(mu.hat))
  mu.hat <- exp(mu.hat)

  S.hat <- S.hat[2:n.doses, 2:n.doses]
  mu.hat <- mu.hat[2:n.doses] - 1
  dose.in <- doses[-1]

  mod.DF <- DoseFinding::MCPMod(dose.in, mu.hat, models = models, S = S.hat,
                  type = "general",
                  placAdj = T, Delta = 0.3)

  mod.pack <- MCPModGen(family = "negative binomial", link = "risk ratio",
                       dose = "dose", resp = "resp", addCovars = ~gender,
                       data = dat, models = models, Delta = 0.3, returnS = T)


  expect_equal(mod.DF$MCTtest, mod.pack$MCPMod$MCTtest, tolerance = 1e-2, scale = 1)

})




## This is the one that doesn't work because it requires a start vector
test_that("MCPModGen Binomial Risk Ratio (addCovars)", {
  ## First generate a negative binomial data set
  doses <- c(0, 1, 2, 3, 4)
  n.doses <- length(doses)
  n.vec <- c(60, 40, 30, 40, 53)
  models <- DoseFinding::Mods(doses = doses, linear = NULL, emax = 1.5,
                             sigEmax = c(2, 8), exponential = 0.7,
                placEff = 0, maxEff = -0.7)
  means <- exp(DoseFinding::getResp(models)[,2]) * 0.8

  dose.vec <- c()
  resp.vec <- c()
  for(i in 1:n.doses){
    dose.vec <- c(dose.vec, rep(doses[i], n.vec[i]))
    resp.vec <- c(resp.vec, rbinom(n.vec[i], size = 1, prob = means[i]))
  }
  dat <- data.frame(dose = dose.vec, resp = resp.vec)
  dat$gender <- rbinom(nrow(dat), 1, prob = 0.3)

  ## Now fit the log risk ratio by hand

  risks.vec <- rep(NA, n.doses)
  k <- 1
  for(i in doses){
    risks.vec[k] <- mean(dat$resp[dat$dose == i])
    k <- k + 1
  }
  emp.rr <- risks.vec/risks.vec[1]

  glm.fit <- glm(resp ~ factor(dose) + gender, data = dat,
                family = binomial(link = "log"),
                start = log(c(mean(dat$resp[dat$dose == 0]), emp.rr[-1], 1)))

  mu.hat <- coef(glm.fit)
  S.hat <- vcov(glm.fit)
  S.hat <- diag(exp(mu.hat))%*%S.hat%*%diag(exp(mu.hat))
  mu.hat <- exp(mu.hat)

  S.hat <- S.hat[2:n.doses,2:n.doses]
  mu.hat <- mu.hat[2:n.doses] - 1
  dose.in <- doses[-1]

  set.seed(51)
  mod.DF.log <- DoseFinding::MCPMod(dose.in, mu.hat, models = models, S = S.hat,
                                   type = "general",
                      placAdj = T, Delta = 0.3)

  ## And now using the package
  set.seed(51)
  mod.pack.log <- MCPModGen(family = "binomial", link = "risk ratio",
                           dose = "dose", resp = "resp", addCovars = ~gender,
                           data = dat, models = models, Delta = 0.3, returnS = T)


  expect_equal(mod.DF.log$MCTtest, mod.pack.log$MCPMod$MCTtest)
})



test_that("MCPModGen Poisson Risk Ratio (addCovars)", {
  doses <- c(0, 1, 2, 3, 4)
  n.doses <- length(doses)
  n.vec <- c(60, 40, 30, 40, 53)
  models <- DoseFinding::Mods(doses = doses, linear = NULL, emax = 1.5,
                             sigEmax = c(2, 8), exponential = 0.7,
                placEff = 1, maxEff = -0.7)
  means <- exp(DoseFinding::getResp(models)[,2]) * 5

  dose.vec <- c()
  resp.vec <- c()
  for(i in 1:n.doses){
    dose.vec <- c(dose.vec, rep(doses[i], n.vec[i]))
    resp.vec <- c(resp.vec, rpois(n.vec[i], lambda = means[i]))
  }
  dat <- data.frame(dose = dose.vec, resp = resp.vec)
  dat$gender <- rbinom(nrow(dat), 1, prob = 0.3)

  ## Now fit the risk ratio by hand

  glm.fit <- glm(resp ~ factor(dose) + gender, data = dat,
                family = poisson(link = "log"))

  mu.hat <- coef(glm.fit)
  S.hat <- vcov(glm.fit)
  S.hat <- diag(exp(mu.hat))%*%S.hat%*%diag(exp(mu.hat))
  mu.hat <- exp(mu.hat)

  S.hat <- S.hat[2:n.doses, 2:n.doses]
  mu.hat <- mu.hat[2:n.doses] - 1
  dose.in <- doses[2:n.doses]

  set.seed(51)
  mod.DF <- DoseFinding::MCPMod(dose.in, mu.hat, models = models,
                               S = S.hat, type = "general",
                  placAdj = T, Delta = 0.3)

  set.seed(51)
  mod.pack <- MCPModGen(family = "poisson", link = "risk ratio",
                       dose = "dose", resp = "resp", addCovars = ~ gender,
                       data = dat, models = models, Delta = 0.3, returnS = T)


  expect_equal(mod.DF$MCTtest, mod.pack$MCPMod$MCTtest)

})




