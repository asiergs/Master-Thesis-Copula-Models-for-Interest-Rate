
mutate <- dplyr::mutate
count <- dplyr::count
filter <- dplyr::filter
group_by <- dplyr::group_by
ungroup <- dplyr::ungroup
wday <- lubridate::wday

diff_narm <- function(ts, lag = 1){na.omit(diff(ts, lag = lag))}

unscale <- function(matrix, scale, center){
  t(apply(matrix,1,function(matrix) matrix * scale + center ))
}

ARIMA_GARCH_test <- function(param, ts, filename, type.model, max_order = 3,
                             start = 1, save_every = 1000){
  
  n <- nrow(param)
  
  if (start != 1){
    results <- as.matrix(read.csv(filename))
  } else{
    results <- matrix(NA, nrow = n, ncol = 7)
    colnames(results) <- c("N","max_p_value", "akaike",
                           "bic","shibata", "hannan_quinn","max_estim")
  }
  
  # do not show warnings
  defaultW <- getOption("warn") 
  options(warn = -1) 
  
  for (i in start:n){
    par <- param[i,]
    par_fixed <- par[,!names(par) %in% c("distr", "mean")]
    par_fixed <- par_fixed[, as.vector(par_fixed[1,] != 1)]
    fixed.pars = as.list(par_fixed)
    mean <- par[,"mean"][[1]]
    distr <- par[,"distr"][[1]]
    
    Garch <- ugarchspec(variance.model = list(model = type.model,
                                              garchOrder = c(max_order,
                                                             max_order,
                                                             max_order)),
                        mean.model = list(armaOrder = c(max_order, max_order),
                                          include.mean = mean),
                        distribution.model = distr, fixed.pars = fixed.pars)
    
    model <- tryCatch(withTimeout(ugarchfit(spec = Garch, data = ts),
                                  timeout = 30, onTimeout = c("error")),
                      error = function(e) 1)
    
    if (typeof(model)!="double") {
      # omega is not checked since it must be there for the model to converge
      model_info <- model@fit$matcoef[,c("Pr(>|t|)"," Estimate")]
      model_info <- model_info[rownames(model_info)!="omega",]
      p.value.max <- max(model_info[,c("Pr(>|t|)")],0, na.rm = TRUE)
      model_info <- model_info[rownames(model_info)!="skew",]
      model_info <- model_info[rownames(model_info)!="shape",]
      max_estim <- max(abs(c(model_info[," Estimate"],0)))
      
      infocriteria <- tryCatch(infocriteria(model), error = function(e) 1)
      
      if (length(infocriteria) == 1){
        akaike <- Inf
        bic <- Inf
        shibata <- Inf
        hannan_quinn <- Inf
      } else {
        akaike <- infocriteria[1]
        bic <- infocriteria[2]
        shibata <- infocriteria[3]
        hannan_quinn <- infocriteria[4]
      }
      results[i,] <- c(i, p.value.max, akaike, bic, shibata, hannan_quinn,
                       max_estim)
    }
    if (i%%save_every == 0) {
      write.csv(as.data.frame(results),filename, row.names = FALSE)
      print(paste0(i, " out of ", n, " - progress ", round(i/n*100, 2), "%"))
      print(Sys.time())
    }
  }
  # show warnings again
  options(warn = defaultW)
  
  write.csv(as.data.frame(results),filename, row.names = FALSE)
  
  as.data.frame(results)
}

testCopula <- function(u, method = "mult", opt.meth = "L-BFGS-B"){
  n <- length(u[,1])
  # Normal copula
  normc <- fitCopula(copula = normalCopula(), data = u, method = "irho")
  resultnorm <- gofCopula(copula = normalCopula(param = normc@estimate), x = u,
                          simulation = method)
  BICnorm <- 1*log(n) - 2*normc@loglik
  AICcnorm <- 2*1 - 2*normc@loglik + (2*1^2+2*1)/(n-1-1)
  # tStudent copula
  tc <- fitCopula(copula = tCopula(), data = u, optim.method = opt.meth)
  # if error change optim.method as opt.meth = "Nelder-Mead"
  tc <- fitCopula(copula = tCopula(df = round(tc@estimate[2]), df.fixed = TRUE),
                  data = u)
  # library only implemented for df = integer
  df = tc@copula@parameters[[2]]
  resultt <- gofCopula(copula = tCopula(param = tc@estimate, df = df,
                                        df.fixed = TRUE),
                       x = u, simulation = method)
  BICt<- 2*log(n) - 2*tc@loglik
  AICct <- 2*2 - 2*tc@loglik + (2*2^2+2*2)/(n-2-1)
  # AMH copula
  amhc <- fitCopula(copula = amhCopula(), data = u)
  resultamh <- gofCopula(copula = amhCopula(param = amhc@estimate), x = u,
                         simulation = method)
  BICamh <- 1*log(n) - 2*amhc@loglik
  AICcamh <- 2*1 - 2*amhc@loglik + (2*1^2+2*1)/(n-1-1)
  # Clayton copula
  cc <- fitCopula(copula = claytonCopula(), data = u)
  resultc <- gofCopula(copula = claytonCopula(param = cc@estimate), x = u,
                       simulation = method)
  BICc <- 1*log(n) - 2*cc@loglik
  AICcc <- 2*1 - 2*cc@loglik + (2*1^2+2*1)/(n-1-1)
  # Frank copula
  fc <- fitCopula(copula = frankCopula(), data = u)
  resultf <- gofCopula(copula = frankCopula(param = fc@estimate), x = u,
                       simulation = method)
  BICf <- 1*log(n) - 2*fc@loglik
  AICcf <- 2*1 - 2*fc@loglik + (2*1^2+2*1)/(n-1-1)
  # Gumbel-Hougaard copula
  ghc <- fitCopula(copula = gumbelCopula(), data = u)
  resultgh <- gofCopula(copula = gumbelCopula(param = ghc@estimate), x = u,
                        simulation = method)
  BICgh <- 1*log(n) - 2*ghc@loglik
  AICcgh <- 2*1 - 2*ghc@loglik + (2*1^2+2*1)/(n-1-1)
  # Joe copula
  jc <- fitCopula(copula = joeCopula(), data = u)
  resultj <- gofCopula(copula = joeCopula(param = jc@estimate), x = u,
                       simulation = method)
  BICj <- 1*log(n) - 2*jc@loglik
  AICcj <- 2*1 - 2*jc@loglik + (2*1^2+2*1)/(n-1-1)
  data.frame("copula" = c("Normal",paste0("t-Student"," (df=",df,")"),"AMH",
                          "Clayton", "Frank","Gumbel-Hougaard","Joe"),
             "log-likelihood" = c(normc@loglik,tc@loglik,amhc@loglik,cc@loglik,
                                  fc@loglik,ghc@loglik,jc@loglik),
             "parameter" = c(normc@estimate,tc@estimate,amhc@estimate,cc@estimate,
                             fc@estimate,ghc@estimate,jc@estimate),
             "s.e." = sqrt(c(normc@var.est,tc@var.est,amhc@var.est,cc@var.est,
                             fc@var.est,ghc@var.est,jc@var.est)),
             "gof pvalue" = c(resultnorm$p.value,resultt$p.value,
                              resultamh$p.value,resultc$p.value,resultf$p.value,
                              resultgh$p.value,resultj$p.value),
             "BIC" = c(BICnorm,BICt,BICamh,BICc,BICf,BICgh,BICj),
             "AICc" = c(AICcnorm,AICct,AICcamh,AICcc,AICcf,AICcgh,AICcj))
}

# a function to fit a GARCH model based in the notation of parameters proposed

fit_GARCH <- function(ts, par, type.model, max_order = 3){
  par_fixed <- par[,!names(par) %in% c("distr", "mean")]
  par_fixed <- par_fixed[, as.vector(par_fixed[1,] != 1)]
  fixed.pars = as.list(par_fixed)
  mean <- par[,"mean"][[1]]
  distr <- par[,"distr"][[1]]
  
  sGarch <- ugarchspec(variance.model = list(model = type.model,
                                             garchOrder = c(max_order,max_order,max_order)),
                       mean.model = list(armaOrder = c(max_order, max_order),
                                         include.mean = mean),
                       distribution.model = distr, fixed.pars = fixed.pars)
  ugarchfit(spec = sGarch, data = dtsPC)
}


fit_GARCH_spec <- function(par, type.model, max_order = 3){
  par_fixed <- par[,!names(par) %in% c("distr", "mean")]
  par_fixed <- par_fixed[, as.vector(par_fixed[1,] != 1)]
  fixed.pars = as.list(par_fixed)
  mean <- par[,"mean"][[1]]
  distr <- par[,"distr"][[1]]
  ugarchspec(variance.model = list(model = type.model,
                                   garchOrder = c(max_order,max_order,max_order)),
             mean.model = list(armaOrder = c(max_order, max_order),
                                   include.mean = mean),
             distribution.model = distr, fixed.pars = fixed.pars)
}

# convert the fitted model to a spec

fit_to_spec <- function(model){
  model.type = model@model$modeldesc$vmodel
  fixed.pars <- as.list(model@model$pars[!is.na(model@model$pars[,"LB"]),"Level"])
  if (!any(names(fixed.pars) == c("gamma1"))){
    fixed.pars$gamma1 <- 0
    fixed.pars$gamma2 <- 0
    fixed.pars$gamma3 <- 0
  }
  maxOrder <- model@model$maxOrder
  distr <- model@model$modeldesc$distribution
  Garch <- ugarchspec(variance.model = list(model = type.model,
                                            garchOrder = c(maxOrder,
                                                           maxOrder,
                                                           maxOrder)),
                      mean.model = list(armaOrder = c(maxOrder, maxOrder)),
                      distribution.model = distr, fixed.pars = fixed.pars)
  Garch
}

# simulate paths providing the residuals as uniform values (p)

sim_ugarch <- function(fit, u){
  # u = random uniform (0.0 to 1.0) vector of residuals in order
  n <- length(u)
  newdates <-tail(fit@model$modeldata$index,1) + 7*seq(1,n)
  coef <- fit@fit$coef
  distr <- fit@model$modeldesc$distribution
  if (is.na(coef["skew"])) skew <- 1 else skew <- coef["skew"]
  if (is.na(coef["shape"])) shape <- 5 else shape <- coef["shape"]
  sres <- matrix(data = qdist(distribution = distr, p = u, shape = shape,
                              skew = skew),
                 nrow = n, ncol = 1)
  ts <- ugarchsim(fit, n.sim = n, m.sim = 1, custom.dist = list(name="custom",
                                                            distfit = sres))
  xts(ts@simulation$seriesSim, order.by = newdates)
}

portfolio_value_calc <- function(portfolio, shock_swaps, times, shock_stock,
                           duration = FALSE){
  # portfolio: list with bond and stock info
  # rates: swap rates in percentage
  # times: the times corresponding with the rates
  # shock_stock: the final value as fraction
  # duration: if duration is wanted as output
  rates <- approx(times, shock_swaps, 1:30)
  value_bond <- 0
  value_stock <- 0
  dur <- 0
  for (i in 1:length(portfolio)) {
    if (names(portfolio[i])=="bond"){
      bond <- portfolio[[i]]
      cashflow <- c(rep(bond$cupon, bond$mat - 1), 1 + bond$cupon)
      PV_rates <- rates$y[1:bond$mat]
      PV_cashflows <- cashflow/(1+PV_rates/100)^(1:bond$mat)
      value_bond <- value_bond + sum(PV_cashflows)*bond$qty
      dur <- dur + sum(PV_cashflows*1:bond$mat)*bond$qty
    } else value_stock <- value_stock + portfolio[[i]]$price*shock_stock
  }
  dur <- dur/value_bond
  if (duration == TRUE) {
    c("bonds" = value_bond, "stocks" = value_stock, "duration" = dur)
  } else c("bonds" = value_bond, "stocks" = value_stock)
}
