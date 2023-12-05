
# clear work space
rm(list = ls())

# 0. Packages and functions ----

library("ggplot2")
library("tidyverse")
library("dplyr")
library("fpp2")
library("cowplot")
library("tseries")
library("copula")
library("zoo")
library("rugarch")
library("xts")
library("goftest")
#install.packages("R.utils")
library("R.utils")
#install.packages("xlsx")
library("xlsx")
# load all user functions
source("functions.R")
#install.packages("latticeExtra")
library("latticeExtra")
#install.packages("corrgram")
library("corrgram")
library("plotly")

old.par <- par(mar = c(0, 0, 0, 0))

# 1. Import and clean data ----

## Import from Github ----

# Define all the bonds and the urls
url <- "https://raw.githubusercontent.com/asiergs/financial_data/main/EUROIS.csv"

## 1.1 Import rates ----
swaps_df <- read.csv(url)
swaps_df2 <- na.omit(swaps_df)

## 1.2 Clean rates data----

# Order the swaps by time
names <- names(swaps_df)
times <- c(0, as.numeric(gsub("[A-Z]", "",names[-1])))

for (i in 2:length(names)){
  if (substring(names[i],nchar(names[i])) != "Y"){
    names[i] <- paste0(names[i],"Y")
  } 
}

names(swaps_df2) <- names

names(times) <- names
times <- sort(times)

swaps_df2 <- swaps_df2[,names(times)]

times <- times[-1]

# Expand the data to work with the data
swaps_df3 <- swaps_df2 %>% pivot_longer(cols = colnames(swaps_df2)[-1],
                                        values_to = "rate")

# Get the time for each swipe and convert date char to date format
swaps_df3 <- swaps_df3 %>% mutate(date = as.Date(Code, format = "%d/%m/%Y"),
                                  time = as.numeric(gsub("[A-Z]", "",name)))

# Get only the columns of interest
swaps_df3 <- swaps_df3 %>% select(date,name,time,rate)

# Check there is no missing data (shall not be with the na.omit())
swaps_df3 %>% group_by(date) %>% count() %>% summarise(max = max(n),
                                                       min = min(n))

## Convert to weekly

start <- min(swaps_df3$date)
start_weekday <- wday(start, week_start = 1)

# Extract the last day from week only
swaps_df3 <- swaps_df3 %>%
  mutate(weekday = wday(date, week_start = 1),
         days_from_start = date - start + start_weekday,
         week = as.numeric(trunc(days_from_start/7))) %>%
  group_by(week) %>% mutate(max_wday = max(weekday)) %>%
  filter(weekday == max_wday) %>% ungroup()

# Convert to matrix

swaps_df4 <- swaps_df3 %>% pivot_wider(id_cols = date, values_from = rate)
dates <- swaps_df4$date

# match all columns in same time

swaps_matrix <- as.matrix(swaps_df4[,-1])
rownames(swaps_matrix) <- as.character(swaps_df4$date)

## 1.3 Import EURSTOXX600 and USDEUR data and join it ----

url <- "https://raw.githubusercontent.com/asiergs/financial_data/main/DJSTOXX.csv"

DJSTOXX <- read.csv(url) %>% mutate(date = as.Date(Date, format = "%d/%m/%Y")) %>%
  select(date, DJSTOXX)

url <- "https://raw.githubusercontent.com/asiergs/financial_data/main/USDEUR.csv"

USDEUR <- read.csv(url) %>% mutate(date = as.Date(Date, format = "%d/%m/%Y")) %>%
  select(date, USDEUR)

market_df <- swaps_df4 %>% left_join(DJSTOXX, join_by(date)) %>%
  left_join(USDEUR, join_by(date))

# 2. PCA & Factorial analysis----

scaled_swaps <- scale(swaps_matrix)

center <- attr(scaled_swaps,"scaled:center")
scale <- attr(scaled_swaps,"scaled:scale")

## 2.1 Get the values of the Principal Components ----

correl_matrix <- cov(scaled_swaps)
eigen_swaps <- eigen(correl_matrix)

PC <- scaled_swaps %*% eigen_swaps$vectors
var_cum_explained <- cumsum(eigen_swaps$values)/sum(eigen_swaps$values)
write.xlsx(as.data.frame(var_cum_explained), file = "tables/var_explained.xlsx")

# undo principal components
PC %*% t(eigen_swaps$vectors)

PC_mat <- t(eigen_swaps$vectors)
colnames(PC_mat) <- names(times)
rownames(PC_mat) <- paste0("PC",1:15)

write.xlsx(as.data.frame(PC_mat), file = "tables/PC_matrix.xlsx")

## 2.2 Fit of the Principal Components to curve ----
dim <- length(times)
curve <- matrix(0,ncol = dim, nrow = dim)

swaps_red_list <- list()
for (i in 1:dim){
  # The not selected PC are set to 0
  PC_red <- PC
  PC_red[,-(1:i)] <- 0
  # The PC are undone and unscaled
  swaps_red_list[[i]]  <- unscale(PC_red %*% t(eigen_swaps$vectors),scale,center)
}

date <- as.character(swaps_df4$date[1220])
ylim <- c(min(swaps_matrix[date,])-0.25,max(swaps_matrix[date,])+0.5)

# Graphs of how PC affect the curve
PC_show <- rep(0,15)

sample_curve <- unscale(PC_show %*% t(eigen_swaps$vectors),scale,center)

for(i in 1:3){
  svg(file = paste0("plots/PCA_",i,"_shock.svg"), width = 5, height = 5)
  par(mfrow=c(1,1), oma = c(0,0,0,0), mar = c(4,4,2,1))
  vector <- c(rep(0,i-1),1,rep(0,15-i))
  plot(times,sample_curve, type="p", lty = 1, pch = 19,
       lwd = 2, ylim = c(0,5), ylab = "swap rate", main = paste0("PC ",i," shock"),
       xlab = "years")
  lines(times,sample_curve, type="l", lty = 1,
        lwd = 2)
  sample_curve_shock <- unscale(vector %*% t(eigen_swaps$vectors),scale,center)
  lines(times,sample_curve_shock, type = "l", lty = 2, lwd = 2, col = "#6A7434")
  sample_curve_shock <- unscale(-vector %*% t(eigen_swaps$vectors),scale,center)
  lines(times,sample_curve_shock, type = "l", lty = 2, lwd = 2, col = "orange")
  legend("bottomright",legend = c("Yield curve", "Positive shock", "Negative shock"),
         lty = c(1,2,2), lwd = c(2,2,2), ncol=1, col= c("black","#6A7434","orange"))
  dev.off()
}


#### INPUT to change graph
n <- 9
#n <- dim
#### INPUT to change graph

svg(file = 'plots/PCA2.svg', width = 10, height = 10)
par(mfrow=c(ceiling(n/3),3), oma = c(2,2,2,2), mar = c(1.5,1.5,4,1.5))
for (i in 1:n) {
  plot(times,swaps_matrix[date,1:dim], type="p", lty = 1, pch = 19,
       lwd = 2, ylim = ylim, ylab = "swap", main = paste0(i," PC"), xlab = "years")
  lines(times,swaps_matrix[date,1:dim], type="l", lty = 1,
        lwd = 2)
  lines(times,swaps_red_list[[i]][date,1:dim], type = "l", lty = 2, lwd = 2, col = "red")
  legend("bottomright",legend = c("Real", "PC approx"), lty = 1:2, lwd = c(2,2), ncol=2,
         col= c("black","red"))
}
dev.off()

## 2.3 PC components interpretation ----
svg(file = 'plots/PCA_inerpret.svg', width = 15/1.2, height = 5/1.2)
par(mfrow=c(1,1), mar = c(5,5,1,1),oma = c(0,0,0,0))
plot(times, eigen_swaps$vectors[,1], type = "l", col = "#6A7434",
     ylim = c(-1,1), lwd = 2, ylab =  "weight", xlab = "years")
lines(times, eigen_swaps$vectors[,2], type = "l", col = "#E0FF00", lwd = 2)
lines(times, eigen_swaps$vectors[,3], type = "l", col = "orange", lwd = 2)
lines(times, eigen_swaps$vectors[,4], type = "l", col = "black", lwd = 2)
legend("bottom", legend = paste0("PC",1:4), lwd = rep(2,4),
       col = c("#6A7434","#E0FF00","orange","black"), ncol = 4)
dev.off()

colors <- rainbow(dim)

svg(file = 'plots/PCA_time.svg', width = 10, height = 6)
par(mfrow=c(1,1), oma = c(1,1,1,1), mar = c(1,1,1,1))
plot(dates,PC[,1],type="l",col=colors[1], ylim = c(-8,6), xlab = "time", ylab = "")
for (i in 2:dim) {
  lines(dates,PC[,i],col=colors[i])
}

legend("bottom",legend = paste("PC",1:dim), col = colors,
       lwd = rep(1,12),ncol=6)
dev.off()

# 3. Time series estimation ----

colnames(PC) <- paste0("PC",seq(1,15))
tsPCX <- xts(PC, order.by = swaps_df4$date)
tsPC <- tsPCX

tsPC1 <- tsPCX[,1]
tsPC2 <- tsPCX[,2]
tsPC3 <- tsPCX[,3]


par(old.par)
svg(file = "plots/ts_PC.svg", width = 16.5/1.6, height = 5/1.6)
par(mfrow=c(1,1), oma = c(1,1,1,1), mar = c(4,1,1,1))
plot(dates, tsPC1, type = "l", col = "#6A7434", xlab = "time", ylab = "")
lines(dates, tsPC2, type = "l",
      col = "#202426")
lines(dates, tsPC3, type = "l",
      col = "orange")
legend("bottomright",legend = paste("PC",1:3), col = c("#6A7434", "#202426", "orange"),
       lwd = rep(1,3),ncol=3)
dev.off()

## 3.1 PC1 ----

### 3.1.1 Stationarity check ----

tsPC <- tsPCX[,1]

svg(file = 'plots/PC1_hist.svg', width = 7, height = 4)
par(mfrow=c(1,1), mar = c(5,5,1,1),oma = c(0,0,0,0))
plot(tsPC, type = "l", main = "PC1 (level)", col = "blue")
dev.off()

ggdraw() +
  draw_plot(ggAcf(tsPC), x = 0, y = .5, width = .5, height = .5) +
  draw_plot(ggPacf(tsPC), x = .5, y = .5, width = .5, height = .5) +
  draw_plot(autoplot(tsPC), x = 0, y = 0, width = 1, height = 0.5)

adf.test(tsPC)
pp.test(tsPC)
kpss.test(tsPC)

dtsPC <- diff_narm(tsPC)

ggdraw() +
  draw_plot(ggAcf(dtsPC), x = 0, y = .5, width = .5, height = .5) +
  draw_plot(ggPacf(dtsPC), x = .5, y = .5, width = .5, height = .5) +
  draw_plot(autoplot(dtsPC), x = 0, y = 0, width = 1, height = 0.5)

station <- list()
station[[1]] <- adf.test(dtsPC)
station[[2]] <- pp.test(dtsPC)
station[[3]] <- kpss.test(dtsPC)

statistic <- c()
pvalue <- c()
alternative <- c()

for (i in 1:3) statistic <- c(statistic, station[[i]]$statistic)
for (i in 1:3) pvalue <- c(pvalue, station[[i]]$p.value)
for (i in 1:3) alternative <- c(alternative, station[[i]]$alternative)
alternative[3] <- "non stationary"

station <- data.frame("statistic" = statistic, "pvalue" = pvalue, "alternative hypothesis"=
                        alternative, row.names = c("Dickey-Fuller", "Phillips-Perron", "KPSS"))

write.xlsx(station, file = "tables/stationarityPC1.xlsx")

### 3.1.2 Estimating the model ARMA GARCH ----

#### 3.1.2.1 ARIMA sGARCH model fitting ----

type.model <- "sGARCH"

filename_results_dtsPC <- paste0("results_",type.model,"_dts",
                                        colnames(tsPC),".csv")
filename_param_dtsPC <- paste0("param_",type.model,"_dts",
                                      colnames(tsPC),".csv")

# tests multiple ARMA order GARCH order models
results <- tryCatch(read.csv(filename_results_dtsPC), error = function(e) 1,
                    warning = function(w) 1)
# check if there is a file with the results and otherwise, calculate them
if (typeof(results) == "double"){
  comb_full <- crossing(ar1 = 0:1, ar2 = 0:1, ar3 = 0:1,
                        ma1 = 0:1, ma2 = 0:1, ma3 = 0:1,
                        alpha1 = 0:1, alpha2 = 0:1, alpha3 = 0:1,
                        beta1 = 0:1, beta2 = 0:1, beta3 = 0:1,
                        mean = c(TRUE, FALSE),
                        distr = c("std","norm","ged","nig"))
  write.csv(comb_full, filename_param_dtsPC, row.names = FALSE)
  results <- ARIMA_GARCH_test(param = comb_full, ts = dtsPC,
                              filename = filename_results_dtsPC,
                              type.model = type.model)
  param <- comb_full
} else param <- read.csv(filename_param_dtsPC)

# filter only best models
alpha <- 1 - 0.995

valid_models <- results %>% mutate(rankAIC = rank(-akaike),
                                   rankBIC = rank(-bic),
                                   rankSHIB = rank(-shibata),
                                   rankHANQ = rank(-hannan_quinn),
                                   score = rankAIC + rankBIC + rankSHIB +
                                     rankHANQ) %>%
  filter(max_p_value < alpha, max_estim < 1) %>% arrange(desc(score)) %>%
  select(-rankAIC, -rankBIC, -rankSHIB, -rankHANQ)

# see models selected
param_valid <- param[valid_models$N,]
valid_models

valid_models$distr <- param_valid$distr
param_valid

# see one particular models of those selected
(model_sGarch_PC1 <- fit_GARCH(dtsPC, param[valid_models$N[95],], type.model = type.model))

#### 3.1.2.2 ARIMA eGARCH model fitting ----

type.model <- "eGARCH"

filename_results_dtsPC <- paste0("results_",type.model,"_dts",
                                 colnames(tsPC),".csv")
filename_param_dtsPC <- paste0("param_",type.model,"_dts",
                               colnames(tsPC),".csv")

# tests multiple ARMA order GARCH order models
results <- tryCatch(read.csv(filename_results_dtsPC), error = function(e) 1,
                    warning = function(w) 1)
# check if there is a file with the results and otherwise, calculate them
if (typeof(results) == "double"){
  comb_full <- crossing(ar1 = 0:1, ar2 = 0:1, ar3 = 0:1,
                        ma1 = 0:1, ma2 = 0:1, ma3 = 0:1,
                        alpha1 = 0:1, alpha2 = 0:1, alpha3 = 0:1,
                        beta1 = 0:1, beta2 = 0:1, beta3 = 0:1,
                        mean = c(TRUE, FALSE),
                        distr = c("std","norm","ged","nig"))
  write.csv(comb_full, filename_param_dtsPC, row.names = FALSE)
  results <- ARIMA_GARCH_test(param = comb_full, ts = dtsPC,
                              filename = filename_results_dtsPC,
                              type.model = type.model)
  param <- comb_full
} else param <- read.csv(filename_param_dtsPC)

# filter only best models
alpha <- 1 - 0.995

valid_models <- results %>% mutate(rankAIC = rank(-akaike),
                                   rankBIC = rank(-bic),
                                   rankSHIB = rank(-shibata),
                                   rankHANQ = rank(-hannan_quinn),
                                   score = rankAIC + rankBIC + rankSHIB +
                                     rankHANQ) %>%
  filter(max_p_value < alpha, max_estim < 1) %>% arrange(desc(score)) %>%
  select(-rankAIC, -rankBIC, -rankSHIB, -rankHANQ)

# see models selected
param_valid <- param[valid_models$N,]
valid_models

valid_models$distr <- param_valid$distr

# see one particular models of those selected
model_eGarch_PC1 <- fit_GARCH(dtsPC, param[valid_models$N[12],], type.model = type.model)

#### 3.1.2.3 ARIMA gjrGARCH model fitting ---- 

type.model <- "gjrGARCH"

filename_results_dtsPC <- paste0("results_",type.model,"_dts",
                                 colnames(tsPC),".csv")
filename_param_dtsPC <- paste0("param_",type.model,"_dts",
                               colnames(tsPC),".csv")

# tests multiple ARMA order GARCH order models
results <- tryCatch(read.csv(filename_results_dtsPC), error = function(e) 1,
                    warning = function(w) 1)
# check if there is a file with the results and otherwise, calculate them
if (typeof(results) == "double"){
  comb_full <- crossing(ar1 = 0:1, ar2 = 0:1, ar3 = 0:1,
                        ma1 = 0:1, ma2 = 0:1, ma3 = 0:1,
                        alpha1 = 0:1, alpha2 = 0:1, alpha3 = 0:1,
                        beta1 = 0:1, beta2 = 0:1, beta3 = 0:1,
                        mean = c(TRUE, FALSE),
                        distr = c("std","norm","ged","nig"))
  write.csv(comb_full, filename_param_dtsPC, row.names = FALSE)
  results <- ARIMA_GARCH_test(param = comb_full, ts = dtsPC,
                              filename = filename_results_dtsPC,
                              type.model = type.model)
  param <- comb_full
} else param <- read.csv(filename_param_dtsPC)

# filter only best models
alpha <- 1 - 0.995

valid_models <- results %>% mutate(rankAIC = rank(-akaike),
                                   rankBIC = rank(-bic),
                                   rankSHIB = rank(-shibata),
                                   rankHANQ = rank(-hannan_quinn),
                                   score = rankAIC + rankBIC + rankSHIB +
                                     rankHANQ) %>%
  filter(max_p_value < alpha, max_estim < 0.9) %>% arrange(desc(score)) %>%
  select(-rankAIC, -rankBIC, -rankSHIB, -rankHANQ)

# see models selected
param_valid <- param[valid_models$N,]
valid_models

valid_models$distr <- param_valid$distr

# see one particular models of those selected
model_gjrGarch_PC1 <- fit_GARCH(dtsPC, param[valid_models$N[3],], type.model = type.model)

### 3.1.3 VaR backtesting ####

model_dtsPC1 <- model_sGarch_PC1
model <- model_sGarch_PC1

spec <- fit_GARCH_spec(param[valid_models$N[9],], type.model = type.model)

modelroll <- ugarchroll(
  spec = spec, data=dtsPC, n.ahead = 1, forecast.length = 1119, refit.every = 100,
  refit.window = c("recursive"), solver = "hybrid", calculate.VaR = TRUE,
  VaR.alpha = c(0.01,0.05,0.95,0.99))

VaR001 <- modelroll@forecast$VaR[,"alpha(1%)"]
VaR005 <- modelroll@forecast$VaR[,"alpha(5%)"]
VaR095 <- modelroll@forecast$VaR[,"alpha(95%)"]
VaR099 <- modelroll@forecast$VaR[,"alpha(99%)"]
real <- modelroll@forecast$VaR[,"realized"]

mean(real<VaR001)
mean(real<VaR005)
mean(real<VaR095)
mean(real<VaR099)

VaR <- list()
VaR[[1]] <- VaRTest(alpha = 0.01, real, VaR001, conf.level = 0.995)
VaR[[2]] <- VaRTest(alpha = 0.05, real, VaR005, conf.level = 0.995)
VaR[[3]] <- VaRTest(alpha = 0.05, -real, -VaR095, conf.level = 0.995)
VaR[[4]] <- VaRTest(alpha = 0.01, -real, -VaR099, conf.level = 0.995)

VaRDur <- list()
VaRDur[[1]] <- VaRDurTest(alpha = 0.01, real, VaR001, conf.level = 0.995)
VaRDur[[2]] <- VaRDurTest(alpha = 0.05, real, VaR005, conf.level = 0.995)
VaRDur[[3]] <- VaRDurTest(alpha = 0.01, -real, -VaR095, conf.level = 0.995)
VaRDur[[4]] <- VaRDurTest(alpha = 0.05, -real, -VaR099, conf.level = 0.995)

DAC <- list()
DAC[[1]] <- DACTest(modelroll@forecast$density$Mu, real, test = "PT", conf.level = 0.95)
DAC[[2]] <- DACTest(modelroll@forecast$density$Mu, real, test = "AG", conf.level = 0.95)

results <- matrix(NA,ncol = 3, nrow = 4)
colnames(results) <- c("unconditional", "conditional", "duration")
rownames(results) <- c(0.01,0.05,0.95,0.99)
for (i in 1:4){
  results[i,] <- c(VaR[[i]]$uc.LRp, VaR[[i]]$cc.LRp, VaRDur[[i]]$LRp) 
}

VaR_down <- quantile(model, probs = c(0.01, 0.05))
VaR[[1]] <- VaRTest(alpha = 0.01, dtsPC, VaR_down[,1], conf.level = 0.95)
VaR[[2]] <- VaRTest(alpha = 0.05, dtsPC, VaR_down[,2], conf.level = 0.95)
VaRDur[[1]] <- VaRDurTest(alpha = 0.01, dtsPC, VaR_down[,1], conf.level = 0.95)
VaRDur[[2]] <- VaRDurTest(alpha = 0.05, dtsPC, VaR_down[,2], conf.level = 0.95)

VaR_up <- -quantile(model, probs = c(0.95, 0.99))
VaR[[3]] <- VaRTest(alpha = 0.05, -dtsPC, VaR_up[,1], conf.level = 0.95)
VaR[[4]] <- VaRTest(alpha = 0.01, -dtsPC, VaR_up[,2], conf.level = 0.95)
VaRDur[[3]] <- VaRDurTest(alpha = 0.01, -dtsPC, VaR_up[,1], conf.level = 0.95)
VaRDur[[4]] <- VaRDurTest(alpha = 0.05, -dtsPC, VaR_up[,2], conf.level = 0.95)

DAC[[1]] <- DACTest(fitted(model), dtsPC, test = "PT", conf.level = 0.95)
DAC[[2]] <- DACTest(fitted(model), dtsPC, test = "AG", conf.level = 0.95)

for (i in 1:4){
  results[i,] <- c(VaR[[i]]$uc.LRp, VaR[[i]]$cc.LRp, VaRDur[[i]]$LRp) 
}


write.xlsx(as.data.frame(results), file = "tables/backtestPC1.xlsx")

## 3.2 PC2 ----

### 3.2.1 Stationarity check ----

tsPC <- tsPCX[,2]

svg(file = 'plots/PC2_hist.svg', width = 7, height = 4)
par(mfrow=c(1,1), mar = c(5,5,1,1),oma = c(0,0,0,0))
plot(tsPC, type = "l", main = "PC1 (level)", col = "blue")
dev.off()

ggdraw() +
  draw_plot(ggAcf(tsPC), x = 0, y = .5, width = .5, height = .5) +
  draw_plot(ggPacf(tsPC), x = .5, y = .5, width = .5, height = .5) +
  draw_plot(autoplot(tsPC), x = 0, y = 0, width = 1, height = 0.5)

adf.test(tsPC)
pp.test(tsPC)
kpss.test(tsPC)

dtsPC <- diff_narm(tsPC)

ggdraw() +
  draw_plot(ggAcf(dtsPC), x = 0, y = .5, width = .5, height = .5) +
  draw_plot(ggPacf(dtsPC), x = .5, y = .5, width = .5, height = .5) +
  draw_plot(autoplot(dtsPC), x = 0, y = 0, width = 1, height = 0.5)

station <- list()
station[[1]] <- adf.test(dtsPC)
station[[2]] <- pp.test(dtsPC)
station[[3]] <- kpss.test(dtsPC)

statistic <- c()
pvalue <- c()
alternative <- c()

for (i in 1:3) statistic <- c(statistic, station[[i]]$statistic)
for (i in 1:3) pvalue <- c(pvalue, station[[i]]$p.value)
for (i in 1:3) alternative <- c(alternative, station[[i]]$alternative)
alternative[3] <- "non stationary"

station <- data.frame("statistic" = statistic, "pvalue" = pvalue, "alternative hypothesis"=
                     alternative, row.names = c("Dickey-Fuller", "Phillips-Perron", "KPSS"))

write.xlsx(station, file = "tables/stationarityPC2.xlsx")

### 3.2.2 Estimating the model ARMA GARCH ----

#### 3.2.2.1 ARIMA sGARCH model fitting ----

type.model <- "sGARCH"

filename_results_dtsPC <- paste0("results_",type.model,"_dts",
                                 colnames(tsPC),".csv")
filename_param_dtsPC <- paste0("param_",type.model,"_dts",
                               colnames(tsPC),".csv")

# tests multiple ARMA order GARCH order models
results <- tryCatch(read.csv(filename_results_dtsPC), error = function(e) 1,
                    warning = function(w) 1)
# check if there is a file with the results and otherwise, calculate them
if (typeof(results) == "double"){
  comb_full <- crossing(ar1 = 0:1, ar2 = 0:1, ar3 = 0:1,
                        ma1 = 0:1, ma2 = 0:1, ma3 = 0:1,
                        alpha1 = 0:1, alpha2 = 0:1, alpha3 = 0:1,
                        beta1 = 0:1, beta2 = 0:1, beta3 = 0:1,
                        mean = c(TRUE, FALSE),
                        distr = c("std","norm","ged","nig"))
  write.csv(comb_full, filename_param_dtsPC, row.names = FALSE)
  results <- ARIMA_GARCH_test(param = comb_full, ts = dtsPC,
                              filename = filename_results_dtsPC,
                              type.model = type.model)
  param <- comb_full
} else param <- read.csv(filename_param_dtsPC)

# filter only best models
alpha <- 1 - 0.995

valid_models <- results %>% mutate(rankAIC = rank(-akaike),
                                   rankBIC = rank(-bic),
                                   rankSHIB = rank(-shibata),
                                   rankHANQ = rank(-hannan_quinn),
                                   score = rankAIC + rankBIC + rankSHIB +
                                     rankHANQ) %>%
  filter(max_p_value < alpha, max_estim < 1) %>% arrange(desc(score)) %>%
  select(-rankAIC, -rankBIC, -rankSHIB, -rankHANQ)

# see models selected
param_valid <- param[valid_models$N,]
valid_models

valid_models$distr <- param_valid$distr

# see one particular models of those selected
(model_sGarch_PC2 <- fit_GARCH(dtsPC, param[valid_models$N[12],], type.model = type.model))


#### 3.2.2.2 ARIMA eGARCH model fitting ----

type.model <- "eGARCH"

filename_results_dtsPC <- paste0("results_",type.model,"_dts",
                                 colnames(tsPC),".csv")
filename_param_dtsPC <- paste0("param_",type.model,"_dts",
                               colnames(tsPC),".csv")

# tests multiple ARMA order GARCH order models
results <- tryCatch(read.csv(filename_results_dtsPC), error = function(e) 1,
                    warning = function(w) 1)
# check if there is a file with the results and otherwise, calculate them
if (typeof(results) == "double"){
  comb_full <- crossing(ar1 = 0:1, ar2 = 0:1, ar3 = 0:1,
                        ma1 = 0:1, ma2 = 0:1, ma3 = 0:1,
                        alpha1 = 0:1, alpha2 = 0:1, alpha3 = 0:1,
                        beta1 = 0:1, beta2 = 0:1, beta3 = 0:1,
                        mean = c(TRUE, FALSE),
                        distr = c("std","norm","ged","nig"))
  write.csv(comb_full, filename_param_dtsPC, row.names = FALSE)
  results <- ARIMA_GARCH_test(param = comb_full, ts = dtsPC,
                              filename = filename_results_dtsPC,
                              type.model = type.model)
  param <- comb_full
} else param <- read.csv(filename_param_dtsPC)

# filter only best models
alpha <- 1 - 0.995

valid_models <- results %>% mutate(rankAIC = rank(-akaike),
                                   rankBIC = rank(-bic),
                                   rankSHIB = rank(-shibata),
                                   rankHANQ = rank(-hannan_quinn),
                                   score = rankAIC + rankBIC + rankSHIB +
                                     rankHANQ) %>%
  filter(max_p_value < alpha, max_estim < 1) %>% arrange(desc(score)) %>%
  select(-rankAIC, -rankBIC, -rankSHIB, -rankHANQ)

# see models selected
param_valid <- param[valid_models$N,]
valid_models

valid_models$distr <- param_valid$distr

# see one particular models of those selected
model_eGarch_PC2 <- fit_GARCH(dtsPC, param[valid_models$N[3],], type.model = type.model)

#### 3.2.2.3 ARIMA gjrGARCH model fitting ---- 

type.model <- "gjrGARCH"

filename_results_dtsPC <- paste0("results_",type.model,"_dts",
                                 colnames(tsPC),".csv")
filename_param_dtsPC <- paste0("param_",type.model,"_dts",
                               colnames(tsPC),".csv")

# tests multiple ARMA order GARCH order models
results <- tryCatch(read.csv(filename_results_dtsPC), error = function(e) 1,
                    warning = function(w) 1)
# check if there is a file with the results and otherwise, calculate them
if (typeof(results) == "double"){
  comb_full <- crossing(ar1 = 0:1, ar2 = 0:1, ar3 = 0:1,
                        ma1 = 0:1, ma2 = 0:1, ma3 = 0:1,
                        alpha1 = 0:1, alpha2 = 0:1, alpha3 = 0:1,
                        beta1 = 0:1, beta2 = 0:1, beta3 = 0:1,
                        mean = c(TRUE, FALSE),
                        distr = c("std","norm","ged","nig"))
  write.csv(comb_full, filename_param_dtsPC, row.names = FALSE)
  results <- ARIMA_GARCH_test(param = comb_full, ts = dtsPC,
                              filename = filename_results_dtsPC,
                              type.model = type.model)
  param <- comb_full
} else param <- read.csv(filename_param_dtsPC)

# filter only best models
alpha <- 1 - 0.995

valid_models <- results %>% mutate(rankAIC = rank(-akaike),
                                   rankBIC = rank(-bic),
                                   rankSHIB = rank(-shibata),
                                   rankHANQ = rank(-hannan_quinn),
                                   score = rankAIC + rankBIC + rankSHIB +
                                     rankHANQ) %>%
  filter(max_p_value < alpha, max_estim < 1) %>% arrange(desc(score)) %>%
  select(-rankAIC, -rankBIC, -rankSHIB, -rankHANQ)

# see models selected
param_valid <- param[valid_models$N,]
valid_models

valid_models$distr <- param_valid$distr

# see one particular models of those selected
model_gjrGarch_PC2 <- fit_GARCH(dtsPC, param[valid_models$N[3],], type.model = type.model)

### 3.2.3 VaR backtesting ####

model_dtsPC2 <- model_sGarch_PC2
model <- model_sGarch_PC2

spec <- fit_GARCH_spec(param[valid_models$N[9],], type.model = type.model)

modelroll <- ugarchroll(
  spec=spec, data=dtsPC, n.ahead = 1, forecast.length = 1119, refit.every = 100,
  refit.window = c("recursive"), solver = "hybrid", calculate.VaR = TRUE,
  VaR.alpha = c(0.01,0.05,0.95,0.99))

VaR001 <- modelroll@forecast$VaR[,"alpha(1%)"]
VaR005 <- modelroll@forecast$VaR[,"alpha(5%)"]
VaR095 <- modelroll@forecast$VaR[,"alpha(95%)"]
VaR099 <- modelroll@forecast$VaR[,"alpha(99%)"]
real <- modelroll@forecast$VaR[,"realized"]

mean(real<VaR001)
mean(real<VaR005)
mean(real<VaR095)
mean(real<VaR099)

VaR <- list()
VaR[[1]] <- VaRTest(alpha = 0.01, real, VaR001, conf.level = 0.995)
VaR[[2]] <- VaRTest(alpha = 0.05, real, VaR005, conf.level = 0.995)
VaR[[3]] <- VaRTest(alpha = 0.05, -real, -VaR095, conf.level = 0.995)
VaR[[4]] <- VaRTest(alpha = 0.01, -real, -VaR099, conf.level = 0.995)

VaRDur <- list()
VaRDur[[1]] <- VaRDurTest(alpha = 0.01, real, VaR001, conf.level = 0.995)
VaRDur[[2]] <- VaRDurTest(alpha = 0.05, real, VaR005, conf.level = 0.995)
VaRDur[[3]] <- VaRDurTest(alpha = 0.01, -real, -VaR095, conf.level = 0.995)
VaRDur[[4]] <- VaRDurTest(alpha = 0.05, -real, -VaR099, conf.level = 0.995)

DAC <- list()
DAC[[1]] <- DACTest(modelroll@forecast$density$Mu, real, test = "PT", conf.level = 0.95)
DAC[[2]] <- DACTest(modelroll@forecast$density$Mu, real, test = "AG", conf.level = 0.95)

results <- matrix(NA,ncol = 3, nrow = 4)
colnames(results) <- c("unconditional", "conditional", "duration")
rownames(results) <- c(0.01,0.05,0.95,0.99)
for (i in 1:4){
  results[i,] <- c(VaR[[i]]$uc.LRp, VaR[[i]]$cc.LRp, VaRDur[[i]]$LRp) 
}

write.xlsx(as.data.frame(results), file = "tables/backtestPC2.xlsx")

VaR_down <- quantile(model, probs = c(0.01, 0.05))
VaR[[1]] <- VaRTest(alpha = 0.01, dtsPC, VaR_down[,1], conf.level = 0.95)
VaR[[2]] <- VaRTest(alpha = 0.05, dtsPC, VaR_down[,2], conf.level = 0.95)
VaRDur[[1]] <- VaRDurTest(alpha = 0.01, dtsPC, VaR_down[,1], conf.level = 0.95)
VaRDur[[2]] <- VaRDurTest(alpha = 0.05, dtsPC, VaR_down[,2], conf.level = 0.95)

VaR_up <- -quantile(model, probs = c(0.95, 0.99))
VaR[[3]] <- VaRTest(alpha = 0.05, -dtsPC, VaR_up[,1], conf.level = 0.95)
VaR[[4]] <- VaRTest(alpha = 0.01, -dtsPC, VaR_up[,2], conf.level = 0.95)
VaRDur[[3]] <- VaRDurTest(alpha = 0.01, -dtsPC, VaR_up[,1], conf.level = 0.95)
VaRDur[[4]] <- VaRDurTest(alpha = 0.05, -dtsPC, VaR_up[,2], conf.level = 0.95)

DAC[[1]] <- DACTest(fitted(model), dtsPC, test = "PT", conf.level = 0.95)
DAC[[2]] <- DACTest(fitted(model), dtsPC, test = "AG", conf.level = 0.95)

for (i in 1:4){
  results[i,] <- c(VaR[[i]]$uc.LRp, VaR[[i]]$cc.LRp, VaRDur[[i]]$LRp) 
}


write.xlsx(as.data.frame(results), file = "tables/backtestPC2.xlsx")

## 3.3 PC3 ----

### 3.3.1 Stationarity check ----

tsPC <- tsPCX[,3]

svg(file = 'plots/PC3_hist.svg', width = 7, height = 4)
par(mfrow=c(1,1), mar = c(5,5,1,1),oma = c(0,0,0,0))
plot(tsPC, type = "l", main = "PC1 (level)", col = "blue")
dev.off()


ggdraw() +
  draw_plot(ggAcf(tsPC), x = 0, y = .5, width = .5, height = .5) +
  draw_plot(ggPacf(tsPC), x = .5, y = .5, width = .5, height = .5) +
  draw_plot(autoplot(tsPC), x = 0, y = 0, width = 1, height = 0.5)

adf.test(tsPC)
pp.test(tsPC)
kpss.test(tsPC)

dtsPC <- diff_narm(tsPC)

ggdraw() +
  draw_plot(ggAcf(dtsPC), x = 0, y = .5, width = .5, height = .5) +
  draw_plot(ggPacf(dtsPC), x = .5, y = .5, width = .5, height = .5) +
  draw_plot(autoplot(dtsPC), x = 0, y = 0, width = 1, height = 0.5)

station <- list()
station[[1]] <- adf.test(dtsPC)
station[[2]] <- pp.test(dtsPC)
station[[3]] <- kpss.test(dtsPC)

statistic <- c()
pvalue <- c()
alternative <- c()

for (i in 1:3) statistic <- c(statistic, station[[i]]$statistic)
for (i in 1:3) pvalue <- c(pvalue, station[[i]]$p.value)
for (i in 1:3) alternative <- c(alternative, station[[i]]$alternative)
alternative[3] <- "non stationary"

station <- data.frame("statistic" = statistic, "pvalue" = pvalue, "alternative hypothesis"=
                        alternative, row.names = c("Dickey-Fuller", "Phillips-Perron", "KPSS"))

write.xlsx(station, file = "tables/stationarityPC3.xlsx")

### 3.3.2 Estimating the model ARMA GARCH ----

#### 3.3.2.1 ARIMA sGARCH model fitting ----

type.model <- "sGARCH"

filename_results_dtsPC <- paste0("results_",type.model,"_dts",
                                 colnames(tsPC),".csv")
filename_param_dtsPC <- paste0("param_",type.model,"_dts",
                               colnames(tsPC),".csv")

# tests multiple ARMA order GARCH order models
results <- tryCatch(read.csv(filename_results_dtsPC), error = function(e) 1,
                    warning = function(w) 1)
# check if there is a file with the results and otherwise, calculate them
if (typeof(results) == "double"){
  comb_full <- crossing(ar1 = 0:1, ar2 = 0:1, ar3 = 0:1,
                        ma1 = 0:1, ma2 = 0:1, ma3 = 0:1,
                        alpha1 = 0:1, alpha2 = 0:1, alpha3 = 0:1,
                        beta1 = 0:1, beta2 = 0:1, beta3 = 0:1,
                        mean = c(TRUE, FALSE),
                        distr = c("std","norm","ged","nig"))
  write.csv(comb_full, filename_param_dtsPC, row.names = FALSE)
  results <- ARIMA_GARCH_test(param = comb_full, ts = dtsPC,
                              filename = filename_results_dtsPC,
                              type.model = type.model)
  param <- comb_full
} else param <- read.csv(filename_param_dtsPC)

# filter only best models
alpha <- 1 - 0.995

valid_models <- results %>% mutate(rankAIC = rank(-akaike),
                                   rankBIC = rank(-bic),
                                   rankSHIB = rank(-shibata),
                                   rankHANQ = rank(-hannan_quinn),
                                   score = rankAIC + rankBIC + rankSHIB +
                                     rankHANQ) %>%
  filter(max_p_value < alpha, max_estim < 1) %>% arrange(desc(score)) %>%
  select(-rankAIC, -rankBIC, -rankSHIB, -rankHANQ)

# see models selected
param_valid <- param[valid_models$N,]
valid_models

valid_models$distr <- param_valid$distr

# see one particular models of those selected
(model_sGarch_PC3 <- fit_GARCH(dtsPC, param[valid_models$N[26],], type.model = type.model))


#### 3.3.2.2 ARIMA eGARCH model fitting ----

type.model <- "eGARCH"

filename_results_dtsPC <- paste0("results_",type.model,"_dts",
                                 colnames(tsPC),".csv")
filename_param_dtsPC <- paste0("param_",type.model,"_dts",
                               colnames(tsPC),".csv")

# tests multiple ARMA order GARCH order models
results <- tryCatch(read.csv(filename_results_dtsPC), error = function(e) 1,
                    warning = function(w) 1)
# check if there is a file with the results and otherwise, calculate them
if (typeof(results) == "double"){
  comb_full <- crossing(ar1 = 0:1, ar2 = 0:1, ar3 = 0:1,
                        ma1 = 0:1, ma2 = 0:1, ma3 = 0:1,
                        alpha1 = 0:1, alpha2 = 0:1, alpha3 = 0:1,
                        beta1 = 0:1, beta2 = 0:1, beta3 = 0:1,
                        mean = c(TRUE, FALSE),
                        distr = c("std","norm","ged","nig"))
  write.csv(comb_full, filename_param_dtsPC, row.names = FALSE)
  results <- ARIMA_GARCH_test(param = comb_full, ts = dtsPC,
                              filename = filename_results_dtsPC,
                              type.model = type.model)
  param <- comb_full
} else param <- read.csv(filename_param_dtsPC)

# filter only best models
alpha <- 1 - 0.995

valid_models <- results %>% mutate(rankAIC = rank(-akaike),
                                   rankBIC = rank(-bic),
                                   rankSHIB = rank(-shibata),
                                   rankHANQ = rank(-hannan_quinn),
                                   score = rankAIC + rankBIC + rankSHIB +
                                     rankHANQ) %>%
  filter(max_p_value < alpha, max_estim < 1) %>% arrange(desc(score)) %>%
  select(-rankAIC, -rankBIC, -rankSHIB, -rankHANQ)

# see models selected
param_valid <- param[valid_models$N,]
valid_models

valid_models$distr <- param_valid$distr

# see one particular models of those selected
model_eGarch_PC3 <- fit_GARCH(dtsPC, param[valid_models$N[3],], type.model = type.model)

#### 3.3.2.3 ARIMA gjrGARCH model fitting ---- 

type.model <- "gjrGARCH"

filename_results_dtsPC <- paste0("results_",type.model,"_dts",
                                 colnames(tsPC),".csv")
filename_param_dtsPC <- paste0("param_",type.model,"_dts",
                               colnames(tsPC),".csv")

# tests multiple ARMA order GARCH order models
results <- tryCatch(read.csv(filename_results_dtsPC), error = function(e) 1,
                    warning = function(w) 1)
# check if there is a file with the results and otherwise, calculate them
if (typeof(results) == "double"){
  comb_full <- crossing(ar1 = 0:1, ar2 = 0:1, ar3 = 0:1,
                        ma1 = 0:1, ma2 = 0:1, ma3 = 0:1,
                        alpha1 = 0:1, alpha2 = 0:1, alpha3 = 0:1,
                        beta1 = 0:1, beta2 = 0:1, beta3 = 0:1,
                        mean = c(TRUE, FALSE),
                        distr = c("std","norm","ged","nig"))
  write.csv(comb_full, filename_param_dtsPC, row.names = FALSE)
  results <- ARIMA_GARCH_test(param = comb_full, ts = dtsPC,
                              filename = filename_results_dtsPC,
                              type.model = type.model)
  param <- comb_full
} else param <- read.csv(filename_param_dtsPC)

# filter only best models
alpha <- 1 - 0.995

valid_models <- results %>% mutate(rankAIC = rank(-akaike),
                                   rankBIC = rank(-bic),
                                   rankSHIB = rank(-shibata),
                                   rankHANQ = rank(-hannan_quinn),
                                   score = rankAIC + rankBIC + rankSHIB +
                                     rankHANQ) %>%
  filter(max_p_value < alpha, max_estim < 1) %>% arrange(desc(score)) %>%
  select(-rankAIC, -rankBIC, -rankSHIB, -rankHANQ)

# see models selected
param_valid <- param[valid_models$N,]
valid_models

valid_models$distr <- param_valid$distr

# see one particular models of those selected
model_gjrGarch_PC3 <- fit_GARCH(dtsPC, param[valid_models$N[3],], type.model = type.model)

### 3.3.3 VaR backtesting ####

model_dtsPC3 <- model_sGarch_PC3
model <- model_sGarch_PC3

spec <- fit_GARCH_spec(param[valid_models$N[9],], type.model = type.model)

modelroll=ugarchroll(
  spec=spec, data=dtsPC, n.ahead = 1, forecast.length = 1119, refit.every = 100,
  refit.window = c("recursive"), solver = "hybrid", calculate.VaR = TRUE,
  VaR.alpha = c(0.01,0.05,0.95,0.99))

VaR001 <- modelroll@forecast$VaR[,"alpha(1%)"]
VaR005 <- modelroll@forecast$VaR[,"alpha(5%)"]
VaR095 <- modelroll@forecast$VaR[,"alpha(95%)"]
VaR099 <- modelroll@forecast$VaR[,"alpha(99%)"]
real <- modelroll@forecast$VaR[,"realized"]

mean(real<VaR001)
mean(real<VaR005)
mean(real<VaR095)
mean(real<VaR099)

VaR <- list()
VaR[[1]] <- VaRTest(alpha = 0.01, real, VaR001, conf.level = 0.995)
VaR[[2]] <- VaRTest(alpha = 0.05, real, VaR005, conf.level = 0.995)
VaR[[3]] <- VaRTest(alpha = 0.05, -real, -VaR095, conf.level = 0.995)
VaR[[4]] <- VaRTest(alpha = 0.01, -real, -VaR099, conf.level = 0.995)

VaRDur <- list()
VaRDur[[1]] <- VaRDurTest(alpha = 0.01, real, VaR001, conf.level = 0.995)
VaRDur[[2]] <- VaRDurTest(alpha = 0.05, real, VaR005, conf.level = 0.995)
VaRDur[[3]] <- VaRDurTest(alpha = 0.01, -real, -VaR095, conf.level = 0.995)
VaRDur[[4]] <- VaRDurTest(alpha = 0.05, -real, -VaR099, conf.level = 0.995)

DAC <- list()
DAC[[1]] <- DACTest(modelroll@forecast$density$Mu, real, test = "PT", conf.level = 0.95)
DAC[[2]] <- DACTest(modelroll@forecast$density$Mu, real, test = "AG", conf.level = 0.95)

results <- matrix(NA,ncol = 3, nrow = 4)
colnames(results) <- c("unconditional", "conditional", "duration")
rownames(results) <- c(0.01,0.05,0.95,0.99)
for (i in 1:4){
  results[i,] <- c(VaR[[i]]$uc.LRp, VaR[[i]]$cc.LRp, VaRDur[[i]]$LRp) 
}

write.xlsx(as.data.frame(results), file = "tables/backtestPC3.xlsx")

VaR_down <- quantile(model, probs = c(0.01, 0.05))
VaR[[1]] <- VaRTest(alpha = 0.01, dtsPC, VaR_down[,1], conf.level = 0.95)
VaR[[2]] <- VaRTest(alpha = 0.05, dtsPC, VaR_down[,2], conf.level = 0.95)
VaRDur[[1]] <- VaRDurTest(alpha = 0.01, dtsPC, VaR_down[,1], conf.level = 0.95)
VaRDur[[2]] <- VaRDurTest(alpha = 0.05, dtsPC, VaR_down[,2], conf.level = 0.95)

VaR_up <- -quantile(model, probs = c(0.95, 0.99))
VaR[[3]] <- VaRTest(alpha = 0.05, -dtsPC, VaR_up[,1], conf.level = 0.95)
VaR[[4]] <- VaRTest(alpha = 0.01, -dtsPC, VaR_up[,2], conf.level = 0.95)
VaRDur[[3]] <- VaRDurTest(alpha = 0.01, -dtsPC, VaR_up[,1], conf.level = 0.95)
VaRDur[[4]] <- VaRDurTest(alpha = 0.05, -dtsPC, VaR_up[,2], conf.level = 0.95)

DAC[[1]] <- DACTest(fitted(model), dtsPC, test = "PT", conf.level = 0.95)
DAC[[2]] <- DACTest(fitted(model), dtsPC, test = "AG", conf.level = 0.95)

for (i in 1:4){
  results[i,] <- c(VaR[[i]]$uc.LRp, VaR[[i]]$cc.LRp, VaRDur[[i]]$LRp) 
}


write.xlsx(as.data.frame(results), file = "tables/backtestPC3.xlsx")


## 3.X STOXX600 ----

DJSTOXX <- xts(market_df$DJSTOXX, order.by = market_df$date)

### 3.X.1 Stationarity check ----

tsPC <- DJSTOXX
colnames(tsPC) <- "DJSTOXX"

svg(file = 'plots/ts_DJSTOXX.svg', width = 16.5/1.6, height = 5/1.6)
par(mfrow=c(1,1), mar = c(5,5,1,1),oma = c(0,0,0,0))
plot(dates, tsPC, type = "l", col = "#6A7434", xlab = "time", ylab = "")
legend("bottomright",legend = paste("EUROSTOXX 600"), col = c("#6A7434"),lwd = 1)
dev.off()

ggdraw() +
  draw_plot(ggAcf(tsPC), x = 0, y = .5, width = .5, height = .5) +
  draw_plot(ggPacf(tsPC), x = .5, y = .5, width = .5, height = .5) +
  draw_plot(autoplot(tsPC), x = 0, y = 0, width = 1, height = 0.5)

adf.test(tsPC)
pp.test(tsPC)
kpss.test(tsPC)

dtsPC <- diff_narm(tsPC)/tsPC

ggdraw() +
  draw_plot(ggAcf(dtsPC), x = 0, y = .5, width = .5, height = .5) +
  draw_plot(ggPacf(dtsPC), x = .5, y = .5, width = .5, height = .5) +
  draw_plot(autoplot(dtsPC), x = 0, y = 0, width = 1, height = 0.5)

station <- list()
station[[1]] <- adf.test(dtsPC)
station[[2]] <- pp.test(dtsPC)
station[[3]] <- kpss.test(dtsPC)

statistic <- c()
pvalue <- c()
alternative <- c()

for (i in 1:3) statistic <- c(statistic, station[[i]]$statistic)
for (i in 1:3) pvalue <- c(pvalue, station[[i]]$p.value)
for (i in 1:3) alternative <- c(alternative, station[[i]]$alternative)
alternative[3] <- "non stationary"

station <- data.frame("statistic" = statistic, "pvalue" = pvalue, "alternative hypothesis"=
                        alternative, row.names = c("Dickey-Fuller", "Phillips-Perron", "KPSS"))

write.xlsx(station, file = "tables/stationarityDJSTOXX.xlsx")

#### 3.X.2.1 ARIMA sGARCH model fitting ----

type.model <- "sGARCH"

filename_results_dtsPC <- paste0("results_",type.model,"_dts",
                                 colnames(tsPC),".csv")
filename_param_dtsPC <- paste0("param_",type.model,"_dts",
                               colnames(tsPC),".csv")

# tests multiple ARMA order GARCH order models
results <- tryCatch(read.csv(filename_results_dtsPC), error = function(e) 1,
                    warning = function(w) 1)
# check if there is a file with the results and otherwise, calculate them
if (typeof(results) == "double"){
  comb_full <- crossing(ar1 = 0:1, ar2 = 0:1, ar3 = 0:1,
                        ma1 = 0:1, ma2 = 0:1, ma3 = 0:1,
                        alpha1 = 0:1, alpha2 = 0:1, alpha3 = 0:1,
                        beta1 = 0:1, beta2 = 0:1, beta3 = 0:1,
                        mean = c(TRUE, FALSE),
                        distr = c("std","norm","ged","nig"))
  write.csv(comb_full, filename_param_dtsPC, row.names = FALSE)
  results <- ARIMA_GARCH_test(param = comb_full, ts = dtsPC,
                              filename = filename_results_dtsPC,
                              type.model = type.model)
  param <- comb_full
} else param <- read.csv(filename_param_dtsPC)

# filter only best models
alpha <- 1 - 0.995

valid_models <- results %>% mutate(rankAIC = rank(-akaike),
                                   rankBIC = rank(-bic),
                                   rankSHIB = rank(-shibata),
                                   rankHANQ = rank(-hannan_quinn),
                                   score = rankAIC + rankBIC + rankSHIB +
                                     rankHANQ) %>%
  filter(max_p_value < alpha, max_estim < 1) %>% arrange(desc(score)) %>%
  select(-rankAIC, -rankBIC, -rankSHIB, -rankHANQ)

# see models selected
param_valid <- param[valid_models$N,]
valid_models

valid_models$distr <- param_valid$distr

# see one particular models of those selected
(model_sGarch_DJSTOXX <- fit_GARCH(dtsPC, param[valid_models$N[8],], type.model = type.model))

### 3.X.3 VaR backtesting ####

model_dtsDJSTOXX <- model_sGarch_DJSTOXX
model <- model_sGarch_DJSTOXX

spec <- fit_GARCH_spec(param[valid_models$N[9],], type.model = type.model)

modelroll=ugarchroll(
  spec=spec, data=dtsPC, n.ahead = 1, forecast.length = 1119, refit.every = 100,
  refit.window = c("recursive"), solver = "hybrid", calculate.VaR = TRUE,
  VaR.alpha = c(0.01,0.05,0.95,0.99))

VaR001 <- modelroll@forecast$VaR[,"alpha(1%)"]
VaR005 <- modelroll@forecast$VaR[,"alpha(5%)"]
VaR095 <- modelroll@forecast$VaR[,"alpha(95%)"]
VaR099 <- modelroll@forecast$VaR[,"alpha(99%)"]
real <- modelroll@forecast$VaR[,"realized"]

mean(real<VaR001)
mean(real<VaR005)
mean(real<VaR095)
mean(real<VaR099)

VaR <- list()
VaR[[1]] <- VaRTest(alpha = 0.01, real, VaR001, conf.level = 0.995)
VaR[[2]] <- VaRTest(alpha = 0.05, real, VaR005, conf.level = 0.995)
VaR[[3]] <- VaRTest(alpha = 0.05, -real, -VaR095, conf.level = 0.995)
VaR[[4]] <- VaRTest(alpha = 0.01, -real, -VaR099, conf.level = 0.995)

VaRDur <- list()
VaRDur[[1]] <- VaRDurTest(alpha = 0.01, real, VaR001, conf.level = 0.995)
VaRDur[[2]] <- VaRDurTest(alpha = 0.05, real, VaR005, conf.level = 0.995)
VaRDur[[3]] <- VaRDurTest(alpha = 0.01, -real, -VaR095, conf.level = 0.995)
VaRDur[[4]] <- VaRDurTest(alpha = 0.05, -real, -VaR099, conf.level = 0.995)

DAC <- list()
DAC[[1]] <- DACTest(modelroll@forecast$density$Mu, real, test = "PT", conf.level = 0.95)
DAC[[2]] <- DACTest(modelroll@forecast$density$Mu, real, test = "AG", conf.level = 0.95)

results <- matrix(NA,ncol = 3, nrow = 4)
colnames(results) <- c("unconditional", "conditional", "duration")
rownames(results) <- c(0.01,0.05,0.95,0.99)
for (i in 1:4){
  results[i,] <- c(VaR[[i]]$uc.LRp, VaR[[i]]$cc.LRp, VaRDur[[i]]$LRp) 
}

write.xlsx(as.data.frame(results), file = "tables/backtestDJSTOXX.xlsx")

VaR_down <- quantile(model, probs = c(0.01, 0.05))
VaR[[1]] <- VaRTest(alpha = 0.01, dtsPC, VaR_down[,1], conf.level = 0.95)
VaR[[2]] <- VaRTest(alpha = 0.05, dtsPC, VaR_down[,2], conf.level = 0.95)
VaRDur[[1]] <- VaRDurTest(alpha = 0.01, dtsPC, VaR_down[,1], conf.level = 0.95)
VaRDur[[2]] <- VaRDurTest(alpha = 0.05, dtsPC, VaR_down[,2], conf.level = 0.95)

VaR_up <- -quantile(model, probs = c(0.95, 0.99))
VaR[[3]] <- VaRTest(alpha = 0.05, -dtsPC, VaR_up[,1], conf.level = 0.95)
VaR[[4]] <- VaRTest(alpha = 0.01, -dtsPC, VaR_up[,2], conf.level = 0.95)
VaRDur[[3]] <- VaRDurTest(alpha = 0.01, -dtsPC, VaR_up[,1], conf.level = 0.95)
VaRDur[[4]] <- VaRDurTest(alpha = 0.05, -dtsPC, VaR_up[,2], conf.level = 0.95)

DAC[[1]] <- DACTest(fitted(model), dtsPC, test = "PT", conf.level = 0.95)
DAC[[2]] <- DACTest(fitted(model), dtsPC, test = "AG", conf.level = 0.95)

for (i in 1:4){
  results[i,] <- c(VaR[[i]]$uc.LRp, VaR[[i]]$cc.LRp, VaRDur[[i]]$LRp) 
}


write.xlsx(as.data.frame(results), file = "tables/backtestDJSTOXX.xlsx")

## 3.Y USDEUR ----

USDEUR <- xts(market_df$USDEUR, order.by = market_df$date)

### 3.Y.1 Stationarity check ----

tsPC <- USDEUR
colnames(tsPC) <- "USDEUR"

svg(file = 'plots/ts_USDEUR.svg', width = 16.5/1.6, height = 5/1.6)
par(mfrow=c(1,1), mar = c(5,5,1,1),oma = c(0,0,0,0))
plot(dates, tsPC, type = "l", col = "#6A7434", xlab = "time", ylab = "")
legend("bottomright",legend = paste("USDEUR"), col = c("#6A7434"),lwd = 1)
dev.off()

ggdraw() +
  draw_plot(ggAcf(tsPC), x = 0, y = .5, width = .5, height = .5) +
  draw_plot(ggPacf(tsPC), x = .5, y = .5, width = .5, height = .5) +
  draw_plot(autoplot(tsPC), x = 0, y = 0, width = 1, height = 0.5)

adf.test(tsPC)
pp.test(tsPC)
kpss.test(tsPC)

dtsPC <- diff_narm(tsPC)/tsPC

ggdraw() +
  draw_plot(ggAcf(dtsPC), x = 0, y = .5, width = .5, height = .5) +
  draw_plot(ggPacf(dtsPC), x = .5, y = .5, width = .5, height = .5) +
  draw_plot(autoplot(dtsPC), x = 0, y = 0, width = 1, height = 0.5)

station <- list()
station[[1]] <- adf.test(dtsPC)
station[[2]] <- pp.test(dtsPC)
station[[3]] <- kpss.test(dtsPC)

statistic <- c()
pvalue <- c()
alternative <- c()

for (i in 1:3) statistic <- c(statistic, station[[i]]$statistic)
for (i in 1:3) pvalue <- c(pvalue, station[[i]]$p.value)
for (i in 1:3) alternative <- c(alternative, station[[i]]$alternative)
alternative[3] <- "non stationary"

station <- data.frame("statistic" = statistic, "pvalue" = pvalue, "alternative hypothesis"=
                        alternative, row.names = c("Dickey-Fuller", "Phillips-Perron", "KPSS"))

write.xlsx(station, file = "tables/stationarityUSDEUR.xlsx")

#### 3.1.2.1 ARIMA sGARCH model fitting ----

type.model <- "sGARCH"

filename_results_dtsPC <- paste0("results_",type.model,"_dts",
                                 colnames(tsPC),".csv")
filename_param_dtsPC <- paste0("param_",type.model,"_dts",
                               colnames(tsPC),".csv")

# tests multiple ARMA order GARCH order models
results <- tryCatch(read.csv(filename_results_dtsPC), error = function(e) 1,
                    warning = function(w) 1)
# check if there is a file with the results and otherwise, calculate them
if (typeof(results) == "double"){
  comb_full <- crossing(ar1 = 0:1, ar2 = 0:1, ar3 = 0:1,
                        ma1 = 0:1, ma2 = 0:1, ma3 = 0:1,
                        alpha1 = 0:1, alpha2 = 0:1, alpha3 = 0:1,
                        beta1 = 0:1, beta2 = 0:1, beta3 = 0:1,
                        mean = c(TRUE, FALSE),
                        distr = c("std","norm","ged","nig"))
  write.csv(comb_full, filename_param_dtsPC, row.names = FALSE)
  results <- ARIMA_GARCH_test(param = comb_full, ts = dtsPC,
                              filename = filename_results_dtsPC,
                              type.model = type.model)
  param <- comb_full
} else param <- read.csv(filename_param_dtsPC)

# filter only best models
alpha <- 1 - 0.995

valid_models <- results %>% mutate(rankAIC = rank(-akaike),
                                   rankBIC = rank(-bic),
                                   rankSHIB = rank(-shibata),
                                   rankHANQ = rank(-hannan_quinn),
                                   score = rankAIC + rankBIC + rankSHIB +
                                     rankHANQ) %>%
  filter(max_p_value < alpha, max_estim < 1) %>% arrange(desc(score)) %>%
  select(-rankAIC, -rankBIC, -rankSHIB, -rankHANQ)

# see models selected
param_valid <- param[valid_models$N,]
valid_models

valid_models$distr <- param_valid$distr

# see one particular models of those selected
(model_sGarch_USDEUR <- fit_GARCH(dtsPC, param[valid_models$N[9],], type.model = type.model))

### 3.Y.3 VaR backtesting ####

model_dtsUSDEUR <- model_sGarch_USDEUR
model <- model_sGarch_USDEUR

spec <- fit_GARCH_spec(param[valid_models$N[9],], type.model = type.model)

modelroll=ugarchroll(
  spec=spec, data=dtsPC, n.ahead = 1, forecast.length = 1119, refit.every = 100,
  refit.window = c("recursive"), solver = "hybrid", calculate.VaR = TRUE,
  VaR.alpha = c(0.01,0.05,0.95,0.99))

VaR001 <- modelroll@forecast$VaR[,"alpha(1%)"]
VaR005 <- modelroll@forecast$VaR[,"alpha(5%)"]
VaR095 <- modelroll@forecast$VaR[,"alpha(95%)"]
VaR099 <- modelroll@forecast$VaR[,"alpha(99%)"]
real <- modelroll@forecast$VaR[,"realized"]

mean(real<VaR001)
mean(real<VaR005)
mean(real<VaR095)
mean(real<VaR099)

VaR <- list()
VaR[[1]] <- VaRTest(alpha = 0.01, real, VaR001, conf.level = 0.995)
VaR[[2]] <- VaRTest(alpha = 0.05, real, VaR005, conf.level = 0.995)
VaR[[3]] <- VaRTest(alpha = 0.05, -real, -VaR095, conf.level = 0.995)
VaR[[4]] <- VaRTest(alpha = 0.01, -real, -VaR099, conf.level = 0.995)

VaRDur <- list()
VaRDur[[1]] <- VaRDurTest(alpha = 0.01, real, VaR001, conf.level = 0.995)
VaRDur[[2]] <- VaRDurTest(alpha = 0.05, real, VaR005, conf.level = 0.995)
VaRDur[[3]] <- VaRDurTest(alpha = 0.01, -real, -VaR095, conf.level = 0.995)
VaRDur[[4]] <- VaRDurTest(alpha = 0.05, -real, -VaR099, conf.level = 0.995)

DAC <- list()
DAC[[1]] <- DACTest(modelroll@forecast$density$Mu, real, test = "PT", conf.level = 0.95)
DAC[[2]] <- DACTest(modelroll@forecast$density$Mu, real, test = "AG", conf.level = 0.95)

results <- matrix(NA,ncol = 3, nrow = 4)
colnames(results) <- c("unconditional", "conditional", "duration")
rownames(results) <- c(0.01,0.05,0.95,0.99)
for (i in 1:4){
  results[i,] <- c(VaR[[i]]$uc.LRp, VaR[[i]]$cc.LRp, VaRDur[[i]]$LRp) 
}

write.xlsx(as.data.frame(results), file = "tables/backtestlUSDEUR.xlsx")

VaR_down <- quantile(model, probs = c(0.01, 0.05))
VaR[[1]] <- VaRTest(alpha = 0.01, dtsPC, VaR_down[,1], conf.level = 0.95)
VaR[[2]] <- VaRTest(alpha = 0.05, dtsPC, VaR_down[,2], conf.level = 0.95)
VaRDur[[1]] <- VaRDurTest(alpha = 0.01, dtsPC, VaR_down[,1], conf.level = 0.95)
VaRDur[[2]] <- VaRDurTest(alpha = 0.05, dtsPC, VaR_down[,2], conf.level = 0.95)

VaR_up <- -quantile(model, probs = c(0.95, 0.99))
VaR[[3]] <- VaRTest(alpha = 0.05, -dtsPC, VaR_up[,1], conf.level = 0.95)
VaR[[4]] <- VaRTest(alpha = 0.01, -dtsPC, VaR_up[,2], conf.level = 0.95)
VaRDur[[3]] <- VaRDurTest(alpha = 0.01, -dtsPC, VaR_up[,1], conf.level = 0.95)
VaRDur[[4]] <- VaRDurTest(alpha = 0.05, -dtsPC, VaR_up[,2], conf.level = 0.95)

DAC[[1]] <- DACTest(fitted(model), dtsPC, test = "PT", conf.level = 0.95)
DAC[[2]] <- DACTest(fitted(model), dtsPC, test = "AG", conf.level = 0.95)

for (i in 1:4){
  results[i,] <- c(VaR[[i]]$uc.LRp, VaR[[i]]$cc.LRp, VaRDur[[i]]$LRp) 
}


write.xlsx(as.data.frame(results), file = "tables/backtestUSDEUR.xlsx")

## 3.0 Summary ----

write.xlsx(as.data.frame(model_dtsPC1@fit$matcoef), file = "tables/modelPC1.xlsx")
write.xlsx(as.data.frame(model_dtsPC2@fit$matcoef), file = "tables/modelPC2.xlsx")
write.xlsx(as.data.frame(model_dtsPC3@fit$matcoef), file = "tables/modelPC3.xlsx")
write.xlsx(as.data.frame(model_dtsDJSTOXX@fit$matcoef), file = "tables/modelDJSTOXX.xlsx")
write.xlsx(as.data.frame(model_dtsUSDEUR@fit$matcoef), file = "tables/modelUSDEUR.xlsx")

write.xlsx(as.data.frame(infocriteria(model_dtsPC1)), file = "tables/infocriteriaPC1.xlsx")
write.xlsx(as.data.frame(infocriteria(model_dtsPC2)), file = "tables/infocriteriaPC2.xlsx")
write.xlsx(as.data.frame(infocriteria(model_dtsPC3)), file = "tables/infocriteriaPC3.xlsx")
write.xlsx(as.data.frame(infocriteria(model_dtsDJSTOXX)), file = "tables/modelDJSTOXX.xlsx")
write.xlsx(as.data.frame(infocriteria(model_dtsUSDEUR)), file = "tables/modelUSDEUR.xlsx")


# 4. Copula estimation and diagnosis ----

## 4.1 Get residuals and analysis ----

model_dtsPC1
model_dtsPC2
model_dtsPC3
model_dtsDJSTOXX
model_dtsUSDEUR

models <- list(model_dtsPC1, model_dtsPC2, model_dtsPC3, model_dtsDJSTOXX, model_dtsUSDEUR)
models_name <- c("dtsPC1", "dtsPC2", "dtsPC3", "dtsDJSTOXX", "dtsUSDEUR")

for (i in 1:length(models)){
  model <- models[[i]]
  
  # get standarized residuals
  s_res <- residuals(model)/sigma(model)
  svg(file = paste0("plots/hist_",models_name[i],".svg"), width = 4, height = 4)
  par(mfrow=c(1,1), mar = c(2,4,1,1),oma = c(0,0,0,0))
  hist(s_res, breaks = 50, main = paste0(gsub("dts","",models_name[i]),
                                         " model residuals"), xlab = "")
  dev.off()
  sd(s_res)
  
  svg(file = paste0("plots/autocorr_",models_name[i],".svg"), width = 4, height = 2.5)
  ggdraw() +
    draw_plot(ggAcf(s_res, main = paste0(gsub("dts","",models_name[i])), xlab = ""),
              x = 0, y = 0, width = .5, height = 1) +
    draw_plot(ggPacf(s_res, main = paste0(gsub("dts","",models_name[i])), xlab = ""),
              x = 0.5, y = 0, width = .5, height = 1)
  dev.off()

  coef <- model@fit$coef
  
  lbtest <- matrix(NA, nrow = 20, ncol = 2)
  for (j in 1:20) {
    lb <- Box.test(s_res,lag = j, "Ljung-Box")
    lbtest[j,1] <- lb$statistic
    lbtest[j,2] <- lb$p.value
  }
  write.xlsx(as.data.frame(lbtest), file = paste0("tables/lbtest_",
                                                  models_name[i],".xlsx"))
  
  
  if (is.na(coef["skew"])) skew <- 1 else skew <- coef["skew"]
  if (is.na(coef["shape"])) shape <- 5 else shape <- coef["shape"]
  
  dist_model <- model@model$modeldesc$distribution
  
  # get cummulative distribution values
  u_res <- pdist(distribution = dist_model, s_res, skew = skew, shape = shape)
  u_res <- abs(u_res)
  hist(u_res, breaks = 40)
  
  # plot it for reference
  svg(file = paste0("plots/qqplot_",models_name[i],".svg"), width = 4, height = 4)
  par(mfrow=c(1,1), mar = c(2,2,2,1),oma = c(0,0,0,0))
  qqplot(rdist(distribution = dist_model,100000,
               skew = skew, shape = shape),coredata(s_res),
         main = paste0(gsub("dts","",models_name[i]), " model residuals"),
         xlab = "")
  qqline(coredata(s_res))
  dev.off()
  
  svg(file = paste0("plots/cumplot_",models_name[i],".svg"), width = 4, height = 4)
  par(mfrow=c(1,1), mar = c(2,4,2,1),oma = c(0,0,0,0))
  plot(ecdf(coredata(s_res)), col = "red", xlim = c(-4,4),
       main = paste0(gsub("dts","",models_name[i]), " model residuals"), xlab = "")
  lines(seq(-8,8,0.025),pdist(seq(-8,8,0.025),distribution = dist_model,
              skew = skew, shape = shape), col = "black")
  legend("bottomright", legend = c("Fitted", "Empirical"), lwd = c(2,2),
         col = c("black", "red"), ncol = 1)
  dev.off()
  
  
  # goodness of fit tests
  jb <- jarque.bera.test(qnorm(u_res))
  ks <- ks.test(rdist(distribution = dist_model,100000,
                skew = skew, shape = shape),coredata(s_res))
  cvm <- cvm.test(qnorm(u_res),"norm")
  t <- t.test(s_res)
  
  goftests <- matrix(NA, ncol = 2, nrow = 4)
  goftests[1,] <- c(jb$statistic, jb$p.value)
  goftests[2,] <- c(ks$statistic, ks$p.value)
  goftests[3,] <- c(cvm$statistic, cvm$p.value)
  goftests[4,] <- c(t$statistic[[1]], t$p.value[[1]])
  rownames(goftests) <- c("jarque-bera", "ks", "cvm", "t")
  
  write.xlsx(as.data.frame(goftests), file = paste0("tables/goftest_",
                                                  models_name[i],".xlsx"))
    
  if(i==1){
    s_res_all <- s_res
    # set u_res_all so it is a xts class object
    u_res_all <- s_res
    u_res_all[,1] <- u_res
  } else {
    s_res_all <- merge(s_res_all, s_res)
    u_res_all <- merge(u_res_all, u_res)
  }
}

colnames(s_res_all) <- models_name
colnames(u_res_all) <- models_name

## 4.2 Two by two copula estimation ----

pairs(coredata(s_res_all))

series_names <- paste0(gsub("dts","",models_name))

svg(file = paste0("plots/copula_all.svg"), width = 12, height = 12)
pairs(coredata(u_res_all),oma = c(2,2,0,0))
dev.off()

colnames(u_res_all) <- series_names 
R <- cor(u_res_all)
svg(file = paste0("plots/corrgram.svg"), width = 4, height = 4)
corrgram(R,oma = c(1,1,1,1))
dev.off()

# independence tests and see what copula
d <- indepTestSim(nrow(u_res_all), p = 2)
combos <- crossing(u1 = 1:5, u2 = 1:5) %>% filter(u1 != u2)
indep_pvalue <- c()
loglike <- c()
pvalue <- c()
copula <- c()
for (i in 1:nrow(combos)){
  u1 <- combos[i,1]
  u2 <- combos[i,2]
  indep_pvalue <- c(indep_pvalue, indepTest(coredata(u_res_all[,c(u1,u2)]), d)$pvalue)
  results <- testCopula(coredata(u_res_all[,c(u1,u2)]))
  results <- na.omit(results)
  loglike <- c(loglike, max(results$log.likelihood, na.rm = TRUE))
  pvalue <- c(pvalue, max(results$gof.pvalue))
  copula <- c(copula, results$copula[max(results$log.likelihood, na.rm = TRUE)==results$log.likelihood])
}

results <- cbind(combos, indep_pvalue, loglike,  pvalue, copula)
results$u1 <- series_names[results$u1]
results$u2 <- series_names[results$u2]

write.xlsx(results, file = "tables/2by2copulas.xlsx")

## 4.3 Unique copula estimation ----

u <- coredata(u_res_all[,1:4])
N <- 1000
opt.meth = "L-BFGS-B"

nc <- fitCopula(copula = normalCopula(dim = 4, dispstr = "un"), data = u,
                optim.method = opt.meth, method = "irho")

nc@loglik <- sum(log(dCopula(u, nc@copula)))

(gofnc <- gofCopula(normalCopula(dim = 4, dispstr = "un"), method = "Sn",
                    estim.method = "mpl", N = N,x = u, simulation = "mult"))

tc <- fitCopula(copula = tCopula(dim = 4, dispstr = "un"), data = u,
                optim.method = opt.meth, method = "itau.mpl")

(goftc <- gofCopula(tCopula(dim = 4, dispstr = "un", df.fixed = TRUE,
                            df = round(tc@estimate[7])),
                    estim.method = "mpl", N = N,x = u, simulation = "mult"))

# if error change optim.method as opt.meth = "Nelder-Mead"
tc <- fitCopula(copula = tCopula(dim = 4, dispstr = "un",
                                 df = round(tail(tc@estimate,1)), df.fixed = TRUE),
                data = u,optim.method = opt.meth, method = "ml")


u_t2 <- rCopula(n = 1219, copula = tc@copula)
u_t <- pCopula(u = u, copula = tc@copula)
pairs(u_t2)
pairs(u)

u_t2 <- pCopula(u = u_t2, copula = tc@copula)
hist(u_t)
hist(u_t2)

fc <- fitCopula(copula = frankCopula(dim = 4), data = u,
                optim.method = opt.meth, method = "ml")

(goffc <- gofCopula(copula = frankCopula(dim = 4), data = u,
                    estim.method = "mpl", N = N,x = u, simulation = "mult"))

jc <- fitCopula(copula = joeCopula(dim = 4), data = u,
                optim.method = opt.meth, method = "ml")

(gofjc <- gofCopula(copula = joeCopula(dim = 4), data = u,
                    estim.method = "mpl", N = N,x = u, simulation = "mult"))

cc <- fitCopula(copula = claytonCopula(dim = 4), data = u,
                optim.method = opt.meth, method = "ml")

(gofcc <- gofCopula(copula = claytonCopula(dim = 4), data = u,
                    estim.method = "mpl", N = N,x = u, simulation = "mult"))

results <- data.frame("copula" = c("Normal","t-Student",
                                  "Clayton", "Frank","Joe"),
                      "log-likelihood" = c(nc@loglik,tc@loglik,cc@loglik,
                                           fc@loglik,jc@loglik),
                      "gof pvalue" = c(gofnc$p.value,goftc$p.value,
                                       gofcc$p.value,goffc$p.value,gofjc$p.value),
                      "parameter" = c(NA, round(tail(tc@copula@parameters,1),0),
                                      cc@estimate, fc@estimate, jc@estimate))

write.xlsx(results, file = "tables/fit_copulas.xlsx")

write.xlsx(as.data.frame(nc@estimate), file = "tables/normal_param.xlsx")
write.xlsx(as.data.frame(tc@estimate), file = "tables/tStud_param.xlsx")



combos <- crossing(u1 = 1:4, u2 = 1:4) %>% filter(u1 != u2) %>% as.matrix()

## 4.4. Copula fit plots ----
# it wont create the graphs properly in the loop so it is created manually

i <- 1
u_graph <- u_res_all[,combos[i,]]
u_rand <- rCopula(n = 1000000, copula = nc@copula)
u <- seq(0, 1, length.out = 50)
grid <- as.matrix(expand.grid(u1 = u, u2 = u))
val <- cbind(grid, z = C.n(grid, X = u_rand[,combos[i,]]))
cpG <- contourplot2(val, region = FALSE,
                    key = list(corner = c(0.01, 0.01),
                               lines = list(col = c(1,4), lwd = 2),
                               text = list(c("Copula",
                                             "Data")), oma =c(0,0,0,0),mar = c(2,2,2,1)),
                    xlab= paste0("u",combos[i,1]),
                    ylab = paste0("u",combos[i,2]))
u <- seq(0, 1, length.out = 50)
grid <- as.matrix(expand.grid(u1 = u, u2 = u))
val <- cbind(grid, z = C.n(grid, X = u_graph))
cpCn <- contourplot2(val, region = FALSE, labels = FALSE, col = 4, oma = c(0,0,0,0),mar = c(2,2,2,1))
file <- paste0("plots/fit_normal_",combos[i,1],"_",combos[i,2],".svg")
svg(file = file, width = 5, height = 5)
cpG + cpCn
dev.off()
i <- i + 1 

i <- 1
u_graph <- u_res_all[,combos[i,]]
u_rand <- rCopula(n = 1000000, copula = tc@copula)
u <- seq(0, 1, length.out = 50)
grid <- as.matrix(expand.grid(u1 = u, u2 = u))
val <- cbind(grid, z = C.n(grid, X = u_rand[,combos[i,]]))
cpG <- contourplot2(val, region = FALSE,
                    key = list(corner = c(0.01, 0.01),
                               lines = list(col = c(1,3), lwd = 2),
                               text = list(c("Copula",
                                             "Data")), oma =c(0,0,0,0),mar = c(2,2,2,1)),
                    xlab= paste0("u",combos[i,1]),
                    ylab = paste0("u",combos[i,2]))
u <- seq(0, 1, length.out = 50)
grid <- as.matrix(expand.grid(u1 = u, u2 = u))
val <- cbind(grid, z = C.n(grid, X = u_graph))
cpCn <- contourplot2(val, region = FALSE, labels = FALSE, col = 3, oma = c(0,0,0,0),mar = c(2,2,2,1))
file <- paste0("plots/fit_t_",combos[i,1],"_",combos[i,2],".svg")
svg(file = file, width = 5, height = 5)
cpG + cpCn
dev.off()
i <- i + 1

i <- 1
u_graph <- u_res_all[,combos[i,]]
u_rand <- rCopula(n = 1000000, copula = cc@copula)
u <- seq(0, 1, length.out = 50)
grid <- as.matrix(expand.grid(u1 = u, u2 = u))
val <- cbind(grid, z = C.n(grid, X = u_rand[,combos[i,]]))
cpG <- contourplot2(val, region = FALSE,
                    key = list(corner = c(0.01, 0.01),
                               lines = list(col = c(1,4), lwd = 2),
                               text = list(c("Copula",
                                             "Data")), oma =c(0,0,0,0),mar = c(2,2,2,1)),
                    xlab= paste0("u",combos[i,1]),
                    ylab = paste0("u",combos[i,2]))
u <- seq(0, 1, length.out = 50)
grid <- as.matrix(expand.grid(u1 = u, u2 = u))
val <- cbind(grid, z = C.n(grid, X = u_graph))
cpCn <- contourplot2(val, region = FALSE, labels = FALSE, col = 4, oma = c(0,0,0,0),mar = c(2,2,2,1))
file <- paste0("plots/fit_clayton_",combos[i,1],"_",combos[i,2],".svg")
svg(file = file, width = 5, height = 5)
cpG + cpCn
dev.off()
i <- i + 1 

# 5. Simulation ####

n_sim <- 1000

u_boot <- rCopula(n = n_sim*52, copula = tc@copula)
colnames(u_boot) <- colnames(u_res_all)[1:4]

## 5.1 Interest rate model ----

### 5.1.1 PC1 ----
nPC <- 1
model <- model_dtsPC1
dtsPC <- diff_narm(tsPCX[,nPC])
pb <- txtProgressBar(1,n_sim)
for (i in 1:n_sim){
  pos <- (i-1)*52+(1:52)
  simulated <- sim_ugarch(model, u_boot[pos,nPC])
  complete <- c(simulated, dtsPC)
  series <- cumsum(complete) + coredata(tsPCX[1,nPC])[[1]]
  if (i == 1) series_sim <- series else series_sim <- cbind(series_sim, series)
  setTxtProgressBar(pb, i)
}

colnames(series_sim) <- paste0("sim_",1:n_sim)
series_sim_PC1 <- series_sim

### 5.1.2 PC2 ----
nPC <- 2
model <- model_dtsPC2
dtsPC <- diff_narm(tsPCX[,nPC])

for (i in 1:n_sim){
  pos <- (i-1)*52+(1:52)
  simulated <- sim_ugarch(model, u_boot[pos,nPC])
  complete <- c(simulated, dtsPC)
  series <- cumsum(complete) + coredata(tsPCX[1,nPC])[[1]]
  if (i == 1) series_sim <- series else series_sim <- cbind(series_sim, series)
  setTxtProgressBar(pb, i)
}

colnames(series_sim) <- paste0("sim_",1:n_sim)
series_sim_PC2 <- series_sim

### 5.1.3 PC3 ----
nPC <- 3
model <- model_dtsPC3
dtsPC <- diff_narm(tsPCX[,nPC])

for (i in 1:n_sim){
  pos <- (i-1)*52+(1:52)
  simulated <- sim_ugarch(model, u_boot[pos,nPC])
  complete <- c(simulated, dtsPC)
  series <- cumsum(complete) + coredata(tsPCX[1,nPC])[[1]]
  if (i == 1) series_sim <- series else series_sim <- cbind(series_sim, series)
  setTxtProgressBar(pb, i)
}

colnames(series_sim) <- paste0("sim_",1:n_sim)
series_sim_PC3 <- series_sim

## 5.2 EuroStoxx600 model ----

nPC <- 4
model <- model_dtsDJSTOXX
dtsPC <- diff_narm(DJSTOXX)/lag(DJSTOXX,1)
init <- coredata(DJSTOXX[-1][1])[[1]]
init <- init/(1+coredata(dtsPC[1])[[1]])

for (i in 1:n_sim){
  pos <- (i-1)*52+(1:52)
  simulated <- sim_ugarch(model, u_boot[pos,nPC])
  complete <- c(simulated, dtsPC)
  # in this case to undo is not simply cumulative sum since it is a percentage
  series <- cumprod(1+complete)*init
  if (i == 1) series_sim <- series else series_sim <- cbind(series_sim, series)
  setTxtProgressBar(pb, i)
}

colnames(series_sim) <- paste0("sim_",1:n_sim)
series_sim_DJSTOXX <- series_sim

## 5.3 Undo PC and organize data ----

market_df_sim <- mutate(market_df, type = "hist")
PC_sim <- matrix(data = 0, nrow = 52, ncol = 15)
market_sub_df <- market_df_sim[1:52,]
dates_sim <- index(tail(series_sim,52))

for (i in 1:n_sim) {
  # undo PC components to swaps curve
  PC1_sim <- series_sim_PC1[,i]
  PC2_sim <- series_sim_PC2[,i]
  PC3_sim <- series_sim_PC3[,i]
  PC_sim[,1:3] <- unlist(coredata(tail(cbind(PC1_sim,PC2_sim,PC3_sim),52)))
  # see PC decomposition code for reference
  swaps_sub <- PC_sim %*% t(eigen_swaps$vectors)
  swaps_sub <- unscale(swaps_sub, scale, center)
  # add also the EUROSTOXX series
  DJSTOXX_sim <- as.vector(unlist(coredata(tail(series_sim_DJSTOXX[,i],52))))
  # bind to already existing data
  market_sub_df2 <- as.data.frame(swaps_sub) %>% mutate(type = paste0("sim_",i))
  market_sub_df2$date <- dates_sim
  market_sub_df2$DJSTOXX <- DJSTOXX_sim
  market_sub_df[,colnames(market_sub_df2)] <- market_sub_df2[,colnames(market_sub_df2)]
  market_df_sim <- rbind(market_df_sim,market_sub_df)
  setTxtProgressBar(pb, i)
}

## 5.4 Plot results ----

### 5.4.1 Yield curve ----

yield <- market_df_sim %>% filter(date == (head(dates_sim,1)-7)) %>%
  subset(select = -c(USDEUR,type,date)) %>% as.matrix()
estim <- market_df_sim %>% filter(date == tail(dates_sim,1)) %>%
  subset(select = -c(USDEUR,type,date)) %>% as.matrix()

quantiles <- c(0.005, 0.05,0.1)

env_top <- apply(estim, 2, quantile, probs = quantiles)
env_bot <- apply(estim, 2, quantile, probs = 1-quantiles)

# art 166 and 167
up <- c(.7,.7,.64,.59,.55,.52,.49,.47,.44,.42,0.37,.33,.26,.26-5/70*.06,
        .26-10/70*.06)
down <- c(.75,.65,.56,.5,.46,.42,.39,.36,.33,.31,.29,.27,.29,.29-5/70*0.09,
          .29-10/70*0.09)

eiopa_up <- (1+up)*yield[1,1:length(times)]
eiopa_down <- (1-down)*yield[1,1:length(times)]

cols <- c("orange", "#6A7434", "#E0FF00")

svg(file = "plots/shocks_interest.svg", width = 17/1.4, height = 6.42/1.4)
par(mfrow=c(1,1), mar = c(4,4,3,1),oma = c(0,0,0,0))
plot(times, yield[1,1:15], type = "l", ylim = c(0,7), col = "#202426", lwd = "2",
     xlab = "maturity", ylab = "rate", main = "")
for (i in 1:nrow(env_top)){
  lines(times, env_top[i,1:15], type = "l", col = cols[i], lwd = 2)
  lines(times, env_bot[i,1:15], type = "l", col = cols[i], lwd = 2)
}

lines(times, eiopa_down, type = "l", lty = 2, col = cols[1], lwd = 3)
lines(times, eiopa_up, type = "l", lty = 2, col = cols[1], lwd = 3)
legend("top", legend = c("Model estimate / Current", "Model VaR10% envelope",
                           "Model VaR5% envelope", "Model VaR0.5% envelope",
                           "EIOPA VaR0.5%"),
       lwd = c(2,2,2,2,3), lty = c(1,1,1,1,3),
       col = c("#202426", cols[3], cols[2], cols[1], cols[1]), ncol = 3)
dev.off()

### 5.4.2 Equity ----

stock <- yield[,"DJSTOXX"]
stock_estim <- estim[,"DJSTOXX"]

stock_top <- quantile(stock_estim, quantiles)/stock
stock_bot <- quantile(stock_estim, 1-quantiles)/stock

# art 169
SA_interval <- DJSTOXX[paste(tail(dates,1)-3*365,tail(dates,1),sep ="::")]
SA <- 1/2*(stock-mean(SA_interval))/mean(SA_interval)
SA <- max(min(0.1,SA),-0.1)
eiopa_stock_up <- (1 + 0.39 + SA)
eiopa_stock_down <- (1 - 0.39 - SA)

svg(file = "plots/shocks_equity.svg", width = 17/1.4, height = 6.42/1.4)
par(mfrow=c(1,1), mar = c(4,4,3,1),oma = c(0,0,0,0))
hist(stock_estim/stock, breaks = 250, freq = FALSE,
     xlab = "Value",
     main = "",
     xlim = c(0.25,1.75), ylim = c(0,3))
lines(c(1,1), c(0,100), type = "l", col = "black", lwd = 2)
for (i in 1:nrow(env_top)){
  lines(c(stock_top[i],stock_top[i]), c(0,100), type = "l", col = cols[i], lwd = 2)
  lines(c(stock_bot[i],stock_bot[i]), c(0,100), type = "l", col = cols[i], lwd = 2)
}

lines(c(eiopa_stock_down,eiopa_stock_down), c(0,100), type = "l", col = cols[1],
      lwd = 3, lty = 2)
legend("right", legend = c("Model estimate", "Model VaR10%", "Model VaR5%",
                           "Model VaR0.5%", "EIOPA VaR0.5%"),
       lwd = c(2,2,2,2,3), lty = c(1,1,1,1,3),
       col = c("#202426", cols[3], cols[2], cols[1], cols[1]), ncol = 1)
dev.off()

### 5.5.1 Portfolio 1 calculation ----

portfolio <- list(bond = list(mat = 8, cupon = 0.02, qty = 1000),
                  bond = list(mat = 10, cupon = 0.04, qty = 1000),
                  bond = list(mat = 5, cupon = 0.03, qty = 1000),
                  bond = list(mat = 3, cupon = 0.01, qty = 1000),
                  stock = list(price = 500))

portfolio_value <- portfolio_value_calc(portfolio, yield[1:length(times)],
                                  times, 1, duration = TRUE)
portfolio_value <- c(portfolio_value, "total" = sum(portfolio_value[1:2]))
portfolio_eiopa_results <- portfolio_value_calc(portfolio, eiopa_up, times,
                                          eiopa_stock_down)

portfolio_eiopa_results2 <- portfolio_value_calc(portfolio, eiopa_down, times,
                                                 eiopa_stock_up)

#### 5.5.1.1 Standard formula ----

eiopa_cor <- matrix(c(1,0.5,0.5,1), nrow = 2, ncol = 2)
eiopa_cor2 <- matrix(c(1,0,0,1), nrow = 2, ncol = 2)

shocks_eiopa <- c(portfolio_eiopa_results["bonds"]-portfolio_value["bonds"],
                  portfolio_eiopa_results["stocks"]-portfolio_value["stocks"])

shocks_eiopa["total"] <- -sqrt((-shocks_eiopa) %*% eiopa_cor %*% (-shocks_eiopa))
names(shocks_eiopa) <- paste0("eiopa_shock_",c("bonds","stocks","total"))

shocks_eiopa2 <- c(portfolio_eiopa_results2["bonds"]-portfolio_value["bonds"],
                  portfolio_eiopa_results2["stocks"]-portfolio_value["stocks"])

shocks_eiopa2["total"] <- sqrt((-shocks_eiopa2) %*% eiopa_cor2 %*% (-shocks_eiopa2))
names(shocks_eiopa2) <- paste0("eiopa_shock_",c("bonds","stocks","total"))

#### 5.5.1.2 Stochastic model ----

for (i in 1:n_sim){
  portfolio_sim <- portfolio_value_calc(portfolio, estim[i,1:length(times)], times,
                                        stock_estim[i]/as.vector(stock))
  if (i == 1) portfolio_sim_results <- portfolio_sim else
  portfolio_sim_results <- rbind(portfolio_sim_results,portfolio_sim)
}

results_sim_portfolio <- market_df_sim %>% filter(date == tail(dates_sim,1)) %>%
  mutate(current_bond_value = portfolio_value["bonds"],
         current_stock_value = portfolio_value["stocks"],
         bonds_duration = portfolio_value["duration"])

results_sim_portfolio$stochastic_bonds_value <- portfolio_sim_results[,"bonds"]
results_sim_portfolio$stochastic_stocks_value <- portfolio_sim_results[,"stocks"]

results_sim_portfolio <- results_sim_portfolio %>%
  mutate(bonds_change = stochastic_bonds_value - current_bond_value,
         stocks_change = stochastic_stocks_value - current_stock_value,
         portfolio_change = bonds_change + stocks_change)

vaR <- quantile(results_sim_portfolio$portfolio_change,c(0.005,0.995))

svg(file = "plots/histogram_portfolio1.svg", width = 10, height = 5)
par(mfrow=c(1,2), mar = c(4.5,4,3,1),oma = c(0,0,0,0))
hist(results_sim_portfolio$bonds_change/portfolio_value["total"],60,
     main = "bonds shock portfolio change", xlab = "portfolio loss")
hist(results_sim_portfolio$stocks_change/portfolio_value["total"],60,
     main = "equity shock portfolio change", xlab = "portfolio loss")
dev.off()
par(mfrow=c(1,1), mar = c(4.5,4,3,1),oma = c(0,0,0,0))

varnames <- c(expression("Shock"["Interest rate"]), expression("Shock"["Equity"]))
u_stock_bonds <- cbind(pobs(results_sim_portfolio$bonds_change),
                       pobs(results_sim_portfolio$stocks_change))

cor(results_sim_portfolio$bonds_change,results_sim_portfolio$stocks_change)

svg(file = "plots/corel_portfolio1.svg", width = 5, height = 5)
pairs(head(u_stock_bonds,1219),
      labels = varnames, cex = 0.5, oma = c(2,2,2,2))
dev.off()

VaR_rates <- quantile(results_sim_portfolio$bonds_change,c(0.005,0.995))
VaR_stock <- quantile(results_sim_portfolio$stocks_change,c(0.005,0.995))

#### 5.5.1.3 Results ----

portfolio1 <- c(portfolio_value,shocks_eiopa,stoch_shock_bonds = VaR_rates[1],
                stoch_shock_stocks = VaR_stock[1], stoch_shock_total = vaR[1])

##### Waterfall diagram Shock Model down ----
range <- min(VaR_rates[1]+VaR_stock[1],shocks_eiopa[1]+shocks_eiopa[2])/portfolio_value["total"]
range <- sign(range)*ceiling(abs(range*10))/10
shocks <- c(VaR_rates[1],VaR_stock[1], vaR[1])/portfolio_value["total"]
x <- list("Interest Shock up", "Equity Shock down", "Diversification", "Aggregated shock")
measure <- c("relative", "relative", "relative", "total")
y <- c(shocks[1], shocks[2], shocks[3]-sum(shocks[1:2]),0)
text <- as.character(abs(round(c(y[1:3],shocks[3]),3)))
signs <- c("-","","+")[sign(c(y[1:3],shocks[3]))+2]
text <- paste(signs, text, sep = "")
data <- data.frame(x=factor(x,levels=x),measure,text,y)
plot_ly(data, name = "", type = "waterfall", measure = ~measure,
  x = ~x, textposition = "outside", y= ~y, text =~text,
  decreasing = list(marker = list(color = "#202426")),
  increasing = list(marker = list(color = "#6A7434")),
  totals = list(marker = list(color = "#8C8C88")),
  connector = list(line = list(color= "rgb(63, 63, 63)"))) %>%
  layout(title = "Model shocks aggregation (Capital loss)",
         xaxis = list(title = ""),
         yaxis = list(title = "",range = list(range, 0)),
         autosize = TRUE,
         showlegend = FALSE)

##### Waterfall diagram Shock EIOPA down ----
shocks <- c(shocks_eiopa[1],shocks_eiopa[2], shocks_eiopa[3])/portfolio_value["total"]
x <- list("Interest Shock up", "Equity Shock down", "Diversification", "Aggregated shock")
measure <- c("relative", "relative", "relative", "total")
y <- c(shocks[1], shocks[2], shocks[3]-sum(shocks[1:2]),0)
text <- as.character(abs(round(c(y[1:3],shocks[3]),3)))
signs <- c("-","","+")[sign(c(y[1:3],shocks[3]))+2]
text <- paste(signs, text, sep = "")
data <- data.frame(x=factor(x,levels=x),measure,text,y)
plot_ly(data, name = "", type = "waterfall", measure = ~measure,
        x = ~x, textposition = "outside", y= ~y, text =~text,
        decreasing = list(marker = list(color = "#202426")),
        increasing = list(marker = list(color = "#6A7434")),
        totals = list(marker = list(color = "#8C8C88")),
        connector = list(line = list(color= "rgb(63, 63, 63)"))) %>%
  layout(title = "EIOPA shocks aggregation (Capital loss)",
         xaxis = list(title = ""),
         yaxis = list(title = "",range = list(range, 0)),
         autosize = TRUE,
         showlegend = FALSE)

##### Waterfall diagram Shock Model up ----
range <- max(VaR_rates[2]+VaR_stock[2],shocks_eiopa2[1]+shocks_eiopa2[2])/portfolio_value["total"]
shocks <- c(VaR_rates[2],VaR_stock[2], vaR[2])/portfolio_value["total"]
x <- list("Interest Shock down", "Equity Shock up", "Diversification", "Aggregated shock")
measure <- c("relative", "relative", "relative", "total")
y <- c(shocks[1], shocks[2], shocks[3]-sum(shocks[1:2]),0)
text <- as.character(abs(round(c(y[1:3],shocks[3]),3)))
signs <- c("-","","+")[sign(c(y[1:3],shocks[3]))+2]
text <- paste(signs, text, sep = "")
data <- data.frame(x=factor(x,levels=x),measure,text,y)
plot_ly(data, name = "", type = "waterfall", measure = ~measure,
        x = ~x, textposition = "outside", y= ~y, text =~text,
        connector = list(line = list(color= "rgb(63, 63, 63)"))) %>%
  layout(title = "Model shocks aggregation (Capital increase)",
         xaxis = list(title = ""),
         yaxis = list(title = "",range = list(0, range)),
         autosize = TRUE,
         showlegend = FALSE)

##### Waterfall diagram Shock EIOPA up ----
shocks <- c(shocks_eiopa2[1],shocks_eiopa2[2], shocks_eiopa2[3])/portfolio_value["total"]
x <- list("Interest Shock down", "Equity Shock up", "Diversification", "Aggregated shock")
measure <- c("relative", "relative", "relative", "total")
y <- c(shocks[1], shocks[2], shocks[3]-sum(shocks[1:2]),0)
text <- as.character(abs(round(c(y[1:3],shocks[3]),3)))
signs <- c("-","","+")[sign(c(y[1:3],shocks[3]))+2]
text <- paste(signs, text, sep = "")
data <- data.frame(x=factor(x,levels=x),measure,text,y)
plot_ly(data, name = "", type = "waterfall", measure = ~measure,
        x = ~x, textposition = "outside", y= ~y, text =~text,
        connector = list(line = list(color= "rgb(63, 63, 63)"))) %>%
  layout(title = "EIOPA shocks aggregation (Capital increase)",
         xaxis = list(title = ""),
         yaxis = list(title = "",range = list(0, range)),
         autosize = TRUE,
         showlegend = FALSE)

### 5.5.2 Portfolio 2 calculation ----

portfolio <- list(bond = list(mat = 10, cupon = 0.02, qty = 1000),
                  bond = list(mat = 15, cupon = 0.04, qty = 1000),
                  bond = list(mat = 20, cupon = 0.03, qty = 1000),
                  bond = list(mat = 5, cupon = 0.01, qty = 1000),
                  stock = list(price = 500))

portfolio_value <- portfolio_value_calc(portfolio, yield[1:length(times)],
                                        times, 1, duration = TRUE)
portfolio_value <- c(portfolio_value, "total" = sum(portfolio_value[1:2]))
portfolio_eiopa_results <- portfolio_value_calc(portfolio, eiopa_up, times,
                                                eiopa_stock_down)

portfolio_eiopa_results2 <- portfolio_value_calc(portfolio, eiopa_down, times,
                                                 eiopa_stock_up)

#### 5.5.1.1 Standard formula ----

eiopa_cor <- matrix(c(1,0.5,0.5,1), nrow = 2, ncol = 2)
eiopa_cor2 <- matrix(c(1,0,0,1), nrow = 2, ncol = 2)

shocks_eiopa <- c(portfolio_eiopa_results["bonds"]-portfolio_value["bonds"],
                  portfolio_eiopa_results["stocks"]-portfolio_value["stocks"])

shocks_eiopa["total"] <- -sqrt((-shocks_eiopa) %*% eiopa_cor %*% (-shocks_eiopa))
names(shocks_eiopa) <- paste0("eiopa_shock_",c("bonds","stocks","total"))

shocks_eiopa2 <- c(portfolio_eiopa_results2["bonds"]-portfolio_value["bonds"],
                   portfolio_eiopa_results2["stocks"]-portfolio_value["stocks"])

shocks_eiopa2["total"] <- sqrt((-shocks_eiopa2) %*% eiopa_cor2 %*% (-shocks_eiopa2))
names(shocks_eiopa2) <- paste0("eiopa_shock_",c("bonds","stocks","total"))

#### 5.5.1.2 Stochastic model ----

for (i in 1:n_sim){
  portfolio_sim <- portfolio_value_calc(portfolio, estim[i,1:length(times)], times,
                                        stock_estim[i]/as.vector(stock))
  if (i == 1) portfolio_sim_results <- portfolio_sim else
    portfolio_sim_results <- rbind(portfolio_sim_results,portfolio_sim)
}

results_sim_portfolio <- market_df_sim %>% filter(date == tail(dates_sim,1)) %>%
  mutate(current_bond_value = portfolio_value["bonds"],
         current_stock_value = portfolio_value["stocks"],
         bonds_duration = portfolio_value["duration"])

results_sim_portfolio$stochastic_bonds_value <- portfolio_sim_results[,"bonds"]
results_sim_portfolio$stochastic_stocks_value <- portfolio_sim_results[,"stocks"]

results_sim_portfolio <- results_sim_portfolio %>%
  mutate(bonds_change = stochastic_bonds_value - current_bond_value,
         stocks_change = stochastic_stocks_value - current_stock_value,
         portfolio_change = bonds_change + stocks_change)

vaR <- quantile(results_sim_portfolio$portfolio_change,c(0.005,0.995))

svg(file = "plots/histogram_portfolio2.svg", width = 10, height = 5)
par(mfrow=c(1,2), mar = c(4.5,4,3,1),oma = c(0,0,0,0))
hist(results_sim_portfolio$bonds_change/portfolio_value["total"],60,
     main = "bonds shock portfolio change", xlab = "portfolio loss")
hist(results_sim_portfolio$stocks_change/portfolio_value["total"],60,
     main = "equity shock portfolio change", xlab = "portfolio loss")
dev.off()
par(mfrow=c(1,1), mar = c(4.5,4,3,1),oma = c(0,0,0,0))

varnames <- c(expression("Shock"["Interest rate"]), expression("Shock"["Equity"]))
u_stock_bonds <- cbind(pobs(results_sim_portfolio$bonds_change),
                       pobs(results_sim_portfolio$stocks_change))

cor(results_sim_portfolio$bonds_change,results_sim_portfolio$stocks_change)

svg(file = "plots/corel_portfolio2.svg", width = 5, height = 5)
pairs(head(u_stock_bonds,1219),
      labels = varnames, cex = 0.5, oma = c(2,2,2,2))
dev.off()

VaR_rates <- quantile(results_sim_portfolio$bonds_change,c(0.005,0.995))
VaR_stock <- quantile(results_sim_portfolio$stocks_change,c(0.005,0.995))

#### 5.5.2.3 Results ----

portfolio2 <- c(portfolio_value,shocks_eiopa,stoch_shock_bonds = VaR_rates[1],
                stoch_shock_stocks = VaR_stock[1], stoch_shock_total = vaR[1])

##### Waterfall diagram Shock Model down ----
range <- min(VaR_rates[1]+VaR_stock[1],shocks_eiopa[1]+shocks_eiopa[2])/portfolio_value["total"]
range <- sign(range)*ceiling(abs(range*10))/10
shocks <- c(VaR_rates[1],VaR_stock[1], vaR[1])/portfolio_value["total"]
x <- list("Interest Shock up", "Equity Shock down", "Diversification", "Aggregated shock")
measure <- c("relative", "relative", "relative", "total")
y <- c(shocks[1], shocks[2], shocks[3]-sum(shocks[1:2]),0)
text <- as.character(abs(round(c(y[1:3],shocks[3]),3)))
signs <- c("-","","+")[sign(c(y[1:3],shocks[3]))+2]
text <- paste(signs, text, sep = "")
data <- data.frame(x=factor(x,levels=x),measure,text,y)
plot_ly(data, name = "", type = "waterfall", measure = ~measure,
        x = ~x, textposition = "outside", y= ~y, text =~text,
        decreasing = list(marker = list(color = "#202426")),
        increasing = list(marker = list(color = "#6A7434")),
        totals = list(marker = list(color = "#8C8C88")),
        connector = list(line = list(color= "rgb(63, 63, 63)"))) %>%
  layout(title = "Model shocks aggregation (Capital loss)",
         xaxis = list(title = ""),
         yaxis = list(title = "",range = list(range, 0)),
         autosize = TRUE,
         showlegend = FALSE)

##### Waterfall diagram Shock EIOPA down ----
shocks <- c(shocks_eiopa[1],shocks_eiopa[2], shocks_eiopa[3])/portfolio_value["total"]
x <- list("Interest Shock up", "Equity Shock down", "Diversification", "Aggregated shock")
measure <- c("relative", "relative", "relative", "total")
y <- c(shocks[1], shocks[2], shocks[3]-sum(shocks[1:2]),0)
text <- as.character(abs(round(c(y[1:3],shocks[3]),3)))
signs <- c("-","","+")[sign(c(y[1:3],shocks[3]))+2]
text <- paste(signs, text, sep = "")
data <- data.frame(x=factor(x,levels=x),measure,text,y)
plot_ly(data, name = "", type = "waterfall", measure = ~measure,
        x = ~x, textposition = "outside", y= ~y, text =~text,
        connector = list(line = list(color= "rgb(63, 63, 63)"))) %>%
  layout(title = "EIOPA shocks aggregation (Capital loss)",
         xaxis = list(title = ""),
         yaxis = list(title = "",range = list(range, 0)),
         autosize = TRUE,
         showlegend = FALSE)

##### Waterfall diagram Shock Model up ----
range <- max(VaR_rates[2]+VaR_stock[2],shocks_eiopa2[1]+shocks_eiopa2[2])/portfolio_value["total"]
range <- sign(range)*ceiling(abs(range*10))/10
shocks <- c(VaR_rates[2],VaR_stock[2], vaR[2])/portfolio_value["total"]
x <- list("Interest Shock down", "Equity Shock up", "Diversification", "Aggregated shock")
measure <- c("relative", "relative", "relative", "total")
y <- c(shocks[1], shocks[2], shocks[3]-sum(shocks[1:2]),0)
text <- as.character(abs(round(c(y[1:3],shocks[3]),3)))
signs <- c("-","","+")[sign(c(y[1:3],shocks[3]))+2]
text <- paste(signs, text, sep = "")
data <- data.frame(x=factor(x,levels=x),measure,text,y)
plot_ly(data, name = "", type = "waterfall", measure = ~measure,
        x = ~x, textposition = "outside", y= ~y, text =~text,
        connector = list(line = list(color= "rgb(63, 63, 63)"))) %>%
  layout(title = "Model shocks aggregation (Capital increase)",
         xaxis = list(title = ""),
         yaxis = list(title = "",range = list(0, range)),
         autosize = TRUE,
         showlegend = FALSE)

##### Waterfall diagram Shock EIOPA up ----
shocks <- c(shocks_eiopa2[1],shocks_eiopa2[2], shocks_eiopa2[3])/portfolio_value["total"]
x <- list("Interest Shock down", "Equity Shock up", "Diversification", "Aggregated shock")
measure <- c("relative", "relative", "relative", "total")
y <- c(shocks[1], shocks[2], shocks[3]-sum(shocks[1:2]),0)
text <- as.character(abs(round(c(y[1:3],shocks[3]),3)))
signs <- c("-","","+")[sign(c(y[1:3],shocks[3]))+2]
text <- paste(signs, text, sep = "")
data <- data.frame(x=factor(x,levels=x),measure,text,y)
plot_ly(data, name = "", type = "waterfall", measure = ~measure,
        x = ~x, textposition = "outside", y= ~y, text =~text,
        connector = list(line = list(color= "rgb(63, 63, 63)"))) %>%
  layout(title = "EIOPA shocks aggregation (Capital increase)",
         xaxis = list(title = ""),
         yaxis = list(title = "",range = list(0, range)),
         autosize = TRUE,
         showlegend = FALSE)

### 5.5.3 Portfolio 3 calculation ----

portfolio <- list(bond = list(mat = 15, cupon = 0.02, qty = 1000),
                  bond = list(mat = 20, cupon = 0.04, qty = 1000),
                  bond = list(mat = 25, cupon = 0.03, qty = 1000),
                  bond = list(mat = 10, cupon = 0.01, qty = 1000),
                  stock = list(price = 500))

portfolio_value <- portfolio_value_calc(portfolio, yield[1:length(times)],
                                        times, 1, duration = TRUE)
portfolio_value <- c(portfolio_value, "total" = sum(portfolio_value[1:2]))
portfolio_eiopa_results <- portfolio_value_calc(portfolio, eiopa_up, times,
                                                eiopa_stock_down)

portfolio_eiopa_results2 <- portfolio_value_calc(portfolio, eiopa_down, times,
                                                 eiopa_stock_up)

#### 5.5.1.1 Standard formula ----

eiopa_cor <- matrix(c(1,0.5,0.5,1), nrow = 2, ncol = 2)
eiopa_cor2 <- matrix(c(1,0,0,1), nrow = 2, ncol = 2)

shocks_eiopa <- c(portfolio_eiopa_results["bonds"]-portfolio_value["bonds"],
                  portfolio_eiopa_results["stocks"]-portfolio_value["stocks"])

shocks_eiopa["total"] <- -sqrt((-shocks_eiopa) %*% eiopa_cor %*% (-shocks_eiopa))
names(shocks_eiopa) <- paste0("eiopa_shock_",c("bonds","stocks","total"))

shocks_eiopa2 <- c(portfolio_eiopa_results2["bonds"]-portfolio_value["bonds"],
                   portfolio_eiopa_results2["stocks"]-portfolio_value["stocks"])

shocks_eiopa2["total"] <- sqrt((-shocks_eiopa2) %*% eiopa_cor2 %*% (-shocks_eiopa2))
names(shocks_eiopa2) <- paste0("eiopa_shock_",c("bonds","stocks","total"))

#### 5.5.1.2 Stochastic model ----

for (i in 1:n_sim){
  portfolio_sim <- portfolio_value_calc(portfolio, estim[i,1:length(times)], times,
                                        stock_estim[i]/as.vector(stock))
  if (i == 1) portfolio_sim_results <- portfolio_sim else
    portfolio_sim_results <- rbind(portfolio_sim_results,portfolio_sim)
}

results_sim_portfolio <- market_df_sim %>% filter(date == tail(dates_sim,1)) %>%
  mutate(current_bond_value = portfolio_value["bonds"],
         current_stock_value = portfolio_value["stocks"],
         bonds_duration = portfolio_value["duration"])

results_sim_portfolio$stochastic_bonds_value <- portfolio_sim_results[,"bonds"]
results_sim_portfolio$stochastic_stocks_value <- portfolio_sim_results[,"stocks"]

results_sim_portfolio <- results_sim_portfolio %>%
  mutate(bonds_change = stochastic_bonds_value - current_bond_value,
         stocks_change = stochastic_stocks_value - current_stock_value,
         portfolio_change = bonds_change + stocks_change)

vaR <- quantile(results_sim_portfolio$portfolio_change,c(0.005,0.995))

svg(file = "plots/histogram_portfolio3.svg", width = 10, height = 5)
par(mfrow=c(1,2), mar = c(4.5,4,3,1),oma = c(0,0,0,0))
hist(results_sim_portfolio$bonds_change/portfolio_value["total"],60,
     main = "bonds shock portfolio change", xlab = "portfolio loss")
hist(results_sim_portfolio$stocks_change/portfolio_value["total"],60,
     main = "equity shock portfolio change", xlab = "portfolio loss")
dev.off()
par(mfrow=c(1,1), mar = c(4.5,4,3,1),oma = c(0,0,0,0))

varnames <- c(expression("Shock"["Interest rate"]), expression("Shock"["Equity"]))
u_stock_bonds <- cbind(pobs(results_sim_portfolio$bonds_change),
                       pobs(results_sim_portfolio$stocks_change))

cor(results_sim_portfolio$bonds_change,results_sim_portfolio$stocks_change)

svg(file = "plots/corel_portfolio3.svg", width = 5, height = 5)
pairs(head(u_stock_bonds,1219),
      labels = varnames, cex = 0.5, oma = c(2,2,2,2))
dev.off()

VaR_rates <- quantile(results_sim_portfolio$bonds_change,c(0.005,0.995))
VaR_stock <- quantile(results_sim_portfolio$stocks_change,c(0.005,0.995))

#### 5.5.3.3 Results ----

portfolio3 <- c(portfolio_value,shocks_eiopa,stoch_shock_bonds = VaR_rates[1],
                stoch_shock_stocks = VaR_stock[1], stoch_shock_total = vaR[1])

##### Waterfall diagram Shock Model down ----
range <- min(VaR_rates[1]+VaR_stock[1],shocks_eiopa[1]+shocks_eiopa[2])/portfolio_value["total"]
range <- sign(range)*ceiling(abs(range*10))/10
shocks <- c(VaR_rates[1],VaR_stock[1], vaR[1])/portfolio_value["total"]
x <- list("Interest Shock up", "Equity Shock down", "Diversification", "Aggregated shock")
measure <- c("relative", "relative", "relative", "total")
y <- c(shocks[1], shocks[2], shocks[3]-sum(shocks[1:2]),0)
text <- as.character(abs(round(c(y[1:3],shocks[3]),3)))
signs <- c("-","","+")[sign(c(y[1:3],shocks[3]))+2]
text <- paste(signs, text, sep = "")
data <- data.frame(x=factor(x,levels=x),measure,text,y)
plot_ly(data, name = "", type = "waterfall", measure = ~measure,
        x = ~x, textposition = "outside", y= ~y, text =~text,
        connector = list(line = list(color= "rgb(63, 63, 63)"))) %>%
  layout(title = "Model shocks aggregation (Capital loss)",
         xaxis = list(title = ""),
         yaxis = list(title = "",range = list(range, 0)),
         autosize = TRUE,
         showlegend = FALSE)

##### Waterfall diagram Shock EIOPA down ----
shocks <- c(shocks_eiopa[1],shocks_eiopa[2], shocks_eiopa[3])/portfolio_value["total"]
x <- list("Interest Shock up", "Equity Shock down", "Diversification", "Aggregated shock")
measure <- c("relative", "relative", "relative", "total")
y <- c(shocks[1], shocks[2], shocks[3]-sum(shocks[1:2]),0)
text <- as.character(abs(round(c(y[1:3],shocks[3]),3)))
signs <- c("-","","+")[sign(c(y[1:3],shocks[3]))+2]
text <- paste(signs, text, sep = "")
data <- data.frame(x=factor(x,levels=x),measure,text,y)
plot_ly(data, name = "", type = "waterfall", measure = ~measure,
        x = ~x, textposition = "outside", y= ~y, text =~text,
        connector = list(line = list(color= "rgb(63, 63, 63)"))) %>%
  layout(title = "EIOPA shocks aggregation (Capital loss)",
         xaxis = list(title = ""),
         yaxis = list(title = "",range = list(range, 0)),
         autosize = TRUE,
         showlegend = FALSE)

##### Waterfall diagram Shock Model up ----
range <- max(VaR_rates[2]+VaR_stock[2],shocks_eiopa2[1]+shocks_eiopa2[2])/portfolio_value["total"]
range <- sign(range)*ceiling(abs(range*10))/10
shocks <- c(VaR_rates[2],VaR_stock[2], vaR[2])/portfolio_value["total"]
x <- list("Interest Shock down", "Equity Shock up", "Diversification", "Aggregated shock")
measure <- c("relative", "relative", "relative", "total")
y <- c(shocks[1], shocks[2], shocks[3]-sum(shocks[1:2]),0)
text <- as.character(abs(round(c(y[1:3],shocks[3]),3)))
signs <- c("-","","+")[sign(c(y[1:3],shocks[3]))+2]
text <- paste(signs, text, sep = "")
data <- data.frame(x=factor(x,levels=x),measure,text,y)
plot_ly(data, name = "", type = "waterfall", measure = ~measure,
        x = ~x, textposition = "outside", y= ~y, text =~text,
        connector = list(line = list(color= "rgb(63, 63, 63)"))) %>%
  layout(title = "Model shocks aggregation (Capital increase)",
         xaxis = list(title = ""),
         yaxis = list(title = "",range = list(0, range)),
         autosize = TRUE,
         showlegend = FALSE)

##### Waterfall diagram Shock EIOPA up ----
shocks <- c(shocks_eiopa2[1],shocks_eiopa2[2], shocks_eiopa2[3])/portfolio_value["total"]
x <- list("Interest Shock down", "Equity Shock up", "Diversification", "Aggregated shock")
measure <- c("relative", "relative", "relative", "total")
y <- c(shocks[1], shocks[2], shocks[3]-sum(shocks[1:2]),0)
text <- as.character(abs(round(c(y[1:3],shocks[3]),3)))
signs <- c("-","","+")[sign(c(y[1:3],shocks[3]))+2]
text <- paste(signs, text, sep = "")
data <- data.frame(x=factor(x,levels=x),measure,text,y)
plot_ly(data, name = "", type = "waterfall", measure = ~measure,
        x = ~x, textposition = "outside", y= ~y, text =~text,
        connector = list(line = list(color= "rgb(63, 63, 63)"))) %>%
  layout(title = "EIOPA shocks aggregation (Capital increase)",
         xaxis = list(title = ""),
         yaxis = list(title = "",range = list(0, range)),
         autosize = TRUE,
         showlegend = FALSE)

## 5.5.4 Summary ----

portfolios <- rbind(portfolio1, portfolio2, portfolio3)

write.xlsx(as.data.frame(rbind(portfolios, portfolios_tstudent)),"data.xlsx")
portfolios_tstudent

cols <- rainbow(5)
for (i in 1:3){
  barplot(abs(portfolios[i,c(1:2,4:10)]), col = cols[c(2,2,1,3,3,1,4,4,1)],
          main = paste0("bond duration = ", round(portfolios[i,"duration"],2)))
}

# Some plots ----
swaps_df3 %>% mutate(coding = fct_reorder(name, sort(time))) %>%
  ggplot() + geom_line(mapping = aes(x = date, y = rate, color = coding))

swaps_df4 %>% ggplot() + geom_line(mapping = aes(x = date, y = TREUR1Y))

bond_all_df %>% filter(date == init_date) %>% ggplot() + 
  geom_line(mapping = aes(x = maturity, y = rate))


bond_all_df$date


dates <- rev((bond_all_df %>% filter(date > init_date, maturity == 1))$date)

for (i in dates){
  sum(bond_all_df$date == i)
}

for(i in dates){
  bond_all_df %>%
    ggplot(aes(times, rate)) +
    geom_line() +
    scale_x_continuous() +
    theme_bw() +
    labs(title = paste0("Year: ",2020), x = 'rate', y = 'maturity') 
}
