# simulate time series for ml forcasting model
# author(s): caitlin allen akselrud
# contact: caitlin.allen_akserud@noaa.gov
# date created: 13.07.2023
# version: 1.0

# libraries ---------------------------------------------------------------
library(tidyverse)
library(here)
library(stats)
library(MARSS)
library(forecast)
library(datasets)
library(TMB)

# functions ---------------------------------------------------------------


# folders -----------------------------------------------------------------

write.dir <- here::here("output") # otherwise will be saved in working directory
dir.create(write.dir, showWarnings = FALSE)
# setwd(write.dir)

# settings and fixed values -----------------------------------------------

t_length <- 100
t <- 1:t_length
dat_length <- 30
burn <- 100-dat_length

# environmental -----------------------------------------------------------
# one: periodic square
set.seed(1011)
period = 8
env1_square <- ifelse(((t %% period) < (0.5*period)),1,0) * rbinom(t_length, 1, 0.9)
# env <- 2*(2*(floor(1/period*t))-floor(2*1/period*t)) +1 #or this- same
plot(env1_square, type = 'l')
plot(env1_square[burn:t_length], type = 'l')
acf(env1_square)
pacf(env1_square)

# two: amplified signal
env <- sin((2*pi*t)/(2/(t)))
plot(env, type = 'l')
plot(env[burn:t_length], type = 'l')
acf(env)
pacf(env)

# three: strong ar (feedback)
set.seed(1234)
AR_lg <- list(order = c(1, 0, 0), ar = 0.9)
AR1_lg <- arima.sim(n = t_length, model = AR_lg, sd = 0.1)
plot(AR1_lg, type = 'l')
plot(AR1_lg[burn:t_length], type = 'l')
acf(AR1_lg)
pacf(AR1_lg)

# four: ar and ma
# cps-like/ prey index that's more env driven
set.seed(9876)
AR_ma <- list(order = c(1, 0, 1), ar = -0.1, ma = -0.1)
AR1_ma <- arima.sim(n = t_length, model = AR_ma, sd = 0.1)
plot(AR1_ma, type = 'l')
plot(AR1_ma[burn:t_length], type = 'l')
acf(AR1_ma)
pacf(AR1_ma)

# five: non-stationary (trend)
set.seed(5678)
ts_ns <- list(order = c(1, 1, 0), ar = -0.5)
ts_ns1 <- arima.sim(n = t_length, model = ts_ns, sd = 0.1) *-1
plot(ts_ns1, type = 'l')
plot(ts_ns1[burn:t_length], type = 'l')
acf(ts_ns1)
pacf(ts_ns1)

env_data <- data.frame(#time = seq(from = 1, to = dat_length, by=1),
  regime = env1_square[(burn + 1):t_length],
  signal = env[(burn + 1):t_length],
  climate = ts_ns1[(burn + 1):t_length],
  pred = AR1_lg[(burn + 1):t_length],
  prey = AR1_ma[(burn + 1):t_length])

m_env <- t(env_data)

colnames(m_env) <- seq(from = lubridate::year(Sys.Date())-dat_length+1,
                       to = lubridate::year(Sys.Date()),
                       by = 1)

n <- nrow(m_env) - 1

# initial state -----------------------------------------------------------
# x_0 = initial state
# x_t = hidden state at each time step
#       includes process error variance
# y_t = current state at each time step (data)
#       includes measurement error variance


# process error -----------------------------------------------------------

## function for simulating various autoregressive process models
proc_sim <- function(t = 30, B = 1, u = 0, CC = matrix(0), cc = matrix(0, ncol = t), Q = 1) {
  ## process model is defined by
  ##
  ## t = 1: x_1 ~ N(0, Q)
  ##
  ## t > 1: x_t = b x_{t-1} + u + CC %*% cc_t + w_t
  ##
  ## `t` is length of desired time series
  ## `B` is the AR(1) coef; default `ar = 1` gives a random walk
  ## `u` is the bias/drift or instantaneous growth coef; default `u = 0` gives unbiased random walk
  ## `CC` is a [1 x p] matrix of effect sizes for p covariates
  ## `cc` is a [p x t] matrix of covariates
  ## `Q` is the process variance; w_t ~ N(0, Q)
  ##
  ## ERROR checks
  ## limits on B
  if(abs(B) > 1) {
    stop("Setting `|B| > 1` will cause extreme booms/busts.")
  }
  ## dims for covariates
  if(ncol(CC) != nrow(cc)) {
    stop("The number of cols in `CC` must equal the number of rows in `cc`.")
  }
  ## positive variance
  if(Q <= 0) {
    stop("The process variance `Q` must be a positive number.")
    }
  ##
  ## initialize states and process errors
  ## t = 1 drawn from normal distribution
  xx <- ww <- rnorm(t, 0, Q)
  ## simulate process for t > 1
  for(i in 2:t) {
    xx[i] <- B * xx[i-1] + u + CC %*% cc[,i] + ww[i]
  }
  return(xx)
}


# observation error -------------------------------------------------------

obs_sim <- function(xx, a = 0, DD = matrix(0), dd = matrix(0, ncol = length(xx)), R = 1) {
  ## observation model is defined by
  ##
  ## y_t = x_t + a + DD %*% dd_t + v_t
  ##
  ## `xx` is simulated process (state)
  ## `a` is an offset
  ## `DD` is a [1 x m] matrix of effect sizes for p covariates
  ## `dd` is a [m x t] matrix of covariates
  ## `R` is the obs variance; v_t ~ N(0, R)
  ##
  ## ERROR checks
  ## dims for covariates
  if(ncol(DD) != nrow(dd)) {
    stop("The number of cols in `CC` must equal the number of rows in `cc`.")
  }
  ## positive variance
  if(R <= 0) {
    stop("The observation variance `R` must be a positive number.")
  }
  ##
  ## create obs errors
  vv <- rnorm(length(xx), 0, R)
  ## add obs error
  yy <- xx + a + DD %*% dd + vv
  return(as.vector(yy))
}


## simple plotting function ------------------------------------------------

plot_ssm <- function(state, obs) {
  plot.ts(obs, ylim = range(state, obs), las = 1, type = "n",
          ylab = expression(paste(italic(x[t])," or ", italic(y[t]))))
  lines(seq(length(state)), state,
        type = "o", pch = 16)
  lines(seq(length(obs)), obs, col = "blue",
        type = "o", pch = 16)
}


## examples ---------------------------------------------------------------

## random walk
xx <- proc_sim()
yy <- obs_sim(xx)
plot_ssm(xx, yy) #black = "state", blue = "observed"

## biased random walk with offset obs
xx <- proc_sim(u = 1)
yy <- obs_sim(xx, a = 10)
plot_ssm(xx, yy)

## stationary AR(1) with Var(v_t) = 2; high obsv variance
xx <- proc_sim(B = 0.5)
yy <- obs_sim(xx, R = 2)
plot_ssm(xx, yy)

## stationary AR(1) with sinusoidal covariate
xx <- proc_sim(B = 0.5,
               CC = matrix(1),
               cc = matrix(sin(2 * pi * seq(30) / 15), nrow = 1))
yy <- obs_sim(xx)
plot_ssm(xx, yy)


# ar coefficient  ---------------------------------------------------------

xx <-  proc_sim(t = dat_length,
                B = 1,                     ## `B` is the AR(1) coef
                u = 0,                     ## `u` is the bias/drift or instantaneous growth coef
                CC = t(matrix(c(1, .5, .25, .7, 2))),            ## `CC` is a [1 x p] matrix of effect sizes for p covariates
                cc = m_env,                ## `cc` is a [p x t] matrix of covariates
                Q = 1)                     ## `Q` is the process variance; w_t ~ N(0, Q)

yy <- obs_sim(xx,
              a = 0,                       ## `a` is an offset
              DD = matrix(0),              ## `DD` is a [1 x m] matrix of effect sizes for p covariates
              dd = matrix(0, ncol = length(xx)),   ## `dd` is a [m x t] matrix of covariates
              R = 1)                       ## `R` is the obs variance; v_t ~ N(0, R)

plot_ssm(xx, yy)


