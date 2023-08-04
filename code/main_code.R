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
# one: periodic square with a little randomness
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


# observation error -------------------------------------------------------


# ar coefficient  ---------------------------------------------------------



