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
burn <- 100-30


# initial state -----------------------------------------------------------


# process error -----------------------------------------------------------


# observation error -------------------------------------------------------


# ar coefficient  ---------------------------------------------------------


# environmental -----------------------------------------------------------

# one: periodic square
period = 8
env1_square <- ifelse(((t %% period) < (0.5*period)),1,0)
# env <- 2*(2*(floor(1/period*t))-floor(2*1/period*t)) +1 #or this- same
plot(env1_square, type = 'l')
acf(env1_square)
pacf(env1_square)

# two: amplified signal
env <- sin((2*pi*t)/(2/(t)))
plot(env, type = 'l')
plot(env[burn:t_length], type = 'l')
acf(env)
pacf(env)

# three: strong ar
AR_lg <- list(order = c(1, 0, 0), ar = 0.9)
AR1_lg <- arima.sim(n = t_length, model = AR_lg, sd = 0.1)
plot(AR1_lg, type = 'l')
plot(AR1_lg[burn:t_length], type = 'l')
acf(AR1_lg)
pacf(AR1_lg)

# four: ar and ma
# cps-like/ prey index that's more env driven
AR_ma <- list(order = c(1, 0, 1), ar = -0.1, ma = -0.1)
AR1_ma <- arima.sim(n = t_length, model = AR_ma, sd = 0.1)
plot(AR1_ma, type = 'l')
plot(AR1_ma[burn:t_length], type = 'l')
acf(AR1_ma)
pacf(AR1_ma)
