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
# x_t+1  = B*x_t + u + C*env_t + w_t
# B = identity matrix
# (m = number of pop trajectories (num rows)) = 1 gor now
# u = scalar on B #if b is not one, you de-mean the data
# B = 1 is random walk setup
# u  = 0.01 = 1% per year pop inc
# C = m rows and # env cols
# C = matrix(0.1, -0.1, 0.2) #effects of environment

# sim data
# 1) use MARSSsimulate but you're locked intp additive effects
# 2) brute force simulation eg

# x_t+1  = B*x_t + u + C*env_t + w_t
# y_t = x_t + v_t
# See and Holms adjusted number of rows in y

r = .1 #obsv variance
q = .1 #process
nsim = 100
m = 1 #number of populations, so x has 1 row

# set initital state
r_err <- rnorm
sim_x[,1] <- rnorm(nsim, 0, 0.1) #0 in log scale; just some variance
sim_y[,1] <- sim_x[,1] + r_err


sim_dat = matrix(NA, nsim, t_length)
for(t in 2:t_length)
{
  r_err <- rnorm(nsim, 0, sd = sqrt(r)) #r is variance, no sd, so put in sd
  q_err <- rnorm(nsim, 0, sd = sqrt(q))
  sim_x[,t] <- B * sim_data[,t-1] + U + q_err
  # this only works if m = 1 (pop data = 1 row)
  sim_y[,t] <- B * sim_data[,t-1] + U + r_err
  # can play with C to get diff env effects
}

# observation error -------------------------------------------------------


# ar coefficient  ---------------------------------------------------------


# notes -------------------------------------------------------------------



