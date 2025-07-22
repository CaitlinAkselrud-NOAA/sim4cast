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
library(cowplot)
library(patchwork)

# functions ---------------------------------------------------------------


# folders -----------------------------------------------------------------

write.dir <- here::here("output") # otherwise will be saved in working directory
dir.create(write.dir, showWarnings = FALSE)
# setwd(write.dir)

# settings and fixed values -----------------------------------------------

t_length <- 100
t <- 1:t_length
dat_length <- 30
burn <- 100-dat_length +1

# environmental -----------------------------------------------------------
make_env_plots <- function(env_index, name = "", burn = burn, t_length = t_length)
{
  index_p <- ggplot() +
    geom_line(aes(x = seq(from = 1, to = length(burn:t_length),
                          by = 1), y = env_index[burn:t_length]))+
    theme_classic() +
    labs(x = '',
         y = paste(name, "index"),
         title = name)

  acf_p <- ggAcf(env_index[burn:t_length], lag.max = length(env_index[burn:t_length])) +
    theme_classic() +
    labs(x = '',
         title = '')

  pacf_p <- ggPacf(env_index[burn:t_length], lag.max = length(env_index[burn:t_length])) +
    theme_classic() +
    labs(x = 'Years',
         title = '')

  p3 <- index_p + acf_p + pacf_p + plot_layout(ncol = 1)

  # ggsave(plot = p3, filename = paste0(name, "index.jpg"), path = here::here("figures"),
  #        width = 6, height = 12)
  return(p3)
}

all_100 <- NULL
all_100_sc <- NULL

for(iter in 1:100)
{
# one: periodic square
# set.seed(1011)
period = 8
env1_square <- ifelse(((t %% period) < (0.5*period)),1,0) * rbinom(t_length, 1, 0.9)
env1_square <- env1_square %>%
  as_tibble() %>%
  mutate(
    regime = if_else(
      value == 0,
      rnorm(n(), mean = 3, sd = 1), # n() gives the number of rows in the current group (here, all rows)
      rnorm(n(), mean = 5, sd = 1)))

# plot(env1_square$value, type = 'l')
# plot(env1_square$regime, type = 'l')
# plot(env1_square$regime[burn:t_length], type = 'l')


# par(mfrow=c(3,1), mar = c(1,4,2,2))
# plot(env1_square[burn:t_length], type = 'l', main = "Regime", ylab = "Regime index", xlab = "", xaxt="n")
# par(mar = c(1,4,1,2))
# acf(env1_square, main = "", lag.max = length(env1_square$regime[burn:t_length]), xlab = "", xaxt="n")
# par(mar = c(4,4,1,2))
# pacf(env1_square, main = "", lag.max = length(env1_square$regime[burn:t_length]), xlab = "Year")
#
# jpeg(here::here("figures","env1_square.jpg"), width = 600, height = 350)
# plot(env1_square[burn:t_length], type = 'l')
# dev.off()
regime_p <- make_env_plots(env_index = env1_square$regime, name = "Regime", burn = burn, t_length = t_length)

# two: amplified signal
# env <- sin((2*pi*t)/(2/(t)))
env <- sin(pi*(t^2))
# plot(env, type = 'l')
# plot(env[burn:t_length], type = 'l')
# acf(env)
# pacf(env)
# plot(diff(env))

# jpeg(here::here("figures","env_amp.jpg"), width = 600, height = 350)
# plot(env[burn:t_length], type = 'l')
# dev.off()
signal_p <- make_env_plots(env_index = env, name = "Signal", burn = burn, t_length = t_length)

# three: strong ar (feedback)
# set.seed(1234)
AR_lg <- list(order = c(1, 0, 0), ar = 0.9)
AR1_lg <- arima.sim(n = t_length, model = AR_lg, sd = 0.1)
# plot(AR1_lg, type = 'l')
# plot(AR1_lg[burn:t_length], type = 'l')
# acf(AR1_lg)
# pacf(AR1_lg)

# jpeg(here::here("figures","AR_lg.jpg"), width = 600, height = 350)
# plot(AR1_lg[burn:t_length], type = 'l')
# dev.off()
pred_p <- make_env_plots(env_index = AR1_lg, name = "Predator", burn = burn, t_length = t_length)


# four: ar and ma
# cps-like/ prey index that's more env driven
# set.seed(9876)
AR_ma <- list(order = c(1, 0, 1), ar = -0.1, ma = -0.1)
AR1_ma <- arima.sim(n = t_length, model = AR_ma, sd = 0.1)
# plot(AR1_ma, type = 'l')
# plot(AR1_ma[burn:t_length], type = 'l')
# acf(AR1_ma)
# pacf(AR1_ma)

# AR1_ma is the simulated prey time series
# add black swan (bs) event(s)
bs_regime <- rbinom(t_length, 1, 0.015)
bs_lower <- rnorm(n = t_length/2, mean = mean(AR1_ma) - sd(AR1_ma)*8, sd = sd(AR1_ma)/2)
bs_upper <- rnorm(n = t_length/2, mean = mean(AR1_ma) + sd(AR1_ma)*8, sd = sd(AR1_ma)/2)
bimodal <- c(bs_lower, bs_upper)
# par(mfrow=c(2,1))
# hist(bimodal)
# hist(AR1_ma)
# hist(c(AR1_ma, bimodal))

prey_bs <- bind_cols(AR1_ma = AR1_ma, bs_regime = bs_regime, bimodal = sample(bimodal)) %>%
  mutate(new_prey = if_else(bs_regime == 1, bimodal, AR1_ma))

# plot(AR1_ma, type = 'l', ylim = c(-0.6, 0.6))
# plot(prey_bs$new_prey, type = 'l')
# plot(AR1_ma[burn:t_length], type = 'l', ylim = c(-0.6, 0.6))
# plot(prey_bs$new_prey[burn:t_length], type = 'l')

# devs <- rnorm(t_length, mean = c(mu1, mu2)[bs_regime],  sd = c(sd1, sd2)[bs_regime])

# jpeg(here::here("figures","AR_ma.jpg"), width = 600, height = 350)
# plot(AR1_ma[burn:t_length], type = 'l')
# dev.off()
# prey_p <- make_env_plots(env_index = AR1_ma, name = "Prey", burn = burn, t_length = t_length)
# prey_p <- make_env_plots(env_index = x, name = "Prey", burn = burn, t_length = t_length)
prey_p <- make_env_plots(env_index = prey_bs$new_prey, name = "Prey", burn = burn, t_length = t_length)

# five: non-stationary (trend)
# set.seed(5678)
ts_ns <- list(order = c(1, 1, 0), ar = -0.5)
ts_ns1 <- arima.sim(n = t_length, model = ts_ns, sd = 0.1) *-1
# plot(ts_ns1, type = 'l')
# plot(ts_ns1[burn:t_length], type = 'l')
# acf(ts_ns1)
# pacf(ts_ns1)

# jpeg(here::here("figures","ts_ns1.jpg"), width = 600, height = 350)
# plot(ts_ns1[burn:t_length], type = 'l')
# dev.off()
climate_p <- make_env_plots(env_index = ts_ns1, name = "Climate", burn = burn, t_length = t_length)


# six: make random walk
rw_t <- rw <- rnorm(n = 100) #process

# rw_targ <- rw_t +rw_t_1

for(t in 2:100) {
  rw_t[t] <- rw_t[t-1] + rw[t]
}

# seven: create lagged rw
rw_t_1 <- lag(rw_t, n = 1) #prev time step

# get correct time for rw's
burn_to <- t_length - burn + 1
random_walk <- tail(rw_t, n = burn_to)
random_walk_lag <- tail(rw_t_1, n = burn_to)
# random_walk_target <- tail(rw_targ, n = burn_to)

rw_p <- make_env_plots(env_index = rw_t, name = "Random walk", burn = burn, t_length = t_length)
rwlag_p <- make_env_plots(env_index = rw_t_1, name = "Random walk lag", burn = burn, t_length = t_length)

# white noise
wn <- rnorm(n=100)
wn_p <- make_env_plots(env_index = wn, name = "White noise", burn = burn, t_length = t_length)

# env indices- full -------------------------------------------------------

env_data_full <- data.frame(
  regime = env1_square$regime,
  signal = env,
  climate = tail(as.vector(ts_ns1), n = length(env1_square$regime)),
  pred = AR1_lg,
  prey = prey_bs$new_prey,
  random_walk = tail(rw_t, n = length(env1_square$regime)),
  random_walk_lag = tail(rw_t_1, n = length(env1_square$regime)),
  white_noise = wn) #%>%
# scale(center = TRUE, scale = TRUE)

m_env_full <- t(env_data_full)

colnames(m_env_full) <- seq(from = lubridate::year(Sys.Date())-100+1,
                       to = lubridate::year(Sys.Date()),
                       by = 1)

n_full <- nrow(m_env_full) - 1

env_data_full <- as_tibble(env_data_full)

env_data_full_sc <- as_tibble(env_data_full) %>%
  scale(center = TRUE, scale = TRUE) %>% as_tibble


# env indices - last 30 ---------------------------------------------------

env_data <- env_data_full %>% slice_tail(n = dat_length)

env_data_sc <- env_data_full_sc %>% slice_tail(n = dat_length)

env_p <- regime_p | signal_p |climate_p |pred_p | prey_p |rw_p |rwlag_p |wn_p
ggsave(plot = env_p, filename = paste0(iter, "_all_env.jpg"), path = here::here("figures"),
       width = 21, height = 9)


# make more plots ---------------------------------------------------------

regime_p_30sc <- make_env_plots(env_index = as.vector(env_data_full_sc$regime), name = "Regime", burn = burn, t_length = t_length)
signal_p_30sc <- make_env_plots(env_index = as.vector(env_data_full_sc$signal), name = "Signal", burn = burn, t_length = t_length)
climate_p_30sc <- make_env_plots(env_index = as.vector(env_data_full_sc$climate), name = "Climate", burn = burn, t_length = t_length)
pred_p_30sc <- make_env_plots(env_index = as.vector(env_data_full_sc$pred), name = "Predator", burn = burn, t_length = t_length)
prey_p_30sc <- make_env_plots(env_index = as.vector(env_data_full_sc$prey), name = "Prey", burn = burn, t_length = t_length)
rw_p_30sc <- make_env_plots(env_index = as.vector(env_data_full_sc$random_walk), name = "Random walk", burn = burn, t_length = t_length)
rwlag_p_30sc <- make_env_plots(env_index = as.vector(env_data_full_sc$random_walk_lag), name = "Random walk lag", burn = burn, t_length = t_length)
wn_p_30sc <- make_env_plots(env_index = as.vector(env_data_full_sc$white_noise), name = "White noise", burn = burn, t_length = t_length)


env_p_30sc <- regime_p_30sc | signal_p_30sc |climate_p_30sc |pred_p_30sc | prey_p_30sc |rw_p_30sc |rwlag_p_30sc | wn_p_30sc
ggsave(plot = env_p_30sc, filename = paste0(iter, "_all_env_sc.jpg"), path = here::here("figures"),
       width = 21, height = 9)


# 100_sc

regime_p_100sc <- make_env_plots(env_index = as.vector(env_data_full_sc$regime), name = "Regime", burn = 1, t_length = t_length)
signal_p_100sc <- make_env_plots(env_index = as.vector(env_data_full_sc$signal), name = "Signal", burn = 1, t_length = t_length)
climate_p_100sc <- make_env_plots(env_index = as.vector(env_data_full_sc$climate), name = "Climate", burn = 1, t_length = t_length)
pred_p_100sc <- make_env_plots(env_index = as.vector(env_data_full_sc$pred), name = "Predator", burn = 1, t_length = t_length)
prey_p_100sc <- make_env_plots(env_index = as.vector(env_data_full_sc$prey), name = "Prey", burn = 1, t_length = t_length)
rw_p_100sc <- make_env_plots(env_index = as.vector(env_data_full_sc$random_walk), name = "Random walk", burn = 1, t_length = t_length)
rwlag_p_100sc <- make_env_plots(env_index = as.vector(env_data_full_sc$random_walk_lag), name = "Random walk lag", burn = 1, t_length = t_length)
wn_p_100sc <- make_env_plots(env_index = as.vector(env_data_full_sc$white_noise), name = "White noise", burn = 1, t_length = t_length)


env_p_100sc <- regime_p_100sc | signal_p_100sc |climate_p_100sc |pred_p_100sc | prey_p_100sc |rw_p_100sc |rwlag_p_100sc |wn_p_100sc
ggsave(plot = env_p_100sc, filename = paste0(iter, "_all_env_100sc.jpg"), path = here::here("figures"),
       width = 21, height = 9)


# 100

regime_p_100 <- make_env_plots(env_index = as.vector(env_data_full$regime), name = "Regime", burn = 1, t_length = t_length)
signal_p_100 <- make_env_plots(env_index = as.vector(env_data_full$signal), name = "Signal", burn = 1, t_length = t_length)
climate_p_100 <- make_env_plots(env_index = as.vector(env_data_full$climate), name = "Climate", burn = 1, t_length = t_length)
pred_p_100 <- make_env_plots(env_index = as.vector(env_data_full$pred), name = "Predator", burn = 1, t_length = t_length)
prey_p_100 <- make_env_plots(env_index = as.vector(env_data_full$prey), name = "Prey", burn = 1, t_length = t_length)
rw_p_100 <- make_env_plots(env_index = as.vector(env_data_full$random_walk), name = "Random walk", burn = 1, t_length = t_length)
rwlag_p_100 <- make_env_plots(env_index = as.vector(env_data_full$random_walk_lag), name = "Random walk", burn = 1, t_length = t_length)
wn_p_100 <- make_env_plots(env_index = as.vector(env_data_full$white_noise), name = "White noise", burn = 1, t_length = t_length)

env_p_100 <- regime_p_100 | signal_p_100 |climate_p_100 |pred_p_100 | prey_p_100 |rw_p_100 |rwlag_p_100 | wn_p_100
ggsave(plot = env_p_100, filename = paste0(iter, "_all_env_100.jpg"), path = here::here("figures"),
       width = 21, height = 9)

# sim notes ---------------------------------------------------------------

# simple: regime + signal + climate + pred + prey
# surplus prod model w/ deviates that are a fxn of env indices
#  + add on fishing
#  + age-struct pop- cohorts

# correlations in space in time
# bad sampling
# finite sample size

# variance-bias trade-off (aggregation level vs number of sample sizes)
# trade-off between sample size and variance

# neural nets- tensorflow, complex, same as trees, harder to fit
# is there only so much information, and you only get so much more from more complex models;
# # OR do we need more coaching on model config?

# start with simple sim and simple model (RF)
# how weird and specific do we have to make things to get it right?

# underperformance vs overperformance errors
# under = lags incorrect, eg
# 0ver = r^2 = .99 without thinking about spatio-temporal autocorr
# the less well specified your model, the porer performance; if you get it wrong and have over perf, bigger prob


# make target data ---------------------------------------------------------

# generate all combos of target data
# target = rw + regime + signal + drift + AR + ARMA + lag
base_simple <- env_data_sc$random_walk +
  env_data_sc$regime +
  env_data_sc$signal +
  env_data_sc$climate  +
  env_data_sc$pred +
  env_data_sc$prey +
  env_data_sc$random_walk_lag +
  env_data_sc$white_noise

base_simple_100 <- rowSums( cbind(env_data_full_sc$random_walk,
  env_data_full_sc$regime,
  env_data_full_sc$signal,
  env_data_full_sc$climate,
  env_data_full_sc$pred,
  env_data_full_sc$prey,
  env_data_full_sc$random_walk_lag,
  env_data_full_sc$white_noise), na.rm = TRUE)

base_complex <- env_data_sc$random_walk *
  env_data_sc$regime *
  env_data_sc$signal *
  env_data_sc$climate  *
  env_data_sc$pred *
  env_data_sc$prey *
  env_data_sc$random_walk_lag *
  env_data_sc$white_noise

base_complex_100 <- env_data_full_sc$random_walk *
  env_data_full_sc$regime *
  env_data_full_sc$signal *
  env_data_full_sc$climate  *
  env_data_full_sc$pred *
  env_data_full_sc$prey *
  env_data_full_sc$random_walk_lag *
  env_data_full_sc$white_noise

base_complex_100[1] <- env_data_full_sc$random_walk[1] *
  env_data_full_sc$regime[1] *
  env_data_full_sc$signal[1] *
  env_data_full_sc$climate[1]  *
  env_data_full_sc$pred[1] *
  env_data_full_sc$prey[1] *
  env_data_full_sc$white_noise[1]


base_complex_log <- log(as.vector(base_complex) + (1 - min(as.vector(base_complex), na.rm = T)))
base_complex_log_100 <- log(as.vector(base_complex_100) + (1 - min(as.vector(base_complex_100), na.rm = T)))

random_walk_target <- env_data_sc$random_walk + env_data_sc$random_walk_lag
random_walk_target_100 <- rowSums( cbind(env_data_full_sc$random_walk, env_data_full_sc$random_walk_lag), na.rm = T)

# input_data_30 <- env_data %>%
#   bind_cols(base_simple = as.vector(base_simple),
#             base_complex = as.vector(base_complex),
#             base_complex_log = base_complex_log,
#             random_walk_target = random_walk_target) %>%
#   # scale(center = TRUE, scale = TRUE) %>%
#   bind_cols(sim_year = seq(from = lubridate::year(Sys.Date())-dat_length+1,
#                            to = lubridate::year(Sys.Date()),
#                            by = 1))

input_data_sc_30 <- env_data_sc %>%
  bind_cols(base_simple = as.vector(base_simple),
            base_complex = as.vector(base_complex),
            base_complex_log = base_complex_log,
            random_walk_target = random_walk_target) %>%
  # scale(center = TRUE, scale = TRUE) %>%
  bind_cols(sim_year = seq(from = lubridate::year(Sys.Date())-dat_length+1,
                           to = lubridate::year(Sys.Date()),
                           by = 1))

input_data_sc_100 <- env_data_full_sc %>%
  bind_cols(base_simple = as.vector(base_simple_100),
            base_complex = as.vector(base_complex_100),
            base_complex_log = base_complex_log_100,
            random_walk_target_100 = random_walk_target_100) %>%
  bind_cols(sim_year = seq(from = lubridate::year(Sys.Date())-100+1,
                           to = lubridate::year(Sys.Date()),
                           by = 1))

input_data_30 <- env_data %>%
  bind_cols(base_simple = as.vector(base_simple),
            base_complex = as.vector(base_complex),
            base_complex_log = base_complex_log,
            random_walk_target = random_walk_target) %>%
  # scale(center = TRUE, scale = TRUE) %>%
  bind_cols(sim_year = seq(from = lubridate::year(Sys.Date())-dat_length+1,
                           to = lubridate::year(Sys.Date()),
                           by = 1))

input_data_100 <- env_data_full %>%
  bind_cols(base_simple = as.vector(base_simple_100),
            base_complex = as.vector(base_complex_100),
            base_complex_log = base_complex_log_100,
            random_walk_target_100 = random_walk_target_100) %>%
  bind_cols(sim_year = seq(from = lubridate::year(Sys.Date())-100+1,
                           to = lubridate::year(Sys.Date()),
                           by = 1))



write_csv(input_data_sc_30, file = here::here("output", paste0(iter, "_sim_input_data_sc_30_", Sys.Date(),".csv")))
write_csv(input_data_sc_100, file = here::here("output", paste0(iter, "_sim_input_data_sc_100_", Sys.Date(),".csv")))
write_csv(input_data_30, file = here::here("output", paste0(iter, "_sim_input_data_unscaled_30_", Sys.Date(),".csv"))) #all targets are sc, but env indices are unscaled here
write_csv(input_data_100, file = here::here("output", paste0(iter, "_sim_input_data_unscaled_100_", Sys.Date(),".csv"))) #all targets are sc, but env indices are unscaled here

all_100 <- bind_rows(all_100, input_data_100)
all_100_sc <- bind_rows(all_100_sc, input_data_sc_100)


}

write_csv(all_100, file = here::here("output", paste0("complete_sim_input_data_unscaled_100_", Sys.Date(),".csv"))) #all targets are sc, but env indices are unscaled here
write_csv(all_100_sc, file = here::here("output", paste0("complete_sim_input_data_scaled_100_", Sys.Date(),".csv"))) #all targets are sc, but env indices are unscaled here

# also create example with multicollinear factors
