# simulate time series for ml forcasting model
# author(s): caitlin allen akselrud
# contact: caitlin.allen_akserud@noaa.gov
# date created: 13.07.2023
# version: 2.0 - standardized process error

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

all_100 <- NULL
all_100_sc <- NULL

for(iter in 1:100)
{
  # generate ts -------------------------------------------------------------

  # * scaled devs -----------------------------------------------------------

  dev <- rnorm(n = t_length, mean = 0, sd = 1) # raw devs
  scaled_dev <-  dev / sd(dev)
  # sd(scaled_dev)

  # hist(scaled_dev)
  # hist(dev)
  # * regime short ----------------------------------------------------------

  # regime short
  # one: periodic square
  # set.seed(1011)
  period = 8
  env1_square <- ifelse(((t %% period) < (0.5*period)),1,0) * rbinom(t_length, 1, 0.9)
  env1_square <- env1_square %>%
    as_tibble() %>%
    bind_cols(scaled_dev = scaled_dev) %>%
    mutate(
      mean = if_else(
        value == 0,3,5),
      regime = mean+scaled_dev)

  # sigma0 <- env1_square %>%
  #   dplyr::filter(value == 0) %>%
  #   summarize(sigma0 = sd(regime))
  # sigma1 <- env1_square %>%
  #   dplyr::filter(value == 1) %>%
  #   summarize(sigma1 = sd(regime))

  # sd ~1 for each regime, using scaled devs

  # * regime long -----------------------------------------------------------

  period = 40
  env2_square <- ifelse(((t %% period) < (0.5*period)),1,0) * rbinom(t_length, 1, 0.9)
  env2_square <- env2_square %>%
    as_tibble() %>%
    bind_cols(scaled_dev = scaled_dev) %>%
    mutate(
      mean = if_else(
        value == 0,3,5),
      regime = mean+scaled_dev)

  # sigma0 <- env2_square %>%
  #   dplyr::filter(value == 0) %>%
  #   summarize(sigma0 = sd(regime))
  # sigma1 <- env2_square %>%
  #   dplyr::filter(value == 1) %>%
  #   summarize(sigma1 = sd(regime))

  # sd ~1 for each regime, using scaled devs

  # * increasing cycle ------------------------------------------------------

  # signal
  # two: amplified signal
  # env <- sin((2*pi*t)/(2/(t)))
  env_base <- (2*(sin(pi*(t^2))))/(sd(sin(pi*(t^2))))
  env <- env_base + scaled_dev

  # sd(env_base)
  # sd(env) #ok b/c most sd comes from underlying process

  # proof of concept
  # check 10k sims
  # n_sim = 10000
  # as_dist <- rep(0,n_sim)
  # as_table <- matrix(nrow = n_sim, ncol = length(env))
  # for(iter in 1:n_sim)
  # {
  #   env_base <- (2*(sin(pi*(t^2))))/(sd(sin(pi*(t^2))))
  #   devs_as <- rnorm(n = length(env_base), mean = 0, sd = 1)
  #   scaled_dev_as <- devs_as / sd(devs_as)
  #   env <- env_base + scaled_dev_as
  #
  #   as_dist[iter] <- sd(env)
  #   as_table[iter,] <- env
  # }
  #
  # annual_sd <- as_table %>% as_tibble %>%
  #   summarise(across(where(is.numeric), sd))

  # range(annual_sd)
  # hist(as.vector(annual_sd) %>% unlist)
  # plot(env, type = 'l')

  # argument: because the point is non-stationary variance,
  #   we are showing that the process error for each year is ~1

  # * environmental trend ---------------------------------------------------

  # environmental trend
  # five: non-stationary (trend)

  phi <-  -0.5
  # corrected_sd <- sqrt(1-(phi^2))

  ts_ns <- list(order = c(1, 1, 0), ar = phi)
  ts_ns1 <- arima.sim(n = t_length, model = ts_ns, innov = scaled_dev) #'innov' setting assigns pre-fix devs

  # * Autocorrelated process ------------------------------------------------

  # autocorrelated process
  # three: strong ar (feedback)

  phi <- 0.9
  # corrected_sd <- sqrt(1-phi^2)
  AR_lg <- list(order = c(1, 0, 0), ar = phi)
  AR1_lg <- arima.sim(n = t_length, model = AR_lg, innov = scaled_dev)

  # * black swan ------------------------------------------------------------

  # black swan
  # four: ar and ma

  phi = -0.1
  theta = -0.1
  # variance <- (1 + theta^2 + 2 * phi * theta) / (1 - phi^2)
  # corrected_sd <- 1 / sqrt(variance)

  AR_ma <- list(order = c(1, 0, 1), ar = phi, ma = theta)
  AR1_ma <- arima.sim(n = t_length, model = AR_ma, innov = scaled_dev)

  # AR1_ma is the simulated prey time series
  # add black swan (bs) event(s)
  bs_regime <- rbinom(t_length, 1, 0.015) # random chance 1-2 times per 100 years is an extreme event
  bs_lower <- rnorm(n = t_length/2, mean = mean(AR1_ma) - sd(AR1_ma)*8, sd = sd(AR1_ma)/2) #lower extremes
  bs_upper <- rnorm(n = t_length/2, mean = mean(AR1_ma) + sd(AR1_ma)*8, sd = sd(AR1_ma)/2) #upper extremes
  bimodal <- c(bs_lower, bs_upper) #put the extreme cases together

  prey_bs <- bind_cols(AR1_ma = AR1_ma, bs_regime = bs_regime, bimodal = sample(bimodal)) %>%
    mutate(new_prey = if_else(bs_regime == 1, bimodal, AR1_ma)) #if bs event, randomly sample an extreme event (high or low)

  # * random walk -----------------------------------------------------------
  rw_t <- rw <- scaled_dev #process

  for(k in 2:100) {
    rw_t[k] <- rw_t[k-1] + rw[k]
  }

  # lagged rw
  rw_t_1 <- lag(rw_t, n = 1) #prev time step


  # * white noise -----------------------------------------------------------

  # wn <- scaled_dev
  # generate a new wn process, so that wn is different from process error in other ts
  # I think it makes sense that wn is the only process without the same process error...
  #      unless we want to explicitly test whether RF can pick up process error?
  wn <- rnorm(n=100, mean = 0, sd = 1)
  wn <- wn / sd(wn)

# env indices- full -------------------------------------------------------

env_data_full <- data.frame(
  regime_short = env1_square$regime,
  regime_long = env2_square$regime,
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
  env_data_sc$regime_short +
  env_data_sc$regime_long +
  env_data_sc$signal +
  env_data_sc$climate  +
  env_data_sc$pred +
  env_data_sc$prey +
  env_data_sc$random_walk_lag +
  env_data_sc$white_noise

base_simple_100 <- rowSums( cbind(env_data_full_sc$random_walk,
  env_data_full_sc$regime_short,
  env_data_full_sc$regime_long,
  env_data_full_sc$signal,
  env_data_full_sc$climate,
  env_data_full_sc$pred,
  env_data_full_sc$prey,
  env_data_full_sc$random_walk_lag,
  env_data_full_sc$white_noise), na.rm = TRUE)

base_complex <- env_data_sc$random_walk *
  env_data_sc$regime_short *
  env_data_sc$regime_long *
  env_data_sc$signal *
  env_data_sc$climate  *
  env_data_sc$pred *
  env_data_sc$prey *
  env_data_sc$random_walk_lag *
  env_data_sc$white_noise

base_complex_100 <- env_data_full_sc$random_walk *
  env_data_full_sc$regime_short *
  env_data_full_sc$regime_long *
  env_data_full_sc$signal *
  env_data_full_sc$climate  *
  env_data_full_sc$pred *
  env_data_full_sc$prey *
  env_data_full_sc$random_walk_lag *
  env_data_full_sc$white_noise

base_complex_100[1] <- env_data_full_sc$random_walk[1] *
  env_data_full_sc$regime_short[1] *
  env_data_full_sc$regime_long[1] *
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

# check results -----------------------------------------------------

all_100 <- read_csv(here::here("output", "complete_sim_input_data_unscaled_100_2025-08-31.csv"))

all_100_sc <- read_csv(here::here("output", "complete_sim_input_data_scaled_100_2025-08-31.csv"))

# make plots ---------------------------------------------------------
regime_p <- make_env_plots(env_index = as.vector(all_100$regime_short), name = "Regime short", burn = burn, t_length = t_length)
regime2_p <- make_env_plots(env_index = as.vector(all_100$regime_long), name = "Regime long", burn = burn, t_length = t_length)
signal_p <- make_env_plots(env_index = as.vector(all_100$signal), name = "Signal", burn = burn, t_length = t_length)
climate_p <- make_env_plots(env_index = as.vector(all_100$climate), name = "Climate", burn = burn, t_length = t_length)
pred_p <- make_env_plots(env_index = as.vector(all_100$pred), name = "Predator", burn = burn, t_length = t_length)
prey_p <- make_env_plots(env_index = as.vector(all_100$prey), name = "Prey", burn = burn, t_length = t_length)
rw_p <- make_env_plots(env_index = as.vector(all_100$random_walk), name = "Random walk", burn = burn, t_length = t_length)
wn_p <- make_env_plots(env_index = as.vector(all_100$white_noise), name = "White noise", burn = burn, t_length = t_length)

env_p <- regime_p | regime2_p | signal_p |climate_p |pred_p | prey_p |rw_p |wn_p #|rwlag_p
ggsave(plot = env_p, filename = paste0(iter, "_all_env.jpg"), path = here::here("figures"),
width = 21, height = 9)

regimes_p_30sc <- make_env_plots(env_index = as.vector(all_100_sc$regime_short), name = "Regime short", burn = burn, t_length = t_length)
regimel_p_30sc <- make_env_plots(env_index = as.vector(all_100_sc$regime_long), name = "Regime long", burn = burn, t_length = t_length)
signal_p_30sc <- make_env_plots(env_index = as.vector(all_100_sc$signal), name = "Signal", burn = burn, t_length = t_length)
climate_p_30sc <- make_env_plots(env_index = as.vector(all_100_sc$climate), name = "Climate", burn = burn, t_length = t_length)
pred_p_30sc <- make_env_plots(env_index = as.vector(all_100_sc$pred), name = "Predator", burn = burn, t_length = t_length)
prey_p_30sc <- make_env_plots(env_index = as.vector(all_100_sc$prey), name = "Prey", burn = burn, t_length = t_length)
rw_p_30sc <- make_env_plots(env_index = as.vector(all_100_sc$random_walk), name = "Random walk", burn = burn, t_length = t_length)
wn_p_30sc <- make_env_plots(env_index = as.vector(all_100_sc$white_noise), name = "White noise", burn = burn, t_length = t_length)


env_p_30sc <- regimes_p_30sc | regimel_p_30sc | signal_p_30sc |climate_p_30sc |pred_p_30sc | prey_p_30sc |rw_p_30sc  | wn_p_30sc #|rwlag_p_30sc
ggsave(plot = env_p_30sc, filename = paste0(iter, "_all_env_sc.jpg"), path = here::here("figures"),
       width = 21, height = 9)


# 100_sc

regimes_p_100sc <- make_env_plots(env_index = as.vector(all_100_sc$regime_short), name = "Regime short", burn = 1, t_length = t_length)
regimel_p_100sc <- make_env_plots(env_index = as.vector(all_100_sc$regime_long), name = "Regime long", burn = 1, t_length = t_length)
signal_p_100sc <- make_env_plots(env_index = as.vector(all_100_sc$signal), name = "Signal", burn = 1, t_length = t_length)
climate_p_100sc <- make_env_plots(env_index = as.vector(all_100_sc$climate), name = "Climate", burn = 1, t_length = t_length)
pred_p_100sc <- make_env_plots(env_index = as.vector(all_100_sc$pred), name = "Predator", burn = 1, t_length = t_length)
prey_p_100sc <- make_env_plots(env_index = as.vector(all_100_sc$prey), name = "Prey", burn = 1, t_length = t_length)
rw_p_100sc <- make_env_plots(env_index = as.vector(all_100_sc$random_walk), name = "Random walk", burn = 1, t_length = t_length)
wn_p_100sc <- make_env_plots(env_index = as.vector(all_100_sc$white_noise), name = "White noise", burn = 1, t_length = t_length)


env_p_100sc <- regimes_p_100sc | regimel_p_100sc | signal_p_100sc |climate_p_100sc |pred_p_100sc | prey_p_100sc |rw_p_100sc  |wn_p_100sc #|rwlag_p_100sc
ggsave(plot = env_p_100sc, filename = paste0(iter, "_all_env_100sc.jpg"), path = here::here("figures"),
       width = 21, height = 9)


# 100

regimes_p_100 <- make_env_plots(env_index = as.vector(all_100$regime_short), name = "Regime short", burn = 1, t_length = t_length)
regimel_p_100 <- make_env_plots(env_index = as.vector(all_100$regime_long), name = "Regime long", burn = 1, t_length = t_length)
signal_p_100 <- make_env_plots(env_index = as.vector(all_100$signal), name = "Signal", burn = 1, t_length = t_length)
climate_p_100 <- make_env_plots(env_index = as.vector(all_100$climate), name = "Climate", burn = 1, t_length = t_length)
pred_p_100 <- make_env_plots(env_index = as.vector(all_100$pred), name = "Predator", burn = 1, t_length = t_length)
prey_p_100 <- make_env_plots(env_index = as.vector(all_100$prey), name = "Prey", burn = 1, t_length = t_length)
rw_p_100 <- make_env_plots(env_index = as.vector(all_100$random_walk), name = "Random walk", burn = 1, t_length = t_length)
wn_p_100 <- make_env_plots(env_index = as.vector(all_100$white_noise), name = "White noise", burn = 1, t_length = t_length)

env_p_100 <- regimes_p_100 | regimel_p_100 | signal_p_100 |climate_p_100 |pred_p_100 | prey_p_100 |rw_p_100  | wn_p_100 #|rwlag_p_100
ggsave(plot = env_p_100, filename = paste0(iter, "_all_env_100.jpg"), path = here::here("figures"),
       width = 21, height = 9)

