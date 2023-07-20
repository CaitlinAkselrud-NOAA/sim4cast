# simulate time series for ml forcasting model
# author(s): caitlin allen akselrud
# contact: caitlin.allen_akserud@noaa.gov
# date created: 13.07.2023
# version: 1.0

# devtools::install_github("timjmiller/wham", dependencies=TRUE)


# libraries ---------------------------------------------------------------
library(tidyverse)
library(wham)
library(here)

# data --------------------------------------------------------------------


# environmental covariates ------------------------------------------------


# fit wham model ----------------------------------------------------------
# base = Ex 1 model 3 (m3)
# use example 2 to set up models with diff cov effects on recr
# try changing sigmaR and sigma_olderages for survival effects

# simulate wham time series -----------------------------------------------


# check variance and autocorrelation --------------------------------------
# remember: you want combos of high/low AC and high/low variance

# save time series output -------------------------------------------------


