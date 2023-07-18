# simulate time series for ml forcasting model
# author(s): caitlin allen akselrud
# contact: caitlin.allen_akserud@noaa.gov
# date created: 13.07.2023
# version: 1.0

devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(tidyverse)
library(wham)
library(here)

# base = Ex 1 model 3 (m3)

# add environmental covariates

# try changing sigmaR and sigma_olderages for survival effects

# remember: you want combos of high/low AC and high/low variance
