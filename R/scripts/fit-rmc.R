library(tidyverse)
library(rutils)
library(ggrepel)
library(grid)
library(gridExtra)
library(wiqid)


home_grown <- c("R/utils/rmc-thm.R")
walk(home_grown, source)

tbl_train <- readRDS(file = "data/infpro_task-cat_beh/tbl_train.RDS")
tbl_transfer <- readRDS(file = "data/infpro_task-cat_beh/tbl_transfer.RDS")


tbl_train_101 <- tbl_train %>% filter(participant == 101) %>%
  rename(x1 = d1i_z, x2 = d2i_z)

# order: coupling, phi, salience_l, a_0, lambda_0
# optimal for 100 trials and 9 categories: 0.24191989 2.62343018 0.01000000 0.04399989 0.01000000
params_init <- c(.1, 3, .01, .02, .01)
bounds_lower <- c(.01, .01, .00001, .00001, .00001)
bounds_upper <- c(.5, 10, 3, 1, 1)

# plan(multisession = min((parallel::detectCores() - 1), length(l_tbl_used)))


wrap_optim <- function(tbl, n_categories, ...) {
  r <- optim(tbl_data = tbl, n_categories = n_categories, ...)
  beepr::beep()
  return(r)
}

# start 18:54
wrap_optim(
  tbl_train_101[, c("x1", "x2", "response_int")] %>% rename(category = response_int) %>% head(4) %>% select(x1, x2, category), 
  3, params_init, wrap_rmc, method = "L-BFGS-B", 
  lower = bounds_lower, upper = bounds_upper, control = list(factr = .001)
  )
