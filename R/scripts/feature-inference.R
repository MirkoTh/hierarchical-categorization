library(tidyverse)
library(cmdstanr)
library(rutils)
library(ggrepel)
library(grid)
library(gridExtra)
library(furrr)
library(loo)

utils_loc <- c("R/utils/plotting-utils.R", "R/utils/utils.R")
walk(utils_loc, source)


# Load Data and Preprocess Them -------------------------------------------

tbl_both <- readRDS(file = "data/infpro_task-cat_beh/tbl_both.RDS")
tbl_train <- readRDS(file = "data/infpro_task-cat_beh/tbl_train.RDS")
tbl_transfer <- readRDS(file = "data/infpro_task-cat_beh/tbl_transfer.RDS")

# params to revert back to untransformed space
mean_d1i <- mean(tbl_both$d1i)
sd_d1i <- sd(tbl_both$d1i)
mean_d2i <- mean(tbl_both$d2i)
sd_d2i <- sd(tbl_both$d2i)
l_pars_tf <- list(
  mean_d1i = mean_d1i, sd_d1i = sd_d1i, mean_d2i = mean_d2i, sd_d2i = sd_d2i
)



# example using category learning training stimuli ------------------------

p_id <- 101
file_loc_gcm <- str_c("data/infpro_task-cat_beh/models/gcm-model-", p_id, ".RDS")
file_loc_gaussian <- str_c(
  "data/infpro_task-cat_beh/models/gaussian-model-", p_id, ".RDS"
)
m_gcm <- readRDS(file_loc_gcm)
m_pt <- readRDS(file_loc_gaussian)
post_c <- m_gcm$draws(variables = "c", format = "df") %>% as_tibble()
post_pts <- m_pt$draws(variables = c("mu1", "mu2"), format = "df")


tbl_train %>% group_by(participant, d1i, d2i, category) %>%
  count()
# all category A and category B stimuli were seen 17 times during training
# we can therefore use all distinct category exemplars only once
# when computing similarities towards within-category exemplars as 17 cancels out

# all the exemplars observed during training that can be referred to in memory
l_tbl_exemplars <- tbl_train %>% 
  mutate(category = fct_recode(category, B = "C", C = "B")) %>%
  filter(participant == p_id) %>%
  group_by(category, d1i_z, d2i_z) %>%
  count() %>% select(-n) %>%
  split(.$category)

# just uses cat learn training data to get a set of inference cues
# has to be replaced once completion data are available
# as exemplars will contain the cue from one dimension
# plus a fine grid over plausible category values from the second dimension

l_tbl_lookup <- map2(
  l_tbl_exemplars, names(l_tbl_exemplars),
  ~ crossing(
    category = ..2,
    d1i_z = c(unique(..1$d1i_z), seq(min(..1$d1i_z), max(..1$d1i_z), by = .1)),
    d2i_z = c(unique(..1$d2i_z), seq(min(..1$d2i_z), max(..1$d2i_z), by = .1))
  )
)

# iterate over A & C categories aka target categories
ids_cat <- c("A", "C")
# iterate over cue dimensions
ids_dim <- c("d1i_z", "d2i_z")
tbl_cat_dim <- crossing(i_cat = ids_cat, i_dim = ids_dim)

l_closest <- pmap(
  tbl_cat_dim, 
  max_sim_responses, 
  post_c = post_c,
  l_tbl_lookup = l_tbl_lookup,
  l_tbl_exemplars = l_tbl_exemplars
)
tbl_closest <- reduce(l_closest, rbind)
tbl_closest$d1i_z * sd_d1i + mean_d1i



# implementation using empirical inference data ---------------------------


tbl_completion_prep <- read_csv(file = "data/infpro_task-cat_beh/sub-all_task-inf_beh-distances.csv")
cols_required <- c(
  "participant", "category", "rep", "cuedim", "cue_val",
  "respdim", "resp_i", "representation", "distance"
)
tbl_completion <- tbl_completion_prep[, cols_required] %>%
  group_by(participant, category, cuedim, respdim, cue_val, rep) %>%
  filter(representation == "prototype_phys") %>% ungroup() %>%
  select(-representation)

n_workers_available <- parallel::detectCores()
plan(multisession, workers = n_workers_available - 2)
safe_distances <- safely(distance_from_model_based_inference)

p_ids <- sort(unique(tbl_completion_prep$participant))


# GCM ---------------------------------------------------------------------


l_results_gcm <- future_map(
  p_ids, 
  safe_distances, 
  tbl_completion = tbl_completion, 
  tbl_train = tbl_train, 
  l_pars_tf = l_pars_tf,
  modeltype = "gcm",
  .progress = TRUE
)

# ok
l_gcm_results <- map(l_results_gcm, "result")
# not ok
map(l_results_gcm, "error") %>% reduce(c)
saveRDS(l_gcm_results, file = "data/infpro_task-cat_beh/inference-gcm-based.RDS")

tbl_gcm_results <- map(l_gcm_results, "tbl_empirical") %>% reduce(rbind)
tbl_gcm_results <- tbl_gcm_results %>%
  arrange(rep, participant, cue_val, cuedim)
tbl_completion <- tbl_completion %>%
  arrange(rep, participant, cue_val, cuedim)

tbl_gcm_results <- tbl_gcm_results %>% 
  left_join(
    tbl_completion %>% rename(distance_pt_phys = distance) %>% select(-resp_i),
    by = c("participant", "category", "rep", "cuedim", "cue_val", "respdim")
    )

tbl_gcm_results %>% 
  mutate(distance = abs(distance)) %>%
  rename(GCM = distance, `Physical PT` = distance_pt_phys) %>%
  pivot_longer(c(GCM, `Physical PT`)) %>% 
  ggplot(aes(value, group = name)) +
  geom_density(aes(color = name)) +
  theme_bw() +
  scale_color_brewer(palette = "Set1", name = "Model") +
  labs(x = "Distance", y = "Density")

tbl_gcm_results %>% 
  rename(GCM = distance) %>%
  ggplot(aes(GCM)) +
  geom_density() +
  theme_bw() +
  scale_color_brewer(palette = "Set1", name = "Model") +
  labs(x = "Distance", y = "Density")

ggplot(tbl_gcm_results %>% filter(category == "A"), aes(resp_i, group = cuedim)) +
  geom_density() +
  facet_grid(cuedim ~ cue_val)


tbl_lookup <- map(l_gcm_results, "tbl_lookup") %>% reduce(rbind)
write.csv(tbl_lookup, "data/infpro_task-cat_beh/gcm-inference-distances.csv")

# plot heat maps for some exemplary participants
# cue x1
ggplot(
  tbl_lookup %>% 
    filter(participant > 140 & participant %% 2 == 0) %>%
    mutate(category = str_c("Category = ", category))
  , aes(d1i, resp_i)
) + geom_tile(aes(fill = distance)) +
  facet_grid(participant ~ category) +
  scale_fill_gradient2(name = "Distance", low = "#FF6666", high = "#339966") +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  theme_bw() +
  labs(x = "Cue", y = "Response")


# cue x2
ggplot(
  tbl_lookup %>% 
    filter(participant > 140 & participant %% 2 == 1) %>%
    mutate(category = str_c("Category = ", category))
  , aes(resp_i, d2i)
  ) + geom_tile(aes(fill = distance)) +
  facet_grid(participant ~ category) +
  scale_fill_gradient2(name = "Distance", low = "#FF6666", high = "#339966") +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  theme_bw() +
  labs(x = "Response", y = "Cue")





l_posteriors <- map(p_ids, load_parameter_posteriors)
l_c <- map(l_posteriors, ~ .x$gcm$c)
tbl_c <- reduce(l_c, cbind) %>% as.data.frame()
colnames(tbl_c) <- p_ids
tbl_c$trial_id <- 1:nrow(tbl_c)
tbl_c %>% pivot_longer(-trial_id) %>%
  ggplot(aes(value)) +
  geom_histogram(color = "white", fill = "dodgerblue") +
  facet_wrap(~ name)





# Gaussian Prototype ------------------------------------------------------

l_results_gaussian <- future_map(
  p_ids, 
  safe_distances, 
  tbl_completion = tbl_completion, 
  tbl_train = tbl_train, 
  l_pars_tf = l_pars_tf,
  modeltype = "gaussian",
  .progress = TRUE
)

# ok
l_gaussian_results <- map(l_results_gaussian, "result")
# not ok
map(l_results_gaussian, "error") %>% reduce(c)
saveRDS(l_gaussian_results, file = "data/infpro_task-cat_beh/inference-gaussian-based.RDS")


# do a few checks
p_id <- 101
l_posteriors <- load_parameter_posteriors(p_id)
post_gaussian <- l_posteriors$gaussian
post_gaussian <- post_gaussian %>% rename(`mu2[1]` = `mu2[3]`, `mu2[3]` = `mu2[1]`) %>% select(-c(.chain, .iteration, .draw, `mu1[2]`, `mu2[2]`))
post_gaussian_mns <- apply(post_gaussian, 2, mean)
post_gaussian_mns * 2.72 + 5.5

tbl_gaussian_results <- map(l_gaussian_results, "tbl_empirical") %>% reduce(rbind)
tbl_gaussian_results <- tbl_gaussian_results %>%
  arrange(rep, participant, cue_val, cuedim)
tbl_completion <- tbl_completion %>%
  arrange(rep, participant, cue_val, cuedim)

tbl_gaussian_results <- tbl_gaussian_results %>% 
  left_join(
    tbl_completion %>% rename(distance_pt_phys = distance) %>% select(-resp_i),
    by = c("participant", "category", "rep", "cuedim", "cue_val", "respdim")
  )



tbl_gaussian_results %>% 
  mutate(distance = abs(distance)) %>%
  rename(Gaussian = distance, `Physical PT` = distance_pt_phys) %>%
  pivot_longer(c(Gaussian, `Physical PT`)) %>% 
  ggplot(aes(value, group = name)) +
  geom_density(aes(color = name)) +
  theme_bw() +
  scale_color_brewer(palette = "Set1", name = "Model") +
  labs(x = "Distance", y = "Density")

tbl_gaussian_results %>% 
  rename(Gaussian = distance) %>%
  ggplot(aes(Gaussian)) +
  geom_density() +
  theme_bw() +
  scale_color_brewer(palette = "Set1", name = "Model") +
  labs(x = "Distance", y = "Density")

ggplot(tbl_gaussian_results %>% filter(category == "A"), aes(resp_i, group = cuedim)) +
  geom_density() +
  facet_grid(cuedim ~ cue_val)


tbl_lookup <- map(l_gaussian_results, "tbl_lookup") %>% reduce(rbind)
write.csv(tbl_lookup, "data/infpro_task-cat_beh/gaussian-inference-distances.csv")
