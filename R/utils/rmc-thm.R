# original code from Badham, Sanborn, & Maylor (2017)
# https://osf.io/m7tck/files/
# date of dump: 2021-09-08
# implementation of the Rational Model of Categorization (Anderson, 1991) for the 
# Shepard, Hovland, and Jenkins (1961) task

predict_rmc_continuous <- function(
  #' predict from Anderson's RMC using continuous features
  #' 
  #' @description predict from Anderson's rational model of categorization (1991)
  #' using the original inference algorithm for the posterior (aka local MAP)
  #' @param stimuli \code{matrix} containing the feature values of one stimulus per row
  #' @param features_cat \code{vector} with names of categorical features as strings
  #' @param features_cont \code{vector} with names of continuous features as strings
  #' @param n_values_cat number of levels for the categorical features
  #' @param n_categories number of categories
  #' @param feedback \code{integer vector} containing category labels as integers
  #' @param params \code{list} containing model parameters 
  #' (salience_f, salience_l, a_0, lambda_0, coupling, phi)
  #' @param previous_learning \code{list} containing cluster_counts vector and
  #' feature_counts array for previously learned categories
  #' @param assignments \code{vector} of integers stating the category
  #' assignments. Defaults to NULL such that inferred categories are saved
  #' @param print_posterior {logical} stating whether the posterior should
  #' be printed while code is running
  #' @return the predicted category probabilities and category assignments
  #' as a \code{list}
  #' 
  stimuli,
  features_cat,
  features_cont,
  n_values_cat,
  n_categories,
  feedback,
  params,
  previous_learning = NULL,
  assignments = NULL, 
  print_posterior = FALSE
) {
  # unpack parameters
  env <- rlang::current_env()
  list2env(params, env)
  
  # stimuli information
  n_stimuli <- nrow(stimuli)
  n_features_cat <- length(features_cat) + 1 # including label feature
  n_features_cont <- length(features_cont)
  assignments <- rep(0, n_stimuli) # assignment of stimuli to clusters
  max_clusters <- nrow(stimuli)
  
  salience <- c(rep(salience_f, times = (n_features_cat - 1)), salience_l)
  
  # cluster_counts: counts of stimuli in each cluster
  # feature_counts: cluster x feature x value array, includes prior pseudocounts
  # n_cluster: # position of the lowest currently empty cluster
  if (!is.null(previous_learning)) {
    cluster_counts <- previous_learning[["cluster_counts"]]
    feature_counts <- previous_learning[["feature_counts"]]
    n_clusters <- length(cluster_counts[cluster_counts != 0]) + 1
  } else {
    cluster_counts <- rep(0, max_clusters)
    feature_counts <- array(
      rep(salience, each = max_clusters),
      dim = c(max_clusters, n_features_cat, max(n_values_cat, n_categories))
    )
    n_clusters <- 1
  }
  
  out <- c()
  out$cat_probs <- matrix(nrow = nrow(stimuli), ncol = n_categories)
  
  for(i in 1:n_stimuli){
    i_prev <- i
    stimuli_concat <- stimuli
    assignments_concat <- assignments
    if (!is.null(previous_learning)) {
      stimuli_concat <- rbind(previous_learning[["stimuli"]], stimuli[i, ])
      assignments_concat <- c(previous_learning[["assignments"]], assignments[i])
      i_prev <- nrow(previous_learning[["stimuli"]]) + 1
    }
    
    # calculate prior, likelihoods, and posterior
    # clusters that already have stimuli
    log_prior <- log(coupling) + log(cluster_counts[1:n_clusters])
    # new cluster prior (in first empty position)
    log_prior[n_clusters] <- log(1 - coupling)
    log_prior <- log_prior - log(
      1 - coupling * (1 - sum(cluster_counts[1:n_clusters]))
    )
    
    # likelihood for categorical features
    pdfs_cat_log <- pdf_cat_log(
      stimuli_concat, i_prev, features_cat, salience, feature_counts, cluster_counts,
      n_values_cat, n_categories, n_clusters, n_features_cat
    )
    log_likelihood_prep <- pdfs_cat_log
    
    # likelihood for continuous features
    if (n_features_cont > 0) {
      pdfs_cont_log <- pdf_cont_log(
        stimuli_concat, assignments_concat, i_prev, features_cont, mu_0, lambda_0, a_0, sigma_sq_0
      )
      pdfs_cont_log <- array(
        rep(pdfs_cont_log, n_categories),
        dim = c(n_features_cont,n_clusters,  n_categories)
      ) %>% aperm(c(2, 1, 3))
      log_likelihood_prep <- abind::abind(
        pdfs_cat_log,
        pdfs_cont_log,
        along = 2
      )
    }
    
    log_likelihood <- apply(
      log_likelihood_prep, MARGIN = c(1, 3), FUN = sum
    )
    log_posterior <- matrix(
      log_prior, nrow = n_clusters, ncol = n_categories
    ) + log_likelihood
    
    if(print_posterior & !is.null(feedback)){
      print(exp(log_posterior[, feedback[i]]))
    }
    
    # compute prediction
    label_posterior <- colSums(exp(log_posterior - max(log_posterior)))
    out$cat_probs[i, ] <- label_posterior^phi / sum(label_posterior^phi)
    
    if (!is.null(feedback)) {
      # update cluster assignment and count variables using Anderson update rule
      assignments[i] <- which.max(log_posterior[, feedback[i]])
      cluster_counts[assignments[i]] <- cluster_counts[assignments[i]] + 1
      feature_index_update <- cbind(
        rep(assignments[i], times=n_features_cat), 
        1:n_features_cat, 
        c(stimuli_concat[i, features_cat] %>% as_vector(), feedback[i])
      )
      feature_counts[feature_index_update] <- (
        feature_counts[feature_index_update] + 1
      )
      n_clusters <- max(assignments) + 1
    }
    
  }
  out$stimuli <- stimuli
  out$assignments <- assignments
  out$feature_counts <- feature_counts
  out$cluster_counts <- cluster_counts
  return(out)
}

pdf_cont_log <- function(
  #' pds of values from categorical feature dimensions
  #' 
  #' @param stimuli \code{tibble} with continuous features as columns
  #' @param assignments \code{integer vector} with previous cluster assignments 
  #' @param i running index of trials in experiment
  #' @param features_cont \code{character vector} with continuous feature names
  #' @param mu_0 mean of prior mean normal distribution
  #' @param lambda_0 confidence in prior mean
  #' @param a_0 df of prior chisq variance distribution
  #' @param sigma_sq_0 expected prior variance
  #' @return 3D \code{array} with the log probabilities of item combinations per category
  #'
  stimuli, assignments, i, features_cont, 
  mu_0, lambda_0, a_0, sigma_sq_0
) {
  tbl_part_1 <- tibble(stimuli[1:(i-1), ], assignments = assignments[1:(i-1)])
  tbl_part_2 <- crossing(
    stimuli[which(assignments[1:i] == 0), ],
    assignments = unique(assignments)
  )
  tbl_part_2$assignments[tbl_part_2$assignments == 0] <- max(assignments) + 1
  tbl_relevant <- rbind(tbl_part_1, tbl_part_2)
  n_clusters <- length(unique(tbl_relevant$assignments))
  tbl_summary <- tbl_relevant %>%
    select(all_of(features_cont), assignments) %>%
    group_by(assignments) %>%
    summarize_if(
      is.numeric, list(length, mean, var)
    )
  tbl_summary[is.na(tbl_summary)] <- 0
  names(tbl_summary) <- str_replace(names(tbl_summary), "fn1$", "n") %>% 
    str_replace("fn2$", "mean") %>%
    str_replace("fn3$", "var")
  tbl_summary$lambda_i <- tbl_summary[, 2] %>% as_vector() + lambda_0
  tbl_summary$a_i <- tbl_summary[, 2] %>% as_vector() + a_0
  col_longer <- names(tbl_summary)[str_detect(names(tbl_summary), "n$|mean$|var$")]
  tbl_summary <- tbl_summary %>%
    pivot_longer(all_of(col_longer)) %>%
    separate(name, c("variable", "aggregation"), sep = "_") %>%
    pivot_wider(names_from = variable, values_from = value)
  agg_stats <- tbl_summary %>%
    split(~ aggregation) %>%
    map(~ select(.x, all_of(features_cont)))
  mu_i <- (
    (lambda_0 * mu_0 + agg_stats[["n"]] * agg_stats[["mean"]]) /
      (lambda_0 + agg_stats[["n"]])
  )
  above_1 <- a_0 * sigma_sq_0
  above_2 <- (agg_stats[["n"]] - 1) * agg_stats[["var"]]
  above_3 <- (
    (lambda_0 * agg_stats[["n"]]) / 
      (lambda_0 + agg_stats[["n"]])
  ) * (mu_0 - agg_stats[["mean"]]) ^ 2
  below <- a_0 + agg_stats[["n"]]
  sigma_sq_i <- (above_1 + above_2 + above_3) / below
  update_i <- tbl_summary %>% 
    group_by(assignments) %>%
    summarize(
      lambda_i = min(lambda_i),
      a_i = min(a_i)
    ) %>% ungroup() %>% select(-assignments)
  scale <- sqrt(sigma_sq_i) * sqrt(1 + 1 / update_i$lambda_i)
  
  v_x <- rep(stimuli[i, features_cont] %>% as_vector(), n_clusters)
  v_mu_i <- as.vector(t(mu_i))
  v_scale <- as.vector(t(scale))
  v_ai <- rep(update_i$a_i, each = length(features_cont))
  
  v_out <- suppressWarnings(
    pmap_dbl(
      list(v_x, v_mu_i, v_scale, v_ai), dt2)
  ) %>% log()
  
  return(v_out)
}


pdf_cat_log <- function(
  #' pds of values from continuous feature dimensions
  #' 
  #' @param stimuli \code{tibble} with continuous features as columns
  #' @param i running index of trials in experiment
  #' @param features_cat \code{character vector} with categorical feature names
  #' @param salience \code {integer vector} categorical and label priors
  #' @param feature_counts feature counts \code{matrix}
  #' @param cluster_counts cluster count \code{integer vector}
  #' @param n_values_cat \code{integer} levels of the categorical features (assuming all the same)
  #' @param n_categories \code{integer} stating the nr of categories
  #' @param n_clusters \code{integer} stating the max nr of clusters
  #' @param n_features_cat \code{integer} stating the number of categorical features
  #' @return 3D \code{array} with the log probabilities of item combinations per category
  #'
  stimuli, i, features_cat, salience, feature_counts, cluster_counts,
  n_values_cat, n_categories, n_clusters, n_features_cat
) {
  # add labels to end of stimuli
  possible_stimuli <- cbind(
    matrix(rep(stimuli[i, features_cat], n_categories), 
           nrow = n_categories, byrow = TRUE),
    seq(1, (n_categories), by = 1)
  )
  
  # cols: cluster nr, feature nr (including cat label), values
  feature_index <- cbind(
    rep(1:n_clusters, times = n_features_cat * n_categories),
    rep(1:n_features_cat, times = n_categories, each = n_clusters),
    rep(t(possible_stimuli) %>% as_vector(), each = n_clusters)
  )
  this_num <- array(
    feature_counts[feature_index], 
    dim = c(n_clusters, n_features_cat, n_categories)
  )
  # what goes into den(ominator) to achieve a uniform prior?
  multiply_salience <- c(
    rep(n_values_cat, (n_features_cat - 1)), 
    n_categories
  )
  this_den <- array(
    outer(cluster_counts[1:n_clusters], multiply_salience * salience, FUN = "+"),
    dim=c(n_clusters, n_features_cat, n_categories)
  )
  return(log(this_num) - log(this_den))
}



predict_rmc_n <- function(
  #' predict n randomly sampled stimuli from stimulus space using Anderson's RMC
  #' 
  #' @param tbl \code{tibble} with feature values and category labels as columns
  #' @param n number of trials to sample with replacement and predict
  #' @param features_cat \code{character vector} with names of categorical features
  #' @param features_cont \code{character vector} with names of continuous features
  #' @param n_values_cat \code{integer} levels of the cateorical features (assuming all the same)
  #' @param max_clusters \code{integer} max nr of clusters to use in the code
  #' @param coupling \code{integer} coupling probability c
  #' @param salience_f \code{integer} feature-salience prior
  #' @param salience_l \code{integer} label-salience prior
  #' @param phi \code{numeric} scaling parameter for response probabilities
  #' @param assignments \code{vector} of integers stating the category
  #' assignments. Defaults to NULL such that inferred categories are saved
  #' @param print_posterior {logical} stating whether the posterior should
  #' be printed while code is running
  #' @return the predicted category probabilities and category assignments
  #' as a \code{list}
  #'
  tbl,
  n,
  features_cat,
  features_cont,
  n_values_cat,
  max_clusters,
  params,
  assignments = NULL
) {
  
  idxs <- 1:nrow(tbl)
  idx_shuffle <- sample(idxs, n, replace = TRUE)
  tbl_used <- tbl[idx_shuffle, ]
  
  stimuli <- tbl_used[, c("x1", "x2")]
  feedback <- tbl_used$category
  
  env <- rlang::current_env()
  list2env(params, env)
  
  l_pred <- predict_rmc_continuous(
    stimuli = stimuli,
    features_cat = features_cat,
    features_cont = features_cont,
    n_values_cat = n_values_cat,
    n_categories = length(unique(tbl_used$category)),
    feedback = feedback,
    params = params,
    previous_learning = NULL,
    assignments = assignments
  )
  tbl_used$preds <- apply(
    l_pred$cat_probs, 1, FUN = function(x) which.max(x)
  )
  tbl_used$n_training <- n
  tbl_used$n_categories <- length(unique(tbl_used$category))
  
  return(tbl_used)
}


summarize_blocks <- function(
  #' summarize categorized stimuli into n_blocks and 
  #' show plot summarized by block if required
  #' 
  #' @param tbl \code{tibble} with category labels and category predictions as columns
  #' @param n_trials_per_block number of traisl per block to summarize
  #' @param show_plot \code{logical} should a raster plot be shown by block
  #' to see categorization predictions along with true categories? defaults to TRUE
  
  #' @return a \code{list} with two tibbles; (a) the summarized results and
  #' (b) the assignments in the originally handed over tibble
  #'
  tbl, 
  n_trials_per_block, 
  show_plot = TRUE
) {
  n_blocks <- ceiling(nrow(tbl)/n_trials_per_block)
  tbl$block_nr <- rep(
    seq(1, n_blocks), 
    each = ceiling(nrow(tbl)/n_blocks)
  )[1:nrow(tbl)]
  # summarize accuracy per block
  tbl_results <- tbl %>%
    group_by(n_categories, n_training, block_nr) %>%
    summarize(
      accuracy = mean(category == preds)
    ) %>% ungroup()
  # prepare tbl for grouped plot
  tbl_long <- tbl %>%
    pivot_longer(
      cols = c(category, preds),
      names_to = "Label",
      values_to = "Value"
    )
  # plots: accuracy over blocks and raster
  pl <- ggplot(tbl_long, aes(x1, x2, group = Value)) +
    geom_raster(aes(fill = Value)) +
    facet_wrap(block_nr ~ Label, ncol = ceiling(n_blocks / 2) * 2)
  if (show_plot) {
    grid.draw(pl)
  }
  return(list(
    tbl_results, # summarized results
    tbl # category assignments
  ))
  
}


plot_block_summary <- function(l) {
  #' summarize categorized stimuli into n_blocks and 
  #' show plot summarized by block if required
  #' 
  #' @param l nested \code{list} containing \code{lists} 
  #' with block summaries and tibbles of all trials
  #' @return nothing, just plots
  #' 
  if (length(l) > 2) {
    l_blocks <- map(l, 1)
    tbl_blocks <- reduce(l_blocks, rbind)
  } else {
    tbl_blocks <- l[[1]]
  }
  tbl_blocks$n_training <- as.factor(tbl_blocks$n_training)
  tbl_summary <- tbl_blocks %>%
    group_by(n_categories, n_training, block_nr) %>%
    summarize(
      accuracy = mean(accuracy)
    )
  
  pl <- ggplot(tbl_summary, aes(block_nr, accuracy, group = n_training)) +
    geom_line(aes(color = n_training)) +
    geom_point(color = "white", size = 3) +
    geom_point(aes(color = n_training)) +
    facet_wrap(~ n_categories, ncol = 3) +
    theme_bw() +
    scale_color_brewer(name = "Nr. Training\nTrials", palette = "Set1") +
    scale_x_continuous(breaks = seq(1, max(tbl_blocks$block_nr))) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = "Block Nr.",
      y = "Proportion Correct"
    )
  grid.draw(pl)
}

wrap_rmc <- function (params, tbl_data, n_categories) { #coupling, phi
  #' wrap prediction function for fitting purposes
  #' 
  params <- list(
    "coupling" = params[[1]],
    "phi" = params[[2]],
    "salience_f" = 1, # not used as all continuous features
    "salience_l" = params[[3]],
    "a_0" = params[[4]],
    "lambda_0" = params[[5]],
    "sigma_sq_0" = (max(tbl_data[, c("x1", "x2")]) / 4) ^ 2,
    "mu_0" = (min(tbl_data[, c("x1", "x2")]) + max(tbl_data[, c("x1", "x2")])) / 2
  )
  l_preds <- predict_rmc_continuous(
    stimuli = tbl_data[, c("x1", "x2")],
    features_cat = c(),
    features_cont = c("x1", "x2"),
    n_values_cat = c(),
    n_categories = n_categories,
    feedback = tbl_data$category,
    params = params,
    print_posterior = FALSE
  )
  neg_ll <- -sum(log(l_preds$cat_probs[cbind(1:nrow(tbl_data), tbl_data$category)]))
  n_cluster <- length(unique(l_preds$assignments))
  
  neg_ll
  #list(neg_ll, n_cluster)
}


summarize_cat_probs <- function(l_pred, tbl_used, n_categories, n_trials) {
  #' summarize categorized stimuli into blocks of n_trials
  #' 
  probs_c <- l_pred$cat_probs[cbind(1:nrow(l_pred$cat_probs), tbl_used$category)]
  tbl_probs <- tibble(
    trial = 1:length(probs_c),
    probability = probs_c
  )
  tbl_probs$block_nr <- cut(
    tbl_probs$trial, c(seq(0, max(tbl_probs$trial), by = n_trials), Inf),
    labels = FALSE
  )
  tbl_probs %>% 
    mutate(
      n_categories = n_categories,
      length_training = max(trial),
      n_clusters = length(unique(l_pred$assignments))
    ) %>%
    group_by(n_categories, length_training, block_nr) %>%
    summarize(
      probability_mn = mean(probability),
      n_clusters = max(n_clusters)
      ) %>%
    ungroup()
}

predict_given_fit <- function(tbl_used, l_fit, n_categories) {
  #' predict from the model given fixed parameters
  #' 
  parms <- l_fit$result$par
  stimuli <- tbl_used[, c("x1", "x2")]
  features_cont <- c("x1", "x2") 
  features_cat <- c() 
  n_values_cat <- NULL
  feedback <- NULL
  print_posterior <- FALSE
  previous_learning <- NULL
  feedback <- tbl_used$category
  params <- list(
    "coupling" = parms[1],
    "phi" = parms[2],
    "salience_f" = 1,
    "salience_l" = parms[3],
    "a_0" = parms[4], #2,
    "lambda_0" = parms[5],
    "sigma_sq_0" = (max(tbl_used[, c("x1", "x2")]) / 4) ^ 2,
    "mu_0" = (min(tbl_used[, c("x1", "x2")]) + max(tbl_used[, c("x1", "x2")])) / 2
  )

  l_pred <- predict_rmc_continuous(
    stimuli = tbl_used[, c("x1", "x2")], 
    features_cat = c(),
    features_cont = c("x1", "x2"),
    n_values_cat = NULL,
    n_categories = n_categories,
    feedback = tbl_used$category,
    params = params,
    previous_learning = NULL, 
    print_posterior = FALSE
  )
  
  return(l_pred)
}



