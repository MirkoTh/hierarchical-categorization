distance_gcm <- function(i, tbl_x, r) {
  N <- nrow(tbl_x)
  m <- matrix(unlist(rep(tbl_x[i, c("x1", "x2")], N)), N, 2, byrow = TRUE)
  colnames(m) <- c("x1", "x2")
  tbl_single <- as_tibble(m)
  ((tbl_single$x1 - tbl_x$x1)^r +
      (tbl_single$x2 - tbl_x$x2)^r
  )^(1 / r)
}

distance_1d <- function(i, tbl_x, c) {
  N <- nrow(tbl_x)
  m <- matrix(unlist(rep(tbl_x[i, c], N)), N, 1, byrow = TRUE)
  colnames(m) <- c
  tbl_single <- as_tibble(m)
  abs(tbl_single[, c] - tbl_x[, c])
}


write_gcm_stan_file <- function() {
  write_stan_file("
data {
  int n_stim;
  int n_cat;
  array[n_stim] int n_trials; // n trials
  array[n_stim] int n_correct; // n correct categorization responses
  array[n_stim] int n_correct_predict; // n correct categorization responses on test set
  array[n_stim] int n_trials_per_item;
  array[n_stim] int cat; // actual category for a given stimulus
  array[n_stim, n_stim] real<lower=0> d1;
  array[n_stim, n_stim] real<lower=0> d2;

}


parameters {
  real <lower=0> c;
  //real <lower=0,upper=1> w;
  simplex[n_cat] bs;
}

transformed parameters {
  array[n_stim, n_stim] real <lower=0,upper=1> s;
  array[n_stim, n_cat] real <lower=0> sumsim;
  array[n_stim] real <lower=0,upper=1> theta;


  // Similarities
  for (i in 1:n_stim){
    for (k in 1:n_cat) {
      sumsim[i, k] = 0;
    }
    for (j in 1:n_stim){
      //s[i, j] = exp(-square(c)*(w*square(d1[i, j])+(1-w)*square(d2[i, j])));
      s[i, j] = exp(-square(c)*(.5*square(d1[i, j])+(.5)*square(d2[i, j])));
      sumsim[i, cat[j]] = sumsim[i, cat[j]] + s[i,j] * bs[cat[j]];
    }
    theta[i] = sumsim[i, cat[i]] / sum(sumsim[i, ]);
  }
}

model {
  n_correct ~ binomial(n_trials, theta);
  c ~ uniform(0, 10);
  //w ~ beta(1, 1);

}

generated quantities {
  array[n_stim] real log_lik_pred;

  for (n in 1:n_stim) {
  log_lik_pred[n] = binomial_lpmf(n_correct_predict[n] | n_trials_per_item[n], theta[n]);
  }
}

")
}


write_gcm_stan_file_predict <- function() {
  write_stan_file("
data {
  int n_stim;
  int n_cat;
  array[n_stim] int n_trials; // n trials
  array[n_stim] int n_correct; // n correct categorization responses
  array[n_stim] int cat; // actual category for a given stimulus
  array[n_stim, n_stim] real<lower=0> d1;
  array[n_stim, n_stim] real<lower=0> d2;

  int n_stim_predict;
  array[n_stim_predict] int n_trials_predict; // n trials on test set
  array[n_stim_predict] int n_correct_predict; // n correct categorization responses on test set
  array[n_stim_predict] int cat_predict; // actual category for a given stimulus
  array[n_stim_predict, n_stim_predict] real<lower=0> d1_predict;
  array[n_stim_predict, n_stim_predict] real<lower=0> d2_predict;

}


parameters {
  real <lower=0> c;
  //real <lower=0,upper=1> w;
  simplex[n_cat] bs;
}

transformed parameters {
  array[n_stim, n_stim] real <lower=0,upper=1> s;
  array[n_stim, n_cat] real <lower=0> sumsim;
  array[n_stim] real <lower=0,upper=1> theta;


  // Similarities
  for (i in 1:n_stim){
    for (k in 1:n_cat) {
      sumsim[i, k] = 0;
    }
    for (j in 1:n_stim){
      //s[i, j] = exp(-square(c)*(w*square(d1[i, j])+(1-w)*square(d2[i, j])));
      s[i, j] = exp(-square(c)*(.5*square(d1[i, j])+(.5)*square(d2[i, j])));
      sumsim[i, cat[j]] = sumsim[i, cat[j]] + s[i,j] * bs[cat[j]];
    }
    theta[i] = sumsim[i, cat[i]] / sum(sumsim[i, ]);
  }
}

model {
  n_correct ~ binomial(n_trials, theta);
  c ~ uniform(0, 10);
  //w ~ beta(1, 1);

}

generated quantities {
  array[n_stim_predict] real log_lik_pred;
  array[n_stim_predict, n_stim_predict] real <lower=0,upper=1> s_predict;
  array[n_stim_predict, n_cat] real <lower=0> sumsim_predict;
  array[n_stim_predict] real <lower=0,upper=1> theta_predict;


  // Similarities
  for (i in 1:n_stim_predict){
    for (k in 1:n_cat) {
      sumsim_predict[i, k] = 0;
    }
    for (j in 1:n_stim_predict){
      //s_predict[i, j] = exp(-square(c)*(w*square(d1_predict[i, j])+(1-w)*square(d2_predict[i, j])));
      s_predict[i, j] = exp(-square(c)*(.5*square(d1_predict[i, j])+(.5)*square(d2_predict[i, j])));
      sumsim_predict[i, cat_predict[j]] = sumsim_predict[i, cat_predict[j]] + s_predict[i,j] * bs[cat_predict[j]];
    }
    theta_predict[i] = sumsim_predict[i, cat_predict[i]] / sum(sumsim_predict[i, ]);
    log_lik_pred[i] = binomial_lpmf(n_correct_predict[i] | n_trials_predict[i], theta_predict[i]);
  }

}

")
}

write_prototype_stan_file_predict <- function() {
  write_stan_file("
data {
  int n_stim;
  int n_cat;
  array[n_stim] int n_trials; // n trials
  array[n_stim] int n_correct; // n correct categorization responses
  array[n_stim] int cat; // actual category for a given stimulus
  array[n_stim, n_cat] real<lower=0> d1;  // distances dimension 1
  array[n_stim, n_cat] real<lower=0> d2;  // distances dimension 2

  int n_stim_predict;
  array[n_stim_predict] int n_trials_predict; // n trials on test set
  array[n_stim_predict] int n_correct_predict; // n correct categorization responses on test set
  array[n_stim_predict] int cat_predict; // actual category for a given stimulus
  array[n_stim_predict, n_cat] real<lower=0> d1_predict;
  array[n_stim_predict, n_cat] real<lower=0> d2_predict;

}


parameters {
  real <lower=0> c;
  //real <lower=0,upper=1> w;
  simplex[n_cat] bs;
}

transformed parameters {
  array[n_stim, n_cat] real <lower=0,upper=1> s;
  array[n_stim, n_cat] real <lower=0> sim_bs;
  array[n_stim] real <lower=0,upper=1> theta;

  // Similarities
  for (i in 1:n_stim) {
    for (k in 1:n_cat) {
      s[i, k] = exp(-square(c)*(.5*square(d1[i, k])+(.5)*square(d2[i, k])));
      sim_bs[i, k] = s[i,k] * bs[k];
    }
    theta[i] = sim_bs[i, cat[i]] / sum(sim_bs[i, ]);
  }
}

model {
  n_correct ~ binomial(n_trials, theta);
  c ~ uniform(0, 10);
  //w ~ beta(1, 1);

}

generated quantities {
  array[n_stim_predict] real log_lik_pred;
  array[n_stim_predict, n_cat] real <lower=0,upper=1> s_predict;
  array[n_stim_predict, n_cat] real <lower=0> sim_bs_predict;
  array[n_stim_predict] real <lower=0,upper=1> theta_predict;


  // Similarities
  for (i in 1:n_stim_predict) {
    for (k in 1:n_cat) {
      s_predict[i, k] = exp(-square(c)*(.5*square(d1_predict[i,k])+(.5)*square(d2_predict[i,k])));
      sim_bs_predict[i, k] = s_predict[i,k] * bs[k];
    }
    theta_predict[i] = sim_bs_predict[i, cat_predict[i]] / sum(sim_bs_predict[i, ]);
    log_lik_pred[i] = binomial_lpmf(n_correct_predict[i] | n_trials_predict[i], theta_predict[i]);
  }

}

")
}


write_flexprototype_stan_file_predict <- function() {
  write_stan_file("
data {
  int n_obs;
  int n_stim;
  int n_cat;
  int n_dim;  // number of dimensions
  array[n_stim] int n_trials; // n trials
  array[n_stim] int n_correct; // n correct categorization responses
  array[n_stim] int cat; // actual category for a given stimulus
  // array[n_stim, n_cat] real<lower=0> d1;  // distances dimension 1
  // array[n_stim, n_cat] real<lower=0> d2;  // distances dimension 2
  matrix[n_obs, n_dim] y;
  matrix[n_stim, n_dim] y_unique;

  
  int n_obs_predict;
  int n_stim_predict;
  array[n_stim_predict] int n_trials_predict; // n trials on test set
  array[n_stim_predict] int n_correct_predict; // n correct categorization responses on test set
  array[n_stim_predict] int cat_predict; // actual category for a given stimulus
  // array[n_stim_predict, n_cat] real<lower=0> d1_predict;
  // array[n_stim_predict, n_cat] real<lower=0> d2_predict;
  matrix[n_obs_predict, n_dim] y_predict;
  matrix[n_stim_predict, n_dim] y_unique_predict;
}


parameters {
  real <lower=0> c;
  //real <lower=0,upper=1> w;
  simplex[n_cat] bs;
  
  ordered[n_cat] mu1; //category means dim1
  ordered[n_cat] mu2; //category means dim2
}

transformed parameters {
  array[n_stim, n_cat] real<lower=0> d1;  // distances dimension 1
  array[n_stim, n_cat] real<lower=0> d2;  // distances dimension 2

  
  array[n_stim, n_cat] real <lower=0,upper=1> s;
  array[n_stim, n_cat] real <lower=0> sim_bs;
  array[n_stim] real <lower=0,upper=1> theta;

  // Distances
  for (st in 1:n_stim) {
    for (k in 1:n_cat) {
      d1[st, k] = 0;
      d2[st, k] = 0;
      int count1;
      int count2;
      count1 = 0;
      count2 = 0;
      for (i in 1:n_obs) {
        if (y[i,1] == y_unique[st,1]) {
          d1[st, k] += abs(y[i,1] - mu1[k]);
          count1 += 1;
        }
        if (y[i,2] == y_unique[st,2]) {
          d2[st, k] += abs(y[i,2] - mu2[k]);
          count2 += 1;
        }
      }
      d1[st, k] = d1[st, k] / count1;
      d2[st, k] = d2[st, k] / count2;
    }
  }

  // Similarities
  for (i in 1:n_stim) {
    for (k in 1:n_cat) {
      s[i, k] = exp(-square(c)*(.5*square(d1[i, k])+(.5)*square(d2[i, k])));
      sim_bs[i, k] = s[i,k] * bs[k];
    }
    theta[i] = sim_bs[i, cat[i]] / sum(sim_bs[i, ]);
  }
}

model {
  n_correct ~ binomial(n_trials, theta);
  c ~ uniform(0, 10);
  //w ~ beta(1, 1);

  mu1[1] ~ normal(-1.5, .5);
  mu1[2] ~ normal(0, .5);
  mu1[3] ~ normal(1.5, .5);
  mu2[1] ~ normal(-1.5, .5);
  mu2[2] ~ normal(0, .5);
  mu2[3] ~ normal(1.5, .5);
  
}

generated quantities {
  array[n_stim_predict, n_cat] real<lower=0> d1_predict;
  array[n_stim_predict, n_cat] real<lower=0> d2_predict;
  //int <lower=0> count1;
  //int <lower=0> count2;
  array[n_stim_predict] real log_lik_pred;
  array[n_stim_predict, n_cat] real <lower=0,upper=1> s_predict;
  array[n_stim_predict, n_cat] real <lower=0> sim_bs_predict;
  array[n_stim_predict] real <lower=0,upper=1> theta_predict;

  //// Distances
  //for (i in 1:n_stim) {
  //  for (k in 1:n_cat) {
  //    d1_predict[i, k] = abs(y_predict[i,1] - mu1[k]);
  //    d1_predict[i, k] = abs(y_predict[i,2] - mu2[k]);
  //  }
  //}
  
  // Distances
  for (st in 1:n_stim_predict) {
    for (k in 1:n_cat) {
      d1_predict[st, k] = 0;
      d2_predict[st, k] = 0;
      int count1;
      int count2;
      count1 = 0;
      count2 = 0;
      for (i in 1:n_obs_predict) {
        if (y_predict[i,1] == y_unique_predict[st,1]) {
          d1_predict[st, k] += abs(y_predict[i,1] - mu1[k]);
          count1 += 1;
        }
        if (y_predict[i,2] == y_unique_predict[st,2]) {
          d2_predict[st, k] += abs(y_predict[i,2] - mu2[k]);
          count2 += 1;
        }
      }
      d1_predict[st, k] = d1_predict[st, k] / count1;
      d2_predict[st, k] = d2_predict[st, k] / count2;
    }
  }

  // Similarities
  for (i in 1:n_stim_predict) {
    for (k in 1:n_cat) {
      s_predict[i, k] = exp(-square(c)*(.5*square(d1_predict[i,k])+(.5)*square(d2_predict[i,k])));
      sim_bs_predict[i, k] = s_predict[i,k] * bs[k];
    }
    theta_predict[i] = sim_bs_predict[i, cat_predict[i]] / sum(sim_bs_predict[i, ]);
    log_lik_pred[i] = binomial_lpmf(n_correct_predict[i] | n_trials_predict[i], theta_predict[i]);
  }

}

")
}

write_gaussian_naive_bayes_stan <- function() {
  write_stan_file("

data {
 int D; //number of dimensions
 int K; //number of categories
 int N; //number of data
 int n_stim; // nr of different stimuli
 array[K] int n_response_per_cat;
 array[n_stim] int cat_true; // true category for a stimulus
 array[N] int cat; //category response for a stimulus
 matrix[N, D] y;
 matrix[n_stim, D] y_unique;
 array[n_stim] int n_correct_predict; // n correct categorization responses on test set
 array[n_stim] int n_trials_per_item;
}

parameters {
 ordered[K] mu1; //category means d1
 ordered[K] mu2; //category means d2
 array[D, K] real<lower=0> sigma; //variance
 simplex[K] cat_prior;
}

transformed parameters {
 vector[n_stim] theta;

 for (n in 1:n_stim){
   vector[K] LL1_unique;
   vector[K] LL2_unique;
   for (k in 1:K) {
     LL1_unique[k] = normal_lpdf(y_unique[n, 1] | mu1[k], sigma[1, k]);
     LL2_unique[k] = normal_lpdf(y_unique[n, 2] | sort_desc(mu2)[k], sigma[2, k]);
   }
   theta[n] = cat_prior[cat_true[n]] * exp(LL1_unique[cat_true[n]] + LL2_unique[cat_true[n]]) / (
   cat_prior[1] * exp(LL1_unique[1] + LL2_unique[1]) +
   cat_prior[2] * exp(LL1_unique[2] + LL2_unique[2]) +
   cat_prior[3] * exp(LL1_unique[3] + LL2_unique[3])
   );
 }
}

model {

 for(k in 1:K){
   for (d in 1:D) {
     sigma[d, k] ~ uniform(0.1, 5);
   }
 }

 n_response_per_cat ~ multinomial(cat_prior);


 mu1[1] ~ normal(-1.5, .5);
 mu1[2] ~ normal(0, .5);
 mu1[3] ~ normal(1.5, .5);
 mu2[1] ~ normal(-1.5, .5);
 mu2[2] ~ normal(0, .5);
 mu2[3] ~ normal(1.5, .5);

 for (n in 1:N){
 vector[K] LL1;
 vector[K] LL2;
   for (k in 1:K) {
     LL1[k] = normal_lpdf(y[n, 1] | mu1[k], sigma[1, k]);
     LL2[k] = normal_lpdf(y[n, 2] | sort_desc(mu2)[k], sigma[2, k]);
   }
  target += LL1[cat[n]] + LL2[cat[n]] + log(cat_prior[cat[n]]);
 }
}

generated quantities {
  array[n_stim] real log_lik_pred;

  for (n in 1:n_stim) {
  log_lik_pred[n] = binomial_lpmf(n_correct_predict[n] | n_trials_per_item[n], theta[n]);
  }
}


")
}


write_gaussian_naive_bayes_stan_recovery <- function() {
  write_stan_file("

data {
 int D; //number of dimensions
 int K; //number of categories
 int N; //number of data
 int n_stim; // nr of different stimuli
 array[n_stim] int cat_true; // true category for a stimulus
 array[N] int cat; //category response for a stimulus
 array[K] int n_response_per_cat;
 matrix[N, D] y;
 matrix[n_stim, D] y_unique;
 array[n_stim] int n_correct_predict; // n correct categorization responses on test set
 array[n_stim] int n_trials_per_item;

}

parameters {
 array[D] ordered[K] mu; //category means
 array[D, K] real<lower=0> sigma; //variance
 simplex[K] cat_prior;
}

transformed parameters {
 vector[n_stim] theta;

 for (n in 1:n_stim){
   vector[K] LL1_unique;
   vector[K] LL2_unique;
   for (k in 1:K) {
     LL1_unique[k] = normal_lpdf(y_unique[n, 1] | mu[1][k], sigma[1, k]);
     LL2_unique[k] = normal_lpdf(y_unique[n, 2] | mu[2][k], sigma[2, k]);
   }
   theta[n] =  cat_prior[cat_true[n]] * exp(LL1_unique[cat_true[n]] + LL2_unique[cat_true[n]]) / (
   cat_prior[1] * exp(LL1_unique[1] + LL2_unique[1]) +
   cat_prior[2] * exp(LL1_unique[2] + LL2_unique[2]) +
   cat_prior[3] * exp(LL1_unique[3] + LL2_unique[3])
   );
 }
}

model {

 for(k in 1:K){
   for (d in 1:D) {
     sigma[d, k] ~ uniform(0.1, 5);
   }
 }
 mu[1][1] ~ normal(-1.5, .5);
 mu[1][2] ~ normal(0, .5);
 mu[1][3] ~ normal(1.5, .5);
 mu[2][1] ~ normal(-1.5, .5);
 mu[2][2] ~ normal(0, .5);
 mu[2][3] ~ normal(1.5, .5);

 for (n in 1:N){
 vector[K] LL1;
 vector[K] LL2;
   for (k in 1:K) {
     LL1[k] = normal_lpdf(y[n, 1] | mu[1][k], sigma[1, k]);
     LL2[k] = normal_lpdf(y[n, 2] | mu[2][k], sigma[2, k]);
   }
  target += LL1[cat[n]] + LL2[cat[n]] + log(cat_prior[cat[n]]);
 }
}


generated quantities {
  array[n_stim] real log_lik_pred;

  for (n in 1:n_stim) {
  log_lik_pred[n] = binomial_lpmf(n_correct_predict[n] | n_trials_per_item[n], theta[n]);
  }
}

")
}




write_gaussian_multi_bayes_stan <- function() {
  write_stan_file("

data {
 int D; //number of dimensions
 int K; //number of categories
 int N; //number of data
 int n_stim; // nr of different stimuli
 array[n_stim] int cat_true; // true category for a stimulus
 array[N] int cat; //category response for a stimulus
 array[K] int n_response_per_cat;
 matrix[N, D] y;
 matrix[n_stim, D] y_unique;
 array[n_stim] int n_correct_predict; // n correct categorization responses on test set
 array[n_stim] int n_trials_per_item;
}

parameters {
 array[K] ordered[D] mu; //category means
 array[K] cholesky_factor_corr[D] L; //cholesky factor of covariance
 // Sqrt of variances for each variate
 array[K] vector<lower=0>[D] L_std;
 simplex[K] cat_prior;
}

transformed parameters {
 vector[n_stim] theta;

 for (n in 1:n_stim){
   vector[K] LL_unique;
   for (k in 1:K) {
     matrix[D, D] L_Sigma_unique = diag_pre_multiply(L_std[k], L[k]);
     LL_unique[k] = multi_normal_cholesky_lpdf(y_unique[n] | mu[k], L_Sigma_unique);
   }
   theta[n] =  cat_prior[cat_true[n]] * exp(LL_unique[cat_true[n]]) / (
   cat_prior[1] * exp(LL_unique[1]) +
   cat_prior[2] * exp(LL_unique[2]) +
   cat_prior[3] * exp(LL_unique[3])
   );
   //theta[n] = exp(LL_unique[cat_true[n]]) / (exp(LL_unique[1]) + exp(LL_unique[2]) + exp(LL_unique[3]));
 }
}

model {

 for(k in 1:K){
   mu[k] ~ normal(0, 1);
   L[k] ~ lkj_corr_cholesky(1);
   L_std[k] ~ normal(0, 2.5);
 }

 for (n in 1:N){
 vector[K] LL;
 for (k in 1:K) {
     matrix[D, D] L_Sigma = diag_pre_multiply(L_std[k], L[k]);
     LL[k] = multi_normal_cholesky_lpdf(y[n] | mu[k], L_Sigma);
   }
  target += LL[cat[n]] + log(cat_prior[cat[n]]);
 }
}

generated quantities {
 array[K] corr_matrix[D] Sigma;
 array[n_stim] real log_lik_pred;

 for (k in 1:K){
   Sigma[k] = multiply_lower_tri_self_transpose(L[k]);
 }

 for (n in 1:n_stim) {
   log_lik_pred[n] = binomial_lpmf(n_correct_predict[n] | n_trials_per_item[n], theta[n]) + log(cat_prior[cat_true[n]]);
 }
}

")
}


gcm_distances_similarities <- function(tbl_df, l_params) {
  #' gcm distances and similarities for all observations
  #'
  #' @description compute pairwise distances and similarities between
  #' all observations (i.e., rows)
  #'
  #' @param tbl_df tbl df with x1 and x2 features
  #' @param l_params gcm model parameters
  #' @return nested list with distances and similarities
  #'
  # compute pairwise distances
  m_distances_x1 <- map(1:nrow(tbl_df), distance_1d, tbl_df, "x1") %>%
    unlist() %>%
    matrix(byrow = TRUE, nrow = nrow(tbl_df), ncol = nrow(tbl_df))
  m_distances_x2 <- map(1:nrow(tbl_df), distance_1d, tbl_df, "x2") %>%
    unlist() %>%
    matrix(byrow = TRUE, nrow = nrow(tbl_df), ncol = nrow(tbl_df))
  m_similarities <- exp(-l_params[["c"]]^2 * (l_params[["w"]] * m_distances_x1^2 + (1 - l_params[["w"]]) * m_distances_x2^2))
  
  l_dist_sim <- list(
    m_distances_x1 = m_distances_x1,
    m_distances_x2 = m_distances_x2,
    m_similarities = m_similarities
  )
  
  return(l_dist_sim)
}



gcm_response_probabilities <- function(i, tbl_df, m_sims, l_params) {
  #' compute gcm response probabilities given pairwise similarities
  #'
  #' @description computes summed similarity of within category items
  #' divided by summed similarity of all items
  #'
  #' @param i idx of the row, for which to compute the response probability
  #' @param tbl_df tbl df with x1 and x2 features and category column
  #' @param m_sims nxn matrix with pairwise similarities
  #' @param l_params list with model parameters
  #' @return dbl with response probability of being correct
  #'
  sum_sim_category <- function(i, tbl_df) {
    cats_distinct <- unique(tbl_df$category)
    l_mask <- map(cats_distinct, ~ tbl_df$category == .x)
    map_dbl(l_mask, ~ sum(m_sims[i, ][.x]))
  }
  
  sum_sim_categories <- sum_sim_category(i, tbl_df)
  sum_sim_biased <- sum_sim_categories * l_params[["b"]]
  
  sum_sim_biased[tbl_df$category[i]] / sum(sum_sim_biased)
}


my_rbinom <- function(n_trials, prob_correct_true) {
  rbinom(n = 1, size = n_trials, prob = prob_correct_true)
}




aggregate_by_stimulus_and_response <- function(tbl_stim_id, tbl_df) {
  #' aggregate responses by participant, stimulus id, category, and response
  #'
  #' @description make sure categories not responded to are filled with 0s
  #'
  #' @param tbl_stim_id tbl_df containing all training stim_id with
  #' associated x values and categories
  #' @param tbl_df tbl df with all category learning training data
  #' @return aggregated tbl df
  #'
  
  session <- tbl_df$session[1]
  tbl_design <- tbl_stim_id %>%
    crossing(
      response = unique(tbl_df$response),
      participant = unique(tbl_df$participant)
    ) %>%
    relocate(stim_id, .before = category)
  
  tbl_train_agg <- tbl_df %>%
    group_by(participant, session, stim_id, d1i, d2i, d1i_z, d2i_z, category, response) %>%
    summarize(
      n_responses = n()
    )
  
  tbl_ <- tbl_design %>%
    left_join(
      tbl_train_agg,
      by = c(
        "participant", "stim_id", "d1i", "d2i", "d1i_z", "d2i_z", "category", "response"
      )
    )
  
  tbl_$n_responses[is.na(tbl_$n_responses)] <- 0
  tbl_$session[is.na(tbl_$session)] <- session
  tbl_ %>%
    group_by(participant, d1i, d2i) %>%
    mutate(
      n_trials = sum(n_responses),
      prop_responses = n_responses / n_trials
    ) %>%
    ungroup()
}


save_my_png <- function(pl, f_name, vals_size) {
  #' save three plots from the categorization model
  #'
  png(filename = f_name, vals_size[1], vals_size[2], "in", res = 200)
  grid.draw(pl)
  dev.off()
}

bayesian_gcm <- function(tbl_participant, l_stan_params, mod_gcm) {
  #' fit by-participant gcm stan model
  #'
  #' @description fit stan model, save three plots, and return loo
  #'
  #' @param tbl_participant by-participant aggregated responses from infpro task
  #' @param l_stan_params list with parameters used for stan model
  #' @param mod_gcm the compiled cmdstanr model
  #' @return loo
  #'
  
  tbl_train <- tbl_participant %>% filter(session == "train")
  tbl_transfer <- tbl_participant %>% filter(session == "transfer")
  participant_sample <- tbl_train$participant[1]
  tbl_gcm_train <- tbl_train %>%
    filter(category == response) %>%
    mutate(prop_correct = prop_responses)
  tbl_gcm_transfer <- tbl_transfer %>%
    filter(category == response) %>%
    mutate(prop_correct = prop_responses)
  
  # compute pairwise distances
  m_distances_x1_train <- map(1:nrow(tbl_gcm_train), distance_1d, tbl_gcm_train, "d1i_z") %>%
    unlist() %>%
    matrix(byrow = TRUE, nrow = nrow(tbl_gcm_train), ncol = nrow(tbl_gcm_train))
  m_distances_x2_train <- map(1:nrow(tbl_gcm_train), distance_1d, tbl_gcm_train, "d2i_z") %>%
    unlist() %>%
    matrix(byrow = TRUE, nrow = nrow(tbl_gcm_train), ncol = nrow(tbl_gcm_train))
  
  m_distances_x1_transfer <- map(1:nrow(tbl_gcm_transfer), distance_1d, tbl_gcm_transfer, "d1i_z") %>%
    unlist() %>%
    matrix(byrow = TRUE, nrow = nrow(tbl_gcm_transfer), ncol = nrow(tbl_gcm_transfer))
  m_distances_x2_transfer <- map(1:nrow(tbl_gcm_transfer), distance_1d, tbl_gcm_transfer, "d2i_z") %>%
    unlist() %>%
    matrix(byrow = TRUE, nrow = nrow(tbl_gcm_transfer), ncol = nrow(tbl_gcm_transfer))
  
  l_data <- list(
    # train
    n_stim = nrow(tbl_gcm_train),
    n_cat = length(unique(tbl_gcm_train$category)),
    n_trials = tbl_gcm_train$n_trials,
    n_correct = tbl_gcm_train$n_responses,
    cat = tbl_gcm_train$category_int,
    d1 = m_distances_x1_train,
    d2 = m_distances_x2_train,
    # transfer / predict
    n_stim_predict = nrow(tbl_gcm_transfer),
    n_trials_predict = tbl_gcm_transfer$n_trials,
    n_correct_predict = tbl_gcm_transfer$n_responses,
    cat_predict = tbl_gcm_transfer$category_int,
    d1_predict = m_distances_x1_transfer,
    d2_predict = m_distances_x2_transfer
  )
  
  fit_gcm <- mod_gcm$sample(
    data = l_data, chains = l_stan_params$n_chains,
    iter_sampling = l_stan_params$n_samples, iter_warmup = l_stan_params$n_warmup
  )
  
  file_loc <- str_c("data/infpro_task-cat_beh/models/gcm-model-", participant_sample, ".RDS")
  fit_gcm$save_object(file = file_loc)
  
  loo_gcm <- fit_gcm$loo(variables = "log_lik_pred")
  
  pars_interest <- c("theta_predict", "bs", "c") # , "w")
  pars_interest_no_theta <- c("bs", "c") # , "w")
  tbl_draws <- fit_gcm$draws(variables = pars_interest, format = "df")
  
  names_thetas <- names(tbl_draws)[startsWith(names(tbl_draws), "theta_predict")]
  tbl_gcm_transfer$pred_theta <- colMeans(tbl_draws[, names_thetas])
  tbl_gcm_transfer$pred_difference <- tbl_gcm_transfer$pred_theta - tbl_gcm_transfer$prop_responses
  
  
  tbl_summary <- fit_gcm$summary(variables = c("theta", "bs", "c"))
  tbl_summary_nok <- tbl_summary %>% filter(rhat > 1.02 | rhat < 0.98)
  if (nrow(tbl_summary_nok) > 0) {
    stop(str_c(
      "participant = ", participant_sample, "; Rhat for some parameters not ok; ",
      "model can be found under: ",
    ))
  }
  tbl_summary$participant <- participant_sample
  file_loc <- str_c("data/infpro_task-cat_beh/models/gcm-summary-", participant_sample, ".RDS")
  saveRDS(tbl_summary, file_loc)
  
  idx_no_theta <- map(pars_interest_no_theta, ~ str_detect(tbl_summary$variable, .x)) %>%
    reduce(rbind) %>%
    colSums()
  tbl_label <- tbl_summary[as.logical(idx_no_theta), ]
  
  tbl_posterior <- tbl_draws %>%
    dplyr::select(starts_with(pars_interest_no_theta), .chain) %>%
    rename(chain = .chain) %>%
    pivot_longer(starts_with(pars_interest_no_theta), names_to = "parameter", values_to = "value") %>%
    filter(parameter != "chain")
  tbl_posterior$parameter <- fct_inorder(tbl_posterior$parameter)
  
  pl_thetas <- plot_item_thetas(tbl_gcm_transfer, str_c("GCM; Participant = ", participant_sample))
  pl_posteriors <- plot_posteriors(tbl_posterior, tbl_label)
  pl_pred_uncertainty <- plot_proportion_responses(tbl_gcm_transfer, participant_sample, color_pred_difference = TRUE)
  
  # save plots
  c_names <- function(x, y) str_c("data/infpro_task-cat_beh/model-plots/", x, y, ".png")
  l_pl_names <- map(c("gcm-thetas-", "gcm-posteriors-", "gcm-uncertainty-"), c_names, y = participant_sample)
  l_pl <- list(pl_thetas, pl_posteriors, pl_pred_uncertainty)
  l_vals_size <- list(c(3.5, 3.5), c(8.5, 3), c(7.5, 7.5))
  pwalk(list(l_pl, l_pl_names, l_vals_size), save_my_png)
  
  return(loo_gcm)
}

bayesian_prototype <- function(tbl_participant, l_stan_params, mod_prototype) {
  #' fit by-participant gcm stan model
  #'
  #' @description fit stan model, save three plots, and return loo
  #'
  #' @param tbl_participant by-participant aggregated responses from infpro task
  #' @param l_stan_params list with parameters used for stan model
  #' @param mod_gcm the compiled cmdstanr model
  #' @return loo
  #'
  
  tbl_train <- tbl_participant %>% filter(session == "train")
  tbl_transfer <- tbl_participant %>% filter(session == "transfer")
  participant_sample <- tbl_train$participant[1]
  tbl_gcm_train <- tbl_train %>%
    filter(category == response) %>%
    mutate(prop_correct = prop_responses)
  tbl_gcm_transfer <- tbl_transfer %>%
    filter(category == response) %>%
    mutate(prop_correct = prop_responses)
  
  tbl_prototypes <- tbl_gcm_train %>%
    group_by(category) %>%
    summarise(
      pt1 = mean(d1i), pt2 = mean(d2i),
      pt1_z = mean(d1i_z), pt2_z = mean(d2i_z)
    ) %>%
    rename(pt_cat = category)
  
  # compute 1D distances between stimuli and prototypes
  tbl_prototype_dist_train <- tbl_gcm_train %>%
    left_join(tbl_prototypes, by = character()) %>%
    mutate(
      dist1 = abs(d1i - pt1), dist1_z = abs(d1i_z - pt1_z),
      dist2 = abs(d2i - pt2), dist2_z = abs(d2i_z - pt2_z)
    ) %>%
    select(stim_id, pt_cat, dist1_z, dist2_z) %>%
    pivot_wider(names_from = "pt_cat", values_from = c(dist1_z, dist2_z))
  
  tbl_prototype_dist_transfer <- tbl_gcm_transfer %>%
    left_join(tbl_prototypes, by = character()) %>%
    mutate(
      dist1 = abs(d1i - pt1), dist1_z = abs(d1i_z - pt1_z),
      dist2 = abs(d2i - pt2), dist2_z = abs(d2i_z - pt2_z)
    ) %>%
    select(stim_id, pt_cat, dist1_z, dist2_z) %>%
    pivot_wider(names_from = "pt_cat", values_from = c(dist1_z, dist2_z))
  
  # Save prototype distances as matrices
  m_distances_x1_train <- tbl_prototype_dist_train %>%
    select(starts_with("dist1_z")) %>%
    as.matrix() %>%
    unname()
  m_distances_x2_train <- tbl_prototype_dist_train %>%
    select(starts_with("dist2_z")) %>%
    as.matrix() %>%
    unname()
  m_distances_x1_transfer <- tbl_prototype_dist_transfer %>%
    select(starts_with("dist1_z")) %>%
    as.matrix() %>%
    unname()
  m_distances_x2_transfer <- tbl_prototype_dist_transfer %>%
    select(starts_with("dist2_z")) %>%
    as.matrix() %>%
    unname()
  
  
  l_data <- list(
    # train
    n_stim = nrow(tbl_gcm_train),
    n_cat = length(unique(tbl_gcm_train$category)),
    n_trials = tbl_gcm_train$n_trials,
    n_correct = tbl_gcm_train$n_responses,
    cat = tbl_gcm_train$category_int,
    d1 = m_distances_x1_train,
    d2 = m_distances_x2_train,
    # transfer / predict
    n_stim_predict = nrow(tbl_gcm_transfer),
    n_trials_predict = tbl_gcm_transfer$n_trials,
    n_correct_predict = tbl_gcm_transfer$n_responses,
    cat_predict = tbl_gcm_transfer$category_int,
    d1_predict = m_distances_x1_transfer,
    d2_predict = m_distances_x2_transfer
  )
  
  fit_gcm <- mod_prototype$sample(
    data = l_data, chains = l_stan_params$n_chains,
    iter_sampling = l_stan_params$n_samples, iter_warmup = l_stan_params$n_warmup
  )
  
  file_loc <- str_c("data/infpro_task-cat_beh/models/prototype-model-", participant_sample, ".RDS")
  fit_gcm$save_object(file = file_loc)
  
  loo_gcm <- fit_gcm$loo(variables = "log_lik_pred")
  
  pars_interest <- c("theta_predict", "bs", "c") # , "w")
  pars_interest_no_theta <- c("bs", "c") # , "w")
  tbl_draws <- fit_gcm$draws(variables = pars_interest, format = "df")
  
  names_thetas <- names(tbl_draws)[startsWith(names(tbl_draws), "theta_predict")]
  tbl_gcm_transfer$pred_theta <- colMeans(tbl_draws[, names_thetas])
  tbl_gcm_transfer$pred_difference <- tbl_gcm_transfer$pred_theta - tbl_gcm_transfer$prop_responses
  
  
  tbl_summary <- fit_gcm$summary(variables = c("theta", "bs", "c"))
  tbl_summary_nok <- tbl_summary %>% filter(rhat > 1.02 | rhat < 0.98)
  if (nrow(tbl_summary_nok) > 0) {
    stop(str_c(
      "participant = ", participant_sample, "; Rhat for some parameters not ok; ",
      "model can be found under: ",
    ))
  }
  tbl_summary$participant <- participant_sample
  file_loc <- str_c("data/infpro_task-cat_beh/models/prototype-summary-", participant_sample, ".RDS")
  saveRDS(tbl_summary, file_loc)
  
  idx_no_theta <- map(pars_interest_no_theta, ~ str_detect(tbl_summary$variable, .x)) %>%
    reduce(rbind) %>%
    colSums()
  tbl_label <- tbl_summary[as.logical(idx_no_theta), ]
  
  tbl_posterior <- tbl_draws %>%
    dplyr::select(starts_with(pars_interest_no_theta), .chain) %>%
    rename(chain = .chain) %>%
    pivot_longer(starts_with(pars_interest_no_theta), names_to = "parameter", values_to = "value") %>%
    filter(parameter != "chain")
  tbl_posterior$parameter <- fct_inorder(tbl_posterior$parameter)
  
  pl_thetas <- plot_item_thetas(tbl_gcm_transfer, str_c("GCM; Participant = ", participant_sample))
  pl_posteriors <- plot_posteriors(tbl_posterior, tbl_label)
  pl_pred_uncertainty <- plot_proportion_responses(tbl_gcm_transfer, participant_sample, color_pred_difference = TRUE)
  
  # save plots
  c_names <- function(x, y) str_c("data/infpro_task-cat_beh/model-plots/", x, y, ".png")
  l_pl_names <- map(c("prototype-thetas-", "prototype-posteriors-", "prototype-uncertainty-"), c_names, y = participant_sample)
  l_pl <- list(pl_thetas, pl_posteriors, pl_pred_uncertainty)
  l_vals_size <- list(c(3.5, 3.5), c(8.5, 3), c(7.5, 7.5))
  pwalk(list(l_pl, l_pl_names, l_vals_size), save_my_png)
  
  return(loo_gcm)
}

bayesian_flexprototype <- function(tbl_participant, tbl_participant_agg, l_stan_params, mod_prototype) {
  #' fit by-participant gcm stan model
  #'
  #' @description fit stan model, save three plots, and return loo
  #'
  #' @param tbl_participant by-participant by-trial responses from infpro task
  #' @param tbl_participant_agg by-participant aggregated responses from infpro task
  #' @param l_stan_params list with parameters used for stan model
  #' @param mod_gcm the compiled cmdstanr model
  #' @return loo
  #'
  
  tbl_train <- tbl_participant %>% filter(session == "train")
  tbl_transfer <- tbl_participant %>% filter(session == "transfer")
  tbl_train_agg <- tbl_participant_agg %>% filter(session == "train")
  tbl_transfer_agg <- tbl_participant_agg %>% filter(session == "transfer")
  participant_sample <- tbl_train$participant[1]
  tbl_gcm_train <- tbl_train_agg %>%
    filter(category == response) %>%
    mutate(prop_correct = prop_responses)
  tbl_gcm_transfer <- tbl_transfer_agg %>%
    filter(category == response) %>%
    mutate(prop_correct = prop_responses)
  
  
  l_data <- list(
    # train
    n_obs = nrow(tbl_train),
    n_stim = nrow(tbl_gcm_train),  # tbl_gcm_train
    n_cat = length(unique(tbl_gcm_train$category)),
    n_dim = 2,
    n_trials = tbl_gcm_train$n_trials,
    n_correct = tbl_gcm_train$n_responses,
    cat = tbl_gcm_train$category_int,
    y = tbl_train[, c("d1i_z", "d2i_z")] %>% as.matrix(),
    y_unique = tbl_gcm_train[, c("d1i_z", "d2i_z")] %>% as.matrix(),
    
    # transfer / predict
    n_obs_predict = nrow(tbl_transfer),
    n_stim_predict = nrow(tbl_gcm_transfer),  
    n_trials_predict = tbl_gcm_transfer$n_trials,
    n_correct_predict = tbl_gcm_transfer$n_responses,
    cat_predict = tbl_gcm_transfer$category_int,
    y_predict = tbl_transfer[, c("d1i_z", "d2i_z")] %>% as.matrix(),
    y_unique_predict = tbl_gcm_transfer[, c("d1i_z", "d2i_z")] %>% as.matrix()
  )
  
  
  fit_gcm <- mod_prototype$sample(
    data = l_data, chains = l_stan_params$n_chains,
    iter_sampling = l_stan_params$n_samples, iter_warmup = l_stan_params$n_warmup
  )
  
  file_loc <- str_c("data/infpro_task-cat_beh/models/prototype-model-", participant_sample, ".RDS")
  fit_gcm$save_object(file = file_loc)
  
  loo_gcm <- fit_gcm$loo(variables = "log_lik_pred")
  
  pars_interest <- c("theta_predict", "bs", "c") # , "w")
  pars_interest_no_theta <- c("bs", "c") # , "w")
  tbl_draws <- fit_gcm$draws(variables = pars_interest, format = "df")
  
  names_thetas <- names(tbl_draws)[startsWith(names(tbl_draws), "theta_predict")]
  tbl_gcm_transfer$pred_theta <- colMeans(tbl_draws[, names_thetas])
  tbl_gcm_transfer$pred_difference <- tbl_gcm_transfer$pred_theta - tbl_gcm_transfer$prop_responses
  
  
  tbl_summary <- fit_gcm$summary(variables = c("theta", "bs", "c"))
  tbl_summary_nok <- tbl_summary %>% filter(rhat > 1.02 | rhat < 0.98)
  if (nrow(tbl_summary_nok) > 0) {
    stop(str_c(
      "participant = ", participant_sample, "; Rhat for some parameters not ok; ",
      "model can be found under: ",
    ))
  }
  tbl_summary$participant <- participant_sample
  file_loc <- str_c("data/infpro_task-cat_beh/models/prototype-summary-", participant_sample, ".RDS")
  saveRDS(tbl_summary, file_loc)
  
  idx_no_theta <- map(pars_interest_no_theta, ~ str_detect(tbl_summary$variable, .x)) %>%
    reduce(rbind) %>%
    colSums()
  tbl_label <- tbl_summary[as.logical(idx_no_theta), ]
  
  tbl_posterior <- tbl_draws %>%
    dplyr::select(starts_with(pars_interest_no_theta), .chain) %>%
    rename(chain = .chain) %>%
    pivot_longer(starts_with(pars_interest_no_theta), names_to = "parameter", values_to = "value") %>%
    filter(parameter != "chain")
  tbl_posterior$parameter <- fct_inorder(tbl_posterior$parameter)
  
  pl_thetas <- plot_item_thetas(tbl_gcm_transfer, str_c("GCM; Participant = ", participant_sample))
  pl_posteriors <- plot_posteriors(tbl_posterior, tbl_label)
  pl_pred_uncertainty <- plot_proportion_responses(tbl_gcm_transfer, participant_sample, color_pred_difference = TRUE)
  
  # save plots
  c_names <- function(x, y) str_c("data/infpro_task-cat_beh/model-plots/", x, y, ".png")
  l_pl_names <- map(c("prototype-thetas-", "prototype-posteriors-", "prototype-uncertainty-"), c_names, y = participant_sample)
  l_pl <- list(pl_thetas, pl_posteriors, pl_pred_uncertainty)
  l_vals_size <- list(c(3.5, 3.5), c(8.5, 3), c(7.5, 7.5))
  pwalk(list(l_pl, l_pl_names, l_vals_size), save_my_png)
  
  return(loo_gcm)
}

bayesian_gaussian_naive_bayes <- function(tbl_participant, tbl_participant_agg, l_stan_params, mod_gaussian) {
  #' fit by-participant 1D Gaussian naive bayes stan model
  #'
  #' @description fit stan model, save three plots, and return loo
  #'
  #' @param tbl_participant by-participant and by-trial responses from infpro task
  #' @param tbl_participant_agg by-participant aggregated responses from infpro task
  #' @param l_stan_params list with parameters used for stan model
  #' @param mod_gaussian the compiled cmdstanr model
  #' @return loo
  #'
  
  tbl_train <- tbl_participant %>% filter(session == "train")
  tbl_participant_agg <- tbl_participant_agg %>%
    filter(response == category & session == "transfer") %>%
    mutate(prop_correct = prop_responses)
  participant_sample <- tbl_train$participant[1]
  
  l_data <- list(
    D = 2, K = length(unique(tbl_train$category)),
    N = nrow(tbl_train),
    y = tbl_train[, c("d1i_z", "d2i_z")] %>% as.matrix(),
    cat = tbl_train$response_int,
    n_response_per_cat = tbl_train %>% group_by(response_int) %>% count() %>% ungroup() %>% dplyr::select(n) %>% as_vector(),
    cat_true = tbl_participant_agg$response_int,
    n_stim = nrow(tbl_participant_agg),
    y_unique = tbl_participant_agg[, c("d1i_z", "d2i_z")] %>% as.matrix(),
    n_correct_predict = tbl_participant_agg$n_responses,
    n_trials_per_item = tbl_participant_agg$n_trials
  )
  fit_gaussian <- mod_gaussian$sample(
    data = l_data, chains = l_stan_params$n_chains,
    iter_sampling = l_stan_params$n_samples, iter_warmup = l_stan_params$n_warmup
  )
  file_loc <- str_c(
    "data/infpro_task-cat_beh/models/gaussian-model-", participant_sample, ".RDS"
  )
  fit_gaussian$save_object(file = file_loc)
  
  # loo
  loo_gaussian <- fit_gaussian$loo(variables = "log_lik_pred")
  
  pars_interest <- c("mu1", "mu2", "sigma", "theta", "cat_prior")
  pars_interest_no_theta <- c("mu1", "mu2", "sigma", "cat_prior")
  tbl_draws <- fit_gaussian$draws(variables = pars_interest, format = "df")
  tbl_draws <- tbl_draws %>% rename(`mu2[1]` = `mu2[3]`, `mu2[3]` = `mu2[1]`)
  
  names_thetas <- names(tbl_draws)[startsWith(names(tbl_draws), "theta")]
  tbl_participant_agg$pred_theta <- colMeans(tbl_draws[, names_thetas])
  tbl_participant_agg$pred_difference <- tbl_participant_agg$pred_theta - tbl_participant_agg$prop_responses
  
  
  tbl_summary <- fit_gaussian$summary(variables = pars_interest)
  # rename to correct for reverse ordering in stan model
  old_3 <- tbl_summary$variable == "mu2[3]"
  old_1 <- tbl_summary$variable == "mu2[1]"
  tbl_summary$variable[old_3] <- "mu2[1]"
  tbl_summary$variable[old_1] <- "mu2[3]"
  
  tbl_summary_nok <- tbl_summary %>% filter(rhat > 1.02 | rhat < 0.98)
  if (nrow(tbl_summary_nok) > 0) {
    stop(str_c(
      "participant = ", participant_sample, "; Rhat for some parameters not ok",
      "model can be found under: ", file_loc
    ))
  }
  
  tbl_summary$participant <- participant_sample
  file_loc <- str_c("data/infpro_task-cat_beh/models/gaussian-summary-", participant_sample, ".RDS")
  saveRDS(tbl_summary, file_loc)
  
  idx_no_theta <- map(pars_interest_no_theta, ~ str_detect(tbl_summary$variable, .x)) %>%
    reduce(rbind) %>%
    colSums()
  tbl_label <- tbl_summary[as.logical(idx_no_theta), ]
  
  tbl_posterior <- tbl_draws %>%
    dplyr::select(starts_with(pars_interest_no_theta), .chain) %>%
    rename(chain = .chain) %>%
    pivot_longer(starts_with(pars_interest_no_theta), names_to = "parameter", values_to = "value") %>%
    filter(parameter != "chain")
  tbl_posterior$parameter <- fct_inorder(tbl_posterior$parameter)
  
  pl_thetas <- plot_item_thetas(tbl_participant_agg, str_c("Gaussian 1D; Participant = ", participant_sample))
  pl_posteriors <- plot_posteriors(tbl_posterior, tbl_label, n_cols = 3)
  pl_pred_uncertainty <- plot_proportion_responses(tbl_participant_agg, participant_sample, color_pred_difference = TRUE)
  
  # save plots
  c_names <- function(x, y) str_c("data/infpro_task-cat_beh/model-plots/", x, y, ".png")
  l_pl_names <- map(
    c("gaussian-1d-thetas-", "gaussian-1d-posteriors-", "gaussian-1d-uncertainty-"),
    c_names,
    y = participant_sample
  )
  l_pl <- list(pl_thetas, pl_posteriors, pl_pred_uncertainty)
  l_vals_size <- list(c(3.5, 3.5), c(7.5, 8.5), c(7.5, 7.5))
  pwalk(list(l_pl, l_pl_names, l_vals_size), save_my_png)
  
  return(loo_gaussian)
}


bayesian_gaussian_multi_bayes <- function(tbl_participant, tbl_participant_agg, l_stan_params, mod_multi) {
  #' fit by-participant multivariate 2D Gaussian stan model
  #'
  #' @description fit stan model, save three plots, and return loo
  #'
  #' @param tbl_participant by-participant and by-trial responses from infpro task
  #' @param tbl_participant_agg by-participant aggregated responses from infpro task
  #' @param l_stan_params list with parameters used for stan model
  #' @param mod_multi the compiled cmdstanr model
  #' @return loo
  #'
  
  tbl_train <- tbl_participant %>% filter(session == "train")
  tbl_participant_agg <- tbl_participant_agg %>%
    filter(response == category & session == "transfer") %>%
    mutate(prop_correct = prop_responses)
  participant_sample <- tbl_train$participant[1]
  
  l_data <- list(
    D = 2, K = length(unique(tbl_train$category)),
    N = nrow(tbl_train),
    y = tbl_train[, c("d1i_z", "d2i_z")] %>% as.matrix(),
    cat = tbl_train$response_int,
    n_response_per_cat = tbl_train %>% group_by(response_int) %>% count() %>% ungroup() %>% dplyr::select(n) %>% as_vector(),
    cat_true = tbl_participant_agg$response_int,
    n_stim = nrow(tbl_participant_agg),
    y_unique = tbl_participant_agg[, c("d1i_z", "d2i_z")] %>% as.matrix(),
    n_correct_predict = tbl_participant_agg$n_responses,
    n_trials_per_item = tbl_participant_agg$n_trials
  )
  fit_multi <- mod_multi$sample(
    data = l_data, chains = l_stan_params$n_chains,
    iter_sampling = l_stan_params$n_samples, iter_warmup = l_stan_params$n_warmup
  )
  file_loc <- str_c(
    "data/infpro_task-cat_beh/models/multi-model-", participant_sample, ".RDS"
  )
  fit_multi$save_object(file = file_loc)
  
  # loo
  loo_multi <- fit_multi$loo(variables = "log_lik_pred")
  
  
  pars_interest <- c("mu", "Sigma", "theta", "cat_prior")
  tbl_draws <- fit_multi$draws(variables = pars_interest, format = "df")
  names_params <- names(tbl_draws)
  names_sigmas <- names_params[startsWith(names_params, "Sigma")]
  sigmas_keep <- map(c("2,2]", "1,1]"), ~ as.integer(endsWith(names_sigmas, .x))) %>%
    reduce(rbind) %>%
    colSums()
  names_sigmas_keep <- names_sigmas[as.logical(abs(sigmas_keep - 1))]
  pars_interest_no_theta <- c("mu", names_sigmas_keep, "cat_prior")
  pars_interest_no_theta_idx <- str_replace(pars_interest_no_theta, "\\[", "\\\\[")
  
  
  names_thetas <- names(tbl_draws)[startsWith(names(tbl_draws), "theta")]
  tbl_participant_agg$pred_theta <- colMeans(tbl_draws[, names_thetas])
  tbl_participant_agg$pred_difference <- tbl_participant_agg$pred_theta - tbl_participant_agg$prop_responses
  
  
  tbl_summary <- fit_multi$summary(variables = pars_interest)
  tbl_summary_nok <- tbl_summary %>% filter(rhat > 1.02 | rhat < 0.98)
  if (nrow(tbl_summary_nok) > 0) {
    stop(str_c(
      "participant = ", participant_sample, "; Rhat for some parameters not ok",
      "model can be found under: ", file_loc
    ))
  }
  
  tbl_summary$participant <- participant_sample
  file_loc <- str_c("data/infpro_task-cat_beh/models/multi-model-summary-", participant_sample, ".RDS")
  saveRDS(tbl_summary, file_loc)
  
  idx_no_theta <- map(pars_interest_no_theta_idx, ~ str_detect(tbl_summary$variable, .x)) %>%
    reduce(rbind) %>%
    colSums()
  tbl_label <- tbl_summary[as.logical(idx_no_theta), ]
  
  tbl_posterior <- tbl_draws %>%
    dplyr::select(starts_with(pars_interest_no_theta), .chain) %>%
    rename(chain = .chain) %>%
    pivot_longer(starts_with(pars_interest_no_theta), names_to = "parameter", values_to = "value") %>%
    filter(parameter != "chain")
  tbl_posterior$parameter <- fct_inorder(tbl_posterior$parameter)
  
  
  pl_thetas <- plot_item_thetas(tbl_participant_agg, str_c("MV Gaussian; Participant = ", participant_sample))
  pl_posteriors <- plot_posteriors(tbl_posterior, tbl_label, n_cols = 5)
  pl_pred_uncertainty <- plot_proportion_responses(tbl_participant_agg, participant_sample, color_pred_difference = TRUE)
  
  
  
  # save plots
  c_names <- function(x, y) str_c("data/infpro_task-cat_beh/model-plots/", x, y, ".png")
  l_pl_names <- map(
    c("multi-gaussian-thetas-", "multi-gaussian-posteriors-", "multi-gaussian-uncertainty-"),
    c_names,
    y = participant_sample
  )
  l_pl <- list(pl_thetas, pl_posteriors, pl_pred_uncertainty)
  l_vals_size <- list(c(3.5, 3.5), c(11, 5.5), c(7.5, 7.5))
  pwalk(list(l_pl, l_pl_names, l_vals_size), save_my_png)
  
  return(loo_multi)
}



max_sim_responses <- function(i_cat, i_dim, post_c, l_tbl_lookup, l_tbl_exemplars, l_tbl_cues = NULL) {
  #' iterates over max_sim_response for all possible cue values
  #' on the given dimension
  #'
  #' @description compute maximally similar responses for given
  #' category and dimension
  #'
  #' @param i_cat category value
  #' @param i_dim dimension value
  #' @param post_c posterior samples of gcm c parameter
  #' @param l_tbl_lookup tbl df containing lookup table with
  #' grid of possible response values on both dimensions
  #' @param l_tbl_exemplars tbl df with all exemplars encountered
  #' during category learning
  #' @param l_tbl_cues list containing tbl dfs with inference cues;
  #' optional, if nothing is provided, just uses all possible cues
  #' from the category learning training data (all available x1s and x2s)
  #' @return tbl df with maximally similar responses for all cues
  #'
  # iterate over all cues from the given dimension
  if (!is.null(l_tbl_cues)) l_tbl_range <- l_tbl_cues
  if (is.null(l_tbl_cues)) l_tbl_range <- l_tbl_exemplars
  i_cue <- seq(1, nrow(l_tbl_range[[i_cat]]), by = 1)
  map(
    i_cue,
    max_sim_response,
    post_c = post_c,
    l_tbl_lookup = l_tbl_lookup,
    l_tbl_exemplars = l_tbl_exemplars,
    i_cat = i_cat,
    i_dim = i_dim,
    l_tbl_cues = l_tbl_cues
  ) %>% reduce(rbind)
}

max_sim_response <- function(i_cue, i_cat, i_dim, post_c, l_tbl_lookup, l_tbl_exemplars, l_tbl_cues = NULL) {
  #' @description compute maximally similar response
  #' for given category, dimension, and cue
  #'
  #' @param i_cue cue value
  #' @param i_cat category value
  #' @param i_dim dimension value
  #' @param post_c posterior samples of gcm c parameter
  #' @param l_tbl_lookup tbl df containing lookup table with
  #' grid of possible response values on both dimensions
  #' @param l_tbl_exemplars tbl df with all exemplars encountered
  #' during category learning
  #' @param l_tbl_cues list containing tbl dfs with inference cues;
  #' optional, if nothing is provided, just uses all possible cues
  #' from the category learning training data (all available x1s and x2s)
  #' @return one row tbl df with maximally similar response
  #'
  if (!is.null(l_tbl_cues)) {
    cue <- l_tbl_cues[[i_cat]][i_cue, i_dim] %>% as_vector()
  } else {
    cue <- l_tbl_exemplars[[i_cat]][i_cue, i_dim] %>% as_vector()
  }
  # then compute the similarities to all exemplars from all available values on the grid
  rows_selected <- as.logical(l_tbl_lookup[[i_cat]][, i_dim] == cue)
  grid_vals <- l_tbl_lookup[[i_cat]][rows_selected, ]
  
  sum_within_category_similarity <- function(x1, x2, c, tbl_df) {
    sum(map_dbl(c, ~ sum(exp(-.x^2 * ((x1 - tbl_df$d1i_z)^2 + (x2 - tbl_df$d2i_z)^2)))))
  }
  
  idx_min <- which.max(map2_dbl(
    grid_vals$d1i_z, grid_vals$d2i_z,
    sum_within_category_similarity,
    c = as.list(post_c$c),
    tbl_df = l_tbl_exemplars[[i_cat]]
  ))
  
  return(tibble(grid_vals[idx_min, ], cue_dim = i_dim))
}


varied_cues <- function(tbl_df) {
  #' @description cross cue-values from one dimension with fine grid
  #' of values from other dimension
  #'
  #' @param tbl_df tbl df with unique cue values from 1D and NAs on other dimension
  #' @return tbl df with seen cues crossed with fine grid from other dimension
  #'
  mask <- which(!is.na(tbl_df[1, ])) %>% as_vector()
  tbl_df_filtered <- tbl_df[, mask]
  dim_cued <- names(tbl_df_filtered)[2]
  if (dim_cued == "d1i_z") dim_response <- "d2i_z"
  if (dim_cued == "d2i_z") dim_response <- "d1i_z"
  
  tbl_cross <- crossing(
    dim_cued = c(unique(tbl_df_filtered[, dim_cued]) %>% as_vector(), seq(min(tbl_df_filtered[, dim_cued]), max(tbl_df_filtered[, dim_cued]), by = .1)),
    dim_response = seq(-1.5, 1.5, by = .05)
  )
  names(tbl_cross) <- c(dim_cued, dim_response)
  return(tbl_cross)
}





load_parameter_posteriors <- function(p_id) {
  #' @description load posterior samples from gcm and gaussian nb for one participant
  #'
  file_loc_gcm <- str_c("data/infpro_task-cat_beh/models/gcm-model-", p_id, ".RDS")
  file_loc_gaussian <- str_c(
    "data/infpro_task-cat_beh/models/gaussian-model-", p_id, ".RDS"
  )
  m_gcm <- readRDS(file_loc_gcm)
  m_pt <- readRDS(file_loc_gaussian)
  post_c <- m_gcm$draws(variables = "c", format = "df") %>% as_tibble()
  post_pts <- m_pt$draws(variables = c("mu1", "mu2"), format = "df")
  
  return(list(gcm = post_c, gaussian = post_pts))
}



model_based_inference_responses_gcm <- function(tbl_completion, tbl_train, p_id, l_pars_tf) {
  #' model-based inferences given gcm model posteriors
  #'
  #' @description computes inference responses that maximize within-category
  #' similarity given posterior samples from c parameter of gcm
  #'
  #' @param tbl_completion empirical completion data
  #' @param tbl_train empirical category learning data from training
  #' @param p_id participant id
  #' @param l_pars_tf means and sds in untransformed space
  #' @return list containing tbl dfs with maximally similar responses
  #' for all categories and cues
  #'
  
  # all the exemplars observed during training that can be referred to in memory
  l_tbl_exemplars <- tbl_train %>%
    mutate(category = fct_recode(category, B = "C", C = "B")) %>%
    filter(participant == p_id) %>%
    group_by(category, d1i_z, d2i_z) %>%
    count() %>%
    select(-n) %>%
    split(.$category)
  
  # load posterior samples
  l_posteriors <- load_parameter_posteriors(p_id)
  post_c <- l_posteriors$gcm
  
  # get unique inference cues and bring them in format to compute model-based responses
  l_tbl_cues <- tbl_completion %>%
    mutate(cue_val_dupl = cue_val, cuedim_dupl = cuedim) %>%
    pivot_wider(
      id_cols = c(participant, category, cuedim_dupl, cue_val_dupl, respdim, rep, resp_i),
      names_from = cuedim, values_from = cue_val, names_glue = "{cuedim}_{.value}"
    ) %>%
    mutate(
      d1i_z = (`1_cue_val` - l_pars_tf$mean_d1i) / l_pars_tf$sd_d1i,
      d2i_z = (`2_cue_val` - l_pars_tf$mean_d2i) / l_pars_tf$sd_d2i,
    ) %>%
    filter(participant == p_id) %>%
    select(-c(
      participant, cuedim_dupl, cue_val_dupl, respdim,
      rep, resp_i, `1_cue_val`, `2_cue_val`
    )) %>%
    group_by(category, d1i_z, d2i_z) %>%
    count() %>%
    select(-n) %>%
    split(.$category)
  
  # expand the not-varied dimension as a grid for a given participant
  l_tbl_lookup <- map(l_tbl_cues, varied_cues)
  
  tbl_cat_dim <- crossing(i_cat = c("A", "B"), i_dim = names(l_tbl_lookup$A)[1])
  
  l_closest <- pmap(
    tbl_cat_dim,
    max_sim_responses,
    post_c = post_c,
    l_tbl_lookup = l_tbl_lookup,
    l_tbl_exemplars = l_tbl_exemplars,
    l_tbl_cues = l_tbl_cues
  )
  names(l_closest) <- tbl_cat_dim$i_cat
  
  return(l_closest)
}



model_based_inference_responses_gaussian <- function(tbl_completion, tbl_train, p_id, l_pars_tf) {
  #' model-based inferences given Gaussian prototype model posteriors
  #'
  #' @description computes inference responses that minimize distances
  #' to representational prototype given posterior distributions of pts
  #'
  #' @param tbl_completion empirical completion data
  #' @param tbl_train empirical category learning data from training
  #' @param p_id participant id
  #' @param l_pars_tf means and sds in untransformed space
  #' @return list containing tbl dfs with maximally similar responses
  #' for all categories and cues
  #'
  
  # all the exemplars observed during training that can be referred to in memory
  l_tbl_exemplars <- tbl_train %>%
    mutate(category = fct_recode(category, B = "C", C = "B")) %>%
    filter(participant == p_id) %>%
    group_by(category, d1i_z, d2i_z) %>%
    count() %>%
    select(-n) %>%
    split(.$category)
  
  # load posterior samples
  l_posteriors <- load_parameter_posteriors(p_id)
  post_gaussian <- l_posteriors$gaussian
  post_gaussian <- post_gaussian %>% rename(`mu2[1]` = `mu2[3]`, `mu2[3]` = `mu2[1]`) %>% select(-c(.chain, .iteration, .draw, `mu1[2]`, `mu2[2]`))
  post_gaussian_mns <- apply(post_gaussian, 2, mean)
  idx_c1 <- str_detect(names(post_gaussian_mns), "\\[1\\]")
  idx_c3 <- str_detect(names(post_gaussian_mns), "\\[3\\]")
  
  
  # get unique inference cues and bring them in format to compute model-based responses
  l_tbl_cues <- tbl_completion %>%
    mutate(cue_val_dupl = cue_val, cuedim_dupl = cuedim) %>%
    pivot_wider(
      id_cols = c(participant, category, cuedim_dupl, cue_val_dupl, respdim, rep, resp_i),
      names_from = cuedim, values_from = cue_val, names_glue = "{cuedim}_{.value}"
    ) %>%
    mutate(
      d1i_z = (`1_cue_val` - l_pars_tf$mean_d1i) / l_pars_tf$sd_d1i,
      d2i_z = (`2_cue_val` - l_pars_tf$mean_d2i) / l_pars_tf$sd_d2i,
    ) %>%
    filter(participant == p_id) %>%
    select(-c(
      participant, cuedim_dupl, cue_val_dupl, respdim,
      rep, resp_i, `1_cue_val`, `2_cue_val`
    )) %>%
    group_by(category, d1i_z, d2i_z) %>%
    count() %>%
    select(-n) %>%
    split(.$category)
  
  # expand the not-varied dimension as a grid for a given participant
  l_tbl_lookup <- map(l_tbl_cues, varied_cues)
  tbl_cat_dim <- crossing(i_cat = c("A", "B"), i_dim = names(l_tbl_lookup$A)[1])
  
  
  tbl_pts <- data.frame(rbind(
    post_gaussian_mns[idx_c1],
    post_gaussian_mns[idx_c3]
  )) %>% as_tibble()
  colnames(tbl_pts) <- c("d1i_z", "d2i_z")
  tbl_pts$category <- c("A", "B")
  
  
  closest_response_to_prototype <- function(i_cat, i_dim, l_tbl_cue, tbl_pts) {
    resp_dim <- c("d1i_z", "d2i_z")[!(c("d1i_z", "d2i_z") == i_dim)]
    response_pt <- tbl_pts %>% 
      filter(category == i_cat) %>%
      select(all_of(resp_dim)) %>%
      as_vector()
    cues <- l_tbl_cues[[i_cat]]
    cues[, resp_dim] <- response_pt
    cues$cue_dim <- i_dim
    return(cues)
  }
  
  l_closest <- pmap(
    tbl_cat_dim,
    closest_response_to_prototype,
    l_tbl_cue = l_tbl_cues,
    tbl_pts = tbl_pts
  )
  names(l_closest) <- tbl_cat_dim$i_cat
  
  return(l_closest)
}


distance_from_model_based_inference <- function(p_id, tbl_completion, tbl_train, l_pars_tf, modeltype) {
  #' compute gcm-based responses and compare them to empirical responses
  #'
  #' @description computes for every 1D inference cue what response GCM would give
  #'
  #' @param p_id participant
  #' @param tbl_completion tbl df with empirical inference data
  #' @param tbl_train tbl df with category learning data from training session
  #' @param l_pars_tf mean and sds in untransformed space
  #' @param modeltype gcm or gaussian category learning model
  #' @return tbl df with distance as additional column
  #'
  if (modeltype == "gcm"){
    l_model_based <- model_based_inference_responses_gcm(tbl_completion, tbl_train, p_id, l_pars_tf)
  } else if (modeltype == "gaussian") {
    l_model_based <- model_based_inference_responses_gaussian(tbl_completion, tbl_train, p_id, l_pars_tf)
  }
  
  add_category <- function(tbl_df, cat) {
    tbl_df$category <- cat
    return(tbl_df)
  }
  # transform back into original 2D space
  tbl_model_based <- map2(l_model_based, names(l_model_based), add_category) %>%
    reduce(rbind) %>%
    mutate(
      d1i = d1i_z * l_pars_tf$sd_d1i + l_pars_tf$mean_d1i,
      d2i = d2i_z * l_pars_tf$sd_d2i + l_pars_tf$mean_d2i,
      cuedim = c(1, 2)[as.integer(cue_dim == "d2i_z") + 1],
      participant = p_id
    ) %>%
    select(-c(d1i_z, d2i_z, cue_dim)) %>%
    relocate(cuedim, .before = d1i) %>%
    relocate(participant, .before = cuedim)
  # select right cue value for given row
  tbl_model_based$cue_val <- pmap_dbl(
    tbl_model_based[, c("d1i", "d2i", "cuedim")], ~ c(..1, ..2)[..3]
  )
  
  # EMPIRICAL
  # join gcm-based responses back into empirical tbl df
  tbl_model <- tbl_model_based %>%
    left_join(tbl_completion, by = c("participant", "category", "cuedim", "cue_val"))
  # compute distance from model-based response
  tbl_model$distance <- pmap_dbl(
    tbl_model[, c("d1i", "d2i", "respdim", "resp_i")],
    ~ c(..1, ..2)[..3] - ..4
  )
  
  # LOOKUP TABLE
  tbl_lookup <- lookup_table_possible_responses(tbl_completion, tbl_model_based, p_id)
  # compute distance from model-based response
  tbl_lookup$distance <- pmap_dbl(
    tbl_lookup[, c("d1i", "d2i", "respdim", "resp_i")],
    ~ c(..1, ..2)[..3] - ..4
  )
  
  # plot and save histograms by category
  pl <- ggplot(
    tbl_model %>% mutate(category = str_c("Category ", category)), aes(distance)
  ) +
    geom_histogram(color = "white", fill = "#66CCFF", bins = 15) +
    facet_wrap(~category) +
    theme_bw() +
    labs(
      x = str_c("Distance from ", modeltype,"-Based Response"),
      y = "Nr. Responses",
      title = str_c("Participant = ", p_id)
    )
  f_name_pl <- str_c("data/infpro_task-cat_beh/figures/ds-", modeltype, "-based-p-", p_id, ".png")
  save_my_png(pl, f_name_pl, c(4, 3))
  
  return(list(tbl_empirical = tbl_model, tbl_lookup = tbl_lookup))
}



lookup_table_possible_responses <- function(tbl_completion, tbl_model_based, p_id) {
  #' expand a grid on all possible completion responses for a given participant
  #'
  #' @description joins fine grid of possible response values
  #' for all cue values
  #'
  #' @param tbl_completion tbl df with empirical inference data
  #' @param tbl_model_based tbl df with model-based responses
  #' @param p_id participant id
  #' @return tbl df with distance as additional column
  #'
  tbl_cues <- tbl_completion %>%
    filter(participant == p_id) %>%
    group_by(participant, category, cuedim, cue_val, respdim) %>%
    count() %>%
    select(cue_val) %>%
    ungroup()
  tbl_design <- crossing(
    cue_val = unique(tbl_cues$cue_val), resp_i = seq(-5, 15, by = .05)
  )
  tbl_design <- tbl_design %>% left_join(tbl_cues, by = "cue_val")
  tbl_lookup <- tbl_model_based %>%
    left_join(tbl_design, by = c("participant", "category", "cuedim", "cue_val"))
  
  return(tbl_lookup)
}
