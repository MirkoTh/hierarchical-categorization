# Modeling Setup

We fit the two models as Bayesian models in the R statistics environment using the cmdstanr package. We compared them using the LOO method implemented in the loo R package.

Two of the used packages cannot be installed from CRAN:

- First, install cmdstanr using these instructions: https://mc-stan.org/cmdstanr/articles/cmdstanr.html
- Second, pull the following github repo: https://github.com/MirkoTh/rutils and install the package from source or use the devtool package and install rtools as follows: devtools::install("your-selected-path\\rutils")

Download the files called "infpro_task-cat_beh" and "infpro_task-cat2_beh" provided by the first author of the associated manuscript. Create a folder called "data" on the root level of the project (i.e., /hierarchical-categorization/data) and copy the downloaded files to that folder.

# Fitting the Two Categorization Models

Run the file R/scripts/fit-all-participants-infpro.R


# Comparing Predictions of the Two Categorization Models for the Feature Inference Task

Run the file R/scripts/feature-inference.R
