# Assignment 03
$$
\DeclareMathOperator{\cor}{cor}
\DeclareMathOperator{\cov}{cov}
\DeclareMathOperator{\quantile}{quantile}
$$

## Instructions

1. [Fork this repository](https://help.github.com/articles/using-pull-requests/) to your GitHub account.
2. Write your solutions in R Markdown in a file named `solutions.Rmd`.
3. When you are ready to submit your assignment, [initiate a pull request](https://help.github.com/articles/using-pull-requests/#initiating-the-pull-request). Title your
pull request "Submission".

To update your fork from the upstream repository:

1. On your fork, e.g. `https://github.com/jrnold/Assignment_03` click on "New Pull reqest"
2. Set your fork `jrnold/Assignment_03` as the base fork on the left, and `UW-POLS503/Assignment_03` as the head fork on the right. In both cases the branch will be master. This means, compare any chanes in the head fork that are not in the base fork. You will see differences between the `US-POLS503` repo and your fork. Click on "Create Pull Request", and if there are no issues, "Click Merge" A quick way is to use this link, but change the `jrnold` to your own username: `https://github.com/jrnold/Assignment_03/compare/master...UW-POLS503:master`.

We'll use these packages,

```r
library("foreign")
library("dplyr")
library("broom")
library("ggplot2")
```
Since we are going to do some simulation, we shoudl set a seed, so the results are exactly replicable.

```r
set.seed(1234)
```
Since some of these computations will take time, we can cache the results so that knitr will
only run code that has changed.

```r
knitr::opts_chunk$set(cache = TRUE, autodep = TRUE)
```

## Problem 1

Let's run some regressions from Nunn and Wantchekon!

```r
nunn <- read.dta("Nunn_Wantchekon_AER_2011.dta") %>% tbl_df()
```

- `trust_neighbors`: Trust of neighbors
- `exports`: Slave exports in 1000s
- `export_pop`: 
- Individual controls: `age`, `age2`, `male`, `urban_dum`, `education`, `occupation`, `religion`, `living_conditions`
- District controls: `district_ethnic_frac`, `frac_ethnicity_in_district`
- Country-fixed effects `isocode`


```r
mod1 <- lm(trust_neighbors ~ exports, data = nunn) 
```

Interpret the coefficient's magnitude and statistical significance.

Some questions on p-values

- What is the null hypothesis of the t-tests
- Explain the meaning of the p-value
- Is the p-value the probability that the null hypothesis is correct?

What other variables does Nunn include? Include those in Table 1, Model 1. Run that regression.

Does the R^2 match that in Table 1?

Do the standard errors match those in Table 1? Any guesses why?

Run the regression in Table 1, model 6. Why is it "log(1 + exports / pop)" instead of "log(exports / pop)"?


```r
resampler_coef <- function(mod, .data, iter = 1) {
  # Remove missing values
  .data <- na.omit(.data)
  # Coefficients
  beta <- coef(mod)
  # mod$terms contains the formula used in the regression
  X <- model.matrix(mod$terms, data = .data)
  # estimate of std. dev. of errors
  sigma <- sqrt(sum(mod$residuals ^ 2) / mod$df.residual)
  # This produces the same result
  # sigma <- summary(mod1)$sigma  
  # Number of observations
  n <- nrow(X)
  # Name of dependent variable
  outcome_var_name <- all.vars(mod$terms)[1]
  # List to save results
  results <- vector(mode = "list", length = iter)
  for (i in seq_len(iter)) {
    # draw errors
    errors <- rnorm(n, mean = 0, sd = sigma)
    # create new outcome variable from errors
    y <- X %*% beta + errors
    # replace outcome variable
    .data[[outcome_var_name]] <- y
    # run regression
    newmod <- lm(mod$terms, data = .data)
    # Save coefficients as a data frame to the list
    results[[i]] <- tidy(newmod) %>% mutate(.iter = i)
  }
  # Convert the list of data frames to a single data frame by stacking the iterations
  bind_rows(results)
}
```

- Plot the distributions of the coefficients
- Calculate the correlation matrix of the coefficients. How similar is it to that from `vcov`?


### F-test example

- Run F-tests of the multiple regression model vs. the model with no controls.
- Run and interpet an F-test on some reasonable group of variables.

F-test simulations


```r
resampler_models <- function(mod, .data, iter = 1) {
  # Remove missing values
  .data <- na.omit(.data)
  # Coefficients
  beta <- coef(mod)
  # mod$terms contains the formula used in the regression
  X <- model.matrix(mod$terms, data = .data)
  # estimate of std. dev. of errors
  sigma <- sqrt(sum(mod$residuals ^ 2) / mod$df.residual)
  # This produces the same result
  # sigma <- summary(mod1)$sigma  
  # Number of observations
  n <- nrow(X)
  # Name of dependent variable
  outcome_var_name <- all.vars(mod$terms)[1]
  # List to save results
  results <- vector(mode = "list", length = iter)
  for (i in seq_len(iter)) {
    # draw errors
    errors <- rnorm(n, mean = 0, sd = sigma)
    # create new outcome variable from errors
    y <- X %*% beta + errors
    # replace outcome variable
    .data[[outcome_var_name]] <- y
    # run regression
    newmod <- lm(mod$terms, data = .data)
    # Save model stats as a data frame to the list
    results[[i]] <- glimpse(newmod) %>% mutate(.iter = i)
  }
  # Convert the list of data frames to a single data frame by stacking the iterations
  bind_rows(results)
}
```

### Bootstrap example

Example of a single bootstrap replication. We draw N observations *with replacement* 
from the original data.

```r
nunn_bootstrapped <- bootstrap(nunn, 1)
```

To get bootstrap standard errors, we draw `m` replications, run the regression, 
and save the estimates. 

```r
bootstrap(nunn, 1024) %>%
  do(tidy(lm(trust_neighbors ~ exports, data = nunn)))
```

There are several ways to calculate standard errors from bootstraped replications.

- Calculate the standard error from these simulations by taking the standard deviation of the estimates.
- Calculate the confidence interval using the 2.5% and 97.5% quantiles in the replications

However, in the bootstrap, we should draw the bootstrap samples the same way the sample
was drawn from the population. Why might this not be the case in what we just did? 

## Multiple comparisons and F-test


```r
noise <- data.frame(matrix(rnorm(2100), nrow = 100, ncol = 21))
summary(lm(noise))
```

- How many variables have t-tests that are significant?
- Is the F-test significant? 
- Explain the difference
