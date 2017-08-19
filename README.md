condensier
==========

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/condensier)](http://cran.r-project.org/package=condensier)
[![](http://cranlogs.r-pkg.org/badges/condensier)](http://cran.rstudio.com/web/packages/condensier/index.html)
[![Travis-CI Build Status](https://travis-ci.org/osofr/condensier.svg?branch=master)](https://travis-ci.org/osofr/condensier)
[![Build status](https://ci.appveyor.com/api/projects/status/8c2xg9pgsohappsu/branch/master?svg=true)](https://ci.appveyor.com/project/osofr/condensier/branch/master)
[![Coverage Status](https://img.shields.io/codecov/c/github/osofr/condensier/master.svg)](https://codecov.io/github/osofr/condensier?branch=master)

### Installation

To install the development version of `condensier` (requires the `devtools` package):

```R
devtools::install_github('osofr/condensier', build_vignettes = FALSE)
```


### Instructions

Simulate some data with continuous outcome ("sA"):

```R
library("simcausal")
D <- DAG.empty()
D <-
D + node("W1", distr = "rbern", prob = 0.5) +
  node("W2", distr = "rbern", prob = 0.3) +
  node("W3", distr = "rbern", prob = 0.3) +
  node("sA.mu", distr = "rconst", const = (0.98 * W1 + 0.58 * W2 + 0.33 * W3)) +
  node("sA", distr = "rnorm", mean = sA.mu, sd = 1)
D <- set.DAG(D, n.test = 10)
datO <- sim(D, n = 10000, rndseed = 12345)
```

Fit conditional density using equal mass bins (same number of obs per bin):

```R
library("condensier")
dens_fit <- fit_density(
    X = c("W1", "W2", "W3"), 
    Y = "sA", 
    input_data = datO, 
    nbins = 20, 
    bin_method = "equal.mass",
    bin_estimator = speedglmR6$new())
```

Wrapper function to predict the conditional probability (likelihood) for new observations:

```R
newdata <- datO[1:5, c("W1", "W2", "W3", "sA"), with = FALSE]
preds <- predict_probability(dens_fit, newdata)
```

Wrapper function to sample the values from the conditional density fit:
```R
sampledY <- sample_value(dens_fit, newdata)
```


### Fitting Super Learner density with `sl3` package

Any binary-outcome regression learner available in `sl3` package can be used as a "drop-in" learner for conditional bin hazard. Below, we use `xgboost` R package implementation of GMB for estimation of the bin hazards (hazards are pooled across all bins into a single dataset, a single regression fit is then performed across all bins).

```R
library(sl3)
lrn <- Lrnr_condensier$new(task, nbins = 35, bin_method = "equal.len", pool = TRUE, bin_estimator = 
  Lrnr_xgboost$new(nrounds = 50, objective = "reg:logistic"))
```

Finally, multiple candidate density estimators can be optimally stacked or combined with a Super Learner. The convex combination of the candidates is found by minimizing the cross-validated negative loglikelihood loss function. In this example we define 3 candidate density learners:

```R
lrn1 <- Lrnr_condensier$new(task, nbins = 25, bin_method = "equal.len", pool = TRUE, 
  bin_estimator = Lrnr_glm_fast$new(family = "binomial"))
lrn2 <- Lrnr_condensier$new(task, nbins = 20, bin_method = "equal.mass", pool = TRUE,
  bin_estimator = Lrnr_xgboost$new(nrounds = 50, objective = "reg:logistic"))
lrn3 <- Lrnr_condensier$new(task, nbins = 35, bin_method = "equal.len", pool = TRUE,
  bin_estimator = Lrnr_xgboost$new(nrounds = 50, objective = "reg:logistic"))
```

We proceed by training the Super Learner (with 10 fold cross-validation) and then finding the optimal convex combination of the candidate densities with the meta-learner `Lrnr_solnp_density`:
```R
sl <- Lrnr_sl$new(learners = list(lrn1, lrn2, lrn3),
                  metalearner = Lrnr_solnp_density$new())
sl_fit <- sl$train(task)
```

To predict for new data, wrap the desired dataset into an `sl3-task` object and call predict on above `sl_fit` object:
```R
newdata <- datO[1:5, c("W1", "W2", "W3", "sA")]
new_task <- sl3_Task$new(newdata, covariates=c("W1", "W2", "W3"),outcome="sA" )
sl_fit$predict(new_task)
```


### Nesting the Super Learner for bin hazards with density Super Learner

Note that `bin_estimator` can be also a Super-Learner object from `sl3`. In this case the bin hazard will be estimated by stacking several candidate estimators. For example, below, we define a single density learner `lrn`,  with the hazard estimator defined by the Super-Learner that stacks two candidates (GLM and `xgboost` GBM). Note that in contrast to the above example, this Super-Learner fit will be optimized for the logistic regression problem (estimating pooled bin hazards), but still using internal 10-fold cross-validation. 

```R
library(sl3)
lrn <- Lrnr_condensier$new(task, nbins = 35, bin_method = "equal.len", pool = TRUE, bin_estimator = 
  Lrnr_sl$new(
    learners = list(
      Lrnr_glm_fast$new(family = "binomial"),
      Lrnr_xgboost$new(nrounds = 50, objective = "reg:logistic")
      ),
    metalearner = Lrnr_glm$new()
    ))
binSL_fit <- lrn$train(task)
```

In prinicple, one can nest the two of the above described types of Super Learners: the Super Learner that fits the bin hazard of each candidate density and the Super Learner that finds the optimal combination of the candidate densities. However, due to potential performance constraints, we currently advise against that. 

### Stacking and cross-validating candidate densities with `sl3` package

One can build a custom version of their own Super Learner by using the stacking and cross-validation procedures availabe in `sl3`. Here we define a stack of 3 learners, then train all 3 and predict for new data (likelihood):
```R
learner_stack <- Stack$new(lrn1, lrn2, lrn3)
stack_fit <- learner_stack$train(task)
preds <- stack_fit$predict(new_task)
```

Here we cross-validate all 3 learners in the stack, using the default 10-fold CV:
```R
cv_stack <- Lrnr_cv$new(learner_stack)
cv_fit <- cv_stack$train(task)
```


### Funding
The development of this package was funded through an NIH grant (R01 AI074345-07).

### Copyright
The contents of this repository are distributed under the MIT license.
```
The MIT License (MIT)

Copyright (c) 2017 Oleg Sofrygin 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
