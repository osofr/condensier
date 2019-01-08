---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# R/`condensier`: Non-parametric Multivariate Conditional Density Estimation with Binned Histograms


[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/condensier)](https://CRAN.R-project.org/package=condensier)
[![](https://cranlogs.r-pkg.org/badges/condensier)](https://CRAN.R-project.org/package=condensier)
[![Travis-CI Build Status](https://travis-ci.org/osofr/condensier.svg?branch=master)](https://travis-ci.org/osofr/condensier)
[![codecov](https://codecov.io/gh/osofr/condensier/branch/master/graph/badge.svg)](https://codecov.io/gh/osofr/condensier)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

> Fit a conditional density `f(A|W)`, where `A` can be continuous and multivariate and `W` is set of predictors.
> This estimator breaks up the support of a continuous `A` into discrete bins and fits the conditional hazard for each bin. By default the logistic regression will be used for fitting each bin hazard. Alternatively, arbitrary machine learning algorithms can be used via learners available in `sl3` R package (see example below).
> Given several competing candidate density estimators, one can find the optimal convex combination these candidate estimators by using Super Learner [`sl3`].
> For detailed description of the estimator implemented in this package see <a name=cite-diaz2011super></a>([Díaz Muñoz and van der
Laan, 2011](#bib-diaz2011super)) and <a name=cite-munoz2012population></a>([Muñoz and van der
Laan, 2012](#bib-munoz2012population)).

__Authors:__ [Oleg Sofrygin](https://github.com/osofr), [Frank Blaauw](https://github.com/frbl), [Antoine Chambaz](https://github.com/achambaz), Mark van der Laan


### Installation

To install the development version of `condensier` (requires the `devtools` package):


```r
devtools::install_github('osofr/condensier', build_vignettes = FALSE)
```


### Instructions

Simulate some data with continuous outcome (`"sA"`):


```r
library("simcausal")
D <- DAG.empty()
D <-
D + node("W1", distr = "rbern", prob = 0.5) +
  node("W2", distr = "rbern", prob = 0.3) +
  node("W3", distr = "rbern", prob = 0.3) +
  node("sA.mu", distr = "rconst", const = (0.98 * W1 + 0.58 * W2 + 0.33 * W3)) +
  node("sA", distr = "rnorm", mean = sA.mu, sd = 1)
D <- set.DAG(D, n.test = 10)
#> ...automatically assigning order attribute to some nodes...
#> node W1, order:1
#> node W2, order:2
#> node W3, order:3
#> node sA.mu, order:4
#> node sA, order:5
datO <- sim(D, n = 10000, rndseed = 12345)
#> simulating observed dataset from the DAG object
```

Fit conditional density using equal mass bins (same number of observations per bin):


```r
library("condensier")
#> condensier
#> The condensier package is still in beta testing. Interpret results with caution.
dens_fit <- fit_density(
    X = c("W1", "W2", "W3"), 
    Y = "sA", 
    input_data = datO, 
    nbins = 20, 
    bin_method = "equal.mass",
    bin_estimator = speedglmR6$new())
```

Wrapper function to predict the conditional probability (likelihood) for new observations:


```r
newdata <- datO[1:5, c("W1", "W2", "W3", "sA"), with = FALSE]
preds <- predict_probability(dens_fit, newdata)
```

Wrapper function to sample the values from the conditional density fit:


```r
sampledY <- sample_value(dens_fit, newdata)
```

Fit conditional density using custom bin definitions (argument `intrvls`):


```r
dens_fit <- fit_density(
    X = c("W1", "W2", "W3"),
    Y = "sA",
    input_data = datO,
    bin_estimator = speedglmR6$new(),
    intrvls = list(sA = seq(-4,4, by = 0.1)))
```

Fit conditional density using custom bin definitions and
pool all bin indicators into a single long-format dataset.
The pooling results in a single regression that is fit for all bin hazards,
with a bin indicator added as an additional covariate.


```r
dens_fit <- fit_density(
    X = c("W1", "W2", "W3"),
    Y = "sA",
    input_data = datO,
    bin_estimator = speedglmR6$new(),
    intrvls = list(sA = seq(-4,4, by = 0.1)),
    pool = TRUE)
```

### Fitting Super Learner density with `sl3` package

Any binary-outcome regression learner available in `sl3` package can be used as a "drop-in" learner for conditional bin hazard. Below, we use `xgboost` R package to define a new estimator of the bin hazard. Note that below, we are setting the tuning parameter `pool` to `TRUE`. This will have an effect of "pooling" all discrete bin indicators into a single dataset (with bin number added as a new covariate). This is followed by a single regression fit that is performed for all bins simultaneously (hence saving a lot of computation time and allowing the algorithm to perform smoothing over the bins).


```r
library("sl3")
#> Error in library("sl3"): there is no package called 'sl3'

task <- sl3_Task$new(datO, covariates=c("W1", "W2", "W3"), outcome="sA")
#> Error in eval(expr, envir, enclos): object 'sl3_Task' not found
lrn <- Lrnr_condensier$new(nbins = 10, bin_method = "equal.len", pool = TRUE, 
  bin_estimator = Lrnr_xgboost$new(nrounds = 5, objective = "reg:logistic"))
#> Error in eval(expr, envir, enclos): object 'Lrnr_condensier' not found

trained_lrn = lrn$train(task)
#> Error in eval(expr, envir, enclos): object 'lrn' not found

newdata <- datO[1:5, c("W1", "W2", "W3", "sA")]
new_task <- sl3_Task$new(newdata, covariates=c("W1", "W2", "W3"),outcome="sA" )
#> Error in eval(expr, envir, enclos): object 'sl3_Task' not found
pred_probs = trained_lrn$predict(new_task)
#> Error in eval(expr, envir, enclos): object 'trained_lrn' not found
pred_probs
#> Error in eval(expr, envir, enclos): object 'pred_probs' not found
```

Now that we have defined the candidate bin hazard estimator, it is time to train the model and obtained predictions (likelihood) based on new observations


```r
trained_lrn = lrn$train(task)
#> Error in eval(expr, envir, enclos): object 'lrn' not found

newdata <- datO[1:5, c("W1", "W2", "W3", "sA")]
new_task <- sl3_Task$new(newdata, covariates=c("W1", "W2", "W3"),outcome="sA" )
#> Error in eval(expr, envir, enclos): object 'sl3_Task' not found
pred_probs = trained_lrn$predict(new_task)
#> Error in eval(expr, envir, enclos): object 'trained_lrn' not found
pred_probs
#> Error in eval(expr, envir, enclos): object 'pred_probs' not found
```

Finally, multiple candidate density estimators can be optimally stacked or combined with a Super Learner. The convex combination of the candidates is found by minimizing the cross-validated negative loglikelihood loss function. In this example we define 3 candidate density learners:


```r
lrn1 <- Lrnr_condensier$new(nbins = 25, bin_method = "equal.len", pool = TRUE, 
  bin_estimator = Lrnr_glm_fast$new(family = "binomial"))
lrn2 <- Lrnr_condensier$new(nbins = 20, bin_method = "equal.mass", pool = TRUE,
  bin_estimator = Lrnr_xgboost$new(nrounds = 50, objective = "reg:logistic"))
lrn3 <- Lrnr_condensier$new(nbins = 35, bin_method = "equal.len", pool = TRUE,
  bin_estimator = Lrnr_xgboost$new(nrounds = 50, objective = "reg:logistic"))
```

We proceed by training the Super Learner (with 10 fold cross-validation) and then finding the optimal convex combination of the candidate densities with the meta-learner `Lrnr_solnp_density`:


```r
sl <- Lrnr_sl$new(learners = list(lrn1, lrn2, lrn3),
                  metalearner = Lrnr_solnp_density$new())
sl_fit <- sl$train(task)
```

To predict for new data, wrap the desired dataset into an `sl3-task` object and call predict on above `sl_fit` object:


```r
newdata <- datO[1:5, c("W1", "W2", "W3", "sA")]
new_task <- sl3_Task$new(newdata, covariates=c("W1", "W2", "W3"),outcome="sA" )
sl_fit$predict(new_task)
```


### Nesting the Super Learner for bin hazards with density Super Learner

Note that `bin_estimator` can be also a Super-Learner object from `sl3`. In this case the bin hazard will be estimated by stacking several candidate estimators. For example, below, we define a single density learner `lrn`,  with the hazard estimator defined by the Super-Learner that stacks two candidates (GLM and `xgboost` GBM). Note that in contrast to the above example, this Super-Learner fit will be optimized for the logistic regression problem (estimating pooled bin hazards), but still using internal 10-fold cross-validation. 


```r
library("sl3")
lrn <- Lrnr_condensier$new(nbins = 35, bin_method = "equal.len", pool = TRUE, bin_estimator = 
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


```r
learner_stack <- Stack$new(lrn1, lrn2, lrn3)
stack_fit <- learner_stack$train(task)
preds <- stack_fit$predict(new_task)
```

Here we cross-validate all 3 learners in the stack, using the default 10-fold CV:


```r
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

### References

<a name=bib-diaz2011super></a>[[1]](#cite-diaz2011super) I. Díaz
Muñoz and M. J. van der Laan. "Super learner based conditional
density estimation with application to marginal structural
models". In: _The international journal of biostatistics_ 7.1
(2011), pp. 1-20.

<a
name=bib-munoz2012population></a>[[2]](#cite-munoz2012population)
I. D. Muñoz and M. van der Laan. "Population intervention causal
effects based on stochastic interventions". In: _Biometrics_ 68.2
(2012), pp. 541-549.
