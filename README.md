condensier
==========

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/condensier)](http://cran.r-project.org/package=condensier)
[![](http://cranlogs.r-pkg.org/badges/condensier)](http://cran.rstudio.com/web/packages/condensier/index.html)
[![Travis-CI Build Status](https://travis-ci.org/osofr/condensier.svg?branch=master)](https://travis-ci.org/osofr/condensier)
[![Coverage Status](https://coveralls.io/repos/osofr/condensier/badge.png?branch=master&service=github)](https://coveralls.io/r/osofr/condensier?branch=master)

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
    bin.method = "equal.mass",
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
