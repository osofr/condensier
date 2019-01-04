### --- Test setup ---
if(FALSE) {
  # # CHECK AND BUILD PACKAGE:
  # # library("condensier")
  # library("roxygen2")
  # library("devtools")
  # library("testthat")
  # library("data.table")
  # setwd(".."); setwd(".."); getwd()
  # document()
  # load_all("./", create = FALSE) # load all R files in /R and datasets in /data. Ignores NAMESPACE:
  # # condensier:::debug_set() # SET TO DEBUG MODE
  # # devtools::use_mit_license()
  
  # setwd("..");
  # install("condensier", build_vignettes = FALSE, dependencies = FALSE) # INSTALL W/ devtools:

  # # system("echo $PATH") # see the current path env var
  # # system("R CMD Rd2pdf condensier")  # just create the pdf manual from help files

  # getwd()
  # # setwd("./condensier"); setwd(".."); getwd()
  # devtools::check(args = "--as-cran")
  # devtools::check(args = c("--no-vignettes"), build_args = c("--no-build-vignettes")) # runs faster
  # # devtools::check() # runs check with devtools
  # # devtools::check(args = c("--no-vignettes"), build_args = c("--no-build-vignettes")) # runs check with devtools
  # # devtools::build_win(args = "--compact-vignettes") # build package on CRAN servers (windows os?)
  # # devtools::build()
  # devtools::build(args = "--compact-vignettes")
  # # devtools::build(args = c("--compact-vignettes", "--resave-data"))
  # # devtools::build(args = "--no-build-vignettes") # build package tarball compacting vignettes
  # # devtools::build() # build package tarball
  # devtools::check_rhub()

  # setwd("..")
  # system("R CMD check --as-cran condensier_0.1.0.tar.gz") # check R package tar ball prior to CRAN submission
  #     ## system("R CMD check --no-manual --no-vignettes condensier") # check without building the pdf manual and not building vignettes
  #     ## system("R CMD build condensier --no-build-vignettes --as-cran")
  #     # system("R CMD build condensier --resave-data")
  # # devtools::use_travis() # SET UP TRAVIS CONFIG FILE
  # # INSTALLING FROM SOURCE:
  # # install.packages("./condensier_0.2.0.tar.gz", repos = NULL, type="source", dependencies=TRUE)
  # # library(condensier)
  # # condensier:::debug_set() # SET TO DEBUG MODE
  # # condensier:::debug_off() # SET DEBUG MODE OFF

  # # To install a specific branch:
  # # devtools::install_github('osofr/simcausal', ref = "simnet", build_vignettes = FALSE)
  # # devtools::install_github('osofr/condensier', ref = "master", build_vignettes = FALSE)

  # # TEST COVERATE:
  # # if your working directory is in the packages base directory
  # # package_coverage()
  # # or a package in another directory
  # # cov <- package_coverage("condensier")
  # # view results as a data.frame
  # # as.data.frame(cov)
  # # zero_coverage() can be used to filter only uncovered lines.
  # # zero_coverage(cov)
}


