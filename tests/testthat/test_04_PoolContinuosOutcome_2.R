library(mockery)
context("BinDat Pooled Continuous Outcome")


# ------------------------------------------------------------------------------------------------
# Test for pooled fitting of the bin indicators
# AN IDEA FOR TESTING pooled regression:
# USE IT TO ESTIMATE POOLED IPTW FOR LONGITUDINAL DATA WITH SEVERAL TIME POINTS (RN's simulation)
# ------------------------------------------------------------------------------------------------
test.PoolContRegression <- function() {
  # require(data.table)
  # gvars <- condensier:::gvars
  # reg_test <- RegressionClass$new(outvar.class = c(gvars$sVartypes$cont),
  #                                 outvar = c("sA"),
  #                                 predvars = c("W1", "W2", "W3"),
  #                                 # subset = list(quote(TRUE)))
  #                                 subset = list("sA"),
  #                                 pool = TRUE, nbins = 10, bin_bymass = FALSE)
  # datO <- get.testDat(nsamp = 1000)
  # nodeobjs <- get.testDatNet(datO)
  # datNetObs <- nodeobjs$datNetObs
  # class(datNetObs) # [1] "DataStore" "DatNet"      "R6"
  # model3 <- SummariesModel$new(reg = reg_test, DataStore.g0 = nodeobjs$datNetObs)
  # # Matrix of all summary measures: (sW,sA)
  # head(nodeobjs$datNetObs$mat.sVar); class(nodeobjs$datNetObs$mat.sVar)
  # head(datNetObs$mat.bin.sVar)
  # binfit_time <- system.time(
  #   model3$fit(data = nodeobjs$datNetObs)
  #   # Error in 1L:self$nbins : argument of length 0
  # )
  # binfit_time
  # binpredict_time <- system.time(
  #   probAeqa <- model3$predictAeqa(newdata = nodeobjs$datNetObs)
  # )
  # binpredict_time
  # [1] "fit (10K)"
  # $coef
  #  Intercept     bin_ID         W1         W2         W3
  # -2.7756215  0.1553186 -1.0014477 -0.5720651 -0.3339728
  # [1] "res_DT: "
  #            ID ProbAeqa_long
  #      1:     1    0.97396036
  #      2:     1    0.96971742
  #      3:     1    0.96480811
  #      4:     1    0.95913645
  #      5:     1    0.95259565
  #     ---
  # 104496:  9998    0.07668105
  # 104497:  9999    0.93215687
  # 104498:  9999    0.92165035
  # 104499:  9999    0.09032560
  # 104500: 10000    0.06784313
  # [1] "res_DT_short: "
  #           ID    cumprob
  #     1:     1 0.06099655
  #     2:     2 0.06145986
  #     3:     3 0.03836225
  #     4:     4 0.05821479
  #     5:     5 0.07303417
  #    ---
  #  9996:  9996 0.05119563
  #  9997:  9997 0.05896735
  #  9998:  9998 0.06414013
  #  9999:  9999 0.07760077
  # 10000: 10000 0.06784313
  # [1] "head(ProbAeqa, 50)"
  #  [1] 0.060996548 0.061459862 0.038362248 0.058214786 0.073034166 0.064140127 0.060658764 0.050023002 0.026039639 0.075325033 0.029168620
  # [12] 0.054538219 0.054031618 0.083549453 0.008653412 0.029594466 0.077600772 0.081220201 0.068319822 0.061459862 0.071357407 0.039453938
  # [23] 0.075325033 0.039007914 0.057871503 0.077600772 0.057871503 0.058967354 0.064140127 0.043973691 0.046655735 0.079794387 0.074434114
  # [34] 0.058967354 0.067843133 0.063492979 0.033237556 0.064138704 0.056974041 0.065426910 0.037236039 0.029168620 0.056974041 0.047226347
  # [45] 0.043973691 0.084256432 0.060173071 0.073034166 0.029168620 0.060183301
}
