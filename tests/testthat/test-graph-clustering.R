test_that("component.cluster works", {
  # load test fixture
  edge.info <- read.csv(test_path("test.tn93.csv"))
  seqs <- ape::read.FASTA(test_path("test.fasta"))
  seq.info <- parse.headers(
    names(seqs), var.names=c('accn', 'coldate', 'subtype'),
    var.transformations=c(as.character, as.Date, as.factor))
  seq.info$colyear <- data.table::year(seq.info$coldate)
  which.new <- (seq.info$colyear >= 2012)
  obj <- read.edges(edge.info, seq.info, which.new)
  
  result <- component.cluster(obj, dist.thresh=0.06)
  expect_true(is.data.table(result))
  
  expected <- c(2, 1, 1, 2)  # (1,2), 3, 4, (5,6)
  expect_equal(result$Size, expected)
  expected <- c(1, 1, 0, 1)
  expect_equal(result$Growth, expected)
  expected <- list(
    c(as.Date("2008-11-04"), as.Date("2009-04-28")),
    as.Date("2011-02-01"),
    as.Date("2008-09-30"),
    c(as.Date("2010-03-02"), as.Date("2008-05-20"))
  )
  expect_equal(result$coldate, expected)
  
  # relaxing threshold to merge nodes 3 and 4
  result <- component.cluster(obj, dist.thresh=0.085)
  expect_equal(result$Size, c(2,2,2))
  expect_equal(result$Growth, c(1,1,1))
  result <- component.cluster(obj, dist.thresh=0.086)
  expect_equal(result$Size, c(4,2))
  expect_equal(result$Growth, c(2,1))
  
  # tightening threshold causes new case to drop out
  result <- component.cluster(obj, dist.thresh=0.05)
  expect_equal(result$Size, c(2, 1, 1, 2))
  expect_equal(result$Growth, c(1, 0, 0, 1))
  
  result <- component.cluster(obj, dist.thresh=0.04)
  expect_equal(result$Size, c(1, 1, 1, 1, 1, 1))
  expect_equal(result$Growth, c(0, 1, 0, 0, 1, 0))
})


test_that("issue 22 is resolved", {
  ## Generate weights using legacy (v1.0) code:
  
  testthat::source_test_helpers()
  iFile <- test_path("test2.tn93.csv")
  iG <- impTN93(iFile, minNS=1)
  subG <- dFilt(iG, 0.02)  # current test case
  
  # step through compAnalyze()
  
  # e is edge list filtered by threshold - only one in-edge per new case
  dMax <- max(subG$e$Distance)  # longest distance under threshold is 0.0199944 
  # used because we have not passed threshold parameter to this function?
  
  # f is unfiltered edge list, excluding edges to new cases and between cases 
  # from the same time period, AND limiting to one in-edge per newer case
  tTab <- table(subG$f$tMax)  # most recent time for every bipartite edge
  
  # v is node list derived from unfiltered edge list
  vTab <- table(subG$v$Time)  # sample times of all nodes
  
  # how many filtered edges associated with each sampling time?
  eTab <- sapply(as.numeric(names(vTab)), function(t){
    nrow(subset(subG$e, (t1==t|t2==t)))
  })
  names(eTab) <- names(vTab)
  
  ageD <- bind_rows(lapply(as.numeric(names(tTab)), function(t) {
    # for every sampling time (not including earliest)
    temp <- subset(subG$f, tMax==t)  # get associated bipartite edges
    dfs <- split(temp, temp$tDiff)  # split by start times
    
    # how many edges in bipartite graph meet distance threshold?
    Positive <- sapply(dfs, function(df){length(which(df$Distance<=dMax))})
    vTotal <- rep((vTab[[as.character(t)]]),length(dfs))  # FIXME: unused?
    tDiff <- as.numeric(names(Positive))
    oeDens <- sapply(tDiff, function(tD){
      # to calculate the average degree (out-edge density TO ANY TIME POINT) of 
      # known cases from the same time point
      oTime <- t-tD  # start time for this bipartite graph
      
      # return the number of filtered edges with a node from time oTime
      # divide by total number of nodes from oTime
      return(eTab[as.character(oTime)]/vTab[as.character(oTime)])
    })
    res <- data.frame(
      time=t,
      Positive=as.numeric(Positive), 
      vTotal=vTab[[as.character(t)]], 
      oeDens=as.numeric(oeDens), 
      tDiff)  # note this can yield multiple rows given tDiff
    return(res)
  }))
  
  mod <- glm(cbind(Positive, vTotal) ~ tDiff+oeDens, data=ageD, family='binomial')
  
  # predicted probability of an edge to a new case given time difference and 
  # mean out-degree for a node from that time point
  subG$v$Weight <- predict(
    mod, type='response',
    data.frame(
      tDiff = max(subG$v$Time)-subG$v$Time, 
      oeDens = as.numeric(eTab[as.character(subG$v$Time)] / 
                            vTab[as.character(subG$v$Time)])
    ))
  
  ## Now attempt to generate the same weights using current code
  
  tn93 <- read.csv(test_path("test2.tn93.csv"))
  seqs <- ape::read.FASTA(test_path("test2.fasta"))
  
  seq.info <- parse.headers(
    names(seqs), var.names=c('accn', 'coldate', 'subtype'),
    var.transformations=c(as.character, as.Date, as.factor))
  seq.info$colyear <- data.table::year(seq.info$coldate)
  
  which.new <- (seq.info$colyear >= 2012)
  obj <- read.edges(tn93, seq.info, which.new)
  times <- obj$seq.info$colyear
  
  obj$seq.info$weights <- fit.decay(obj, times, dist.thresh=0.02)
  
  ## COMPARE RESULTS
  idx <- match(subG$v$ID, obj$seq.info$accn)
  
  expect_equal(obj$seq.info$accn[idx], subG$v$ID)
  expect_equal(obj$seq.info$colyear[idx], subG$v$Time)
  expect_equal(obj$seq.info$weights[idx], subG$v$Weight, tolerance=1e-6)
  
  
  ## compare model fits
  # old method
  subG <- simGrow(subG)
  cPred <- subset(subG$v, Time<max(Time))[,c("Weight", "Cluster")]
  
  #Create two data frames from two predictive models, one based on absolute size (NULL) and our date-informed model
  df1 <- data.frame(Growth = as.numeric(subG$g), 
                    Pred = sapply(names(subG$c), function(x) { 
                      sum(subset(cPred, Cluster==as.numeric(x))$Weight) 
                    }))
  df2 <- data.frame(Growth = as.numeric(subG$g), 
                    Pred = as.numeric(subG$c) * (sum(as.numeric(subG$g)) /
                                                   sum(as.numeric(subG$c))))
  fit1 <- glm(Growth ~ Pred, data = df1, family = "poisson")
  fit2 <- glm(Growth ~ Pred, data = df2, family = "poisson")
  subG$gaic <- fit1$aic-fit2$aic
  
  # new method
  cluster.data <- component.cluster(obj, dist.thresh=0.02, setID=0)
  cluster.data$RangeID <- 0
  expect_equal(as.integer(table(cPred$Cluster)), cluster.data$Size)
  result <- sapply(cluster.data$weights, sum)
  expect_equal(sort(result), sort(df1$Pred), tolerance=1e-6)
  
  ptrans <- list("weights"=sum)
  pmods <- list(
    "NullModel"=function(x) glm(Growth~Size, data=x, family="poisson"),
    "AltModel"=function(x) glm(Growth~weights, data=x, family="poisson")
  )
  result <- fit.analysis(cluster.data, transforms=ptrans, models=pmods)
  expect_equal(subG$gaic, result$AltModel[[1]]$aic - result$NullModel[[1]]$aic)
})


