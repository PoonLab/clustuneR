test_that("multi.cluster works for components", {
  # load test fixture
  edge.info <- read.csv(test_path("test.tn93.csv"))
  seqs <- ape::read.FASTA(test_path("test.fasta"))
  seq.info <- parse.headers(
    names(seqs), var.names=c('accn', 'coldate', 'subtype'),
    var.transformations=c(as.character, as.Date, as.factor))
  seq.info$colyear <- data.table::year(seq.info$coldate)
  which.new <- (seq.info$colyear >= 2012)
  obj <- read.edges(edge.info, seq.info, which.new)
  
  param.list <- lapply(c(0.04, 0.05, 0.06, 0.085, 0.086), function(x) 
    list(dist.thresh=x))
  result <- multi.cluster(obj, param.list, component.cluster)
  expect_true(is.data.table(result))
  
  # check SetID
  expected <- c(rep(1,6), rep(2,4), rep(3,4), rep(4,3), rep(5,2))
  expect_equal(result$SetID, expected)
  expected <- c(rep(0.04,6), rep(0.05,4), rep(0.06,4), rep(0.085,3), 
                rep(0.086,2))
  expect_equal(result$DistThresh, expected)
  
  # check cluster sizes
  expected <- c(rep(1,6), c(2,1,1,2), c(2,1,1,2), c(2,2,2), c(4,2))
  expect_equal(result$Size, expected)
  expected <- c(c(0,1,0,0,1,0), c(1,0,0,1), c(1,1,0,1), c(1,1,1), c(2,1))
  expect_equal(result$Growth, expected)
  
  # TODO: test that this works with reduced edge list
  
  # check that function rejects incompoatible parameter list
  expect_error(multi.cluster(obj, param.list, step.cluster))
})


test_that("fit.analysis works", {
  set.seed(1)
  lo.risk <- rpois(100, 0.1)
  hi.risk <- rpois(100, 0.2)
  
  # partition into clusters
  sizes <- c(50, 20, 10, 5, 5, 2, 2, 2, 1, 1, 1, 1)
  temp <- sapply(cumsum(sizes), function(i) sum(lo.risk[1:i]))
  lo.grow <- c(temp[1], diff(temp))
  temp <- sapply(cumsum(sizes), function(i) sum(hi.risk[1:i]))
  hi.grow <- c(temp[1], diff(temp))
  
  cluster.data <- data.table(
    SetID = 1,
    RangeID=1,
    Header = NA,
    Size = rep(sizes, 2),
    Growth = c(lo.grow, hi.grow),
    Predictor = rep(c(0, 1), each=12)  # collapse
  )
  fit0 <- glm(Growth ~ Size, data=cluster.data, family='poisson')
  fit1 <- glm(Growth ~ Size + Predictor, data=cluster.data, family='poisson')

  ptrans <- list("Predictor"=identity)  
  pmods <- list(
    "NullModel"=function(x) glm(Growth~Size, data=x, family="poisson"),
    "AltModel"=function(x) glm(Growth~Size+Predictor, data=x, family="poisson")
    )
  result <- fit.analysis(cluster.data, transforms=ptrans, models=pmods)
  
  expect_equal(result$NullModel[[1]]$aic, fit0$aic)
  expect_equal(result$AltModel[[1]]$aic, fit1$aic)
  
  # scrambling predictor values among clusters should reduce dAIC
  cluster.data$Predictor <- 
    lapply(cluster.data$Size, function(x) sample(c(0,1), size=x, replace=T))
  ptrans <- list("Predictor"=mean)
  result <- fit.analysis(cluster.data, transforms=ptrans, models=pmods)
  
})
