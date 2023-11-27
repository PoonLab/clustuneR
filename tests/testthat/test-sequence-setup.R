# excerpt from NA dataset (Genbank)
headers <- c("KU190029.1_2013-09-10_subA", "KU190030.1_2013-11-21_subA", 
             "KU190730.1_2013-01-02_subB", "KU190731.1_2013-01-02_subB", 
             "KU190732.1_2013-01-04_subB", "KU190733.1_2013-01-08_subB")

test_that("parse headers", {
  result <- parse.headers(headers, var.names=c("accn", "coldate", "subtype"),
                          var.transformations=c(as.character, as.Date, as.factor),
                          sep="_")
  expected <- data.table(
    Header=headers,
    accn=c("KU190029.1", "KU190030.1", "KU190730.1", 
           "KU190731.1", "KU190732.1", "KU190733.1"),
    coldate=c(as.Date("2013-09-10"), as.Date("2013-11-21"), 
              as.Date("2013-01-02"), as.Date("2013-01-02"),
              as.Date("2013-01-04"), as.Date("2013-01-08")),
    subtype=as.factor(c("subA", "subA", "subB", "subB", "subB", "subB"))
  )
  expect_equal(result, expected)
})

test_that("Empty var.names raises warning, default values", {
  fcall <- function() parse.headers(
    headers, var.transformations=c(as.character, as.Date, as.factor), sep='_'
  )
  expect_warning(result <- fcall())
  expected <- c("Header", "V1", "V2", "V3")
  expect_equal(names(result), expected)
})

test_that("var.transformations defaults to character values", {
  expect_warning(
    result <- parse.headers(
      headers, var.names=c("accn", "coldate", "subtype"), sep="_")
  )
  expected <- c("2013-09-10", "2013-11-21", "2013-01-02", "2013-01-02",
                "2013-01-04", "2013-01-08")
  expect_equal(result$coldate, expected)
})
