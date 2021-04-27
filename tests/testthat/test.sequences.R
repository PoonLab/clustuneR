test_that("Sequence Headers will be overwritten by full headers if named by User", {
  expect_warning(x <- pull.headers(alignment.ex,var.names = c("Header", "CollectionDate", "Subtype"),
                              var.transformations =list(as.character, as.Date, as.factor)),
                 "'Header' is contained within var.names, this will be overwritten")
  expect_true(all(x$Header%in%names(alignment.ex)))
})



