seq.info <- parse.headers(names(alignment.ex), var.names = c("ID", "CollectionDate", "Subtype"),
                         var.transformations = list(as.character, as.Date, as.factor))
