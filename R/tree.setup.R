#' Prepare a tree for clustering
#'
#' Extends an ape tree object to include annotated node and path information.
#' By default this will include some additional info like paths, bootstraps and growth.
#'
#' @param t: An inputted tree using ape's tree handling
#' @param seq.info: A data frame or data.table object containing the sequences.
#' By default, new sequences are assigned as those not in the tree.
#' @param log.file: A path to the logfile from a tree construction run
#' @param seqs.full: The full alignment. Including sequences excluded from the tree.
#' @param mc.cores: Passed to annotate.nodes as a parallel option
#' @return The tree annotated with node information and seq.info
extend.tree <- function(t, seq.info, mc.cores = 1, log.file=NA, full.align=character(0), locus = "LOCUS") {

  # Root the tree (if unrooted) and resolve multichotomies
  if (!ape::is.rooted(t)) {
    t <- phangorn::midpoint(t)
  }
  t <- ape::multi2di(t)

  # Check Sequence names inputs
  var.names <- colnames(seq.info)
  if (!("Header" %in% var.names)) {
    stop("Header must be contained within var.names")
  } else {
    if (length(unique(seq.info[, "Header"])) != length(seq.info[, "Header"])) {
      seq.info <- seq.info[!duplicated(Header), ]
      warning("Duplicate Headers have been arbitrarily removed")
    }
    if (!all(t$tip.label %in% seq.info$Header)) {
      stop("At least 1 tip labels in tree is not represent seq.labels")
    }
  }
  seq.info <- annotate.new(seq.info, which(!(seq.info$Header %in% t$tip.label)))
  t$seq.info <- seq.info

  t$node.info <- annotate.nodes(t, mc.cores)
  t$path.info <- annotate.paths(t)

  #Extend with growth_info
  if(is.na(log.file)|(length(full.align)==0)) {
    warning("Ignoring growth information, path to logfile and full sequence alignment required.")
    t$growth.info <- data.table::data.table(
      "Header" = character(0), "NeighbourDes" = numeric(0), "Bootstrap" = numeric(0),
      "TermDistance" = numeric(0), "PendantDistance" = numeric(0), "Terminal" = logical(0)
    )
  } else {

    if(!all(names(full.align)%in%seq.info$Header)){
      stop("Headers in full Sequence alignment do not correspond to seq.info data")
    }

    stats.json <- translate.log(log.file)
    refpkg <- taxit.create(t, full.align, stats.json)
    growth.info.trees  <- run.pplacer_guppy(refpkg)
    t$growth.info <- annotate.growth(t, growth.info.trees, mc.cores = mc.cores)
  }

  return(t)
}

##-TO-DO: Mask from user
#' Prepare nodes in tree for clustering
#'
#' Called by extend.tree. Adds additional node info to the tree. Required for mono.cluster()
#' NOTE: Unlabeled nodes are defaulted to a bootstrap value of 1
#'
#' @param t: An inputted tree using ape's tree handling. This must be annotated with seq.info
#' @param mc.cores: A parallel option
#' @return A slightly more detailed data table to replace "node.info" to a given tree.
annotate.nodes <- function(t, mc.cores = 1) {

  # Unlabelled nodes are defaulted to a bootstrap value of 1
  t$node.label[which(t$node.label %in% "")] <- 1

  nodes <- 1:(2 * length(t$tip.label) - 1)

  # Store node info in data.table
  node.info <- data.table::data.table()
  node.info[, "NodeID" := nodes]
  node.info[, "Bootstrap" := c(rep(100, length(t$tip.label)), as.numeric(t$node.label))]

  node.info[
    is.na(Bootstrap),
    "Bootstrap" := 10^ceiling(log10(max(node.info$Bootstrap[!is.na(node.info$Bootstrap)])))
  ]
  node.info[, "Bootstrap" := (node.info$Bootstrap) / max(node.info$Bootstrap)]

  # Get descendant information for each node
  des <- phangorn::Descendants(t, type = "all")
  node.info[, "Descendants" := des]

  # Get pairwise patristic distance info
  patristic.dists <- ape::dist.nodes(t)
  dists.by.des <- lapply(des, function(x) {
    patristic.dists[x[x <= length(t$tip.label)], x[x <= length(t$tip.label)]]
  })
  node.info$max.patristic.dist <- sapply(dists.by.des, function(x) {
    max(x)
  })
  node.info$mean.patristic.dist <- sapply(dists.by.des, function(x) {
    mean(x[x > 0])
  })

  return(node.info)
}

##-TO-DO: Mask from user
#' Get paths relative to starting nodes
#'
#' Called by extend.tree. Adds path info to the tree. Required for step.cluster()
#'
#' @param t: An inputted tree using ape's tree handling
#' @return A matrix labelled "path.info" to attach to a given tree.
#' For each node in the path the branch lengths (below node) and bootstraps are given
#' For the terminal node, no branch length is given below the node and the bootstrap is 1
annotate.paths <- function(t) {

  # Get paths and length information from terminal nodes
  lens <- ape::node.depth.edgelength(t)
  paths <- ape::nodepath(t)
  names(paths) <- sapply(paths, function(p) {
    p[length(p)]
  })

  # Extend apes nodepath() function to internal nodes
  i <- 1
  while (length(paths) < nrow(t$node.info) + 1) {
    new.paths <- sapply(paths[i:length(paths)], function(x) {
      x[-length(x)]
    })
    new.paths <- new.paths[!duplicated(new.paths)]
    names(new.paths) <- sapply(new.paths, function(p) {
      p[length(p)]
    })

    i <- length(paths) + 1
    paths <- c(paths, new.paths[which(!(names(new.paths) %in% names(paths)))])
  }
  paths <- paths[which(!sapply(paths, function(p) {
    length(p)
  }) == 0)]
  paths <- paths[order(as.numeric(names(paths)))]
  lens <- lapply(paths, function(x) {
    c((lens[x[-1]] - lens[x[-length(x)]]), NA)
  })

  # Obtain bootstrap branchlength and node number information for all paths
  boots <- lapply(paths, function(x) {
    t$node.info$Bootstrap[x]
  })
  path.info <- lapply(1:nrow(t$node.info), function(i) {
    m <- matrix(ncol = length(paths[[i]]), nrow = 3)
    rownames(m) <- c("Node", "Boot", "BranchLength")
    m[1, ] <- rev(paths[[i]])
    m[2, ] <- rev(boots[[i]])
    m[3, ] <- rev(lens[[i]])
    return(m)
  })

  return(path.info)
}


#' Add the growth information onto the known tree.
#'
#' This uses pplacer to ensure that previously defined clusters remain the same.
#' New tips are assigned a set of possible placements in the tree.
#'
#' @param t: An inputted tree using ape's tree handling
#' @param t.growth: A set of trees from pplacer. Each with a newly added tip.
#' @param mc.cores: A parallel option
#' @return The input tree annotated with growth information stored as growth.info.
annotate.growth <- function(t, t.grown, mc.cores = 1) {

  # Obtain placement information from trees.
  # New IDs, bootstap + branch lengths of new node
  growth.info <- dplyr::bind_rows(
    parallel::mclapply(t.grown, function(x) {
      new.ids <- setdiff(x$tip.label, t$tip.label)
      new.tips <- which(x$tip.label %in% new.ids)

      h <- gsub("_#.*", "", new.ids[1])

      new.edges <- which(x$edge[, 2] %in% new.tips)
      term.dists <- x$edge.length[new.edges]

      new.nodes <- x$edge[new.edges, 1]
      new.nodes.des <- which(x$edge[, 1] %in% new.nodes)
      node.boots <- sapply(new.ids, function(x) {
        l <- strsplit(x, "=")[[1]]
        as.numeric(l[length(l)])
      })

      pen.edges <- dplyr::setdiff(new.nodes.des, new.edges)
      pen.dists <- x$edge.length[pen.edges]

      neighbour.des <- lapply(new.nodes, function(n) {
        phangorn::Descendants(x, n, "tips")[[1]]
      })
      old.tips <- lapply(neighbour.des, function(des) {
        which(t$tip.label %in% x$tip.label[des])
      })
      terminal <- sapply(old.tips, function(x) {
        length(x) == 1
      })

      DT <- data.table::data.table(
        "Header" = h, "NeighbourDes" = old.tips, "Bootstrap" = node.boots,
        "TermDistance" = term.dists, "PendantDistance" = pen.dists, "Terminal" = terminal
      )

      return(DT)
    }, mc.cores = mc.cores)
  )

  # Collapse neighbour node descendant tips to their MRCA
  collapsed.neighbours <- sapply(growth.info[!(Terminal), NeighbourDes], function(des) {
    ape::getMRCA(t, des)
  })
  temp <- growth.info[, NeighbourDes]
  temp[which(!growth.info$Terminal)] <- collapsed.neighbours

  suppressWarnings(growth.info[, "NeighbourNode" := unlist(temp)])
  growth.info[, NeighbourDes := NULL]

  if (!all(growth.info$Header %in% t$seq.info[(New), Header])) {
    warning("Not all newly added sequences are noted in the seq.info of the tree")
  }

  return(growth.info)
}

##- TO-DO: Make this masked from user -##
##- TO-DO: Automatically identify program -##
#' Build a stats.json
#'
#' Parses the logfile output of tree building software to find information relevant to pplacer.
#' Prints stats to a temporary ref.pkg file and returns the path to that file.
#' NOTE: Currently compatible with GTR substitution model and either FastTree, RAxML or IQ-TREE logfiles.
#' Working to extend this to PhyML logfiles, and then to different models of evolution
#'
#' @param log.file: A path to the logfile from a tree construction run
#' @param substitution.model: The substitution model used. Currently only accepts "GTR"
#' @return A json output to be written to a given stats.file for pplacer
translate.log <- function(log.file) {

  # Open connection to log.file
  con <- file(log.file)
  lns <- readLines(con)
  close(con)

  #Identify program with tell-tale string
  if(any(grepl("FastTree", lns))){
    program <- "FastTree"
  } else if(any(grepl("RAxML", lns))) {
    program <- "RAxML"
  } else if(any(grepl("IQ-TREE", lns))) {
    program <- "IQ-TREE"
  } else{
    stop("Unrecognized log file. Note that currently only FastTree,
         IQ-TREE and RAxML logfiles are recognized")
  }

  # Extracts and normalizes list of frequencies
  if (program %in% "FastTree") {
    p <- "GTRRates"
    s <- lns[grep(p, lns)]
    s <- strsplit(s, "\t")[[1]]
    s <- as.numeric(s[c(2, 3, 4, 5, 6, 7)])
  }
  if (program %in% "IQ-TREE") {
    p <- "Rate parameters"
    s <- lns[grep(p, lns)]
    s <- strsplit(s[length(s)], " ")[[1]]
    s <- as.numeric(s[c(5, 20, 11, 8, 14, 17)])
  }
  if (program %in% "RAxML") {
    p <- " rates"
    s <- lns[grep(p, lns)]
    s <- strsplit(s[length(s)], " ")[[1]]
    s <- as.numeric(s[c(10, 15, 12, 11, 13, 14)])
  }

  # Write stats information to .json
  stats.json <- jsonlite::toJSON(list(
    "empirical_frequencies" = TRUE,
    "datatype" = "DNA",
    "subs_model" = "GTR",
    "program" = program,
    ## -TO-DO: Test Correctness of gamma assumption -##
    "ras_model" = "gamma",
    "gamma" = list(
      "alpha" = 1.0,
      "n_cats" = as.integer(20)
    ),
    "subs_rates" = list(
      "ac" = s[1],
      "gt" = s[2],
      "at" = s[3],
      "ag" = s[4],
      "cg" = s[5],
      "ct" = s[6]
    )
  ), pretty = T, always_decimal = T, auto_unbox = T)

  return(stats.json)
}

##- TO-DO: Less direct user interaction? -##
#' Create a refpkg
#'
#' A wrapper for the taxit create function used by pplacer. This will generate a summary json
#' See pplacer's basic function regarding alignments and tree function
#' NOTE: Creates temporary files. These are only deleted with the end of the session
#'
#' @param t: The tree (made on a subset of the full alignment)
#' @param seqs.full: The full alignment. Including sequences excluded from the tree.
#' @param stats.json: Path to the full alignment file (new + old seqs)
#' @param locus: Extra information required for the summary json
#' @return Path to a temporary refpkg directory
taxit.create <- function(t, seqs.full, stats.json, locus = "LOCUS") {

  # Set up and populate temporary file system
  refpkg <- tempdir()

  seq.file <- paste0(refpkg, "/seq.fasta")
  seq.file.name <- tail(strsplit(seq.file, "/")[[1]], 1)
  ape::write.FASTA(seqs.full, seq.file)

  stats.file <- paste0(refpkg, "/stats.json")
  stats.file.name <- tail(strsplit(stats.file, "/")[[1]], 1)
  conn <- file(stats.file)
  write(stats.json, conn)
  close(conn)

  tree.file <- paste0(refpkg, "/tree.nwk")
  tree.file.name <- tail(strsplit(tree.file, "/")[[1]], 1)
  ape::write.tree(t, tree.file)

  log.file <- paste0(refpkg, "/log.txt")
  log.file.name <- tail(strsplit(log.file, "/")[[1]], 1)
  conn <- file(log.file)
  write("Sample log file. Created Using taxit.create() wrapper",conn)
  close(conn)

  # Create generic JSON summary for refpkg.
  summary.json <- jsonlite::toJSON(list(
    "files" = list(
      "aln_fasta" = seq.file.name,
      "phylo_model" = stats.file.name,
      "tree" = tree.file.name,
      "tree_stats" = log.file.name
    ),
    "rollback" = NULL,
    "log" = c(
      "Stripped refpkg (removed 0 files)",
      "Loaded initial files into empty refpkg"
    ),
    "metadata" = list(
      "create_date" = as.character(Sys.Date()),
      "format_version" = "1.1",
      "locus" = locus
    ),
    "rollforward" = NULL,
    "md5" = list(
      "aln_fasta" = digest::digest(seq.file.name, algo = "md5"),
      "phylo_model" = digest::digest(stats.file.name, algo = "md5"),
      "tree" = digest::digest(tree.file.name, algo = "md5"),
      "tree_stats" = digest::digest(log.file.name, algo = "md5")
    )
  ), pretty = T, always_decimal = T, auto_unbox = T)
  summary.file <- paste0(refpkg, "/CONTENTS.json")
  con <- file(summary.file)
  write(summary.json, summary.file)
  close(con)

  return(refpkg)
}

##- TO-DO: Include package binaries such that it is not a requirement to install both pplacer and guppy -##
#' Get grown phylogenies
#'
#' A wrapper for pplacer's basic run function coupled with guppy's sing function.
#' Together, these extend fixed trees with most likely placement locations.
#' @param refpkg: A reference package to use as input for pplacer
#' @return A set of trees, each containing 1 new sequence.
run.pplacer_guppy <- function(refpkg) {


  # Run pplacer to obtain placements
  system(paste0(
    "export LC_ALL=C ; ", "pplacer -c ", refpkg,
    " -o ", refpkg, "/placements.jplace",
    " --verbosity 0",
    " ", refpkg, "/seq.fasta"
  ))

  # Run guppy to obtain trees
  system(paste0(
    "guppy sing ", refpkg, "/placements.jplace",
    " -o ", refpkg, "/growth.tre"
  ))

  ts <- ape::read.tree(paste0(refpkg, "/growth.tre"))

  return(ts)
}
