require(phytools)
require(ape)
require(digest)

#' Prepare a tree for clustering
#'
#' Extends an ape tree object to include annotated node and path information.
#' By default this will include some additional info like paths, bootstraps and 
#' a table of cluster growth. This function makes calls to numerous helpers,
#' including the binaries of pplacer and guppy software installed with this package
#' (Matsen, 2010)
#'
#' @param phy:  ape::phylo object.  Un-extended.
#' @param seq.info: A data frame or data.table object containing the sequences.
#' New sequences are assigned as those not in the tree.
#' @param log.file: A path to the logfile from a tree construction run
#' @param seqs.full: The full alignment. Including sequences excluded from the tree.
#' @return The tree annotated with node information and seq.info
#' @export
#' @example examples/extend.tree_ex.R
extend.tree <- function(phy, seq.info=data.frame(), full.align=character(0), 
                        log.file=NA, locus = "LOCUS") {

  # Midpoint root for consistency and resolve multichotomies
  phy <- phytools::midpoint.root(phy)
  phy <- ape::multi2di(phy)

  # Check Sequence names inputs
  if (nrow(seq.info) == 0){
    if(length(full.align) != 0){
      warning("No sequence meta-data included, creating default seq.info input",
              "from headers in alignment")
      seq.info <- data.frame("Header"=names(full.align)) 
    } else {
      warning("No sequence meta-data included, creating default seq.info input", 
              "from tree tip.labels")
      seq.info <- data.frame("Header"=phy$tip.label) 
    }
  }
  
  var.names <- colnames(seq.info)
  if (!("Header" %in% var.names)) {
    stop("Data frame `var.name` must contain Header")
  } 
  else {
    if (length(unique(seq.info[, "Header"])) != length(seq.info[, "Header"])) {
      seq.info <- seq.info[!duplicated(Header), ]
      warning("Duplicate Headers have been arbitrarily removed")
    }
    if (!all(phy$tip.label %in% seq.info$Header)) {
      stop("At least 1 tip labels in tree is not represent seq.labels")
    }
  }
  
  # label new sequences
  seq.info$New <- FALSE
  which.new <- which(!(seq.info$Header %in% phy$tip.label))
  seq.info$New[which.new] <- TRUE

  phy$seq.info <- seq.info

  phy$node.info <- annotate.nodes(phy, mc.cores)
  phy$path.info <- annotate.paths(phy)

  # Extend with growth_info
  if(is.na(log.file) | (length(full.align)==0)) {
    warning("Ignoring growth information, path to logfile and full sequence ", 
            "alignment required.")
    phy$growth.info <- data.frame(
      "Header"=character(0), "NeighbourDes"=numeric(0), "Bootstrap"=numeric(0),
      "TermDistance"=numeric(0), "PendantDistance"=numeric(0), "Terminal"=logical(0)
    )
  } 
  else {
    if(!all(names(full.align) %in% seq.info$Header)){
      stop("Headers in full sequence alignment do not correspond to seq.info data")
    }

    stats.json <- translate.log(log.file)
    refpkg <- taxit.create(phy, full.align, stats.json)
    growth.info.trees  <- run.pplacer_guppy(refpkg)
    
    phy$growth.info <- annotate.growth(phy, growth.info.trees)
  }

  return(phy)
}


#' Prepare nodes in tree for clustering
#'
#' Groups meta data by descendant. This also includes bootstrap certainty referenced
#' by node labels in the ape object. Unlabeled nodes are defaulted to a bootstrap 
#' value of 1 (maximum certainty), including the root.
#'
#' @param phy: An inputted tree passed by extend.tree. Should contain seq.info
#' @param max.boot: numeric, user specifies maximum bootstrap value (usually 
#'                  1 or 100).
#' @return:  data frame, "node.info" to a annotate the information under a given
#' node at a tree.
annotate.nodes <- function(phy, max.boot=NA) {

  # We assume that internal node labels represent bootstrap support values
  # Assign maximum bootstrap value to unlabelled internal nodes
  phy$node.label <- as.numeric(phy$node.label)
  
  # are bootstraps scaled to 1 or 100?
  max.boot <- 1
  if (sum(phy$node.label > 1, na.rm=T) / phy$Nnode > 0.5) max.boot <- 100
  
  # assign max value to nodes without bootstrap support
  phy$node.label[which(is.na(phy$node.label))] <- max.boot

  # number of nodes in a rooted binary tree is 2n-1
  # note tips are numbered first (1..n), followed by internal nodes
  node.info <- data.frame(
    NodeID = 1:(2*Ntip(phy)-1),
    Bootstrap = c(rep(max.boot, Ntip(phy)), phy$node.label)
    )
  
  # FIXME: why are we normalizing bootstrap values?
  node.info$Bootstrap <- node.info$Bootstrap / max(node.info$Bootstrap)
  
  # Get descendant information for each node
  des <- phangorn::Descendants(phy, type="all")
  node.info$Descendants <- des

  # Get pairwise patristic distance info
  patristic.dists <- ape::dist.nodes(phy)
  dists.by.des <- lapply(des, function(d) {
    idx <- d[d<=Ntip(phy)]
    patristic.dists[idx, idx]
  })
  node.info$max.patristic.dist <- sapply(dists.by.des, max)
  node.info$mean.patristic.dist <- sapply(
    dists.by.des, function(x) { 
      xx <- x[x > 0]
      if (length(xx) > 0) { return (mean(xx)) }
      else { return (NA) }
      }
    )

  return(node.info)
}

#' Get paths relative to starting nodes
#'
#' Called by extend.tree. Adds path info to the tree, ie. The path of branches from
#' each node to the root. Required for step.cluster()
#'
#' @param phy:  ape::phylo object. An inputted tree An inputted tree passed by 
#'              extend.tree.
#' @return A matrix labelled "path.info" to attach to a given tree. For each node 
#' in the path the branch lengths (below node) and bootstraps are given. For terminal 
#' nodes, no branch length is given below the node and the bootstrap is 1
annotate.paths <- function(phy) {

  # Get paths and length information from terminal nodes
  lens <- ape::node.depth.edgelength(phy)
  paths <- ape::nodepath(phy)
  names(paths) <- sapply(paths, function(p) {
    p[length(p)]
  })

  # Extend apes nodepath() function to internal nodes
  i <- 1
  while (length(paths) < nrow(phy$node.info) + 1) {
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
    c(NA, (lens[x[-1]] - lens[x[-length(x)]]))
  })

  # Obtain bootstrap branch length and node number information for all paths
  boots <- lapply(paths, function(x) {
    phy$node.info$Bootstrap[x]
  })
  
  path.info <- lapply(1:nrow(phy$node.info), function(i) {
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
#' This uses pplacer software to ensure that previously defined clusters remain 
#' the same. New tips are assigned a set of possible placements in the tree.
#'
#' @param t: An inputted tree using ape's tree handling passed by extend.tree
#' @param t.growth: A set of trees from pplacer. Each with a newly added tip.
#' @return The growth.info for a given tree. This includes the branch length to the 
#' newly added terminal node, whether or not the neighbour to the newly added terminal
#' node is also terminal, and the branch length to that neighbour, from the newly 
#' placed internal node (the "PendantLength")
annotate.growth <- function(phy, phy.grown) {

  growth.info <- lapply(phy.grown, function(x) {
    # each x is a tree with a new tip grafted
    
    # extract set of possible placements of new tip
    new.ids <- setdiff(x$tip.label, phy$tip.label)
    new.tips <- which(x$tip.label %in% new.ids)
    
    # extract tip label prefix (without M value from pplacer)
    prefix <- gsub("_#.*", "", new.ids[1])
    
    # retrieve indices of in-edges to new tips
    new.edges <- which(x$edge[, 2] %in% new.tips)
    term.dists <- x$edge.length[new.edges]
    
    new.nodes <- x$edge[new.edges, 1]
    new.nodes.des <- which(x$edge[, 1] %in% new.nodes)
    node.boots <- sapply(new.ids, function(x) {
      l <- strsplit(x, "_M=")[[1]]
      as.numeric(l[length(l)])
    })
    
    # pendant edges 
    pen.edges <- setdiff(new.nodes.des, new.edges)
    pen.dists <- x$edge.length[pen.edges]
    
    neighbour.des <- lapply(new.nodes, function(n) {
      phangorn::Descendants(x, n, "tips")[[1]]
    })
    old.tips <- lapply(neighbour.des, function(des) {
      which(phy$tip.label %in% x$tip.label[des])
    })
    terminal <- sapply(old.tips, function(x) {
      length(x) == 1
    })
    
    list(
      Header=prefix,
      NeighborDes=old.tips,
      Bootstrap=node.boots,
      TermDistance=term.dists,
      PendantDistance=pen.dists,
      Terminal=terminal
    )
  })
  
  # Obtain placement information from trees.
  # New IDs, bootstap + branch lengths of new node

  # Collapse neighbour node descendant tips to their MRCA
  for (i in 1:length(growth.info)) {
    collapsed.neighbours <- sapply(
      which(!growth.info[[i]]$Terminal), function(j) {
        ape::getMRCA(phy, growth.info[[i]]$NeighborDes[[j]])
        })
    nd <- growth.info[[i]]$NeighborDes
    nd[which(!growth.info[[i]]$Terminal)] <- collapsed.neighbours
    suppressWarnings(growth.info[[i]]$NeighbourNode <- unlist(nd))
    growth.info[[i]]$NeighborDes <- NULL
  }
  
  # stretch out list to data frame
  df <- data.frame(
    Header=unlist(sapply(growth.info, 
                         function(x) rep(x$Header, length(x$Bootstrap))))
    )
  for (field in names(growth.info[[1]])) {
    if (field == "Header") next;
    df[[field]] <- unlist(sapply(growth.info, function(x) x[[field]]))
  }
  row.names(df) <- NULL

  if (!all(growth.info$Header %in% phy$seq.info$Header[phy$seq.info$New])) {
    # FIXME: sometimes tip labels are merged in pplacer output
    warning("Not all newly added sequences are noted in the seq.info of the tree")
  }

  return(growth.info)
}

#' Build a stats.json
#'
#' Parses the logfile output of tree building software to find information relevant 
#' to pplacer. Prints stats to a temporary ref.pkg file and returns the file path
#' NOTE: Currently compatible with GTR substitution model and either FastTree, 
#' RAxML or IQ-TREE logfiles. Working to extend this to PhyML logfiles, and to 
#' different models of evolution
#'
#' @param log.file: A path to the logfile from a tree construction run.
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

#' Create a refpkg
#'
#' A wrapper for the taxit.create function used by pplacer. This will generate a 
#' temporary refpkg file based on passed tree and alignment data summary json. 
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

#' obtain placement information to identify grown phylogenies
#'
#' A wrapper for pplacer's basic run function coupled with guppy's sing function.
#' Together, these extend fixed trees with most likely placement locations.
#' 
#' @param refpkg: A reference package to use as input for pplacer
#' @return A set of trees, each containing 1 new sequence as a placement
run.pplacer_guppy <- function(refpkg, quiet=FALSE) {
  if (!quiet) {
    print("Running pplacer...")
  }
  pplacer <- system.file("pplacer", package="clustuneR")
  if (pplacer == "") {
    # in case user has not installed package, but running from pkg directory
    pplacer <- "inst/pplacer"
  }
  guppy <- system.file("guppy", package="clustuneR")
  if (guppy == "") {
    guppy <- "inst/guppy"
  }
  
  # Run pplacer to obtain placements
  system(paste0(
    "export LC_ALL=C; ", pplacer, " -c ", refpkg,
    " -o ", refpkg, "/placements.jplace",
    " --verbosity 0 ",
    refpkg, "/seq.fasta"
  ))

  # Run guppy to obtain trees
  system(paste0(
    guppy, " sing ", refpkg, "/placements.jplace",
    " -o ", refpkg, "/growth.tre"
  ))

  ts <- ape::read.tree(paste0(refpkg, "/growth.tre"))

  return(ts)
}
