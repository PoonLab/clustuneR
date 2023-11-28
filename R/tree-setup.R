require(ape)
require(digest)  # md5 checksum
require(data.table)

#' Prepare a tree for clustering analysis
#'
#' @param phy:  ape::phylo object.  The tree should only contain "known" 
#'              sequences, i.e., none of the "New" sequences we want to predict.
#' @param seq.info:  data.table object.  Contains metadata for both known and 
#'                   new sequences, where the latter are identified by not 
#'                   appearing in the input tree `phy`.  Must contain sequence
#'                   labels under column label `Header`.
#' @param keep_root:  logical, override midpoint rooting.  Valid only if input
#'                    tree is already rooted.
#' @return The tree annotated with node information and seq.info
#' @export
import.tree <- function(phy, seq.info=data.table(), keep_root=FALSE) {
  # Midpoint root for consistency and resolve multichotomies
  if (is.rooted(phy)) {
    if (!keep_root) {
      warning("Re-rooting tree at midpoint. To retain original root, re-run",
              "with keep_root=TRUE.")
      phy <- phytools::midpoint.root(phy)  
    }
  } else {
    phy <- phytools::midpoint.root(phy)  
  }
  phy <- ape::multi2di(phy)
  cat(paste("\nRead in tree with", Ntip(phy), "tips\n"))

  # Check Sequence names inputs
  if (nrow(seq.info) == 0) {
      warning("No sequence meta-data included, creating default seq.info input", 
              "from tree tip.labels")
      seq.info <- data.table("Header"=phy$tip.label) 
  }
  
  # check sequence metadata (seq.info)
  var.names <- colnames(seq.info)
  if (!("Header" %in% var.names)) {
    stop("`var.name` must contain Header")
  } 
  else {
    if (length(unique(seq.info[, "Header"])) != length(seq.info[, "Header"])) {
      seq.info <- seq.info[!duplicated(Header), ]
      warning("Duplicate Headers have been arbitrarily removed")
    }
    if (!all(phy$tip.label %in% seq.info$Header)) {
      stop("At least one tip label in tree is not represented seq.info")
    }
  }
  
  # label new sequences
  seq.info$New <- !(seq.info$Header %in% phy$tip.label)
  cat(paste("\nIdentified", sum(seq.info$New), "new sequences\n"))
  phy$seq.info <- seq.info

  # parse bootstrap values, find descendants, calculate patristic distances
  phy$node.info <- annotate.nodes(phy)
  
  # get path lengths (from node to root)
  phy$path.info <- annotate.paths(phy)
  
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
  if (is.na(max.boot)) {
    max.boot <- 1
    if (sum(phy$node.label > 1, na.rm=T) / phy$Nnode > 0.5) {
      max.boot <- 100  
    }
  }
  
  # assign max value to nodes without bootstrap support
  phy$node.label[which(is.na(phy$node.label))] <- max.boot

  # number of nodes in a rooted binary tree is 2n-1
  # note tips are numbered first (1..n), followed by internal nodes
  node.info <- data.table(
    NodeID = 1:(2*Ntip(phy)-1),
    Bootstrap = c(rep(max.boot, Ntip(phy)), phy$node.label)
    )
  
  # normalize so we are always working on a scale of (0,1)
  node.info$Bootstrap <- node.info$Bootstrap / max.boot
  
  # Get descendant information for each node
  des <- phangorn::Descendants(phy, type="all")
  node.info[, "Descendants" := des]

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
      else { return (NA) }  # avoid NaN values
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
#' @return list of N matrices, where N is the number of nodes in the input tree.
#'         Each matrix has dimension 3 x P, where P is the path length from 
#'         the node to the root (inclusive).  The three rows are:
#'         - Node: integer index of node in path.
#'         - Boot: bootstrap support associated with the node
#'         - BranchLength: length of branch associated with the node
#'         For the root, no branch length is given (NA) and Boot is set to 1.
annotate.paths <- function(phy) {
  # Get node indices for all paths from terminal nodes to root
  paths <- ape::nodepath(phy)
  names(paths) <- sapply(paths, function(p) {
    p[length(p)]  # retrieve node index of tip
  })

  # Append path indices from internal nodes to root (not including root itself)
  i <- 1
  while (length(paths) < nrow(phy$node.info) + 1) {
    new.paths <- sapply(paths[i:length(paths)], function(x) {
      x[-length(x)]
    })
    new.paths <- new.paths[!duplicated(new.paths)]
    names(new.paths) <- sapply(new.paths, function(p) {
      p[length(p)]  # label with last node index
    })
    # FIXME: is this really the most efficient way?
    i <- length(paths) + 1
    paths <- c(paths, new.paths[which(!(names(new.paths) %in% names(paths)))])
  }
  
  # drop paths of no nodes and reorder
  path.len <- sapply(paths, length) 
  paths <- paths[which(path.len > 0)]
  paths <- paths[order(as.numeric(names(paths)))]
  
  # retrieve branch lengths for all nodes in paths
  lens <- lapply(paths, function(p) phy$edge.length[match(p, phy$edge[,2])])

  # Obtain bootstrap supports for all nodes in paths
  boots <- lapply(paths, function(p) phy$node.info$Bootstrap[p])
  
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


#' Place new sequences on the tree.
#' This function makes calls to pplacer and guppy binaries (Matsen, 2010) 
#' provided with this package.  pplacer ensures that previously defined clusters 
#' remain the same. New tips are assigned a set of possible placements in the 
#' tree.
#' 
#' @param phy:  ape::phylo object.  Must be annotated with seq.info, node.info,
#'              and path.info by calling import.tree().
#' @param seqs.full:  ape::DNAbin or AAbin object.  This sequence alignment 
#'                    must contain all known sequences in the tree `phy`, as 
#'                    well as all new sequences we want to graft to the tree
#'                    as potential cluster growth.
#' @param log.file:  character.  A path to the logfile from a tree reconstruction 
#'                   run.  This file can be produced by IQTREE, FastTree or RAxML.
#' @return  ape::phylo object with growth.info field
#' @export
extend.tree <- function (phy, seqs.full, log.file=NA, locus = "LOCUS") {
  # Extend with growth_info
  if (is.na(log.file)) {
    warning("Ignoring growth information, path to logfile and full sequence ", 
            "alignment required.")
    phy$growth.info <- data.table(
      "Header"=character(0), "NeighbourDes"=numeric(0), "Bootstrap"=numeric(0),
      "TermDistance"=numeric(0), "PendantDistance"=numeric(0), "Terminal"=logical(0)
    )
    return(phy)
  } 
  
  if(!all(names(full.align) %in% seq.info$Header)){
    stop("Headers in full sequence alignment do not correspond to seq.info data")
  }
  
  # call pplacer to graft new sequences to the tree
  stats.json <- translate.log(log.file)
  refpkg <- taxit.create(phy, full.align, stats.json)
  ptrees  <- run.pplacer_guppy(refpkg)
  
  # process trees with new tips
  growth.info <- lapply(ptrees, function(x) annotate.growth(phy, x))
  growth.info <- do.call(rbind, growth.info)
  
  # Collapse neighbour node descendant tips to their MRCA
  collapsed.neighbours <- sapply(
    growth.info[!(Terminal), NeighbourDes], function(des) {
      ape::getMRCA(phy, des)
    })
  
  temp <- growth.info[, NeighbourDes]
  temp[which(!growth.info$Terminal)] <- collapsed.neighbours
  
  suppressWarnings(growth.info[, "NeighbourNode" := unlist(temp)])
  growth.info[, NeighbourDes := NULL]
  
  if (!all(growth.info$Header %in% phy$seq.info$Header[phy$seq.info$New])) {
    # FIXME: sometimes tip labels are merged in pplacer output
    warning("Not all newly added sequences are noted in the seq.info of the tree")
  }
  phy$growth.info <- growth.info  # attach to tree
  
  return(phy)
}


#' Add the growth information onto the known tree.
#'
#' @param phy:  ape::phylo object.  Return value of extend.tree().
#' @param ptree:  ape::phylo object.  A tree generated by pplacer with a 
#'                newly added tip placed by maximum likelihood.
#' @return data.table. The growth.info for a given tree. This includes the 
#'         branch length to the newly added terminal node, whether or not the 
#'         neighbour to the newly added terminal node is also terminal, and the 
#'         branch length to that neighbour, from the newly placed internal node 
#'         (the "PendantLength")
annotate.growth <- function(phy, ptree) {
  # extract set of possible placements of new tip
  new.ids <- setdiff(ptree$tip.label, phy$tip.label)
  new.tips <- which(ptree$tip.label %in% new.ids)  # node indices
  
  # extract tip label prefix (without M value from pplacer)
  header <- gsub("_#.*", "", new.ids[1])
  
  # retrieve indices of in-edges to new tips
  new.edges <- which(ptree$edge[, 2] %in% new.tips)
  term.dists <- ptree$edge.length[new.edges]
  
  # get index of parent node created by grafting new tip to a branch
  new.nodes <- ptree$edge[new.edges, 1]
  new.nodes.des <- which(ptree$edge[, 1] %in% new.nodes)  # all descendants
  node.boots <- sapply(new.ids, function(x) {
    l <- strsplit(x, "_M=")[[1]]
    as.numeric(l[length(l)])
  })

  # a pendant edge is a branch from the parent node to a descendant 
  # in the original tree
  pen.edges <- setdiff(new.nodes.des, new.edges)
  pen.dists <- ptree$edge.length[pen.edges]
  
  neighbour.des <- lapply(new.nodes, function(n) {
    phangorn::Descendants(ptree, n, "tips")[[1]]
  })
  old.tips <- lapply(neighbour.des, function(des) {
    which(phy$tip.label %in% ptree$tip.label[des])
  })
  terminal <- sapply(old.tips, function(x) {
    length(x) == 1
  })
  
  growth.info <- data.table::data.table(
    "Header" = header, "NeighbourDes" = old.tips, 
    "Bootstrap" = node.boots, "TermDistance" = term.dists, 
    "PendantDistance" = pen.dists, "Terminal" = terminal
  )
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
