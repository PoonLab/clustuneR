# A library of functions 
library(dplyr, quietly=TRUE, verbose = FALSE)

#Creates a set of data-frames representing a graph of sequences, with the edges between those sequences representing the TN93 Distance.
#Sequences must be dated with the date separated from the id by '_'. 
impTN93 <- function(iFile, minNS=63){
  #@param iFile: The name/path of the input file (expecting tn93 output csv)
  #@param minNS: The minimum number of acceptible new Sequences. By default we keep this high.
  #@return: A list of 3 Data frames. An edge list (weighted by TN93 genetic distance), a vertex list, 
  #         and a list of minimum edges, for the future establishment of a timepoint-based model
  
  #Reading input file as a data frame. This will essentially act as an edgelist
  idf <- read.csv(iFile, stringsAsFactors = F)
  
  # accession
  temp1 <- sapply(idf$ID1, function(x) (strsplit(x,'_')[[1]])[[1]])
  # collection date
  temp2 <- sapply(idf$ID1, function(x) (strsplit(x,'_')[[1]])[[2]])
  
  temp3 <- sapply(idf$ID2, function(x) (strsplit(x,'_')[[1]])[[1]])
  temp4 <- sapply(idf$ID2, function(x) (strsplit(x,'_')[[1]])[[2]])
  
  #Create a data frame from the imported edge list. 
  el <- data.frame(ID1=as.character(temp1), t1=year(as.Date(temp2)), 
                   ID2=as.character(temp3), t2=year(as.Date(temp4)), 
                   Distance = as.numeric(idf$Distance), stringsAsFactors= F)
  
  #Obtain the maximum time and time difference between the head and tail of each edge
  el$tMax <- pmax(el$t1, el$t2)
  el$tDiff <- abs(el$t1-el$t2)
  
  #Create a list of vertices based on the edge list. Original sequences with no edges will not be considered.
  vl <- unique(data.frame(ID = c(el$ID1, el$ID2), Time = c(el$t1, el$t2), stringsAsFactors=F))
  
  #Order edges and vertices by time point and sort them into a single, larger list
  #If the newest timepoint contains a small number of sequences, we remove the newest year from consideration
  g <- list(v=vl[order(vl$Time),], e=el[order(el$tMax),], f=el[order(el$tMax),])
  while(nrow(subset(g$v,Time==max(Time)))<=minNS) {
    keepT <- head(as.numeric(names(table(g$v$Time))),-1)
    g <- tFilt(g, keepT)
  }
  
  #Permanently remove edges from the new year such that only the closest edge between new vertices and old vertices remains
  #Obtain new vertices and remove any internal edges within the new vertices 
  #We are not interested in completely new clustering
  nV <- subset(g$v, Time==max(Time))
  
  #The subset of vertices excluding those at the oldest year
  subV <- subset(g$v, Time>min(Time))
  
  #Obtain the closest retrospective edge of every vertex beyond the oldest year
  minE <- bind_rows(lapply(1:nrow(subV), function(i){
    v <- subV[i,]  # for each node in subV
    incE <- subset(g$e, (ID1%in%v$ID)|(ID2%in%v$ID))
    retE <- subset(incE, (tMax==v$Time)&(tDiff>0))  # retrospective, bipartite edges
    retE[which(retE$Distance==min(retE$Distance))[[1]],]  # reduce to shortest edge
  }))
  
  #Only closest retrospective edges are kept for edges from new cases. 
  g$e <- subset(g$e, tMax!=max(tMax))
  g$e <- rbind(g$e, subset(minE, tMax==max(tMax)))
  g$f <- subset(minE, tMax < max(tMax))

  return(g)
}

#Create clusters based on component clustering by some measure of genetic distance
compClu <- function(iG) {
  #@param iG: The inputted graph. Expecting all vertices, but some edges filtered by distance.
  #@return: The inputted graph, annotated with a cluster size summary and case membership in the vertices section
  
  #Simplify the list of unsorted vertices (just id's) and edges (just head and tail id's)
  vid <- iG$v[,"ID"]
  adj <- iG$e[,c("ID1","ID2")]
  
  #Initialize the first cluster name and a column for cluster membership. c0 will be reserved for all singletons if sing=T
  iG$v$Cluster <- vector(mode="numeric", length=nrow(iG$v))
  
  #The search vertex becomes the first member of the first cluster and is removed from the searchable set of cluster names
  i <- 1
  srchV <- vid[1]
  memV <- srchV
  vid <- setdiff(vid, memV)
  
  #Assigning Cluster Membership
  repeat {
    
    #Remove edges internal to search query and list outgoing edges
    adj <- subset(adj, !(ID1%in%srchV & ID2%in%srchV))
    exE <- subset(adj, ID1%in%srchV | ID2%in%srchV)

    #Find all neighbouring vertices to the search vertex (or search vertices) through external edges
    #These are then added to the list of member vertices and removed from the list of searchable vertices
    nbV <- setdiff(c(exE$ID1,exE$ID2), srchV)
    memV <- c(memV, nbV) 
    vid <- setdiff(vid, nbV)

    #If there are no more neigbours to the search vertices, the cluster is completed and we reset the search parameters
    if (length(nbV)==0) {
      
      iG$v$Cluster[iG$v$ID%in%memV] <- i
      
      #The end condition, catching the event that there are no vertices to assign to clusters
      if (length(vid)==0) {break}
      
      #Reset search parameters
      i <- i+1
      srchV <- vid[1]
      memV <- srchV
      vid <- setdiff(vid, memV)

      next
    }
    
    #Remove all edges within the current cluster from the adjacency list
    adj <- subset(adj, !(ID1%in%srchV | ID2%in%srchV))
    srchV <- nbV
  }

  #Add some summary information regarding clusters
  iG$c <- table(iG$v$Cluster)
  
  return(iG)
}

#A simple function, removing edges that sit above a maximum reporting distance (@param:maxD).
dFilt <- function(iG, maxD) {
  iG$e <- subset(iG$e, Distance<=maxD)
  return(iG)
}

#A simple function, removing vertices that sit above a maximum time point (@param: maxT)
tFilt <- function(iG, keepT) {
  iG$v <- subset(iG$v, Time%in%keepT)
  iG$e <- subset(iG$e, tMax%in%keepT)
  return(iG)
}


#' Simulate the growth of clusters, showing the difference in cluster size between 
#' the newest and the penultimate time point
#' The frame of reference for clusters is the penultimate year, simulating one 
#' making forcasting decisions based on one time point and validating them with 
#' the next
#' @param: The inputted graph. Expecting all vertices, but some edges filtered 
#' by distance.
#' @return: The same cluster annotated with the actual growth and cluster 
#' information
simGrow <- function(iG) {
  #Obtain clusters at the new time point, after removing singletons
  nG <- iG
  
  # subset of new nodes that have no edges to known cases
  nSing <- subset(nG$v, (!ID%in%c(nG$e$ID1, nG$e$ID2) & Time==max(Time)))
  
  # subset of new nodes that DO have edges to known cases, AND all known cases
  nG$v <- subset(nG$v, !(!ID%in%c(nG$e$ID1,nG$e$ID2) & Time==max(Time)))
  nG <- compClu(nG)  # assign nodes to connected components
  
  # obtain clusters at an old time point
  keepT <- head(as.numeric(names(table(iG$v$Time))),-1)
  oG <- compClu(tFilt(iG, keepT))
  
  # Define growth as the difference in cluster size between new and old graphs 
  # After clsFilter(), nG will have the same number of clusters as oG, and 
  # similar membership
  iG$g <- nG$c-oG$c
  iG$c <- oG$c
  
  #Re-Add the singletons, citing new singletons as members of the cluster 0
  if (nrow(nSing)>0){
    nSing$Cluster <- 0
  }
  iG$v <- rbind(nG$v, nSing)
  
  return(iG)
}

#Obtains some likelihood data in order to weight cases based on their recency  
likData <- function(iG) {
  #@param iG: The inputted graph. Expecting the entire Graph with new year included.
  #@return: A data frame of "Positives" (related cases) between one time point and another, annotated with the number of possible total positives.
  #         Each positive also carries it's Genetic Distance measurement and time point difference (between time points)
  
  #Take in total graph without the newest time point (otherwise, we include the validation set in or data set)
  keepT <- tail(head(as.numeric(names(table(iG$v$Time))),-1),-1)
  subG <- tFilt(iG, keepT)
  
  #Obtain the closest retrospective edge of every vertex
  f <- bind_rows(lapply(2:nrow(subG$v), function(i){
    
    v <- subG$v[i,]
    
    incE <- subset(iG$e, (ID1%in%v$ID)|(ID2%in%v$ID))
    retE <- subset(incE, (tMax==v$Time)&(tDiff>0))
    minE <- retE[which(retE$Distance==min(retE$Distance))[[1]],]
    df <- data.frame(Distance=minE$Distance, tMax=minE$tMax, tDiff=minE$tDiff, vTotal=nrow(subset(iG$v, Time==v$Time)))
    
    return(df)
    
  }))
  
  return(f)
}

#Analyze a given Clustered Graph to establish the difference between the performance of two different models
#Performance is defined as the ability for cluster growth to fit a predictive model.
compAnalyze <- function(subG) {
  #@param subG: A subGraph cut based on a threshold distance, expecting a member of the multiGraph set
  #@return: A graph annotated with growth, cluster info and level of predictive performance (measured through GAIC)
  
  #Obtain some summarized information from the sub-Graph
  dMax <- max(subG$e$Distance)  # e is edge list filtered by threshold

  # f is unfiltered edge list, excluding edges to new cases and between cases from time
  tTab <- table(subG$f$tMax)  # table of sampling times of most recent node for every edge in $f
  
  # v is node list derived from unfiltered edge list (excludes singletons if input is not complete graph)
  vTab <- table(subG$v$Time)  # table of sampling times, including most recent nodes
  eTab <- sapply(as.numeric(names(vTab)), function(t){
    nrow(subset(subG$e, (t1==t|t2==t)))  # how many filtered edges associated with each sampling time?
  })
  names(eTab) <- names(vTab)
  
  #Take the total edge frequency data from the graph and format this information into successes and attempts
  #An edge to the newest year falling below the max distance is considered a success
  ageD <- bind_rows(lapply(as.numeric(names(tTab)), function(t) {
    # for every sampling time (not including earliest)
    temp <- subset(subG$f, tMax==t)  # get associated bipartite edges
    dfs <- split(temp, temp$tDiff)  # split by start times
    
    # how many edges in bipartite graph meet distance threshold?
    Positive <- sapply(dfs, function(df){length(which(df$Distance<=dMax))})
    vTotal <- rep((vTab[[as.character(t)]]),length(dfs))  # FIXME: unused?
    tDiff <- as.numeric(names(Positive))
    oeDens <- sapply(tDiff, function(tD){
      # to calculate the average degree (out-edge density) of known cases from the same time point
      oTime <- t-tD  # start time for this bipartite graph
      
      # return the number of filtered edges with a node from time oTime
      # divide by total number of nodes from oTime
      return(eTab[as.character(oTime)]/vTab[as.character(oTime)])
    })
    
    res <- data.frame(Positive=as.numeric(Positive), vTotal=vTab[[as.character(t)]], oeDens=oeDens, tDiff)
    return(res)
  }))
  
  #Obtain a model of case connection frequency to new cases as predicted by individual case age
  #Use this to weight cases by age
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
  # subG$v$Weight <- sapply(subG$v$ID, function(id) {as.numeric(substr(id,8,8))})
  
  #Create clusters for this subgraph and measure growth
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
  
  #Save, gaic, model and age data as part of the output
  subG$gaic <- fit1$aic-fit2$aic
  subG$ageMod <- mod
  subG$ageFit <- fit1
  subG$nullFit <- fit2
  subG$f <- ageD

  return(subG)
}

#Run across a set of several subGraphs created at various filters, analyzing GAIC at each with clusterAnalyze
gaicRun <- function(iG, cutoffs=NA) {
  #@param iG: Expecting the entire Graph, but in some cases may take a subset  
  #@return: A data frame of each runs cluster information (clusterAnalyze output)
  
  #Initialize a set of cutoffs to observe (based on the genetic distance distribution)
  if (all(is.na(cutoffs))) {
    steps <- head(hist(subset(iG$e, Distance<0.05)$Distance, plot=FALSE)$breaks, -5)
    cutoffs <- seq(0 , max(steps), max(steps)/50) 
  }

  #A set of several graphs created at different cutoffs
  gs <- lapply(cutoffs, function(d) {dFilt(iG, d)})
  
  #Generate cluster data for each subGraph in gs
  res <- lapply(gs, function(subG) {compAnalyze(subG)})
  names(res) <- cutoffs
  
  return(res)
}
