library(igraph)
library(tidyverse)
library(data.table)

calculate_rd <- function(g, alpha, beta){
  g <- as.undirected(g)
  u <- rep(NA, vcount(g))
  V(g)$tag <- 0
  u[1] <- which(closeness(g) == max(closeness(g)))
  V(g)$tag[u[1]] <- 1
  tmp <- NULL
  for(i in na.omit(u)){
    tmp <- c(tmp, neighbors(g, V(g)[i]))
  }
  tmp <- subset(tmp, tmp %in% V(g)[V(g)$tag==0])
  for(i in 2:vcount(g)){
    xc <- as.data.frame(table(tmp))
    cand <- as.numeric(levels(droplevels(xc$tmp[xc$Freq==max(xc$Freq)])))
    if(length(cand) > 1){
      degsmp <- degree(g)[cand]
      cand <- cand[which.max(degsmp)]
      if(length(cand) > 1){
        cand <- sample(cand, 1)
      }
    }
    u[i] <- cand
    V(g)$tag[cand] <- 1
    tmp <- c(tmp, neighbors(g, V(g)[u[i]]))
    tmp <- subset(tmp, tmp %in% V(g)[V(g)$tag==0])
  }
  rd <- 0
  for(i in 2:length(u)){
    lower <- i - (alpha-1)
    if(lower < 1){
      lower <- 1
    }
    s <- induced_subgraph(g, u[lower:i])
    rd[i] <- (2*ecount(s))/(vcount(s)*(vcount(s)-1))
  }
  return(rd)
}
cascade_capacity <- function(g, adopt){
  g <- as.undirected(g)
  qlist <- seq(0,0.5,0.01)
  prlist <- function(x1){
    rnp <- V(g)$info[unlist(x1)]
    pt <- sum(rnp)/length(rnp)
  }
  len <- length(V(g))
  for(q in qlist){
    V(g)$info <- 0
    V(g)$info[adopt] <- 1
    for(i in 1:(vcount(g)-1)){
      fr <- sum(V(g)$info)/len
      uninit <- V(g)[V(g)$info==0]  # uninitiated nodes
      tmp <- adjacent_vertices(g, uninit)
      p <- unlist(lapply(tmp, prlist))
      V(g)$info[uninit[(is.finite(p)) & (p >= q)]] <- 1
      dfr <- (sum(V(g)$info)/len) - fr
      if(dfr==0) break
    }
    if(fr < 1) break
  }
  if(q > 0){q <- (q - 0.01)}
  return(q)
}
rumor_transmit <- function(g, infect, tsteps=100, adopt){
  g <- as.undirected(g)
  V(g)$info <- 0
  V(g)$info[adopt] <- 1
  for(t in 1:tsteps){
    tmp <- unlist(adjacent_vertices(g, V(g)[V(g)$info >= 1]))
    prst <- sample(c(0,1), length(tmp), 
                   replace=TRUE, prob =c((1-infect), infect))
    V(g)$info[tmp[prst==1]] <- V(g)$info[tmp[prst==1]] + 1
  }
  return(g)
}
graph_char <- function(gr){
  directed <- is_directed(gr)
  if(directed==TRUE){
    directed <- 'Directed'
    links <- sum(degree(gr, mode='out'))
  }else{
    directed <- 'Undirected'
    links <- (sum(degree(gr)))/2
  }
  numNodes <- length(V(gr))
  
  
  #links <- sum(degree(gr, mode='out'))
  if(is_directed(gr)==TRUE){
    comp <- paste(count_components(gr),
                  count_components(gr,mode='strong'),sep=', ')
  }else{
    comp <- count_components(gr)
  }
  
  dmax <- max(degree(gr))
  avgpath <- mean_distance(gr, directed=is_directed(gr))#,
  #unconnected= !(is_connected(gr)))
  diam <- diameter(gr)
  ccl <- transitivity(gr, type='localaverage')
  ccg <- transitivity(gr, type='global')
  
  grlist <- c(directed, numNodes, links, comp,
              dmax, signif(avgpath,4), diam, 
              signif(ccl,4), signif(ccg,4))
  return(grlist)
}
node.centrality <- function(gr){
  if(is_directed(gr)==TRUE){
    degcent <- which.max(degree(gr, mode='in'))
    degcent[2] <- which.max(degree(gr, mode='out'))
  }else{
    degcent <- which.max(degree(gr))
  }
  tmp <- eccentricity(gr)
  tmp[tmp==0] <- NA
  ec.cent <- which(tmp == min(tmp, na.rm=TRUE))
  close.cent <- which(closeness(gr) == max(closeness(gr)))
  btwn.cent <- which(betweenness(gr) == max(betweenness(gr)))
  pr.cent <- which(page_rank(gr)$vector == max(page_rank(gr)$vector))
  k.auth <- which(authority_score(gr)$vector == max(authority_score(gr)$vector))
  k.hub <- which(hub_score(gr)$vector == max(hub_score(gr)$vector))
  attrlist <- list('degree'=degcent, 'eccentricity'=ec.cent,
                   'closeness' = close.cent, 'between'=btwn.cent,
                   'pagerank'=pr.cent, 'authority'=k.auth, 'hub'=k.hub)
  return(attrlist)
}
centrality.measure <- function(gr){
  if(is_directed(gr)==TRUE){
    degcent <- max(degree(gr, mode='in'))
    degcent[2] <- max(degree(gr, mode='out'))
  }else{
    degcent <- max(degree(gr))
  }
  eccent <- 1 / (min(eccentricity(gr)[eccentricity(gr) != 0]))
  closecent <- 1/max(closeness(gr))
  btwncent <- max(betweenness(gr))
  pagerank <- signif(max(page_rank(gr)$vector),4)
  k.auth <- max(authority_score(gr)$vector)
  k.hub <- max(hub_score(gr)$vector)
  attrlist <- list('degree'=degcent, 'eccentricity'=eccent,
                   'closeness' = closecent, 'between'=btwncent,
                   'pagerank'=pagerank, 'authority'=k.auth, 'hub'=k.hub)
  return(attrlist)
}


