library(igraph)
library(tidyverse)
library(data.table)

# who has the most neighbors in different communities?
cl <- cluster_louvain(g)
V(g)$cl <- cl$membership
V(g)$num_cl <- 0  # number of connected communities

calc_numcl <- function(x){return(length(table(V(g)$cl[neighbors(g,x)])))}
V(g)$num_cl <- unlist(lapply(V(g), calc_numcl))

cpset <- read_csv('CPSet.csv', col_names = FALSE)
cpset <- cpset$X2
V(g)$cp <- cpset

#calc_numcp
# write this function




# add assignment tables -------------------------------------------------------------
library(knitr)
#gr1 <- read_graph('polblogs/polblogs.gml', format='gml')
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

df <- NULL
df <- rbind(c('Facebook',graph_char(g)))
colnames(df) <- c('Network','Type','n','m','c (weak, strong)','d','l','L','ccl','ccg')
rownames(df) <- c()

df_networkchar <- kable(df, 
caption = 'Table 1. FB network and properties. Type indicates directed/undirected, 
n = number of nodes, m = number of links, c = number of connected components 
(weak/strong separate in the case of a directed network), d = maximum degree, 
l = average path length, L = diameter, ccl = average local clustering coefficient, 
and ccg = global clustering coefficient (3 times number of triangles/ number of 
connected triplets). (See appendix for description/rationale behind procedure.)')


#  (i) Degree, (ii) Eccentricity, (iii) Closeness,
# (iv) Betweenness, (v) PageRank, (vi) Kleinberg’s Authority score, and (vii) Kleinberg’s Hub score
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
  close.cent <- which(closeness(gr) == min(closeness(gr)))
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
  closecent <- 1 / (min(closeness(gr)))
  btwncent <- max(betweenness(gr))
  pagerank <- signif(max(page_rank(gr)$vector),4)
  k.auth <- max(authority_score(gr)$vector)
  k.hub <- max(hub_score(gr)$vector)
  attrlist <- list('degree'=degcent, 'eccentricity'=eccent,
                   'closeness' = closecent, 'between'=btwncent,
                   'pagerank'=pagerank, 'authority'=k.auth, 'hub'=k.hub)
  return(attrlist)
}

nodecentr <- node.centrality(g)
centr_measure <- centrality.measure(g)
labelrows <- c('degree','eccentricity','closeness','between','pagerank',
               'authority','hub')
df <- data.table(labelrows,nodecentr, centr_measure)
colnames(df) <- c('centrality measure','FB (max. node)', 'FB (values)')
df_centrality <- kable(df, 
                       caption = 'Table 1. Node(s) in each network with the 
                       highest scores in terms of centrality measures (listed 
                       in the leftmost column). Column labels include network 
                       labels, (max. node) indicates the node with the highest 
                       centrality measure in the network, and (values) indicate 
                       the maximized centrality measure itself for that node. For 
                       entries ending in "...", see the appendix for a complete 
                       list of values.')



  
grm <- laplacian_matrix(gr)
eig_fb <- eigen(grm)
#eig <- eigen(grm)$values
lambda1 <- max(eig)
lambda2 <- eig[which.min(eig) - 1]

ev1 <- eig_fb$vector[,1]
ev2 <- eig_fb$vector[,(vcount(gr)-1)]
df <- data.frame(ev1,ev2); df$x <- 1:vcount(gr)

ggplot(df, aes(x=x,y=ev1)) +
  geom_line(colour='red') + #geom_point(colour='red') +
  geom_line(aes(y=ev2), colour='blue') + #geom_point(aes(y=ev2), colour='blue') + 
  labs(x='vertex id number', y='eigenvector value')

ggplot(data.frame(x=1:7, y=c(2456,564,425,221,210,111,52)),
       aes(x=x, y=y)) + geom_point() + scale_x_log10() +
  scale_y_log10()

mean(degree(g))
mean(degree(g, v=V(g)[V(g) %in% coreset$X1]))
mean(degree(g, v=V(g)[!(V(g) %in% coreset$X1)]))

median(degree(g))
median(degree(g, v=V(g)[V(g) %in% coreset$X1]))
median(degree(g, v=V(g)[!(V(g) %in% coreset$X1)]))

sd(degree(g))
sd(degree(g, v=V(g)[V(g) %in% coreset$X1]))
sd(degree(g, v=V(g)[!(V(g) %in% coreset$X1)]))

IQR(degree(g))
IQR(degree(g, v=V(g)[V(g) %in% coreset$X1]))
IQR(degree(g, v=V(g)[!(V(g) %in% coreset$X1)]))



