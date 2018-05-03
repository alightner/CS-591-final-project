library(tidyverse)
library(igraph)
library(data.table)

g <- read_graph('facebook_combined.txt', format='edgelist')
g <- as.undirected(g)

df <- data.frame(matrix(nrow=2,ncol=5))
len <- length(V(g))

page <- page_rank(g)
adopt <- c(which.max(page$vector), which.max(degree(g)))
#adopt <- which(page$vector >= 0.004)
#adopt <- which.max(closeness(g))
#adopt <- which.min(closeness(g))
#adopt <- which.max(degree(g)) == max close
#adopt <- which.max(betweenness(g)) == max between

# visualization
E(g)$width <- 0.05
V(g)$frame.color <- NA
V(g)$color <- 'black'
V(g)$color[which.max(page_rank(g)$vector)] <- 'yellow'
V(g)$color[which.max(degree(g))] <- 'blue'
plot(g, vertex.size=1, vertex.label=NA,
     layout=layout.lgl(g))

plot(g, vertex.size=1, vertex.label=NA,
     layout=layout.fruchterman.reingold(g))

#plot(g, vertex.size=1, vertex.label=NA, 
 #    edge.color='blue')

#plot(g, vertex.size=1, vertex.label=NA,
 #    layout=layout.lgl(g),vertex.color=graphCol)

#plot(g, vertex.size=1, vertex.label=NA,
 #    layout=layout.fruchterman.reingold(g),vertex.color=graphCol)

# start simulation
finfr <- 1
q <- 0

starttime <- Sys.time()

while(finfr==1){
  
  V(g)$info <- 0
  V(g)$info[adopt] <- 1
  
  t <- 0
  dfr <- 1
  infr <- sum(V(g)$info)/len
  
  while(dfr > 0){
    fr <- sum(V(g)$info)/len
    #uninit <- V(g)[V(g)$info==0]
    tmp <- adjacent_vertices(g, V(g)[V(g)$info==0])
    p <- rep(NA, length(tmp))
    for(i in 1:length(p)){
      rnp <- V(g)$info[unlist(tmp[i])]
      p[i] <- sum(rnp)/length(rnp)
    }
    V(g)$info[uninit[(is.finite(p)) & (p >= q)]] <- 1
    t <- t+1
    dfr <- (sum(V(g)$info)/len) - fr
  }
  
  finfr <- fr
  
  df[1,] <- df[2,]
  df[2,] <- c(q,t,infr,fr,'max_PR')
  #df <- rbind(df, c(q,t,infr,fr,'max_PR'))
  
  q <- q+0.01
  
}

endtime <- Sys.time()

colnames(df) <- c('q','time','init_fr', 'end_fr','pr')

endtime - starttime

#df_maxPR <- df
#df_leq002 <- df
#df_leq004 <- df
#df_maxclose <- df
#df_minclose <- df
#df_both <- df



# ----- Epidemiological -----

g <- read_graph('facebook_combined.txt', format='edgelist')
g <- as.undirected(g)

E(g)$width <- 0.1
V(g)$frame.color <- NA
#plot(g, vertex.size=1, vertex.label=NA, 
 #    vertex.color='black',edge.color='blue', layout=layout.lgl(g))
V(g)$info <- 0

# set parameters
infect <- 0.1
tsteps <- 100

V(g)$info <- sample(c(0,1), size=length(V(g)), 
                    replace=TRUE, prob=c((1-pr), pr))
V(g)[info==0]$color <- 'black'
V(g)[info==1]$color <- 'yellow'

plot(g, vertex.size=1, vertex.label=NA, 
     edge.color='blue', layout=layout.lgl(g))

s <- Sys.time()
for(t in 1:tsteps){
  tmp <- unlist(adjacent_vertices(g, V(g)[V(g)$info >= 1]))
  prst <- sample(c(0,1), length(tmp), 
                 replace=TRUE, prob =c((1-infect), infect))
  V(g)$info[tmp[prst==1]] <- V(g)$info[tmp[prst==1]] + 1
}
e <- Sys.time()
e-s
# coloring nodes by salience
E(g)$width <- 0.05
heard <- V(g)$info
fine <- 500 # this will adjust the resolving power
pal <- colorRampPalette(c('darkgreen','lightgreen'))
graphCol = pal(fine)[as.numeric(cut(btw,breaks = fine))]
plot(g, vertex.size=1, vertex.label=NA,
     layout=layout.lgl(g),vertex.color=graphCol)

plot(g, vertex.size=1, vertex.label=NA,
     layout=layout.fruchterman.reingold(g),vertex.color=graphCol)

plot(g, vertex.size=1, vertex.label=NA,
     layout=layout.kamada.kawai(g),vertex.color=graphCol)

pagerank_distr <- ggplot(data.frame(x=page_rank(g)$vector), aes(x=x)) +
  geom_density() + 
  geom_vline(xintercept = median(page_rank(g)$vector),linetype=3) +
  geom_vline(xintercept = (median(page_rank(g)$vector) + 
                             IQR(page_rank(g)$vector)),linetype=3) +
  geom_vline(xintercept = 0.002,linetype=3) + 
  geom_vline(xintercept = 0.004,linetype=3) +
  geom_vline(xintercept = max(page_rank(g)$vector),linetype=3) +
  theme_bw() +
  labs(x='PageRank')
ggsave('pagerank_distr.jpg', pagerank_distr)
  






































# ----- Cascade capacity -----

# set approach (test)

g <- read_graph('karate/karate.gml', format='gml')
df <- data.frame(matrix(nrow=2,ncol=5))
len <- length(V(g))

page <- page_rank(g)
adopt <- which.max(page$vector)

finfr <- 1
q <- 0

while(finfr==1){
  
  V(g)$info <- 0
  V(g)$info[adopt] <- 1
  
  t <- 0
  dfr <- 1
  infr <- sum(V(g)$info)/len
  
  while(dfr > 0){
    fr <- sum(V(g)$info)/len
    uninit <- V(g)[V(g)$info==0]
    tmp <- adjacent_vertices(g, uninit)
    p <- rep(NA, length(tmp))
    for(i in 1:length(tmp)){
      rnp <- V(g)$info[unlist(tmp[i])]
      p[i] <- sum(rnp)/length(rnp)
    }
    V(g)$info[uninit[(is.finite(p)) & (p >= q)]] <- 1
    t <- t+1
    dfr <- (sum(V(g)$info)/len) - fr
  }
  
  finfr <- fr
  
  df[1,] <- df[2,]
  df[2,] <- c(q,t,infr,fr,'max_PR')
  #df <- rbind(df, c(q,t,infr,fr,'max_PR'))
  
  q <- q+0.01
  
}

colnames(df) <- c('q','time','init_fr', 'end_fr','pr')



# probabilistic approach (Facebook data)

g <- read_graph('facebook_combined.txt', format='edgelist')
g <- as.undirected(g)

#g <- read_graph('karate/karate.gml', format='gml')
# set parameters
E(g)$width <- 0.2
V(g)$frame.color <- NA
plot(g, vertex.size=5, vertex.label=NA, 
     vertex.color='black',edge.color='blue')

pr <- 0.1  # probability of initiating rumor, story, etc.
qlist <- seq(0,0.5,by=0.01)
#df <- data.table()
#df <- data.frame()
df <- data.frame(matrix(NA, nrow = length(qlist), ncol = 5))
colnames(df) <- c('q','time','init_fr', 'end_fr','pr')
#q,t,infr,fr,pr
len <- length(V(g))

starttime <- Sys.time()

trial <- round(runif(1,min=1,max=100000))

# change qlist to while < delta qlist? or something like that?

for(q in qlist){
  set.seed(trial)
  V(g)$info <- sample(c(0,1), size=length(V(g)), 
                      replace=TRUE, prob=c((1-pr), pr))
  t <- 0
  dfr <- 1
  infr <- sum(V(g)$info)/len
  while(dfr > 0){
    fr <- sum(V(g)$info)/len
    uninit <- V(g)[V(g)$info==0]
    tmp <- adjacent_vertices(g, uninit)
    p <- rep(NA, length(tmp))
    for(i in 1:length(tmp)){
      rnp <- V(g)$info[unlist(tmp[i])]
      p[i] <- sum(rnp)/length(rnp)
    }
    V(g)$info[uninit[(is.finite(p)) & (p >= q)]] <- 1
    t <- t+1
    dfr <- (sum(V(g)$info)/len) - fr
  }
  
  df[(q*100),] <- c(q,t,infr,fr,pr)
  #df <- rbind(df, c(q,t,infr,fr,pr))
}

endtime <- Sys.time()

#colnames(df) <- c('q','time','init_fr', 'end_fr','pr')

endtime - starttime

#Time difference of 1.627798 secs! add multiple seed trials?

ggplot(df, aes(x=q, y=end_fr)) + geom_point() + geom_smooth(method='lm') +
  geom_jitter()