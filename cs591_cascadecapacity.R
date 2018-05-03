library(tidyverse)
library(igraph)

# ---- rumor or trend simulation / cascade (probabilistic approach) --------------------
#g <- read_graph('facebook_combined.txt', format='edgelist')
#g <- as.undirected(g)
#plot(g, vertex.size=1, vertex.label=NA, 
 #    layout=layout.fruchterman.reingold(g))
#plot(g, vertex.size=1, vertex.label=NA, 
 #    layout=layout.kamada.kawai(g))


set.seed(192837)
g <- erdos.renyi.game(100,0.1)
E(g)$width <- 0.5
V(g)$frame.color <- NA
plot(g, vertex.size=5, vertex.label=NA, 
     vertex.color='black',edge.color='blue')
V(g)$info <- 0


# set parameters
dtab <- data.frame()
set.seed(34564)
pr <- 0.01  # probability of initiating rumor, story, etc.

starttime <- Sys.time()

qlist <- seq(0,0.5,by=0.01)
for(q in qlist){
  V(g)$info <- sample(c(0,1), size=length(V(g)), 
                      replace=TRUE, prob=c((1-pr), pr))
  V(g)[info==0]$color <- 'blue'
  V(g)[info==1]$color <- 'red'
  #plot(g, vertex.size=3, vertex.label=NA)
  #should this be more systematic, rather than random (e.g., iteratively start from each node once)?
  # df <- data.frame()
  
  # start while loop here
  fin <- FALSE
  t <- 0
  
  while(fin == FALSE){
    frac <- sum(V(g)$info)/length(V(g))
    # df <- rbind(df, c(t, frac))
    for(i in V(g)[V(g)$info==0]){
      tmp <- V(g)$info[neighbors(g, i, mode='all')]
      p <- sum(tmp)/length(tmp)
      if(is.finite(p) & p >= q){
        V(g)[i]$info <- 1
      }
    }
    frac2 <- sum(V(g)$info)/length(V(g))
    
    t <- t+1
    fin <- frac == frac2
  }
  # end while loop here
  
  dtab <- rbind(dtab, c(q,t,frac,pr))
  
  # colnames(df) <- c('time','freq_a')
  #  V(g)[info==0]$color <- 'blue'
  # V(g)[info==1]$color <- 'red'
  #plot(g, vertex.size=3, vertex.label=NA)
  
  # plot(df$freq_a ~ df$time)
}
colnames(dtab) <- c('threshold_q','time', 'freq_a','pr')

endtime <- Sys.time()

endtime - starttime

ggplot(dtab, aes(x=threshold_q, y=freq_a)) + geom_point()
write.table(dtab, "fbprobcascade_pr001.txt", sep=",")
write.table(dtab, file="fbprobcascade_pr001.csv",sep=",")


# ---- rumor or trend simulation / cascade (degree centrality) --------------------
g <- read_graph('facebook_combined.txt', format='edgelist')
g <- as.undirected(g)

#plot(g, vertex.size=1, vertex.label=NA, 
#     layout=layout.fruchterman.reingold(g))
#plot(g, vertex.size=1, vertex.label=NA, 
 #    layout=layout.kamada.kawai(g))


# set parameters
#dtab <- data.frame()
#set.seed(34564)
#pr <- 0.01  # probability of initiating rumor, story, etc.

starttime <- Sys.time()

qlist <- seq(0,0.5,by=0.01)
for(q in qlist){
  
  V(g)$info <- 0
  V(g)[which.max(degree(g))]$info <- 1
 
  fin <- FALSE
  t <- 0
  
  while(fin == FALSE){
    frac <- sum(V(g)$info)/length(V(g))
    for(i in V(g)[V(g)$info==0]){
      tmp <- V(g)$info[neighbors(g, i, mode='all')]
      p <- sum(tmp)/length(tmp)
      if(is.finite(p) & p >= q){
        V(g)[i]$info <- 1
      }
    }
    frac2 <- sum(V(g)$info)/length(V(g))
    
    t <- t+1
    fin <- frac == frac2
  }
  # end while loop here
  
  dtab <- rbind(dtab, c(q,t,frac,'degree'))
  
  # colnames(df) <- c('time','freq_a')
  #  V(g)[info==0]$color <- 'blue'
  # V(g)[info==1]$color <- 'red'
  #plot(g, vertex.size=3, vertex.label=NA)
  
  # plot(df$freq_a ~ df$time)
}
#colnames(dtab) <- c('threshold_q','time', 'freq_a','pr')

endtime <- Sys.time()

endtime - starttime

dtab$threshold_q <- as.numeric(dtab$threshold_q)
dtab$freq_a <- as.numeric(dtab$freq_a)

ggplot(dtab, aes(x=threshold_q, y=freq_a, colour=pr)) + 
  geom_point() + geom_line()
#write.table(dtab, "fbprobcascade2.txt", sep=",")
#write.table(dtab, file="fbprobcascade2.csv",sep=",")

#
# -------- Cascade capacity -----------
finfr <- 1 # finished frequency always starts at 1
q <- 0 # q (threshold) starts at 0, CC at 1

s1 <- Sys.time()

# adopt, len

# func
prlist <- function(x1){
  rnp <- V(g)$info[unlist(x1)]
  pt <- sum(rnp)/length(rnp)
}

while(finfr==1){
  
  V(g)$info <- 0
  V(g)$info[adopt] <- 1
  
  t <- 0
  dfr <- 1
  infr <- sum(V(g)$info)/len
  
  while(dfr > 0){
    fr <- sum(V(g)$info)/len
    uninit <- V(g)[V(g)$info==0]
    tmp <- adjacent_vertices(g, V(g)[V(g)$info==0])
    # p <- rep(NA, length(tmp))
    
    p <- lapply(tmp, prlist)
    p <- unlist(p)
    
    #    for(i in 1:length(p)){
    #     rnp <- V(g)$info[unlist(tmp[i])]
    #    p[i] <- sum(rnp)/length(rnp)
    # }
    
    
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

e <- Sys.time()

colnames(df) <- c('q','time','init_fr', 'end_fr','pr')

e - s1












