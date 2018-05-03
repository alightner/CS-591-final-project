library(igraph)
library(tidyverse)

# ---------- CP structures (Xiang et al. 2018) ---------------
# karate club test data
g <- read_graph('karate/karate.gml', format='gml')
#g <- read_graph('dolphins/dolphins.gml', format='gml')
#g <- read_graph('polblogs/polblogs.gml', format='gml')
g <- read_graph('facebook_combined.txt', format='edgelist')
#plot(g, vertex.size=3)
g <- as.undirected(g)

# set parameters
alpha <- round(mean(degree(g)))
#alpha <- 5
beta <- 0.75
#randomtick <- 0

s1 <- Sys.time()

V(g)$tag <- 0
u <- which(closeness(g) == max(closeness(g)))
V(g)$tag[u[1]] <- 1
for(i in 2:vcount(g)){
  tmp <- NULL
  for(k in u){
    tmp <- c(tmp, neighbors(g, V(g)[k]))
  }
  tmp <- subset(tmp, tmp %in% V(g)[V(g)$tag==0])
  #neighbors(g, V(g)[k])[V(g)$tag==1]
  xc <- as.data.frame(table(tmp))
  cand <- as.numeric(levels(droplevels(xc$tmp[xc$Freq==max(xc$Freq)])))
  # test case, add 24 and 47
  #cand <- c(24, 25, 47)
  if(length(cand) > 1){
    degsmp <- degree(g)[cand]
    cand <- cand[which.max(degsmp)]
    if(length(cand) > 1){
      cand <- sample(cand, 1)
      #randomtick <- randomtick + 1
    }
  }
  u[i] <- cand
  V(g)$tag[cand] <- 1
}
e <- Sys.time()

rd <- 0
for(i in 2:length(u)){
  lower <- i - (alpha-1)
  if(lower < 1){
    lower <- 1
  }
  s <- induced_subgraph(g, u[lower:i])
  rd[i] <- (2*ecount(s))/(vcount(s)*(vcount(s)-1))
}
e - s1

karate_rd <- ggplot(data.frame(rank=1:length(rd), rd_1=rd, node=u), 
       aes(x=rank, y=rd_1)) + 
  geom_point() + 
  #geom_text(aes(label=node), nudge_x = 0.5, nudge_y = 0.05,size=2.5) +
  geom_line() + 
  #geom_hline(yintercept = beta, linetype=3) + 
  theme_bw() +
  labs(x='Rank', y='RD') + geom_hline(yintercept = 1, linetype=3)

write.table(u, file="FBsetU3395_2.csv",sep=",", row.names = FALSE)
write.table(u, file="FBsetU3395_2.txt", sep=",", row.names = FALSE)
ggsave('karate_rd.jpg',karate_rd)

V(g)$cp[u] <- rd
V(g)$color[V(g)$cp >= beta] <- 'red'
V(g)$color[V(g)$cp < beta] <- 'blue'

plot(g, vertex.size=5)

# ---------- communities ------------

louvain <- cluster_louvain(g)
labelprop <- cluster_louvain(g)

# ---- rumor or trend simulation / cascade (probabilistic approach) --------------------
g <- read_graph('karate/karate.gml', format='gml')
# set parameters
set.seed(34564)
E(g)$width <- 0.5
V(g)$frame.color <- NA
plot(g, vertex.size=5, vertex.label=NA, 
     vertex.color='black',edge.color='blue')
V(g)$info <- 0

pr <- 0.1  # probability of initiating rumor, story, etc.
qlist <- seq(0,0.5,by=0.01)
df <- data.frame()

# qlist here

V(g)$info <- sample(c(0,1), size=length(V(g)), 
                    replace=TRUE, prob=c((1-pr), pr))

fr <- sum(V(g)$info)/length(V(g))


for(i in V(g)[V(g)$info==0]){
  tmp <- V(g)$info[neighbors(g, i, mode='all')]
  p <- sum(tmp)/length(tmp)
  if(is.finite(p) & p >= q){
    V(g)[i]$info <- 1
  }
}



for(q in qlist){
  
  
  
  
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


# ---- rumor or trend simulation / cascade (set nodes approach) --------------------

g <- read_graph('karate/karate.gml', format='gml')
# set parameters
qlist <- seq(0,0.5,by=0.01)

stime <- Sys.time()

dtab <- data.frame()
for(q in qlist){
  for(nd in 1:vcount(g)){
    V(g)$info <- 0
    V(g)[nd]$info <- 1
    
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
    dtab <- rbind(dtab, c(q,t,frac,nd))
  }
  
  
  
  
  
  
  
  # colnames(df) <- c('time','freq_a')
  #  V(g)[info==0]$color <- 'blue'
  # V(g)[info==1]$color <- 'red'
  #plot(g, vertex.size=3, vertex.label=NA)
  
  # plot(df$freq_a ~ df$time)
  
}
colnames(dtab) <- c('threshold_q','time', 'freq_a','start_node')

endtime <- Sys.time()

endtime - stime  # Time difference of 2.457374 mins

# ---------- adjacency matrix plot ----------------
#amg <- get.adjacency(g)
library(reshape2)

Am <- as.matrix(as_adjacency_matrix(g, type='both'))
#Am <- as.matrix(g[u])
longData<-melt(Am)

adjm_plot <- ggplot(longData, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="white", high="black") +
  #scale_y_reverse() + 
  theme_bw() + guides(fill=FALSE)

# ------- epidemic models for salience (probabilistic) ------------
# set parameters
pr <- 0.1
# karate club test data
g <- read_graph('karate/karate.gml', format='gml')
plot(g, vertex.size=3)

infect <- 0.95
max_salience <- 10

V(g)$info <- sample(c(0,1), size=length(V(g)), 
                    replace=TRUE, prob=c((1-pr), pr))
V(g)[info==0]$color <- 'blue'
V(g)[info==1]$color <- 'red'

plot(g, vertex.size=3, vertex.label=NA)
df <- data.frame()

while(max(V(g)$info) < max_salience){
  for(i in V(g)[V(g)$info >= 1]){
    tmp <- V(g)$info[neighbors(g, i, mode='all')]
    for(k in 1:length(tmp)){
      if((sample(c(0,1), 1, prob =c((1-infect), infect)))==1){
        V(g)$info[neighbors(g, i, mode='all')][k] <- 
          V(g)$info[neighbors(g, i, mode='all')][k] + 1
      }
    }
  } 
}

V(g)[info==0]$color <- 'blue'
V(g)[info > 1]$color <- 'red'

plot(g, vertex.size=3, vertex.label=NA)

#sum(V(g)$info)/length(V(g))
plot(density(V(g)$info))
#max(V(g)$info)
#V(g)$info
#  }

length(V(g)$info[V(g)$info > 0])/
  length(V(g)$info)

table(V(g)$info)

#median(V(g)$info)  
#IQR(V(g)$info)

ggplot(data.frame(v=as.vector(V(g)),deg=degree(g),salience=V(g)$info),
       aes(x=deg,y=salience)) + geom_point()


# -------- data analyses and experimentation ----------------

V(g)$color[u[1:5]] <- 'red'
V(g)$color[u[6]] <- 'blue'
V(g)$color[u[7:10]] <- 'red'
V(g)$color[u[11:length(u)]] <- 'blue'

plot(g, vertex.size=5, vertex.label=V(g)$info)

V(g)$color[V(g) %in% cset] <- 'yellow'
V(g)$color[!(V(g) %in% cset)] <- 'black'
plot(g, vertex.size=3, vertex.label=NA, layout=layout_nicely(g))







