library(igraph)
library(tidyverse)
library(data.table)

#
# -------- Xiang et al. (2018) RD calculation -----------

#g <- read_graph('karate/karate.gml', format='gml')
#g <- as.undirected(g)

g <- read_graph('facebook_combined.txt', format='edgelist')
g <- as.undirected(g)

# set parameters
alpha <- floor(mean(degree(g)))
beta <- 0.75

u <- rep(NA, vcount(g))
V(g)$tag <- 0

u[1] <- which(closeness(g) == max(closeness(g)))
V(g)$tag[u[1]] <- 1

tmp <- NULL
for(i in na.omit(u)){
  tmp <- c(tmp, neighbors(g, V(g)[i]))
}

tmp <- subset(tmp, tmp %in% V(g)[V(g)$tag==0])

s1 <- Sys.time()

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

#Time difference of 47.34357 secs on FB network

ggplot(data.frame(rank=1:length(rd), rd_1=rd, node=u), 
       aes(x=rank, y=rd_1)) + 
  geom_point(size=0.25) + 
  #geom_text(aes(label=node), nudge_x = 0.5, nudge_y = 0.05,size=2.5) +
  geom_line() + 
  #geom_hline(yintercept = beta, linetype=3) + 
  theme_bw() +
  labs(x='Rank', y='RD') + geom_hline(yintercept = 0.85, linetype=3) +
  geom_hline(yintercept = 1, linetype=3)





#
# -------- Cascade capacity -----------

#g <- read_graph('karate/karate.gml', format='gml')
#g <- as.undirected(g)

g <- read_graph('facebook_combined.txt', format='edgelist')
g <- as.undirected(g)

finfr <- 1 # finished frequency always starts at 1
q <- 0 # q (threshold) starts at 0, CC at 1

# probability of initiated (behavior A) neighbors function
prlist <- function(x1){
  rnp <- V(g)$info[unlist(x1)]
  pt <- sum(rnp)/length(rnp)
}

df <- data.frame(matrix(nrow=2,ncol=5))
len <- length(V(g))

state <- 'max_PR'
page <- page_rank(g)
adopt <- c(which.max(page$vector), which.max(degree(g)))
#adopt <- which(page$vector >= 0.004)
#adopt <- which.max(closeness(g))
#adopt <- which.min(closeness(g))
#adopt <- which.max(degree(g)) == max close
#adopt <- which.max(betweenness(g)) == max between

s1 <- Sys.time()

while(finfr==1){
  
  V(g)$info <- 0
  V(g)$info[adopt] <- 1
  
  t <- 0
  dfr <- 1
  infr <- sum(V(g)$info)/len
  
  while(dfr > 0){
    fr <- sum(V(g)$info)/len
    uninit <- V(g)[V(g)$info==0]  # uninitiated nodes
    tmp <- adjacent_vertices(g, uninit)
    p <- unlist(lapply(tmp, prlist))
    V(g)$info[uninit[(is.finite(p)) & (p >= q)]] <- 1
    t <- t+1
    dfr <- (sum(V(g)$info)/len) - fr
  }
  
  finfr <- fr
  
  df[1,] <- df[2,]
  df[2,] <- c(q,t,infr,fr,state)
  #df <- rbind(df, c(q,t,infr,fr,'max_PR'))
  q <- q+0.01
  
}

e <- Sys.time()

colnames(df) <- c('q','time','init_fr', 'end_fr','pr')

e - s1

#
# -------- Epidemiological model idea transmission -----------

g <- read_graph('karate/karate.gml', format='gml')
g <- as.undirected(g)

#g <- read_graph('facebook_combined.txt', format='edgelist')
#g <- as.undirected(g)

# set parameters
infect <- 0.1
tsteps <- 100

state <- 'maxpagerank'

if(state=='maxpagerank'){
  page <- page_rank(g)
  adopt <- c(which.max(page$vector), which.max(degree(g)))
}
if(state=='maxdegree'){
  adopt <- which.max(degree(g))
}

V(g)$info <- 0
V(g)$info[adopt] <- 1

s1 <- Sys.time()

for(t in 1:tsteps){
  tmp <- unlist(adjacent_vertices(g, V(g)[V(g)$info >= 1]))
  prst <- sample(c(0,1), length(tmp), 
                 replace=TRUE, prob =c((1-infect), infect))
  V(g)$info[tmp[prst==1]] <- V(g)$info[tmp[prst==1]] + 1
}

e <- Sys.time()

e-s1

# coloring nodes by salience
E(g)$width <- 0.5
V(g)$frame.color <- NA
heard <- V(g)$info
fine <- 500 # this will adjust the resolving power
pal <- colorRampPalette(c('black','white'))
graphCol = pal(fine)[as.numeric(cut(btw,breaks = fine))]
plot(g, vertex.size=5, vertex.label=NA,
     layout=layout.fruchterman.reingold(g),vertex.color=graphCol)







