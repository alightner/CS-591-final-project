library(igraph)
library(tidyverse)

g <- read_graph('facebook_combined.txt', format='edgelist')
g <- as.undirected(g)

# set parameters
alpha <- round(mean(degree(g)))
beta <- 0.75

u <- read_csv('FBsetU3395_2.csv')
u <- as.vector(u$x)

#V(g)$tag <- 0
#u <- which(closeness(g) == max(closeness(g)))
#V(g)$tag[u[1]] <- 1
#u[3396:4039] <- NA

tmp <- NULL
for(i in u){
  tmp <- c(tmp, neighbors(g, V(g)[i]))
}

# tmp <- rep(NA, (ecount(g)*2))
#ct <- 1
#tmp <- c(tmp,neighbors(g, V(g)[i]))
#length(subset(tmp, tmp %in% V(g)[V(g)$tag==0]))
tmp <- subset(tmp, tmp %in% V(g)[V(g)$tag==0])

for(i in 3397:vcount(g)){
  
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

FBgraph_cp <- ggplot(data.frame(rank=1:length(rd), rd_1=rd, node=u), 
                     aes(x=rank, y=rd_1)) + 
  geom_point() + 
  #geom_text(aes(label=node), nudge_x = 0.5, nudge_y = 0.05,size=2.5) +
  geom_line() + 
  #geom_hline(yintercept = beta, linetype=3) + 
  theme_bw() +
  labs(x='Rank', y='RD')

write.table(u, file="FBsetU3395_final.csv",sep=",", row.names = FALSE)
write.table(u, file="FBsetU3395_final.txt", sep=",", row.names = FALSE)
write.table(rd, file='FBrd_data.csv', sep=',', row.names = FALSE)
write.table(rd, file='FBrd_data.txt', sep=',', row.names = FALSE)

ggsave('FBgraph_cp3.jpg',FBgraph_cp)

V(g)$cp[u] <- rd
V(g)$color[V(g)$cp >= beta] <- 'red'
V(g)$color[V(g)$cp < beta] <- 'blue'

plot(g, vertex.size=5)


# ------- algo. 3 and 4 ---------

cset <- NULL
numcore <- 1

for(i in alpha:vcount(g)){
  if(rd[u[i]] >= beta){
    cset <- c(cset, u[(i-alpha+1):i])
    if((rd[u[i-1]] >= beta) & (i > alpha)){
      numcore <- numcore + 1
    }
  }
}

if(rd[u[vcount(g)]] < beta){
  numcore <- numcore - 1
}

cset
numcore

  
  
  