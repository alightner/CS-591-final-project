library(tidyverse)
library(igraph)
library(data.table)

g <- read_graph('facebook_combined.txt', format='edgelist')
g <- as.undirected(g)

# ----- Xiang et al. algorithm -----
alpha <- round(mean(degree(g)))
beta <- 0.75
#fun_cp on g after consolidating in fun_script
u <- read_csv('FBsetU_final.csv')
u <- as.vector(u$x)

rd <- read_csv('FBrd_data.csv')
rd <- as.vector(rd$x)

fbCPplot1 <- ggplot(data.frame(rank=1:length(rd), rd_1=rd, node=u), 
       aes(x=rank, y=rd_1)) + 
  geom_point(size=0.25) + 
  #geom_text(aes(label=node), nudge_x = 0.5, nudge_y = 0.05,size=2.5) +
  geom_line() + 
  #geom_hline(yintercept = beta, linetype=3) + 
  theme_bw() +
  labs(x='Rank', y='RD') + geom_hline(yintercept = 0.85, linetype=3) +
  geom_hline(yintercept = 1, linetype=3)

ggsave('fbCPplot1.jpg', fbCPplot1)

# ----- Centrality set -----

page <- page_rank(g)
which.max(page$vector)
which(page$vector >= 0.002)
plot(density(page$vector))


# ----- Community coloring and graph visual -----

g <- read_graph('facebook_combined.txt', format='edgelist')
g <- as.undirected(g)

#cm <- cluster_label_prop(g)
cm <- cluster_louvain(g)
#cm <- fastgreedy.community(g)
V(g)$color <- cm$membership
V(g)$frame.color <- NA
E(g)$width <- 0.2
weights <- ifelse(crossing(cm, g), 1, 100)
layout <- layout_with_lgl(g)#, weights=weights)
plot(g, vertex.size=1, vertex.label=NA, layout=layout)

#, layout=layout)

set.seed(192837)
g <- grg.game(100, 0.2)
E(g)$width <- 0.2
plot(g, vertex.size=2, vertex.label=NA)
#g <- erdos.renyi.game(100,0.1)
#cl <- fastgreedy.community(g)
cl <- cluster_louvain(g)
weights <- ifelse(crossing(cl, g), 1, 100)
layout <- layout_with_kk(g, weights=weights)
plot(g, vertex.size=1, vertex.label=NA, 
     vertex.color=cl$membership, layout=layout)

cl <- cluster_louvain(g)
members_cl <- ggplot(data.frame(x=cl$membership), aes(x=x)) + geom_histogram(binwidth = 1) +
  theme_bw() + labs(x='membership') 

ggsave('members_cl.jpg', members_cl)

# ----- Cascade capacity -----

# probabilistic approach (test)

g <- read_graph('karate/karate.gml', format='gml')
# set parameters
E(g)$width <- 0.5
V(g)$frame.color <- NA
plot(g, vertex.size=5, vertex.label=NA, 
     vertex.color='black',edge.color='blue')

pr <- 0.1  # probability of initiating rumor, story, etc.
qlist <- seq(0,0.5,by=0.01)
#df <- data.table()
df <- data.frame()
len <- length(V(g))

starttime <- Sys.time()

trial <- round(runif(1,min=1,max=100000))

# start while here
finfr <- 1
while(finfr==1){
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
  
  df <- rbind(df, c(q,t,infr,fr,pr))
  q <- q+0.01
  finfr <- fr
}

endtime <- Sys.time()

colnames(df) <- c('q','time','init_fr', 'end_fr','pr')

endtime - starttime

#Time difference of 1.627798 secs! add multiple seed trials?

ggplot(df, aes(x=q, y=end_fr)) + geom_point() + geom_smooth(method='lm') +
  geom_jitter()


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

# ----- Epidemiological -----

g <- read_graph('facebook_combined.txt', format='edgelist')
g <- as.undirected(g)

E(g)$width <- 0.1
V(g)$frame.color <- NA
plot(g, vertex.size=1, vertex.label=NA, 
     vertex.color='black',edge.color='blue', layout=layout.lgl(g))
V(g)$info <- 0

# set parameters
pr <- 0.01
infect <- 0.01
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
E(g)$width <- 0.1
heard <- V(g)$info
fine <- 500 # this will adjust the resolving power
pal <- colorRampPalette(c('black','white'))
graphCol = pal(fine)[as.numeric(cut(btw,breaks = fine))]
plot(g, vertex.size=1, vertex.label=NA,
     layout=layout.lgl(g),edge.color='blue',vertex.color=graphCol)

tt1 <- table(V(g)$info)
pp1 <- ggplot(data.frame(heard=V(g)$info, deg=degree(g), 
                  close=closeness(g), member=cluster_louvain(g)$membership),
       aes(x=heard)) + geom_histogram(binwidth=5)

pp2 <- ggplot(data.frame(heard=V(g)$info, deg=degree(g), 
                  close=closeness(g), member=cluster_louvain(g)$membership),
       aes(x=heard, y=deg)) + geom_point(alpha=0.1) #+ scale_x_log10() + scale_y_log10()

pp3 <- ggplot(data.frame(heard=V(g)$info, deg=degree(g), 
                  close=page_rank(g)$vector, member=cluster_louvain(g)$membership),
       aes(x=close, y=heard)) + geom_point(alpha=0.1) + scale_x_log10() #+ scale_y_log10()

pp4 <- ggplot(data.frame(heard=V(g)$info, deg=degree(g), 
                  close=closeness(g), member=cluster_louvain(g)$membership),
       aes(x=factor(member), y=heard)) + geom_violin()


tt2 <- table(V(g)$info)
pp5 <- ggplot(data.frame(heard=V(g)$info, deg=degree(g), 
                         close=closeness(g), member=cluster_louvain(g)$membership),
              aes(x=heard)) + geom_histogram(binwidth=5)

pp6 <- ggplot(data.frame(heard=V(g)$info, deg=degree(g), 
                         close=closeness(g), member=cluster_louvain(g)$membership),
              aes(x=heard, y=deg)) + geom_point(alpha=0.1) #+ scale_x_log10() + scale_y_log10()

pp7 <- ggplot(data.frame(heard=V(g)$info, deg=degree(g), 
                         close=page_rank(g)$vector, member=cluster_louvain(g)$membership),
              aes(x=close, y=heard)) + geom_point(alpha=0.1) + scale_x_log10() #+ scale_y_log10()

pp8 <- ggplot(data.frame(heard=V(g)$info, deg=degree(g), 
                         close=closeness(g), member=cluster_louvain(g)$membership),
              aes(x=factor(member), y=heard)) + geom_violin()

PR1 <- pp3 + labs(x='PageRank') + theme_bw()
PR2 <- pp7 + labs(x='PageRank') + theme_bw()
ggsave('heard_pr1.jpg', PR1)
ggsave('heard_pr01.jpg', PR2)


plot(g, vertex.size=1, vertex.label=NA,
     layout=layout.fruchterman.reingold(g),vertex.color=graphCol)

plot(g, vertex.size=1, vertex.label=NA,
     layout=layout.kamada.kawai(g),vertex.color=graphCol)


