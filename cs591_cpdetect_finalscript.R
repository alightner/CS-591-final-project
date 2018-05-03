library(igraph)
library(tidyverse)
library(data.table)
library(knitr)
library(gridExtra)

source('~/Desktop/Docs/Spring 2018 Classes/CS 591/igraph_projects/cs591_cpdetect_funclist.R')

g <- read_graph('facebook_combined.txt', format='edgelist')
g <- as.undirected(g)

# ------------- RD curve plot ------------------------------

#rd <- read_csv('FBrd_data.csv')
#u <- read_csv('FBsetU_final.csv')

ggplot(data.frame(rank=1:length(rd), rd_1=rd, node=u), 
       aes(x=rank, y=rd_1)) + 
  geom_point(size=0.25) + 
  #geom_text(aes(label=node), nudge_x = 0.5, nudge_y = 0.05,size=2.5) +
  geom_line() + 
  #geom_hline(yintercept = beta, linetype=3) + 
  theme_bw() +
  labs(x='Rank', y='RD') + geom_hline(yintercept = 0.75, linetype=3) +
  geom_hline(yintercept = 1, linetype=3)



# ------------- general network property (table) ------------------------------
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


# ------------- centrality nodes and measures (table) ------------------------------
nodecentr <- node.centrality(g)
centr_measure <- centrality.measure(g)
labelrows <- c('degree','eccentricity','closeness','between','pagerank',
               'authority','hub')
df <- data.table(labelrows,nodecentr, centr_measure)

cclist_val1 <- read_csv('cclist_values.csv')
cclist1 <- read_csv('cclist.csv')

nodelist <- unlist(df$`FB (max. node)`)
nlcs <- cclist_val1$x[cclist1$x %in% nodelist]

# c() manual can be re-coded to obviate copying, time permitting...
df$cascade <- c(0.07,0.07,0.07,0.07,0.07,0.1,0.1)

colnames(df) <- c('centrality measure','FB (max. node)', 'FB (values)','cascade capacity')

df_centrality <- kable(df, 
                       caption = 'Table 1. Node(s) in each network with the 
                       highest scores in terms of centrality measures (listed 
                       in the leftmost column). Column labels include network 
                       labels, (max. node) indicates the node with the highest 
                       centrality measure in the network, and (values) indicate 
                       the maximized centrality measure itself for that node. For 
                       entries ending in "...", see the appendix for a complete 
                       list of values.')


# ------- spectral graph properties (plot and eigenvalues) ---------------------

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

# ------------- degree distribution ------------------------------

ggplot(data.frame(x=1:length(degree_distribution(g)), 
                  y=degree_distribution(g)), aes(x=x,y=y)) + 
  geom_point(alpha=0.3) + scale_x_log10() + scale_y_log10()


# ------------- cascade capacity from central nodes? ------------------------------

xc <- as.data.frame(table(as.numeric(unlist(nodecentr))))
centr_nodes <- as.numeric(levels(droplevels(xc$Var1)))

#cc_fromcentrnodes <- cascade_capacity(g)

# ------------- cascade capacity from community nodes? ------------------------------

# who has the most neighbors in different communities?
cl <- cluster_louvain(g)
V(g)$cl <- cl$membership
V(g)$num_cl <- 0  # number of connected communities?
calc_numcl <- function(x){return(length(table(V(g)$cl[neighbors(g,x)])))}
V(g)$num_cl <- unlist(lapply(V(g), calc_numcl))

table(unlist(V(g)$num_cl))
max_fromcccl <- NULL
max_fromcccl[1] <- V(g)[unlist(V(g)$num_cl) == max(unlist(V(g)$num_cl))] # max = 9
max_fromcccl[2:10] <- V(g)[unlist(V(g)$num_cl) == 6]
max_fromcccl[11:27] <- V(g)[unlist(V(g)$num_cl) == 5]

# ------------- cascade capacity from CP structure nodes ------------------------------

cpset <- read_csv('CPSet.csv', col_names = FALSE)
cpset <- cpset$X2
V(g)$cp <- cpset
V(g)$num_cp <- 0  # number of connected communities?
calc_numcp <- function(x){return(length(table(V(g)$cp[neighbors(g,x)])))}
V(g)$num_cp <- unlist(lapply(V(g), calc_numcp))

table(unlist(V(g)$num_cp))
max_fromcccp <- NULL
max_fromcccp[1:4] <- V(g)[unlist(V(g)$num_cp) == max(unlist(V(g)$num_cp))] # max = 5
max_fromcccp[5:37] <- V(g)[unlist(V(g)$num_cp) == 4]

# ------------- run cascade capacity from all nodes ------------------------------
cclist <- c(centr_nodes,max_fromcccl,max_fromcccp)
xc <- as.data.frame(table(as.numeric(unlist(cclist))))
cclist <- as.numeric(levels(droplevels(xc$Var1)))  # removing redundant nodes

# uncomment below to run - takes some time to run at this stage, use read_csv above

#cclist_val <- rep(NA, length(cclist))
#for(i in 1:length(cclist)){
 # cclist_val[i] <- cascade_capacity(g, cclist[i])
#}

#write.table(cclist, file='cclist.csv', sep=',', row.names = FALSE)
#write.table(cclist_val, file='cclist_values.csv', sep=',', row.names = FALSE)


cclist_val1$x[cclist1$x %in% max_fromcccl]
cclist_val1$x[cclist1$x %in% max_fromcccp]


# ------------- core set identification ------------------------------
coreset <- read_csv('CoreSet.csv', col_names = FALSE)
V(g)$color <- 'blue'
V(g)$color[V(g) %in% coreset$X1] <- 'yellow'
E(g)$width <- 0.05
V(g)$frame.color <- NA
plot(g, vertex.size=1, vertex.label=NA,
     layout=layout.fruchterman.reingold(g))


# ------------- CP set identification ------------------------------
cpset <- read_csv('CPSet.csv', col_names = FALSE)
cpset <- cpset$X2
V(g)$cp <- cpset
V(g)$color[V(g)$cp==0] <- 'green'
V(g)$color[V(g)$cp==1] <- 'red'
V(g)$color[V(g)$cp==2] <- 'blue'
V(g)$color[V(g)$cp==3] <- 'pink'
V(g)$color[V(g)$cp==4] <- 'gray'
V(g)$color[V(g)$cp==5] <- 'orange'
V(g)$color[V(g)$cp==6] <- 'black'
V(g)$color[V(g) %in% coreset$X1] <- 'yellow'
E(g)$width <- 0.05
V(g)$frame.color <- NA

plot(g, vertex.size=1, vertex.label=NA,
     layout=layout.fruchterman.reingold(g))

plot(g, vertex.size=1, vertex.label=NA,
     layout=layout.kamada.kawai(g))

plot(g, vertex.size=1, vertex.label=NA,
     layout=layout.lgl(g))

plot(g, vertex.size=1, vertex.label=NA)

# ------------- statistics ------------------------------

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

mean(cclist_val1$x); sd(cclist_val1$x)
median(cclist_val1$x); IQR(cclist_val1$x)


cls <- cclist_val1$x[cclist1$x %in% max_fromcccl]
cps <- cclist_val1$x[cclist1$x %in% max_fromcccp]
range(cls);mean(cls); sd(cls); median(cls); IQR(cls)
range(cps);mean(cps); sd(cps); median(cps); IQR(cps)

ggplot(data.frame(x=c(cls,cps,nlcs),
                  type=c(rep('LC',length(cls)),
                  rep('CPC', length(cps)), 
                  rep('max. centrality', length(nlcs)))),
       aes(x=x, color=factor(type))) + geom_density() + 
  theme_bw() + labs(x='cascade capacity', color='init. adopters')

# ------------- Epidemic model, or rumor transmission ------------------------------

pr_infect <- seq(0,0.2,by=0.01)
s <- Sys.time()
rmt <- rumor_transmit(g, infect = 0.01, tsteps = 100, 108)
dist108 <- V(rmt)$info

rmt <- rumor_transmit(g, infect = 0.01, tsteps = 100, 693)
V(rmt)$info[693]
dist693 <- V(rmt)$info

# LC and CPC = 1
rmt <- rumor_transmit(g, infect = 0.01, tsteps = 100, 18)
V(rmt)$info[18]
dist18 <- V(rmt)$info

ggplot(data.frame(x=cl$membership, y=dist108), aes(x=factor(x), y=y)) + 
  geom_boxplot() + theme_bw()+ 
  geom_hline(yintercept = median(dist108), linetype=3) +
  labs(x='community',y='times heard') + geom_vline(xintercept = 8, linetype=3)#+ 
 # geom_vline(xintercept = 8.5, linetype=3)

p1 <- ggplot(data.frame(x=c(dist108,dist18,dist693),
                  type=c(rep(108,length(dist108)),
                         rep(18, length(dist18)),
                         rep(693, length(dist693)))),
       aes(x=x, color=factor(type))) + geom_density() +
  theme_bw() + labs(x='number of times heard', color='node') +
  ggtitle('Pr=0.01')

p2 <- ggplot(data.frame(x=c(dist108,dist18,dist693),
                        type=c(rep(108,length(dist108)),
                               rep(18, length(dist18)),
                               rep(693, length(dist693)))),
             aes(x=x, color=factor(type))) + geom_density() +
  theme_bw() + labs(x='number of times heard', color='node') +
  ggtitle('Pr=0.05')

p3 <- ggplot(data.frame(x=c(dist108,dist18,dist693),
                        type=c(rep(108,length(dist108)),
                               rep(18, length(dist18)),
                               rep(693, length(dist693)))),
             aes(x=x, color=factor(type))) + geom_density() +
  theme_bw() + labs(x='number of times heard', color='node') +
  ggtitle('Pr=0.10')

grid.arrange(p1,p2,p3)


# coloring nodes by salience
#V(g)$info <- dist18
#V(g)$info <- dist693
#V(g)$info <- dist108
E(g)$width <- 0.05
V(g)$frame.color <- NA
heard <- V(g)$info
fine <- 500 # this will adjust the resolving power
pal <- colorRampPalette(c('black','yellow'))
graphCol = pal(fine)[as.numeric(cut(heard,breaks = fine))]

plot(g, vertex.size=2, vertex.label=NA,
     layout=layout.fruchterman.reingold(g),vertex.color=graphCol)

plot(g, vertex.size=2, vertex.label=NA,
     layout=layout.lgl(g),vertex.color=graphCol)

clcent <- closeness(g)
degcent <- degree(g)
prank <- page_rank(g)
prank <- prank$vector
#ggplot(data.frame(x=clcent, y=dist108), aes(x=x,y=y)) + geom_point()
#ggplot(data.frame(x=degcent, y=dist108), aes(x=x,y=y)) + geom_point()
#ggplot(data.frame(x=prank, y=dist108), aes(x=x,y=y)) + geom_point()


