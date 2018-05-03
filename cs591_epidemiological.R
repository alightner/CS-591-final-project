library(tidyverse)
library(igraph)

#g <- read_graph('facebook_combined.txt', format='edgelist')
#g <- as.undirected(g)

# test case
set.seed(192837)
g <- erdos.renyi.game(100,0.1)
E(g)$width <- 0.5
V(g)$frame.color <- NA
plot(g, vertex.size=5, vertex.label=NA, 
     vertex.color='black',edge.color='blue')
V(g)$info <- 0

# set parameters
pr <- 0.1

infect <- 0.05
tsteps <- 100

V(g)$info <- sample(c(0,1), size=length(V(g)), 
                    replace=TRUE, prob=c((1-pr), pr))
V(g)[info==0]$color <- 'black'
V(g)[info==1]$color <- 'yellow'

plot(g, vertex.size=5, vertex.label=NA, 
     edge.color='blue', layout=layout.kamada.kawai(g))

for(t in 1:tsteps){
  tmp <- unlist(adjacent_vertices(g, V(g)[V(g)$info >= 1]))
  prst <- sample(c(0,1), length(tmp), 
                 replace=TRUE, prob =c((1-infect), infect))
  V(g)$info[tmp[prst==1]] <- V(g)$info[tmp[prst==1]] + 1
}

# coloring nodes by salience
heard <- V(g)$info
fine <- 500 # this will adjust the resolving power
pal <- colorRampPalette(c('darkblue','white'))
graphCol = pal(fine)[as.numeric(cut(btw,breaks = fine))]
plot(g, vertex.size=5, vertex.label=NA,
     vertex.color=graphCol)




#V(g)[info==0]$color <- 'blue'
#V(g)[info > 1]$color <- 'red'

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

# coloring nodes by salience
heard <- V(g)$info

fine <- 500 # this will adjust the resolving power.
pal <- colorRampPalette(c('black','skyblue'))

#this gives you the colors you want for every point
graphCol = pal(fine)[as.numeric(cut(btw,breaks = fine))]

# now you just need to plot it with those colors
plot(adj, vertex.size=5, vertex.label=NA,
     vertex.color=graphCol)

ggplot(data.frame(x=degree(g), y=heard), aes(x=x,y=y)) + 
  geom_point()

ggplot(data.frame(x=page_rank(g)$vector, y=heard), aes(x=x,y=y)) + 
  geom_point() + geom_smooth(method='lm')

dff <- data.frame(x=page_rank(g)$vector, y=heard)

m <- lm(y ~ x, data = dff)
summary(m)

# -------- color gradient code ----------------

#Some sample data
x <- runif(100)
dat <- data.frame(x = x,y = x^2 + 1)

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))

#This adds a column of color values
# based on the y values
dat$Col <- rbPal(10)[as.numeric(cut(dat$y,breaks = 10))]

plot(dat$x,dat$y,pch = 20,col = dat$Col)

# example 2

btw <- degree(g)

fine <- 500 # this will adjust the resolving power.
pal <- colorRampPalette(c('black','skyblue'))

#this gives you the colors you want for every point
graphCol = pal(fine)[as.numeric(cut(btw,breaks = fine))]

# now you just need to plot it with those colors
plot(adj, vertex.size=5, vertex.label=NA,
     vertex.color=graphCol)

