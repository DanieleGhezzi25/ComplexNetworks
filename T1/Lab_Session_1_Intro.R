library(igraph)
library(ggplot2)
library(RColorBrewer)
library(rgl)


################
# CREATE NETWORKS
################

# creating a simple empty network by scratch

g <- make_empty_graph()
summary(g)

# creating a simple directed network with a few edges by scratch
g <- make_graph(edges = c(1,2, 3,4, 5,6), n=8, directed = TRUE)
summary(g)

# check if it is directed
is.directed(g)


# creating a simple undirected network with a few edges by scratch: way 1
g <- make_graph(edges = c(1,2, 3,4, 5,6), n=8, directed = FALSE)
summary(g)

# check if it is directed
is.directed(g)



# creating a simple undirected network with a few edges by scratch: way 2
g <- make_graph(~ 1--2, 2--3, 3--4, 4--5, 5--6, 6--1, 4, 5, 6, 7, 8)
summary(g)


# see edge list
get.edgelist(g)

# see adj matrix (default: sparse matrix format)
get.adjacency(g)

# see thecombinatorial Laplacian matrix (default: sparse matrix format)
laplacian_matrix(g)

# see the normalized Laplacian (ie, the one describing a RW) (default: sparse matrix format)
laplacian_matrix(g, normalized=TRUE)

# check if the latest sums to zero along -> error because of the matrix format
colSums(laplacian_matrix(g, normalized=TRUE, sparse=T))

# check if the latest sums to zero along -> no error
colSums(laplacian_matrix(g, normalized=TRUE, sparse=F))

################
# NODES
################

# view and count nodes
V(g)
vcount(g)

# node attributes
V(g)$name
V(g)$label
V(g)$label <- letters[1:vcount(g)]

V(g)$color
V(g)$color <- rainbow(vcount(g))

# add nodes
g <- add_vertices(g, 3)

# check their names
V(g)$name

# fix their names
V(g)$name <- 1:vcount(g)

# delete a node with invalid ID
g <- delete_vertices(g, 10)

# delete a node with valid ID
g <- delete_vertices(g, 8)


################
# EDGES
################

# view and count edges
E(g)
ecount(g)

# check their list
E(g)

# get their weights
E(g)$weight 

# set their weights
E(g)$weight <- rpois(ecount(g),2)+1


# add edges
g <- add_edges(g, edges=c(2,7, 3,1))

# error if you have invalid node IDs
g + edges(c(45,46))

# delete an edge 
idx <- get.edge.ids(g, c(2,7))
g <- delete_edges(g, idx)



##########################
# CREATE FROM DATA & PLOT
##########################

set.seed(12345)
my_edges <- data.frame(from=sample(letters, 10), to=sample(letters, 10), weight=runif(10))
my_edges

# creating a simple graph from an edges list
g <- graph_from_data_frame(my_edges, directed=TRUE)

# same, but customizing the ordering of nodes
nodes_ordering <- sort(union(my_edges$from, my_edges$to))
nodes_ordering
g <- graph_from_data_frame(my_edges, directed=TRUE, vertices=nodes_ordering)

plot(g)


# more complex
n <- 200
my_edges <- data.frame(from=sample(letters, n, replace=TRUE), to=sample(letters, n, replace=TRUE), weight=runif(n))

g <- graph_from_data_frame(my_edges, directed=TRUE)
plot(g)

# use also weights and change arrow size
plot(g, edge.width=E(g)$weight, edge.arrow.size=0.4)


##########################
# LAYOUT
##########################

g <- make_graph('Zachary')
plot(g)

# remove labels and change node colors
plot(g, vertex.color="steelblue", vertex.label=NA)

# remove labels and change node colors, curve edges a little
plot(g, vertex.color="steelblue", vertex.label=NA, edge.curved=0.2)

lay <- layout_in_circle(g)
plot(g, layout=lay, vertex.color="steelblue", vertex.label=NA, edge.curved=0.2)

lay <- layout_with_fr(g)
plot(g, layout=lay, vertex.color="steelblue", vertex.label=NA, edge.curved=0.2)

lay <- layout_with_kk(g)
plot(g, layout=lay, vertex.color="steelblue", vertex.label=NA, edge.curved=0.2)



##########################
# BASIC MODELS
##########################

N <- 100

# ERDOS-RENYI
# at the percolation point
g <- erdos.renyi.game(N, p=1/N, directed=F)
lay <- layout_with_fr(g)
plot(g, layout=lay, vertex.color="steelblue", vertex.label=NA, edge.curved=0.2)

# fix the viz
plot(g, layout=lay, vertex.color="steelblue", vertex.label=NA, edge.color="gray80", vertex.size=2.5, vertex.frame.color=NA)


# above the percolation point
g <- erdos.renyi.game(N, p=log(N)/N, directed=F)
lay <- layout_with_fr(g)
plot(g, layout=lay, vertex.color="steelblue", vertex.label=NA, edge.color="gray80", vertex.size=2.5, vertex.frame.color=NA)



# GEOMETRIC RANDOM GRAPH
g <- sample_grg(N, radius=0.2)
lay <- layout_with_fr(g)
plot(g, layout=lay, vertex.color="steelblue", vertex.label=NA, edge.color="gray80", vertex.size=2.5, vertex.frame.color=NA)



# BARABASI-ALBERT
N <- 500
g <- barabasi.game(N, m=1, directed=FALSE)
lay <- layout_with_fr(g)
plot(g, layout=lay, vertex.color="steelblue", vertex.label=NA, edge.color="gray80", vertex.size=2.5, vertex.frame.color=NA)

# highlight degree of nodes
plot(g, vertex.color="steelblue", vertex.label=NA, edge.curved=0.2, layout=lay, vertex.size=4*log10(degree(g)), edge.color="gray80", vertex.frame.color=NA)


# SBM with 4 groups
N <- 1000
OM <- matrix(0.0003, ncol=4, nrow=4)
diag(OM) <- 0.03
g <- sample_sbm(N, pref.matrix=OM, block.sizes=c(200, 200, 300, 300))
lay <- layout_with_fr(g)
plot(g, layout=lay, vertex.color="steelblue", vertex.label=NA, edge.color="gray80", vertex.size=2.5, vertex.frame.color=NA)

# more complex
OM <- matrix(runif(4*4,0.3/N,0.5/N), ncol=4, nrow=4)
diag(OM) <- runif(4, 3*log(N)/N, 5*log(N)/N)
# symmetrize it
OM <- (OM + t(OM))/2
g <- sample_sbm(N, pref.matrix=OM, block.sizes=c(100, 200, 300, 400))
lay <- layout_with_fr(g)
plot(g, layout=lay, vertex.color="steelblue", vertex.label=NA, edge.color="gray80", vertex.size=2.5, vertex.frame.color=NA)


