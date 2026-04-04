library(igraph) # is written in C and Fortran
library(ggplot2)
library(RColorBrewer)
library(rgl)

g <- make_empty_graph()
g
summary(g)

g <- make_graph(edges = c(1,2, 3,4, 5,6), n=8, directed=TRUE)
g # 1 to 2, 3 to 4, 5 to 6, while 7 and 8 are isolated vertices
summary(g)
is_directed(g)

# NOTE: name the nodes ALWAYS starting from 1 (not 0), ALWAYS use integers !!!

g <- make_graph(edges = c(1,2, 3,4, 5,6), n=8, directed=FALSE)
g # 1--2, 3--4, 5--6 (undirected), while 7 and 8 are isolated vertices
summary(g)
is_directed(g)

g <- make_graph(~ 1--2, 2--3, 3--4, 4--5, 5--6, 6--1, 4, 5, 6, 7, 8)
g 
get.edgelist(g)
as_edgelist(g) # is a DataFrame
get.adjacency(g) # sparse matrix format (igraph works with this format)

laplacian_matrix(g) # equivalent to L = D^-1 A
laplacian_matrix(g, normalized=TRUE)
laplacian_matrix(g, normalized=TRUE, sparse=FALSE) # write all zeros
colSums(laplacian_matrix(g, normalized=TRUE, sparse=FALSE)) # sparse=TRUE does not work with colSums !


V(g) # vertices
vcount(g) # number of vertices
V(g)$name # names of vertices
V(g)$label # labels of vertices (metadata)
letters # built-in R vector of letters
V(g)$label <- letters[1:vcount(g)] # assign names to label's vertices
V(g)$color
rainbow(n=4) # first 4 colors of the rainbow (built-in R function)
V(g)$color <- rainbow(n=vcount(g)) # assign different colors to vertices
g # now we have more attributes (name, label, color)

# how to add nodes?
g <- add_vertices(g, 3)
g
V(g)$name # we have NA names for new nodes
V(g)$name <- 1:vcount(g)
V(g)$name

g <- delete_vertices(g, 10)
g
V(g)$name

# for edges
E(g)
ecount(g) # number of edges
E(g)$weight # weights of edges
E(g)$weight <- rpois(n=ecount(g), lambda=2) + 1 # random Poisson weights for edges + 1
E(g)$weight
# NOTE: we want to avoid 0 weights, because they are equivalent to no edge at all
# to avoid this, we sum 1 to the random weights
g <- add_edges(g, edges=c(2,7, 3,1))
# or:  g + edges(c(2,7, 3,1))

g + edges(c(45,46)) # gives an error

idx <- get_edge_ids(g, c(2,7))
g <- delete_edges(g, idx)
g


# network of random list of random letters
set.seed(12345)
letters
sample(letters, size=10)
# dataframe with 3 columns and 10 rows
my_edges <- data.frame(from=sample(letters, 10), to=sample(letters, 10), weight=runif(10))
my_edges
g <- graph_from_data_frame(my_edges, directed=TRUE)
g
V(g)$name
nodes_ordering <- sort(union(my_edges$from, my_edges$to)) # union gives the same V(g)$name object
g <- graph_from_data_frame(my_edges, directed=TRUE, vertices=data.frame(name=nodes_ordering)) 
V(g)$name
plot(g)

n <- 200
my_edges <- data.frame(from=sample(letters, n, replace=TRUE), to=sample(letters, n, replace=TRUE), weight=runif(n))
my_edges
g <- graph_from_data_frame(my_edges, directed=TRUE)
plot(g, edge.width=E(g)$weight, edge.arrow.size=0.4)


# Zachary's karate club network
g <- make_graph('Zachary')
plot(g, vertex.color='steelblue', vertex.label=NA, edge.curved=0.2)
layout_with_fr(g) # Fruchterman-Reingold layout => it tries to place the nodes in a nice way (see on internet how - opt problem)
# layout_with_kk(g) # Kamada-Kawai layout
# layout_nicely(g) # it is lighter than the previous two (which tends to show the hubs nicely)
plot(g, vertex.color='steelblue', vertex.label=NA, edge.curved=0.2, layout=layout_with_kk(g))
# NOTE: a graph plot is not so much informative, because we can place the node in a infinite number of ways
# it is not as informative as a cartesian plot


# Erdos-Renyi random graph
# remember that p=1/N is the critical point for the emergence of a giant component
N <- 100
g <- erdos.renyi.game(n=N, p=1/N, directed=FALSE) # p or m => canonical or microcanonical ensemble
lay <- layout_with_fr(g)
plot(g, layout=lay, vertex.color='steelblue', vertex.label=NA, edge.curved=0.2, vertex.size=2.5, vertex.frame.color="white")
# if p=log(N)/N => we are in the supercritical regime, where we have a giant component
g <- erdos.renyi.game(n=N, p=log(N)/N, directed=FALSE)
plot(g, layout=lay, vertex.color='steelblue', vertex.label=NA, edge.curved=0.2, vertex.size=2.5, vertex.frame.color="white")



# Geometric random graph
g <- sample_grg(n=100, radius=0.2)
lay <- layout_with_fr(g)
plot(g, layout=lay, vertex.color='steelblue', vertex.label=NA, edge.curved=0.2, vertex.size=2.5, vertex.frame.color="white")


# Barabasi-Albert random graph
N <- 500
g <- barabasi.game(N, m=1, directed=F)
lay <- layout_with_fr(g)
plot(g, layout=lay, vertex.color='steelblue', vertex.label=NA, edge.curved=0.2, vertex.size=2.5, vertex.frame.color="white")

degree(g)
plot(g, layout=lay, vertex.color='steelblue', vertex.label=NA, edge.curved=0.2, vertex.size=1.5*log(degree(g)), vertex.frame.color="white")



# Stochastic block model
N <- 1000
OM <- matrix(0.0003, ncol=4, nrow=4) # block matrix
diag(OM) <- 0.03 # higher probability of connection within the same block
g <- sample_sbm(n=N, pref.matrix=OM, block.sizes=c(200,200,300,300))
lay <- layout_with_fr(g)
plot(g, layout=lay, vertex.color='steelblue', vertex.label=NA, edge.curved=0.2, vertex.size=2.5, vertex.frame.color="white")

'''
QM <- matrix(runif(4*4, 0.3/N, 0.5/N), ncol=4, nrow=4)
QM <- (QM + t(QM))/2 # make it symmetric
diag(QM) <- diag(QM)*10 # higher probability of connection within the same block
g <- sample_sbm(n=N, pref.matrix=QM, block.sizes=c(200,200,300,300))
lay <- layout_with_fr(g)
plot(g, layout=lay, vertex.color='steelblue', vertex.label=NA, edge.curved=0.2, vertex.size=2.5, vertex.frame.color="white")
'''