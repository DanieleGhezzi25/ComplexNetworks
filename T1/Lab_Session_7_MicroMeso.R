library(ggplot2)
library(RColorBrewer)
library(rgl)

source("common.R")

##########################
# CENTRALITY
##########################

# set and view a palette
#mypal <- colorRampPalette(rev(c(brewer.pal(9, "YlGnBu"))))(20)
mypal <- viridis::viridis_pal()(20)

plot(1, 1, type = "n", xlim = c(1, length(mypal)), ylim = c(1, 2), xlab = "", ylab = "")
for (i in 1:length(mypal)) {
  rect(i - 0.5, 1, i + 0.5, 2, col = mypal[i])
}
legend("topright", legend = mypal, fill = mypal, bg = "white")


set.seed(12345)

N <- 1000
#g <- barabasi.game(N, m=1, directed=TRUE)
g <- watts.strogatz.game(dim=2, size=20, p=0.05, nei=1)
lay3D <- layout_with_fr(g, dim=3)
lay2D <- layout_with_fr(g, dim=2)


# degree distribution
deg <- degree(g)
#hist(degree(g))
centr_deg <- deg

# normalize and convert to colors
cols_deg <- vec2pal(centr_deg, mypal)

# color nodes by centrality
rgl.clear()
rglplot(g, layout=lay3D, vertex.color=cols_deg, vertex.size=deg, vertex.label=NA, edge.color="gray80", vertex.frame.color=NA, edge.width=0.5)


# betweenness
centr_betw <- betweenness(g, normalized=TRUE)

# normalize and convert to colors
cols_betw <- vec2pal(log(1+centr_betw), mypal)

# color nodes by centrality
rgl.clear()
rglplot(g, layout=lay3D, vertex.color=cols_betw, vertex.size=sqrt(deg), vertex.label=NA, edge.color="gray80", vertex.frame.color=NA)


# let's calculate some centralities now and we will compare them later at once

# closeness
centr_clos <- closeness(g, normalized=TRUE)
cols_clos <- vec2pal(centr_clos, mypal)

# pagerank
centr_pr <- page_rank(g)$vector
cols_pr <- vec2pal(centr_pr, mypal)

# eigenvector
centr_eig <- eigen_centrality(g)$vector
cols_eig <- vec2pal(centr_eig, mypal)

# hub score
centr_hub <- hub_score(g)$vector
cols_hub <- vec2pal(centr_hub, mypal)

# authority score
centr_auth <- authority_score(g)$vector
cols_auth <- vec2pal(centr_auth, mypal)



# build the dataframe

# Create the dataframe
df <- data.frame(
  Degree = centr_deg,
  Betweenness = centr_betw,
  Closeness = centr_clos,
  PageRank = centr_pr,
  Eigenvector = centr_eig,
  Hub = centr_hub,
  Authority = centr_auth
)

# scatter matrix to see the patterns:
pairs(df, pch=19, cex=0.5)

# now let's plot the results on the network and compare:
par(mfrow=c(3, 2), mar=c(2, 2, 1, 1))

plot(g, layout=lay2D, vertex.color=cols_deg, vertex.size=deg, vertex.label=NA, edge.color="gray80", vertex.frame.color=NA, edge.width=0.3, main="Degree")
plot(g, layout=lay2D, vertex.color=cols_betw, vertex.size=deg, vertex.label=NA, edge.color="gray80", vertex.frame.color=NA, edge.width=0.3, main="Betweenness")
plot(g, layout=lay2D, vertex.color=cols_clos, vertex.size=deg, vertex.label=NA, edge.color="gray80", vertex.frame.color=NA, edge.width=0.3, main="Closeness")
plot(g, layout=lay2D, vertex.color=cols_pr, vertex.size=deg, vertex.label=NA, edge.color="gray80", vertex.frame.color=NA, edge.width=0.3, main="PageRank")
plot(g, layout=lay2D, vertex.color=cols_eig, vertex.size=deg, vertex.label=NA, edge.color="gray80", vertex.frame.color=NA, edge.width=0.3, main="Eigenvector")
plot(g, layout=lay2D, vertex.color=cols_hub, vertex.size=deg, vertex.label=NA, edge.color="gray80", vertex.frame.color=NA, edge.width=0.3, main="Hub")




##########################
# COMMUNITY DETECTION
##########################

# reset the last par
dev.off()

# Zachary
g <- make_graph('Zachary')
plot(g)
comms_lou <- cluster_louvain(g, resolution=1)

comms_lou$membership

plot(g, vertex.color=comms_lou$membership, vertex.size=degree(g), vertex.label=NA, edge.color="gray80", vertex.frame.color="white")

# test the 3D version!
# rgl.clear()
# rglplot(g, layout=layout_with_fr(g, dim=3), 
#      vertex.label=NA, 
#      vertex.color=comms_lou$membership,
#      vertex.size=degree(g),
#      edge.color="gray70",
#      edge.width=0.25)

# use rgl.clear() to clean
# use rgl.bg(color="black") to change background color

comms_spin <- cluster_spinglass(g, gamma=1)
comms_info <- cluster_infomap(g)

modularity(comms_lou)
modularity(comms_spin)
modularity(comms_info)

lay <- layout_with_fr(g, dim=2)

par(mfrow=c(1, 3), mar=c(0, 0, 1, 1))

plot(g, layout=lay, vertex.color=comms_lou$membership, vertex.size=degree(g), vertex.label=NA, edge.color="gray80", vertex.frame.color="white",  main=paste("Louvain - ", round(modularity(comms_lou),4)) ) 
plot(g, layout=lay, vertex.color=comms_spin$membership, vertex.size=degree(g), vertex.label=NA, edge.color="gray80", vertex.frame.color="white", main=paste("Spin - ", round(modularity(comms_spin),4)) )
plot(g, layout=lay, vertex.color=comms_info$membership, vertex.size=degree(g), vertex.label=NA, edge.color="gray80", vertex.frame.color="white", main=paste("Infomap - ", round(modularity(comms_info),4)) )




# SBM with 4 groups
N <- 1000
OM <- matrix(0.001, ncol=4, nrow=4)
diag(OM) <- 0.03
g <- sample_sbm(N, pref.matrix=OM, block.sizes=c(200, 200, 300, 300))

# verify that there is only one CC
clusters(g)

lay <- layout_with_fr(g)
planted_partition <- c(rep(1,200), rep(2,200), rep(3,300), rep(4,300))

comms_lou <- cluster_louvain(g, resolution=1)
comms_spin <- cluster_spinglass(g, gamma=1)
comms_info <- cluster_infomap(g)

modularity(comms_lou)
modularity(comms_spin)
modularity(comms_info)

par(mfrow=c(1, 3), mar=c(0, 0, 1, 1))

plot(g, layout=lay, vertex.color=comms_lou$membership, vertex.size=sqrt(degree(g)), vertex.label=NA, edge.color="gray80", vertex.frame.color="white",  main=paste("Louvain - ", round(modularity(comms_lou),4)) ) 
plot(g, layout=lay, vertex.color=comms_spin$membership, vertex.size=sqrt(degree(g)), vertex.label=NA, edge.color="gray80", vertex.frame.color="white", main=paste("Spin - ", round(modularity(comms_spin),4)) )
plot(g, layout=lay, vertex.color=comms_info$membership, vertex.size=sqrt(degree(g)), vertex.label=NA, edge.color="gray80", vertex.frame.color="white", main=paste("Infomap - ", round(modularity(comms_info),4)) )



# understand who's doing better, compare the planted partition vs the result
compare(planted_partition, comms_lou$membership, method="nmi")

compare(planted_partition, comms_spin$membership, method="nmi")

compare(planted_partition, comms_info$membership, method="nmi")


# test the 3D version!
# rglplot(g, layout=layout_with_fr(g, dim=3), 
#      vertex.label=NA, 
#      vertex.color=comms_spin$membership,
#      vertex.size=degree(g)*0.5,
#      edge.color="gray70",
#      edge.width=0.25)

# export with: rgl.snapshot('3dplot.png', fmt = 'png')


# let's explore the impact of the resolution parameter

df_mod <- data.frame()
for (gam in 10^seq(-2, 1, 0.25)) {
    comms_lou <- cluster_louvain(g, resolution = gam)
    modularity <- modularity(comms_lou)
    groups <- length(unique(membership(comms_lou)))

    cat(paste(" Groups:", groups,"\n"))
    df_mod <- rbind(df_mod, data.frame(gamma = gam, modularity = modularity, groups=groups))
}

# plot it
ggplot(df_mod, aes(x = gamma, y = modularity)) + theme_bw() + geom_line()


