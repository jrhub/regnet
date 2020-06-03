#' plot a regnet object
#'
#' plot the network structures of the identified genetic variants.
#'
#' @param x regnet object.
#' @param subnetworks whether to plot sub-networks
#' @param vsize the size of the vertex
#' @param labelDist the distance of the label from the center of the vertex.
#' @param minVertices the minimum number of vertices a sub-network should contain.
#' @param theta the multiplier fro the width of the edge. Specifically, \eqn{edge.width=\theta\times adjacency}. The defualt is 1.
#' @param ... other plot arguments
#'
#' @details This function depends on the "igraph" package in generating the network graphs. It returns a (list of) igraph object(s),
#' on which users can do further modification on the network graphs.
#'
#' @return an object of class "igraph" is returned in default.
#' When \emph{subnetworks=TRUE}, a list of "igraph" objects (sub-networks) is returned.
#'
#' @usage \method{plot}{regnet}(x, subnetworks=FALSE, vsize=10, labelDist=2, minVertices=2, theta=1, \dots)
#' @seealso \code{\link{regnet}}
#'
#' @examples
#' data(ContExample)
#' X = rgn.tcga$X
#' Y = rgn.tcga$Y
#' clv = (1:2)
#' fit = regnet(X, Y, "continuous", "network", rgn.tcga$lamb1, rgn.tcga$lamb2, clv =clv, alpha.i=0.5)
#'
#' plot(fit)
#' plot(fit, subnetworks = TRUE, vsize=20, labelDist = 3, theta = 5)
#'
#'
#'@export
plot.regnet=function(x, subnetworks=FALSE, vsize=10, labelDist=2, minVertices=2, theta=1, ...){

  adjacency = x$Adj
  if(nrow(adjacency)==0) return(NULL)
  net0 <- igraph::graph.adjacency(adjacency, mode="undirected", weighted=TRUE, diag=FALSE)
  igraph::V(net0)$color = "skyblue"

  if(!subnetworks){
    plot(net0, vertex.size=vsize, edge.color="gray40", vertex.label=NA)
    net0
  }else{
    nets = igraph::decompose.graph(net0, mode="weak", min.vertices = minVertices)
    largest = 0
    if(length(nets)==0) return(NULL)
    for(i in 1: length(nets)){
        plot(nets[[i]], vertex.size=vsize, vertex.label.dist=labelDist, edge.width	= theta*(igraph::E(nets[[i]])$weight),
             edge.color="gray75", main = paste("sub-network ",i))
    }
    nets
  }
}
