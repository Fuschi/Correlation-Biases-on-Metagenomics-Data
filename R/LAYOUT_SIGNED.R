#' Signed-Weighted Graph Layout
#'
#' @description It elaborates the coordinates for the representation of
#' the vertices of the graph considering only the links with a positive sign.
LAYOUT_SIGNED <- function(g){
  g.sub <- igraph::subgraph.edges(graph=g,
                              eids=which(igraph::E(g)$weight>0),
                              delete.vertices=FALSE)
  
  layout <- igraph::layout.fruchterman.reingold(g.sub)
  return(layout)
}
