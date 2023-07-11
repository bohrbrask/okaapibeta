
#' Create colours based on values
#'
#' `contincols()` creates a set of colours based on a set of values and any colour scale, such that the distance between the colours correspond to the distance between the values.
#' The colours can be used to colour the nodes of an `igraph` network plot, based on a node attribute.
#'
#' For valid formats of cols see the col argument of the `col2rgb()` function.
#' Values in `colvals` must be between 0 and 1.
#' If the colours are used to colour the nodes of a network, then the length of `colvals` should equal the number of nodes in that network.
#'
#'
#' @param cols A vector with two or more colours from which the colour scale is created by extrapolation between the colours
#' @param colvals A vector with values (for example values of a node attribute)
#'
#' @return A vector with colours that can be used directly to set the node colours of an `igraph` network plot.
#' @export
#'
#' @author
#' Josefine Bohr Brask, \email{bohrbrask@@gmail.com}
#'
#' @seealso
#' \code{\link{igraph}}
#'
#'
#' @examples
#' # An example of colouring an igraph network based on node values:
#'
#'
#' # set colours for colour scale:
#' cols <- c('#CCFFFF','#003366')
#'
#' # Make simulated node attribute values:
#' colvals <- c(rep(0.1,30), rep(0.4,30), rep(0.9,40))
#'
#' # create colours for the network nodes, based on their attribute values:
#' nodecols <- contincols(cols, colvals)
#'
#'# make a network and colour it:
#' library(igraph)
#' net <- sample_gnp(100,0.1)
#' V(net)$color <- nodecols
#' plot.igraph(net)






contincols <- function(cols, colvals){
  if(any(is.na(cols))) {
    stop("cols must not contain NA's")
  }
  if(any(colvals>1) || any(colvals<0)) {
    stop("colvals must be between 0 and 1")
  }
  colfunction = colorRamp(cols)
  nodecolsraw = colfunction(colvals)
  nodecols = rgb(nodecolsraw[,1], nodecolsraw[,2], nodecolsraw[,3], maxColorValue = 255)
  return(nodecols)
}
