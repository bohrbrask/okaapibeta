

#' Measure network metrics on a single network
#'
#'`netmetrics()` measures a set of network metrics on an input network.
#'
#'
#'The function uses the `igraph` package for some of the metric calculations.
#'
#'
#' About `net`:
#' * The input network must be undirected, and the network matrix must be symmetric with zeros in the diagonal.
#'
#' About `metrics`:
#'
#' * 'clust' gives the global clustering coefficient, also called transitivity (measured as the ratio of triangles to connected triples).
#' * 'clustw' gives the weighted clustering coefficient (as defined in Barrat et al. 2004, see References). Nodes with one link get a weighted clustering of zero.
#'
#' * 'path' gives the average path length, i.e. the average unweighted distance between nodes.
#' * 'pathw' gives the average weighted path length, i.e. the average weighted distance between nodes (with smaller distances for larger weights).
#'   - For networks with more than one component: nodes in different network components are considered to have no path between them and these dyads therefore do not contribute to path length calculations.
#'
#' * 'avdeg' gives the average degree.
#' * 'avdegw' gives the average weighted degree (also called strength).
#'
#' * 'degassort' gives the degree assortativity.
#' * 'degassortw' gives the weighted degree assortativity.
#'
#' * 'degvar' gives the degree variation (variance).
#' * 'degvarw' gives the weighted degree variation (variance).
#'
#'
#' @param net A network (adjacency matrix)
#' @param metrics A vector of strings indicating the metrics that should be measured (see Details)
#'
#' @return A vector with a value for each of the metrics specified in the `metrics` argument. The columns are named with the metric names given in the `metrics` argument.

#' @export

#' @references
#' Alain Barrat, Marc Barthelemy, Romualdo Pastor-Satorras, Alessandro Vespignani: The architecture of complex weighted networks, Proc. Natl. Acad. Sci. USA 101, 3747 (2004)
#'
#' @author
#' Josefine Bohr Brask, \email{bohrbrask@@gmail.com}
#'
#' @seealso  \code{\link{traitnet()}}
#' @seealso  \code{\link{traitnetsmetrics()}}
#' @seealso  \code{\link{igraph}}
#'
#' @examples
#'
#' # An example that uses a network created with the traitnet function:
#'
#' traitnetoutput <- traitnet(n=100, k=10, traittypes = c('cate', 'tnorm'),
#' allncats = c(5, NA), preftypes = c('sim', 'pop'), linktype = 'stow',
#' wvals = c(0.8, 0.1), onecomp = TRUE, visnet = TRUE)
#'
#' net <- traitnetoutput$net
#'
#' metrics <- c('clust', 'clustw', 'path', 'pathw', 'avdeg', 'avdegw',
#' 'degassort', 'degassortw', 'degvar', 'degvarw')
#'
#' netmetrics(net, metrics)
#'
#'
#' # An example of calculating clustering and path length on
#' two common network types (made with igraph):
#'
#'library(igraph)
#'
#' # A random network
#' randnet <-  as_adjacency_matrix(sample_gnp(100, 0.1), type = 'both',
#' names = FALSE, sparse = FALSE)
#'
#' # A scale-free (Barabasi-Albert) network
#' scalenet <-  as_adjacency_matrix(sample_pa(100, directed = FALSE),
#' type = 'both', names = FALSE, sparse = FALSE)
#'
#' # Calculate clustering and path length for the two networks
#' metrics <- c('clust', 'path')
#' randnetmetrics <- netmetrics(randnet, metrics)
#' scalenetmetrics <- netmetrics(scalenet, metrics)









netmetrics <- function(net, metrics){

  if(any(diag(net)!=0)) {
    stop("The input network must have only zeros in the diagonal.")
  }
  if(!isSymmetric(net)) {
    stop("The input network must be symmetric.")
  }
  if(any(is.na(net))) {
    stop("The network must not contain NA's.")
  }
  if (!all(metrics %in% c('clust', 'clustw', 'path', 'pathw', 'avdeg', 'avdegw', 'degassort', 'degassortw', 'degvar', 'degvarw')) ) {
    stop('input metrics are incorrect.')
  }

  metricvals <- rep(NA, length(metrics))
  names(metricvals) <- metrics
  netigraph <- igraph::graph_from_adjacency_matrix(net, mode = "undirected", diag = FALSE, weighted = TRUE) # get igraph version of network
  degs <- colSums(net>0)
  weidegs <- colSums(net)

  if ('clust' %in% metrics){
    metricvals['clust'] <- igraph::transitivity(netigraph, type= "global") # calculate and store unweighted transitivity (global clustering = ratio of triangles and connected triples)
  }
  if ('clustw' %in% metrics){
    transes <- igraph::transitivity(netigraph, type= "weighted")
    transesnonan <- transes
    transesnonan[transesnonan =='NaN'] <- 0 #we set NaN values to zero (individuals with only one link get NaN but we give them zero instead)
    metricvals['clustw'] <- mean(transesnonan)     # weighted clustering
  }
  if ('path' %in% metrics){
    metricvals['path'] <- igraph::mean_distance(netigraph, directed = FALSE, weights = NA) # calculate unweighted mean distance (average path length)
  }
  if ('pathw' %in% metrics){
    nonzeroweights <- net[net>0] # get nonzero links
    invnonzeroweights <- 1/nonzeroweights # get the inverse of the nonzero edge weights
    nonzeroplaces <- which(net>0) # find places of nonzero links
    invnet <- net
    invnet[nonzeroplaces] <- invnonzeroweights #insert the inverted link strengths in place of the originals
    invnetigraph <- igraph::graph_from_adjacency_matrix(invnet, mode = "undirected", diag = FALSE, weighted = TRUE)
    metricvals['pathw'] <- igraph::mean_distance(invnetigraph, directed = F)
  }
  if ('avdeg' %in% metrics){
    metricvals['avdeg'] <- mean(degs) # calculate mean degree (unweighted)
  }
  if ('avdegw' %in% metrics){
    metricvals['avdegw'] <- mean(weidegs) # calculate mean weighted degree
  }
  if ('degassort' %in% metrics){
    metricvals['degassort'] <- igraph::assortativity_degree(netigraph, directed = F)# calculate degree assortativity
  }
  if ('degassortw' %in% metrics){
    metricvals['degassortw'] <- igraph::assortativity(netigraph, weidegs, directed = F)# calculate weighted degree assortativity
  }
  if ('degvar' %in% metrics){
    metricvals['degvar'] <- var(degs)   # calculate degree variation
  }
  if ('degvarw' %in% metrics){
    metricvals['degvarw'] <- var(weidegs)   # calculate weighted degree variation
  }

  return(metricvals)

}


