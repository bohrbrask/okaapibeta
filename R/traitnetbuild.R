
#' Build a network matrix based on social attraction values
#'
#' `traitnetbuild()` generates a network (an adjacency matrix) based on social attraction values,
#' such that pairs of individuals with higher social attraction are more likely to get a link.
#'
#' The number of links is calculated based on the set average degree and number of individuals.
#' Then the dyads to be linked are found by a weighted draw, where dyads with higher social attraction have a higher chance of being picked.
#' The link weights are determined based on the set link type and added to the existing links.
#'
#' To construct networks based on the trait preference model, it is recommended to use the `traitnet()` and `traitnetsmetrics()` functions.
#' Those functions call this one, and there is thus no need to use this one directly.
#'
#' The number of individuals `n` must fit with the length of the `sociatvals` vector (which equals the number of dyads).
#'
#' If there are not enough dyads with social attraction values above zero to give the input average degree,
#' then random links will be added among the remaining dyads to get the input average degree.
#'
#' About `linktype`:
#' * 'unw' gives unweighted (binary) links. Existing links are indicated with the value 1, non-existing links are indicated with 0.
#' * 'stow' gives weighted links with stochasticity, where the weights are drawn from normal distributions with the social attraction values as means (sd = 0.05).
#' * 'detw' gives deterministic weighted links, where the weight equals the social attraction of the dyad.
#' * For weighted links, the link weights are in the interval \[0,1], and existing links have a minimum weight of 0.001.
#'
#'
#' @param n An integer indicating the number of individuals (network nodes)
#' @param k An integer indicating the average degree (average number of links connected to each node)
#' @param linktype A string indicating the type of link ('stow', 'detw', or 'unw'; see Details)
#' @param sociatvals A vector with a social attraction value (a number between 0 and 1) for each pair of individuals (each dyad of nodes)
#'
#' @return An undirected network matrix (adjacency matrix). The matrix is symmetric, with zeros in the diagonal.
#' @export
#'
#'@references
#' \url{https://arxiv.org/abs/2303.08107}
#'
#' @author
#' Josefine Bohr Brask, \email{bohrbrask@@gmail.com}
#'
#' @seealso
#' \code{\link{traitnet()}}
#'
#'
#' @examples
#'
#' n <- 100
#' k <- 10
#' linktype <- 'stow'
#'
#' # For the example we use random values as social attraction values
#' # (see the examples in the `traitnetsociat()`  help page to generate
#' social attraction values based on the trait preference model):
#' ndyads <- (n*n-n)/2
#' sociatvals <-  runif(ndyads)
#'
#' # Make the network matrix:
#' net <- traitnetbuild(n, k, linktype, sociatvals)



traitnetbuild <- function(n, k, linktype, sociatvals){

  if (n <= 0 || is.na(n)) {
    stop('n must be positive.')
  }
  if (k <= 0 || is.na(k)) {
    stop('k must be positive.')
  }
  if (!(linktype %in% c('stow','detw','unw')) ) {
    stop('input linktype is incorrect')
  }
  if (length(sociatvals) != ((n*n-n)/2)) {
    stop('length of sociatvals does not fit with n')
  }

  ndyads <- (n*n-n)/2
  inames <- rep(1:(n-1),(n-1):1)
  jnames <- sequence((n-1):1, from = 2:n)
  L <- (k*n)/2
  unwnet <- matrix(0,n,n)
  linkdyadnums = wrswoR::sample_int_crank(ndyads, L, sociatvals)
  linkinames = inames[linkdyadnums]
  linkjnames = jnames[linkdyadnums]
  linkinamesboth <- c(linkinames, linkjnames)
  linkjnamesboth <- c(linkjnames, linkinames)
  unwnet[cbind(linkinamesboth, linkjnamesboth)] <- 1

  if (linktype =='detw'){
    net <- unwnet
    linkweightsmod <- sociatvals[linkdyadnums]
    linkweightsmod[linkweightsmod==0]<- 0.001
    linkweights <- linkweightsmod
    net[cbind(linkinamesboth, linkjnamesboth)] <- linkweights
  } else if (linktype =='stow'){
    net <- unwnet
    linkwmin <- 0.001
    linkwvar <- 0.05
    rawlinkw <- sociatvals[linkdyadnums]
    rawlinkdraw <- rnorm(length(rawlinkw), mean = rawlinkw, sd = linkwvar)
    linkweightsmod <- rawlinkdraw
    linkweightsmod[linkweightsmod <= 0] <- linkwmin
    linkweightsmod[linkweightsmod > 1] <- 1
    linkweights <- linkweightsmod
    net[cbind(linkinamesboth, linkjnamesboth)] <- linkweights
  }
  else if (linktype =='unw'){
    net <- unwnet
  }
  return(net)
}

