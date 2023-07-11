
#' Generate a network based on the trait preference model
#'
#'`traitnet()` makes a network matrix (and if desired a plot of the network), based on the trait preference model.
#'
#'
#' This function can be considered both as a tool to generate networks based on social preferences for traits (node attributes),
#' And as a tool to generate networks with specific structural properties (potentially disregarding the preference interpretation).
#'
#' The function can create networks going from randomly structured to strongly influenced by trait preferences.
#' One or more traits, and different types of social preferences, can influence the network structure.
#' The trait values of individuals can be simulated (based on different trait types), or given as input by the user.
#' The traits (and thereby also the preference types) can vary in their importance (their influence on the structure).
#'
#' Network structures that can be generated include random, modular, assorted (with or without modularity), centralized, and mixes of these in any ratio.
#' * Popularity preferences give centralization.
#' * Similarity preferences with categorical traits give modularity and trait assortativity.
#' * Similarity preferences with continuous traits give assortativity (without modularity).
#'
#' For details about how the trait preference model generates networks, see Brask et al. in References below.
#' Briefly, the model calculates a social attraction value for each pair of individuals, based on their trait values
#' and the social preference type and importance (weight) given for each trait. From this, a network is constructed where pairs
#' with higher social attraction have a higher chance of getting a link.
#'
#'
#' The `traittypes`, `allncats`, `preftypes`, and `wvals` arguments must fit together in terms of order (and have the same length),
#' because the function assumes that the entries in the same place in these vectors fit together.
#' For example, the first trait type in the `traittypes` vector will be used with the first preference type in the `preftypes` vector,
#' the first importance value in the `wvals` vector, and the first category number (which may be NA) in the `allncats` vector.
#'
#' About `traittypes`:
#' * 'cate'  gives categorical trait values. The values consist of the integers from 1 to `ncats`, which each indicates a category. Each category is drawn with the same probability.
#' * 'tnorm' gives normally distributed trait values (truncated). The values are drawn from a truncated normal distribution that lie in the interval \[0,1], with mean = 0.5 and standard deviation = 0.25.
#' * 'ranks' gives trait values indicating ranks (e.g. social rank). The values are equally distributed in the interval \[0,1]. This corresponds to ranks indicated by integers from 1 to `n`, transformed to the interval \[0,1].
#' * 'own' must be given for any trait where trait values are given as input by the user. For such traits, the class needs to be specified in the `owntraitclasses` argument.
#'
#' About `preftypes`:
#' * There are two preference types: 'similarity' and 'popularity' (as described in Brask et al., see References).
#' * 'sim' indicates that individuals prefer others that are similar to themselves in terms of the trait.
#' * 'pop' indicates that individuals with high trait values are preferred (popular).
#'
#' About `linktype`:
#' * 'unw' gives unweighted (binary) links. Existing links are indicated with the value 1, non-existing links are indicated with 0.
#' * 'stow' gives weighted links with stochasticity, where the weights are drawn from distributions with the social attraction values as means.
#' * 'detw' gives deterministic weighted links, where the weight equals the social attraction of the dyad.
#' * For weighted links, the link weights are in the interval \[0,1], and existing links have a minimum weight of 0.001.
#'
#' About `wvals`:
#' * The sum of the `wval` vector, and each value in the vector, must lie in the interval \[0,1], where 0 indicates no importance and 1 indicates full importance.
#'
#' About `visnet`:
#' * The node size is based on the popularity trait (if any) with highest w value (importance). The node colour is based on the similarity trait (if any) with highest w value (importance).
#'
#' About `owntraitvals` and `owntraitclasses`:
#' * These arguments should only be given if the user wants to provide trait values for one or more of the traits (for example real-world trait values).
#'
#' About `owntraitvals`:
#' * If the user provides trait values for more than one trait, then `owntraitvals` must have one column for each of these traits, and `n` rows.
#' If the user provides trait values for a single trait then `owntraitvals` can be a vector of length `n` or a matrix with `n` rows and 1 column.
#' Categorical trait values should consist of integers that each indicates a category.
#' Continuous trait values should consist of numbers between 0 and 1, and empirical trait values on other scales should be transformed to this scale before being used.
#' The trait values must fit with the trait classes given in the `owntraitclasses` argument.
#'
#' About `owntraitclasses`:
#' * Give 'cate' for categorical traits.
#' * Give 'cont' for continuous traits.


#' @param n An integer indicating the number of individuals (network nodes)
#' @param k An integer indicating the average degree (average number of links connected to each node)
#' @param traittypes A vector of strings indicating the trait type for each trait ('cate', 'tnorm', 'ranks', or 'own'; see Details)
#' @param allncats A vector of integers indicating the number of categories for each trait (numbers for non-categorical traits are not used and may be set to NA, but must be included in the vector)
#' @param preftypes A vector of strings indicating the preference type to be used with each trait ('sim', or 'pop'; see Details)
#' @param linktype A string indicating the type of link ('stow', 'detw', or 'unw'; see Details)
#' @param wvals A vector with the importance (weight) of each trait-preference combination
#' @param onecomp Logical. If TRUE then the generated network will have all the nodes in a single component (i.e. at least indirectly connected). If FALSE then the generated network may consist of more than one component
#' @param visnet Logical. If TRUE then a plot of the network is made. If false then no plot is made (see also Value below)
#' @param owntraitclasses A vector of strings indicating the trait class for each trait for which trait values are provided by the user ('cate' or 'cont'; see Details)
#' @param owntraitvals A matrix with a column with trait values for traits where trait values are provided by the user (can be a vector if there is only a single trait)
#'
#' @return A list containing the following:
#' * `net`: the network matrix
#' * `traitnetplot`: a plot of the network (only if `visnet` = TRUE)
#' * `traitvals`: a matrix with all the trait values (both generated and user-supplied)
#' * `sociatvals`: a vector with a social attraction value for each dyad


#' @export

#' @references
#' \url{https://arxiv.org/abs/2303.08107}
#'
#' @author
#' Josefine Bohr Brask, \email{bohrbrask@@gmail.com}
#'
#' @seealso
#' \code{\link{traitnetsmetrics()}}
#'
#'
#' @examples
#'
#' # A network based on two traits with popularity and similarity preferences,
#' # respectively (this corresponds to the minimal version of the trait preference
#' # model, see References):
#'
#' traitnet(n=100, k=10, traittypes = c('cate', 'tnorm'), allncats = c(5, NA),
#' preftypes = c('sim', 'pop'), linktype = 'stow', wvals = c(0.5, 0.5),
#' onecomp = TRUE, visnet = TRUE)
#'
#'
#' # The same as the above example, but where the categorical similarity trait is of
#' # high importance and the other trait (with popularity) of no importance:
#'
#' traitnet(n=100, k=10, traittypes = c('cate', 'tnorm'), allncats = c(5, NA),
#' preftypes = c('sim', 'pop'), linktype = 'stow', wvals = c(0.9, 0),
#' onecomp = TRUE, visnet = TRUE)
#'
#'
#' # An example with three traits, where we specify the arguments before the call,
#' and we store the network matrix afterwards:
#'
#' n <- 100
#' k <- 10
#' linktype <- 'stow'
#' traittypes <- c('cate', 'tnorm', 'ranks')
#' allncats <- c(5, NA, NA)
#' preftypes <- c('sim', 'pop', 'sim')
#' wvals <- c(0.5, 0.1, 0.3)
#'
#' traitnetoutput <- traitnet(n=100, k=10, traittypes = traittypes, allncats = allncats,
#' preftypes = preftypes, linktype = linktype, wvals = wvals, onecomp = TRUE,
#' visnet = TRUE)
#'
#' thenet <- traitnetoutput$net
#'
#'
#' # An example with two traits, where the trait values are input by the user.
#' # (for the example, the trait values are simulated (before the call),
#' # but empirical trait data can be input in the same way):
#'
#' n <- 100
#' k <- 10
#' linktype <- 'stow'
#' traittypes <- c('own', 'own')
#' allncats <- c(NA, 2)
#' preftypes <- c('pop', 'sim')
#' wvals <- c(0.1, 0.85)
#' owntraitvals <- cbind(sample(100)/100, sample(1:2,100, replace = T))
#' owntraitclasses <- c('cont', 'cate')
#'
#' traitnet(n=100, k=10, traittypes = traittypes, allncats = allncats,
#' preftypes = preftypes, linktype = linktype, wvals = wvals, onecomp = TRUE,
#' visnet = TRUE, owntraitclasses = owntraitclasses, owntraitvals = owntraitvals)








traitnet <- function(n, k, traittypes, allncats, preftypes, linktype, wvals, onecomp, visnet, owntraitclasses = NULL, owntraitvals = NULL){

  if (n <= 0 || is.na(n)) {
    stop('n must be positive.')
  }
  if (k <= 0 || is.na(k)) {
    stop('k must be positive.')
  }
  if(sum(wvals) > 1) {
    stop("The sum of wvals must not exceed 1.")
  }
  if(any(wvals <0)||any(wvals >1)) {
    stop("all entries in wvals must be between 0 and 1.")
  }
  inputlengths <- c(length(traittypes), length(allncats), length(preftypes), length(wvals))
  if(any(inputlengths != inputlengths[1])) {
    stop("lengths of input arguments do not fit.")
  }
  if (!(linktype %in% c('stow','detw','unw')) ) {
    stop('input linktype is incorrect.')
  }
  if (!all(traittypes %in% c('cate','tnorm','ranks', 'own')) ) {
    stop('input traittypes is incorrect.')
  }
  if(any(is.na(owntraitvals))) {
    stop("owntraitvals must not contain NA's.")
  }
  if (!is.null(owntraitclasses)) {
    if (!all(owntraitclasses %in% c('cate','cont')) ) {
      stop('owntraitclasses are incorrect.')
    }
    inputlengthsown <- c(length(owntraitclasses), ncol(as.matrix(owntraitvals)), length(which(traittypes == 'own')) )
    if(any(inputlengthsown != inputlengthsown[1])) {
      stop("lengths of input arguments do not fit.")
    }
  }
  if (!all(preftypes %in% c('sim','pop')) ) {
    stop('input preftypes are incorrect.')
  }
  if (any(wvals[traittypes=='cate' & preftypes == 'sim'] >0.99) & onecomp == TRUE){
    warning(prompt = 'using sim with wval > 0.99 is not recommended as it may be slow or impossible to generate the network.', immediate. = TRUE)
  }
  if (!is.null(owntraitvals)) {
    owntmat <- as.matrix(owntraitvals)
    owntcate <- owntmat[,owntraitclasses == 'cate']
    owntcont <- owntmat[,owntraitclasses == 'cont']
    owntcatediv <- owntcate/round(owntcate)
    if (any(owntcatediv != 1)) {
      stop('categorical trait values must be integers.')
    }
    if (any(owntcont<0)||any(owntcont>1)) {
      stop('continuous trait values must be in the range from 0 to 1.')
    }
  }

  traitvals <-  matrix(NA, nrow = n, ncol = length(traittypes))
  if (any(traittypes != 'own')){
    nonownplaces <- which(traittypes != 'own')
    for (curtraitnum in seq_along(nonownplaces)){
      curtraitplace <- nonownplaces[curtraitnum]
      traitvals[, curtraitplace] <- traitvalues(n=n, traittype = traittypes[curtraitplace], ncats = allncats[curtraitplace]) # make the traitvals for this trait and add them to the traitvals matrix
    }
  }
  if (any(traittypes == 'own')){
    ownplaces <- which(traittypes == 'own')
    traitvals[,ownplaces] <- owntraitvals
  }

  traitclassesraw <- traittypes
  traitclassesraw[traitclassesraw %in% c('tnorm', 'ranks')] <- 'cont'  # set the relevant types from the types vector to cont
  if (any(traittypes == 'own')){ # add the classes from the classvector for traits that are inputted from the user (if any)
    ownplaces <- which(traittypes == 'own')
    traitclassesraw[ownplaces] <- owntraitclasses
  }
  traitclasses <- traitclassesraw

  sociatvals <- traitnetsociat(traitclasses = traitclasses, preftypes = preftypes, wvals = wvals, traitvals = traitvals)

  if (onecomp == FALSE){
    net <- traitnetbuild(n = n, k = k, linktype = linktype, sociatvals = sociatvals)
  } else if (onecomp == TRUE){ # if set to only one component, recreate networks until getting one with a single component
    componentnum <- 0
    while (round(componentnum) != 1) {  #to make sure the network is in one component
      net <- traitnetbuild(n = n, k = k, sociatvals = sociatvals, linktype = linktype)
      componentnum <- sna::components(net) # number of components of the current network. Note: igraph also has a function with this name.
    }
  }

  if (visnet == TRUE){
    traitnetplot <- traitnetvisual(net = net, traitclasses = traitclasses, preftypes = preftypes, linktype = linktype, wvals = wvals, traitvals = traitvals)
  }

  if (visnet == TRUE){
    return(list(net=net, traitnetplot=traitnetplot, traitvals = traitvals, sociatvals = sociatvals))
  } else if (visnet == FALSE){
    return(list(net = net, traitvals = traitvals, sociatvals = sociatvals))
  }

}



