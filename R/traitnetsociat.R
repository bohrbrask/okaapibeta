

#' Calculate social attraction values
#'
#' `traitnetsociat()` calculates social attraction values between individuals, based on social preferences for traits (attributes).
#' The social attraction can be based on one or multiple traits, which may be used with different preference types.
#'
#'
#'
#'
#'
#' The `traitnetsociat()` function calculates the social attraction for each pair of individuals, based on the two individuals'
#' trait values (attributes), and the preference type and importance (weight) set for each trait (such that traits with higher weight have more influence on the social attraction).
#' The calculation follows the preference model (see References).
#'
#' To construct networks based on the trait preference model, it is recommended to use the `traitnet()` and `traitnetsmetrics()` functions.
#' Those functions call this one, and there is thus no need to use this one directly.
#'
#' The function assumes that each input is ordered in the same way. For example, the first trait class given in the
#' `traitclasses` vector will be matched with the first preference type given in the `preftype` vector, the first
#' importance value given in the `wvals` vector, and the first column of trait values in the `traitvals` matrix. This also
#' means that the trait values for a given trait must be of the class given for that trait.
#'
#' About `traitclasses`:
#' * 'cate' indicates a categorical trait.
#' * 'cont' indicates a continuous trait.
#'
#' About `preftypes`:
#' * There are two preference types: 'similarity' and 'popularity' (as described in Brask et al., see References).
#' * 'sim' indicates that individuals prefer others that are similar to themselves in terms of the trait.
#' * 'pop' indicates that individuals with high trait values are preferred (popular).
#'
#' About `wvals`:
#' * The sum of the `wval` vector, and each value in the vector, must lie in the interval \[0,1], where 0 indicates no importance and 1 indicates full importance.
#'
#' About `traitvals`:
#' * The trait values must fit with the trait classes given in the `traitclasses` argument.
#' Categorical trait values should consist of integers that each indicates a category.
#' Continuous trait values should consist of numbers between 0 and 1. Empirical trait values on other scales should be transformed to this scale before being used.
#'
#'
#'
#'
#' @param traitclasses A vector of strings indicating the trait class for each trait ('cate' or 'cont'; see Details)
#' @param preftypes A vector of strings indicating the preference type to be used with each trait ('sim', or 'pop'; see Details)
#' @param wvals A vector with the importance (weight) of each trait-preference combination
#' @param traitvals A matrix with a column of trait values for each trait, and a row for each individual. Can be a vector if there is only a single trait
#'
#' @return A vector of social attraction values (one value for each pair of individuals).
#' @export

#' @references
#' \url{https://arxiv.org/abs/2303.08107}
#'
#' @author
#' Josefine Bohr Brask, \email{bohrbrask@@gmail.com}
#'
#' @seealso
#' \code{\link{traitnet()}}
#'

#' @examples
#'
#' # An example of social attraction values that are based on the
#' # combined preferences for three traits of different type:
#'
#' # Create trait values based on some trait types:
#' traittypes <- c('cate', 'tnorm', 'ranks')
#' n <- 100
#' tval1 <- traitvalues(n = n, traittype = traittypes[1], ncats = 2)
#' tval2 <- traitvalues(n = n, traittype = traittypes[2])
#' tval3 <- traitvalues(n = n, traittype = traittypes[3])
#' traitvals <- cbind(tval1, tval2, tval3)
#'
#' # Set the class, preference type and importance for each trait
#' (The trait classes must fit with the trait values):
#' traitclasses <- c('cate', 'cont', 'cont')
#' preftypes <- c('sim', 'pop', 'sim')
#' wvals <- c(0.7, 0.1, 0.1)
#'
#' # Generate the social attraction values:
#' sociatvals <- traitnetsociat(traitclasses, preftypes, wvals, traitvals)
#'
#'
#'




traitnetsociat <- function(traitclasses, preftypes, wvals, traitvals){

  if(sum(wvals) > 1) {
    stop("The sum of wvals must not exceed 1.")
  }
  if(any(wvals <0)||any(wvals >1)) {
    stop("all entries in wvals must be between 0 and 1.")
  }
  if (!all(traitclasses %in% c('cate','cont')) ) {
    stop('input traitclasses are incorrect.')
  }
  inputlengths <- c(length(traitclasses), length(preftypes), length(wvals), ncol(traitvals))
  if(any(inputlengths != length(traitclasses))) {
    stop("inputs must have the same lengths.")
  }
  if (!all(preftypes %in% c('sim','pop')) ) {
    stop('input preftypes are incorrect.')
  }
  if(any(is.na(traitvals))) {
    stop("traitvals must not contain NA's.")
  }
  tmat <- as.matrix(traitvals)
  tcate <- tmat[,traitclasses == 'cate']
  tcont <- tmat[,traitclasses == 'cont']
  tcatediv <- tcate/round(tcate)
  if (any(tcatediv != 1)) {
    stop('categorical trait values must be integers.')
  }
  if (any(tcont<0)||any(tcont>1)) {
    stop('continuous trait values must be in the range from 0 to 1.')
  }

  simprefs <- function(simtype, tsim){
    n <- length(tsim)
    inames <- rep(1:(n-1),(n-1):1)
    jnames <- sequence((n-1):1, from = 2:n)
    tsimalli <- tsim[inames]
    tsimallj <- tsim[jnames]
    if (simtype == 'cate') {
      psim  <- ifelse(tsimalli-tsimallj ==0, 1,0)
    } else if (simtype == 'cont'){
      psim <- 1 - abs(tsimalli-tsimallj)
    }
  }

  popprefs <- function(poptype, tpop){
    n <- length(tpop)
    inames <- rep(1:(n-1),(n-1):1)
    jnames <- sequence((n-1):1, from = 2:n)
    tpopalli <- tpop[inames]
    tpopallj <- tpop[jnames]
    if (poptype == 'cate') {
      npopcats <- length(unique(tpop))
      popcatvals = seq(0,1,length.out = npopcats)
      ppop <- (popcatvals[tpopalli]+popcatvals[tpopallj])/2
    } else if (poptype == 'cont'){
      ppop <- (tpopalli+tpopallj)/2
    }
  }

  traitvals <- as.matrix(traitvals)
  n <- nrow(traitvals)
  ndyads <- (n*n-n)/2
  ptraits <- matrix(NA, nrow = ndyads, ncol = length(traitclasses))
  for (traitnum in seq_along(traitclasses)) {
    if (preftypes[traitnum]=='sim'){
      ptraits[,traitnum] <- simprefs(simtype = traitclasses[traitnum], tsim = traitvals[,traitnum])
    } else if (preftypes[traitnum]=='pop'){
      ptraits[,traitnum] <- popprefs(poptype = traitclasses[traitnum], tpop = traitvals[,traitnum])
    }
  }
  prand <- runif(ndyads)
  wrand <- 1-sum(wvals)
  ptws <- t(t(ptraits)*wvals)
  prandtws <- wrand*prand
  ptwsums <- rowSums(ptws)
  sociatvals <- ptwsums+prandtws

  return(sociatvals)

}

