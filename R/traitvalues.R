

#' Generate trait values
#'
#' `traitvalues()` generates a vector of trait values based on a given trait type.
#'
#'
#' To generate trait values, the function draws values from a distribution that corresponds to the given trait type.
#'
#' To make networks based on trait preferences, it is recommended to use the `traitnet()` or `traitnetsmetrics()` functions.
#' Those functions call this function and there is no need to use this function directly.
#'
#' About `traittype`:
#' * 'cate'  gives categorical trait values. The values consist of the integers from 1 to `ncats`, which each indicates a category. Each category is drawn with the same probability.
#' * 'tnorm' gives normally distributed trait values (truncated). The values are drawn from a truncated normal distribution that lie in the interval \[0,1], with mean = 0.5 and standard deviation = 0.25.
#' * 'ranks' gives trait values indicating ranks (e.g. social rank). The values are homogeneously distributed in the interval \[0,1]. This corresponds to ranks indicated by integers from 1 to `n`, transformed to the interval \[0,1].


#' @param n An integer indicating the number of trait values to be made (i.e. number of individuals for which the trait values are made)
#' @param traittype A string indicating the trait type ('cate', 'tnorm', or 'ranks'; see Details)
#' @param ncats An integer indicating the number of categories (only used if `traittype` is set to categorical)
#'
#' @return A vector of trait values. The vector has length `n`. The values are in random order. Non-categorical values are in the interval \[0,1]; categorical values are given by integers.

#' @export
#'
#' @references
#' \url{https://arxiv.org/abs/2303.08107}
#'
#' @author
#' Josefine Bohr Brask, \email{bohrbrask@@gmail.com}
#'
#' @seealso \code{\link{traitnet()}}
#' @seealso \code{\link{traitnetsmetrics()}}
#'
#' @examples
#' # Categorical trait values with two categories, for 100 individuals:
#' traitvalues(n = 100, traittype = 'cate', ncats = 2)
#'
#' @examples
#' # Continuous trait values (from a truncated normal distribution), for 100 individuals:
#' traitvalues(n = 100, traittype = 'tnorm')
#'
#' @examples
#' # Rank trait values, for 50 individuals:
#' traitvalues(n = 50, traittype = 'ranks')



traitvalues <- function(n, traittype, ncats) {
  if (n <= 0 || is.na(n)) {
    stop('n must be positive.')
  }
  if (traittype == 'cate') {
    if (ncats <= 0 || length(ncats)==0 || is.na(ncats)) {
      stop('ncats must be a positive integer.')
    }
    traitvals <- sample(ncats,n, replace = TRUE)
  } else if (traittype == 'tnorm') {
    traitvals <- truncnorm::rtruncnorm(n, a=0,b=1,mean=0.5, 0.25)
  }else if (traittype == 'ranks') {
    traitvals <- sample(n,n, replace = FALSE)/n
  }
  return(traitvals)
}

