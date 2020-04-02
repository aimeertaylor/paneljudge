compute_cardinalities <- function(frequencies){
  cardinalities <- apply(frequencies, 1, function(x){sum(x > 0)})
  return(cardinalities)
}

###########################################################################
#' Compute marker effective cardinalities
#'
#' Given a matrix frequencies \code{compute_eff_cardinalities()} returns the
#' effective cardinalites, \eqn{K_t^{'}}, for \eqn{t = 1,\ldots,m} markers. Each
#' calculated as followed (described in detail in Taylor et al. 2009)
#' without accounting for sample bias:
#' \deqn{K_t^{'} = \dfrac{1}{h_t}}
#'
#' @param frequencies ndata (marker count) by Kmax (max cardinality over m
#'   markers) matrix if Kt < Kmax for any t in 1:m, then fs[t,1:Kt] in (0,1) &
#'   fs[t,(Kt+1):Kmax] = 0 Example: if Kt = 2 < Kmax = 4 then fs[t,] might look
#'   like [0.2, 0.7, 0, 0].
#'
#' @return Effective cardinalities, \eqn{K_t^{'}}, for \eqn{t = 1,\ldots,m}
#'   markers.
#' @examples
#' compute_eff_cardinalities(fs = frequencies$Colombia)
#'
#' @references Taylor, A.R., Jacob, P.E., Neafsey, D.E. and Buckee, C.O., 2019.
#'   Estimating relatedness between malaria parasites. Genetics, 212(4),
#'   pp.1337-1351.
###########################################################################
compute_eff_cardinalities <- function(frequencies){
  eff_cardinalities = 1/rowSums(frequencies^2)
  return(eff_cardinalities)
}


###########################################################################
#' Compute marker diversities
#'
#' Given a matrix frequencies \code{compute_diversities()} returns the
#' diversities, \eqn{h_t}, for \eqn{t = 1,\ldots,m} markers. Each
#' calculated as described in Taylor et al. 2009
#' without accounting for sample bias:
#' \deqn{h_t^{'} = \dfrac{1}{h_t}}
#'
#' @param frequencies ndata (marker count) by Kmax (max cardinality over m
#'   markers) matrix if Kt < Kmax for any t in 1:m, then fs[t,1:Kt] in (0,1) &
#'   fs[t,(Kt+1):Kmax] = 0 Example: if Kt = 2 < Kmax = 4 then fs[t,] might look
#'   like [0.2, 0.7, 0, 0].
#'
#' @return Diversities, \eqn{h_t}, for \eqn{t = 1,\ldots,m}
#'   markers.
#' @examples
#' compute_diversities(fs = frequencies$Colombia)
#'
#' @references Taylor, A.R., Jacob, P.E., Neafsey, D.E. and Buckee, C.O., 2019.
#'   Estimating relatedness between malaria parasites. Genetics, 212(4),
#'   pp.1337-1351.
###########################################################################
compute_diversities <- function(frequencies){
  diversities = 1-rowSums(frequencies^2)
  return(diversities)
}


