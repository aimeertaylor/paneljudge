% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_marker_summaries.R
\name{compute_diversities}
\alias{compute_diversities}
\title{Function to compute marker diversities}
\usage{
compute_diversities(fs, warn_fs = TRUE)
}
\arguments{
\item{fs}{Matrix of marker allele frequencies, i.e. the \eqn{ft}s in [1].
Specifically, a \eqn{m} by \eqn{Kmax} matrix, where \eqn{m} is the marker
count and \eqn{Kmax} is the maximum cardinality (per-marker allele count)
observed over all \eqn{m} markers. If, for any \eqn{t = 1,...,m}, the
maximum cardinality exceeds that of the \eqn{t}-th marker (i.e. if
\eqn{Kmax > Kt}), then all \code{fs[t,1:Kt]} are in (0,1] and all
\code{fs[t,(Kt+1):Kmax]} are zero. For example, if \eqn{Kt = 2} and
\eqn{Kmax = 4} then \code{fs[t,]} might look like \code{[0.3, 0.7, 0, 0]}.}

\item{warn_fs}{Logical indicating if the function should return warnings
following allele frequency checks.}
}
\value{
Diversities for \eqn{t = 1,\ldots,m} markers.
}
\description{
Given a matrix of marker allele frequencies, \code{compute_diversities}
returns the diversities of \eqn{t = 1,...,m} markers, where \eqn{m} is the
marker count. Each diversity is calculated as described in [1], i.e. without
correcting for finite sample sizes or considering uncertainty.
}
\examples{
compute_diversities(fs = frequencies$Colombia)

}
\references{
\enumerate{ \item Taylor, A.R., Jacob, P.E., Neafsey, D.E. and
  Buckee, C.O., 2019. Estimating relatedness between malaria parasites.
  Genetics, 212(4), pp.1337-1351.}
}
