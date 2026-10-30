#' Aggregated Cauchy association test
#'
#' Combines P values that may be arbitrarily correlated, which is what makes it suitable for
#' combining the tissues in which the same gene is tested.
#'
#' @param p a vector of P values.
#' @param weights optional non-negative weights.
#' @return a single P value.
#' @export
acat <- function(p, weights = NULL) {
  p <- p[is.finite(p)]
  if (!length(p)) return(NA_real_)
  p <- pmin(pmax(p, 1e-300), 1 - 1e-16)
  w <- if (is.null(weights)) rep(1 / length(p), length(p)) else weights[is.finite(p)] / sum(weights[is.finite(p)])
  t <- sum(w * tan((0.5 - p) * pi))
  0.5 - atan(t) / pi
}
