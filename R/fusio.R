#' Joint effects of genetically regulated expression on two traits
#'
#' For every gene in \code{instruments}, the instruments are aligned to both GWAS and the effects of
#' expression on the two traits are estimated jointly with the single-exposure, multiple-outcome
#' model of FusioMR. Estimating both traits in one model, rather than one at a time, keeps the two
#' effects on a common scale and lets the sampler share information about pleiotropy between them.
#'
#' @param gwas1,gwas2 tables from \code{\link{read_gwas}}.
#' @param instruments output of \code{\link{build_instruments}} (one or several tissues).
#' @param control a \code{FusioMR::parameter_control} list; the defaults follow the published
#'   analysis (hybrid prior off, 20,000 iterations).
#' @param min_iv minimum instruments per gene after alignment to both GWAS.
#' @param verbose print each gene as it is fitted.
#' @return a data.table with one row per gene and tissue: b1, se1, p1, b2, se2, p2, n_iv.
#' @export
fusio_pair <- function(gwas1, gwas2, instruments, control = NULL, min_iv = 5, verbose = FALSE) {
  if (!requireNamespace("FusioMR", quietly = TRUE))
    stop("fusio_pair needs FusioMR: remotes::install_github(\"kangbw702/FusioMR\")")
  if (is.null(control)) control <- FusioMR::parameter_control(niter = 20000, c_gamma = 4, c_theta = 6,
                                                              kappa_gamma = 1, kappa_theta = 2)
  iv <- data.table::as.data.table(instruments)
  g1 <- data.table::as.data.table(gwas1); g2 <- data.table::as.data.table(gwas2)
  data.table::setkey(g1, SNP); data.table::setkey(g2, SNP)
  ## align an outcome to the instrument's effect allele; drop variants that do not match either way
  align <- function(d, g) {
    m <- g[d$SNP, nomatch = NA]
    s <- data.table::fifelse(m$EA == d$EA & m$OA == d$OA, 1,
         data.table::fifelse(m$EA == d$OA & m$OA == d$EA, -1, NA_real_))
    list(b = m$BETA * s, se = m$SE)
  }
  res <- iv[, {
    a1 <- align(.SD, g1); a2 <- align(.SD, g2)
    ok <- is.finite(a1$b) & is.finite(a2$b) & is.finite(a1$se) & is.finite(a2$se) &
          a1$se > 0 & a2$se > 0 & is.finite(BETA) & SE > 0
    if (sum(ok) < min_iv) {
      data.table::data.table(b1 = NA_real_, se1 = NA_real_, p1 = NA_real_,
                             b2 = NA_real_, se2 = NA_real_, p2 = NA_real_, n_iv = sum(ok))
    } else {
      if (verbose) cat("gene", .BY[[1]], "tissue", .BY[[2]], "n_iv", sum(ok), "\n")
      fit <- try(FusioMR::fusiomr(BETA[ok], SE[ok],
                                  cbind(a1$b[ok], a2$b[ok]), cbind(a1$se[ok], a2$se[ok]),
                                  model = "semo", control = control), silent = TRUE)
      if (inherits(fit, "try-error")) {
        data.table::data.table(b1 = NA_real_, se1 = NA_real_, p1 = NA_real_,
                               b2 = NA_real_, se2 = NA_real_, p2 = NA_real_, n_iv = sum(ok))
      } else {
        e <- unlist(fit$est); s <- unlist(fit$se); p <- unlist(fit$pval)
        data.table::data.table(b1 = e[1], se1 = s[1], p1 = p[1],
                               b2 = e[2], se2 = s[2], p2 = p[2], n_iv = sum(ok))
      }
    }
  }, by = .(gene, tissue)]
  res[is.finite(b1)][]
}
