#' Combine the tissues of each gene
#'
#' A gene is usually testable in several tissues, and which tissue carries the regulatory effect is
#' rarely known in advance. The tissue-level P values are therefore combined with the aggregated
#' Cauchy test, which is valid however strongly the tissues are correlated, and the direction is
#' taken from the tissue with the smallest P value. Because ACAT returns one P value per gene, the
#' number of tissues does not multiply the number of tests.
#'
#' @param fit output of \code{\link{fusio_pair}}.
#' @param z_max drop gene-tissue pairs whose |z| exceeds this on either trait.
#' @return one row per gene: effects and P values for both traits, and the tissue they came from.
#' @export
aggregate_genes <- function(fit, z_max = 9) {
  d <- data.table::as.data.table(fit)
  n0 <- nrow(d)
  d <- d[is.finite(b1) & is.finite(b2) & abs(b1 / se1) < z_max & abs(b2 / se2) < z_max]
  if (!nrow(d)) {
    warning("every gene was removed by z_max = ", z_max,
            "; raise it if the estimates are legitimately this precise (", n0, " gene-tissue pairs in)")
    return(d[0])
  }
  d[, `:=`(pmin_ = pmin(p1, p2))]
  d[order(pmin_), .(
      b1 = b1[1], se1 = se1[1], p1_acat = acat(p1),
      b2 = b2[1], se2 = se2[1], p2_acat = acat(p2),
      tissue = tissue[1], n_tissue = .N, n_iv = n_iv[1]), by = gene]
}

#' Expression-effect genetic correlation between the two traits
#'
#' rho_GE is the correlation, across genes, of the standardised effect of genetically regulated
#' expression on one trait with that on the other. The correlation is taken on z scores
#' (effect divided by its standard error) rather than on the effects themselves: genes differ
#' greatly in how precisely their expression effect is estimated, and correlating raw effects lets
#' the imprecise ones dominate. Genes in the extended MHC are dropped because its long-range linkage
#' disequilibrium breaks the independence of gene-level estimates, and \code{z_max} removes the
#' outliers that single poorly instrumented genes can produce.
#'
#' @param genes output of \code{\link{aggregate_genes}}.
#' @param block optional block label per gene (for example a 5-Mb bin); the confidence interval is
#'   then a jackknife over blocks rather than the parametric interval, because neighbouring genes
#'   share instruments and are not independent.
#' @param chr,pos optional gene coordinates, used only to drop the extended MHC.
#' @param z_max drop genes whose |z| exceeds this on either trait.
#' @param method "pearson" or "spearman".
#' @return a list with the correlation, its standard error, confidence interval, P value and n.
#' @export
rho_ge <- function(genes, block = NULL, chr = NULL, pos = NULL, z_max = 9,
                   method = c("pearson", "spearman")) {
  method <- match.arg(method)
  d <- data.table::as.data.table(genes)
  d[, `:=`(z1 = b1 / se1, z2 = b2 / se2)]
  keep <- is.finite(d$z1) & is.finite(d$z2) & abs(d$z1) < z_max & abs(d$z2) < z_max
  if (!is.null(chr) && !is.null(pos))            # extended MHC, GRCh37/38 chr6:25-34 Mb
    keep <- keep & !(as.character(chr) %in% c("6", "chr6") & pos >= 25e6 & pos <= 34e6)
  d <- d[keep]; blk <- if (is.null(block)) NULL else block[keep]
  if (nrow(d) < 3) stop("rho_GE needs at least three genes passing the filters; got ", nrow(d))
  r <- stats::cor(d$z1, d$z2, method = method)
  ct <- stats::cor.test(d$z1, d$z2, method = method)
  out <- list(rho = r, n_genes = nrow(d), P = ct$p.value, method = method,
              se = NA_real_, ci = as.numeric(ct$conf.int))
  if (!is.null(blk)) {
    u <- unique(blk); k <- length(u)
    loo <- vapply(u, function(x) stats::cor(d$z1[blk != x], d$z2[blk != x], method = method), numeric(1))
    out$se <- sqrt((k - 1) / k * sum((loo - mean(loo))^2))
    out$ci <- r + c(-1.96, 1.96) * out$se
    out$n_blocks <- k
  }
  out
}

#' Directional classes of pleiotropy
#'
#' Classifies each gene by the signs of its two effects, after selecting the genes that are
#' associated with both traits. Same-signed effects are antagonistic in the sense used for disease
#' and fitness, where a liability-raising change is paired with a fitness benefit; opposite-signed
#' effects are synergistic. Rename them if that reading does not fit your traits.
#'
#' @param genes output of \code{\link{aggregate_genes}}.
#' @param fdr Benjamini-Hochberg threshold applied to each trait separately.
#' @return a list with \code{genes} (the selected genes, each labelled with its class and with
#'   \code{type}, antagonistic or synergistic), \code{shares} (the count and share of each class),
#'   the two overall shares, the number of classified genes, and \code{P}, a two-sided binomial test
#'   of whether antagonistic and synergistic genes are equally common.
#' @export
direction_classes <- function(genes, fdr = 0.05) {
  d <- data.table::as.data.table(genes)
  d[, `:=`(q1 = stats::p.adjust(p1_acat, "BH"), q2 = stats::p.adjust(p2_acat, "BH"))]
  s <- d[q1 < fdr & q2 < fdr]
  if (!nrow(s)) return(list(genes = s, shares = NULL, n_genes = 0L))
  s[, class := paste0(data.table::fifelse(b1 > 0, "d+", "d-"), data.table::fifelse(b2 > 0, "f+", "f-"))]
  ## the label that matters biologically: a liability-raising change paired with a fitness benefit
  s[, type := data.table::fifelse(class %in% c("d+f+", "d-f-"), "antagonistic", "synergistic")]
  sh <- s[, .(n = .N, share = .N / nrow(s)), by = .(type, class)][order(type, class)]
  n_ant <- sum(s$type == "antagonistic")
  list(genes = s[], shares = sh[],
       antagonistic = n_ant / nrow(s),
       synergistic  = 1 - n_ant / nrow(s),
       n_genes = nrow(s),
       P = stats::binom.test(n_ant, nrow(s), 0.5)$p.value)
}
