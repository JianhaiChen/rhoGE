#' Build cis-eQTL instruments for one tissue
#'
#' Selects the cis-eQTLs of each gene below \code{p_thresh} and prunes them by linkage disequilibrium
#' in a reference panel, keeping the most strongly associated variant of each clump. Genes left with
#' fewer than \code{min_iv} instruments are dropped, because the Mendelian randomization model is not
#' identified with too few. This is the step that takes the time; its result depends only on the eQTL
#' resource and the reference panel, so build it once and reuse it for every pair of traits.
#'
#' @param eqtl a table with columns gene, SNP, EA, OA, BETA, SE, P (or a file readable by
#'   \code{data.table::fread}), restricted to cis variants.
#' @param bfile PLINK prefix of the reference panel (same build as the eQTL positions).
#' @param p_thresh,r2,window instrument P value, clumping r-squared and window in base pairs.
#' @param min_iv minimum number of instruments a gene must retain.
#' @param tissue label written into the returned table, so instrument sets of several tissues can be
#'   combined and \code{fusio_pair} can keep them apart.
#' @return a data.table of retained instruments, with the same columns plus the tissue label.
#' @export
build_instruments <- function(eqtl, bfile, p_thresh = 1e-3, r2 = 0.1, window = 1e5,
                              min_iv = 5, tissue = NA_character_) {
  if (!requireNamespace("BEDMatrix", quietly = TRUE))
    stop("instrument building needs the BEDMatrix package")
  if (is.character(eqtl)) eqtl <- data.table::fread(eqtl, showProgress = FALSE)
  d <- data.table::as.data.table(eqtl)[P < p_thresh & is.finite(BETA) & SE > 0]
  if (!nrow(d)) return(d[0])
  keep <- d[, .N, by = gene][N >= min_iv, gene]
  d <- d[gene %in% keep]
  if (!nrow(d)) return(d[0])
  bim <- data.table::fread(paste0(bfile, ".bim"), header = FALSE,
                           col.names = c("chr", "SNP", "cm", "pos", "a1", "a2"), showProgress = FALSE)
  G <- BEDMatrix::BEDMatrix(bfile)
  shared <- intersect(unique(d$SNP), bim$SNP)
  if (length(shared) < min_iv) return(d[0])
  idx <- match(shared, bim$SNP)
  X <- G[, idx, drop = FALSE]; colnames(X) <- shared
  pos <- stats::setNames(bim$pos[idx], shared)
  clump <- function(snp, p) {
    ## keep p aligned with snp when variants absent from the reference are dropped, otherwise the
    ## clump is seeded by the wrong variant and the retained set is not LD-pruned at all
    ok <- snp %in% shared
    snp <- snp[ok]; p <- p[ok]
    if (!length(snp)) return(character(0))
    o <- order(p); snp <- snp[o]
    kept <- character(0); rest <- snp
    while (length(rest)) {
      top <- rest[1]; kept <- c(kept, top); rest <- rest[-1]
      if (!length(rest)) break
      near <- abs(pos[rest] - pos[top]) <= window
      if (any(near)) {
        r <- suppressWarnings(stats::cor(X[, top], X[, rest[near], drop = FALSE],
                                         use = "pairwise.complete.obs"))
        r[is.na(r)] <- 0
        rest <- setdiff(rest, rest[near][as.numeric(r)^2 > r2])
      }
    }
    kept
  }
  out <- d[, {
    k <- clump(SNP, P)
    if (length(k) >= min_iv) .SD[SNP %in% k] else .SD[0]
  }, by = gene]
  out[, tissue := tissue][]
}
