#' Read and harmonise a GWAS summary-statistics file
#'
#' The file needs an rsID, both alleles, an effect estimate and its standard error. Column names are
#' matched case-insensitively against the usual spellings, so most public files are read without
#' renaming anything. Effects are returned for \code{ea}; nothing is flipped here, because alignment
#' happens against the eQTL instrument in \code{\link{fusio_pair}}.
#'
#' @param file path to a (optionally gzipped) table.
#' @param snp,ea,oa,beta,se,pval column names, if the automatic match fails.
#' @return a data.table with columns SNP, EA, OA, BETA, SE and (if present) P.
#' @export
read_gwas <- function(file, snp = NULL, ea = NULL, oa = NULL, beta = NULL, se = NULL, pval = NULL) {
  d <- data.table::fread(file, showProgress = FALSE)
  pick <- function(user, alts) {
    if (!is.null(user)) return(user)
    hit <- names(d)[tolower(names(d)) %in% alts]
    if (!length(hit)) stop("could not find a column for ", alts[1],
                           "; name it explicitly. Columns are: ", paste(names(d), collapse = ", "))
    hit[1]
  }
  cols <- c(SNP  = pick(snp,  c("snp", "rsid", "rsids", "variant_id", "markername")),
            EA   = pick(ea,   c("ea", "a1", "alt", "effect_allele", "allele1")),
            OA   = pick(oa,   c("oa", "a2", "ref", "other_allele", "allele0", "allele2")),
            BETA = pick(beta, c("beta", "b", "effect", "slope")),
            SE   = pick(se,   c("se", "sebeta", "standard_error", "slope_se")))
  out <- d[, cols, with = FALSE]
  data.table::setnames(out, names(cols))
  pcol <- tryCatch(pick(pval, c("p", "pval", "pvalue", "p_value", "p.value")), error = function(e) NA_character_)
  if (!is.na(pcol)) out[, P := d[[pcol]]]
  out[, SNP := sub(",.*", "", SNP)]
  out <- out[SNP != "" & is.finite(BETA) & is.finite(SE) & SE > 0]
  unique(out, by = "SNP")
}
