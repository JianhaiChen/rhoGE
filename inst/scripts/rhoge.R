#!/usr/bin/env Rscript
## Command-line entry point: two GWAS in, rho_GE and the directional classes out.
suppressMessages({library(data.table); library(rhoGE)})
a <- commandArgs(TRUE)
get <- function(flag, default = NULL) {
  i <- match(flag, a)
  if (is.na(i)) { if (is.null(default)) stop("missing ", flag) else return(default) }
  a[i + 1]
}
if (!length(a) || "--help" %in% a) {
  cat("Usage: rhoge.R --gwas1 F --gwas2 F --instruments F --out PREFIX [--niter 20000] [--zmax 9]\n")
  quit(save = "no")
}
g1 <- read_gwas(get("--gwas1")); g2 <- read_gwas(get("--gwas2"))
iv <- fread(get("--instruments"), showProgress = FALSE)
ctl <- FusioMR::parameter_control(niter = as.integer(get("--niter", "20000")),
                                  c_gamma = 4, c_theta = 6, kappa_gamma = 1, kappa_theta = 2)
fit <- fusio_pair(g1, g2, iv, control = ctl)
gen <- aggregate_genes(fit, z_max = as.numeric(get("--zmax", "9")))
out <- get("--out")
fwrite(fit, paste0(out, ".gene_tissue.tsv.gz"), sep = "\t")
fwrite(gen, paste0(out, ".gene.tsv.gz"), sep = "\t")
r <- rho_ge(gen); cl <- direction_classes(gen)
cat(sprintf("rho_GE = %.4f (95%% CI %.4f to %.4f, P = %.3g, %d genes)\n",
            r$rho, r$ci[1], r$ci[2], r$P, r$n_genes))
if (!is.null(cl$shares)) { cat("directional classes:\n"); print(cl$shares) }
fwrite(data.table(rho = r$rho, ci_low = r$ci[1], ci_high = r$ci[2], P = r$P, n_genes = r$n_genes,
                  antagonistic = cl$antagonistic, synergistic = cl$synergistic),
       paste0(out, ".rho_ge.tsv"), sep = "\t")
