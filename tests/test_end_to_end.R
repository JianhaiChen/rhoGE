## Minimal end-to-end check with simulated data: no eQTL files or reference panel needed.
## Twenty genes, ten instruments each; the true expression effects on the two traits are correlated,
## so rho_GE should recover that correlation and the directional classes should be mostly same-signed.
suppressMessages({library(data.table); library(rhoGE)})
set.seed(42)
G <- 20; K <- 10
true <- data.table(gene = sprintf("G%02d", 1:G), a1 = rnorm(G, 0.3, 0.2))
true[, a2 := 0.6 * a1 + rnorm(G, 0, 0.12)]          # correlated effects on the two traits
iv <- rbindlist(lapply(1:G, function(i) data.table(
  gene = true$gene[i], tissue = "sim",
  SNP = sprintf("rs%04d", (i - 1) * K + 1:K),
  EA = "A", OA = "G",
  BETA = rnorm(K, 0.25, 0.05), SE = 0.02, P = 1e-8)))
mk <- function(a) {
  b <- iv$BETA * a[match(iv$gene, true$gene)] + rnorm(nrow(iv), 0, 0.01)
  data.table(SNP = iv$SNP, EA = "A", OA = "G", BETA = b, SE = 0.01)
}
g1 <- mk(true$a1); g2 <- mk(true$a2)
fit <- fusio_pair(g1, g2, iv, control = FusioMR::parameter_control(niter = 3000))
stopifnot(nrow(fit) == G)
## the simulated effects are estimated far more precisely than real ones, so the |z| outlier
## filter of the published analysis would remove every gene here
gen <- aggregate_genes(fit, z_max = Inf)
## simulated effects are estimated far more precisely than real ones, so the |z| filter is
## relaxed here for the same reason as in aggregate_genes above
r <- rho_ge(gen, z_max = Inf)
cat(sprintf("genes fitted: %d\n", nrow(gen)))
cat(sprintf("recovered effects correlate with the truth: trait 1 r = %.3f, trait 2 r = %.3f\n",
            cor(gen$b1, true$a1[match(gen$gene, true$gene)]),
            cor(gen$b2, true$a2[match(gen$gene, true$gene)])))
cat(sprintf("rho_GE (z scores) = %.3f (P = %.2g); correlation of the simulated effects = %.3f\n",
            r$rho, r$P, cor(true$a1, true$a2)))
cl <- direction_classes(gen, fdr = 0.05)
cat("directional classes:\n"); print(cl$shares)
stopifnot(r$rho > 0.5, nrow(gen) >= 0.8 * G)
cat("PASS\n")
