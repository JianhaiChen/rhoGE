## Does rho_ge() reproduce a published number?
## Takes the gene-level ACAT summaries of the disease-fertility study and recomputes rho_GE for
## Alzheimer's disease through the package. Needs the study's intermediate files, so it is a
## validation script rather than a unit test; it is skipped when they are not present.
suppressMessages({library(data.table); library(rhoGE)})
M <- Sys.getenv("RHOGE_VALIDATION_DIR", "/Volumes/X10Pro/data/synergistic/mr_res")
need <- file.path(M, c("gene_level_b1_disease_acat_summary.dropbox.csv",
                       "gene_level_b2_fitness_acat_summary.dropbox.csv", "gene_coords_v19.tsv"))
if (!all(file.exists(need))) { cat("validation data not present, skipping\n"); quit(save = "no") }
strip <- function(x) sub("[.].*", "", x)
co <- fread(need[3], col.names = c("g", "chr", "pos"))[, g := strip(g)]
b1 <- fread(need[1], select = c("gene", "disease", "top_beta", "top_se"))[, g := strip(gene)]
b2 <- fread(need[2], select = c("gene", "fitness", "top_beta", "top_se"))[, g := strip(gene)]
published <- list(male = -0.057, female = 0.016)
for (sx in c("male", "female")) {
  m <- merge(b1[disease == "Alzheimer2022", .(gene = g, b1 = top_beta, se1 = top_se)],
             b2[fitness == sx,              .(gene = g, b2 = top_beta, se2 = top_se)], by = "gene")
  m <- merge(m, co[, .(gene = g, chr, pos)], by = "gene")
  r <- rho_ge(m, chr = m$chr, pos = m$pos)        # z scores, MHC dropped, |z| < 9
  cat(sprintf("%-7s rho_GE = %+.4f (P = %.2g, n = %d); published %+.3f\n",
              sx, r$rho, r$P, r$n_genes, published[[sx]]))
  stopifnot(abs(r$rho - published[[sx]]) < 0.001)
}
cat("PASS: the package reproduces the published rho_GE\n")
