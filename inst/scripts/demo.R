## End-to-end example on the data shipped with the package: 25 chr22 genes with their real
## LD-pruned skeletal-muscle eQTL instruments and two small simulated GWAS. Runs in a few minutes.
library(rhoGE)
set.seed(1)   # FusioMR is a sampler; fix the seed so the example reproduces
x <- function(f) system.file("extdata", f, package = "rhoGE")
iv <- data.table::fread(x("instruments_chr22_muscle.tsv.gz"))
g1 <- read_gwas(x("gwas_trait1.tsv.gz")); g2 <- read_gwas(x("gwas_trait2.tsv.gz"))
fit <- fusio_pair(g1, g2, iv, control = FusioMR::parameter_control(niter = 4000))
genes <- aggregate_genes(fit)
co <- data.table::fread(x("gene_coords_v19.tsv.gz"), col.names = c("gene", "chr", "pos"))
genes <- merge(genes, co, by = "gene")
print(rho_ge(genes, block = genes$pos %/% 5e6, chr = genes$chr, pos = genes$pos))
cls <- direction_classes(genes)
print(cls$shares)
cat(sprintf("antagonistic %.1f%% of %d classified genes, binomial P = %.3g\n",
            100 * cls$antagonistic, cls$n_genes, cls$P))
