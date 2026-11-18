# rhoGE

Expression-effect genetic correlation (ρ_GE) between **any two traits**, from GWAS summary
statistics and cis-eQTLs.

For each gene and tissue, cis-eQTLs are used as instruments and the effects of genetically
regulated expression on both traits are estimated **jointly** with
[FusioMR](https://github.com/kangbw702/FusioMR) (`model = "semo"`). Tissue-level evidence is
combined per gene with the aggregated Cauchy association test (ACAT), and ρ_GE is the correlation
of the gene-level effects across genes. The same per-gene estimates give the four directional
classes (d+f+, d+f−, d−f+, d−f−) that separate antagonistic from synergistic pleiotropy.

## Install

```r
remotes::install_github("kangbw702/FusioMR")   # the MR engine
remotes::install_github("JianhaiChen/rhoGE")
```

Tested on macOS (R 4.3, Apple silicon) and Linux (R 4.2/4.3). Instrument building additionally
needs `BEDMatrix` and a PLINK reference panel.

## Use

```r
library(rhoGE)

g1 <- read_gwas("trait1.gz")          # column names are matched automatically
g2 <- read_gwas("trait2.gz")

## once per eQTL resource; see inst/scripts for building the full GTEx set
iv <- build_instruments(eqtl_table, bfile = "1kg_EUR", p_thresh = 1e-3, r2 = 0.1, min_iv = 5)

fit   <- fusio_pair(g1, g2, iv)       # per gene and tissue
genes <- aggregate_genes(fit)         # ACAT across tissues, one row per gene
rho   <- rho_ge(genes, block = genes$block)   # correlation, jackknife CI over genomic blocks
cls   <- direction_classes(genes)     # per gene: class and antagonistic/synergistic; shares and a binomial test
```

Command line:

```sh
Rscript inst/scripts/rhoge.R --gwas1 trait1.gz --gwas2 trait2.gz \
        --instruments gtex_v10_instruments.tsv.gz --out mypair
```

## Example

The package ships 25 chromosome-22 genes with their real skeletal-muscle instruments and two small
simulated GWAS, so the whole chain runs in a few minutes without any download:

```sh
Rscript -e 'source(system.file("scripts/demo.R", package = "rhoGE"))'
```

Gene coordinates for the jackknife blocks are included as
`system.file("extdata/gene_coords_v19.tsv.gz", package = "rhoGE")`.

## Choices worth knowing

- **Instruments** default to cis-eQTLs at *P* < 10⁻³, LD-pruned at r² < 0.1 within 100 kb, and a
  gene is kept only with ≥5 instruments. Building them is the slow step and depends only on the
  eQTL resource, so build once and reuse.
- **All tissues are used.** Regulatory effects are often tissue-specific and are frequently missed
  in the tissue one would pick a priori. ACAT returns one valid P value per gene under arbitrary
  correlation between tissues, so adding tissues does not multiply the number of tests.
- **Jackknife over blocks, not genes.** Neighbouring genes share instruments; pass `block` (for
  example 5-Mb bins) so the confidence interval reflects that.
- **The class names come from the disease–fitness setting**, where a liability-raising expression
  change paired with a fitness benefit is antagonistic. Rename them if that reading does not fit
  your traits.

## Prebuilt instruments

`gtex_v10_instruments.tsv.gz` (165 MB, release page) holds the instruments used in the paper:
5,351,666 cis-eQTLs over 50 GTEx v10 tissues, selected at *P* < 10⁻³ and LD-pruned at r² < 0.1
within 100 kb against a 1000 Genomes EUR reference, genes with ≥5 instruments only. Columns are
`gene SNP EA OA BETA SE P tissue`, ready to pass to `fusio_pair()`. Variants are keyed by rsID and
the alleles follow GRCh38; align your GWAS by rsID rather than by position.

## Citation

If you use this package, cite FusioMR (Kang et al., *Am J Hum Genet*, 2026) and the paper that
introduced ρ_GE in this form.
