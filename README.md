# protivity

**protivity** is an R package for running **COMPASS**: context-matched protein activity inference from transcriptomic data.

COMPASS applies cancer-context-specific perturbation signatures to user-provided transcriptomic data. The main function, `compass()`, supports two workflows:

- **GSVA-based scoring** for gene-by-sample expression matrices
- **FGSEA-based enrichment** for ranked gene-level statistic vectors

## Installation

```r
# install.packages("remotes")
remotes::install_github("landgrebe-a/protivity")
```

```r
library(protivity)
```

## Optional dependencies

Some workflows require suggested Bioconductor packages such as `Biobase`, `GSVA`, `fgsea`, and `GSEABase`.

If needed, install them with:

```r
# install.packages("BiocManager")
BiocManager::install(c("Biobase", "GSVA", "fgsea", "GSEABase"))
```

## Quick start: GSVA workflow

The package includes `kebir_gb`, a small glioblastoma example dataset stored as a `Biobase::ExpressionSet`.

```r
data(kebir_gb, package = "protivity")

expr_mat <- Biobase::exprs(kebir_gb)
sample_metadata <- Biobase::pData(kebir_gb)

gsva_result <- compass(
  input = expr_mat,
  context = "glioma",
  algorithm = "gsva",
  unique = TRUE,
  verbose = FALSE
)

head(gsva_result)
```

For the GSVA workflow, `compass()` returns a `protivity_result` data frame. Rows correspond to targets/proteins when `unique = TRUE` and to individual COMPASS signatures when `unique = FALSE`. Score columns correspond to samples.

Optional GSVA score scaling modes include `"sample_z"`, `"sample_pop_sd"`, and `"signature_z"`.

## FGSEA workflow

For FGSEA, provide a named numeric vector of gene-level statistics or rankings.

```r
relapsed <- sample_metadata$relapse_TYPE != "n"
naive <- sample_metadata$relapse_TYPE == "n"

stats_vec <- rowMeans(expr_mat[, relapsed, drop = FALSE]) -
  rowMeans(expr_mat[, naive, drop = FALSE])

stats_vec <- stats_vec[is.finite(stats_vec) & !is.na(stats_vec)]
stats_vec <- sort(stats_vec, decreasing = TRUE)

fgsea_result <- compass(
  input = stats_vec,
  context = "glioma",
  algorithm = "fgsea",
  unique = TRUE,
  verbose = FALSE
)

head(fgsea_result)
```

For the FGSEA workflow, `compass()` returns a `protivity_result` data frame containing COMPASS metadata together with FGSEA enrichment statistics such as `target`, `conf`, `pathway`, `signature`, `NES`, `pval`, `padj`, `ES`, and `size`.

## Unique target-level output

COMPASS reference data may contain multiple perturbation signatures for the same target/protein. With `unique = TRUE`, protivity groups signatures by their target/protein annotation.

If validation metadata are available, validated signatures are prioritized. Among the remaining candidates, signatures with the highest `cps_conf_total` are retained.

For GSVA, all equally prioritized signatures for a target are scored separately and then averaged sample-wise, resulting in one row per target/protein.

For FGSEA, retained signatures are tested separately. If multiple equally prioritized signatures remain for a target/protein, they are returned as separate rows rather than being averaged or arbitrarily collapsed.

## Returning gene sets and metadata

By default, `compass()` returns the main COMPASS result object.

If you want to inspect the gene sets and signature metadata used for a run, set `return_gene_sets = TRUE`:

```r
gsva_full <- compass(
  input = expr_mat,
  context = "glioma",
  algorithm = "gsva",
  unique = TRUE,
  return_gene_sets = TRUE,
  verbose = FALSE
)

names(gsva_full)
```

This returns a list containing the COMPASS result, the gene sets used for scoring, signature metadata, and unique-selection metadata when applicable.

## Reference gene sets

To inspect or export the underlying COMPASS reference gene sets directly, use `compass_gsc()`:

```r
gene_sets <- compass_gsc(
  context = "glioma",
  output = "list",
  verbose = FALSE
)

length(gene_sets)
names(gene_sets)[1:5]
```

`compass_gsc()` can return gene sets as a named list, a data frame, or a `GSEABase::GeneSetCollection`.

```r
gene_sets_df <- compass_gsc(
  context = "glioma",
  output = "df",
  verbose = FALSE
)

if (requireNamespace("GSEABase", quietly = TRUE)) {
  gene_sets_gsc <- compass_gsc(
    context = "glioma",
    output = "gsc",
    verbose = FALSE
  )
}
```

## Main functions

- `compass()` runs COMPASS on a user-provided expression matrix or ranked gene-level vector.
- `compass_gsc()` retrieves COMPASS reference gene sets for a selected cancer context.

## Available contexts

Currently supported contexts are:

```r
c(
  "glioma",
  "melanoma",
  "nsclc",
  "gastric",
  "ovarian",
  "crc",
  "breast",
  "prostate",
  "pdac",
  "headneck"
)
```

## Citation

If you use **protivity** or COMPASS-derived scores in scientific work, please cite the COMPASS paper.

## Status

This package is currently under active development.
