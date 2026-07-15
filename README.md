# vmpaR

**vmpaR** is an R package for running **VMPA**: context-matched protein activity inference from transcriptomic data.

VMPA applies cancer-context-specific perturbation signatures to user-provided transcriptomic data. The main function, `vmpa()`, supports two workflows:

- **GSVA-based scoring** for gene-by-sample expression matrices
- **FGSEA-based enrichment** for ranked gene-level statistic vectors

## Installation

```r
# install.packages("remotes")
remotes::install_github("landgrebe-a/vmpaR")
```

```r
library(vmpaR)
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
data(kebir_gb, package = "vmpaR")

expr_mat <- Biobase::exprs(kebir_gb)
sample_metadata <- Biobase::pData(kebir_gb)

gsva_result <- vmpa(
  input = expr_mat,
  context = "glioma",
  algorithm = "gsva",
  unique = TRUE,
  verbose = FALSE
)

head(gsva_result)
```

For the GSVA workflow, `vmpa()` returns a `vmpa_result` data frame. Rows correspond to targets/proteins when `unique = TRUE` and to individual VMPA signatures when `unique = FALSE`. Score columns correspond to samples.

Optional GSVA score scaling modes include `"sample_z"`, `"sample_pop_sd"`, and `"signature_z"`.

## Gene identifiers and overlap

Input identifiers must be unique, case-sensitive human gene symbols matching
the symbols used by the VMPA/LINCS reference, for example `TP53`, `EGFR`, or
`AKT1`. Ensembl and Entrez identifiers are not converted automatically.
Missing, blank, whitespace-padded, and duplicated input symbols are rejected.

VMPA calculates the input overlap separately for every selected gene set. The
default thresholds preserve the original backend-specific workflows:

- `gsva_min_size = 1L`
- `fgsea_min_size = 10L`

With `verbose = TRUE`, `vmpa()` reports the number of matching input genes,
the number of gene sets retained by the relevant threshold, and the number
actually used by the backend. The analysis stops before the backend call if
there is no matching gene symbol or no selected gene set reaches the threshold.

If a symbol occurs more than once in a VMPA reference matrix, the occurrence
with the lowest reference-signature value is retained. Each symbol counts at
most once toward the requested gene-set size `n`.

## FGSEA workflow

For FGSEA, provide a named numeric vector of gene-level statistics or rankings.

```r
relapsed <- sample_metadata$relapse_TYPE != "n"
naive <- sample_metadata$relapse_TYPE == "n"

stats_vec <- rowMeans(expr_mat[, relapsed, drop = FALSE]) -
  rowMeans(expr_mat[, naive, drop = FALSE])

stats_vec <- stats_vec[is.finite(stats_vec) & !is.na(stats_vec)]
stats_vec <- sort(stats_vec, decreasing = TRUE)

fgsea_result <- vmpa(
  input = stats_vec,
  context = "glioma",
  algorithm = "fgsea",
  unique = TRUE,
  verbose = FALSE
)

head(fgsea_result)
```

For the FGSEA workflow, `vmpa()` returns a `vmpa_result` data frame containing VMPA metadata together with FGSEA enrichment statistics such as `target`, `conf`, `pathway`, `signature`, `NES`, `pval`, `padj`, `ES`, and `size`.

## Unique target-level output

VMPA reference data may contain multiple perturbation signatures for the same target/protein. With `unique = TRUE`, vmpaR groups signatures by their target/protein annotation.

If validation metadata are available, validated signatures are prioritized. Among the remaining candidates, signatures with the highest `cps_conf_total` are retained.

For GSVA, all equally prioritized signatures for a target are scored separately and then averaged sample-wise, resulting in one row per target/protein.

For FGSEA, retained signatures are tested separately. If multiple equally prioritized signatures remain for a target/protein, they are returned as separate rows rather than being averaged or arbitrarily collapsed.

## Returning gene sets and metadata

By default, `vmpa()` returns the main VMPA result object.

If you want to inspect the gene sets and signature metadata used for a run, set `return_gene_sets = TRUE`:

```r
gsva_full <- vmpa(
  input = expr_mat,
  context = "glioma",
  algorithm = "gsva",
  unique = TRUE,
  return_gene_sets = TRUE,
  verbose = FALSE
)

names(gsva_full)
```

This returns a list containing the VMPA result, the gene sets used for scoring, signature metadata, and unique-selection metadata when applicable.

## Reference gene sets

To inspect or export the underlying VMPA reference gene sets directly, use `vmpa_gsc()`:

The context-specific VMPA reference subsets are published on Figshare:
<https://doi.org/10.6084/m9.figshare.32060643>.

```r
gene_sets <- vmpa_gsc(
  context = "glioma",
  output = "list",
  verbose = FALSE
)

length(gene_sets)
names(gene_sets)[1:5]
```

`vmpa_gsc()` can return gene sets as a named list, a data frame, or a `GSEABase::GeneSetCollection`.

```r
gene_sets_df <- vmpa_gsc(
  context = "glioma",
  output = "df",
  verbose = FALSE
)

if (requireNamespace("GSEABase", quietly = TRUE)) {
  gene_sets_gsc <- vmpa_gsc(
    context = "glioma",
    output = "gsc",
    verbose = FALSE
  )
}
```

## Main functions

- `vmpa()` runs VMPA on a user-provided expression matrix or ranked gene-level vector.
- `vmpa_gsc()` retrieves VMPA reference gene sets for a selected cancer context.

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

If you use **vmpaR** or VMPA-derived scores in scientific work, please cite the VMPA paper.

## Status

This package is currently under active development.
