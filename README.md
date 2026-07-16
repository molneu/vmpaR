# vmpaR

`vmpaR` implements **VMPA (Virtual Mapping of Active Proteomes)**, a method for context-specific, transcriptome-based inference of target-associated protein activity.

The package applies curated CRISPR loss-of-function signatures derived from the CMap/LINCS L1000 resource to either an expression matrix or a ranked gene-level vector. It provides two analysis backends:

- **GSVA** for sample-wise scoring of gene-by-sample expression data
- **FGSEA** for enrichment analysis of ranked gene-level statistics

> [!NOTE]
> `vmpaR` is under active development and is not yet available from CRAN or Bioconductor. For a reproducible analysis, install a specific package commit and report the corresponding Figshare dataset version.

## What VMPA does

For a selected cancer context, `vmpaR`:

1. loads the corresponding curated LINCS reference subset;
2. filters perturbation signatures by confidence, target, and optional cancer-driver annotation;
3. represents each retained CRISPR perturbation by its `n` most downregulated unique genes; and
4. scores those target-associated gene sets with GSVA or FGSEA.

The default gene-set size is `n = 250`. Multiple LINCS signatures can exist for the same target. With `unique = TRUE`, the package prioritizes signatures using the available reference metadata and confidence score. Equally prioritized GSVA signatures are averaged at target level; equally prioritized FGSEA signatures remain separate rows.

### Interpretation

VMPA scores are **context-conditioned activity proxies derived from transcriptional patterns**. They are not direct measurements of protein abundance, phosphorylation, enzymatic activity, or causal regulation.

- In the GSVA workflow, a higher score indicates stronger relative enrichment of the target-associated transcriptional program in a sample.
- In the FGSEA workflow, a positive normalized enrichment score (`NES`) indicates enrichment toward the positive end of the supplied ranking; a negative `NES` indicates enrichment toward the negative end. Its biological meaning therefore depends on how the ranking statistic was defined.

Scores should primarily be interpreted within a consistently processed dataset and with regard to the selected reference context.

## Installation

`vmpaR` requires R 4.3 or later.

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install the core data dependency, analysis backends, and example dependencies.
BiocManager::install(c("cmapR", "Biobase", "GSVA", "fgsea", "GSEABase"))

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("landgrebe-a/vmpaR")
```

Then load the package:

```r
library(vmpaR)
```

`GSVA` is required only for the GSVA backend, `fgsea` only for the FGSEA backend, and `GSEABase` only for exporting a `GeneSetCollection`. `Biobase` is used by the bundled example dataset.

The first analysis for a context downloads its reference subset from Figshare. Subsequent calls reuse the copy under `file.path(tools::R_user_dir("vmpaR", "cache"), "subsets")`.

## Quick start: GSVA

The package includes `kebir_gb`, a glioblastoma expression dataset stored as a `Biobase::ExpressionSet`.

```r
data("kebir_gb", package = "vmpaR")

expr_mat <- Biobase::exprs(kebir_gb)

gsva_result <- vmpa(
  input = expr_mat,
  context = "glioma",
  algorithm = "gsva",
  unique = TRUE,
  verbose = FALSE
)

gsva_result
```

The result is a `vmpa_result` data frame. With `unique = TRUE`, each row represents one target and the remaining columns contain one score per sample. Set `unique = FALSE` to retain individual reference signatures.

Optional post-GSVA scaling is available through `gsva_score_scaling`:

```r
gsva_scaled <- vmpa(
  input = expr_mat,
  context = "glioma",
  algorithm = "gsva",
  gsva_score_scaling = "sample_z",
  verbose = FALSE
)
```

Available modes are `"none"` (default), `"sample_z"`, `"sample_pop_sd"`, and `"signature_z"`.

## Quick start: FGSEA

FGSEA requires a named numeric vector of gene-level ranking statistics. The following mean-difference ranking is a minimal runnable example using the included data:

```r
sample_metadata <- Biobase::pData(kebir_gb)

relapsed <- sample_metadata$relapse_TYPE == "cR"
naive <- sample_metadata$relapse_TYPE == "n"

stopifnot(any(relapsed), any(naive))

stats_vec <- rowMeans(expr_mat[, relapsed, drop = FALSE]) -
  rowMeans(expr_mat[, naive, drop = FALSE])

stats_vec <- sort(stats_vec[is.finite(stats_vec)], decreasing = TRUE)

fgsea_result <- vmpa(
  input = stats_vec,
  context = "glioma",
  algorithm = "fgsea",
  unique = TRUE,
  seed = 123,
  verbose = FALSE
)

fgsea_result
```

This simple ranking demonstrates the API; it is not a substitute for a design-aware differential-expression model. For scientific analyses, supply an appropriate statistic that accounts for pairing, covariates, and the experimental design.

The FGSEA result contains target and signature metadata together with enrichment statistics such as `NES`, `pval`, `padj`, `ES`, `size`, and `log2err` when available. The default seed is applied locally around the FGSEA backend, so the caller's global random-number state is restored after the run. Set `seed = NULL` to disable package-controlled seeding.

## Input requirements

VMPA reference sets use human gene symbols. Input identifiers must therefore be unique, case-sensitive human gene symbols such as `TP53`, `EGFR`, or `AKT1`.

| Workflow | Required input |
|---|---|
| GSVA | Numeric matrix or data frame with genes in rows, samples in columns, gene symbols in row names, and sample identifiers in column names |
| FGSEA | Named numeric vector with gene symbols as names and gene-level ranking statistics as values |

The package does not automatically convert Ensembl or Entrez identifiers. Missing, blank, whitespace-padded, or duplicated input symbols are rejected. GSVA input must not contain `NA`, `NaN`, or infinite values; non-finite FGSEA statistics should be removed before analysis.

Gene-set overlap is evaluated after all reference filters have been applied. By default, a set must contain at least:

- one represented gene for GSVA (`gsva_min_size = 1L`); or
- ten represented genes for FGSEA (`fgsea_min_size = 10L`).

Sets below the relevant threshold are excluded. The analysis stops before calling the backend if no input gene matches the reference or no set reaches the selected minimum. With `verbose = TRUE`, `vmpa()` reports global overlap and the number of sets retained.

## Important arguments

| Argument | Default | Purpose |
|---|---:|---|
| `context` | required | Selects one of the supported context-specific reference subsets |
| `algorithm` | `"gsva"` | Chooses the GSVA or FGSEA workflow |
| `n` | `250L` | Number of bottom-ranked unique genes used per reference signature |
| `min_conf` | `1L` | Minimum rule-based confidence score; must be `1`, `2`, or `3` |
| `targets` | `NULL` | Restricts the analysis to selected target symbols |
| `driver_filter` | `FALSE` | Retains targets with any non-`"None"` cancer-driver annotation; this is not a canonical-driver-only filter |
| `unique` | `TRUE` | Reduces repeated reference signatures to target-level candidates |
| `gsva_score_scaling` | `"none"` | Controls optional post-GSVA score scaling |
| `gsva_min_size` | `1L` | Minimum overlapping genes per GSVA gene set |
| `fgsea_min_size` | `10L` | Minimum finite, ranked genes per FGSEA gene set |
| `seed` | `123L` | Locally controls FGSEA randomness; has no effect on GSVA |
| `return_gene_sets` | `FALSE` | Returns the result together with the gene sets and selection metadata actually used |

See `?vmpa` for the complete argument and return-value documentation.

## Inspect the gene sets used in a run

Set `return_gene_sets = TRUE` to obtain the exact post-filtering and post-overlap gene sets passed to the selected backend:

```r
full_result <- vmpa(
  input = expr_mat,
  context = "glioma",
  algorithm = "gsva",
  return_gene_sets = TRUE,
  verbose = FALSE
)

names(full_result)
head(full_result$vmpa_result)
head(names(full_result$gene_sets))
head(full_result$signature_metadata)
```

The returned list also records the selected context, algorithm, `unique` setting, and unique-selection metadata.

## Export reference gene sets

Use `vmpa_gsc()` to retrieve reference gene sets independently of an analysis run:

```r
gene_sets <- vmpa_gsc(
  context = "glioma",
  n = 250,
  min_conf = 1,
  unique = FALSE,
  output = "list",
  verbose = FALSE
)

length(gene_sets)
head(names(gene_sets))
```

Supported output formats are:

- `"list"`: named list of gene vectors;
- `"df"`: long-format data frame; and
- `"gsc"`: `GSEABase::GeneSetCollection`.

## Supported contexts

Each context currently maps to one LINCS reference cell line. The package does not infer the context from the query data.

| `context` | Reference cell line | Biological context |
|---|---|---|
| `glioma` | U251MG | Glioma / glioblastoma |
| `melanoma` | A375 | Melanoma |
| `nsclc` | A549 | Non-small cell lung cancer |
| `gastric` | AGS | Gastric cancer |
| `ovarian` | ES2 | Ovarian cancer |
| `crc` | HT29 | Colorectal cancer |
| `breast` | MCF7 | Breast cancer |
| `prostate` | PC3 | Prostate cancer |
| `pdac` | YAPC | Pancreatic ductal adenocarcinoma |
| `headneck` | BICR6 | Head and neck cancer |

Select the context on biological grounds and report it explicitly when presenting results. A context label should not be interpreted as covering all cell states or molecular subtypes of that cancer.

## Reference data and provenance

The context-specific VMPA reference subsets are downloaded automatically on first use and stored in the package cache. The curated dataset is published on Figshare:

- [VMPA signatures - Figshare dataset](https://doi.org/10.6084/m9.figshare.32060643)

The reference catalogue is derived from publicly accessible CMap/LINCS L1000 Level 5 CRISPR perturbation signatures and includes VMPA-specific selection and annotation. The original data were generated by the Connectivity Map at the Broad Institute as part of the NIH LINCS program.

Original source and policy information:

- [CMap/LINCS public BigQuery dataset](https://console.cloud.google.com/bigquery?p=cmap-big-table&d=cmap_lincs_public_views&page=dataset)
- [NIH LINCS data-release policy](https://github.com/dhimmel/integrate/blob/3b16651051ae12129ddc2250e8d3e6d4050dd349/licenses/custom/LINCS.md)

For a reproducible report, record the package commit, the version-specific Figshare DOI, the selected context, and all non-default arguments.

## Citation

If you use `vmpaR` or the VMPA reference signatures, cite the Figshare dataset and the original CMap/LINCS publication:

> Cima I, Landgrebe A, Sure U, Scheffler B. VMPA signatures. Figshare. Dataset. 2026. <https://doi.org/10.6084/m9.figshare.32060643.v4>

> Subramanian A, et al. A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles. *Cell*. 2017;171(6):1437-1452.e17. <https://doi.org/10.1016/j.cell.2017.10.049>

Please use the version-specific citation displayed on Figshare for a published analysis. The VMPA manuscript should be cited in addition once its final citation is available.

The package-generated software citation can be displayed with:

```r
citation("vmpaR")
```

## License

The `vmpaR` source code is licensed under the [GNU General Public License version 3](https://www.gnu.org/licenses/gpl-3.0.html). The reference data are distributed separately through Figshare under CC BY 4.0; consult the dataset record for provenance and reuse information.
