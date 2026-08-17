# vmpaR

`vmpaR` implements **VMPA (Virtual Mapping of Proteome Activity)**, a method for
context-specific, transcriptome-based inference of protein activity.

## Desktop application

VMPA is also available as a standalone desktop application for macOS Apple
Silicon and Windows 64-bit.

[Download VMPA desktop applications](https://doi.org/10.6084/m9.figshare.33028712)

The R package and desktop application use the same VMPA scoring engine. The
desktop application additionally provides graphical input preprocessing,
gene-identifier annotation, visualization, and downstream analysis workflows.

## Installation

`vmpaR` requires R 4.3 or later.

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("molneu/vmpaR")
```

Then load the package:

```r
library(vmpaR)
```

The first analysis for a context downloads its reference subset from Figshare.
Subsequent calls reuse the copy under:

```r
file.path(tools::R_user_dir("vmpaR", "cache"), "subsets")
```

## Citation

If you use vmpaR, either as R package or as desktop application, please cite our manuscript:

> Cima I, et al. Modeling the active proteome by context-matched perturbation
> signatures. Manuscript in preparation.

The package-generated citation can be displayed with:

```r
citation("vmpaR")
```

## What VMPA does

For a selected cancer context, `vmpaR`:

1. loads a curated perturbation reference subset;
2. filters the perturbation reference by confidence, target, and optional
   cancer-driver annotation;
3. represents each retained perturbation by its `n` most downregulated unique
   genes; and
4. scores those target-associated gene sets using either GSVA or FGSEA.

The default gene-set size is `n = 250`. Multiple signatures can exist for the
same target. With `unique = TRUE`, the package prioritizes signatures using the
available reference metadata and confidence score. Equally prioritized GSVA
signatures are averaged at target level; equally prioritized FGSEA signatures
remain separate rows.

## Quick start

The main function is `vmpa()`. The package includes an example dataset:

```r
data("kebir_gb", package = "vmpaR")

expr_mat <- Biobase::exprs(kebir_gb)

result <- vmpa(
  input = expr_mat,
  context = "glioma",
  algorithm = "gsva",
  unique = FALSE,
  driver_filter = TRUE,
  gsva_score_scaling = "sample_pop_sd",
  verbose = FALSE
)

result
```

The result is a `vmpa_result` data frame. With `unique = TRUE`, each row
represents one target and the remaining columns contain one score per sample.
Set `unique = FALSE` to retain individual reference signatures.

Available `gsva_score_scaling` modes are `"none"` (default), `"sample_z"`,
`"sample_pop_sd"`, and `"signature_z"`.

## Important arguments for `vmpa()`

| Argument | Default | Purpose |
|---|---:|---|
| `context` | required | Selects one of the supported context-specific reference subsets; see "Supported contexts" below |
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

Set `return_gene_sets = TRUE` to obtain the exact post-filtering and post-overlap
gene sets passed to the selected backend:

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

The returned list also records the selected context, algorithm, `unique`
setting, and unique-selection metadata.

## Export reference gene sets

Use `vmpa_gsc()` to retrieve reference gene sets independently of an analysis
run:

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

## Supported biological contexts

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

## Reference data and provenance

The context-specific VMPA reference subsets are downloaded automatically on
first use and stored in the package cache. The curated dataset is published on
Figshare:

- [VMPA signatures - Figshare dataset](https://doi.org/10.6084/m9.figshare.32060643)

The reference catalogue is derived from publicly accessible CMap/LINCS L1000
Level 5 CRISPR perturbation signatures and includes VMPA-specific selection and
annotation. The original data were generated by the Connectivity Map at the
Broad Institute as part of the NIH LINCS program.

Original source and policy information:

- [CMap/LINCS public BigQuery dataset](https://console.cloud.google.com/bigquery?p=cmap-big-table&d=cmap_lincs_public_views&page=dataset)
- [NIH LINCS data-release policy](https://github.com/dhimmel/integrate/blob/3b16651051ae12129ddc2250e8d3e6d4050dd349/licenses/custom/LINCS.md)

## License

The `vmpaR` source code is licensed under the
[GNU General Public License version 3](https://www.gnu.org/licenses/gpl-3.0.html).
The reference data are distributed separately through Figshare under CC BY 4.0;
consult the dataset record for provenance and reuse information.
