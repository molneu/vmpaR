# protivity

**protivity** is an R package that implements the **COMPASS** method for context-matched protein activity inference from transcriptomic data.

At its core, COMPASS treats the transcriptome as a sensor of perturbation state. A user-provided transcriptomic query is compared against a context-specific library of perturbation-derived reference signatures, and the resulting compatibility landscape is summarized at the level of targeted proteins and signatures.

This repository currently contains the development version of the package, including:

- internal subset resolution and cache scaffolding for context-specific reference databases
- core COMPASS scoring in `gsva` and `fgsea` modes
- structured result summarization via `compass_summary()`
- downstream annotation via `compass_annotate()`
- explicit comparison utilities via `compass_compare()` and `compass_compare_wide()`

---

## Conceptual overview

The current package structure follows three layers.

**Layer 1: reference resolution and gene-set construction**  
`compass()` resolves the context-specific subset database internally, reads the corresponding subset object, and builds one reference gene set per perturbation signature using the most downregulated genes in that signature.

**Layer 2: scoring**  
The user supplies either an expression matrix or a ranked gene-level vector, and `compass()` scores that query against the selected COMPASS reference library.

**Layer 3: interpretation and downstream analysis**  
The raw output can then be summarized, annotated, and compared using dedicated downstream functions.

---

## Current user-facing API

### `compass()`
Main user-facing scoring function.

This is the primary package entry point. It:

1. resolves the context-specific subset database
2. builds COMPASS reference gene sets internally
3. runs the selected analysis mode

The current supported modes are:

- `mode = "gsva"`
- `mode = "fgsea"`

### `compass_summary()`
Structured interpretation layer for COMPASS result objects.

This function summarizes COMPASS outputs into readable positive/negative activity blocks and optional downstream diagnostics such as RNA-activity discordance, variability, consistency, and shortlist-style prioritization.

### `compass_annotate()`
Optional downstream annotation utility.

This function augments COMPASS outputs with GO-based and CancerMine-based annotation layers.

### `compass_compare()`
Explicit group-wise comparison utility for COMPASS score matrices, typically after `mode = "gsva"`.

### `compass_compare_wide()`
Utility to pivot `compass_compare()` output into a wide comparison table for practical downstream inspection.

---

## Internal package structure

The current implementation uses internal helpers that are not intended as part of the public API.

In particular:

- `load_compass_db()` resolves local or cached subset files and contains the current cache/download scaffold
- internal gene-set construction is implemented in `R/compass_gene_set_building.R`
- internal mode-specific scoring logic is implemented in `R/compass_modes.R`

The former standalone `compass_gsc()` function is no longer part of the intended user-facing API. Gene-set construction is now handled internally by `compass()`.

---

## Current data model

### Context-specific subset files

The package currently works with one subset database per cancer context, stored as `.rds` files such as:

- `glioma_subset.rds`
- `breast_subset.rds`
- `nsclc_subset.rds`
- `headneck_subset.rds`

Supported contexts in the current code are:

- `breast`
- `crc`
- `gastric`
- `glioma`
- `headneck`
- `melanoma`
- `nsclc`
- `ovarian`
- `pdac`
- `prostate`

### Subset resolution and cache behavior

`compass()` resolves subset files internally through the package database-loading layer.

The current lookup order is:

1. an explicitly supplied `subset_dir`
2. the package cache directory
3. optional download logic if enabled

The cache location is based on:

```r
tools::R_user_dir("protivity", which = "cache")
```

and subset files are stored in the `subsets/` subdirectory of that cache.

At the moment, the download/caching architecture is scaffolded in the package, but remote URLs for actual subset downloads still need to be configured.

---

## How `compass()` currently builds reference gene sets

For the selected context, `compass()` currently applies the following internal logic:

1. load the context-specific subset object
2. filter signatures by minimum `cps_conf_total`
3. optionally restrict to selected targets
4. optionally restrict to signatures with non-`"None"` `cancer_driver_summary`
5. rank genes within each retained signature
6. take the bottom `n` genes per signature
7. use these gene sets as the reference library for the requested mode

The current default reference size is:

```r
n = 250
```

This means the current implementation uses the `250` most downregulated genes per retained reference signature by default.

---

## `compass()` modes

### `mode = "gsva"`

Input:

- numeric matrix or data frame
- genes in rows
- samples in columns

Output:

- numeric matrix of COMPASS scores
- reference signatures in rows
- samples in columns

Interpretation:

Each value reflects how compatible a sample is with a given perturbation-derived COMPASS reference signature.

### `mode = "fgsea"`

Input:

- named numeric vector
- names = gene symbols
- values = ranked gene-level statistics

Output:

- FGSEA-style result data frame
- including pathway, NES, p-value, adjusted p-value, and parsed confidence metadata

Interpretation:

Each COMPASS reference signature is tested as an enrichment set against the ranked query vector.

---

## What a COMPASS score means

A COMPASS score is not automatically a direct treatment-vs-control contrast.

Instead, it represents the compatibility of the observed transcriptomic profile with a selected perturbation-derived reference signature.

That means interpretation remains context-dependent:

- in `gsva` mode, scores are sample-level compatibility scores across the reference library
- in `fgsea` mode, each reference signature is tested against a ranked query vector

---

## Example usage

### GSVA mode

```r
compass_res <- compass(
  input = expr_mat,
  context = "glioma",
  mode = "gsva",
  subset_dir = "path/to/subsets",
  min_conf = 1,
  n = 250
)
```

### FGSEA mode

```r
compass_res <- compass(
  input = stats_vec,
  context = "glioma",
  mode = "fgsea",
  subset_dir = "path/to/subsets",
  min_conf = 1,
  n = 250
)
```

### Downstream summary

```r
summary_df <- compass_summary(compass_res)
```

### Downstream annotation

```r
annotated_df <- compass_annotate(summary_df)
```

### Group-wise comparison

```r
cmp <- compass_compare(
  compass_res = compass_res_matrix,
  sample_info = sample_info,
  sample_col = "sample",
  group_col = "group",
  mode = "pairwise",
  contrast = c("treated", "control")
)
```

---

## Development status

This repository is still in active development.

In particular, the following parts should be treated as evolving:

- final package metadata and documentation surface
- remote subset download configuration
- cache/distribution workflow for context-specific subset databases
- downstream interpretation heuristics and reporting style
