# protivity

**protivity** is an R package that implements the **COMPASS** method for context-matched protein activity inference from transcriptomic data.

## Installation

```r
# install.packages("remotes")
remotes::install_github("landgrebe-a/protivity")
library(protivity)
```

## Quick start

```r
data(kebir_gb, package = "protivity")

df <- Biobase::exprs(kebir_gb)
sample_metadata <- Biobase::pData(kebir_gb)

gsva_result_example <- compass(
  input = df,
  context = "glioma",
  algorithm = "gsva",
  unique = TRUE
)
```

Reference subsets are downloaded automatically on first use and cached locally for reuse.

`kebir_gb` is provided as a Biobase `ExpressionSet`. Use `Biobase::exprs(kebir_gb)` for the expression matrix and `Biobase::pData(kebir_gb)` for sample metadata.

Optional GSVA score normalization modes include `sample_z`, `sample_pop_sd`, and `signature_z`.

See `?compass` for an additional FGSEA example using a ranked gene-level vector.

## Main functions

- `compass()` runs COMPASS scoring on a user-provided expression matrix or ranked gene-level vector.
- `compass_gsc()` extracts COMPASS reference gene sets for a selected cancer context.

## Reference

If you use **protivity**, please cite the COMPASS paper.

## Status

This package is currently under active development.
