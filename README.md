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

gsva_result_example <- compass(
  input = df,
  context = "glioma",
  algorithm = "gsva"
)
```

Reference subsets are downloaded automatically on first use and cached locally for reuse.

## Main functions

- `compass()` runs COMPASS scoring on a user-provided expression matrix or ranked gene-level vector.
- `compass_gsc()` extracts COMPASS reference gene sets for a selected cancer context.

## Reference

If you use **protivity**, please cite the COMPASS paper.

## Status

This package is currently under active development.
