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
data(example_gbm_expr, package = "protivity")

gsva_result_example <- compass(
  input = example_gbm_expr,
  context = "glioma",
  algorithm = "gsva"
)
```

```r
data(example_u251_azd8055_rank, package = "protivity")

fgsea_result_example <- compass(
  input = example_u251_azd8055_rank,
  context = "glioma",
  algorithm = "fgsea",
  targets = "MTOR"
)
```

Reference subsets are downloaded automatically on first use and cached locally for reuse.

## Main functions

- `compass()` runs COMPASS scoring on a user-provided expression matrix or ranked gene-level vector.
- `compass_gsc()` extracts COMPASS reference gene sets for a selected cancer context.
- `compass_summary()` summarizes COMPASS results for downstream interpretation.
- `compass_annotate()` adds optional annotation layers to COMPASS outputs.

## Reference

If you use **protivity**, please cite the COMPASS paper.

## Status

This package is currently under active development.
