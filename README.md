# protivity

**protivity** is an R package that implements the **COMPASS** method for context-matched protein activity inference from transcriptomic data.

## Installation

```r
# install.packages("remotes")
remotes::install_github("landgrebe-a/protivity")
library(protivity)
```

## Main functions

- `compass()` runs COMPASS scoring on a user-provided expression matrix or ranked gene-level vector.
- `compass_gsc()` extracts COMPASS reference gene sets for a selected cancer context.
- `compass_summary()` summarizes COMPASS results for downstream interpretation.
- `compass_annotate()` adds optional annotation layers to COMPASS outputs.

## Reference

If you use **protivity**, please cite the COMPASS paper.

## Status

This package is currently under active development.
