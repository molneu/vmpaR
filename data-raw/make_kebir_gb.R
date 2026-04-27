# Create package example dataset: kebir_gb
# Source object: data-raw/ExpressionSet_GSE145128.rds

kebir_gb <- readRDS("data-raw/ExpressionSet_GSE145128.rds")

class(kebir_gb)
dim(Biobase::exprs(kebir_gb))
dim(Biobase::pData(kebir_gb))

usethis::use_data(
  kebir_gb,
  overwrite = TRUE,
  compress = "xz"
)