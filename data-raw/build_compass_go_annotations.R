# build_compass_go_annotations.R
# Build GO annotations for COMPASS targets only
# Output:
#   annotation/compass_go_annotations.rds
#   annotation/compass_go_unmatched_targets.csv
#
# Expected project structure:
#   COMPASS R Package/
#     build_compass_go_annotations.R
#     annotation/
#       annotation_raw/
#         goa_human.gaf.gz
#         go-basic.obo
#     subsets/                 OR
#     Data/subsets/

# -------------------------------------------------------------------------
# 0) Checks and paths
# -------------------------------------------------------------------------

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

annotation_dir <- file.path(project_root, "annotation")
annotation_raw_dir <- file.path(annotation_dir, "annotation_raw")

gaf_path <- file.path(annotation_raw_dir, "goa_human.gaf.gz")
obo_path <- file.path(annotation_raw_dir, "go-basic.obo")
out_path <- file.path(annotation_dir, "compass_go_annotations.rds")
unmatched_path <- file.path(annotation_dir, "compass_go_unmatched_targets.csv")

subset_dir_candidates <- c(
  file.path(project_root, "subsets"),
  file.path(project_root, "Data", "subsets")
)

subset_dir <- subset_dir_candidates[file.exists(subset_dir_candidates)][1]

if (is.na(subset_dir) || length(subset_dir) == 0L) {
  stop(
    "Could not find subset directory.\n",
    "Expected either:\n",
    "  - ", file.path(project_root, "subsets"), "\n",
    "  - ", file.path(project_root, "Data", "subsets"),
    call. = FALSE
  )
}

if (!file.exists(gaf_path)) {
  stop("Missing GAF file: ", gaf_path, call. = FALSE)
}

if (!file.exists(obo_path)) {
  stop("Missing OBO file: ", obo_path, call. = FALSE)
}

if (!requireNamespace("cmapR", quietly = TRUE)) {
  stop("Package `cmapR` must be installed to read COMPASS subset objects.", call. = FALSE)
}

subset_files <- list.files(
  subset_dir,
  pattern = "_subset\\.rds$",
  full.names = TRUE
)

if (length(subset_files) == 0L) {
  stop("No *_subset.rds files found in: ", subset_dir, call. = FALSE)
}

message("Project root: ", project_root)
message("Subset dir:   ", subset_dir)
message("Found ", length(subset_files), " subset files.")

# -------------------------------------------------------------------------
# 1) Extract COMPASS targets from subset RDS files
#    IMPORTANT:
#    Use cdesc$cmap_name as biological target label.
#    Do NOT parse cdesc$id blindly because technical IDs may contain
#    constructs / controls / non-gene labels.
# -------------------------------------------------------------------------

.clean_target_label <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- x[!is.na(x) & x != ""]
  x
}

extract_compass_targets <- function(subset_file) {
  gct <- readRDS(subset_file)
  
  if (!methods::is(gct, "GCT")) {
    stop("Object in file is not a GCT: ", subset_file, call. = FALSE)
  }
  
  cdesc <- gct@cdesc
  
  if (!"cmap_name" %in% colnames(cdesc)) {
    stop("Missing `cmap_name` in cdesc for file: ", subset_file, call. = FALSE)
  }
  
  targets <- .clean_target_label(cdesc$cmap_name)
  unique(targets)
}

compass_targets <- sort(unique(unlist(lapply(subset_files, extract_compass_targets), use.names = FALSE)))

message("Unique raw COMPASS targets from cmap_name: ", length(compass_targets))

# -------------------------------------------------------------------------
# 2) Read GOA human GAF
# -------------------------------------------------------------------------

gaf_colnames <- c(
  "DB",
  "DB_Object_ID",
  "DB_Object_Symbol",
  "Qualifier",
  "GO_ID",
  "DB_Reference",
  "Evidence_Code",
  "With_From",
  "Aspect",
  "DB_Object_Name",
  "DB_Object_Synonym",
  "DB_Object_Type",
  "Taxon",
  "Date",
  "Assigned_By",
  "Annotation_Extension",
  "Gene_Product_Form_ID"
)

go_gaf <- utils::read.delim(
  file = gzfile(gaf_path, open = "rt"),
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "!",
  fill = TRUE,
  stringsAsFactors = FALSE,
  col.names = gaf_colnames
)

message("Raw GOA rows: ", nrow(go_gaf))

# keep only protein entries
go_gaf <- go_gaf[go_gaf$DB_Object_Type == "protein", , drop = FALSE]

# remove explicit negation annotations
is_not <- grepl("(^|\\|)NOT($|\\|)", go_gaf$Qualifier)
go_gaf <- go_gaf[!is_not, , drop = FALSE]

# keep only rows with required fields
go_gaf <- go_gaf[
  !is.na(go_gaf$DB_Object_Symbol) & go_gaf$DB_Object_Symbol != "" &
    !is.na(go_gaf$GO_ID) & go_gaf$GO_ID != "",
  ,
  drop = FALSE
]

message("GOA rows after basic filtering: ", nrow(go_gaf))

# -------------------------------------------------------------------------
# 3) Restrict GOA to COMPASS targets
# -------------------------------------------------------------------------

matched_targets <- sort(intersect(compass_targets, unique(go_gaf$DB_Object_Symbol)))
unmatched_targets <- sort(setdiff(compass_targets, matched_targets))

message("Matched COMPASS targets in GOA:   ", length(matched_targets))
message("Unmatched COMPASS targets in GOA: ", length(unmatched_targets))

if (length(unmatched_targets) > 0L) {
  utils::write.csv(
    data.frame(target = unmatched_targets, stringsAsFactors = FALSE),
    file = unmatched_path,
    row.names = FALSE
  )
  message("Saved unmatched targets to: ", unmatched_path)
  message("Examples unmatched: ", paste(utils::head(unmatched_targets, 25), collapse = ", "))
}

go_gaf <- go_gaf[go_gaf$DB_Object_Symbol %in% matched_targets, , drop = FALSE]

message("Filtered GOA rows for matched COMPASS targets: ", nrow(go_gaf))

# -------------------------------------------------------------------------
# 4) Parse OBO: GO ID -> term name + namespace
# -------------------------------------------------------------------------

parse_go_obo <- function(obo_file) {
  lines <- readLines(obo_file, warn = FALSE, encoding = "UTF-8")
  
  term_starts <- which(lines == "[Term]")
  
  terms <- vector("list", length(term_starts))
  
  for (i in seq_along(term_starts)) {
    start <- term_starts[i]
    end <- if (i < length(term_starts)) term_starts[i + 1L] - 1L else length(lines)
    
    block <- lines[start:end]
    
    id_line <- grep("^id: GO:", block, value = TRUE)
    name_line <- grep("^name: ", block, value = TRUE)
    namespace_line <- grep("^namespace: ", block, value = TRUE)
    obsolete_line <- grep("^is_obsolete: true$", block, value = TRUE)
    
    if (length(id_line) == 0L || length(name_line) == 0L || length(namespace_line) == 0L) {
      next
    }
    
    if (length(obsolete_line) > 0L) {
      next
    }
    
    terms[[i]] <- data.frame(
      go_id = sub("^id: ", "", id_line[1]),
      go_term = sub("^name: ", "", name_line[1]),
      namespace = sub("^namespace: ", "", namespace_line[1]),
      stringsAsFactors = FALSE
    )
  }
  
  terms <- terms[!vapply(terms, is.null, logical(1))]
  
  out <- do.call(rbind, terms)
  out <- unique(out)
  rownames(out) <- NULL
  out
}

go_terms <- parse_go_obo(obo_path)

message("Parsed GO terms from OBO: ", nrow(go_terms))

# -------------------------------------------------------------------------
# 5) Build annotation table
# -------------------------------------------------------------------------

go_ann <- data.frame(
  protein = go_gaf$DB_Object_Symbol,
  go_id = go_gaf$GO_ID,
  aspect = go_gaf$Aspect,
  evidence_code = go_gaf$Evidence_Code,
  db_object_id = go_gaf$DB_Object_ID,
  source = "GOA_HUMAN",
  stringsAsFactors = FALSE
)

go_ann <- merge(
  x = go_ann,
  y = go_terms,
  by = "go_id",
  all.x = TRUE,
  sort = FALSE
)

go_ann <- go_ann[, c(
  "protein",
  "go_id",
  "go_term",
  "namespace",
  "aspect",
  "evidence_code",
  "db_object_id",
  "source"
), drop = FALSE]

go_ann <- go_ann[!is.na(go_ann$go_term) & go_ann$go_term != "", , drop = FALSE]
go_ann <- unique(go_ann)
go_ann <- go_ann[order(go_ann$protein, go_ann$namespace, go_ann$go_term, go_ann$evidence_code), , drop = FALSE]

rownames(go_ann) <- NULL

message("Final annotation rows: ", nrow(go_ann))
message("Annotated matched COMPASS targets: ", length(unique(go_ann$protein)))

# -------------------------------------------------------------------------
# 6) Save
# -------------------------------------------------------------------------

saveRDS(go_ann, file = out_path)

message("Saved: ", out_path)

# Optional quick preview
print(utils::head(go_ann, 10))