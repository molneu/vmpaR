#' Annotate COMPASS outputs with GO terms and CancerMine roles
#'
#' Adds downstream annotation layers to COMPASS-derived objects.
#'
#' Supported input types:
#' - data.frame-like COMPASS outputs with a `protein` column
#' - data.frame-like COMPASS outputs with a `pathway` column
#' - matrix-like COMPASS outputs with pathway names as rownames
#'
#' @param compass_obj A COMPASS result object:
#'   - data.frame with `protein` and/or `pathway`
#'   - matrix with pathway names as rownames
#' @param annotation_path Optional path to `compass_go_annotations.rds`.
#'   If `NULL`, the function checks `getOption("compass.go_annotation_path")`
#'   and otherwise falls back to `annotation/compass_go_annotations.rds`
#'   in the current working directory.
#' @param namespaces Optional character vector to filter GO namespaces.
#'   Allowed values:
#'   - `"biological_process"`
#'   - `"molecular_function"`
#'   - `"cellular_component"`
#'   Default: `NULL` (keep all).
#' @param evidence_codes Optional character vector of GO evidence codes to keep.
#'   Default: `NULL` (keep all).
#' @param collapse_terms Logical. If `TRUE`, collapse GO terms per protein into
#'   compact summary columns. Default: `TRUE`.
#' @param max_terms_per_namespace Integer. Maximum number of GO terms shown per
#'   namespace when `collapse_terms = TRUE`. Default: 8.
#' @param use_go Logical. If `TRUE`, attach GO-based annotation. Default: `TRUE`.
#' @param use_cancermine Logical. If `TRUE`, attach CancerMine-based annotation.
#'   Default: `TRUE`.
#' @param cancermine_path Optional path to `cancermine_collated.tsv`.
#'   If `NULL`, the function checks `getOption("compass.cancermine_path")`
#'   and otherwise falls back to
#'   `annotation/annotation_raw/cancermine_collated.tsv`.
#' @param context Optional COMPASS context. One of:
#'   `"breast"`, `"crc"`, `"gastric"`, `"glioma"`, `"headneck"`,
#'   `"melanoma"`, `"nsclc"`, `"ovarian"`, `"pdac"`, `"prostate"`.
#'   Used only for context-filtered CancerMine annotation.
#' @param cancermine_scope One of `"all"`, `"context"`, `"both"`.
#'   - `"all"`: pan-cancer CancerMine summary
#'   - `"context"`: only context-filtered CancerMine summary
#'   - `"both"`: both pan-cancer and context-filtered summaries
#'   Default: `"all"`.
#' @param cancermine_cancer Optional character vector. If provided, this
#'   overrides automatic context mapping and keeps only these
#'   `cancer_normalized` values.
#' @param cancermine_roles Optional character vector. If provided, keep only
#'   selected roles such as `"Oncogene"`, `"Tumor_Suppressor"`, `"Driver"`.
#'   Default: `NULL` (keep all roles).
#' @param cancermine_min_citations Integer. Minimum `citation_count` required to
#'   keep a CancerMine row. Default: 1.
#' @param max_cancermine_items Integer. Maximum number of roles/cancers shown
#'   in collapsed CancerMine annotation strings. Default: 5.
#'
#' @return
#' If `compass_obj` is a data.frame and `collapse_terms = TRUE`:
#'   the same data.frame plus annotation columns.
#'
#' If `compass_obj` is a matrix and `collapse_terms = TRUE`:
#'   a list with:
#'   - `data`: the original matrix
#'   - `row_annotation`: one row per matrix row/pathway
#'
#' If `collapse_terms = FALSE`:
#'   a list with:
#'   - `data`: the original object
#'   - `go_annotations`: the filtered long-format GO table (if requested)
#'   - `cancermine_annotations_all`: pan-cancer CancerMine rows (if requested)
#'   - `cancermine_annotations_context`: context-filtered CancerMine rows
#'     (if requested)
#'
#' @export
compass_annotate <- function(compass_obj,
                             annotation_path = NULL,
                             namespaces = NULL,
                             evidence_codes = NULL,
                             collapse_terms = TRUE,
                             max_terms_per_namespace = 8L,
                             use_go = TRUE,
                             use_cancermine = TRUE,
                             cancermine_path = NULL,
                             context = NULL,
                             cancermine_scope = c("all", "context", "both"),
                             cancermine_cancer = NULL,
                             cancermine_roles = NULL,
                             cancermine_min_citations = 1L,
                             max_cancermine_items = 5L) {
  cancermine_scope <- match.arg(cancermine_scope)
  
  if (!is.logical(collapse_terms) || length(collapse_terms) != 1L || is.na(collapse_terms)) {
    stop("`collapse_terms` must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.logical(use_go) || length(use_go) != 1L || is.na(use_go)) {
    stop("`use_go` must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.logical(use_cancermine) || length(use_cancermine) != 1L || is.na(use_cancermine)) {
    stop("`use_cancermine` must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!isTRUE(use_go) && !isTRUE(use_cancermine)) {
    stop("At least one of `use_go` or `use_cancermine` must be TRUE.", call. = FALSE)
  }
  
  if (!is.numeric(max_terms_per_namespace) ||
      length(max_terms_per_namespace) != 1L ||
      is.na(max_terms_per_namespace) ||
      max_terms_per_namespace <= 0) {
    stop("`max_terms_per_namespace` must be a single positive number.", call. = FALSE)
  }
  max_terms_per_namespace <- as.integer(max_terms_per_namespace)
  
  if (!is.numeric(max_cancermine_items) ||
      length(max_cancermine_items) != 1L ||
      is.na(max_cancermine_items) ||
      max_cancermine_items <= 0) {
    stop("`max_cancermine_items` must be a single positive number.", call. = FALSE)
  }
  max_cancermine_items <- as.integer(max_cancermine_items)
  
  if (!is.numeric(cancermine_min_citations) ||
      length(cancermine_min_citations) != 1L ||
      is.na(cancermine_min_citations) ||
      cancermine_min_citations < 0) {
    stop("`cancermine_min_citations` must be a single non-negative number.", call. = FALSE)
  }
  cancermine_min_citations <- as.integer(cancermine_min_citations)
  
  valid_namespaces <- c(
    "biological_process",
    "molecular_function",
    "cellular_component"
  )
  
  if (!is.null(namespaces)) {
    if (!is.character(namespaces)) {
      stop("`namespaces` must be NULL or a character vector.", call. = FALSE)
    }
    bad_ns <- setdiff(namespaces, valid_namespaces)
    if (length(bad_ns) > 0L) {
      stop(
        "Invalid `namespaces`: ",
        paste(bad_ns, collapse = ", "),
        "\nAllowed values: ",
        paste(valid_namespaces, collapse = ", "),
        call. = FALSE
      )
    }
  }
  
  if (!is.null(evidence_codes) && !is.character(evidence_codes)) {
    stop("`evidence_codes` must be NULL or a character vector.", call. = FALSE)
  }
  
  valid_contexts <- c(
    "breast", "crc", "gastric", "glioma", "headneck",
    "melanoma", "nsclc", "ovarian", "pdac", "prostate"
  )
  
  if (!is.null(context)) {
    if (!is.character(context) || length(context) != 1L || is.na(context) || !context %in% valid_contexts) {
      stop(
        "`context` must be NULL or one of: ",
        paste(valid_contexts, collapse = ", "),
        call. = FALSE
      )
    }
  }
  
  if (!is.null(cancermine_cancer) && !is.character(cancermine_cancer)) {
    stop("`cancermine_cancer` must be NULL or a character vector.", call. = FALSE)
  }
  
  if (!is.null(cancermine_roles) && !is.character(cancermine_roles)) {
    stop("`cancermine_roles` must be NULL or a character vector.", call. = FALSE)
  }
  
  if (cancermine_scope %in% c("context", "both") &&
      is.null(context) &&
      is.null(cancermine_cancer) &&
      isTRUE(use_cancermine)) {
    stop(
      "For `cancermine_scope = \"context\"` or `\"both\"`, please provide either `context` or `cancermine_cancer`.",
      call. = FALSE
    )
  }
  
  if (is.data.frame(compass_obj)) {
    df <- compass_obj
    
    if (!"protein" %in% colnames(df)) {
      if ("pathway" %in% colnames(df)) {
        df$protein <- .compass_annotate_derive_protein(df$pathway)
      } else {
        stop(
          "For data.frame input, `compass_obj` must contain either `protein` or `pathway`.",
          call. = FALSE
        )
      }
    }
    
    protein_vec <- df$protein
    base_df <- df
  } else if (is.matrix(compass_obj)) {
    if (is.null(rownames(compass_obj))) {
      stop(
        "For matrix input, `compass_obj` must have rownames representing pathways.",
        call. = FALSE
      )
    }
    
    protein_vec <- .compass_annotate_derive_protein(rownames(compass_obj))
    base_df <- data.frame(
      pathway = rownames(compass_obj),
      protein = protein_vec,
      stringsAsFactors = FALSE
    )
  } else {
    stop(
      "`compass_obj` must be either a data.frame or a matrix.",
      call. = FALSE
    )
  }
  
  go_ann <- NULL
  go_ann_collapsed <- NULL
  
  if (isTRUE(use_go)) {
    annotation_path <- .compass_resolve_go_annotation_path(annotation_path)
    
    if (!file.exists(annotation_path)) {
      stop(
        "GO annotation file not found: ", annotation_path,
        "\nPlease run `build_compass_go_annotations.R` first.",
        call. = FALSE
      )
    }
    
    go_ann <- readRDS(annotation_path)
    
    required_go_cols <- c(
      "protein", "go_id", "go_term", "namespace",
      "aspect", "evidence_code", "db_object_id", "source"
    )
    
    missing_go_cols <- setdiff(required_go_cols, colnames(go_ann))
    if (length(missing_go_cols) > 0L) {
      stop(
        "GO annotation file is missing required columns: ",
        paste(missing_go_cols, collapse = ", "),
        call. = FALSE
      )
    }
    
    if (!is.null(namespaces)) {
      go_ann <- go_ann[go_ann$namespace %in% namespaces, , drop = FALSE]
    }
    
    if (!is.null(evidence_codes)) {
      go_ann <- go_ann[go_ann$evidence_code %in% evidence_codes, , drop = FALSE]
    }
    
    rownames(go_ann) <- NULL
    
    if (isTRUE(collapse_terms)) {
      go_ann_collapsed <- .compass_collapse_go_annotations(
        go_ann = go_ann,
        max_terms_per_namespace = max_terms_per_namespace
      )
    }
  }
  
  cm_all_ann <- NULL
  cm_ctx_ann <- NULL
  cm_all_collapsed <- NULL
  cm_ctx_collapsed <- NULL
  
  if (isTRUE(use_cancermine)) {
    cancermine_path <- .compass_resolve_cancermine_path(cancermine_path)
    
    if (!file.exists(cancermine_path)) {
      stop(
        "CancerMine file not found: ", cancermine_path,
        call. = FALSE
      )
    }
    
    cm_raw <- .compass_read_cancermine_annotations(
      cancermine_path = cancermine_path,
      cancermine_roles = cancermine_roles,
      cancermine_min_citations = cancermine_min_citations
    )
    
    if (cancermine_scope %in% c("all", "both")) {
      cm_all_ann <- cm_raw
      if (isTRUE(collapse_terms)) {
        cm_all_collapsed <- .compass_collapse_cancermine_annotations(
          cm_ann = cm_all_ann,
          max_cancermine_items = max_cancermine_items
        )
        cm_all_collapsed <- .compass_prefix_annotation_columns(cm_all_collapsed, prefix = "cm_all")
      }
    }
    
    if (cancermine_scope %in% c("context", "both")) {
      ctx_cancers <- cancermine_cancer
      
      if (is.null(ctx_cancers)) {
        ctx_cancers <- .compass_context_to_cancermine_cancers(
          context = context,
          available_cancers = unique(cm_raw$cancer_normalized)
        )
      }
      
      cm_ctx_ann <- cm_raw[cm_raw$cancer_normalized %in% ctx_cancers, , drop = FALSE]
      rownames(cm_ctx_ann) <- NULL
      
      if (isTRUE(collapse_terms)) {
        cm_ctx_collapsed <- .compass_collapse_cancermine_annotations(
          cm_ann = cm_ctx_ann,
          max_cancermine_items = max_cancermine_items
        )
        cm_ctx_collapsed <- .compass_prefix_annotation_columns(cm_ctx_collapsed, prefix = "cm_ctx")
      }
    }
  }
  
  if (!isTRUE(collapse_terms)) {
    return(list(
      data = compass_obj,
      go_annotations = go_ann,
      cancermine_annotations_all = cm_all_ann,
      cancermine_annotations_context = cm_ctx_ann
    ))
  }
  
  annotation_tables <- list()
  
  if (!is.null(go_ann_collapsed)) {
    annotation_tables[[length(annotation_tables) + 1L]] <- go_ann_collapsed
  }
  
  if (!is.null(cm_all_collapsed)) {
    annotation_tables[[length(annotation_tables) + 1L]] <- cm_all_collapsed
  }
  
  if (!is.null(cm_ctx_collapsed)) {
    annotation_tables[[length(annotation_tables) + 1L]] <- cm_ctx_collapsed
  }
  
  anno_match <- .compass_build_annotation_match(
    protein_vec = protein_vec,
    annotation_tables = annotation_tables
  )
  
  if (is.data.frame(compass_obj)) {
    out <- cbind(base_df, anno_match)
    rownames(out) <- NULL
    return(out)
  }
  
  row_annotation <- cbind(base_df, anno_match)
  rownames(row_annotation) <- NULL
  
  list(
    data = compass_obj,
    row_annotation = row_annotation
  )
}

# Internal helpers -----------------------------------------------------------

.compass_resolve_go_annotation_path <- function(annotation_path = NULL) {
  if (!is.null(annotation_path)) {
    return(annotation_path)
  }
  
  opt_path <- getOption("compass.go_annotation_path", default = NULL)
  if (!is.null(opt_path)) {
    return(opt_path)
  }
  
  file.path(getwd(), "annotation", "compass_go_annotations.rds")
}

.compass_resolve_cancermine_path <- function(cancermine_path = NULL) {
  if (!is.null(cancermine_path)) {
    return(cancermine_path)
  }
  
  opt_path <- getOption("compass.cancermine_path", default = NULL)
  if (!is.null(opt_path)) {
    return(opt_path)
  }
  
  file.path(getwd(), "annotation", "annotation_raw", "cancermine_collated.tsv")
}

.compass_annotate_derive_protein <- function(pathways) {
  out <- ifelse(
    grepl("_XPR", pathways),
    sub("_XPR.*$", "", pathways),
    sub("_c[0-9]+$", "", pathways)
  )
  
  out[is.na(out)] <- ""
  out
}

.compass_build_annotation_match <- function(protein_vec,
                                            annotation_tables = list()) {
  if (length(annotation_tables) == 0L) {
    return(data.frame())
  }
  
  out_list <- lapply(annotation_tables, function(tbl) {
    matched <- tbl[
      match(protein_vec, tbl$protein),
      setdiff(colnames(tbl), "protein"),
      drop = FALSE
    ]
    rownames(matched) <- NULL
    matched
  })
  
  out <- do.call(cbind, out_list)
  rownames(out) <- NULL
  out
}

.compass_prefix_annotation_columns <- function(df, prefix) {
  if (nrow(df) == 0L) {
    return(df)
  }
  
  keep_names <- colnames(df)
  keep_names[keep_names != "protein"] <- paste0(prefix, "_", keep_names[keep_names != "protein"])
  colnames(df) <- keep_names
  df
}

.compass_collapse_go_annotations <- function(go_ann, max_terms_per_namespace = 8L) {
  if (nrow(go_ann) == 0L) {
    return(data.frame(
      protein = character(),
      n_go_total = integer(),
      n_go_bp = integer(),
      n_go_mf = integer(),
      n_go_cc = integer(),
      go_bp_terms = character(),
      go_mf_terms = character(),
      go_cc_terms = character(),
      go_all_terms = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  split_by_protein <- split(go_ann, go_ann$protein)
  
  collapsed_list <- lapply(names(split_by_protein), function(prot) {
    sub_df <- split_by_protein[[prot]]
    
    bp_terms <- unique(sub_df$go_term[sub_df$namespace == "biological_process"])
    mf_terms <- unique(sub_df$go_term[sub_df$namespace == "molecular_function"])
    cc_terms <- unique(sub_df$go_term[sub_df$namespace == "cellular_component"])
    
    bp_terms <- sort(bp_terms)
    mf_terms <- sort(mf_terms)
    cc_terms <- sort(cc_terms)
    
    bp_show <- utils::head(bp_terms, max_terms_per_namespace)
    mf_show <- utils::head(mf_terms, max_terms_per_namespace)
    cc_show <- utils::head(cc_terms, max_terms_per_namespace)
    
    all_terms <- unique(c(bp_terms, mf_terms, cc_terms))
    all_terms <- sort(all_terms)
    all_show <- utils::head(all_terms, max_terms_per_namespace)
    
    data.frame(
      protein = prot,
      n_go_total = length(unique(sub_df$go_id)),
      n_go_bp = length(bp_terms),
      n_go_mf = length(mf_terms),
      n_go_cc = length(cc_terms),
      go_bp_terms = if (length(bp_show) > 0L) paste(bp_show, collapse = "; ") else NA_character_,
      go_mf_terms = if (length(mf_show) > 0L) paste(mf_show, collapse = "; ") else NA_character_,
      go_cc_terms = if (length(cc_show) > 0L) paste(cc_show, collapse = "; ") else NA_character_,
      go_all_terms = if (length(all_show) > 0L) paste(all_show, collapse = "; ") else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, collapsed_list)
  out <- out[order(out$protein), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.compass_read_cancermine_annotations <- function(cancermine_path,
                                                 cancermine_roles = NULL,
                                                 cancermine_min_citations = 1L) {
  cm <- utils::read.delim(
    file = cancermine_path,
    sep = "\t",
    header = TRUE,
    quote = "",
    stringsAsFactors = FALSE
  )
  
  required_cols <- c(
    "matching_id",
    "role",
    "cancer_id",
    "cancer_normalized",
    "gene_hugo_id",
    "gene_entrez_id",
    "gene_normalized",
    "citation_count"
  )
  
  missing_cols <- setdiff(required_cols, colnames(cm))
  if (length(missing_cols) > 0L) {
    stop(
      "CancerMine file is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  cm$protein <- as.character(cm$gene_normalized)
  cm$role <- as.character(cm$role)
  cm$cancer_normalized <- as.character(cm$cancer_normalized)
  cm$citation_count <- suppressWarnings(as.numeric(cm$citation_count))
  
  cm <- cm[
    !is.na(cm$protein) & cm$protein != "" &
      !is.na(cm$role) & cm$role != "" &
      !is.na(cm$cancer_normalized) & cm$cancer_normalized != "" &
      is.finite(cm$citation_count),
    ,
    drop = FALSE
  ]
  
  cm <- cm[cm$citation_count >= cancermine_min_citations, , drop = FALSE]
  
  if (!is.null(cancermine_roles)) {
    cm <- cm[cm$role %in% cancermine_roles, , drop = FALSE]
  }
  
  cm <- cm[, c(
    "protein",
    "role",
    "cancer_id",
    "cancer_normalized",
    "gene_hugo_id",
    "gene_entrez_id",
    "citation_count",
    "matching_id"
  ), drop = FALSE]
  
  rownames(cm) <- NULL
  cm
}

.compass_collapse_cancermine_annotations <- function(cm_ann,
                                                     max_cancermine_items = 5L) {
  if (nrow(cm_ann) == 0L) {
    return(data.frame(
      protein = character(),
      cancermine_n_records = integer(),
      cancermine_total_citations = numeric(),
      cancermine_roles = character(),
      cancermine_top_role = character(),
      cancermine_top_role_citations = numeric(),
      cancermine_cancers = character(),
      cancermine_top_cancer = character(),
      cancermine_top_cancer_citations = numeric(),
      cancermine_oncogene_citations = numeric(),
      cancermine_tsg_citations = numeric(),
      cancermine_driver_citations = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  
  split_by_protein <- split(cm_ann, cm_ann$protein)
  
  collapsed_list <- lapply(names(split_by_protein), function(prot) {
    sub_df <- split_by_protein[[prot]]
    
    role_totals <- stats::aggregate(
      citation_count ~ role,
      data = sub_df,
      FUN = sum
    )
    role_totals <- role_totals[order(-role_totals$citation_count, role_totals$role), , drop = FALSE]
    
    cancer_totals <- stats::aggregate(
      citation_count ~ cancer_normalized,
      data = sub_df,
      FUN = sum
    )
    cancer_totals <- cancer_totals[order(-cancer_totals$citation_count, cancer_totals$cancer_normalized), , drop = FALSE]
    
    role_show <- utils::head(
      paste0(role_totals$role, "(", role_totals$citation_count, ")"),
      max_cancermine_items
    )
    
    cancer_show <- utils::head(
      paste0(cancer_totals$cancer_normalized, "(", cancer_totals$citation_count, ")"),
      max_cancermine_items
    )
    
    og_cit <- sum(sub_df$citation_count[sub_df$role == "Oncogene"], na.rm = TRUE)
    tsg_cit <- sum(sub_df$citation_count[sub_df$role == "Tumor_Suppressor"], na.rm = TRUE)
    drv_cit <- sum(sub_df$citation_count[sub_df$role == "Driver"], na.rm = TRUE)
    
    data.frame(
      protein = prot,
      cancermine_n_records = nrow(sub_df),
      cancermine_total_citations = sum(sub_df$citation_count, na.rm = TRUE),
      cancermine_roles = if (length(role_show) > 0L) paste(role_show, collapse = "; ") else NA_character_,
      cancermine_top_role = if (nrow(role_totals) > 0L) role_totals$role[1] else NA_character_,
      cancermine_top_role_citations = if (nrow(role_totals) > 0L) role_totals$citation_count[1] else NA_real_,
      cancermine_cancers = if (length(cancer_show) > 0L) paste(cancer_show, collapse = "; ") else NA_character_,
      cancermine_top_cancer = if (nrow(cancer_totals) > 0L) cancer_totals$cancer_normalized[1] else NA_character_,
      cancermine_top_cancer_citations = if (nrow(cancer_totals) > 0L) cancer_totals$citation_count[1] else NA_real_,
      cancermine_oncogene_citations = og_cit,
      cancermine_tsg_citations = tsg_cit,
      cancermine_driver_citations = drv_cit,
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, collapsed_list)
  out <- out[order(out$protein), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.compass_context_to_cancermine_cancers <- function(context, available_cancers) {
  patterns <- switch(
    context,
    breast = c("breast"),
    crc = c("colorectal", "colon", "rectal"),
    gastric = c("gastric", "stomach"),
    glioma = c(
      "glioma",
      "glioblastoma",
      "gliosarcoma",
      "oligodendroglioma",
      "brain cancer",
      "brain glioma",
      "brain stem glioma",
      "diffuse midline glioma",
      "high grade glioma",
      "low grade glioma",
      "optic nerve glioma",
      "mixed glioma"
    ),
    headneck = c(
      "head and neck",
      "oral",
      "orophary",
      "laryn",
      "hypophary",
      "nasophary"
    ),
    melanoma = c("melanoma"),
    nsclc = c("lung non-small cell carcinoma", "non-small cell lung"),
    ovarian = c("ovarian"),
    pdac = c("pancreatic ductal adenocarcinoma", "pancreatic cancer"),
    prostate = c("prostate"),
    NULL
  )
  
  if (is.null(patterns)) {
    return(character())
  }
  
  matched <- unique(unlist(lapply(patterns, function(pat) {
    grep(pat, available_cancers, value = TRUE, ignore.case = TRUE)
  }), use.names = FALSE))
  
  sort(unique(matched))
}