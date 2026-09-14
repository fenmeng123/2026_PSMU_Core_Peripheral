# S6.6: Bootstrap intervals for S6.1 community-level bridge profiles
# Author: Kunru Song

options(stringsAsFactors = FALSE, scipen = 999)
try(Sys.setlocale("LC_CTYPE", ".UTF-8"), silent = TRUE)

required_packages <- c("bruceR", "dplyr", "openxlsx", "readxl", "tidyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing required R package(s): ", paste(missing_packages, collapse = ", "),
    ". Install them under R 4.5.2 before rerunning S6.6."
  )
}

resolve_script_file <- function(expected_basename) {
  candidates <- character(0)
  script_argument <- grep(
    "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
  )
  if (length(script_argument) == 1L) {
    candidates <- c(candidates, sub("^--file=", "", script_argument[[1L]]))
  }
  frame_files <- vapply(sys.frames(), function(frame) {
    frame_file <- frame[["ofile"]]
    if (is.null(frame_file) || !length(frame_file)) return(NA_character_)
    as.character(frame_file[[1L]])
  }, character(1))
  candidates <- c(candidates, frame_files[!is.na(frame_files)])
  candidates <- unique(candidates[nzchar(candidates)])
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates)) {
    candidates <- normalizePath(candidates, winslash = "/", mustWork = TRUE)
    candidates <- candidates[basename(candidates) == expected_basename]
  }
  # A same-named working-directory file must not override the invoked script.
  if (!length(candidates) && file.exists(expected_basename)) {
    candidates <- normalizePath(expected_basename, winslash = "/", mustWork = TRUE)
  }
  if (!length(candidates)) {
    stop(
      "Cannot locate ", expected_basename, ". Run it with Rscript or call ",
      "source() with the script path."
    )
  }
  candidates[[length(candidates)]]
}

script_file <- resolve_script_file("S6_6_BridgeProfileBootstrapIntervals.R")
script_dir <- dirname(script_file)
setwd(script_dir)
message("S6.6 working directory: ", normalizePath(".", winslash = "/"))

analysis_dir <- "../Res_2_Data_for_Analysis"
results_dir <- "../Res_3_Results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

network_file <- file.path(analysis_dir, "results_s5_3_mlvar_networks.rda")
s6_1_file <- file.path(
  analysis_dir, "results_s6_1_normalized_bridge_heatmaps.rda"
)
full_cluster_file <- file.path(
  analysis_dir, "checkpoint_s5_6_cluster_bootstrap.rda"
)
sparse_checkpoint_file <- file.path(
  analysis_dir, "checkpoint_s6_4_nonparametric_bootstrap_accuracy.rda"
)
result_file <- file.path(
  analysis_dir, "results_s6_6_bridge_profile_bootstrap_intervals.rda"
)
table_file <- file.path(
  results_dir, "Table_S6_6_BridgeProfileBootstrapIntervals.xlsx"
)

load_single_object <- function(path, expected_name) {
  if (!file.exists(path)) stop("Required input does not exist: ", path)
  input_environment <- new.env(parent = emptyenv())
  loaded_names <- load(path, envir = input_environment)
  if (!identical(loaded_names, expected_name)) {
    stop(
      "Expected sole object `", expected_name, "` in ", basename(path),
      "; found: ", paste(loaded_names, collapse = ", ")
    )
  }
  input_environment[[expected_name]]
}

network_results <- load_single_object(
  network_file, "results_s5_3_mlvar_networks"
)
s6_1_results <- load_single_object(
  s6_1_file, "results_s6_1_normalized_bridge_heatmaps"
)
full_cluster <- load_single_object(
  full_cluster_file, "checkpoint_s5_6_cluster_bootstrap"
)
sparse_checkpoint <- load_single_object(
  sparse_checkpoint_file, "checkpoint_s6_4_nonparametric_bootstrap_accuracy"
)

network_names <- c("Temporal", "Contemporaneous", "Between-person")
network_ids <- c(
  Temporal = "temporal", Contemporaneous = "contemporaneous",
  `Between-person` = "between"
)
expected_edge_counts <- c(
  Temporal = 144L, Contemporaneous = 66L, `Between-person` = 66L
)
group_levels <- c("Engagement Time", "Peripheral", "Core", "Distress")
group_display <- c(
  `Engagement Time` = "Engagement", Peripheral = "Peripheral",
  Core = "Core", Distress = "Distress"
)

node_map <- network_results$node_map |>
  dplyr::arrange(.data$order) |>
  dplyr::mutate(
    profile_group = dplyr::case_when(
      .data$community == "Engagement Time" ~ "Engagement Time",
      .data$community == "PSMU" & .data$subgroup == "Peripheral" ~ "Peripheral",
      .data$community == "PSMU" & .data$subgroup == "Core" ~ "Core",
      .data$community == "Distress" ~ "Distress",
      TRUE ~ NA_character_
    )
  )
node_order <- node_map$node
group_nodes <- stats::setNames(
  lapply(group_levels, function(group_name) {
    node_map$node[node_map$profile_group == group_name]
  }),
  group_levels
)

stopifnot(
  nrow(node_map) == 12L,
  !anyNA(node_map$profile_group),
  full_cluster$replications == 1000L,
  identical(dim(full_cluster$edge_estimates), c(1000L, 276L)),
  sparse_checkpoint$replications_per_resampling == 1000L,
  identical(sparse_checkpoint$raw_p_threshold, .05),
  all(unlist(sparse_checkpoint$validation)),
  identical(
    colnames(full_cluster$edge_estimates),
    names(sparse_checkpoint$baseline_full_edges)
  ),
  max(abs(full_cluster$baseline_edges - sparse_checkpoint$baseline_full_edges)) < 1e-12
)

edge_vector_to_matrix <- function(edge_values, network_name) {
  edge_keys <- names(edge_values)
  if (is.null(edge_keys)) stop("Edge values must retain their canonical names.")
  edge_parts <- strsplit(edge_keys, "|", fixed = TRUE)
  key_network <- vapply(edge_parts, `[[`, character(1), 1L)
  row_node <- vapply(edge_parts, `[[`, character(1), 2L)
  column_node <- vapply(edge_parts, `[[`, character(1), 3L)
  if (!all(key_network == network_name)) {
    stop("An edge vector contains keys from the wrong network: ", network_name)
  }
  matrix_object <- matrix(
    0, nrow = 12L, ncol = 12L,
    dimnames = list(node_order, node_order)
  )
  indices <- cbind(match(row_node, node_order), match(column_node, node_order))
  if (anyNA(indices)) stop("An edge key contains an unknown node.")
  matrix_object[indices] <- as.numeric(edge_values)
  if (network_name != "Temporal") {
    matrix_object[cbind(indices[, 2L], indices[, 1L])] <- as.numeric(edge_values)
  }
  matrix_object
}

# This exactly reproduces S6.1: signed arithmetic means, no absolute values,
# and no self-loops. Temporal rows are predictors at t-1 and columns are
# outcomes at t; undirected diagonal cells use each unique within-group edge.
aggregate_group_profile <- function(matrix_object, directed) {
  output <- matrix(
    NA_real_, nrow = 4L, ncol = 4L,
    dimnames = list(group_levels, group_levels)
  )
  for (row_group in group_levels) {
    for (column_group in group_levels) {
      row_nodes <- group_nodes[[row_group]]
      column_nodes <- group_nodes[[column_group]]
      block <- matrix_object[row_nodes, column_nodes, drop = FALSE]
      if (identical(row_group, column_group)) {
        if (directed) {
          values <- block[outer(row_nodes, column_nodes, FUN = "!=")]
        } else {
          values <- block[lower.tri(block, diag = FALSE)]
        }
      } else {
        values <- as.numeric(block)
      }
      output[row_group, column_group] <- mean(values)
    }
  }
  output
}

profiles_from_edge_estimates <- function(edge_estimates, network_name) {
  directed <- identical(network_name, "Temporal")
  output <- matrix(
    NA_real_, nrow = nrow(edge_estimates), ncol = 16L,
    dimnames = list(NULL, as.vector(outer(group_levels, group_levels, paste, sep = "|")))
  )
  for (replication in seq_len(nrow(edge_estimates))) {
    edge_vector <- edge_estimates[replication, ]
    names(edge_vector) <- colnames(edge_estimates)
    profile <- aggregate_group_profile(
      edge_vector_to_matrix(edge_vector, network_name), directed = directed
    )
    output[replication, ] <- as.vector(profile)
  }
  output
}

profile_from_edge_vector <- function(edge_vector, network_name) {
  aggregate_group_profile(
    edge_vector_to_matrix(edge_vector, network_name),
    directed = identical(network_name, "Temporal")
  )
}

original_support_from_s5_3 <- function(network_name, full_keys) {
  network_id <- unname(network_ids[[network_name]])
  edge_table <- network_results$edge_tables[[network_id]]
  if (network_name == "Temporal") {
    support_keys <- paste(
      network_name, edge_table$from, edge_table$to, sep = "|"
    )
  } else {
    support_keys <- paste(
      network_name, edge_table$node1, edge_table$node2, sep = "|"
    )
  }
  support <- stats::setNames(edge_table$supported_p05, support_keys)
  matched_support <- support[full_keys]
  if (anyNA(matched_support)) {
    stop("S5.3 original-sample support is incomplete for ", network_name, ".")
  }
  as.numeric(matched_support)
}

summarize_variant <- function(
    edge_estimates, point_edge_vector, network_name, method_id, method_label
) {
  bootstrap_profiles <- profiles_from_edge_estimates(edge_estimates, network_name)
  point_profile <- profile_from_edge_vector(point_edge_vector, network_name)
  cell_grid <- expand.grid(
    row_group = group_levels,
    column_group = group_levels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  data.frame(
    network = network_name,
    method_id = method_id,
    method = method_label,
    row_group = cell_grid$row_group,
    column_group = cell_grid$column_group,
    estimate = as.vector(point_profile),
    ci_lower = apply(
      bootstrap_profiles, 2L, stats::quantile,
      probs = .025, names = FALSE, type = 7L
    ),
    ci_upper = apply(
      bootstrap_profiles, 2L, stats::quantile,
      probs = .975, names = FALSE, type = 7L
    ),
    bootstrap_replications = nrow(bootstrap_profiles),
    stringsAsFactors = FALSE
  )
}

method_labels <- c(
  fixed_original_mask = "Raw-p thresholded (fixed original-sample mask)",
  replication_specific = "Replication-specific raw-p thresholded",
  full_unthresholded = "Full unthresholded"
)

summary_rows <- list()
bootstrap_method_audit <- list()
row_counter <- 0L
for (network_name in network_names) {
  full_keys <- names(sparse_checkpoint$baseline_full_edge_vectors[[network_name]])
  full_columns <- match(full_keys, colnames(full_cluster$edge_estimates))
  if (anyNA(full_columns)) stop("Full bootstrap edge columns are incomplete.")
  full_estimates <- full_cluster$edge_estimates[, full_columns, drop = FALSE]
  colnames(full_estimates) <- full_keys
  full_point <- sparse_checkpoint$baseline_full_edge_vectors[[network_name]]
  sparse_point <- sparse_checkpoint$baseline_sparse_edge_vectors[[network_name]]
  original_support <- original_support_from_s5_3(network_name, full_keys)

  fixed_estimates <- sweep(full_estimates, 2L, original_support, `*`)
  replication_specific_estimates <-
    sparse_checkpoint$cluster$sparse_edges_by_network[[network_name]]

  stopifnot(
    identical(dim(full_estimates), c(1000L, expected_edge_counts[[network_name]])),
    identical(colnames(replication_specific_estimates), full_keys),
    max(abs(sparse_point - full_point * original_support)) < 1e-12,
    all(fixed_estimates[, original_support == 0, drop = FALSE] == 0),
    all(is.finite(replication_specific_estimates))
  )

  variants <- list(
    fixed_original_mask = list(estimates = fixed_estimates, point = sparse_point),
    replication_specific = list(
      estimates = replication_specific_estimates, point = sparse_point
    ),
    full_unthresholded = list(estimates = full_estimates, point = full_point)
  )
  for (method_id in names(variants)) {
    row_counter <- row_counter + 1L
    summary_rows[[row_counter]] <- summarize_variant(
      variants[[method_id]]$estimates,
      variants[[method_id]]$point,
      network_name,
      method_id,
      unname(method_labels[[method_id]])
    )
  }
  bootstrap_method_audit[[network_name]] <- data.frame(
    network = network_name,
    edge_count = length(full_keys),
    original_supported_edges = sum(original_support),
    fixed_mask_bootstrap_replications = nrow(fixed_estimates),
    replication_specific_bootstrap_replications =
      nrow(replication_specific_estimates),
    full_bootstrap_replications = nrow(full_estimates),
    stringsAsFactors = FALSE
  )
}

intervals_long <- dplyr::bind_rows(summary_rows) |>
  dplyr::mutate(
    row_index = match(.data$row_group, group_levels),
    column_index = match(.data$column_group, group_levels),
    visible = .data$network == "Temporal" | .data$row_index >= .data$column_index,
    row_label = unname(group_display[.data$row_group]),
    column_label = unname(group_display[.data$column_group]),
    cell_label = dplyr::if_else(
      .data$network == "Temporal",
      paste0(.data$row_label, "\u2192", .data$column_label),
      paste0(.data$column_label, "-", .data$row_label)
    )
  ) |>
  dplyr::filter(.data$visible) |>
  dplyr::arrange(
    match(.data$network, network_names), .data$row_index, .data$column_index,
    match(.data$method_id, names(method_labels))
  )

fixed_point_profiles <- intervals_long |>
  dplyr::filter(.data$method_id == "fixed_original_mask")
s6_1_point_profiles <- s6_1_results$plot_data |>
  dplyr::filter(.data$visible) |>
  dplyr::mutate(
    network = dplyr::recode(
      .data$network,
      temporal = "Temporal", contemporaneous = "Contemporaneous",
      between = "Between-person"
    ),
    row_group = as.character(.data$row_group),
    column_group = as.character(.data$column_group)
  ) |>
  dplyr::select("network", "row_group", "column_group", s6_1_value = "value")
point_reproduction_audit <- fixed_point_profiles |>
  dplyr::select("network", "row_group", "column_group", recomputed = "estimate") |>
  dplyr::left_join(
    s6_1_point_profiles,
    by = c("network", "row_group", "column_group")
  ) |>
  dplyr::mutate(
    absolute_difference = abs(.data$recomputed - .data$s6_1_value),
    passed = .data$absolute_difference < 1e-12
  )
stopifnot(
  nrow(intervals_long) == 108L,
  nrow(point_reproduction_audit) == 36L,
  !anyNA(point_reproduction_audit$s6_1_value),
  all(point_reproduction_audit$passed),
  all(intervals_long$ci_lower <= intervals_long$ci_upper),
  all(intervals_long$bootstrap_replications == 1000L)
)

format_number <- function(value) {
  value[abs(value) < .0005] <- 0
  sprintf("%.3f", value)
}
format_estimate_ci <- function(estimate, lower, upper) {
  paste0(
    format_number(estimate), " [", format_number(lower), ", ",
    format_number(upper), "]"
  )
}

workbook_sheets <- stats::setNames(lapply(network_names, function(network_name) {
  network_data <- intervals_long |>
    dplyr::filter(.data$network == network_name) |>
    dplyr::mutate(
      formatted = format_estimate_ci(.data$estimate, .data$ci_lower, .data$ci_upper)
    ) |>
    dplyr::select("cell_label", "method_id", "formatted") |>
    tidyr::pivot_wider(names_from = "method_id", values_from = "formatted")
  network_data <- network_data[, c(
    "cell_label", "fixed_original_mask", "replication_specific",
    "full_unthresholded"
  )]
  names(network_data) <- c(
    "Cell",
    "Raw-p thresholded (fixed original mask)",
    "Replication-specific thresholded",
    "Full unthresholded"
  )
  as.data.frame(network_data, stringsAsFactors = FALSE)
}), network_names)

bruceR::export(workbook_sheets, file = table_file, overwrite = TRUE)

style_workbook <- function(path, sheet_data) {
  workbook <- openxlsx::loadWorkbook(path)
  header_style <- openxlsx::createStyle(
    fontName = "Aptos", fontSize = 10, textDecoration = "bold",
    fontColour = "#FFFFFF", fgFill = "#1F4E78",
    halign = "center", valign = "center", wrapText = TRUE,
    border = "bottom", borderColour = "#163A5C"
  )
  body_style <- openxlsx::createStyle(
    fontName = "Aptos", fontSize = 10, valign = "center", wrapText = TRUE
  )
  for (sheet_name in names(sheet_data)) {
    data <- sheet_data[[sheet_name]]
    rows <- nrow(data) + 1L
    columns <- ncol(data)
    openxlsx::showGridLines(workbook, sheet_name, showGridLines = FALSE)
    openxlsx::freezePane(workbook, sheet_name, firstRow = TRUE)
    openxlsx::setRowHeights(workbook, sheet_name, rows = 1L, heights = 34)
    openxlsx::addStyle(
      workbook, sheet_name, header_style, rows = 1L,
      cols = seq_len(columns), gridExpand = TRUE, stack = TRUE
    )
    openxlsx::addStyle(
      workbook, sheet_name, body_style, rows = 2L:rows,
      cols = seq_len(columns), gridExpand = TRUE, stack = TRUE
    )
    openxlsx::setColWidths(
      workbook, sheet_name, cols = 1L:columns,
      widths = c(24, 40, 36, 30)
    )
  }
  # bruceR/openxlsx can retain empty drawing relationships in data-only
  # workbooks. Remove those orphan parts so strict OOXML readers can import it.
  workbook$worksheets_rels <- lapply(
    workbook$worksheets_rels,
    function(relationships) {
      relationships[grepl("printerSettings", relationships, fixed = TRUE)]
    }
  )
  workbook$Content_Types <- workbook$Content_Types[
    !grepl("/xl/drawings/", workbook$Content_Types, fixed = TRUE)
  ]
  openxlsx::saveWorkbook(workbook, path, overwrite = TRUE)
}
style_workbook(table_file, workbook_sheets)

reloaded_sheet_names <- readxl::excel_sheets(table_file)
reloaded_sheets <- stats::setNames(lapply(reloaded_sheet_names, function(sheet) {
  as.data.frame(readxl::read_excel(table_file, sheet = sheet), stringsAsFactors = FALSE)
}), reloaded_sheet_names)
workbook_audit <- data.frame(
  check = c(
    "Exact sheet names and order", "Expected visible-cell row counts",
    "Exact four-column structure", "No missing cells", "No Excel formula errors"
  ),
  passed = c(
    identical(reloaded_sheet_names, network_names),
    identical(vapply(reloaded_sheets, nrow, integer(1)), c(
      Temporal = 16L, Contemporaneous = 10L, `Between-person` = 10L
    )),
    all(vapply(reloaded_sheets, ncol, integer(1)) == 4L),
    all(vapply(reloaded_sheets, function(data) !anyNA(data), logical(1))),
    !any(vapply(reloaded_sheets, function(data) {
      any(grepl("^#(REF!|VALUE!|NAME\\?|DIV/0!|N/A)", as.matrix(data)))
    }, logical(1)))
  ),
  stringsAsFactors = FALSE
)
stopifnot(all(workbook_audit$passed), file.info(table_file)$size > 0)

results_s6_6_bridge_profile_bootstrap_intervals <- list(
  schema_version = "1.0.0",
  created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  analysis_rule = list(
    cell_estimand = paste(
      "Signed arithmetic mean per possible edge; self-loops excluded;",
      "undirected diagonal cells use unique within-group edges"
    ),
    fixed_original_mask = paste(
      "Original-sample raw-p p < .05 support applied as a fixed mask to",
      "all 1,000 full participant-cluster bootstrap edge vectors"
    ),
    replication_specific = paste(
      "Each of 1,000 participant-cluster bootstrap replications thresholded",
      "by its own raw-p p < .05 support; undirected networks use the AND rule"
    ),
    full_unthresholded =
      "All 1,000 participant-cluster bootstrap edge vectors retained in full",
    interval = "Percentile 95% interval (2.5th and 97.5th percentiles; type 7)"
  ),
  intervals_long = intervals_long,
  workbook_sheets = workbook_sheets,
  bootstrap_method_audit = dplyr::bind_rows(bootstrap_method_audit),
  point_reproduction_audit = point_reproduction_audit,
  workbook_audit = workbook_audit,
  provenance = list(
    network_input = paste0("Res_2_Data_for_Analysis/", basename(network_file)),
    network_md5 = unname(tools::md5sum(network_file)),
    s6_1_input = paste0("Res_2_Data_for_Analysis/", basename(s6_1_file)),
    s6_1_md5 = unname(tools::md5sum(s6_1_file)),
    full_cluster_input = paste0(
      "Res_2_Data_for_Analysis/", basename(full_cluster_file)
    ),
    full_cluster_md5 = unname(tools::md5sum(full_cluster_file)),
    sparse_checkpoint_input = paste0(
      "Res_2_Data_for_Analysis/", basename(sparse_checkpoint_file)
    ),
    sparse_checkpoint_md5 = unname(tools::md5sum(sparse_checkpoint_file)),
    workbook_output = paste0("Res_3_Results/", basename(table_file)),
    workbook_md5 = unname(tools::md5sum(table_file)),
    generating_script = paste0(
      "Step_6_mlVAR_Reproduce_Centrality/", basename(script_file)
    ),
    r_version = R.version.string,
    package_versions = stats::setNames(
      vapply(required_packages, function(package_name) {
        as.character(utils::packageVersion(package_name))
      }, character(1)),
      required_packages
    )
  )
)
save(
  results_s6_6_bridge_profile_bootstrap_intervals,
  file = result_file, version = 2
)

cat("S6.6 completed successfully.\n")
cat("Saved workbook: ", table_file, "\n", sep = "")
cat("Saved result object: ", result_file, "\n", sep = "")
print(workbook_audit, row.names = FALSE)
