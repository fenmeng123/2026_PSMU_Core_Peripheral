# S6.8: A priori paired-bootstrap difference tests
# Author: Kunru Song

options(stringsAsFactors = FALSE, scipen = 999)
try(Sys.setlocale("LC_CTYPE", ".UTF-8"), silent = TRUE)

required_packages <- c("bruceR", "dplyr", "openxlsx", "readxl")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing required R package(s): ", paste(missing_packages, collapse = ", "),
    ". Install them under R 4.5.2 before rerunning S6.8."
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

script_file <- resolve_script_file("S6_8_APrioriBootstrapDifferenceTests.R")
script_dir <- dirname(script_file)
setwd(script_dir)
message("S6.8 working directory: ", normalizePath(".", winslash = "/"))

analysis_dir <- "../Res_2_Data_for_Analysis"
results_dir <- "../Res_3_Results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

s6_2_file <- file.path(
  analysis_dir, "results_s6_2_node_bridge_expected_influence.rda"
)
s5_3_file <- file.path(
  analysis_dir, "results_s5_3_mlvar_networks.rda"
)
cluster_checkpoint_file <- file.path(
  analysis_dir, "checkpoint_s5_6_cluster_bootstrap.rda"
)
result_file <- file.path(
  analysis_dir, "results_s6_8_apriori_bootstrap_difference_tests.rda"
)
table_file <- file.path(
  results_dir, "Table_S6_8_APrioriBootstrapDifferenceTests.xlsx"
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

s6_2 <- load_single_object(
  s6_2_file, "results_s6_2_node_bridge_expected_influence"
)
s5_3 <- load_single_object(
  s5_3_file, "results_s5_3_mlvar_networks"
)
cluster_checkpoint <- load_single_object(
  cluster_checkpoint_file, "checkpoint_s5_6_cluster_bootstrap"
)

point_estimates <- s6_2$point_estimates
bridge_bootstrap <- s6_2$bootstrap_estimates
edge_bootstrap <- cluster_checkpoint$edge_estimates
planned_edge_contrasts <- s5_3$planned_contrasts

required_point_columns <- c(
  "panel_id", "network", "metric", "target_id", "node", "node_display",
  "key", "estimate"
)
stopifnot(
  nrow(point_estimates) == 48L,
  all(required_point_columns %in% names(point_estimates)),
  identical(dim(bridge_bootstrap), c(1000L, 48L)),
  identical(dim(edge_bootstrap), c(1000L, 276L)),
  nrow(planned_edge_contrasts) == 6L,
  all(is.finite(point_estimates$estimate)),
  all(is.finite(bridge_bootstrap)),
  all(is.finite(edge_bootstrap)),
  !anyDuplicated(point_estimates$key),
  setequal(point_estimates$key, colnames(bridge_bootstrap))
)

current_s5_3_md5 <- unname(tools::md5sum(s5_3_file))
current_cluster_md5 <- unname(tools::md5sum(cluster_checkpoint_file))
stopifnot(
  identical(current_s5_3_md5, s6_2$provenance$network_md5),
  identical(current_cluster_md5, s6_2$provenance$cluster_checkpoint_md5),
  identical(colnames(edge_bootstrap), names(cluster_checkpoint$baseline_edges))
)

panel_order <- c(
  "temporal_in", "temporal_out", "contemporaneous", "between"
)
target_order <- c("engagement", "distress")
node_order <- paste0("PMU", 1:6)
node_display_map <- stats::setNames(paste0("P", 1:6), node_order)
panel_labels <- data.frame(
  panel_id = panel_order,
  network = c("Temporal", "Temporal", "Contemporaneous", "Between-person"),
  metric = c("Bridge In-EI", "Bridge Out-EI", "Bridge EI", "Bridge EI"),
  stringsAsFactors = FALSE
)
target_labels <- c(
  engagement = "Engagement Community", distress = "Distress Community"
)

point_lookup <- stats::setNames(point_estimates$estimate, point_estimates$key)
bridge_bootstrap <- bridge_bootstrap[, point_estimates$key, drop = FALSE]

confidence_level <- 0.999
alpha <- 1 - confidence_level
interval_probabilities <- c(alpha / 2, 1 - alpha / 2)
quantile_type <- 6L
bootstrap_replications <- nrow(bridge_bootstrap)
maximum_directly_resolved_confidence <- 1 - 2 / bootstrap_replications
tail_resolution_sufficient <- confidence_level <=
  maximum_directly_resolved_confidence

summarize_contrast_matrix <- function(contrast_matrix, specification) {
  stopifnot(
    is.matrix(contrast_matrix),
    nrow(contrast_matrix) == bootstrap_replications,
    ncol(contrast_matrix) == nrow(specification),
    identical(colnames(contrast_matrix), specification$contrast_id),
    all(is.finite(contrast_matrix)),
    all(is.finite(specification$original_difference))
  )
  intervals <- t(apply(
    contrast_matrix,
    2L,
    stats::quantile,
    probs = interval_probabilities,
    type = quantile_type,
    names = FALSE,
    na.rm = FALSE
  ))
  colnames(intervals) <- c("ci_lower_99_9", "ci_upper_99_9")
  minimums <- apply(contrast_matrix, 2L, min)
  maximums <- apply(contrast_matrix, 2L, max)

  specification |>
    dplyr::mutate(
      bootstrap_median = apply(contrast_matrix, 2L, stats::median),
      ci_lower_99_9 = intervals[.data$contrast_id, "ci_lower_99_9"],
      ci_upper_99_9 = intervals[.data$contrast_id, "ci_upper_99_9"],
      ci_equals_observed_envelope =
        abs(.data$ci_lower_99_9 - minimums[.data$contrast_id]) < 1e-15 &
        abs(.data$ci_upper_99_9 - maximums[.data$contrast_id]) < 1e-15,
      excludes_zero_99_9 =
        .data$ci_lower_99_9 > 0 | .data$ci_upper_99_9 < 0,
      exploratory_decision = dplyr::if_else(
        .data$excludes_zero_99_9,
        "Excludes zero: strong exploratory evidence",
        "Includes zero: insufficient evidence of a difference"
      )
    )
}

# -----------------------------------------------------------------------------
# Family 1: 24 within-node target contrasts (Distress minus Engagement)
# -----------------------------------------------------------------------------
node_target_specification <- expand.grid(
  panel_id = panel_order,
  node = node_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) |>
  dplyr::left_join(panel_labels, by = "panel_id") |>
  dplyr::arrange(match(.data$panel_id, panel_order), match(.data$node, node_order)) |>
  dplyr::mutate(
    contrast_id = sprintf("NT%02d", dplyr::row_number()),
    node_label = unname(node_display_map[.data$node]),
    contrast_direction = "Distress target - Engagement target",
    positive_term = paste(.data$node_label, "Distress-target", .data$metric),
    negative_term = paste(.data$node_label, "Engagement-target", .data$metric),
    contrast_label = paste0(
      .data$node_label, ": Distress Community - Engagement Community"
    ),
    distress_key = paste(.data$panel_id, "distress", .data$node, sep = "|"),
    engagement_key = paste(.data$panel_id, "engagement", .data$node, sep = "|"),
    original_difference =
      unname(point_lookup[.data$distress_key]) -
      unname(point_lookup[.data$engagement_key])
  )

node_target_bootstrap <- matrix(
  NA_real_, nrow = bootstrap_replications,
  ncol = nrow(node_target_specification),
  dimnames = list(rownames(bridge_bootstrap), node_target_specification$contrast_id)
)
for (index in seq_len(nrow(node_target_specification))) {
  specification_row <- node_target_specification[index, ]
  node_target_bootstrap[, specification_row$contrast_id] <-
    bridge_bootstrap[, specification_row$distress_key] -
    bridge_bootstrap[, specification_row$engagement_key]
}
node_target_results <- summarize_contrast_matrix(
  node_target_bootstrap, node_target_specification
)

# -----------------------------------------------------------------------------
# Family 2: 64 pairwise Core-minus-Peripheral node contrasts
# -----------------------------------------------------------------------------
peripheral_nodes <- paste0("PMU", 1:2)
core_nodes <- paste0("PMU", 3:6)
core_peripheral_specification <- expand.grid(
  panel_id = panel_order,
  target_id = target_order,
  peripheral_node = peripheral_nodes,
  core_node = core_nodes,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) |>
  dplyr::left_join(panel_labels, by = "panel_id") |>
  dplyr::arrange(
    match(.data$panel_id, panel_order),
    match(.data$target_id, target_order),
    match(.data$peripheral_node, peripheral_nodes),
    match(.data$core_node, core_nodes)
  ) |>
  dplyr::mutate(
    contrast_id = sprintf("CP%02d", dplyr::row_number()),
    target_community = unname(target_labels[.data$target_id]),
    peripheral_label = unname(node_display_map[.data$peripheral_node]),
    core_label = unname(node_display_map[.data$core_node]),
    contrast_direction = "Core node - Peripheral node",
    positive_term = paste(.data$core_label, .data$target_community, .data$metric),
    negative_term = paste(
      .data$peripheral_label, .data$target_community, .data$metric
    ),
    contrast_label = paste0(
      .data$core_label, " - ", .data$peripheral_label,
      " (", .data$target_community, ")"
    ),
    core_key = paste(
      .data$panel_id, .data$target_id, .data$core_node, sep = "|"
    ),
    peripheral_key = paste(
      .data$panel_id, .data$target_id, .data$peripheral_node, sep = "|"
    ),
    original_difference =
      unname(point_lookup[.data$core_key]) -
      unname(point_lookup[.data$peripheral_key])
  )

core_peripheral_bootstrap <- matrix(
  NA_real_, nrow = bootstrap_replications,
  ncol = nrow(core_peripheral_specification),
  dimnames = list(
    rownames(bridge_bootstrap), core_peripheral_specification$contrast_id
  )
)
for (index in seq_len(nrow(core_peripheral_specification))) {
  specification_row <- core_peripheral_specification[index, ]
  core_peripheral_bootstrap[, specification_row$contrast_id] <-
    bridge_bootstrap[, specification_row$core_key] -
    bridge_bootstrap[, specification_row$peripheral_key]
}
core_peripheral_results <- summarize_contrast_matrix(
  core_peripheral_bootstrap, core_peripheral_specification
)

# -----------------------------------------------------------------------------
# Family 3: the six original S5.3 Temporal edge-family contrasts
# -----------------------------------------------------------------------------
split_nodes <- function(value) strsplit(value, ", ", fixed = TRUE)[[1L]]
node_label_lookup <- stats::setNames(s5_3$node_map$label, s5_3$node_map$node)
node_short_lookup <- stats::setNames(
  c("E1", "E2", "E3", paste0("P", 1:6), "D3", "D2", "D1"),
  s5_3$node_map$node
)

temporal_edge_specification <- planned_edge_contrasts |>
  dplyr::mutate(
    contrast_id = sprintf("TE%02d", dplyr::row_number()),
    network = "Temporal",
    metric = "Mean edge-family contrast",
    outcome_label = unname(node_label_lookup[.data$outcome]),
    outcome_short = unname(node_short_lookup[.data$outcome]),
    contrast_direction = .data$contrast,
    original_difference = .data$estimate
  ) |>
  dplyr::select(
    "contrast_id", "network", "metric", "family", "outcome",
    "outcome_label", "outcome_short", "contrast_direction",
    "positive_group", "negative_group", "original_difference"
  )

temporal_edge_bootstrap <- matrix(
  NA_real_, nrow = nrow(edge_bootstrap),
  ncol = nrow(temporal_edge_specification),
  dimnames = list(rownames(edge_bootstrap), temporal_edge_specification$contrast_id)
)
edge_original_reproduction <- numeric(nrow(temporal_edge_specification))
for (index in seq_len(nrow(temporal_edge_specification))) {
  specification_row <- temporal_edge_specification[index, ]
  positive_nodes <- split_nodes(specification_row$positive_group)
  negative_nodes <- split_nodes(specification_row$negative_group)
  positive_edges <- paste(
    "Temporal", positive_nodes, specification_row$outcome, sep = "|"
  )
  negative_edges <- paste(
    "Temporal", negative_nodes, specification_row$outcome, sep = "|"
  )
  if (!all(c(positive_edges, negative_edges) %in% colnames(edge_bootstrap))) {
    stop("Required Temporal edge columns are missing for ", specification_row$contrast_id)
  }
  temporal_edge_bootstrap[, specification_row$contrast_id] <-
    rowMeans(edge_bootstrap[, positive_edges, drop = FALSE]) -
    rowMeans(edge_bootstrap[, negative_edges, drop = FALSE])
  edge_original_reproduction[[index]] <-
    mean(cluster_checkpoint$baseline_edges[positive_edges]) -
    mean(cluster_checkpoint$baseline_edges[negative_edges])
}
edge_original_max_abs_difference <- max(abs(
  edge_original_reproduction - temporal_edge_specification$original_difference
))
if (edge_original_max_abs_difference > 1e-12) {
  stop("S6.8 failed to reproduce one or more original S5.3 Temporal contrasts.")
}
temporal_edge_results <- summarize_contrast_matrix(
  temporal_edge_bootstrap, temporal_edge_specification
)

stopifnot(
  nrow(node_target_results) == 24L,
  nrow(core_peripheral_results) == 64L,
  nrow(temporal_edge_results) == 6L,
  identical(dim(node_target_bootstrap), c(1000L, 24L)),
  identical(dim(core_peripheral_bootstrap), c(1000L, 64L)),
  identical(dim(temporal_edge_bootstrap), c(1000L, 6L)),
  all(node_target_results$ci_equals_observed_envelope),
  all(core_peripheral_results$ci_equals_observed_envelope),
  all(temporal_edge_results$ci_equals_observed_envelope)
)

format_number <- function(value) {
  value[abs(value) < .00005] <- 0
  sprintf("%.4f", value)
}
format_estimate_ci <- function(estimate, lower, upper) {
  paste0(
    format_number(estimate), " [", format_number(lower), ", ",
    format_number(upper), "]"
  )
}
yes_no <- function(value) ifelse(value, "Yes", "No")

node_target_sheet <- node_target_results |>
  dplyr::transmute(
    `Contrast ID` = .data$contrast_id,
    Network = .data$network,
    `Bridge EI metric` = .data$metric,
    Node = .data$node_label,
    `Contrast direction` = .data$contrast_direction,
    `Original difference` = .data$original_difference,
    `Bootstrap median` = .data$bootstrap_median,
    `99.9% CI lower` = .data$ci_lower_99_9,
    `99.9% CI upper` = .data$ci_upper_99_9,
    `Estimate [99.9% CI]` = format_estimate_ci(
      .data$original_difference, .data$ci_lower_99_9, .data$ci_upper_99_9
    ),
    `Excludes zero` = yes_no(.data$excludes_zero_99_9),
    `Exploratory decision` = .data$exploratory_decision
  )

core_peripheral_sheet <- core_peripheral_results |>
  dplyr::transmute(
    `Contrast ID` = .data$contrast_id,
    Network = .data$network,
    `Bridge EI metric` = .data$metric,
    `Target community` = .data$target_community,
    `Peripheral node` = .data$peripheral_label,
    `Core node` = .data$core_label,
    `Contrast direction` = .data$contrast_direction,
    `Original difference` = .data$original_difference,
    `Bootstrap median` = .data$bootstrap_median,
    `99.9% CI lower` = .data$ci_lower_99_9,
    `99.9% CI upper` = .data$ci_upper_99_9,
    `Estimate [99.9% CI]` = format_estimate_ci(
      .data$original_difference, .data$ci_lower_99_9, .data$ci_upper_99_9
    ),
    `Excludes zero` = yes_no(.data$excludes_zero_99_9),
    `Exploratory decision` = .data$exploratory_decision
  )

temporal_edge_sheet <- temporal_edge_results |>
  dplyr::transmute(
    `Contrast ID` = .data$contrast_id,
    Network = .data$network,
    `Contrast family` = .data$family,
    Outcome = paste0(.data$outcome_short, ": ", .data$outcome_label),
    `Contrast direction` = .data$contrast_direction,
    `Positive group` = .data$positive_group,
    `Negative group` = .data$negative_group,
    `Original difference` = .data$original_difference,
    `Bootstrap median` = .data$bootstrap_median,
    `99.9% CI lower` = .data$ci_lower_99_9,
    `99.9% CI upper` = .data$ci_upper_99_9,
    `Estimate [99.9% CI]` = format_estimate_ci(
      .data$original_difference, .data$ci_lower_99_9, .data$ci_upper_99_9
    ),
    `Excludes zero` = yes_no(.data$excludes_zero_99_9),
    `Exploratory decision` = .data$exploratory_decision
  )

workbook_sheets <- list(
  Node_Target = as.data.frame(node_target_sheet, stringsAsFactors = FALSE),
  Core_Peripheral = as.data.frame(
    core_peripheral_sheet, stringsAsFactors = FALSE
  ),
  Temporal_Edge_Family = as.data.frame(
    temporal_edge_sheet, stringsAsFactors = FALSE
  )
)

common_note <- paste(
  "Note. All contrasts are paired within the same 1,000 participant-cluster",
  "nonparametric bootstrap replications and use full unthresholded networks.",
  "Node-target contrasts are Distress-target minus Engagement-target Bridge EI;",
  "Core-Peripheral node contrasts are Core node minus Peripheral node; the six",
  "Temporal edge-family contrast directions reproduce S5.3. The exploratory",
  "decision uses a nominal two-sided 99.9% percentile interval (quantile type 6)",
  "and asks whether it excludes zero. With B = 1,000, the 0.05th and 99.95th",
  "percentiles are not directly tail-resolved and reduce to the observed minimum",
  "and maximum. The decision is therefore a conservative all-replication",
  "sign-stability envelope, not an FDR-adjusted result or p value. Temporal",
  "associations are conditional predictive associations and are not causal effects."
)
workbook_notes <- stats::setNames(
  rep(common_note, length(workbook_sheets)), names(workbook_sheets)
)

bruceR::export(workbook_sheets, file = table_file, overwrite = TRUE)

style_workbook <- function(path, sheet_data, notes) {
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
  number_style <- openxlsx::createStyle(numFmt = "0.0000")
  evidence_style <- openxlsx::createStyle(
    fontName = "Aptos", fontSize = 10, fontColour = "#006100",
    fgFill = "#C6EFCE", valign = "center", wrapText = TRUE
  )
  note_style <- openxlsx::createStyle(
    fontName = "Aptos", fontSize = 9, fontColour = "#404040",
    fgFill = "#F2F2F2", valign = "top", wrapText = TRUE,
    border = "top", borderColour = "#A6A6A6"
  )

  for (sheet_name in names(sheet_data)) {
    data <- sheet_data[[sheet_name]]
    last_data_row <- nrow(data) + 1L
    last_column <- ncol(data)
    note_row <- last_data_row + 2L
    numeric_columns <- which(names(data) %in% c(
      "Original difference", "Bootstrap median",
      "99.9% CI lower", "99.9% CI upper"
    ))
    evidence_rows <- which(data[["Excludes zero"]] == "Yes") + 1L

    openxlsx::showGridLines(workbook, sheet_name, showGridLines = FALSE)
    openxlsx::freezePane(workbook, sheet_name, firstRow = TRUE)
    openxlsx::pageSetup(
      workbook, sheet_name, orientation = "landscape", paperSize = 9,
      fitToWidth = 1, fitToHeight = 0, printTitleRows = 1L
    )
    openxlsx::addFilter(
      workbook, sheet_name, row = 1L, cols = seq_len(last_column)
    )
    openxlsx::setRowHeights(workbook, sheet_name, rows = 1L, heights = 34)
    openxlsx::addStyle(
      workbook, sheet_name, header_style, rows = 1L,
      cols = seq_len(last_column), gridExpand = TRUE, stack = TRUE
    )
    openxlsx::addStyle(
      workbook, sheet_name, body_style, rows = 2L:last_data_row,
      cols = seq_len(last_column), gridExpand = TRUE, stack = TRUE
    )
    openxlsx::addStyle(
      workbook, sheet_name, number_style, rows = 2L:last_data_row,
      cols = numeric_columns, gridExpand = TRUE, stack = TRUE
    )
    if (length(evidence_rows)) {
      openxlsx::addStyle(
        workbook, sheet_name, evidence_style, rows = evidence_rows,
        cols = which(names(data) %in% c(
          "Excludes zero", "Exploratory decision"
        )), gridExpand = TRUE, stack = TRUE
      )
    }
    openxlsx::writeData(
      workbook, sheet_name, notes[[sheet_name]],
      startRow = note_row, startCol = 1L, colNames = FALSE
    )
    openxlsx::mergeCells(
      workbook, sheet_name, cols = seq_len(last_column), rows = note_row
    )
    openxlsx::addStyle(
      workbook, sheet_name, note_style, rows = note_row,
      cols = seq_len(last_column), gridExpand = TRUE, stack = TRUE
    )
    openxlsx::setRowHeights(workbook, sheet_name, rows = note_row, heights = 72)

    widths <- rep(14, last_column)
    names(widths) <- names(data)
    widths[intersect(names(widths), c("Contrast ID", "Node"))] <- 11
    widths[intersect(names(widths), c("Network", "Bridge EI metric"))] <- 18
    widths[intersect(names(widths), c(
      "Target community", "Peripheral node", "Core node"
    ))] <- 18
    widths[intersect(names(widths), c(
      "Contrast direction", "Positive group", "Negative group"
    ))] <- 34
    widths[intersect(names(widths), "Outcome")] <- 22
    widths[intersect(names(widths), "Estimate [99.9% CI]")] <- 28
    widths[intersect(names(widths), "Exploratory decision")] <- 42
    openxlsx::setColWidths(
      workbook, sheet_name, cols = seq_len(last_column), widths = widths
    )
  }

  # Remove orphan drawing relationships retained by bruceR/openxlsx in a
  # data-only workbook so strict OOXML readers can import all three sheets.
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
style_workbook(table_file, workbook_sheets, workbook_notes)

expected_sheet_names <- names(workbook_sheets)
reloaded_sheet_names <- readxl::excel_sheets(table_file)
reloaded_sheets <- stats::setNames(lapply(expected_sheet_names, function(sheet) {
  expected_data <- workbook_sheets[[sheet]]
  as.data.frame(readxl::read_excel(
    table_file,
    sheet = sheet,
    range = readxl::cell_limits(
      c(1L, 1L), c(nrow(expected_data) + 1L, ncol(expected_data))
    )
  ), stringsAsFactors = FALSE)
}), expected_sheet_names)

formula_error_pattern <- "^#(REF!|VALUE!|NAME\\?|DIV/0!|N/A)"
workbook_audit <- data.frame(
  check = c(
    "Exact sheet names and order",
    "Exact 24/64/6 contrast row counts",
    "Expected column counts",
    "No missing table cells",
    "No Excel formula errors",
    "Nonempty workbook"
  ),
  passed = c(
    identical(reloaded_sheet_names, expected_sheet_names),
    identical(
      unname(vapply(reloaded_sheets, nrow, integer(1))),
      c(24L, 64L, 6L)
    ),
    identical(
      unname(vapply(reloaded_sheets, ncol, integer(1))),
      unname(vapply(workbook_sheets, ncol, integer(1)))
    ),
    all(vapply(reloaded_sheets, function(data) !anyNA(data), logical(1))),
    !any(vapply(reloaded_sheets, function(data) {
      any(grepl(formula_error_pattern, as.matrix(data)))
    }, logical(1))),
    file.exists(table_file) && file.info(table_file)$size > 0
  ),
  stringsAsFactors = FALSE
)

validation <- data.frame(
  check = c(
    "S6.2 full-unthresholded bootstrap input matches current S5.3/S5.6 files",
    "Node-target contrast count and bootstrap dimensions",
    "Core-Peripheral contrast count and bootstrap dimensions",
    "Temporal edge-family contrast count and bootstrap dimensions",
    "Original S5.3 Temporal contrasts exactly reproduced",
    "All contrast estimates and intervals finite",
    "99.9% type-6 endpoints equal the B=1000 observed envelopes",
    "Workbook reload and structural audit"
  ),
  passed = c(
    identical(current_s5_3_md5, s6_2$provenance$network_md5) &&
      identical(current_cluster_md5, s6_2$provenance$cluster_checkpoint_md5),
    nrow(node_target_results) == 24L &&
      identical(dim(node_target_bootstrap), c(1000L, 24L)),
    nrow(core_peripheral_results) == 64L &&
      identical(dim(core_peripheral_bootstrap), c(1000L, 64L)),
    nrow(temporal_edge_results) == 6L &&
      identical(dim(temporal_edge_bootstrap), c(1000L, 6L)),
    edge_original_max_abs_difference <= 1e-12,
    all(is.finite(c(
      node_target_results$original_difference,
      node_target_results$ci_lower_99_9,
      node_target_results$ci_upper_99_9,
      core_peripheral_results$original_difference,
      core_peripheral_results$ci_lower_99_9,
      core_peripheral_results$ci_upper_99_9,
      temporal_edge_results$original_difference,
      temporal_edge_results$ci_lower_99_9,
      temporal_edge_results$ci_upper_99_9
    ))),
    all(c(
      node_target_results$ci_equals_observed_envelope,
      core_peripheral_results$ci_equals_observed_envelope,
      temporal_edge_results$ci_equals_observed_envelope
    )),
    all(workbook_audit$passed)
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) {
  print(validation, row.names = FALSE)
  stop("One or more S6.8 validation checks failed.")
}

evidence_counts <- data.frame(
  contrast_family = c(
    "Node-target", "Core-Peripheral node", "Temporal edge-family"
  ),
  total = c(24L, 64L, 6L),
  excludes_zero_99_9 = c(
    sum(node_target_results$excludes_zero_99_9),
    sum(core_peripheral_results$excludes_zero_99_9),
    sum(temporal_edge_results$excludes_zero_99_9)
  ),
  stringsAsFactors = FALSE
)

results_s6_8_apriori_bootstrap_difference_tests <- list(
  schema_version = "1.0.0",
  created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  analysis_rule = list(
    framework = "Engagement-Harm Differentiation",
    resampling = paste(
      "Paired values from the same 1,000 participant-cluster nonparametric",
      "bootstrap replications"
    ),
    network_basis = "Full unthresholded original and bootstrap networks",
    node_target_direction = "Distress-target Bridge EI minus Engagement-target Bridge EI",
    core_peripheral_direction = "Core-node Bridge EI minus Peripheral-node Bridge EI",
    temporal_edge_direction = paste(
      "S5.3 directions: Core minus Peripheral for Distress outcomes;",
      "Peripheral minus Core for Engagement outcomes"
    ),
    confidence_level = confidence_level,
    interval_probabilities = interval_probabilities,
    quantile_type = quantile_type,
    decision = "The nominal 99.9% percentile interval excludes zero",
    p_values = "Not calculated",
    fdr = "Not applied"
  ),
  bootstrap_resolution = list(
    replications = bootstrap_replications,
    requested_confidence = confidence_level,
    maximum_directly_resolved_two_sided_confidence =
      maximum_directly_resolved_confidence,
    tail_resolution_sufficient = tail_resolution_sufficient,
    type6_endpoint_behavior = paste(
      "At B=1000, the nominal 0.05th and 99.95th percentile endpoints",
      "equal the observed minimum and maximum; exclusion of zero therefore",
      "requires all 1,000 bootstrap contrasts to have the same sign."
    )
  ),
  node_target = list(
    results = node_target_results,
    bootstrap_contrasts = node_target_bootstrap
  ),
  core_peripheral = list(
    results = core_peripheral_results,
    bootstrap_contrasts = core_peripheral_bootstrap
  ),
  temporal_edge_family = list(
    results = temporal_edge_results,
    bootstrap_contrasts = temporal_edge_bootstrap,
    original_reproduction_max_abs_difference = edge_original_max_abs_difference
  ),
  evidence_counts = evidence_counts,
  workbook_sheets = workbook_sheets,
  workbook_notes = workbook_notes,
  workbook_audit = workbook_audit,
  validation = validation,
  provenance = list(
    s6_2_input = paste0("Res_2_Data_for_Analysis/", basename(s6_2_file)),
    s6_2_md5 = unname(tools::md5sum(s6_2_file)),
    s5_3_input = paste0("Res_2_Data_for_Analysis/", basename(s5_3_file)),
    s5_3_md5 = current_s5_3_md5,
    cluster_checkpoint_input = paste0(
      "Res_2_Data_for_Analysis/", basename(cluster_checkpoint_file)
    ),
    cluster_checkpoint_md5 = current_cluster_md5,
    workbook_output = paste0("Res_3_Results/", basename(table_file)),
    workbook_md5 = unname(tools::md5sum(table_file)),
    generating_script = paste0(
      "Step_6_mlVAR_Reproduce_Centrality/", basename(script_file)
    ),
    script_md5 = unname(tools::md5sum(script_file)),
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
  results_s6_8_apriori_bootstrap_difference_tests,
  file = result_file,
  version = 2
)

reload_environment <- new.env(parent = emptyenv())
reloaded_names <- load(result_file, envir = reload_environment)
if (!identical(
  reloaded_names, "results_s6_8_apriori_bootstrap_difference_tests"
)) {
  stop("Saved S6.8 RDA does not contain its sole expected object.")
}
reloaded_result <- reload_environment[[reloaded_names]]
stopifnot(
  all(reloaded_result$validation$passed),
  all(reloaded_result$workbook_audit$passed),
  identical(
    reloaded_result$provenance$workbook_md5,
    unname(tools::md5sum(table_file))
  )
)

cat("S6.8 completed successfully.\n")
cat("Saved workbook: ", table_file, "\n", sep = "")
cat("Saved result object: ", result_file, "\n", sep = "")
cat(
  "Nominal confidence: 99.9%; B=1000 directly resolves at most 99.8%; ",
  "type-6 endpoints equal min/max.\n",
  sep = ""
)
print(evidence_counts, row.names = FALSE)
print(validation, row.names = FALSE)
