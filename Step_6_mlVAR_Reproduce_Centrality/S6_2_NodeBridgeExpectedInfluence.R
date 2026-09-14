# S6.2: Target-specific raw one-step Bridge Expected Influence
# Author: Kunru Song

options(stringsAsFactors = FALSE, scipen = 999)
try(Sys.setlocale("LC_CTYPE", ".UTF-8"), silent = TRUE)

required_packages <- c(
  "dplyr", "ggplot2", "networktools", "patchwork", "scales", "svglite"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required R package(s): ", paste(missing_packages, collapse = ", "),
    ". Install them under R 4.5.2 before rerunning S6.2."
  )
}

# Locate this file under both command-line Rscript execution and interactive
# source(). This avoids a machine-specific absolute project-root path.
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
    if (is.null(frame_file) || length(frame_file) == 0L) {
      return(NA_character_)
    }
    as.character(frame_file[[1L]])
  }, character(1))
  candidates <- c(candidates, frame_files[!is.na(frame_files)])

  candidates <- unique(candidates[nzchar(candidates)])
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) > 0L) {
    candidates <- normalizePath(candidates, winslash = "/", mustWork = TRUE)
    candidates <- candidates[basename(candidates) == expected_basename]
  }
  # Use the working directory only when neither invocation mechanism identifies
  # this script; a same-named file there must not override the sourced file.
  if (length(candidates) == 0L && file.exists(expected_basename)) {
    candidates <- normalizePath(expected_basename, winslash = "/", mustWork = TRUE)
  }
  if (length(candidates) == 0L) {
    stop(
      "Cannot locate ", expected_basename, ". Run it with Rscript or call ",
      "source() with the script path."
    )
  }
  candidates[[length(candidates)]]
}

script_file <- resolve_script_file("S6_2_NodeBridgeExpectedInfluence.R")
script_dir <- dirname(script_file)
setwd(script_dir)
message("S6.2 working directory: ", normalizePath(".", winslash = "/"))

analysis_dir <- "../Res_2_Data_for_Analysis"
results_dir <- "../Res_3_Results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

network_file <- file.path(analysis_dir, "results_s5_3_mlvar_networks.rda")
cluster_checkpoint_file <- file.path(
  analysis_dir, "checkpoint_s5_6_cluster_bootstrap.rda"
)
result_file <- file.path(
  analysis_dir, "results_s6_2_node_bridge_expected_influence.rda"
)
figure_paths <- c(
  temporal_svg = file.path(
    results_dir, "Figure_S6_2_NodeBridgeExpectedInfluence_Temporal.svg"
  ),
  undirected_svg = file.path(
    results_dir, "Figure_S6_2_NodeBridgeExpectedInfluence_Undirected.svg"
  )
)

load_single_object <- function(path, expected_name) {
  if (!file.exists(path)) {
    stop("Required input does not exist: ", path)
  }
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
cluster_checkpoint <- load_single_object(
  cluster_checkpoint_file, "checkpoint_s5_6_cluster_bootstrap"
)

required_network_components <- c("full_matrices", "node_map", "provenance")
if (!all(required_network_components %in% names(network_results))) {
  stop("The S5.3 network object is missing required components.")
}
required_checkpoint_components <- c(
  "type", "replications", "network_md5", "baseline_edges", "edge_estimates"
)
if (!all(required_checkpoint_components %in% names(cluster_checkpoint))) {
  stop("The S5.6 cluster-bootstrap checkpoint is missing required components.")
}

node_map <- network_results$node_map |>
  dplyr::arrange(.data$order)
node_names <- node_map$node
communities <- stats::setNames(node_map$community, node_names)
psmu_nodes <- node_map$node[node_map$community == "PSMU"]
engagement_nodes <- node_map$node[node_map$community == "Engagement Time"]
distress_nodes <- node_map$node[node_map$community == "Distress"]
psmu_display <- stats::setNames(paste0("P", seq_along(psmu_nodes)), psmu_nodes)

stopifnot(
  nrow(node_map) == 12L,
  identical(node_names, rownames(network_results$full_matrices$temporal)),
  identical(node_names, colnames(network_results$full_matrices$temporal)),
  identical(length(psmu_nodes), 6L),
  identical(length(engagement_nodes), 3L),
  identical(length(distress_nodes), 3L),
  identical(unname(communities[psmu_nodes]), rep("PSMU", 6L))
)

full_matrices <- lapply(
  network_results$full_matrices,
  function(matrix_object) matrix_object[node_names, node_names, drop = FALSE]
)
if (!identical(names(full_matrices), c("temporal", "contemporaneous", "between"))) {
  stop("Unexpected full-matrix names in the S5.3 object.")
}
for (network_name in names(full_matrices)) {
  matrix_object <- full_matrices[[network_name]]
  if (!is.matrix(matrix_object) || !identical(dim(matrix_object), c(12L, 12L))) {
    stop(network_name, " full matrix must be 12 x 12.")
  }
  if (anyNA(matrix_object) || any(!is.finite(matrix_object))) {
    stop(network_name, " full matrix contains missing or non-finite values.")
  }
  if (network_name != "temporal") {
    if (!isTRUE(all.equal(matrix_object, t(matrix_object), tolerance = 1e-12))) {
      stop(network_name, " full matrix must be symmetric.")
    }
    if (max(abs(diag(matrix_object))) > 1e-12) {
      stop(network_name, " full matrix must have a structural-zero diagonal.")
    }
  }
}

network_md5 <- unname(tools::md5sum(network_file))
cluster_checkpoint_md5 <- unname(tools::md5sum(cluster_checkpoint_file))
checkpoint_network_md5 <- unname(as.character(cluster_checkpoint$network_md5))
network_file_md5_matches_checkpoint <- identical(checkpoint_network_md5, network_md5)

edge_estimates <- cluster_checkpoint$edge_estimates
if (!is.matrix(edge_estimates) || !identical(dim(edge_estimates), c(1000L, 276L))) {
  stop("The cluster-bootstrap edge matrix must be exactly 1000 x 276.")
}
if (is.null(colnames(edge_estimates)) || anyDuplicated(colnames(edge_estimates))) {
  stop("The cluster-bootstrap edge matrix needs unique edge names.")
}
if (anyNA(edge_estimates) || any(!is.finite(edge_estimates))) {
  stop("The cluster-bootstrap edge matrix contains missing or non-finite values.")
}

edge_vectors_from_matrices <- function(matrices) {
  temporal_indices <- expand.grid(
    from = node_names,
    to = node_names,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  temporal_values <- matrices$temporal[cbind(
    match(temporal_indices$from, node_names),
    match(temporal_indices$to, node_names)
  )]
  names(temporal_values) <- paste(
    "Temporal", temporal_indices$from, temporal_indices$to, sep = "|"
  )

  undirected_values <- function(matrix_object, network_label) {
    pair_indices <- which(upper.tri(matrix_object), arr.ind = TRUE)
    values <- matrix_object[pair_indices]
    names(values) <- paste(
      network_label,
      rownames(matrix_object)[pair_indices[, 1L]],
      colnames(matrix_object)[pair_indices[, 2L]],
      sep = "|"
    )
    values
  }

  c(
    temporal_values,
    undirected_values(matrices$contemporaneous, "Contemporaneous"),
    undirected_values(matrices$between, "Between-person")
  )
}

baseline_edges_from_s5_3 <- edge_vectors_from_matrices(full_matrices)
if (!identical(names(baseline_edges_from_s5_3), colnames(edge_estimates))) {
  stop("S5.3 baseline-edge names do not match the S5.6 bootstrap edge columns.")
}
if (!identical(names(cluster_checkpoint$baseline_edges), colnames(edge_estimates))) {
  stop("S5.6 stored baseline-edge names do not match its bootstrap edge columns.")
}
baseline_edge_max_abs_difference <- max(abs(
  baseline_edges_from_s5_3 - cluster_checkpoint$baseline_edges
))
if (baseline_edge_max_abs_difference > 1e-12) {
  stop("S5.3 and S5.6 baseline full edges are not numerically identical.")
}

panel_specification <- data.frame(
  panel = LETTERS[1:4],
  panel_id = c("temporal_in", "temporal_out", "contemporaneous", "between"),
  network = c("Temporal", "Temporal", "Contemporaneous", "Between-person"),
  metric = c(
    "Bridge In-EI", "Bridge Out-EI", "Bridge EI", "Bridge EI"
  ),
  title = c(
    "Bridge In-EI",
    "Bridge Out-EI",
    "Contemporaneous",
    "Between-person"
  ),
  directed = c(TRUE, TRUE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

target_specification <- data.frame(
  target_id = c("engagement", "distress"),
  target_community = c("Engagement Time", "Distress"),
  target_label = c("Engagement Community", "Distress Community"),
  target_color = c("#D99A28", "#319B55"),
  stringsAsFactors = FALSE
)

estimand_specification <- expand.grid(
  panel_id = panel_specification$panel_id,
  target_id = target_specification$target_id,
  node = psmu_nodes,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) |>
  dplyr::left_join(panel_specification, by = "panel_id") |>
  dplyr::left_join(target_specification, by = "target_id") |>
  dplyr::mutate(
    node_display = unname(psmu_display[.data$node]),
    key = paste(.data$panel_id, .data$target_id, .data$node, sep = "|")
  ) |>
  dplyr::arrange(
    match(.data$panel_id, panel_specification$panel_id),
    match(.data$node, psmu_nodes),
    match(.data$target_id, target_specification$target_id)
  )

target_nodes_by_id <- list(
  engagement = engagement_nodes,
  distress = distress_nodes
)

raw_bridge_vector_networktools <- function(matrices) {
  output <- stats::setNames(rep(NA_real_, nrow(estimand_specification)), estimand_specification$key)

  for (panel_id in panel_specification$panel_id) {
    matrix_object <- switch(
      panel_id,
      temporal_in = t(matrices$temporal),
      temporal_out = matrices$temporal,
      contemporaneous = matrices$contemporaneous,
      between = matrices$between,
      stop("Unknown panel_id: ", panel_id)
    )
    directed <- panel_specification$directed[
      match(panel_id, panel_specification$panel_id)
    ]

    for (target_id in target_specification$target_id) {
      target_community <- target_specification$target_community[
        match(target_id, target_specification$target_id)
      ]
      bridge_result <- suppressMessages(networktools::bridge(
        network = matrix_object,
        communities = unname(communities),
        useCommunities = c("PSMU", target_community),
        directed = directed,
        normalize = FALSE
      ))
      values <- bridge_result[["Bridge Expected Influence (1-step)"]][psmu_nodes]
      keys <- paste(panel_id, target_id, psmu_nodes, sep = "|")
      output[keys] <- unname(values)
    }
  }

  if (anyNA(output) || any(!is.finite(output))) {
    stop("networktools returned missing or non-finite S6.2 values.")
  }
  output[estimand_specification$key]
}

raw_bridge_vector_manual <- function(matrices) {
  output <- stats::setNames(rep(NA_real_, nrow(estimand_specification)), estimand_specification$key)

  for (row_index in seq_len(nrow(estimand_specification))) {
    specification_row <- estimand_specification[row_index, ]
    target_nodes <- target_nodes_by_id[[specification_row$target_id]]
    node <- specification_row$node
    value <- switch(
      specification_row$panel_id,
      temporal_in = sum(matrices$temporal[target_nodes, node]),
      temporal_out = sum(matrices$temporal[node, target_nodes]),
      contemporaneous = sum(matrices$contemporaneous[node, target_nodes]),
      between = sum(matrices$between[node, target_nodes]),
      stop("Unknown panel_id: ", specification_row$panel_id)
    )
    output[[specification_row$key]] <- value
  }
  output[estimand_specification$key]
}

original_networktools <- raw_bridge_vector_networktools(full_matrices)
original_manual <- raw_bridge_vector_manual(full_matrices)
original_networktools_manual_max_abs_difference <- max(abs(
  original_networktools - original_manual
))
if (original_networktools_manual_max_abs_difference > 1e-12) {
  stop("Original-sample networktools values do not equal manual raw signed sums.")
}

undirected_edge_name <- function(network_label, node_a, node_b) {
  index_a <- match(node_a, node_names)
  index_b <- match(node_b, node_names)
  if (is.na(index_a) || is.na(index_b) || index_a == index_b) {
    stop("Invalid undirected edge request: ", node_a, " / ", node_b)
  }
  if (index_a < index_b) {
    paste(network_label, node_a, node_b, sep = "|")
  } else {
    paste(network_label, node_b, node_a, sep = "|")
  }
}

bootstrap_raw_bridge_matrix <- matrix(
  NA_real_,
  nrow = nrow(edge_estimates),
  ncol = nrow(estimand_specification),
  dimnames = list(rownames(edge_estimates), estimand_specification$key)
)

for (row_index in seq_len(nrow(estimand_specification))) {
  specification_row <- estimand_specification[row_index, ]
  node <- specification_row$node
  target_nodes <- target_nodes_by_id[[specification_row$target_id]]

  edge_names <- switch(
    specification_row$panel_id,
    temporal_in = paste("Temporal", target_nodes, node, sep = "|"),
    temporal_out = paste("Temporal", node, target_nodes, sep = "|"),
    contemporaneous = vapply(
      target_nodes,
      function(target_node) {
        undirected_edge_name("Contemporaneous", node, target_node)
      },
      character(1)
    ),
    between = vapply(
      target_nodes,
      function(target_node) {
        undirected_edge_name("Between-person", node, target_node)
      },
      character(1)
    ),
    stop("Unknown panel_id: ", specification_row$panel_id)
  )

  if (!all(edge_names %in% colnames(edge_estimates))) {
    stop("One or more required bootstrap edge columns are missing.")
  }
  bootstrap_raw_bridge_matrix[, specification_row$key] <- rowSums(
    edge_estimates[, edge_names, drop = FALSE]
  )
}

if (anyNA(bootstrap_raw_bridge_matrix) ||
    any(!is.finite(bootstrap_raw_bridge_matrix))) {
  stop("The derived 1000 x 48 bootstrap Bridge EI matrix is not fully finite.")
}

edge_vector_to_matrices <- function(edge_vector) {
  if (!identical(names(edge_vector), colnames(edge_estimates))) {
    stop("Replication edge vector does not follow the canonical edge order.")
  }
  parsed_names <- strsplit(names(edge_vector), "|", fixed = TRUE)
  parsed_matrix <- do.call(rbind, parsed_names)

  temporal <- matrix(
    0, nrow = 12L, ncol = 12L,
    dimnames = list(node_names, node_names)
  )
  contemporaneous <- temporal
  between <- temporal

  temporal_rows <- parsed_matrix[, 1L] == "Temporal"
  temporal[cbind(
    match(parsed_matrix[temporal_rows, 2L], node_names),
    match(parsed_matrix[temporal_rows, 3L], node_names)
  )] <- edge_vector[temporal_rows]

  fill_undirected <- function(network_label, output_matrix) {
    selected <- parsed_matrix[, 1L] == network_label
    indices <- cbind(
      match(parsed_matrix[selected, 2L], node_names),
      match(parsed_matrix[selected, 3L], node_names)
    )
    output_matrix[indices] <- edge_vector[selected]
    output_matrix[cbind(indices[, 2L], indices[, 1L])] <- edge_vector[selected]
    diag(output_matrix) <- 0
    output_matrix
  }

  list(
    temporal = temporal,
    contemporaneous = fill_undirected("Contemporaneous", contemporaneous),
    between = fill_undirected("Between-person", between)
  )
}

# A deterministic spread of replications is recalculated with networktools.
# Increase audit_replications below (for example to 25 or 100) for a slower,
# denser implementation audit; this does not change the 1,000-replication CIs.
audit_replications <- unique(as.integer(round(seq(1, 1000, length.out = 12L))))
bootstrap_networktools_audit_max_abs_difference <- 0
for (replication in audit_replications) {
  replication_matrices <- edge_vector_to_matrices(edge_estimates[replication, ])
  networktools_values <- raw_bridge_vector_networktools(replication_matrices)
  stored_values <- bootstrap_raw_bridge_matrix[replication, ]
  bootstrap_networktools_audit_max_abs_difference <- max(
    bootstrap_networktools_audit_max_abs_difference,
    abs(networktools_values - stored_values)
  )
}
if (bootstrap_networktools_audit_max_abs_difference > 1e-12) {
  stop("Bootstrap batch sums do not equal the audited networktools values.")
}

bootstrap_intervals <- t(apply(
  bootstrap_raw_bridge_matrix,
  2L,
  stats::quantile,
  probs = c(.025, .975),
  na.rm = FALSE,
  names = FALSE
))
colnames(bootstrap_intervals) <- c("ci_lower", "ci_upper")

point_estimates <- estimand_specification |>
  dplyr::mutate(
    estimate = unname(original_networktools[.data$key]),
    ci_lower = bootstrap_intervals[.data$key, "ci_lower"],
    ci_upper = bootstrap_intervals[.data$key, "ci_upper"],
    includes_zero = .data$ci_lower <= 0 & .data$ci_upper >= 0,
    ci_exclusion_marker = dplyr::if_else(
      !.data$includes_zero,
      "*",
      NA_character_
    ),
    node_display = factor(.data$node_display, levels = rev(unname(psmu_display))),
    target_label = factor(
      .data$target_label,
      levels = target_specification$target_label
    )
  )

if (any(point_estimates$ci_lower > point_estimates$ci_upper)) {
  stop("Internal interval-order validation failed.")
}
if (any(!point_estimates$includes_zero & point_estimates$estimate == 0)) {
  stop("A CI-excluding estimate is exactly zero, so its annotation side is undefined.")
}

# -----------------------------------------------------------------------------
# EDITABLE FIGURE AESTHETICS
# -----------------------------------------------------------------------------
# Every object in this block may be edited without changing the statistical
# results. Suggested alternatives are given in each comment.

engagement_color <- "#D99A28" # EDIT color: use any hex code or named R color, e.g., "#E69F00", "orange3", or "#CC79A7".
distress_color <- "#319B55" # EDIT color: use any hex code or named R color, e.g., "#009E73", "forestgreen", or "#0072B2".
zero_line_color <- "#4D4D4D" # EDIT zero line: alternatives include "black", "grey40", or a hex code.
major_grid_color <- "#E6E6E6" # EDIT grid color: alternatives include "grey90", "#F0F0F0", or NA to remove it in theme().
node_separator_color <- "#E1E1E1" # EDIT node separators: alternatives include "grey88", "#ECECEC", or "#D9D9D9".
node_separator_linewidth <- 0.38 # EDIT separator width: useful range 0.2-0.7; larger values make row divisions more prominent.
node_separator_linetype <- "solid" # EDIT separator type: choose "solid", "dashed", "dotted", "dotdash", or numeric 1-6.
panel_border_color <- "black" # EDIT frame: alternatives include "#222222", "grey30", or another named/hex color.
bar_outline_linewidth <- 0.75 # EDIT bar outline width: useful range 0.4-1.4; larger values produce heavier outlines.
bar_width <- 0.40 # EDIT bar thickness: choose 0-1; typical values are 0.4, 0.5, 0.6, or 0.7.
dodge_width <- 0.64 # EDIT separation: typical values 0.55-0.85; keep it at least as large as bar_width.
errorbar_linewidth <- 1 # EDIT CI line width: useful range 0.35-1.2.
errorbar_cap_height <- 0.3 # EDIT CI cap height: use 0 for no caps, or about 0.08-0.25 for visible caps.
endpoint_point_size <- 3 # EDIT endpoint size: typical values 1.5-3.5.
endpoint_point_shape <- 16 # EDIT point shape: 16 or 19 = solid circle; 15 = square; 17 = triangle; 21 = bordered circle.
ci_exclusion_star_size <- 6 # EDIT asterisk size: typical values are 2.5-5; larger values make CI-exclusion markers more prominent.
ci_exclusion_star_offset_fraction <- 0.06 # EDIT outward offset as a fraction of each panel's CI span; typical values are 0.01-0.04.
zero_line_linetype <- "solid" # EDIT linetype: choose "solid", "dashed", "dotted", "dotdash", or numeric 1-6.
zero_line_linewidth <- 0.5 # EDIT zero-line width: useful range 0.25-0.9.
base_font_size <- 9 # EDIT global font size: typical manuscript values are 8-12 points.
axis_text_size <- 10 # EDIT tick-label size: typical values are 7-11 points.
axis_title_size <- 5 # EDIT x-axis-title size: typical values are 8-12 points.
panel_title_size <- 10 # EDIT panel-title size: typical values are 9-13 points.
legend_text_size <- 8.2 # EDIT legend-text size: typical values are 7-11 points.
panel_tag_size <- 11 # EDIT A-D tag size: typical values are 9-14 points.
panel_spacing_lines <- 0.9 # EDIT panel spacing: typical values are 0.5-1.5 text lines.
x_axis_expansion <- c(.06, .10) # EDIT x margins: use two nonnegative multipliers, e.g., c(.03,.06), c(.05,.10), or c(.10,.10).

# EDIT THE TEMPORAL FIGURE WIDTH AND HEIGHT HERE (inches).
# Width alternatives: 7-8.27 for manuscript/A4 width, or 10-12 for a wide supplement.
# Height alternatives: 3.5-5.5 for two side-by-side six-row panels.
temporal_figure_width <- 6
temporal_figure_height <- 5
# EDIT THE UNDIRECTED FIGURE WIDTH AND HEIGHT HERE (inches).
# Width alternatives: 7-8.27 for manuscript/A4 width, or 10-12 for a wide supplement.
# Height alternatives: 3.5-5.5 for two side-by-side six-row panels.
undirected_figure_width <- 6
undirected_figure_height <- 5

plot_colors <- c(
  "Engagement Community" = engagement_color,
  "Distress Community" = distress_color
)

dodge_position <- ggplot2::position_dodge(width = dodge_width) # EDIT `width`: optional values are usually 0.6-0.9; larger values separate paired bars more.

base_panel_theme <- ggplot2::theme_minimal(base_size = base_font_size) + # EDIT base theme: alternatives are theme_classic(), theme_bw(), theme_light(), or theme_minimal().
  ggplot2::theme( # EDIT any element below; use element_blank() to remove it or element_text()/element_line() to show it.
    panel.grid.minor = ggplot2::element_blank(), # EDIT minor grid: replace with element_line(color="grey95") to display it.
    panel.grid.major.y = ggplot2::element_blank(), # LOCKED off because explicit geom_hline separators are drawn between nodes; default major-y lines would cross bar centers.
    panel.grid.major.x = ggplot2::element_line(color = major_grid_color, linewidth = .32), # EDIT vertical grid color/width; use element_blank() to remove.
    panel.border = ggplot2::element_rect(color = panel_border_color, fill = NA, linewidth = .65), # EDIT border color/width; keep fill=NA so data remain visible.
    axis.text = ggplot2::element_text(size = axis_text_size, color = "black"), # EDIT size/color; optional colors include "black", "grey20", or hex values.
    axis.title.x = ggplot2::element_text(size = axis_title_size, color = "black", margin = ggplot2::margin(t = 5)), # EDIT size/color/top margin; typical margins are 3-8 points.
    axis.title.y = ggplot2::element_blank(), # EDIT y title: replace with element_text(size=axis_title_size) if a y-axis title is wanted.
    plot.title = ggplot2::element_text(size = panel_title_size, face = "bold", hjust = .5, color = "black"), # EDIT face: "plain", "bold", "italic"; hjust: 0 left, .5 center, 1 right.
    legend.position = "bottom", # EDIT legend position: choose "bottom", "top", "left", "right", "inside", or "none".
    legend.direction = "horizontal", # EDIT legend direction: choose "horizontal" or "vertical".
    legend.text = ggplot2::element_text(size = legend_text_size, color = "black"), # EDIT legend size/color; typical size is 7-11.
    legend.title = ggplot2::element_blank(), # EDIT legend title: replace with element_text() and set labs(colour="Title") if desired.
    plot.margin = ggplot2::margin(5, 7, 5, 6), # EDIT margins in points: order is top, right, bottom, left.
    panel.spacing = grid::unit(panel_spacing_lines, "lines") # EDIT spacing unit: "lines", "pt", "cm", or "in" are valid grid units.
  )

make_panel <- function(panel_id) {
  panel_data <- point_estimates |>
    dplyr::filter(.data$panel_id == .env$panel_id)
  panel_ci_span <- diff(range(c(panel_data$ci_lower, panel_data$ci_upper, 0)))
  if (!is.finite(panel_ci_span) || panel_ci_span <= 0) {
    stop("Cannot calculate a finite positive panel span for CI annotations.")
  }
  panel_star_offset <- ci_exclusion_star_offset_fraction * panel_ci_span
  panel_data <- panel_data |>
    dplyr::mutate(
      ci_exclusion_star_x = dplyr::case_when(
        !.data$includes_zero & .data$estimate < 0 ~
          .data$ci_lower - .env$panel_star_offset,
        !.data$includes_zero & .data$estimate > 0 ~
          .data$ci_upper + .env$panel_star_offset,
        TRUE ~ NA_real_
      )
    )
  panel_title <- panel_specification$title[
    match(panel_id, panel_specification$panel_id)
  ]

  ggplot2::ggplot( # EDIT plotting object only if changing the panel data source or global mappings.
    data = panel_data, # EDIT data: normally keep panel_data; a filtered data frame can be supplied for a subset.
    mapping = ggplot2::aes( # EDIT global aesthetics: optional mappings include alpha, linetype, shape, or fill.
      y = .data$node_display, # LOCKED row mapping P1-P6; reverse factor levels above to reverse vertical order.
      group = .data$target_label # LOCKED dodge grouping; optional grouping variables must still distinguish the two communities.
    )
  ) +
    ggplot2::geom_hline( # EDIT node-row separators: explicit lines lie between P1-P6 rather than crossing their bar centers.
      yintercept = seq(1.5, 5.5, by = 1), # LOCKED separator positions for six discrete rows; edit only if the number/order of nodes changes.
      color = node_separator_color, # EDIT separator color above; named R colors and hex codes are valid.
      linetype = node_separator_linetype, # EDIT separator type above: "solid", "dashed", "dotted", "dotdash", or numeric 1-6.
      linewidth = node_separator_linewidth # EDIT separator width above: useful range 0.2-0.7.
    ) +
    ggplot2::geom_vline( # EDIT zero reference layer: remove this entire layer if no reference line is wanted.
      xintercept = 0, # EDIT reference location: use another numeric value only for a different reference hypothesis.
      color = zero_line_color, # EDIT line color: named R colors or hex codes are valid.
      linetype = zero_line_linetype, # EDIT line type: "solid", "dashed", "dotted", "dotdash", or numeric 1-6.
      linewidth = zero_line_linewidth # EDIT line width: typical values are 0.25-0.9.
    ) +
    ggplot2::geom_col( # EDIT bar layer: geom_col() uses the observed estimate as bar length; geom_bar() is not appropriate here.
      mapping = ggplot2::aes( # EDIT bar aesthetics: fill, alpha, or linetype could be mapped here if the design changes.
        x = .data$estimate, # LOCKED bar length: raw signed one-step Bridge EI.
        colour = .data$target_label # LOCKED community-colored outline, retained even when the 95% CI includes zero.
      ),
      fill = "white", # LOCKED requested white interior; optional alternatives are NA, "transparent", a named color, or a hex code.
      width = bar_width, # EDIT thickness: choose 0-1; typical compact values are 0.35-0.60.
      linewidth = bar_outline_linewidth, # EDIT outline width: typical values are 0.4-1.4.
      position = dodge_position, # EDIT paired-bar spacing through dodge_width above; alternatives include "identity" or position_dodge2().
      show.legend = TRUE # EDIT legend inclusion: TRUE shows color keys; FALSE suppresses this layer from the legend.
    ) +
    ggplot2::geom_errorbar( # EDIT CI layer: horizontal intervals use xmin/xmax plus orientation="y".
      mapping = ggplot2::aes( # EDIT CI aesthetics: alpha or linetype may be added; keep xmin/xmax for percentile limits.
        xmin = .data$ci_lower, # LOCKED lower endpoint: 2.5th percentile of 1,000 full-network bootstrap values.
        xmax = .data$ci_upper, # LOCKED upper endpoint: 97.5th percentile of 1,000 full-network bootstrap values.
        colour = .data$target_label # LOCKED community color for every interval, regardless of whether its 95% CI includes zero.
      ),
      orientation = "y", # LOCKED horizontal orientation; use "x" only for vertical error bars with different mappings.
      width = errorbar_cap_height, # EDIT cap height: 0 removes caps; typical visible values are 0.08-0.25.
      linewidth = errorbar_linewidth, # EDIT line width: typical values are 0.35-1.2.
      position = dodge_position, # LOCKED paired alignment; edit dodge_width above if changing spacing.
      show.legend = FALSE # EDIT legend inclusion: FALSE avoids duplicate legend keys; TRUE would add this layer.
    ) +
    ggplot2::geom_point( # EDIT endpoint layer: this marks the signed terminal of every bar.
      mapping = ggplot2::aes( # EDIT point aesthetics: shape, size, fill, or alpha can be mapped or fixed.
        x = .data$estimate, # LOCKED x position: original-sample raw one-step Bridge EI.
        colour = .data$target_label # LOCKED community color for every point, regardless of whether its 95% CI includes zero.
      ),
      shape = endpoint_point_shape, # EDIT shape: 16/19 circle, 15 square, 17 triangle, 21 bordered circle.
      size = endpoint_point_size, # EDIT size: typical values are 1.5-3.5.
      position = dodge_position, # LOCKED paired alignment; edit dodge_width above if changing spacing.
      show.legend = FALSE # EDIT legend inclusion: FALSE avoids duplicate keys; TRUE would add point keys.
    ) +
    ggplot2::geom_text( # EDIT CI-exclusion annotation layer: remove this layer to suppress all asterisks without changing estimates.
      mapping = ggplot2::aes( # EDIT asterisk aesthetics: fontface, alpha, or family can be added while retaining x/label/color mappings.
        x = .data$ci_exclusion_star_x, # LOCKED side rule: just left of ci_lower for negative estimates and just right of ci_upper for positive estimates.
        label = .data$ci_exclusion_marker, # LOCKED marker rule: "*" only when the unrounded percentile 95% CI excludes zero.
        colour = .data$target_label # LOCKED asterisk color: the corresponding Engagement or Distress community color.
      ),
      size = ci_exclusion_star_size, # EDIT text size above; typical values are 2.5-5.
      hjust = 0.5, # EDIT horizontal alignment: 0 left, 0.5 centered at the offset location, or 1 right.
      vjust = 0.75, # EDIT vertical alignment: 0 bottom, 0.5 centered, or 1 top; values outside 0-1 add further nudging.
      position = dodge_position, # LOCKED paired-row alignment; edit dodge_width above only if all paired geoms should move together.
      na.rm = TRUE, # LOCKED omission: rows whose CI includes zero carry NA and therefore receive no annotation or warning.
      show.legend = FALSE # EDIT legend inclusion: FALSE keeps the shared legend limited to the two communities.
    ) +
    ggplot2::scale_colour_manual( # EDIT discrete palette; scale_color_manual() is the equivalent US-spelling alias.
      values = plot_colors, # EDIT colors in plot_colors above; every named level needs one color.
      breaks = target_specification$target_label, # LOCKED legend to the two community colors shared by bars, points, intervals, and asterisks.
      drop = FALSE # EDIT unused levels: FALSE retains both community keys; TRUE drops absent target communities.
    ) +
    ggplot2::scale_x_continuous( # EDIT x scale; alternatives include scale_x_reverse() or custom limits/breaks.
      breaks = scales::breaks_pretty(n = 5), # EDIT tick count: typical n values are 4-7; explicit numeric breaks are also allowed.
      labels = scales::label_number(accuracy = .01), # EDIT display precision: use .001 for three decimals or .1 for one decimal.
      expand = ggplot2::expansion(mult = x_axis_expansion) # EDIT margins via x_axis_expansion above; add=c(0,.01) is another option.
    ) +
    ggplot2::labs( # EDIT visible text labels; set any item to NULL to remove it.
      title = panel_title, # EDIT panel title in panel_specification above; NULL removes it.
      x = NULL, # EDIT x-axis title; alternatives include "Bridge Expected Influence" or NULL.
      y = NULL, # EDIT y-axis title: use "PSMU node" or another string if desired.
      colour = NULL # EDIT legend title: use "Encoding" or "Community / CI status" if a title is wanted.
    ) +
    base_panel_theme # EDIT shared theme in base_panel_theme above rather than duplicating changes across panels.
}

panel_list <- lapply(panel_specification$panel_id, make_panel)

compose_two_panel_figure <- function(selected_panels) {
  patchwork::wrap_plots( # EDIT layout constructor: plot_spacer() can be inserted for intentional blank cells.
    selected_panels, # EDIT panel order by changing the selected panel indices supplied below.
    ncol = 2L, # EDIT columns: 2 gives the requested side-by-side layout; 1 creates a vertical stack.
    nrow = 1L, # EDIT rows: 1 gives the requested side-by-side layout; 2 works with ncol=1 for a vertical stack.
    byrow = TRUE, # EDIT fill order: TRUE fills left-to-right; FALSE matters only for layouts with multiple rows and columns.
    guides = "collect" # EDIT legend handling: "collect" makes one shared legend; "keep" leaves panel-specific legends.
  ) +
    patchwork::plot_annotation( # EDIT shared annotations; title/subtitle/caption may be added here if needed.
      tag_levels = "A", # EDIT tag sequence within each figure: "A" uppercase, "a" lowercase, "1" numeric, or NULL for no tags.
      tag_suffix = ".", # EDIT tag suffix: ".", ")", ":", or "" are common choices.
      theme = ggplot2::theme( # EDIT the A-B tag appearance and outer margins.
        plot.tag = ggplot2::element_text(face = "bold", size = panel_tag_size, color = "black"), # EDIT tag face/size/color; face options are "plain", "bold", or "italic".
        plot.tag.position = c(.008, .992), # EDIT tag position with x/y values from 0 to 1; c(0,1) is the upper-left corner.
        plot.margin = ggplot2::margin(5, 6, 5, 5) # EDIT outer margins in points: top, right, bottom, left.
      )
    ) &
    ggplot2::theme( # EDIT shared collected-legend placement after patchwork composition.
      legend.position = "bottom", # EDIT shared legend: "bottom", "top", "left", "right", or "none".
      legend.box = "horizontal" # EDIT legend box: "horizontal" or "vertical".
    )
}

temporal_figure <- compose_two_panel_figure(panel_list[1:2])
undirected_figure <- compose_two_panel_figure(panel_list[3:4])

# EDIT WIDTH/HEIGHT at temporal_figure_width and temporal_figure_height in the aesthetic block above.
ggplot2::ggsave( # EDIT temporal SVG export: keep ggsave() so the figure remains reproducible from this script.
  filename = figure_paths[["temporal_svg"]], # EDIT filename only if intentionally changing the canonical temporal output name/path.
  plot = temporal_figure, # LOCKED temporal composition: Bridge In-EI left and Bridge Out-EI right.
  device = svglite::svglite, # EDIT SVG device: svglite is preferred; alternatives include "svg" or grDevices::svg.
  width = temporal_figure_width, # EDIT temporal width in inches in the aesthetic block above.
  height = temporal_figure_height, # EDIT temporal height in inches in the aesthetic block above.
  units = "in", # EDIT units: "in", "cm", or "mm"; width/height must use the same unit.
  bg = "white" # EDIT background: "white", "transparent", or another named/hex color.
)

# EDIT WIDTH/HEIGHT at undirected_figure_width and undirected_figure_height in the aesthetic block above.
ggplot2::ggsave( # EDIT undirected SVG export: keep ggsave() so the figure remains reproducible from this script.
  filename = figure_paths[["undirected_svg"]], # EDIT filename only if intentionally changing the canonical undirected output name/path.
  plot = undirected_figure, # LOCKED undirected composition: Contemporaneous left and Between-person right.
  device = svglite::svglite, # EDIT SVG device: svglite is preferred; alternatives include "svg" or grDevices::svg.
  width = undirected_figure_width, # EDIT undirected width in inches in the aesthetic block above.
  height = undirected_figure_height, # EDIT undirected height in inches in the aesthetic block above.
  units = "in", # EDIT units: "in", "cm", or "mm"; width/height must use the same unit.
  bg = "white" # EDIT background: "white", "transparent", or another named/hex color.
)

if (!all(file.exists(figure_paths)) || any(file.info(figure_paths)$size <= 0)) {
  stop("One or more S6.2 figure outputs are missing or empty.")
}

manuscript_note <- paste(
  paste(
    "Note. Bars show target-specific raw signed one-step Bridge Expected",
    "Influence for P1-P6 from the original-sample full unthresholded networks."
  ),
  paste(
    "Each value is the signed sum of the three edges linking a PSMU node to",
    "the Engagement Time or Distress community; negative and positive edges",
    "can cancel."
  ),
  paste(
    "Points mark original-sample estimates and horizontal error bars show",
    "percentile 95% confidence intervals from 1,000 participant-cluster",
    "nonparametric bootstrap full-network edge vectors."
  ),
  paste(
    "Bar outlines, points, and intervals use their target-community colors.",
    "An asterisk just beyond the left interval endpoint for a negative estimate",
    "or the right endpoint for a positive estimate indicates that the unrounded",
    "95% bootstrap interval excludes zero; intervals including zero are not",
    "annotated. This is not a multiplicity-adjusted significance test or a",
    "direct comparison between communities or nodes."
  ),
  paste(
    "Temporal Bridge In-EI summarizes incoming lagged cross-community",
    "prediction and Temporal Bridge Out-EI summarizes outgoing lagged",
    "cross-community prediction. Undirected Bridge EI is net signed",
    "conditional connectivity, not directional or causal influence."
  )
)

numerical_audit <- data.frame(
  check = c(
    "checkpoint_baseline_reproduces_current_full_edges",
    "baseline_edges_match",
    "cluster_edge_matrix_is_1000_by_276",
    "bootstrap_bridge_matrix_is_1000_by_48",
    "all_bootstrap_bridge_values_finite",
    "original_networktools_equals_manual",
    "audited_bootstrap_networktools_equals_batch_sums",
    "all_intervals_ordered",
    "all_outputs_nonempty"
  ),
  passed = c(
    baseline_edge_max_abs_difference <= 1e-12,
    baseline_edge_max_abs_difference <= 1e-12,
    identical(dim(edge_estimates), c(1000L, 276L)),
    identical(dim(bootstrap_raw_bridge_matrix), c(1000L, 48L)),
    all(is.finite(bootstrap_raw_bridge_matrix)),
    original_networktools_manual_max_abs_difference <= 1e-12,
    bootstrap_networktools_audit_max_abs_difference <= 1e-12,
    all(point_estimates$ci_lower <= point_estimates$ci_upper),
    all(file.info(figure_paths)$size > 0)
  ),
  stringsAsFactors = FALSE
)
if (!all(numerical_audit$passed)) {
  stop("One or more S6.2 numerical/output validation checks failed.")
}

results_s6_2_node_bridge_expected_influence <- list(
  schema_version = "1.2.0",
  analysis_rule = list(
    framework = "Engagement-Harm Differentiation",
    estimand = "Target-specific signed raw one-step Bridge Expected Influence",
    normalization = FALSE,
    network_basis = "Original-sample and bootstrap full unthresholded networks",
    temporal_in = "Signed sum of target-community predictors into each PSMU outcome",
    temporal_out = "Signed sum from each PSMU predictor to target-community outcomes",
    undirected = "Signed sum of PSMU-to-target-community conditional edges",
    self_loops = "Not part of any cross-community Bridge EI sum",
    bootstrap = paste(
      "1,000 participant-cluster nonparametric bootstrap full edge vectors;",
      "percentile 95% CI using R default quantile type"
    ),
    ci_includes_zero = "Unrounded ci_lower <= 0 and ci_upper >= 0",
    threshold = "None; no p-value support mask is used"
  ),
  panel_specification = panel_specification,
  target_specification = target_specification,
  node_map = node_map,
  point_estimates = point_estimates,
  bootstrap_estimates = bootstrap_raw_bridge_matrix,
  bootstrap_intervals = bootstrap_intervals,
  networktools_audit_replications = audit_replications,
  numerical_audit = numerical_audit,
  audit_values = list(
    baseline_edge_max_abs_difference = baseline_edge_max_abs_difference,
    original_networktools_manual_max_abs_difference =
      original_networktools_manual_max_abs_difference,
    bootstrap_networktools_audit_max_abs_difference =
      bootstrap_networktools_audit_max_abs_difference,
    intervals_including_zero = sum(point_estimates$includes_zero),
    intervals_excluding_zero = sum(!point_estimates$includes_zero)
  ),
  figure_specification = list(
    layout = "Two separate 1 x 2 figures; panels A-B within each figure",
    colors = plot_colors,
    bar_fill = "white",
    bar_width = bar_width,
    bar_outline_rule = "Always target-community color",
    point_and_errorbar_rule = "Always target-community color",
    ci_exclusion_star_rule = paste(
      "A community-colored asterisk appears just outside ci_lower for negative",
      "estimates or ci_upper for positive estimates only when the unrounded",
      "percentile 95% CI excludes zero"
    ),
    ci_exclusion_star_size = ci_exclusion_star_size,
    ci_exclusion_star_offset_fraction = ci_exclusion_star_offset_fraction,
    node_separator_rule = paste(
      "Explicit horizontal geom_hline separators at y = 1.5, 2.5, ..., 5.5;",
      "not default major or minor grid lines"
    ),
    point_shape = endpoint_point_shape,
    panel_border = panel_border_color,
    temporal = list(
      width_inches = temporal_figure_width,
      height_inches = temporal_figure_height,
      svg_path = figure_paths[["temporal_svg"]],
      svg_md5 = unname(tools::md5sum(figure_paths[["temporal_svg"]]))
    ),
    undirected = list(
      width_inches = undirected_figure_width,
      height_inches = undirected_figure_height,
      svg_path = figure_paths[["undirected_svg"]],
      svg_md5 = unname(tools::md5sum(figure_paths[["undirected_svg"]]))
    )
  ),
  manuscript_note = manuscript_note,
  provenance = list(
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    script = basename(script_file),
    script_md5 = unname(tools::md5sum(script_file)),
    network_file = network_file,
    network_md5 = network_md5,
    checkpoint_stored_network_md5 = checkpoint_network_md5,
    network_file_md5_matches_checkpoint = network_file_md5_matches_checkpoint,
    network_compatibility_rule = paste(
      "S5.3 file MD5 drift is allowed only because all 276 current full",
      "baseline edges exactly reproduce the checkpoint baseline; this follows",
      "the accepted S5.6 checkpoint-reuse rule"
    ),
    cluster_checkpoint_file = cluster_checkpoint_file,
    cluster_checkpoint_md5 = cluster_checkpoint_md5,
    cluster_replications = nrow(edge_estimates),
    R_version = R.version.string,
    package_versions = vapply(
      required_packages,
      function(package_name) as.character(utils::packageVersion(package_name)),
      character(1)
    )
  )
)

save(
  results_s6_2_node_bridge_expected_influence,
  file = result_file,
  compress = "xz"
)

reload_environment <- new.env(parent = emptyenv())
reloaded_names <- load(result_file, envir = reload_environment)
if (!identical(reloaded_names, "results_s6_2_node_bridge_expected_influence")) {
  stop("The saved S6.2 RDA does not reload as its sole expected object.")
}
reloaded_result <- reload_environment[[reloaded_names]]
if (!all(reloaded_result$numerical_audit$passed) ||
    !identical(dim(reloaded_result$bootstrap_estimates), c(1000L, 48L))) {
  stop("The reloaded S6.2 result object failed its numerical audit.")
}

message("S6.2 completed successfully.")
message("Result RDA: ", normalizePath(result_file, winslash = "/"))
message("Temporal SVG: ", normalizePath(figure_paths[["temporal_svg"]], winslash = "/"))
message("Undirected SVG: ", normalizePath(figure_paths[["undirected_svg"]], winslash = "/"))
message(
  "95% CIs including zero: ", sum(point_estimates$includes_zero),
  "; excluding zero: ", sum(!point_estimates$includes_zero)
)
