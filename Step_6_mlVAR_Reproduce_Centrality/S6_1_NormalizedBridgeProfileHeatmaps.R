# S6.1: Normalized bridge-profile heatmaps from raw-p sparse mlVAR networks
# Author: Kunru Song

options(stringsAsFactors = FALSE, scipen = 999)
try(Sys.setlocale("LC_CTYPE", ".UTF-8"), silent = TRUE)

required_packages <- c(
  "dplyr", "ggplot2", "patchwork", "scales", "svglite", "tidyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required R package(s): ", paste(missing_packages, collapse = ", "),
    ". Install them under R 4.5.2 before rerunning S6.1."
  )
}

# Locate this file under both command-line Rscript execution and interactive
# source(). This deliberately avoids a machine-specific project-root path.
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

  if (file.exists(expected_basename)) {
    candidates <- c(candidates, expected_basename)
  }
  candidates <- unique(candidates[nzchar(candidates)])
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) > 0L) {
    candidates <- normalizePath(
      candidates, winslash = "/", mustWork = TRUE
    )
    candidates <- candidates[basename(candidates) == expected_basename]
  }
  if (length(candidates) == 0L) {
    stop(
      "Cannot locate ", expected_basename, ". Run it with Rscript or call ",
      "source() with the script path."
    )
  }
  candidates[[length(candidates)]]
}

script_file <- resolve_script_file("S6_1_NormalizedBridgeProfileHeatmaps.R")
script_dir <- dirname(script_file)
setwd(script_dir)
message("S6.1 working directory: ", normalizePath(".", winslash = "/"))

analysis_dir <- "../Res_2_Data_for_Analysis"
results_dir <- "../Res_3_Results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

network_input_file <- file.path(analysis_dir, "results_s5_3_mlvar_networks.rda")
result_file <- file.path(
  analysis_dir, "results_s6_1_normalized_bridge_heatmaps.rda"
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
  network_input_file, "results_s5_3_mlvar_networks"
)

required_network_components <- c(
  "full_matrices", "sparse_matrices", "edge_tables", "node_map", "provenance"
)
if (!all(required_network_components %in% names(network_results))) {
  stop("The S5.3 object is missing one or more required components.")
}
if (!identical(
  names(network_results$sparse_matrices),
  c("temporal", "contemporaneous", "between")
)) {
  stop("Unexpected sparse-matrix names in the S5.3 object.")
}

network_input_md5 <- unname(tools::md5sum(network_input_file))
# S5.3 is the sole numerical input. Its exact file hash is recorded below;
# matrix dimensions, canonical node order and raw-p support are validated here.

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

group_levels <- c("Engagement Time", "Peripheral", "Core", "Distress")
group_axis_labels <- c(
  "Engagement Time" = "Engagement\nTime",
  "Peripheral" = "Peripheral",
  "Core" = "Core",
  "Distress" = "Distress"
)
expected_group_sizes <- c(
  "Engagement Time" = 3L, "Peripheral" = 2L, "Core" = 4L, "Distress" = 3L
)
observed_group_sizes <- table(factor(node_map$profile_group, levels = group_levels))

stopifnot(
  nrow(node_map) == 12L,
  !anyNA(node_map$profile_group),
  identical(as.integer(observed_group_sizes), unname(expected_group_sizes)),
  identical(as.character(names(observed_group_sizes)), names(expected_group_sizes))
)

group_nodes <- stats::setNames(
  lapply(group_levels, function(group_name) {
    node_map$node[node_map$profile_group == group_name]
  }),
  group_levels
)

validate_sparse_matrix <- function(matrix_object, network_name) {
  if (!is.matrix(matrix_object) || !identical(dim(matrix_object), c(12L, 12L))) {
    stop(network_name, " sparse matrix must be 12 x 12.")
  }
  if (!identical(rownames(matrix_object), node_map$node) ||
      !identical(colnames(matrix_object), node_map$node)) {
    stop(network_name, " sparse matrix does not follow the canonical node order.")
  }
  if (anyNA(matrix_object) || any(!is.finite(matrix_object))) {
    stop(network_name, " sparse matrix contains a missing or non-finite value.")
  }
  if (network_name != "temporal") {
    if (!isTRUE(all.equal(matrix_object, t(matrix_object), tolerance = 1e-12))) {
      stop(network_name, " matrix must be symmetric.")
    }
    if (max(abs(diag(matrix_object))) > 1e-12) {
      stop(network_name, " matrix must have a structural-zero diagonal.")
    }
  }
  invisible(TRUE)
}

invisible(Map(
  validate_sparse_matrix,
  network_results$sparse_matrices,
  names(network_results$sparse_matrices)
))

expected_sparse_from_support <- function(network_name, directed) {
  full_matrix <- network_results$full_matrices[[network_name]][
    node_map$node, node_map$node, drop = FALSE
  ]
  edge_table <- network_results$edge_tables[[network_name]]
  expected <- matrix(
    0, nrow = 12L, ncol = 12L,
    dimnames = list(node_map$node, node_map$node)
  )

  if (directed) {
    indices <- cbind(
      match(edge_table$from, node_map$node),
      match(edge_table$to, node_map$node)
    )
    expected[indices] <- full_matrix[indices] * as.numeric(edge_table$supported_p05)
  } else {
    indices <- cbind(
      match(edge_table$node1, node_map$node),
      match(edge_table$node2, node_map$node)
    )
    values <- full_matrix[indices] * as.numeric(edge_table$supported_p05)
    expected[indices] <- values
    expected[cbind(indices[, 2L], indices[, 1L])] <- values
  }
  expected
}

expected_sparse_matrices <- list(
  temporal = expected_sparse_from_support("temporal", directed = TRUE),
  contemporaneous = expected_sparse_from_support(
    "contemporaneous", directed = FALSE
  ),
  between = expected_sparse_from_support("between", directed = FALSE)
)
sparse_support_verified <- all(vapply(names(expected_sparse_matrices), function(
  network_name
) {
  isTRUE(all.equal(
    network_results$sparse_matrices[[network_name]],
    expected_sparse_matrices[[network_name]],
    tolerance = 1e-12, check.attributes = TRUE
  ))
}, logical(1)))
if (!sparse_support_verified) {
  stop("S5.3 sparse matrices do not exactly match their raw-p support flags.")
}

# Compute a signed mean per possible edge. No absolute values enter any mean.
# For group diagonals, summarize edges between distinct nodes only: temporal
# autoregressive self-loops and undirected structural-zero diagonals are excluded.
aggregate_group_matrix <- function(matrix_object, directed) {
  group_matrix <- matrix(
    NA_real_, nrow = length(group_levels), ncol = length(group_levels),
    dimnames = list(group_levels, group_levels)
  )
  possible_edges <- matrix(
    NA_integer_, nrow = length(group_levels), ncol = length(group_levels),
    dimnames = list(group_levels, group_levels)
  )

  for (row_group in group_levels) {
    for (column_group in group_levels) {
      row_nodes <- group_nodes[[row_group]]
      column_nodes <- group_nodes[[column_group]]
      block <- matrix_object[row_nodes, column_nodes, drop = FALSE]

      if (identical(row_group, column_group)) {
        if (directed) {
          keep <- outer(row_nodes, column_nodes, FUN = "!=")
          values <- block[keep]
        } else {
          values <- block[lower.tri(block, diag = FALSE)]
        }
      } else {
        values <- as.numeric(block)
      }

      if (length(values) == 0L || any(!is.finite(values))) {
        stop("A group block has no valid edge values: ", row_group, " / ", column_group)
      }
      group_matrix[row_group, column_group] <- mean(values)
      possible_edges[row_group, column_group] <- length(values)
    }
  }

  if (!directed && !isTRUE(all.equal(group_matrix, t(group_matrix), tolerance = 1e-12))) {
    stop("Aggregated undirected group matrix is not symmetric.")
  }

  list(values = group_matrix, possible_edges = possible_edges)
}

profile_objects <- list(
  temporal = aggregate_group_matrix(
    network_results$sparse_matrices$temporal, directed = TRUE
  ),
  contemporaneous = aggregate_group_matrix(
    network_results$sparse_matrices$contemporaneous, directed = FALSE
  ),
  between = aggregate_group_matrix(
    network_results$sparse_matrices$between, directed = FALSE
  )
)

make_plot_data <- function(network_name, lower_triangle_only) {
  value_matrix <- profile_objects[[network_name]]$values
  edge_count_matrix <- profile_objects[[network_name]]$possible_edges
  plot_data <- tidyr::expand_grid(
    row_group = group_levels,
    column_group = group_levels
  ) |>
    dplyr::mutate(
      row_index = match(.data$row_group, group_levels),
      column_index = match(.data$column_group, group_levels),
      value = value_matrix[cbind(.data$row_index, .data$column_index)],
      possible_edges = edge_count_matrix[cbind(.data$row_index, .data$column_index)],
      visible = !lower_triangle_only | .data$row_index >= .data$column_index,
      network = network_name,
      row_group = factor(.data$row_group, levels = rev(group_levels)),
      column_group = factor(.data$column_group, levels = group_levels)
    )
  plot_data
}

profile_plot_data <- dplyr::bind_rows(
  make_plot_data("temporal", lower_triangle_only = FALSE),
  make_plot_data("contemporaneous", lower_triangle_only = TRUE),
  make_plot_data("between", lower_triangle_only = TRUE)
)

# abs() is used only to set a separate symmetric range for each network;
# plotted/annotated values remain signed means, so positive and negative edges
# can cancel. The three color scales are independent and their colorbars hidden.
network_order <- c("temporal", "contemporaneous", "between")
color_limits <- tapply(
  abs(profile_plot_data$value[profile_plot_data$visible]),
  profile_plot_data$network[profile_plot_data$visible],
  max,
  na.rm = TRUE
)[network_order]
color_limits <- stats::setNames(as.numeric(color_limits), network_order)
if (any(!is.finite(color_limits)) || any(color_limits <= 0)) {
  stop("Invalid network-specific color limit.")
}

format_three_decimals <- function(value) {
  rounded_value <- round(value, 3L)
  rounded_value[rounded_value == 0] <- 0
  sprintf("%.3f", rounded_value)
}

profile_plot_data <- profile_plot_data |>
  dplyr::mutate(
    panel_color_limit = unname(color_limits[.data$network]),
    value_label = dplyr::if_else(
      .data$visible, format_three_decimals(.data$value), ""
    ),
    text_color = dplyr::if_else(
      .data$visible & abs(.data$value) / .data$panel_color_limit >= 0.52,
      "#FFFFFF", "#202020"
    )
  )

make_network_fill_scale <- function(network_name) {
  panel_limit <- unname(color_limits[[network_name]])
  ggplot2::scale_fill_gradient2(
    low = "#2166AC", mid = "#FFFFFF", high = "#B2182B",
    midpoint = 0, limits = c(-panel_limit, panel_limit),
    oob = scales::squish, na.value = "#FFFFFF",
    name = NULL, guide = "none"
  )
}

network_titles <- c(
  temporal = "Temporal network",
  contemporaneous = "Contemporaneous network",
  between = "Between-person network"
)

standalone_cell_text_size <- 4.1
combined_cell_text_size <- 5.2

make_heatmap_plot <- function(
  network_name, temporal = FALSE,
  cell_text_size = standalone_cell_text_size
) {
  panel_data <- profile_plot_data |>
    dplyr::filter(.data$network == network_name, .data$visible)

  ggplot2::ggplot(
    panel_data,
    ggplot2::aes(x = .data$column_group, y = .data$row_group, fill = .data$value)
  ) +
    ggplot2::geom_tile(color = "#D9D9D9", linewidth = 0.55) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$value_label, color = .data$text_color),
      size = cell_text_size, fontface = "plain", show.legend = FALSE
    ) +
    ggplot2::scale_color_identity() +
    make_network_fill_scale(network_name) +
    ggplot2::scale_x_discrete(
      drop = FALSE, labels = group_axis_labels, expand = c(0, 0)
    ) +
    ggplot2::scale_y_discrete(
      drop = FALSE, labels = group_axis_labels, expand = c(0, 0)
    ) +
    ggplot2::coord_fixed(clip = "off") +
    ggplot2::labs(
      title = unname(network_titles[[network_name]]),
      x = if (temporal) "Outcome (t)" else NULL,
      y = if (temporal) "Predictor (t - 1)" else NULL
    ) +
    ggplot2::theme_minimal(base_size = 11.5) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        face = "bold", color = "#202020", hjust = 0.5, size = 12.5
      ),
      axis.title = ggplot2::element_text(color = "#202020", size = 11),
      axis.text = ggplot2::element_text(color = "#202020", size = 10.2),
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 5)),
      axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 5)),
      axis.ticks = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )
}

temporal_plot <- make_heatmap_plot("temporal", temporal = TRUE)
contemporaneous_plot <- make_heatmap_plot("contemporaneous", temporal = FALSE)
between_plot <- make_heatmap_plot("between", temporal = FALSE)

# The combined figure intentionally retains only the three network panel titles.
# Overall title/subtitle, axis titles, tick labels, and tick marks are suppressed.
# Equal patchwork weights plus identical fixed-coordinate scales keep all three
# 4 x 4 heatmap plotting regions the same size.
make_combined_heatmap_plot <- function(network_name, temporal = FALSE) {
  make_heatmap_plot(
    network_name = network_name,
    temporal = temporal,
    cell_text_size = combined_cell_text_size
  ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )
}

combined_heatmap_plot_list <- list(
  temporal = make_combined_heatmap_plot("temporal", temporal = TRUE),
  contemporaneous = make_combined_heatmap_plot(
    "contemporaneous", temporal = FALSE
  ),
  between = make_combined_heatmap_plot("between", temporal = FALSE)
)
combined_panel_widths <- rep(1, length(combined_heatmap_plot_list))
combined_plot <- patchwork::wrap_plots(
  combined_heatmap_plot_list,
  nrow = 1L, widths = combined_panel_widths
) &
  ggplot2::theme(legend.position = "none")

temporal_panel_plot <- temporal_plot
contemporaneous_panel_plot <- contemporaneous_plot
between_panel_plot <- between_plot

figure_paths <- list(
  temporal = c(
    svg = file.path(results_dir, "Figure_S6_1_NormalizedBridgeProfiles_Temporal.svg")
  ),
  contemporaneous = c(
    svg = file.path(results_dir, "Figure_S6_1_NormalizedBridgeProfiles_Contemporaneous.svg")
  ),
  between = c(
    svg = file.path(results_dir, "Figure_S6_1_NormalizedBridgeProfiles_BetweenPerson.svg")
  ),
  combined = c(
    svg = file.path(results_dir, "Figure_S6_1_NormalizedBridgeProfiles_Combined.svg")
  )
)

# EDIT THE TEMPORAL HEATMAP WIDTH AND HEIGHT HERE (inches).
temporal_width <- 5.9
temporal_height <- 5.7
ggplot2::ggsave(
  figure_paths$temporal[["svg"]], temporal_panel_plot,
  device = svglite::svglite, width = temporal_width, height = temporal_height,
  units = "in", bg = "white"
)

# EDIT THE CONTEMPORANEOUS HEATMAP WIDTH AND HEIGHT HERE (inches).
contemporaneous_width <- 5.6
contemporaneous_height <- 5.4
ggplot2::ggsave(
  figure_paths$contemporaneous[["svg"]], contemporaneous_panel_plot,
  device = svglite::svglite, width = contemporaneous_width,
  height = contemporaneous_height, units = "in", bg = "white"
)

# EDIT THE BETWEEN-PERSON HEATMAP WIDTH AND HEIGHT HERE (inches).
between_width <- 5.6
between_height <- 5.4
ggplot2::ggsave(
  figure_paths$between[["svg"]], between_panel_plot,
  device = svglite::svglite, width = between_width, height = between_height,
  units = "in", bg = "white"
)

# EDIT THE COMBINED THREE-HEATMAP WIDTH AND HEIGHT HERE (inches).
combined_width <- 16.2
combined_height <- 6.7
ggplot2::ggsave(
  figure_paths$combined[["svg"]], combined_plot,
  device = svglite::svglite, width = combined_width, height = combined_height,
  units = "in", bg = "white"
)

figure_dimensions <- data.frame(
  figure = names(figure_paths),
  width_in = c(
    temporal_width, contemporaneous_width, between_width, combined_width
  ),
  height_in = c(
    temporal_height, contemporaneous_height, between_height, combined_height
  ),
  stringsAsFactors = FALSE
)

figure_audit <- lapply(names(figure_paths), function(figure_name) {
  paths <- figure_paths[[figure_name]]
  data.frame(
    figure = figure_name,
    svg_file = basename(paths[["svg"]]),
    svg_bytes = as.numeric(file.info(paths[["svg"]])$size),
    svg_md5 = unname(tools::md5sum(paths[["svg"]])),
    saved_by_ggsave = TRUE,
    stringsAsFactors = FALSE
  )
}) |>
  dplyr::bind_rows() |>
  dplyr::left_join(figure_dimensions, by = "figure")

heatmap_plot_list <- list(
  temporal = temporal_plot,
  contemporaneous = contemporaneous_plot,
  between = between_plot
)
plotted_color_limits <- vapply(heatmap_plot_list, function(plot_object) {
  max(abs(plot_object$scales$get_scales("fill")$limits))
}, numeric(1))
network_specific_scales_verified <- isTRUE(all.equal(
  unname(plotted_color_limits), unname(color_limits), tolerance = 1e-12
))

has_visible_guide <- function(grob_object) {
  guide_indices <- grep("^guide-box", grob_object$layout$name)
  if (length(guide_indices) == 0L) return(FALSE)
  any(vapply(guide_indices, function(index) {
    !inherits(grob_object$grobs[[index]], "zeroGrob")
  }, logical(1)))
}

guide_visibility_audit <- local({
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  c(
    standalone_hidden = all(vapply(heatmap_plot_list, function(plot_object) {
      identical(plot_object$scales$get_scales("fill")$guide, "none") &&
        identical(plot_object$theme$legend.position, "none") &&
        !has_visible_guide(ggplot2::ggplotGrob(plot_object))
    }, logical(1))),
    combined_hidden = !has_visible_guide(patchwork::patchworkGrob(combined_plot))
  )
})
standalone_colorbars_hidden <- unname(guide_visibility_audit[["standalone_hidden"]])
combined_colorbar_hidden <- unname(guide_visibility_audit[["combined_hidden"]])

combined_axis_theme_elements <- c(
  "axis.title.x", "axis.title.y", "axis.text.x", "axis.text.y",
  "axis.ticks.x", "axis.ticks.y"
)
combined_axes_hidden_verified <- all(vapply(
  combined_heatmap_plot_list,
  function(plot_object) {
    all(vapply(combined_axis_theme_elements, function(element_name) {
      inherits(plot_object$theme[[element_name]], "element_blank")
    }, logical(1))) &&
      is.null(plot_object$labels$x) && is.null(plot_object$labels$y)
  },
  logical(1)
))
combined_panel_titles_verified <- identical(
  unname(vapply(
    combined_heatmap_plot_list,
    function(plot_object) plot_object$labels$title,
    character(1)
  )),
  unname(network_titles)
)
combined_overall_annotation_absent <-
  is.null(combined_plot$patches$annotation$title) &&
  is.null(combined_plot$patches$annotation$subtitle)
combined_equal_panel_sizes_verified <-
  length(combined_panel_widths) == 3L &&
  length(unique(combined_panel_widths)) == 1L
combined_cell_text_size_verified <- all(vapply(
  combined_heatmap_plot_list,
  function(plot_object) {
    isTRUE(all.equal(
      plot_object$layers[[2L]]$aes_params$size,
      combined_cell_text_size,
      tolerance = 0
    ))
  },
  logical(1)
)) && combined_cell_text_size > standalone_cell_text_size
combined_canvas_unchanged_verified <-
  isTRUE(all.equal(combined_width, 16.2, tolerance = 0)) &&
  isTRUE(all.equal(combined_height, 6.7, tolerance = 0))

visible_labels <- profile_plot_data$value_label[profile_plot_data$visible]
validation <- data.frame(
  criterion = c(
    "Three raw-p sparse matrices are 12 x 12 in canonical node order",
    "Sparse matrices exactly apply the S5.3 original-sample raw-p support flags",
    "Four display groups have sizes 3, 2, 4, and 3",
    "Temporal plot displays all 16 cells",
    "Each undirected plot displays the 10-cell lower triangle",
    "All displayed labels retain exactly three decimals",
    "Signed means are used without absolute-value aggregation",
    "Within-group diagonals exclude node self-loops/structural zeros",
    "Each network uses its own symmetric blue-white-red color scale",
    "Colorbars are hidden in standalone and combined figures",
    "Combined figure retains panel titles but has no overall title or subtitle",
    "Combined figure hides all axis titles, tick labels, and tick marks",
    "Combined figure uses equal panel sizes and larger cell text",
    "Combined figure retains the locked 16.2 x 6.7 inch canvas",
    "All four SVGs were saved by ggsave",
    "No Step 5 artifact is overwritten"
  ),
  passed = c(
    TRUE,
    sparse_support_verified,
    identical(as.integer(observed_group_sizes), c(3L, 2L, 4L, 3L)),
    sum(profile_plot_data$network == "temporal" & profile_plot_data$visible) == 16L,
    all(vapply(c("contemporaneous", "between"), function(network_name) {
      sum(profile_plot_data$network == network_name & profile_plot_data$visible) == 10L
    }, logical(1))),
    all(grepl("^-?[0-9]+\\.[0-9]{3}$", visible_labels)),
    TRUE,
    all(diag(profile_objects$temporal$possible_edges) == c(6L, 2L, 12L, 6L)) &&
      all(diag(profile_objects$contemporaneous$possible_edges) == c(3L, 1L, 6L, 3L)) &&
      all(diag(profile_objects$between$possible_edges) == c(3L, 1L, 6L, 3L)),
    network_specific_scales_verified,
    standalone_colorbars_hidden && combined_colorbar_hidden,
    combined_panel_titles_verified && combined_overall_annotation_absent,
    combined_axes_hidden_verified,
    combined_equal_panel_sizes_verified && combined_cell_text_size_verified,
    combined_canvas_unchanged_verified,
    all(figure_audit$svg_bytes > 1000) && all(figure_audit$saved_by_ggsave),
    all(grepl("^Figure_S6_1_", figure_audit$svg_file))
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) {
  stop(
    "S6.1 validation failed: ",
    paste(validation$criterion[!validation$passed], collapse = "; ")
  )
}

results_s6_1_normalized_bridge_heatmaps <- list(
  profile_matrices = lapply(profile_objects, `[[`, "values"),
  possible_edge_counts = lapply(profile_objects, `[[`, "possible_edges"),
  plot_data = profile_plot_data |>
    dplyr::mutate(
      row_group = as.character(.data$row_group),
      column_group = as.character(.data$column_group)
    ),
  group_map = node_map |>
    dplyr::select("order", "node", "short_label", "profile_group"),
  aggregation_rule = list(
    statistic = "Signed arithmetic mean per possible edge in each group block",
    network_version = "Raw-p sparse",
    support_rule = paste(
      "Original-sample unadjusted two-sided p < .05; undirected networks use",
      "the S5.3 conservative node-wise AND rule"
    ),
    unsupported_edge_value = 0,
    absolute_value_aggregation = FALSE,
    temporal_orientation = "Rows are predictors at t - 1; columns are outcomes at t",
    temporal_diagonal = "Signed mean of within-group directed edges excluding node self-loops",
    undirected_diagonal = "Signed mean of unique within-group edges excluding structural zeros",
    undirected_display = "Lower triangle including the group diagonal"
  ),
  color_scale = list(
    shared_across_networks = FALSE,
    limits = data.frame(
      network = names(color_limits),
      lower_limit = -unname(color_limits),
      upper_limit = unname(color_limits),
      stringsAsFactors = FALSE
    ),
    colorbar_visible = FALSE,
    negative_color = "#2166AC",
    zero_color = "#FFFFFF",
    positive_color = "#B2182B"
  ),
  combined_figure_specification = list(
    panel_titles_visible = TRUE,
    overall_title_visible = FALSE,
    overall_subtitle_visible = FALSE,
    axis_titles_visible = FALSE,
    axis_tick_labels_visible = FALSE,
    axis_tick_marks_visible = FALSE,
    equal_panel_widths = combined_panel_widths,
    standalone_cell_text_size = standalone_cell_text_size,
    combined_cell_text_size = combined_cell_text_size,
    width_in = combined_width,
    height_in = combined_height
  ),
  figure_audit = figure_audit,
  validation = validation,
  provenance = list(
    network_input = "Res_2_Data_for_Analysis/results_s5_3_mlvar_networks.rda",
    network_input_md5 = network_input_md5,
    network_version = "Raw-p sparse",
    supported_edges = vapply(
      network_results$edge_tables,
      function(edge_table) sum(edge_table$supported_p05),
      integer(1)
    ),
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

save(results_s6_1_normalized_bridge_heatmaps, file = result_file, version = 2)
rm(results_s6_1_normalized_bridge_heatmaps)
results_s6_1_normalized_bridge_heatmaps <- load_single_object(
  result_file, "results_s6_1_normalized_bridge_heatmaps"
)
stopifnot(
  all(results_s6_1_normalized_bridge_heatmaps$validation$passed),
  identical(
    results_s6_1_normalized_bridge_heatmaps$provenance$network_input_md5,
    network_input_md5
  ),
  nrow(results_s6_1_normalized_bridge_heatmaps$plot_data) == 48L
)

cat("S6.1 completed successfully.\n")
cat("Network-specific signed color limits:\n")
print(data.frame(
  network = names(color_limits),
  limit = unname(color_limits),
  row.names = NULL
))
cat("All colorbars hidden: TRUE\n")
cat("Saved result object: ", result_file, "\n", sep = "")
cat("Saved figures:\n")
cat(paste0("- ", unlist(figure_paths), collapse = "\n"), "\n")
