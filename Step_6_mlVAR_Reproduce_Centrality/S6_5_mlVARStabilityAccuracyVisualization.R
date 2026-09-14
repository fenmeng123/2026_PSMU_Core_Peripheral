# S6.5: mlVAR stability and edge-accuracy visualization
# Author: Kunru Song

options(stringsAsFactors = FALSE, scipen = 999)
try(Sys.setlocale("LC_CTYPE", ".UTF-8"), silent = TRUE)

required_packages <- c(
  "dplyr", "ggplot2", "patchwork", "scales", "svglite"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing required R package(s): ", paste(missing_packages, collapse = ", "),
    ". Install them under R 4.5.2 before rerunning S6.5."
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
    if (is.null(frame_file) || !length(frame_file)) return(NA_character_)
    as.character(frame_file[[1L]])
  }, character(1))
  candidates <- c(candidates, frame_files[!is.na(frame_files)])

  if (file.exists(expected_basename)) candidates <- c(candidates, expected_basename)
  candidates <- unique(candidates[nzchar(candidates)])
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates)) {
    candidates <- normalizePath(candidates, winslash = "/", mustWork = TRUE)
    candidates <- candidates[basename(candidates) == expected_basename]
  }
  if (!length(candidates)) {
    stop(
      "Cannot locate ", expected_basename, ". Run it with Rscript or call ",
      "source() with the script path."
    )
  }
  candidates[[length(candidates)]]
}

script_file <- resolve_script_file("S6_5_mlVARStabilityAccuracyVisualization.R")
script_dir <- dirname(script_file)
setwd(script_dir)
message("S6.5 working directory: ", normalizePath(".", winslash = "/"))

analysis_dir <- "../Res_2_Data_for_Analysis"
results_dir <- "../Res_3_Results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

network_file <- file.path(analysis_dir, "results_s5_3_mlvar_networks.rda")
sparse_checkpoint_file <- file.path(
  analysis_dir, "checkpoint_s6_4_nonparametric_bootstrap_accuracy.rda"
)
case_drop_edge_checkpoint_file <- file.path(
  analysis_dir, "checkpoint_s6_3_case_drop_edges.rda"
)
result_file <- file.path(
  analysis_dir, "results_s6_5_mlvar_stability_accuracy.rda"
)
figure_paths <- c(
  svg = file.path(results_dir, "Figure_S6_5_mlVAR_StabilityAccuracy_A4.svg")
)

load_single_object <- function(path, expected_name) {
  if (!file.exists(path)) {
    stop(
      "Required input does not exist: ", path, ". Run ",
      "S6_4_NonParametricBootstrapAccuracy.R first."
    )
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
sparse_checkpoint <- load_single_object(
  sparse_checkpoint_file, "checkpoint_s6_4_nonparametric_bootstrap_accuracy"
)
case_drop_edge_checkpoint <- load_single_object(
  case_drop_edge_checkpoint_file, "checkpoint_s6_3_case_drop_edges"
)

network_md5 <- unname(tools::md5sum(network_file))
sparse_checkpoint_md5 <- unname(tools::md5sum(sparse_checkpoint_file))
case_drop_edge_checkpoint_md5 <-
  unname(tools::md5sum(case_drop_edge_checkpoint_file))
network_order <- c("Temporal", "Contemporaneous", "Between-person")
expected_edge_counts <- c(
  Temporal = 144L, Contemporaneous = 66L, `Between-person` = 66L
)
drop_proportions <- c(.10, .20, .30, .40, .50, .60, .70, .75)

required_checkpoint_components <- c(
  "schema_version", "replications_per_resampling", "parallel_workers",
  "internal_mlvar_cores", "raw_p_threshold", "temporal_support_rule",
  "undirected_support_rule", "network_md5", "baseline_sparse_edge_vectors",
  "baseline_sparse_edges", "cluster",
  "case_drop", "validation"
)
stopifnot(
  all(required_checkpoint_components %in% names(sparse_checkpoint)),
  sparse_checkpoint$schema_version %in% c("2.0.0", "2.1.0"),
  sparse_checkpoint$replications_per_resampling == 1000L,
  sparse_checkpoint$parallel_workers == 6L,
  sparse_checkpoint$internal_mlvar_cores == 1L,
  identical(sparse_checkpoint$raw_p_threshold, .05),
  identical(sparse_checkpoint$network_md5, network_md5),
  all(unlist(sparse_checkpoint$validation))
)

stopifnot(
  case_drop_edge_checkpoint$schema_version %in% c("1.1.0", "1.2.0"),
  case_drop_edge_checkpoint$replications == 1000L,
  identical(dim(case_drop_edge_checkpoint$edge_estimates), c(1000L, 276L)),
  identical(
    colnames(case_drop_edge_checkpoint$edge_estimates),
    names(case_drop_edge_checkpoint$baseline_edges)
  ),
  identical(case_drop_edge_checkpoint$network_md5, network_md5),
  identical(
    case_drop_edge_checkpoint$input_md5, sparse_checkpoint$input_md5
  ),
  identical(
    case_drop_edge_checkpoint$specification, sparse_checkpoint$specification
  ),
  all(is.finite(case_drop_edge_checkpoint$edge_estimates)),
  all(is.finite(case_drop_edge_checkpoint$baseline_edges)),
  nrow(case_drop_edge_checkpoint$status) == 1000L,
  all(case_drop_edge_checkpoint$status$status == "success"),
  identical(
    as.numeric(case_drop_edge_checkpoint$status$drop_proportion),
    as.numeric(rep(drop_proportions, each = 125L))
  ),
  all(unlist(case_drop_edge_checkpoint$validation))
)

for (resampling in c("cluster", "case_drop")) {
  object <- sparse_checkpoint[[resampling]]
  stopifnot(
    identical(dim(object$p_values), c(1000L, 276L)),
    identical(dim(object$support), c(1000L, 276L)),
    identical(dim(object$sparse_edges), c(1000L, 276L)),
    identical(object$support, is.finite(object$p_values) & object$p_values < .05),
    all(object$sparse_edges[!object$support] == 0),
    all(is.finite(object$sparse_edges)),
    nrow(object$status) == 1000L,
    all(object$status$status == "success")
  )
  for (network_name in network_order) {
    stopifnot(
      identical(
        dim(object$sparse_edges_by_network[[network_name]]),
        c(1000L, unname(expected_edge_counts[[network_name]]))
      ),
      identical(
        colnames(object$sparse_edges_by_network[[network_name]]),
        names(sparse_checkpoint$baseline_sparse_edge_vectors[[network_name]])
      )
    )
  }
}
stopifnot(
  identical(
    as.numeric(sparse_checkpoint$case_drop$status$drop_proportion),
    as.numeric(rep(drop_proportions, each = 125L))
  ),
  identical(
    names(sparse_checkpoint$baseline_sparse_edges),
    colnames(sparse_checkpoint$cluster$sparse_edges)
  ),
  identical(
    names(sparse_checkpoint$baseline_sparse_edges),
    colnames(sparse_checkpoint$case_drop$sparse_edges)
  )
)

constant_row_mask <- function(estimate_matrix) {
  apply(estimate_matrix, 1L, function(values) {
    all(values == values[[1L]])
  })
}

safe_row_correlations <- function(estimate_matrix, baseline_vector, label) {
  if (all(baseline_vector == baseline_vector[[1L]])) {
    stop("The original-sample vector has zero variance for ", label, ".")
  }
  zero_variance <- constant_row_mask(estimate_matrix)
  correlations <- numeric(nrow(estimate_matrix))
  if (any(!zero_variance)) {
    correlations[!zero_variance] <- apply(
      estimate_matrix[!zero_variance, , drop = FALSE],
      1L, stats::cor,
      y = baseline_vector, use = "everything"
    )
  }
  # A completely constant bootstrap vector has an undefined Pearson
  # correlation. Code it conservatively as zero rather than dropping that
  # replication and thereby inflating stability.
  correlations[zero_variance] <- 0
  if (length(correlations) != nrow(estimate_matrix) || any(!is.finite(correlations))) {
    stop("Non-finite sparse-network correlation encountered for ", label, ".")
  }
  as.numeric(correlations)
}

summarize_stability <- function(data) {
  data |>
    dplyr::group_by(
      .data$panel_id, .data$network, .data$metric,
      .data$series, .data$drop_proportion
    ) |>
    dplyr::summarise(
      retained_participants = dplyr::first(.data$retained_participants),
      successful_replications = dplyr::n(),
      mean_correlation = mean(.data$correlation),
      median_correlation = stats::median(.data$correlation),
      correlation_q025 = stats::quantile(
        .data$correlation, .025, names = FALSE, type = 7L
      ),
      correlation_q05 = stats::quantile(
        .data$correlation, .05, names = FALSE, type = 7L
      ),
      correlation_q95 = stats::quantile(
        .data$correlation, .95, names = FALSE, type = 7L
      ),
      correlation_q975 = stats::quantile(
        .data$correlation, .975, names = FALSE, type = 7L
      ),
      probability_at_least_0_70 = mean(.data$correlation >= .70),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      .data$panel_id, .data$series, .data$drop_proportion
    )
}

calculate_cs_0_70 <- function(summary_data) {
  summary_data |>
    dplyr::group_by(
      .data$panel_id, .data$network, .data$metric, .data$series
    ) |>
    dplyr::arrange(.data$drop_proportion, .by_group = TRUE) |>
    dplyr::summarise(
      cs_coefficient = {
        contiguous_pass <- cumprod(.data$correlation_q05 >= .70) == 1L
        if (any(contiguous_pass)) {
          max(.data$drop_proportion[contiguous_pass])
        } else {
          0
        }
      },
      .groups = "drop"
    )
}

edge_panel_map <- data.frame(
  panel_id = c("D", "E", "F"),
  network = network_order,
  metric = "Edge weights",
  series = "Edge weights",
  stringsAsFactors = FALSE
)
edge_stability_correlations <- lapply(seq_len(nrow(edge_panel_map)), function(i) {
  network_name <- edge_panel_map$network[[i]]
  edge_estimates <-
    sparse_checkpoint$case_drop$sparse_edges_by_network[[network_name]]
  data.frame(
    replication = sparse_checkpoint$case_drop$status$replication,
    drop_proportion = sparse_checkpoint$case_drop$status$drop_proportion,
    retained_participants =
      sparse_checkpoint$case_drop$status$retained_participants,
    panel_id = edge_panel_map$panel_id[[i]],
    network = network_name,
    metric = "Edge weights",
    series = "Edge weights",
    network_basis = "Replication-specific raw-p sparse",
    zero_variance_input = constant_row_mask(edge_estimates),
    correlation = safe_row_correlations(
      edge_estimates,
      sparse_checkpoint$baseline_sparse_edge_vectors[[network_name]],
      paste(network_name, "Edge Stability")
    ),
    stringsAsFactors = FALSE
  )
}) |>
  dplyr::bind_rows()
edge_stability_summary <- summarize_stability(edge_stability_correlations)
edge_stability_cs <- calculate_cs_0_70(edge_stability_summary)

selected_centrality_metrics <- data.frame(
  panel_id = c("A", "A", "B", "C"),
  network = c("Temporal", "Temporal", "Contemporaneous", "Between-person"),
  metric = c(
    "Overall bridge incoming expected influence (one-step)",
    "Overall bridge outgoing expected influence (one-step)",
    "Overall bridge expected influence (one-step)",
    "Overall bridge expected influence (one-step)"
  ),
  direction = c("incoming", "outgoing", "undirected", "undirected"),
  series = c("Bridge In-EI", "Bridge Out-EI", "Bridge EI", "Bridge EI"),
  stringsAsFactors = FALSE
)

node_map <- network_results$node_map |>
  dplyr::arrange(.data$order)
node_order <- node_map$node
node_community <- stats::setNames(node_map$community, node_map$node)

# Build a node-by-edge incidence matrix for signed one-step bridge EI. Each
# retained cross-community edge contributes its signed weight to the source
# node (outgoing), target node (incoming), or both endpoints (undirected).
build_bridge_incidence <- function(edge_keys, network_name, direction) {
  edge_parts <- strsplit(edge_keys, "|", fixed = TRUE)
  key_network <- vapply(edge_parts, `[[`, character(1), 1L)
  source_node <- vapply(edge_parts, `[[`, character(1), 2L)
  target_node <- vapply(edge_parts, `[[`, character(1), 3L)
  if (!all(key_network == network_name)) {
    stop("Edge keys do not match the requested network: ", network_name)
  }
  cross_community <-
    node_community[source_node] != node_community[target_node]
  incidence <- matrix(
    0, nrow = length(node_order), ncol = length(edge_keys),
    dimnames = list(node_order, edge_keys)
  )
  edge_index <- seq_along(edge_keys)
  if (direction == "outgoing") {
    incidence[cbind(match(source_node, node_order), edge_index)] <-
      as.numeric(cross_community)
  } else if (direction == "incoming") {
    incidence[cbind(match(target_node, node_order), edge_index)] <-
      as.numeric(cross_community)
  } else if (direction == "undirected") {
    incidence[cbind(match(source_node, node_order), edge_index)] <-
      as.numeric(cross_community)
    incidence[cbind(match(target_node, node_order), edge_index)] <-
      as.numeric(cross_community)
  } else {
    stop("Unknown bridge-EI direction: ", direction)
  }
  incidence
}

centrality_estimate_objects <- lapply(
  seq_len(nrow(selected_centrality_metrics)), function(i) {
    specification <- selected_centrality_metrics[i, ]
    edge_estimates <-
      case_drop_edge_checkpoint$edge_estimates_by_network[[specification$network]]
    baseline_edges <-
      case_drop_edge_checkpoint$baseline_edge_vectors[[specification$network]]
    incidence <- build_bridge_incidence(
      colnames(edge_estimates), specification$network, specification$direction
    )
    centrality_estimates <- edge_estimates %*% t(incidence)
    baseline_centrality <- as.numeric(incidence %*% baseline_edges)
    names(baseline_centrality) <- rownames(incidence)
    colnames(centrality_estimates) <- rownames(incidence)
    list(
      specification = specification,
      incidence = incidence,
      estimates = centrality_estimates,
      baseline = baseline_centrality
    )
  }
)

centrality_stability_correlations <- lapply(
  centrality_estimate_objects, function(centrality_object) {
    specification <- centrality_object$specification
    data.frame(
      replication = case_drop_edge_checkpoint$status$replication,
      drop_proportion = case_drop_edge_checkpoint$status$drop_proportion,
      retained_participants =
        case_drop_edge_checkpoint$status$retained_participants,
      panel_id = specification$panel_id,
      network = specification$network,
      metric = specification$metric,
      series = specification$series,
      network_basis = "Full unthresholded network",
      zero_variance_input = constant_row_mask(centrality_object$estimates),
      correlation = safe_row_correlations(
        centrality_object$estimates,
        centrality_object$baseline,
        paste(specification$network, specification$metric)
      ),
      stringsAsFactors = FALSE
    )
  }
) |>
  dplyr::bind_rows()
centrality_stability_summary <- summarize_stability(
  centrality_stability_correlations
)
centrality_cs <- calculate_cs_0_70(centrality_stability_summary)

# The temporal overall one-step formulas must reproduce the independently summed
# temporal bridge EI stored by S6.3. This audit fixes the incoming/outgoing
# orientation while B-C add the signed undirected analogue requested here.
centrality_definition_audit <- dplyr::bind_rows(lapply(
  centrality_estimate_objects[1:2], function(centrality_object) {
    direction <- centrality_object$specification$direction
    reference_metric <- if (direction == "incoming") {
      "Bridge incoming expected influence"
    } else {
      "Bridge outgoing expected influence"
    }
    reference_keys <- paste(
      "Temporal", reference_metric, node_order, sep = "|"
    )
    reference_columns <- match(
      reference_keys, colnames(case_drop_edge_checkpoint$centrality_reproduction)
    )
    if (anyNA(reference_columns)) {
      stop("S6.3 temporal centrality audit columns are incomplete.")
    }
    maximum_difference <- max(abs(
      centrality_object$estimates -
        case_drop_edge_checkpoint$centrality_reproduction[, reference_columns]
    ))
    data.frame(
      network = "Temporal",
      direction = direction,
      comparison = "Incidence-matrix product versus S6.3 direct signed edge-sum",
      maximum_absolute_difference = maximum_difference,
      tolerance = 1e-12,
      passed = maximum_difference < 1e-12,
      stringsAsFactors = FALSE
    )
  }
))
stopifnot(
  nrow(centrality_definition_audit) == 2L,
  all(centrality_definition_audit$passed)
)

zero_variance_correlation_summary <- dplyr::bind_rows(
  edge_stability_correlations,
  centrality_stability_correlations
) |>
  dplyr::filter(.data$zero_variance_input) |>
  dplyr::count(
    .data$panel_id, .data$network, .data$metric,
    .data$series, .data$drop_proportion,
    name = "replications"
  ) |>
  dplyr::arrange(.data$panel_id, .data$series, .data$drop_proportion)

edge_accuracy_plot_data <- lapply(network_order, function(network_name) {
  estimates <- sparse_checkpoint$cluster$sparse_edges_by_network[[network_name]]
  baseline <- sparse_checkpoint$baseline_sparse_edge_vectors[[network_name]]
  keys <- colnames(estimates)
  parts <- strsplit(keys, "|", fixed = TRUE)
  result <- data.frame(
    network = network_name,
    edge_key = keys,
    node1 = vapply(parts, `[[`, character(1), 2L),
    node2 = vapply(parts, `[[`, character(1), 3L),
    original_weight = as.numeric(baseline),
    cluster_bootstrap_median = apply(estimates, 2L, stats::median),
    cluster_ci_lower = apply(
      estimates, 2L, stats::quantile,
      probs = .025, names = FALSE, type = 7L
    ),
    cluster_ci_upper = apply(
      estimates, 2L, stats::quantile,
      probs = .975, names = FALSE, type = 7L
    ),
    successful_replications = nrow(estimates),
    stringsAsFactors = FALSE
  )
  # Case-insensitive, bytewise tie-breaking preserves the frozen alphabetical
  # order of zero edges independently of operating-system collation settings.
  result <- result[order(
    result$original_weight, tolower(result$edge_key), result$edge_key,
    method = "radix"
  ), ]
  result$edge_rank <- seq_len(nrow(result))
  result$edge_count <- nrow(result)
  result
}) |>
  dplyr::bind_rows()

replication_support_summary <- dplyr::bind_rows(lapply(
  c("cluster", "case_drop"), function(resampling) {
    counts <- sparse_checkpoint[[resampling]]$support_counts
    data.frame(
      resampling = resampling,
      replication = seq_len(nrow(counts)),
      Temporal = counts[, "Temporal"],
      Contemporaneous = counts[, "Contemporaneous"],
      `Between-person` = counts[, "Between-person"],
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
))

network_colors <- c(
  Temporal = "#3C5488",
  Contemporaneous = "#00A087",
  `Between-person` = "#E64B35"
)
temporal_bridge_colors <- c(
  `Bridge In-EI` = "#3C5488",
  `Bridge Out-EI` = "#E64B35"
)
accuracy_original_color <- "#E64B35"
accuracy_bootstrap_median_color <- "#111111"
temporal_legend_inside <- c(.975, .025)
default_cs_annotation_y <- .055
temporal_cs_annotation_y <- .160
reference_lines <- data.frame(
  threshold = c(.30, .50, .70),
  label = c("r = 0.3", "r = 0.5", "r = 0.7"),
  linetype = factor(c("0.3", "0.5", "0.7"), levels = c("0.3", "0.5", "0.7")),
  x_label = .745,
  stringsAsFactors = FALSE
)
panel_border_color <- "#111111"
panel_border_linewidth <- .55

base_panel_theme <- function() {
  ggplot2::theme_minimal(base_size = 7.6) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        color = "#E7E7E7", linewidth = .28
      ),
      panel.border = ggplot2::element_rect(
        color = panel_border_color, fill = NA,
        linewidth = panel_border_linewidth
      ),
      plot.title = ggplot2::element_text(
        face = "bold", color = "#202020", hjust = .5,
        size = 8.4, margin = ggplot2::margin(b = 2.5)
      ),
      axis.title = ggplot2::element_text(color = "#202020", size = 6.9),
      axis.text = ggplot2::element_text(color = "#202020", size = 6.1),
      axis.ticks = ggplot2::element_line(color = "#202020", linewidth = .3),
      axis.ticks.length = grid::unit(1.3, "mm"),
      plot.margin = ggplot2::margin(4.5, 5, 4, 4.5)
    )
}

make_stability_panel <- function(
    panel_data, title_text, color_values, cs_annotation,
    show_legend = FALSE, cs_annotation_y = default_cs_annotation_y) {
  panel_data$series <- factor(panel_data$series, levels = names(color_values))
  observed_minimum <- min(panel_data$correlation_q025)
  y_minimum <- max(-1, min(0, floor(observed_minimum * 10) / 10))

  panel <- ggplot2::ggplot(
    panel_data,
    ggplot2::aes(
      x = .data$drop_proportion,
      y = .data$mean_correlation,
      group = .data$series,
      color = .data$series,
      fill = .data$series
    )
  ) +
    ggplot2::geom_hline(
      data = reference_lines,
      ggplot2::aes(yintercept = .data$threshold, linetype = .data$linetype),
      inherit.aes = FALSE, color = "#777777", linewidth = .38
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = .data$correlation_q025,
        ymax = .data$correlation_q975
      ),
      alpha = .14, color = NA, show.legend = FALSE
    ) +
    ggplot2::geom_line(linewidth = .75, show.legend = show_legend) +
    ggplot2::geom_point(size = 1.55, stroke = .2, show.legend = show_legend) +
    ggplot2::geom_text(
      data = reference_lines,
      ggplot2::aes(x = .data$x_label, y = .data$threshold, label = .data$label),
      inherit.aes = FALSE, hjust = 1, vjust = -.45,
      size = 1.85, color = "#555555"
    ) +
    ggplot2::annotate(
      "label", x = .105, y = cs_annotation_y, label = cs_annotation,
      hjust = 0, vjust = 0, size = 1.85, lineheight = .90,
      color = "#202020", fill = scales::alpha("white", .88),
      linewidth = 0, label.padding = grid::unit(.7, "mm")
    ) +
    ggplot2::scale_color_manual(values = color_values, drop = FALSE) +
    ggplot2::scale_fill_manual(values = color_values, drop = FALSE) +
    ggplot2::scale_linetype_manual(
      values = c("0.3" = "dotted", "0.5" = "dashed", "0.7" = "dotdash"),
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(
      breaks = drop_proportions,
      labels = scales::label_percent(accuracy = 1),
      limits = c(.095, .755),
      expand = ggplot2::expansion(mult = c(.01, .01))
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::breaks_pretty(n = 5),
      labels = scales::label_number(accuracy = .01),
      limits = c(y_minimum, 1.02),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
      title = title_text,
      x = "Participants dropped",
      y = "Average correlation with original sample",
      color = NULL, fill = NULL
    ) +
    base_panel_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45, hjust = 1, vjust = 1, size = 5.8
      )
    )

  if (show_legend) {
    panel +
      ggplot2::theme(
        legend.position = "inside",
        legend.position.inside = temporal_legend_inside,
        legend.justification = c(1, 0),
        legend.direction = "vertical",
        legend.background = ggplot2::element_rect(
          fill = scales::alpha("white", .82), color = NA
        ),
        legend.key.height = grid::unit(2.3, "mm"),
        legend.key.width = grid::unit(3.5, "mm"),
        legend.text = ggplot2::element_text(size = 5.6),
        legend.margin = ggplot2::margin(1, 1, 1, 1)
      ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(
          override.aes = list(linewidth = .75, size = 1.4)
        )
      )
  } else {
    panel + ggplot2::theme(legend.position = "none")
  }
}

cs_text <- function(cs_data, network_name, metric_name = NULL) {
  selected <- cs_data[cs_data$network == network_name, ]
  if (!is.null(metric_name)) selected <- selected[selected$metric == metric_name, ]
  if (nrow(selected) != 1L) stop("Expected one CS value for ", network_name)
  sprintf("CS (0.7) = %.3f", selected$cs_coefficient)
}

temporal_cs <- centrality_cs |>
  dplyr::filter(.data$network == "Temporal") |>
  dplyr::arrange(match(.data$series, names(temporal_bridge_colors)))
temporal_cs_text <- paste(
  sprintf("%s: CS (0.7) = %.3f", temporal_cs$series, temporal_cs$cs_coefficient),
  collapse = "\n"
)

panel_a <- make_stability_panel(
  dplyr::filter(centrality_stability_summary, .data$panel_id == "A"),
  "Temporal Bridge EI\nStability",
  temporal_bridge_colors,
  temporal_cs_text,
  show_legend = TRUE,
  cs_annotation_y = temporal_cs_annotation_y
)
panel_b <- make_stability_panel(
  dplyr::filter(centrality_stability_summary, .data$panel_id == "B"),
  "Contemporaneous Bridge EI\nStability",
  c(`Bridge EI` = unname(network_colors[["Contemporaneous"]])),
  cs_text(centrality_cs, "Contemporaneous")
)
panel_c <- make_stability_panel(
  dplyr::filter(centrality_stability_summary, .data$panel_id == "C"),
  "Between-person Bridge EI\nStability",
  c(`Bridge EI` = unname(network_colors[["Between-person"]])),
  cs_text(centrality_cs, "Between-person")
)
panel_d <- make_stability_panel(
  dplyr::filter(edge_stability_summary, .data$panel_id == "D"),
  "Temporal Edge\nStability",
  c(`Edge weights` = unname(network_colors[["Temporal"]])),
  cs_text(edge_stability_cs, "Temporal")
)
panel_e <- make_stability_panel(
  dplyr::filter(edge_stability_summary, .data$panel_id == "E"),
  "Contemporaneous Edge\nStability",
  c(`Edge weights` = unname(network_colors[["Contemporaneous"]])),
  cs_text(edge_stability_cs, "Contemporaneous")
)
panel_f <- make_stability_panel(
  dplyr::filter(edge_stability_summary, .data$panel_id == "F"),
  "Between-person Edge\nStability",
  c(`Edge weights` = unname(network_colors[["Between-person"]])),
  cs_text(edge_stability_cs, "Between-person")
)

make_accuracy_panel <- function(network_name, title_text) {
  panel_data <- edge_accuracy_plot_data |>
    dplyr::filter(.data$network == network_name) |>
    dplyr::arrange(.data$edge_rank)
  edge_count <- unique(panel_data$edge_count)
  y_range <- range(c(panel_data$cluster_ci_lower, panel_data$cluster_ci_upper))
  y_span <- diff(y_range)
  annotation_y <- y_range[[2L]] - .035 * y_span

  ggplot2::ggplot(panel_data, ggplot2::aes(x = .data$edge_rank)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$cluster_ci_lower, ymax = .data$cluster_ci_upper),
      fill = "#BDBDBD", alpha = .58, color = NA
    ) +
    ggplot2::geom_hline(yintercept = 0, color = "#777777", linewidth = .35) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$original_weight),
      color = accuracy_original_color, linewidth = .68, linetype = "solid"
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$cluster_bootstrap_median),
      color = accuracy_bootstrap_median_color,
      linewidth = .62, linetype = "dashed"
    ) +
    ggplot2::annotate(
      "label", x = max(1, round(edge_count * .025)), y = annotation_y,
      label = sprintf("n = %d edges", edge_count),
      hjust = 0, vjust = 1, size = 1.9,
      color = "#202020", fill = scales::alpha("white", .88),
      linewidth = 0, label.padding = grid::unit(.7, "mm")
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::breaks_pretty(n = 5),
      expand = ggplot2::expansion(mult = c(.01, .01))
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::breaks_pretty(n = 5),
      expand = ggplot2::expansion(mult = c(.03, .04))
    ) +
    ggplot2::labs(
      title = title_text,
      x = "Edges ordered by original weight",
      y = "Edge weight"
    ) +
    base_panel_theme() +
    ggplot2::theme(legend.position = "none")
}

panel_g <- make_accuracy_panel("Temporal", "Temporal Edge\nAccuracy")
panel_h <- make_accuracy_panel(
  "Contemporaneous", "Contemporaneous Edge\nAccuracy"
)
panel_i <- make_accuracy_panel(
  "Between-person", "Between-person Edge\nAccuracy"
)

manuscript_note <- paste(
  paste(
    "Note. Panels A-C show overall one-step bridge expected-influence",
    "stability from 1,000 participant case-drop bootstrap refits estimated",
    "on full unthresholded networks. Temporal results retain separate",
    "incoming and outgoing directions; the undirected networks use signed",
    "overall Bridge EI without community-specific series. Panels D-F show",
    "edge stability from the same case-drop",
    "refits after thresholding each replication by its own raw p results",
    "(125 replications at each dropped proportion). Lines represent mean",
    "correlations with the corresponding original-sample network, and shaded",
    "bands represent central 95% intervals."
  ),
  paste(
    "CS (0.7) is the maximum proportion of participants that can be dropped",
    "while retaining a 95% probability that the correlation with the",
    "original sample is at least .70."
  ),
  paste(
    "Panels G-I show edge accuracy from 1,000 nonparametric",
    "participant-cluster bootstrap refits; gray areas indicate percentile",
    "95% intervals, red solid lines indicate original-sample sparse edge",
    "weights, and black dashed lines indicate bootstrap medians."
  ),
  paste(
    "For panels D-I, every bootstrap replication is thresholded by its own",
    "unadjusted raw p < .05 results before any correlation, median, or",
    "interval is calculated; the two undirected networks use the",
    "conservative two-direction AND rule."
  ),
  "The Temporal edge panels include all 12 modeled autoregressive edges."
)

combined_figure <- patchwork::wrap_plots(
  panel_a, panel_b, panel_c,
  panel_d, panel_e, panel_f,
  panel_g, panel_h, panel_i,
  ncol = 3L, nrow = 3L, byrow = TRUE,
  heights = c(1, 1, 1)
) +
  patchwork::plot_annotation(
    tag_levels = "A",
    tag_suffix = ".",
    theme = ggplot2::theme(
      plot.tag = ggplot2::element_text(
        face = "bold", size = 9.5, color = "#111111"
      ),
      plot.tag.position = c(.008, .992),
      plot.margin = ggplot2::margin(4, 5, 4, 4)
    )
  )

# EDIT THE COMBINED A4 FIGURE WIDTH AND HEIGHT HERE (inches).
combined_width <- 6
combined_height <- 7
ggplot2::ggsave(
  figure_paths[["svg"]], combined_figure,
  device = svglite::svglite,
  width = combined_width, height = combined_height,
  units = "in", bg = "white"
)

panel_specification <- data.frame(
  panel = LETTERS[1:9],
  row = rep(1:3, each = 3),
  column = rep(1:3, times = 3),
  title = c(
    "Temporal Bridge EI Stability",
    "Contemporaneous Bridge EI Stability",
    "Between-person Bridge EI Stability",
    "Temporal Edge Stability",
    "Contemporaneous Edge Stability",
    "Between-person Edge Stability",
    "Temporal Edge Accuracy",
    "Contemporaneous Edge Accuracy",
    "Between-person Edge Accuracy"
  ),
  resampling_source = c(
    rep("Participant case-drop bootstrap", 6L),
    rep("Nonparametric participant-cluster bootstrap", 3L)
  ),
  network_within_replication = c(
    rep("Unthresholded", 3L),
    rep("Replication-specific raw-p sparse", 6L)
  ),
  stringsAsFactors = FALSE
)

results_s6_5_mlvar_stability_accuracy <- list(
  schema_version = "3.1.0",
  analysis_rule = list(
    bridge_stability = paste(
      "Panels A-C use overall signed one-step Bridge EI recomputed from",
      "the full unthresholded S6.3 participant case-drop edge vectors;",
      "A retains incoming/outgoing direction and B-C use one undirected series"
    ),
    edge_stability_and_accuracy = paste(
      "Panels D-I use replication-specific raw-p sparse networks from",
      "the S6.4 nonparametric bootstrap accuracy checkpoint"
    ),
    sparse_temporal = sparse_checkpoint$temporal_support_rule,
    sparse_undirected = sparse_checkpoint$undirected_support_rule,
    sparse_unsupported_edges =
      "Set to exactly zero before every D-I calculation",
    original_sample = paste(
      "Full unthresholded S6.3 edge-derived centrality for A-C;",
      "S5.3 raw-p sparse edges",
      "for D-I"
    )
  ),
  replication_support_summary = replication_support_summary,
  edge_stability_correlations = edge_stability_correlations,
  edge_stability_summary = edge_stability_summary,
  edge_stability_cs = edge_stability_cs,
  centrality_stability_correlations = centrality_stability_correlations,
  centrality_stability_summary = centrality_stability_summary,
  centrality_cs = centrality_cs,
  centrality_definition_audit = centrality_definition_audit,
  zero_variance_correlation_summary = zero_variance_correlation_summary,
  edge_accuracy_plot_data = edge_accuracy_plot_data,
  panel_specification = panel_specification,
  manuscript_note = manuscript_note,
  figure_specification = list(
    layout = "3 x 3, A-I, left to right and top to bottom",
    width_in = combined_width,
    height_in = combined_height,
    orientation = "Compact A4-compatible Word insertion",
    title = NULL,
    subtitle = NULL,
    embedded_caption = NULL,
    stability_y_axis_title = "Average correlation with original sample",
    stability_x_axis_title = "Participants dropped",
    reference_correlations = c(.30, .50, .70),
    cs_threshold = .70,
    cs_quantile = .05,
    bridge_stability_network_basis = "Unthresholded",
    bridge_stability_definition =
      "Overall signed one-step Bridge EI across the three fitted communities",
    edge_stability_accuracy_network_basis =
      "Replication-specific raw-p sparse",
    panel_border_color = panel_border_color,
    panel_border_linewidth = panel_border_linewidth,
    network_colors = network_colors,
    temporal_bridge_colors = temporal_bridge_colors,
    temporal_bridge_legend_position = "inside lower-right",
    temporal_bridge_legend_position_inside = temporal_legend_inside,
    default_cs_annotation_y = default_cs_annotation_y,
    temporal_bridge_cs_annotation_y = temporal_cs_annotation_y,
    edge_accuracy_original_color = accuracy_original_color,
    edge_accuracy_original_linetype = "solid",
    edge_accuracy_bootstrap_median_color = accuracy_bootstrap_median_color,
    edge_accuracy_bootstrap_median_linetype = "dashed",
    edge_accuracy_band_color = "#BDBDBD",
    temporal_accuracy_includes_autoregressive_edges = TRUE
  ),
  provenance = list(
    sparse_checkpoint_input = paste0(
      "Res_2_Data_for_Analysis/", basename(sparse_checkpoint_file)
    ),
    sparse_checkpoint_md5 = sparse_checkpoint_md5,
    case_drop_edge_checkpoint_input = paste0(
      "Res_2_Data_for_Analysis/", basename(case_drop_edge_checkpoint_file)
    ),
    case_drop_edge_checkpoint_md5 = case_drop_edge_checkpoint_md5,
    network_input = paste0("Res_2_Data_for_Analysis/", basename(network_file)),
    network_md5 = network_md5,
    cluster_replications = 1000L,
    case_drop_replications = 1000L,
    generating_script = paste0(
      "Step_6_mlVAR_Reproduce_Centrality/", basename(script_file)
    ),
    recovery_script = paste0(
      "Step_6_mlVAR_Reproduce_Centrality/",
      "S6_4_NonParametricBootstrapAccuracy.R"
    ),
    bridge_stability_source = paste(
      "S6.3 full unthresholded case-drop edge vectors, converted to overall",
      "signed one-step Bridge EI"
    ),
    edge_stability_accuracy_source =
      "S6.4 replication-specific sparse checkpoint",
    no_mlvar_refit_in_visualizer = TRUE,
    r_version = R.version.string,
    package_versions = stats::setNames(
      vapply(required_packages, function(package_name) {
        as.character(utils::packageVersion(package_name))
      }, character(1)),
      required_packages
    )
  )
)

save(results_s6_5_mlvar_stability_accuracy, file = result_file, version = 2)

cat("S6.5 visualization completed with overall unthresholded Bridge EI in A-C and sparse D-I.\n")
cat("Edge Stability CS (0.7):\n")
print(edge_stability_cs, row.names = FALSE)
cat("Bridge-centrality Stability CS (0.7):\n")
print(centrality_cs, row.names = FALSE)
cat("Manuscript figure note:\n", manuscript_note, "\n", sep = "")
cat("Saved result object: ", result_file, "\n", sep = "")
cat("Saved figures:\n")
cat(paste0("- ", unname(figure_paths)), sep = "\n")
cat("\n")
