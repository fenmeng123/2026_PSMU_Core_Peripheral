# Author: Kunru Song

locate_project_root <- function(start = getwd()) {
  # Public code folders remain available even when private data/docs are absent.
  markers <- c("Step_1_Descriptive_and_Reliability", "Step_5_mlVAR_Network_Analysis")
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  if (length(current) != 1L) stop("Supply one project location")
  if (!dir.exists(current)) current <- dirname(current)
  repeat {
    if (all(dir.exists(file.path(current, markers)))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  stop("Could not locate the public analysis project from: ", start)
}

load_single_object <- function(path, expected_name) {
  environment <- new.env(parent = emptyenv())
  loaded <- load(path, envir = environment)
  if (!identical(loaded, expected_name)) {
    stop("Expected only ", expected_name, "; found: ", paste(loaded, collapse = ", "))
  }
  environment[[expected_name]]
}

select_named_network_values <- function(values, expected_names) {
  value_names <- if (is.matrix(values)) colnames(values) else names(values)
  if (!is.numeric(values) || is.null(value_names) || anyNA(value_names) ||
      anyDuplicated(value_names) || !length(expected_names) ||
      anyNA(expected_names) || anyDuplicated(expected_names) ||
      !all(expected_names %in% value_names)) {
    stop("Missing, duplicated or invalid named network summary fields")
  }
  # Select the currently specified fields by name, never by column position.
  # This is a lossless view of retained values, not a new estimate or a refit.
  if (is.matrix(values)) values[, expected_names, drop = FALSE] else values[expected_names]
}

retain_case_drop_centrality <- function(checkpoint) {
  parts <- strsplit(names(checkpoint$baseline_centrality), "|", fixed = TRUE)
  if (!length(parts) || any(lengths(parts) != 3L)) {
    stop("Invalid case-drop centrality keys")
  }
  nodes <- unique(vapply(parts, `[[`, character(1), 3L))
  temporal_metrics <- c(
    "Outgoing expected influence", "Incoming expected influence",
    "Bridge outgoing expected influence", "Bridge incoming expected influence"
  )
  keys <- c(
    unlist(lapply(nodes, function(node) paste("Temporal", temporal_metrics, node, sep = "|")),
           use.names = FALSE),
    paste("Contemporaneous", "Strength", nodes, sep = "|"),
    paste("Between-person", "Strength", nodes, sep = "|")
  )
  checkpoint$baseline_centrality <- select_named_network_values(checkpoint$baseline_centrality, keys)
  checkpoint$centrality_estimates <- select_named_network_values(checkpoint$centrality_estimates, keys)
  checkpoint
}

collapse_messages <- function(x) {
  if (is.null(x) || !length(x)) return("")
  paste(unique(unlist(x, recursive = TRUE, use.names = FALSE)), collapse = " | ")
}

max_gradient <- function(model) {
  gradient <- tryCatch(model@optinfo$derivs$gradient, error = function(e) NULL)
  if (is.null(gradient) || !length(gradient) || any(!is.finite(gradient))) return(NA_real_)
  max(abs(gradient))
}

lmer_diagnostic_table <- function(models, node_names, stage) {
  if (length(models) != length(node_names)) {
    stop("Expected ", length(node_names), " models for ", stage, "; found ", length(models))
  }
  dplyr::bind_rows(lapply(seq_along(models), function(index) {
    model <- models[[index]]
    variance <- as.data.frame(lme4::VarCorr(model))
    diagonal_variance <- variance[is.na(variance$var2), , drop = FALSE]
    messages <- collapse_messages(model@optinfo$conv$lme4$messages)
    data.frame(
      stage = stage,
      outcome = node_names[[index]],
      observations = as.integer(stats::nobs(model)),
      groups = as.integer(length(unique(lme4::getME(model, "flist")[[1L]]))),
      fixed_effects = as.integer(length(lme4::fixef(model))),
      random_variance_terms = as.integer(nrow(diagonal_variance)),
      zero_random_variances = as.integer(sum(
        is.finite(diagonal_variance$vcov) & diagonal_variance$vcov < 1e-8
      )),
      singular_1e4 = lme4::isSingular(model, tol = 1e-4),
      singular_1e6 = lme4::isSingular(model, tol = 1e-6),
      optimizer_code = as.integer(model@optinfo$conv$opt %||% 0L),
      max_abs_gradient = max_gradient(model),
      convergence_message = messages,
      aic = stats::AIC(model),
      bic = stats::BIC(model),
      finite_fixed_effects = all(is.finite(lme4::fixef(model))),
      stringsAsFactors = FALSE
    )
  }))
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

extract_mlvar_matrices <- function(fit) {
  node_names <- fit$input$vars
  temporal <- t(fit$results$Beta$mean[, , 1L])
  temporal_se <- t(fit$results$Beta$SE[, , 1L])
  temporal_p <- t(fit$results$Beta$P[, , 1L])
  contemporaneous <- fit$results$Theta$pcor$mean
  contemporaneous_gamma_p <- fit$results$Gamma_Theta$P
  between <- fit$results$Omega_mu$pcor$mean
  between_gamma_p <- fit$results$Gamma_Omega_mu$P

  matrices <- list(
    temporal = temporal,
    temporal_se = temporal_se,
    temporal_p = temporal_p,
    contemporaneous = contemporaneous,
    contemporaneous_gamma_p = contemporaneous_gamma_p,
    between = between,
    between_gamma_p = between_gamma_p
  )
  matrices <- lapply(matrices, function(matrix) {
    if (!identical(dim(matrix), c(length(node_names), length(node_names)))) {
      stop("Unexpected mlVAR matrix dimensions")
    }
    dimnames(matrix) <- list(node_names, node_names)
    matrix
  })

  matrices
}

validate_network_matrix <- function(matrix, node_names, symmetric = FALSE) {
  stopifnot(
    is.matrix(matrix),
    identical(dim(matrix), c(length(node_names), length(node_names))),
    identical(rownames(matrix), node_names),
    identical(colnames(matrix), node_names),
    all(is.finite(matrix))
  )
  if (symmetric && !isTRUE(all.equal(matrix, t(matrix), tolerance = 1e-8))) {
    stop("Expected symmetric network matrix")
  }
  TRUE
}

matrix_spectral_radius <- function(matrix) {
  max(Mod(eigen(matrix, only.values = TRUE)$values))
}

edge_table_directed <- function(weight, se, p, node_map) {
  indices <- expand.grid(
    from_index = seq_len(nrow(weight)),
    to_index = seq_len(ncol(weight)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  indices |>
    dplyr::mutate(
      from = rownames(weight)[from_index],
      to = colnames(weight)[to_index],
      estimate = weight[cbind(from_index, to_index)],
      se = se[cbind(from_index, to_index)],
      p = p[cbind(from_index, to_index)],
      ci_lower = estimate - stats::qnorm(0.975) * se,
      ci_upper = estimate + stats::qnorm(0.975) * se,
      self_loop = from == to
    ) |>
    dplyr::left_join(
      node_map |>
        dplyr::select(from = node, from_label = label,
          from_community = community, from_subgroup = subgroup),
      by = "from"
    ) |>
    dplyr::left_join(
      node_map |>
        dplyr::select(to = node, to_label = label,
          to_community = community, to_subgroup = subgroup),
      by = "to"
    )
}

edge_table_undirected <- function(weight, gamma_p, node_map) {
  pair <- which(upper.tri(weight), arr.ind = TRUE)
  data.frame(
    node1 = rownames(weight)[pair[, 1L]],
    node2 = colnames(weight)[pair[, 2L]],
    estimate = weight[pair],
    p_node1_on_node2 = gamma_p[pair],
    p_node2_on_node1 = t(gamma_p)[pair],
    stringsAsFactors = FALSE
  ) |>
    dplyr::mutate(
      p_and = pmax(p_node1_on_node2, p_node2_on_node1),
      p_or = pmin(p_node1_on_node2, p_node2_on_node1)
    ) |>
    dplyr::left_join(
      node_map |>
        dplyr::select(node1 = node, node1_label = label,
          node1_community = community, node1_subgroup = subgroup),
      by = "node1"
    ) |>
    dplyr::left_join(
      node_map |>
        dplyr::select(node2 = node, node2_label = label,
          node2_community = community, node2_subgroup = subgroup),
      by = "node2"
    )
}

round_numeric_columns <- function(data, digits = 4L) {
  data |>
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, digits)))
}

style_and_sanitize_workbook <- function(path, sheet_data) {
  workbook <- openxlsx::loadWorkbook(path)
  header_style <- openxlsx::createStyle(
    fontName = "Aptos", fontSize = 10, textDecoration = "bold",
    fontColour = "#FFFFFF", fgFill = "#1F4E78",
    halign = "center", valign = "center", wrapText = TRUE,
    border = "bottom", borderColour = "#163A5C"
  )
  body_style <- openxlsx::createStyle(
    fontName = "Aptos", fontSize = 10, valign = "top", wrapText = TRUE
  )
  integer_style <- openxlsx::createStyle(numFmt = "0")
  decimal_style <- openxlsx::createStyle(numFmt = "0.0000")

  for (sheet_name in names(sheet_data)) {
    data <- sheet_data[[sheet_name]]
    row_count <- nrow(data) + 1L
    column_count <- ncol(data)
    openxlsx::showGridLines(workbook, sheet = sheet_name, showGridLines = FALSE)
    openxlsx::freezePane(workbook, sheet = sheet_name, firstRow = TRUE)
    openxlsx::setRowHeights(workbook, sheet = sheet_name, rows = 1L, heights = 30)
    openxlsx::addStyle(
      workbook, sheet = sheet_name, style = header_style,
      rows = 1L, cols = seq_len(column_count), gridExpand = TRUE, stack = TRUE
    )
    if (row_count > 1L) {
      openxlsx::addStyle(
        workbook, sheet = sheet_name, style = body_style,
        rows = 2L:row_count, cols = seq_len(column_count),
        gridExpand = TRUE, stack = TRUE
      )
    }
    integer_columns <- which(vapply(data, is.integer, logical(1)))
    numeric_columns <- setdiff(which(vapply(data, is.numeric, logical(1))), integer_columns)
    if (row_count > 1L && length(integer_columns)) {
      openxlsx::addStyle(
        workbook, sheet = sheet_name, style = integer_style,
        rows = 2L:row_count, cols = integer_columns, gridExpand = TRUE, stack = TRUE
      )
    }
    if (row_count > 1L && length(numeric_columns)) {
      openxlsx::addStyle(
        workbook, sheet = sheet_name, style = decimal_style,
        rows = 2L:row_count, cols = numeric_columns, gridExpand = TRUE, stack = TRUE
      )
    }
    widths <- vapply(seq_len(column_count), function(index) {
      values <- c(names(data)[index], as.character(data[[index]]))
      values[is.na(values)] <- ""
      min(42, max(11, max(nchar(values, type = "width"), na.rm = TRUE) + 2L))
    }, numeric(1))
    openxlsx::setColWidths(
      workbook, sheet = sheet_name, cols = seq_len(column_count), widths = widths
    )
  }

  workbook$worksheets_rels <- lapply(
    workbook$worksheets_rels,
    function(x) x[grepl("printerSettings", x, fixed = TRUE)]
  )
  workbook$Content_Types <- workbook$Content_Types[
    !grepl("/xl/drawings/", workbook$Content_Types, fixed = TRUE)
  ]
  openxlsx::saveWorkbook(workbook, path, overwrite = TRUE)
  invisible(path)
}
