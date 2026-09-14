# S6.7: Results-writing table for S6.2 node bridge expected influence
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
    ". Install them under R 4.5.2 before rerunning S6.7."
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

script_file <- resolve_script_file("S6_7_NodeBridgeExpectedInfluenceTable.R")
script_dir <- dirname(script_file)
setwd(script_dir)
message("S6.7 working directory: ", normalizePath(".", winslash = "/"))

analysis_dir <- "../Res_2_Data_for_Analysis"
results_dir <- "../Res_3_Results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

s6_2_file <- file.path(
  analysis_dir, "results_s6_2_node_bridge_expected_influence.rda"
)
result_file <- file.path(
  analysis_dir, "results_s6_7_node_bridge_expected_influence_table.rda"
)
table_file <- file.path(
  results_dir, "Table_S6_7_NodeBridgeExpectedInfluence.xlsx"
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

s6_2_results <- load_single_object(
  s6_2_file, "results_s6_2_node_bridge_expected_influence"
)
point_estimates <- s6_2_results$point_estimates
required_columns <- c(
  "panel_id", "target_id", "node", "node_display", "estimate",
  "ci_lower", "ci_upper", "includes_zero"
)
stopifnot(
  nrow(point_estimates) == 48L,
  all(required_columns %in% names(point_estimates)),
  all(is.finite(point_estimates$estimate)),
  all(is.finite(point_estimates$ci_lower)),
  all(is.finite(point_estimates$ci_upper)),
  all(point_estimates$ci_lower <= point_estimates$ci_upper),
  all(point_estimates$includes_zero ==
        (point_estimates$ci_lower <= 0 & point_estimates$ci_upper >= 0))
)

panel_map <- data.frame(
  panel_id = c("temporal_in", "temporal_out", "contemporaneous", "between"),
  sheet_name = c(
    "Temporal_In", "Temporal_Out", "Contemporaneous", "Between-person"
  ),
  stringsAsFactors = FALSE
)
node_order <- paste0("P", 1:6)
target_order <- c("engagement", "distress")

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

table_long <- point_estimates |>
  dplyr::mutate(
    panel_id = factor(.data$panel_id, levels = panel_map$panel_id),
    target_id = factor(.data$target_id, levels = target_order),
    node_display = factor(.data$node_display, levels = node_order),
    estimate_ci = format_estimate_ci(
      .data$estimate, .data$ci_lower, .data$ci_upper
    )
  ) |>
  dplyr::arrange(.data$panel_id, .data$node_display, .data$target_id) |>
  dplyr::mutate(
    panel_id = as.character(.data$panel_id),
    target_id = as.character(.data$target_id),
    node_display = as.character(.data$node_display)
  )

workbook_sheets <- stats::setNames(lapply(panel_map$panel_id, function(current_panel_id) {
  panel_data <- table_long |>
    dplyr::filter(.data$panel_id == .env$current_panel_id) |>
    dplyr::select("node_display", "target_id", "estimate_ci") |>
    tidyr::pivot_wider(names_from = "target_id", values_from = "estimate_ci") |>
    dplyr::arrange(match(.data$node_display, node_order))
  panel_data <- panel_data[, c("node_display", "engagement", "distress")]
  names(panel_data) <- c(
    "Node", "Engagement Community", "Distress Community"
  )
  as.data.frame(panel_data, stringsAsFactors = FALSE)
}), panel_map$sheet_name)

stopifnot(
  nrow(dplyr::distinct(
    table_long, .data$panel_id, .data$node_display, .data$target_id
  )) == 48L,
  identical(vapply(workbook_sheets, nrow, integer(1)), stats::setNames(
    rep(6L, 4L), panel_map$sheet_name
  )),
  all(vapply(workbook_sheets, function(data) !anyNA(data), logical(1)))
)

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
    openxlsx::setRowHeights(workbook, sheet_name, rows = 1L, heights = 32)
    openxlsx::addStyle(
      workbook, sheet_name, header_style, rows = 1L,
      cols = seq_len(columns), gridExpand = TRUE, stack = TRUE
    )
    openxlsx::addStyle(
      workbook, sheet_name, body_style, rows = 2L:rows,
      cols = seq_len(columns), gridExpand = TRUE, stack = TRUE
    )
    openxlsx::setColWidths(
      workbook, sheet_name, cols = 1L:columns, widths = c(12, 30, 30)
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
    "Exact sheet names and order", "Six PSMU nodes per sheet",
    "Both target communities per node", "No missing cells",
    "No Excel formula errors"
  ),
  passed = c(
    identical(reloaded_sheet_names, panel_map$sheet_name),
    all(vapply(reloaded_sheets, nrow, integer(1)) == 6L),
    all(vapply(reloaded_sheets, ncol, integer(1)) == 3L),
    all(vapply(reloaded_sheets, function(data) !anyNA(data), logical(1))),
    !any(vapply(reloaded_sheets, function(data) {
      any(grepl("^#(REF!|VALUE!|NAME\\?|DIV/0!|N/A)", as.matrix(data)))
    }, logical(1)))
  ),
  stringsAsFactors = FALSE
)
stopifnot(all(workbook_audit$passed), file.info(table_file)$size > 0)

results_s6_7_node_bridge_expected_influence_table <- list(
  schema_version = "1.0.0",
  created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  analysis_rule = paste(
    "Direct transcription of all 48 S6.2 raw one-step target-community",
    "Bridge EI estimates and full unthresholded participant-cluster",
    "bootstrap percentile 95% intervals; no centrality is recomputed in S6.7"
  ),
  estimates_long = table_long,
  workbook_sheets = workbook_sheets,
  workbook_audit = workbook_audit,
  provenance = list(
    s6_2_input = paste0("Res_2_Data_for_Analysis/", basename(s6_2_file)),
    s6_2_md5 = unname(tools::md5sum(s6_2_file)),
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
  results_s6_7_node_bridge_expected_influence_table,
  file = result_file, version = 2
)

cat("S6.7 completed successfully.\n")
cat("Saved workbook: ", table_file, "\n", sep = "")
cat("Saved result object: ", result_file, "\n", sep = "")
print(workbook_audit, row.names = FALSE)
