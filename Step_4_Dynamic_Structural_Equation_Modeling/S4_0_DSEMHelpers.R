# Author: Kunru Song

# Generic Step 4 utilities only. Scientific model syntax remains in each owning Rmd.

locate_project_root <- function(start = getwd()) {
  candidate <- normalizePath(start, winslash = "/", mustWork = TRUE)
  if (!dir.exists(candidate)) candidate <- dirname(candidate)
  repeat {
    if (all(dir.exists(file.path(candidate, c(
      "Step_1_Descriptive_and_Reliability", "Step_4_Dynamic_Structural_Equation_Modeling"
    ))))) return(candidate)
    parent <- dirname(candidate)
    if (identical(parent, candidate)) break
    candidate <- parent
  }
  stop("Could not identify the project root from the public Step folders")
}

load_single_object <- function(path, expected_name) {
  environment <- new.env(parent = emptyenv())
  loaded <- load(path, envir = environment)
  if (!identical(loaded, expected_name)) {
    stop("Expected only ", expected_name, "; found: ", paste(loaded, collapse = ", "))
  }
  environment[[expected_name]]
}

mplus_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

summary_value <- function(summary_row, candidates, type = c("numeric", "character")) {
  type <- match.arg(type)
  hit <- intersect(candidates, names(summary_row))
  if (!length(hit)) return(if (type == "numeric") NA_real_ else NA_character_)
  value <- summary_row[[hit[[1L]]]][[1L]]
  if (type == "numeric") mplus_numeric(value) else as.character(value)
}

flatten_messages <- function(x) {
  if (is.null(x) || !length(x)) return("")
  paste(unlist(x, recursive = TRUE, use.names = FALSE), collapse = " | ")
}

parse_output_integer <- function(lines, label) {
  hit <- grep(label, lines, ignore.case = TRUE, value = TRUE, fixed = TRUE)
  if (!length(hit)) return(NA_integer_)
  numbers <- regmatches(hit[[1L]], gregexpr("[0-9]+", hit[[1L]]))[[1L]]
  if (!length(numbers)) NA_integer_ else as.integer(tail(numbers, 1L))
}

parse_tech8_psr <- function(lines) {
  start <- grep("TECHNICAL 8 OUTPUT FOR BAYES ESTIMATION", lines, fixed = TRUE)
  if (!length(start)) return(data.frame())
  section <- lines[seq.int(start[[length(start)]], length(lines))]
  values <- lapply(section, function(line) {
    match <- regexec(
      "^\\s*([0-9]+)\\s+([0-9]+\\.?[0-9]*)\\s+([0-9]+)\\s*$", line,
      perl = TRUE
    )
    captured <- regmatches(line, match)[[1L]]
    if (length(captured) != 4L) return(NULL)
    data.frame(
      iteration = as.integer(captured[[2L]]),
      psr = as.numeric(captured[[3L]]),
      parameter_highest_psr = as.integer(captured[[4L]]),
      stringsAsFactors = FALSE
    )
  })
  values <- values[!vapply(values, is.null, logical(1))]
  if (!length(values)) data.frame() else dplyr::bind_rows(values)
}

numeric_parameter_table <- function(parsed, table = "unstandardized") {
  parameters <- parsed[["parameters"]][[table]]
  if (is.null(parameters) || !nrow(parameters)) {
    stop("No parameters parsed from Mplus table: ", table)
  }
  numeric_names <- intersect(
    c("est", "posterior_sd", "pval", "lower_2.5ci", "upper_2.5ci"),
    names(parameters)
  )
  parameters |>
    dplyr::mutate(dplyr::across(dplyr::all_of(numeric_names), mplus_numeric)) |>
    as.data.frame()
}

dsem_data_matches <- function(path, data, usevariables) {
  if (!file.exists(path) || !all(usevariables %in% names(data)) ||
      anyDuplicated(usevariables)) return(FALSE)
  expected <- as.matrix(data[usevariables])
  if (!is.numeric(expected) || any(!is.finite(expected) & !is.na(expected))) return(FALSE)
  observed <- tryCatch(
    as.matrix(utils::read.table(path, header = FALSE, na.strings = ".",
                               colClasses = "numeric", comment.char = "", quote = "")),
    error = function(e) NULL
  )
  if (is.null(observed) || !identical(dim(observed), dim(expected)) ||
      !identical(is.na(unname(observed)), is.na(unname(expected)))) return(FALSE)
  keep <- !is.na(expected)
  # Allow decimal serialization roundoff only; check each cell and its position.
  all(is.finite(observed[keep])) &&
    all(abs(observed[keep] - expected[keep]) <= 1e-12 * pmax(1, abs(expected[keep])))
}

run_bayes_dsem <- function(
  data, variable_block, model_block, constraint_block, usevariables,
  mplus_directory, stem, title, iterations, seed, expected_clusters,
  mplus_command, expected_mplus_version = "9", reuse_existing = TRUE,
  parameter_scale = c("stdyx", "unstandardized")
) {
  parameter_scale <- match.arg(parameter_scale)
  stopifnot(
    identical(iterations, 80000L), file.exists(mplus_command),
    is.character(expected_mplus_version), length(expected_mplus_version) == 1L,
    nzchar(expected_mplus_version)
  )
  suffix <- "80k"
  run_stem <- paste0(stem, "_", suffix)
  data_path <- file.path(mplus_directory, paste0(run_stem, ".dat"))
  input_path <- file.path(mplus_directory, paste0(run_stem, ".inp"))
  output_path <- file.path(mplus_directory, paste0(run_stem, ".out"))

  analysis_block <- paste(
    "TYPE = TWOLEVEL;",
    "ESTIMATOR = BAYES;",
    "CHAINS = 2;",
    paste0("FBITERATIONS = ", iterations, ";"),
    paste0("BSEED = ", seed, ";"),
    "PROCESSORS = 2;",
    sep = "\n"
  )

  output_block <- if (parameter_scale == "stdyx") {
    paste("TECH1;", "TECH8;", "STANDARDIZED(STDYX);", "CINTERVAL(HPD);", sep = "\n")
  } else {
    paste("TECH1;", "TECH8;", "CINTERVAL(HPD);", sep = "\n")
  }

  object <- MplusAutomation::mplusObject(
    TITLE = paste("! Author: Kunru Song", title, suffix, sep = "\n"),
    VARIABLE = variable_block,
    ANALYSIS = analysis_block,
    MODEL = model_block,
    MODELCONSTRAINT = constraint_block,
    OUTPUT = output_block,
    PLOT = "TYPE = PLOT2;",
    usevariables = usevariables,
    rdata = data,
    autov = FALSE,
    quiet = TRUE
  )

  existing_input <- if (file.exists(input_path)) readLines(input_path, warn = FALSE) else ""
  existing_output <- if (file.exists(output_path)) readLines(output_path, warn = FALSE) else ""
  expected_model_lines <- trimws(c(
    strsplit(variable_block, "\n", fixed = TRUE)[[1L]],
    strsplit(model_block, "\n", fixed = TRUE)[[1L]],
    strsplit(constraint_block, "\n", fixed = TRUE)[[1L]]
  ))
  existing_trimmed <- trimws(existing_input)
  standardized_input_match <- if (parameter_scale == "stdyx") {
    any(grepl(
      "STANDARDIZED\\s*\\(STDYX\\)\\s*;", existing_input,
      ignore.case = TRUE, perl = TRUE
    ))
  } else {
    !any(grepl("STANDARDIZED", existing_input, ignore.case = TRUE))
  }
  standardized_output_match <- if (parameter_scale == "stdyx") {
    any(grepl("STDYX Standardization", existing_output, fixed = TRUE))
  } else {
    TRUE
  }
  reusable <- isTRUE(reuse_existing) && file.exists(data_path) &&
    file.exists(input_path) && file.exists(output_path) &&
    dsem_data_matches(data_path, data, usevariables) &&
    any(trimws(existing_input) == paste0("NAMES = ", paste(usevariables, collapse = " "), ";")) &&
    any(trimws(existing_input) == "MISSING=.;") &&
    any(trimws(existing_input) == paste0("BSEED = ", seed, ";")) &&
    any(trimws(existing_output) == paste("Mplus VERSION", expected_mplus_version)) &&
    any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", existing_output, fixed = TRUE)) &&
    any(grepl(paste0("FBITERATIONS\\s*=\\s*", iterations, ";"),
              existing_input, ignore.case = TRUE, perl = TRUE)) &&
    any(grepl("CHAINS\\s*=\\s*2;", existing_input, ignore.case = TRUE, perl = TRUE)) &&
    standardized_input_match && standardized_output_match &&
    all(expected_model_lines %in% existing_trimmed)

  if (reusable) {
    modeler_result <- NULL
  } else {
    modeler_result <- MplusAutomation::mplusModeler(
      object = object,
      dataout = data_path,
      modelout = input_path,
      run = 1L,
      check = FALSE,
      varwarnings = TRUE,
      Mplus_command = mplus_command,
      writeData = "always",
      hashfilename = FALSE,
      killOnFail = TRUE,
      quiet = TRUE
    )
  }

  stopifnot(file.exists(data_path), file.exists(input_path), file.exists(output_path))
  input_lines <- readLines(input_path, warn = FALSE)
  output_lines <- readLines(output_path, warn = FALSE)
  parsed <- MplusAutomation::readModels(output_path, what = "all", quiet = TRUE)
  unstandardized_parameters <- numeric_parameter_table(parsed, "unstandardized")
  stdyx_parameters <- if (parameter_scale == "stdyx") {
    numeric_parameter_table(parsed, "stdyx.standardized")
  } else {
    data.frame()
  }
  reported_parameters <- if (parameter_scale == "stdyx") {
    stdyx_parameters
  } else {
    unstandardized_parameters
  }
  psr_history <- parse_tech8_psr(output_lines)
  if (!nrow(psr_history)) stop("TECH8 PSR history was not parsed from ", output_path)

  second_half <- psr_history[psr_history$iteration >= iterations / 2, , drop = FALSE]
  final_psr <- psr_history$psr[[nrow(psr_history)]]
  max_second_half_psr <- max(second_half$psr, na.rm = TRUE)
  summary_row <- parsed[["summaries"]][1, , drop = FALSE]
  parser_error_text <- flatten_messages(parsed[["errors"]])
  fbit_advisory <- grepl(
    "USE THE FBITERATIONS OPTION TO INCREASE THE NUMBER OF ITERATIONS",
    parser_error_text, fixed = TRUE
  )
  substantive_parser_error <- nzchar(parser_error_text) &&
    !(fbit_advisory && !grepl(
      "FATAL ERROR|DID NOT TERMINATE|MODEL ESTIMATION DID NOT",
      parser_error_text, ignore.case = TRUE, perl = TRUE
    ))
  explicit_failure <- any(grepl(
    "DID NOT TERMINATE NORMALLY|MODEL ESTIMATION DID NOT TERMINATE|NOT CONVERGED",
    output_lines, ignore.case = TRUE, perl = TRUE
  ))

  checks <- list(
    author_present = any(trimws(input_lines) == "! Author: Kunru Song"),
    twolevel_declared = any(grepl(
      "TYPE\\s*=\\s*TWOLEVEL;", input_lines, ignore.case = TRUE, perl = TRUE
    )),
    bayes_declared = any(grepl(
      "ESTIMATOR\\s*=\\s*BAYES;", input_lines, ignore.case = TRUE, perl = TRUE
    )),
    chains_two_declared = any(grepl(
      "CHAINS\\s*=\\s*2;", input_lines, ignore.case = TRUE, perl = TRUE
    )),
    fixed_iterations_declared = any(grepl(
      paste0("FBITERATIONS\\s*=\\s*", iterations, ";"),
      input_lines, ignore.case = TRUE, perl = TRUE
    )),
    point_mean_absent = !any(grepl("POINT\\s*=\\s*MEAN", input_lines, ignore.case = TRUE)),
    stdyx_declared = if (parameter_scale == "stdyx") any(grepl(
      "STANDARDIZED\\s*\\(STDYX\\)\\s*;", input_lines,
      ignore.case = TRUE, perl = TRUE
    )) else !any(grepl("STANDARDIZED", input_lines, ignore.case = TRUE)),
    stdyx_output_present = if (parameter_scale == "stdyx") any(grepl(
      "STDYX Standardization", output_lines, fixed = TRUE
    )) else TRUE,
    hpd_declared = any(grepl("CINTERVAL\\s*\\(HPD\\)", input_lines, ignore.case = TRUE)),
    random_type_absent = !any(grepl(
      "TYPE\\s*=\\s*TWOLEVEL\\s+RANDOM", input_lines, ignore.case = TRUE, perl = TRUE
    )),
    normal_termination = any(grepl(
      "THE MODEL ESTIMATION TERMINATED NORMALLY", output_lines, fixed = TRUE
    )),
    no_explicit_failure = !explicit_failure,
    no_substantive_parser_error = !substantive_parser_error,
    mplus_version_match = identical(
      summary_value(summary_row, "Mplus.version", "character"),
      expected_mplus_version
    ),
    estimator_bayes = identical(summary_value(summary_row, "Estimator", "character"), "BAYES"),
    chains_two_output = parse_output_integer(
      output_lines, "Number of Markov chain Monte Carlo (MCMC) chains"
    ) == 2L,
    clusters_match = parse_output_integer(output_lines, "Number of clusters") == expected_clusters,
    completed_fixed_iterations = max(psr_history$iteration) == iterations,
    final_psr_le_1_05 = is.finite(final_psr) && final_psr <= 1.05,
    second_half_psr_le_1_05 = is.finite(max_second_half_psr) && max_second_half_psr <= 1.05,
    focal_output_finite = all(is.finite(reported_parameters$est)) &&
      all(is.finite(reported_parameters$posterior_sd)) &&
      all(is.finite(reported_parameters$lower_2.5ci)) &&
      all(is.finite(reported_parameters$upper_2.5ci))
  )
  passed <- all(unlist(checks))

  list(
    iteration_label = suffix,
    requested_iterations = iterations,
    seed = seed,
    expected_mplus_version = expected_mplus_version,
    reused_existing_output = reusable,
    paths = list(data = data_path, input = input_path, output = output_path),
    file_md5 = unname(tools::md5sum(c(data_path, input_path, output_path))),
    syntax = list(
      variable = variable_block, analysis = analysis_block,
      model = model_block, constraint = constraint_block
    ),
    convergence = list(
      passed = passed,
      checks = checks,
      final_psr = final_psr,
      max_second_half_psr = max_second_half_psr,
      psr_history = psr_history
    ),
    parameters = reported_parameters,
    parameters_stdyx = stdyx_parameters,
    parameters_unstandardized = unstandardized_parameters,
    reporting_scale = if (parameter_scale == "stdyx") {
      "Mplus STDYX standardization"
    } else {
      "Mplus unstandardized"
    },
    parsed_summary = summary_row,
    warnings = flatten_messages(parsed[["warnings"]]),
    advisories = if (fbit_advisory) parser_error_text else "",
    errors = if (substantive_parser_error) parser_error_text else ""
  )
}

extract_focal_parameters <- function(parameters, dictionary) {
  rows <- lapply(seq_len(nrow(dictionary)), function(index) {
    definition <- dictionary[index, , drop = FALSE]
    if (definition$type == "path") {
      hit <- parameters[
        toupper(parameters$paramHeader) == toupper(definition$paramHeader) &
          toupper(parameters$param) == toupper(definition$param) &
          toupper(parameters$BetweenWithin) == toupper(definition$BetweenWithin),
        , drop = FALSE
      ]
    } else {
      hit <- parameters[toupper(parameters$param) == toupper(definition$param), , drop = FALSE]
    }
    if (nrow(hit) != 1L) {
      stop("Expected one row for focal parameter ", definition$key, "; found ", nrow(hit))
    }
    data.frame(
      key = definition$key,
      role = definition$role,
      mode = definition$mode,
      outcome = definition$outcome,
      predictor = definition$predictor,
      lag = definition$lag,
      estimate = hit$est,
      posterior_sd = hit$posterior_sd,
      posterior_p_one_tailed = hit$pval,
      lower_hpd = hit$lower_2.5ci,
      upper_hpd = hit$upper_2.5ci,
      excludes_zero = hit$lower_2.5ci > 0 | hit$upper_2.5ci < 0,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

scale_duration_effects <- function(focal, scaling, level) {
  mode_sd <- scaling |>
    dplyr::filter(.data$level == .env$level) |>
    dplyr::select(mode, scaling_variable = variable, scaling_constant = within_sd)
  applicable <- focal |>
    dplyr::filter(!is.na(mode), outcome %in% c("Core", "Peripheral", "Distress")) |>
    dplyr::left_join(mode_sd, by = "mode")
  if (anyNA(applicable$scaling_constant)) stop("Missing duration scaling constant")
  applicable |>
    dplyr::mutate(
      estimate_1sd = estimate * scaling_constant,
      posterior_sd_1sd = posterior_sd * scaling_constant,
      lower_hpd_1sd = lower_hpd * scaling_constant,
      upper_hpd_1sd = upper_hpd * scaling_constant,
      outcome_unit = "Raw 1-5 construct-mean points",
      interpretation_metric = "Effect per 1 pooled within-person SD of log1p duration"
    )
}

convergence_table <- function(run_80k) {
  data.frame(
    iteration_run = run_80k$iteration_label,
    requested_iterations = run_80k$requested_iterations,
    chains = 2L,
    final_psr = run_80k$convergence$final_psr,
    max_second_half_psr = run_80k$convergence$max_second_half_psr,
    normal_termination = run_80k$convergence$checks$normal_termination,
    all_run_checks_passed = run_80k$convergence$passed,
    warnings = run_80k$warnings,
    advisories = run_80k$advisories,
    errors = run_80k$errors,
    canonical = run_80k$convergence$passed,
    stringsAsFactors = FALSE
  )
}

round_numeric_columns <- function(data, digits = 5L) {
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
  decimal_style <- openxlsx::createStyle(numFmt = "0.00000")
  for (sheet_name in names(sheet_data)) {
    data <- sheet_data[[sheet_name]]
    rows <- nrow(data) + 1L
    columns <- ncol(data)
    openxlsx::showGridLines(workbook, sheet_name, showGridLines = FALSE)
    openxlsx::freezePane(workbook, sheet_name, firstRow = TRUE)
    openxlsx::setRowHeights(workbook, sheet_name, rows = 1L, heights = 30)
    openxlsx::addStyle(
      workbook, sheet_name, header_style, rows = 1L, cols = seq_len(columns),
      gridExpand = TRUE, stack = TRUE
    )
    if (rows > 1L) openxlsx::addStyle(
      workbook, sheet_name, body_style, rows = 2L:rows, cols = seq_len(columns),
      gridExpand = TRUE, stack = TRUE
    )
    integer_columns <- which(vapply(data, is.integer, logical(1)))
    numeric_columns <- setdiff(which(vapply(data, is.numeric, logical(1))), integer_columns)
    if (rows > 1L && length(integer_columns)) openxlsx::addStyle(
      workbook, sheet_name, integer_style, rows = 2L:rows, cols = integer_columns,
      gridExpand = TRUE, stack = TRUE
    )
    if (rows > 1L && length(numeric_columns)) openxlsx::addStyle(
      workbook, sheet_name, decimal_style, rows = 2L:rows, cols = numeric_columns,
      gridExpand = TRUE, stack = TRUE
    )
    widths <- vapply(seq_len(columns), function(column) {
      values <- c(names(data)[column], as.character(data[[column]]))
      values[is.na(values)] <- ""
      min(45, max(11, max(nchar(values, type = "width"), na.rm = TRUE) + 2L))
    }, numeric(1))
    openxlsx::setColWidths(workbook, sheet_name, seq_len(columns), widths)
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

validate_workbook <- function(path, expected_sheets, expected_data = NULL) {
  stopifnot(file.exists(path), identical(readxl::excel_sheets(path), expected_sheets))
  formula_errors <- unlist(lapply(expected_sheets, function(sheet) {
    data <- readxl::read_excel(path, sheet = sheet, col_types = "text")
    values <- unlist(data, use.names = FALSE)
    values[grepl("#REF!|#VALUE!|#DIV/0!|#NAME\\?|#N/A", values)]
  }))
  stopifnot(!length(formula_errors))
  if (!is.null(expected_data)) {
    stopifnot(identical(names(expected_data), expected_sheets))
    for (sheet in expected_sheets) {
      observed <- readxl::read_excel(path, sheet = sheet, .name_repair = "minimal")
      expected <- expected_data[[sheet]]
      stopifnot(
        nrow(observed) == nrow(expected),
        ncol(observed) == ncol(expected),
        identical(names(observed), names(expected))
      )
      for (column in seq_along(expected)) {
        expected_column <- expected[[column]]
        observed_column <- observed[[column]]
        if (is.numeric(expected_column)) {
          observed_numeric <- suppressWarnings(as.numeric(observed_column))
          equal <- (is.na(expected_column) & is.na(observed_numeric)) |
            (!is.na(expected_column) & !is.na(observed_numeric) &
               abs(expected_column - observed_numeric) <= 5e-6)
        } else {
          expected_text <- as.character(expected_column)
          observed_text <- as.character(observed_column)
          expected_text[!is.na(expected_text) & expected_text == ""] <- NA_character_
          observed_text[!is.na(observed_text) & observed_text == ""] <- NA_character_
          equal <- (is.na(expected_text) & is.na(observed_text)) |
            (!is.na(expected_text) & !is.na(observed_text) &
               expected_text == observed_text)
        }
        if (!all(equal)) {
          bad_rows <- which(!equal)
          stop(
            "Workbook value mismatch in sheet '", sheet,
            "', column '", names(expected)[[column]],
            "', row(s) ", paste(head(bad_rows + 1L, 10L), collapse = ", "),
            call. = FALSE
          )
        }
      }
    }
  }
  TRUE
}

six_variable_model_block <- function(engagement = c("li", "la", "lp"), controls) {
  stopifnot(length(engagement) == 3L, length(controls) >= 1L)
  e1 <- engagement[[1L]]
  e2 <- engagement[[2L]]
  e3 <- engagement[[3L]]
  controls_text <- paste(controls, collapse = " ")
  paste(
    "%WITHIN%",
    paste0(e1, " ON ", e1, "&1 (arli);"),
    paste0(e1, " ON ", e2, "&1 ", e3, "&1 core&1 peri&1 dist&1;"),
    paste0(e1, " ON ", controls_text, ";"),
    paste0(e2, " ON ", e1, "&1;"),
    paste0(e2, " ON ", e2, "&1 (arla);"),
    paste0(e2, " ON ", e3, "&1 core&1 peri&1 dist&1;"),
    paste0(e2, " ON ", controls_text, ";"),
    paste0(e3, " ON ", e1, "&1 ", e2, "&1;"),
    paste0(e3, " ON ", e3, "&1 (arlp);"),
    paste0(e3, " ON core&1 peri&1 dist&1;"),
    paste0(e3, " ON ", controls_text, ";"),
    paste0("core ON ", e1, "&1 (aic);"),
    paste0("core ON ", e2, "&1 (aac);"),
    paste0("core ON ", e3, "&1 (apc);"),
    "core ON core&1 (arco);",
    "core ON peri&1 dist&1;",
    paste0("core ON ", controls_text, ";"),
    paste0("peri ON ", e1, "&1 (aip);"),
    paste0("peri ON ", e2, "&1 (aap);"),
    paste0("peri ON ", e3, "&1 (app);"),
    "peri ON core&1;",
    "peri ON peri&1 (arpe);",
    "peri ON dist&1;",
    paste0("peri ON ", controls_text, ";"),
    paste0("dist ON ", e1, "&1 (d1i);"),
    paste0("dist ON ", e2, "&1 (d1a);"),
    paste0("dist ON ", e3, "&1 (d1p);"),
    "dist ON core&1 (bc);",
    "dist ON peri&1 (bp);",
    "dist ON dist&1 (ardi);",
    paste0("dist ON ", e1, "&2 (cdi);"),
    paste0("dist ON ", e2, "&2 (cda);"),
    paste0("dist ON ", e3, "&2 (cdp);"),
    paste0("dist ON ", controls_text, ";"),
    paste0(e1, " WITH ", e2, " ", e3, " core peri dist;"),
    paste0(e2, " WITH ", e3, " core peri dist;"),
    paste0(e3, " WITH core peri dist;"),
    "core WITH peri dist;",
    "peri WITH dist;",
    "%BETWEEN%",
    paste0(e1, " WITH ", e2, " ", e3, " core peri dist;"),
    paste0(e2, " WITH ", e3, " core peri dist;"),
    paste0(e3, " WITH core peri dist;"),
    "core WITH peri dist;",
    "peri WITH dist;",
    sep = "\n"
  )
}

assert_free_engagement_lag_paths <- function(
  model_block, engagement = c("li", "la", "lp")
) {
  stopifnot(length(engagement) == 3L)
  e1 <- engagement[[1L]]
  e2 <- engagement[[2L]]
  e3 <- engagement[[3L]]
  model_lines <- trimws(strsplit(model_block, "\n", fixed = TRUE)[[1L]])
  required_lines <- c(
    paste0(e1, " ON ", e1, "&1 (arli);"),
    paste0(e2, " ON ", e1, "&1;"),
    paste0(e2, " ON ", e2, "&1 (arla);"),
    paste0(e3, " ON ", e1, "&1 ", e2, "&1;"),
    paste0(e3, " ON ", e3, "&1 (arlp);")
  )
  forbidden_lines <- c(
    paste0(e2, " ON ", e1, "&1 ", e2, "&1 (arla);"),
    paste0(e3, " ON ", e1, "&1 ", e2, "&1 ", e3, "&1 (arlp);")
  )
  stopifnot(
    all(required_lines %in% model_lines),
    !any(forbidden_lines %in% model_lines)
  )
  invisible(TRUE)
}

six_variable_constraint_block <- function() {
  paste(
    "NEW(iic iip iac iap ipc ipp dic dia dip bdiff h3i h3a h3p);",
    "iic = aic * bc;",
    "iip = aip * bp;",
    "iac = aac * bc;",
    "iap = aap * bp;",
    "ipc = apc * bc;",
    "ipp = app * bp;",
    "dic = iic - iip;",
    "dia = iac - iap;",
    "dip = ipc - ipp;",
    "bdiff = bc - bp;",
    "h3i = aip - aic;",
    "h3a = aap - aac;",
    "h3p = app - apc;",
    sep = "\n"
  )
}

h3_contrast_dictionary <- function(predictor_names) {
  stopifnot(length(predictor_names) == 3L)
  data.frame(
    key = c("h3i", "h3a", "h3p"),
    type = "constraint", paramHeader = NA_character_,
    param = c("H3I", "H3A", "H3P"), BetweenWithin = NA_character_,
    role = "H3 Engagement-to-symptoms Peripheral-minus-Core contrast",
    mode = c("Interactive/message", "Active/posting", "Passive/browsing"),
    outcome = "Peripheral minus Core", predictor = predictor_names,
    lag = "t-1 to t", stringsAsFactors = FALSE
  )
}

six_variable_focal_dictionary <- function(
  engagement = c("li", "la", "lp"),
  predictor_names = c("InteractiveLog", "ActiveLog", "PassiveLog")
) {
  stopifnot(length(engagement) == 3L, length(predictor_names) == 3L)
  e <- toupper(engagement)
  modes <- c("Interactive/message", "Active/posting", "Passive/browsing")
  dictionary <- data.frame(
    key = c(
      "aic", "aac", "apc", "aip", "aap", "app",
      "bc", "bp", "d1i", "d1a", "d1p", "cdi", "cda", "cdp",
      "iic", "iip", "iac", "iap", "ipc", "ipp", "dic", "dia", "dip", "bdiff"
    ),
    type = c(rep("path", 14), rep("constraint", 10)),
    paramHeader = c(
      rep("CORE.ON", 3), rep("PERI.ON", 3), rep("DIST.ON", 8),
      rep(NA_character_, 10)
    ),
    param = c(
      paste0(e, "&1"), paste0(e, "&1"),
      "CORE&1", "PERI&1", paste0(e, "&1"), paste0(e, "&2"),
      "IIC", "IIP", "IAC", "IAP", "IPC", "IPP", "DIC", "DIA", "DIP", "BDIFF"
    ),
    BetweenWithin = c(rep("Within", 14), rep(NA_character_, 10)),
    role = c(
      rep("Engagement to mediator", 6), rep("Mediator to Distress", 2),
      rep("Lag-1 Engagement to Distress adjustment", 3), rep("Lag-2 direct effect", 3),
      rep("Specific indirect effect", 6), rep("Indirect Core-minus-Peripheral contrast", 3),
      "Core-minus-Peripheral component-path contrast"
    ),
    mode = c(
      modes, modes, NA, NA, modes, modes,
      modes[[1L]], modes[[1L]], modes[[2L]], modes[[2L]], modes[[3L]], modes[[3L]],
      modes, NA
    ),
    outcome = c(rep("Core", 3), rep("Peripheral", 3), rep("Distress", 18)),
    predictor = c(
      predictor_names, predictor_names, "Core", "Peripheral", predictor_names,
      predictor_names,
      paste0(predictor_names[[1L]], " via ", c("Core", "Peripheral")),
      paste0(predictor_names[[2L]], " via ", c("Core", "Peripheral")),
      paste0(predictor_names[[3L]], " via ", c("Core", "Peripheral")),
      "Interactive Core minus Peripheral", "Active Core minus Peripheral",
      "Passive Core minus Peripheral", "Core minus Peripheral"
    ),
    lag = c(
      rep("t-1 to t", 11), rep("t-2 to t", 3),
      rep("t-2 through t-1 to t", 9), "t-1 to t"
    ),
    stringsAsFactors = FALSE
  )
  dplyr::bind_rows(dictionary, h3_contrast_dictionary(predictor_names))
}

compare_sensitivity_with_primary <- function(sensitivity, primary, sensitivity_name) {
  sensitivity |>
    dplyr::select(
      key, role, mode, outcome, predictor, lag,
      sensitivity_estimate = estimate, sensitivity_lower_hpd = lower_hpd,
      sensitivity_upper_hpd = upper_hpd,
      sensitivity_excludes_zero = excludes_zero
    ) |>
    dplyr::inner_join(
      primary |>
        dplyr::select(
          key, primary_estimate = estimate, primary_lower_hpd = lower_hpd,
          primary_upper_hpd = upper_hpd, primary_excludes_zero = excludes_zero
        ),
      by = "key"
    ) |>
    dplyr::mutate(
      sensitivity = sensitivity_name,
      direction_consistent = sign(sensitivity_estimate) == sign(primary_estimate),
      zero_classification_consistent = sensitivity_excludes_zero == primary_excludes_zero
    )
}
