# Author: Kunru Song
# Fresh reproduction: source("Reproduce.R") or Rscript Reproduce.R run
# Read-only preflight: Rscript Reproduce.R check
# Definitions only: options(psmu.reproduce.mode = "load"); source("Reproduce.R")

psmu_root <- function() {
  frames <- lapply(sys.frames(), function(frame) frame[["ofile"]])
  frames <- Filter(function(x) is.character(x) && length(x) == 1L &&
    !is.na(x) && basename(x) == "Reproduce.R", frames)
  args <- commandArgs(trailingOnly = FALSE)
  cli <- sub("^--file=", "", args[grepl("^--file=", args)])
  candidates <- c(rev(unlist(frames)), cli[basename(cli) == "Reproduce.R"])
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) stop("Cannot locate Reproduce.R. Invoke its full path.")
  dirname(normalizePath(candidates[1L], winslash = "/", mustWork = TRUE))
}

psmu_manifest <- function() {
  dirs <- c("Step_1_Descriptive_and_Reliability", "Step_2_Measurement_Models",
    "Step_3_Concurrent_Multilevel_Associations", "Step_4_Dynamic_Structural_Equation_Modeling",
    "Step_5_mlVAR_Network_Analysis", "Step_6_mlVAR_Reproduce_Centrality")
  data <- function(x) file.path("Res_2_Data_for_Analysis", x)
  result <- function(x) file.path("Res_3_Results", x)
  entry <- "dat_focal_analysis.rda"
  measure <- "dat_step2_measurement.rda"
  concurrent <- "dat_step3_concurrent.rda"
  dsem <- "dat_step4_dsem.rda"
  mlvar <- "dat_step5_mlvar.rda"
  fit <- "fit_s5_2_primary_mlvar.rda"
  networks <- "results_s5_3_mlvar_networks.rda"
  cluster <- "checkpoint_s5_6_cluster_bootstrap.rda"
  cases <- "checkpoint_s5_6_case_drop.rda"
  edges <- "checkpoint_s6_3_case_drop_edges.rda"
  sparse <- "checkpoint_s6_4_nonparametric_bootstrap_accuracy.rda"
  bridge <- "results_s6_2_node_bridge_expected_influence.rda"
  profiles <- "results_s6_1_normalized_bridge_heatmaps.rda"
  stage <- function(step, id, stem, inputs, rdas = character(), results = character(),
                    extension = "Rmd", args = character(), output_dirs = character()) {
    source <- file.path(dirs[step], paste0(id, "_", stem, ".", extension))
    outputs <- c(if (length(rdas)) data(rdas), if (length(results)) result(results),
      if (extension == "Rmd") sub("[.]Rmd$", ".html", source))
    list(id = id, source = source, inputs = data(inputs), outputs = outputs,
         args = args, output_dirs = output_dirs,
         mplus_dir = if (step %in% 2:4 && !id %in% c("S2_1", "S3_1", "S4_1")) {
           file.path(dirs[step], "mplus", tolower(id))
         } else character())
  }
  stages <- list(
    stage(1, "S1_3", "Reliability", entry, results = "Table_S1_3_Reliability.xlsx"),
    stage(1, "S1_4", "DescriptiveAnalysis", entry, results = "Table_S1_4_DescriptiveAnalysis.xlsx"),
    stage(1, "S1_5", "MissingValueAnalysis", entry, results = "Table_S1_5_MVA_Report.xlsx"),
    stage(2, "S2_1", "Prepare", entry, measure),
    stage(2, "S2_2_1", "BaselineCFA_PSMU18items_WLSMV", measure,
      "fit_s2_2_1_baseline_cfa_wlsmv.rda",
      c("Table_S2_2_1_BaselineCFA_PSMU18items_WLSMV.xlsx",
        "Figure_S2_2_1_BaselineCFA_PSMU18items_WLSMV_Fit.svg",
        "Figure_S2_2_1_BaselineCFA_PSMU18items_WLSMV_Loadings.svg")),
    stage(2, "S2_2_2", "BaselineCFA_PSMU18items_MLR", measure,
      "fit_s2_2_2_baseline_cfa_mlr.rda",
      c("Table_S2_2_2_BaselineCFA_PSMU18items_MLR.xlsx",
        "Figure_S2_2_2_BaselineCFA_PSMU18items_MLR_Fit.svg",
        "Figure_S2_2_2_BaselineCFA_PSMU18items_MLR_Loadings.svg")),
    stage(2, "S2_3", "BaselineCFA_ModelComparisonTest",
      c(measure, "fit_s2_2_1_baseline_cfa_wlsmv.rda")),
    stage(2, "S2_4", "BaselineSEM_PSMU_PsychologicalDistress", measure,
      "fit_s2_8_baseline_sem.rda", c("Table_S2_4_BaselineSEM_PSMU_PsychologicalDistress.xlsx",
        "Figure_S2_4_BaselineSEM_PSMU_PsychologicalDistress.svg")),
    stage(2, "S2_5", "EMA_MultilevelCFA_ModelComparisonTest", measure,
      "fit_s2_9_multilevel_cfa_new_model_order.rda", "Table_S2_5_EMA_MultilevelCFA_ModelComparisonTest.xlsx"),
    stage(2, "S2_6", "EMA_MultilevelCFA_ModelComparisonTest_HeywoodCaseFixed", measure,
      "fit_s2_10_multilevel_cfa_new_model_order_p2var0.rda",
      "Table_S2_6_EMA_MultilevelCFA_ModelComparisonTest_HeywoodCaseFixed.xlsx"),
    stage(3, "S3_1", "PrepareConcurrentData",
      c(entry, "fit_s2_10_multilevel_cfa_new_model_order_p2var0.rda"), concurrent),
    stage(3, "S3_2", "EngagementTimeCorePeripheral", concurrent,
      "fit_s3_2_engagement_core_peripheral.rda",
      c("Table_S3_2_EngagementTimeCorePeripheral.xlsx", "Figure_S3_2_EngagementTimeCorePeripheral.svg")),
    stage(3, "S3_3", "UnadjustedDistressModel", concurrent, "fit_s3_3_unadjusted_distress.rda",
      c("Table_S3_3_UnadjustedDistress.xlsx", "Figure_S3_3_UnadjustedDistress.svg")),
    stage(3, "S3_4", "EngagementAdjustedDistressModel",
      c(concurrent, "fit_s3_3_unadjusted_distress.rda"), "fit_s3_4_adjusted_distress.rda",
      c("Table_S3_4_EngagementAdjustedDistress.xlsx", "Figure_S3_4_EngagementAdjustedDistress.svg")),
    stage(4, "S4_1", "PrepareDSEMData", entry, dsem),
    stage(4, "S4_2", "PrimaryPromptDSEM", dsem, "fit_s4_2_primary_prompt_dsem.rda",
      c("Table_S4_2_PrimaryPromptDSEM.xlsx", "Figure_S4_2_PrimaryPromptDSEM.svg")),
    stage(4, "S4_3", "DailyAggregationDSEM", c(dsem, "fit_s4_2_primary_prompt_dsem.rda"),
      "fit_s4_3_daily_aggregation_dsem.rda",
      c("Table_S4_3_DailyAggregationDSEM.xlsx", "Figure_S4_3_DailyAggregationDSEM.svg")),
    stage(4, "S4_4", "FrequencySensitivityDSEM", c(dsem, "fit_s4_2_primary_prompt_dsem.rda"),
      "fit_s4_4_frequency_sensitivity_dsem.rda",
      c("Table_S4_4_FrequencySensitivityDSEM.xlsx", "Figure_S4_4_FrequencySensitivityDSEM.svg")),
    stage(5, "S5_1", "PreparemlVARData", entry, mlvar),
    stage(5, "S5_2", "PrimarymlVAR", mlvar, fit),
    stage(5, "S5_3", "ExtractmlVARNetworks", c(mlvar, fit), networks, "Table_S5_3_mlVARNetworks.xlsx"),
    stage(5, "S5_4", "VisualizemlVARNetworks", c(mlvar, fit, networks), "results_s5_4_mlvar_figures.rda",
      c("Table_S5_4_mlVAR_AdjacencyMatrices.xlsx", paste0("Figure_S5_4_mlVAR_",
        c("Temporal", "Temporal_WithAutoregression", "Contemporaneous", "BetweenPerson",
          "ThreeNetworks", "Temporal_Heatmap", "Contemporaneous_Heatmap", "BetweenPerson_Heatmap"), ".svg"))),
    stage(5, "S5_5", "CharacterizemlVARNetworks", c(mlvar, fit, networks),
      "results_s5_5_mlvar_characteristics.rda", c("Table_S5_5_mlVARCharacteristics.xlsx",
        paste0("Figure_S5_5_mlVAR_", c("CentralityFull", "CentralitySparse", "Predictability", "BridgeProfiles"), ".svg"))),
    stage(5, "S5_6", "mlVARBootstrapWorker", c(mlvar, networks), c(cluster, cases),
      extension = "R", args = c("all", "1", "1")),
    stage(5, "S5_7", "mlVAR_full_detrend", c(mlvar, fit),
      c("fit_s5_7_mlvar_full_detrend.rda", "results_s5_7_mlvar_full_detrend_comparison.rda"),
      "Table_S5_7_mlVAR_FullDetrendEdgeComparison.xlsx"),
    stage(6, "S6_1", "NormalizedBridgeProfileHeatmaps", networks, profiles,
      paste0("Figure_S6_1_NormalizedBridgeProfiles_", c("Temporal", "Contemporaneous", "BetweenPerson", "Combined"), ".svg"),
      extension = "R"),
    stage(6, "S6_2", "NodeBridgeExpectedInfluence", c(networks, cluster), bridge,
      paste0("Figure_S6_2_NodeBridgeExpectedInfluence_", c("Temporal", "Undirected"), ".svg"), extension = "R"),
    stage(6, "S6_3", "CaseDropEdgeStabilityBootstrap", c(mlvar, networks, cluster, cases), edges,
      extension = "R", args = "run",
      output_dirs = data("checkpoints_s6_3_case_drop_edges/replications")),
    stage(6, "S6_4", "NonParametricBootstrapAccuracy", c(mlvar, networks, cluster, cases, edges), sparse,
      extension = "R", args = "run", output_dirs = data(c(
        "checkpoints_s6_4_nonparametric_bootstrap_accuracy/cluster/replications",
        "checkpoints_s6_4_nonparametric_bootstrap_accuracy/case_drop/replications"))),
    stage(6, "S6_5", "mlVARStabilityAccuracyVisualization", c(networks, edges, sparse),
      "results_s6_5_mlvar_stability_accuracy.rda", "Figure_S6_5_mlVAR_StabilityAccuracy_A4.svg", extension = "R"),
    stage(6, "S6_6", "BridgeProfileBootstrapIntervals", c(networks, profiles, cluster, sparse),
      "results_s6_6_bridge_profile_bootstrap_intervals.rda", "Table_S6_6_BridgeProfileBootstrapIntervals.xlsx", extension = "R"),
    stage(6, "S6_7", "NodeBridgeExpectedInfluenceTable", bridge,
      "results_s6_7_node_bridge_expected_influence_table.rda", "Table_S6_7_NodeBridgeExpectedInfluence.xlsx", extension = "R"),
    stage(6, "S6_8", "APrioriBootstrapDifferenceTests", c(networks, bridge, cluster),
      "results_s6_8_apriori_bootstrap_difference_tests.rda", "Table_S6_8_APrioriBootstrapDifferenceTests.xlsx")
  )
  names(stages) <- vapply(stages, `[[`, character(1), "id")
  # Rmd render executes S6.8's R implementation once. Helpers are definitions,
  # not independent analysis stages. The shard merger is unnecessary with one shard.
  sources <- c(vapply(stages, `[[`, character(1), "source"),
    file.path(dirs[4], "S4_0_DSEMHelpers.R"), file.path(dirs[5], "S5_0_mlVARHelpers.R"),
    file.path(dirs[6], "S6_8_APrioriBootstrapDifferenceTests.R"))
  list(entry = data(entry), stages = stages, sources = sources)
}

psmu_parse_source <- function(path) {
  if (endsWith(path, ".R")) return(parse(path, keep.source = FALSE, encoding = "UTF-8"))
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  starts <- grep("^```\\{r[ ,}]", lines)
  expressions <- lapply(starts, function(i) {
    end <- which(seq_along(lines) > i & grepl("^```[ ]*$", lines))[1L]
    if (is.na(end)) stop("Unclosed R chunk in ", basename(path))
    parse(text = if (end == i + 1L) character() else lines[seq.int(i + 1L, end - 1L)], keep.source = FALSE)
  })
  as.expression(unlist(lapply(expressions, as.list), recursive = FALSE))
}

psmu_packages <- function(expressions) {
  packages <- c("processx", "rmarkdown", "knitr")
  walk <- function(x) {
    if (missing(x)) return(invisible(NULL))
    if (is.call(x)) {
      head <- x[[1L]]
      if (is.symbol(head) && as.character(head) %in% c("::", ":::")) {
        packages <<- c(packages, as.character(x[[2L]]))
      }
      if (is.symbol(head) && as.character(head) %in% c("library", "require", "requireNamespace") && length(x) >= 2L) {
        if (is.character(x[[2L]])) packages <<- c(packages, x[[2L]])
        if (as.character(head) %in% c("library", "require") && is.symbol(x[[2L]]) &&
            !"character.only" %in% names(as.list(x))) packages <<- c(packages, as.character(x[[2L]]))
      }
      if (identical(head, as.name("<-")) && identical(x[[2L]], as.name("required_packages")) &&
          is.call(x[[3L]]) && identical(x[[3L]][[1L]], as.name("c"))) {
        values <- as.list(x[[3L]])[-1L]
        packages <<- c(packages, unlist(Filter(is.character, values), use.names = FALSE))
      }
      # Defaults may contain missing arguments; package calls occur in bodies.
      for (i in seq_along(x)) if (!is.pairlist(x[[i]])) walk(x[[i]])
    } else if (is.expression(x)) for (item in as.list(x)) walk(item)
  }
  for (expression in expressions) walk(expression)
  sort(unique(packages))
}

psmu_validate_graph <- function(manifest) {
  produced <- manifest$entry
  for (stage in manifest$stages) {
    missing <- setdiff(stage$inputs, produced)
    if (length(missing)) stop(stage$id, " has unresolved upstream inputs: ", paste(missing, collapse = ", "))
    if (any(stage$outputs %in% produced)) stop("Duplicate output owner: ", stage$id)
    produced <- c(produced, stage$outputs)
  }
  paths <- c(manifest$entry, manifest$sources, produced,
    unlist(lapply(manifest$stages, `[[`, "output_dirs"), use.names = FALSE),
    unlist(lapply(manifest$stages, `[[`, "mplus_dir"), use.names = FALSE))
  if (any(grepl("(^/|^[A-Za-z]:|(^|/)\\.\\.(/|$)|\\\\)", paths))) stop("Manifest paths must be project-relative.")
  invisible(TRUE)
}

psmu_preflight <- function(root, manifest = psmu_manifest()) {
  psmu_validate_graph(manifest)
  missing_sources <- manifest$sources[!file.exists(file.path(root, manifest$sources))]
  if (length(missing_sources)) stop("Missing public source files: ", paste(missing_sources, collapse = ", "))
  expressions <- lapply(file.path(root, manifest$sources), psmu_parse_source)
  packages <- psmu_packages(expressions)
  available <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  missing_packages <- packages[!available]
  problems <- character()
  if (length(missing_packages)) problems <- c(problems,
    paste("Install required R packages:", paste(missing_packages, collapse = ", ")))
  if (as.character(getRversion()) != "4.5.2") problems <- c(problems, "Use R 4.5.2 for the validated environment.")
  versions <- c(mlVAR = "0.6.1", qgraph = "1.9.8")
  for (package in names(versions)) {
    if (package %in% packages[available] && as.character(utils::packageVersion(package)) != versions[[package]]) {
      problems <- c(problems, paste("Validated version required:", package, versions[[package]]))
    }
  }
  entry_path <- file.path(root, manifest$entry)
  if (!file.exists(entry_path)) {
    problems <- c(problems, paste("Place the author-provided prepared entry at", manifest$entry))
  } else {
    env <- new.env(parent = emptyenv())
    objects <- load(entry_path, envir = env)
    required <- c("Data_ID", "Name_ID", "Day", "Time", "AgeInYears", "Gender", "Grade",
      "FamilyResidence", "ParentsMaritalStatus", "ParentsHighestEdu", "TotalFamilyIncome", "SiblingStructure",
      paste0("V_BSMAS_", 1:18), paste0("K_CESD_", 1:7), paste0("L_DASS_", 1:7), paste0("M_PSS_", 1:10),
      paste0("mu_time", 1:3), paste0("mu_frequency", 1:3), paste0("PMU", 1:6), paste0("DASS", 1:6))
    valid <- identical(objects, "dat_focal_analysis") && is.data.frame(env$dat_focal_analysis) &&
      setequal(names(env$dat_focal_analysis), required) && ncol(env$dat_focal_analysis) == length(required) &&
      !anyDuplicated(names(env$dat_focal_analysis)) && nrow(env$dat_focal_analysis) == 18810L &&
      length(unique(env$dat_focal_analysis$Name_ID)) == 209L
    if (valid) {
      keys <- env$dat_focal_analysis[c("Name_ID", "Day", "Time")]
      valid <- !anyNA(keys) && !anyDuplicated(keys) &&
        all(keys$Day %in% 1:30) && all(keys$Time %in% 1:3) && all(table(keys$Name_ID) == 90L)
    }
    if (!valid) problems <- c(problems, "Prepared-entry object/schema/dimensions do not match the analysis contract.")
  }
  mplus <- Sys.getenv("MPLUS_COMMAND", unset = "")
  if (!nzchar(mplus)) mplus <- unname(Sys.which("Mplus"))
  if (!nzchar(mplus) || !file.exists(mplus)) problems <- c(problems,
    "Set MPLUS_COMMAND to your licensed Mplus 8.7 executable. Preflight does not run Mplus.")
  if (nzchar(mplus) && file.exists(mplus)) mplus <- normalizePath(mplus, winslash = "/", mustWork = TRUE)
  pandoc <- if ("rmarkdown" %in% packages[available]) rmarkdown::find_pandoc()$dir else ""
  if (!nzchar(pandoc)) problems <- c(problems, "Install Pandoc or set RSTUDIO_PANDOC to its executable directory.")
  runtime <- file.path(root, "Res_2_Data_for_Analysis", ".reproduction")
  temporary <- Sys.getenv("PSMU_TMPDIR", unset = file.path(runtime, "tmp"))
  if (grepl("[^!-~]", temporary)) problems <- c(problems,
    "Set PSMU_TMPDIR to an ASCII-only, space-free temporary directory for Windows R/Excel/ZIP compatibility.")
  if (!grepl("^([A-Za-z]:[/\\\\]|/)", temporary)) problems <- c(problems,
    "PSMU_TMPDIR must be an absolute path.")
  if (file.exists(temporary) && !dir.exists(temporary)) problems <- c(problems,
    "PSMU_TMPDIR points to a file, not a directory.")
  cat(length(manifest$stages), "ordered stages;", length(manifest$sources), "sources parsed;",
      length(packages), "R package requirements. No analysis executed.\n")
  if (length(problems)) cat(paste0("- ", problems, collapse = "\n"), "\n")
  list(ok = !length(problems), problems = problems, root = root, manifest = manifest,
       packages = packages, mplus = mplus, pandoc = pandoc, runtime = runtime, temporary = temporary)
}

psmu_hashes <- function(root, paths) {
  paths <- sort(unique(paths))
  absolute <- file.path(root, paths)
  if (any(!file.exists(absolute))) stop("Missing required files: ", paste(paths[!file.exists(absolute)], collapse = ", "))
  hashes <- unname(tools::md5sum(absolute))
  if (anyNA(hashes)) stop("Could not fingerprint required files.")
  stats::setNames(hashes, paths)
}

psmu_quote_library_paths <- function(paths) {
  # R expands library paths as globs at startup. Treat bracketed paths literally.
  vapply(strsplit(paths, "", fixed = TRUE), function(chars) {
    replacements <- c("[" = "[[]", "]" = "[]]", "*" = "[*]", "?" = "[?]")
    special <- chars %in% names(replacements)
    chars[special] <- unname(replacements[chars[special]])
    paste0(chars, collapse = "")
  }, character(1), USE.NAMES = FALSE)
}

psmu_receipt <- function(directory) {
  files <- list.files(directory, pattern = "^[0-9]{6}[.]rds$", full.names = TRUE)
  if (!length(files)) return(NULL)
  readRDS(utils::tail(sort(files), 1L))
}

psmu_record <- function(value, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  number <- length(list.files(directory, pattern = "^[0-9]{6}[.]rds$")) + 1L
  path <- file.path(directory, sprintf("%06d.rds", number))
  if (file.exists(path)) stop("Receipt collision; no existing receipt was replaced.")
  temporary <- tempfile("receipt-", tmpdir = directory)
  saveRDS(value, temporary, version = 2)
  if (!file.rename(temporary, path)) stop("Could not install stage receipt.")
  invisible(path)
}

psmu_stage_files <- function(root, stage) {
  files <- stage$outputs
  for (directory in stage$output_dirs) {
    replicas <- list.files(file.path(root, directory), pattern = "[.]rda$", full.names = FALSE)
    if (length(replicas) != 1000L) stop(stage$id, " requires 1,000 completed files in ", directory)
    files <- c(files, file.path(directory, replicas))
  }
  for (directory in stage$mplus_dir) {
    audits <- list.files(file.path(root, directory), pattern = "[.](inp|out|dat|gh5)$",
                        full.names = FALSE, recursive = TRUE)
    if (!any(endsWith(audits, ".inp")) || !any(endsWith(audits, ".out"))) {
      stop(stage$id, " is missing its Mplus syntax/output files.")
    }
    files <- c(files, file.path(directory, audits))
  }
  files
}

psmu_validate_bootstrap <- function(root, stage) {
  if (stage$id != "S5_6") return(invisible(TRUE))
  for (relative in stage$outputs) {
    e <- new.env(parent = emptyenv())
    name <- sub("[.]rda$", "", basename(relative))
    stopifnot(identical(load(file.path(root, relative), envir = e), name))
    checkpoint <- e[[name]]
    if (nrow(checkpoint$status) != 1000L || !all(checkpoint$status$status == "success")) {
      stop("S5.6 is incomplete or contains failed fits. Downstream stages were not started.")
    }
  }
  invisible(TRUE)
}

psmu_execute <- function(stage, config) {
  script <- file.path(config$root, stage$source)
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  arguments <- if (endsWith(script, ".Rmd")) {
    c("--vanilla", file.path(config$root, "Reproduce.R"), "--render", stage$source)
  } else c("--vanilla", script, stage$args)
  child_env <- c(R_LIBS = paste(psmu_quote_library_paths(.libPaths()), collapse = .Platform$path.sep),
    MPLUS_COMMAND = config$mplus, RSTUDIO_PANDOC = config$pandoc,
    TMPDIR = config$temporary, TMP = config$temporary, TEMP = config$temporary)
  processx::run(rscript, arguments, wd = config$root, env = c("current", child_env),
    stdout = file.path(config$runtime, "logs", paste0(stage$id, ".stdout.log")),
    stderr = file.path(config$runtime, "logs", paste0(stage$id, ".stderr.log")),
    error_on_status = FALSE, echo = FALSE, windows_hide_window = TRUE)
}

psmu_run <- function(config, executor = psmu_execute) {
  if (!isTRUE(config$ok)) stop("Preflight failed. No stages started.")
  root <- config$root
  manifest <- config$manifest
  signature <- list(sources = psmu_hashes(root, c("Reproduce.R", manifest$sources)),
    entry = psmu_hashes(root, manifest$entry), R = as.character(getRversion()),
    packages = stats::setNames(vapply(config$packages, function(p) as.character(utils::packageVersion(p)), character(1)), config$packages),
    mplus = unname(tools::md5sum(config$mplus)))
  # Fail before the first stage if any accepted local results lack our receipts.
  for (stage in manifest$stages) {
    receipt_dir <- file.path(config$runtime, "receipts", stage$id)
    record <- psmu_receipt(receipt_dir)
    local_outputs <- stage$outputs[grepl("^Res_[23]_", stage$outputs)]
    if (is.null(record) && (any(file.exists(file.path(root, local_outputs))) ||
                           any(dir.exists(file.path(root, stage$output_dirs))))) {
      stop("Unmanaged results for ", stage$id, ". Use a fresh checkout plus only the prepared entry. Nothing was rerun.")
    }
    if (!is.null(record) && !identical(record$signature, signature)) {
      stop("Source, data or environment changed since ", stage$id, ". Use a separate fresh reproduction directory.")
    }
  }
  dir.create(config$runtime, recursive = TRUE, showWarnings = FALSE)
  lock <- file.path(config$runtime, "run.lock")
  if (!dir.create(lock, showWarnings = FALSE)) stop("Another run or a stale run.lock exists; inspect it before resuming.")
  on.exit({
    if (identical(normalizePath(lock, winslash = "/", mustWork = FALSE),
                  paste0(normalizePath(config$runtime, winslash = "/", mustWork = TRUE), "/run.lock"))) {
      unlink(lock, recursive = TRUE)
    }
  }, add = TRUE)
  saveRDS(list(pid = Sys.getpid(), started = Sys.time()), file.path(lock, "owner.rds"))
  dir.create(config$temporary, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(config$runtime, "logs"), showWarnings = FALSE)
  for (stage in manifest$stages) {
    receipt_dir <- file.path(config$runtime, "receipts", stage$id)
    record <- psmu_receipt(receipt_dir)
    inputs <- psmu_hashes(root, stage$inputs)
    if (!is.null(record) && !identical(record$inputs, inputs)) stop("Upstream inputs changed for ", stage$id)
    if (!is.null(record) && identical(record$status, "complete")) {
      if (!identical(record$outputs, psmu_hashes(root, psmu_stage_files(root, stage)))) {
        stop("Completed output changed or is missing for ", stage$id, "; it was not regenerated.")
      }
      cat(stage$id, "unchanged; skipped.\n")
      next
    }
    record <- list(signature = signature, inputs = inputs, status = "running", started = Sys.time())
    psmu_record(record, receipt_dir)
    cat(stage$id, "started. Logs are in Res_2_Data_for_Analysis/.reproduction/logs/.\n")
    result <- executor(stage, config)
    if (!identical(as.integer(result$status), 0L)) stop(stage$id, " failed; inspect its log before resuming.")
    psmu_validate_bootstrap(root, stage)
    record$outputs <- psmu_hashes(root, psmu_stage_files(root, stage))
    record$status <- "complete"
    record$completed <- Sys.time()
    psmu_record(record, receipt_dir)
    cat(stage$id, "completed.\n")
  }
  invisible(TRUE)
}

psmu_main <- function(root, mode, args = character()) {
  if (mode == "load") return(invisible(NULL))
  library <- Sys.getenv("PSMU_R_LIBRARY", unset = "")
  if (nzchar(library)) {
    previous <- .libPaths()
    .libPaths(psmu_quote_library_paths(c(normalizePath(library, winslash = "/", mustWork = TRUE), previous)))
    on.exit(.libPaths(psmu_quote_library_paths(previous)), add = TRUE)
  }
  if (mode == "--render") {
    manifest <- psmu_manifest()
    reports <- manifest$sources[endsWith(manifest$sources, ".Rmd")]
    if (length(args) != 1L || !args %in% reports) stop("Unknown report requested.")
    rmarkdown::render(file.path(root, args), envir = new.env(parent = globalenv()), quiet = TRUE)
    return(invisible(TRUE))
  }
  if (!mode %in% c("check", "run")) stop("Use check or run.")
  config <- psmu_preflight(root)
  if (mode == "run") return(psmu_run(config))
  if (!config$ok) stop("Preflight failed; no stages started.")
  invisible(config)
}

.psmu_arguments <- commandArgs(trailingOnly = TRUE)
.psmu_direct <- any(grepl("^--file=.*Reproduce[.]R$", commandArgs(trailingOnly = FALSE)))
.psmu_mode <- if (.psmu_direct && length(.psmu_arguments)) .psmu_arguments[1L] else getOption("psmu.reproduce.mode", "run")
psmu_main(psmu_root(), .psmu_mode, if (.psmu_direct) .psmu_arguments[-1L] else character())
