# Author: Kunru Song

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
  if (!length(candidates) && file.exists(expected_basename)) {
    candidates <- normalizePath(expected_basename, winslash = "/", mustWork = TRUE)
  }
  if (!length(candidates)) {
    stop("Cannot locate ", expected_basename, ". Use Rscript or source() with its path.")
  }
  candidates[[length(candidates)]]
}

resolve_s6_3_mode <- function(script_path, args_all, args, override = NULL) {
  script_argument <- grep("^--file=", args_all, value = TRUE)
  invoked_path <- if (length(script_argument) == 1L) {
    sub("^--file=", "", script_argument[[1L]])
  } else ""
  direct <- nzchar(invoked_path) && file.exists(invoked_path) && identical(
    normalizePath(invoked_path, winslash = "/", mustWork = TRUE), script_path
  )
  # A source() caller's command-line options belong to that caller, not S6.3.
  mode <- if (!is.null(override)) override else if (direct && length(args)) args[[1L]] else "run"
  if (!is.character(mode) || length(mode) != 1L || is.na(mode) ||
      !tolower(mode) %in% c("run", "audit", "merge")) {
    stop("S6.3 mode must be one of: run, audit, merge")
  }
  tolower(mode)
}

script_path <- resolve_script_file("S6_3_CaseDropEdgeStabilityBootstrap.R")
mode <- resolve_s6_3_mode(
  script_path, commandArgs(trailingOnly = FALSE), commandArgs(trailingOnly = TRUE),
  override = getOption("s6_3.mode", NULL)
)
project_root <- normalizePath(
  file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE
)
helper_file <- file.path(
  project_root, "Step_5_mlVAR_Network_Analysis", "S5_0_mlVARHelpers.R"
)
source(helper_file, local = TRUE, encoding = "UTF-8")

Sys.setenv(
  LANG = "C", LC_ALL = "C",
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)
try(Sys.setlocale("LC_ALL", "C"), silent = TRUE)
options(warn = 1, mc.cores = 1L)

# The user-approved memory ceiling is six simultaneous mlVAR fits. Each fit is
# single-core internally, so this is six total model workers, not 6 x 6 cores.
parallel_workers <- 6L
internal_mlvar_cores <- 1L
n_replications <- 1000L
expected_drop_proportions <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.75)
expected_replications_per_proportion <- 125L
# S5.6 used 12 internal mlVAR cores; S6.3 is intentionally limited to one
# internal core in each of six simultaneous fits. Numerical equivalence is
# therefore audited by both a tight absolute tolerance and vector correlation.
centrality_absolute_tolerance <- 1e-2
centrality_correlation_threshold <- 0.99999
stability_q05_absolute_tolerance <- 1e-3

required_packages <- c("lme4", "mlVAR")
available <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
if (!all(available)) {
  stop(
    "Missing required R packages: ",
    paste(names(available)[!available], collapse = ", ")
  )
}

input_file <- file.path(
  project_root, "Res_2_Data_for_Analysis", "dat_step5_mlvar.rda"
)
network_file <- file.path(
  project_root, "Res_2_Data_for_Analysis", "results_s5_3_mlvar_networks.rda"
)
s5_6_cluster_file <- file.path(
  project_root, "Res_2_Data_for_Analysis", "checkpoint_s5_6_cluster_bootstrap.rda"
)
s5_6_case_drop_file <- file.path(
  project_root, "Res_2_Data_for_Analysis", "checkpoint_s5_6_case_drop.rda"
)
checkpoint_root <- file.path(
  project_root, "Res_2_Data_for_Analysis", "checkpoints_s6_3_case_drop_edges"
)
success_dir <- file.path(checkpoint_root, "replications")
failure_dir <- file.path(checkpoint_root, "failed_attempts")
invalid_dir <- file.path(checkpoint_root, "invalid_files")
log_dir <- file.path(checkpoint_root, "logs")
lock_dir <- file.path(checkpoint_root, "run.lock")
manifest_file <- file.path(checkpoint_root, "manifest_s6_3_case_drop_edges.rda")
canonical_checkpoint_file <- file.path(
  project_root, "Res_2_Data_for_Analysis", "checkpoint_s6_3_case_drop_edges.rda"
)

stopifnot(all(file.exists(c(
  helper_file, input_file, network_file, s5_6_cluster_file, s5_6_case_drop_file
))))

dat_step5_mlvar <- load_single_object(input_file, "dat_step5_mlvar")
results_s5_3_mlvar_networks <- load_single_object(
  network_file, "results_s5_3_mlvar_networks"
)
checkpoint_s5_6_cluster_bootstrap <- load_single_object(
  s5_6_cluster_file, "checkpoint_s5_6_cluster_bootstrap"
)
checkpoint_s5_6_case_drop <- load_single_object(
  s5_6_case_drop_file, "checkpoint_s5_6_case_drop"
)

node_map <- dat_step5_mlvar$node_map
node_names <- node_map$node
model_data <- dat_step5_mlvar$prompt[c(node_names, "Name_ID", "occasion")]
participant_series <- split(model_data, model_data$Name_ID)
participant_ids <- names(participant_series)

stopifnot(
  length(participant_ids) == 209L,
  nrow(model_data) == 18810L,
  all(vapply(participant_series, nrow, integer(1)) == 90L),
  all(vapply(participant_series, function(x) identical(x$occasion, 1:90), logical(1)))
)

matrix_list <- list(
  temporal = results_s5_3_mlvar_networks$full_matrices$temporal,
  contemporaneous = results_s5_3_mlvar_networks$full_matrices$contemporaneous,
  between = results_s5_3_mlvar_networks$full_matrices$between
)

edge_vectors_by_network <- function(matrices) {
  temporal_indices <- expand.grid(
    from = node_names, to = node_names,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  temporal_values <- matrices$temporal[cbind(
    match(temporal_indices$from, node_names),
    match(temporal_indices$to, node_names)
  )]
  names(temporal_values) <- paste(
    "Temporal", temporal_indices$from, temporal_indices$to, sep = "|"
  )

  undirected_values <- function(matrix, network) {
    pair <- which(upper.tri(matrix), arr.ind = TRUE)
    value <- matrix[pair]
    names(value) <- paste(
      network, rownames(matrix)[pair[, 1L]], colnames(matrix)[pair[, 2L]], sep = "|"
    )
    value
  }

  list(
    Temporal = temporal_values,
    Contemporaneous = undirected_values(
      matrices$contemporaneous, "Contemporaneous"
    ),
    `Between-person` = undirected_values(matrices$between, "Between-person")
  )
}

flatten_edge_vectors <- function(edge_vectors) {
  c(
    edge_vectors[["Temporal"]],
    edge_vectors[["Contemporaneous"]],
    edge_vectors[["Between-person"]]
  )
}

centrality_vector <- function(matrices) {
  community <- stats::setNames(node_map$community, node_names)
  temporal <- matrices$temporal[node_names, node_names, drop = FALSE]
  contemporaneous <- matrices$contemporaneous[node_names, node_names, drop = FALSE]
  between <- matrices$between[node_names, node_names, drop = FALSE]
  diag(temporal) <- 0
  diag(contemporaneous) <- 0
  diag(between) <- 0

  values <- numeric(0)
  for (node in node_names) {
    other <- community != community[[node]]
    values[paste("Temporal", "Outgoing expected influence", node, sep = "|")] <-
      sum(temporal[node, ])
    values[paste("Temporal", "Incoming expected influence", node, sep = "|")] <-
      sum(temporal[, node])
    values[paste(
      "Temporal", "Bridge outgoing expected influence", node, sep = "|"
    )] <- sum(temporal[node, other])
    values[paste(
      "Temporal", "Bridge incoming expected influence", node, sep = "|"
    )] <- sum(temporal[other, node])
  }
  for (network in c("Contemporaneous", "Between-person")) {
    matrix <- if (network == "Contemporaneous") contemporaneous else between
    for (node in node_names) {
      values[paste(network, "Strength", node, sep = "|")] <-
        sum(abs(matrix[node, ]))
    }
  }
  values
}

baseline_edge_vectors <- edge_vectors_by_network(matrix_list)
baseline_edges <- flatten_edge_vectors(baseline_edge_vectors)
baseline_centrality <- centrality_vector(matrix_list)
checkpoint_s5_6_case_drop$baseline_centrality <- select_named_network_values(
  checkpoint_s5_6_case_drop$baseline_centrality, names(baseline_centrality)
)
checkpoint_s5_6_case_drop$centrality_estimates <- select_named_network_values(
  checkpoint_s5_6_case_drop$centrality_estimates, names(baseline_centrality)
)
input_md5 <- unname(tools::md5sum(input_file))
network_md5 <- unname(tools::md5sum(network_file))
s5_6_cluster_md5 <- unname(tools::md5sum(s5_6_cluster_file))
s5_6_case_drop_md5 <- unname(tools::md5sum(s5_6_case_drop_file))

case_drop_schedule <- checkpoint_s5_6_case_drop$status$drop_proportion
retained_schedule <- checkpoint_s5_6_case_drop$status$retained_participants
base_seed <- checkpoint_s5_6_case_drop$base_seed
model_specification <- checkpoint_s5_6_case_drop$specification
upstream_hash_audit <- data.frame(
  source = c("Current S5.1 input", "Current S5.3 network RDA"),
  current_md5 = c(input_md5, network_md5),
  s5_6_cluster_recorded_md5 = c(
    checkpoint_s5_6_cluster_bootstrap$input_md5,
    checkpoint_s5_6_cluster_bootstrap$network_md5
  ),
  s5_6_case_drop_recorded_md5 = c(
    checkpoint_s5_6_case_drop$input_md5,
    checkpoint_s5_6_case_drop$network_md5
  ),
  exact_hash_match = c(
    identical(checkpoint_s5_6_case_drop$input_md5, input_md5) &&
      identical(checkpoint_s5_6_cluster_bootstrap$input_md5, input_md5),
    identical(checkpoint_s5_6_case_drop$network_md5, network_md5) &&
      identical(checkpoint_s5_6_cluster_bootstrap$network_md5, network_md5)
  ),
  stringsAsFactors = FALSE
)

stopifnot(
  checkpoint_s5_6_cluster_bootstrap$replications == n_replications,
  checkpoint_s5_6_case_drop$replications == n_replications,
  checkpoint_s5_6_case_drop$replications_per_proportion ==
    expected_replications_per_proportion,
  identical(
    checkpoint_s5_6_case_drop$drop_proportions,
    expected_drop_proportions
  ),
  identical(checkpoint_s5_6_case_drop$input_md5, input_md5),
  identical(checkpoint_s5_6_cluster_bootstrap$input_md5, input_md5),
  identical(checkpoint_s5_6_case_drop$status$status, rep("success", n_replications)),
  all(table(case_drop_schedule) == expected_replications_per_proportion),
  length(baseline_edges) == 276L,
  length(baseline_edge_vectors[["Temporal"]]) == 144L,
  length(baseline_edge_vectors[["Contemporaneous"]]) == 66L,
  length(baseline_edge_vectors[["Between-person"]]) == 66L,
  identical(names(baseline_edges), names(checkpoint_s5_6_cluster_bootstrap$baseline_edges)),
  isTRUE(all.equal(
    baseline_edges, checkpoint_s5_6_cluster_bootstrap$baseline_edges,
    tolerance = 1e-12, check.attributes = TRUE
  )),
  isTRUE(all.equal(
    baseline_centrality, checkpoint_s5_6_case_drop$baseline_centrality,
    tolerance = 1e-12, check.attributes = TRUE
  )),
  identical(
    names(baseline_centrality),
    colnames(checkpoint_s5_6_case_drop$centrality_estimates)
  ),
  identical(
    checkpoint_s5_6_cluster_bootstrap$specification,
    checkpoint_s5_6_case_drop$specification
  )
)

replication_file <- function(replication) {
  file.path(success_dir, sprintf("rep_%04d_case_drop_edges.rda", replication))
}

failure_file <- function(replication) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  file.path(
    failure_dir,
    sprintf(
      "rep_%04d_failed_%s_pid%d.rda", replication, timestamp, Sys.getpid()
    )
  )
}

save_single_object_atomic <- function(object, object_name, path) {
  environment <- new.env(parent = emptyenv())
  environment[[object_name]] <- object
  temporary_path <- paste0(path, ".pid", Sys.getpid(), ".new")
  if (file.exists(temporary_path)) unlink(temporary_path)
  save(list = object_name, envir = environment, file = temporary_path, version = 2)
  if (!file.rename(temporary_path, path)) {
    copied <- file.copy(temporary_path, path, overwrite = TRUE)
    unlink(temporary_path)
    if (!copied) stop("Could not atomically write: ", path)
  }
  invisible(path)
}

validate_replication_object <- function(object, replication) {
  is.list(object) &&
    object$schema_version %in% c("1.0.0", "1.1.0") &&
    identical(object$status, "success") &&
    identical(object$replication, as.integer(replication)) &&
    identical(object$seed, as.integer(base_seed + replication)) &&
    identical(object$drop_proportion, case_drop_schedule[[replication]]) &&
    identical(object$retained_participants, retained_schedule[[replication]]) &&
    identical(object$input_md5, input_md5) &&
    identical(object$network_md5, network_md5) &&
    identical(object$s5_6_case_drop_md5, s5_6_case_drop_md5) &&
    identical(object$internal_mlvar_cores, internal_mlvar_cores) &&
    length(object$selected_participant_positions) == retained_schedule[[replication]] &&
    !anyDuplicated(object$selected_participant_positions) &&
    all(object$selected_participant_positions %in% seq_along(participant_ids)) &&
    is.list(object$edge_vectors) &&
    identical(names(object$edge_vectors), names(baseline_edge_vectors)) &&
    identical(
      vapply(object$edge_vectors, length, integer(1)),
      c(Temporal = 144L, Contemporaneous = 66L, `Between-person` = 66L)
    ) &&
    identical(names(object$edges), names(baseline_edges)) &&
    length(object$edges) == 276L &&
    all(is.finite(object$edges)) &&
    identical(names(object$centrality), names(baseline_centrality)) &&
    length(object$centrality) == 72L &&
    all(is.finite(object$centrality))
}

load_valid_replication <- function(replication, path = replication_file(replication)) {
  if (!file.exists(path)) return(NULL)
  object <- tryCatch(
    {
      value <- load_single_object(path, "replication_s6_3_case_drop_edges")
      value$centrality <- select_named_network_values(
        value$centrality, names(baseline_centrality)
      )
      value
    },
    error = function(error) NULL
  )
  if (is.null(object) || !validate_replication_object(object, replication)) {
    return(NULL)
  }
  object
}

archive_invalid_file <- function(path) {
  if (!file.exists(path)) return(invisible(FALSE))
  destination <- file.path(
    invalid_dir,
    paste0(
      tools::file_path_sans_ext(basename(path)), "_",
      format(Sys.time(), "%Y%m%d_%H%M%S"), ".rda"
    )
  )
  moved <- file.rename(path, destination)
  if (!moved) {
    copied <- file.copy(path, destination, overwrite = FALSE)
    if (copied) unlink(path)
    moved <- copied
  }
  if (!moved) stop("Could not archive invalid checkpoint: ", path)
  invisible(TRUE)
}

scan_replications <- function(archive_invalid = FALSE) {
  valid <- logical(n_replications)
  for (replication in seq_len(n_replications)) {
    path <- replication_file(replication)
    if (!file.exists(path)) next
    object <- load_valid_replication(replication, path)
    valid[[replication]] <- !is.null(object)
    if (is.null(object) && archive_invalid) archive_invalid_file(path)
  }
  data.frame(
    replication = seq_len(n_replications),
    drop_proportion = case_drop_schedule,
    retained_participants = retained_schedule,
    valid = valid,
    path = vapply(seq_len(n_replications), replication_file, character(1)),
    stringsAsFactors = FALSE
  )
}

partial_centrality_reproduction_audit <- function(scan) {
  replications <- scan$replication[scan$valid]
  if (!length(replications)) {
    return(list(
      replications_checked = 0L,
      maximum_absolute_difference = NA_real_,
      minimum_vector_correlation = NA_real_,
      absolute_tolerance = centrality_absolute_tolerance,
      correlation_threshold = centrality_correlation_threshold,
      passed = TRUE
    ))
  }
  comparison_by_replication <- vapply(replications, function(replication) {
    object <- load_valid_replication(replication)
    accepted <- checkpoint_s5_6_case_drop$centrality_estimates[replication, ]
    c(
      maximum_absolute_difference = max(abs(object$centrality - accepted)),
      vector_correlation = stats::cor(object$centrality, accepted)
    )
  }, numeric(2))
  maximum_by_replication <- comparison_by_replication[
    "maximum_absolute_difference", 
  ]
  correlation_by_replication <- comparison_by_replication[
    "vector_correlation", 
  ]
  list(
    replications_checked = length(replications),
    maximum_absolute_difference = max(maximum_by_replication),
    minimum_vector_correlation = min(correlation_by_replication),
    absolute_tolerance = centrality_absolute_tolerance,
    correlation_threshold = centrality_correlation_threshold,
    passed = all(is.finite(c(maximum_by_replication, correlation_by_replication))) &&
      max(maximum_by_replication) <= centrality_absolute_tolerance &&
      min(correlation_by_replication) >= centrality_correlation_threshold
  )
}

write_manifest <- function(scan, state) {
  partial_reproduction <- partial_centrality_reproduction_audit(scan)
  if (!partial_reproduction$passed) {
    stop(
      "Existing S6.3 files fail the S5.6 centrality reproduction audit; ",
      "maximum absolute difference = ",
      format(partial_reproduction$maximum_absolute_difference, scientific = TRUE),
      "; minimum vector correlation = ",
      format(partial_reproduction$minimum_vector_correlation, digits = 10)
    )
  }
  manifest_s6_3_case_drop_edges <- list(
    schema_version = "1.1.0",
    state = state,
    updated = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    parallel_workers = parallel_workers,
    internal_mlvar_cores = internal_mlvar_cores,
    replications = n_replications,
    successful_files = sum(scan$valid),
    pending_files = sum(!scan$valid),
    balance = aggregate(
      valid ~ drop_proportion + retained_participants,
      data = scan, FUN = sum
    ),
    input_md5 = input_md5,
    network_md5 = network_md5,
    s5_6_cluster_md5 = s5_6_cluster_md5,
    s5_6_case_drop_md5 = s5_6_case_drop_md5,
    upstream_hash_audit = upstream_hash_audit,
    partial_centrality_reproduction_audit = partial_reproduction,
    checkpoint_root = checkpoint_root,
    canonical_checkpoint_file = canonical_checkpoint_file
  )
  save_single_object_atomic(
    manifest_s6_3_case_drop_edges,
    "manifest_s6_3_case_drop_edges", manifest_file
  )
  invisible(manifest_s6_3_case_drop_edges)
}

make_resample <- function(selected_positions) {
  result <- do.call(rbind, lapply(seq_along(selected_positions), function(index) {
    series <- participant_series[[selected_positions[[index]]]]
    series$Name_ID <- index
    series
  }))
  rownames(result) <- NULL
  result
}

fit_resample_values <- function(data) {
  captured_warnings <- character(0)
  started <- Sys.time()
  fit <- withCallingHandlers(
    mlVAR::mlVAR(
      data = data,
      vars = node_names,
      idvar = "Name_ID",
      beepvar = "occasion",
      lags = 1,
      estimator = "lmer",
      temporal = "orthogonal",
      contemporaneous = "orthogonal",
      nCores = internal_mlvar_cores,
      verbose = FALSE,
      scale = TRUE,
      scaleWithin = FALSE,
      na.rm = TRUE,
      full_detrend = FALSE
    ),
    warning = function(warning) {
      captured_warnings <<- unique(c(captured_warnings, conditionMessage(warning)))
      invokeRestart("muffleWarning")
    }
  )
  matrices <- extract_mlvar_matrices(fit)[c("temporal", "contemporaneous", "between")]
  stopifnot(all(vapply(matrices, function(x) {
    identical(dim(x), c(12L, 12L)) && all(is.finite(x))
  }, logical(1))))

  all_models <- c(fit$output$temporal, fit$output$contemporaneous)
  convergence_messages <- vapply(all_models, function(model) {
    collapse_messages(model@optinfo$conv$lme4$messages)
  }, character(1))
  optimizer_codes <- vapply(all_models, function(model) {
    as.integer(model@optinfo$conv$opt %||% 0L)
  }, integer(1))
  singular <- vapply(all_models, lme4::isSingular, logical(1), tol = 1e-4)
  nonboundary_messages <- nzchar(convergence_messages) &
    !grepl("boundary (singular) fit", convergence_messages, fixed = TRUE)

  edge_vectors <- edge_vectors_by_network(matrices)
  result <- list(
    edge_vectors = edge_vectors,
    edges = flatten_edge_vectors(edge_vectors),
    centrality = centrality_vector(matrices),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    optimizer_nonzero = sum(optimizer_codes != 0L),
    singular_models = sum(singular),
    nonboundary_messages = sum(nonboundary_messages),
    warning_text = paste(captured_warnings, collapse = " | "),
    temporal_spectral_radius = matrix_spectral_radius(matrices$temporal)
  )
  rm(fit, all_models, matrices)
  invisible(gc(verbose = FALSE))
  result
}

write_worker_log <- function(text) {
  path <- file.path(log_dir, sprintf("worker_pid_%d.log", Sys.getpid()))
  cat(text, "\n", file = path, append = TRUE)
  invisible(path)
}

fit_one_replication <- function(replication) {
  started <- Sys.time()
  set.seed(base_seed + replication)
  selected_positions <- sample(
    seq_along(participant_ids), retained_schedule[[replication]], replace = FALSE
  )
  data <- make_resample(selected_positions)
  result <- tryCatch(fit_resample_values(data), error = function(error) error)
  finished <- Sys.time()

  if (inherits(result, "error")) {
    failure_s6_3_case_drop_edges <- list(
      schema_version = "1.1.0", status = "failed",
      replication = as.integer(replication),
      seed = as.integer(base_seed + replication),
      drop_proportion = case_drop_schedule[[replication]],
      retained_participants = retained_schedule[[replication]],
      selected_participant_positions = as.integer(selected_positions),
      input_md5 = input_md5, network_md5 = network_md5,
      s5_6_case_drop_md5 = s5_6_case_drop_md5,
      internal_mlvar_cores = internal_mlvar_cores,
      started = format(started, "%Y-%m-%d %H:%M:%S %Z"),
      finished = format(finished, "%Y-%m-%d %H:%M:%S %Z"),
      elapsed_seconds = as.numeric(difftime(finished, started, units = "secs")),
      error_text = conditionMessage(result)
    )
    path <- failure_file(replication)
    save_single_object_atomic(
      failure_s6_3_case_drop_edges,
      "failure_s6_3_case_drop_edges", path
    )
    write_worker_log(sprintf(
      "rep=%04d drop=%.2f retained=%d FAILED elapsed=%.1fs error=%s",
      replication, case_drop_schedule[[replication]],
      retained_schedule[[replication]],
      failure_s6_3_case_drop_edges$elapsed_seconds,
      failure_s6_3_case_drop_edges$error_text
    ))
    rm(data, result)
    invisible(gc(verbose = FALSE))
    return(list(replication = replication, status = "failed", path = path))
  }

  replication_s6_3_case_drop_edges <- list(
    schema_version = "1.1.0", status = "success",
    replication = as.integer(replication),
    seed = as.integer(base_seed + replication),
    drop_proportion = case_drop_schedule[[replication]],
    retained_participants = retained_schedule[[replication]],
    selected_participant_positions = as.integer(selected_positions),
    input_md5 = input_md5, network_md5 = network_md5,
    s5_6_case_drop_md5 = s5_6_case_drop_md5,
    model_specification = model_specification,
    internal_mlvar_cores = internal_mlvar_cores,
    edge_vectors = result$edge_vectors,
    edges = result$edges,
    centrality = result$centrality,
    diagnostics = list(
      elapsed_seconds = result$elapsed_seconds,
      optimizer_nonzero = result$optimizer_nonzero,
      singular_models = result$singular_models,
      nonboundary_messages = result$nonboundary_messages,
      temporal_spectral_radius = result$temporal_spectral_radius,
      warning_text = result$warning_text
    ),
    started = format(started, "%Y-%m-%d %H:%M:%S %Z"),
    finished = format(finished, "%Y-%m-%d %H:%M:%S %Z"),
    r_version = R.version.string,
    package_versions = vapply(
      c("mlVAR", "lme4"), function(package) as.character(packageVersion(package)),
      character(1)
    )
  )
  stopifnot(validate_replication_object(
    replication_s6_3_case_drop_edges, replication
  ))
  path <- replication_file(replication)
  save_single_object_atomic(
    replication_s6_3_case_drop_edges,
    "replication_s6_3_case_drop_edges", path
  )
  stopifnot(!is.null(load_valid_replication(replication, path)))
  write_worker_log(sprintf(
    paste0(
      "rep=%04d drop=%.2f retained=%d SUCCESS elapsed=%.1fs ",
      "optimizer_nonzero=%d singular=%d nonboundary=%d"
    ),
    replication, case_drop_schedule[[replication]],
    retained_schedule[[replication]], result$elapsed_seconds,
    result$optimizer_nonzero, result$singular_models,
    result$nonboundary_messages
  ))
  rm(data, result, replication_s6_3_case_drop_edges)
  invisible(gc(verbose = FALSE))
  list(replication = replication, status = "success", path = path)
}

process_is_running <- function(pid) {
  if (.Platform$OS.type != "windows") return(FALSE)
  output <- tryCatch(
    suppressWarnings(system2(
      "tasklist",
      c("/FI", shQuote(sprintf("PID eq %d", as.integer(pid))), "/NH"),
      stdout = TRUE, stderr = TRUE
    )),
    error = function(error) character(0)
  )
  any(grepl(sprintf("(^|[^0-9])%d([^0-9]|$)", as.integer(pid)), output))
}

acquire_run_lock <- function() {
  if (dir.exists(lock_dir)) {
    lock_file <- file.path(lock_dir, "lock.rda")
    lock <- tryCatch(
      load_single_object(lock_file, "s6_3_run_lock"),
      error = function(error) NULL
    )
    if (!is.null(lock) && process_is_running(lock$pid)) {
      stop("An S6.3 run is already active with PID ", lock$pid)
    }
    stale <- paste0(
      lock_dir, "_stale_", format(Sys.time(), "%Y%m%d_%H%M%S")
    )
    if (!file.rename(lock_dir, stale)) {
      stop("A stale run lock exists and could not be archived: ", lock_dir)
    }
  }
  if (!dir.create(lock_dir, recursive = FALSE, showWarnings = FALSE)) {
    stop("Could not acquire S6.3 run lock")
  }
  s6_3_run_lock <- list(
    pid = Sys.getpid(),
    started = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    script_path = script_path,
    parallel_workers = parallel_workers,
    internal_mlvar_cores = internal_mlvar_cores
  )
  save_single_object_atomic(s6_3_run_lock, "s6_3_run_lock", file.path(lock_dir, "lock.rda"))
  invisible(s6_3_run_lock)
}

release_run_lock <- function() {
  if (dir.exists(lock_dir)) unlink(lock_dir, recursive = TRUE, force = TRUE)
  invisible(TRUE)
}

merge_completed_replications <- function() {
  scan <- scan_replications(archive_invalid = FALSE)
  if (!all(scan$valid)) {
    stop("Cannot merge: ", sum(!scan$valid), " replication files are missing or invalid")
  }

  edge_estimates <- matrix(
    NA_real_, nrow = n_replications, ncol = length(baseline_edges),
    dimnames = list(seq_len(n_replications), names(baseline_edges))
  )
  centrality_reproduction <- matrix(
    NA_real_, nrow = n_replications, ncol = length(baseline_centrality),
    dimnames = list(seq_len(n_replications), names(baseline_centrality))
  )
  status <- data.frame(
    replication = seq_len(n_replications),
    drop_proportion = case_drop_schedule,
    retained_participants = retained_schedule,
    status = "success", elapsed_seconds = NA_real_,
    optimizer_nonzero = NA_integer_, singular_models = NA_integer_,
    nonboundary_messages = NA_integer_, temporal_spectral_radius = NA_real_,
    warning_text = NA_character_, started = NA_character_, finished = NA_character_,
    file = scan$path, file_md5 = NA_character_, stringsAsFactors = FALSE
  )

  for (replication in seq_len(n_replications)) {
    object <- load_valid_replication(replication)
    stopifnot(!is.null(object))
    edge_estimates[replication, ] <- object$edges
    centrality_reproduction[replication, ] <- object$centrality
    status$elapsed_seconds[[replication]] <- object$diagnostics$elapsed_seconds
    status$optimizer_nonzero[[replication]] <- object$diagnostics$optimizer_nonzero
    status$singular_models[[replication]] <- object$diagnostics$singular_models
    status$nonboundary_messages[[replication]] <-
      object$diagnostics$nonboundary_messages
    status$temporal_spectral_radius[[replication]] <-
      object$diagnostics$temporal_spectral_radius
    status$warning_text[[replication]] <- object$diagnostics$warning_text
    status$started[[replication]] <- object$started
    status$finished[[replication]] <- object$finished
  }
  status$file_md5 <- unname(tools::md5sum(status$file))

  accepted_centrality <- checkpoint_s5_6_case_drop$centrality_estimates
  centrality_difference <- centrality_reproduction - accepted_centrality
  max_abs_difference <- max(abs(centrality_difference), na.rm = TRUE)
  reproduction_correlations <- vapply(seq_len(n_replications), function(replication) {
    stats::cor(
      centrality_reproduction[replication, ],
      accepted_centrality[replication, ]
    )
  }, numeric(1))
  min_reproduction_correlation <- min(reproduction_correlations)
  vector_reproduction_passed <-
    is.finite(max_abs_difference) &&
    is.finite(min_reproduction_correlation) &&
    max_abs_difference <= centrality_absolute_tolerance &&
    min_reproduction_correlation >= centrality_correlation_threshold

  centrality_keys <- strsplit(colnames(accepted_centrality), "|", fixed = TRUE)
  metric_group <- vapply(
    centrality_keys,
    function(parts) paste(parts[[1L]], parts[[2L]], sep = "|"),
    character(1)
  )
  metric_names <- unique(metric_group)
  stability_correlations <- function(estimates) {
    result <- matrix(
      NA_real_, nrow = n_replications, ncol = length(metric_names),
      dimnames = list(seq_len(n_replications), metric_names)
    )
    for (metric in metric_names) {
      index <- which(metric_group == metric)
      result[, metric] <- apply(
        estimates[, index, drop = FALSE], 1L, stats::cor,
        y = baseline_centrality[index]
      )
    }
    result
  }
  stability_q05 <- function(correlations) {
    result <- matrix(
      NA_real_, nrow = length(expected_drop_proportions),
      ncol = length(metric_names),
      dimnames = list(as.character(expected_drop_proportions), metric_names)
    )
    for (drop in expected_drop_proportions) {
      rows <- which(case_drop_schedule == drop)
      result[as.character(drop), ] <- apply(
        correlations[rows, , drop = FALSE], 2L,
        stats::quantile, probs = 0.05, names = FALSE, type = 7
      )
    }
    result
  }
  stability_cs <- function(q05) {
    vapply(seq_len(ncol(q05)), function(column) {
      contiguous_pass <- cumprod(q05[, column] >= 0.70) == 1L
      if (!any(contiguous_pass)) 0 else {
        max(expected_drop_proportions[contiguous_pass])
      }
    }, numeric(1), USE.NAMES = FALSE) |>
      stats::setNames(colnames(q05))
  }

  accepted_stability_correlations <- stability_correlations(accepted_centrality)
  reproduced_stability_correlations <- stability_correlations(
    centrality_reproduction
  )
  stopifnot(all(is.finite(c(
    accepted_stability_correlations, reproduced_stability_correlations
  ))))
  accepted_stability_q05 <- stability_q05(accepted_stability_correlations)
  reproduced_stability_q05 <- stability_q05(reproduced_stability_correlations)
  accepted_cs <- stability_cs(accepted_stability_q05)
  reproduced_cs <- stability_cs(reproduced_stability_q05)
  q05_max_abs_difference <- max(abs(
    reproduced_stability_q05 - accepted_stability_q05
  ))
  cs_identical <- identical(accepted_cs, reproduced_cs)
  stability_consequence_passed <-
    cs_identical &&
    is.finite(q05_max_abs_difference) &&
    q05_max_abs_difference <= stability_q05_absolute_tolerance
  reproduction_passed <-
    vector_reproduction_passed && stability_consequence_passed

  if (!reproduction_passed) {
    stop(
      "The S6.3 refits did not reproduce accepted S5.6 centrality within tolerance; ",
      "maximum absolute difference = ", format(max_abs_difference, scientific = TRUE),
      "; minimum vector correlation = ",
      format(min_reproduction_correlation, digits = 10),
      "; CS coefficients identical = ", cs_identical,
      "; maximum q05 difference = ",
      format(q05_max_abs_difference, scientific = TRUE)
    )
  }

  cs_audit <- data.frame(
    metric = metric_names,
    accepted_cs_0_70 = unname(accepted_cs),
    reproduced_cs_0_70 = unname(reproduced_cs),
    cs_identical = unname(accepted_cs == reproduced_cs),
    q05_maximum_absolute_difference = vapply(metric_names, function(metric) {
      max(abs(
        reproduced_stability_q05[, metric] -
          accepted_stability_q05[, metric]
      ))
    }, numeric(1)),
    stringsAsFactors = FALSE
  )

  network_columns <- list(
    Temporal = grep("^Temporal\\|", colnames(edge_estimates)),
    Contemporaneous = grep("^Contemporaneous\\|", colnames(edge_estimates)),
    `Between-person` = grep("^Between-person\\|", colnames(edge_estimates))
  )
  edge_estimates_by_network <- lapply(network_columns, function(index) {
    edge_estimates[, index, drop = FALSE]
  })

  checkpoint_s6_3_case_drop_edges <- list(
    schema_version = "1.2.0",
    type = paste(
      "participant case-dropping bootstrap without replacement;",
      "one durable file per replication"
    ),
    created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    replications = n_replications,
    base_seed = base_seed,
    drop_proportions = expected_drop_proportions,
    replications_per_proportion = expected_replications_per_proportion,
    parallel_workers = parallel_workers,
    internal_mlvar_cores = internal_mlvar_cores,
    specification = model_specification,
    input_md5 = input_md5,
    network_md5 = network_md5,
    s5_6_cluster_md5 = s5_6_cluster_md5,
    s5_6_case_drop_md5 = s5_6_case_drop_md5,
    upstream_hash_audit = upstream_hash_audit,
    baseline_edge_vectors = baseline_edge_vectors,
    baseline_edges = baseline_edges,
    edge_estimates = edge_estimates,
    edge_estimates_by_network = edge_estimates_by_network,
    centrality_reproduction = centrality_reproduction,
    centrality_reproduction_audit = list(
      accepted_source = s5_6_case_drop_file,
      accepted_s5_6_internal_mlvar_cores =
        checkpoint_s5_6_case_drop$worker_cores,
      s6_3_internal_mlvar_cores = internal_mlvar_cores,
      absolute_tolerance = centrality_absolute_tolerance,
      correlation_threshold = centrality_correlation_threshold,
      maximum_absolute_difference = max_abs_difference,
      minimum_vector_correlation = min_reproduction_correlation,
      vector_reproduction_passed = vector_reproduction_passed,
      stability_q05_absolute_tolerance = stability_q05_absolute_tolerance,
      stability_q05_maximum_absolute_difference = q05_max_abs_difference,
      cs_coefficients_identical = cs_identical,
      passed = reproduction_passed
    ),
    centrality_stability_reproduction = list(
      metric_names = metric_names,
      accepted_correlations = accepted_stability_correlations,
      reproduced_correlations = reproduced_stability_correlations,
      accepted_q05 = accepted_stability_q05,
      reproduced_q05 = reproduced_stability_q05,
      cs_audit = cs_audit
    ),
    status = status,
    balance = aggregate(
      replication ~ drop_proportion + retained_participants,
      data = status, FUN = length
    ),
    checkpoint_directory = checkpoint_root,
    package_versions = vapply(
      c("mlVAR", "lme4"), function(package) as.character(packageVersion(package)),
      character(1)
    ),
    r_version = R.version.string,
    validation = list(
      all_1000_replication_files_valid = all(scan$valid),
      every_edge_vector_finite = all(is.finite(edge_estimates)),
      exact_edge_dimensions = identical(dim(edge_estimates), c(1000L, 276L)),
      exact_network_dimensions = identical(
        vapply(edge_estimates_by_network, ncol, integer(1)),
        c(Temporal = 144L, Contemporaneous = 66L, `Between-person` = 66L)
      ),
      six_workers_one_internal_core =
        parallel_workers == 6L && internal_mlvar_cores == 1L,
      no_failed_replications = all(status$status == "success"),
      centrality_reproduction_passed = reproduction_passed
    )
  )
  stopifnot(all(unlist(checkpoint_s6_3_case_drop_edges$validation)))
  save_single_object_atomic(
    checkpoint_s6_3_case_drop_edges,
    "checkpoint_s6_3_case_drop_edges", canonical_checkpoint_file
  )
  reloaded <- load_single_object(
    canonical_checkpoint_file, "checkpoint_s6_3_case_drop_edges"
  )
  stopifnot(
    identical(dim(reloaded$edge_estimates), c(1000L, 276L)),
    all(unlist(reloaded$validation))
  )
  invisible(checkpoint_s6_3_case_drop_edges)
}

run_s6_3 <- function(mode) {
# Acquire the controller lock before any checkpoint scan can move/write files.
# Audit remains read-only, including when output directories do not yet exist.
if (mode != "audit") {
  dir.create(checkpoint_root, recursive = TRUE, showWarnings = FALSE)
  acquire_run_lock()
  on.exit(release_run_lock(), add = TRUE)
  for (path in c(success_dir, failure_dir, invalid_dir, log_dir)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}
initial_scan <- scan_replications(archive_invalid = identical(mode, "run"))
if (mode != "audit") {
  write_manifest(initial_scan, if (all(initial_scan$valid)) "complete" else "incomplete")
}
cat(sprintf(
  paste0(
    "S6.3 audit: valid=%d pending=%d workers=%d internal_mlvar_cores=%d ",
    "mode=%s time=%s\n"
  ),
  sum(initial_scan$valid), sum(!initial_scan$valid),
  parallel_workers, internal_mlvar_cores, mode,
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
))
partial_audit <- partial_centrality_reproduction_audit(initial_scan)
if (partial_audit$replications_checked > 0L) {
  cat(sprintf(
    paste0(
      "S6.3 partial centrality reproduction: n=%d max_abs_diff=%.3e ",
      "min_cor=%.9f abs_tolerance=%.1e cor_threshold=%.6f passed=%s\n"
    ),
    partial_audit$replications_checked,
    partial_audit$maximum_absolute_difference,
    partial_audit$minimum_vector_correlation,
    partial_audit$absolute_tolerance,
    partial_audit$correlation_threshold,
    partial_audit$passed
  ))
}

if (mode == "audit") {
  return(invisible(initial_scan))
}
if (mode == "merge") {
  return(invisible(merge_completed_replications()))
}

pending <- initial_scan$replication[!initial_scan$valid]
if (!length(pending)) {
  merge_completed_replications()
  write_manifest(scan_replications(), "complete")
  cat("S6.3 is already complete; canonical checkpoint was revalidated.\n")
  return(invisible(initial_scan))
}

cat(sprintf(
  "S6.3 run starting %d pending replications with %d single-core workers.\n",
  length(pending), parallel_workers
))
flush.console()

worker_console_file <- file.path(log_dir, "psock_worker_console.log")
cluster <- parallel::makePSOCKcluster(
  parallel_workers, outfile = worker_console_file
)
on.exit({
  if (!is.null(cluster)) try(parallel::stopCluster(cluster), silent = TRUE)
}, add = TRUE, after = FALSE)
parallel::clusterEvalQ(cluster, {
  Sys.setenv(
    LANG = "C", LC_ALL = "C",
    OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )
  try(Sys.setlocale("LC_ALL", "C"), silent = TRUE)
  options(warn = 1, mc.cores = 1L)
  NULL
})
parallel::clusterExport(
  cluster,
  varlist = c(
    "helper_file", "node_names", "node_map", "participant_series",
    "participant_ids", "base_seed", "case_drop_schedule", "retained_schedule",
    "input_md5", "network_md5", "s5_6_case_drop_md5",
    "internal_mlvar_cores", "model_specification", "baseline_edge_vectors",
    "baseline_edges", "baseline_centrality", "success_dir", "failure_dir",
    "log_dir", "edge_vectors_by_network", "flatten_edge_vectors",
    "centrality_vector", "replication_file", "failure_file",
    "save_single_object_atomic", "validate_replication_object",
    "load_valid_replication", "make_resample", "fit_resample_values",
    "write_worker_log", "fit_one_replication"
  ),
  envir = environment()
)
parallel::clusterEvalQ(cluster, {
  source(helper_file, encoding = "UTF-8")
  stopifnot(
    requireNamespace("mlVAR", quietly = TRUE),
    requireNamespace("lme4", quietly = TRUE)
  )
  NULL
})

run_results <- parallel::parLapplyLB(cluster, pending, fit_one_replication)
parallel::stopCluster(cluster)
cluster <- NULL

final_scan <- scan_replications(archive_invalid = TRUE)
write_manifest(
  final_scan,
  if (all(final_scan$valid)) "complete" else "incomplete_after_run"
)
run_status <- vapply(run_results, `[[`, character(1), "status")
cat(sprintf(
  "S6.3 run returned success=%d failed=%d; valid files=%d pending=%d.\n",
  sum(run_status == "success"), sum(run_status == "failed"),
  sum(final_scan$valid), sum(!final_scan$valid)
))

if (!all(final_scan$valid)) {
  stop(
    "S6.3 remains incomplete. Re-run the same script; only missing/invalid ",
    "replications will be fitted."
  )
}

merge_completed_replications()
write_manifest(scan_replications(), "complete")
cat(sprintf(
  "S6.3 COMPLETE: 1000 durable edge-vector files and %s written at %s.\n",
  basename(canonical_checkpoint_file),
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
))
invisible(final_scan)
}

run_s6_3(mode)
