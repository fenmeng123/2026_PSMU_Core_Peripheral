# S6.4: Nonparametric bootstrap edge-accuracy data preparation
# Author: Kunru Song

options(stringsAsFactors = FALSE, scipen = 999, warn = 1, mc.cores = 1L)
Sys.setenv(
  LANG = "C", LC_ALL = "C",
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)
try(Sys.setlocale("LC_ALL", "C"), silent = TRUE)

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
    stop(
      "Cannot locate ", expected_basename, ". Run it with Rscript or call ",
      "source() with the script path."
    )
  }
  candidates[[length(candidates)]]
}

resolve_s6_4_mode <- function(script_path, args_all, args, override = NULL) {
  script_argument <- grep("^--file=", args_all, value = TRUE)
  invoked_path <- if (length(script_argument) == 1L) {
    sub("^--file=", "", script_argument[[1L]])
  } else ""
  direct <- nzchar(invoked_path) && file.exists(invoked_path) && identical(
    normalizePath(invoked_path, winslash = "/", mustWork = TRUE), script_path
  )
  # A source() caller's command-line options belong to that caller, not S6.4.
  mode <- if (!is.null(override)) override else if (direct && length(args)) args[[1L]] else "run"
  if (!is.character(mode) || length(mode) != 1L || is.na(mode) ||
      !tolower(mode) %in% c("run", "audit", "merge")) {
    stop("S6.4 mode must be one of: run, audit, merge")
  }
  tolower(mode)
}

script_file <- resolve_script_file("S6_4_NonParametricBootstrapAccuracy.R")
mode <- resolve_s6_4_mode(
  script_file, commandArgs(trailingOnly = FALSE), commandArgs(trailingOnly = TRUE),
  override = getOption("s6_4.mode", NULL)
)
script_dir <- dirname(script_file)
setwd(script_dir)
message(
  "S6.4 nonparametric bootstrap accuracy working directory: ",
  normalizePath(".", winslash = "/")
)

analysis_dir <- "../Res_2_Data_for_Analysis"
step5_dir <- "../Step_5_mlVAR_Network_Analysis"
helper_file <- file.path(step5_dir, "S5_0_mlVARHelpers.R")
source(helper_file, local = TRUE, encoding = "UTF-8")

# EDIT ONLY IF THE USER APPROVES A DIFFERENT MEMORY/CPU CEILING.
# Six workers are started; every worker fits mlVAR with one internal core.
parallel_workers <- 6L
internal_mlvar_cores <- 1L
n_replications <- 1000L
raw_p_threshold <- .05
expected_drop_proportions <- c(.10, .20, .30, .40, .50, .60, .70, .75)
expected_replications_per_proportion <- 125L
edge_absolute_tolerance <- 1e-2
edge_correlation_threshold <- .99999

required_packages <- c("lme4", "mlVAR")
available <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
if (!all(available)) {
  stop(
    "Missing required R packages: ",
    paste(names(available)[!available], collapse = ", ")
  )
}

input_file <- file.path(analysis_dir, "dat_step5_mlvar.rda")
network_file <- file.path(analysis_dir, "results_s5_3_mlvar_networks.rda")
cluster_source_file <- file.path(
  analysis_dir, "checkpoint_s5_6_cluster_bootstrap.rda"
)
case_source_file <- file.path(
  analysis_dir, "checkpoint_s6_3_case_drop_edges.rda"
)
s5_6_case_source_file <- file.path(
  analysis_dir, "checkpoint_s5_6_case_drop.rda"
)
s6_3_replication_dir <- file.path(
  analysis_dir, "checkpoints_s6_3_case_drop_edges", "replications"
)
checkpoint_root <- file.path(
  analysis_dir, "checkpoints_s6_4_nonparametric_bootstrap_accuracy"
)
cluster_success_dir <- file.path(checkpoint_root, "cluster", "replications")
case_success_dir <- file.path(checkpoint_root, "case_drop", "replications")
failure_dir <- file.path(checkpoint_root, "failed_attempts")
invalid_dir <- file.path(checkpoint_root, "invalid_files")
log_dir <- file.path(checkpoint_root, "logs")
lock_dir <- file.path(checkpoint_root, "run.lock")
manifest_file <- file.path(
  checkpoint_root, "manifest_s6_4_nonparametric_bootstrap_accuracy.rda"
)
canonical_checkpoint_file <- file.path(
  analysis_dir, "checkpoint_s6_4_nonparametric_bootstrap_accuracy.rda"
)

# Existing completed replication files predate this script rename. Keep their
# historical sole-object name readable, while all newly fitted files use the
# new S6.4 object name. This avoids rewriting 2,000 valid durable checkpoints.
replication_object_name <-
  "replication_s6_4_nonparametric_bootstrap_accuracy"
legacy_replication_object_name <- "replication_s6_4_sparse_bootstrap"

stopifnot(all(file.exists(c(
  helper_file, input_file, network_file, cluster_source_file,
  case_source_file, s5_6_case_source_file
))))

dat_step5_mlvar <- load_single_object(input_file, "dat_step5_mlvar")
results_s5_3_mlvar_networks <- load_single_object(
  network_file, "results_s5_3_mlvar_networks"
)
checkpoint_s5_6_cluster_bootstrap <- load_single_object(
  cluster_source_file, "checkpoint_s5_6_cluster_bootstrap"
)
checkpoint_s6_3_case_drop_edges <- load_single_object(
  case_source_file, "checkpoint_s6_3_case_drop_edges"
)
checkpoint_s5_6_case_drop <- load_single_object(
  s5_6_case_source_file, "checkpoint_s5_6_case_drop"
)

node_map <- dat_step5_mlvar$node_map
node_names <- node_map$node
model_data <- dat_step5_mlvar$prompt[c(node_names, "Name_ID", "occasion")]
participant_series <- split(model_data, model_data$Name_ID)
participant_ids <- names(participant_series)
network_order <- c("Temporal", "Contemporaneous", "Between-person")
expected_edge_counts <- c(
  Temporal = 144L, Contemporaneous = 66L, `Between-person` = 66L
)

stopifnot(
  length(participant_ids) == 209L,
  nrow(model_data) == 18810L,
  all(vapply(participant_series, nrow, integer(1)) == 90L),
  all(vapply(participant_series, function(x) identical(x$occasion, 1:90), logical(1)))
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
      network, rownames(matrix)[pair[, 1L]], colnames(matrix)[pair[, 2L]],
      sep = "|"
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

p_vectors_by_network <- function(extracted) {
  temporal_indices <- expand.grid(
    from = node_names, to = node_names,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  temporal_p <- extracted$temporal_p[cbind(
    match(temporal_indices$from, node_names),
    match(temporal_indices$to, node_names)
  )]
  names(temporal_p) <- paste(
    "Temporal", temporal_indices$from, temporal_indices$to, sep = "|"
  )

  undirected_p_and <- function(gamma_p, network) {
    pair <- which(upper.tri(gamma_p), arr.ind = TRUE)
    value <- pmax(gamma_p[pair], t(gamma_p)[pair])
    names(value) <- paste(
      network, rownames(gamma_p)[pair[, 1L]], colnames(gamma_p)[pair[, 2L]],
      sep = "|"
    )
    value
  }
  list(
    Temporal = temporal_p,
    Contemporaneous = undirected_p_and(
      extracted$contemporaneous_gamma_p, "Contemporaneous"
    ),
    `Between-person` = undirected_p_and(
      extracted$between_gamma_p, "Between-person"
    )
  )
}

flatten_vectors <- function(vectors) {
  c(vectors[["Temporal"]], vectors[["Contemporaneous"]], vectors[["Between-person"]])
}

sparsify_extracted_matrices <- function(extracted, support_by_network) {
  temporal <- extracted$temporal[node_names, node_names, drop = FALSE]
  temporal[] <- 0
  temporal_indices <- expand.grid(
    from = seq_along(node_names), to = seq_along(node_names),
    KEEP.OUT.ATTRS = FALSE
  )
  temporal[cbind(temporal_indices$from, temporal_indices$to)] <-
    ifelse(
      support_by_network$Temporal,
      extracted$temporal[cbind(temporal_indices$from, temporal_indices$to)],
      0
    )

  make_undirected <- function(weight, support) {
    sparse <- matrix(
      0, nrow = length(node_names), ncol = length(node_names),
      dimnames = list(node_names, node_names)
    )
    pair <- which(upper.tri(sparse), arr.ind = TRUE)
    sparse[pair] <- ifelse(support, weight[pair], 0)
    sparse <- sparse + t(sparse)
    diag(sparse) <- 0
    sparse
  }
  list(
    temporal = temporal,
    contemporaneous = make_undirected(
      extracted$contemporaneous, support_by_network$Contemporaneous
    ),
    between = make_undirected(
      extracted$between, support_by_network[["Between-person"]]
    )
  )
}

baseline_sparse_matrices <- results_s5_3_mlvar_networks$sparse_matrices[
  c("temporal", "contemporaneous", "between")
]
baseline_full_matrices <- results_s5_3_mlvar_networks$full_matrices[
  c("temporal", "contemporaneous", "between")
]
baseline_full_edge_vectors <- edge_vectors_by_network(baseline_full_matrices)
baseline_full_edges <- flatten_vectors(baseline_full_edge_vectors)
baseline_sparse_edge_vectors <- edge_vectors_by_network(baseline_sparse_matrices)
baseline_sparse_edges <- flatten_vectors(baseline_sparse_edge_vectors)

input_md5 <- unname(tools::md5sum(input_file))
network_md5 <- unname(tools::md5sum(network_file))
cluster_source_md5 <- unname(tools::md5sum(cluster_source_file))
case_source_md5 <- unname(tools::md5sum(case_source_file))
s5_6_case_source_md5 <- unname(tools::md5sum(s5_6_case_source_file))
cluster_base_seed <- checkpoint_s5_6_cluster_bootstrap$base_seed
case_base_seed <- checkpoint_s5_6_case_drop$base_seed
case_drop_schedule <- checkpoint_s5_6_case_drop$status$drop_proportion
case_retained_schedule <- checkpoint_s5_6_case_drop$status$retained_participants
model_specification <- checkpoint_s5_6_cluster_bootstrap$specification

stopifnot(
  checkpoint_s5_6_cluster_bootstrap$replications == n_replications,
  checkpoint_s6_3_case_drop_edges$replications == n_replications,
  checkpoint_s5_6_case_drop$replications == n_replications,
  identical(checkpoint_s5_6_cluster_bootstrap$input_md5, input_md5),
  identical(checkpoint_s6_3_case_drop_edges$input_md5, input_md5),
  identical(checkpoint_s5_6_case_drop$input_md5, input_md5),
  identical(checkpoint_s6_3_case_drop_edges$network_md5, network_md5),
  identical(checkpoint_s5_6_cluster_bootstrap$specification, model_specification),
  identical(checkpoint_s5_6_case_drop$specification, model_specification),
  identical(checkpoint_s6_3_case_drop_edges$specification, model_specification),
  identical(dim(checkpoint_s5_6_cluster_bootstrap$edge_estimates), c(1000L, 276L)),
  identical(dim(checkpoint_s6_3_case_drop_edges$edge_estimates), c(1000L, 276L)),
  identical(names(baseline_sparse_edges), colnames(
    checkpoint_s5_6_cluster_bootstrap$edge_estimates
  )),
  identical(names(baseline_sparse_edges), colnames(
    checkpoint_s6_3_case_drop_edges$edge_estimates
  )),
  isTRUE(all.equal(
    baseline_full_edges, checkpoint_s5_6_cluster_bootstrap$baseline_edges,
    tolerance = 1e-12, check.attributes = TRUE
  )),
  isTRUE(all.equal(
    baseline_full_edges, checkpoint_s6_3_case_drop_edges$baseline_edges,
    tolerance = 1e-12, check.attributes = TRUE
  )),
  identical(as.numeric(expected_edge_counts), c(144, 66, 66)),
  all(table(case_drop_schedule) == expected_replications_per_proportion)
)

s6_3_replication_file <- function(replication) {
  file.path(
    s6_3_replication_dir,
    sprintf("rep_%04d_case_drop_edges.rda", replication)
  )
}

case_selected_positions <- lapply(seq_len(n_replications), function(replication) {
  object <- load_single_object(
    s6_3_replication_file(replication), "replication_s6_3_case_drop_edges"
  )
  positions <- as.integer(object$selected_participant_positions)
  set.seed(case_base_seed + replication)
  regenerated <- sample(
    seq_along(participant_ids), case_retained_schedule[[replication]],
    replace = FALSE
  )
  stopifnot(
    identical(positions, as.integer(regenerated)),
    length(positions) == case_retained_schedule[[replication]],
    !anyDuplicated(positions)
  )
  positions
})

replication_file <- function(resampling, replication) {
  directory <- if (resampling == "cluster") cluster_success_dir else case_success_dir
  file.path(
    directory,
    sprintf("rep_%04d_%s_sparse.rda", replication, resampling)
  )
}

failure_file <- function(resampling, replication) {
  file.path(
    failure_dir,
    sprintf(
      "rep_%04d_%s_failed_%s_pid%d.rda", replication, resampling,
      format(Sys.time(), "%Y%m%d_%H%M%S"), Sys.getpid()
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

expected_seed <- function(resampling, replication) {
  as.integer((if (resampling == "cluster") cluster_base_seed else case_base_seed) + replication)
}

expected_positions <- function(resampling, replication) {
  if (resampling == "case_drop") return(case_selected_positions[[replication]])
  set.seed(expected_seed(resampling, replication))
  as.integer(sample(seq_along(participant_ids), length(participant_ids), replace = TRUE))
}

accepted_full_edges <- function(resampling, replication) {
  source <- if (resampling == "cluster") {
    checkpoint_s5_6_cluster_bootstrap$edge_estimates
  } else {
    checkpoint_s6_3_case_drop_edges$edge_estimates
  }
  source[replication, ]
}

validate_replication_object <- function(object, resampling, replication) {
  expected_selected <- expected_positions(resampling, replication)
  expected_drop <- if (resampling == "case_drop") {
    case_drop_schedule[[replication]]
  } else {
    NA_real_
  }
  expected_retained <- if (resampling == "case_drop") {
    case_retained_schedule[[replication]]
  } else {
    length(participant_ids)
  }
  expected_source_md5 <- if (resampling == "cluster") {
    cluster_source_md5
  } else {
    case_source_md5
  }
  p_rule <- is.finite(object$p_values) & object$p_values < raw_p_threshold
  sparse_rule <- ifelse(object$support, object$full_edges, 0)
  computed_support_counts <- vapply(network_order, function(network) {
    sum(object$support[startsWith(names(object$support), paste0(network, "|"))])
  }, integer(1))

  is.list(object) &&
    object$schema_version %in% c("2.0.0", "2.1.0") &&
    identical(object$status, "success") &&
    identical(object$resampling, resampling) &&
    identical(object$replication, as.integer(replication)) &&
    identical(object$seed, expected_seed(resampling, replication)) &&
    isTRUE(all.equal(object$drop_proportion, expected_drop, check.attributes = FALSE)) &&
    identical(object$retained_participants, as.integer(expected_retained)) &&
    identical(object$selected_participant_positions, expected_selected) &&
    identical(object$input_md5, input_md5) &&
    identical(object$network_md5, network_md5) &&
    identical(object$source_checkpoint_md5, expected_source_md5) &&
    identical(object$internal_mlvar_cores, internal_mlvar_cores) &&
    identical(names(object$full_edges), names(baseline_sparse_edges)) &&
    identical(names(object$p_values), names(baseline_sparse_edges)) &&
    identical(names(object$support), names(baseline_sparse_edges)) &&
    identical(names(object$sparse_edges), names(baseline_sparse_edges)) &&
    length(object$full_edges) == 276L &&
    all(is.finite(object$full_edges)) &&
    all(is.finite(object$p_values) | is.na(object$p_values)) &&
    identical(as.logical(object$support), as.logical(p_rule)) &&
    isTRUE(all.equal(
      as.numeric(object$sparse_edges), as.numeric(sparse_rule),
      tolerance = 0, check.attributes = FALSE
    )) &&
    identical(
      as.integer(object$support_counts[network_order]),
      unname(computed_support_counts)
    ) &&
    isTRUE(object$reproduction_audit$passed)
}

load_valid_replication <- function(
    resampling, replication,
    path = replication_file(resampling, replication)) {
  if (!file.exists(path)) return(NULL)
  object <- tryCatch(
    {
      replication_environment <- new.env(parent = emptyenv())
      loaded_names <- load(path, envir = replication_environment)
      permitted_names <- c(
        replication_object_name, legacy_replication_object_name
      )
      if (length(loaded_names) != 1L || !loaded_names %in% permitted_names) {
        stop("Unexpected replication object in ", basename(path))
      }
      value <- replication_environment[[loaded_names]]
      value$sparse_centrality <- NULL
      value
    },
    error = function(error) NULL
  )
  if (is.null(object) || !validate_replication_object(
    object, resampling, replication
  )) return(NULL)
  object
}

archive_invalid_file <- function(path, resampling) {
  if (!file.exists(path)) return(invisible(FALSE))
  destination <- file.path(
    invalid_dir,
    paste0(
      resampling, "_", tools::file_path_sans_ext(basename(path)), "_",
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

task_table <- rbind(
  data.frame(resampling = "cluster", replication = seq_len(n_replications)),
  data.frame(resampling = "case_drop", replication = seq_len(n_replications))
)
task_table$task_id <- seq_len(nrow(task_table))

scan_replications <- function(archive_invalid = FALSE) {
  valid <- logical(nrow(task_table))
  paths <- character(nrow(task_table))
  for (index in seq_len(nrow(task_table))) {
    resampling <- task_table$resampling[[index]]
    replication <- task_table$replication[[index]]
    path <- replication_file(resampling, replication)
    paths[[index]] <- path
    if (!file.exists(path)) next
    object <- load_valid_replication(resampling, replication, path)
    valid[[index]] <- !is.null(object)
    if (is.null(object) && archive_invalid) archive_invalid_file(path, resampling)
  }
  data.frame(task_table, valid = valid, path = paths, stringsAsFactors = FALSE)
}

write_manifest <- function(scan, state) {
  manifest_s6_4_nonparametric_bootstrap_accuracy <- list(
    schema_version = "2.1.0",
    state = state,
    updated = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    parallel_workers = parallel_workers,
    internal_mlvar_cores = internal_mlvar_cores,
    replications_per_resampling = n_replications,
    successful_files = sum(scan$valid),
    pending_files = sum(!scan$valid),
    completion_by_resampling = aggregate(valid ~ resampling, scan, sum),
    raw_p_threshold = raw_p_threshold,
    temporal_rule = "two-sided unadjusted p < .05 per replication",
    undirected_rule = "p_AND = max(two node-wise p values) < .05 per replication",
    input_md5 = input_md5,
    network_md5 = network_md5,
    historical_cluster_network_md5 =
      checkpoint_s5_6_cluster_bootstrap$network_md5,
    historical_case_drop_network_md5 = checkpoint_s5_6_case_drop$network_md5,
    cluster_source_md5 = cluster_source_md5,
    case_source_md5 = case_source_md5,
    canonical_checkpoint_file = canonical_checkpoint_file
  )
  save_single_object_atomic(
    manifest_s6_4_nonparametric_bootstrap_accuracy,
    "manifest_s6_4_nonparametric_bootstrap_accuracy", manifest_file
  )
  invisible(manifest_s6_4_nonparametric_bootstrap_accuracy)
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
  extracted <- extract_mlvar_matrices(fit)
  weights <- extracted[c("temporal", "contemporaneous", "between")]
  stopifnot(all(vapply(weights, function(x) {
    identical(dim(x), c(12L, 12L)) && all(is.finite(x))
  }, logical(1))))

  edge_vectors <- edge_vectors_by_network(weights)
  p_vectors <- p_vectors_by_network(extracted)
  support_by_network <- lapply(
    p_vectors, function(p) is.finite(p) & p < raw_p_threshold
  )
  sparse_matrices <- sparsify_extracted_matrices(extracted, support_by_network)
  sparse_edge_vectors <- edge_vectors_by_network(sparse_matrices)

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

  result <- list(
    full_edge_vectors = edge_vectors,
    full_edges = flatten_vectors(edge_vectors),
    p_values_by_network = p_vectors,
    p_values = flatten_vectors(p_vectors),
    support_by_network = support_by_network,
    support = flatten_vectors(support_by_network),
    sparse_edge_vectors = sparse_edge_vectors,
    sparse_edges = flatten_vectors(sparse_edge_vectors),
    support_counts = vapply(support_by_network, sum, integer(1)),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    optimizer_nonzero = sum(optimizer_codes != 0L),
    singular_models = sum(singular),
    nonboundary_messages = sum(nonboundary_messages),
    warning_text = paste(captured_warnings, collapse = " | "),
    temporal_spectral_radius = matrix_spectral_radius(weights$temporal)
  )
  rm(fit, all_models, extracted, weights, sparse_matrices)
  invisible(gc(verbose = FALSE))
  result
}

write_worker_log <- function(text) {
  path <- file.path(log_dir, sprintf("worker_pid_%d.log", Sys.getpid()))
  cat(text, "\n", file = path, append = TRUE)
  invisible(path)
}

fit_one_task <- function(task_id) {
  resampling <- task_table$resampling[[task_id]]
  replication <- task_table$replication[[task_id]]
  seed <- expected_seed(resampling, replication)
  selected_positions <- expected_positions(resampling, replication)
  drop_proportion <- if (resampling == "case_drop") {
    case_drop_schedule[[replication]]
  } else {
    NA_real_
  }
  retained_participants <- length(selected_positions)
  started <- Sys.time()
  data <- make_resample(selected_positions)
  result <- tryCatch(fit_resample_values(data), error = function(error) error)
  finished <- Sys.time()

  if (!inherits(result, "error")) {
    accepted <- accepted_full_edges(resampling, replication)
    maximum_absolute_difference <- max(abs(result$full_edges - accepted))
    vector_correlation <- stats::cor(result$full_edges, accepted)
    reproduction_passed <-
      is.finite(maximum_absolute_difference) &&
      is.finite(vector_correlation) &&
      maximum_absolute_difference <= edge_absolute_tolerance &&
      vector_correlation >= edge_correlation_threshold
    if (!reproduction_passed) {
      result <- simpleError(sprintf(
        paste0(
          "Accepted full-edge reproduction failed: max_abs_diff=%.6g; ",
          "correlation=%.10f"
        ), maximum_absolute_difference, vector_correlation
      ))
    }
  }

  if (inherits(result, "error")) {
    failure_s6_4_nonparametric_bootstrap_accuracy <- list(
      schema_version = "2.1.0", status = "failed",
      resampling = resampling, replication = as.integer(replication),
      seed = seed, drop_proportion = drop_proportion,
      retained_participants = as.integer(retained_participants),
      selected_participant_positions = as.integer(selected_positions),
      input_md5 = input_md5, network_md5 = network_md5,
      source_checkpoint_md5 = if (resampling == "cluster") {
        cluster_source_md5
      } else {
        case_source_md5
      },
      internal_mlvar_cores = internal_mlvar_cores,
      started = format(started, "%Y-%m-%d %H:%M:%S %Z"),
      finished = format(finished, "%Y-%m-%d %H:%M:%S %Z"),
      elapsed_seconds = as.numeric(difftime(finished, started, units = "secs")),
      error_text = conditionMessage(result)
    )
    path <- failure_file(resampling, replication)
    save_single_object_atomic(
      failure_s6_4_nonparametric_bootstrap_accuracy,
      "failure_s6_4_nonparametric_bootstrap_accuracy", path
    )
    write_worker_log(sprintf(
      "%s rep=%04d FAILED elapsed=%.1fs error=%s",
      resampling, replication,
      failure_s6_4_nonparametric_bootstrap_accuracy$elapsed_seconds,
      failure_s6_4_nonparametric_bootstrap_accuracy$error_text
    ))
    rm(data, result)
    invisible(gc(verbose = FALSE))
    return(list(task_id = task_id, status = "failed", path = path))
  }

  accepted <- accepted_full_edges(resampling, replication)
  reproduction_audit <- list(
    accepted_checkpoint = if (resampling == "cluster") {
      basename(cluster_source_file)
    } else {
      basename(case_source_file)
    },
    maximum_absolute_difference = max(abs(result$full_edges - accepted)),
    vector_correlation = stats::cor(result$full_edges, accepted),
    absolute_tolerance = edge_absolute_tolerance,
    correlation_threshold = edge_correlation_threshold,
    passed = TRUE
  )
  replication_s6_4_nonparametric_bootstrap_accuracy <- list(
    schema_version = "2.1.0", status = "success",
    resampling = resampling, replication = as.integer(replication),
    seed = seed, drop_proportion = drop_proportion,
    retained_participants = as.integer(retained_participants),
    selected_participant_positions = as.integer(selected_positions),
    input_md5 = input_md5, network_md5 = network_md5,
    source_checkpoint_md5 = if (resampling == "cluster") {
      cluster_source_md5
    } else {
      case_source_md5
    },
    model_specification = model_specification,
    raw_p_threshold = raw_p_threshold,
    temporal_support_rule = "two-sided unadjusted p < .05",
    undirected_support_rule = "p_AND = max(two node-wise p values) < .05",
    internal_mlvar_cores = internal_mlvar_cores,
    full_edge_vectors = result$full_edge_vectors,
    full_edges = result$full_edges,
    p_values_by_network = result$p_values_by_network,
    p_values = result$p_values,
    support_by_network = result$support_by_network,
    support = result$support,
    sparse_edge_vectors = result$sparse_edge_vectors,
    sparse_edges = result$sparse_edges,
    support_counts = result$support_counts,
    reproduction_audit = reproduction_audit,
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
      c("mlVAR", "lme4"),
      function(package) as.character(packageVersion(package)), character(1)
    )
  )
  validation_passed <- tryCatch(
    validate_replication_object(
      replication_s6_4_nonparametric_bootstrap_accuracy,
      resampling, replication
    ),
    error = function(error) FALSE
  )
  if (!isTRUE(validation_passed)) {
    failure_s6_4_nonparametric_bootstrap_accuracy <- list(
      schema_version = "2.1.0", status = "failed_validation",
      resampling = resampling, replication = as.integer(replication),
      seed = seed,
      candidate = replication_s6_4_nonparametric_bootstrap_accuracy,
      error_text = "Candidate replication object failed structural validation"
    )
    path <- failure_file(resampling, replication)
    save_single_object_atomic(
      failure_s6_4_nonparametric_bootstrap_accuracy,
      "failure_s6_4_nonparametric_bootstrap_accuracy", path
    )
    write_worker_log(sprintf(
      "%s rep=%04d FAILED_VALIDATION path=%s",
      resampling, replication, basename(path)
    ))
    rm(data, result, replication_s6_4_nonparametric_bootstrap_accuracy)
    invisible(gc(verbose = FALSE))
    return(list(task_id = task_id, status = "failed", path = path))
  }
  path <- replication_file(resampling, replication)
  save_single_object_atomic(
    replication_s6_4_nonparametric_bootstrap_accuracy,
    replication_object_name, path
  )
  stopifnot(!is.null(load_valid_replication(resampling, replication, path)))
  write_worker_log(sprintf(
    paste0(
      "%s rep=%04d SUCCESS elapsed=%.1fs support=%d/%d/%d ",
      "max_abs_diff=%.3e cor=%.9f"
    ),
    resampling, replication, result$elapsed_seconds,
    result$support_counts[["Temporal"]],
    result$support_counts[["Contemporaneous"]],
    result$support_counts[["Between-person"]],
    reproduction_audit$maximum_absolute_difference,
    reproduction_audit$vector_correlation
  ))
  rm(data, result, replication_s6_4_nonparametric_bootstrap_accuracy)
  invisible(gc(verbose = FALSE))
  list(task_id = task_id, status = "success", path = path)
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
      load_single_object(
        lock_file, "s6_4_nonparametric_bootstrap_accuracy_run_lock"
      ),
      error = function(error) NULL
    )
    if (!is.null(lock) && process_is_running(lock$pid)) {
      stop("An S6.4 nonparametric bootstrap run is active with PID ", lock$pid)
    }
    stale <- paste0(lock_dir, "_stale_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    if (!file.rename(lock_dir, stale)) {
      stop("A stale run lock exists and could not be archived: ", lock_dir)
    }
  }
  if (!dir.create(lock_dir, recursive = FALSE, showWarnings = FALSE)) {
    stop("Could not acquire S6.4 nonparametric bootstrap run lock")
  }
  s6_4_nonparametric_bootstrap_accuracy_run_lock <- list(
    pid = Sys.getpid(),
    started = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    script_file = script_file,
    parallel_workers = parallel_workers,
    internal_mlvar_cores = internal_mlvar_cores
  )
  save_single_object_atomic(
    s6_4_nonparametric_bootstrap_accuracy_run_lock,
    "s6_4_nonparametric_bootstrap_accuracy_run_lock",
    file.path(lock_dir, "lock.rda")
  )
  invisible(s6_4_nonparametric_bootstrap_accuracy_run_lock)
}

release_run_lock <- function() {
  if (dir.exists(lock_dir)) unlink(lock_dir, recursive = TRUE, force = TRUE)
  invisible(TRUE)
}

assemble_resampling <- function(resampling) {
  p_values <- matrix(
    NA_real_, n_replications, length(baseline_sparse_edges),
    dimnames = list(seq_len(n_replications), names(baseline_sparse_edges))
  )
  support <- matrix(
    FALSE, n_replications, length(baseline_sparse_edges),
    dimnames = dimnames(p_values)
  )
  sparse_edges <- matrix(
    NA_real_, n_replications, length(baseline_sparse_edges),
    dimnames = dimnames(p_values)
  )
  status <- data.frame(
    replication = seq_len(n_replications),
    resampling = resampling,
    drop_proportion = if (resampling == "case_drop") case_drop_schedule else NA_real_,
    retained_participants = if (resampling == "case_drop") {
      case_retained_schedule
    } else {
      rep(length(participant_ids), n_replications)
    },
    status = "success", elapsed_seconds = NA_real_,
    optimizer_nonzero = NA_integer_, singular_models = NA_integer_,
    nonboundary_messages = NA_integer_, temporal_spectral_radius = NA_real_,
    warning_text = NA_character_, started = NA_character_, finished = NA_character_,
    file = vapply(
      seq_len(n_replications),
      function(replication) replication_file(resampling, replication),
      character(1)
    ),
    file_md5 = NA_character_, reproduction_max_abs_difference = NA_real_,
    reproduction_correlation = NA_real_, stringsAsFactors = FALSE
  )

  for (replication in seq_len(n_replications)) {
    object <- load_valid_replication(resampling, replication)
    if (is.null(object)) stop("Invalid file during merge: ", resampling, " ", replication)
    p_values[replication, ] <- object$p_values
    support[replication, ] <- object$support
    sparse_edges[replication, ] <- object$sparse_edges
    status$elapsed_seconds[[replication]] <- object$diagnostics$elapsed_seconds
    status$optimizer_nonzero[[replication]] <- object$diagnostics$optimizer_nonzero
    status$singular_models[[replication]] <- object$diagnostics$singular_models
    status$nonboundary_messages[[replication]] <- object$diagnostics$nonboundary_messages
    status$temporal_spectral_radius[[replication]] <- object$diagnostics$temporal_spectral_radius
    status$warning_text[[replication]] <- object$diagnostics$warning_text
    status$started[[replication]] <- object$started
    status$finished[[replication]] <- object$finished
    status$reproduction_max_abs_difference[[replication]] <-
      object$reproduction_audit$maximum_absolute_difference
    status$reproduction_correlation[[replication]] <-
      object$reproduction_audit$vector_correlation
  }
  status$file_md5 <- unname(tools::md5sum(status$file))
  network_columns <- lapply(network_order, function(network) {
    startsWith(colnames(sparse_edges), paste0(network, "|"))
  })
  names(network_columns) <- network_order
  sparse_edges_by_network <- lapply(network_columns, function(index) {
    sparse_edges[, index, drop = FALSE]
  })
  p_values_by_network <- lapply(network_columns, function(index) {
    p_values[, index, drop = FALSE]
  })
  support_by_network <- lapply(network_columns, function(index) {
    support[, index, drop = FALSE]
  })
  support_counts <- vapply(
    support_by_network, rowSums, numeric(n_replications)
  )
  list(
    p_values = p_values,
    p_values_by_network = p_values_by_network,
    support = support,
    support_by_network = support_by_network,
    sparse_edges = sparse_edges,
    sparse_edges_by_network = sparse_edges_by_network,
    support_counts = support_counts,
    status = status,
    reproduction_audit = list(
      maximum_absolute_difference = max(status$reproduction_max_abs_difference),
      minimum_vector_correlation = min(status$reproduction_correlation),
      absolute_tolerance = edge_absolute_tolerance,
      correlation_threshold = edge_correlation_threshold,
      passed = all(status$reproduction_max_abs_difference <= edge_absolute_tolerance) &&
        all(status$reproduction_correlation >= edge_correlation_threshold)
    )
  )
}

merge_completed_replications <- function() {
  scan <- scan_replications(archive_invalid = FALSE)
  if (!all(scan$valid)) {
    stop("Cannot merge: ", sum(!scan$valid), " replication files are missing or invalid")
  }
  cluster <- assemble_resampling("cluster")
  case_drop <- assemble_resampling("case_drop")
  validation <- list(
    all_2000_replication_files_valid = all(scan$valid),
    exact_sparse_edge_dimensions =
      identical(dim(cluster$sparse_edges), c(1000L, 276L)) &&
      identical(dim(case_drop$sparse_edges), c(1000L, 276L)),
    replicate_specific_support_rule =
      identical(cluster$support, is.finite(cluster$p_values) & cluster$p_values < .05) &&
      identical(case_drop$support, is.finite(case_drop$p_values) & case_drop$p_values < .05),
    unsupported_edges_are_exact_zero =
      all(cluster$sparse_edges[!cluster$support] == 0) &&
      all(case_drop$sparse_edges[!case_drop$support] == 0),
    every_sparse_value_finite =
      all(is.finite(cluster$sparse_edges)) &&
      all(is.finite(case_drop$sparse_edges)),
    exact_network_dimensions = all(vapply(network_order, function(network) {
      identical(
        dim(cluster$sparse_edges_by_network[[network]]),
        c(1000L, unname(expected_edge_counts[[network]]))
      ) && identical(
        dim(case_drop$sparse_edges_by_network[[network]]),
        c(1000L, unname(expected_edge_counts[[network]]))
      )
    }, logical(1))),
    exact_case_drop_balance = identical(
      as.integer(table(case_drop$status$drop_proportion)), rep(125L, 8L)
    ),
    six_workers_one_internal_core =
      parallel_workers == 6L && internal_mlvar_cores == 1L,
    full_edge_reproduction_passed =
      isTRUE(cluster$reproduction_audit$passed) &&
      isTRUE(case_drop$reproduction_audit$passed)
  )
  stopifnot(all(unlist(validation)))

  checkpoint_s6_4_nonparametric_bootstrap_accuracy <- list(
    schema_version = "2.1.0",
    created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    replications_per_resampling = n_replications,
    parallel_workers = parallel_workers,
    internal_mlvar_cores = internal_mlvar_cores,
    raw_p_threshold = raw_p_threshold,
    temporal_support_rule = "two-sided unadjusted p < .05 within each replication",
    undirected_support_rule = paste(
      "p_AND = max(two node-wise p values) < .05 within each replication"
    ),
    specification = model_specification,
    input_md5 = input_md5,
    network_md5 = network_md5,
    historical_cluster_network_md5 =
      checkpoint_s5_6_cluster_bootstrap$network_md5,
    historical_case_drop_network_md5 = checkpoint_s5_6_case_drop$network_md5,
    cluster_source_md5 = cluster_source_md5,
    case_source_md5 = case_source_md5,
    s5_6_case_source_md5 = s5_6_case_source_md5,
    baseline_sparse_matrices = baseline_sparse_matrices,
    baseline_full_edge_vectors = baseline_full_edge_vectors,
    baseline_full_edges = baseline_full_edges,
    baseline_sparse_edge_vectors = baseline_sparse_edge_vectors,
    baseline_sparse_edges = baseline_sparse_edges,
    baseline_support_counts = vapply(
      baseline_sparse_edge_vectors, function(x) sum(x != 0), integer(1)
    ),
    cluster = cluster,
    case_drop = case_drop,
    checkpoint_directory = checkpoint_root,
    r_version = R.version.string,
    package_versions = vapply(
      c("mlVAR", "lme4"),
      function(package) as.character(packageVersion(package)), character(1)
    ),
    validation = validation
  )
  save_single_object_atomic(
    checkpoint_s6_4_nonparametric_bootstrap_accuracy,
    "checkpoint_s6_4_nonparametric_bootstrap_accuracy",
    canonical_checkpoint_file
  )
  reloaded <- load_single_object(
    canonical_checkpoint_file,
    "checkpoint_s6_4_nonparametric_bootstrap_accuracy"
  )
  stopifnot(all(unlist(reloaded$validation)))
  invisible(checkpoint_s6_4_nonparametric_bootstrap_accuracy)
}

main <- function(mode) {
# Read-only audit does not create directories. All checkpoint mutations require
# ownership of the controller lock, including invalid-file archival at startup.
if (mode != "audit") {
  dir.create(checkpoint_root, recursive = TRUE, showWarnings = FALSE)
  acquire_run_lock()
  on.exit(release_run_lock(), add = TRUE)
  for (path in c(cluster_success_dir, case_success_dir, failure_dir, invalid_dir, log_dir)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}
initial_scan <- scan_replications(archive_invalid = identical(mode, "run"))
if (mode != "audit") {
  write_manifest(initial_scan, if (all(initial_scan$valid)) "complete" else "incomplete")
}
cat(sprintf(
  paste0(
      "S6.4 nonparametric bootstrap audit: valid=%d pending=%d ",
    "(cluster=%d/%d, case_drop=%d/%d) workers=%d internal_cores=%d mode=%s\n"
  ),
  sum(initial_scan$valid), sum(!initial_scan$valid),
  sum(initial_scan$valid & initial_scan$resampling == "cluster"), n_replications,
  sum(initial_scan$valid & initial_scan$resampling == "case_drop"), n_replications,
  parallel_workers, internal_mlvar_cores, mode
))

if (mode == "audit") {
  invisible(initial_scan)
} else if (mode == "merge") {
  merge_completed_replications()
    cat("S6.4 nonparametric bootstrap merge completed.\n")
} else {
  pending <- initial_scan$task_id[!initial_scan$valid]
  if (length(pending)) {
    cat(sprintf(
      paste0(
        "S6.4 nonparametric bootstrap fitting %d pending tasks with ",
        "six single-core workers.\n"
      ),
      length(pending)
    ))
    flush.console()
    worker_console_file <- file.path(log_dir, "psock_worker_console.log")
    cluster_handle <- parallel::makePSOCKcluster(
      parallel_workers, outfile = worker_console_file
    )
    on.exit({
      if (!is.null(cluster_handle)) try(parallel::stopCluster(cluster_handle), silent = TRUE)
    }, add = TRUE, after = FALSE)
    parallel::clusterEvalQ(cluster_handle, {
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
      cluster_handle,
      varlist = c(
        "helper_file", "node_names", "node_map", "participant_series",
        "participant_ids", "network_order", "expected_edge_counts",
        "raw_p_threshold", "internal_mlvar_cores", "input_md5", "network_md5",
        "cluster_source_md5", "case_source_md5", "cluster_base_seed",
        "case_base_seed", "case_drop_schedule", "case_retained_schedule",
        "case_selected_positions", "model_specification", "task_table",
        "baseline_sparse_edges", "replication_object_name",
        "legacy_replication_object_name",
        "checkpoint_s5_6_cluster_bootstrap", "checkpoint_s6_3_case_drop_edges",
        "cluster_source_file", "case_source_file", "edge_absolute_tolerance",
        "edge_correlation_threshold", "cluster_success_dir", "case_success_dir",
        "failure_dir", "log_dir", "edge_vectors_by_network",
        "p_vectors_by_network", "flatten_vectors", "sparsify_extracted_matrices",
        "replication_file", "failure_file",
        "save_single_object_atomic", "expected_seed", "expected_positions",
        "accepted_full_edges", "validate_replication_object",
        "load_valid_replication", "make_resample", "fit_resample_values",
        "write_worker_log", "fit_one_task"
      ),
      envir = environment()
    )
    parallel::clusterEvalQ(cluster_handle, {
      source(helper_file, encoding = "UTF-8")
      stopifnot(
        requireNamespace("mlVAR", quietly = TRUE),
        requireNamespace("lme4", quietly = TRUE)
      )
      NULL
    })
    run_results <- parallel::parLapplyLB(cluster_handle, pending, fit_one_task)
    parallel::stopCluster(cluster_handle)
    cluster_handle <- NULL
    run_status <- vapply(run_results, `[[`, character(1), "status")
    cat(sprintf(
      "S6.4 nonparametric bootstrap returned success=%d failed=%d.\n",
      sum(run_status == "success"), sum(run_status == "failed")
    ))
  }
  final_scan <- scan_replications(archive_invalid = TRUE)
  write_manifest(
    final_scan,
    if (all(final_scan$valid)) "complete" else "incomplete_after_run"
  )
  if (!all(final_scan$valid)) {
    stop(
      "S6.4 nonparametric bootstrap remains incomplete. Re-run this script; only ",
      "missing or invalid replications will be fitted."
    )
  }
  merge_completed_replications()
  write_manifest(scan_replications(), "complete")
  cat(sprintf(
    paste0(
      "S6.4 nonparametric bootstrap COMPLETE: 2,000 durable files and ",
      "%s written.\n"
    ),
    basename(canonical_checkpoint_file)
  ))
}
}

main(mode)
