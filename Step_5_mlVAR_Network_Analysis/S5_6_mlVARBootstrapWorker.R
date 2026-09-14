# Author: Kunru Song

args_all <- commandArgs(trailingOnly = FALSE)
script_argument <- args_all[grepl("^--file=", args_all)]
if (length(script_argument) != 1L) stop("Could not identify worker script path")
script_path <- normalizePath(
  sub("^--file=", "", script_argument), winslash = "/", mustWork = TRUE
)
project_root <- normalizePath(
  file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE
)
source(
  file.path(dirname(script_path), "S5_0_mlVARHelpers.R"),
  encoding = "UTF-8"
)

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) tolower(args[[1L]]) else "all"
if (!mode %in% c("cluster", "case-drop", "all")) {
  stop("Mode must be cluster, case-drop, or all")
}
shard_index <- if (length(args) >= 2L) as.integer(args[[2L]]) else 1L
shard_count <- if (length(args) >= 3L) as.integer(args[[3L]]) else 1L
requested_worker_cores <- if (length(args) >= 4L) as.integer(args[[4L]]) else 12L
if (
  anyNA(c(shard_index, shard_count, requested_worker_cores)) ||
    shard_count < 1L || shard_index < 1L || shard_index > shard_count ||
    requested_worker_cores < 1L
) {
  stop("Invalid shard_index, shard_count, or requested_worker_cores")
}
shard_suffix <- if (shard_count == 1L) "" else {
  sprintf("_shard%02dof%02d", shard_index, shard_count)
}

Sys.setenv(LANG = "C", LC_ALL = "C")
try(Sys.setlocale("LC_ALL", "C"), silent = TRUE)
options(warn = 1)

input_file <- file.path(
  project_root, "Res_2_Data_for_Analysis", "dat_step5_mlvar.rda"
)
network_file <- file.path(
  project_root, "Res_2_Data_for_Analysis", "results_s5_3_mlvar_networks.rda"
)
canonical_cluster_checkpoint_file <- file.path(
  project_root, "Res_2_Data_for_Analysis", "checkpoint_s5_6_cluster_bootstrap.rda"
)
canonical_case_drop_checkpoint_file <- file.path(
  project_root, "Res_2_Data_for_Analysis", "checkpoint_s5_6_case_drop.rda"
)
cluster_checkpoint_file <- sub(
  "\\.rda$", paste0(shard_suffix, ".rda"), canonical_cluster_checkpoint_file
)
case_drop_checkpoint_file <- sub(
  "\\.rda$", paste0(shard_suffix, ".rda"), canonical_case_drop_checkpoint_file
)

dat_step5_mlvar <- load_single_object(input_file, "dat_step5_mlvar")
results_s5_3_mlvar_networks <- load_single_object(
  network_file, "results_s5_3_mlvar_networks"
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

n_cluster_replications <- 1000L
n_case_drop_replications <- 1000L
drop_proportions <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.75)
case_drop_schedule <- rep(drop_proportions, each = n_case_drop_replications / length(drop_proportions))
worker_cores <- min(requested_worker_cores, parallel::detectCores(logical = FALSE))
assigned_cluster_replications <- which(
  (seq_len(n_cluster_replications) - 1L) %% shard_count == shard_index - 1L
)
assigned_case_drop_replications <- which(
  (seq_len(n_case_drop_replications) - 1L) %% shard_count == shard_index - 1L
)

matrix_list <- list(
  temporal = results_s5_3_mlvar_networks$full_matrices$temporal,
  contemporaneous = results_s5_3_mlvar_networks$full_matrices$contemporaneous,
  between = results_s5_3_mlvar_networks$full_matrices$between
)

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
    values[paste("Temporal", "Bridge outgoing expected influence", node, sep = "|")] <-
      sum(temporal[node, other])
    values[paste("Temporal", "Bridge incoming expected influence", node, sep = "|")] <-
      sum(temporal[other, node])
  }
  for (network in c("Contemporaneous", "Between-person")) {
    matrix <- if (network == "Contemporaneous") contemporaneous else between
    for (node in node_names) {
      values[paste(network, "Strength", node, sep = "|")] <- sum(abs(matrix[node, ]))
    }
  }
  values
}

edge_vector <- function(matrices) {
  temporal_indices <- expand.grid(
    from = node_names, to = node_names,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  temporal_values <- matrices$temporal[cbind(
    match(temporal_indices$from, node_names), match(temporal_indices$to, node_names)
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
  c(
    temporal_values,
    undirected_values(matrices$contemporaneous, "Contemporaneous"),
    undirected_values(matrices$between, "Between-person")
  )
}

baseline_edges <- edge_vector(matrix_list)
baseline_centrality <- if (mode %in% c("case-drop", "all")) {
  centrality_vector(matrix_list)
} else NULL

make_resample <- function(selected_ids) {
  result <- do.call(rbind, lapply(seq_along(selected_ids), function(index) {
    series <- participant_series[[selected_ids[[index]]]]
    series$Name_ID <- index
    series
  }))
  rownames(result) <- NULL
  result
}

fit_resample <- function(data, compute_centrality = TRUE) {
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
      nCores = worker_cores,
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
  stopifnot(
    all(vapply(matrices, function(x) {
      identical(dim(x), c(12L, 12L)) && all(is.finite(x))
    }, logical(1)))
  )

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
    edges = edge_vector(matrices),
    centrality = if (compute_centrality) centrality_vector(matrices) else NULL,
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

save_checkpoint <- function(object, object_name, path) {
  environment <- new.env(parent = emptyenv())
  environment[[object_name]] <- object
  new_path <- paste0(path, ".new")
  previous_path <- paste0(path, ".previous")
  save(list = object_name, envir = environment, file = new_path, version = 2)
  if (file.exists(path)) file.copy(path, previous_path, overwrite = TRUE)
  copied <- file.copy(new_path, path, overwrite = TRUE)
  unlink(new_path)
  if (!copied) stop("Could not update checkpoint: ", path)
  invisible(path)
}

load_checkpoint <- function(path, object_name) {
  candidates <- c(path, paste0(path, ".previous"))
  for (candidate in candidates[file.exists(candidates)]) {
    value <- tryCatch(
      load_single_object(candidate, object_name), error = function(error) NULL
    )
    if (!is.null(value)) {
      if (identical(object_name, "checkpoint_s5_6_case_drop")) {
        value$baseline_centrality <- select_named_network_values(
          value$baseline_centrality, names(baseline_centrality)
        )
        value$centrality_estimates <- select_named_network_values(
          value$centrality_estimates, names(baseline_centrality)
        )
      }
      return(value)
    }
  }
  NULL
}

new_status <- function(replications) {
  data.frame(
    replication = seq_len(replications), status = "pending",
    drop_proportion = NA_real_, retained_participants = NA_integer_,
    started = NA_character_, finished = NA_character_, elapsed_seconds = NA_real_,
    optimizer_nonzero = NA_integer_, singular_models = NA_integer_,
    nonboundary_messages = NA_integer_, temporal_spectral_radius = NA_real_,
    warning_text = NA_character_, error_text = NA_character_,
    stringsAsFactors = FALSE
  )
}

run_cluster <- function() {
  object_name <- "checkpoint_s5_6_cluster_bootstrap"
  checkpoint <- load_checkpoint(cluster_checkpoint_file, object_name)
  # Downstream S6 consumers use full edges, schedules and fit diagnostics only.
  # Discard unused summary fields in memory when resuming an older checkpoint.
  if (!is.null(checkpoint)) {
    checkpoint$baseline_centrality <- NULL
    checkpoint$centrality_estimates <- NULL
  }
  if (is.null(checkpoint)) {
    checkpoint <- list(
      type = "participant-cluster nonparametric bootstrap with replacement",
      replications = n_cluster_replications,
      base_seed = 5261000L,
      worker_cores = worker_cores,
      shard_index = shard_index,
      shard_count = shard_count,
      assigned_replications = assigned_cluster_replications,
      specification = list(
        estimator = "lmer", lags = 1L, temporal = "orthogonal",
        contemporaneous = "orthogonal", scale = TRUE, scaleWithin = FALSE,
        na.rm = TRUE, full_detrend = FALSE, beepvar = "occasion", dayvar = NULL
      ),
      input_md5 = unname(tools::md5sum(input_file)),
      network_md5 = unname(tools::md5sum(network_file)),
      baseline_edges = baseline_edges,
      edge_estimates = matrix(
        NA_real_, nrow = n_cluster_replications, ncol = length(baseline_edges),
        dimnames = list(seq_len(n_cluster_replications), names(baseline_edges))
      ),
      status = new_status(n_cluster_replications),
      created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    )
    if (shard_count > 1L && file.exists(canonical_cluster_checkpoint_file)) {
      inherited <- load_checkpoint(
        canonical_cluster_checkpoint_file, object_name
      )
      if (!is.null(inherited)) {
        inherited_rows <- intersect(
          assigned_cluster_replications,
          which(inherited$status$status != "pending")
        )
        inherited_compatible <-
          identical(inherited$input_md5, unname(tools::md5sum(input_file))) &&
          isTRUE(all.equal(
            inherited$baseline_edges, baseline_edges,
            tolerance = 1e-12, check.attributes = TRUE
          )) && identical(inherited$specification, checkpoint$specification)
        if (length(inherited_rows) && inherited_compatible) {
          checkpoint$status[inherited_rows, ] <- inherited$status[inherited_rows, ]
          checkpoint$edge_estimates[inherited_rows, ] <-
            inherited$edge_estimates[inherited_rows, ]
        }
      }
    }
  }
  stopifnot(
    checkpoint$replications == n_cluster_replications,
    identical(colnames(checkpoint$edge_estimates), names(baseline_edges)),
    identical(checkpoint$input_md5, unname(tools::md5sum(input_file))),
    isTRUE(all.equal(
      checkpoint$baseline_edges, baseline_edges,
      tolerance = 1e-12, check.attributes = TRUE
    )),
    identical(checkpoint$specification, list(
      estimator = "lmer", lags = 1L, temporal = "orthogonal",
      contemporaneous = "orthogonal", scale = TRUE, scaleWithin = FALSE,
      na.rm = TRUE, full_detrend = FALSE, beepvar = "occasion", dayvar = NULL
    ))
  )

  pending <- intersect(
    assigned_cluster_replications,
    which(checkpoint$status$status == "pending")
  )
  for (replication in pending) {
    set.seed(checkpoint$base_seed + replication)
    selected <- sample(participant_ids, length(participant_ids), replace = TRUE)
    data <- make_resample(selected)
    checkpoint$status$started[[replication]] <- format(
      Sys.time(), "%Y-%m-%d %H:%M:%S %Z"
    )
    result <- tryCatch(
      fit_resample(data, compute_centrality = FALSE), error = function(error) error
    )
    checkpoint$status$finished[[replication]] <- format(
      Sys.time(), "%Y-%m-%d %H:%M:%S %Z"
    )
    checkpoint$status$retained_participants[[replication]] <- length(selected)
    if (inherits(result, "error")) {
      checkpoint$status$status[[replication]] <- "failed"
      checkpoint$status$error_text[[replication]] <- conditionMessage(result)
    } else {
      checkpoint$status$status[[replication]] <- "success"
      checkpoint$status$elapsed_seconds[[replication]] <- result$elapsed_seconds
      checkpoint$status$optimizer_nonzero[[replication]] <- result$optimizer_nonzero
      checkpoint$status$singular_models[[replication]] <- result$singular_models
      checkpoint$status$nonboundary_messages[[replication]] <- result$nonboundary_messages
      checkpoint$status$temporal_spectral_radius[[replication]] <-
        result$temporal_spectral_radius
      checkpoint$status$warning_text[[replication]] <- result$warning_text
      checkpoint$edge_estimates[replication, ] <- result$edges
    }
    save_checkpoint(checkpoint, object_name, cluster_checkpoint_file)
    cat(sprintf(
      paste0(
        "CLUSTER shard=%d/%d rep=%04d assigned=%04d/%04d %s elapsed=%.1fs ",
        "success=%d failed=%d time=%s\n"
      ),
      shard_index, shard_count,
      replication, match(replication, assigned_cluster_replications),
      length(assigned_cluster_replications), checkpoint$status$status[[replication]],
      checkpoint$status$elapsed_seconds[[replication]] %||% NA_real_,
      sum(checkpoint$status$status[assigned_cluster_replications] == "success"),
      sum(checkpoint$status$status[assigned_cluster_replications] == "failed"),
      format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ))
    flush.console()
    rm(data, result)
    invisible(gc(verbose = FALSE))
  }
  stopifnot(!any(
    checkpoint$status$status[assigned_cluster_replications] == "pending"
  ))
  checkpoint$assigned_completed <- length(assigned_cluster_replications)
  checkpoint$completed <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  save_checkpoint(checkpoint, object_name, cluster_checkpoint_file)
  invisible(checkpoint)
}

run_case_drop <- function() {
  object_name <- "checkpoint_s5_6_case_drop"
  checkpoint <- load_checkpoint(case_drop_checkpoint_file, object_name)
  if (is.null(checkpoint)) {
    status <- new_status(n_case_drop_replications)
    status$drop_proportion <- case_drop_schedule
    status$retained_participants <- pmax(
      2L, as.integer(round(length(participant_ids) * (1 - case_drop_schedule)))
    )
    checkpoint <- list(
      type = "participant case-dropping bootstrap without replacement",
      replications = n_case_drop_replications,
      base_seed = 5262000L,
      drop_proportions = drop_proportions,
      replications_per_proportion = n_case_drop_replications / length(drop_proportions),
      correlation_threshold = 0.70,
      worker_cores = worker_cores,
      shard_index = shard_index,
      shard_count = shard_count,
      assigned_replications = assigned_case_drop_replications,
      specification = list(
        estimator = "lmer", lags = 1L, temporal = "orthogonal",
        contemporaneous = "orthogonal", scale = TRUE, scaleWithin = FALSE,
        na.rm = TRUE, full_detrend = FALSE, beepvar = "occasion", dayvar = NULL
      ),
      input_md5 = unname(tools::md5sum(input_file)),
      network_md5 = unname(tools::md5sum(network_file)),
      baseline_centrality = baseline_centrality,
      centrality_estimates = matrix(
        NA_real_, nrow = n_case_drop_replications, ncol = length(baseline_centrality),
        dimnames = list(seq_len(n_case_drop_replications), names(baseline_centrality))
      ),
      status = status,
      created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    )
    if (shard_count > 1L && file.exists(canonical_case_drop_checkpoint_file)) {
      inherited <- load_checkpoint(
        canonical_case_drop_checkpoint_file, object_name
      )
      if (!is.null(inherited)) {
        inherited_rows <- intersect(
          assigned_case_drop_replications,
          which(inherited$status$status != "pending")
        )
        inherited_compatible <-
          identical(inherited$input_md5, unname(tools::md5sum(input_file))) &&
          isTRUE(all.equal(
            inherited$baseline_centrality, baseline_centrality,
            tolerance = 1e-12, check.attributes = TRUE
          ))
        if (length(inherited_rows) && inherited_compatible) {
          checkpoint$status[inherited_rows, ] <- inherited$status[inherited_rows, ]
          checkpoint$centrality_estimates[inherited_rows, ] <-
            inherited$centrality_estimates[inherited_rows, ]
        }
      }
    }
  }
  stopifnot(
    checkpoint$replications == n_case_drop_replications,
    identical(checkpoint$drop_proportions, drop_proportions),
    identical(colnames(checkpoint$centrality_estimates), names(baseline_centrality)),
    identical(checkpoint$input_md5, unname(tools::md5sum(input_file))),
    isTRUE(all.equal(
      checkpoint$baseline_centrality, baseline_centrality,
      tolerance = 1e-12, check.attributes = TRUE
    ))
  )

  pending <- intersect(
    assigned_case_drop_replications,
    which(checkpoint$status$status == "pending")
  )
  for (replication in pending) {
    set.seed(checkpoint$base_seed + replication)
    retained <- checkpoint$status$retained_participants[[replication]]
    selected <- sample(participant_ids, retained, replace = FALSE)
    data <- make_resample(selected)
    checkpoint$status$started[[replication]] <- format(
      Sys.time(), "%Y-%m-%d %H:%M:%S %Z"
    )
    result <- tryCatch(fit_resample(data), error = function(error) error)
    checkpoint$status$finished[[replication]] <- format(
      Sys.time(), "%Y-%m-%d %H:%M:%S %Z"
    )
    if (inherits(result, "error")) {
      checkpoint$status$status[[replication]] <- "failed"
      checkpoint$status$error_text[[replication]] <- conditionMessage(result)
    } else {
      checkpoint$status$status[[replication]] <- "success"
      checkpoint$status$elapsed_seconds[[replication]] <- result$elapsed_seconds
      checkpoint$status$optimizer_nonzero[[replication]] <- result$optimizer_nonzero
      checkpoint$status$singular_models[[replication]] <- result$singular_models
      checkpoint$status$nonboundary_messages[[replication]] <- result$nonboundary_messages
      checkpoint$status$temporal_spectral_radius[[replication]] <-
        result$temporal_spectral_radius
      checkpoint$status$warning_text[[replication]] <- result$warning_text
      checkpoint$centrality_estimates[replication, ] <- result$centrality
    }
    save_checkpoint(checkpoint, object_name, case_drop_checkpoint_file)
    cat(sprintf(
      paste0(
        "CASE_DROP shard=%d/%d rep=%04d assigned=%04d/%04d ",
        "drop=%.2f retained=%d %s elapsed=%.1fs ",
        "success=%d failed=%d time=%s\n"
      ),
      shard_index, shard_count,
      replication, match(replication, assigned_case_drop_replications),
      length(assigned_case_drop_replications),
      checkpoint$status$drop_proportion[[replication]], retained,
      checkpoint$status$status[[replication]],
      checkpoint$status$elapsed_seconds[[replication]] %||% NA_real_,
      sum(checkpoint$status$status[assigned_case_drop_replications] == "success"),
      sum(checkpoint$status$status[assigned_case_drop_replications] == "failed"),
      format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ))
    flush.console()
    rm(data, result)
    invisible(gc(verbose = FALSE))
  }
  stopifnot(!any(
    checkpoint$status$status[assigned_case_drop_replications] == "pending"
  ))
  checkpoint$assigned_completed <- length(assigned_case_drop_replications)
  checkpoint$completed <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  save_checkpoint(checkpoint, object_name, case_drop_checkpoint_file)
  invisible(checkpoint)
}

cat(sprintf(
  "S5.6 worker started mode=%s shard=%d/%d assigned=%d cores=%d time=%s\n",
  mode, shard_index, shard_count, length(assigned_cluster_replications),
  worker_cores, format(Sys.time(), "%Y-%m-%d %H:%M:%S")
))
flush.console()

if (mode %in% c("cluster", "all")) run_cluster()
if (mode %in% c("case-drop", "all")) run_case_drop()

cat(sprintf(
  "S5.6 worker completed mode=%s shard=%d/%d time=%s\n",
  mode, shard_index, shard_count, format(Sys.time(), "%Y-%m-%d %H:%M:%S")
))
