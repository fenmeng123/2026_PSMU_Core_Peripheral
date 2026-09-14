# Author: Kunru Song

args_all <- commandArgs(trailingOnly = FALSE)
script_argument <- args_all[grepl("^--file=", args_all)]
if (length(script_argument) != 1L) stop("Could not identify merger script path")
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
shard_count <- if (length(args)) as.integer(args[[1L]]) else 2L
poll_seconds <- if (length(args) >= 2L) as.numeric(args[[2L]]) else 60
if (
  is.na(shard_count) || shard_count < 2L ||
    is.na(poll_seconds) || poll_seconds < 1
) {
  stop("shard_count must be at least 2 and poll_seconds must be positive")
}

analysis_directory <- file.path(project_root, "Res_2_Data_for_Analysis")
canonical_paths <- c(
  cluster = file.path(
    analysis_directory, "checkpoint_s5_6_cluster_bootstrap.rda"
  ),
  case_drop = file.path(
    analysis_directory, "checkpoint_s5_6_case_drop.rda"
  )
)
object_names <- c(
  cluster = "checkpoint_s5_6_cluster_bootstrap",
  case_drop = "checkpoint_s5_6_case_drop"
)

shard_path <- function(type, shard_index) {
  sub(
    "\\.rda$",
    sprintf("_shard%02dof%02d.rda", shard_index, shard_count),
    canonical_paths[[type]]
  )
}

load_with_previous_fallback <- function(path, object_name) {
  candidates <- c(path, paste0(path, ".previous"))
  for (candidate in candidates[file.exists(candidates)]) {
    value <- tryCatch(
      load_single_object(candidate, object_name),
      error = function(error) NULL
    )
    if (!is.null(value)) return(value)
  }
  NULL
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
  if (!copied) stop("Could not update canonical checkpoint: ", path)
  invisible(path)
}

assigned_replications <- function(shard_index, replications = 1000L) {
  which(
    (seq_len(replications) - 1L) %% shard_count == shard_index - 1L
  )
}

canonical_complete <- function(value) {
  !is.null(value) &&
    identical(value$replications, 1000L) &&
    nrow(value$status) == 1000L &&
    !any(value$status$status == "pending") &&
    !is.null(value$completed)
}

shards_ready <- function(shards) {
  if (any(vapply(shards, is.null, logical(1)))) return(FALSE)
  all(vapply(seq_along(shards), function(index) {
    shard <- shards[[index]]
    expected <- assigned_replications(index)
    identical(shard$replications, 1000L) &&
      identical(shard$shard_index, as.integer(index)) &&
      identical(shard$shard_count, shard_count) &&
      identical(shard$assigned_replications, expected) &&
      identical(shard$assigned_completed, length(expected)) &&
      !is.null(shard$completed) &&
      !any(shard$status$status[expected] == "pending")
  }, logical(1)))
}

validate_shards <- function(shards, type) {
  replication_owners <- do.call(rbind, lapply(seq_along(shards), function(index) {
    data.frame(
      replication = shards[[index]]$assigned_replications,
      shard_index = index
    )
  }))
  owner_table <- table(replication_owners$replication)
  stopifnot(
    length(owner_table) == 1000L,
    identical(as.integer(names(owner_table)), seq_len(1000L)),
    all(owner_table == 1L)
  )

  reference <- shards[[1L]]
  for (shard in shards[-1L]) {
    stopifnot(
      identical(shard$replications, reference$replications),
      identical(shard$base_seed, reference$base_seed),
      identical(shard$specification, reference$specification),
      identical(shard$input_md5, reference$input_md5),
      identical(shard$network_md5, reference$network_md5)
    )
    if (type == "cluster") {
      stopifnot(
        identical(shard$baseline_edges, reference$baseline_edges),
        identical(
          colnames(shard$edge_estimates),
          colnames(reference$edge_estimates)
        )
      )
    } else {
      stopifnot(
        identical(shard$baseline_centrality, reference$baseline_centrality),
        identical(
          colnames(shard$centrality_estimates),
          colnames(reference$centrality_estimates)
        ),
        identical(shard$drop_proportions, reference$drop_proportions),
        identical(
          shard$replications_per_proportion,
          reference$replications_per_proportion
        ),
        identical(
          shard$status$drop_proportion,
          reference$status$drop_proportion
        )
      )
    }
  }
  invisible(TRUE)
}

merge_shards <- function(shards, type) {
  if (type == "case_drop") shards <- lapply(shards, retain_case_drop_centrality)
  validate_shards(shards, type)
  merged <- shards[[1L]]
  merged$status[,] <- shards[[1L]]$status[,]
  merged$status$status <- "pending"
  merged$status$started <- NA_character_
  merged$status$finished <- NA_character_
  merged$status$elapsed_seconds <- NA_real_
  merged$status$optimizer_nonzero <- NA_integer_
  merged$status$singular_models <- NA_integer_
  merged$status$nonboundary_messages <- NA_integer_
  merged$status$temporal_spectral_radius <- NA_real_
  merged$status$warning_text <- NA_character_
  merged$status$error_text <- NA_character_
  if (type == "cluster") {
    merged$baseline_centrality <- NULL
    merged$centrality_estimates <- NULL
    merged$edge_estimates[,] <- NA_real_
  } else {
    merged$centrality_estimates[,] <- NA_real_
  }

  for (index in seq_along(shards)) {
    rows <- assigned_replications(index)
    merged$status[rows, ] <- shards[[index]]$status[rows, ]
    if (type == "cluster") {
      merged$edge_estimates[rows, ] <- shards[[index]]$edge_estimates[rows, ]
    } else {
      merged$centrality_estimates[rows, ] <-
        shards[[index]]$centrality_estimates[rows, ]
    }
  }

  merged$worker_cores <- sum(vapply(shards, function(x) x$worker_cores, numeric(1)))
  merged$shard_index <- NULL
  merged$shard_count <- shard_count
  merged$assigned_replications <- seq_len(1000L)
  merged$assigned_completed <- 1000L
  merged$execution <- list(
    strategy = sprintf("%d deterministic disjoint replication shards", shard_count),
    ownership_rule = "(replication - 1) modulo shard_count",
    shard_worker_cores = vapply(shards, function(x) x$worker_cores, numeric(1)),
    shard_files = vapply(
      seq_len(shard_count), function(index) {
        normalizePath(shard_path(type, index), winslash = "/", mustWork = TRUE)
      }, character(1)
    ),
    shard_md5 = unname(vapply(
      seq_len(shard_count), function(index) {
        tools::md5sum(shard_path(type, index))
      }, character(1)
    )),
    merged = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
  merged$completed <- merged$execution$merged

  stopifnot(
    !any(merged$status$status == "pending"),
    sum(merged$status$status %in% c("success", "failed")) == 1000L,
    all(merged$status$replication == seq_len(1000L)),
    sum(merged$status$status == "failed") <= 50L
  )
  if (type == "case_drop") {
    drop_attempts <- table(merged$status$drop_proportion)
    stopifnot(length(drop_attempts) == 8L, all(drop_attempts == 125L))
  }
  merged
}

merge_one_if_ready <- function(type) {
  canonical <- load_with_previous_fallback(
    canonical_paths[[type]], object_names[[type]]
  )
  if (canonical_complete(canonical)) return(TRUE)

  shards <- lapply(seq_len(shard_count), function(index) {
    load_with_previous_fallback(
      shard_path(type, index), object_names[[type]]
    )
  })
  if (!shards_ready(shards)) return(FALSE)

  merged <- merge_shards(shards, type)
  backup_path <- paste0(canonical_paths[[type]], ".pre_shard")
  if (file.exists(canonical_paths[[type]]) && !file.exists(backup_path)) {
    file.copy(canonical_paths[[type]], backup_path, overwrite = FALSE)
  }
  save_checkpoint(merged, object_names[[type]], canonical_paths[[type]])
  reloaded <- load_single_object(canonical_paths[[type]], object_names[[type]])
  stopifnot(canonical_complete(reloaded), identical(reloaded$status, merged$status))
  cat(sprintf(
    "S5.6 merged %s shards: success=%d failed=%d time=%s\n",
    type, sum(merged$status$status == "success"),
    sum(merged$status$status == "failed"),
    format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ))
  flush.console()
  TRUE
}

cat(sprintf(
  "S5.6 shard merger started shards=%d time=%s\n",
  shard_count, format(Sys.time(), "%Y-%m-%d %H:%M:%S")
))
flush.console()

last_counts <- NULL
repeat {
  cluster_done <- merge_one_if_ready("cluster")
  case_drop_done <- merge_one_if_ready("case_drop")
  if (cluster_done && case_drop_done) break

  counts <- unlist(lapply(c("cluster", "case_drop"), function(type) {
    vapply(seq_len(shard_count), function(index) {
      value <- load_with_previous_fallback(
        shard_path(type, index), object_names[[type]]
      )
      if (is.null(value)) 0L else {
        rows <- assigned_replications(index)
        sum(value$status$status[rows] != "pending")
      }
    }, integer(1))
  }))
  if (!identical(counts, last_counts)) {
    cat(sprintf(
      paste0(
        "S5.6 shard merger waiting cluster=%s case_drop=%s time=%s\n"
      ),
      paste(counts[seq_len(shard_count)], collapse = "/"),
      paste(counts[shard_count + seq_len(shard_count)], collapse = "/"),
      format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ))
    flush.console()
    last_counts <- counts
  }
  Sys.sleep(poll_seconds)
}

cat(sprintf(
  "S5.6 shard merger completed successfully time=%s\n",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
))
