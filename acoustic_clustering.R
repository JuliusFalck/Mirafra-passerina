# Audio-only clustering; irruption labels are used only for external validation.
acoustic_adjusted_rand <- function(x, y) {
  if (length(x) != length(y) || length(x) < 2L || anyNA(x) || anyNA(y)) {
    stop("ARI requires two complete label vectors of equal length >= 2.")
  }
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  joint <- sum(choose2(tab))
  a <- sum(choose2(rowSums(tab)))
  b <- sum(choose2(colSums(tab)))
  expected <- a * b / choose2(length(x))
  denominator <- (a + b) / 2 - expected
  # Identical trivial partitions (all one group or all singletons).
  if (denominator == 0) return(1)
  (joint - expected) / denominator
}

acoustic_cluster_robustness <- function(pair_data, metadata,
                                       nperm = 9999L, seed = 37L) {
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("Install the recommended R package 'cluster' for silhouette widths.")
  }
  if (!all(c("id_i", "id_j", "combo_sim") %in% names(pair_data)) ||
      !all(c("name", "irruption", "season") %in% names(metadata))) {
    stop("Need pair IDs/combo_sim and metadata name/irruption/season.")
  }
  if (length(nperm) != 1L || !is.finite(nperm) ||
      nperm < 1 || nperm != as.integer(nperm)) stop("nperm must be a positive integer.")
  ids <- as.character(metadata$name)
  n <- length(ids)
  if (n < 3L || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) ||
      anyNA(metadata$irruption) || anyNA(metadata$season)) {
    stop("Need >= 3 unique recording IDs with complete irruption and season labels.")
  }
  # Reconstruct by recording ID, never by an assumed metadata row order.
  i <- match(as.character(pair_data$id_i), ids)
  j <- match(as.character(pair_data$id_j), ids)
  similarity <- pair_data$combo_sim
  if (anyNA(i) || anyNA(j) || any(i == j) || !is.numeric(similarity) ||
      any(!is.finite(similarity)) || any(similarity < 0 | similarity > 1)) {
    stop("Pairs must reference known distinct IDs and finite similarities in [0, 1].")
  }
  keys <- paste(pmin(i, j), pmax(i, j), sep = ":")
  if (nrow(pair_data) != choose(n, 2) || anyDuplicated(keys)) {
    stop("Need exactly one composite similarity for every unordered recording pair.")
  }
  similarity_matrix <- diag(1, n)
  dimnames(similarity_matrix) <- list(ids, ids)
  similarity_matrix[cbind(i, j)] <- similarity
  similarity_matrix[cbind(j, i)] <- similarity
  acoustic_distance <- stats::as.dist(1 - similarity_matrix)
  if (diff(range(acoustic_distance)) == 0) {
    stop("All composite dissimilarities are equal; an acoustic partition is not identifiable.")
  }
  # Average linkage accepts general dissimilarities; no Euclidean assumption.
  tree <- stats::hclust(acoustic_distance, method = "average")
  candidate_k <- seq.int(2L, n - 1L)
  partitions <- matrix(stats::cutree(tree, k = candidate_k), nrow = n)
  mean_silhouette <- vapply(seq_along(candidate_k), function(z) {
    mean(cluster::silhouette(partitions[, z], acoustic_distance)[, "sil_width"])
  }, numeric(1))
  # Select solely on audio, with ties resolved toward fewer clusters.
  best <- which.max(mean_silhouette)
  selected_k <- candidate_k[best]
  selected <- partitions[, best]

  reference <- as.character(metadata$irruption)
  season <- as.character(metadata$season)
  if (any(vapply(split(season, reference), function(z) length(unique(z)),
                 integer(1)) != 1L)) {
    stop("Each irruption must belong to exactly one season.")
  }
  reference_k <- length(unique(reference))
  # Secondary, externally constrained resolution; not the audio-selected result.
  matched <- stats::cutree(tree, k = reference_k)
  profile <- data.frame(
    k = candidate_k, mean_silhouette = mean_silhouette,
    adjusted_rand = vapply(seq_along(candidate_k), function(z) {
      acoustic_adjusted_rand(partitions[, z], reference)
    }, numeric(1)),
    selected_by_audio = candidate_k == selected_k
  )
  solutions <- list(audio_selected = selected, matched_irruption_k = matched)
  summary <- do.call(rbind, lapply(names(solutions), function(label) {
    membership <- solutions[[label]]
    k <- length(unique(membership))
    data.frame(solution = label, n_recordings = n, n_irruptions = reference_k,
               k = k, mean_silhouette = if (k > 1L && k < n) {
                 mean(cluster::silhouette(membership, acoustic_distance)[, "sil_width"])
               } else NA_real_,
               adjusted_rand = acoustic_adjusted_rand(membership, reference),
               singleton_clusters = sum(table(membership) == 1L))
  }))
  rownames(summary) <- NULL

  # Local seed: reproducible without changing subsequent analyses' RNG state.
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  season_indices <- split(seq_len(n), season)
  can_permute_within <- any(vapply(season_indices, function(idx) {
    length(unique(reference[idx])) > 1L
  }, logical(1)))
  null_global <- matrix(NA_real_, nperm, length(solutions),
                        dimnames = list(NULL, names(solutions)))
  null_within <- null_global
  for (b in seq_len(nperm)) {
    # Both fixed partitions use the same recording-label shuffles.
    # The matched K remains valid because permutations preserve irruption counts.
    permuted_global <- reference[sample.int(n)]
    null_global[b, ] <- vapply(solutions, function(membership) {
      acoustic_adjusted_rand(membership, permuted_global)
    }, numeric(1))
    if (can_permute_within) {
      permuted <- reference
      for (idx in season_indices) {
        permuted[idx] <- reference[idx[sample.int(length(idx))]]
      }
      null_within[b, ] <- vapply(solutions, function(membership) {
        acoustic_adjusted_rand(membership, permuted)
      }, numeric(1))
    }
  }
  permutation_summary <- do.call(rbind, lapply(names(solutions), function(label) {
    observed <- summary$adjusted_rand[match(label, summary$solution)]
    do.call(rbind, lapply(c("global", "within_season"), function(scheme) {
      null <- if (scheme == "global") null_global[, label] else null_within[, label]
      available <- any(is.finite(null))
      data.frame(solution = label, permutation = scheme,
                 adjusted_rand = observed,
                 null_mean = if (available) mean(null) else NA_real_,
                 p_greater = if (available) {
                   (1 + sum(null >= observed - 1e-12)) / (nperm + 1)
                 } else NA_real_,
                 nperm = if (available) nperm else 0L,
                 note = if (available) "" else "No season contains multiple irruptions")
    }))
  }))
  assignments <- data.frame(name = ids, irruption = reference, season = season,
                            acoustic_cluster = unname(selected),
                            acoustic_cluster_matched_k = unname(matched),
                            stringsAsFactors = FALSE)
  contingency <- as.data.frame(table(
    irruption = factor(reference, levels = unique(reference)),
    acoustic_cluster = factor(selected, levels = sort(unique(selected)))
  ), responseName = "n_recordings")
  contingency$within_irruption_fraction <- contingency$n_recordings /
    as.numeric(table(reference)[as.character(contingency$irruption)])
  list(tree = tree, distance = acoustic_distance, selected_k = selected_k,
       profile = profile, summary = summary, assignments = assignments,
       contingency = contingency, permutation_summary = permutation_summary,
       permutation_null = data.frame(
         global = null_global[, "audio_selected"],
         within_season = null_within[, "audio_selected"],
         matched_irruption_k_global = null_global[, "matched_irruption_k"],
         matched_irruption_k_within_season = null_within[, "matched_irruption_k"]))
}
