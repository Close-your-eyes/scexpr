# =============================================================================
# EDUCATIONAL SOUPX-LIKE AMBIENT-RNA CORRECTION
# =============================================================================
#
# Purpose
# -------
# This script exposes the statistical ideas behind SoupX using ordinary R
# matrices. It intentionally does NOT call the SoupX package.
#
# The four ideas demonstrated are:
#   1. Learn the ambient-RNA composition from empty droplets.
#   2. Find marker genes that should be absent from selected cell clusters.
#   3. Estimate the contamination fraction from those "negative" observations.
#   4. Allocate and subtract the expected ambient counts without producing
#      negative expression values.
#
# This is a teaching implementation, not a drop-in replacement for SoupX.
# Production SoupX has more careful marker selection, pruning, uncertainty
# handling, cluster-to-cell expansion, sparse-matrix code, and diagnostics.
#
# Matrix convention used throughout:
#   rows    = genes
#   columns = droplets or cells
# =============================================================================


# =============================================================================
# PART 0 — A SMALL, RUNNABLE EXAMPLE
# =============================================================================
#
# The raw matrix contains two empty droplets plus six cell-containing droplets.
# HBB is abundant in the empty droplets, representing RNA released by lysed red
# blood cells. Some HBB molecules therefore appear in every cell type.

genes <- c("HBB", "LST1", "CD3D", "MS4A1", "MALAT1")
droplets <- c("empty_1", "empty_2", "T_1", "T_2",
              "Mono_1", "Mono_2", "B_1", "B_2")

raw_umi <- matrix(
  c(
    # empty_1 empty_2 T_1 T_2 Mono_1 Mono_2 B_1 B_2
          20,      15,   2,  3,     2,     3,  2,  3,  # HBB
           1,       1,   0,  1,    14,    16,  1,  0,  # LST1
           1,       0,  12, 14,     0,     1,  1,  0,  # CD3D
           0,       1,   1,  0,     0,     1, 11, 13,  # MS4A1
           8,      10,  25, 27,    24,    26, 23, 25   # MALAT1
  ),
  nrow = length(genes),
  byrow = TRUE,
  dimnames = list(genes, droplets)
)

cell_barcodes <- c("T_1", "T_2", "Mono_1", "Mono_2", "B_1", "B_2")
filtered_umi <- raw_umi[, cell_barcodes, drop = FALSE]

# SoupX relies on preliminary clusters. These do not need to be final cell-type
# annotations; they only need to pool transcriptionally similar cells.
clusters <- c(
  T_1 = "T", T_2 = "T",
  Mono_1 = "Mono", Mono_2 = "Mono",
  B_1 = "B", B_2 = "B"
)


# =============================================================================
# PART 1 — ESTIMATE THE AMBIENT "SOUP" PROFILE
# =============================================================================
#
# Empty droplets contain no intact cell, so their UMIs approximate the ambient
# RNA pool. If S_g is the number of empty-droplet UMIs for gene g, then
#
#                         p_g = S_g / sum_g(S_g)
#
# is the probability that a randomly sampled ambient molecule belongs to gene g.

estimate_soup_profile <- function(raw_counts, empty_barcodes) {
  stopifnot(is.matrix(raw_counts))
  stopifnot(all(empty_barcodes %in% colnames(raw_counts)))

  ambient_counts <- rowSums(raw_counts[, empty_barcodes, drop = FALSE])
  if (sum(ambient_counts) == 0) {
    stop("The selected empty droplets contain no UMIs.")
  }

  ambient_counts / sum(ambient_counts)
}

# Here the empty droplets are known. With real data, SoupX usually selects
# low-UMI barcodes from the unfiltered/raw droplet matrix.
empty_barcodes <- setdiff(colnames(raw_umi), colnames(filtered_umi))
soup_profile <- estimate_soup_profile(raw_umi, empty_barcodes)

cat("\nSTEP 1: ambient RNA proportions\n")
print(round(sort(soup_profile, decreasing = TRUE), 3))


# =============================================================================
# PART 2 — UNDERSTAND TF-IDF MARKER SELECTION
# =============================================================================
#
# SoupX's automated procedure seeks genes that are specific to one cluster.
# Conceptually, its TF-IDF marker score rewards a gene when:
#   - it is detected frequently inside one cluster (term frequency), and
#   - it is detected infrequently across all cells (inverse document frequency).
#
# This compact implementation is for illustration. Production SoupX also uses
# statistical filtering and additional thresholds.

calculate_tfidf_markers <- function(cell_counts, cluster_labels) {
  stopifnot(identical(colnames(cell_counts), names(cluster_labels)))

  detected <- cell_counts > 0
  global_frequency <- rowMeans(detected)
  inverse_document_frequency <- -log(pmax(global_frequency, 1e-12))

  cluster_names <- unique(cluster_labels)
  answer <- lapply(cluster_names, function(k) {
    cells_in_cluster <- names(cluster_labels)[cluster_labels == k]
    within_cluster_frequency <- rowMeans(
      detected[, cells_in_cluster, drop = FALSE]
    )

    data.frame(
      gene = rownames(cell_counts),
      cluster = k,
      within_cluster_frequency = within_cluster_frequency,
      global_frequency = global_frequency,
      tfidf = within_cluster_frequency * inverse_document_frequency,
      row.names = NULL
    )
  })

  answer <- do.call(rbind, answer)
  answer[order(answer$cluster, -answer$tfidf), ]
}

tfidf_table <- calculate_tfidf_markers(filtered_umi, clusters)

cat("\nSTEP 2: strongest illustrative TF-IDF marker for each cluster\n")
print(do.call(
  rbind,
  lapply(split(tfidf_table, tfidf_table$cluster), head, n = 2)
))


# =============================================================================
# PART 3 — DEFINE WHERE MARKERS SHOULD NOT BE EXPRESSED
# =============================================================================
#
# The essential SoupX insight is not merely that a gene is a marker. It is that
# a marker should be absent from the other clusters.
#
# marker_owner maps each gene to the cluster that genuinely expresses it.
# An NA owner means the gene should be endogenous in none of the sampled cells.
# HBB is the clearest example here: erythroid cells are absent from the dataset,
# so HBB observed in T, Mono, or B cells is treated as ambient evidence.

marker_owner <- c(
  HBB = NA_character_,
  LST1 = "Mono",
  CD3D = "T",
  MS4A1 = "B"
)

make_negative_pairs <- function(marker_owner, cluster_names) {
  pairs <- lapply(names(marker_owner), function(gene) {
    owner <- unname(marker_owner[gene])
    negative_clusters <- if (is.na(owner)) {
      cluster_names
    } else {
      setdiff(cluster_names, owner)
    }

    data.frame(
      gene = gene,
      cluster = negative_clusters,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, pairs)
}

negative_pairs <- make_negative_pairs(marker_owner, unique(clusters))


# =============================================================================
# PART 4 — ESTIMATE THE CONTAMINATION FRACTION, rho
# =============================================================================
#
# Aggregate similar cells first. For gene g in cluster k:
#
#   O_gk = observed UMI count
#   E_gk = N_k * p_g
#
# where N_k is the cluster's total UMI count and p_g is the soup proportion.
# E_gk is the expected count if 100% of the cluster library were soup.
# If the contamination fraction is rho, SoupX uses the approximation
#
#                   O_gk ~ Poisson(rho * E_gk).
#
# Therefore O_gk / E_gk is an intuitive, noisy estimate of rho.

aggregate_by_cluster <- function(cell_counts, cluster_labels) {
  cluster_names <- unique(cluster_labels)
  aggregated <- sapply(cluster_names, function(k) {
    cells <- names(cluster_labels)[cluster_labels == k]
    rowSums(cell_counts[, cells, drop = FALSE])
  })

  rownames(aggregated) <- rownames(cell_counts)
  aggregated
}

build_rho_evidence <- function(cell_counts,
                               cluster_labels,
                               soup_profile,
                               negative_pairs,
                               maximum_plausible_rho = 0.50,
                               expression_fdr = 0.05,
                               minimum_exposure = 0.25) {
  aggregated <- aggregate_by_cluster(cell_counts, cluster_labels)
  cluster_umis <- colSums(aggregated)

  evidence <- negative_pairs
  evidence$observed <- mapply(
    function(g, k) aggregated[g, k],
    evidence$gene,
    evidence$cluster
  )
  evidence$exposure <- mapply(
    function(g, k) cluster_umis[k] * soup_profile[g],
    evidence$gene,
    evidence$cluster
  )
  evidence$naive_rho <- evidence$observed / evidence$exposure

  # Conservative pruning principle:
  # Test whether expression is greater than even a generously contaminated
  # droplet could explain. A small upper-tail probability suggests genuine
  # endogenous expression, so that marker/cluster pair should not estimate rho.
  evidence$p_if_maximally_contaminated <- stats::ppois(
    evidence$observed - 1,
    lambda = maximum_plausible_rho * evidence$exposure,
    lower.tail = FALSE
  )
  evidence$q_value <- stats::p.adjust(
    evidence$p_if_maximally_contaminated,
    method = "BH"
  )
  evidence$use_for_rho <-
    evidence$exposure >= minimum_exposure &
    evidence$q_value >= expression_fdr

  evidence
}

rho_evidence <- build_rho_evidence(
  filtered_umi,
  clusters,
  soup_profile,
  negative_pairs
)

cat("\nSTEP 3: negative marker/cluster evidence\n")
print(rho_evidence)


# Combine the retained Poisson likelihoods with a weak gamma prior. SoupX uses
# gamma-based posterior information as well; evaluating a grid makes the idea
# especially visible.
estimate_rho_on_grid <- function(evidence,
                                 rho_grid = seq(0.001, 0.50, by = 0.001),
                                 prior_shape = 2,
                                 prior_rate = 20) {
  use <- evidence$use_for_rho
  if (!any(use)) {
    stop("No negative marker/cluster pairs survived pruning.")
  }

  observed <- evidence$observed[use]
  exposure <- evidence$exposure[use]

  log_prior <- stats::dgamma(
    rho_grid,
    shape = prior_shape,
    rate = prior_rate,
    log = TRUE
  )

  log_likelihood <- vapply(rho_grid, function(rho) {
    sum(stats::dpois(observed, lambda = rho * exposure, log = TRUE))
  }, numeric(1))

  log_posterior <- log_prior + log_likelihood
  posterior <- exp(log_posterior - max(log_posterior))
  posterior <- posterior / sum(posterior)

  list(
    rho = rho_grid[which.max(posterior)],
    curve = data.frame(
      rho = rho_grid,
      log_likelihood = log_likelihood,
      posterior_probability = posterior
    )
  )
}

rho_fit <- estimate_rho_on_grid(rho_evidence)
rho_hat <- rho_fit$rho

cat(sprintf("\nSTEP 4: estimated contamination fraction rho = %.3f\n", rho_hat))

# Uncomment in an interactive R session to visualize uncertainty:
plot(
  rho_fit$curve$rho,
  rho_fit$curve$posterior_probability,
  type = "l",
  xlab = "Contamination fraction (rho)",
  ylab = "Posterior probability",
  main = "Ambient-RNA contamination estimate"
)
abline(v = rho_hat, col = "red", lty = 2)


# =============================================================================
# PART 5 — CONSTRAINED SUBTRACTION OF AMBIENT COUNTS
# =============================================================================
#
# For a cell or cluster containing N UMIs, the expected number of ambient UMIs is
#
#                              M = rho * N.
#
# The first guess for gene g is M * p_g. Direct subtraction can fail because a
# sparse count vector may contain fewer observed UMIs than this expectation.
# The allocator below therefore:
#   - never assigns more ambient UMIs to a gene than were observed;
#   - redistributes unused ambient mass to other genes according to p_g; and
#   - stops after approximately M molecules have been assigned.

allocate_ambient_counts <- function(observed,
                                    target_ambient,
                                    soup_profile,
                                    tolerance = 1e-10,
                                    max_iterations = 1000L) {
  stopifnot(identical(names(observed), names(soup_profile)))

  observed <- as.numeric(observed)
  names(observed) <- names(soup_profile)
  target_ambient <- min(target_ambient, sum(observed))

  allocation <- setNames(numeric(length(observed)), names(observed))
  remaining <- target_ambient

  for (iteration in seq_len(max_iterations)) {
    capacity <- observed - allocation
    eligible <- capacity > tolerance & soup_profile > 0

    if (remaining <= tolerance || !any(eligible)) {
      break
    }

    proposal <- remaining *
      soup_profile[eligible] / sum(soup_profile[eligible])
    addition <- pmin(proposal, capacity[eligible])

    if (sum(addition) <= tolerance) {
      break
    }

    allocation[eligible] <- allocation[eligible] + addition
    remaining <- remaining - sum(addition)
  }

  attr(allocation, "unallocated_target") <- remaining
  allocation
}


# SoupX performs the subtraction on aggregated clusters to reduce ambiguity from
# single-cell sparsity. This teaching version follows that principle, then
# distributes each gene's inferred cluster-level contamination back to cells in
# proportion to their observed counts for that gene.
correct_counts_by_cluster <- function(cell_counts,
                                      cluster_labels,
                                      soup_profile,
                                      rho) {
  ambient_matrix <- matrix(
    0,
    nrow = nrow(cell_counts),
    ncol = ncol(cell_counts),
    dimnames = dimnames(cell_counts)
  )

  for (k in unique(cluster_labels)) {
    cells <- names(cluster_labels)[cluster_labels == k]
    cluster_observed <- rowSums(cell_counts[, cells, drop = FALSE])
    cluster_target <- rho * sum(cluster_observed)

    cluster_ambient <- allocate_ambient_counts(
      observed = cluster_observed,
      target_ambient = cluster_target,
      soup_profile = soup_profile
    )

    # Expand each gene's cluster-level ambient estimate back to cells. The
    # production SoupX expansion is more sophisticated; proportional expansion
    # is used here because every arithmetic step remains easy to inspect.
    for (gene in rownames(cell_counts)) {
      gene_total <- sum(cell_counts[gene, cells])
      if (gene_total > 0 && cluster_ambient[gene] > 0) {
        ambient_matrix[gene, cells] <-
          cluster_ambient[gene] * cell_counts[gene, cells] / gene_total
      }
    }
  }

  corrected <- cell_counts - ambient_matrix
  corrected[corrected < 1e-12] <- 0

  list(
    corrected = corrected,
    ambient = ambient_matrix
  )
}

correction <- correct_counts_by_cluster(
  filtered_umi,
  clusters,
  soup_profile,
  rho_hat
)

cat("\nSTEP 5: inferred ambient counts removed from each cell\n")
print(round(correction$ambient, 2))

cat("\nSTEP 6: corrected expression matrix\n")
print(round(correction$corrected, 2))


# =============================================================================
# PART 6 — SIMPLE DIAGNOSTICS
# =============================================================================

removed_by_gene <- sort(rowSums(correction$ambient), decreasing = TRUE)
fraction_removed <- sum(correction$ambient) / sum(filtered_umi)

cat("\nDIAGNOSTICS\n")
cat(sprintf("Fraction of cell-containing UMIs removed: %.3f\n", fraction_removed))
cat("Total inferred ambient counts by gene:\n")
print(round(removed_by_gene, 2))


# =============================================================================
# USING YOUR OWN MATRICES
# =============================================================================
#
# Replace raw_umi and filtered_umi at the top of this script. Then provide:
#
#   clusters <- setNames(your_cluster_vector, colnames(filtered_umi))
#
# and a marker-owner map such as:
#
#   marker_owner <- c(
#     HBB   = NA_character_,  # absent from every sampled cluster
#     CD3D  = "T",
#     MS4A1 = "B",
#     LST1  = "Mono"
#   )
#
# Requirements:
#   - raw_umi and filtered_umi must have the same gene rows;
#   - raw_umi must include low-UMI/empty droplets;
#   - filtered_umi columns must be the called cells;
#   - names(clusters) must equal colnames(filtered_umi);
#   - each sample or capture channel should be treated separately.
#
# For a production analysis, use the actual SoupX package after understanding
# these steps. Its optimized implementation is safer for sparse, large datasets.
