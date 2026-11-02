# =============================================================================
# 02_GGMs.R
#
# Cross-sectional Gaussian Graphical Models at T0, T1 and T2: estimation,
# stability, recovery simulations, strength and bridge strength, and paired
# Network Comparison Tests between waves.
#
#   source("00_setup.R")
#   source("02_GGMs.R")
#
# Inputs
#   data/residuals.csv
#
# Outputs (-> manuscript)
#   results/figures/ggm_T0.png, ggm_T1.png, ggm_T2.png        -> Figure 1
#   results/nct_global_invariance.csv                         -> Figure 1 (S, M, Holm p)
#   results/figures/strength_T0/T1/T2.png                     -> Figure 2a-c
#   results/figures/bridge_strength_T0/T1/T2.png              -> Figure 2d-f
#   results/centrality_T0_T1_T2.csv                           -> Supplementary Table 2
#   results/bridge_centrality_T0_T1_T2.csv                    -> Supplementary Table 3
#   results/figures/ggm_<wave>_case_stability.png             -> Supplementary Figures 2, 4, 6
#   results/figures/ggm_<wave>_simulation.png                 -> Supplementary Figures 3, 5, 7
#   results/ggm_cor_stability.csv                             CS-coefficients
#   results/ggm_density.csv                                   non-zero edges per wave
#   results/nct_edge_invariance.csv                           edge-level NCT (not reported)
#   results/nct_centrality_invariance.csv                     centrality NCT (not reported)
#
# Runtime: the bootstrap and simulation steps dominate (about 1-2 hours with
# N_BOOT_NET = 1000 and N_CORES = 4).
# =============================================================================

stopifnot(exists("NODES"))   # 00_setup.R must be sourced first


# -----------------------------------------------------------------------------
# 1. Wave-specific data sets
# -----------------------------------------------------------------------------
net_impute  <- read.csv(RESIDUALS_FILE, row.names = 1, check.names = FALSE)
analytic_id <- rownames(net_impute)

df_T0 <- net_impute[, NODES, drop = FALSE]
df_T1 <- net_impute[, paste0("T1", NODES), drop = FALSE]
df_T2 <- net_impute[, paste0("T2", NODES), drop = FALSE]

colnames(df_T1) <- NODES
colnames(df_T2) <- NODES

rownames(df_T0) <- rownames(df_T1) <- rownames(df_T2) <- analytic_id

# Same participants, same order, same variables: precondition for paired NCT.
stopifnot(
  identical(rownames(df_T0), rownames(df_T1)),
  identical(rownames(df_T1), rownames(df_T2)),
  identical(colnames(df_T0), colnames(df_T1)),
  identical(colnames(df_T1), colnames(df_T2))
)

waves_data <- list(T0 = df_T0, T1 = df_T1, T2 = df_T2)


# -----------------------------------------------------------------------------
# 2. Network estimation
# -----------------------------------------------------------------------------
# Graphical LASSO with EBIC model selection (EBICglasso, gamma = 0.5).
nwk   <- lapply(waves_data, function(d) bootnet::estimateNetwork(d, default = "EBICglasso"))
edges <- lapply(nwk, bootnet::getWmat)

for (w in names(edges)) {
  write.csv(edges[[w]], file.path(RESULTS_DIR, paste0("ggm_edges_", w, ".csv")))
}

density_tbl <- purrr::imap_dfr(edges, function(m, w) {
  network_density(m, directed = FALSE) %>% mutate(Wave = w)
}) %>% select(Wave, everything())

write.csv(density_tbl, file.path(RESULTS_DIR, "ggm_density.csv"), row.names = FALSE)
print(density_tbl)


# -----------------------------------------------------------------------------
# 3. Network figures (Figure 1)
# -----------------------------------------------------------------------------
# T0 fixes the layout; T1 and T2 reuse it so that the three graphs are directly
# comparable.
plot_ggm <- function(net, title, layout_in) {
  qgraph::qgraph(
    net$graph,
    layout       = layout_in,
    labels       = LABELS_SHORT[NODES],
    nodeNames    = LABELS_FULL[NODES],
    color        = NODE_COLORS[NODES],
    border.color = NODE_BORDERS[NODES],
    border.width = 2,
    label.color  = NODE_LABEL_COLORS[NODES],
    groups       = NODE_GROUPS,
    legend       = TRUE,
    legend.mode  = "names",
    legend.cex   = 0.6,
    GLratio      = 3.2,
    label.cex    = 0.9,
    label.prop   = 1,
    negDashed    = FALSE,
    curveAll     = TRUE,
    curveDefault = 0.35,
    posCol       = "blue",
    negCol       = "red",
    vsize        = 9,
    font         = 2,
    title        = title
  )
}

png(file.path(FIG_DIR, "ggm_T0.png"), width = 12, height = 8, units = "in", res = 300)
plot_T0 <- plot_ggm(nwk$T0, "T0 - GGM (EBICglasso)", "spring")
dev.off()

png(file.path(FIG_DIR, "ggm_T1.png"), width = 12, height = 8, units = "in", res = 300)
plot_T1 <- plot_ggm(nwk$T1, "T1 (12 months) - GGM (EBICglasso)", plot_T0$layout)
dev.off()

png(file.path(FIG_DIR, "ggm_T2.png"), width = 12, height = 8, units = "in", res = 300)
plot_T2 <- plot_ggm(nwk$T2, "T2 (24 months) - GGM (EBICglasso)", plot_T0$layout)
dev.off()

plots <- list(T0 = plot_T0, T1 = plot_T1, T2 = plot_T2)


# -----------------------------------------------------------------------------
# 4. Stability and recovery
# -----------------------------------------------------------------------------
# Case-dropping bootstrap -> CS-coefficients. All indices are computed so that
# their CS-coefficients are on record; only strength and bridge strength are
# plotted and reported.
b_case <- lapply(names(nwk), function(w) {
  set.seed(SEED)
  bootnet::bootnet(
    nwk[[w]],
    nBoots      = N_BOOT_NET,
    type        = "case",
    nCores      = N_CORES,
    statistics  = c("strength", "closeness", "betweenness", "expectedInfluence",
                    "bridgeStrength", "bridgeExpectedInfluence"),
    communities = BRIDGE_COMMUNITIES
  )
})
names(b_case) <- names(nwk)

cs_tbl <- purrr::imap_dfr(b_case, function(b, w) {
  cs <- bootnet::corStability(b, verbose = FALSE)
  tibble::tibble(Wave = w, Statistic = names(cs), CS = as.numeric(cs))
})

write.csv(cs_tbl, file.path(RESULTS_DIR, "ggm_cor_stability.csv"), row.names = FALSE)
print(cs_tbl)

# Recovery simulations across sample sizes (netSimulator)
sim <- lapply(names(nwk), function(w) {
  set.seed(SEED)
  bootnet::netSimulator(
    input   = nwk[[w]],
    default = "EBICglasso",
    nCases  = SIM_CASES,
    nReps   = SIM_REPS
  )
})
names(sim) <- names(nwk)

saveRDS(list(case = b_case, simulation = sim),
        file.path(RESULTS_DIR, "ggm_stability_objects.rds"))

for (w in names(nwk)) {

  # Supplementary Figures 2, 4, 6
  png(file.path(FIG_DIR, paste0("ggm_", w, "_case_stability.png")),
      width = 9, height = 5, units = "in", res = 300)
  print(plot(b_case[[w]], statistics = c("strength", "bridgeStrength")))
  dev.off()

  # Supplementary Figures 3, 5, 7
  png(file.path(FIG_DIR, paste0("ggm_", w, "_simulation.png")),
      width = 10, height = 5, units = "in", res = 300)
  print(plot(sim[[w]]))
  dev.off()
}


# -----------------------------------------------------------------------------
# 5. Strength and bridge strength (Figure 2, Supplementary Tables 2-3)
# -----------------------------------------------------------------------------
centrality_all <- purrr::imap_dfr(nwk, function(net, w) centrality_table(net, w))

write.csv(centrality_all,
          file.path(RESULTS_DIR, "centrality_T0_T1_T2.csv"),
          row.names = FALSE)

bridges <- lapply(plots, function(p) {
  networktools::bridge(
    p,
    communities    = BRIDGE_COMMUNITIES,
    useCommunities = "all",
    directed       = NULL,
    nodes          = NULL
  )
})

bridge_all <- purrr::imap_dfr(bridges, function(br, w) bridge_table(br, w))

write.csv(bridge_all,
          file.path(RESULTS_DIR, "bridge_centrality_T0_T1_T2.csv"),
          row.names = FALSE)

# One panel per wave and index. Values are z-standardised WITHIN each wave, so
# each panel shows the relative ranking of nodes at that assessment; absolute
# values are in the CSV files above.
for (w in names(nwk)) {

  p_str <- centrality_all %>%
    filter(Wave == w) %>%
    transmute(Node, Index = "Strength", Value = z_score(Strength)) %>%
    centrality_panel()

  ggsave(file.path(FIG_DIR, paste0("strength_", w, ".png")),
         p_str, width = 5, height = 5, dpi = 300)

  p_brg <- bridge_all %>%
    filter(Wave == w) %>%
    transmute(Node, Index = "Bridge strength", Value = z_score(BridgeStrength)) %>%
    centrality_panel()

  ggsave(file.path(FIG_DIR, paste0("bridge_strength_", w, ".png")),
         p_brg, width = 5, height = 5, dpi = 300)
}

# Cross-wave overlays, same style
p_str_all <- centrality_all %>%
  group_by(Wave) %>%
  transmute(Wave, Node, Index = "Strength", Value = z_score(Strength)) %>%
  ungroup() %>%
  centrality_panel(colour_by = "Wave")

ggsave(file.path(FIG_DIR, "strength_across_waves.png"),
       p_str_all, width = 6, height = 5, dpi = 300)

p_brg_all <- bridge_all %>%
  group_by(Wave) %>%
  transmute(Wave, Node, Index = "Bridge strength", Value = z_score(BridgeStrength)) %>%
  ungroup() %>%
  centrality_panel(colour_by = "Wave")

ggsave(file.path(FIG_DIR, "bridge_strength_across_waves.png"),
       p_brg_all, width = 6, height = 5, dpi = 300)


# -----------------------------------------------------------------------------
# 6. Paired Network Comparison Tests
# -----------------------------------------------------------------------------
# T0, T1 and T2 are the same participants, so relabelling is performed within
# each pair of observations (paired = TRUE). The package is asked for RAW
# p-values (p.adjust.methods = "none"); multiplicity corrections are applied
# once, below: Holm across the three temporal contrasts for the global tests,
# Benjamini-Hochberg within each contrast for edges and centrality.

nct_centrality_indices <- c("strength", "expectedInfluence",
                            "bridgeStrength", "bridgeExpectedInfluence")

run_nct <- function(net_a, net_b) {
  set.seed(SEED)
  NetworkComparisonTest::NCT(
    net_a, net_b,
    it               = N_NCT,
    paired           = TRUE,
    test.edges       = TRUE,
    edges            = "all",
    p.adjust.methods = "none",
    test.centrality  = TRUE,
    centrality       = nct_centrality_indices,
    nodes            = "all",
    communities      = BRIDGE_COMMUNITIES,
    progressbar      = FALSE,
    verbose          = FALSE
  )
}

nct_list <- list(
  "T0 vs T1" = run_nct(nwk$T0, nwk$T1),
  "T1 vs T2" = run_nct(nwk$T1, nwk$T2),
  "T0 vs T2" = run_nct(nwk$T0, nwk$T2)
)

saveRDS(nct_list, file.path(RESULTS_DIR, "nct_paired_GGM_T0_T1_T2.rds"))


## 6a. Global invariance (reported in Figure 1 and Results) --------------------
# S (global strength difference) and M (maximum edge difference) are
# magnitudes, so their permutation p-values are one-sided. The package values
# are cross-checked against the stored permutation distributions.
perm_p_onesided <- function(observed, permuted) {
  permuted <- permuted[is.finite(permuted)]
  if (!length(permuted) || !is.finite(observed)) return(NA_real_)
  (sum(permuted >= observed) + 1) / (length(permuted) + 1)
}

nct_global <- purrr::imap_dfr(nct_list, function(res, cmp) {
  tibble::tibble(
    Comparison          = cmp,
    GlobalStrength_net1 = res$glstrinv.sep[1],
    GlobalStrength_net2 = res$glstrinv.sep[2],
    S                   = res$glstrinv.real,
    p_S                 = res$glstrinv.pval,
    p_S_recomputed      = perm_p_onesided(res$glstrinv.real, as.numeric(res$glstrinv.perm)),
    M                   = res$nwinv.real,
    p_M                 = res$nwinv.pval,
    p_M_recomputed      = perm_p_onesided(res$nwinv.real, as.numeric(res$nwinv.perm))
  )
}) %>%
  mutate(
    # Two families of three hypotheses each (S and M), Holm within each family
    p_S_holm = p.adjust(p_S, method = "holm"),
    p_M_holm = p.adjust(p_M, method = "holm")
  )

write.csv(nct_global,
          file.path(RESULTS_DIR, "nct_global_invariance.csv"),
          row.names = FALSE)
print(nct_global)

# With N_NCT permutations the smallest attainable p-value is 1 / (N_NCT + 1);
# report anything at that floor as p < 0.001 rather than p = 0.


## 6b. Edge-level invariance ----------------------------------------------------
# Not reported in the manuscript: following the omnibus logic of the NCT, edge-
# level differences are not interpreted when the structure test (M) is not
# significant after Holm correction. The table is written for completeness.
#
# The edge p-values are read from einv.pvals BY COLUMN NAME and matched BY NODE
# PAIR. Selecting the last column by position is unsafe: in the installed
# version the last column of einv.pvals is the edge-invariance statistic E, not
# the p-value.

edge_index <- function(labs) {
  k   <- length(labs)
  idx <- which(upper.tri(matrix(NA_real_, k, k)), arr.ind = TRUE)
  tibble::tibble(
    i     = as.integer(idx[, "row"]),
    j     = as.integer(idx[, "col"]),
    Node1 = labs[idx[, "row"]],
    Node2 = labs[idx[, "col"]]
  ) %>%
    mutate(key = paste(pmin(Node1, Node2), pmax(Node1, Node2), sep = "|"))
}

observed_edge_diff <- function(res, labs) {
  ei <- edge_index(labs)
  er <- res$einv.real
  if ((is.matrix(er) || is.data.frame(er)) && all(dim(as.matrix(er)) == length(labs))) {
    ei$Difference <- as.matrix(er)[cbind(ei$i, ei$j)]
  } else if (!is.null(er) && length(as.numeric(er)) == nrow(ei)) {
    ei$Difference <- as.numeric(er)
  } else {
    ei$Difference <- NA_real_
  }
  ei
}

# Orient einv.perm as edges x permutations, if its shape is recognised
align_perm <- function(perm, n_edges, n_it) {
  if (is.null(perm)) return(NULL)
  if (is.list(perm) && !is.data.frame(perm)) {
    perm <- tryCatch(do.call(cbind, lapply(perm, as.numeric)), error = function(e) NULL)
    if (is.null(perm)) return(NULL)
  }
  p <- tryCatch(as.matrix(perm), error = function(e) NULL)
  if (is.null(p) || !is.numeric(p) || length(dim(p)) != 2) return(NULL)
  if (nrow(p) == n_edges && ncol(p) != n_edges) return(p)
  if (ncol(p) == n_edges && nrow(p) != n_edges) return(t(p))
  if (nrow(p) == n_it) return(t(p))
  if (ncol(p) == n_it) return(p)
  NULL
}

# Edge differences are signed, so the edge-level test is two-sided
perm_p_twosided <- function(observed, permuted) {
  permuted <- permuted[is.finite(permuted)]
  if (!length(permuted) || !is.finite(observed)) return(NA_real_)
  (sum(abs(permuted) >= abs(observed)) + 1) / (length(permuted) + 1)
}

package_edge_pvals <- function(res, labs) {
  k     <- length(labs)
  empty <- tibble::tibble(key = character(), p_package = numeric())
  ep    <- res$einv.pvals
  if (is.null(ep)) return(empty)
  ep <- as.data.frame(ep)
  if (ncol(ep) < 3) return(empty)

  to_label <- function(v) {
    if (is.factor(v)) v <- as.character(v)
    if (is.numeric(v) && all(v %in% seq_len(k), na.rm = TRUE)) labs[v] else as.character(v)
  }

  num_cols <- setdiff(which(vapply(ep, is.numeric, logical(1))), 1:2)
  p_col    <- num_cols[grepl("p", tolower(names(ep)[num_cols]), fixed = TRUE)]
  if (!length(p_col)) return(empty)

  tibble::tibble(
    Node1     = to_label(ep[[1]]),
    Node2     = to_label(ep[[2]]),
    p_package = as.numeric(ep[[p_col[1]]])
  ) %>%
    mutate(key = paste(pmin(Node1, Node2), pmax(Node1, Node2), sep = "|")) %>%
    select(key, p_package)
}

nct_edge_table <- function(res, comparison, labs, n_it) {
  ei <- observed_edge_diff(res, labs)
  pm <- align_perm(res$einv.perm, nrow(ei), n_it)

  ei$p_recomputed <- if (!is.null(pm)) {
    vapply(seq_len(nrow(ei)),
           function(e) perm_p_twosided(ei$Difference[e], pm[e, ]),
           numeric(1))
  } else NA_real_

  left_join(ei, package_edge_pvals(res, labs), by = "key") %>%
    mutate(
      p_raw      = dplyr::coalesce(p_recomputed, p_package),
      Comparison = comparison
    )
}

nct_edges <- purrr::imap_dfr(nct_list, ~ nct_edge_table(.x, .y, colnames(df_T0), N_NCT)) %>%
  group_by(Comparison) %>%
  mutate(p_BH = if (any(is.finite(p_raw))) p.adjust(p_raw, method = "BH") else NA_real_) %>%
  ungroup() %>%
  mutate(
    Edge       = paste(LABELS_SHORT[Node1], "--", LABELS_SHORT[Node2]),
    Comparison = factor(Comparison, levels = names(nct_list))
  ) %>%
  select(Comparison, Edge, Difference, p_recomputed, p_package, p_raw, p_BH) %>%
  arrange(Comparison, p_BH)

write.csv(nct_edges,
          file.path(RESULTS_DIR, "nct_edge_invariance.csv"),
          row.names = FALSE)

# Cross-check: if both sources are available they must agree
edge_check <- nct_edges %>%
  group_by(Comparison) %>%
  summarise(
    n_recomputed = sum(is.finite(p_recomputed)),
    n_package    = sum(is.finite(p_package)),
    max_abs_diff = safe_max(abs(p_recomputed - p_package)),
    spearman_r   = safe_cor(p_recomputed, p_package),
    .groups = "drop"
  )
print(edge_check)


## 6c. Centrality invariance ---------------------------------------------------
# Not reported in the manuscript; written for completeness. diffcen.pval is raw
# (p.adjust.methods = "none"); BH is applied within each contrast and index.
nct_cen_table <- function(res, comparison) {
  real <- as.data.frame(res$diffcen.real)
  pval <- as.data.frame(res$diffcen.pval)
  stopifnot(identical(dim(real), dim(pval)))
  real$Node <- rownames(real)
  pval$Node <- rownames(pval)
  left_join(
    pivot_longer(real, -Node, names_to = "Measure", values_to = "Difference"),
    pivot_longer(pval, -Node, names_to = "Measure", values_to = "p_raw"),
    by = c("Node", "Measure")
  ) %>%
    mutate(Comparison = comparison)
}

nct_centrality <- purrr::imap_dfr(nct_list, ~ nct_cen_table(.x, .y)) %>%
  group_by(Comparison, Measure) %>%
  mutate(p_BH = if (any(is.finite(p_raw))) p.adjust(p_raw, method = "BH") else NA_real_) %>%
  ungroup() %>%
  mutate(
    Node       = recode(as.character(Node), !!!as.list(LABELS_SHORT)),
    Measure    = recode(Measure,
                        strength                = "Strength",
                        expectedInfluence       = "Expected influence",
                        bridgeStrength          = "Bridge strength",
                        bridgeExpectedInfluence = "Bridge expected influence",
                        .default = Measure),
    Comparison = factor(Comparison, levels = names(nct_list))
  ) %>%
  select(Comparison, Node, Measure, Difference, p_raw, p_BH) %>%
  arrange(Comparison, Measure, p_BH)

write.csv(nct_centrality,
          file.path(RESULTS_DIR, "nct_centrality_invariance.csv"),
          row.names = FALSE)


## 6d. Permutation distributions ------------------------------------------------
png(file.path(FIG_DIR, "nct_permutation_strength.png"),
    width = 12, height = 4, units = "in", res = 300)
par(mfrow = c(1, 3))
for (cmp in names(nct_list)) {
  plot(nct_list[[cmp]], what = "strength")
  title(sub = cmp)
}
par(mfrow = c(1, 1))
dev.off()

png(file.path(FIG_DIR, "nct_permutation_network.png"),
    width = 12, height = 4, units = "in", res = 300)
par(mfrow = c(1, 3))
for (cmp in names(nct_list)) {
  plot(nct_list[[cmp]], what = "network")
  title(sub = cmp)
}
par(mfrow = c(1, 1))
dev.off()


message("02_GGMs.R complete.")
