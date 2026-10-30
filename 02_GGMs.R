# =============================================================================
# 02_GGMs.R
#
# Cross-sectional Gaussian Graphical Models at T0, T1 and T2, with accuracy and
# stability analyses, standard and bridge centrality, and paired Network
# Comparison Tests between waves.
#
#   source("R/00_setup.R")
#   source("R/02_GGMs.R")
#
# Inputs
#   data/residuals.csv
#
# Outputs
#   results/ggm_edges_<wave>.csv            weight matrices
#   results/ggm_density.csv
#   results/ggm_cor_stability.csv           CS-coefficients
#   results/centrality_T0_T1_T2.csv
#   results/bridge_centrality_T0_T1_T2.csv
#   results/nct_global_invariance.csv
#   results/nct_edge_invariance.csv
#   results/nct_centrality_invariance.csv
#   results/nct_paired_GGM_T0_T1_T2.rds
#   results/figures/ggm_<wave>.png and accompanying diagnostic figures
#
# Runtime note: the bootstrap and simulation steps dominate. With N_BOOT_NET =
# 1000 and N_CORES = 4 expect roughly 1-2 hours in total on a laptop.
# =============================================================================

stopifnot(exists("NODES"))


# -----------------------------------------------------------------------------
# 1. Wave-specific data sets
# -----------------------------------------------------------------------------
net_impute <- read.csv(RESIDUALS_FILE, row.names = 1, check.names = FALSE)

analytic_id <- rownames(net_impute)

df_T0 <- net_impute[, NODES, drop = FALSE]
df_T1 <- net_impute[, paste0("T1", NODES), drop = FALSE]
df_T2 <- net_impute[, paste0("T2", NODES), drop = FALSE]

colnames(df_T1) <- NODES
colnames(df_T2) <- NODES

rownames(df_T0) <- rownames(df_T1) <- rownames(df_T2) <- analytic_id

# The three data sets must contain the same participants, in the same order,
# with the same variables. This is a precondition for paired = TRUE in NCT().
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
# Regularised partial-correlation networks: graphical LASSO with EBIC model
# selection (EBICglasso, tuning gamma = 0.5, the bootnet default).

nwk <- lapply(waves_data, function(d) {
  bootnet::estimateNetwork(d, default = "EBICglasso")
})

edges <- lapply(nwk, bootnet::getWmat)

for (w in names(edges)) {
  write.csv(edges[[w]], file.path(RESULTS_DIR, paste0("ggm_edges_", w, ".csv")))
}

density_tbl <- purrr::imap_dfr(edges, function(m, w) {
  network_density(m, directed = FALSE) %>% mutate(Wave = w)
}) %>%
  select(Wave, everything())

write.csv(density_tbl, file.path(RESULTS_DIR, "ggm_density.csv"), row.names = FALSE)
print(density_tbl)


# -----------------------------------------------------------------------------
# 3. Network figures
# -----------------------------------------------------------------------------
# T0 fixes the layout; T1 and T2 reuse it so that the three graphs are visually
# comparable.

plot_ggm <- function(net, title, layout_in) {
  qgraph::qgraph(
    net$graph,
    layout       = layout_in,
    labels       = LABELS_SHORT[NODES],
    nodeNames    = LABELS_FULL[NODES],
    color        = NODE_COLORS[NODES],
    groups       = NODE_GROUPS,
    legend       = TRUE,
    legend.mode  = "names",
    legend.cex   = 0.6,
    GLratio      = 3.2,
    label.cex    = 0.9,
    label.color  = "black",
    label.prop   = 1,
    border.width = 1,
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
# 4. Accuracy and stability
# -----------------------------------------------------------------------------
# (a) Case-dropping subset bootstrap -> correlation-stability coefficient
# (b) Non-parametric bootstrap       -> edge-weight accuracy, difference tests
# (c) netSimulator                   -> structure recovery across sample sizes

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

b_np <- lapply(names(nwk), function(w) {
  set.seed(SEED)
  bootnet::bootnet(
    nwk[[w]],
    nBoots     = N_BOOT_NET,
    type       = "nonparametric",
    nCores     = N_CORES,
    statistics = c("edge", "strength", "closeness",
                   "betweenness", "expectedInfluence")
  )
})
names(b_np) <- names(nwk)

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

saveRDS(list(case = b_case, nonparametric = b_np, simulation = sim),
        file.path(RESULTS_DIR, "ggm_bootstrap_objects.rds"))

# Diagnostic figures
for (w in names(nwk)) {

  png(file.path(FIG_DIR, paste0("ggm_", w, "_case_stability.png")),
      width = 9, height = 5, units = "in", res = 300)
  print(plot(b_case[[w]], statistics = c("strength", "closeness",
                                         "betweenness", "expectedInfluence")))
  dev.off()

  png(file.path(FIG_DIR, paste0("ggm_", w, "_edge_accuracy.png")),
      width = 8, height = 8, units = "in", res = 300)
  print(plot(b_np[[w]], labels = TRUE, order = "sample"))
  dev.off()

  png(file.path(FIG_DIR, paste0("ggm_", w, "_simulation.png")),
      width = 10, height = 5, units = "in", res = 300)
  print(plot(sim[[w]]))
  dev.off()

  # Bootstrapped centrality difference tests. These, not the point estimates,
  # are the basis for any claim that one node is more central than another.
  for (stat in c("strength", "expectedInfluence")) {
    png(file.path(FIG_DIR, paste0("ggm_", w, "_difftest_", stat, ".png")),
        width = 7, height = 6, units = "in", res = 300)
    print(plot(b_np[[w]], statistics = stat, plot = "difference",
               order = "sample", labels = TRUE))
    dev.off()
  }
}


# -----------------------------------------------------------------------------
# 5. Standard centrality
# -----------------------------------------------------------------------------
centrality_all <- purrr::imap_dfr(nwk, function(net, w) centrality_table(net, w))

write.csv(centrality_all,
          file.path(RESULTS_DIR, "centrality_T0_T1_T2.csv"),
          row.names = FALSE)

png(file.path(FIG_DIR, "centrality_across_waves.png"),
    width = 10, height = 7, units = "in", res = 300)
print(qgraph::centralityPlot(
  lapply(nwk, net_labelled),
  include = c("Strength", "Closeness", "Betweenness", "ExpectedInfluence"),
  scale   = "z-scores",
  orderBy = "Strength"
))
dev.off()

# Reporting note. Strength and expected influence are sums of edge weights, so
# their stability follows directly from edge stability; in a network containing
# negative edges, expected influence is the more interpretable of the two.
# Closeness and betweenness rely on shortest paths, which have no clear meaning
# in a partial-correlation network, and are the least stable indices: report
# them only if their CS-coefficients in ggm_cor_stability.csv reach 0.25 at a
# minimum, 0.50 preferably.


# -----------------------------------------------------------------------------
# 6. Bridge centrality
# -----------------------------------------------------------------------------
# Two communities: psychopathology (5 PANSS factors) and functioning (4 HoNOS
# items). Bridge indices isolate the portion of a node's connectivity that
# crosses that boundary, and are the indices most directly tied to the aim of
# the study.

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

bridge_fig <- bridge_all %>%
  select(Wave, Node, BridgeStrength, BridgeEI1) %>%
  pivot_longer(c(BridgeStrength, BridgeEI1),
               names_to = "Index", values_to = "Value") %>%
  mutate(
    Wave  = factor(Wave, levels = c("T0", "T1", "T2")),
    Node  = factor(Node, levels = LABELS_SHORT[NODES]),
    Index = recode(Index,
                   BridgeStrength = "Bridge strength",
                   BridgeEI1      = "Bridge expected influence (1-step)")
  ) %>%
  ggplot(aes(x = Node, y = Value, colour = Wave, group = Wave)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  facet_wrap(~ Index, scales = "free_y") +
  labs(title = "Bridge centrality across waves",
       x = "Node", y = "Value", colour = "Wave") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "bridge_centrality_across_waves.png"),
       bridge_fig, width = 11, height = 6, dpi = 300)


# -----------------------------------------------------------------------------
# 7. Paired Network Comparison Tests
# -----------------------------------------------------------------------------
# T0, T1 and T2 are the same participants measured three times. The default
# permutation scheme of NCT (paired = FALSE) reshuffles observations freely
# between samples and is invalid for this design: it would treat repeated
# measurements as if they came from different people, and the resulting
# p-values would be anticonservative. With paired = TRUE, relabelling is
# performed within each pair of observations.
#
# The package documentation states that dependent data are currently assumed to
# entail one group measured twice. Three pairwise comparisons across three waves
# is a reasonable extension, but the three tests are not independent: they are
# adjusted jointly below and interpreted as a set.

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
    p.adjust.methods = "fdr",
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


## 7a. Global invariance -------------------------------------------------------
nct_global <- purrr::imap_dfr(nct_list, function(res, cmp) {
  tibble::tibble(
    Comparison             = cmp,
    GlobalStrength_net1    = res$glstrinv.sep[1],
    GlobalStrength_net2    = res$glstrinv.sep[2],
    S_globalStrengthDiff   = res$glstrinv.real,
    p_globalStrength       = res$glstrinv.pval,
    M_maxEdgeDiff          = res$nwinv.real,
    p_structure            = res$nwinv.pval
  )
}) %>%
  mutate(
    p_globalStrength_holm = p.adjust(p_globalStrength, method = "holm"),
    p_structure_holm      = p.adjust(p_structure,      method = "holm")
  )

write.csv(nct_global,
          file.path(RESULTS_DIR, "nct_global_invariance.csv"),
          row.names = FALSE)
print(nct_global)

# A non-significant M means no single edge differs beyond what permutation would
# produce. It does NOT establish that the two networks are identical: absence of
# evidence for a difference is not evidence of invariance, and NCT has limited
# power at these sample sizes. A significant S with a non-significant M points to
# a diffuse change in overall connectivity rather than a localised one.


## 7b. Edge-level invariance ---------------------------------------------------
nct_edge_table <- function(res, comparison, labs) {

  tryCatch({

    ep <- as.data.frame(res$einv.pvals)
    names(ep)[ncol(ep)] <- "p_fdr"
    names(ep)[1:2] <- c("Var1", "Var2")

    # Coerce node identifiers to character unconditionally. If they stayed
    # factors, LABELS_SHORT[Var1] would index by level code and silently
    # produce wrong labels.
    ep <- ep %>% mutate(Var1 = as.character(Var1), Var2 = as.character(Var2))

    er <- res$einv.real

    if (is.matrix(er) || is.data.frame(er)) {
      er <- as.matrix(er)
      dimnames(er) <- list(labs, labs)
      er_long <- er %>%
        as.data.frame() %>%
        tibble::rownames_to_column("Var1") %>%
        pivot_longer(-Var1, names_to = "Var2", values_to = "Difference")
      ep <- left_join(ep, er_long, by = c("Var1", "Var2"))
    }

    ep %>% mutate(Comparison = comparison)

  }, error = function(e) {
    message("Edge extraction failed for ", comparison, ": ", e$message)
    NULL
  })
}

nct_edges <- purrr::imap_dfr(nct_list, ~ nct_edge_table(.x, .y, colnames(df_T0)))

if (nrow(nct_edges) > 0) {
  nct_edges <- nct_edges %>%
    mutate(Edge = paste(LABELS_SHORT[Var1], "--", LABELS_SHORT[Var2])) %>%
    select(Comparison, Edge, any_of("Difference"), p_fdr) %>%
    arrange(p_fdr)

  write.csv(nct_edges,
            file.path(RESULTS_DIR, "nct_edge_invariance.csv"),
            row.names = FALSE)
}


## 7c. Centrality invariance, standard and bridge ------------------------------
nct_cen_table <- function(res, comparison) {

  tryCatch({

    real <- as.data.frame(res$diffcen.real)
    pval <- as.data.frame(res$diffcen.pval)

    real$Node <- rownames(real)
    pval$Node <- rownames(pval)

    left_join(
      pivot_longer(real, -Node, names_to = "Measure", values_to = "Difference"),
      pivot_longer(pval, -Node, names_to = "Measure", values_to = "p_fdr"),
      by = c("Node", "Measure")
    ) %>%
      mutate(Comparison = comparison)

  }, error = function(e) {
    message("Centrality extraction failed for ", comparison, ": ", e$message)
    NULL
  })
}

nct_centrality <- purrr::imap_dfr(nct_list, ~ nct_cen_table(.x, .y))

if (nrow(nct_centrality) > 0) {

  nct_centrality <- nct_centrality %>%
    mutate(
      Node     = as.character(Node),
      Node     = recode(Node, !!!as.list(LABELS_SHORT)),
      IsBridge = grepl("^bridge", Measure, ignore.case = TRUE),
      Measure  = recode(
        Measure,
        strength                = "Strength",
        expectedInfluence       = "Expected influence",
        bridgeStrength          = "Bridge strength",
        bridgeExpectedInfluence = "Bridge expected influence",
        .default = Measure
      ),
      Comparison = factor(Comparison, levels = names(nct_list))
    ) %>%
    select(Comparison, Node, Measure, IsBridge, Difference, p_fdr) %>%
    arrange(p_fdr)

  write.csv(select(nct_centrality, -IsBridge),
            file.path(RESULTS_DIR, "nct_centrality_invariance.csv"),
            row.names = FALSE)

  # Differences are computed as first network minus second: in "T0 vs T1" a
  # positive value means the node bridged the two communities more strongly at
  # T0 than at T1.
  bridge_diff_plot <- nct_centrality %>%
    filter(IsBridge) %>%
    mutate(
      Node        = factor(Node, levels = LABELS_SHORT[NODES]),
      Significant = if_else(p_fdr < .05, "FDR p < .05", "n.s.")
    ) %>%
    ggplot(aes(x = Node, y = Difference, fill = Significant)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    facet_grid(Measure ~ Comparison, scales = "free_y") +
    scale_fill_manual(values = c("FDR p < .05" = "#C51B7D", "n.s." = "grey70")) +
    labs(title = "Between-wave differences in bridge centrality",
         x = "Node", y = "Difference (net1 - net2)", fill = "") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")

  ggsave(file.path(FIG_DIR, "nct_bridge_centrality_differences.png"),
         bridge_diff_plot, width = 11, height = 6, dpi = 300)
}


## 7d. Permutation distributions -----------------------------------------------
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

# The NCT compares structure, not severity: because EBICglasso operates on the
# correlation matrix, a wave in which every domain improves in absolute terms can
# still produce an identical network. Symptomatic change over follow-up is
# described by the Friedman tests in 01_preprocessing.R, not here.


message("02_GGMs.R complete.")
