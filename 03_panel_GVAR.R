# =============================================================================
# 03_panel_GVAR.R
#
# Panel graphical vector-autoregressive model over the three waves: model
# estimation, model selection, extraction of the temporal / contemporaneous /
# between-subjects networks, and temporal centrality.
#
#   source("R/00_setup.R")
#   source("R/03_panel_GVAR.R")
#
# Inputs
#   data/residuals.csv
#
# Outputs
#   results/panelgvar_fit.csv
#   results/panelgvar_model_selection.csv
#   results/panelgvar_model_final.rds
#   results/panelgvar_temporal_PDC_full.csv
#   results/panelgvar_temporal_PDC_crosslagged.csv
#   results/panelgvar_autoregressions.csv
#   results/panelgvar_temporal_edges_retained.csv
#   results/panelgvar_temporal_edge_types.csv
#   results/panelgvar_temporal_strength.csv
#   results/panelgvar_contemporaneous.csv
#   results/panelgvar_between.csv
#   results/panelgvar_undirected_edges.csv
#   results/panelgvar_density.csv
#   results/figures/panelgvar_*.png
# =============================================================================

stopifnot(exists("NODES"))


# -----------------------------------------------------------------------------
# 1. Data
# -----------------------------------------------------------------------------
net_impute <- read.csv(RESIDUALS_FILE, row.names = 1, check.names = FALSE)

data_pgvar <- net_impute[, NODE_VARS, drop = FALSE] %>% as_tibble()

stopifnot(ncol(data_pgvar) == length(NODE_VARS))

message("Panel-GVAR input: ", nrow(data_pgvar), " participants x ",
        ncol(data_pgvar), " variables")


# -----------------------------------------------------------------------------
# 2. Model estimation
# -----------------------------------------------------------------------------
# psychonetrics::dlvm1 with FIML, a within-person GGM and a Cholesky
# parameterisation of the between-person structure (Epskamp, 2020,
# Psychometrika). The model decomposes the covariance structure into:
#   - a temporal network        : directed lag-1 effects across adjacent waves
#   - a contemporaneous network : within-person partial correlations at a wave,
#                                 after accounting for temporal effects
#   - a between-subjects network: stable individual differences in mean level

model <- psychonetrics::dlvm1(
  data           = data_pgvar,
  vars           = DESIGN_MAT,
  estimator      = "FIML",
  within_latent  = "ggm",
  between_latent = "chol"
)

set.seed(SEED)

model_unconstrained <- psychonetrics::runmodel(model)

model_pruned <- model_unconstrained %>%
  psychonetrics::prune(alpha = PRUNE_ALPHA, recursive = TRUE) %>%
  psychonetrics::stepup(alpha = STEPUP_ALPHA, verbose = FALSE)

# Terminology. The unconstrained model is NOT saturated: dlvm1 imposes
# stationarity (the same beta and the same within-person GGM for the T0->T1 and
# T1->T2 transitions), which consumes degrees of freedom. The remaining df are
# reported in panelgvar_fit.csv. "Unconstrained" or "fully connected" is the
# correct label throughout.


# -----------------------------------------------------------------------------
# 3. Fit and model selection
# -----------------------------------------------------------------------------
fit_unc <- psychonetrics::fit(model_unconstrained)
fit_spa <- psychonetrics::fit(model_pruned)

fit_tbl <- as.data.frame(fit_unc) %>%
  rename(Unconstrained = Value) %>%
  left_join(as.data.frame(fit_spa) %>% rename(Sparse = Value), by = "Measure")

write.csv(fit_tbl, file.path(RESULTS_DIR, "panelgvar_fit.csv"), row.names = FALSE)
print(fit_tbl)

npar_unc <- fit_value(fit_unc, "npar");   df_unc  <- fit_value(fit_unc, "df")
npar_spa <- fit_value(fit_spa, "npar");   df_spa  <- fit_value(fit_spa, "df")
bic_unc  <- fit_value(fit_unc, "bic");    bic_spa <- fit_value(fit_spa, "bic")
aic_unc  <- fit_value(fit_unc, "aic.ll"); aic_spa <- fit_value(fit_spa, "aic.ll")

cmp <- psychonetrics::compare(
  unconstrained = model_unconstrained,
  sparse        = model_pruned
)
print(cmp)

cmp_df       <- as.data.frame(cmp)
p_col        <- grep("^p", names(cmp_df), value = TRUE, ignore.case = TRUE)[1]
chisq_diff_p <- suppressWarnings(as.numeric(cmp_df[[p_col]]))
chisq_diff_p <- chisq_diff_p[!is.na(chisq_diff_p)]
chisq_diff_p <- if (length(chisq_diff_p) == 0) NA_real_ else chisq_diff_p[1]

# Selection rule. The sparse model is retained when BIC favours it. BIC
# penalises the parameters that pruning removes and is the criterion least
# affected by the fact that pruning was carried out on these same data. The
# chi-square difference test is reported descriptively alongside it, and can be
# significant even when every parsimony-adjusted criterion favours the sparse
# model; because pruning and step-up were data-driven on the same data, that
# test is biased IN FAVOUR of the sparse model, so a significant result is
# conservative evidence of some loss of fit.

use_sparse <- isTRUE(bic_spa < bic_unc)

model_final       <- if (use_sparse) model_pruned else model_unconstrained
model_final_label <- if (use_sparse) "sparse (pruned + step-up)" else "unconstrained"

selection_tbl <- tibble::tibble(
  Quantity = c("Retained model",
               "npar unconstrained", "df unconstrained",
               "npar sparse", "df sparse",
               "BIC unconstrained", "BIC sparse", "BIC difference",
               "AIC unconstrained", "AIC sparse",
               "Chi-square difference test p"),
  Value = c(model_final_label,
            npar_unc, df_unc, npar_spa, df_spa,
            round(bic_unc, 2), round(bic_spa, 2), round(bic_unc - bic_spa, 2),
            round(aic_unc, 2), round(aic_spa, 2),
            format.pval(chisq_diff_p, digits = 3, eps = .001))
)

write.csv(selection_tbl,
          file.path(RESULTS_DIR, "panelgvar_model_selection.csv"),
          row.names = FALSE)
print(selection_tbl)

saveRDS(model_final, file.path(RESULTS_DIR, "panelgvar_model_final.rds"))


# -----------------------------------------------------------------------------
# 4. Network extraction
# -----------------------------------------------------------------------------
# Temporal edges are extracted as PDC (partial directed correlations:
# standardised lag-1 coefficients), NOT as beta. The two differ in orientation:
# getmatrix(m, "beta") is column -> row, getmatrix(m, "PDC") is row -> column,
# which is the convention qgraph expects. Extracting beta and plotting it as if
# it were PDC reverses every arrow and swaps in- and out-strength.

features_temporal        <- psychonetrics::getmatrix(model_final, "PDC")
features_contemporaneous <- psychonetrics::getmatrix(model_final, "omega_zeta_within")
features_between         <- psychonetrics::getmatrix(model_final, "omega_zeta_between")

labels_pgvar <- LABELS_SHORT[rownames(DESIGN_MAT)]

dimnames(features_temporal)        <- list(labels_pgvar, labels_pgvar)
dimnames(features_contemporaneous) <- list(labels_pgvar, labels_pgvar)
dimnames(features_between)         <- list(labels_pgvar, labels_pgvar)

# Autoregressive effects are removed before plotting and before centrality and
# community computations, and are reported separately.
features_temporal_cross <- features_temporal
diag(features_temporal_cross) <- 0

write.csv(features_temporal,        file.path(RESULTS_DIR, "panelgvar_temporal_PDC_full.csv"))
write.csv(features_temporal_cross,  file.path(RESULTS_DIR, "panelgvar_temporal_PDC_crosslagged.csv"))
write.csv(features_contemporaneous, file.path(RESULTS_DIR, "panelgvar_contemporaneous.csv"))
write.csv(features_between,         file.path(RESULTS_DIR, "panelgvar_between.csv"))

autoreg_tbl <- tibble::tibble(
  Node           = unname(LABELS_FULL[rownames(DESIGN_MAT)]),
  Autoregression = diag(features_temporal)
) %>% arrange(desc(Autoregression))

write.csv(autoreg_tbl,
          file.path(RESULTS_DIR, "panelgvar_autoregressions.csv"),
          row.names = FALSE)
print(autoreg_tbl)


# -----------------------------------------------------------------------------
# 5. Temporal network figure
# -----------------------------------------------------------------------------
# curve is set to a non-zero value so that reciprocal pairs (A -> B and B -> A)
# remain visually distinguishable. No display threshold is applied: every edge
# surviving pruning is shown.

png(file.path(FIG_DIR, "panelgvar_temporal.png"),
    width = 10, height = 9, units = "in", res = 300)
g_temp <- qgraph::qgraph(
  features_temporal_cross,
  layout      = "spring",
  directed    = TRUE,
  diag        = FALSE,
  labels      = labels_pgvar,
  nodeNames   = LABELS_FULL[rownames(DESIGN_MAT)],
  color       = NODE_COLORS[rownames(DESIGN_MAT)],
  groups      = NODE_GROUPS,
  legend      = TRUE,
  legend.mode = "names",
  legend.cex  = 0.55,
  GLratio     = 3,
  curveAll    = TRUE,
  curve       = 0.3,
  vsize       = 9,
  asize       = 4,
  edge.labels = FALSE,
  theme       = "colorblind",
  title       = "Temporal network (lag-1, cross-lagged PDC)"
)
dev.off()


# -----------------------------------------------------------------------------
# 6. Retained temporal edges and their directionality
# -----------------------------------------------------------------------------
temporal_edges <- make_long(features_temporal_cross, "PDC") %>%
  filter(PDC != 0) %>%
  mutate(
    Sign   = if_else(PDC > 0, "positive", "negative"),
    AbsPDC = abs(PDC),
    Type   = edge_type(From, To)
  ) %>%
  arrange(desc(AbsPDC))

write.csv(temporal_edges,
          file.path(RESULTS_DIR, "panelgvar_temporal_edges_retained.csv"),
          row.names = FALSE)
print(temporal_edges)

# Central descriptive result: how many temporal edges run from functioning to
# psychopathology versus the reverse. A non-trivial number of
# Functioning -> Psychopathology edges is what supports the claim that
# functional domains are not merely downstream consequences of symptoms.
edge_types <- temporal_edges %>%
  group_by(Type) %>%
  summarise(
    n_edges     = dplyr::n(),
    n_positive  = sum(Sign == "positive"),
    n_negative  = sum(Sign == "negative"),
    mean_absPDC = mean(AbsPDC),
    max_absPDC  = max(AbsPDC),
    .groups = "drop"
  ) %>%
  arrange(desc(n_edges))

write.csv(edge_types,
          file.path(RESULTS_DIR, "panelgvar_temporal_edge_types.csv"),
          row.names = FALSE)
print(edge_types)


# -----------------------------------------------------------------------------
# 7. Temporal centrality: in- and out-strength
# -----------------------------------------------------------------------------
# PDC is oriented row -> column, so rows send and columns receive.
in_strength  <- colSums(abs(features_temporal_cross))
out_strength <- rowSums(abs(features_temporal_cross))

strength_tbl <- tibble::tibble(
  Node        = unname(LABELS_FULL[rownames(DESIGN_MAT)]),
  Short       = unname(labels_pgvar),
  InStrength  = as.numeric(in_strength),
  OutStrength = as.numeric(out_strength)
) %>% arrange(desc(OutStrength))

write.csv(strength_tbl,
          file.path(RESULTS_DIR, "panelgvar_temporal_strength.csv"),
          row.names = FALSE)
print(strength_tbl)

# Nodes with high out-strength send the largest amount of standardised lag-1
# influence to the rest of the system: they precede subsequent change in other
# domains and are the candidate leverage points of the network. Nodes with high
# in-strength are predicted by the previous state of the system and are better
# read as downstream markers of network activation than as intervention targets.

radar_data <- rbind(
  max         = rep(max(c(in_strength, out_strength)), length(in_strength)),
  min         = rep(0, length(in_strength)),
  InStrength  = in_strength,
  OutStrength = out_strength
) %>% as.data.frame()
colnames(radar_data) <- labels_pgvar

radar_cols <- c(rgb(0.2, 0.5, 0.5, 0.4), rgb(0.8, 0.2, 0.5, 0.4))

png(file.path(FIG_DIR, "panelgvar_temporal_strength_radar.png"),
    width = 3000, height = 2000, res = 300)
fmsb::radarchart(
  radar_data, axistype = 1,
  pcol = radar_cols, pfcol = scales::alpha(radar_cols, 0.5),
  plwd = 2, plty = 1, cglcol = "black", cglty = 1, cglwd = 0.5,
  seg = 4, axislabcol = "black", vlcex = 0.9,
  caxislabels = round(seq(0, max(c(in_strength, out_strength)), length.out = 5), 2),
  title = "Temporal network strength (In vs Out)"
)
legend(x = 1, y = 1.2, legend = c("InStrength", "OutStrength"),
       bty = "n", pch = 20, col = radar_cols, cex = 1.1)
dev.off()


# -----------------------------------------------------------------------------
# 8. Contemporaneous and between-subjects networks
# -----------------------------------------------------------------------------
plot_undirected <- function(m, title, file) {
  png(file, width = 10, height = 9, units = "in", res = 300)
  g <- qgraph::qgraph(
    m,
    layout      = g_temp$layout,
    directed    = FALSE,
    diag        = FALSE,
    labels      = labels_pgvar,
    nodeNames   = LABELS_FULL[rownames(DESIGN_MAT)],
    color       = NODE_COLORS[rownames(DESIGN_MAT)],
    groups      = NODE_GROUPS,
    legend      = TRUE,
    legend.mode = "names",
    legend.cex  = 0.55,
    GLratio     = 3,
    curveAll    = FALSE,
    vsize       = 9,
    edge.labels = FALSE,
    theme       = "colorblind",
    title       = title
  )
  dev.off()
  invisible(g)
}

plot_undirected(features_contemporaneous,
                "Contemporaneous network (within-person)",
                file.path(FIG_DIR, "panelgvar_contemporaneous.png"))

plot_undirected(features_between,
                "Between-subjects network",
                file.path(FIG_DIR, "panelgvar_between.png"))

# Contemporaneous edges are within-person partial correlations at the same wave,
# after accounting for lag-1 effects: co-occurrence, not directional effect.
# Between-subjects edges are partial correlations between participants' mean
# levels: trait-like covariation, not within-person process.

undirected_edges <- function(m, label) {
  make_long(m, "PartialCor") %>%
    filter(From != To, PartialCor != 0) %>%
    rowwise() %>%
    mutate(Pair = paste(sort(c(From, To)), collapse = " -- ")) %>%
    ungroup() %>%
    distinct(Pair, .keep_all = TRUE) %>%
    arrange(desc(abs(PartialCor))) %>%
    transmute(Network = label, Pair, PartialCor)
}

write.csv(
  bind_rows(
    undirected_edges(features_contemporaneous, "Contemporaneous"),
    undirected_edges(features_between,         "Between-subjects")
  ),
  file.path(RESULTS_DIR, "panelgvar_undirected_edges.csv"),
  row.names = FALSE
)

density_pgvar <- bind_rows(
  network_density(features_temporal_cross, directed = TRUE) %>%
    mutate(Network = "Temporal (cross-lagged)"),
  network_density(features_contemporaneous) %>% mutate(Network = "Contemporaneous"),
  network_density(features_between)         %>% mutate(Network = "Between-subjects")
) %>% select(Network, everything())

write.csv(density_pgvar,
          file.path(RESULTS_DIR, "panelgvar_density.csv"),
          row.names = FALSE)
print(density_pgvar)


message("03_panel_GVAR.R complete. Retained model: ", model_final_label)
