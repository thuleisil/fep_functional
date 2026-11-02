# =============================================================================
# 00_setup.R
#
# The Interplay between Psychopathology and Functioning in First-Episode
# Psychosis: a two-year longitudinal network study.
#
# Parma Early Psychosis (Pr-EP) programme - Tor Vergata / UniBo collaboration.
#
# Paths, packages, node structure, labels, colours, analysis constants and
# shared helper functions. Every other script starts by sourcing this file
# from the repository root:
#
#   source("00_setup.R")
#
# Nothing here reads data or produces output.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------
# All paths are relative to the repository root. Open the .Rproj file, or set
# the working directory to the repository root, before sourcing. No setwd() is
# used anywhere in the pipeline.

if (!dir.exists("data")) {
  stop(
    "Could not find the 'data/' directory.\n",
    "  Current working directory: ", getwd(), "\n",
    "  Set the working directory to the repository root and create data/."
  )
}

DATA_DIR    <- "data"
RESULTS_DIR <- "results"
FIG_DIR     <- file.path(RESULTS_DIR, "figures")
COMM_DIR    <- file.path(RESULTS_DIR, "community")

for (d in c(RESULTS_DIR, FIG_DIR, COMM_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Raw SPSS file. Not distributed with this repository (see README).
SAV_FILE <- file.path(DATA_DIR, "FEP-PrEP-Roma.sav")

# Intermediate objects written by 01_preprocessing.R
RESIDUALS_FILE <- file.path(DATA_DIR, "residuals.csv")
SAMPLE_FILE    <- file.path(DATA_DIR, "analytic_sample.rds")


# -----------------------------------------------------------------------------
# 2. Packages
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(haven)                   # SPSS import
  library(tidyverse)               # data manipulation and plotting
  library(labelled)                # labelled vectors from SPSS
  library(rstatix)                 # Friedman test

  library(qgraph)                  # network visualisation and centrality
  library(bootnet)                 # GGM estimation, stability, simulation
  library(networktools)            # bridge centrality
  library(NetworkComparisonTest)   # permutation comparison of networks
  library(psychonetrics)           # panel-GVAR (dlvm1)
  library(igraph)                  # graph objects and community detection
  library(reshape2)                # melt() for the heatmap
  library(fmsb)                    # radar chart
  library(scales)                  # colour transparency
})


# -----------------------------------------------------------------------------
# 3. Analysis constants
# -----------------------------------------------------------------------------
# Every stochastic step sets this seed immediately before running.
SEED <- 1234

# Cross-sectional GGMs
N_BOOT_NET <- 1000                 # case-dropping and non-parametric bootstrap
N_CORES    <- 4
SIM_CASES  <- seq(200, 1000, 200)  # netSimulator sample sizes
SIM_REPS   <- 100                  # replications per sample size

# Network Comparison Tests
N_NCT <- 1000                      # permutations per pairwise comparison

# Panel-GVAR model selection
PRUNE_ALPHA  <- 0.05               # prune(alpha = , recursive = TRUE)
STEPUP_ALPHA <- 0.05               # stepup(alpha = )

# Spinglass community detection
N_SPIN_ITER <- 1000                # random seeds
N_SPINS     <- 6                   # maximum number of communities per run


# -----------------------------------------------------------------------------
# 4. Covariates
# -----------------------------------------------------------------------------
# Each of the 27 node variables is residualised on these five time-invariant
# baseline covariates before any network is estimated.
COVARIATES <- c(
  "ETA",            # age
  "GENDER",         # sex
  "ETNIA",          # ethnicity
  "DUP",            # duration of untreated psychosis
  "T0_DIAGNOSI"     # baseline DSM-5 diagnosis
)


# -----------------------------------------------------------------------------
# 5. Nodes and design matrix
# -----------------------------------------------------------------------------
# Nine nodes: five PANSS factor scores plus four HoNOS items, measured at
# T0 (baseline), T1 (12 months) and T2 (24 months). 9 x 3 = 27 variables.

NODES <- c(
  "PANSS_positive",
  "PANSS_negative",
  "PANSS_disorganization",
  "PANSS_excitement",
  "PANSS_emotional",
  "HONOS9",
  "HONOS10",
  "HONOS11",
  "HONOS12"
)

WAVES <- c("", "T1", "T2")   # "" = baseline

DESIGN_MAT <- sapply(
  WAVES,
  function(w) if (w == "") NODES else paste0(w, NODES)
)
rownames(DESIGN_MAT) <- NODES

NODE_VARS <- unique(as.vector(DESIGN_MAT))

# PANSS five-factor composition (unweighted sums of the items listed)
PANSS_FACTORS <- list(
  PANSS_positive        = c("PANSSP1", "PANSSP3", "PANSSP5", "PANSSP6", "PANSSG9"),
  PANSS_negative        = c("PANSSN1", "PANSSN2", "PANSSN3", "PANSSN4", "PANSSN6", "PANSSG7"),
  PANSS_disorganization = c("PANSSP2", "PANSSN5", "PANSSN7", "PANSSG11"),
  PANSS_excitement      = c("PANSSP4", "PANSSP7", "PANSSG8", "PANSSG14"),
  PANSS_emotional       = c("PANSSG1", "PANSSG2", "PANSSG3", "PANSSG4", "PANSSG6")
)

HONOS_ITEMS <- c("HONOS9", "HONOS10", "HONOS11", "HONOS12")


# -----------------------------------------------------------------------------
# 6. Labels, colours, communities
# -----------------------------------------------------------------------------
LABELS_SHORT <- c(
  PANSS_positive        = "Pos",
  PANSS_negative        = "Neg",
  PANSS_disorganization = "Dis",
  PANSS_excitement      = "Exc",
  PANSS_emotional       = "Emo",
  HONOS9                = "H9",
  HONOS10               = "H10",
  HONOS11               = "H11",
  HONOS12               = "H12"
)

LABELS_FULL <- c(
  PANSS_positive        = "Pos - PANSS positive",
  PANSS_negative        = "Neg - PANSS negative",
  PANSS_disorganization = "Dis - PANSS disorganisation",
  PANSS_excitement      = "Exc - PANSS excitement",
  PANSS_emotional       = "Emo - PANSS emotional distress",
  HONOS9                = "H9 - Relationships / social functioning",
  HONOS10               = "H10 - Activities of daily living",
  HONOS11               = "H11 - Living conditions",
  HONOS12               = "H12 - Occupational & recreational activities"
)

# Node palette, chosen for print:
#   * the two groups differ in hue AND lightness, so they remain
#     distinguishable in greyscale and for deuteranopia / protanopia;
#   * neither hue is blue or red, so node fills never compete with the blue
#     (positive) / red (negative) edge colours;
#   * each node has a darker border of the same hue;
#   * label colour is set per group for adequate contrast on both fills.
PAL_PSY <- c(fill = "#E9C46A", border = "#A97E16", label = "#000000")  # amber
PAL_FUN <- c(fill = "#2A9D8F", border = "#1A6B61", label = "#FFFFFF")  # teal

NODE_COLORS       <- unname(c(rep(PAL_PSY["fill"],   5), rep(PAL_FUN["fill"],   4)))
NODE_BORDERS      <- unname(c(rep(PAL_PSY["border"], 5), rep(PAL_FUN["border"], 4)))
NODE_LABEL_COLORS <- unname(c(rep(PAL_PSY["label"],  5), rep(PAL_FUN["label"],  4)))

names(NODE_COLORS) <- names(NODE_BORDERS) <- names(NODE_LABEL_COLORS) <- NODES

# Accent colour for ggplot figures
ACCENT_COLOR <- unname(PAL_FUN["fill"])

NODE_GROUPS <- list(
  "PANSS - psychopathology" = 1:5,
  "HoNOS - functioning"     = 6:9
)

# Two-community vector used by networktools::bridge() and by NCT()
BRIDGE_COMMUNITIES <- c(rep("Psychopathology", 5), rep("Functioning", 4))


# -----------------------------------------------------------------------------
# 7. Helper functions
# -----------------------------------------------------------------------------

#' Mean with a normal-theory 95% confidence interval (no Hmisc dependency).
mean_ci <- function(x, conf = 0.95) {
  x <- x[is.finite(x)]
  n <- length(x)
  m <- mean(x)
  if (n < 2) return(data.frame(y = m, ymin = m, ymax = m))
  se    <- sd(x) / sqrt(n)
  tcrit <- qt(1 - (1 - conf) / 2, df = n - 1)
  data.frame(y = m, ymin = m - tcrit * se, ymax = m + tcrit * se)
}

#' z-standardisation.
z_score <- function(x) as.numeric(scale(x))

#' Tidy table of the standard centrality indices for one wave.
centrality_table <- function(net, wave) {
  qgraph::centrality_auto(net$graph)$node.centrality %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Var") %>%
    dplyr::mutate(
      Wave  = wave,
      Node  = unname(LABELS_SHORT[Var]),
      Label = unname(LABELS_FULL[Var])
    ) %>%
    dplyr::select(
      Wave, Node, Label,
      dplyr::any_of(c("Strength", "Closeness", "Betweenness", "ExpectedInfluence"))
    )
}

#' Tidy table of bridge centrality for one wave.
bridge_table <- function(br, wave) {
  tibble::tibble(
    Wave              = wave,
    Node              = unname(LABELS_SHORT[NODES]),
    Label             = unname(LABELS_FULL[NODES]),
    Community         = BRIDGE_COMMUNITIES,
    BridgeStrength    = as.numeric(br$`Bridge Strength`),
    BridgeBetweenness = as.numeric(br$`Bridge Betweenness`),
    BridgeCloseness   = as.numeric(br$`Bridge Closeness`),
    BridgeEI1         = as.numeric(br$`Bridge Expected Influence (1-step)`),
    BridgeEI2         = as.numeric(br$`Bridge Expected Influence (2-step)`)
  )
}

#' One styling function for every centrality figure, so that the strength and
#' bridge-strength panels are identical in theme, labels and scale.
#'
#' dat: data frame with columns Node | Index | Value (Value already
#' standardised). With colour_by = NULL a single black series is drawn; with
#' colour_by = "Wave" one series per wave.
centrality_panel <- function(dat, order_index = NULL,
                             colour_by = NULL, xlab = "z-score") {

  if (is.null(order_index)) order_index <- as.character(dat$Index[1])

  ord <- dat %>%
    dplyr::filter(Index == order_index) %>%
    dplyr::group_by(Node) %>%
    dplyr::summarise(v = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(v) %>%
    dplyr::pull(Node)

  dat <- dat %>%
    dplyr::mutate(
      Node  = factor(as.character(Node), levels = as.character(ord)),
      Index = factor(Index, levels = unique(Index))
    )

  p <- if (is.null(colour_by)) {
    ggplot(dat, aes(x = Value, y = Node, group = 1)) +
      geom_line(linewidth = 0.5, colour = "black") +
      geom_point(size = 2, colour = "black")
  } else {
    ggplot(dat, aes(x = Value, y = Node,
                    group  = .data[[colour_by]],
                    colour = .data[[colour_by]])) +
      geom_line(linewidth = 0.5) +
      geom_point(size = 2)
  }

  p +
    facet_wrap(~ Index, nrow = 1, scales = "free_x") +
    labs(x = xlab, y = NULL, colour = NULL) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90", colour = "grey20"),
      strip.text       = element_text(size = 11),
      axis.text.y      = element_text(size = 10),
      legend.position  = if (is.null(colour_by)) "none" else "bottom"
    )
}

#' Extract a single fit measure from a psychonetrics fit object.
fit_value <- function(f, measure) {
  f <- as.data.frame(f)
  v <- suppressWarnings(as.numeric(f$Value[f$Measure == measure]))
  if (length(v) == 0) NA_real_ else v[1]
}

#' Convert a square matrix into a long From / To / value data frame.
make_long <- function(mat, value_name) {
  mat %>%
    as.data.frame() %>%
    tibble::rownames_to_column("From") %>%
    tidyr::pivot_longer(cols = -From, names_to = "To", values_to = value_name)
}

#' Classify a directed edge by the domains it connects.
edge_type <- function(from, to) {
  psy <- LABELS_SHORT[1:5]
  fun <- LABELS_SHORT[6:9]
  dplyr::case_when(
    from %in% psy & to %in% fun ~ "Psychopathology -> Functioning",
    from %in% fun & to %in% psy ~ "Functioning -> Psychopathology",
    from %in% fun & to %in% fun ~ "Functioning -> Functioning",
    TRUE                        ~ "Psychopathology -> Psychopathology"
  )
}

#' Density of a weight matrix.
network_density <- function(m, directed = FALSE) {
  mm <- m
  diag(mm) <- 0
  k <- nrow(mm)
  n_possible <- if (directed) k * (k - 1) else choose(k, 2)
  n_edges    <- if (directed) sum(mm != 0) else sum(mm[upper.tri(mm)] != 0)
  tibble::tibble(
    `Non-zero edges` = n_edges,
    `Possible edges` = n_possible,
    Density          = n_edges / n_possible
  )
}

#' Guarded utilities: return NA instead of failing on empty input.
safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  v <- suppressWarnings(stats::cor(x[ok], y[ok], method = method))
  if (is.finite(v)) v else NA_real_
}

safe_max <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else max(x) }


# -----------------------------------------------------------------------------
# 8. Session report
# -----------------------------------------------------------------------------
message("Setup complete.")
message("  nodes        : ", length(NODES), " x ", length(WAVES),
        " waves = ", length(NODE_VARS), " variables")
message("  covariates   : ", paste(COVARIATES, collapse = ", "))
message("  results dir  : ", normalizePath(RESULTS_DIR, mustWork = FALSE))
