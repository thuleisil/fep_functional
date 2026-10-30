# =============================================================================
# 00_setup.R
#
# Functional domains and psychopathological dimensions in first-episode
# psychosis: cross-sectional GGMs and a two-year panel-GVAR network analysis.
#
# Parma Early Psychosis (Pr-EP) programme - Tor Vergata / UniBo collaboration.
#
# This script defines paths, packages, node structure, labels, analysis
# constants and shared helper functions. Every other script in the pipeline
# starts by sourcing it:
#
#   source("R/00_setup.R")
#
# Nothing here reads data or produces output.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Paths
# -----------------------------------------------------------------------------
# Paths are relative to the repository root. Open the .Rproj file, or set the
# working directory to the repository root before sourcing. No setwd() is used
# anywhere in the pipeline: inside a knitr chunk it is reverted at the end of
# the chunk, which silently breaks every relative path downstream.

if (!dir.exists("data")) {
  stop(
    "Could not find the 'data/' directory.\n",
    "  Current working directory: ", getwd(), "\n",
    "  Set the working directory to the repository root before sourcing."
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
# Edit the file name here if your local copy is named differently.
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
  library(gt)                      # tables
  library(rstatix)                 # Friedman test

  library(qgraph)                  # network visualisation and centrality
  library(bootnet)                 # GGM estimation, bootstrap, simulation
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
N_BOOT_NET <- 1000    # bootnet iterations (case-dropping and non-parametric)
N_CORES    <- 4       # cores for bootnet
SIM_CASES  <- seq(200, 1000, 200)
SIM_REPS   <- 100

# Network Comparison Tests
N_NCT <- 1000         # permutations per pairwise comparison

# Panel-GVAR
PRUNE_ALPHA  <- 0.05  # prune(alpha = , recursive = TRUE)
STEPUP_ALPHA <- 0.05  # stepup(alpha = )

# Spinglass community detection
N_SPIN_ITER <- 1000   # random seeds
N_SPINS     <- 6      # maximum number of communities per run


# -----------------------------------------------------------------------------
# 4. Covariates
# -----------------------------------------------------------------------------
# Each of the 27 node variables is residualised on these five covariates before
# any network is estimated. Pharmacological and psychosocial-treatment
# indicators are deliberately NOT included; see the README for the rationale
# and for the sensitivity analysis that reintroduces them.
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
# Nine nodes: the five PANSS factors of the van der Gaag five-factor model plus
# four HoNOS items indexing social and functional impairment. Measured at three
# waves: T0 (baseline), T1 (12 months), T2 (24 months). 9 x 3 = 27 variables.

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

# PANSS five-factor composition (van der Gaag et al., 2006)
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

NODE_COLORS <- c(
  rep("#C51B7D", 5),   # PANSS - psychopathology
  rep("#1B9E77", 4)    # HoNOS - functioning
)
names(NODE_COLORS) <- NODES

NODE_GROUPS <- list(
  "PANSS - psychopathology" = 1:5,
  "HoNOS - functioning"     = 6:9
)

# Two-community vector used by networktools::bridge() and by NCT()
BRIDGE_COMMUNITIES <- c(rep("Psychopathology", 5), rep("Functioning", 4))


# -----------------------------------------------------------------------------
# 7. Helper functions
# -----------------------------------------------------------------------------

#' Mean with a normal-theory 95% confidence interval.
#'
#' Replacement for ggplot2::mean_cl_normal, which requires Hmisc to be
#' installed. Used by the residual violin plots in 01_preprocessing.R.
mean_ci <- function(x, conf = 0.95) {
  x <- x[is.finite(x)]
  n <- length(x)
  m <- mean(x)
  if (n < 2) return(data.frame(y = m, ymin = m, ymax = m))
  se    <- sd(x) / sqrt(n)
  tcrit <- qt(1 - (1 - conf) / 2, df = n - 1)
  data.frame(y = m, ymin = m - tcrit * se, ymax = m + tcrit * se)
}

#' Weight matrix of a bootnet network with short node labels.
net_labelled <- function(net) {
  m <- net$graph
  dimnames(m) <- list(LABELS_SHORT[NODES], LABELS_SHORT[NODES])
  m
}

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
    Wave      = wave,
    Node      = unname(LABELS_SHORT[NODES]),
    Label     = unname(LABELS_FULL[NODES]),
    Community = BRIDGE_COMMUNITIES,
    BridgeStrength     = as.numeric(br$`Bridge Strength`),
    BridgeBetweenness  = as.numeric(br$`Bridge Betweenness`),
    BridgeCloseness    = as.numeric(br$`Bridge Closeness`),
    BridgeEI1          = as.numeric(br$`Bridge Expected Influence (1-step)`),
    BridgeEI2          = as.numeric(br$`Bridge Expected Influence (2-step)`)
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


# -----------------------------------------------------------------------------
# 8. Session report
# -----------------------------------------------------------------------------
message("Setup complete.")
message("  nodes        : ", length(NODES), " x ", length(WAVES),
        " waves = ", length(NODE_VARS), " variables")
message("  covariates   : ", paste(COVARIATES, collapse = ", "))
message("  results dir  : ", normalizePath(RESULTS_DIR, mustWork = FALSE))
