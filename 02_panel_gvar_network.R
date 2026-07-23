# =============================================================================
# 02_panel_gvar_network.R
# Panel-GVAR (psychonetrics::dlvm1) on the residualised nodes. Estimates the
# saturated model, prunes it, extracts the temporal / contemporaneous /
# between-subjects networks, plots them, and computes temporal centrality.
#
# Requires data/residuals.csv (produced by 01_preprocessing_covariates.R).
# =============================================================================

library(tidyverse)
library(psychonetrics)
library(qgraph)
library(fmsb)
library(scales)

# Node-label helper (drops "PANSS_" prefix, capitalises). Also defined in
# 00_setup.R; redefined here so this script is runnable standalone.
clean_labels <- function(x) {
  x <- gsub("^PANSS_", "", x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# -----------------------------------------------------------------------------
# 1. Load residuals and build the design matrix (node x wave)
# -----------------------------------------------------------------------------
data <- read_csv("data/residuals.csv") %>%
  dplyr::select(contains(c("PANSS", "HONOS"))) %>%
  as_tibble()

all_names <- colnames(data)
symptoms  <- all_names[!grepl("^T[0-9]", all_names)]   # 9 baseline nodes
waves     <- c("", "T1", "T2")

design_mat <- sapply(
  waves,
  function(w) if (w == "") symptoms else paste0(w, symptoms)
)
rownames(design_mat) <- symptoms
dim(design_mat)   # 9 x 3

# -----------------------------------------------------------------------------
# 2. Panel-GVAR estimation (dlvm1): saturated then pruned (sparse) model
# -----------------------------------------------------------------------------
model <- dlvm1(
  data          = data,
  vars          = design_mat,
  estimator     = "FIML",
  within_latent = "ggm",
  between_latent = "chol"
)

set.seed(1234)

model_sat    <- model %>% runmodel()
model_pruned <- model_sat %>% prune(alpha = 0.05, adjust = "none")
# NOTE: pruning alpha here (0.5, FDR) differs from 04_community_spinglass.R
# (0.05, recursive). Keep consistent with the manuscript / reconcile if needed.

# -----------------------------------------------------------------------------
# 3. Fit indices (descriptive) and model comparison
# -----------------------------------------------------------------------------
fit_sat    <- fit(model_sat)
fit_pruned <- fit(model_pruned)

psychonetrics::compare(
  saturate = model_sat,
  sparse   = model_pruned
)

# -----------------------------------------------------------------------------
# 4. Extract the three networks
# -----------------------------------------------------------------------------
features_temporal        <- getmatrix(model_pruned, "beta")               # lag-1 (unstandardised)
features_contemporaneous <- getmatrix(model_pruned, "omega_zeta_within")  # within-person
features_between         <- getmatrix(model_pruned, "omega_zeta_between") # between-subjects

# -----------------------------------------------------------------------------
# 5. Clean node labels for plotting (Positive, Negative, ... ; HoNOS9-12)
# -----------------------------------------------------------------------------
labels_clean <- clean_labels(rownames(design_mat))

# -----------------------------------------------------------------------------
# 6. Network plots (shared layout from the temporal network)
# -----------------------------------------------------------------------------
g_temp <- qgraph(
  features_temporal,
  layout = "spring", directed = TRUE, diag = TRUE,
  labels = labels_clean, curveAll = TRUE, curve = 0,
  vsize = 6, asize = 4, edge.labels = FALSE, theme = "colorblind",
  title = "Temporal network (lag-1)"
)

g_cont <- qgraph(
  features_contemporaneous,
  layout = g_temp$layout, directed = FALSE, diag = FALSE,
  labels = labels_clean, curveAll = FALSE, curve = 1,
  vsize = 6, edge.labels = FALSE, theme = "colorblind",
  title = "Contemporaneous network (within-person)"
)

g_between <- qgraph(
  features_between,
  layout = g_temp$layout, directed = FALSE, diag = FALSE,
  labels = labels_clean, curveAll = FALSE, curve = 1,
  vsize = 6, edge.labels = FALSE, theme = "colorblind",
  title = "Between-subjects network"
)

# -----------------------------------------------------------------------------
# 7. Temporal centrality: in-strength (received) and out-strength (sent)
# -----------------------------------------------------------------------------
in_strength  <- colSums(abs(features_temporal))
out_strength <- rowSums(abs(features_temporal))
names(in_strength)  <- labels_clean
names(out_strength) <- labels_clean

in_strength
out_strength

# -----------------------------------------------------------------------------
# 8. Radar plot of temporal strength (In vs Out)
# -----------------------------------------------------------------------------
global_strength <- c(in_strength, out_strength)

max_min_data <- matrix(NA, nrow = 4, ncol = length(in_strength))
max_min_data[1, ] <- max(global_strength)
max_min_data[2, ] <- 0
max_min_data[3, ] <- in_strength
max_min_data[4, ] <- out_strength
rownames(max_min_data) <- c("max", "min", "InStrength", "OutStrength")
colnames(max_min_data) <- labels_clean
max_min_data <- as.data.frame(max_min_data)

colors_in   <- c(rgb(0.2, 0.5, 0.5, 0.4), rgb(0.8, 0.2, 0.5, 0.4))
caxislabels <- seq(0, round(max(global_strength), 2), length.out = 5)

draw_radar <- function() {
  radarchart(
    max_min_data, axistype = 1,
    pcol = colors_in, pfcol = alpha(colors_in, 0.5), plwd = 2,
    cglcol = "black", cglty = 1, cglwd = 0.5, seg = 4,
    axislabcol = "black", vlcex = 0.7,
    caxislabels = round(caxislabels, 2),
    title = "Temporal Network Strength (In vs Out)"
  )
  legend(x = 1, y = 1.2, legend = c("InStrength", "OutStrength"),
         bty = "n", pch = 20, col = colors_in, cex = 1.2)
}

draw_radar()

# -----------------------------------------------------------------------------
# 9. Save radar plot
# -----------------------------------------------------------------------------
if (!dir.exists("results")) dir.create("results")
tiff("results/radar_temporal_strength.tiff", width = 3000, height = 2000, res = 300)
draw_radar()
dev.off()

# Save the pruned model so downstream scripts can reuse it if desired
saveRDS(model_pruned, "results/model_pruned_dlvm1.RDS")
