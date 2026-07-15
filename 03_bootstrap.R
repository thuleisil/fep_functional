# =============================================================================
# 03_bootstrap.R
# Bootstrap stability of the temporal (lag-1) network. Re-estimates the dlvm1
# model on 1000 subsamples (75%, without replacement), summarises each edge
# (mean, SE, 95% CI), and compares bootstrap vs original pruned estimates.
#
# Requires data/residuals.csv (produced by 01_preprocessing_covariates.R).
# =============================================================================

library(tidyverse)
library(psychonetrics)
library(qgraph)
library(ggplot2)
library(tibble)

clean_labels <- function(x) {
  x <- gsub("^PANSS_", "", x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# -----------------------------------------------------------------------------
# 1. Load data and build the design matrix
# -----------------------------------------------------------------------------
data_full <- read_csv("data/residuals.csv") %>%
  dplyr::select(contains(c("PANSS", "HONOS"))) %>%
  as_tibble()

all_names <- colnames(data_full)
nodes     <- all_names[!grepl("^T[0-9]", all_names)]
waves     <- c("", "T1", "T2")

design_mat <- sapply(
  waves,
  function(w) if (w == "") nodes else paste0(w, nodes)
)
rownames(design_mat) <- nodes

labelsD <- clean_labels(nodes)   # cleaned labels for plots / matrices

# -----------------------------------------------------------------------------
# 2. Original pruned model (for the bootstrap-vs-original comparison)
#    Same specification and pruning as 02_panel_gvar_network.R.
# -----------------------------------------------------------------------------
set.seed(1234)
model_pruned <- dlvm1(
  data = data_full, vars = design_mat, estimator = "FIML",
  within_latent = "ggm", between_latent = "chol"
) %>%
  runmodel() %>%
  prune(alpha = 0.5, adjust = "fdr")

# -----------------------------------------------------------------------------
# 3. Bootstrap loop (1000 iterations, 75% subsample without replacement)
# -----------------------------------------------------------------------------
features_boot_temporal <- list()
features_boot_contempo <- list()
features_boot_between  <- list()
fit_boot               <- list()

set.seed(1234)
reps       <- 1000
start_time <- Sys.time()

for (i in 1:reps) {
  message("Running bootstrap iteration ", i, " of ", reps)

  rowstouse <- sort(sample(1:nrow(data_full),
                           size = floor(nrow(data_full) * 0.75),
                           replace = FALSE))
  boot_data <- data_full[rowstouse, ]

  model_boot <- dlvm1(
    data = boot_data, vars = design_mat, estimator = "FIML",
    within_latent = "ggm", between_latent = "chol"
  )

  fit_i <- tryCatch(runmodel(model_boot),
                    error = function(e) { message("Model failed at iteration ", i); NULL })

  if (!is.null(fit_i)) {
    features_boot_temporal[[i]] <- getmatrix(fit_i, "beta")
    features_boot_contempo[[i]] <- getmatrix(fit_i, "omega_zeta_within")
    features_boot_between[[i]]  <- getmatrix(fit_i, "omega_zeta_between")
    fit_boot[[i]] <- tryCatch(fit(fit_i) %>% as.data.frame(), error = function(e) NULL)
  }
}

end_time <- Sys.time()
print(end_time - start_time)

# -----------------------------------------------------------------------------
# 4. Drop failed iterations
# -----------------------------------------------------------------------------
features_boot_temporal <- features_boot_temporal[!sapply(features_boot_temporal, is.null)]
features_boot_contempo <- features_boot_contempo[!sapply(features_boot_contempo, is.null)]
features_boot_between  <- features_boot_between[!sapply(features_boot_between, is.null)]
length(features_boot_temporal)

# -----------------------------------------------------------------------------
# 5. Save raw bootstrap results
# -----------------------------------------------------------------------------
if (!dir.exists("results")) dir.create("results")
saveRDS(features_boot_temporal, "results/features_boot_temporal_dlvm1.RDS")
saveRDS(features_boot_contempo, "results/features_boot_contempo_dlvm1.RDS")
saveRDS(features_boot_between,  "results/features_boot_between_dlvm1.RDS")
saveRDS(fit_boot,               "results/fit_boot_dlvm1.RDS")

# -----------------------------------------------------------------------------
# 6. Summarise bootstrap matrices (mean, SD, SE, 95% CI, significant mean)
# -----------------------------------------------------------------------------
summarise_boot <- function(x) {
  boot_array <- array(unlist(x), dim = c(dim(x[[1]]), length(x)))

  boot_mean <- apply(boot_array, c(1, 2), mean, na.rm = TRUE)
  boot_sd   <- apply(boot_array, c(1, 2), sd,   na.rm = TRUE)
  boot_se   <- boot_sd / sqrt(length(x))
  boot_lower <- boot_mean - 1.96 * boot_se
  boot_upper <- boot_mean + 1.96 * boot_se

  boot_sig <- boot_mean
  boot_sig[boot_lower < 0 & boot_upper > 0] <- 0   # CI crosses zero -> set to 0

  dimnames(boot_mean)  <- list(labelsD, labelsD)
  dimnames(boot_sd)    <- list(labelsD, labelsD)
  dimnames(boot_se)    <- list(labelsD, labelsD)
  dimnames(boot_lower) <- list(labelsD, labelsD)
  dimnames(boot_upper) <- list(labelsD, labelsD)
  dimnames(boot_sig)   <- list(labelsD, labelsD)

  list(mean = boot_mean, sd = boot_sd, se = boot_se,
       lowerCI = boot_lower, upperCI = boot_upper, sigmean = boot_sig)
}

av_temporal <- summarise_boot(features_boot_temporal)
av_contempo <- summarise_boot(features_boot_contempo)
av_between  <- summarise_boot(features_boot_between)

write.csv(av_temporal$mean,    "results/bootstrap_temporal_mean.csv")
write.csv(av_temporal$se,      "results/bootstrap_temporal_SE.csv")
write.csv(av_temporal$lowerCI, "results/bootstrap_temporal_lowerCI.csv")
write.csv(av_temporal$upperCI, "results/bootstrap_temporal_upperCI.csv")
write.csv(av_temporal$sigmean, "results/bootstrap_temporal_sigmean.csv")

# -----------------------------------------------------------------------------
# 7. Bootstrapped temporal network plot
# -----------------------------------------------------------------------------
qgraph(
  av_temporal$sigmean,
  directed = TRUE, diag = TRUE, labels = labelsD,
  theme = "colorblind", threshold = 0.1,
  title = "Bootstrapped temporal network"
)

# -----------------------------------------------------------------------------
# 8. Bootstrap vs original estimates (per-edge CI plot)
# -----------------------------------------------------------------------------
real_temporal <- getmatrix(model_pruned, "beta")

diag(real_temporal)       <- 0
diag(av_temporal$mean)    <- 0
diag(av_temporal$lowerCI) <- 0
diag(av_temporal$upperCI) <- 0

dimnames(real_temporal)      <- list(labelsD, labelsD)
dimnames(av_temporal$mean)   <- list(labelsD, labelsD)
dimnames(av_temporal$lowerCI) <- list(labelsD, labelsD)
dimnames(av_temporal$upperCI) <- list(labelsD, labelsD)

make_long <- function(mat, value_name) {
  mat %>%
    as.data.frame() %>%
    rownames_to_column("Node2") %>%
    pivot_longer(cols = -Node2, names_to = "Node1", values_to = value_name)
}

boot_df  <- make_long(av_temporal$mean,    "BootstrapMean")
lower_df <- make_long(av_temporal$lowerCI, "LowerCI")
upper_df <- make_long(av_temporal$upperCI, "UpperCI")
real_df  <- make_long(real_temporal,       "RealEstimate")

intervals <- boot_df %>%
  left_join(lower_df, by = c("Node1", "Node2")) %>%
  left_join(upper_df, by = c("Node1", "Node2")) %>%
  left_join(real_df,  by = c("Node1", "Node2")) %>%
  mutate(
    Edge        = paste(Node2, "→", Node1),
    SelfLoop    = Node1 == Node2,
    Significant = !(LowerCI < 0 & UpperCI > 0),
    AbsEstimate = abs(BootstrapMean)
  ) %>%
  filter(!SelfLoop) %>%
  arrange(desc(AbsEstimate))

intervals$Edge <- factor(intervals$Edge, levels = rev(intervals$Edge))
sum(is.na(intervals$RealEstimate))   # check: should be 0

p_boot_original <- ggplot(intervals, aes(y = Edge)) +
  geom_errorbarh(aes(xmin = LowerCI, xmax = UpperCI),
                 height = 0.18, color = "grey60", linewidth = 3, alpha = 0.8) +
  geom_point(aes(x = RealEstimate, color = "Original model"), size = 3, alpha = 0.85) +
  geom_point(aes(x = BootstrapMean, color = "Bootstrap"), size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8) +
  scale_color_manual(values = c("Original model" = "#D64933", "Bootstrap" = "black"),
                     breaks = c("Original model", "Bootstrap")) +
  labs(title = "Bootstrap temporal edge estimates", x = "Edge estimate", y = "", color = "") +
  theme_bw() +
  theme(
    axis.text.y     = element_text(size = 8),
    axis.text.x     = element_text(size = 10),
    plot.title      = element_text(size = 15, face = "bold"),
    legend.position = "right",
    legend.text     = element_text(size = 11)
  )

print(p_boot_original)

ggsave("results/bootstrap_vs_original_temporal.png",
       plot = p_boot_original, width = 7, height = 8, dpi = 300)
