# =============================================================================
# 01_preprocessing_covariates.R
# Residualise each PANSS factor / HoNOS item (at each wave) on demographic and
# clinical covariates, z-standardise the residuals, and write data/residuals.csv
# (the input to all network scripts). Also produces the residual violin plots.
# =============================================================================

source("scripts/00_setup.R")

# -----------------------------------------------------------------------------
# 1. Node variables (PANSS factors + HoNOS) and design matrix (node x wave)
# -----------------------------------------------------------------------------
net_impute_vars <- prep_df_drop %>%
  dplyr::select(contains(c("PANSS", "HONOS")))

all_names <- colnames(net_impute_vars)

# Baseline nodes: those NOT starting with a wave prefix (T1/T2)
symptoms <- all_names[!grepl("^T[0-9]", all_names)]
length(symptoms)   # 9 nodes (5 PANSS factors + 4 HoNOS)
symptoms

waves <- c("", "T1", "T2")   # "" = baseline

design_mat <- sapply(
  waves,
  function(w) if (w == "") symptoms else paste0(w, symptoms)
)
rownames(design_mat) <- symptoms

# -----------------------------------------------------------------------------
# 2. Covariate regression: residualise every node x wave on the covariates
# -----------------------------------------------------------------------------
node_vars <- unique(as.vector(design_mat))   # all node columns across waves

covars <- c(
  "GENDER", "ETA", "ETNIA", "ANTIPSIC1", "ANTIDEP1", "STABILIZ1", "BENZO",
  "T0_DIAGNOSI", "T1_CBT", "T1_PSICED", "T1_CMREC"
)

net_impute <- prep_df_drop %>%
  dplyr::select(all_of(c(covars, node_vars))) %>%
  as.data.frame()

n_v <- nrow(design_mat)   # nodes
n_t <- ncol(design_mat)   # waves

set.seed(1234)

for (k in 1:n_v) {
  for (i in 1:n_t) {
    x_var <- design_mat[k, i]                              # e.g. "PANSS_positive", "T1PANSS_positive"
    f     <- as.formula(paste(x_var, "~", paste(covars, collapse = " + ")))
    fit   <- lm(f, data = net_impute)
    net_impute[[x_var]] <- residuals(fit)                  # replace with residuals
  }
}

# -----------------------------------------------------------------------------
# 3. Center & z-standardise the residualised node variables (mean 0, SD 1)
# -----------------------------------------------------------------------------
net_impute <- net_impute %>%
  mutate(across(all_of(node_vars),
                ~ as.numeric(scale(.x, center = TRUE, scale = TRUE))))

dim(net_impute)

# -----------------------------------------------------------------------------
# 4. Save residuals (input for the network / bootstrap / community scripts)
# -----------------------------------------------------------------------------
if (!dir.exists("data")) dir.create("data")
write.csv(net_impute, "data/residuals.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. Visual inspection: violin plots of standardised residuals by node x wave
# -----------------------------------------------------------------------------
node_long <- net_impute %>%
  dplyr::select(dplyr::all_of(c(symptoms, paste0("T1", symptoms), paste0("T2", symptoms)))) %>%
  pivot_longer(cols = everything(), names_to = "Var", values_to = "Value") %>%
  mutate(
    Time = dplyr::case_when(
      grepl("^T1", Var) ~ "T1",
      grepl("^T2", Var) ~ "T2",
      TRUE              ~ "T0"
    ),
    Item = gsub("^T[0-9]", "", Var),
    Time = factor(Time, levels = c("T0", "T1", "T2"))
  )

p_residuals <- ggplot(node_long, aes(x = Time, y = Value)) +
  geom_violin(trim = FALSE) +
  stat_summary(color = "red", size = 0.3, fun.data = mean_cl_normal) +
  facet_wrap(~ Item, ncol = 3) +
  labs(
    title = "Standardised residuals by node and timepoint",
    x     = "Time point",
    y     = "Standardised residuals"
  ) +
  theme_bw() +
  theme(
    text            = element_text(size = 9),
    strip.text      = element_text(size = 7),
    axis.title      = element_text(size = 10),
    axis.text.x     = element_text(size = 7),
    axis.text.y     = element_text(size = 7),
    plot.title      = element_text(size = 12, face = "bold"),
    panel.spacing.x = unit(0.15, "lines"),
    panel.spacing.y = unit(0.15, "lines")
  )

print(p_residuals)

if (!dir.exists("results")) dir.create("results")
ggsave("results/residuals_violin.png", p_residuals, width = 10, height = 9, dpi = 300)
