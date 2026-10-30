# =============================================================================
# 01_preprocessing.R
#
# Builds the analytic sample, computes the PANSS five-factor scores, residualises
# every node on the covariate set, and z-standardises the residuals.
#
#   source("R/00_setup.R")
#   source("R/01_preprocessing.R")
#
# Inputs
#   data/FEP-PrEP-Roma 2.sav        raw SPSS file (not distributed, see README)
#
# Outputs
#   data/residuals.csv              27 residualised, z-standardised node variables
#   data/analytic_sample.rds        analytic sample with raw scores and covariates
#   results/table1a_whole_cohort.csv
#   results/table1b_analytic_sample.csv
#   results/node_descriptives_by_wave.csv
#   results/figures/node_trajectories.png
#   results/figures/residual_distributions.png
# =============================================================================

stopifnot(exists("NODES"))   # 00_setup.R must be sourced first


# -----------------------------------------------------------------------------
# 1. Load raw SPSS data
# -----------------------------------------------------------------------------
if (!file.exists(SAV_FILE)) {
  stop("Raw data file not found:\n  ", SAV_FILE,
       "\nSee the Data availability section of the README.")
}

prep_df <- haven::read_sav(SAV_FILE)

message("Raw cohort: ", nrow(prep_df), " participants")


# -----------------------------------------------------------------------------
# 2. Recode categorical variables
# -----------------------------------------------------------------------------
prep_df_mod <- prep_df %>%
  mutate(across(c(GENDER, ETNIA, T0_DIAGNOSI), ~ as.factor(.x)))


# -----------------------------------------------------------------------------
# 3. Variable selection and complete cases
# -----------------------------------------------------------------------------
# HoNOS variables are named explicitly rather than selected with a `:` range,
# so that the selection does not depend on the column order of the .sav file.

honos_t0 <- HONOS_ITEMS
honos_t1 <- paste0("T1", HONOS_ITEMS)
honos_t2 <- paste0("T2", HONOS_ITEMS)

stopifnot(all(c(honos_t0, honos_t1, honos_t2) %in% names(prep_df_mod)))

#' Row sums for one PANSS factor at one wave.
#'
#' na.rm = FALSE throughout: a participant with a missing item contributes NA
#' to the factor rather than an artificially low sum. This is also what makes
#' the whole-cohort and analytic-sample tables comparable.
panss_factor <- function(data, items, prefix = "") {
  cols <- paste0(prefix, items)
  rowSums(data[, cols, drop = FALSE], na.rm = FALSE)
}

prep_df_mod <- prep_df_mod %>%
  select(
    contains("PANSS"),
    all_of(c(honos_t0, honos_t1, honos_t2)),
    all_of(COVARIATES)
  ) %>%
  drop_na()

for (w in WAVES) {
  for (fac in names(PANSS_FACTORS)) {
    prep_df_mod[[paste0(w, fac)]] <-
      panss_factor(prep_df_mod, PANSS_FACTORS[[fac]], prefix = w)
  }
}


# -----------------------------------------------------------------------------
# 4. Analytic data frame
# -----------------------------------------------------------------------------
prep_df_drop <- prep_df_mod %>%
  select(all_of(c(COVARIATES, NODE_VARS))) %>%
  mutate(ID = seq_len(n())) %>%
  select(ID, everything())

n_analytic  <- nrow(prep_df_drop)
analytic_id <- as.character(prep_df_drop$ID)

message("Analytic sample (complete cases at T0, T1, T2): n = ", n_analytic)

saveRDS(prep_df_drop, SAMPLE_FILE)


# -----------------------------------------------------------------------------
# 5. Descriptive tables
# -----------------------------------------------------------------------------
# Whole cohort, baseline characteristics. PANSS factors computed with
# na.rm = FALSE, as in the analytic sample.
whole_cohort <- prep_df
for (fac in names(PANSS_FACTORS)) {
  whole_cohort[[fac]] <- panss_factor(whole_cohort, PANSS_FACTORS[[fac]])
}

describe_vars <- c(COVARIATES, names(PANSS_FACTORS), HONOS_ITEMS)

#' Compact descriptive summary: median (Q1-Q3) for numeric variables,
#' n (%) per level for factors, plus the number of missing values.
describe_sample <- function(data, vars) {
  purrr::map_dfr(vars, function(v) {
    x <- data[[v]]
    if (is.null(x)) return(NULL)
    x <- haven::zap_labels(x)
    n_miss <- sum(is.na(x))
    if (is.numeric(x)) {
      q <- quantile(x, c(.25, .5, .75), na.rm = TRUE)
      tibble::tibble(
        Variable = v, Level = NA_character_,
        Summary  = sprintf("%.2f (%.2f-%.2f)", q[2], q[1], q[3]),
        Missing  = n_miss
      )
    } else {
      tab <- table(as.factor(x), useNA = "no")
      tibble::tibble(
        Variable = v,
        Level    = names(tab),
        Summary  = sprintf("%d (%.1f%%)", as.integer(tab),
                           100 * as.integer(tab) / sum(tab)),
        Missing  = n_miss
      )
    }
  })
}

table1a <- describe_sample(whole_cohort, describe_vars)
table1b <- describe_sample(prep_df_drop, describe_vars)

write.csv(table1a, file.path(RESULTS_DIR, "table1a_whole_cohort.csv"), row.names = FALSE)
write.csv(table1b, file.path(RESULTS_DIR, "table1b_analytic_sample.csv"), row.names = FALSE)


# -----------------------------------------------------------------------------
# 6. Node descriptives across waves + Friedman tests
# -----------------------------------------------------------------------------
nodes_long <- prep_df_drop %>%
  select(ID, all_of(NODE_VARS)) %>%
  mutate(across(all_of(NODE_VARS), ~ as.numeric(.x))) %>%
  pivot_longer(cols = all_of(NODE_VARS), names_to = "Var", values_to = "Score") %>%
  mutate(
    Time = case_when(
      grepl("^T1", Var) ~ "T1",
      grepl("^T2", Var) ~ "T2",
      TRUE              ~ "T0"
    ),
    Node = gsub("^T[12]", "", Var),
    Time = factor(Time, levels = c("T0", "T1", "T2")),
    Node = factor(Node, levels = NODES, labels = LABELS_FULL[NODES])
  )

node_summary <- nodes_long %>%
  group_by(Node, Time) %>%
  summarise(
    Median = median(Score, na.rm = TRUE),
    Q1     = quantile(Score, 0.25, na.rm = TRUE),
    Q3     = quantile(Score, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Median_IQR = sprintf("%.2f (%.2f-%.2f)", Median, Q1, Q3)) %>%
  select(Node, Time, Median_IQR) %>%
  pivot_wider(names_from = Time, values_from = Median_IQR)

friedman_tbl <- nodes_long %>%
  group_by(Node) %>%
  rstatix::friedman_test(Score ~ Time | ID) %>%
  select(Node, friedman_statistic = statistic, df, p)

node_descriptives <- left_join(node_summary, friedman_tbl, by = "Node")

write.csv(node_descriptives,
          file.path(RESULTS_DIR, "node_descriptives_by_wave.csv"),
          row.names = FALSE)


# -----------------------------------------------------------------------------
# 7. Trajectory figure
# -----------------------------------------------------------------------------
traj_plot <- nodes_long %>%
  group_by(Node, Time) %>%
  summarise(MeanScore = mean(Score, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Domain = if_else(grepl("^H", Node), "HoNOS - functioning", "PANSS - psychopathology"),
    Node   = factor(as.character(Node),
                    levels = LABELS_FULL[NODES],
                    labels = LABELS_SHORT[NODES])
  ) %>%
  ggplot(aes(x = Time, y = MeanScore, colour = Node, group = Node)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.6) +
  facet_wrap(~ Domain, scales = "free_y") +
  labs(title = "Raw node scores across timepoints",
       x = "Time point", y = "Mean score", colour = "") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 3))

ggsave(file.path(FIG_DIR, "node_trajectories.png"),
       traj_plot, width = 11, height = 6, dpi = 300)


# -----------------------------------------------------------------------------
# 8. Covariate regression and standardisation
# -----------------------------------------------------------------------------
# Each of the 27 node variables is regressed on the five covariates and the
# residuals are retained, then centred and z-standardised. The residuals index
# node severity independent of age, sex, ethnicity, DUP and baseline diagnosis.
#
# No detrending is applied: systematic change in mean level across waves is left
# in the data and is described by the Friedman tests above.

net_impute <- prep_df_drop %>%
  dplyr::select(all_of(c(COVARIATES, NODE_VARS))) %>%
  mutate(across(all_of(NODE_VARS), ~ as.numeric(.x))) %>%
  mutate(ETA = as.numeric(ETA), DUP = as.numeric(DUP)) %>%
  as.data.frame()

for (k in seq_len(nrow(DESIGN_MAT))) {
  for (i in seq_len(ncol(DESIGN_MAT))) {

    x_var <- DESIGN_MAT[k, i]

    f <- as.formula(paste(x_var, "~", paste(COVARIATES, collapse = " + ")))

    net_impute[[x_var]] <- residuals(lm(f, data = net_impute))
  }
}

net_impute <- net_impute %>%
  mutate(across(all_of(NODE_VARS),
                ~ as.numeric(scale(.x, center = TRUE, scale = TRUE))))

# Participant IDs are carried as row names so that the wave-specific data frames
# built in 02_GGMs.R can be verified to contain the same people in the same
# order - a precondition for the paired Network Comparison Tests.
rownames(net_impute) <- analytic_id

write.csv(net_impute[, NODE_VARS], RESIDUALS_FILE, row.names = TRUE)

message("Residuals written to ", RESIDUALS_FILE,
        "  (", nrow(net_impute), " x ", length(NODE_VARS), ")")


# -----------------------------------------------------------------------------
# 9. Residual distributions
# -----------------------------------------------------------------------------
resid_plot <- net_impute %>%
  dplyr::select(all_of(NODE_VARS)) %>%
  pivot_longer(cols = everything(), names_to = "Var", values_to = "Value") %>%
  mutate(
    Time = case_when(
      grepl("^T1", Var) ~ "T1",
      grepl("^T2", Var) ~ "T2",
      TRUE              ~ "T0"
    ),
    Node = gsub("^T[12]", "", Var),
    Time = factor(Time, levels = c("T0", "T1", "T2")),
    Node = factor(Node, levels = NODES, labels = LABELS_FULL[NODES])
  ) %>%
  ggplot(aes(x = Time, y = Value)) +
  geom_violin(trim = FALSE) +
  stat_summary(colour = "red", size = 0.3, fun.data = mean_ci) +
  facet_wrap(~ Node, ncol = 3) +
  labs(title = "Standardised residuals by node and timepoint",
       x = "Time point", y = "Standardised residuals") +
  theme_bw() +
  theme(text = element_text(size = 10), strip.text = element_text(size = 8))

ggsave(file.path(FIG_DIR, "residual_distributions.png"),
       resid_plot, width = 12, height = 12, dpi = 300)


message("01_preprocessing.R complete.")
