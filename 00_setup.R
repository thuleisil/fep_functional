# =============================================================================
# 00_setup.R
# PrEP - TV / UniBo - FEP panel-GVAR network analysis
# Load raw data, build the PANSS five-factor scores (van der Gaag) and HoNOS
# items across the three waves (T0, T1, T2), and produce the analytic data frame
# `prep_df_drop` (complete cases, n = 368).
# =============================================================================

library(haven)
library(tidyverse)
library(summarytools)
library(networktools)
library(qgraph)
library(bootnet)
library(psychonetrics)

# -----------------------------------------------------------------------------
# 1. Load raw SPSS data
# -----------------------------------------------------------------------------
prep_df <- read_sav(file = "data/FEP-PrEP-Roma.sav")

# -----------------------------------------------------------------------------
# 2. Recode sociodemographic / clinical variables
# -----------------------------------------------------------------------------
prep_df_mod <- prep_df %>%
  mutate(gender = as.factor(GENDER)) %>%
  mutate(
    age                            = ETA,
    race                           = as.factor(ETNIA),
    religion                       = as.factor(Regione_Nascita),
    ses                            = as.factor(STATOCIV),
    housing_condition              = as.factor(DOMICIL),
    school_years                   = ANNISCOL,
    occupation                     = as.factor(Occupazione),
    referral                       = as.factor(INVIO),
    previous_psychiatrist_contact  = as.factor(PRE_CONTATTO),
    previous_suicide_attempt       = as.factor(PRE_TS),
    diagnosis_t0                   = as.factor(T0_DIAGNOSI),
    dui                            = DUI,
    substance_abuse                = as.factor(ABUSO_Sostanze),
    antipsychotics_t0              = as.factor(ANTIPSIC1)
  ) %>%
  mutate(across(c(ANTIDEP1, ANTIDEP2, STABILIZ1, STABILIZ2, BENZO), ~ as.factor(.x)))

# -----------------------------------------------------------------------------
# 3. PANSS five-factor model (van der Gaag) at T0, T1, T2 via rowSums
#    Keep only PANSS + HoNOS9-12 + covariates, then complete cases (n = 368)
# -----------------------------------------------------------------------------
prep_df_mod <- prep_df_mod %>%
  select(
    contains("PANSS"),
    HONOS9:HONOS12,
    T1HONOS9:T1HONOS12,
    T2HONOS9:T2HONOS12,
    ETA, GENDER, ETNIA,
    ANTIPSIC1, ANTIDEP1, STABILIZ1, BENZO,
    T0_DIAGNOSI, T1_CBT, T1_PSICED, T1_CMREC
  ) %>%
  drop_na() %>%
  mutate(
    # ---- T0 ----
    PANSS_positive        = rowSums(select(., PANSSP1, PANSSP3, PANSSP5, PANSSP6, PANSSG9),        na.rm = TRUE),
    PANSS_negative        = rowSums(select(., PANSSN1, PANSSN2, PANSSN3, PANSSN4, PANSSN6, PANSSG7), na.rm = FALSE),
    PANSS_disorganization = rowSums(select(., PANSSP2, PANSSN5, PANSSN7, PANSSG11),                na.rm = FALSE),
    PANSS_excitement      = rowSums(select(., PANSSP4, PANSSP7, PANSSG8, PANSSG14),                na.rm = FALSE),
    PANSS_emotional       = rowSums(select(., PANSSG1, PANSSG2, PANSSG3, PANSSG4, PANSSG6),        na.rm = FALSE),

    # ---- T1 ----
    T1PANSS_positive        = rowSums(select(., T1PANSSP1, T1PANSSP3, T1PANSSP5, T1PANSSP6, T1PANSSG9),        na.rm = FALSE),
    T1PANSS_negative        = rowSums(select(., T1PANSSN1, T1PANSSN2, T1PANSSN3, T1PANSSN4, T1PANSSN6, T1PANSSG7), na.rm = FALSE),
    T1PANSS_disorganization = rowSums(select(., T1PANSSP2, T1PANSSN5, T1PANSSN7, T1PANSSG11),                  na.rm = FALSE),
    T1PANSS_excitement      = rowSums(select(., T1PANSSP4, T1PANSSP7, T1PANSSG8, T1PANSSG14),                  na.rm = FALSE),
    T1PANSS_emotional       = rowSums(select(., T1PANSSG1, T1PANSSG2, T1PANSSG3, T1PANSSG4, T1PANSSG6),        na.rm = FALSE),

    # ---- T2 ----
    T2PANSS_positive        = rowSums(select(., T2PANSSP1, T2PANSSP3, T2PANSSP5, T2PANSSP6, T2PANSSG9),        na.rm = FALSE),
    T2PANSS_negative        = rowSums(select(., T2PANSSN1, T2PANSSN2, T2PANSSN3, T2PANSSN4, T2PANSSN6, T2PANSSG7), na.rm = FALSE),
    T2PANSS_disorganization = rowSums(select(., T2PANSSP2, T2PANSSN5, T2PANSSN7, T2PANSSG11),                  na.rm = FALSE),
    T2PANSS_excitement      = rowSums(select(., T2PANSSP4, T2PANSSP7, T2PANSSG8, T2PANSSG14),                  na.rm = FALSE),
    T2PANSS_emotional       = rowSums(select(., T2PANSSG1, T2PANSSG2, T2PANSSG3, T2PANSSG4, T2PANSSG6),        na.rm = FALSE)
  )

# -----------------------------------------------------------------------------
# 4. Final analytic data frame: covariates + 5 PANSS factors x 3 waves + HoNOS
# -----------------------------------------------------------------------------
prep_df_factors <- prep_df_mod %>%
  select(
    # Sociodemographic / clinical covariates
    ETA, GENDER, ETNIA,
    ANTIPSIC1, ANTIDEP1, STABILIZ1, BENZO,
    T0_DIAGNOSI, T1_CBT, T1_PSICED, T1_CMREC,

    # PANSS factors - T0
    PANSS_positive, PANSS_negative, PANSS_disorganization, PANSS_excitement, PANSS_emotional,
    # PANSS factors - T1
    T1PANSS_positive, T1PANSS_negative, T1PANSS_disorganization, T1PANSS_excitement, T1PANSS_emotional,
    # PANSS factors - T2
    T2PANSS_positive, T2PANSS_negative, T2PANSS_disorganization, T2PANSS_excitement, T2PANSS_emotional,

    # HoNOS items - T0, T1, T2
    c(HONOS9:HONOS12),
    c(T1HONOS9:T1HONOS12),
    c(T2HONOS9:T2HONOS12)
  )

# Analytic data frame used by all downstream scripts
prep_df_drop <- prep_df_factors

# -----------------------------------------------------------------------------
# Helper: clean node labels for plotting (drop the "PANSS_" prefix, capitalise).
# Shared by the network scripts so that nodes read "Excitement", "Negative",
# ... instead of "PANSS_excitement". HoNOS9-12 are left unchanged.
# -----------------------------------------------------------------------------
clean_labels <- function(x) {
  x <- gsub("^PANSS_", "", x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
