# The Interplay between Psychopathology and Functioning in First-Episode Psychosis

Analysis code for a two-year longitudinal network study of psychopathological
dimensions and functional domains in first-episode psychosis (FEP), using data
from the **Parma Early Psychosis (Pr-EP) programme**.

The pipeline estimates cross-sectional Gaussian Graphical Models at each wave,
compares them with paired Network Comparison Tests, fits a panel graphical
vector-autoregressive model across the three waves, and examines the community
structure of the resulting temporal network.

---

## Study design

- **Design**: prospective cohort of FEP patients enrolled in the Pr-EP
  programme, assessed at baseline (T0), 12 months (T1) and 24 months (T2).
- **Analytic sample**: participants with complete PANSS, HoNOS and covariate
  data at all three waves.
- **Nodes (9)**: five PANSS factor scores and four HoNOS items.

  | Node | Label | Content |
  |------|-------|---------|
  | Pos  | PANSS positive | P1, P3, P5, P6, G9 |
  | Neg  | PANSS negative | N1, N2, N3, N4, N6, G7 |
  | Dis  | PANSS disorganisation | P2, N5, N7, G11 |
  | Exc  | PANSS excitement | P4, P7, G8, G14 |
  | Emo  | PANSS emotional distress | G1, G2, G3, G4, G6 |
  | H9   | HoNOS 9 | Problems with relationships |
  | H10  | HoNOS 10 | Problems with activities of daily living |
  | H11  | HoNOS 11 | Problems with living conditions |
  | H12  | HoNOS 12 | Problems with occupational and recreational activities |

  9 nodes × 3 waves = 27 variables.

- **Covariates (5)**: age, sex, ethnicity, duration of untreated psychosis
  (DUP) and baseline DSM-5 diagnosis. Every node variable is residualised on
  these covariates and z-standardised before any network is estimated.

---

## Repository structure

```
.
├── 00_setup.R                          paths, packages, nodes, labels, colours, constants, helpers
├── 01_preprocessing.R                  analytic sample, PANSS factors, residualisation
├── 02_GGMs.R                           cross-sectional networks, stability, strength, bridge strength, NCT
├── 03_panel_GVAR.R                     panel-GVAR, model selection, temporal / contemporaneous / between networks
├── 04_community_analysis_spinglass.R   iterated Spinglass community detection on the temporal network
└── README.md
```

Create a `data/` folder in the repository root and place the raw SPSS file
there (see *Data availability*). Output is written to `results/`, which the
scripts create.

Run the scripts in order, from the repository root:

```r
source("00_setup.R")
source("01_preprocessing.R")
source("02_GGMs.R")
source("03_panel_GVAR.R")
source("04_community_analysis_spinglass.R")
```

`00_setup.R` must be sourced first in every session. Each later script reads the
intermediate files written by the previous one, so scripts `02`–`04` can be run
independently once `01_preprocessing.R` has produced `data/residuals.csv`.
No script calls `setwd()`.

---

## Where each manuscript element comes from

| Manuscript | Script | Output |
|---|---|---|
| Table 1 | `01` | `results/table1b_analytic_sample.csv` |
| Figure 1 (networks) | `02` | `results/figures/ggm_T0.png`, `ggm_T1.png`, `ggm_T2.png` |
| Figure 1 (S, M, Holm-adjusted p) | `02` | `results/nct_global_invariance.csv` |
| Figure 2a–c | `02` | `results/figures/strength_T0.png`, `_T1`, `_T2` |
| Figure 2d–f | `02` | `results/figures/bridge_strength_T0.png`, `_T1`, `_T2` |
| Figure 3a | `03` | `results/figures/panelgvar_temporal.png` |
| Figure 3b | `03` | `results/figures/panelgvar_temporal_strength_radar.png` |
| Figure 4 | `04` | `results/figures/community_co_occurrence_heatmap.png` |
| Model fit, AIC / BIC | `03` | `results/panelgvar_fit.csv`, `panelgvar_model_selection.csv` |
| Supplementary Table 2 | `02` | `results/centrality_T0_T1_T2.csv` |
| Supplementary Table 3 | `02` | `results/bridge_centrality_T0_T1_T2.csv` |
| Supplementary Table 4 | `03` | `results/panelgvar_temporal_edges_retained.csv` |
| Supplementary Table 5 | `03` | `results/panelgvar_temporal_strength.csv` |
| Supplementary Table 6 | `03` | `results/panelgvar_contemporaneous.csv` |
| Supplementary Table 7 | `03` | `results/panelgvar_between.csv` |
| Supplementary Table 8 | `04` | `results/community/community_co_occurrence_matrix.csv` |
| Supplementary Figure 1 | `01` | `results/figures/residual_distributions.png` |
| Supplementary Figures 2, 4, 6 | `02` | `results/figures/ggm_T0_case_stability.png`, `_T1_`, `_T2_` |
| Supplementary Figures 3, 5, 7 | `02` | `results/figures/ggm_T0_simulation.png`, `_T1_`, `_T2_` |
| Supplementary Figure 8 | `03` | `results/figures/panelgvar_contemporaneous.png` |
| Supplementary Figure 9 | `03` | `results/figures/panelgvar_between.png` |
| CS-coefficients | `02` | `results/ggm_cor_stability.csv` |

---

## Methods summary

### Preprocessing (`01_preprocessing.R`)

Complete cases on all PANSS items, HoNOS 9–12 and the five covariates at all
three waves. PANSS factor scores are unweighted item sums computed with
`na.rm = FALSE`. Each of the 27 node variables is regressed on the five
covariates; residuals are centred and z-standardised. Participant IDs are
carried through `data/residuals.csv` so that the wave-specific data sets used
for the paired NCT can be verified to contain the same people in the same order.

### Cross-sectional networks (`02_GGMs.R`)

Regularised partial-correlation networks at T0, T1 and T2 (graphical LASSO with
EBIC model selection, `EBICglasso`, γ = 0.5) via `bootnet::estimateNetwork`. The
T0 layout is reused at T1 and T2.

- **Stability**: case-dropping subset bootstrap (1,000 iterations) and
  correlation-stability (CS) coefficient. CS-coefficients are computed for all
  standard and bridge indices; strength and bridge strength are the indices
  plotted and reported.
- **Recovery**: `netSimulator`, n = 200–1,000 in steps of 200, 100 replications
  per sample size.
- **Strength** (`qgraph::centrality_auto`) and **bridge strength**
  (`networktools::bridge`, two communities: psychopathology and functioning).
  In the figures, values are z-standardised **within each wave**: each panel
  shows the relative ranking of nodes at that assessment, and absolute values
  are reported in the CSV files.

**Network Comparison Tests.** `NetworkComparisonTest::NCT` with
`paired = TRUE` and 1,000 permutations for each pairwise contrast
(T0–T1, T1–T2, T0–T2). The three waves are the same participants, so
relabelling is performed within each pair of observations. The package is
called with `p.adjust.methods = "none"`; corrections are applied once in the
script:

- **global strength (S)** and **network structure (M)**: Holm across the three
  temporal contrasts, separately for each of the two families of tests. These
  are the results reported in the manuscript. Package p-values are cross-checked
  against the stored permutation distributions;
- **edge-level** and **centrality** invariance: Benjamini–Hochberg within each
  contrast. These are written for completeness but not interpreted, because the
  structure test was not significant after Holm correction.

Edge p-values are read from `einv.pvals` by **column name** and matched by
**node pair**. In NetworkComparisonTest 2.2.3 the last column of `einv.pvals`
is the edge-invariance statistic E, not the p-value, so selecting it by position
would mislabel a statistic as a p-value.

### Panel-GVAR (`03_panel_GVAR.R`)

`psychonetrics::dlvm1` with `estimator = "FIML"`, `within_latent = "ggm"` and
`between_latent = "chol"` (Epskamp, 2020). Three networks are extracted:
temporal (directed, lag-1), contemporaneous (within-person, same wave) and
between-subjects (stable individual differences).

An **unconstrained** model is fitted first, then a sparse model with
`prune(alpha = 0.05, recursive = TRUE)` followed by `stepup(alpha = 0.05)`.
The sparse model is retained when BIC favours it; AIC, BIC, the χ² difference
test and all fit indices are saved. The unconstrained model is not saturated:
`dlvm1` imposes stationarity (the same temporal and contemporaneous parameters
for T0→T1 and T1→T2), which leaves degrees of freedom.

**Temporal edges are extracted as PDC** (partial directed correlations), not as
`beta`:

```r
getmatrix(model, "beta")   # oriented column -> row
getmatrix(model, "PDC")    # oriented row -> column, the qgraph convention
```

Reading `beta` as if it were PDC reverses every arrow and swaps in- and
out-strength. Here in-strength is `colSums(abs(PDC))` and out-strength is
`rowSums(abs(PDC))`, both computed with autoregressive effects removed.

### Community analysis (`04_community_analysis_spinglass.R`)

On the pruned cross-lagged temporal matrix: diagonal set to zero, reciprocal
edges averaged into an undirected weighted graph, zero-weight edges removed.
Spinglass (`implementation = "neg"` for signed edges, at most six communities)
is run within each connected component and repeated over 1,000 random seeds. A
node × node co-occurrence matrix records the proportion of runs in which two
nodes share a community; nodes are ordered by hierarchical clustering on
1 − co-occurrence. Because the graph is collapsed to undirected, community
results do not depend on the PDC / `beta` orientation.

The script also writes the distribution of the number of communities per run
and a table pairing each co-occurrence value with the sign of the direct edge
(if any) and the geodesic distance. High co-occurrence does not imply a direct
edge: partitions are transitive, and under the signed implementation negative
edges push nodes into different communities.

---

## Software

R 4.6.0. Package versions used for the published analysis:

| Package | Version | Purpose |
|---|---|---|
| tidyverse | 2.0.0 | data manipulation and plotting |
| haven, labelled | 2.5.5, 2.16.0 | SPSS import |
| bootnet | 1.8 | GGM estimation, stability, simulation |
| qgraph | 1.9.8 | network visualisation, strength |
| networktools | 1.6.0 | bridge strength |
| NetworkComparisonTest | 2.2.3 | permutation comparison of networks |
| psychonetrics | 0.15 | panel-GVAR (`dlvm1`) |
| igraph | 2.3.2 | Spinglass community detection |
| rstatix, reshape2, fmsb, scales | 0.7.3, 1.4.5, 0.7.6, 1.4.0 | tests and figures |

All stochastic steps set `SEED <- 1234` immediately beforehand.

**Runtime.** The stability bootstraps and recovery simulations in `02_GGMs.R`
dominate (about 1–2 hours with `N_BOOT_NET = 1000` and `N_CORES = 4`). Reduce
`N_BOOT_NET` and `SIM_REPS` in `00_setup.R` for a quick test run.

---

## Data availability

The raw SPSS file (`data/FEP-PrEP-Roma.sav`) contains individual-level clinical
data from the Pr-EP programme and is **not distributed**. Requests for access
should be addressed to the Pr-EP programme, subject to the applicable
data-protection and ethics approvals. `data/residuals.csv`, written by
`01_preprocessing.R`, is derived from individual-level data and is likewise not
included.

---

## References

- Epskamp S. Psychometric network models from time-series and panel data.
  *Psychometrika* 2020;85:206–231.
- Epskamp S, Borsboom D, Fried EI. Estimating psychological networks and their
  accuracy: a tutorial paper. *Behav Res Methods* 2018;50:195–212.
- Jones PJ, Ma R, McNally RJ. Bridge centrality: a network approach to
  understanding comorbidity. *Multivariate Behav Res* 2021;56:353–367.
- van Borkulo CD, van Bork R, Boschloo L, et al. Comparing network structures on
  three aspects: a permutation test. *Psychol Methods* 2023;28:1273–1285.
- Blanken TF, Isvoranu AM, Epskamp S. Estimating network structures using model
  selection. In: *Network Psychometrics with R*. Routledge; 2022:111–132.
- Traag VA, Bruggeman J. Community detection in networks with positive and
  negative links. *Phys Rev E* 2009;80:036115.
- Werner MA, de Ron J, Fried EI, Robinaugh DJ. Iterated community detection in
  psychological networks. *Psychol Methods* 2025.
- Arribas M, Barnby JM, Patel R, et al. Longitudinal evolution of the
  transdiagnostic prodrome to severe mental disorders: a dynamic temporal
  network analysis. *Mol Psychiatry* 2025;30:2931–2942.

---

## Citation

If you use this code, please cite the accompanying article. The full reference
will be added on publication.

## Contact

Tommaso B. Jannini — Department of Experimental Medicine, University of Rome
Tor Vergata.
