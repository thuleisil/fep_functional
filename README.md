# Functional domains and psychopathological dimensions in first-episode psychosis

Analysis code for a two-year longitudinal network study of the relationship between
psychopathological dimensions and functional impairment in first-episode psychosis
(FEP), using data from the **Parma Early Psychosis (Pr-EP) programme**.

The pipeline estimates cross-sectional Gaussian Graphical Models at each wave,
compares them with paired Network Comparison Tests, fits a panel graphical
vector-autoregressive model across the three waves, and examines the community
structure of the resulting temporal network.

---

## Study design

- **Design**: prospective cohort of FEP patients enrolled in the Pr-EP programme,
  assessed at three waves — T0 (baseline), T1 (12 months), T2 (24 months).
- **Nodes (9)**: the five PANSS factor scores of the van der Gaag five-factor model
  (positive, negative, disorganisation, excitement, emotional distress) plus four
  HoNOS items indexing social and functional impairment:

  | Node | Label | Content |
  |------|-------|---------|
  | Pos  | PANSS positive | P1, P3, P5, P6, G9 |
  | Neg  | PANSS negative | N1, N2, N3, N4, N6, G7 |
  | Dis  | PANSS disorganisation | P2, N5, N7, G11 |
  | Exc  | PANSS excitement | P4, P7, G8, G14 |
  | Emo  | PANSS emotional distress | G1, G2, G3, G4, G6 |
  | H9   | HoNOS 9 | Relationships / social functioning |
  | H10  | HoNOS 10 | Activities of daily living |
  | H11  | HoNOS 11 | Living conditions |
  | H12  | HoNOS 12 | Occupational and recreational activities |

  9 nodes × 3 waves = **27 variables**.

- **Covariates (5)**: age, sex, ethnicity, duration of untreated psychosis (DUP) and
  baseline DSM-5 diagnosis. Every node variable is residualised on these five
  covariates and z-standardised before any network is estimated.

---

## Repository structure

```
.
├── R/
│   ├── 00_setup.R                        paths, packages, nodes, labels, constants, helpers
│   ├── 01_preprocessing.R                analytic sample, PANSS factors, residualisation
│   ├── 02_GGMs.R                         cross-sectional networks, centrality, bridge, NCT
│   ├── 03_panel_GVAR.R                   panel-GVAR, model selection, temporal centrality
│   └── 04_community_analysis_spinglass.R iterated Spinglass, co-occurrence
├── data/                                 raw SPSS file (not distributed) + residuals
├── results/
│   ├── figures/
│   └── community/
└── README.md
```

Run the scripts in order, from the repository root:

```r
source("R/00_setup.R")
source("R/01_preprocessing.R")
source("R/02_GGMs.R")
source("R/03_panel_GVAR.R")
source("R/04_community_analysis_spinglass.R")
```

`00_setup.R` must be sourced first in every session; each subsequent script reads
the intermediate files written by the previous one, so they can also be run
independently once `01_preprocessing.R` has produced `data/residuals.csv`.

No script calls `setwd()`. All paths are relative to the repository root.

---

## Methods summary

### Preprocessing (`01_preprocessing.R`)

Complete cases on all PANSS items, HoNOS 9–12 and the five covariates at all three
waves. PANSS factor scores are computed with `na.rm = FALSE`, so a participant with
a missing item contributes `NA` to that factor rather than an artificially low sum;
this also keeps the whole-cohort and analytic-sample tables comparable.

Each of the 27 node variables is then regressed on the five covariates, and the
residuals are centred and z-standardised. No detrending is applied: systematic
change in mean level across waves is left in the data and described by Friedman
tests.

Participant IDs are carried as row names from the analytic sample through
`data/residuals.csv`, so that the wave-specific data frames used later can be
verified to contain the same people in the same order.

### Cross-sectional networks (`02_GGMs.R`)

Regularised partial-correlation networks estimated separately at T0, T1 and T2 with
the graphical LASSO and EBIC model selection (EBICglasso, tuning γ = 0.5), via
`bootnet::estimateNetwork`.

Accuracy and stability: case-dropping subset bootstrap (correlation-stability
coefficient), non-parametric bootstrap of edge weights with centrality difference
tests, and `netSimulator` recovery across sample sizes. 1000 bootstrap iterations
throughout.

Centrality: strength, closeness, betweenness and expected influence, plus **bridge
strength and bridge expected influence** (`networktools::bridge`) computed over two
communities — psychopathology (5 PANSS factors) and functioning (4 HoNOS items).

Two reporting notes carried in the script:

- Report strength and expected influence. Both are sums of edge weights, so their
  stability follows from edge stability, and in a network containing negative edges
  expected influence is the more interpretable of the two. Closeness and betweenness
  rely on shortest paths, which have no clear meaning in a partial-correlation
  network, and are the least stable indices.
- A difference in centrality between two nodes should be claimed only when the
  bootstrapped difference test supports it, never from point estimates alone.

Comparison between waves: `NetworkComparisonTest::NCT` with **`paired = TRUE`**,
1000 permutations, testing network-structure invariance (*M*), global-strength
invariance (*S*), edge-level invariance (FDR-corrected across the 36 possible
edges) and centrality invariance for the standard and bridge indices.

> `paired = TRUE` is not optional here. The three waves are the same participants
> measured three times; the default permutation scheme reshuffles observations
> freely between samples and would treat repeated measurements as independent
> groups, yielding anticonservative p-values. The script verifies that the three
> wave-specific data frames contain the same participant IDs in the same order
> before running the tests.

A non-significant *M* means no single edge differs beyond what permutation would
produce. It does **not** establish that the networks are identical: absence of
evidence for a difference is not evidence of invariance. The NCT also compares
*structure*, not *severity* — a wave in which every domain improves in absolute
terms can still produce an identical network.

### Panel-GVAR (`03_panel_GVAR.R`)

`psychonetrics::dlvm1` with `estimator = "FIML"`, `within_latent = "ggm"` and
`between_latent = "chol"` (Epskamp, 2020, *Psychometrika*). Three networks are
extracted: temporal (directed lag-1), contemporaneous (within-person, same wave)
and between-subjects (stable individual differences).

An **unconstrained** model is fitted first, then a sparse model via
`prune(alpha = 0.05, recursive = TRUE)` followed by `stepup(alpha = 0.05)`. The
sparse model is retained when BIC favours it; the choice is made in code
(`use_sparse`), and every downstream network is extracted from `model_final`.

Two points worth stating explicitly in any write-up:

- **The unconstrained model is not saturated.** `dlvm1` imposes stationarity — the
  same β and the same within-person GGM for the T0→T1 and T1→T2 transitions — which
  consumes degrees of freedom. The remaining df are in `panelgvar_fit.csv`. Use
  "unconstrained" or "fully connected", not "saturated".
- **The model comparison is post-hoc.** Pruning and step-up are data-driven on the
  same data, so the χ² difference test is biased in favour of the sparse model and
  can still be significant when every parsimony-adjusted criterion prefers it. The
  BIC difference is the more defensible basis for the decision.

`panelgvar_temporal_edge_types.csv` is the central descriptive result: it counts
how many temporal edges run from functioning to psychopathology versus the reverse.

### Temporal edge orientation

Temporal edges are extracted as **PDC** (partial directed correlations: standardised
lag-1 coefficients), never as `beta`. This matters:

```r
getmatrix(model, "beta")   # oriented column -> row
getmatrix(model, "PDC")    # oriented row -> column, the qgraph convention
```

Extracting `beta` and plotting or reading it as though it were PDC **reverses every
arrow** and swaps in-strength with out-strength. Throughout this pipeline,
in-strength is `colSums(abs(PDC))` and out-strength is `rowSums(abs(PDC))`, and
directed figures are drawn with a non-zero `curve` so that reciprocal pairs
(A → B and B → A) remain distinguishable. No display threshold is applied: every
edge surviving pruning is shown.

### Community analysis (`04_community_analysis_spinglass.R`)

On the pruned cross-lagged temporal matrix: diagonal set to zero, directed matrix
collapsed to an undirected weighted graph by averaging reciprocal edges, zero-weight
edges removed, then Spinglass (`implementation = "neg"`, `spins = 6`) run separately
within each connected component. Repeated over 1000 random seeds; nodes ordered by
hierarchical clustering on 1 − co-occurrence.

Because the graph is collapsed to undirected, the community results are invariant to
the PDC/`beta` orientation question above.

Because the network is sparse, the graph may split into several connected
components. The distribution of the number of communities per iteration
(`community_n_communities_per_iteration.csv`) should be reported alongside the
heatmap, so readers can judge how much of the partition is driven by disconnection
rather than by modularity.

The script also writes `community_co_occurrence_vs_direct_edges.csv`, which pairs
each co-occurrence value with the corresponding edge sign and geodesic distance. It
exists because co-occurrence and edges diverge for two reasons that change how
results should be worded:

1. **A partition is transitive; an edge is not.** Community detection assigns nodes,
   not edges. If A–B and B–C are connected, all three are typically placed in one
   community, after which A and C co-occur in every iteration despite having no
   direct edge.
2. **Negative edges repel.** Under `implementation = "neg"` (Traag & Bruggeman,
   2009), positive edges pull nodes together while negative edges push them into
   *different* communities. A node joined to a community only by a negative edge is
   assigned away from it, and may end up grouped with nodes to which it has no
   positive connection at all.

A node placed in a community by attraction and one placed there by repulsion occupy
the same cell of the heatmap but mean different things. Before writing that a
functional domain "belongs to" a psychopathological community, check the
`DirectEdge` column. Intermediate co-occurrence values (roughly 0.3–0.7) indicate an
unstable assignment, not a partial association, and such nodes are better reported
as unassigned than allocated to the nearer cluster.

Communities here are defined on the **temporal** structure, so a community is a set
of nodes that tend to influence each other across waves — a dynamic subsystem, not a
cross-sectional symptom cluster.

---

## Software

Developed under R 4.6.0. Core packages:

| Package | Purpose |
|---------|---------|
| `haven`, `tidyverse`, `labelled` | data import and manipulation |
| `bootnet` | GGM estimation, bootstrap, simulation |
| `qgraph` | network visualisation and centrality |
| `networktools` | bridge centrality |
| `NetworkComparisonTest` | permutation comparison of networks |
| `psychonetrics` | panel-GVAR (`dlvm1`) |
| `igraph` | graph objects, Spinglass community detection |
| `rstatix`, `gt`, `reshape2`, `fmsb`, `scales` | tests, tables, figures |

All stochastic steps set `SEED <- 1234` immediately beforehand.

**Runtime.** The GGM bootstraps and recovery simulations in `02_GGMs.R` dominate and
take roughly one to two hours with `N_BOOT_NET = 1000` and `N_CORES = 4`. Reduce
`N_BOOT_NET` in `00_setup.R` for a quick test run, and restore it before producing
final results. The panel-GVAR and the community analysis are comparatively fast.

---

## Data availability

The raw SPSS file (`data/FEP-PrEP-Roma.sav`) contains individual-level clinical data
from the Pr-EP programme and is **not distributed with this repository**. Requests
for access should be directed to the Pr-EP programme, subject to the applicable
data-protection and ethics approvals.

`01_preprocessing.R` writes `data/residuals.csv` (residualised, z-standardised node
scores). This file is derived from individual-level data and is likewise not
included here.

Scripts `02`–`04` run from `data/residuals.csv` alone, so the analysis can be
reproduced without access to the raw file if the residuals are obtained separately.

---

## References

- Epskamp, S. (2020). Psychometric network models from time-series and panel data.
  *Psychometrika*, 85, 206–231.
- Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating psychological networks
  and their accuracy: A tutorial paper. *Behavior Research Methods*, 50, 195–212.
- Jones, P. J., Ma, R., & McNally, R. J. (2021). Bridge centrality: A network
  approach to understanding comorbidity. *Multivariate Behavioral Research*, 56,
  353–367.
- Traag, V. A., & Bruggeman, J. (2009). Community detection in networks with positive
  and negative links. *Physical Review E*, 80, 036115.
- van Borkulo, C. D., van Bork, R., Boschloo, L., et al. (2023). Comparing network
  structures on three aspects: A permutation test. *Psychological Methods*, 28,
  1273–1285.
- van der Gaag, M., Hoffman, T., Remijsen, M., et al. (2006). The five-factor model
  of the Positive and Negative Syndrome Scale II. *Schizophrenia Research*, 85,
  280–287.
- Arribas, M., et al. (2025). *Molecular Psychiatry* — temporal network and
  community-analysis template followed here.
- Blanken, T. F., et al. (2022) — pruning of non-significant edges.
- Werner, M., de Ron, J., Fried, E. I., & Robinaugh, D. J. (2025). *Psychological
  Methods* — iterated community detection and co-occurrence.

> Verify page numbers and volumes against your reference manager before submission;
> the entries above are provided for orientation, not as a formatted bibliography.

---

## Citation

If you use this code, please cite the accompanying manuscript. A full citation will
be added here on publication.

## Contact

Tommaso B. Jannini — Department of Experimental Medicine, University of Rome Tor
Vergata.
