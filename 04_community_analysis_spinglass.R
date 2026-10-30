# =============================================================================
# 04_community_analysis_spinglass.R
#
# Iterated Spinglass community detection on the pruned cross-lagged temporal
# matrix, with a node x node co-occurrence matrix and a diagnostic table
# relating co-occurrence to the underlying temporal edges.
#
#   source("R/00_setup.R")
#   source("R/04_community_analysis_spinglass.R")
#
# Inputs
#   results/panelgvar_temporal_PDC_crosslagged.csv   (from 03_panel_GVAR.R)
#
# Outputs
#   results/community/community_spinglass_memberships.csv
#   results/community/community_co_occurrence_matrix.csv
#   results/community/community_pairwise_co_occurrences.csv
#   results/community/community_n_communities_per_iteration.csv
#   results/community/community_co_occurrence_vs_direct_edges.csv
#   results/community/community_undirected_edges_used.csv
#   results/figures/community_co_occurrence_heatmap.png / .tiff
#
# Procedure, following Arribas et al. (2025, Mol Psychiatry) and the iterated
# community-detection framework of Werner, de Ron, Fried & Robinaugh (2025,
# Psychological Methods):
#   1. set the diagonal (autoregressions) to zero
#   2. collapse the directed matrix into an undirected weighted graph by
#      averaging reciprocal edges
#   3. remove zero-weight edges
#   4. run Spinglass with implementation = "neg" separately within each
#      connected component; "neg" accommodates signed edges (Traag &
#      Bruggeman, 2009, Phys Rev E)
#   5. repeat over N_SPIN_ITER random seeds and build a co-occurrence matrix
#   6. order nodes by hierarchical clustering on 1 - co-occurrence
# =============================================================================

stopifnot(exists("NODES"))


# -----------------------------------------------------------------------------
# 1. Temporal matrix and undirected graph
# -----------------------------------------------------------------------------
temporal_file <- file.path(RESULTS_DIR, "panelgvar_temporal_PDC_crosslagged.csv")

if (!file.exists(temporal_file)) {
  stop("Run 03_panel_GVAR.R first: ", temporal_file, " not found.")
}

features_temporal_cross <- as.matrix(
  read.csv(temporal_file, row.names = 1, check.names = FALSE)
)

labels_pgvar <- colnames(features_temporal_cross)

diag(features_temporal_cross) <- 0

g_directed <- igraph::graph_from_adjacency_matrix(
  features_temporal_cross,
  mode     = "directed",
  weighted = TRUE,
  diag     = FALSE
)

g_undirected <- igraph::as_undirected(
  g_directed,
  mode           = "collapse",
  edge.attr.comb = list(weight = "mean")
)

g_undirected <- igraph::delete_edges(
  g_undirected,
  igraph::E(g_undirected)[igraph::E(g_undirected)$weight == 0]
)

message("Undirected graph: ", igraph::vcount(g_undirected), " nodes, ",
        igraph::ecount(g_undirected), " edges, ",
        igraph::components(g_undirected)$no, " connected component(s)")

write.csv(igraph::as_data_frame(g_undirected, what = "edges"),
          file.path(COMM_DIR, "community_undirected_edges_used.csv"),
          row.names = FALSE)

# Note. Because the graph is collapsed to undirected, the community results do
# not depend on the row/column orientation of the temporal matrix: they would be
# identical whether edges were extracted as PDC or as beta. Orientation matters
# only for the directed figures and for in-/out-strength in 03_panel_GVAR.R.


# -----------------------------------------------------------------------------
# 2. Spinglass within connected components
# -----------------------------------------------------------------------------
# Spinglass requires a connected graph, so it is run separately within each
# component. Isolated nodes and edgeless components become singletons.

run_spinglass_components <- function(graph, spins = N_SPINS, seed = 1) {

  set.seed(seed)

  comp <- igraph::components(graph)
  membership_global <- rep(NA_integer_, igraph::vcount(graph))
  current_comm <- 1L

  for (cc in unique(comp$membership)) {

    node_ids <- which(comp$membership == cc)

    if (length(node_ids) == 1) {
      membership_global[node_ids] <- current_comm
      current_comm <- current_comm + 1L
      next
    }

    subg <- igraph::induced_subgraph(graph, vids = node_ids)

    if (igraph::ecount(subg) == 0) {
      for (n in node_ids) {
        membership_global[n] <- current_comm
        current_comm <- current_comm + 1L
      }
      next
    }

    comm <- tryCatch(
      igraph::cluster_spinglass(
        subg,
        spins          = spins,
        implementation = "neg",
        weights        = igraph::E(subg)$weight
      ),
      error = function(e) NULL
    )

    if (is.null(comm)) {
      for (n in node_ids) {
        membership_global[n] <- current_comm
        current_comm <- current_comm + 1L
      }
    } else {
      sub_membership <- igraph::membership(comm)
      for (m in unique(sub_membership)) {
        membership_global[node_ids[sub_membership == m]] <- current_comm
        current_comm <- current_comm + 1L
      }
    }
  }

  membership_global
}


# -----------------------------------------------------------------------------
# 3. Iterate over random seeds
# -----------------------------------------------------------------------------
node_names_graph <- igraph::V(g_undirected)$name
n_nodes <- length(node_names_graph)

spin_memberships <- matrix(NA_integer_, nrow = N_SPIN_ITER, ncol = n_nodes)
colnames(spin_memberships) <- node_names_graph

for (i in seq_len(N_SPIN_ITER)) {
  spin_memberships[i, ] <- run_spinglass_components(g_undirected, seed = i)
}

write.csv(spin_memberships,
          file.path(COMM_DIR, "community_spinglass_memberships.csv"),
          row.names = FALSE)

n_comm_per_iter <- apply(spin_memberships, 1,
                         function(x) length(unique(na.omit(x))))

n_comm_tbl <- tibble::tibble(k = n_comm_per_iter) %>%
  count(k, name = "n_iterations") %>%
  mutate(Proportion = n_iterations / N_SPIN_ITER)

write.csv(n_comm_tbl,
          file.path(COMM_DIR, "community_n_communities_per_iteration.csv"),
          row.names = FALSE)
print(n_comm_tbl)

# Because the network is sparse, the graph may split into several connected
# components. The distribution above should be reported alongside the heatmap,
# so that readers can judge how much of the partition is driven by disconnection
# rather than by modularity.


# -----------------------------------------------------------------------------
# 4. Co-occurrence matrix
# -----------------------------------------------------------------------------
co_occ_matrix <- matrix(0, n_nodes, n_nodes)

for (r in seq_len(N_SPIN_ITER)) {
  memb <- spin_memberships[r, ]
  same <- outer(memb, memb, FUN = function(a, b) !is.na(a) & !is.na(b) & a == b)
  co_occ_matrix <- co_occ_matrix + same
}

co_occ_matrix <- co_occ_matrix / N_SPIN_ITER
dimnames(co_occ_matrix) <- list(node_names_graph, node_names_graph)

write.csv(co_occ_matrix,
          file.path(COMM_DIR, "community_co_occurrence_matrix.csv"))

co_occ_pairs <- co_occ_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Node1") %>%
  pivot_longer(-Node1, names_to = "Node2", values_to = "CoOccurrence") %>%
  filter(Node1 != Node2) %>%
  rowwise() %>%
  mutate(Pair = paste(sort(c(Node1, Node2)), collapse = " -- ")) %>%
  ungroup() %>%
  distinct(Pair, .keep_all = TRUE) %>%
  arrange(desc(CoOccurrence))

write.csv(co_occ_pairs,
          file.path(COMM_DIR, "community_pairwise_co_occurrences.csv"),
          row.names = FALSE)


# -----------------------------------------------------------------------------
# 5. Heatmap
# -----------------------------------------------------------------------------
hc <- hclust(as.dist(1 - co_occ_matrix))
node_order <- order.dendrogram(as.dendrogram(hc))

reordered <- co_occ_matrix[node_order, node_order]
diag(reordered) <- NA

heat_df <- reshape2::melt(reordered)
heat_df$value <- as.numeric(heat_df$value)

community_heatmap <- ggplot(heat_df, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  coord_fixed() +
  scale_fill_distiller(name = "Co-occurrence", palette = "Blues",
                       direction = 1, limits = c(0, 1), na.value = "white") +
  labs(title = "Community co-occurrence heatmap", x = "Node", y = "Node") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 10))

ggsave(file.path(FIG_DIR, "community_co_occurrence_heatmap.png"),
       community_heatmap, width = 8, height = 7, dpi = 300)
ggsave(file.path(FIG_DIR, "community_co_occurrence_heatmap.tiff"),
       community_heatmap, width = 8, height = 7, dpi = 300)


# -----------------------------------------------------------------------------
# 6. Co-occurrence versus direct edges
# -----------------------------------------------------------------------------
# High co-occurrence between two nodes does NOT imply a direct temporal edge.
# Two mechanisms produce divergence, and both change how results are worded:
#
#   (a) A partition is transitive; an edge is not. Community detection assigns
#       nodes, not edges. If A-B and B-C are connected, all three are typically
#       placed in one community, after which A and C co-occur in every
#       iteration even though no A-C edge exists. Such pairs appear below with
#       DirectEdge = "none" and Geodesic = 2.
#
#   (b) Negative edges repel. With implementation = "neg", positive edges pull
#       nodes into the same community while negative edges push them into
#       DIFFERENT ones. A node joined to a community only by a negative edge is
#       assigned away from it, and ends up grouped with nodes to which it may
#       have no positive connection at all.
#
# Before writing that a functional domain "belongs to" a psychopathological
# community, check the DirectEdge column: if the only relevant edge is negative,
# the accurate statement is that the domain showed an inverse temporal
# association with the other subsystem and was consequently assigned elsewhere.
#
# Intermediate co-occurrence (roughly 0.3-0.7) marks an unstable assignment, not
# a partial association: such nodes should be reported as unassigned rather than
# allocated to the nearer cluster.

und_w <- igraph::as_adjacency_matrix(g_undirected, attr = "weight", sparse = FALSE)
geo   <- igraph::distances(g_undirected, weights = NA)

co_occ_vs_edges <- co_occ_pairs %>%
  mutate(
    EdgeWeight = mapply(function(a, b) und_w[a, b], Node1, Node2),
    DirectEdge = case_when(
      EdgeWeight > 0 ~ "positive",
      EdgeWeight < 0 ~ "negative",
      TRUE           ~ "none"
    ),
    Geodesic = mapply(function(a, b) geo[a, b], Node1, Node2),
    Geodesic = if_else(is.infinite(Geodesic), NA_real_, Geodesic)
  ) %>%
  select(Pair, CoOccurrence, DirectEdge, EdgeWeight, Geodesic) %>%
  arrange(desc(CoOccurrence))

write.csv(co_occ_vs_edges,
          file.path(COMM_DIR, "community_co_occurrence_vs_direct_edges.csv"),
          row.names = FALSE)
print(co_occ_vs_edges, n = Inf)


# Communities here are defined on the TEMPORAL structure, so a community is a
# set of nodes that tend to influence each other across waves - a dynamic
# subsystem, not a cross-sectional symptom cluster. The two results of primary
# interest are (i) whether the psychopathological dimensions form a coherent
# cluster, and (ii) whether any HoNOS functional domain is assigned to the same
# community as a PANSS dimension, and by attraction rather than repulsion.

message("04_community_analysis_spinglass.R complete.")
