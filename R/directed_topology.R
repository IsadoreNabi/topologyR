# ----------------------------------------------------------------------
# Internal helpers for the generalized Alexandrov construction (C-40)
# ----------------------------------------------------------------------

# Neighbors the C++ engine would use: valid 1-based indices only.
# The engine silently ignores out-of-range entries; the R-side analysis
# must see exactly the digraph the engine sees.
alexandrov_valid_neighbors <- function(out_adjacency, n_elements) {
  lapply(out_adjacency, function(nb) {
    nb <- as.integer(nb)
    nb[!is.na(nb) & nb >= 1L & nb <= n_elements]
  })
}

# Does every edge go from a lower to a higher index? That is the contract
# the reverse-order propagation of the engine actually relies on -- it is
# strictly stronger than acyclicity.
alexandrov_index_ordered <- function(adj) {
  for (v in seq_along(adj)) {
    if (any(adj[[v]] <= v)) return(FALSE)
  }
  TRUE
}

# Iterative Tarjan: strongly connected components of the digraph.
# Returns comp_of (1-based membership, components emitted sinks-first,
# i.e. in reverse topological order of the condensation) and n_comp.
alexandrov_tarjan_scc <- function(adj, n) {
  index_of <- integer(n)      # 0 = unvisited
  lowlink <- integer(n)
  on_stack <- logical(n)
  scc_stack <- integer(0)
  comp_of <- integer(n)
  n_comp <- 0L
  counter <- 0L

  for (root in seq_len(n)) {
    if (index_of[root] != 0L) next
    # Explicit DFS stack: each frame is (vertex, next-neighbor position)
    stack_v <- integer(0)
    stack_i <- integer(0)
    stack_v <- c(stack_v, root)
    stack_i <- c(stack_i, 1L)
    counter <- counter + 1L
    index_of[root] <- counter
    lowlink[root] <- counter
    scc_stack <- c(scc_stack, root)
    on_stack[root] <- TRUE

    while (length(stack_v) > 0L) {
      top <- length(stack_v)
      v <- stack_v[top]
      i <- stack_i[top]
      nb <- adj[[v]]
      if (i <= length(nb)) {
        stack_i[top] <- i + 1L
        w <- nb[i]
        if (index_of[w] == 0L) {
          counter <- counter + 1L
          index_of[w] <- counter
          lowlink[w] <- counter
          scc_stack <- c(scc_stack, w)
          on_stack[w] <- TRUE
          stack_v <- c(stack_v, w)
          stack_i <- c(stack_i, 1L)
        } else if (on_stack[w]) {
          lowlink[v] <- min(lowlink[v], index_of[w])
        }
      } else {
        # Post-order: fold lowlink into the parent, pop SCC if root
        stack_v <- stack_v[-top]
        stack_i <- stack_i[-top]
        if (length(stack_v) > 0L) {
          p <- stack_v[length(stack_v)]
          lowlink[p] <- min(lowlink[p], lowlink[v])
        }
        if (lowlink[v] == index_of[v]) {
          n_comp <- n_comp + 1L
          repeat {
            w <- scc_stack[length(scc_stack)]
            scc_stack <- scc_stack[-length(scc_stack)]
            on_stack[w] <- FALSE
            comp_of[w] <- n_comp
            if (w == v) break
          }
        }
      }
    }
  }
  list(comp_of = comp_of, n_comp = n_comp)
}

#' Generate the Alexandrov topology of the reachability preorder of a digraph
#'
#' @description
#' Computes the Alexandrov topology induced by a directed graph via the
#' reachability preorder. The open sets are the upsets of reachability:
#' vertex \eqn{v}'s minimal open set is the set of all vertices reachable
#' from \eqn{v} via directed paths (including \eqn{v} itself).
#'
#' The input is \strong{any} directed graph. Three regimes are handled,
#' and the result reports which one applied:
#' \enumerate{
#'   \item \strong{Index-ordered DAG} (every edge goes from a lower to a
#'     higher index -- what directed visibility graphs produce): the
#'     original single-pass bitset engine runs unchanged, bit-identical
#'     to versions up to 0.2.0.
#'   \item \strong{DAG in arbitrary vertex order}: the vertices are
#'     re-indexed along a topological order first. Versions up to 0.2.0
#'     required index order but did not check it, and returned wrong
#'     upsets on such inputs without any warning.
#'   \item \strong{Digraph with directed cycles}: mutually reachable
#'     vertices have identical upsets by definition of the preorder, so
#'     the topology is computed on the condensation (Tarjan, 1972) and
#'     expanded back to the vertices. Versions up to 0.2.0 returned a
#'     valid-looking but wrong topology on cyclic inputs.
#' }
#'
#' The Alexandrov topology is a \strong{lower bound} on the resolution of
#' the Nada et al. (2018) topology: \eqn{\tau_A \subseteq \tau^+_{\mathrm{Nada}}}.
#' The difference in base sizes quantifies how much additional structure
#' the intersection/union closure of Nada et al. captures beyond the pure
#' order structure.
#'
#' @param out_adjacency A list of integer vectors, where element \code{i}
#'   contains the 1-based indices of the direct successors of node
#'   \code{i}. Any directed graph is accepted; see \emph{Description} for
#'   the three regimes. Typically obtained from
#'   \code{\link{horizontal_visibility_graph}} or
#'   \code{\link{natural_visibility_graph}} with \code{directed = TRUE}.
#' @param n_elements Integer. Total number of elements in the ground set V.
#' @param max_open_sets Integer. Maximum number of open sets to enumerate
#'   (default: 0, meaning skip enumeration). The Alexandrov topology can have
#'   up to \eqn{2^n} open sets, so enumeration is only feasible for small
#'   \eqn{n}. Connectivity is always computed from the base without
#'   enumeration.
#' @param check_connected Logical. Whether to compute exact topological
#'   connectivity via the specialization preorder (default: \code{TRUE}).
#' @param verify_axioms Logical. Whether to verify topology axioms after
#'   enumeration (default: \code{FALSE}). Only meaningful when the topology
#'   is fully enumerated.
#' @param expect Character. What the caller asserts about the input:
#'   \code{"any"} (default) accepts every digraph and dispatches to the
#'   right regime; \code{"index_ordered_dag"} \emph{verifies} the strict
#'   contract of the visibility-graph case -- every edge from a lower to a
#'   higher index -- and stops with an error if it does not hold. The
#'   check is \eqn{O(m)}. Use it when the caller's own semantics require
#'   the strict case (e.g. time-indexed series), so that a malformed
#'   input fails loudly instead of being silently generalized.
#' @return A \code{list} with the same structure as
#'   \code{\link{generate_topology}}:
#'   \describe{
#'     \item{subbase}{List of integer vectors. The distinct upsets
#'       (= base for Alexandrov).}
#'     \item{base}{List of integer vectors. Same as subbase (upsets are already
#'       closed under intersection).}
#'     \item{base_complete}{Always \code{TRUE} (no iteration needed).}
#'     \item{connected}{Logical. Exact topological connectivity.}
#'     \item{components}{List of integer vectors. Exact connected components.}
#'     \item{topology}{List of integer vectors, or \code{NULL}. Open sets if
#'       enumerated.}
#'     \item{n_open_sets}{Integer, or \code{NA}. Number of open sets if enumerated.}
#'     \item{complete}{Logical. Whether enumeration finished.}
#'     \item{input_index_ordered}{Logical. Whether every edge of the input
#'       went from a lower to a higher index (regime 1).}
#'     \item{input_acyclic}{Logical. Whether the input had no directed
#'       cycle.}
#'     \item{collapsed_classes}{List of integer vectors: the non-trivial
#'       strongly connected components, each a class of topologically
#'       indistinguishable vertices. \code{list()} when there are none.
#'       A non-empty value is a finding about the system, not an error:
#'       within each class, the topology cannot distinguish.}
#'   }
#'
#' @details
#' For regime 1 the algorithm computes reachability in
#' \eqn{O(n \cdot m / 64)} time by processing vertices in reverse index
#' order (= reverse time order for visibility graphs) and propagating
#' reachability sets as multi-word bitsets with OR operations. For
#' regimes 2 and 3 the strongly connected components are found with
#' Tarjan's algorithm (\eqn{O(n + m)}), the condensation is re-indexed
#' along a topological order, the same bitset engine runs on the
#' condensation, and the resulting sets are expanded back to the original
#' vertices. Since vertices in one strongly connected component are
#' topologically indistinguishable, the expansion is a bijection between
#' the open sets of the condensation and the open sets of the input, so
#' every reported field (connectivity, components, enumeration,
#' axiom checks) transfers exactly.
#'
#' The Alexandrov topology on a finite set equipped with a preorder is
#' canonical: there is a bijection between preorders and Alexandrov
#' topologies (Alexandrov, 1937). For directed visibility graphs, the
#' reachability preorder captures the temporal ordering structure of the
#' time series.
#'
#' @references
#' Alexandrov, P. (1937). Diskrete Raume. \emph{Matematicheskii Sbornik},
#' 2(44), 501-519.
#'
#' Tarjan, R. E. (1972). Depth-first search and linear graph algorithms.
#' \emph{SIAM Journal on Computing}, 1(2), 146-160.
#'
#' @seealso \code{\link{generate_topology}} for the Nada et al. topology,
#'   \code{\link{generate_bitopology}} for the complete bitopological analysis.
#'
#' @examples
#' series <- c(3, 1, 4, 1, 5)
#' g <- horizontal_visibility_graph(series, directed = TRUE)
#' alex <- generate_alexandrov_topology(g$out_adjacency, g$n)
#' alex$connected
#' length(alex$base)
#'
#' # A directed cycle: all three vertices are mutually reachable, so the
#' # topology collapses them into one class
#' cyc <- generate_alexandrov_topology(list(2L, 3L, 1L), 3L)
#' cyc$collapsed_classes
#'
#' @export
generate_alexandrov_topology <- function(out_adjacency, n_elements,
                                         max_open_sets = 0L,
                                         check_connected = TRUE,
                                         verify_axioms = FALSE,
                                         expect = c("any", "index_ordered_dag")) {
  expect <- match.arg(expect)
  if (!is.list(out_adjacency)) {
    stop("'out_adjacency' must be a list of integer vectors.", call. = FALSE)
  }
  n_elements <- as.integer(n_elements)
  if (n_elements < 1L) {
    stop("'n_elements' must be >= 1.", call. = FALSE)
  }
  if (length(out_adjacency) != n_elements) {
    stop("Length of 'out_adjacency' must equal 'n_elements'.", call. = FALSE)
  }
  max_open_sets <- as.integer(max_open_sets)
  enumerate <- max_open_sets > 0L

  adj <- alexandrov_valid_neighbors(out_adjacency, n_elements)
  index_ordered <- alexandrov_index_ordered(adj)

  if (expect == "index_ordered_dag" && !index_ordered) {
    offender <- NULL
    for (v in seq_along(adj)) {
      bad <- adj[[v]][adj[[v]] <= v]
      if (length(bad)) { offender <- c(v, bad[1L]); break }
    }
    stop("'expect = \"index_ordered_dag\"' but the input has the edge ",
         offender[1L], " -> ", offender[2L],
         ", which does not go from a lower to a higher index. ",
         "The strict visibility-graph contract does not hold; pass ",
         "'expect = \"any\"' if generalization to the reachability ",
         "preorder is intended.", call. = FALSE)
  }

  if (index_ordered) {
    # Regime 1: the original engine, unchanged and bit-identical.
    result <- alexandrov_topology_engine(
      out_adjacency = out_adjacency,
      n_elements = n_elements,
      max_open_sets = max_open_sets,
      enumerate_topology = enumerate,
      check_connected = check_connected,
      verify_axioms = verify_axioms
    )
    result$input_index_ordered <- TRUE
    result$input_acyclic <- TRUE
    result$collapsed_classes <- list()
    return(result)
  }

  # Regimes 2 and 3: condense, re-index topologically, run the same
  # engine on the condensation, expand back to the vertices.
  scc <- alexandrov_tarjan_scc(adj, n_elements)
  n_comp <- scc$n_comp
  # Tarjan emits components sinks-first (reverse topological order), so
  # reversing the emission index yields a topological order in which
  # every condensation edge goes from a lower to a higher class index.
  topo_id <- n_comp + 1L - scc$comp_of

  members <- vector("list", n_comp)
  for (v in seq_len(n_elements)) {
    cid <- topo_id[v]
    members[[cid]] <- c(members[[cid]], v)
  }
  members <- lapply(members, sort)

  cond_adj <- vector("list", n_comp)
  for (v in seq_len(n_elements)) {
    cv <- topo_id[v]
    cw <- topo_id[adj[[v]]]
    cw <- cw[cw != cv]
    if (length(cw)) {
      cond_adj[[cv]] <- c(cond_adj[[cv]], cw)
    }
  }
  cond_adj <- lapply(cond_adj, function(x) {
    if (is.null(x)) integer(0) else sort(unique(as.integer(x)))
  })

  cond <- alexandrov_topology_engine(
    out_adjacency = cond_adj,
    n_elements = n_comp,
    max_open_sets = max_open_sets,
    enumerate_topology = enumerate,
    check_connected = check_connected,
    verify_axioms = verify_axioms
  )

  # Expansion: a set of classes becomes the sorted set of their vertices.
  # It is a bijection between the condensation's open sets and the
  # input's open sets, so counts and flags transfer unchanged.
  expand <- function(class_set) {
    if (length(class_set) == 0L) return(integer(0))
    sort(unlist(members[class_set], use.names = FALSE))
  }

  result <- cond
  result$subbase <- lapply(cond$subbase, expand)
  result$base <- lapply(cond$base, expand)
  result$components <- lapply(cond$components, expand)
  if (!is.null(cond$topology)) {
    result$topology <- lapply(cond$topology, expand)
  }

  collapsed <- members[lengths(members) > 1L]
  if (length(collapsed)) {
    collapsed <- collapsed[order(vapply(collapsed, min, integer(1)))]
  }
  result$input_index_ordered <- FALSE
  result$input_acyclic <- (n_comp == n_elements)
  result$collapsed_classes <- collapsed
  result
}


#' Generate a bitopological analysis from a time series
#'
#' @description
#' Performs a complete bitopological analysis of a time series by constructing
#' directed visibility graphs and computing four topological spaces:
#' \enumerate{
#'   \item \strong{Undirected Nada topology} \eqn{\tau}: from the symmetric
#'     (undirected) visibility graph (standard Nada et al. 2018).
#'   \item \strong{Forward Nada topology} \eqn{\tau^+}: from closed forward
#'     neighborhoods \eqn{N^+[v] = \{v\} \cup \{w : v \to w\}}.
#'   \item \strong{Backward Nada topology} \eqn{\tau^-}: from closed backward
#'     neighborhoods \eqn{N^-[v] = \{v\} \cup \{w : w \to v\}}.
#'   \item \strong{Alexandrov topology} \eqn{\tau_A}: from reachability upsets
#'     of the directed graph.
#' }
#'
#' The pair \eqn{(X, \tau^+, \tau^-)} forms a \strong{bitopological space}
#' (Kelly, 1963). The divergence between \eqn{\tau^+} and \eqn{\tau^-}
#' quantifies temporal irreversibility: in a reversible process,
#' \eqn{\tau^+ \cong \tau^-}; in an irreversible process (e.g., asymmetric
#' business cycles), they diverge.
#'
#' @param series Numeric vector representing the time series.
#' @param graph_type Character. Either \code{"hvg"} (Horizontal Visibility
#'   Graph) or \code{"nvg"} (Natural Visibility Graph). Default: \code{"hvg"}.
#' @param max_open_sets Integer. Maximum open sets to enumerate per topology
#'   (default: 0, skip enumeration). Required for pairwise connectedness check.
#' @param max_base_sets Integer. Maximum base elements during intersection
#'   closure for the Nada topologies (default: 100,000).
#' @param alexandrov Logical. Whether to also compute the Alexandrov topology
#'   for resolution comparison (default: \code{TRUE}).
#' @return A \code{list} with:
#'   \describe{
#'     \item{graph}{The directed visibility graph (output of
#'       \code{\link{horizontal_visibility_graph}} or
#'       \code{\link{natural_visibility_graph}} with \code{directed = TRUE}).}
#'     \item{undirected}{The undirected Nada topology
#'       (output of \code{\link{generate_topology}}).}
#'     \item{forward}{The forward Nada topology \eqn{\tau^+}.}
#'     \item{backward}{The backward Nada topology \eqn{\tau^-}.}
#'     \item{alexandrov}{The Alexandrov topology \eqn{\tau_A} (if requested).}
#'     \item{invariants}{Bitopological invariants (output of
#'       \code{\link{bitopology_invariants}}).}
#'   }
#'
#' @details
#' The existing undirected engine is called three times (undirected, forward,
#' backward) and the Alexandrov engine once. No new C++ code is needed for the
#' forward and backward topologies: passing directed adjacency lists to the
#' existing \code{\link{generate_topology}} produces the correct closed
#' neighborhoods \eqn{N^+[v]} or \eqn{N^-[v]} automatically (the engine
#' always adds the self-loop \eqn{v \in N[v]}).
#'
#' \strong{The Alexandrov branch is constant for this input class.} Because
#' consecutive observations always see each other, every directed
#' visibility graph contains the path \eqn{1 \to 2 \to \dots \to n};
#' reachability is therefore the total order, and the Alexandrov topology
#' is the same chain of upsets \eqn{\{i, \dots, n\}} for \emph{every}
#' series. It carries no information about the particular series analyzed
#' -- it is a property of the visibility-graph class (see the
#' corresponding section in \code{\link{bitopology_invariants}}). Callers
#' that consume only the Nada side can set \code{alexandrov = FALSE} and
#' skip that cost; the Alexandrov branch earns its keep on general
#' digraphs passed directly to
#' \code{\link{generate_alexandrov_topology}}, where it does vary with
#' the input.
#'
#' @references
#' Kelly, J. C. (1963). Bitopological spaces. \emph{Proceedings of the
#' London Mathematical Society}, 3(1), 71-89.
#'
#' Lacasa, L., & Toral, R. (2010). Description of stochastic and chaotic
#' series using visibility graphs. \emph{Physical Review E}, 82(3), 036120.
#'
#' @seealso \code{\link{bitopology_invariants}} for computing invariants
#'   from pre-computed topologies.
#'
#' @examples
#' series <- c(3, 1, 4, 1, 5, 9, 2, 6)
#' bt <- generate_bitopology(series, graph_type = "hvg")
#' bt$invariants$forward_components
#' bt$invariants$backward_components
#' bt$invariants$irreversibility_components
#'
#' @export
generate_bitopology <- function(series,
                                graph_type = c("hvg", "nvg"),
                                max_open_sets = 0L,
                                max_base_sets = 100000L,
                                alexandrov = TRUE) {
  graph_type <- match.arg(graph_type)

  if (!is.numeric(series)) {
    stop("'series' must be a numeric vector.", call. = FALSE)
  }
  if (anyNA(series)) {
    stop("'series' must not contain NA values.", call. = FALSE)
  }
  if (!all(is.finite(series))) {
    stop("'series' must not contain Inf or NaN values.", call. = FALSE)
  }
  if (length(series) < 2L) {
    stop("'series' must have at least 2 observations.", call. = FALSE)
  }

  # Build directed visibility graph
  if (graph_type == "hvg") {
    g <- horizontal_visibility_graph(series, directed = TRUE)
  } else {
    g <- natural_visibility_graph(series, directed = TRUE)
  }

  n <- g$n

  # Undirected topology (standard Nada et al.)
  tau_undirected <- generate_topology(g$adjacency, n,
                                      max_open_sets = max_open_sets,
                                      max_base_sets = max_base_sets)

  # Forward topology τ+ (from N+[v])
  tau_forward <- generate_topology(g$out_adjacency, n,
                                   max_open_sets = max_open_sets,
                                   max_base_sets = max_base_sets)

  # Backward topology τ- (from N-[v])
  tau_backward <- generate_topology(g$in_adjacency, n,
                                    max_open_sets = max_open_sets,
                                    max_base_sets = max_base_sets)

  # Alexandrov topology (optional)
  tau_alexandrov <- NULL
  if (alexandrov) {
    tau_alexandrov <- generate_alexandrov_topology(
      g$out_adjacency, n,
      max_open_sets = max_open_sets
    )
  }

  # Compute invariants
  inv <- bitopology_invariants(tau_forward, tau_backward, n,
                               tau_undirected, tau_alexandrov)

  list(
    graph = g,
    undirected = tau_undirected,
    forward = tau_forward,
    backward = tau_backward,
    alexandrov = tau_alexandrov,
    invariants = inv
  )
}


#' Compute bitopological invariants from forward and backward topologies
#'
#' @description
#' Given the forward topology \eqn{\tau^+} and backward topology \eqn{\tau^-}
#' (as returned by \code{\link{generate_topology}}), computes bitopological
#' invariants that quantify temporal irreversibility and structural asymmetry.
#'
#' @param forward The forward Nada topology \eqn{\tau^+} (output of
#'   \code{\link{generate_topology}} called with forward adjacency).
#' @param backward The backward Nada topology \eqn{\tau^-} (output of
#'   \code{\link{generate_topology}} called with backward adjacency).
#' @param n_elements Integer. Total number of elements in the ground set V.
#' @param undirected The undirected Nada topology (optional, for comparison).
#' @param alexandrov The Alexandrov topology (optional, for resolution
#'   comparison).
#' @return A \code{list} with:
#'   \describe{
#'     \item{forward_components}{Integer. Number of connected components in
#'       \eqn{\tau^+}.}
#'     \item{backward_components}{Integer. Number of connected components in
#'       \eqn{\tau^-}.}
#'     \item{forward_connected}{Logical. Whether \eqn{\tau^+} is connected.}
#'     \item{backward_connected}{Logical. Whether \eqn{\tau^-} is connected.}
#'     \item{forward_base_size}{Integer. Number of base elements in
#'       \eqn{\tau^+}.}
#'     \item{backward_base_size}{Integer. Number of base elements in
#'       \eqn{\tau^-}.}
#'     \item{irreversibility_components}{Numeric in the range 0 to 1. Normalized
#'       asymmetry of component counts:
#'       \eqn{|C^+ - C^-| / \max(C^+, C^-)}.
#'       Zero means equal component counts (necessary for reversibility).
#'       Note this is a coarse reading: it compares \emph{counts}, and two
#'       partitions can differ while their counts agree.}
#'     \item{asymmetry_direction}{Integer. \eqn{C^- - C^+}. Positive means
#'       the forward topology is more connected (fewer components) than the
#'       backward topology, consistent with gradual expansions (forward
#'       visibility) and abrupt contractions (backward obstruction).}
#'     \item{pairwise}{A list with pairwise connectedness results. Contains
#'       \code{pairwise_connected} (logical or \code{NA}),
#'       \code{clopen_forward} (integer vector, the forward-open / backward-
#'       closed set if disconnected), and \code{clopen_backward} (its
#'       complement). Requires both topologies to be fully enumerated
#'       (\code{max_open_sets > 0} and \code{complete = TRUE}).}
#'     \item{resolution}{A list comparing the Nada and Alexandrov topologies
#'       (if \code{alexandrov} is provided). Contains base size comparisons
#'       and relative resolution gains.}
#'     \item{undirected_info}{A list with undirected topology summary (if
#'       \code{undirected} is provided).}
#'   }
#'
#' @details
#' \strong{Why there is no \code{irreversibility_base} field.} Versions up
#' to 0.2.0 reported the normalized asymmetry of the two base sizes,
#' \eqn{||B^+| - |B^-|| / \max}. That quantity is \strong{identically zero
#' by a theorem of the construction}, for every digraph -- not only for
#' visibility graphs -- so it measures nothing and was removed. The two
#' base sizes remain reported, since each is informative on its own.
#'
#' \emph{Theorem.} Let \eqn{R} be any binary relation on a finite set
#' \eqn{V}. For \eqn{v \in V} let \eqn{F(v) = \{w : v R w\}} (row) and
#' \eqn{C(v) = \{w : w R v\}} (column), and let \eqn{\mathcal{B}(F)} be
#' the set of non-empty intersections of non-empty subfamilies of
#' \eqn{\{F(v)\}}, likewise \eqn{\mathcal{B}(C)}. Then
#' \eqn{|\mathcal{B}(F)| = |\mathcal{B}(C)|}.
#'
#' \emph{Proof sketch.} Consider the formal context \eqn{K = (V, V, R)}.
#' The intents of \eqn{K} are the intersections \eqn{\cap_{v \in A} F(v)}
#' over all \eqn{A \subseteq V} (with \eqn{A = \emptyset} giving \eqn{V});
#' the extents are the same over columns. The Galois connection induced by
#' \eqn{R} gives the classical bijection between intents and extents
#' (Ganter & Wille, 1999; the duality is Birkhoff's, 1940), so
#' \eqn{|\mathrm{Int}| = |\mathrm{Ext}|}. Adjusting the two boundary cases
#' -- whether \eqn{V} itself is an intersection of rows, and whether some
#' non-empty subfamily intersects to the empty set -- contributes the same
#' correction on both sides, because a full row on one side is exactly what
#' prevents an empty intersection on the other. Subtracting,
#' \eqn{|\mathcal{B}(F)| = |\mathcal{B}(C)|}.
#'
#' The forward and backward closed neighborhoods are rows and columns of
#' the reflexive closure of the visibility relation, and the Nada closure
#' generates all finite intersections, so
#' \code{forward_base_size == backward_base_size} whenever both closures
#' ran to completion (\code{base_complete}). The only exception is
#' truncation: if \code{max_base_sets} cuts one side, both sides report
#' \code{base_complete = FALSE} and the sizes are not comparable.
#'
#' \strong{The Alexandrov branch is constant on visibility graphs.} In a
#' visibility graph -- horizontal or natural -- consecutive observations
#' always see each other, so the digraph contains every edge
#' \eqn{i \to i+1}, reachability is the total order \eqn{1 < 2 < \dots < n},
#' and the Alexandrov topology is always the same chain of upsets
#' \eqn{\{i, \dots, n\}}, independent of the series. Consequently, for
#' visibility-graph input: \code{alexandrov_base_size} \eqn{= n},
#' \code{alexandrov_components} \eqn{= 1}, \code{alexandrov_connected}
#' \eqn{=} \code{TRUE}, and the resolution-gain fields reduce to
#' \code{fwd/bwd_base_size - n}: all of their information comes from the
#' Nada side. These fields are genuinely informative for \emph{general}
#' digraphs, where the Alexandrov topology varies with the input; for
#' visibility graphs they are properties of the graph class, not
#' measurements of the series.
#'
#' \strong{Pairwise connectedness} (Kelly, 1963): the bitopological space
#' \eqn{(X, \tau^+, \tau^-)} is pairwise disconnected if there exists a
#' proper non-empty subset \eqn{A \subset X} that is simultaneously
#' \eqn{\tau^+}-open and \eqn{\tau^-}-closed (its complement is
#' \eqn{\tau^-}-open). This is checked by iterating over all
#' \eqn{\tau^+} open sets and testing whether their complement is
#' \eqn{\tau^-}-open.
#'
#' \strong{Irreversibility prediction} (Lacasa & Toral, 2010): for a time
#' series with asymmetric dynamics (e.g., gradual expansions and abrupt
#' recessions), the forward topology \eqn{\tau^+} should be more connected
#' than \eqn{\tau^-}, because forward visibility is less obstructed during
#' gradual rises than backward visibility after abrupt drops.
#'
#' \strong{Resolution gain}: the difference
#' \eqn{|\mathcal{B}_{\mathrm{Nada}}| - |\mathcal{B}_A|} quantifies how
#' much additional structure the Nada intersection/union closure captures
#' beyond the pure Alexandrov order structure. A large gain means the
#' Nada pipeline is extracting genuinely new topological information.
#'
#' @references
#' Kelly, J. C. (1963). Bitopological spaces. \emph{Proceedings of the
#' London Mathematical Society}, 3(1), 71-89.
#'
#' Ganter, B., & Wille, R. (1999). \emph{Formal Concept Analysis:
#' Mathematical Foundations}. Springer.
#'
#' Birkhoff, G. (1940). \emph{Lattice Theory}. American Mathematical
#' Society.
#'
#' @seealso \code{\link{generate_bitopology}} for the complete pipeline.
#'
#' @examples
#' series <- c(3, 1, 4, 1, 5, 9, 2, 6)
#' g <- horizontal_visibility_graph(series, directed = TRUE)
#' tf <- generate_topology(g$out_adjacency, g$n, max_open_sets = 0L)
#' tb <- generate_topology(g$in_adjacency, g$n, max_open_sets = 0L)
#' inv <- bitopology_invariants(tf, tb, g$n)
#' inv$asymmetry_direction
#'
#' @export
bitopology_invariants <- function(forward, backward, n_elements,
                                  undirected = NULL, alexandrov = NULL) {
  if (!is.list(forward) || is.null(forward$components)) {
    stop("'forward' must be a topology result from generate_topology().",
         call. = FALSE)
  }
  if (!is.list(backward) || is.null(backward$components)) {
    stop("'backward' must be a topology result from generate_topology().",
         call. = FALSE)
  }
  n_elements <- as.integer(n_elements)

  # --- Component analysis ---
  fwd_n_comp <- length(forward$components)
  bwd_n_comp <- length(backward$components)
  fwd_connected <- isTRUE(forward$connected)
  bwd_connected <- isTRUE(backward$connected)

  # --- Base analysis ---
  fwd_base_size <- length(forward$base)
  bwd_base_size <- length(backward$base)

  # --- Irreversibility indices ---
  max_comp <- max(fwd_n_comp, bwd_n_comp)
  irreversibility_components <- if (max_comp > 0L) {
    abs(fwd_n_comp - bwd_n_comp) / max_comp
  } else {
    0
  }

  # No base-size asymmetry index is computed: ||B+|-|B-||/max is
  # identically zero by the intents-extents duality (see @details), so
  # it does not measure. The two sizes are reported separately.

  # Positive = forward more connected (fewer components)
  asymmetry_direction <- bwd_n_comp - fwd_n_comp

  # --- Pairwise connectedness (Kelly, 1963) ---
  pairwise <- list(
    pairwise_connected = NA,
    clopen_forward = NULL,
    clopen_backward = NULL
  )

  fwd_enumerated <- !is.null(forward$topology) && isTRUE(forward$complete)
  bwd_enumerated <- !is.null(backward$topology) && isTRUE(backward$complete)

  if (fwd_enumerated && bwd_enumerated) {
    # Encode τ- open sets as character keys for O(1) lookup
    bwd_env <- new.env(hash = TRUE, parent = emptyenv())
    for (s in backward$topology) {
      if (length(s) == 0L || length(s) == n_elements) next
      key <- paste(sort(s), collapse = ",")
      bwd_env[[key]] <- TRUE
    }

    all_elements <- seq_len(n_elements)
    found_clopen <- FALSE

    for (U in forward$topology) {
      if (length(U) == 0L || length(U) == n_elements) next
      complement <- setdiff(all_elements, U)
      comp_key <- paste(sort(complement), collapse = ",")
      if (!is.null(bwd_env[[comp_key]])) {
        pairwise$pairwise_connected <- FALSE
        pairwise$clopen_forward <- sort(U)
        pairwise$clopen_backward <- sort(complement)
        found_clopen <- TRUE
        break
      }
    }

    if (!found_clopen) {
      pairwise$pairwise_connected <- TRUE
    }
  }

  # --- Resolution comparison with Alexandrov ---
  resolution <- NULL
  if (!is.null(alexandrov) && is.list(alexandrov)) {
    alex_base_size <- length(alexandrov$base)
    alex_n_comp <- length(alexandrov$components)

    resolution <- list(
      alexandrov_base_size = alex_base_size,
      alexandrov_components = alex_n_comp,
      alexandrov_connected = isTRUE(alexandrov$connected),
      nada_forward_base_gain = fwd_base_size - alex_base_size,
      nada_backward_base_gain = bwd_base_size - alex_base_size,
      nada_forward_relative_gain = (fwd_base_size - alex_base_size) /
        max(alex_base_size, 1L),
      nada_backward_relative_gain = (bwd_base_size - alex_base_size) /
        max(alex_base_size, 1L)
    )
  }

  # --- Undirected comparison ---
  undirected_info <- NULL
  if (!is.null(undirected) && is.list(undirected)) {
    undirected_info <- list(
      undirected_components = length(undirected$components),
      undirected_connected = isTRUE(undirected$connected),
      undirected_base_size = length(undirected$base)
    )
  }

  list(
    forward_components = fwd_n_comp,
    backward_components = bwd_n_comp,
    forward_connected = fwd_connected,
    backward_connected = bwd_connected,
    forward_base_size = fwd_base_size,
    backward_base_size = bwd_base_size,
    irreversibility_components = irreversibility_components,
    asymmetry_direction = asymmetry_direction,
    pairwise = pairwise,
    resolution = resolution,
    undirected_info = undirected_info
  )
}
