# ======================================================================
# Tests for directed visibility graphs, Alexandrov topology,
# bitopological analysis, and invariants
# ======================================================================

# ------------------------------------------------------------------
# 1. Directed visibility graphs
# ------------------------------------------------------------------

test_that("HVG directed=FALSE returns no directed adjacency", {
  g <- horizontal_visibility_graph(c(3, 1, 4, 1, 5))
  expect_null(g$out_adjacency)
  expect_null(g$in_adjacency)
})

test_that("HVG directed=TRUE returns out_adjacency and in_adjacency", {
  g <- horizontal_visibility_graph(c(3, 1, 4, 1, 5), directed = TRUE)
  expect_true(!is.null(g$out_adjacency))
  expect_true(!is.null(g$in_adjacency))
  expect_length(g$out_adjacency, g$n)
  expect_length(g$in_adjacency, g$n)
})

test_that("HVG directed edges satisfy from < to (DAG property)", {
  g <- horizontal_visibility_graph(c(3, 1, 4, 1, 5), directed = TRUE)
  # All out-neighbors of v should have index > v
  for (v in seq_len(g$n)) {
    out_nb <- g$out_adjacency[[v]]
    if (length(out_nb) > 0) {
      expect_true(all(out_nb > v),
                  info = paste("out_adjacency of", v, "has backward edges"))
    }
    in_nb <- g$in_adjacency[[v]]
    if (length(in_nb) > 0) {
      expect_true(all(in_nb < v),
                  info = paste("in_adjacency of", v, "has forward edges"))
    }
  }
})

test_that("NVG directed=TRUE returns valid directed adjacency", {
  g <- natural_visibility_graph(c(3, 1, 4, 1, 5), directed = TRUE)
  expect_true(!is.null(g$out_adjacency))
  expect_true(!is.null(g$in_adjacency))
  # DAG property
  for (v in seq_len(g$n)) {
    if (length(g$out_adjacency[[v]]) > 0)
      expect_true(all(g$out_adjacency[[v]] > v))
    if (length(g$in_adjacency[[v]]) > 0)
      expect_true(all(g$in_adjacency[[v]] < v))
  }
})

test_that("directed adjacency is consistent with undirected", {
  g <- horizontal_visibility_graph(c(3, 1, 4, 1, 5), directed = TRUE)
  # For each v, sort(union(out_adj[v], in_adj[v])) == sort(adj[v])
  for (v in seq_len(g$n)) {
    directed_union <- sort(c(g$out_adjacency[[v]], g$in_adjacency[[v]]))
    undirected <- sort(g$adjacency[[v]])
    expect_equal(directed_union, undirected,
                 info = paste("vertex", v))
  }
})

test_that("directed graph of single observation", {
  g <- horizontal_visibility_graph(c(42), directed = TRUE)
  expect_length(g$out_adjacency, 1L)
  expect_equal(g$out_adjacency[[1]], integer(0))
  expect_equal(g$in_adjacency[[1]], integer(0))
})

test_that("HVG directed: known edge structure for c(3, 1, 4, 1, 5)", {
  # Manual computation:
  # Edges: (1,2), (1,3), (2,3), (3,4), (3,5), (4,5)
  g <- horizontal_visibility_graph(c(3, 1, 4, 1, 5), directed = TRUE)
  expect_equal(sort(g$out_adjacency[[1]]), c(2L, 3L))
  expect_equal(g$out_adjacency[[2]], 3L)
  expect_equal(sort(g$out_adjacency[[3]]), c(4L, 5L))
  expect_equal(g$out_adjacency[[4]], 5L)
  expect_equal(g$out_adjacency[[5]], integer(0))
  expect_equal(g$in_adjacency[[1]], integer(0))
  expect_equal(g$in_adjacency[[2]], 1L)
  expect_equal(sort(g$in_adjacency[[3]]), c(1L, 2L))
  expect_equal(g$in_adjacency[[4]], 3L)
  expect_equal(sort(g$in_adjacency[[5]]), c(3L, 4L))
})


# ------------------------------------------------------------------
# 2. Alexandrov topology
# ------------------------------------------------------------------

test_that("Alexandrov topology: path graph 1->2->3", {
  # out_adjacency for path: 1->[2], 2->[3], 3->[]
  out_adj <- list(2L, 3L, integer(0))
  alex <- generate_alexandrov_topology(out_adj, 3L)

  # reach[3]={3}, reach[2]={2,3}, reach[1]={1,2,3}
  # Base = {{3}, {2,3}, {1,2,3}} = 3 distinct sets
  expect_equal(length(alex$base), 3L)
  expect_true(alex$connected)
  expect_equal(length(alex$components), 1L)
  expect_true(alex$base_complete)
})

test_that("Alexandrov topology: enumeration produces correct open sets", {
  out_adj <- list(2L, 3L, integer(0))
  alex <- generate_alexandrov_topology(out_adj, 3L, max_open_sets = 100L)

  # Topology = {empty, {3}, {2,3}, {1,2,3}} = 4 open sets
  expect_equal(alex$n_open_sets, 4L)
  expect_true(alex$topology_complete)
})

test_that("Alexandrov topology: disconnected DAG", {
  # Two isolated components: 1->2, 3->4
  out_adj <- list(2L, integer(0), 4L, integer(0))
  alex <- generate_alexandrov_topology(out_adj, 4L)

  expect_false(alex$connected)
  expect_equal(length(alex$components), 2L)
})

test_that("Alexandrov topology: single vertex", {
  out_adj <- list(integer(0))
  alex <- generate_alexandrov_topology(out_adj, 1L)

  expect_equal(length(alex$base), 1L)
  expect_true(alex$connected)
})

test_that("Alexandrov topology: the transitive tournament on three vertices has three base sets", {
  out_adj <- list(c(2L, 3L), 3L, integer(0))
  alex <- generate_alexandrov_topology(out_adj, 3L)

  # reach[3]={3}, reach[2]={2,3}, reach[1]={1,2,3} → 3 base elements
  expect_equal(length(alex$base), 3L)
  expect_true(alex$connected)
})

test_that("Alexandrov topology: star DAG (hub reaches all)", {
  # 1->{2,3,4,5}, 2->{}, 3->{}, 4->{}, 5->{}
  out_adj <- list(c(2L, 3L, 4L, 5L), integer(0), integer(0),
                  integer(0), integer(0))
  alex <- generate_alexandrov_topology(out_adj, 5L)

  # reach[2]={2}, reach[3]={3}, reach[4]={4}, reach[5]={5}, reach[1]={1,2,3,4,5}
  # Base has 5 elements. Singletons are pairwise incomparable.
  # S(1)={V} ⊆ S(k) for any k? S(1) has only base[4]={1,2,3,4,5}.
  # S(2) has base[0]={2} and base[4]={1,2,3,4,5}. S(1) ⊆ S(2)? Yes.
  # So 1 is comparable with everyone → connected.
  expect_equal(length(alex$base), 5L)
  expect_true(alex$connected)
})

test_that("every reachable set is open in the forward topology (tau_A within tau+)", {
  set.seed(4301)
  for (type in c("hvg", "nvg")) {
    for (k in 1:15) {
      n <- sample(3:25, 1)
      g <- if (type == "hvg") horizontal_visibility_graph(stats::rnorm(n), directed = TRUE)
           else natural_visibility_graph(stats::rnorm(n), directed = TRUE)
      alex <- generate_alexandrov_topology(g$out_adjacency, g$n)
      nada_fwd <- generate_topology(g$out_adjacency, g$n, max_open_sets = 0L)
      for (U in alex$base) {
        inside <- Filter(function(B) all(B %in% U), nada_fwd$base)
        expect_setequal(unique(unlist(inside)), U)
      }
      expect_gte(length(nada_fwd$base), length(alex$base))
    }
  }
})

test_that("Alexandrov verify_axioms passes for small example", {
  out_adj <- list(2L, 3L, integer(0))
  alex <- generate_alexandrov_topology(out_adj, 3L,
                                       max_open_sets = 100L,
                                       verify_axioms = TRUE)
  expect_true(alex$axioms_ok)
})


# ------------------------------------------------------------------
# 3. Forward and backward Nada topologies
# ------------------------------------------------------------------

test_that("Nada-forward and Nada-backward from directed graph", {
  g <- horizontal_visibility_graph(c(3, 1, 4, 1, 5), directed = TRUE)

  tau_fwd <- generate_topology(g$out_adjacency, g$n, max_open_sets = 0L)
  tau_bwd <- generate_topology(g$in_adjacency, g$n, max_open_sets = 0L)

  # Both should have valid results
  expect_true(!is.null(tau_fwd$connected))
  expect_true(!is.null(tau_bwd$connected))
  expect_true(length(tau_fwd$base) > 0)
  expect_true(length(tau_bwd$base) > 0)
})

test_that("forward and backward base sizes are equal on an asymmetric series (duality)", {
  series <- c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5)
  for (type in c("hvg", "nvg")) {
    inv <- generate_bitopology(series, graph_type = type)$invariants
    expect_true(inv$forward_base_complete && inv$backward_base_complete)
    expect_identical(inv$forward_base_size, inv$backward_base_size)
    expect_identical(inv$resolution$nada_forward_base_gain,
                     inv$resolution$nada_backward_base_gain)
  }
})


# ------------------------------------------------------------------
# 4. Bitopological analysis (full pipeline)
# ------------------------------------------------------------------

test_that("generate_bitopology returns all expected components", {
  bt <- generate_bitopology(c(3, 1, 4, 1, 5, 9, 2, 6))
  expect_true(!is.null(bt$graph))
  expect_true(!is.null(bt$undirected))
  expect_true(!is.null(bt$forward))
  expect_true(!is.null(bt$backward))
  expect_true(!is.null(bt$alexandrov))
  expect_true(!is.null(bt$invariants))
})

test_that("generate_bitopology with NVG", {
  bt <- generate_bitopology(c(3, 1, 4, 1, 5), graph_type = "nvg")
  expect_true(!is.null(bt$graph))
  expect_true(!is.null(bt$invariants))
})

test_that("generate_bitopology without Alexandrov", {
  bt <- generate_bitopology(c(3, 1, 4, 1, 5), alexandrov = FALSE)
  expect_null(bt$alexandrov)
  expect_null(bt$invariants$resolution)
})

test_that("generate_bitopology input validation", {
  expect_error(generate_bitopology("not numeric"))
  expect_error(generate_bitopology(c(1, NA, 3)))
  expect_error(generate_bitopology(c(1, Inf)))
  expect_error(generate_bitopology(42))  # length 1
})


# ------------------------------------------------------------------
# 5. Bitopological invariants
# ------------------------------------------------------------------

test_that("bitopology_invariants returns all expected fields", {
  g <- horizontal_visibility_graph(c(3, 1, 4, 1, 5), directed = TRUE)
  tf <- generate_topology(g$out_adjacency, g$n, max_open_sets = 0L)
  tb <- generate_topology(g$in_adjacency, g$n, max_open_sets = 0L)
  inv <- bitopology_invariants(tf, tb, g$n)

  expect_true(!is.null(inv$forward_components))
  expect_true(!is.null(inv$backward_components))
  expect_true(!is.null(inv$forward_connected))
  expect_true(!is.null(inv$backward_connected))
  expect_true(!is.null(inv$irreversibility_components))
  expect_true(!is.null(inv$asymmetry_direction))
  expect_true(!is.null(inv$pairwise))
  # Retired in 0.3.0: identically zero by the intents-extents duality
  expect_null(inv$irreversibility_base)
})

test_that("irreversibility index is in valid range", {
  bt <- generate_bitopology(c(3, 1, 4, 1, 5, 9, 2, 6))
  inv <- bt$invariants
  expect_true(inv$irreversibility_components >= 0)
  expect_true(inv$irreversibility_components <= 1)
})

test_that("a palindrome has zero asymmetry and zero I_C, exactly", {
  for (type in c("hvg", "nvg")) {
    inv <- generate_bitopology(c(1, 3, 5, 3, 1), graph_type = type)$invariants
    expect_identical(inv$asymmetry_direction, 0L)
    expect_identical(inv$irreversibility_components, 0)
  }
})

test_that("asymmetry_direction is C- minus C+", {
  bt <- generate_bitopology(c(3, 1, 4, 1, 5, 9, 2, 6))
  inv <- bt$invariants
  # asymmetry_direction = bwd_comp - fwd_comp
  expect_equal(inv$asymmetry_direction,
               inv$backward_components - inv$forward_components)
})

test_that("resolution gain is non-negative", {
  bt <- generate_bitopology(c(3, 1, 4, 1, 5, 9, 2, 6))
  res <- bt$invariants$resolution
  expect_true(!is.null(res))
  expect_true(res$nada_forward_base_gain >= 0,
              info = "Nada base must be >= Alexandrov base")
  expect_true(res$nada_backward_base_gain >= 0)
})


# ------------------------------------------------------------------
# 6. Pairwise connectedness
# ------------------------------------------------------------------

# The engine warns when it enumerates a space it already found disconnected;
# the tests below enumerate on purpose, so only that warning is muffled.
enumerate_quietly <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("already determined to be disconnected", conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  })
}

test_that("pairwise connectedness is decided without enumeration", {
  y <- c(3, 1, 4, 1, 5)
  pw <- generate_bitopology(y, max_open_sets = 0L)$invariants$pairwise
  expect_false(pw$pairwise_connected)
  # The witness is checked against the complete enumerations.
  full <- enumerate_quietly(generate_bitopology(y, max_open_sets = 10000L))
  key <- function(s) paste(sort(s), collapse = ",")
  expect_true(key(pw$clopen_forward) %in% vapply(full$forward$topology, key, ""))
  expect_true(key(pw$clopen_backward) %in% vapply(full$backward$topology, key, ""))
  expect_setequal(c(pw$clopen_forward, pw$clopen_backward), seq_along(y))
  expect_identical(full$invariants$pairwise, pw)
})

# Algorithm 5 of the accompanying article, kept as the second oracle of the
# exact criterion: it looks up the complement of each enumerated proper
# forward-open set among the enumerated proper backward-open sets, answers
# FALSE with the first witness found, TRUE when there is none and both
# enumerations are certified complete, and NA otherwise.
algorithm5 <- function(tf, tb, n) {
  key <- function(s) paste(sort(s), collapse = ",")
  proper_bwd <- tb$topology[lengths(tb$topology) > 0L & lengths(tb$topology) < n]
  bwd <- vapply(proper_bwd, key, "")
  for (U in tf$topology) {
    if (length(U) == 0L || length(U) == n) next
    W <- setdiff(seq_len(n), U)
    if (key(W) %in% bwd) {
      return(list(pairwise_connected = FALSE, clopen_forward = sort(U),
                  clopen_backward = W))
    }
  }
  complete <- isTRUE(tf$topology_complete) && isTRUE(tb$topology_complete)
  list(pairwise_connected = if (complete) TRUE else NA)
}

# The exact criterion on one digraph against algorithm 5 on its complete
# enumerations: the same answer, and a witness made of enumerated open sets
# that partition the vertices.
expect_algorithm5_agrees <- function(out_adj, in_adj, n, info) {
  key <- function(s) paste(sort(s), collapse = ",")
  tf <- enumerate_quietly(generate_topology(out_adj, n,
                                            max_open_sets = as.integer(2^n)))
  tb <- enumerate_quietly(generate_topology(in_adj, n,
                                            max_open_sets = as.integer(2^n)))
  expect_true(isTRUE(tf$topology_complete) && isTRUE(tb$topology_complete),
              info = info)
  pw <- bitopology_invariants(tf, tb, n)$pairwise
  expect_identical(pw$pairwise_connected, algorithm5(tf, tb, n)$pairwise_connected,
                   info = info)
  if (isFALSE(pw$pairwise_connected)) {
    expect_true(key(pw$clopen_forward) %in% vapply(tf$topology, key, ""), info = info)
    expect_true(key(pw$clopen_backward) %in% vapply(tb$topology, key, ""), info = info)
    expect_setequal(c(pw$clopen_forward, pw$clopen_backward), seq_len(n))
    expect_length(intersect(pw$clopen_forward, pw$clopen_backward), 0L)
  }
  pw$pairwise_connected
}

test_that("the exact criterion agrees with algorithm 5 on every digraph with at most three vertices", {
  seen <- c(conn = 0L, disc = 0L)
  for (n in 1:3) {
    arcs <- expand.grid(j = seq_len(n), i = seq_len(n))
    arcs <- arcs[arcs$i != arcs$j, ]
    for (code in seq_len(2^nrow(arcs)) - 1L) {
      on <- bitwAnd(code, as.integer(2^(seq_len(nrow(arcs)) - 1L))) > 0L
      out_adj <- lapply(seq_len(n), function(v) arcs$j[on & arcs$i == v])
      in_adj <- lapply(seq_len(n), function(v) arcs$i[on & arcs$j == v])
      ans <- expect_algorithm5_agrees(out_adj, in_adj, n, sprintf("n = %d, code %d", n, code))
      if (isTRUE(ans)) seen[["conn"]] <- seen[["conn"]] + 1L else seen[["disc"]] <- seen[["disc"]] + 1L
    }
  }
  expect_identical(sum(seen), 69L)
  expect_gt(seen[["conn"]], 0L)
  expect_gt(seen[["disc"]], 0L)
})

test_that("the exact criterion agrees with algorithm 5 on random digraphs with cycles", {
  set.seed(4203)
  seen <- c(conn = 0L, disc = 0L)
  for (r in 1:300) {
    n <- sample(2:7, 1)
    A <- matrix(stats::runif(n * n) < stats::runif(1), n)
    diag(A) <- FALSE
    out_adj <- lapply(seq_len(n), function(v) which(A[v, ]))
    in_adj <- lapply(seq_len(n), function(v) which(A[, v]))
    ans <- expect_algorithm5_agrees(out_adj, in_adj, n, sprintf("case %d", r))
    if (isTRUE(ans)) seen[["conn"]] <- seen[["conn"]] + 1L else seen[["disc"]] <- seen[["disc"]] + 1L
  }
  expect_gt(seen[["conn"]], 0L)
  expect_gt(seen[["disc"]], 0L)
})

test_that("a directed visibility graph is pairwise disconnected by every split into segments", {
  # Referent: the theorem in ?bitopology_invariants. Every arc goes forward in
  # time, so each final segment is forward-open and each initial segment is
  # backward-open; the segments are looked up in the complete enumerations.
  set.seed(4201)
  key <- function(s) paste(sort(s), collapse = ",")
  for (r in 1:40) {
    n <- sample(2:8, 1)
    y <- sample(0:4, n, TRUE) + 0
    bt <- enumerate_quietly(generate_bitopology(
      y, graph_type = if (r %% 2) "hvg" else "nvg",
      max_open_sets = 100000L, alexandrov = FALSE))
    expect_true(isTRUE(bt$forward$topology_complete) &&
                  isTRUE(bt$backward$topology_complete))
    fwd <- vapply(bt$forward$topology, key, "")
    bwd <- vapply(bt$backward$topology, key, "")
    for (k in seq_len(n - 1L)) {
      expect_true(key((k + 1L):n) %in% fwd, info = sprintf("case %d, k = %d", r, k))
      expect_true(key(seq_len(k)) %in% bwd, info = sprintf("case %d, k = %d", r, k))
    }
    pw <- bt$invariants$pairwise
    expect_false(pw$pairwise_connected)
    expect_true(key(pw$clopen_forward) %in% fwd, info = sprintf("case %d", r))
    expect_true(key(pw$clopen_backward) %in% bwd, info = sprintf("case %d", r))
    # {n} is forward-open and not backward-open, so the two topologies always
    # differ (see ?generate_bitopology).
    expect_true(key(n) %in% fwd)
    expect_false(key(n) %in% bwd)
  }
})

test_that("a digraph with cycles can be pairwise connected, with or without enumeration", {
  # Arcs 1->2, 2->1, 2->3, 3->1, 3->2: the forward topology is
  # {{}, {1,2}, V}, the backward one {{}, {2,3}, V}, and {3} is not
  # backward-open, so no witness exists.
  out_adj <- list(2L, c(1L, 3L), c(1L, 2L))
  in_adj <- list(c(2L, 3L), c(1L, 3L), 2L)
  tf <- generate_topology(out_adj, 3L, max_open_sets = 100L)
  tb <- generate_topology(in_adj, 3L, max_open_sets = 100L)
  key <- function(s) paste(sort(s), collapse = ",")
  expect_setequal(vapply(tf$topology, key, ""), c("", "1,2", "1,2,3"))
  expect_setequal(vapply(tb$topology, key, ""), c("", "2,3", "1,2,3"))
  expect_true(bitopology_invariants(tf, tb, 3L)$pairwise$pairwise_connected)
  tf0 <- generate_topology(out_adj, 3L, max_open_sets = 0L)
  tb0 <- generate_topology(in_adj, 3L, max_open_sets = 0L)
  expect_true(bitopology_invariants(tf0, tb0, 3L)$pairwise$pairwise_connected)
})


# ------------------------------------------------------------------
# 7. Nada-forward vs Alexandrov resolution (concrete example)
# ------------------------------------------------------------------

test_that("Nada-forward has strictly more base elements than Alexandrov for path 1->2->3", {
  # Nada-forward: N+[1]={1,2}, N+[2]={2,3}, N+[3]={3}
  # Base = {{1,2}, {2,3}, {3}, {2}} (intersection {1,2}∩{2,3}={2})
  # Alexandrov: reach[1]={1,2,3}, reach[2]={2,3}, reach[3]={3}
  # Base = {{1,2,3}, {2,3}, {3}}
  out_adj <- list(2L, 3L, integer(0))

  nada <- generate_topology(out_adj, 3L, max_open_sets = 0L)
  alex <- generate_alexandrov_topology(out_adj, 3L)

  # Nada should have strictly more base elements
  expect_true(length(nada$base) > length(alex$base),
              info = paste("Nada:", length(nada$base),
                           "Alex:", length(alex$base)))
})


# ------------------------------------------------------------------
# 8. Existing tests still pass (regression)
# ------------------------------------------------------------------

test_that("existing undirected topology still works after refactor", {
  series <- c(3, 1, 4, 1, 5)
  g <- horizontal_visibility_graph(series)
  topo <- generate_topology(g$adjacency, g$n)
  expect_true(!is.null(topo$connected))
  expect_true(length(topo$base) > 0)
})

test_that("K6 complete graph still gives indiscrete topology (N[v] check)", {
  n <- 6
  adj <- lapply(seq_len(n), function(v) setdiff(seq_len(n), v))
  topo <- generate_topology(adj, n, max_open_sets = 100L)
  # With N[v], complete graph gives indiscrete topology
  expect_equal(topo$n_open_sets, 2L)
  expect_true(topo$connected)
})
