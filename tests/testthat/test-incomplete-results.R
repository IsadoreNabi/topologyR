# ======================================================================
# Incomplete results (0.4.0): the six states of a computation
#
# base complete or truncated, times enumeration skipped, complete or
# truncated. The referent for the components is the computation with the
# complete base; the referent for the enumeration is the complete
# enumeration. A truncated base must never certify a topology, and the
# components must not change with the truncation.
# ======================================================================

nvg_example <- function() natural_visibility_graph(c(3, 1, 4, 1, 5, 9, 2, 6))
disconnected_warning <- "already determined to be disconnected"
canon_family <- function(sets) {
  sort(vapply(sets, function(s) paste(sort(as.integer(s)), collapse = ","),
              character(1)))
}

test_that("the reference computation is complete and certified", {
  g <- nvg_example()
  full <- generate_topology(g$adjacency, g$n, max_open_sets = 0L)
  expect_true(full$base_complete)
  expect_identical(length(full$base), 15L)
  expect_identical(full$topology_complete, NA)
  expect_null(full$topology)
  expect_identical(full$n_open_sets, NA_integer_)
  expect_warning(enum <- generate_topology(g$adjacency, g$n, max_open_sets = 100000L,
                                           verify_axioms = TRUE),
                 disconnected_warning)
  expect_true(enum$topology_complete)
  expect_true(enum$axioms_ok)
  expect_identical(enum$n_open_sets, 30L)
})

test_that("a truncated enumeration is not certified", {
  g <- nvg_example()
  expect_warning(cut <- generate_topology(g$adjacency, g$n, max_open_sets = 5L,
                                          verify_axioms = TRUE),
                 disconnected_warning)
  expect_true(cut$base_complete)
  expect_false(cut$topology_complete)
  expect_null(cut$axioms_ok)
})

test_that("a truncated base keeps the exact components and never certifies a topology", {
  g <- nvg_example()
  full <- generate_topology(g$adjacency, g$n, max_open_sets = 0L)
  for (mb in c(8L, 10L, 12L)) {
    cut <- generate_topology(g$adjacency, g$n, max_open_sets = 0L, max_base_sets = mb)
    expect_false(cut$base_complete)
    expect_identical(length(cut$base), mb + 1L)
    expect_identical(cut$topology_complete, NA)
    expect_identical(canon_family(cut$components), canon_family(full$components))
    for (mo in c(20L, 100000L)) {
      expect_warning(enum <- generate_topology(g$adjacency, g$n, max_open_sets = mo,
                                               max_base_sets = mb, verify_axioms = TRUE),
                     disconnected_warning)
      expect_false(enum$topology_complete)
      expect_null(enum$axioms_ok)
    }
  }
  expect_warning(partial <- generate_topology(g$adjacency, g$n, max_open_sets = 100000L,
                                              max_base_sets = 12L),
                 disconnected_warning)
  expect_lt(partial$n_open_sets, 30L)
})

test_that("components stay exact under truncation on random visibility graphs", {
  set.seed(4401)
  truncated <- 0L
  for (k in 1:20) {
    n <- sample(8:30, 1)
    g <- natural_visibility_graph(stats::rnorm(n), directed = TRUE)
    for (adj in list(g$adjacency, g$out_adjacency, g$in_adjacency)) {
      full <- generate_topology(adj, n, max_open_sets = 0L)
      cut <- generate_topology(adj, n, max_open_sets = 0L, max_base_sets = n)
      truncated <- truncated + !cut$base_complete
      expect_identical(canon_family(cut$components), canon_family(full$components))
    }
  }
  expect_gt(truncated, 10L)
})

test_that("the invariants report NA for truncated bases and stop without components", {
  g <- natural_visibility_graph(c(3, 1, 4, 1, 5, 9, 2, 6), directed = TRUE)
  tf <- generate_topology(g$out_adjacency, g$n, max_open_sets = 0L)
  tb <- generate_topology(g$in_adjacency, g$n, max_open_sets = 0L)
  tf_cut <- generate_topology(g$out_adjacency, g$n, max_open_sets = 0L, max_base_sets = 8L)
  alex <- generate_alexandrov_topology(g$out_adjacency, g$n)
  ref <- bitopology_invariants(tf, tb, g$n, alexandrov = alex)
  inv <- bitopology_invariants(tf_cut, tb, g$n, alexandrov = alex)
  expect_false(inv$forward_base_complete)
  expect_true(inv$backward_base_complete)
  expect_identical(inv$forward_base_size, NA_integer_)
  expect_identical(inv$resolution$nada_forward_base_gain, NA_integer_)
  expect_identical(inv$forward_components, ref$forward_components)
  expect_identical(inv$asymmetry_direction, ref$asymmetry_direction)
  expect_identical(ref$forward_base_size, ref$backward_base_size)
  expect_identical(ref$resolution$nada_forward_base_gain,
                   ref$resolution$nada_backward_base_gain)
  no_comp <- generate_topology(g$out_adjacency, g$n, max_open_sets = 0L,
                               check_connected = FALSE)
  expect_error(bitopology_invariants(no_comp, tb, g$n), "check_connected = FALSE")
  expect_error(bitopology_invariants(tf, tb, g$n + 1L), "does not match")
})

test_that("the generator rejects invalid limits", {
  g <- nvg_example()
  expect_error(generate_topology(g$adjacency, g$n, max_open_sets = -1L), "non-negative")
  expect_error(generate_topology(g$adjacency, g$n, max_base_sets = 0L), "positive")
})

test_that("count arguments are validated before any conversion, never truncated", {
  g <- nvg_example()
  tf <- generate_topology(g$adjacency, g$n, max_open_sets = 0L)
  # Fractions, strings, logical values, NA and vectors are rejected.
  expect_error(generate_topology(g$adjacency, 8.5, max_open_sets = 0L), "n_elements")
  expect_error(generate_topology(g$adjacency, g$n, max_open_sets = 2.7), "max_open_sets")
  expect_error(generate_topology(g$adjacency, g$n, max_open_sets = "100"), "max_open_sets")
  expect_error(generate_topology(g$adjacency, g$n, max_open_sets = TRUE), "max_open_sets")
  expect_error(generate_topology(g$adjacency, g$n, max_open_sets = NA), "max_open_sets")
  expect_error(generate_topology(g$adjacency, g$n, max_open_sets = c(1, 2)), "max_open_sets")
  expect_error(generate_topology(g$adjacency, g$n, max_open_sets = Inf), "max_open_sets")
  expect_error(generate_topology(g$adjacency, g$n, max_open_sets = 0L, max_base_sets = 3.9),
               "max_base_sets")
  expect_error(is_topology_connected_exact(list(integer(0), 1:3), 3.2), "n_elements")
  expect_error(bitopology_invariants(tf, tf, 8.4), "n_elements")
  expect_error(generate_alexandrov_topology(g$adjacency, 8.5), "n_elements")
  expect_error(generate_alexandrov_topology(g$adjacency, g$n, max_open_sets = -1),
               "max_open_sets")
  expect_error(generate_bitopology(c(3, 1, 4), max_open_sets = 10.5), "max_open_sets")
  # Whole numbers stored as doubles are accepted and give the same result.
  expect_identical(generate_topology(g$adjacency, 8, max_open_sets = 0)$components,
                   tf$components)
})

test_that("the flags are logical, NA included, on every engine", {
  g <- nvg_example()
  a <- generate_topology(g$adjacency, g$n, max_open_sets = 0L, check_connected = FALSE)
  expect_identical(a$connected, NA)
  expect_identical(a$topology_complete, NA)
  expect_identical(a$components, list())
  n <- 200L
  adj <- lapply(seq_len(n), function(v) setdiff(seq_len(n), v))
  b <- generate_topology(adj, n, max_open_sets = 0L, check_connected = FALSE)
  expect_identical(b$connected, NA)
  expect_identical(b$topology_complete, NA)
  al <- generate_alexandrov_topology(g$adjacency, g$n,
                                     check_connected = FALSE)
  expect_identical(al$connected, NA)
  expect_identical(al$topology_complete, NA)
  expect_identical(generate_topology(g$adjacency, g$n, max_open_sets = 0L)$connected,
                   FALSE)
})
