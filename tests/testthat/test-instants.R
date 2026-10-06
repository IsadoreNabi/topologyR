# ======================================================================
# Unequally spaced instants (0.4.0)
#
# (a) times = NULL is the equispaced case: identical results for 1..n,
#     for Date instants one day apart, and for any exact affine image
#     of 1..n.
# (b) Theorem T1 with instants: reversing the values and reflecting the
#     instants (t -> -t, read backwards) exchanges the two Nada
#     topologies through r(i) = n + 1 - i, on random series with random
#     irregular instants.
# (c) Negative control, the instance of the formal statement
#     nEdgeT_not_equivariant: instants 0, 1, 4 and values 0, 2, 4.
# (d) Exact affine changes of the time scale and of the values leave
#     the natural graph unchanged; the horizontal graph ignores the
#     instants and is invariant under strictly increasing maps of the
#     values.
# (e) Input validation.
# ======================================================================

edge_set <- function(g) sort(paste(g$edges$from, g$edges$to, sep = "-"))
reflect_set <- function(s, n) sort(n + 1L - as.integer(s))
canon_family <- function(sets) {
  sort(vapply(sets, function(s) paste(sort(as.integer(s)), collapse = ","),
              character(1)))
}
visibility <- function(type, x, times) {
  if (type == "hvg") horizontal_visibility_graph(x, directed = TRUE, times = times)
  else natural_visibility_graph(x, directed = TRUE, times = times)
}

same_bitopology <- function(a, b, where) {
  expect_identical(edge_set(a$graph), edge_set(b$graph), info = where)
  for (side in c("undirected", "forward", "backward")) {
    expect_identical(canon_family(a[[side]]$base), canon_family(b[[side]]$base),
                     info = paste(where, side, "base"))
    expect_identical(canon_family(a[[side]]$components),
                     canon_family(b[[side]]$components),
                     info = paste(where, side, "components"))
  }
  expect_identical(a$invariants$asymmetry_direction,
                   b$invariants$asymmetry_direction, info = where)
}

test_that("times = NULL is the equispaced case (a)", {
  set.seed(4201)
  series <- c(list(c(1, 1, 1), c(1, 0.999, 1), c(0, 1, 2, 3)),
              lapply(1:12, function(k) stats::rnorm(sample(5:30, 1))))
  for (type in c("hvg", "nvg")) {
    for (k in seq_along(series)) {
      x <- series[[k]]
      n <- length(x)
      ref <- generate_bitopology(x, graph_type = type)
      where <- sprintf("%s series %d", type, k)
      same_bitopology(ref, generate_bitopology(x, graph_type = type, times = 1:n), where)
      same_bitopology(ref, generate_bitopology(x, graph_type = type,
                                               times = as.Date("2020-01-01") + 0:(n - 1)),
                      paste(where, "Date"))
      same_bitopology(ref, generate_bitopology(x, graph_type = type,
                                               times = 4 * (1:n) + 1024),
                      paste(where, "affine"))
    }
  }
})

test_that("reversal with reflected instants exchanges the two topologies (b, T1)", {
  set.seed(4202)
  for (type in c("hvg", "nvg")) {
    for (k in 1:25) {
      n <- sample(4:40, 1)
      x <- stats::rnorm(n)
      t <- cumsum(stats::rexp(n))
      rv <- time_reverse(x, t)
      expect_identical(rv$series, rev(x))
      expect_identical(rv$times, -rev(t))
      g <- visibility(type, x, t)
      gr <- visibility(type, rv$series, rv$times)
      where <- sprintf("%s case %d (n = %d)", type, k, n)
      reflected_edges <- unlist(lapply(seq_len(n), function(i) {
        if (!length(gr$out_adjacency[[i]])) character(0)
        else paste(reflect_set(gr$out_adjacency[[i]], n), n + 1L - i, sep = "-")
      }))
      expect_identical(sort(reflected_edges), edge_set(g), info = where)
      for (i in seq_len(n)) {
        fwd_rev <- c(i, gr$out_adjacency[[i]])
        bwd_orig <- c(n + 1L - i, g$in_adjacency[[n + 1L - i]])
        expect_identical(reflect_set(fwd_rev, n), sort(as.integer(bwd_orig)),
                         info = paste(where, "N+ onto N-, vertex", i))
      }
      bt <- generate_bitopology(x, graph_type = type, times = t, alexandrov = FALSE)
      br <- generate_bitopology(rv$series, graph_type = type, times = rv$times,
                                alexandrov = FALSE)
      expect_identical(canon_family(lapply(br$forward$base, reflect_set, n = n)),
                       canon_family(bt$backward$base), info = paste(where, "B+ onto B-"))
      expect_identical(canon_family(lapply(br$forward$components, reflect_set, n = n)),
                       canon_family(bt$backward$components),
                       info = paste(where, "components"))
      expect_identical(br$invariants$forward_components, bt$invariants$backward_components)
      expect_identical(br$invariants$backward_components, bt$invariants$forward_components)
      expect_identical(br$invariants$asymmetry_direction, -bt$invariants$asymmetry_direction)
      expect_identical(br$invariants$irreversibility_components,
                       bt$invariants$irreversibility_components)
    }
  }
})

test_that("reversing the values alone is not equivariant (c, negative control)", {
  t <- c(0, 1, 4)
  x <- c(0, 2, 4)
  original <- natural_visibility_graph(x, times = t)
  values_only <- natural_visibility_graph(rev(x), times = t)
  reflected <- natural_visibility_graph(rev(x), times = -rev(t))
  expect_identical(edge_set(original), c("1-2", "2-3"))
  expect_identical(edge_set(values_only), c("1-2", "1-3", "2-3"))
  expect_false(identical(edge_set(values_only), edge_set(original)))
  expect_identical(edge_set(reflected), edge_set(original))
  expect_identical(edge_set(natural_visibility_graph(x)), c("1-2", "2-3"))
  expect_identical(edge_set(natural_visibility_graph(rev(x))), c("1-2", "2-3"))
})

test_that("exact affine maps leave the graphs unchanged (d)", {
  set.seed(4203)
  for (k in 1:20) {
    n <- sample(4:30, 1)
    t <- cumsum(sample(1:6, n, TRUE)) / 4
    x <- sample(-40:40, n, TRUE) / 8
    ref_n <- edge_set(natural_visibility_graph(x, times = t))
    expect_identical(edge_set(natural_visibility_graph(x, times = 8 * t - 512)), ref_n)
    expect_identical(edge_set(natural_visibility_graph(16 * x + 32, times = t)), ref_n)
    ref_h <- edge_set(horizontal_visibility_graph(x))
    expect_identical(edge_set(horizontal_visibility_graph(x, times = t)), ref_h)
    expect_identical(edge_set(horizontal_visibility_graph(exp(x))), ref_h)
  }
})

test_that("the instants are validated (e)", {
  x <- c(1, 3, 2, 5)
  expect_error(natural_visibility_graph(x, times = c(1, 2, 2, 3)), "strictly increasing")
  expect_error(natural_visibility_graph(x, times = c(1, 3, 2, 4)), "strictly increasing")
  expect_error(natural_visibility_graph(x, times = 1:3), "same length")
  expect_error(natural_visibility_graph(x, times = c(1, 2, NA, 4)), "finite")
  expect_error(natural_visibility_graph(x, times = c(1, 2, Inf, 4)), "finite")
  expect_error(natural_visibility_graph(x, times = letters[1:4]), "numeric")
  expect_error(horizontal_visibility_graph(x, times = c(4, 3, 2, 1)), "strictly increasing")
  expect_error(generate_bitopology(x, times = c(0, 1)), "same length")
  expect_identical(natural_visibility_graph(x)$times, as.double(1:4))
  expect_identical(natural_visibility_graph(x, times = as.Date("2024-03-01") + c(0, 2, 3, 9))$times,
                   as.numeric(as.Date("2024-03-01") + c(0, 2, 3, 9)))
})
