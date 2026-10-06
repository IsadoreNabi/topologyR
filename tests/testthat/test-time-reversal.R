# ======================================================================
# T-1 gate battery (0.4.0): the time-reversal theorem
#
# The referent is EXTERNAL to the engine: the index reflection
# r(i) = n + 1 - i, applied by hand to the objects the engine returns.
# The theorem says that reading a series backwards exchanges the two
# Nada topologies through r, so on a palindrome -- a series equal to
# its own reversal -- r has to be a symmetry carrying the forward
# structure onto the backward one edge by edge, neighbourhood by
# neighbourhood, base set by base set and component by component. The
# asymmetry direction is then zero exactly, since it is a difference of
# two integers, and no tolerance is involved anywhere in the battery.
#
# Only palindromes appear here, and none of them is a repeated ramp
# with an abrupt fall. The point of the restriction is that on a
# palindrome the theorem has a consequence checkable within a single
# call, and that the consequence is checked as an identity of families
# of sets rather than as an equality of counts: a change that broke the
# reflection equivariance of either visibility criterion, of the
# intersection closure or of the component search would surface as a
# mismatch of families even where the counts happened to agree.
# ======================================================================

reflect <- function(idx, n) sort(n + 1L - as.integer(idx))

# Canonical form for comparing families of sets regardless of order.
canon_family <- function(sets) {
  sort(vapply(sets,
              function(s) paste(sort(as.integer(s)), collapse = ","),
              character(1)))
}

# Palindromes: three written by hand, seven mirrored from continuous
# noise, of both parities of length.
palindrome_cases <- function() {
  by_hand <- list(
    c(1, 4, 2, 5, 2, 4, 1),
    c(3, 1, 2, 7, 2, 1, 3),
    c(0, 2, 9, 9, 2, 0)
  )
  set.seed(4075)
  mirrored <- lapply(seq_len(7), function(k) {
    half <- stats::rnorm(sample(3:20, 1))
    if (stats::runif(1) < 0.5) c(half, rev(half)) else c(half, stats::rnorm(1), rev(half))
  })
  c(by_hand, mirrored)
}

PALINDROMES <- palindrome_cases()

visibility_of <- function(series, graph_type) {
  if (graph_type == "hvg") {
    horizontal_visibility_graph(series, directed = TRUE)
  } else {
    natural_visibility_graph(series, directed = TRUE)
  }
}

# ----------------------------------------------------------------------
# The graph level: edges and closed neighbourhoods
# ----------------------------------------------------------------------

test_that("reflection is a symmetry of the visibility graph of a palindrome (T-1)", {
  for (graph_type in c("hvg", "nvg")) {
    for (case in seq_along(PALINDROMES)) {
      series <- PALINDROMES[[case]]
      expect_identical(series, rev(series))
      g <- visibility_of(series, graph_type)
      n <- g$n
      where <- sprintf("%s case %d (n=%d)", graph_type, case, n)

      edges <- unlist(lapply(seq_len(n), function(i) {
        if (length(g$out_adjacency[[i]]) == 0L) character(0)
        else paste(i, g$out_adjacency[[i]], sep = "-")
      }))
      reflected <- unlist(lapply(seq_len(n), function(i) {
        if (length(g$out_adjacency[[i]]) == 0L) character(0)
        else paste(reflect(g$out_adjacency[[i]], n), n + 1L - i, sep = "-")
      }))
      expect_identical(sort(reflected), sort(edges),
                       info = paste(where, "-- the edge set is not r-invariant"))

      for (i in seq_len(n)) {
        forward_nb <- sort(c(i, as.integer(g$out_adjacency[[i]])))
        backward_nb <- sort(c(n + 1L - i, as.integer(g$in_adjacency[[n + 1L - i]])))
        expect_identical(reflect(forward_nb, n), backward_nb,
                         info = sprintf("%s -- r(N+[%d]) is not N-[%d]", where, i, n + 1L - i))
      }
    }
  }
})

# ----------------------------------------------------------------------
# The topological level: bases and components
# ----------------------------------------------------------------------

test_that("reflection exchanges the two Nada topologies of a palindrome (T-1)", {
  for (graph_type in c("hvg", "nvg")) {
    for (case in seq_along(PALINDROMES)) {
      series <- PALINDROMES[[case]]
      bt <- generate_bitopology(series, graph_type = graph_type,
                                max_open_sets = 0L, alexandrov = FALSE)
      n <- length(series)
      where <- sprintf("%s case %d (n=%d)", graph_type, case, n)

      if (isTRUE(bt$forward$base_complete) && isTRUE(bt$backward$base_complete)) {
        expect_identical(canon_family(lapply(bt$forward$base, reflect, n = n)),
                         canon_family(bt$backward$base),
                         info = paste(where, "-- r does not carry B+ onto B-"))
      }
      expect_identical(canon_family(lapply(bt$forward$components, reflect, n = n)),
                       canon_family(bt$backward$components),
                       info = paste(where, "-- r does not carry the components onto each other"))
    }
  }
})

# ----------------------------------------------------------------------
# The reported invariants
# ----------------------------------------------------------------------

test_that("the asymmetry direction of a palindrome is exactly zero (T-1)", {
  for (graph_type in c("hvg", "nvg")) {
    for (case in seq_along(PALINDROMES)) {
      series <- PALINDROMES[[case]]
      inv <- generate_bitopology(series, graph_type = graph_type,
                                 max_open_sets = 0L, alexandrov = FALSE)$invariants
      where <- sprintf("%s case %d (n=%d)", graph_type, case, length(series))

      expect_identical(inv$forward_components, inv$backward_components,
                       info = paste(where, "-- C+ != C-"))
      expect_identical(inv$asymmetry_direction, 0L,
                       info = paste(where, "-- D != 0"))
      expect_identical(inv$irreversibility_components, 0,
                       info = paste(where, "-- I_C != 0"))
      expect_identical(inv$forward_base_size, inv$backward_base_size,
                       info = paste(where, "-- |B+| != |B-|"))
    }
  }
})
