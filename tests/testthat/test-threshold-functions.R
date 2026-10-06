# ---- Threshold function tests ----

test_that("calculate_thresholds returns all five methods", {
  th <- calculate_thresholds(rnorm(50))
  expect_named(th, c("mean_diff", "median_diff", "sd", "iqr", "dbscan"))
  expect_true(all(vapply(th, is.numeric, logical(1))))
  expect_true(all(vapply(th, function(x) x > 0, logical(1))))
})

test_that("calculate_thresholds validates input", {
  expect_error(calculate_thresholds("text"), "numeric")
  expect_error(calculate_thresholds(1), "at least 2")
  expect_error(calculate_thresholds(c(1, NA)), "NA")
})

test_that("calculate_topology returns positive integer", {
  bs <- calculate_topology(rnorm(20), threshold = 0.5)
  expect_true(is.numeric(bs))
  expect_true(bs >= 2L)  # at least {} and V
})

test_that("calculate_topology validates input", {
  expect_error(calculate_topology(1:5, -1), "non-negative")
  expect_error(calculate_topology("text", 1), "numeric")
})

test_that("analyze_topology_factors returns data.frame with plot attr", {
  set.seed(42)
  result <- analyze_topology_factors(rnorm(30), plot = TRUE)
  expect_s3_class(result, "data.frame")
  expect_true("factor" %in% names(result))
  expect_true("base_size" %in% names(result))
  expect_false(is.null(attr(result, "plot")))
})

test_that("analyze_topology_factors validates input", {
  expect_error(analyze_topology_factors("text"), "numeric")
  expect_error(analyze_topology_factors(c(1, NA)), "NA")
})

test_that("visualize_topology_thresholds returns data.frame", {
  set.seed(42)
  result <- visualize_topology_thresholds(rnorm(20), plot = TRUE)
  expect_s3_class(result, "data.frame")
  expect_true("method" %in% names(result))
  plots <- attr(result, "plots")
  expect_true(is.list(plots))
  expect_equal(length(plots), 3L)
})

# ---- Legacy connectivity function tests ----

test_that("is_topology_connected works on basic examples", {
  expect_true(is_topology_connected(list(c(1, 2, 3))))
  expect_false(is_topology_connected(list()))
  expect_true(is_topology_connected(list(c(1, 2), c(2, 3))))
})

test_that("is_topology_connected2 works on basic examples", {
  expect_true(is_topology_connected2(list(c(1, 2, 3))))
  expect_false(is_topology_connected2(list()))
})

test_that("is_topology_connected_manual checks coverage", {
  expect_true(is_topology_connected_manual(list(c(1, 2, 3), c(3, 4, 5))))
  expect_false(is_topology_connected_manual(list(c(1, 2), c(4, 5))))
})

# ---- 0.4.0: what the legacy functions compute, against independent referents ----

pairwise_intersections <- function(x, h) {
  nb <- lapply(seq_along(x), function(i) which(abs(x - x[i]) <= h))
  out <- list()
  for (i in seq_along(x)) for (j in i:length(x)) {
    s <- intersect(nb[[i]], nb[[j]])
    if (length(s)) out <- c(out, list(sort(s)))
  }
  unique(out)
}

test_that("the size columns of analyze_topology_factors describe the neighbourhood intersections", {
  set.seed(4104)
  x <- stats::rnorm(25)
  res <- analyze_topology_factors(x, factors = c(1, 3, 9), plot = FALSE)
  for (r in seq_len(nrow(res))) {
    inter <- pairwise_intersections(x, stats::IQR(x) / res$factor[r])
    expect_identical(res$max_set_size[r], max(lengths(inter)))
    expect_identical(res$min_set_size[r], min(lengths(inter)))
    expect_identical(res$base_size[r],
                     length(unique(c(list(integer(0), seq_along(x)), inter))))
    expect_gte(res$min_set_size[r], 1L)
  }
})

test_that("calculate_topology counts the empty set, the whole set and the intersections", {
  set.seed(4105)
  x <- stats::rnorm(15)
  inter <- pairwise_intersections(x, 0.4)
  expect_identical(calculate_topology(x, 0.4),
                   length(unique(c(list(integer(0), seq_along(x)), inter))))
  expect_error(calculate_topology(c(1, NA, 3), 0.4), "NA")
  expect_error(analyze_topology_factors(1:5, factors = c(1, -2)), "positive")
})

test_that("the documented counterexamples of the legacy connectivity checks hold", {
  connected_tau <- list(integer(0), 3L, c(1L, 3L), c(2L, 3L), 1:3)
  expect_true(is_topology_connected_exact(connected_tau, 3L)$connected)
  expect_false(is_topology_connected2(list(c(1L, 3L), c(2L, 3L))))
  discrete <- list(integer(0), 1L, 2L, c(1L, 2L))
  expect_false(is_topology_connected_exact(discrete, 2L)$connected)
  expect_true(is_topology_connected(discrete))
  expect_true(is_topology_connected2(discrete))
})
