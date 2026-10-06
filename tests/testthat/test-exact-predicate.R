# ======================================================================
# The exact chord predicate of the natural visibility criterion (0.4.0)
#
# Referents EXTERNAL to the engine:
#   (i)  fixtures/exact_predicate_cases.csv: 1,800 sextuples in nine
#        regimes (normal, collinear dyadic, decimal collinear, decimal
#        instants and values, extreme exponents, subnormal, near-collinear
#        within a few ulps, large integers, beyond the bound of the
#        floating-point filter), written as hexadecimal doubles, with the
#        sign of D computed in exact rational arithmetic (Python fractions)
#        outside the package by dev/make_exact_predicate_fixture.py, which
#        regenerates the file deterministically;
#   (ii) an R reference on small integers and dyadic rationals, where double
#        arithmetic is exact, so the plain formula is the truth;
#   (iii) the literal O(n^3) definition of the natural visibility graph in R,
#        on integer and dyadic series, where the same holds.
# ======================================================================

chord_sign <- function(ta, xa, tb, xb, tk, xk) {
  topologyR:::chord_below_sign_cpp(ta, xa, tb, xb, tk, xk)
}

test_that("the predicate reproduces the exact signs of the external fixture", {
  cases <- utils::read.csv(test_path("fixtures", "exact_predicate_cases.csv"),
                           colClasses = "character")
  expect_identical(nrow(cases), 1800L)
  num <- function(v) as.numeric(v)
  got <- chord_sign(num(cases$ta), num(cases$xa), num(cases$tb),
                    num(cases$xb), num(cases$tk), num(cases$xk))
  expected <- as.integer(cases$sign)
  for (kind in unique(cases$kind)) {
    sel <- cases$kind == kind
    expect_identical(got[sel], expected[sel], info = kind)
  }
  expect_true(any(expected == 0L))
})

test_that("the predicate equals the plain formula where double arithmetic is exact", {
  set.seed(4101)
  m <- 20000
  ta <- sample(-1000:1000, m, TRUE); tk <- sample(-1000:1000, m, TRUE)
  tb <- sample(-1000:1000, m, TRUE)
  xa <- sample(-1000:1000, m, TRUE) / 8; xb <- sample(-1000:1000, m, TRUE) / 8
  xk <- sample(-1000:1000, m, TRUE) / 8
  # Exactness bound: every value is a multiple of 1/8, every difference of
  # instants an integer below 2^11 in magnitude and every value below 2^7, so
  # each product and each partial sum is a multiple of 1/8 below 2^20 in
  # magnitude and needs at most 23 significant bits, fewer than the 53 of a
  # double.
  expect_true(all(abs(c(ta, tb, tk)) <= 1000) && all(abs(c(xa, xb, xk)) <= 125))
  expect_lt(3 * 2000 * 125, 2^20)
  truth <- sign(xa * (tb - tk) + xb * (tk - ta) - xk * (tb - ta))
  expect_identical(chord_sign(ta, xa, tb, xb, tk, xk), as.integer(truth))
})

test_that("the predicate is invariant under the reflection of the instants", {
  set.seed(4102)
  m <- 5000
  ta <- stats::rnorm(m); tb <- stats::rnorm(m); tk <- stats::rnorm(m)
  xa <- stats::rnorm(m); xb <- stats::rnorm(m); xk <- stats::rnorm(m)
  expect_identical(chord_sign(-tb, xb, -ta, xa, -tk, xk),
                   chord_sign(ta, xa, tb, xb, tk, xk))
  expect_identical(chord_sign(4 * ta, 8 * xa, 4 * tb, 8 * xb, 4 * tk, 8 * xk),
                   chord_sign(ta, xa, tb, xb, tk, xk))
})

nvg_literal <- function(y, t = seq_along(y)) {
  n <- length(y)
  e <- NULL
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      k <- if (j > i + 1L) (i + 1L):(j - 1L) else integer(0)
      if (all(y[k] * (t[j] - t[i]) < y[i] * (t[j] - t[k]) + y[j] * (t[k] - t[i]))) {
        e <- rbind(e, c(i, j))
      }
    }
  }
  e
}

edge_keys <- function(e) sort(paste(e[, 1], e[, 2], sep = "-"))

test_that("the natural visibility graph equals the literal definition on exact data", {
  set.seed(4103)
  for (r in 1:150) {
    n <- sample(3:25, 1)
    y <- switch(r %% 3 + 1,
                as.numeric(sample(0:4, n, TRUE)),
                sample(-20:20, n, TRUE) / 4,
                as.numeric(seq_len(n) * sample(-2:2, 1) + sample(0:1, n, TRUE)))
    t <- if (r %% 2 == 0) cumsum(sample(1:4, n, TRUE)) / 2 else seq_len(n)
    # Exactness bound of the literal reference: values are multiples of 1/4
    # below 2^6, instants multiples of 1/2 below 2^7, so every product is a
    # multiple of 1/8 below 2^13 and every comparison is exact in double.
    expect_true(max(abs(y)) < 2^6 && max(abs(t)) < 2^7)
    g <- natural_visibility_graph(y, times = if (r %% 2 == 0) t else NULL)
    expect_identical(edge_keys(cbind(g$edges$from, g$edges$to)),
                     edge_keys(nvg_literal(y, t)),
                     info = sprintf("case %d (n = %d)", r, n))
  }
})

test_that("a straight line has only the consecutive edges, exactly", {
  for (n in c(3L, 10L, 60L)) {
    g <- natural_visibility_graph(3 * seq_len(n) - 7)
    expect_identical(g$n_edges, n - 1L)
    g <- natural_visibility_graph(round(10 * (0.1 * seq_len(n) + 0.3)))
    expect_identical(g$n_edges, n - 1L)
  }
})

test_that("the decimal recipe recovers the integers of the recorded decimals", {
  # The proof in ?natural_visibility_graph needs 10^d to be an exact double
  # for d <= 22; sprintf("%.0f") prints the exact value of a double.
  for (d in 0:22) {
    expect_identical(sprintf("%.0f", 10^d), paste0("1", strrep("0", d)))
  }
  # y is built as the double nearest to m / 10^d (one correctly rounded
  # division), the hypothesis of the recipe, including integers next to the
  # bound 2^50.
  set.seed(4104)
  for (d in c(0:6, 10, 15, 22)) {
    m_true <- c(sample(-10^6:10^6, 300), 2^50 - sample(0:1000, 50),
                -(2^50 - sample(0:1000, 50)))
    y <- m_true / 10^d
    m <- round(10^d * y)
    expect_identical(m, m_true, info = sprintf("d = %d", d))
    expect_true(all(abs(m) <= 2^50) && all(m / 10^d == y), info = sprintf("d = %d", d))
  }
})

test_that("the ordered sextuple that overflows without the bound has the exact sign", {
  # t_a = -M < t_k = 2^1023 < t_b = M, with M the largest double, and
  # x_a = 5/4, x_k = 0, x_b = -1/2. Exactly,
  # D = (5/4)(M - 2^1023) - (1/2)(2^1023 + M) = -2^1021 - 3 * 2^969 < 0.
  # The inputs exceed 2^510, so the exact stage decides; a filter without
  # that bound returns +1 under rounding toward zero (see
  # ?natural_visibility_exactness).
  M <- .Machine$double.xmax
  expect_identical(M, 2^1023 * (2 - 2^-52))
  expect_identical(chord_sign(-M, 1.25, M, -0.5, 2^1023, 0), -1L)
})
