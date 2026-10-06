#' Horizontal visibility graph of a time series
#'
#' @description
#' Builds the horizontal visibility graph of a series (Luque et al., 2009): two
#' observations are adjacent when every observation strictly between them is
#' strictly lower than both. With \code{directed = TRUE} it also returns the
#' arcs from the earlier to the later observation of each edge, which
#' [generate_bitopology()] turns into the forward and backward topologies.
#' The instants of observation enter only through their order, so
#' \code{times} is validated and recorded but does not change the graph.
#'
#' @param series Numeric vector without \code{NA}, \code{NaN} or infinite
#'   values. A series of length 0 or 1 gives a graph without edges. No default.
#' @param directed Logical scalar (default \code{FALSE}): whether to return the
#'   directed adjacency lists as well.
#' @param times \code{NULL} (default), meaning the instants
#'   \eqn{1, \ldots, n}, or a numeric, \code{Date} or date-time vector of
#'   finite, strictly increasing instants of the same length as \code{series};
#'   dates count in days and date-times in seconds.
#' @return A \code{list} with:
#'   \describe{
#'     \item{edges}{\code{data.frame} with integer columns \code{from} and
#'       \code{to}, \code{from < to}, one row per edge.}
#'     \item{n}{Integer: the number of observations.}
#'     \item{n_edges}{Integer: the number of edges.}
#'     \item{adjacency}{List of integer vectors: element \code{i} holds the
#'       neighbours of observation \code{i}, its open neighbourhood
#'       \eqn{N(i)}; [generate_topology()] adds \code{i} itself.}
#'     \item{times}{Numeric vector: the instants used.}
#'     \item{out_adjacency}{Only with \code{directed = TRUE}: list whose
#'       element \code{i} holds the later neighbours of \code{i}.}
#'     \item{in_adjacency}{Only with \code{directed = TRUE}: list whose
#'       element \code{i} holds the earlier neighbours of \code{i}.}
#'   }
#'
#' @details
#' The series is scanned from left to right with a stack whose values decrease
#' strictly from bottom to top, in \eqn{O(n)} time. When observation \eqn{i}
#' arrives, every stacked index \eqn{j} with \eqn{y_j \le y_i} is popped and
#' joined to \eqn{i}: nothing between \eqn{j} and \eqn{i} reaches \eqn{y_j},
#' or it would have popped \eqn{j} earlier. The index left on top, if any, is
#' higher than \eqn{y_i} and is joined to \eqn{i} unless a popped value equals
#' \eqn{y_i}, which blocks at height \eqn{\min = y_i}. Every other earlier index
#' is blocked: one deeper in the stack by that top, which is higher than
#' \eqn{y_i}; one popped before by the index that popped it, which lies in
#' between and is at least as high. Each edge is thus reported exactly once,
#' always with \code{from < to}.
#'
#' @section Methodological notes:
#' The criterion uses comparisons of values and nothing else, so the graph is
#' computed exactly on the numbers received, is invariant under any strictly
#' increasing transformation of the values, and does not depend on the
#' spacing of the instants. Decimal data are converted to binary doubles
#' before they reach the function, and the graph is that of the doubles. A
#' conversion that rounds each decimal to the nearest double keeps the order
#' of the decimals, and keeps distinct any two with at most fifteen
#' significant digits in the normal range (magnitudes from about
#' \eqn{2.2 \times 10^{-308}} to \eqn{1.8 \times 10^{308}}), because their
#' relative gap, at least \eqn{10^{-15}}, exceeds the relative spacing of
#' doubles, at most \eqn{2^{-52}}; outside that range distinct decimals can
#' collapse:
#' \eqn{2 \times 10^{-324}}, \eqn{10^{-324}} and \eqn{2 \times 10^{-324}} all
#' become zero, and the edge between the first and the third disappears.
#'
#' The function stops with an error when the floating-point environment does
#' not provide IEEE 754 double arithmetic (see [natural_visibility_exactness]
#' for the test). Under the denormals-are-zero mode, which some libraries
#' switch on, the comparisons would treat every subnormal value as zero; the
#' flush-to-zero mode, which changes results of arithmetic and not these
#' comparisons, is refused as well, because the same test guards the natural
#' criterion, which depends on it.
#'
#' @section Dependencies:
#' The compiled engine (C++ through 'Rcpp').
#'
#' @references
#'   Luque, B., Lacasa, L., Ballesteros, F., & Luque, J. (2009). Horizontal
#'   visibility graphs: Exact results for random time series. Physical Review
#'   E, 80(4), Article 046103. https://doi.org/10.1103/PhysRevE.80.046103
#'
#' @seealso [natural_visibility_graph()], [generate_bitopology()],
#'   [time_reverse()].
#'
#' @examples
#' series <- c(3, 1, 4, 1, 5, 9, 2, 6)
#' g <- horizontal_visibility_graph(series)
#' g$n_edges
#' head(g$edges)
#'
#' gd <- horizontal_visibility_graph(series, directed = TRUE)
#' gd$out_adjacency[[1]]
#' gd$in_adjacency[[5]]
#' @export
horizontal_visibility_graph <- function(series, directed = FALSE, times = NULL) {
  validate_series(series)
  n <- length(series)
  times <- validate_times(times, n)
  raw <- if (n < 2L) list(from = integer(0), to = integer(0)) else hvg_cpp(as.double(series))
  visibility_result(raw, n, times, directed)
}


#' Natural visibility graph of a time series
#'
#' @description
#' Builds the natural visibility graph of a series (Lacasa et al., 2008): two
#' observations are adjacent when every observation strictly between them lies
#' strictly below the straight segment that joins them, drawn at the instants
#' of observation, equally spaced unless \code{times} says otherwise. The
#' criterion is decided exactly on the numbers received, under the hypotheses
#' on the floating-point arithmetic stated in [natural_visibility_exactness].
#' With
#' \code{directed = TRUE} it also returns the arcs from the earlier to the
#' later observation of each edge, which [generate_bitopology()] turns into the
#' forward and backward topologies.
#'
#' @inheritParams horizontal_visibility_graph
#' @return A \code{list} with the same elements as
#'   [horizontal_visibility_graph()].
#'
#' @details
#' Observations \eqn{a < b} are adjacent when, for every \eqn{k} between them,
#' \deqn{D(a, b, k) = y_a (t_b - t_k) + y_b (t_k - t_a) - y_k (t_b - t_a) > 0.}
#' The graph is built by divide and conquer on the maximum. Each level of the
#' recursion makes at most one decision per point of its intervals, so the
#' number of decisions is \eqn{O(n^2)} in the worst case (monotone data, whose
#' maximum is always at an end), \eqn{O(n \log n)} when the maximum of every
#' interval splits it in parts of bounded proportion, and \eqn{O(n \log n)} in
#' expectation when the values are exchangeable and without ties, because the
#' position of each maximum is then uniform, as the pivot of randomized
#' quicksort. Three facts make the construction correct for any strictly
#' increasing instants. (i) If \eqn{p} maximizes \eqn{y} on an interval and
#' \eqn{a < p < b}, the chord from \eqn{a} to \eqn{b} takes at \eqn{t_p} a
#' convex combination of \eqn{y_a} and \eqn{y_b}, at most \eqn{y_p}, so \eqn{p}
#' blocks every pair across it. (ii) For \eqn{j < p}, every \eqn{k} between
#' them is below the chord from \eqn{j} to \eqn{p} if and only if the slope
#' from \eqn{j} to \eqn{p} is smaller than the slope from each such \eqn{k} to
#' \eqn{p}; scanning leftwards from \eqn{p}, the smallest of those slopes
#' belongs to the last index found visible, \eqn{m}, so \eqn{j} is visible if
#' and only if \eqn{D(j, p, m) > 0}, and symmetrically on the right. (iii) The
#' edges inside each half depend only on the points of that half.
#'
#' Each decision is therefore the sign of one quantity \eqn{D}, and it is
#' computed exactly under those hypotheses. When the six inputs of a decision are at most
#' \eqn{2^{510}} in magnitude, so that no intermediate value can overflow, the
#' engine evaluates in floating
#' point the orientation determinant \eqn{O = -D}, as the difference of two
#' rounded products \eqn{l} and \eqn{r}, and accepts the sign of the result
#' when its magnitude exceeds \eqn{8uS}, where \eqn{u = 2^{-53}} and \eqn{S}
#' is the computed value of \eqn{|l| + |r|}: above that threshold the sign is
#' proved correct. Otherwise, near a tie, or when \eqn{S} does not exceed the
#' \eqn{2^{-930}}, about \eqn{1.1 \times 10^{-280}}, \eqn{D} is evaluated exactly in integer
#' arithmetic on the binary mantissas and exponents of the six products it
#' contains. This is the adaptive scheme of Shewchuk (1997). The proofs of
#' both stages, and the hypotheses on the arithmetic under which they hold
#' (every IEEE 754 rounding mode; evaluation in binary64, in a format with at
#' least 64 significand bits, or in the x87 format at 53 bits), with what is
#' checked of them, are written in [natural_visibility_exactness]; every
#' intermediate value of the floating-point evaluation is stored, so that the
#' compiler cannot fuse a product into the subtraction or reorder the
#' operations. A floating-point decision takes a fixed number of operations;
#' an exact one allocates and compares two accumulators of \eqn{L} words of 32
#' bits, where \eqn{L} grows with the spread of the binary exponents of the six
#' products and is at most 137, so its number of word operations is
#' proportional to \eqn{L}.
#' The function stops with an error when the floating-point environment does
#' not provide IEEE 754 double arithmetic (subnormal numbers flushed to zero,
#' or operations with fewer than 53 significant bits), a condition that some
#' libraries create and under which the proofs do not apply.
#'
#' @section Methodological notes:
#' \emph{What exact means here.} Under the hypotheses of
#' [natural_visibility_exactness], the graph is the natural visibility graph of
#' the doubles received, independent of the order of floating-point
#' operations. Consequences, under the same hypotheses: reversing the values
#' and reflecting the instants (see [time_reverse()]) exchanges the forward and
#' backward topologies exactly, for every input; and the graph does not change
#' when the values are replaced by \eqn{\alpha y + \beta} and the instants by
#' \eqn{\gamma t + \delta}, with \eqn{\alpha, \gamma > 0}, provided that the new
#' values and instants are finite doubles equal to these exact images (no
#' rounding), so that the instants stay strictly increasing. Proof: each
#' \eqn{D} computed for the images equals \eqn{\alpha \gamma D}, because the
#' coefficients of \eqn{y_a}, \eqn{y_b} and \eqn{y_k} in \eqn{D} add up to zero,
#' so \eqn{\beta} cancels, and every difference of instants is multiplied by
#' \eqn{\gamma}, so \eqn{\delta} cancels.
#' Versions up to 0.3.0 accepted only equally spaced instants and compared
#' the rounded slopes \eqn{\mathrm{fl}(\mathrm{fl}(y_p - y_j) / (p - j))}. When
#' no difference of values overflows and no slope is non-zero and below
#' \eqn{2^{-1022}} in magnitude, each rounded slope is within relative error
#' \eqn{2u + u^2} of the exact one (\eqn{u = 2^{-53}}, rounding to nearest), so
#' the two engines can disagree only on two slopes whose difference is at most
#' about \eqn{2u} times the sum of their magnitudes, that is, when a point lies
#' on a chord or very near it. Outside that range of magnitudes, with values
#' near the largest doubles or slopes in the subnormal range, the rounded
#' slopes of 0.3.0 could also err far from any tie.
#'
#' \emph{Decimal data.} A decimal such as 0.1 has no exact binary
#' representation, so points that are collinear in decimal can be slightly
#' below or above the line in binary, and the exact criterion then follows the
#' binary numbers. Suppose the values were recorded with \eqn{d} decimal
#' places, \eqn{q_i = m_i / 10^d} with a non-negative integer \eqn{d \le 22}
#' and integers \eqn{|m_i| \le 2^{50}}, and that each was stored as the double
#' \eqn{y_i} nearest to it. Then \code{m <- round(10^d * y)} returns the
#' integers \eqn{m_i}, and the graph of the recorded decimals is the graph of
#' \code{m}, computed exactly; the same applies to the instants. The check
#' \code{all(abs(m) <= 2^50) && all(m / 10^d == y)} passes in that case, fails
#' whenever some \eqn{y_i} is not the double nearest to a decimal with
#' \eqn{d} places, and, when it passes, \eqn{m_i / 10^d} is the only decimal
#' with \eqn{d} places whose nearest double is \eqn{y_i}. Proof, in the
#' rounding to nearest in which R computes: \eqn{10^d} is an exact double for
#' \eqn{d \le 22}, because \eqn{5^{22} < 2^{53}}. Since \eqn{y_i} is nearest
#' to \eqn{q_i}, \eqn{y_i = q_i (1 + \delta)} with \eqn{|\delta| < u = 2^{-53}};
#' the computed product \code{10^d * y[i]} equals
#' \eqn{m_i (1 + \delta)(1 + \delta')} with \eqn{|\delta'| < u}, so it differs
#' from \eqn{m_i} by less than \eqn{2^{50} (2u + u^2)}, which is below
#' \eqn{1/2}, and \code{round()} returns \eqn{m_i}; \code{m[i] / 10^d}, a single
#' correctly rounded division, returns the double nearest to \eqn{q_i}, that
#' is \eqn{y_i}. When the check passes, \eqn{y_i} is the double nearest to
#' \eqn{m_i / 10^d}, and no other decimal with \eqn{d} places shares it: two
#' such decimals differ by at least \eqn{10^{-d}}, while two numbers rounded to
#' the same \eqn{y_i} differ by at most \eqn{2u |y_i|}, which is below
#' \eqn{10^{-d} / 2}. Finally, multiplying every value by \eqn{10^d > 0}
#' multiplies every \eqn{D} by \eqn{10^d}, which keeps its sign, and the
#' integers \eqn{m_i} are exact doubles. Without the bound the recipe can
#' fail: the decimals
#' \eqn{2^{53}}, \eqn{2^{53} + 1} and \eqn{2^{53} + 2} are stored as
#' \eqn{2^{53}}, \eqn{2^{53}} and \eqn{2^{53} + 2}, and rounding cannot recover
#' the middle one.
#'
#' \emph{What the graph keeps of the series.} Lacasa et al. (2008) describe the
#' graph as inheriting properties of the series: in the words of their
#' abstract, periodic series convert into regular graphs, random series into
#' random graphs, and fractal series into scale-free networks. These are
#' statements of the abstract, which does not say for which series they hold;
#' they do not hold for every finite series: the periodic series
#' 0, 1, 0, 1, 0 at the equally spaced instants 1, ..., 5 has, under the
#' natural criterion, the edges 1-2, 2-3, 2-4, 3-4 and 4-5, and the degrees
#' 1, 3, 2, 3, 1.
#'
#' @section Dependencies:
#' The compiled engine (C++ through 'Rcpp').
#'
#' @references
#'   Lacasa, L., Luque, B., Ballesteros, F., Luque, J., & Nuño, J. C. (2008).
#'   From time series to complex networks: The visibility graph. Proceedings
#'   of the National Academy of Sciences, 105(13), 4972-4975.
#'   https://doi.org/10.1073/pnas.0709247105
#'
#'   Shewchuk, J. R. (1997). Adaptive precision floating-point arithmetic and
#'   fast robust geometric predicates. Discrete & Computational Geometry,
#'   18(3), 305-363. https://doi.org/10.1007/PL00009321
#'
#' @seealso [horizontal_visibility_graph()], [generate_bitopology()],
#'   [time_reverse()].
#'
#' @examples
#' series <- c(3, 1, 4, 1, 5, 9, 2, 6)
#' g <- natural_visibility_graph(series)
#' g$n_edges
#' head(g$edges)
#'
#' natural_visibility_graph(c(0, 2, 4), times = c(0, 1, 4))$edges
#'
#' y <- round(0.1 * (1:40) + 0.3, 1)
#' natural_visibility_graph(y)$n_edges
#' m <- round(10 * y)
#' all(abs(m) <= 2^50) && all(m / 10 == y)
#' natural_visibility_graph(m)$n_edges
#' @export
natural_visibility_graph <- function(series, directed = FALSE, times = NULL) {
  validate_series(series)
  n <- length(series)
  times <- validate_times(times, n)
  raw <- if (n < 2L) list(from = integer(0), to = integer(0)) else nvg_cpp(as.double(series), times)
  visibility_result(raw, n, times, directed)
}


#' Reverse a time series together with its instants
#'
#' @description
#' Returns the series read backwards together with the reflection of its
#' instants, the operation under which the forward and backward topologies of
#' a visibility graph are exchanged (theorem T1, stated and proved in
#' [generate_bitopology()]). Use it instead of reversing the values alone,
#' which is a different operation as soon as the instants are not equally
#' spaced.
#'
#' @inheritParams horizontal_visibility_graph
#' @return A \code{list} with \code{series}, the reversed values
#'   \eqn{\tilde{y}_i = y_{n+1-i}}, and \code{times}: \code{NULL} when
#'   \code{times} is \code{NULL}, and otherwise the numeric vector
#'   \eqn{\tilde{t}_i = -t_{n+1-i}}, again strictly increasing.
#'
#' @details
#' Reflecting the instants keeps every gap between consecutive observations
#' and reverses their order. Any reflection \eqn{c - t_{n+1-i}} gives the same
#' visibility graphs in exact arithmetic, because both criteria are invariant
#' under translations of the time axis, and \eqn{c = t_1 + t_n} keeps the
#' observation window. The function takes \eqn{c = 0} because negation is exact
#' in floating point, whereas \eqn{c - t} can round. With \code{times = NULL}
#' the reflection of \eqn{1, \ldots, n} is an exact translation of the same
#' equally spaced grid, so \code{NULL} is returned.
#'
#' @section Methodological notes:
#' Reversing the values over the same instants coincides with this operation
#' only when the gaps between consecutive instants read the same in both
#' directions. With the instants 0, 1, 4 and the values 0, 2, 4, the reversed
#' values 4, 2, 0 over the same instants make the first and the third
#' observation visible to each other, while in the original series they are
#' not; with the reflected instants -4, -1, 0 the exchange of the two
#' topologies holds.
#'
#' @section Dependencies:
#' Base R only.
#'
#' @seealso [generate_bitopology()] for the theorem and its consequences,
#'   [natural_visibility_graph()] and [horizontal_visibility_graph()].
#'
#' @examples
#' rv <- time_reverse(c(0, 2, 4), times = c(0, 1, 4))
#' rv
#' natural_visibility_graph(rv$series, times = rv$times)$edges
#' natural_visibility_graph(c(0, 2, 4), times = c(0, 1, 4))$edges
#' natural_visibility_graph(c(4, 2, 0), times = c(0, 1, 4))$edges
#' @export
time_reverse <- function(series, times = NULL) {
  validate_series(series)
  n <- length(series)
  list(series = rev(series),
       times = if (is.null(times)) NULL else -rev(validate_times(times, n)))
}


# ----------------------------------------------------------------------
# Internal helpers
# ----------------------------------------------------------------------

# A numeric series without missing or infinite values (any length).
validate_series <- function(series) {
  if (!is.numeric(series)) {
    stop("'series' must be a numeric vector.", call. = FALSE)
  }
  if (anyNA(series)) {
    stop("'series' must not contain NA values.", call. = FALSE)
  }
  if (!all(is.finite(series))) {
    stop("'series' must not contain Inf or NaN values.", call. = FALSE)
  }
  invisible(TRUE)
}

# Instants of observation: NULL means 1, ..., n. Numeric, Date and date-time
# vectors are accepted and converted to numbers (days for Date, seconds for
# date-times); the converted numbers must be finite and strictly increasing.
validate_times <- function(times, n) {
  if (is.null(times)) return(as.double(seq_len(n)))
  if (inherits(times, "POSIXlt")) times <- as.POSIXct(times)
  if (inherits(times, "Date") || inherits(times, "POSIXct")) {
    times <- as.numeric(times)
  }
  if (!is.numeric(times)) {
    stop("'times' must be NULL, a numeric vector, a Date vector or a ",
         "date-time vector.", call. = FALSE)
  }
  times <- as.double(times)
  if (length(times) != n) {
    stop("'times' must have the same length as 'series' (", n, ").",
         call. = FALSE)
  }
  if (anyNA(times) || !all(is.finite(times))) {
    stop("'times' must contain finite values only.", call. = FALSE)
  }
  if (n >= 2L) {
    bad <- which(times[-1L] <= times[-n])
    if (length(bad)) {
      stop("'times' must be strictly increasing; times[", bad[1L] + 1L,
           "] is not greater than times[", bad[1L], "].", call. = FALSE)
    }
  }
  times
}

# Assemble the result of a visibility-graph constructor from the raw edge
# list returned by the engine (every edge with from < to).
visibility_result <- function(raw, n, times, directed) {
  from <- as.integer(raw$from)
  to <- as.integer(raw$to)
  adjacency <- lapply(seq_len(n), function(i) integer(0))
  for (k in seq_along(from)) {
    adjacency[[from[k]]] <- c(adjacency[[from[k]]], to[k])
    adjacency[[to[k]]] <- c(adjacency[[to[k]]], from[k])
  }
  result <- list(
    edges = data.frame(from = from, to = to),
    n = n,
    n_edges = length(from),
    adjacency = adjacency,
    times = times
  )
  if (directed) {
    out_adj <- lapply(seq_len(n), function(i) integer(0))
    in_adj <- lapply(seq_len(n), function(i) integer(0))
    for (k in seq_along(from)) {
      out_adj[[from[k]]] <- c(out_adj[[from[k]]], to[k])
      in_adj[[to[k]]] <- c(in_adj[[to[k]]], from[k])
    }
    result$out_adjacency <- out_adj
    result$in_adjacency <- in_adj
  }
  result
}
