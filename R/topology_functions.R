#' Connectivity of the co-occurrence graph of a family of sets
#'
#' @description
#' Builds the co-occurrence graph of a family of sets, in which two elements
#' are adjacent when some member of the family contains both, and reports
#' whether that graph is connected. It is retained from version 0.1.0 as a
#' quick necessary-condition screen; the exact test of topological
#' connectedness is [generate_topology()] for a topology given by a subbase and
#' [is_topology_connected_exact()] for a fully enumerated topology.
#'
#' @param topology A list of integer vectors (positive labels), read as a
#'   family of subsets of the labels that appear in it. A full topology, a base
#'   or a subbase are all admissible. No default.
#' @return A \code{logical} scalar: \code{TRUE} if the co-occurrence graph on
#'   the labels that appear in the family is connected, \code{FALSE} if it is
#'   not or if the family is empty.
#'
#' @details
#' The ground set is taken to be the set of labels that appear in some member,
#' and the graph is traversed from the first of them. Memory grows with the
#' square of the largest label, because the graph is stored as a dense
#' adjacency matrix.
#'
#' @section Methodological notes:
#' For any family of open sets that covers the ground set, a topologically
#' connected space has a connected co-occurrence graph: if the graph split into
#' two non-empty parts, every member would lie inside one part, each part would
#' be a union of members and hence open, and the two parts would disconnect the
#' space. A disconnected graph therefore proves that the space is disconnected.
#' The converse fails: the discrete topology on \{1, 2\}, given by its open sets
#' \{\}, \{1\}, \{2\} and \{1, 2\}, is disconnected while its co-occurrence
#' graph is connected, and any list that contains the whole ground set makes
#' the graph complete. A \code{TRUE} answer is thus inconclusive.
#'
#' @section Dependencies:
#' Base R only.
#'
#' @seealso [is_topology_connected_exact()] for the exact test on an enumerated
#'   topology, [generate_topology()] for exact components from a subbase,
#'   [is_topology_connected2()] and [is_topology_connected_manual()] for the
#'   other legacy diagnostics.
#'
#' @examples
#' is_topology_connected(list(c(1, 2, 3), c(3, 4, 5)))
#' is_topology_connected(list(c(1, 2), c(3, 4)))
#' is_topology_connected(list(integer(0), 1L, 2L, c(1L, 2L)))
#' @export
is_topology_connected <- function(topology) {
  if (!length(topology)) return(FALSE)

  elements <- unique(unlist(topology))
  if (length(elements) <= 1L) return(TRUE)

  n <- max(elements)
  edges <- matrix(0L, nrow = n, ncol = n)
  for (set in topology) {
    for (i in set) {
      for (j in set) {
        if (i != j && i <= n && j <= n) {
          edges[i, j] <- 1L
          edges[j, i] <- 1L
        }
      }
    }
  }

  visited <- rep(FALSE, n)
  stack <- elements[1L]
  while (length(stack) > 0L) {
    v <- stack[1L]
    stack <- stack[-1L]
    if (!visited[v]) {
      visited[v] <- TRUE
      neighbors <- which(edges[v, ] == 1L)
      stack <- c(stack, neighbors[!visited[neighbors]])
    }
  }

  all(visited[elements])
}


#' Label-order reachability within a family of sets
#'
#' @description
#' Builds a directed graph with an arc between consecutive labels, in
#' increasing order, inside each member of a family of sets, and reports
#' whether every label that appears in the family is reachable from the
#' smallest one. It is retained from version 0.1.0 for backward compatibility
#' and is not a test of topological connectedness; use [generate_topology()] or
#' [is_topology_connected_exact()] for that.
#'
#' @param topology A list of integer vectors (positive labels), read as a
#'   family of subsets. No default.
#' @return A \code{logical} scalar: \code{TRUE} if every label that appears in
#'   the family is reachable from the smallest label along the arcs described
#'   above, \code{FALSE} otherwise or if the family is empty.
#'
#' @details
#' Memory grows with the square of the largest label, because the arcs are
#' stored as a dense matrix.
#'
#' @section Methodological notes:
#' The answer depends on how the elements are labelled, so it is not a
#' topological invariant, and it is neither a necessary nor a sufficient
#' condition for connectedness. The base \{1, 3\}, \{2, 3\} generates the
#' connected topology whose open sets are \{\}, \{3\}, \{1, 3\}, \{2, 3\} and
#' \{1, 2, 3\} (no proper non-empty open set has an open complement), yet label
#' 2 is not reachable from label 1 and the function returns \code{FALSE}. The
#' discrete topology on \{1, 2\}, listed as \{\}, \{1\}, \{2\}, \{1, 2\}, is
#' disconnected, yet the function returns \code{TRUE}.
#'
#' @section Dependencies:
#' Base R only.
#'
#' @seealso [generate_topology()] and [is_topology_connected_exact()] for exact
#'   connectedness, [is_topology_connected()] for a necessary condition.
#'
#' @examples
#' is_topology_connected2(list(c(1, 2, 3), c(3, 4, 5)))
#' is_topology_connected2(list(c(1L, 3L), c(2L, 3L)))
#' @export
is_topology_connected2 <- function(topology) {
  if (!length(topology)) return(FALSE)

  elements <- unique(unlist(topology))
  if (length(elements) <= 1L) return(TRUE)

  n <- max(elements)
  edges <- matrix(0L, nrow = n, ncol = n)

  for (set in topology) {
    if (length(set) > 1L) {
      values <- sort(set)
      for (i in seq_len(length(values) - 1L)) {
        edges[values[i], values[i + 1L]] <- 1L
      }
    }
  }

  visited <- rep(FALSE, n)
  start <- min(elements)
  stack <- start
  while (length(stack) > 0L) {
    current <- stack[1L]
    stack <- stack[-1L]
    if (!visited[current]) {
      visited[current] <- TRUE
      neighbors <- which(edges[current, ] == 1L)
      stack <- c(stack, neighbors[!visited[neighbors]])
    }
  }

  all(visited[elements])
}


#' Check whether a family of sets covers the labels 1 to its maximum
#'
#' @description
#' Checks whether every integer from 1 to the largest label of a family of
#' sets appears in at least one member, that is, whether the family covers the
#' ground set \{1, ..., max\}. Despite its name, retained from version 0.1.0,
#' it does not test connectedness.
#'
#' @param topology A list of integer vectors (positive labels). No default.
#' @return A \code{logical} scalar: \code{TRUE} if every integer from 1 to the
#'   largest label appears in some member, \code{FALSE} otherwise or if the
#'   family is empty.
#'
#' @section Methodological notes:
#' Coverage is necessary for a family to be a base or a subbase of a topology
#' on \{1, ..., max\}, but it says nothing about connectedness: the discrete
#' topology covers its ground set and is maximally disconnected.
#'
#' @section Dependencies:
#' Base R only.
#'
#' @seealso [is_topology_connected_exact()] and [generate_topology()] for
#'   connectedness.
#'
#' @examples
#' is_topology_connected_manual(list(c(1, 2, 3), c(3, 4, 5)))
#' is_topology_connected_manual(list(c(1, 2), c(4, 5)))
#' @export
is_topology_connected_manual <- function(topology) {
  all_elements <- unique(unlist(topology))
  if (length(all_elements) == 0L) return(FALSE)
  n <- max(all_elements)
  all(seq_len(n) %in% all_elements)
}


#' Size of the threshold-neighbourhood base for several IQR factors
#'
#' @description
#' For each factor \eqn{f}, builds the threshold neighbourhoods of the data
#' with radius \eqn{h = \mathrm{IQR}/f} and reports the size of the base they
#' generate and the sizes of its largest and smallest members. It shows how the
#' choice of the factor changes the granularity of a metric-threshold topology;
#' it does not select a factor, because no criterion of optimality is defined.
#'
#' @param data Numeric vector with at least two elements and no missing
#'   values. No default.
#' @param factors Positive numeric vector of divisors of the interquartile
#'   range. Default \code{NULL}, which means \code{c(1, 2, 4, 8, 16)}.
#' @param plot Logical scalar; if \code{TRUE} (default) the result carries a
#'   \code{ggplot} object in its \code{"plot"} attribute.
#' @return A \code{data.frame} with one row per factor and the columns
#'   \code{factor} (numeric), \code{threshold} (numeric, \eqn{\mathrm{IQR}/f}),
#'   \code{base_size} (integer, number of sets in the family described in
#'   Details, the empty set and the whole index set included),
#'   \code{max_set_size} and \code{min_set_size} (integers, largest and
#'   smallest cardinality among the non-empty intersections of two
#'   neighbourhoods). If \code{plot = TRUE}, the data frame carries an
#'   attribute \code{"plot"} with a \code{ggplot} object.
#'
#' @details
#' The neighbourhood of observation \eqn{i} is the set of indices \eqn{j} with
#' \eqn{|x_j - x_i| \le h}. The family counted in \code{base_size} consists of
#' the empty set, the whole index set and every non-empty intersection of two
#' neighbourhoods, deduplicated. Each neighbourhood collects the observations
#' whose values lie in a closed interval, and an intersection of several such
#' sets is the intersection of the two whose intervals have the largest left
#' and the smallest right endpoint; the pairwise intersections are therefore
#' all the non-empty finite intersections, and together with the whole index
#' set they form the base generated by the neighbourhoods.
#'
#' @section Methodological notes:
#' These are metric-threshold neighbourhoods on the values, not the
#' graph-induced neighbourhoods of the visibility-graph pipeline; the radius
#' is a tuning choice of the user. Up to version 0.3.0 the two size columns
#' were taken over the whole family, so \code{min_set_size} was always 0 (the
#' empty set) and \code{max_set_size} always \eqn{n} (the whole index set);
#' they are now taken over the neighbourhood intersections, which is where the
#' factor acts.
#'
#' @section Dependencies:
#' \code{stats::IQR()} for the radius; 'ggplot2' (Imports) for the optional
#' plot.
#'
#' @seealso [calculate_topology()] for a single threshold,
#'   [calculate_thresholds()] for data-driven radii, [generate_topology()] for
#'   graph-induced topologies.
#'
#' @examples
#' set.seed(1)
#' results <- analyze_topology_factors(rnorm(50), plot = FALSE)
#' results
#' @export
analyze_topology_factors <- function(data, factors = NULL, plot = TRUE) {
  if (!is.numeric(data) || length(data) < 2L) {
    stop("'data' must be a numeric vector with at least 2 elements.",
         call. = FALSE)
  }
  if (anyNA(data)) {
    stop("'data' must not contain NA values.", call. = FALSE)
  }
  if (is.null(factors)) {
    factors <- c(1, 2, 4, 8, 16)
  }
  if (!is.numeric(factors) || !length(factors) || anyNA(factors) ||
      any(factors <= 0)) {
    stop("'factors' must be a numeric vector of positive numbers.",
         call. = FALSE)
  }

  results <- lapply(factors, function(f) {
    threshold <- stats::IQR(data) / f
    inter <- threshold_intersections(data, threshold)
    sizes <- lengths(inter)
    list(
      factor = f,
      threshold = threshold,
      base_size = length(unique(c(list(integer(0), seq_along(data)), inter))),
      max_set_size = max(sizes),
      min_set_size = min(sizes)
    )
  })

  results_df <- do.call(rbind, lapply(results, data.frame))

  if (plot) {
    p <- ggplot2::ggplot(results_df, ggplot2::aes(x = factor)) +
      ggplot2::geom_line(ggplot2::aes(y = base_size, color = "Base Size")) +
      ggplot2::geom_line(ggplot2::aes(y = max_set_size,
                                      color = "Maximum Set Size")) +
      ggplot2::geom_line(ggplot2::aes(y = min_set_size,
                                      color = "Minimum Set Size")) +
      ggplot2::scale_x_log10() +
      ggplot2::labs(title = "Effect of IQR Factor on Topology",
                    x = "IQR Factor", y = "Size") +
      ggplot2::theme_minimal()
    attr(results_df, "plot") <- p
  }

  results_df
}


#' Candidate radii for threshold neighbourhoods
#'
#' @description
#' Computes five scales of the data that can serve as radii for the
#' threshold neighbourhoods of [calculate_topology()] and
#' [visualize_topology_thresholds()]. They are descriptive scales, not
#' estimates of an optimal radius.
#'
#' @param data Numeric vector with at least two elements and no missing
#'   values. No default.
#' @return A named \code{list} of five numeric scalars:
#'   \describe{
#'     \item{mean_diff}{Mean of the gaps between consecutive sorted values.}
#'     \item{median_diff}{Median of the gaps between consecutive sorted values.}
#'     \item{sd}{Standard deviation of the data.}
#'     \item{iqr}{Interquartile range divided by 4.}
#'     \item{dbscan}{The \eqn{(k n)}-th smallest of the \eqn{n(n-1)/2}
#'       pairwise distances, with \eqn{k = \lceil \log n \rceil} and the
#'       position capped at the number of distances.}
#'   }
#'
#' @section Methodological notes:
#' The \code{dbscan} entry keeps its historical name, but it is a single
#' order statistic of all pairwise distances; it is not the distance from any
#' observation to its \eqn{k}-th nearest neighbour.
#'
#' @section Dependencies:
#' \code{stats::dist()}, \code{stats::median()}, \code{stats::sd()} and
#' \code{stats::IQR()}.
#'
#' @seealso [calculate_topology()], [visualize_topology_thresholds()].
#'
#' @examples
#' set.seed(1)
#' calculate_thresholds(rnorm(100))
#' @export
calculate_thresholds <- function(data) {
  if (!is.numeric(data) || length(data) < 2L) {
    stop("'data' must be a numeric vector with at least 2 elements.",
         call. = FALSE)
  }
  if (anyNA(data)) {
    stop("'data' must not contain NA values.", call. = FALSE)
  }
  sorted_data <- sort(data)
  diffs <- abs(diff(sorted_data))
  k <- ceiling(log(length(data)))
  dist_sorted <- sort(as.numeric(stats::dist(matrix(data, ncol = 1))))

  list(
    mean_diff = mean(diffs),
    median_diff = stats::median(diffs),
    sd = stats::sd(data),
    iqr = stats::IQR(data) / 4,
    dbscan = dist_sorted[min(k * length(data), length(dist_sorted))]
  )
}


#' Size of the threshold-neighbourhood base for one radius
#'
#' @description
#' Builds the threshold neighbourhoods of the data with a given radius and
#' returns the number of sets in the base they generate, counted as in
#' [analyze_topology_factors()].
#'
#' @param data Numeric vector with at least two elements and no missing
#'   values. No default.
#' @param threshold Non-negative numeric scalar, the radius. No default.
#' @return An \code{integer} scalar: the number of distinct sets among the
#'   empty set, the whole index set and the non-empty intersections of two
#'   neighbourhoods (see [analyze_topology_factors()] for why these are all the
#'   finite intersections).
#'
#' @section Methodological notes:
#' The count includes the empty set and the whole index set. Missing values
#' are rejected; up to version 0.3.0 they were silently dropped from every
#' neighbourhood.
#'
#' @section Dependencies:
#' Base R only.
#'
#' @seealso [analyze_topology_factors()], [calculate_thresholds()].
#'
#' @examples
#' set.seed(1)
#' calculate_topology(rnorm(30), threshold = 0.5)
#' @export
calculate_topology <- function(data, threshold) {
  if (!is.numeric(data) || length(data) < 2L) {
    stop("'data' must be a numeric vector with at least 2 elements.",
         call. = FALSE)
  }
  if (anyNA(data)) {
    stop("'data' must not contain NA values.", call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold) ||
      threshold < 0) {
    stop("'threshold' must be a single non-negative number.", call. = FALSE)
  }
  length(unique(c(list(integer(0), seq_along(data)),
                  threshold_intersections(data, threshold))))
}


# The distinct non-empty intersections of two threshold neighbourhoods (a
# neighbourhood with itself included). The base counted by
# analyze_topology_factors() and calculate_topology() adds the empty set and
# the whole index set to these.
threshold_intersections <- function(data, threshold) {
  n <- length(data)
  subbase <- lapply(seq_len(n), function(i) {
    which(abs(data - data[i]) <= threshold)
  })
  inter <- list()
  for (i in seq_len(n)) {
    for (j in seq(i, n)) {
      s <- intersect(subbase[[i]], subbase[[j]])
      if (length(s) > 0L) {
        inter <- c(inter, list(s))
      }
    }
  }
  unique(inter)
}


#' Compare the candidate radii of calculate_thresholds()
#'
#' @description
#' Computes the five candidate radii of [calculate_thresholds()] and the size
#' of the threshold-neighbourhood base that each one generates, and optionally
#' plots them side by side.
#'
#' @param data Numeric vector with at least two elements and no missing
#'   values. No default.
#' @param plot Logical scalar; if \code{TRUE} (default) the result carries
#'   three \code{ggplot} objects in its \code{"plots"} attribute.
#' @return A \code{data.frame} with the columns \code{method} (character),
#'   \code{threshold} (numeric) and \code{base_size} (integer, counted as in
#'   [calculate_topology()]). If \code{plot = TRUE}, it carries an attribute
#'   \code{"plots"}: a list with the elements \code{threshold},
#'   \code{base_size} and \code{scatter}.
#'
#' @section Methodological notes:
#' The comparison is descriptive: it shows how the radius, and with it the
#' granularity of the base, depends on the scale chosen.
#'
#' @section Dependencies:
#' 'ggplot2' (Imports) for the plots; [calculate_thresholds()] and
#' [calculate_topology()] for the numbers.
#'
#' @seealso [calculate_thresholds()], [calculate_topology()],
#'   [analyze_topology_factors()].
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_point geom_text theme_minimal labs
#' @examples
#' \donttest{
#' set.seed(1)
#' results <- visualize_topology_thresholds(rnorm(50))
#' results
#' }
#' @export
visualize_topology_thresholds <- function(data, plot = TRUE) {
  if (!is.numeric(data) || length(data) < 2L) {
    stop("'data' must be a numeric vector with at least 2 elements.",
         call. = FALSE)
  }
  if (anyNA(data)) {
    stop("'data' must not contain NA values.", call. = FALSE)
  }

  thresholds <- calculate_thresholds(data)
  base_sizes <- vapply(thresholds, function(t) calculate_topology(data, t),
                       integer(1))

  df <- data.frame(
    method = names(thresholds),
    threshold = unlist(thresholds),
    base_size = base_sizes,
    stringsAsFactors = FALSE
  )

  if (plot) {
    p1 <- ggplot2::ggplot(df, ggplot2::aes(x = method, y = threshold)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Threshold comparison by method",
                    x = "Method", y = "Threshold")

    p2 <- ggplot2::ggplot(df, ggplot2::aes(x = method, y = base_size)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Base size comparison by method",
                    x = "Method", y = "Base size")

    p3 <- ggplot2::ggplot(df, ggplot2::aes(x = threshold, y = base_size,
                                           label = method)) +
      ggplot2::geom_point() +
      ggplot2::geom_text(hjust = -0.1, vjust = 0) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Relationship between threshold and base size",
                    x = "Threshold", y = "Base size")

    attr(df, "plots") <- list(threshold = p1, base_size = p2, scatter = p3)
  }

  df
}


#' The discrete topology on n points
#'
#' @description
#' Returns the discrete topology on the indices of the data, the finest
#' topology on a set, in which every subset is open. Only the singleton base
#' is returned, because the full topology has \eqn{2^n} open sets.
#'
#' @param data Vector with at least one element; only its length is used.
#'   No default.
#' @return A \code{list} with:
#'   \describe{
#'     \item{subbase}{List of the \eqn{n} singletons.}
#'     \item{base}{The same list: the singletons form the minimal base.}
#'     \item{topology_type}{Character: \code{"discrete"}.}
#'     \item{n}{Integer, the number of points.}
#'     \item{connected}{Logical: \code{TRUE} only when \eqn{n = 1}; every
#'       singleton is open and closed, so two or more points are
#'       disconnected.}
#'   }
#'
#' @section Methodological notes:
#' The discrete topology is a reference point: every topology produced by
#' the package is coarser than or equal to it.
#'
#' @section Dependencies:
#' Base R only.
#'
#' @seealso [complete_topology()] for the opposite extreme, the indiscrete
#'   topology.
#'
#' @examples
#' result <- simplest_topology(c(1, 2, 3, 4, 5))
#' result$connected
#' @export
simplest_topology <- function(data) {
  if (!is.numeric(data) || length(data) < 1L) {
    stop("'data' must be a numeric vector with at least 1 element.",
         call. = FALSE)
  }
  n <- length(data)
  singletons <- lapply(seq_len(n), function(i) i)

  list(
    subbase = singletons,
    base = singletons,
    topology_type = "discrete",
    n = n,
    connected = n <= 1L
  )
}


#' The topology induced by the complete graph
#'
#' @description
#' Runs [generate_topology()] on the complete graph over the indices of the
#' data. Every closed neighbourhood of the complete graph is the whole index
#' set, so the result is the indiscrete topology, whose only open sets are the
#' empty set and the whole set, and the space is connected.
#'
#' @param data Numeric vector with at least two elements and no missing
#'   values; only its length is used. No default.
#' @param verify_axioms Logical scalar passed to [generate_topology()].
#'   Default \code{FALSE}.
#' @return The \code{list} returned by [generate_topology()]: base \{V\},
#'   two open sets when enumerated, one component.
#'
#' @section Methodological notes:
#' The construction uses closed neighbourhoods \eqn{N[v] = \{v\} \cup N(v)},
#' and in the complete graph \eqn{N[v] = V} for every vertex. The result is a
#' reference point: the coarsest topology on the set.
#'
#' @section Dependencies:
#' The compiled engine of [generate_topology()] (via 'Rcpp').
#'
#' @seealso [simplest_topology()] for the opposite extreme, the discrete
#'   topology.
#'
#' @examples
#' result <- complete_topology(c(1, 2, 3, 4, 5))
#' result$n_open_sets
#' result$connected
#' @export
complete_topology <- function(data, verify_axioms = FALSE) {
  if (!is.numeric(data) || length(data) < 2L) {
    stop("'data' must be a numeric vector with at least 2 elements.",
         call. = FALSE)
  }
  if (anyNA(data)) {
    stop("'data' must not contain NA values.", call. = FALSE)
  }
  n <- length(data)

  vertices <- seq_len(n)
  adjacency <- lapply(vertices, function(v) setdiff(vertices, v))

  generate_topology(
    adjacency = adjacency,
    n_elements = n,
    check_connected = TRUE,
    verify_axioms = verify_axioms
  )
}
