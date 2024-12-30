#' @title Topology Analysis Functions for Time Series Data
#' @description Functions for analyzing topological properties of time series data
#' @author José Mauricio Gómez Julián
#' @details
#' This module provides three main functions for analyzing the connectivity
#' of topological structures, particularly focused on economic time series:
#' - is_topology_connected: Uses an undirected graph approach
#' - is_topology_connected2: Uses a directed graph approach
#' - is_topology_connected_manual: Uses a manual checking approach
#'
#' Created as part of research on economic cycle analysis.
#' Check if a topology is connected using undirected graph approach
#'
#' @param topology A list of sets representing the topology
#' @return logical TRUE if topology is connected, FALSE otherwise
#' @export
#' @examples
#' topology <- list(c(1,2,3), c(3,4,5))
#' is_topology_connected(topology)
is_topology_connected <- function(topology) {
  if (!length(topology)) return(FALSE)

  elements <- unique(unlist(topology))
  n <- max(elements)  # Changed from length(elements)

  # Create undirected graph from topology
  edges <- matrix(0, nrow = n, ncol = n)
  for (set in topology) {
    for (i in set) {
      for (j in set) {
        if (i != j && i <= n && j <= n) {  # Added boundary check
          edges[i, j] <- 1
          edges[j, i] <- 1
        }
      }
    }
  }

  # Perform DFS to check connectivity
  visited <- rep(FALSE, n)
  stack <- c(elements[1])

  while (length(stack) > 0) {
    v <- stack[1]
    stack <- stack[-1]
    if (!visited[v]) {
      visited[v] <- TRUE
      neighbors <- which(edges[v, ] == 1)
      stack <- c(stack, elements[neighbors])
    }
  }

  all(visited[elements])  # Only check elements that are present
}

#' Check if a topology is connected using directed graph approach
#'
#' @param topology A list of sets representing the topology
#' @return logical TRUE if topology is connected, FALSE otherwise
#' @export
#' @examples
#' topology <- list(c(1,2,3), c(3,4,5))
#' is_topology_connected2(topology)
is_topology_connected2 <- function(topology) {
  # Handle empty topology case
  if (!length(topology)) return(FALSE)

  # Get unique elements and determine matrix size
  elements <- unique(unlist(topology))
  n <- max(elements)

  # Initialize adjacency matrix
  edges <- matrix(0, nrow = n, ncol = n)

  # Build directed edges
  for (set in topology) {
    if (length(set) > 1) {
      # Get values in current set
      values <- sort(set)
      # Create edges between consecutive elements
      for (i in 1:(length(values)-1)) {
        current_val <- values[i]
        next_val <- values[i+1]
        edges[current_val, next_val] <- 1
      }
    }
  }

  # Initialize DFS variables
  visited <- rep(FALSE, n)
  start <- min(elements)
  stack <- c(start)

  # Perform DFS
  while (length(stack) > 0) {
    current <- stack[1]
    stack <- stack[-1]

    if (!visited[current]) {
      visited[current] <- TRUE
      # Find all unvisited neighbors
      neighbors <- which(edges[current,] == 1)
      stack <- c(stack, neighbors[!visited[neighbors]])
    }
  }

  # Check if all elements in topology are visited
  return(all(visited[elements]))
}

#' Check if a topology is connected using manual checking approach
#'
#' @param topology A list of sets representing the topology
#' @return logical TRUE if topology is connected, FALSE otherwise
#' @export
#' @examples
#' topology <- list(c(1,2,3), c(3,4,5))
#' is_topology_connected_manual(topology)
is_topology_connected_manual <- function(topology) {
  # Get all unique elements in the topology
  all_elements <- unique(unlist(topology))

  # Get elements from original set (assuming they go from 1 to n)
  n <- max(all_elements)
  original_elements <- 1:n

  # Check if each element appears in at least one set
  for (element in original_elements) {
    if (!any(sapply(topology, function(set) element %in% set))) {
      return(FALSE)
    }
  }

  return(TRUE)
}

#' Calculate topology characteristics for different IQR factors
#'
#' @description
#' This function analyzes how different IQR (Interquartile Range) factors affect
#' the topology's characteristics. It helps users determine the optimal factor
#' for their specific data by showing how the factor choice impacts the base size
#' and set sizes in the resulting topology.
#'
#' @param data Numeric vector containing the data to analyze
#' @param factors Numeric vector of factors to test (default: c(1, 2, 4, 8, 16))
#' @return A data frame containing the following columns for each factor:
#'   \itemize{
#'     \item factor: The IQR division factor used
#'     \item threshold: The calculated threshold (IQR/factor)
#'     \item base_size: Number of sets in the topological base
#'     \item max_set_size: Size of the largest set
#'     \item min_set_size: Size of the smallest set
#'   }
#'
#' @details
#' The function works by:
#' 1. Calculating different thresholds using IQR/factor
#' 2. Creating a subbase using these thresholds
#' 3. Generating the base from intersections of subbase elements
#' 4. Analyzing the resulting topology's characteristics
#'
#' A larger factor results in a smaller threshold, which typically leads to
#' a finer topology with more distinct sets but smaller set sizes.
#'
#' @examples
#' # Generate sample data
#' data <- rnorm(100)
#'
#' # Analyze topology with default factors
#' results <- analyze_topology_factors(data)
#' print(results)
#'
#' # Use custom factors
#' custom_results <- analyze_topology_factors(data, factors = c(2, 4, 6))
#' print(custom_results)
#'
#' @export
analyze_topology_factors <- function(data, factors = NULL, plot = TRUE) {
  library(ggplot2)
  # Set default factors if not provided
  if (is.null(factors)) {
    factors <- c(1, 2, 4, 8, 16)
  }

  # Inner function to calculate topology for a single factor
  calculate_topology_with_factor <- function(data, factor) {
    threshold <- IQR(data) / factor
    n <- length(data)

    # Calculate subbase using threshold
    subbase <- lapply(1:n, function(i) {
      which(abs(data - data[i]) <= threshold)
    })

    # Generate base (including empty set and full set)
    base <- list(integer(0), 1:n)
    for (i in 1:n) {
      for (j in i:n) {
        intersection <- intersect(subbase[[i]], subbase[[j]])
        if (length(intersection) > 0) {
          base <- c(base, list(intersection))
        }
      }
    }
    unique(base)
  }

  # Calculate results for each factor
  results <- lapply(factors, function(f) {
    topology <- calculate_topology_with_factor(data, f)
    list(
      factor = f,
      threshold = IQR(data) / f,
      base_size = length(topology),
      max_set_size = max(sapply(topology, length)),
      min_set_size = min(sapply(topology, length))
    )
  })

  # Convert results to data frame
  results_df <- do.call(rbind, lapply(results, data.frame))

  # Optional plotting
  if (plot) {
    p <- ggplot(results_df, aes(x = factor)) +
      geom_line(aes(y = base_size, color = "Tamaño de la Base")) +
      geom_line(aes(y = max_set_size, color = "Tamaño Máximo del Conjunto")) +
      geom_line(aes(y = min_set_size, color = "Tamaño Mínimo del Conjunto")) +
      scale_x_log10() +
      labs(title = "Efecto del Factor IQR en la Topología",
           x = "Factor IQR", y = "Tamaño") +
      theme_minimal()

    # Print the plot
    print(p)
  }

  # Return the results data frame
  return(results_df)
}

#' Calculate multiple threshold methods for topology analysis
#'
#' @param data Numeric vector to calculate thresholds for
#' @return A list of threshold values using different methods
#' @export
calculate_thresholds <- function(data) {
  list(
    mean_diff = mean(abs(diff(sort(data)))),
    median_diff = median(abs(diff(sort(data)))),
    sd = sd(data),
    iqr = IQR(data) / 4,
    dbscan = {
      k <- ceiling(log(length(data)))
      sort(dist(matrix(data, ncol=1)))[k * length(data)]
    }
  )
}

#' Visualize thresholds and base sizes for different methods
#'
#' @param data Numeric vector to analyze
#' @export
calculate_topology <- function(data, threshold) {
  n <- length(data)
  subbase <- lapply(1:n, function(i) {
    which(abs(data - data[i]) <= threshold)
  })

  base <- list(integer(0), 1:n)
  for (i in 1:n) {
    for (j in i:n) {
      intersection <- intersect(subbase[[i]], subbase[[j]])
      if (length(intersection) > 0) {
        base <- c(base, list(intersection))
      }
    }
  }
  base <- unique(base)

  length(base)  # Retorna el tamaño de la base como medida de complejidad
}

#' @param data Numeric vector to analyze
#' @export
visualize_topology_thresholds <- function(data, plot = TRUE) {
  library(ggplot2)

  # Calculate thresholds
  thresholds <- calculate_thresholds(data)

  # Calculate base sizes (using the existing calculate_topology function)
  base_sizes <- sapply(thresholds, function(t) calculate_topology(data, t))

  # Create data frame
  df <- data.frame(
    method = names(thresholds),
    threshold = unlist(thresholds),
    base_size = base_sizes
  )

  # Threshold comparison plot
  p1 <- ggplot(df, aes(x = method, y = threshold)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = "Comparación de umbrales por método",
         x = "Método", y = "Umbral")

  # Base size comparison plot
  p2 <- ggplot(df, aes(x = method, y = base_size)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = "Comparación de tamaños de base por método",
         x = "Método", y = "Tamaño de la base")

  # Threshold vs base size plot
  p3 <- ggplot(df, aes(x = threshold, y = base_size, label = method)) +
    geom_point() +
    geom_text(hjust = -0.1, vjust = 0) +
    theme_minimal() +
    labs(title = "Relación entre umbral y tamaño de base",
         x = "Umbral", y = "Tamaño de la base")

  # Print plots
  print(p1)
  print(p2)
  print(p3)

  # Return data frame for further analysis
  return(df)
}
