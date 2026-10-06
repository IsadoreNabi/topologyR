# Generates the 2,800 test series of the comparison between the natural
# visibility engines of topologyR 0.3.0 and 0.4.0, with the edges that 0.3.0
# returns. topologyR 0.3.0 (CRAN, 2026-08-20) must be installed in the
# library given as the first argument. Output: cases_030.txt, one series per
# line: id|class|values as hexadecimal doubles|edges "from-to" of 0.3.0.
# Usage: Rscript make_series.R <library with topologyR 0.3.0>
args <- commandArgs(TRUE)
.libPaths(c(path.expand(args[1]), .libPaths()))
suppressMessages(library(topologyR)); stopifnot(packageVersion("topologyR") == "0.3.0")
set.seed(20261005)
cases <- list()
add <- function(kind, y) cases[[length(cases) + 1L]] <<- list(kind = kind, y = y)
for (n in c(10, 50, 200)) for (r in 1:300) add(paste0("rnorm_n", n), rnorm(n))
for (r in 1:500) add("integers_0_10", as.numeric(sample(0:10, sample(5:40, 1), replace = TRUE)))
for (r in 1:500) {
  n <- sample(5:40, 1); a <- round(runif(1, -5, 5), 1); b <- round(runif(1, -1, 1), 1)
  y <- round(a + b * seq_len(n) + sample(c(0, 0, 0, 0.1, -0.1), n, replace = TRUE), 1)
  add("decimals_near_lines", y)
}
for (r in 1:500) add("decimals_2_places", round(runif(sample(5:60, 1), 0, 100), 2))
for (r in 1:200) add("large_integers", 1e10 + as.numeric(sample(0:20, sample(5:30, 1), replace = TRUE)))
for (r in 1:200) { n <- sample(5:40, 1); y <- round(cumsum(rnorm(n)), 3); add("walk_3_places", y) }
out <- file("cases_030.txt", "w")
for (k in seq_along(cases)) {
  y <- cases[[k]]$y; g <- natural_visibility_graph(y)
  e <- if (nrow(g$edges)) paste(g$edges$from, g$edges$to, sep = "-", collapse = " ") else ""
  writeLines(paste(k, cases[[k]]$kind, paste(sprintf("%a", y), collapse = ","), e, sep = "|"), out)
}
close(out); cat("series:", length(cases), "\n")
