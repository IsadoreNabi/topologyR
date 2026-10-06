# Recomputes, with the topologyR installed in the library given as the first
# argument, the natural visibility edges of the series of cases_030.txt, and
# writes them in the same format to the file given as the second argument.
# Usage: Rscript edges_installed.R <library> <output file>
args <- commandArgs(TRUE)
.libPaths(c(path.expand(args[1]), .libPaths()))
suppressMessages(library(topologyR))
cat("topologyR", as.character(packageVersion("topologyR")), "\n")
lines <- readLines("cases_030.txt")
out <- file(args[2], "w")
for (line in lines) {
  p <- strsplit(line, "|", fixed = TRUE)[[1]]
  y <- as.numeric(strsplit(p[3], ",", fixed = TRUE)[[1]])
  g <- natural_visibility_graph(y)
  e <- if (nrow(g$edges)) paste(g$edges$from, g$edges$to, sep = "-", collapse = " ") else ""
  writeLines(paste(p[1], p[2], p[3], e, sep = "|"), out)
}
close(out)
