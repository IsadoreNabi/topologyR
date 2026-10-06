# ----------------------------------------------------------------------
# Validation of count arguments
# ----------------------------------------------------------------------

# A count argument (a number of points or a limit on a number of sets) must
# be a single finite whole number in [lower, .Machine$integer.max]. It is
# checked before the conversion to integer, so that a fraction, a string, a
# logical value, NA or a vector is rejected instead of being truncated or
# coerced; the value is returned as an integer.
check_count <- function(x, name, lower) {
  ok <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) &&
    x == trunc(x) && x >= lower && x <= .Machine$integer.max
  if (!ok) {
    what <- switch(as.character(lower),
                   "0" = "a non-negative whole number",
                   "1" = "a positive whole number",
                   sprintf("a whole number >= %d", lower))
    stop(sprintf("'%s' must be %s (a single finite value, not a fraction).",
                 name, what), call. = FALSE)
  }
  as.integer(x)
}
