## topologyR 0.4.0

CRAN hosts 0.3.0. The changes are listed in NEWS.md; in short:

* The licence changes from MIT to GPL (>= 3). Every commit of the package is
  by its author and maintainer.
* The visibility-graph constructors accept unequally spaced instants
  (`times`), and a new function `time_reverse()` returns a series read
  backwards with its instants reflected.
* The natural visibility criterion is decided exactly on the input doubles,
  under IEEE 754 double arithmetic (a floating-point filter whose sign is
  certified by a proven threshold, and an exact integer evaluation when the
  threshold does not decide); the proofs and their hypotheses are in the new
  help page `natural_visibility_exactness`.
* The completeness flag `complete` is renamed `topology_complete` and is `NA`
  when no enumeration was requested; incomplete computations no longer report
  undue exactness.
* The documentation and the user manual were corrected and completed; the
  manual's code is now run by the test suite.

There are no reverse dependencies on CRAN (`tools::package_dependencies()`
with `reverse = TRUE` over all dependency fields, 2026-10-05).

## Test environments

* local: Fedora Linux 44, R 4.6.1, g++ (GCC) 16.2.1.

## R CMD check results

With the distribution's compiler flags: 0 errors | 0 warnings | 1 note.

The note reports non-portable compilation flags. All of them come from the
system R configuration (`R CMD config CXXFLAGS`); the package's `Makevars`
only sets `CXX_STD = CXX17`. With `CXXFLAGS = -g -O2 -Wall -pedantic` in a
user Makevars the check returns `Status: OK`, with no compiler warnings.

The test suite passes with no failures, warnings or skips, including
batteries whose referents are external to the engine: exact signs computed
outside the package in rational arithmetic, the 0.2.0 goldens of the
Alexandrov engine, a pure-R reachability closure, and literal definitions on
data where double arithmetic is exact.
