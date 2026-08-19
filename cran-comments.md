## Submission of topologyR 0.3.0 (bug-fix release)

CRAN currently hosts 0.2.0. This release is confined to the
directed-topology module (`R/directed_topology.R` and its tests);
everything else is untouched. Details in NEWS.md:

* `generate_alexandrov_topology()` assumed -- without checking -- that
  every edge of the input runs from a lower to a higher vertex index
  (what directed visibility graphs produce), and returned a valid-looking
  but wrong topology on any other digraph, silently. It now handles
  arbitrary digraphs correctly: the original engine runs unchanged on its
  original input class (verified bit-identical against stored 0.2.0
  goldens), other DAGs are re-indexed along a topological order, and
  cyclic digraphs are computed on the Tarjan condensation and expanded
  back. One additive argument (`expect`) turns the old implicit
  precondition into a verified contract; three new result fields report
  which case applied.

* `bitopology_invariants()` no longer returns `irreversibility_base`:
  the quantity is identically zero by a theorem of the construction
  (the intents-extents duality of formal concept analysis), now proven
  and documented on the manual page. The two base sizes it compared
  remain reported. A property test keeps the theorem in the suite.

No exported function was removed or renamed.

## Test environments

* local: Fedora Linux 44, R 4.6.0

## R CMD check results

0 errors | 0 warnings | 1 note

The note reports non-portable compilation flags; they are injected by the
distribution's system toolchain, none is set by the package's Makevars,
and the note should not reproduce on CRAN's builders.

The test suite (1,465 assertions) passes with no failures, warnings or
skips, including gate batteries with referents external to the engine:
stored 0.2.0 goldens (bit identity on the original input class), a pure-R
Warshall closure, and hand-built condensations.
