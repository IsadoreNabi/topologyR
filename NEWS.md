# topologyR 0.4.0

## The package is now licensed under GPL (>= 3)

`DESCRIPTION` declares `License: GPL (>= 3)`, and `LICENSE.md` carries the full
text of version 3 of the GNU General Public License, excluded from the built
package. Every commit in the history of the repository is by the author, so the
change of licence needs no other consent. The two-line file required by the
previous MIT declaration is gone.

## Unequally spaced instants

`horizontal_visibility_graph()`, `natural_visibility_graph()` and
`generate_bitopology()` gain an argument `times`: `NULL` (the instants
1, ..., n, as before) or a numeric, `Date` or date-time vector of finite,
strictly increasing instants. The natural criterion uses them as abscissae;
the horizontal criterion depends only on their order, so for it they are
validated and recorded but change nothing. The graph objects gain a field
`times` with the instants used.

The new function `time_reverse()` returns a series read backwards together
with its instants reflected (`-rev(times)`), the operation under which the
forward and backward topologies are exchanged (theorem T1). Reversing the
values alone over the same instants is a different operation when the gaps
between consecutive instants are not palindromic, and for the natural
criterion it can change the graph: with the instants 0, 1, 4 and the values 0,
2, 4 it creates an edge that the original series does not have.

## The natural visibility criterion is decided exactly

Each visibility decision of `natural_visibility_graph()` is the sign of
`y_a (t_b - t_k) + y_b (t_k - t_a) - y_k (t_b - t_a)`, and it is now computed
exactly on the doubles received, under the hypotheses on the floating-point
arithmetic stated in the new help page `?natural_visibility_exactness`: the sign of a floating-point evaluation is
accepted only when the result exceeds a threshold above which that sign is
proved correct, and otherwise the quantity is evaluated in integer arithmetic
on the mantissas and exponents of its six products (Shewchuk, 1997, for the
adaptive scheme). Up to 0.3.0 the engine
compared rounded slopes. On 2,800 test series the two engines return
identical edge lists in 2,558 cases. The other 242, all among the 500 series
of rounded decimals lying almost on a straight line, differ from the
definition applied to the stored numbers: across them 0.3.0 misses 596 edges
that the definition gives and adds 67 that it does not, and 194 of them only
miss edges. The comparison is reproducible with the scripts in
`dev/compare_with_0.3.0`, which also find no discrepancy with 0.4.0. For decimal data the
documentation of `natural_visibility_graph()` explains when and how to obtain
the graph of the recorded decimals instead.

Consequences: theorem T1 holds exactly in the computation for every input,
and the graph does not change under exact positive affine changes of the
values or of the time scale, when the new values are finite doubles equal to the
exact images. The exact evaluation runs when the threshold
does not certify the sign, which happens near a tie (for instance on data
with many collinear points), when an input exceeds 2^510 in magnitude, or
when the computed sum of the magnitudes of the two products does not exceed
2^-930, about 1.1e-280. It costs more than the floating-point
evaluation, so data that trigger it often take longer.

The threshold is valid in every IEEE 754 rounding mode, with evaluation in
binary64, in a format with at least 64 significand bits or in the x87 format at
53 bits. The floating-point stage is used only when the six inputs of a
decision are at most 2^510 in magnitude, so that no intermediate value
overflows (a directed rounding mode can turn an overflow into the largest
finite double instead of an infinity, with a relative error no longer bounded
by 2u), and each of its intermediate values
is stored, so that every operation is rounded to double before the next one
uses it and the compiler cannot fuse a product into the subtraction or
reorder the operations. The build stops with an error when doubles are not
binary64, when operations on doubles are evaluated in another format than
binary64 or one with at least 64 significand bits (`FLT_EVAL_METHOD`), with
`-ffast-math` or `-ffinite-math-only`, and, under GCC, with any option that
makes the compiler set `__GCC_IEC_559` to 0 (among them
`-funsafe-math-optimizations`, `-fassociative-math`, `-freciprocal-math` and
`-fno-signed-zeros`); Clang does not reveal its reassociation options, which
the stored intermediates make harmless. `horizontal_visibility_graph()` and
`natural_visibility_graph()` compute in the non-stop mode of IEEE 754, with
every floating-point trap disabled and the caller's environment restored at
exit, and stop with an error when the floating-point environment does not
provide IEEE 754 double arithmetic: subnormal numbers flushed to zero or read
as zero, or operations with fewer than 53 significant bits (as with the x87
precision control set to 24 bits), conditions that some libraries create.
The proofs of both stages, their hypotheses on the arithmetic and what is
checked of them are written in the new help page
`?natural_visibility_exactness`.

## Incomplete computations are flagged, and only what is exact is reported

* The flag `complete` of `generate_topology()`, `generate_alexandrov_topology()`
  and `complete_topology()` is renamed `topology_complete`, the name the
  accompanying article uses. It is `NA` when no enumeration was requested
  (it used to be `FALSE`, which also meant "truncated"), `FALSE` when the
  enumeration was truncated or the base was truncated, and `TRUE` only when
  the enumerated family is the whole topology. Up to 0.3.0 a truncated base
  could be enumerated to completion and flagged complete although the union
  closure of a truncated base is not the topology. `verify_axioms` now runs
  only when `topology_complete` is `TRUE`.
* The connected components are exact even when `max_base_sets` truncates the
  base, because the engine always keeps the whole subbase; the documentation
  said they could be approximate, and now gives the proof.
* `bitopology_invariants()` reports the two base sizes, and the gains with
  respect to the Alexandrov topology, as `NA` when the corresponding base was
  truncated, and adds the fields `forward_base_complete` and
  `backward_base_complete`. It stops with an error when a topology was
  computed with `check_connected = FALSE` (it used to count zero components)
  or when `n_elements` does not match the points of the components.
* `bitopology_invariants()` decides pairwise connectedness exactly, from the
  two subbases and without enumerating either topology: the space is pairwise
  connected if and only if the digraph with an arc x -> y when y lies in the
  smallest forward-open set of x or x in the smallest backward-open set of y
  is strongly connected, and a strongly connected component that no arc
  leaves is a witness. The answer is never `NA`. Up to 0.3.0 the check needed
  both topologies enumerated, answered `NA` otherwise, and ran only when both
  reported `complete`, a flag that a truncated base did not clear.
* Count arguments (`n_elements`, `max_open_sets`, `max_base_sets`) are
  validated before they are converted to integers: a fraction, a string, a
  logical value, a missing value or a vector is rejected with an error,
  instead of being truncated or coerced (`max_open_sets = 2.7` used to mean
  2). `generate_topology()` rejects a negative `max_open_sets` or a
  non-positive `max_base_sets`, and so does `generate_alexandrov_topology()`
  with a negative `max_open_sets`.
* With `check_connected = FALSE`, `connected` is a logical `NA`, as
  documented; it used to be an integer `NA`.
* `verify_axioms = TRUE` checks every axiom of a topology on a finite set: the
  empty set and the whole set are members and the family is closed under the
  union and the intersection of any two members. It used to check pairwise
  intersections only, although its result is called `axioms_ok`.

## Documentation corrected

* The pair of topologies does not "quantify temporal irreversibility", and a
  reversible process does not make them homeomorphic: theorem T1 has no
  converse. The manual pages now state T1 with its proof and say what it does
  and does not imply; `irreversibility_components` is no longer described as a
  necessary condition for reversibility, and `asymmetry_direction` is no
  longer tied to "gradual expansions and abrupt contractions". The paragraph
  attributing an irreversibility prediction to Lacasa and Toral (2010), whose
  article does not discuss irreversibility, is removed; Lacasa et al. (2012)
  is cited for the visibility-graph approach to irreversibility.
* The base of `generate_alexandrov_topology()` is the minimal base of the
  topology, which is not closed under intersection in general; the page said
  it was.
* The construction is a closed-neighbourhood variant of that of Nada, El Atik
  and Atef (2018), who use the open post classes of the adjacency relation;
  the documentation no longer attributes the closed neighbourhoods to them.
* Kelly (1963) is cited for bitopological spaces only, with its correct
  volume, and no longer next to "temporal irreversibility" in the package
  description or as the source of pairwise connectedness.
* The bitopological space of a directed visibility graph is pairwise
  disconnected for every series with at least two observations: every final
  segment of indices is open in the forward topology and every initial
  segment in the backward one. The documentation says so, and that the
  pairwise check therefore carries no information about a series.
* The base-size gains are counts of the construction, not measures of how
  much finer the Nada topology is; the two gains are always equal when
  defined, by the theorem that makes the two base sizes equal.
* The example of `is_topology_connected_exact()` labelled "A connected
  topology" was not a topology; both examples now are.
* Every exported function documents its parameters, return value, design
  decisions, methodological notes and dependencies, and every reference is
  checked against the list in `dev/REFERENCIAS_APA7.md`.
* The user manual is rewritten to match the package, and its code is run.

## Legacy functions

* `analyze_topology_factors()`: `max_set_size` and `min_set_size` are taken
  over the intersections of the threshold neighbourhoods; they used to include
  the empty set and the whole set and were always 0 and n. `factors` must be
  positive. The documentation no longer promises an optimal factor.
* `calculate_topology()` rejects missing values instead of dropping them.
* `calculate_thresholds()`: the `dbscan` entry is documented as what it is, an
  order statistic of all pairwise distances, not a nearest-neighbour distance.
* `is_topology_connected2()` is documented as neither a necessary nor a
  sufficient condition for connectedness, with a counterexample to each;
  `is_topology_connected()` as a necessary condition only.
* `complete_topology()` no longer claims a length limit of 64, and its page
  says that the result is the indiscrete topology.

## Build

* Compiled objects and an R history file are no longer part of the source
  tree, and the OpenMP flags, unused by the code, are removed from
  `Makevars`.
* The two manuals of version 0.1.0, `inst/rmd/USERMANUAL.Rmd` and its
  Spanish version `inst/rmd/MANUALDEUSUARIO.Rmd`, are deleted: they repeated
  errors that the current manual corrects. The manual is
  `inst/manuals/user_manual.Rmd`, whose code the tests run; the deleted files
  remain in the history of the repository and in version 0.3.0 on CRAN.
* A `.zenodo.json` file carries the metadata of the archived release (title,
  author with ORCID, version and licence); it is excluded from the built
  package.

## Tests

New batteries check the exact predicate against 1,800 cases whose signs are
computed in exact rational arithmetic outside the package, by the
deterministic generator `dev/make_exact_predicate_fixture.py` (the numerical
inequalities of the error bound are checked in exact arithmetic by
`dev/check_filter_constants.py`), the decimal recipe, the natural graph
against its literal definition on data where double arithmetic is exact,
the equality of `times = NULL` with every exact affine image of 1, ..., n,
theorem T1 with random unequally spaced instants, the negative control of
the instants 0, 1, 4, the six states of a truncated computation, the
invariants of incomplete computations, the validation of count arguments,
and pairwise connectedness: the exact criterion against algorithm 5 of the
accompanying article, applied to the complete enumerations of every digraph
with at most three vertices and of 300 random digraphs with cycles, its
decision without enumeration, the segment witnesses of every directed
visibility graph, and a digraph with cycles that is pairwise connected with
or without enumeration. The time-reversal property test over
palindromes added in the development of this version is kept.

# topologyR 0.3.0

## `generate_alexandrov_topology()` accepts any digraph and says which case it got

Up to 0.2.0 the reachability propagation assumed -- without checking --
that every edge goes from a lower to a higher vertex index, which is what
directed visibility graphs produce but not what the signature promised. On
any other input it returned a valid-looking and wrong topology, silently:
a directed cycle, or a DAG whose vertices arrive in a different order,
both produced upsets computed from incomplete reachability. Measured
before the repair, 145 of 300 random digraphs with feedback came back
wrong, and the reversed chain `3 -> 2 -> 1` -- an acyclic input -- came
back wrong too.

The construction is now correct for every directed graph, in three
regimes reported by the result itself:

* **Index-ordered DAG** (the visibility-graph case): the original
  single-pass bitset engine runs unchanged. Verified bit-identical to
  0.2.0 against stored goldens on 64 real visibility graphs.
* **DAG in arbitrary vertex order**: vertices are re-indexed along a
  topological order before propagation.
* **Digraph with directed cycles**: mutually reachable vertices are
  topologically indistinguishable, so the topology is computed on the
  condensation (Tarjan's algorithm) and expanded back; the expansion is a
  bijection on open sets, so connectivity, components, enumeration and
  axiom checks all transfer exactly.

Three new result fields declare what happened: `input_index_ordered`,
`input_acyclic`, and `collapsed_classes` (the non-trivial strongly
connected components -- classes the topology cannot tell apart; a finding
about the system, not an error). A new argument
`expect = c("any", "index_ordered_dag")` lets a caller whose semantics
require the strict visibility contract get an error instead of a
generalization: with `"index_ordered_dag"` the old precondition is
*verified* in O(m) and a violating edge is named in the error message.

The gate batteries run three independent referents: stored 0.2.0 goldens
(bit identity on regime 1), a pure-R Warshall closure (permuted DAGs and
cyclic digraphs -- where the tests also demonstrate the raw 0.2.0 motor
disagreeing with the truth), and hand-built condensations.

## `bitopology_invariants()` no longer reports `irreversibility_base`

The field was the normalized asymmetry of the forward and backward base
sizes -- and it is **identically zero by a theorem of the construction**,
for every digraph, not only for visibility graphs: the two base families
are the intents and extents of the formal context of the visibility
relation, and the Galois connection makes them equinumerous (Ganter &
Wille 1999; the duality is Birkhoff's). An identically null index does
not measure, and its presence invited reading it, so it was removed
rather than kept as a decoy. The two base sizes remain reported -- each
is informative on its own -- and the theorem, with its proof sketch and
the one exception (truncated closures, flagged by `base_complete =
FALSE`), is documented on the manual page. A property test over random
digraphs keeps the theorem in the suite as a regression watch on the
engine.

`irreversibility_components` stays, documented as the coarse
count-comparison it is.

## The Alexandrov branch is documented as constant on visibility graphs

In any visibility graph consecutive observations see each other, so
reachability is the total order and the Alexandrov topology is the same
chain of upsets for every series: `alexandrov_base_size = n`,
`alexandrov_components = 1`, and the resolution-gain fields reduce to
`base_size - n`. The manual pages of `generate_bitopology()` and
`bitopology_invariants()` now say so, and say where the branch is
genuinely informative (general digraphs). Callers that consume only the
Nada side can pass `alexandrov = FALSE`. A property test asserts the
chain on random series for both graph types.
