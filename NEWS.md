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
