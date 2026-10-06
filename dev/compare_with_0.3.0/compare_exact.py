"""Compares the natural visibility edges recorded in a cases file with the
definition applied, in exact rational arithmetic, to the stored doubles, and
classifies the discrepancies: edges the engine adds (extra) and edges it
misses (missing). Two exact references are used: the literal O(n^3)
definition, and the maximum-slope form O(n^2) (j is visible from i when its
slope from i exceeds that of every k between them); they are checked against
each other for n <= 40. Instants: 1..n.

Usage: python3 compare_exact.py cases_030.txt
"""
import sys, collections
from fractions import Fraction as F


def nvg_literal(Y, T):
    n = len(Y); E = set()
    for i in range(n):
        for j in range(i + 1, n):
            if all(Y[k] * (T[j] - T[i]) < Y[i] * (T[j] - T[k]) + Y[j] * (T[k] - T[i])
                   for k in range(i + 1, j)):
                E.add((i + 1, j + 1))
    return E


def nvg_max_slope(Y, T):
    n = len(Y); E = set()
    for i in range(n):
        m = None
        for j in range(i + 1, n):
            s = (Y[j] - Y[i]) / (T[j] - T[i])
            if m is None or s > m:
                E.add((i + 1, j + 1))
                m = s
    return E


def read_cases(path):
    for line in open(path):
        k, kind, ys, es = line.rstrip("\n").split("|")[:4]
        y = [float.fromhex(v) for v in ys.split(",")]
        got = {tuple(map(int, e.split("-"))) for e in es.split()} if es else set()
        yield k, kind, y, got


if __name__ == "__main__":
    total, disc, only_missing = collections.Counter(), collections.Counter(), collections.Counter()
    extra, missing = collections.Counter(), collections.Counter()
    crossed = 0
    for k, kind, y, got in read_cases(sys.argv[1]):
        Y = [F(v) for v in y]; T = [F(i) for i in range(1, len(y) + 1)]
        ref = nvg_max_slope(Y, T)
        if len(y) <= 40:
            assert nvg_literal(Y, T) == ref, ("the two exact references disagree", k)
            crossed += 1
        total[kind] += 1
        if got != ref:
            disc[kind] += 1
            extra[kind] += len(got - ref); missing[kind] += len(ref - got)
            if not (got - ref):
                only_missing[kind] += 1
    for kind in sorted(total):
        print(f"{kind:20s} series={total[kind]:4d} discrepant={disc[kind]:4d} only_missing={only_missing[kind]:4d} "
              f"extra_edges={extra[kind]:4d} missing_edges={missing[kind]:4d}")
    print("total:", sum(total.values()), "| discrepant:", sum(disc.values()),
          "| extra edges:", sum(extra.values()), "| missing edges:", sum(missing.values()),
          "| literal = max-slope checks:", crossed)
