# Comparison of the natural visibility engines of topologyR 0.3.0 and 0.4.0

`make_series.R` generates 2,800 series in eight classes (seed 20261005) and
records the edges that topologyR 0.3.0 returns; `compare_exact.py` compares
them with the definition applied, in exact rational arithmetic, to the stored
doubles, and counts the edges that 0.3.0 adds and misses. `edges_installed.R`
recomputes the edges of the same series with another installed version; on the
edges of 0.4.0 the comparison finds no discrepancy. This directory is excluded
from the built package.

    Rscript make_series.R <library with topologyR 0.3.0>
    python3 compare_exact.py cases_030.txt
    Rscript edges_installed.R <library with topologyR 0.4.0> cases_040.txt
    python3 compare_exact.py cases_040.txt
