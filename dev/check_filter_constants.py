"""Check, in exact rational arithmetic, the numerical facts of the proof
of the error bound of the floating-point filter, written in the help page
natural_visibility_exactness (R/exactness.R).

Run from the root of the package:

    python3 dev/check_filter_constants.py

Every line must end in True.
"""

from fractions import Fraction as F

u = F(1, 2 ** 53)                       # unit roundoff of binary64
e = 2 * u * (1 + F(1, 2 ** 10))         # relative error of one operation
c = F(1, 2 ** 930)                      # the threshold of S, 0x1p-930 in the code
absolute = F(1, 2 ** 1071)              # absolute term of step (b)

checks = [
    ("(i) one directed rounding to the 64-bit extended format, then one to "
     "double: (1 + 2^-63)(1 + 2u) - 1 <= e",
     (1 + F(1, 2 ** 63)) * (1 + 2 * u) - 1 <= e),
    ("(b) g = (1 + e)^3 - 1 < 6.01 u",
     (1 + e) ** 3 - 1 < F(601, 100) * u),
    ("(b) g / (1 - g) < 6.02 u",
     ((1 + e) ** 3 - 1) / (1 - ((1 + e) ** 3 - 1)) < F(602, 100) * u),
    ("(b) g / (1 - g) 2^-1072 + 2^-1072 <= 2^-1071",
     ((1 + e) ** 3 - 1) / (1 - ((1 + e) ** 3 - 1)) * F(1, 2 ** 1072)
     + F(1, 2 ** 1072) <= absolute),
    ("(d) 8 u (1 - e) / (1 + e) > 7.99 u",
     8 * u * (1 - e) / (1 + e) > F(799, 100) * u),
    ("(d) S > 2^-930 implies |l| + |r| > 2^-930 / (1 + e) > 1.1e-280",
     c / (1 + e) > F(11, 10) * F(10) ** -280),
    ("(d) 2^-1071 < 1.97 u 1.1e-280",
     absolute < F(197, 100) * u * F(11, 10) * F(10) ** -280),
    ("(d) the scaling 2^-50 S is exact: 2^-50 2^-930 is above 2^-1022",
     F(1, 2 ** 50) * c > F(1, 2 ** 1022)),
    ("(d) 6.02 u + 1.97 u <= 7.99 u",
     F(602, 100) + F(197, 100) <= F(799, 100)),
    ("(o) inputs at most 2^510: differences <= 2^511, products <= 2^1022, "
     "sum <= 2^1023 < largest double",
     2 * F(2) ** 510 == F(2) ** 511 and F(2) ** 511 * F(2) ** 511 == F(2) ** 1022
     and 2 * F(2) ** 1022 < (2 - F(1, 2 ** 52)) * F(2) ** 1023),
    ("2^-50 = 8 u", F(1, 2 ** 50) == 8 * u),
]

for label, ok in checks:
    print(f"{label}: {ok}")
