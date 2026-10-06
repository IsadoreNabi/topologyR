"""Generate tests/testthat/fixtures/exact_predicate_cases.csv.

The fixture is the external referent of the exact chord predicate of
natural_visibility_graph(): 1,800 sextuples (t_a, x_a, t_b, x_b, t_k, x_k) of
doubles in nine regimes, 200 each, with the sign of

    D(a, b, k) = x_a (t_b - t_k) + x_b (t_k - t_a) - x_k (t_b - t_a)

computed in exact rational arithmetic (fractions.Fraction), outside the
package and independently of it. The ninth regime holds 200 more. Doubles are written in hexadecimal notation,
so that R reads back exactly the same values.

The generator is deterministic (seeded) and uses only the Python standard
library. Run it from the root of the package:

    python3 dev/make_exact_predicate_fixture.py

Regimes:
  normal                    six standard normal draws
  collinear_dyadic          three points on a line with dyadic coefficients,
                            exactly representable: D = 0
  decimal_collinear         decimal values of a line, rounded to one place
  decimal_times_and_values  decimal instants and values of a line
  extreme_exponents         random significands with exponents from -1074
                            to 1023
  subnormal                 small multiples of 2^-1074
  near_collinear_ulp        a point on a line moved by up to two ulps
  large_integers            instants near 2^60 and values that are large
                            multiples of 2^40
  beyond_filter_bound       inputs above 2^510 in magnitude, which skip the
                            floating-point filter, half of them in the
                            pattern in which a filter without that bound
                            fails under directed rounding (an overflow
                            turned into the largest finite double)
"""

import math
import random
import sys
from fractions import Fraction

SEED = 20261005
PER_KIND = 200
OUT = "tests/testthat/fixtures/exact_predicate_cases.csv"


def exact_sign(ta, xa, tb, xb, tk, xk):
    ta, xa, tb, xb, tk, xk = (Fraction(v) for v in (ta, xa, tb, xb, tk, xk))
    d = xa * (tb - tk) + xb * (tk - ta) - xk * (tb - ta)
    return (d > 0) - (d < 0)


def random_double(rng, emin, emax):
    """A double with a uniformly drawn 53-bit significand and exponent."""
    m = (1 << 52) + rng.getrandbits(52)
    e = rng.randint(emin, emax)
    return rng.choice((-1.0, 1.0)) * math.ldexp(m, e - 52)


def normal(rng):
    return tuple(rng.gauss(0.0, 1.0) for _ in range(6))


def collinear_dyadic(rng):
    ta = rng.randint(-50, 50) / 8
    tk = ta + rng.randint(1, 40) / 8
    tb = tk + rng.randint(1, 40) / 8
    a = rng.randint(-100, 100) / 16
    b = rng.randint(-100, 100) / 16
    return (ta, a + b * ta, tb, a + b * tb, tk, a + b * tk)


def decimal_collinear(rng):
    ta = float(rng.randint(0, 30))
    tk = ta + rng.randint(1, 10)
    tb = tk + rng.randint(1, 10)
    a = round(rng.uniform(-5, 5), 1)
    b = round(rng.uniform(-1, 1), 1)
    return (ta, round(a + b * ta, 1), tb, round(a + b * tb, 1),
            tk, round(a + b * tk, 1))


def decimal_times_and_values(rng):
    tt = round(rng.uniform(0, 10), 2)
    ta, tb, tk = tt, tt + 0.7, tt + 0.3
    return (ta, round(0.3 * ta + 0.1, 2), tb, round(0.3 * tb + 0.1, 2),
            tk, round(0.3 * tk + 0.1, 2))


def extreme_exponents(rng):
    return tuple(random_double(rng, -1074, 1023) for _ in range(6))


def subnormal(rng):
    tiny = math.ldexp(1.0, -1074)
    return tuple(rng.randint(-50, 50) * tiny for _ in range(6))


def near_collinear_ulp(rng):
    ta = rng.random()
    tk = ta + rng.random()
    tb = tk + rng.random()
    a, b = rng.gauss(0.0, 1.0), rng.gauss(0.0, 1.0)
    xk = a + b * tk
    d = rng.choice((-2, -1, 0, 1, 2))
    xk = xk + d * sys.float_info.epsilon * abs(xk)
    return (ta, a + b * ta, tb, a + b * tb, tk, xk)


def large_integers(rng):
    big = 2.0 ** 60
    def value():
        return rng.randint(-10 ** 6, 10 ** 6) * 2.0 ** 40
    ta = big + rng.randint(0, 1000)
    tb = big + 2000 + rng.randint(0, 1000)
    tk = big + 1000.5
    return (ta, value(), tb, value(), tk, value())


def beyond_filter_bound(rng):
    if rng.random() < 0.5:
        return tuple(random_double(rng, 1000, 1023) for _ in range(6))
    big = math.ldexp(1.5 + 0.25 * rng.random(), 1023)
    ta, tk = big, -big * (0.9 + 0.1 * rng.random())
    tb = rng.uniform(-1, 1) * big
    return (ta, rng.uniform(-2, 2), tb, rng.uniform(-2, 2), tk, rng.uniform(-2, 2))


REGIMES = [
    ("normal", normal),
    ("collinear_dyadic", collinear_dyadic),
    ("decimal_collinear", decimal_collinear),
    ("decimal_times_and_values", decimal_times_and_values),
    ("extreme_exponents", extreme_exponents),
    ("subnormal", subnormal),
    ("near_collinear_ulp", near_collinear_ulp),
    ("large_integers", large_integers),
    ("beyond_filter_bound", beyond_filter_bound),
]


def main():
    rng = random.Random(SEED)
    lines = ["kind,ta,xa,tb,xb,tk,xk,sign"]
    for kind, draw in REGIMES:
        for _ in range(PER_KIND):
            v = draw(rng)
            assert all(math.isfinite(z) for z in v)
            lines.append(",".join([kind] + [z.hex() for z in v] +
                                  [str(exact_sign(*v))]))
    with open(OUT, "w", newline="\n") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
