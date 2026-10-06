// ======================================================================
// exact_predicate.h -- exact sign of the chord condition of the natural
// visibility criterion, evaluated on the doubles the caller passes
//
// For three observations a, k, b with t_a < t_k < t_b, observation k lies
// strictly below the chord joining a to b if and only if
//
//   D(a, b, k) = x_a (t_b - t_k) + x_b (t_k - t_a) - x_k (t_b - t_a) > 0,
//
// which is the inequality chordBelowT of the formal development, with both
// sides multiplied by t_b - t_a > 0. chord_below_sign() returns the sign of
// D (+1, 0 or -1) computed exactly for every finite input under the
// hypotheses on the arithmetic stated below: the answer is a function of the
// input doubles alone and does not depend on the order in which
// floating-point operations happen to round.
//
// Stage 1, floating-point filter. When the six inputs are at most 2^510 in
// magnitude, the filter computes, one operation per statement and each result
// stored in a volatile variable, l = fl(fl(t_a - t_k) fl(x_b - x_k)),
// r = fl(fl(x_a - x_k) fl(t_b - t_k)), S = fl(|l| + |r|) and det = fl(l - r),
// an evaluation of the orientation determinant O = L - R = -D. If S exceeds
// 2^-930 (about 1.1e-280, written 0x1p-930 and therefore exact in binary64
// and in every wider evaluation format) and |det| > 8 u S, with u = 2^-53, the
// sign of det is the sign of O,
// and the function returns its opposite. Stage 2, exact
// evaluation: the six products of D as integers times powers of two, added in
// fixed-size multiprecision integers aligned at the smallest exponent, and
// compared.
//
// The proofs of both stages, including every bound the code below relies on,
// are written in the help page natural_visibility_exactness (R/exactness.R),
// the single place where they are maintained; the numerical facts of the
// proof of stage 1 are checked in exact rational arithmetic by
// dev/check_filter_constants.py. The help page also states the hypotheses:
// IEEE 754 binary64 doubles in any rounding mode; each operation rounded to
// binary64 directly, or first to a format with at least 64 significand bits,
// or to the x87 extended format with its precision control at 53 bits, and
// then to binary64 on the store; and subnormal numbers handled as IEEE 754
// prescribes; and no floating-point trap enabled during the computation. The
// guards below check at compile time the format of doubles, the evaluation
// method and the options that GCC reveals; ieee_double_environment() checks
// at run time that sums and products of doubles keep 53 significant bits and
// that subnormal numbers are neither flushed nor read as zero; and
// NonStopFloatingPoint, which the callers hold while they compute, disables
// every trap and restores the caller's environment at exit. A run-time
// precision strictly between 53 and 64 bits would be neither covered by the
// proof nor detected. The scheme, a cheap
// floating-point evaluation with a sign certified by a proven threshold and an
// exact evaluation when the threshold does not decide, is the adaptive one of
// Shewchuk (1997).
//
// The code implements the published description of the filter and plain
// schoolbook multiprecision addition; it contains no third-party code.
//
// Reference: Shewchuk, J. R. (1997). Adaptive precision floating-point
// arithmetic and fast robust geometric predicates. Discrete & Computational
// Geometry, 18(3), 305-363. https://doi.org/10.1007/PL00009321
// ======================================================================

#ifndef TOPOLOGYR_EXACT_PREDICATE_H_
#define TOPOLOGYR_EXACT_PREDICATE_H_

// Hypotheses of the filter (proved sufficient in the help page
// natural_visibility_exactness): IEEE 754 binary64 doubles, in any rounding
// mode, each operation rounded to binary64 directly, or first to a format with
// at least 64 significand bits, or to the x87 extended format with its
// precision control at 53 bits; and subnormal numbers handled as IEEE 754
// prescribes. The guards below refuse a build without binary64 doubles, with
// an evaluation method other than binary64 or a format of at least 64
// significand bits, with -ffast-math or -ffinite-math-only, or, under GCC,
// with any option that clears __GCC_IEC_559. Among the options the compilers
// accept, contraction and the reassociation and reciprocal options of Clang
// cannot change the filter: its intermediate values are volatile objects,
// which the C++ standard requires to be written and read back as the program
// says, each statement performs one operation, and there is no floating-point
// division in this file; an accepted
// option can still act through the floating-point environment it may set at
// start-up, of which
// ieee_double_environment() checks the precision of sums and products and the
// handling of subnormal numbers, and NonStopFloatingPoint disables the traps.
// The exact stage uses only integer arithmetic on the decomposition returned
// by frexp(), which is exact.

#include <cfenv>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#if defined(__FAST_MATH__)
#error "topologyR's exact predicate requires IEEE floating-point semantics; do not compile it with -ffast-math."
#endif
#if defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__
#error "topologyR's exact predicate requires IEEE floating-point semantics; do not compile it with -ffinite-math-only."
#endif
#if !defined(FLT_EVAL_METHOD) || \
    !(FLT_EVAL_METHOD == 0 || FLT_EVAL_METHOD == 1 || \
      (FLT_EVAL_METHOD == 2 && LDBL_MANT_DIG >= 64))
#error "topologyR's exact predicate requires operations on doubles evaluated in binary64 (FLT_EVAL_METHOD 0 or 1) or in a format with at least 64 significand bits (FLT_EVAL_METHOD 2)."
#endif
#if defined(__GCC_IEC_559) && __GCC_IEC_559 == 0
#error "topologyR's exact predicate requires IEEE floating-point semantics; an option in use (such as -funsafe-math-optimizations, -fassociative-math, -freciprocal-math or -fno-signed-zeros) breaks it."
#endif

static_assert(std::numeric_limits<double>::is_iec559 &&
              std::numeric_limits<double>::radix == 2 &&
              std::numeric_limits<double>::digits == 53 &&
              std::numeric_limits<double>::min_exponent == -1021 &&
              std::numeric_limits<double>::max_exponent == 1024,
              "topologyR's exact predicate requires IEEE 754 binary64 doubles.");

namespace exactpred {

inline int sign_of(double v) { return (v > 0.0) - (v < 0.0); }

// Whether the floating-point environment meets the two hypotheses of the
// filter that the compiler cannot check: (1) addition and multiplication of
// doubles keep 53 significant bits, which excludes the x87 unit with its
// precision control set to 24 bits; (2) subnormal numbers are neither flushed
// to zero nor read as zero (flush-to-zero and denormals-are-zero modes). The
// tests of (1) use results that are exact in binary64, 1 + 2^-52 and
// (1 + 2^-26)^2 = 1 + 2^-25 + 2^-52, so they do not depend on the rounding
// mode; under flush-to-zero, half of the smallest normal double becomes 0, and
// under denormals-are-zero the smallest subnormal compares equal to 0.
inline bool ieee_double_environment() {
  volatile double one = 1.0;
  volatile double ulp = 0x1p-52;
  volatile double sum = one + ulp;
  volatile double back = sum - one;
  volatile double near_one = 0x1p0 + 0x1p-26;
  volatile double square = near_one * near_one;
  volatile double smallest_normal = std::numeric_limits<double>::min();
  volatile double half = 0.5;
  volatile double halved = smallest_normal * half;
  volatile double tiny = std::numeric_limits<double>::denorm_min();
  const double expected_square = 0x1p0 + 0x1p-25 + 0x1p-52;
  return back == ulp && square == expected_square && halved > 0.0 && tiny > 0.0;
}

// Holds the IEEE 754 non-stop mode while it lives: feholdexcept() saves the
// caller's floating-point environment, clears the exception flags and disables
// every trap, so that an inexact or underflow trap that some library enabled
// cannot interrupt the filter; the destructor restores the saved environment
// with fesetenv(), without raising the exceptions the computation signalled.
class NonStopFloatingPoint {
 public:
  NonStopFloatingPoint() : held_(std::feholdexcept(&saved_) == 0) {}
  ~NonStopFloatingPoint() {
    if (held_) std::fesetenv(&saved_);
  }
  bool held() const { return held_; }
  NonStopFloatingPoint(const NonStopFloatingPoint&) = delete;
  NonStopFloatingPoint& operator=(const NonStopFloatingPoint&) = delete;

 private:
  std::fenv_t saved_;
  bool held_;
};

// Finite double -> (negative?, mantissa m < 2^53, exponent e) with
// |v| = m * 2^e. Zero gives m = 0.
struct Split {
  bool neg;
  uint64_t m;
  int e;
};

inline Split split_double(double v) {
  Split s;
  s.neg = v < 0.0;
  if (v == 0.0) {
    s.m = 0; s.e = 0;
    return s;
  }
  int E = 0;
  double f = std::frexp(std::fabs(v), &E);           // |v| = f 2^E, f in [0.5, 1)
  s.m = static_cast<uint64_t>(std::ldexp(f, 53));    // exact: f has <= 53 bits
  s.e = E - 53;
  return s;
}

// a * b for a, b < 2^53, as four 32-bit limbs (little endian).
inline void mul53(uint64_t a, uint64_t b, uint32_t out[4]) {
  const uint64_t mask = 0xFFFFFFFFULL;
  uint64_t a0 = a & mask, a1 = a >> 32;              // a1 < 2^21
  uint64_t b0 = b & mask, b1 = b >> 32;              // b1 < 2^21
  uint64_t p00 = a0 * b0;                            // < 2^64
  uint64_t p01 = a0 * b1;                            // < 2^53
  uint64_t p10 = a1 * b0;                            // < 2^53
  uint64_t p11 = a1 * b1;                            // < 2^42
  uint64_t w0 = p00 & mask;
  uint64_t c = p00 >> 32;
  uint64_t mid = c + (p01 & mask) + (p10 & mask);    // < 3 * 2^32
  uint64_t w1 = mid & mask;
  c = mid >> 32;
  uint64_t hi = c + (p01 >> 32) + (p10 >> 32) + (p11 & mask);
  uint64_t w2 = hi & mask;
  c = hi >> 32;
  uint64_t w3 = c + (p11 >> 32);
  out[0] = static_cast<uint32_t>(w0);
  out[1] = static_cast<uint32_t>(w1);
  out[2] = static_cast<uint32_t>(w2);
  out[3] = static_cast<uint32_t>(w3);
}

// acc += (4-limb value) * 2^shift, shift >= 0.
inline void add_shifted(std::vector<uint32_t>& acc, const uint32_t v[4], int shift) {
  int limb = shift / 32, bit = shift % 32;
  uint32_t s[5];
  if (bit == 0) {
    s[0] = v[0]; s[1] = v[1]; s[2] = v[2]; s[3] = v[3]; s[4] = 0;
  } else {
    s[0] = v[0] << bit;
    s[1] = (v[1] << bit) | (v[0] >> (32 - bit));
    s[2] = (v[2] << bit) | (v[1] >> (32 - bit));
    s[3] = (v[3] << bit) | (v[2] >> (32 - bit));
    s[4] = v[3] >> (32 - bit);
  }
  uint64_t carry = 0;
  size_t i = static_cast<size_t>(limb);
  for (int k = 0; k < 5; k++, i++) {
    uint64_t sum = static_cast<uint64_t>(acc[i]) + s[k] + carry;
    acc[i] = static_cast<uint32_t>(sum);
    carry = sum >> 32;
  }
  while (carry) {
    uint64_t sum = static_cast<uint64_t>(acc[i]) + carry;
    acc[i] = static_cast<uint32_t>(sum);
    carry = sum >> 32;
    i++;
  }
}

// Exact sign of sum_k c_k x_k y_k over six terms whose coefficients c_k are
// +1 or -1: negative[k] is true when c_k = -1. Passing the coefficients as
// booleans makes any other coefficient impossible to pass.
inline int exact_sign6(const double x[6], const double y[6], const bool negative[6]) {
  Split sx[6], sy[6];
  bool nonzero[6];
  int e[6];
  int emin = 0, emax = 0;
  bool any = false;
  for (int k = 0; k < 6; k++) {
    sx[k] = split_double(x[k]);
    sy[k] = split_double(y[k]);
    nonzero[k] = sx[k].m != 0 && sy[k].m != 0;
    e[k] = sx[k].e + sy[k].e;
    if (nonzero[k]) {
      if (!any || e[k] < emin) emin = e[k];
      if (!any || e[k] > emax) emax = e[k];
      any = true;
    }
  }
  if (!any) return 0;
  // Each product is below 2^106 and is shifted by at most emax - emin bits;
  // the sum of six such values needs three more bits. Two spare limbs absorb
  // the final carries.
  size_t limbs = static_cast<size_t>((emax - emin + 106 + 3) / 32 + 3);
  std::vector<uint32_t> pos(limbs, 0), neg(limbs, 0);
  for (int k = 0; k < 6; k++) {
    if (!nonzero[k]) continue;
    uint32_t prod[4];
    mul53(sx[k].m, sy[k].m, prod);
    bool term_negative = (sx[k].neg != sy[k].neg) != negative[k];
    add_shifted(term_negative ? neg : pos, prod, e[k] - emin);
  }
  for (size_t i = limbs; i-- > 0;) {
    if (pos[i] != neg[i]) return pos[i] > neg[i] ? 1 : -1;
  }
  return 0;
}

// Sign of D(a, b, k) = x_a (t_b - t_k) + x_b (t_k - t_a) - x_k (t_b - t_a):
// +1 if k is strictly below the chord from a to b, 0 if it is on it, -1 if
// it is above. Exact for every finite input under the hypotheses stated at
// the top of this file and in the help page natural_visibility_exactness.
inline int chord_below_sign(double ta, double xa, double tb, double xb,
                            double tk, double xk) {
  const double bound = 0x1p510;                 // no intermediate can overflow
  if (std::fabs(ta) <= bound && std::fabs(xa) <= bound &&
      std::fabs(tb) <= bound && std::fabs(xb) <= bound &&
      std::fabs(tk) <= bound && std::fabs(xk) <= bound) {
    // One operation per statement, each result stored and rounded to double.
    volatile double d1 = ta - tk;
    volatile double d2 = xb - xk;
    volatile double d3 = xa - xk;
    volatile double d4 = tb - tk;
    volatile double l = d1 * d2;
    volatile double r = d3 * d4;
    volatile double al = std::fabs(l);
    volatile double ar = std::fabs(r);
    volatile double S = al + ar;
    volatile double det = l - r;                            // orientation O
    const double s = S, o = det;
    if (s > 0x1p-930 && std::fabs(o) > std::ldexp(s, -50)) {  // 8 u S
      return -sign_of(o);
    }
  }
  // D = x_a t_b - x_a t_k + x_b t_k - x_b t_a - x_k t_b + x_k t_a.
  const double x[6] = {xa, xa, xb, xb, xk, xk};
  const double y[6] = {tb, tk, tk, ta, tb, ta};
  const bool negative[6] = {false, true, false, true, true, false};
  return exact_sign6(x, y, negative);
}

} // namespace exactpred

#endif // TOPOLOGYR_EXACT_PREDICATE_H_
