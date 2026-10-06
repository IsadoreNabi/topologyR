#include <Rcpp.h>
#include <vector>
#include <stack>
#include "exact_predicate.h"

// The comparisons of both criteria and the filter of the exact predicate
// assume IEEE 754 double arithmetic: 53 significant bits, subnormal numbers
// handled as the standard prescribes, and no floating-point trap. Each entry
// point holds the non-stop mode while it computes (NonStopFloatingPoint, which
// restores the caller's environment at exit) and refuses to run when the mode
// cannot be set or when a flush-to-zero or denormals-are-zero mode, or the x87
// precision control set to 24 bits, which some libraries switch on, would
// silently change the graphs.
static void require_ieee_environment(const exactpred::NonStopFloatingPoint& non_stop) {
  if (!non_stop.held()) {
    Rcpp::stop("The floating-point environment could not be set to the IEEE "
               "754 non-stop mode, without traps; topologyR's visibility "
               "graphs require it.");
  }
  if (!exactpred::ieee_double_environment()) {
    Rcpp::stop("The floating-point environment does not provide IEEE 754 "
               "double arithmetic (subnormal numbers are flushed to zero, or "
               "operations keep fewer than 53 significant bits, a mode that "
               "some library has switched on); topologyR's visibility graphs "
               "require IEEE 754 double arithmetic.");
  }
}

// ==========================================================================
// Horizontal Visibility Graph (HVG) — O(n) stack-based algorithm
//
// Two observations (t_a, y_a) and (t_b, y_b) with t_a < t_b are connected
// iff every intermediate observation (t_c, y_c) satisfies:
//   y_c < min(y_a, y_b)
//
// Algorithm: modified "next greater or equal element" using a monotone stack.
// The stack maintains a strictly decreasing sequence of y-values.
//
// Key correctness detail: when popping with <=, if ANY popped element has
// y == y[i], then an equal-valued point exists between the remaining stack
// top and i, blocking that edge (since y_c = y_i = min(y_top, y_i) violates
// the strict inequality). The `had_equal` flag handles this.
//
// Reference: Luque, B., Lacasa, L., Ballesteros, F., & Luque, J. (2009).
// Horizontal visibility graphs: exact results for random time series.
// Physical Review E, 80(4), 046103.
// ==========================================================================

// [[Rcpp::export]]
Rcpp::List hvg_cpp(Rcpp::NumericVector y) {
  exactpred::NonStopFloatingPoint non_stop;
  require_ieee_environment(non_stop);
  int n = y.size();
  if (n < 2) {
    return Rcpp::List::create(
      Rcpp::Named("from") = Rcpp::IntegerVector(),
      Rcpp::Named("to")   = Rcpp::IntegerVector()
    );
  }

  std::vector<int> from_vec;
  std::vector<int> to_vec;
  from_vec.reserve(n * 2);
  to_vec.reserve(n * 2);

  std::stack<int> stk;

  for (int i = 0; i < n; i++) {
    bool had_equal = false;

    // Pop all elements with y <= y[i]. Each popped element j has its
    // "next greater or equal to the right" at i, meaning all values
    // between j and i are < y[j] <= y[i], so (j, i) is an HVG edge.
    while (!stk.empty() && y[stk.top()] <= y[i]) {
      int j = stk.top();
      stk.pop();
      if (y[j] == y[i]) had_equal = true;
      from_vec.push_back(j + 1);  // 1-based for R
      to_vec.push_back(i + 1);
    }

    // The remaining stack top (if any) has y > y[i]. It is connected
    // to i iff no intermediate point has y == y[i] (which would block
    // the horizontal line at height y[i] = min(y_top, y_i)).
    if (!stk.empty() && !had_equal) {
      from_vec.push_back(stk.top() + 1);
      to_vec.push_back(i + 1);
    }

    stk.push(i);
  }

  return Rcpp::List::create(
    Rcpp::Named("from") = Rcpp::wrap(from_vec),
    Rcpp::Named("to")   = Rcpp::wrap(to_vec)
  );
}


// ==========================================================================
// Natural Visibility Graph (NVG) -- O(n^2) decisions in the worst case,
// O(n log n) when the maxima split the intervals in bounded proportions
//
// Two observations (t_a, y_a) and (t_b, y_b) with t_a < t_b are connected
// iff every intermediate observation (t_c, y_c) lies strictly below the
// straight segment joining them:
//   y_c (t_b - t_a) < y_a (t_b - t_c) + y_b (t_c - t_a).
// The instants t are any strictly increasing doubles; the R wrapper passes
// 1, ..., n when the caller gives none.
//
// Algorithm: divide and conquer on the global maximum. The three facts it
// relies on hold for arbitrary strictly increasing instants:
//   (i)  If p maximizes y on [lo, hi] and a < p < b, the chord from a to b
//        takes at t_p the convex combination
//        (y_a (t_b - t_p) + y_b (t_p - t_a)) / (t_b - t_a) <= max(y_a, y_b)
//        <= y_p, so p is not strictly below it: no edge crosses p.
//   (ii) For j < p, every k in (j, p) is strictly below the chord from j to p
//        iff the slope from j to p is smaller than the slope from every such k
//        to p. The smallest of those slopes belongs to the last vertex found
//        visible, m, so j is visible iff there is no intermediate vertex or m
//        lies strictly below the chord from j to p. Symmetrically on the
//        right: j > p is visible iff m lies strictly below the chord from p
//        to j.
//   (iii) Edges inside each half depend only on the points of that half.
// Each visibility decision is therefore the sign of one chord condition,
// evaluated exactly on the input doubles by exactpred::chord_below_sign(), so
// the edge set is the natural visibility graph of the numbers received,
// independent of rounding. Ties are handled by the strict inequality: a
// point on the chord blocks.
//
// Reference: Lacasa, L., Luque, B., Ballesteros, F., Luque, J., & Nuño, J. C.
// (2008). From time series to complex networks: The visibility graph.
// Proceedings of the National Academy of Sciences, 105(13), 4972-4975.
// ==========================================================================

namespace {

// Index of the leftmost maximum of y[lo..hi] (inclusive).
int find_max_idx(const double* y, int lo, int hi) {
  int idx = lo;
  for (int i = lo + 1; i <= hi; i++) {
    if (y[i] > y[idx]) {
      idx = i;
    }
  }
  return idx;
}

void nvg_recursive(const double* t, const double* y, int lo, int hi,
                   std::vector<int>& from_vec, std::vector<int>& to_vec) {
  if (lo >= hi) return;

  // Two adjacent points are always mutually visible
  if (hi - lo == 1) {
    from_vec.push_back(lo + 1);
    to_vec.push_back(hi + 1);
    return;
  }

  int p = find_max_idx(y, lo, hi);

  // Left scan: m is the last vertex found visible from p
  int m = -1;
  for (int j = p - 1; j >= lo; j--) {
    if (m < 0 ||
        exactpred::chord_below_sign(t[j], y[j], t[p], y[p], t[m], y[m]) > 0) {
      from_vec.push_back(j + 1);
      to_vec.push_back(p + 1);
      m = j;
    }
  }

  // Right scan, symmetric
  m = -1;
  for (int j = p + 1; j <= hi; j++) {
    if (m < 0 ||
        exactpred::chord_below_sign(t[p], y[p], t[j], y[j], t[m], y[m]) > 0) {
      from_vec.push_back(p + 1);
      to_vec.push_back(j + 1);
      m = j;
    }
  }

  nvg_recursive(t, y, lo, p - 1, from_vec, to_vec);
  nvg_recursive(t, y, p + 1, hi, from_vec, to_vec);
}

} // anonymous namespace


// [[Rcpp::export]]
Rcpp::List nvg_cpp(Rcpp::NumericVector y, Rcpp::NumericVector t) {
  exactpred::NonStopFloatingPoint non_stop;
  require_ieee_environment(non_stop);
  int n = y.size();
  if (t.size() != n) Rcpp::stop("'y' and 't' must have the same length.");
  if (n < 2) {
    return Rcpp::List::create(
      Rcpp::Named("from") = Rcpp::IntegerVector(),
      Rcpp::Named("to")   = Rcpp::IntegerVector()
    );
  }

  std::vector<int> from_vec;
  std::vector<int> to_vec;
  from_vec.reserve(n * 3);
  to_vec.reserve(n * 3);

  nvg_recursive(t.begin(), y.begin(), 0, n - 1, from_vec, to_vec);

  return Rcpp::List::create(
    Rcpp::Named("from") = Rcpp::wrap(from_vec),
    Rcpp::Named("to")   = Rcpp::wrap(to_vec)
  );
}


// Vectorized access to the exact chord predicate, for the test suite and for
// external verification. Not exported to users.
// [[Rcpp::export]]
Rcpp::IntegerVector chord_below_sign_cpp(Rcpp::NumericVector ta,
                                         Rcpp::NumericVector xa,
                                         Rcpp::NumericVector tb,
                                         Rcpp::NumericVector xb,
                                         Rcpp::NumericVector tk,
                                         Rcpp::NumericVector xk) {
  exactpred::NonStopFloatingPoint non_stop;
  require_ieee_environment(non_stop);
  R_xlen_t n = ta.size();
  if (xa.size() != n || tb.size() != n || xb.size() != n ||
      tk.size() != n || xk.size() != n) {
    Rcpp::stop("all arguments must have the same length.");
  }
  Rcpp::IntegerVector out(n);
  for (R_xlen_t i = 0; i < n; i++) {
    if (!std::isfinite(ta[i]) || !std::isfinite(xa[i]) ||
        !std::isfinite(tb[i]) || !std::isfinite(xb[i]) ||
        !std::isfinite(tk[i]) || !std::isfinite(xk[i])) {
      Rcpp::stop("all arguments must be finite.");
    }
  }
  for (R_xlen_t i = 0; i < n; i++) {
    out[i] = exactpred::chord_below_sign(ta[i], xa[i], tb[i], xb[i],
                                         tk[i], xk[i]);
  }
  return out;
}
