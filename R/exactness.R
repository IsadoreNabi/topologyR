#' Exactness of the natural visibility criterion
#'
#' @description
#' [natural_visibility_graph()] decides every pair by the sign of
#' \deqn{D(a, b, k) = x_a (t_b - t_k) + x_b (t_k - t_a) - x_k (t_b - t_a)}
#' for an observation \eqn{k} between \eqn{a} and \eqn{b}: \eqn{k} lies
#' strictly below the chord joining \eqn{a} to \eqn{b} if and only if
#' \eqn{D > 0}. That sign is computed exactly for every finite input, in a
#' floating-point environment that provides IEEE 754 double arithmetic,
#' under the hypotheses stated below; the graph is therefore the natural
#' visibility graph of the
#' doubles received, whatever the order in which floating-point operations
#' round. This page proves it. The computation has two stages, the adaptive
#' scheme of Shewchuk (1997): a floating-point evaluation whose sign is
#' accepted only when the result exceeds a threshold above which that sign is
#' proved correct, and an exact evaluation in integer arithmetic otherwise. The code is in the file
#' \code{src/exact_predicate.h} of the package sources.
#'
#' @section The floating-point stage:
#' \eqn{D} equals \eqn{-O}, where \eqn{O = L - R} is the orientation
#' determinant of the points \eqn{(t_a, x_a)}, \eqn{(t_b, x_b)},
#' \eqn{(t_k, x_k)}, with \eqn{L = (t_a - t_k)(x_b - x_k)} and
#' \eqn{R = (x_a - x_k)(t_b - t_k)}. When the six inputs have magnitude at most
#' \eqn{2^{510}}, the code computes the four differences, the products
#' \eqn{l = \mathrm{fl}(\mathrm{fl}(t_a - t_k)\,\mathrm{fl}(x_b - x_k))} and
#' \eqn{r} likewise, \eqn{S = \mathrm{fl}(|l| + |r|)} and
#' \eqn{\mathrm{det} = \mathrm{fl}(l - r)}, an evaluation of \eqn{O}, where
#' \eqn{\mathrm{fl}(\cdot)} denotes the double stored after one operation: its
#' exact result rounded to binary64, either directly or after a first rounding
#' to the evaluation format (the proof bounds both roundings together). When
#' \eqn{S > 2^{-930}} (about \eqn{1.1 \times 10^{-280}}; the code writes it as
#' the hexadecimal literal \code{0x1p-930}, which binary64 represents exactly,
#' and so does every wider format in which the comparison may be evaluated, so
#' its value does not depend on the compiler) and
#' \eqn{|\mathrm{det}| > 8uS}, with \eqn{u = 2^{-53}}, it accepts the sign of
#' \eqn{\mathrm{det}} as the sign of \eqn{O} and returns its opposite as the
#' sign of \eqn{D}; otherwise it goes to the exact stage. Each of the ten values
#' above (four differences, two products, two absolute values, \eqn{S} and
#' \eqn{\mathrm{det}}) is stored in a volatile variable, one operation per
#' statement, so it is the result of a single operation rounded to double: the
#' compiler can neither fuse a product into the subtraction nor reorder the
#' operations, and a value computed in a wider format is rounded to double
#' when it is stored.
#'
#' \emph{Proof that the accepted sign is the sign of} \eqn{O}. Assume IEEE 754
#' binary64 doubles in any of the four rounding modes, each operation rounded
#' to binary64 directly, or first to a format with at least 64 significand
#' bits (such as the x87 extended format) or to the x87 extended format with
#' its precision control at 53 bits, and then to binary64 on the store. Under
#' any other evaluation the threshold is not established; the next section
#' says which part of these hypotheses is checked.
#' (o) No operation overflows: the differences have magnitude at most
#' \eqn{2^{511}}, the products at most \eqn{2^{1022}} and \eqn{S} at most
#' \eqn{2^{1023}}, and rounding is monotone, so every rounded value obeys the
#' same bounds.
#' (i) Hence every operation whose exact result has magnitude at least
#' \eqn{2^{-1022}} returns it with relative error at most
#' \eqn{e = 2u(1 + 2^{-10})}, which bounds one directed rounding to double
#' (\eqn{2u}) after one rounding to a 64-bit significand (\eqn{2^{-63}}); a
#' 53-bit significand gives an error of at most \eqn{2u} that the store, exact
#' in that range, does not increase; an
#' operation whose exact result is smaller is exact if it is a subtraction of
#' two doubles (the difference is a multiple of \eqn{2^{-1074}} below
#' \eqn{2^{-1022}}, hence a double), and otherwise has absolute error below
#' \eqn{2^{-1073}}, the sum of at most \eqn{2^{-1074}} from the store and less
#' than \eqn{2^{-1074}} from a previous rounding in a wider format.
#' (a) Each difference is exact or within relative error \eqn{e}.
#' (b) Each product is \eqn{l = L(1 + \theta) + a_l} with
#' \eqn{|\theta| \le g = (1 + e)^3 - 1 < 6.01u} and \eqn{|a_l| < 2^{-1073}},
#' the absolute term covering an underflow; likewise \eqn{r}. Hence
#' \eqn{|(l - r) - O| \le g(|L| + |R|) + 2^{-1072}}, and since
#' \eqn{|L| \le (|l| + 2^{-1073}) / (1 - g)} and \eqn{g / (1 - g) < 6.02u},
#' \eqn{|(l - r) - O| < 6.02u(|l| + |r|) + 2^{-1071}}.
#' (c) \eqn{\mathrm{det} = \mathrm{fl}(l - r)} has the sign of \eqn{l - r}: a
#' non-zero \eqn{l - r} is at least \eqn{2^{-1074}} in magnitude, a double that
#' every rounding mode keeps on the same side of zero, and
#' \eqn{\mathrm{det} = 0} only if \eqn{l = r}; and
#' \eqn{|l - r| \ge |\mathrm{det}| / (1 + e)} (the subtraction is exact when
#' its result is below \eqn{2^{-1022}}).
#' (d) \eqn{S > 2^{-930}} implies that the exact sum \eqn{|l| + |r|} is at
#' least \eqn{2^{-1022}} (a smaller sum of two doubles would be exact, and
#' \eqn{S} would equal it), hence at least \eqn{S / (1 + e)}, and
#' \eqn{S \ge (|l| + |r|)(1 - e)}. If \eqn{|\mathrm{det}| > 8uS}, then
#' \eqn{|l - r| > 8u(1 - e)(|l| + |r|) / (1 + e) > 7.99u(|l| + |r|)}, and with
#' \eqn{|l| + |r| > 2^{-930} / (1 + e) > 1.1 \times 10^{-280}} the term
#' \eqn{2^{-1071}} is below
#' \eqn{1.97u(|l| + |r|)}; so \eqn{|l - r| > |(l - r) - O|}, \eqn{l - r} and
#' \eqn{O} have the same sign, and so does \eqn{\mathrm{det}}. \eqn{\square}
#'
#' Inputs above \eqn{2^{510}} in magnitude go straight to the exact stage.
#' Without that bound, a directed rounding mode can turn an overflow into the
#' largest finite double instead of an infinity (under rounding toward zero,
#' every overflow; under rounding downward or upward, the overflows of one
#' sign), and the relative error, which approaches 1 as the exact result grows,
#' is no longer bounded by \eqn{2u}. With \eqn{M} the
#' largest double, \eqn{2^{1024} - 2^{971}}, the instants
#' \eqn{t_a = -M < t_k = 2^{1023} < t_b = M} and the values
#' \eqn{x_a = 5/4}, \eqn{x_k = 0}, \eqn{x_b = -1/2} give
#' \eqn{D = -2^{1021} - 3 \cdot 2^{969} < 0}, but under rounding toward zero a
#' filter without the bound computes the difference \eqn{t_a - t_k}, of
#' magnitude \eqn{M + 2^{1023}}, as \eqn{-M}, finds \eqn{\mathrm{det} < 0} far
#' above the threshold, and returns \eqn{+1}.
#'
#' @section The exact stage:
#' \emph{Proof that it returns the sign of} \eqn{D}.
#' (1) Decomposition. For a finite non-zero double \eqn{v}, the C function
#' \code{frexp()} returns, exactly, \eqn{f \in [1/2, 1)} and \eqn{E} with
#' \eqn{|v| = f\,2^E}, \eqn{E} from \eqn{-1073} to \eqn{1024}. \eqn{f} has at
#' most 53 significant bits (a subnormal \eqn{v} has fewer, and \code{frexp()}
#' normalizes it), so \eqn{m = f\,2^{53}}, computed exactly by \code{ldexp()},
#' is an integer with \eqn{2^{52} \le m < 2^{53}}, and \eqn{|v| = m\,2^{e}}
#' with \eqn{e = E - 53}, from \eqn{-1126} to \eqn{971}. Zero gives
#' \eqn{m = 0}.
#' (2) Products. \eqn{D} is the sum of the six terms
#' \eqn{x_a t_b}, \eqn{-x_a t_k}, \eqn{x_b t_k}, \eqn{-x_b t_a},
#' \eqn{-x_k t_b}, \eqn{x_k t_a}, each a product of two inputs with the
#' coefficient \eqn{+1} or \eqn{-1}; the function that adds them receives,
#' for each term, only whether its coefficient is \eqn{-1}, so no other
#' coefficient can be passed. Each non-zero term is
#' \eqn{\pm (m_x m_y)\,2^{e_x + e_y}} with \eqn{m_x m_y < 2^{106}}, and the
#' product \eqn{m_x m_y} is computed exactly in four words of 32 bits:
#' writing \eqn{a = a_1 2^{32} + a_0} and \eqn{b = b_1 2^{32} + b_0} with
#' \eqn{a_1, b_1 < 2^{21}}, the four partial products are below \eqn{2^{64}},
#' \eqn{2^{53}}, \eqn{2^{53}} and \eqn{2^{42}}, the middle sum is below
#' \eqn{3 \cdot 2^{32}} and the high sum below \eqn{2^{33}}, so no operation on
#' 64-bit unsigned integers overflows.
#' (3) Alignment. Let \eqn{e_{\min}} and \eqn{e_{\max}} be the smallest and
#' the largest exponent of the non-zero terms and
#' \eqn{\Delta = e_{\max} - e_{\min}}. Then \eqn{2^{-e_{\min}} D = P - N},
#' where \eqn{P} and \eqn{N} are the sums of the terms of each sign, each
#' term \eqn{\pm (m_x m_y) 2^{e_k}}, with \eqn{e_k = e_x + e_y}, contributing
#' the integer \eqn{m_x m_y 2^{e_k - e_{\min}}}; \eqn{2^{-e_{\min}} D} has the
#' sign of \eqn{D}. Each shifted term is below
#' \eqn{2^{106 + \Delta}}, so \eqn{P} and \eqn{N} are below
#' \eqn{6 \cdot 2^{106 + \Delta} < 2^{109 + \Delta}}.
#' (4) Capacity. Each accumulator has
#' \eqn{L = \lfloor (\Delta + 109) / 32 \rfloor + 3} words of 32 bits, at least
#' \eqn{\Delta + 174} bits, so no carry leaves it; a term shifted by
#' \eqn{s \le \Delta} bits is written into the words
#' \eqn{\lfloor s / 32 \rfloor} to \eqn{\lfloor s / 32 \rfloor + 4}, and
#' \eqn{\lfloor \Delta / 32 \rfloor + 4 < L}. Since \eqn{\Delta} is at most
#' \eqn{1942 + 2252 = 4194}, \eqn{L} is at most 137.
#' (5) Comparison. \eqn{P} and \eqn{N} are compared word by word from the
#' most significant one, which orders non-negative integers written in base
#' \eqn{2^{32}}, and \eqn{P = N} gives \eqn{D = 0}. After \code{frexp()} and
#' \code{ldexp()}, which are exact, the stage uses integer arithmetic only, so
#' no rounding, overflow or underflow can occur, and neither the rounding mode
#' nor the compiler options matter. \eqn{\square}
#'
#' @section Compiler options and run-time environment:
#' At compile time the build stops with an error unless doubles have the
#' binary64 format (radix 2, 53 significant digits, its exponent range) and
#' operations on doubles are evaluated in binary64 (\code{FLT_EVAL_METHOD} 0
#' or 1) or in a format with at least 64 significand bits
#' (\code{FLT_EVAL_METHOD} 2 with \code{LDBL_MANT_DIG} at least 64). It also
#' stops with \code{-ffast-math} or \code{-ffinite-math-only}, which let the
#' compiler assume that no infinity or NaN occurs, and, under GCC, with every
#' option that makes the compiler set \code{__GCC_IEC_559} to 0 (among them
#' \code{-funsafe-math-optimizations}, \code{-fassociative-math},
#' \code{-freciprocal-math} and \code{-fno-signed-zeros}). Among the options
#' that the compilers accept, contraction (with any compiler) and the
#' reassociation and reciprocal options of Clang, which does not reveal them
#' through its predefined macros, cannot change the floating-point stage.
#' Each of the ten intermediate values of that stage is a volatile object, and
#' the C++ standard requires every access to a volatile object to be performed
#' as the program writes it: each result is stored into its object, and every
#' operation that uses it reads it from that object, so the operand is the
#' double stored there and not a fused or wider intermediate; each statement
#' that computes one of them performs a single operation; the final test
#' reads \eqn{S} and \eqn{\mathrm{det}} once and compares
#' \eqn{|\mathrm{det}|} with \eqn{2^{-50} S} through an absolute value, a
#' scaling by a power of two (exact, since \eqn{S > 2^{-930}}) and
#' comparisons,
#' which are all exact; and
#' the stage contains no division (its operations are four subtractions, two
#' multiplications, one addition, absolute values, comparisons and a scaling
#' by \eqn{2^{-50}}), so a reciprocal option has nothing to transform. An
#' accepted option can still act through the floating-point environment it
#' sets at start-up, which the run-time tests below cover. The
#' proof assumes, beyond these checks, a compilation that performs each
#' operation as written, as the C++ standard prescribes for operations on
#' volatile objects. What remains is the floating-point environment at run time,
#' which some libraries change: the flush-to-zero and denormals-are-zero modes (Clang
#' links such a start-up routine into executables built with
#' \code{-funsafe-math-optimizations}) make subnormal arithmetic inexact and
#' would also change the comparisons of both visibility criteria, and the
#' precision control of the x87 unit can reduce the significand to 24 bits.
#' Before building a graph, the constructors test, with results that are
#' exact in binary64 and therefore independent of the rounding mode, that
#' \eqn{(1 + 2^{-52}) - 1 = 2^{-52}} and
#' \eqn{(1 + 2^{-26})^2 = 1 + 2^{-25} + 2^{-52}} (53 significant bits), that
#' half of the smallest normal double is still positive and that the smallest
#' subnormal double compares greater than zero; they stop with an error when
#' any test fails. They also compute in the non-stop mode of IEEE 754: they
#' save the floating-point environment, disable every trap (an enabled trap on
#' inexact results or on underflow would otherwise interrupt the filter, whose
#' products can be inexact or underflow) and restore the environment when
#' they return. With the compile-time checks, these tests cover the
#' hypotheses of the proof except one case: a run-time precision strictly
#' between 53 and 64 bits would be neither covered by the proof nor detected,
#' and under it the exactness is not established.
#'
#' @section Reproducible checks:
#' In the source repository, \code{dev/check_filter_constants.py} checks the
#' numerical facts used in the floating-point stage (the constants 6.01,
#' 6.02, 7.99 and 1.97, the bound \eqn{e}, the bound on \eqn{|l| + |r|} from
#' \eqn{2^{-930}}, the absolute terms, the overflow bounds and
#' \eqn{2^{-50} = 8u}) in exact rational arithmetic, and
#' \code{dev/make_exact_predicate_fixture.py} regenerates, deterministically,
#' the 1,800 sextuples of the test suite in nine regimes (among them extreme
#' exponents, subnormal numbers, points within two units of rounding of a
#' line, and inputs beyond the bound of the floating-point stage), with signs
#' computed in exact rational arithmetic outside the package. The test suite
#' compares the computed sign with each of them.
#'
#' @references
#'   Shewchuk, J. R. (1997). Adaptive precision floating-point arithmetic and
#'   fast robust geometric predicates. Discrete & Computational Geometry,
#'   18(3), 305-363. https://doi.org/10.1007/PL00009321
#'
#' @seealso [natural_visibility_graph()], [horizontal_visibility_graph()].
#'
#' @name natural_visibility_exactness
NULL
