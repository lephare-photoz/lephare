/*
  Absolute-magnitude prior: log10-free model rejection in onesource::fit().

  In onesource::fit() every call evaluates, for each of the ~1e6 library
  models, the rest-frame absolute magnitude
      abs_mag = mag0 - 2.5*log10(dm)   [- funz0 if z < 1e-10]
  and rejects the model if abs_mag <= lo or abs_mag >= hi.
  abs_mag is strictly decreasing in dm, so the test is equivalent to
      dm >= thrLo = 10^((mag0' - lo)/2.5)   or   dm <= thrHi = 10^((mag0' - hi)/2.5)
  where mag0' = mag0 (- funz0).  The thresholds only depend on the library and
  on the prior, not on the object: they are computed once per library.

  To keep the decisions *identical* to the original expression, the
  threshold test is only trusted when dm is farther than a relative guard band
  (kGuard, ~1e9 times larger than the rounding error of log10/pow) from a
  threshold; otherwise (and for dm <= 0, NaN, inf, e.g. the -999 sentinel, or
  non finite thresholds) the original log10 expression is evaluated.
*/
#ifndef ABSMAG_PRIOR_H
#define ABSMAG_PRIOR_H

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace absmag {

/// relative guard band around the thresholds (rounding errors are ~1e-14)
constexpr double kGuard = 1.0e-9;

/// Reference (original) abs-mag computation, kept verbatim.
inline double exact_abs_mag(double mag0, double dm, bool z_is_zero,
                            double funz0) {
  double abs_mag = mag0 - 2.5 * log10(dm);
  if (z_is_zero) abs_mag = abs_mag - funz0;
  return abs_mag;
}

/// Original rejection test (reference implementation).
inline bool rejected_exact(double dm, double mag0, bool z_is_zero,
                           double funz0, double lo, double hi) {
  const double a = exact_abs_mag(mag0, dm, z_is_zero, funz0);
  return (a <= lo || a >= hi);
}

/// Fast rejection test; returns exactly rejected_exact(...)
inline bool rejected(double dm, double mag0, bool z_is_zero, double funz0,
                     double lo, double hi, double thr_lo, double thr_hi) {
  // finite, strictly positive dm only; NaN thresholds make every test false
  if (dm > 0 && dm < std::numeric_limits<double>::infinity()) {
    if (dm >= thr_lo * (1. + kGuard) || dm <= thr_hi * (1. - kGuard))
      return true;
    if (dm < thr_lo * (1. - kGuard) && dm > thr_hi * (1. + kGuard))
      return false;
  }
  return rejected_exact(dm, mag0, z_is_zero, funz0, lo, hi);
}

/// Per-library thresholds, rebuilt only when library or prior change.
struct Cache {
  std::vector<double> thr_lo, thr_hi;  ///< NaN where no prior applies / invalid
  size_t n = 0;
  double prior[4] = {0, 0, 0, 0};
  double funz0 = 0;
  const double *mag0_ptr = nullptr, *red_ptr = nullptr;
  double sig = 0;  ///< cheap sampled signature of mag0/red/nlib

  template <class V, class R, class T>
  static double signature(const V& mag0, const R& red, const T& nlib) {
    double s = 0;
    const size_t n = mag0.size();
    const size_t step = n / 61 + 1;
    for (size_t i = 0; i < n; i += step)
      s = s * 1.000001 + mag0[i] + 3. * red[i] + static_cast<int>(nlib[i]);
    return s + (n ? mag0[n - 1] + red[n - 1] : 0.);
  }

  /// (Re)build if needed.  Must be called outside any parallel region.
  /// nlib: 0 = GAL (priors [0],[1]), 1 = QSO (priors [2],[3]).
  template <class V, class R, class T, class P>
  void ensure(const T& nlib, const V& mag0, const R& red, const P& priorLib,
              double funz0_in) {
    const size_t N = mag0.size();
    const double sg = signature(mag0, red, nlib);
    if (N == n && thr_lo.size() == N && mag0.data() == mag0_ptr &&
        red.data() == red_ptr && funz0_in == funz0 && sg == sig &&
        prior[0] == priorLib[0] && prior[1] == priorLib[1] &&
        prior[2] == priorLib[2] && prior[3] == priorLib[3])
      return;
    n = N;
    mag0_ptr = mag0.data();
    red_ptr = red.data();
    funz0 = funz0_in;
    sig = sg;
    for (int k = 0; k < 4; k++) prior[k] = priorLib[k];
    const double nan = std::numeric_limits<double>::quiet_NaN();
    thr_lo.assign(N, nan);
    thr_hi.assign(N, nan);
    const bool gal = (priorLib[0] < 0 && priorLib[1] < 0);
    const bool agn = (priorLib[2] < 0 && priorLib[3] < 0);
    for (size_t i = 0; i < N; i++) {
      const int t = static_cast<int>(nlib[i]);
      double lo, hi;
      if (gal && t == 0) {
        lo = priorLib[0];
        hi = priorLib[1];
      } else if (agn && t == 1) {
        lo = priorLib[2];
        hi = priorLib[3];
      } else {
        continue;
      }
      double m = mag0[i];
      if (red[i] < 1.e-10) m = m - funz0_in;
      const double a = std::pow(10., (m - lo) / 2.5);
      const double b = std::pow(10., (m - hi) / 2.5);
      // keep only usable (finite, > 0) thresholds, else stay NaN => exact path
      if (std::isfinite(a) && a > 0 && std::isfinite(b) && b > 0) {
        thr_lo[i] = a;
        thr_hi[i] = b;
      }
    }
  }
};

}  // namespace absmag

#endif  // ABSMAG_PRIOR_H
