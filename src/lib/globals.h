#ifndef GLOBAL_H
#define GLOBAL_H

#include <bitset>
#include <cmath>
#include <string>
#include <vector>

#define MAX_CONTEXT 1024
#define INVALID_FLUX -9999.0
#define INVALID_VAL -9999.0
#define INVALID_MAG -9999.0
#define INVALID_PHYS -999.0
#define INVALID_Z -99.9
#define INVALID_INDEX -99
#define NULL_FLUX 0.0
#define HIGH_MAG 1000.0
#define HIGH_CHI2 1.0e9
#define EPS_Z 1.0e-10

using namespace std;

// extern tells the compiler this variable is declared elsewhere
extern string lepharedir, lepharework;
extern string fakeString;
extern const string nonestring;
extern int fakeInt;
extern const double c, ckms;
extern const double hplanck;
extern const double kboltzmann;
extern const double Lsol;
extern const double pc;
extern const double pi;
extern const double hc;
extern const double f_ga;

void get_lephare_env();

bool check_first_char(const string &maligne);

double blackbody(double T, double lambda);

int bdincl(int n, long cont, int max);

//! Return the vector of indexes of values in `vec` that match `value` to the
//! required `precision`.
/*!
  \param value input to be compared to the content of `vector`
  \param vec the input vector of values to be compared to `value`
  \param precision the precision with which `value` must compare to the entries
  in `vector`

  \return the vector of indexes in `vector` of the matching values.
*/
vector<size_t> indexes_in_vec(const double &value, const vector<double> &vec,
                              const float &precision);

/*! For a curve defined by vectors (x,y), return the values interpolated
 * at points z
 * @param x : sorted x vector of the curve
 * @param y : y vector of the curve
 * @param z : sorted vector of values for which to get the interpolated values
 * from (x,y) curve, if interpolation is possible. Else, return d
 * @param d : default value in case of extrapolation
 * This implementation is single pass : O(n + m) time, where n = x.size() and m
 * = z.size(), instead of O(m*log(n)), owing to the assumption that x and z are
 * sorted.
 */
vector<double> fast_interpolate(const std::vector<double> &x,
                                const std::vector<double> &y,
                                const std::vector<double> &z, double d);

inline string bool2string(const bool &b) {
  string sb;
  b ? sb = "YES" : sb = "NO";
  return sb;
}

inline bool CHECK_CONTEXT_BIT(unsigned long context, unsigned int n) {
  return std::bitset<MAX_CONTEXT>(context).test(n);
}
inline double POW10D_FAST(double x) { return exp(2.302585092994046 * x); }
inline double POW10D_SLOW(double x) { return pow(10.0, x); }

// Marginal improvement in speed
// In [3]: %timeit lp.POW10DV(np.random.random(100000000))
// 6.38 s ± 24.7 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
//
// In [4]: %timeit lp.POW10DSLOWV(np.random.random(100000000))
// 6.94 s ± 15.8 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

inline vector<double> POW10D_FASTV(vector<double> x) {
  vector<double> res;
  res.reserve(x.size());
  for (auto xx : x) res.push_back(exp(2.302585092994046 * xx));
  return res;
}
inline vector<double> POW10D_SLOWV(vector<double> x) {
  vector<double> res;
  res.reserve(x.size());
  for (auto xx : x) res.push_back(pow(10.0, xx));
  return res;
}
// This is a fast approximation to log2()
// Y = C[0]*F*F*F + C[1]*F*F + C[2]*F + C[3] + E;
inline double log2f_approx(double X) {
  double Y, F;
  int E;
  F = frexp(fabs(X), &E);
  Y = 1.23149591368684;
  Y *= F;
  Y += -4.11852516267426;
  Y *= F;
  Y += 6.02197014179219;
  Y *= F;
  Y += -3.13396450166353;
  Y += E;
  return (Y);
}

inline double LOG10D_SLOW(double x) { return log10(x); }
inline double LOG10D_FAST(double x) {
  return log2f_approx(x) * 0.3010299956639812f;
}

// Same : marginal; ~0.5s for 10M evaluations
// In [6]: %timeit lp.LOG10D_SLOWV(np.random.random(100000000))
// 6.67 s ± 3.43 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
//
// In [7]: %timeit lp.LOG10D_FASTV(np.random.random(100000000))
// 6.2 s ± 21.5 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
inline vector<double> LOG10D_SLOWV(vector<double> x) {
  vector<double> res;
  res.reserve(x.size());
  for (auto xx : x) res.push_back(log10(xx));
  return res;
}
inline vector<double> LOG10D_FASTV(vector<double> x) {
  vector<double> res;
  res.reserve(x.size());
  for (auto xx : x) res.push_back(log2f_approx(xx) * 0.3010299956639812f);
  return res;
}

#define POW10D POW10D_SLOW
#define LOG10D LOG10D_SLOW

inline double mag2flux(double x, double zp = 0.0) {
  return POW10D(-0.4 * (x + 48.6 - zp));
}
inline double flux2mag(double x, double offset = 0.0) {
  return -2.5 * LOG10D(x) - 48.6 + offset;
}
#endif
