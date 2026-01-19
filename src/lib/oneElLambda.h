/*
  17/11/2014
  Class to store one element with lambda, value, origin (used for create vector
  in filter and SED)
*/

// avoid multiple def of the same class
#ifndef ONEEL_H  // check that this keyword has been set already
#define ONEEL_H  // define the keyword to be checked

#include <algorithm>
#include <tuple>
#include <vector>

using namespace std;

/// Small class structure to store one extinction element, including lambda,
/// value, and origin
class oneElLambda {
 public:
  double lamb,  ///< lambda value
      val;      ///< extinction value at lambda

  /*! Constructor
   * @param lambda: the \f$\lambda\f$ value
   * @param value: the value at `lambda`
   !*/
  oneElLambda(double lambda, double value) {
    lamb = lambda;
    val = value;
  }

  /// Copy constructor
  /// @param obj: the object to copy from
  oneElLambda(const oneElLambda& obj) {
    lamb = obj.lamb;
    val = obj.val;
  }

  /// empty destructor
  ~oneElLambda() {}

  /// Check that current lamb < rhs.lamb
  /// @param rhs: the other object to compare lambda ordering to
  inline bool operator<(const oneElLambda& rhs) const {
    return lamb < rhs.lamb;
  }

  inline bool operator<(const double rhs) const { return lamb < rhs; }

  inline bool operator==(const oneElLambda& rhs) const {
    return lamb == rhs.lamb;
  }
};

typedef vector<oneElLambda> oneElVector;
inline std::tuple<std::vector<double>, std::vector<double>> to_tuple(
    const oneElVector& v) {
  std::vector<double> lambdas;
  std::vector<double> vals;
  lambdas.reserve(v.size());
  vals.reserve(v.size());

  for (const auto& el : v) {
    lambdas.push_back(el.lamb);
    vals.push_back(el.val);
  }

  return {std::move(lambdas), std::move(vals)};
}

/*! Interpolate a curve (x,y) at positions given by a vector
 * @param x : x-values of the curve to interpolate
 * @param y : y-values of the curve to interpolate
 * @param xi : x-value at which to interpolate
 * @return: the y-value from (x,y) interpolated at xi.
 * Important : x is expected to be sorted in increasing order.
 * if xi is outside the range of x, 0 is returned.
 */
double interp_linear_point(const std::vector<double>& x,
                           const std::vector<double>& y, double xi);

/*! Interpolate a curve (x,y) at positions given by a vector
 * @param x : x-values of the curve to interpolate
 * @param y : y-values of the curve to interpolate
 * @param q : x-values at which to interpolate
 * @return: the vector of y-values frm (x,y) interpolated at q.
 * Important : x is expected to be sorted in increasing order.
 */
inline std::vector<double> interp_linear_vec(const std::vector<double>& x,
                                             const std::vector<double>& y,
                                             const std::vector<double>& q) {
  std::vector<double> out(q.size());
  // not worth accelerating, due to OMP overhead
  // #pragma omp parallel for
  for (long long i = 0; i < (long long)q.size(); ++i)
    out[i] = interp_linear_point(x, y, q[i]);
  return out;
}

/*! Create a regular grid of double
 * @param lo : start of the grid
 * @param hi : stop of the grid
 * @param dx : interval (throw if non strictly positive)
 */
std::vector<double> make_regular_grid(double lo, double hi, double dx);

/*! Create an irregular grid based on the unions of two vectors of double
 * @param x1 : first vector
 * @param x2 : second vector
 * @param lo : start of the grid
 * @param hi : stop of the grid
 */
std::vector<double> make_union_grid(const std::vector<double>& x1,
                                    const std::vector<double>& x2, double lo,
                                    double hi);

/*! Provide the cross interpolation of (x1, y1) and (x2, y2) in the
 * intersection of x1 and x2.
 * @param x1 : vector x of first function
 * @param y1 : vector y of first function
 * @param x2 : vector x of second function
 * @param y2 : vector y of second function
 * @param dx : grid interval if positive, else indicates that the grid is
 * defined by the union of x1 and x2
 * @return x, y1', y2' : the outputs of the restricted cross interpolation;
 * x is the common vector of values from x1 and x2, and y1' and y2' are
 * respectively the interpolated values of  (x1, y1) and (x2, y2) over x.
 *
 */
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
common_interpolate_combined(const std::vector<double>& x1,
                            const std::vector<double>& y1,
                            const std::vector<double>& x2,
                            const std::vector<double>& y2, double dx);

/*! Helper function to take two vectors of oneElLambda objects and to
 * return the cross interpolation of them restricted to their intersection
 * @param v1 : first vector of oneElLambda
 * @param v2 :  second vector of oneElLambda
 * @return x, y1, y2 : x is the restricted union of the x vectors of v1 and v2
 * and y1, y2 are respectively the interpolated vector of v1 and v2 values at x.
 */
inline std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
restricted_resampling(const oneElVector& v1, const oneElVector& v2, double dx) {
  auto [v1_l, v1_v] = to_tuple(v1);
  auto [v2_l, v2_v] = to_tuple(v2);
  return common_interpolate_combined(v1_l, v1_v, v2_l, v2_v, dx);
}

/*! merge two oneElLambda vectors into a sorted one
 *
 * @param v1: first vector of oneElLambda
 * @param v2: second vector of oneElLambda
 *
 * @return vector merging v1 and v2 into a single vector,
 * ordered by increasing lambda. Each oneElLambda keeps its
 * original attribute `ori`, which allows to keep track of the
 * origin of this element as coming from v1 or v2
 */
inline vector<oneElLambda> concatenate_and_sort(const vector<oneElLambda>& v1,
                                                const vector<oneElLambda>& v2) {
  vector<oneElLambda> res = v1;
  res.insert(res.end(), v2.begin(), v2.end());
  sort(res.begin(), res.end());
  return res;
}

#endif
