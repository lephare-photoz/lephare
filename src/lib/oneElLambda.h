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
  short ori;    ///< index : index referring to filter transmission (ori=0), or
                ///< SED flux (ori=1), or attenuation law (ori=2), or opacity
                ///< (ori=3), or emission line (ori=5)

  /*! Constructor
   * @param lambda: the \f$\lambda\f$ value
   * @param value: the value at `lambda`
   * @param origin: an arbitrary integer that flags the provenance
   * of the oneElLambda object, when part of a vector (used in SED::resample
   * for instance
   !*/
  oneElLambda(double lambda, double value, short origin) {
    lamb = lambda;
    val = value;
    ori = origin;
  }

  /// Copy constructor
  /// @param obj: the object to copy from
  oneElLambda(const oneElLambda& obj) {
    lamb = obj.lamb;
    val = obj.val;
    ori = obj.ori;
  }

  /// empty destructor
  ~oneElLambda() {}

  /*! interpolate linearly between two adjacent lambda values
   * (expected positive and ordered)
   * @param previous: oneEleLambda before the current lambda
   * @param next: oneElLambda after current lambda
   * @return If correctly ordered, set current value to the linear
   * interpolation at lambda between `previous` and `next`
   !*/
  void interp(const oneElLambda& previous, const oneElLambda& next);

  /// Check that current lamb < rhs.lamb
  /// @param rhs: the other object to compare lambda ordering to
  inline bool operator<(const oneElLambda& rhs) const {
    return lamb < rhs.lamb;
  }

  inline bool operator<(const double rhs) const { return lamb < rhs; }
};

typedef vector<oneElLambda> oneElVector;

inline pair<vector<double>, vector<double>> to_pairs(
    const vector<oneElLambda>& v) {
  vector<double> lambdas;
  vector<double> vals;
  for (auto it = v.begin(); it != v.end(); it++) {
    lambdas.push_back(it->lamb);
    vals.push_back(it->val);
  }
  return make_pair(lambdas, vals);
}

inline std::tuple<std::vector<double>, std::vector<double>> to_tuple(
    const oneElVector& v) {
  auto p = to_pairs(v);
  return {p.first, p.second};
}

double interp_linear_point(const std::vector<double>& x,
                           const std::vector<double>& y, double xi);

std::vector<double> interp_linear_vec(const std::vector<double>& x,
                                      const std::vector<double>& y,
                                      const std::vector<double>& q);

std::vector<double> make_regular_grid(double lo, double hi, double dx);

std::vector<double> make_union_grid(const std::vector<double>& x1,
                                    const std::vector<double>& x2, double lo,
                                    double hi);

/*! Provide the cross interpolation of (x1, y1) and (x2, y2) in the
 * intersection of x1 and x2.
 * @param x1 : vector x of first function
 * @param y1 : vector y of first function
 * @param x2 : vector x of second function
 * @param y2 : vector y of second function
 * @output : x, y1', y2' the outputs of the restricted cross interpolation;
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
 * @output x, y1, y2 : x is the restricted union of the x vectors of v1 and v2
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
