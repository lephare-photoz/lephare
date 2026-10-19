/*
  17/11/2014
  Class to store one element with lambda, value, origin (used for create vector
  in filter and SED)
*/

// avoid multiple def of the same class
#ifndef ONEEL_H  // check that this keyword has been set already
#define ONEEL_H  // define the keyword to be checked

#include <algorithm>
#include <cassert>
#include <iostream>
#include <tuple>
#include <utility>
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

  /// Check that current lamb < rhs (comparison against a bare wavelength)
  /// @param rhs: the wavelength value to compare to
  inline bool operator<(const double rhs) const { return lamb < rhs; }

  /// Check that current lamb == rhs.lamb
  /// @param rhs: the other object to compare lambda equality to
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

/*! Create a regular grid of double
 * @param lo : start of the grid
 * @param hi : stop of the grid
 * @param dx : interval (throw if non strictly positive)
 */
std::vector<double> make_regular_grid(double lo, double hi, double dx);

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
