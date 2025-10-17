/*
  17/11/2014
  Class to store one element with lambda, value, origin (used for create vector
  in filter and SED)
*/

// avoid multiple def of the same class
#ifndef ONEEL_H  // check that this keyword has been set already
#define ONEEL_H  // define the keyword to be checked

/// Small class structure to store one extinction element, including lambda,
/// value, and origin
class oneElLambda {
 public:
  double lamb,  ///< lambda value
      val,      ///< extinction value at lambda
      ori;  ///< index : index referring to filter transmission (ori=0), or SED
            ///< flux (ori=1), or attenuation law (ori=2), or opacity (ori=3),
            ///< or emission line (ori=5)

  /*! Constructor
   * @param lambda: the \f$\lambda\f$ value
   * @param value: the value at `lambda`
   * @param origin: an arbitrary integer that flags the provenance
   * of the oneElLambda object, when part of a vector (used in SED::resample
   * for instance
   !*/
  oneElLambda(double lambda, double value, int origin) {
    lamb = lambda;
    val = value;
    ori = origin;
  }

  /// Copy constructor
  /// @param obj: the object to copy from
  oneElLambda(const oneElLambda &obj) {
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
  void interp(const oneElLambda &previous, const oneElLambda &next);

  /// Check that current lamb < rhs.lamb
  /// @param rhs: the other object to compare lambda ordering to
  inline bool operator<(const oneElLambda &rhs) const {
    return lamb < rhs.lamb;
  }
};

#endif
