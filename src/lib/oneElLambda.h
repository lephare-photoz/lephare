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

  // Constructor, passing the values as double
  oneElLambda(double lambin, double valin, int oriin) {
    lamb = lambin;
    val = valin;
    ori = oriin;
  }

  // Copy constructor
  oneElLambda(const oneElLambda &elIn) {
    lamb = elIn.lamb;
    val = elIn.val;
    ori = elIn.ori;
  }

  /* MP: added destructor*/
  ~oneElLambda() {}

  /// interpolate linearly between two adjacent lambda values (expected positive
  /// and ordered)
  void interp(const oneElLambda &previousEl, const oneElLambda &nextEl);

  /// comparison operator for sorting
  inline bool operator<(const oneElLambda &rhs) const {
    return lamb < rhs.lamb;
  }
};

#endif
