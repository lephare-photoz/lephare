/*
02/01/2014
Class to store the LF non parametric data
*/

// avoid multiple def of the same class
#ifndef LF_H  // check that this keyword has been set already
#define LF_H  // define the keyword to be checked

class LFdata {
 private:
  double zmoy;

 public:
  double M, value, errp, errm, z;
  double nor, norerr;
  int bin;

  // constructor of the LFdata class
  // overload constructor, without the redshift bin
  LFdata() {
    M = -999.;
    value = -999.;
    errp = -999.;
    errm = -999.;
    z = -999.;
    bin = -10;
    nor = -999.;
    norerr = -999.;
  }

  // keep it in the same file, since really simple
  LFdata(double zmoy, int k, double norma, double normaerr) {
    M = -999.;
    value = -999.;
    errp = -999.;
    errm = -999.;
    z = zmoy;
    bin = k;
    nor = norma;
    norerr = normaerr;
  }
};

#endif
