/*

  29/05/14
  Implementation of the operator () for the FCN minuit class

*/

#include "ssfrparaFCN.h"

#include <fstream>   // print output file
#include <iostream>  // print standard file
#include <sstream>   // print output file

#include "LFdata.h"
#include "math.h"

using namespace std;
using namespace ROOT::Minuit2;

double ssfrparaFCN::operator()(const vector<double> &par) const {
  double chi2 = 0;

  // Loop over the data to be fitted
  for (unsigned i = 0; i < fssfrData.size(); ++i) {
    // Keep only z>0.2
    if ((fssfrData[i]).z > 0.2 && (fssfrData[i]).value > 0) {
      // model
      double model = par[0] +
                     par[1] * pow(10., (fssfrData[i]).M) / pow(10., 10.5) +
                     par[2] * log10(1 + (fssfrData[i]).z);
      // try with a fit in log(ssfr)
      // double model = par[0] + par[1]*(fssfrData[i]).M +
      // par[2]*log10(1+(fssfrData[i]).z);
      double diff = pow(10., model) - pow(10., (fssfrData[i]).value);
      // error logerr=log(1+err/v)
      double erreur =
          (pow(10., (fssfrData[i]).errp) - 1) * pow(10., (fssfrData[i]).value);
      chi2 += pow(diff / erreur, 2.);
    }
  }

  return chi2;
}
