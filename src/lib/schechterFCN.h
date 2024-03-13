/*

  08/01/2014
  FCN function that we want to minimize for minuit
  This function is a chi2 with the difference between the LFdata and a double schechter function

*/

#ifndef MN_SCHECHTER_FCN_H_
#define MN_SCHECHTER_FCN_H_

#include <sstream> // print output file
#include <fstream> // print output file
#include <iostream> // print standard file

#include "Minuit2/FCNBase.h"
#include <vector>

using namespace std;
using namespace ROOT::Minuit2;

class schechterFCN : public FCNBase {

public:

  //constructor
 schechterFCN(const std::vector<double>& meas, // density measurements
	      const std::vector<double>& pos,  // M positions
	      const std::vector<double>& mvar, // error on the  measurements  
	      const std::vector<int>& lowlim,  // Is the point a lower limit lowlim=1
              const int type,                  //type of schecheter function to use (0:LF,1:MF,2SFRF)
              const double normal,             //Use the normalisation to constrain the fit 
              const double normalerr,          //and its associated error
              const double sigconv)            //sigma of the gaussian used for the convolution
   {   fMeasurements=meas;
       fPositions=pos;
       fMVariances=mvar;
       fLim=lowlim;
       fErrorDef=1.;
       ftype=type ;
       fnormal=normal;
       fnormalerr=normalerr;
       fsigconv=sigconv;
  }

  ~schechterFCN() {}

  virtual double Up() const {return fErrorDef;}
  virtual double operator()(const vector<double>&) const;
  
  vector<double> Measurements() const {return fMeasurements;}
  vector<double> Positions() const {return fPositions;}
  vector<double> Variances() const {return fMVariances;}

  void SetErrorDef(double def) {fErrorDef = def;}

private:

  vector<double> fMeasurements, fPositions, fMVariances;
  vector<int> fLim;
  double fErrorDef,fnormal,fnormalerr,fsigconv;
  int ftype;

};

#endif // MN_SCHECHTER_FCN_H_
