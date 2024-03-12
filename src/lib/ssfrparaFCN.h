/*

  28/05/2014
  FCN function that we want to minimize for minuit to fit the evolution of the sSFR

*/

#ifndef MN_SSFR_FCN_H_
#define MN_SSFR_FCN_H_

#include <sstream> // print output file
#include <fstream> // print output file
#include <iostream> // print standard file

#include "Minuit2/FCNBase.h"
#include <vector>

#include "LFdata.h"    

using namespace std;
using namespace ROOT::Minuit2;

class ssfrparaFCN : public FCNBase {

public:


  //constructor
  ssfrparaFCN(const std::vector<LFdata>& ssfrData)    
   { fssfrData=ssfrData;
     fErrorDef=1.;
   }

  ~ssfrparaFCN() {}
 
  virtual double Up() const {return fErrorDef;}
  virtual double operator()(const vector<double>&) const;
   
  void SetErrorDef(double def) {fErrorDef = def;}

  

private:

  vector<LFdata>   fssfrData;
 double fErrorDef;

};

#endif // MN_SSFR_FCN_H_
