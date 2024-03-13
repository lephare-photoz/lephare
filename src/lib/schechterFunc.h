/*

  03/01/14
  class with the Schechter function

*/

#ifndef MN_SchechterFunc_H_
#define MN_SchechterFunc_H_

#include <math.h>
#include <vector>
#include <fstream> // print output file
#include <iostream> // print stadard file

using namespace std;


class SchechterFunc {

private:
  int ftype;
  double falpha1, falpha2,fmstar,fphi1,fphi2;
  double errna1,errna2,errnmstar,errnp1,errnp2;
  double errpa1,errpa2,errpmstar,errpp1,errpp2;

public:

  double prenorma,premed;
  
  //Constructor of the schechter function, 
  // arguments -> type(LF/MF/SFRF), 5 schechter parameters
  SchechterFunc(const int type, const double alpha1, const double alpha2, const double mstar, const double phi1, const double phi2) : 
    ftype(type), falpha1(alpha1), falpha2(alpha2),fmstar(mstar),fphi1(phi1),fphi2(phi2),
    errna1(-999.), errna2(-999.),errnmstar(-999.),errnp1(-999999.), errnp2(-999999.),
      errpa1(-999.), errpa2(-999.),errpmstar(-999.),errpp1(-999999.), errpp2(-999999.), prenorma(0.), premed(-999) {}

  //Constructor of the schechter function, 
  //arguments -> type(LF/MF/SFRF), a vector of double including the 5 schechter parameters
  SchechterFunc(const int type,  const vector<double> para) : 
    ftype(type), falpha1(para[0]), falpha2(para[1]),fmstar(para[2]),fphi1(para[3]),fphi2(para[4]),
    errna1(-999.), errna2(-999.),errnmstar(-999.),errnp1(-999999.), errnp2(-999999.),
    errpa1(-999.), errpa2(-999.),errpmstar(-999.),errpp1(-999999.), errpp2(-999999.), prenorma(0.), premed(-999)  {}

  //Constructor of the schechter function, 
  //arguments -> type(LF/MF/SFRF), 3 vectors of double, including the 5 schechter parameters, the negative and positive errors associated
  SchechterFunc(const int type,  const  vector<double> para, const vector<double> errn, const vector<double> errp) : 
    ftype(type), falpha1(para[0]), falpha2(para[1]),fmstar(para[2]),fphi1(para[3]),fphi2(para[4]),
    errna1(errn[0]), errna2(errn[1]),errnmstar(errn[2]),errnp1(errn[3]), errnp2(errn[4]),
    errpa1(errp[0]), errpa2(errp[1]),errpmstar(errp[2]),errpp1(errp[3]), errpp2(errp[4]), prenorma(0.), premed(-999) {}

  //Constructor of the schechter function, 
  //arguments -> type(LF/MF/SFRF), vector of double[3], including the 5 schechter parameters, the negative and positive errors associated
  SchechterFunc(const int type,  const  vector< vector <double> > para) : 
    ftype(type), 
    falpha1(para[0][0]), falpha2(para[1][0]),   fmstar(para[2][0]), fphi1(para[3][0]),  fphi2(para[4][0]),
     errna1(para[0][1]),  errna2(para[1][1]),errnmstar(para[2][1]),errnp1(para[3][1]), errnp2(para[4][1]),
     errpa1(para[0][2]),  errpa2(para[1][2]),errpmstar(para[2][2]),errpp1(para[3][2]), errpp2(para[4][2]), prenorma(0.), premed(-999) {}

 ~SchechterFunc() {}


  double a1() const {return falpha1;}
  double a2() const {return falpha2;}
  double mstar() const {return fmstar;}
  double p1() const {return fphi1;}
  double p2() const {return fphi2;}
  double chtype(const int intype) {ftype=intype; return 0;}

  double operator()(const double M, const double sigconv);
  void writeSchech(const string outfile, const double sigconv);
  double schechInt(const double Mmin, const double Mmax);
  double schechIntLD(const double Mmin, const double Mmax);
  double Mmed(const double Mmin, const double Mmax);


};

#endif // MN_SchechterFunc_H_
