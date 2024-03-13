/*

  06/01/14
  functions of the class Schechter function
  1) compute the double schecheter in the case of a MF or LF
  2) write it

*/

#include "schechterFunc.h"

#include <stdlib.h> /* abs */

#include <cmath>
#include <fstream>   // print output file
#include <iostream>  // print stadard file

using namespace std;

// define the operator () which return the value of the double Schechter at a
// given M
double SchechterFunc::operator()(const double M, const double sigconv) {
  double dm4, phi1inter, phi2inter;
  double pi = 4 * atan(1);
  double gaussconv, funcori, funcfin, funcunconv;

  double stepconv;

  stepconv = 6. * sigconv / 20.;  // 30 steps between the -5 sigma and +5 sigma

  // Case were the distance between M and Mstar reasonable
  if (abs(M - fmstar) < 20.) {
    switch (ftype) {
      case 0:  // case LF with Schechter
        funcfin = 0.;
        if (sigconv <= 0.) {
          dm4 = pow(10., -0.4 * (M - fmstar));
          phi1inter = fphi1 * pow(dm4, falpha1);
          phi2inter = fphi2 * pow(dm4, falpha2);
          funcfin = 0.4 * (phi1inter + phi2inter) * log(10.) * exp(-dm4) * dm4;
        } else {
          for (double Mc = M - 3. * sigconv; Mc < M + 3. * sigconv;
               Mc = Mc + stepconv) {
            gaussconv = 1. / (sigconv * sqrt(2. * pi)) *
                        exp(-0.5 * pow((M - Mc) / sigconv, 2.));
            dm4 = pow(10., -0.4 * (Mc - fmstar));
            phi1inter = fphi1 * pow(dm4, falpha1);
            phi2inter = fphi2 * pow(dm4, falpha2);
            funcori =
                0.4 * (phi1inter + phi2inter) * log(10.) * exp(-dm4) * dm4;
            funcfin = funcfin + funcori * gaussconv * stepconv;
          }
        }
        return funcfin;
        break;

      case 1:  // case MF or SFRF with Schechter
        funcfin = 0.;
        if (sigconv <= 0.) {
          dm4 = pow(10., (M - fmstar));
          phi1inter = fphi1 * pow(dm4, falpha1);
          phi2inter = fphi2 * pow(dm4, falpha2);
          funcfin = (phi1inter + phi2inter) * log(10.) * exp(-dm4) * dm4;
        } else {
          for (double Mc = M - 3. * sigconv; Mc < M + 3. * sigconv;
               Mc = Mc + stepconv) {
            gaussconv = 1. / (sigconv * sqrt(2. * pi)) *
                        exp(-0.5 * pow((M - Mc) / sigconv, 2.));
            dm4 = pow(10., (Mc - fmstar));
            phi1inter = fphi1 * pow(dm4, falpha1);
            phi2inter = fphi2 * pow(dm4, falpha2);
            funcori = (phi1inter + phi2inter) * log(10.) * exp(-dm4) * dm4;
            funcfin = funcfin + funcori * gaussconv * stepconv;
          }
        }
        return funcfin;
        break;

      case 2:  // case MF or SFRF with gaussian
        funcfin = 0.;
        if (sigconv <= 0.) {
          funcfin = fphi1 / (falpha1 * sqrt(2. * pi)) *
                    exp(-0.5 * pow((M - fmstar) / falpha1, 2.));
        } else {
          for (double Mc = M - 3. * sigconv; Mc < M + 3. * sigconv;
               Mc = Mc + stepconv) {
            gaussconv = 1. / (sigconv * sqrt(2. * pi)) *
                        exp(-0.5 * pow((M - Mc) / sigconv, 2.));
            funcori = fphi1 / (falpha1 * sqrt(2. * pi)) *
                      exp(-0.5 * pow((Mc - fmstar) / falpha1, 2.));
            funcfin = funcfin + funcori * gaussconv * stepconv;
          }
        }
        return funcfin;
        break;

      case 3:  // case MF or SFRF with double-exponential profile, as Le Floch
               // 05
        funcfin = 0.;
        if (sigconv <= 0.) {
          dm4 = pow(10., (M - fmstar));
          phi1inter =
              pow(dm4, -2.0 - falpha2);  // Change la definition 1-alpha2
                                         // pour que alpha2<alpha1
          funcfin = fphi1 * phi1inter *
                    exp(-0.5 / pow(falpha1, 2.) * pow(log10(1 + dm4), 2.));
        } else {
          for (double Mc = M - 3. * sigconv; Mc < M + 3. * sigconv;
               Mc = Mc + stepconv) {
            gaussconv = 1. / (sigconv * sqrt(2. * pi)) *
                        exp(-0.5 * pow((M - Mc) / sigconv, 2.));
            dm4 = pow(10., (Mc - fmstar));
            phi1inter =
                pow(dm4, -2.0 - falpha2);  // Change la definition 1-alpha2
                                           // pour que alpha2<alpha1
            funcori = fphi1 * phi1inter *
                      exp(-0.5 / pow(falpha1, 2.) * pow(log10(1 + dm4), 2.));
            funcfin = funcfin + funcori * gaussconv * stepconv;
          }
        }
        return funcfin;
        break;

      case 4:  // case sSFRF with gaussian + starburst
        funcfin = 0.;
        if (sigconv <= 0.) {
          funcfin = fphi1 / (falpha1 * sqrt(2. * pi)) *
                        exp(-0.5 * pow((M - fmstar) / falpha1, 2.)) +
                    fphi1 * fphi2 / (0.25 * sqrt(2. * pi)) *
                        exp(-0.5 * pow((M - fmstar - 0.6) / 0.25, 2.));
        } else {
          for (double Mc = M - 3. * sigconv; Mc < M + 3. * sigconv;
               Mc = Mc + stepconv) {
            gaussconv = 1. / (sigconv * sqrt(2. * pi)) *
                        exp(-0.5 * pow((M - Mc) / sigconv, 2.));
            funcori = fphi1 / (falpha1 * sqrt(2. * pi)) *
                          exp(-0.5 * pow((Mc - fmstar) / falpha1, 2.)) +
                      fphi1 * fphi2 / (0.25 * sqrt(2. * pi)) *
                          exp(-0.5 * pow((Mc - fmstar - 0.6) / 0.25, 2.));
            funcfin = funcfin + funcori * gaussconv * stepconv;
          }
        }
        return funcfin;
        break;

      case 5:  // case MF or SFRF with double-exponential profile, as Le Floch
               // 05 + starburst component
        funcfin = 0.;
        if (sigconv <= 0.) {
          dm4 = pow(10., (M - fmstar));
          phi1inter =
              pow(dm4, -2.0 - falpha2);  // Change la definition 1-alpha2
                                         // pour que alpha2<alpha1
          funcfin = fphi1 * phi1inter *
                        exp(-0.5 / pow(falpha1, 2.) * pow(log10(1 + dm4), 2.)) +
                    prenorma * fphi2 / (0.25 * sqrt(2. * pi)) *
                        exp(-0.5 * pow((M - premed - 0.6) / 0.25, 2.));

        } else {
          for (double Mc = M - 3. * sigconv; Mc < M + 3. * sigconv;
               Mc = Mc + stepconv) {
            gaussconv = 1. / (sigconv * sqrt(2. * pi)) *
                        exp(-0.5 * pow((M - Mc) / sigconv, 2.));
            dm4 = pow(10., (Mc - fmstar));
            phi1inter =
                pow(dm4, -2.0 - falpha2);  // Change la definition 1-alpha2
                                           // pour que alpha2<alpha1
            funcori =
                fphi1 * phi1inter *
                    exp(-0.5 / pow(falpha1, 2.) * pow(log10(1 + dm4), 2.)) +
                prenorma * fphi2 / (0.25 * sqrt(2. * pi)) *
                    exp(-0.5 * pow((Mc - premed - 0.6) / 0.25, 2.));
            funcfin = funcfin + funcori * gaussconv * stepconv;
          }
        }
        return funcfin;
        break;

      default:
        cout << "type of fit unknown. Check FIT_TYPE option" << endl;
    }

  } else {
    return (10.e10);
  }
}

// define the function which write the double Schechter function into a file
void SchechterFunc::writeSchech(const string outfile, const double sigconv) {
  ofstream outst;
  double Mmini, Mmaxi, D, D1, D2, D3, D4, D5, lgD1, lgD2, lgD3, lgD4, lgD5;

  // Decide the range to adopt
  switch (ftype) {
    case 0:  // LF
      Mmini = -30.;
      Mmaxi = -5.;
      break;
    case 1:  // MF
      Mmini = 5.;
      Mmaxi = 13.;
      break;
    case 2:  // MF with gauss
      Mmini = 5.;
      Mmaxi = 13.;
      break;
    case 3:  // MF with double expo
      Mmini = 5.;
      Mmaxi = 13.;
      break;
    case 4:  // MF with gauss
      Mmini = 5.;
      Mmaxi = 13.;
      break;
    case 5:  // MF with gauss
      Mmini = 5.;
      Mmaxi = 13.;
      break;
    default:  // SFRF
      Mmini = -7;
      Mmaxi = 5.;
  }

  outst.open(outfile.c_str());
  outst << "# alpha1: " << falpha1 << " err neg: " << errna1
        << " err pos: " << errpa1 << endl;
  outst << "# alpha2: " << falpha2 << " err neg: " << errna2
        << " err pos: " << errpa2 << endl;
  outst << "# mstar: " << fmstar << " err neg: " << errnmstar
        << " err pos: " << errpmstar << endl;
  outst << "# phi1: " << fphi1 / 1000. << " err neg: " << errnp1 / 1000.
        << " err pos: " << errpp1 / 1000. << endl;
  outst << "# phi2: " << fphi2 / 1000. << " err neg: " << errnp2 / 1000.
        << " err pos: " << errpp2 / 1000. << endl;
  // Compute the value of the function at each step and write it
  for (double M = Mmini; M < Mmaxi; M = M + 0.1) {
    // density
    D = this->operator()(M, 0.);        // no convol
    D1 = this->operator()(M, sigconv);  // convol with input
    D2 = this->operator()(M, 0.2);      // extreme conv
    D3 = this->operator()(M, 0.25);     // extreme conv
    D4 = this->operator()(M, 0.3);      // extreme conv
    D5 = this->operator()(M, 0.35);     // extreme conv
    lgD1 = -999.;
    lgD2 = -999.;
    lgD3 = -999.;
    lgD4 = -999.;
    lgD5 = -999.;
    // check if the log makes sens before writing it
    if (D1 > 1.e-100) {
      lgD1 = log10(D1);
    }
    if (D2 > 1.e-100) {
      lgD2 = log10(D2);
    }
    if (D3 > 1.e-100) {
      lgD3 = log10(D3);
    }
    if (D4 > 1.e-100) {
      lgD4 = log10(D4);
    }
    if (D5 > 1.e-100) {
      lgD5 = log10(D5);
    }
    if (D > 1.e-100) {
      outst << M << " " << log10(D) - 3 << " " << lgD1 - 3 << " " << lgD2 - 3
            << " " << lgD3 - 3 << " " << lgD4 - 3 << " " << lgD5 - 3 << endl;
    }
  }
  outst.close();
}

// define the function which compute the integral of the Schechter function
double SchechterFunc::schechInt(const double Mmin, const double Mmax) {
  double Da, Db, Dab;
  double step = (Mmax - Mmin) / 10000.;
  double inte = 0.;

  // Compute the numerical integral of the schechter function between Mmin-Mmax
  for (double M = Mmin; M < Mmax; M = M + step) {
    // density
    Da = this->operator()(M, 0.);
    Db = this->operator()(M + step, 0.);
    Dab = this->operator()(M + step / 2., 0.);
    // sum
    if (Dab > 0.) {
      inte = inte + (Da + 4 * Dab + Db) * step / 6.;
    }
  }
  return (inte);
}

// define the function which compute the LD: integral of the Schechter function
// * L
double SchechterFunc::schechIntLD(const double Mmin, const double Mmax) {
  double Da, Db, Dab;
  double step = (Mmax - Mmin) / 10000.;
  double inte = 0.;

  // Compute the numerical integral of the schechter function between Mmin-Mmax
  for (double M = Mmin; M < Mmax; M = M + step) {
    // density
    Da = this->operator()(M, 0.);
    Db = this->operator()(M + step, 0.);
    Dab = this->operator()(M + step / 2., 0.);
    // sum
    if (Dab > 0.) {
      inte = inte + pow(10., M + step / 2.) * (Da + 4 * Dab + Db) * step / 6.;
    }
  }
  return (inte);
}

// define the SFR median
double SchechterFunc::Mmed(const double Mmin, const double Mmax) {
  double integral, Da, Db, Dab, M;
  double step = (Mmax - Mmin) / 10000.;
  double inte = 0.;

  integral = this->schechInt(Mmin, Mmax);

  // Compute the numerical integral of the schechter function between Mmin-Mmax
  M = Mmin;
  while (inte < (integral / 2.) && M < Mmax) {
    // density
    Da = this->operator()(M, 0.);
    Db = this->operator()(M + step, 0.);
    Dab = this->operator()(M + step / 2., 0.);
    // sum
    if (Dab > 0) {
      inte = inte + (Da + 4 * Dab + Db) * step / 6.;
    }
    // Increment
    M = M + step;
  }
  return (M - step / 2.);
}
