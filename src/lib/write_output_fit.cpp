/*
original: 09/05/2014
updated:

Write the output (latex, parameters, screen)

input: stream from the outputs, parameters, fittype
output: the files
*/
#include <fstream>
#include <iomanip>  // std::setprecision
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "schechterFunc.h"

using namespace std;

void write_output_fit(ofstream &outlatex, ofstream &outpara, const string out,
                      const int fittype, const vector<vector<double>> allPara,
                      const double zmins, const double zmaxs, const int ngals,
                      const double sigconv) {
  // Write the output file for each function
  double prenorma = 0;
  double premed = -99;
  if (fittype == 4 || fittype == 5) {
    SchechterFunc schechpre(3, allPara);
    if (fittype == 4) {
      schechpre.chtype(2);
    }
    prenorma = schechpre.schechInt(0., 20.);
    premed =
        log10(schechpre.schechIntLD(0., 20.) / schechpre.schechInt(0., 20.));
  }
  SchechterFunc finalSchech(fittype, allPara);
  if (fittype == 4 || fittype == 5) {
    finalSchech.prenorma = prenorma;
    finalSchech.premed = premed;
  }
  finalSchech.writeSchech(out, sigconv);

  // OUTPUT SCREEN

  cout << endl
       << "Best fit parameters for the redshift bin:" << zmins << " " << zmaxs
       << endl;
  cout << "alpha 1:   " << allPara[0][0] << "  " << allPara[0][1] << "  "
       << allPara[0][2] << endl;
  cout << "alpha 2:   " << allPara[1][0] << "  " << allPara[1][1] << "  "
       << allPara[1][2] << endl;
  cout << "Mstar:     " << allPara[2][0] << "  " << allPara[2][1] << "  "
       << allPara[2][2] << endl;
  cout << "phistar 1: " << allPara[3][0] << "  " << allPara[3][1] << "  "
       << allPara[3][2] << endl;
  cout << "phistar 2: " << allPara[4][0] << "  " << allPara[4][1] << "  "
       << allPara[4][2] << endl;

  // OUTPUT LATEX
  outlatex << std::fixed;
  outlatex << " & " << setprecision(1) << zmins << "-" << setprecision(1)
           << zmaxs << " & " << ngals << " & ";

  // alpha
  outlatex << setprecision(3) << allPara[0][0] << "$^{{\\rm +";
  outlatex << setprecision(3) << allPara[0][2] << "}}_{{\\rm ";
  outlatex << setprecision(3) << allPara[0][1] << "}}$ & ";

  // In the case double-exponential (sigma) or schechter (alpha2)
  if (fittype == 0 || fittype == 1 || fittype == 3) {
    if (fittype != 3) {
      outlatex << setprecision(3) << allPara[1][0] << "$^{{\\rm +";
    } else {
      // Case double-exponential, to fall back on the right definition
      outlatex << setprecision(3) << allPara[1][0] + 3 << "$^{{\\rm +";
    }
    outlatex << setprecision(3) << allPara[1][2] << "}}_{{\\rm ";
    outlatex << setprecision(3) << allPara[1][1] << "}}$ & ";
  }

  // Mstar
  outlatex << setprecision(3) << allPara[2][0] - 11 << "$^{{\\rm +";
  outlatex << setprecision(3) << allPara[2][2] << "}}_{{\\rm ";
  outlatex << setprecision(3) << allPara[2][1] << "}}$ & ";

  // Phistar 1
  outlatex << setprecision(3) << allPara[3][0] << "$^{{\\rm +";
  outlatex << setprecision(3) << allPara[3][2] << "}}_{{\\rm ";
  outlatex << setprecision(3) << allPara[3][1] << "}}$ & ";

  // In the case schechter (phistar2)
  if (fittype == 0 || fittype == 1) {
    outlatex << setprecision(3) << allPara[4][0] << "$^{{\\rm +";
    outlatex << setprecision(3) << allPara[4][2] << "}}_{{\\rm ";
    outlatex << setprecision(3) << allPara[4][1] << "}}$ & ";
  }

  // Median sSFR
  outlatex << setprecision(3) << allPara[5][0] - 11 << "$^{{\\rm +";
  outlatex << setprecision(3) << allPara[5][1] << "}}_{{\\rm ";
  outlatex << setprecision(3) << allPara[5][2] << "}}$ & ";

  // Average sSFR
  outlatex << setprecision(3) << allPara[6][0] - 11 << "$^{{\\rm +";
  outlatex << setprecision(3) << allPara[6][1] << "}}_{{\\rm ";
  outlatex << setprecision(3) << allPara[6][2] << "}}$   & ";

  // chi2
  outlatex << setprecision(3) << allPara[7][0] << "  \\\\" << endl;

  // OUTPUT PARAMETERS

  outpara << std::fixed;
  outpara << setprecision(1) << zmins << setprecision(1) << " " << zmaxs << " ";

  outpara << setprecision(5) << allPara[0][0] << " ";
  outpara << setprecision(5) << allPara[0][2] << " ";
  outpara << setprecision(5) << allPara[0][1] << " ";

  if (fittype != 3) {
    outpara << setprecision(3) << allPara[1][0] << " ";
  } else {
    // Case double-exponential, to fall back on the right definition
    outpara << setprecision(3) << allPara[1][0] + 3 << " ";
  }
  outpara << setprecision(5) << allPara[1][2] << " ";
  outpara << setprecision(5) << allPara[1][1] << " ";

  outpara << setprecision(5) << allPara[2][0] - 11 << " ";
  outpara << setprecision(5) << allPara[2][2] << " ";
  outpara << setprecision(5) << allPara[2][1] << " ";

  outpara << setprecision(5) << allPara[3][0] << " ";
  outpara << setprecision(5) << allPara[3][2] << " ";
  outpara << setprecision(5) << allPara[3][1] << " ";

  outpara << setprecision(5) << allPara[4][0] << " ";
  outpara << setprecision(5) << allPara[4][2] << " ";
  outpara << setprecision(5) << allPara[4][1] << " ";

  outpara << setprecision(5) << allPara[5][0] << " ";
  outpara << setprecision(5) << allPara[5][1] << " ";
  outpara << setprecision(5) << allPara[5][2] << " ";

  outpara << setprecision(5) << allPara[6][0] << " ";
  outpara << setprecision(5) << allPara[6][1] << " ";
  outpara << setprecision(5) << allPara[6][2] << " " << endl;

  return;
}
