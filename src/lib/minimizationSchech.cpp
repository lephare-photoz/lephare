/*

  08/01/2014
  Core of the minimization
  Use minuit intensively in this function

  output: A vector with the free parameters to be minimized, each element
  including a vector of values, negative and positive errors

  input LFall: The LF data to be fitted
  input zbin : The redshift bin that is considered
  input upar : The set of parameters which are minimized by minuit
  input fittype: The type of function to minimize (LF or MF)
  input eddi=sigma of the gaussian to convolve the schechter func

*/

#include <math.h>    // power
#include <stdlib.h>  // abs, exit

#include <cstring>
#include <fstream>   // print output file
#include <iostream>  // print standard file
#include <sstream>
#include <string>
#include <vector>

#include "LFdata.h"
#include "Minuit2/ContoursError.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MinosError.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnScan.h"
#include "Minuit2/MnSimplex.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnUserParameters.h"
#include "schechterFCN.h"   //Function to be minimized in minuit
#include "schechterFunc.h"  //Function to be minimized in minuit

using namespace std;
using namespace ROOT::Minuit2;

vector<vector<double>> minimizationSchech(vector<LFdata> LFall, const int zbin,
                                          MnUserParameters upar,
                                          const int fittype,
                                          const double sigconv) {
  vector<vector<double>> allPara;
  vector<double> interPara, parafin, paraDupli, meas, pos, err;
  vector<int> lowLim;
  double normal, normalerr;

  // Put all the mesurements in a vector used to initialize the FCN function
  for (vector<LFdata>::iterator it = LFall.begin(); it < LFall.end(); ++it) {
    if (zbin == (*it).bin) {
      meas.push_back(pow(10., (*it).value));
      pos.push_back((*it).M);
      // If the positive error is negative, we face a lower limit
      if ((*it).errp > 0) {
        err.push_back(pow(10., (*it).value) * (pow(10., (*it).errp) - 1.));
        lowLim.push_back(0);
      } else {
        err.push_back(pow(10., (*it).value) *
                      (pow(10., -1. * (*it).errp) - 1.));
        lowLim.push_back(1);
      }
      // keep the normalisation if needed to be fixed
      normal = (*it).nor;
      normalerr = (*it).norerr;
    }
  }

  // Construct the function that will be minimize by initializing the data to be
  // fit
  schechterFCN theFCN(meas, pos, err, lowLim, fittype, normal, normalerr,
                      sigconv);

  // simplex, then migrad, redo if migrad not successfull,  maximum 4 calls
  int ncall = 0;
  bool successmini = false;
  while (ncall < 10) {
    // SIMPLEX minimizer
    MnSimplex simplex(theFCN, upar, 2);
    FunctionMinimum simplexmin = simplex();
    // cout << "Simplex " << ncall << "  " << simplexmin  << endl;

    // Set the result of simplex to migrad input
    for (int k = 0; k < 5; k++) {
      upar.SetValue(k, simplex.Value(k));
    }

    // MIGRAD minimizer
    MnMigrad migrad(theFCN, upar, 2);
    FunctionMinimum min = migrad();
    successmini = min.IsValid();
    cout << "Min " << ncall << "  " << min << endl;

    // In case MIGRAD was sucessfull
    if (successmini) {
      // MINOS ERRORS
      MnMinos minos(theFCN, min);
      cout << "Minos " << ncall << "  " << min << endl;

      // for each parameter
      for (int k = 0; k < 5; k++) {
        // Store the best fit parameters
        interPara.clear();
        interPara.push_back(min.UserState().Value(k));
        if (!min.UserState().Parameter(k).IsFixed()) {
          // Compute minos errors when parameter not fixed
          pair<double, double> errminos = minos(k);
          interPara.push_back(errminos.first);
          interPara.push_back(errminos.second);
          // Errorbars from migrad
          // interPara.push_back(-1*min.UserState().Error(k));
          // interPara.push_back(min.UserState().Error(k));
        } else {
          // Error at -999. for the parameter which is fixed
          interPara.push_back(-999.);
          interPara.push_back(-999.);
        }
        allPara.push_back(interPara);
      }

      // compute the median and averaged SFR
      // Keep the best fit parameters
      parafin.push_back(allPara[0][0]);
      parafin.push_back(allPara[1][0]);
      parafin.push_back(allPara[2][0]);
      parafin.push_back(allPara[3][0]);
      parafin.push_back(allPara[4][0]);
      // Compute the chi2, SFR med, SFR moy for the original parameters
      double prenorma = 0;
      double premed = -99;
      if (fittype == 4 || fittype == 5) {
        SchechterFunc schechpre(3, parafin);
        if (fittype == 4) {
          schechpre.chtype(2);
        }
        prenorma = schechpre.schechInt(0., 20.);
        premed = log10(schechpre.schechIntLD(0., 20.) /
                       schechpre.schechInt(0., 20.));
      }
      SchechterFunc errSchech(fittype, parafin);
      if (fittype == 4 || fittype == 5) {
        errSchech.prenorma = prenorma;
        errSchech.premed = premed;
      }
      double schechMedOri = errSchech.Mmed(0., 20.);
      double schechMoyOri =
          log10(errSchech.schechIntLD(0., 20.) / errSchech.schechInt(0., 20.));

      // create contours factory with FCN and Minimum
      int nbsample = 50;
      vector<double> a1cont;
      vector<double> mcont;
      bool successcont = false;
      // check that mstar is not fixed
      if (allPara[0][2] > -990) {
        // generate the contour
        MnContours contours(theFCN, min);
        // delta chi2+1 contour
        theFCN.SetErrorDef(1);
        vector<pair<double, double>> cont = contours(0, 2, nbsample);
        // plot the contours
        MnPlot plot;
        plot(min.UserState().Value(0), min.UserState().Value(2), cont);
        // put the contours into two vectors alpha1-mstar
        for (vector<pair<double, double>>::iterator it = cont.begin();
             it < cont.end(); ++it) {
          a1cont.push_back((*it).first);
          mcont.push_back((*it).second);
        }
        successcont = 'true';
      }

      // Case in which the parameter alpha 1 is set, or if the full contour is
      // not spanned Span the range with the mstar error
      if (!successcont || a1cont.size() < nbsample) {
        for (int iter = 0; iter < nbsample; ++iter) {
          double step = (allPara[2][2] - allPara[2][1]) / double(nbsample);
          a1cont.push_back(allPara[0][0]);
          mcont.push_back(allPara[2][0] +
                          (iter - double(nbsample / 2.)) * step);
        }
      }

      // section to find the error on the median/moy
      double schechMoyErrn = 999.9;
      double schechMoyErrp = -999.9;
      double schechMedErrn = 999.9;
      double schechMedErrp = -999.9;
      int niter = 0;
      // Keep the best chi2
      double bestchi2 = theFCN(parafin);
      // Loop over all the points of the contours
      for (int iter = 1; iter < nbsample; iter++) {
        double a1b = a1cont[iter];
        double mb = mcont[iter];
        // duplicate the paramters to modify them
        paraDupli.clear();
        paraDupli.push_back(a1b);
        paraDupli.push_back(allPara[1][0]);
        paraDupli.push_back(mb);
        paraDupli.push_back(allPara[3][0]);
        paraDupli.push_back(allPara[4][0]);
        cout << " Iter " << niter << " Para: " << a1b << "  " << mb << " \r "
             << flush;
        cout << " Iter " << niter << " Para: " << a1b << "  " << mb << "  "
             << "  " << bestchi2 << " " << theFCN(paraDupli) << " "
             << bestchi2 + 1 << endl;

        // Compute the chi2, SFR med, SFR moy for the original parameters
        double prenorma = 0;
        double premed = -99;
        if (fittype == 4 || fittype == 5) {
          SchechterFunc schechpre(3, paraDupli);
          if (fittype == 4) {
            schechpre.chtype(2);
          }
          prenorma = schechpre.schechInt(0., 20.);
          premed = log10(schechpre.schechIntLD(0., 20.) /
                         schechpre.schechInt(0., 20.));
        }
        SchechterFunc errSchech(fittype, paraDupli);
        if (fittype == 4 || fittype == 5) {
          errSchech.prenorma = prenorma;
          errSchech.premed = premed;
        }
        double schechMed = errSchech.Mmed(0., 20.);
        if ((schechMed - schechMedOri) > schechMedErrp) {
          schechMedErrp = (schechMed - schechMedOri);
        }
        if ((schechMed - schechMedOri) < schechMedErrn) {
          schechMedErrn = (schechMed - schechMedOri);
        }
        double schechMoy = log10(errSchech.schechIntLD(0., 20.) /
                                 errSchech.schechInt(0., 20.));
        if ((schechMoy - schechMoyOri) > schechMoyErrp) {
          schechMoyErrp = (schechMoy - schechMoyOri);
        }
        if ((schechMoy - schechMoyOri) < schechMoyErrn) {
          schechMoyErrn = (schechMoy - schechMoyOri);
        }
        niter++;
      }

      // Store the averaged function as another parameter
      interPara.clear();
      interPara.push_back(schechMedOri);
      interPara.push_back(schechMedErrp);
      interPara.push_back(schechMedErrn);
      allPara.push_back(interPara);
      // Store the median function as another parameter
      interPara.clear();
      interPara.push_back(schechMoyOri);
      interPara.push_back(schechMoyErrp);
      interPara.push_back(schechMoyErrn);
      allPara.push_back(interPara);

      // Add for the last parameters the best-chi2
      interPara.clear();
      interPara.push_back(min.UserState().Fval());
      interPara.push_back(-999);
      interPara.push_back(-999);
      allPara.push_back(interPara);

      // Return the best fit parameters and stop here the function
      return allPara;

    } else {
      // Migrad failed, Restart with the migrad values in input of simplex
      for (int k = 0; k < 5; k++) {
        upar.SetValue(k, migrad.Value(k));
      }
    }

    ncall++;
  }

  cout << endl
       << "Never successful in the minimization for zbin " << zbin << endl;
  return allPara;
}
