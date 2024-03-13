/*

  03/01/14
  Implementation of the operator () for the FCN minuit class

*/

#include "schechterFCN.h"
#include "schechterFunc.h"
#include "math.h"
#include <sstream> // print output file
#include <fstream> // print output file
#include <iostream> // print standard file
#include <stdlib.h>     /* abs */

#include <cassert>

using namespace std;
using namespace ROOT::Minuit2;

double schechterFCN::operator()(const vector<double>& par) const {

  assert(par.size() == 5);
  double chi2 = 0.;

  //alpha2 need to be lower than alpha1
  //if(par[0]>par[1] && par[3]>0.1){
  if(par[0]>par[1]){

    //Calcule d'abord la mediane et le norma dans le cas ou on ajoute un starburst
    double prenorma=0;
    double premed=-99 ;
    if(ftype==4 || ftype==5){
      SchechterFunc schechpre(3,par);
      if(ftype==4){schechpre.chtype(2);}
      prenorma=schechpre.schechInt(0.,20.);
      premed=log10(schechpre.schechIntLD(0.,20.)/schechpre.schechInt(0.,20.)) ;    //schechpre.Mmed(0.,20.) ;
      //prenorma=par[3];
      //premed=par[2];
      /*if(fabs(premed-par[2])>0.3){
        cout << "Prob " << premed-par[2] <<endl;
        cout << " Ben, c'est la med ? " << premed << " Mstar " << par[2]  << endl;
        cout << " Ben, c'est la norma ? " << prenorma <<  " p1 " <<  par[3] << endl;
        cout << " Reste " <<  " a1 " << par[0] << " a2 " << par[1] <<  " p2 " << par[4] << endl;
	}*/
    }

   //Generate a Schechter function with the considered free parameters
    SchechterFunc schech(ftype,par);
    if(ftype==4 || ftype==5){     
      schech.prenorma=prenorma;
      schech.premed=premed;
    }
  

   //Loop over the data to be fitted
   for(unsigned int n = 0; n < fMeasurements.size(); n++) {
     double diff=schech(fPositions[n],fsigconv) - fMeasurements[n];
       // Normal computation of the chi2 if not in lower limit
       // or if prediction are below the measurement
       if(fLim[n]==0 || diff<0){
	 chi2 += pow( diff/fMVariances[n] ,2.);}
   }

   //Add the constrain on the normalisation if used
   if(fnormal>0 && fnormalerr>0){
    //compute the integral of the Schechter function
    double inte=schech.schechInt(0.,20.);
    chi2 += pow( (inte-fnormal)/fnormalerr ,2.);

   }
   

   //Put the chi2 at huge value if ftype==4 SB mode and phi2>phi1
   if((ftype==4 || ftype==5) && par[4]>0.2){chi2=10e10;}

  }else chi2=10e10;

  return chi2;

}
