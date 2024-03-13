/*

29/05/2014
Program to fit the sSFR parameters

*/

#include <fstream> // print output file
#include <iostream> // print standard file
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <algorithm>    // sort
#include <stdlib.h>     // abs, exit 
#include <cstring>
#include <iomanip>      // std::setprecision


//MINUIT
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnSimplex.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnScan.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/MinosError.h"
#include "Minuit2/ContoursError.h"

#include "LFdata.h"    
#include "ssfrparaFCN.h" //Function to be minimized in minuit


using namespace std;
using namespace ROOT::Minuit2;

int main (int argc, char *argv[])
{

 //Stock the ssfr value in this (even if it was not the original design of the object
 LFdata SSFRone; 
 vector<LFdata> ssfrData;
 ifstream ssfrfiles,onessfrst;
 string lit,name; 
 double a,zmin,zmax,massmed;


 // to estimate the times needed to run the code
 clock_t t1,t2;
 t1=clock();

 //READ THE SSFR IN THE PARAMETERS FILE

 //list of the ssfr files
 ssfrfiles.open("parameter_files.dat");
 //Take the stream line by line
 while(getline( ssfrfiles, lit)){

   //put the line into the stream ss again 
   stringstream ss(lit);
   ss >> name >> massmed;
   //cout << name << " " << massmed << endl;

   //Open the ssfr parameter file
   onessfrst.open(name.c_str());
   while(getline( onessfrst, lit)){

     //put the line into the stream ss again 
     stringstream ss(lit);
     ss >> zmin >> zmax >> a >> a >> a >> a >> a >> a >> a >> a >> a >> a >> a >> a >> a >> a >> a >> SSFRone.value  >> SSFRone.errp >> SSFRone.errm;
     SSFRone.z=(zmin+zmax)/2.;
     SSFRone.M=massmed;
     //cout << " ssfr " <<  SSFRone.z << " " << SSFRone.value  << "  " << SSFRone.errp << " " <<  SSFRone.errm << endl;

     //Everything in one vector
     ssfrData.push_back(SSFRone);
   }
   
 onessfrst.close();

 }


  /*
  START THE FIT
  */

  //Set the parameters
  MnUserParameters upar;
  upar.Add("Norma",10 ,0.5);     //sSFR noralisation
  upar.Add("amass",-0.2,0.1);    //dependency with the mass
  upar.Add("zmass" ,3.4,0.2);    //dependency with z
  
  //Remove the dependency with the mass
  upar.SetValue(1,0);     
  upar.Fix(1);

  //Construct the function that will be minimize by initializing the data to be fit
  ssfrparaFCN theFCN(ssfrData);

  // simplex, then migrad, redo if migrad not successfull,  maximum 4 calls
  int ncall=0;
  bool successmini=false;
  while (ncall<10) {
 
    //SIMPLEX minimizer
    MnSimplex simplex(theFCN, upar,2);
    FunctionMinimum simplexmin = simplex();
    
    
    //Set the result of simplex to migrad input
    for (int k=0 ; k<3 ; k++){upar.SetValue(k,simplex.Value(k));}

    //MIGRAD minimizer
    MnMigrad migrad(theFCN, upar,2);
    FunctionMinimum min = migrad();

    successmini = min.IsValid();
    //cout << "Min " << ncall << "  " << min << endl;

     //In case MIGRAD was sucessfull
    if(successmini){

     //MINOS ERRORS
     MnMinos minos(theFCN,min);
     //for each parameter
     for (int k=0 ; k<3 ; k++){
      //Compute minos errors when parameter not fixed
       if(!min.UserState().Parameter(k).IsFixed()){
        //Compute minos errors when parameter not fixed
         pair<double,double> errminos=minos(k);
         cout << "Parameter: "<< k << " Value: " << min.UserState().Value(k) << " Negative err: " << errminos.first  << " Positive err: " << errminos.second << endl;
       }else{
        cout << "Parameter: "<< k << " Value: " << min.UserState().Value(k) << endl;
       }
     }
     break;

    }else{ 
      // Migrad failed, Restart with the migrad values in input of simplex
      for (int k=0 ; k<5 ; k++){upar.SetValue(k,migrad.Value(k));} 
    }
      
    ncall++;
 }

 

 // display the time needed to run the code
 t2=clock();
 float diff (((float)t2-(float)t1) / CLOCKS_PER_SEC);
 cout << "Times to run the code in sec: " << diff <<endl;


 return 0;

}
 
