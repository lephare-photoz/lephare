/*
    08/01/2014
    Program to fit a LF or a luminosity/mass function to non parametric data

*/

#include <fstream> // print output file
#include <iostream> // print standard file
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <algorithm>    // sort
#include <stdlib.h>     // abs, exit 
#include <getopt.h>     // get option line commands getopt_long
#include <typeinfo>
#include <cstring>
#include <iomanip>      // std::setprecision

//Le Phare
#include "keyword.h"    //our own class to read the keywords
#include "LFdata.h"    

//MINUIT
#include "schechterFunc.h" 
#include "Minuit2/MnUserParameters.h"

using namespace std;
using namespace ROOT::Minuit2;

static const string LEPHAREWORK = "LEPHAREWORK";
static const string LEPHAREDIR = "LEPHAREDIR";

//Declare function prototypes
double strtodouble(const string & inputstring);
bool test_first_char(string maligne);
vector<keyword> read_config  (istream & configst,int & nb_keywords);
vector<keyword> read_command  (int argc, char *argv[], int & nb_keywords, string & configfile);
int find_keyword (const keyword onekey, const string * list_keywords, const int nb_ref_key);
vector<double> split_double(const string s, char delim);
void read_LFdata  (istream & LFst, vector<LFdata> &  LFall, const double zmoy, const int k, const double norma, const double normaerr);
vector<double> opt_demulti(const vector<double>  optin, const int Nzbin) ;
vector< vector<double> > minimizationSchech(const vector<LFdata> LFall, const int k,  const MnUserParameters upar, const int fittype, const double sigconv);
void write_output_fit  (ofstream & outlatex, ofstream & outpara, const string out, const int fittype, const  vector < vector<double> > allPara, const double zmins, const double zmaxs, const int ngals, const double sigconv);


int main (int argc, char *argv[])
{

 //List of the keyword to be found in the config file/command line
 string list_keywords[] =  {"FIT_INPUT", "FIT_OUTPUT", "FIT_FIX_ALPHA1", "FIT_FIX_ALPHA2","FIT_FIX_MSTAR", "FIT_FIX_PHI1", "FIT_FIX_PHI2", "EDDINGTON","FIT_TYPE","FIX_NORMA","FIX_NORMAERR"};
 //Number of keywords
 int nb_ref_key= (int)(sizeof(list_keywords)/sizeof(list_keywords[0]));

 //List of the keywords to be found in the info file
 string list_info[] =  {"Bin", "Number", "Output"};
 //Number of keywords
 int nb_info= (int)(sizeof(list_info)/sizeof(list_info[0]));

 //parameters read in option line if necessary
 string inputname, outputname,lepharedir,lit,configfile;
 vector<double> a1fix,a2fix,phi1fix,phi2fix,mstarfix,zbin,zbmin,zbmax,normafix,normafixerr,eddi;
 vector<int> NbGal;
 vector<keyword> key_config,key_command,key_all,key_info;
 vector<string> inputNoPara;
 int nb_keywords,fittype;
 string lepharework = "./";
 char const* temp;
 ifstream configst, ifstream,inputst,LFst;
 ofstream outputst,outlatex,outpara;

 //Minuit
 vector<double> init_par;

 //All  
 vector<LFdata> LFall;
 vector<double> bestPara, errnPara, errpPara;
 vector < vector<double> > allPara;
 double Mini;

 // to estimate the times needed to run the code
 clock_t t1,t2;
 t1=clock();

 /*
 ENVIRONMENT VARIABLES
 */

 //Check LEPHAREDIR and LEPHAREWORK. Stop if LEPHAREDIR not defined.
 temp = getenv(LEPHAREDIR.c_str());
 if(temp != NULL){lepharedir = string(temp);}else{cout << "Environment variable LEPHAREDIR not defined, need to stop." << endl; exit(1);}   
 temp = getenv(LEPHAREWORK.c_str());
 if(temp != NULL){lepharework = string(temp);}else{cout << "Environment variable LEPHAREWORK not defined, used local directory." << endl;}   

  /*
  ANALYSE KEYWORDS
  */

 //Read the command line and print the keywords
 key_command=read_command( argc, argv,nb_keywords,configfile);
 cout << "Number of keywords read in the command line: " << nb_keywords << endl;

 //stop if no config file in the command line
 if(configfile.length()==0){ 
     cout << "No configuration file indicated with the -c option. Need to stop" << endl; 
     exit(2);
  }else{
     configst.open(configfile.c_str());
  }

 //Read the configuration file and print the keywords
 key_config=read_config(configst,nb_keywords);
 cout << "Number of keywords read in the config file: " << nb_keywords << endl;
 configst.close();

 //concatanate all the keywords from config and from command line
 key_all.reserve( key_config.size() + key_command.size()); 
 key_all.insert( key_all.end(), key_config.begin(), key_config.end() );
 key_all.insert( key_all.end(), key_command.begin(), key_command.end() );

 //Loop over all the keywords
 for (vector<keyword>::iterator it=key_all.begin(); it<key_all.end(); ++it){

   //Name of the keyword which match
    int iarg=find_keyword(*it,list_keywords, nb_ref_key);

    //FIT_INPUT:0  FIT_OUTPUT:1 FIT_FIX_ALPHA1:2 FIT_FIX_ALPHA2=3 FIT_FIX_MSTAR=4 FIT_FIX_PHI1=5 FIT_FIX_PHI2=6 EDDINGTON=7 FIT_TYPE 8 FIX_NORMA 9 FIX_NORMAERR 10
    switch (iarg)
    {
    case -1:
      //cout << " No rules for keyword " << (*it).name << endl;
    case 0:
      inputname= (*it).value;
      break;
    case 1:
      outputname= (*it).value;
      break;
    case 2:
      a1fix=split_double((*it).value,',');
      break;
    case 3:
      a2fix=split_double((*it).value,',');
      break;
    case 4:
      mstarfix=split_double((*it).value,',');
      break;
    case 5:
      phi1fix=split_double((*it).value,',');
      break;
    case 6:
      phi2fix=split_double((*it).value,',');
      break;
    case 7:
      eddi=split_double((*it).value,',');
      break;
    case 8:
      fittype = atoi(((*it).value).c_str());
      break;
    case 9:
      normafix=split_double((*it).value,',');
      break;
    case 10:
      normafixerr=split_double((*it).value,',');
      break;
    }
 }

 //Decide what is the initial Mstar
 switch (fittype){                   
   case 0:
     Mini=-20.;	     //LF	 
    break;
   case 1:
     Mini=10.;	     //MF	 
    break;
   default:
     Mini=10.;	     //SFRF	 
 }

  
 /* 
 READ THE INFO FILE 
 */

 //open the info file 
 inputst.open(inputname.c_str());
 //Read the info file and print the keywords
 key_info=read_config(inputst,nb_keywords);
 inputst.close();
 
 //Name of the output file for the latex table
 string outlat="default.out";
 outlat= outputname +  ".latex";
 outlatex.open(outlat.c_str());
 string outpar="default.out";
 outpar= outputname +  ".para";
 outpara.open(outpar.c_str());

 //Loop over all the keywords
 for (vector<keyword>::iterator it=key_info.begin(); it<key_info.end(); ++it){

    //Name of the keyword which match
    // case 0:"Bin" 1:"Number" 2:"Output"
    int iarg=find_keyword(*it,list_info, nb_info);

    switch (iarg)
    {
    case 0:
      zbin=split_double((*it).value,',');
      zbmin.push_back(zbin[0]);  // minimum value of the redshift bin
      zbmax.push_back(zbin[1]);  // maximum value of the redshift bin
      break;
    case 1:
      NbGal.push_back(atoi(((*it).value).c_str()));  // Number of galaxies used to computed the LF
      break;
    case 2:
      inputNoPara.push_back((*it).value); // File in which the LF data are stored for this redshift bin
      break;
    }
 } 
  

 /*
  CHECK THAT THE KEYWORDS HAVE THE RIGHT DIMENSIONS
 */
 
 unsigned ref_size=zbmin.size();
 vector<double> a1fixe      =opt_demulti(a1fix,int(ref_size));
 vector<double> a2fixe      =opt_demulti(a2fix,int(ref_size));
 vector<double> mstarfixe   =opt_demulti(mstarfix,int(ref_size));
 vector<double> phi1fixe    =opt_demulti(phi1fix,int(ref_size));
 vector<double> phi2fixe    =opt_demulti(phi2fix,int(ref_size));
 vector<double> normafixe   =opt_demulti(normafix,int(ref_size));
 vector<double> normafixeerr=opt_demulti(normafixerr,int(ref_size));
 //Basic check on the fixed parameters
 for (int k=0 ; k<(int)(phi1fixe.size()) ; k++){
   //If phistar1 =0, alpha1 must be fixed too
   if(phi1fixe[k]==0)a1fixe[k]=5.;
   //If simple gaussian fit, phi*_2=0
   if(fittype==2)phi2fixe[k]=0; 
   //If phistar2 =0, alpha2 must be fixed too
   if(phi2fixe[k]==0 && fittype!=3 && fittype!=5)a2fixe[k]=-5.; //dont fix alpha 2 if double-epo
   //If alpha2 fixed for a double expo, change the value to be consistent with our definion
   if(a2fixe[k]>-90 && a2fixe[k]<90 && (fittype==3 || fittype==5))a2fixe[k]=a2fixe[k]-3;
   //sigma2 fixed for double gaussian with SB mode
   if(fittype==4)a2fixe[k]=-5.; //dont fix alpha 2 if double-epo
 }


 // Check that all the vectors read in the .info file have the same size 
 if( (ref_size != zbmax.size()) || (ref_size != inputNoPara.size()) || (ref_size != NbGal.size())){
	cout << "All the  vectors read in the .info file do NOT have the same size ";
	cout << "Sizes: " << zbmin.size() << " " <<  zbmax.size() << " " << inputNoPara.size() << " " << NbGal.size() << endl;
 }



 /*
 INFO PARAMETERS ON SCREEN 
 */

 cout <<  "#######################################" << endl;
 cout <<  "# Fit schechter with OPTIONS   #" << endl;
 cout <<  "# CONFIG FILE  : " << configfile << endl;
 cout <<  "# INFO  FILE, NON PARARAMETRIC ESTIMATOR : " << inputname << endl;
 cout <<  "# OUTPUT FILE  : " << outputname << endl;
 for (int k=0 ; k<(int)zbmin.size() ; k++){
     cout << "# Redshift bin: " << zbmin[k] << " " <<  zbmax[k] << " File: " << inputNoPara[k] << " Nb.Gal: " <<  NbGal[k] << endl;
 }
 cout <<  "#######################################" << endl;

 /*
 READ THE LF DATA FILE
 */

 // Loop on each redshift bin
 for (int k=0 ; k<(int)zbmin.size() ; k++){

  // Open LF/MF file, read it, close it
  LFst.open((inputNoPara[k]).c_str());
  read_LFdata(LFst,LFall,((zbmin[k]+zbmax[k])/2.),k,normafixe[k],normafixeerr[k]);
  LFst.close();

 } 
 

 /*
 START THE FIT
 */

 // Loop on each redshift bin
 // Parallelize the code for this loop !!
 //#pragma omp parallel for ordered schedule(dynamic) private(allPara)
 for (int k=0 ; k<(int)zbmin.size() ; k++){

  //Set the parameters
  MnUserParameters upar;
  upar.Add("alpha1",0.5 ,0.25);       //alpha 1	 
  upar.Add("alpha2",-1.,0.5);	     //alpha 2	 
  upar.Add("mstar" ,Mini,0.25);       //Mstar	 
  upar.Add("phi1",1.000000,0.5);     //phi star 1
  upar.Add("phi2",0.10000000,0.1);     //phi star 2
  upar.SetLowerLimit(3,0.001);     //only positive phistar
  upar.SetLowerLimit(4,-1.e-40);     //only positive phistar
  if(fittype>=2 && fittype<=5)upar.SetLowerLimit(0,0.1); //only positive sigma for double-expo et gaussienne
  if(fittype==2 || fittype==4)upar.SetUpperLimit(0,2); //upper limit for sigma



  //Fix parameters if required
  if(a1fixe[k]>-90 && a1fixe[k]<90){
     upar.SetValue(0,a1fixe[k]);     
     upar.Fix(0);     
  }
  if(a2fixe[k]>-90 && a2fixe[k]<90){
     upar.SetValue(1,a2fixe[k]);     
     upar.Fix(1);     
  }
  if(mstarfixe[k]>-90 && mstarfixe[k]<90){
     upar.SetValue(2,mstarfixe[k]);     
     upar.Fix(2);     
  }
  if(phi1fixe[k]>-90 && phi1fixe[k]<90){
     upar.SetValue(3,phi1fixe[k]);     
     upar.Fix(3);     
  }
  if(phi2fixe[k]>-90 && phi2fixe[k]<90){
     upar.SetValue(4,phi2fixe[k]);     
     upar.Fix(4);     
  }

  //Perform the full minimization
  cout << "Eddington " << eddi[k] << endl;
  allPara=minimizationSchech(LFall,k,upar,fittype,eddi[k]);


  //Name of the output file
  string out="default.out";
  if(k<10){out= outputname + char(k+48) + ".dat";}
    else if(k>=10 && k<20){out= outputname + "1" + char(k+38) + ".dat";}
    else{cout << "Too many redshift bins (>20) for a good naming " << endl;}


  //Check that the output makes sense (converged)
  if(allPara.size()>0){

   //Display in the right order, even when the code is parrallelized
   //#pragma omp ordered
   {

    //write latex, parameter files
     write_output_fit(outlatex, outpara, out, fittype, allPara,zbmin[k],zbmax[k],NbGal[k],eddi[k]);

   }

  }else{

    cout <<  "No best fit parameters for the redshift bin:" << zbmin[k] << " " <<  zbmax[k]   << endl;
    ofstream outst;
    outst.open(out.c_str());
    outst << "# No best fit parameters for the redshift bin " << endl << " -999 -999 " << endl;
    outst.close();
   }



 }//loop on the redshift bins
 
 outlatex.close();
 outpara.close();

 // display the time needed to run the code
 t2=clock();
 float diff (((float)t2-(float)t1) / CLOCKS_PER_SEC);
 cout << "Times to run the code in sec: " << diff <<endl;

 return 0;
}
