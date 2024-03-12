/*

  06/01/14
  function to have the same vector size as the number of redshift bins

*/


#include <math.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;




vector<double> opt_demulti(const vector<double> optin, const int Nzbin) {

  vector<double> optout;
  int nvec=(int)optin.size();

  switch (nvec)
    {
    case 0:
      for (int k=0 ; k<Nzbin ; k++){optout.push_back(-999.9);};
      break;
    case 1:
      for (int k=0 ; k<Nzbin ; k++){optout.push_back(optin[0]);};
      break;
    default :
      if(nvec!= Nzbin){
        cout << "Wrong number of arguments for option " << endl;
        for (int k=0 ; k<Nzbin ; k++){optout.push_back(-999.9);};}else
        optout=optin;
     }
    return optout;
}
