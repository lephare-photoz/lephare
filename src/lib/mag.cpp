/*
 *  Class with all the tools used to create the grid of SED in e(b-v), redshift,
 * etc Create the library of modeled magnitude and k-corrections
 */

#include "mag.h"

#include <algorithm>
#include <cmath>
#include <fstream>  // print output file
#include <functional>
#include <iomanip>   // std::set precisionls
#include <iostream>  // print standard file
#include <sstream>
#include <string>
#include <vector>

// Constructor of the basis class which read the keywords common to the
// QSO/STARS/GAL
Mag::Mag(keymap &key_analysed) {
  /*
    ENVIRONMENT VARIABLES LEPHAREDIR and LEPHAREWORK
  */
  get_lephare_env();

  // Configuration file
  config = key_analysed["c"].value;

  // type of source which is read (Galaxy G, QSO Q, Star S)
  typ = key_analysed["t"].value;

  // Instantiate cosmology
  double h0 = (key_analysed["COSMOLOGY"].split_double("70", 3))[0];
  double om0 = (key_analysed["COSMOLOGY"].split_double("0.3", 3))[1];
  double l0 = (key_analysed["COSMOLOGY"].split_double("0.7", 3))[2];
  // Define one object with this cosmology
  lcdm = cosmo(h0, om0, l0);

  // Load filters
  filtFile = ((key_analysed["FILTER_FILE"]).split_string("filters", 1))[0];
  // Full path to the input filter file
  string fltFile = lepharework + "/filt/" + filtFile + ".dat";
  read_flt(fltFile);

  // mag type AB/VEGA
  magtyp = ((key_analysed["MAGTYPE"]).split_string("AB", 1))[0];

  // extinction laws, multiple laws are possible, number of expected laws
  // unknown in advance -> -1
  extlaw = (key_analysed["EXTINC_LAW"]).split_string("calzetti.dat", -1);
  nextlaw = int(extlaw.size());
  read_ext();

  // possible E(B-V) values, multiple values are possible, number of expected
  // values unknown in advance -> -1
  ebv = (key_analysed["EB_V"]).split_double("0", -1);
  nebv = int(ebv.size());
  // model ranges for each extinction curve
  modext = (key_analysed["MOD_EXTINC"]).split_int("0,0", nextlaw * 2);

  // type of the grid in redshift
  gridType = ((key_analysed["ZGRID_TYPE"]).split_int("0", 1))[0];

  // define the grid in redshift
  dz = ((key_analysed["Z_STEP"]).split_double("0.04", 3))[0];
  zmin = ((key_analysed["Z_STEP"]).split_double("0.", 3))[1];
  zmax = ((key_analysed["Z_STEP"]).split_double("6.", 3))[2];
  if (zmax < zmin) {
    throw runtime_error(
        "You are probably using the old parametrisation of "
        "Z_STEP since Z MIN > Z MAX in Z_STEP. Stop here. ");
  }
  // Output file in ascii ?
  outasc = ((key_analysed["LIB_ASCII"]).split_bool("NO", 1))[0];

  // keyword to add the LDUST component to the stellar component (e.g. in BC03)
  addDust = key_analysed["ADD_DUSTEM"].value;

  // Want to display the template number on the screen
  // VERBOSE output  file -  YES default
  verbose = key_analysed["VERBOSE"].split_bool("YES", 1)[0];

  // Read the extragalactic opacity files into a vector
  opaAll = read_opa();
}

// destructor of the class Mag cleaning all the vectors
Mag::~Mag() {
  close_files();

  extlaw.clear();
  ebv.clear();
  modext.clear();

  extAll.clear();
  opaAll.clear();
  allFlt.clear();
  gridz.clear();
  gridT.clear();
  gridDM.clear();
}

// opens the common files for GAL/QSO/STARS  that will be used
void Mag::open_files() {
  // Name of the input binary SED file
  sedlibFile = lepharework + "/lib_bin/" + lib + ".bin";
  ssedIn.open(sedlibFile.c_str(), ios::binary);
  // Check if file is opened
  if (!ssedIn) {
    throw invalid_argument("Can't open file " + sedlibFile);
  }

  // Name of the input SED doc file
  sedlibFile = lepharework + "/lib_bin/" + lib + ".doc";

  // Name of the output doc file
  docFile = lepharework + "/lib_mag/" + colib + ".doc";
  sdocOut.open(docFile.c_str());
  // Check if file is opened
  if (!sdocOut) {
    throw invalid_argument("Can't open file " + docFile);
  }

  // Name of the output binary and doc files, open it in binary
  binOutFile = lepharework + "/lib_mag/" + colib + ".bin";
  sbinOut.open(binOutFile.c_str(), ios::binary | ios::out);
  // Check if file is opened
  if (!sbinOut) {
    throw invalid_argument("Can't open file " + binOutFile);
  }

  // Open an ascii file to store the predicted magnitudes
  if (outasc) {
    // Name of the output ascii file
    datFile = lepharework + "/lib_mag/" + colib + ".dat";
    sdatOut.open(datFile.c_str());
    // Check if file is opened
    if (!sdatOut) {
      throw invalid_argument("Can't open file " + datFile);
    }
    // header of the .dat file
    sdatOut
        << "# model ext_law E(B-V) L_T(IR) redshift dist_modulus age record "
           "N_filt magnitude_vector kcorr_vector em_lines_fluxes_vector "
        << endl;
  }

  cout << " All files opened " << endl;
}

// open the opacity files
ifstream Mag::open_opa_files() {
  ifstream stream;
  // open the ascii file with all the opacity file listed
  string opaListFile = lepharedir + "/opa/OPACITY.dat";
  stream.open(opaListFile.c_str());
  // Check if file is opened
  if (!stream) {
    throw invalid_argument("Can't open file " + opaListFile);
  }
  return stream;
}

void Mag::close_files() {
  sbinOut.close();
  sdocOut.close();
  ssedIn.close();
  sdatOut.close();
}

// Function of the basis class which read the extinction laws
void Mag::read_ext() {
  // Loop over the possible extinction laws
  for (int k = 0; k < nextlaw; k++) {
    // Instance one ext object
    ext oneext(extlaw[k], k);
    // Name of the extinction law file
    string extFile = lepharedir + "/ext/" + extlaw[k];
    // read the extinction law file
    oneext.read(extFile);
    // store it into the vector of exction laws
    extAll.push_back(oneext);
  }
  // Read the MW extinction curve and store it into the last item
  // Do not increment nextlaw
  ext oneext("MW_seaton.dat", nextlaw);
  string extFile = lepharedir + "/ext/MW_seaton.dat";
  oneext.read(extFile);
  extAll.push_back(oneext);
}

// Function of the basis class which read the IGM opacity
vector<opa> Mag::read_opa() {
  string name;
  double red;

  // In oder to fill the two last elements around Lyman alpha
  // Put 1 for the last element
  // Put the last value of the opa below 1215.67 just before
  oneElLambda beflastOpa(1215.66, 1., 3);
  oneElLambda lastOpa(1215.67, 1., 3);

  ifstream stream = Mag::open_opa_files();
  vector<opa> result;

  // Take the stream line by line: list of each opa file
  for (int i = 0; i < 81; i++) {
    stream >> red >> name;
    opa oneOpa(red, name);
    oneOpa.read();
    // Put as last element a lambda at the Lyman-alpha wavelength with
    // transmission=1 Meiksin case : remove the last element which is after the
    // Lya line
    if (oneOpa.lamb_opa.back().lamb > 1215.66) oneOpa.lamb_opa.pop_back();
    // Put the last transmission value very close to Lyman alpha
    beflastOpa.val = oneOpa.lamb_opa.back().val;
    // Add the two last values close to Lyman alpha
    oneOpa.lamb_opa.push_back(beflastOpa);
    oneOpa.lamb_opa.push_back(lastOpa);
    oneOpa.lmax = 1215.67;
    // Add to the list of opacity
    result.push_back(oneOpa);
  }
  return result;
}

// Function of the basis class which read all the filters
void Mag::read_flt(const string &inputfile) {
  ifstream sfiltIn;
  sfiltIn.open(inputfile.c_str());
  // Check if file is opened
  if (!sfiltIn) {
    throw invalid_argument("Can't open file " + inputfile);
  }

  string dummy;
  int imag;
  // read the number of filter
  sfiltIn >> dummy >> imag;

  // Loop over each filter
  for (int k = 0; k < imag; k++) {
    // Generate one object "flt" and read it
    flt oneFilt(k, sfiltIn, 0, 0);
    // store all filters in a vector
    allFlt.push_back(oneFilt);
  }

  sfiltIn.close();
}

// Read the long wavelength Bethermin+2012 templates to add the dust emission to
// the BC03 templates Associate a b12 SED to each redshift of the grid in
// redshift
void Mag::read_B12() {
  /*
  IMPORTANT NOTE
  There is one limitation with current implementation of the code:
  If several templates from B12 are used (not the first one by default), the fit
  and predicted magnitudes will be correct. But the best-fit template in the
  .spec file will be off in FIR since it is based on the z=0 full template to be
  reconstructed No easy fix yet.
  */

  // Open the file with the list of B12 templates
  ifstream b12mod;
  string b12List = lepharedir + "/sed/GAL/BETHERMIN12/BETHERMIN12_MOD.list";
  b12mod.open(b12List.c_str());
  if (!b12mod)
    throw invalid_argument("Can't open Bethermin+12 list " + b12List);

  // Create a list of SED with the B12 templates. Need one SED for each redshift
  // of gridz.
  string lit, nameSED, bid;
  double b12z;
  size_t gr = 0;
  while (getline(b12mod, lit)) {
    // If the first character of the line is not #
    if (check_first_char(lit)) {
      // Read the SED name + information
      stringstream ss(lit);
      ss >> nameSED >> bid >> b12z;
      string sedpath = lepharedir + "/sed/GAL/" + nameSED;

      // Use as template number the value of the redshift x 100
      GalSED oneSEDFIR(sedpath, int(b12z * 100.));
      oneSEDFIR.read(sedpath);

      // First step
      if (gr == 0) {
        B12SED.push_back(oneSEDFIR);
        gr++;
      }

      // We need one SED for each step of the gridz
      while (((gridz[gr] + gridz[gr + 1]) / 2.) <= b12z && gr > 0 &&
             gr < gridz.size() - 1) {
        B12SED.push_back(oneSEDFIR);
        gr++;
      }
    }
  }
  // Add the SED to go from the last redshift in B12 template to the last
  // redshift of the grid
  for (size_t k = gr; k < gridz.size(); k++) B12SED.push_back(B12SED.back());

  // Close the stream
  b12mod.close();
}

// Define the redshift grid
// Associate it to a grid in age and distance modulus
void Mag::def_zgrid() {
  // redshift grid, depending on the method
  gridz = zgrid(gridType, dz, zmin, zmax);

  // Loop over the redshift grid and measure the age of the Universe and the
  // distance modulus
  for (size_t k = 0; k < gridz.size(); k++) {
    gridT.push_back(lcdm.time(gridz[k]));
    gridDM.push_back(lcdm.distMod(gridz[k]));
  }
}

// General information on the screen valid for GAL/QSO/Stars
void Mag::print_info() {
  cout << "#######################################" << endl;
  cout << "# It s computing the SYNTHETIC MAGNITUDES #" << endl;
  cout << "# For Gal/QSO libraries with these OPTIONS #" << endl;
  cout << "# with the following options :           " << endl;
  cout << "# Config file     : " << config << endl;
  cout << "# Filter file     : " << filtFile << endl;
  cout << "# Magnitude type     : " << magtyp << endl;
}

// Write the documentation in the GALAXY/QSO/STAR case
void Mag::write_doc() {
  sdocOut << "CONFIG_FILE    " << config << endl;
  if (typ[0] == 'G' || typ[0] == 'g') {
    sdocOut << "LIB_TYPE  GALAXY" << endl;
  } else if (typ[0] == 'Q' || typ[0] == 'q') {
    sdocOut << "LIB_TYPE  QSO" << endl;
  } else if (typ[0] == 'S' || typ[0] == 's') {
    sdocOut << "LIB_TYPE  STAR" << endl;
  }
  sdocOut << "LIB_NAME      " << lib << endl;
  sdocOut << "FILTER_FILE     " << filtFile << endl;
  sdocOut << "FILTERS   ";
  for (vector<flt>::iterator itf = allFlt.begin(); itf < allFlt.end(); ++itf)
    sdocOut << itf->name << ",";
  sdocOut << endl << "MAG_TYPE      " << magtyp << endl;
  sdocOut << "AB_COR   ";
  for (vector<flt>::iterator itf = allFlt.begin(); itf < allFlt.end(); ++itf)
    sdocOut << itf->abcorr() << ",";
  sdocOut << endl << "FLUX_COR   ";
  for (vector<flt>::iterator itf = allFlt.begin(); itf < allFlt.end(); ++itf)
    sdocOut << itf->fcorrec() << ",";
  sdocOut << endl << "ZGRID_TYPE   " << gridType << endl;
  sdocOut << "Z_STEP   " << dz << "," << zmin << "," << zmax << endl;
  sdocOut << "COSMOLOGY   " << lcdm << endl;
  sdocOut << "EXTINC_LAW   ";
  for (int k = 0; k < nextlaw; k++) {
    sdocOut << extlaw[k] << ",";
  };
  sdocOut << endl << "MOD_EXTINC   ";
  for (int k = 0; k < nextlaw; k++) {
    sdocOut << modext[k * 2] << "," << modext[2 * k + 1] << ",";
  };
  sdocOut << endl << "EB_V   ";
  for (int k = 0; k < nebv; k++) {
    sdocOut << ebv[k] << ",";
  };
  sdocOut << endl << "EM_LINES   " << emlines << endl;
  sdocOut << "LIB_ASCII   " << (outasc ? "YES" : "NO") << endl;
  time_t result = time(nullptr);
  sdocOut << "CREATION_DATE " << asctime(std::localtime(&result));
}

/*
  FUNCTIONS SPECIFIC TO THE GALAXY
  They inherit from the basis "Mag" class
*/

// constructure of the Galaxy case
// read the keywords missed by the constructor of the basis class
GalMag::GalMag(keymap &key_analysed) : Mag(key_analysed) {
  // Name of the input file, default value "SED"
  lib = ((key_analysed["GAL_LIB_IN"]).split_string("SED", 1))[0];
  // Name of the output file, default value "LIB"
  colib = ((key_analysed["GAL_LIB_OUT"]).split_string("LIB", 1))[0];

  // Emission lines in output
  emlines = ((key_analysed["EM_LINES"]).split_string("EMP_UV", 1))[0];
  // If 'yes' as the old keyword, swich to EMP_UV which should be always working
  if (emlines[0] == 'y' || emlines[0] == 'Y') emlines = "EMP_UV";
  // Check the the keyword has an expected value, otherwise stop
  if (emlines.substr(0, 6).compare("EMP_UV") != 0 &&
      emlines.substr(0, 7).compare("EMP_SFR") != 0 &&
      emlines.substr(0, 4).compare("PHYS") != 0 &&
      emlines.substr(0, 2).compare("NO") != 0) {
    throw invalid_argument(
        "Not the right value for EM_LINES (NO, EMP_SFR, EMP_UV, PHYS)");
  }
  // possible dispersion within the emission lines
  // indicate how many step do you want
  // Don't fill it if EM_LINES NO
  if (emlines.substr(0, 2).compare("NO") != 0) {
    fracEm = (key_analysed["EM_DISPERSION"]).split_double("1", -1);
  } else {
    fracEm.push_back(1.);
  }
}

// Only for galaxies
// Read the binary file containing the SEDs (file created by sedtolib)
// Redshift the SED and apply the extinction and opacity
// Create the emission lines
// Compute the modeled magnitude. Done here to optimize the library size
void GalMag::read_SED() {
  // Obtain system statistics. Check the memory.
#ifdef _linux_
  struct sysinfo si;
  const double megabyte = 1024 * 1024;
#endif

  // find the end of the stream for the binary input file, then go back at the
  // beginning
  ssedIn.seekg(0, ssedIn.end);
  int length = ssedIn.tellg();
  ssedIn.seekg(0, ssedIn.beg);

  // Read until the end of the stream -> all the SED read
  while (ssedIn.tellg() < length) {
    GalSED oneSED("");
    // read one SED in the binary file
    oneSED.readSEDBin(ssedIn);
    //// Check that the lambda coverage is correct
    // oneSED.warning_integrateSED(allFlt, verbose);
    // build the library of SEDs that modify the initial template SED
    vector<GalSED> seds = make_maglib(oneSED);
    // write the result in file
    write_mag(seds);

    // Stop the code if only 2% of the memory left
#ifdef _linux
    {
      sysinfo(&si);
      if (double(si.freeram) / double(si.totalram) < 0.02) {
        cout << "Need to stop the process. Not enough memory.";
        cout << "Free RAM (MegaB) " << si.freeram / megabyte << endl;
        cout << "Total RAM (MegaB) " << si.totalram / megabyte << endl;
        cout << "Possible to subdivide the library if necessary, or reduce the "
                "parameter space."
             << endl;
        throw runtime_error();
      }
    }
#endif
  }  // end of while loop
}

vector<GalSED> GalMag::make_maglib(GalSED &oneSED) {
  vector<GalSED> allSED;
  // build the emission line SED. This changes the state of oneSED
  GalSED oneEm = oneSED.generateEmSED(emlines);

// PARALLELIZE all the 4 loops  [Iary, 12 March 2018]
#pragma omp parallel for ordered schedule(dynamic) collapse(4)
  // Loop over each extinction law
  for (int i = 0; i < nextlaw; i++) {
    // loop over each E(B-V)
    for (int j = 0; j < nebv; j++) {
      // loop over each fraction of emission line flux (add a dispersion in
      // emission lines as a new template)
      for (size_t l = 0; l < fracEm.size(); l++) {
        // Loop over the redshift grid
        for (size_t k = 0; k < gridz.size(); k++) {
          // Select case which need to be considered (no extinction or
          // extinction in the right model range) Remove all cases with
          // extinction not in the right model range The condition i==0 only
          // means that for null extinction the templates are computed only
          // once, using the first extinction law, and it does not matter which
          // one this extinction law is.
          if ((ebv[j] < 1.e-10 && i == 0) ||
              (ebv[j] > 0 && oneSED.nummod >= modext[i * 2] &&
               oneSED.nummod <= modext[i * 2 + 1])) {
            // Generate intermediate Continuum SED, since original one must not
            // change
            GalSED oneSEDInt(oneSED);

            // galaxy redshift/distance in the grid
            oneSEDInt.red = gridz[k];
            oneSEDInt.distMod = gridDM[k];

            // Check that the lambda coverage is correct
            oneSEDInt.warning_integrateSED(allFlt, verbose);

            // Not older than the age of the Universe
            if (gridT[k] > oneSEDInt.age) {
              double LbeforeExt = oneSEDInt.trapzd();

              // product of the SED with the extinction law
              oneSEDInt.applyExt(ebv[j], extAll[i]);

              // Difference between the integrated flux with and without
              // extinction (without is computed just above) flux integrate of
              // the lambda range -> erg/s/cm2. It was for the source at 10pc,
              // and then convert erg/s in Lsol
              double dL = (LbeforeExt - oneSEDInt.trapzd()) / Lsol *
                          (4 * pi * 100 * pow(pc, 2));
              if (oneSEDInt.ltir < 0 && dL > 0) oneSEDInt.ltir = log10(dL);
              // Rescale the B12 to the right dust luminosity (with energy
              // balance) and sum to the stellar continuum is option on.
              if (addDust[0] == 'Y' || addDust[0] == 'y') {
                oneSEDInt.sumSpectra(B12SED[k], dL);
              }

              // Opacity applied in rest-frame, depending on the redshift of the
              // source
              oneSEDInt.applyOpa(opaAll);

              // redshift the SED, and restrict it to the union of the filters
              // support.
              oneSEDInt.redshift();
              oneSEDInt.reduce_memory(allFlt);
              // Compute magnitude
              // Loop over the filters
              oneSEDInt.compute_magnitudes(allFlt);
              // If z>0, no need to keep the spectra
              if (oneSEDInt.red > 1.e-10) oneSEDInt.lamb_flux.clear();

              // Derive the emission line flux in each filter
              if (emlines[0] == 'E' || emlines[0] == 'P') {
                // Generate intermediate EM SED, since original one must not
                // change
                GalSED oneEmInt(oneEm);
                oneEmInt.ebv = ebv[j];
                // set the value of fracEm
                oneEmInt.fracEm = fracEm[l];
                oneEmInt.red = gridz[k];
                // For the emission lines, use only the MW. Change fac_line
                oneEmInt.applyExtLines(extAll[nextlaw]);
                // rescale the lines as a free parameter
                oneEmInt.rescaleEmLines();
                /*
                // Decide to not applied.
                // apply a z dependence of the emission line ratio for OIII
                oneEmInt.zdepEmLines(1);
                */
                // Generate the spectra with the emission lines
                oneEmInt.generateEmSpectra(40);
                // Opacity applied in rest-frame, depending on the redshift of
                // the source
                oneEmInt.applyOpa(opaAll);
                // Save the emission lines rest-frame in the continuum SED
                oneSEDInt.fac_line = oneEmInt.fac_line;
                //
                oneEmInt.redshift();
                oneEmInt.rescale(pow(10., -0.4 * oneSEDInt.distMod));
                oneEmInt.reduce_memory(allFlt);
                if (oneEmInt.lamb_flux.size() > 0) {
                  oneSEDInt.flEm = oneEmInt.compute_fluxes(allFlt);
                } else {
                  oneSEDInt.flEm.assign(allFlt.size(), 0.);
                }
                // indicate that the emission lines have been computed
                oneSEDInt.has_emlines = true;
                if (oneSEDInt.red > 1.e-10) oneEmInt.lamb_flux.clear();
                oneEmInt.clean();
              }

#pragma omp ordered
              {
                // add to all SED (one time if ebv==0)
                allSED.push_back(oneSEDInt);

                // Display in the right order, even when the code is
                // parrallelized
                if (verbose) {
                  cout << "SED " << oneSEDInt.name << " z " << setw(6)
                       << oneSEDInt.red;
                  cout << " Ext law " << extlaw[i] << "  E(B-V) " << ebv[j]
                       << "  Age " << oneSEDInt.age << "  \r " << flush;
                }
                // Cleaning
                oneSEDInt.clean();
              }
            }  // close age condition
          }
        }
      }
    }
  }
  // Now take all the SED for the current initial template
  // Compute the K-correction, and save to file.
  for (size_t k = 0; k < allSED.size(); k++) {
    // compute k-correction
    if (allSED[k].red < 1.e-5) {
      // keep the magnitude at z=0 and put the k-correction at 0
      magko = allSED[k].mag;
      allSED[k].kcorr.assign(allFlt.size(), 0.);
    } else {
      // compute the k-correction mag(z)-mag(z=0)
      for (size_t itf = 0; itf < allFlt.size(); itf++) {
        allSED[k].kcorr.push_back(allSED[k].mag[itf] - allSED[k].distMod -
                                  magko[itf]);
      }
    }
  }
  return allSED;
}

void GalMag::write_mag(const vector<GalSED> &seds) {
  // write the output files
  for (const auto &sed : seds) {
    sed.writeMag(outasc, sbinOut, sdatOut, allFlt, magtyp);
  }
}

/*
  Information in the screen related to the galaxies
*/
void GalMag::print_info() {
  Mag::print_info();

  cout << "# GAL_LIB_IN    :"
       << lepharework + "/lib_bin/" + lib + "(.doc & .bin)" << endl;
  cout << "# GAL_LIB_OUT   :"
       << lepharework + "/lib_mag/" + colib + "(.doc & .bin)" << endl;
  cout << "# ZGRID_TYPE   :" << gridType << endl;
  cout << "# Z_STEP   :" << dz << " " << zmin << " " << zmax << endl;
  cout << "# COSMOLOGY   :" << lcdm << endl;
  cout << "# EXTINC_LAW   :";
  for (int k = 0; k < nextlaw; k++) {
    cout << extlaw[k] << " ";
  };
  cout << endl << "# MOD_EXTINC   :";
  for (int k = 0; k < nextlaw; k++) {
    cout << modext[k * 2] << " " << modext[2 * k + 1] << " ";
  };
  cout << endl << "# EB_V   :";
  for (int k = 0; k < nebv; k++) {
    cout << ebv[k] << " ";
  };
  cout << endl << "# EM_LINES   " << emlines << endl;
  cout << "# EM_DISPERSION   ";
  for (size_t k = 0; k < fracEm.size(); k++) {
    cout << fracEm[k] << ",";
  };
  cout << endl << "# LIB_ASCII   " << (outasc ? "YES" : "NO") << endl;
  time_t result = time(nullptr);
  cout << "# CREATION_DATE " << asctime(std::localtime(&result));
  cout << "#############################################" << endl;
}

/*
  SPECIFIC TO THE QSOs
*/

// constructor for the QSO adding the missing keywords from the basis
// constructor
QSOMag::QSOMag(keymap &key_analysed) : Mag(key_analysed) {
  // Name of the input file, default value "SED"
  lib = ((key_analysed["QSO_LIB_IN"]).split_string("SED", 1))[0];
  // Name of the output file, default value "LIB"
  colib = ((key_analysed["QSO_LIB_OUT"]).split_string("LIB", 1))[0];
}

/*
  Information in the screen related to the QSO
*/
void QSOMag::print_info() {
  Mag::print_info();

  cout << "# QSO_LIB_IN    :"
       << lepharework + "/lib_bin/" + lib + "(.doc & .bin)" << endl;
  cout << "# QSO_LIB_OUT   :"
       << lepharework + "/lib_mag/" + colib + "(.doc & .bin)" << endl;
  cout << "# ZGRID_TYPE   :" << gridType << endl;
  cout << "# Z_STEP   :" << dz << " " << zmin << " " << zmax << endl;
  cout << "# COSMOLOGY   :" << lcdm << endl;
  cout << "# EXTINC_LAW   :";
  for (int k = 0; k < nextlaw; k++) {
    cout << extlaw[k] << " ";
  };
  cout << endl << "# MOD_EXTINC   :";
  for (int k = 0; k < nextlaw; k++) {
    cout << modext[k * 2] << " " << modext[k + 1] << " ";
  };
  cout << endl << "# EB_V   :";
  for (int k = 0; k < nebv; k++) {
    cout << ebv[k] << " ";
  };
  cout << "# LIB_ASCII   " << (outasc ? "YES" : "NO") << endl;
  time_t result = time(nullptr);
  cout << "# CREATION_DATE " << asctime(std::localtime(&result));
  cout << "#############################################" << endl;
}

// Only for QSO
// Read the binary file containing the SEDs (file created by sedtolib)
// Redshift the SED and apply the extinction and opacity
void QSOMag::read_SED() {
  // Obtain system statistics. Check the memory.
#ifdef _linux_
  struct sysinfo si;
  const double megabyte = 1024 * 1024;
#endif

  // find the end of the stream for the binary input file, then go back at the
  // beginning
  ssedIn.seekg(0, ssedIn.end);
  int length = ssedIn.tellg();
  ssedIn.seekg(0, ssedIn.beg);

  // Read until the end of the stream -> all the SED read
  while (ssedIn.tellg() < length) {
    // create initial template SED from current position in SED binary file
    QSOSED oneSED("");
    oneSED.clean();
    oneSED.readSEDBin(ssedIn);
    //// Check that the lambda coverage is correct
    // oneSED.warning_integrateSED(allFlt, verbose);
    // build the library of SEDs that modify the initial template SED
    vector<QSOSED> seds = make_maglib(oneSED);
    // write the result in file
    write_mag(seds);

    // Stop the code if only 2% of the memory left
#ifdef _linux
    {
      sysinfo(&si);
      if (double(si.freeram) / double(si.totalram) < 0.02) {
        cout << "Need to stop the process. Not enough memory.";
        cout << "Free RAM (MegaB) " << si.freeram / megabyte << endl;
        cout << "Total RAM (MegaB) " << si.totalram / megabyte << endl;
        cout << "Possible to subdivide the library if necessary, or reduce the "
                "parameter space."
             << endl;
        throw runtime_error();
      }
    }
#endif
  }
  return;
}

vector<QSOSED> QSOMag::make_maglib(const QSOSED &oneSED) {
  vector<QSOSED> allSED;
#pragma omp parallel for ordered schedule(dynamic) collapse(3)
  // Loop over each extinction law
  for (int i = 0; i < nextlaw; i++) {
    // loop over each E(B-V)
    for (int j = 0; j < nebv; j++) {
      // Loop over the redshift grid
      for (size_t k = 0; k < gridz.size(); k++) {
        // Select case which need to be considered (no extinction or extinction
        // in the right model range) Remove all cases with extinction not in the
        // right model range
        if ((ebv[j] < 1.e-10 && i == 0) ||
            (ebv[j] > 0 && oneSED.nummod >= (modext[i * 2]) &&
             oneSED.nummod <= (modext[i * 2 + 1]))) {
          // Generate intermediate Continuum SED, since original one must not
          // change
          QSOSED oneSEDInt(oneSED);

          // galaxy redshift/distance in the grid
          oneSEDInt.red = gridz[k];
          oneSEDInt.distMod = gridDM[k];

          // Check that the lambda coverage is correct
          oneSEDInt.warning_integrateSED(allFlt, verbose);

          // product of the SED with the extinction law
          oneSEDInt.applyExt(ebv[j], extAll[i]);

          // Opacity applied in rest-frame, depending on the redshift of the
          // source
          oneSEDInt.applyOpa(opaAll);

          // redshift the SED
          oneSEDInt.redshift();

          // Compute magnitude
          oneSEDInt.compute_magnitudes(allFlt);

          // If z>0, no need to keep the spectra
          if (oneSEDInt.red > 1.e-10) oneSEDInt.lamb_flux.clear();

#pragma omp ordered
          {
            // add to all SED (one time if ebv==0)
            allSED.push_back(oneSEDInt);

            // Info screen
            // Display in the right order, even when the code is parrallelized
            if (verbose) {
              cout << "SED " << oneSEDInt.name << " z " << setw(6)
                   << oneSEDInt.red;
              cout << " Ext law " << extlaw[i] << "  E(B-V) " << ebv[j]
                   << "  \r " << flush;
            }
            // Cleaning
            oneSEDInt.clean();
          }
        }
      }
    }
  }
  // Now take all the SED for the current initial template
  // Compute the K-correction
  for (size_t k = 0; k < allSED.size(); k++) {
    // compute k-correction
    if (allSED[k].red < 1.e-5) {
      // keep the magnitude at z=0 and put the k-correction at 0
      magko = allSED[k].mag;
      allSED[k].kcorr.assign(allFlt.size(), 0.);
    } else {
      for (size_t itf = 0; itf < allFlt.size(); itf++) {
        allSED[k].kcorr.push_back(allSED[k].mag[itf] - allSED[k].distMod -
                                  magko[itf]);
      }
    }
  }
  return allSED;
}

void QSOMag::write_mag(const vector<QSOSED> &seds) {
  // write the output files
  for (const auto &sed : seds) {
    sed.writeMag(outasc, sbinOut, sdatOut, allFlt, magtyp);
  }
}

/*
  SPECIFIC TO THE STARS
*/

// Constructor of the stars adding the keywords missing in the basis class
// constructor
StarMag::StarMag(keymap &key_analysed) : Mag(key_analysed) {
  // Name of the input file, default value "SED"
  lib = ((key_analysed["STAR_LIB_IN"]).split_string("SED", 1))[0];
  // Name of the output file, default value "LIB"
  colib = ((key_analysed["STAR_LIB_OUT"]).split_string("LIB", 1))[0];
}

// Information in the screen related to the QSO
void StarMag::print_info() {
  Mag::print_info();
  cout << "# COSMOLOGY   :" << lcdm << endl;

  cout << "# STAR_LIB_IN    :"
       << lepharework + "/lib_bin/" + lib + "(.doc & .bin)" << endl;
  cout << "# STAR_LIB_OUT   :"
       << lepharework + "/lib_mag/" + colib + "(.doc & .bin)" << endl;
  cout << "# LIB_ASCII   " << (outasc ? "YES" : "NO") << endl;
  time_t result = time(nullptr);
  cout << "# CREATION_DATE " << asctime(std::localtime(&result));
  cout << "#############################################" << endl;
}

// Read the binary file containing the SEDs (file created by sedtolib)
// and compute the magnitudes
void StarMag::read_SED() {
  // find the end of the stream for the binary input file
  ssedIn.seekg(0, ssedIn.end);
  int length = ssedIn.tellg();
  ssedIn.seekg(0, ssedIn.beg);

  // read until the end of the stream -> all the SED read
  while (ssedIn.tellg() < length) {
    // Instance one SED to be read
    StarSED oneSED("");
    // read one SED in the binary file
    oneSED.readSEDBin(ssedIn);
    // Check that the lambda coverage is correct
    oneSED.warning_integrateSED(allFlt, verbose);

    vector<StarSED> seds = make_maglib(oneSED);
    write_mag(seds);

    // Info screen
    if (verbose)
      cout << "Read SED " << oneSED.name << " z " << setw(6) << oneSED.red
           << "  \r " << flush;
  }
}

vector<StarSED> StarMag::make_maglib(const StarSED &sed) {
  vector<StarSED> allSED;
  StarSED newsed(sed);
  // compute magnitude for the template directly,
  // as for a star no other extinction or redshifting is applied
  newsed.compute_magnitudes(allFlt);
  // return singleton vector in order to have the same structure as for QSO and
  // Gal
  allSED.push_back(newsed);
  return allSED;
}

void StarMag::write_mag(const vector<StarSED> &seds) {
  // write the output files
  for (const auto &sed : seds) {
    sed.writeMag(outasc, sbinOut, sdatOut, allFlt, magtyp);
  }
}
