/*
 * SEDLib.h
 *
 */

#ifndef SEDLIB_H_
#define SEDLIB_H_

#include <ctime>

#include "SED.h"
#include "globals.h"
#include "keyword.h"

vector<GalSED> readBC03(string sedFile, int nummod, string type,
                        vector<double> &ageSel);
vector<GalSED> readPEGASE(string sedFile, int nummod, string type,
                          vector<double> &ageSel);
vector<bool> closeAge(vector<double> ageSel, vector<double> age);

/*!
 * \brief class for a general SED (Star, QSO, or Galaxy) library
 */
template <class T>
class SEDLib {
 private:
  string config,  ///< configuration file (arg -c at the command line)
      typ;  ///< type of the SED instance (arg -t at the command line), either
            ///< S, G, or Q for a star, galaxy, or QSO template, respectively

 protected:
  ofstream sdocOut,  ///< output stream to store the doc in \a #docFile
      sbinOut,       ///< output stream to store the binary SED library in \a
                     ///< #binFile
      sphysOut;      ///< output stream to store the binary SED library in \a
                     ///< #physFile (only for typ=GAL)
  string ageFile;  ///< optional file providing the selection of ages to keep in
                   ///< the SED library;
  double agemin, agemax;
  vector<double> ageSel;
  string physFile;  ///< output file name, locate in $LEPHAREWORK/lib_bin/; used
                    ///< only for typ GAL

 public:
  /*
   * vector, containing the SEDs for the particular type, that will be unique in
   * the derived class
   */
  vector<T> allSED, resultSED;

  string docFile,  ///< output file for the doc file
      binFile;     ///< output file for the binary SED library
  string modList, libOut, path;
  double fscale;
  SEDLib(string config, string typ);
  SEDLib(keymap &key_analysed, string config, string typ);
  virtual ~SEDLib();

  /// write time of creation and the number of SED recorded to the doc file
  void print_time_tofile(time_t result) {
    sdocOut << "CREATION_DATE " << asctime(std::localtime(&result));
  }
  /// print config and type onscreen and in the doc file
  virtual void print_info();
  /// open the output files in $LEPHAREWORK/lib_bin
  virtual void open_output_files();
  // close the output files
  virtual void close_output_files();
  // read the SEDs from the files
  void read_model_list();
  // writes the SED library to the output files (bin and doc)
  // template class<T>
  void write_SED_lib();
  /// \brief read content of one SED file into a SED vector
  ///
  /// @param sedFile the file to read the SED from
  /// @param sedFormat format of \a sedFile : can be B(C03), P or F for PEGASE
  /// type, or else plain ASCII
  /// @param nummod index of the SED; see SED
  /// @param type type of the SED S|Q|G for star|qso|galaxy; see SED
  // template class<T>
  virtual void readSED(string sedFile, string sedFormat, int nummod,
                       string type);
  /// For GAL, read the file with the selected galaxy ages, provided as kw
  /// SEL_AGE
  void read_age(string ageFich);
};

template <class T>
SEDLib<T>::SEDLib(string conf, string t) {
  config = conf;
  typ = t;
  // redefine typ with the right spelling so that it can
  // be used as conditioner in member functions
  if (typ[0] == 'S' || typ[0] == 's') {
    typ = "STAR";
  } else if (typ[0] == 'Q' || typ[0] == 'q') {
    typ = "QSO";
  } else {
    typ = "GAL";
  }

  // ENVIRONMENT VARIABLES LEPHAREDIR and LEPHAREWORK
  get_lephare_env();
}

template <class T>
SEDLib<T>::SEDLib(keymap &key_analysed, string config, string t)
    : SEDLib(config, t) {
  path = "/sed/" + typ + "/";
  // cout<<typ+"_SED"<< " "<<key_analysed[typ+"_SED"]<<endl;
  modList = ((key_analysed[typ + "_SED"]).split_string("SED.list", 1))[0];
  libOut = ((key_analysed[typ + "_LIB"]).split_string("SED.bin", 1))[0];
  fscale = ((key_analysed[typ + "_FSCALE"]).split_double("1", 1))[0];
  if (typ == "GAL") {
    ageFile = ((key_analysed["SEL_AGE"]).split_string("none", 1))[0];
    // Range of ages to be considered for galaxies
    agemin = ((key_analysed["AGE_RANGE"])
                  .split_double(to_string(INVALID_PHYS), 2))[0];
    agemax = ((key_analysed["AGE_RANGE"])
                  .split_double(to_string(INVALID_PHYS), 2))[1];
  }

  open_output_files();

  // if (typ=="GAL" && ageFile!="none"){
  if (typ == "GAL") {
    read_age(ageFile);
  }
}

template <class T>
void SEDLib<T>::write_SED_lib() {
  // Loop over the SED of QSO
  for (typename vector<T>::iterator it = allSED.begin(); it < allSED.end();
       ++it) {
    // Rescale the flux of each SED according to the factor given in keyword
    it->rescale(fscale);
    if (typ == "GAL") {
      // Compute some SED properties
      it->SEDproperties();
      // Add by Cedric to derive the density of ionizing photons
      it->calc_ph();
    }
    // Write the SED in the output binary file
    it->writeSED(sbinOut, sphysOut, sdocOut);
  }
}

template <class T>
void SEDLib<T>::readSED(string sedFile, string sedFormat, int nummod,
                        string type) {
  // Create one object "SED" and fill it with one ascii file
  T oneSEDascii(sedFile, nummod);
  oneSEDascii.read(sedFile);
  // The vector has only one element SED since only one age
  allSED.push_back(oneSEDascii);
}

template <class T>
void SEDLib<T>::print_info() {
  cout << "#######################################" << endl;
  cout << "# It s translating SEDs to binary library #" << endl;
  cout << "# with the following options :           " << endl;
  cout << "# Config file     : " << config << endl;
  cout << "# Library type     : " << typ << endl;
  sdocOut << "CONFIG_FILE " << config << endl;
  sdocOut << "LIB_TYPE    " << typ << endl;

  cout << "# " + typ + "_SED    :" << modList << endl;
  cout << "# " + typ + "_LIB    :" << libOut << endl;
  cout << "# " + typ + "_LIB doc:" << docFile << endl;
  if (typ == "GAL") {
    cout << "# GAL_LIB phys:" << physFile << endl;
    cout << "# SEL_AGE    :" << ageFile << endl;
  }
  cout << "# " + typ + "_FSCALE :" << fscale << endl;
  if (typ == "GAL") {
    cout << "# AGE_RANGE   " << agemin << " " << agemax << endl;
  }
  sdocOut << "" + typ + "_SED    " << modList << endl;
  sdocOut << "" + typ + "_LIB    " << libOut << endl;
  if (typ == "GAL") {
    sdocOut << "# SEL_AGE    :" << ageFile << endl;
  }
  sdocOut << "" + typ + "_FSCALE " << fscale << endl;
  if (typ == "GAL") {
    sdocOut << "# AGE_RANGE   " << agemin << " " << agemax << endl;
  }
  cout << "#######################################" << endl;
}

template <class T>
SEDLib<T>::~SEDLib() {
  close_output_files();
  if (typ == "AGE") {
    ageSel.clear();
  }
}

template <class T>
void SEDLib<T>::open_output_files() {
  docFile = lepharework + "/lib_bin/" + libOut + ".doc";
  sdocOut.open(docFile.c_str());
  if (!sdocOut) {
    throw invalid_argument("Can't open doc file " + docFile);
  }

  binFile = lepharework + "/lib_bin/" + libOut + ".bin";
  sbinOut.open(binFile.c_str(), ios::binary | ios::out);
  if (!sbinOut) {
    throw invalid_argument("Can't open bin file " + binFile);
  }

  if (typ == "GAL") {
    physFile = lepharework + "/lib_bin/" + libOut + ".phys";
    sphysOut.open(physFile.c_str());
    if (!sphysOut) {
      throw invalid_argument("Can't open phys file " + physFile);
    }
  }
}

template <class T>
void SEDLib<T>::close_output_files() {
  sbinOut.close();
  sdocOut.close();
  if (typ == "GAL") {
    sphysOut.close();
  }
}

template <class T>
void SEDLib<T>::read_model_list() {
  string nameSED, formatSED, lit;
  int nbSED = 0;
  ifstream smod;

  // open the template list file into a stream
  smod.open(modList.c_str());
  if (!smod) {
    throw invalid_argument("Can't open mod file " + modList);
  }

  // Take the template list line by line
  while (getline(smod, lit)) {
    // If the first character of the line is not #
    if (check_first_char(lit)) {
      // put the line into the stream ss again
      stringstream ss(lit);

      // Read the name of the SED and the format which need to be used to read
      // it
      ss >> nameSED;
      formatSED = 'A';  // Default: ascii
      if (!ss.eof()) ss >> formatSED;
      string list = lepharedir + path + nameSED;
      // Read the file and output a vector of SED
      // (in some file, you have several SEDs with different ages)
      readSED(list, formatSED, nbSED + 1, typ);
      nbSED++;
    }
  }

  // Close the stream
  smod.close();

  cout << "Number of SED in the list " << nbSED << endl;
  // Write the documentation
  sdocOut << "NUMBER_SED " << nbSED << endl;
}

/*
 * reads the age of the Galaxy, only for type GAL
 */
template <class T>
void SEDLib<T>::read_age(string ageFich) {
  ifstream sage;
  double dage;
  string lit;

  // Put agemin and agemax in the two first elements of ageSel
  ageSel.clear();
  ageSel.push_back(agemin);
  ageSel.push_back(agemax);

  // If the file with the ages exists
  if (ageFich != "none") {
    // Take the stream line by line
    sage.open(ageFile.c_str());
    // Check if file has opened properly
    if (!sage) {
      cerr << "Can't open age file " << ageFile << endl;
      cerr << "No selection by age. " << endl;
      // throw "Failing opening ",ageFile.c_str();
    }

    while (getline(sage, lit)) {
      // If the first character of the line is not #
      if (check_first_char(lit)) {
        // put the line into the stream ss again
        stringstream ss(lit);
        ss >> dage;

        // fill the age vector in Gyr.
        ageSel.push_back(dage * 1.e9);
      }
    }

    sage.close();
  }
}

template <>
void SEDLib<GalSED>::readSED(string sedFile, string sedFormat, int nummod,
                             string type) {
  resultSED.clear();

  if (sedFormat[0] == 'B') {
    // BC03 case
    resultSED = readBC03(sedFile, nummod, type, ageSel);

  } else if (sedFormat[0] == 'P' || sedFormat[0] == 'F') {
    // PEGASE
    resultSED = readPEGASE(sedFile, nummod, type, ageSel);

  } else {
    // ASCII by default
    GalSED oneSEDascii(sedFile, nummod);
    oneSEDascii.read(sedFile);
    resultSED.push_back(oneSEDascii);
  }

  // Add these SEDs to the ones already read
  allSED.insert(allSED.end(), resultSED.begin(), resultSED.end());
}

typedef SEDLib<StarSED> StarSEDLib;
typedef SEDLib<QSOSED> QSOSEDLib;
typedef SEDLib<GalSED> GalSEDLib;

#endif /* SEDLIB_H_ */
