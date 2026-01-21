/*
 * SEDLib.h
 *
 */

#ifndef SEDLIB_H_
#define SEDLIB_H_

#include <ctime>
#include <filesystem>

#include "SED.h"
#include "globals.h"
#include "keyword.h"

vector<GalSED> readBC03(string sedFile, int nummod, vector<double> &ageSel);
vector<GalSED> readPEGASE(string sedFile, int nummod, vector<double> &ageSel);

/*!
 * @brief Identifies which ages in a library match a set of target selection
 * ages.
 * * This function creates a boolean mask for a vector of ages. If no specific
 * target ages are provided (only the bounds), all ages are selected. Otherwise,
 * it finds the closest match in the library for each target age, provided
 * the match is within a tolerance (1e-5) and falls within the specified
 * [agemin, agemax] range.
 * * @param ageSel A vector where:
 * - ageSel[0] is the minimum allowed age (agemin).
 * - ageSel[1] is the maximum allowed age (agemax).
 * - ageSel[2...n] are the specific target ages to select.
 * @param age    The library of available ages to be filtered.
 * * @return std::vector<bool> A mask of the same size as the 'age' parameter.
 * True if the age is selected, False otherwise.
 * * @note If agemax <= 0, the upper bound check is ignored.
 * @note The matching tolerance is strictly hardcoded at 1e-5.
 */
vector<bool> closeAge(vector<double> ageSel, vector<double> age);

/*!
 * @brief Reads galaxy ages from a file and returns them in a vector.
 *
 * This function creates a vector where the first two elements are agemin and
 * agemax. If a valid filename is provided, it parses the file line-by-line,
 * skipping comments (lines starting with #), and converts the values from Gyr
 * to years.
 *
 * @param filename Name of the file to read. Use "none" to skip file reading.
 * @param agemin   Minimum age boundary (stored at index 0).
 * @param agemax   Maximum age boundary (stored at index 1).
 *
 * @return std::vector<double> A vector containing [agemin, agemax, age1*1e9,
 * age2*1e9, ...].
 *
 * @note Ages in the file are assumed to be in Gyr and are converted to years.
 * @warning If the file cannot be opened, an error is printed to cerr but the
 * function returns the vector containing only agemin and agemax.
 */
vector<double> read_ages_from_file(string filename, double agemin,
                                   double agemax) {
  ifstream sage;
  double dage;
  string lit;

  // Put agemin and agemax in the two first elements of ageSel
  vector<double> ages;
  ages.push_back(agemin);
  ages.push_back(agemax);

  // If the file with the ages exists
  if (filename != "none") {
    // Take the stream line by line
    sage.open(filename.c_str());
    // Check if file has opened properly
    if (!sage) {
      cerr << "Can't open file with the ages to be selected " << filename
           << endl;
      cerr << "No selection by age. " << endl;
      // throw "Failing opening ",filename.c_str();
    }

    while (getline(sage, lit)) {
      // If the first character of the line is not #
      if (check_first_char(lit)) {
        // put the line into the stream ss again
        stringstream ss(lit);
        ss >> dage;

        // fill the age vector, converting to yr.
        ages.push_back(dage * 1.e9);
      }
    }

    sage.close();
  }
  return ages;
}

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

  /*! \brief read content of one SED file into a SED vector
   *
   * @param sedFile the file to read the SED from
   * @param sedFormat format of \a sedFile : can be B(C03), P or F for PEGASE
   * type, or else plain ASCII
   * @param nummod index of the SED; see SED
   * @param type type of the SED S|Q|G for star|qso|galaxy; see SED
   !*/
  virtual void readSED(string sedFile, string sedFormat, int nummod);
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
  modList = ((key_analysed[typ + "_SED"]).split_string("SED.list", 1))[0];
  libOut = ((key_analysed[typ + "_LIB"]).split_string("SED.bin", 1))[0];
  fscale = ((key_analysed[typ + "_FSCALE"]).split_double("1", 1))[0];

  open_output_files();
}

template <>
SEDLib<GalSED>::SEDLib(keymap &key_analysed, string config, string t)
    : SEDLib(config, t) {
  path = "/sed/" + typ + "/";
  modList = ((key_analysed[typ + "_SED"]).split_string("SED.list", 1))[0];
  libOut = ((key_analysed[typ + "_LIB"]).split_string("SED.bin", 1))[0];
  fscale = ((key_analysed[typ + "_FSCALE"]).split_double("1", 1))[0];

  open_output_files();

  ageFile = ((key_analysed["SEL_AGE"]).split_string("none", 1))[0];
  // Range of ages to be considered for galaxies
  agemin =
      ((key_analysed["AGE_RANGE"]).split_double(to_string(INVALID_PHYS), 2))[0];
  agemax =
      ((key_analysed["AGE_RANGE"]).split_double(to_string(INVALID_PHYS), 2))[1];

  ageSel = read_ages_from_file(ageFile, agemin, agemax);
}

template <class T>
void SEDLib<T>::write_SED_lib() {
  // Loop over the SED of QSO
  for (typename vector<T>::iterator it = allSED.begin(); it < allSED.end();
       ++it) {
    // Rescale the flux of each SED according to the factor given in keyword
    it->rescale(fscale);
    // Compute some SED properties: implemented for GalSED only
    it->compute_luminosities();
    // Write the SED in the output binary file
    it->writeSED(sbinOut, sphysOut, sdocOut);
  }
}

template <class T>
void SEDLib<T>::readSED(string sedFile, string sedFormat, int nummod) {
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
    throw invalid_argument("Can't open doc file of the SED library in " +
                           docFile);
  }

  binFile = lepharework + "/lib_bin/" + libOut + ".bin";
  sbinOut.open(binFile.c_str(), ios::binary | ios::out);
  if (!sbinOut) {
    throw invalid_argument("Can't open binary file of the SED library in " +
                           binFile);
  }

  if (typ == "GAL") {
    physFile = lepharework + "/lib_bin/" + libOut + ".phys";
    sphysOut.open(physFile.c_str());
    if (!sphysOut) {
      throw invalid_argument("Can't open phys file of the SED library in " +
                             physFile);
    }
    // header
    sphysOut << "# age luv lopt lnir  ltir mass  sfr  zmet  tau  d4000 qi"
             << endl;
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
    throw invalid_argument("Can't open file with the list of SED to be used " +
                           modList);
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
      // check if nameSED is an absolute path
      if (std::filesystem::path(nameSED).is_relative())
        nameSED = lepharedir + path + nameSED;

      formatSED = 'A';  // Default: ascii
      if (!ss.eof()) ss >> formatSED;
      // Read the file and output a vector of SED
      // (in some file, you have several SEDs with different ages)
      readSED(nameSED, formatSED, nbSED + 1);
      nbSED++;
    }
  }

  // Close the stream
  smod.close();

  cout << "Number of templates in the list " << nbSED << endl;
  cout << "Number of SED in the list (including different ages) "
       << allSED.size() << endl;
  // Write the documentation
  sdocOut << "NUMBER_SED " << nbSED << endl;
}

template <>
void SEDLib<GalSED>::readSED(string sedFile, string sedFormat, int nummod) {
  resultSED.clear();

  if (sedFormat[0] == 'B') {
    // BC03 case
    resultSED = readBC03(sedFile, nummod, ageSel);

  } else if (sedFormat[0] == 'P' || sedFormat[0] == 'F') {
    // PEGASE
    resultSED = readPEGASE(sedFile, nummod, ageSel);

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
