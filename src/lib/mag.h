/*
 *  Class with all the tools used to create the grid of SED in e(b-v), redshift,
 * etc and create the library of modeled magnitude and k-corrections
 */

#ifndef MAG_H_
#define MAG_H_

#include "SED.h"
#include "cosmology.h"
#include "ext.h"
#include "flt.h"
#include "globals.h"
#include "keyword.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/// abstract class to handle everything related to the computation of expected
/// magnitudes
class Mag {
 private:
 protected:
  object_type
      object;       ///< type of object (GAL/QSO/STAR) handled by this instance
  string config;    ///< path to the configuration file (arg -c)
  cosmo lcdm;       ///< fiducial cosmology used to build the redshift grid
  string filtFile,  ///< path to the filter list file (keyword FILTER_FILE)
      magtyp;       ///< magnitude system, "AB" or "VEGA" (keyword MAGTYPE)
  bool outasc,      ///< whether to also write the library in ASCII format
      verbose,      ///< verbosity flag
      add_dust;     ///< whether dust attenuation should be applied
  vector<string> extlaw;  ///< list of extinction law file names to apply
  vector<double> ebv,     ///< grid of E(B-V) values to apply
      magko;              ///< scratch buffer holding the k-corrected
                          ///< magnitudes of the SED currently being written

  vector<int> modext;  ///< (min,max) SED model index ranges to which each
                       ///< extinction law is restricted (keyword MOD_EXTINC),
                       ///< stored as consecutive pairs
  double dz,           ///< redshift step of the grid
      zmin,            ///< minimum redshift of the grid
      zmax;            ///< maximum redshift of the grid
  string lib,          ///< name of the input SED library
      colib;  ///< base name (without extension) of the output magnitude
              ///< library (keyword GAL_LIB_OUT/QSO_LIB_OUT/STAR_LIB_OUT)

  // only for the galaxy, but much easier to keep them here
  string emlines = "NO";  ///< emission line treatment (e.g. "NO", "EMP",
                          ///< "PHYS"), only relevant for galaxies

  string sedlibFile,  ///< path to the binary SED library to read
      docFile,        ///< path to the output .doc file describing the library
      binOutFile,     ///< path to the output binary magnitude library file
      datFile;        ///< path to the output .dat physical parameter file
  ifstream ssedIn;    ///< input stream for the binary SED library
  ofstream sdocOut,   ///< output stream for #docFile
      sbinOut,        ///< output stream for #binOutFile
      sdatOut;        ///< output stream for #datFile

  ext milkyWayExtinction;        ///< Milky Way extinction law, applied when
                                 ///< #applyMilkyWayExtinction is true
  bool applyMilkyWayExtinction;  ///< whether to apply a Milky Way extinction
                                 ///< correction to the library

 public:
  /// Build a Mag instance from the keywords common to GAL/QSO/STAR runs
  /// (config file, object type, cosmology, IGM opacities)
  /// @param key_analysed: map of keyword/value pairs parsed from the
  /// configuration file and command line
  Mag(keymap& key_analysed);
  Mag(){};
  virtual ~Mag();

  /// read the extinction laws into attribute extAll (vector of vectors of type
  /// ext)
  void read_ext();

  /// Read the long wavelength Bethermin+2012 templates
  /// to add the dust emission to the BC03 templates
  void read_B12();

  /// define the vector of redshifts, and associate to it vectors of age and
  /// distance modulus, based on the lcdm attribute
  void def_zgrid();

  /// helper function to set the grid of redshift
  /// @param dz : step in z in the grid
  /// @param zmin : minimum z
  /// @param zmax : maximum z
  inline void set_zgrid(double dz, double zmin, double zmax) {
    gridz = zgrid(dz, zmin, zmax);
  }

  /// Write in file sdocOut the documentation for the GALAXY/QSO/STAR case
  void write_doc();

  /// print general information onscreen, valid for GAL/QSO/STAR objects
  virtual void print_info();
  /// open all the input and output streams needed for the computations
  void open_files();
  /// close all opened files
  void close_files();
  /// read SED files, apply extinction corrections, and store into a vector of
  /// instances of class SED
  virtual void read_SED() = 0;

  vector<ext> extAll;     ///< extinction laws read from #extlaw
  vector<flt> allFlt;     ///< filters read from #filtFile
  vector<GalSED> B12SED;  ///< Bethermin et al. (2012) dust emission templates
  vector<double> gridz;   ///< redshift grid values
  vector<double> gridT,   ///< age of the Universe at each #gridz step
      gridDM;             ///< distance modulus at each #gridz step

  vector<opa> opaAll;  ///< extragalactic opacity curves used along the line
                       ///< of sight
};

/// inherited class handling expected magnitudes from star SED
class StarMag : public Mag {
 public:
  /// Build a StarMag instance and read the star-specific keywords
  /// @param key_analysed: map of keyword/value pairs
  StarMag(keymap& key_analysed);
  StarMag(){};  // LCOV_EXCL_LINE
  ~StarMag(){};

  void print_info();
  void read_SED();
  /// Build the library of synthetic magnitudes for a single star SED, by
  /// looping over the filter set (no extinction or E(B-V) grid for stars)
  /// @param sed: the star SED to process
  /// @return vector of StarSED, one per output entry of the library
  vector<StarSED> make_maglib(const StarSED& sed);
  /// Write the magnitudes of the input SEDs to the output library files
  /// @param seds: SEDs (with their magnitudes already computed) to write out
  void write_mag(const vector<StarSED>& seds);
};

/// inherited class handling expected magnitudes from QSO SED
class QSOMag : public Mag {
 public:
  /// Build a QSOMag instance and read the QSO-specific keywords
  /// @param key_analysed: map of keyword/value pairs
  QSOMag(keymap& key_analysed);
  QSOMag(){};  // LCOV_EXCL_LINE
  ~QSOMag(){};

  void print_info();
  void read_SED();
  /// Build the library of synthetic magnitudes for a single QSO SED, looping
  /// over the E(B-V) and extinction law grids and the filter set
  /// @param oneSED: the QSO SED to process
  /// @return vector of QSOSED, one per output entry of the library
  vector<QSOSED> make_maglib(const QSOSED& oneSED);
  /// Write the magnitudes of the input SEDs to the output library files
  /// @param seds: SEDs (with their magnitudes already computed) to write out
  void write_mag(const vector<QSOSED>& seds);
};

/// inherited class handling expected magnitudes from galaxy SED
class GalMag : public Mag {
 private:
  vector<double> fracEm;

 public:
  /// Build a GalMag instance and read the galaxy-specific keywords
  /// @param key_analysed: map of keyword/value pairs
  GalMag(keymap& key_analysed);
  GalMag(){};
  ~GalMag(){};

  void print_info();
  void read_SED();
  /// Build the library of synthetic magnitudes for a single galaxy SED,
  /// looping over the E(B-V) and extinction law grids and the filter set.
  /// This also generates the emission-line SED (see GalSED::generateEmSED)
  /// @param oneSED: the galaxy SED to process; its state is modified in
  /// place (emission lines are added to it)
  /// @return vector of GalSED, one per output entry of the library
  vector<GalSED> make_maglib(GalSED& oneSED);
  /// Write the magnitudes of the input SEDs to the output library files
  /// @param seds: SEDs (with their magnitudes already computed) to write out
  void write_mag(const vector<GalSED>& seds);
};

#endif /* MAG_H_ */
