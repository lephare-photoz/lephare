/*
 *  Class with all the tools used to create the grid of SED in e(b-v), redshift,
 * etc and create the library of modeled magnitude and k-corrections
 */

#ifndef MAG_H_
#define MAG_H_

#include "SED.h"
#include "cosmology.h"
#include "ext.h"  // extinction law class
#include "flt.h"
#include "globals.h"
#include "keyword.h"
#include "opa.h"  // IGM opacity class
#ifdef _OPENMP
#include <omp.h>
#endif

/// abstract class to handle everything related to the computation of expected
/// magnitudes
class Mag {
 protected:
  string config, typ;
  cosmo lcdm;
  string filtFile, magtyp;
  bool outasc, verbose;
  vector<string> extlaw;
  int nextlaw;
  vector<double> ebv;
  int nebv;
  vector<int> modext;
  int gridType;
  double dz, zmin, zmax;
  string lib, colib, addDust;

  // only for the galaxy, but much easier to keep them here
  string emlines = "NO";

  string sedlibFile, docFile, binOutFile, datFile;
  ifstream ssedIn;
  ofstream sdocOut, sbinOut, sdatOut;

  vector<ext> extAll;
  vector<flt> allFlt;
  vector<GalSED> B12SED;
  vector<double> gridz;
  vector<double> gridT, gridDM;

  vector<double> magko;

  vector<opa> opaAll;

 public:
  Mag(keymap &key_analysed);
  Mag(){};
  virtual ~Mag();

  /// open the IGM opacity files
  static ifstream open_opa_files();
  /// read the extinction laws into attribute extAll (vector of vectors of type
  /// ext)
  void read_ext();
  /// read the IGM opacities into attribute opaAll (vector of vectors of type
  /// opa)
  static vector<opa> read_opa();
  /// read the filter curve and build the corresponding vectors stored in
  /// attribute allFlt
  void read_flt(const string &);
  // Read the long wavelength Bethermin+2012 templates to add the dust emission
  // to the BC03 templates
  void read_B12();
  /// define the vector of redshifts, and associate to it vectors of age and
  /// distance modulus, based on the lcdm attribute
  void def_zgrid();

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
};

/// inherited class handling expected magnitudes from star SED
class StarMag : public Mag {
 public:
  StarMag(keymap &key_analysed);
  StarMag(){};
  ~StarMag(){};

  void print_info();
  void read_SED();
  vector<StarSED> make_maglib(const StarSED &);
  void write_mag(const vector<StarSED> &);
};

/// inherited class handling expected magnitudes from QSO SED
class QSOMag : public Mag {
 public:
  QSOMag(keymap &key_analysed);
  QSOMag(){};
  ~QSOMag(){};

  void print_info();
  void read_SED();
  vector<QSOSED> make_maglib(const QSOSED &);
  void write_mag(const vector<QSOSED> &);
};

/// inherited class handling expected magnitudes from galaxy SED
class GalMag : public Mag {
 private:
  vector<double> fracEm;

 public:
  GalMag(keymap &key_analysed);
  GalMag(){};
  ~GalMag(){};

  void print_info();
  void read_SED();
  vector<GalSED> make_maglib(GalSED &);
  void write_mag(const vector<GalSED> &);
};

#endif /* MAG_H_ */
