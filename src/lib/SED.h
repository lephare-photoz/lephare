/*
  17/11/2014
  Class to store one SED
*/

// avoid multiple def of the same class
#ifndef SED_H  // check that this keyword has been set already
#define SED_H  // define the keyword to be checked

#include <string>
#include <vector>

#include "ext.h"
#include "flt.h"
#include "globals.h"
#include "oneElLambda.h"
#include "onesource.h"
#include "opa.h"

enum object_type { GAL, QSO, STAR };

/*
 * General SED base class
 */
class SED {
 protected:
  int idAge;  ///< index of the age for this SED object

 public:
  vector<oneElLambda> lamb_flux;  ///< vector of oneElLambda elements, capturing
                                  ///< the SED template
  vector<double> kcorr,           ///< k-correction term
      mag_z0;                     ///< magnitude at z=0
  vector<double> mag;             ///< magnitude and flux of the model
  string name;
  bool has_emlines;  ///< True if the emission lines have been computed, false
                     ///< if not
  int nummod;
  object_type nlib;
  int index, index_z0, posz = -1;
  double red, chi2 = 1.e9, dm, lnir, luv, lopt, inter;
  double mass, age, sfr, ssfr,
      ltir;  // need to put it out of GalSED since used in the PDF without
             // knowing that it's a gal.
  double ebv, mag0, distMod;
  int extlawId;
  double qi[4];  ///< Store the number flux (phot/cm\f$^{-2}\f$s\f$^{-1}\f$) of
                 ///< ionizing photons for HeII, HeI, H, and H2. See
                 ///< SED::calc_ph. In practice, qi[2] only is used, and only
                 ///< for the physical modeling of emission lines
                 ///< (EM_LINES="PHYS", see GalMag::read_SED)
  vector<oneElLambda> fac_line;

  // Constructors defined in SED.cpp
  SED(const string nameC, int nummodC = 0, string typeC = "G");
  SED(const string nameC, double tauC, double ageC, int nummodC, string typeC,
      int idAgeC);
  SED(SED const &p) {
    idAge = p.idAge;
    lamb_flux = p.lamb_flux;
    kcorr = p.kcorr;
    mag = p.mag;
    name = p.name;
    has_emlines = p.has_emlines;
    nummod = p.nummod;
    nlib = p.nlib;
    index = p.index;
    index_z0 = p.index_z0;
    red = p.red;
    chi2 = p.chi2;
    dm = p.dm;
    mass = p.mass;
    age = p.age;
    sfr = p.sfr;
    ltir = p.ltir;
    ebv = p.ebv;
    extlawId = p.extlawId;
    fac_line = p.fac_line;
    distMod = p.distMod;
  };
  bool is_gal() { return nlib == GAL ? true : false; }
  bool is_qso() { return nlib == QSO ? true : false; }
  bool is_star() { return nlib == STAR ? true : false; }
  virtual ~SED();
  ///\brief Read sedFile assumed to be ASCII and build the #lamb_flux vector of
  /// oneElLambda elements.
  ///
  /// Negative flux values are set to 0, oneElLambda::ori=1 to indicate that it
  /// is a SED. The #lamb_flux vector is filled iteratively with each line of
  /// the file, and is finally sorted by ascending lambda. More complex input
  /// types are treated in inherited class methods.
  void read(const string &sedFile);
  void warning_integrateSED(const vector<flt> &filters, bool verbose = false);
  vector<double> integrateSED(const flt &filter);
  void resample(vector<oneElLambda> &lamb_all, vector<oneElLambda> &lamb_new,
                const int origine, const double lmin, const double lmax);

  ///\brief Generate a calibration SED based on the argument calib
  ///
  ///@param lmin start of the lambda vector
  ///@param lmax end of the lambda vector
  ///@param Nsteps number of intervals between $lambda values (hence there are
  /// Nsteps+1 values of \lambda)
  ///@param calib: parameter FILTER_CALIB passed as argument to define the
  /// calibration function \f$C(\lambda)\f$
  /// - calib=0 : \f$C(\lambda)=\lambda^{-2}\f$
  /// - calib=1 : \f$C(\lambda)=\lambda^{-1}\f$
  /// - calib=2 : \f$C(\lambda)=\lambda^{-3}\f$
  /// - calib=3 : \f$C(\lambda)=Blackbody(\lambda, T=10000K)\f$
  /// - calib=4 : \f$C(\lambda)=Blackbody(\lambda, T=10000K)\f$
  /// - calib=5 : \f$C(\lambda)=\lambda^{-3}\f$
  void generateCalib(double lmin, double lmax, int Nsteps, int calib);
  /// return the size of the internal vector #lamb_flux
  int size() { return lamb_flux.size(); }
  /// rescale the lamb_flux.val as val *= scaleFac
  void rescale(double scaleFac);
  /// apply  \a #shifts to the magnitude : mag += shifts
  void applyShift(const vector<double> &shifts, const int imagm);
  /// compute magnitude from filters
  void compute_magnitudes(const vector<flt> &filters);
  /// compute fluxed from filters
  vector<double> compute_fluxes(const vector<flt> &filters);
  double trapzd();
  void sumSpectra(SED addSED, const double rescal);
  void reduce_memory(vector<flt> allFlt);

  /*
   * These functions are different depending on the type of SED
   */
  virtual void writeSED(ofstream &ofs, ofstream &ofsPhys, ofstream &ofsDoc);
  inline void writeSED(const string &binFile, const string &physFile,
                       const string &docFile) {
    ofstream sdocOut, sphysOut, sbinOut;
    sdocOut.open(docFile.c_str());
    if (!sdocOut) {
      throw invalid_argument("Can't open doc file " + docFile);
    }
    sbinOut.open(binFile.c_str(), ios::binary | ios::out);
    if (!sbinOut) {
      throw invalid_argument("Can't open bin file " + binFile);
    }
    if (nlib == GAL) {
      sphysOut.open(physFile.c_str());
      if (!sphysOut) {
        throw invalid_argument("Can't open phys file " + physFile);
      }
    }
    writeSED(sbinOut, sphysOut, sdocOut);
  };
  virtual void writeMag(bool outasc, ofstream &ofsBin, ofstream &ofsDat,
                        vector<flt> allFilters, string magtyp) {};

  inline void readSEDBin(const string &fname) {
    ifstream sbinIn;
    sbinIn.open(fname.c_str(), ios::binary);
    if (!sbinIn) {
      throw invalid_argument("Can't open file " + fname);
    }
    readSEDBin(sbinIn);
  }
  /// read the SED library when it is in binary format
  virtual void readSEDBin(ifstream &ins);
  virtual void readMagBin(ifstream &ins) {};
  virtual void sumEmLines() {};
  /// for each magnitude \a #mag[k] compute kcorr = mag[k] - mag_z0[k] - distMod
  virtual void kcorrec(const vector<double> &magz0) {};
  virtual void add_neb_cont() {};  // Add continuum
  /*!
   * Compute the number flux of photons able to ionize HeII, HeI, H, and H2
   * For a given SED, this amounts to compute the integral
   * \f$\int_0^{w_i} SED(\lambda)\cdot \frac{\lambda}{hc}\,d\lambda\quad,\f$
   * where \f$w_i\f$=54.42, 24.52, 13.60, and 1108.7 A for HeII, HeI, H, and H2
   respectively,
   * and where \f$hc\f$ is in ergs.A. This normalization assumes that the SED
   are provided in args/cm2/s/A.
   * In practice the integral is approximated by :
   \f$\sum_{\lambda_{min}}^{w_k}\frac{SED_{j-1}+SED_j}{2}\cdot(\lambda_j-\lambda_{j-1})\cdot\frac{\lambda_j}{hc}\f$.
   *
   * Results are stored in the q_i array member of size 4 of the SED instance.
   */
  virtual void calc_ph() {};  // Number of inoizing photons
  virtual void SEDproperties() {};
  void generate_spectra(double zin = 0.0, double dmin = 1.0,
                        vector<opa> opaAll = {});

  // in read_lib for zphota
  virtual void clean() {
    lamb_flux.clear();
    mag.clear();
    kcorr.clear();
    fac_line.clear();
  };
  virtual void setOthers() {};
  void fit_normalization(const onesource &source, const int imagm);
  inline bool is_same_model(const SED &other) {
    return ((*this).nummod == other.nummod && (*this).ebv == other.ebv &&
            (*this).age == other.age);
  }
  pair<vector<double>, vector<double>> get_data_vector(double, double, bool,
                                                       double offset = 0.0);

  void redshift();
  void applyExt(const double ebv, const ext &oneext);
  void applyExtLines(const ext &oneext);
  void applyOpa(const vector<opa> &opaAll);

  // helper function for python binding
  inline void emplace_back(const double lambda, const double value) {
    lamb_flux.emplace_back(lambda, value, 1);
  }
  inline void set_vector(const vector<double> &x, const vector<double> &y) {
    if (x.size() != y.size()) throw runtime_error("vector sizes are different");
    for (size_t k = 0; k < x.size(); k++) {
      emplace_back(x[k], y[k]);
    }
  }
};

/// concrete SED implementation for galaxy objects
class GalSED : public SED {
 public:
  vector<double> flEm;
  string format;
  double tau, zmet, d4000, fracEm;

  GalSED(SED const &p) : SED(p) { nlib = GAL; };
  GalSED(GalSED const &p) : SED(p) {
    flEm = p.flEm;
    format = p.format;
    tau = p.tau;
    zmet = p.zmet;
    lnir = p.lnir;
    luv = p.luv;
    lopt = p.lopt;
    d4000 = p.d4000;
    fracEm = p.fracEm;
  };

  // Constructors defined in SED.cpp
  GalSED(const string nameC, int nummodC = 0);
  GalSED(const string nameC, double tauC, double ageC, string formatC,
         int nummodC, string typeC, int idAgeC);
  ~GalSED() { flEm.clear(); }

  void SEDproperties();
  void add_neb_cont();
  GalSED generateEmSED(const string &emtype);
  void generateEmEmpUV(double MNUV_int, double NUVR);
  void generateEmEmpSFR(double MNUV_int, double NUVR);
  void generateEmPhys(double zmet, double qi);
  void generateEmSpectra(int nstep);
  void sumEmLines();
  void kcorrec(const vector<double> &magz0);
  void rescaleEmLines();
  void zdepEmLines(int flag);
  void calc_ph();

  void writeSED(ofstream &ofs, ofstream &ofsPhys, ofstream &ofsDoc);
  void readSEDBin(ifstream &ins);

  void writeMag(bool outasc, ofstream &ofsBin, ofstream &ofsDat,
                vector<flt> allFilters, string magtyp) const;
  void readMagBin(ifstream &ins);
  void clean() {
    SED::clean();
    flEm.clear();
  }
};

/// concrete SED implementation for QSO objects
class QSOSED : public SED {
 public:
  QSOSED(SED const &p) : SED(p) { nlib = QSO; };
  QSOSED(QSOSED const &p) : SED(p){};
  QSOSED(const string nameC, int nummodC = 0) : SED(nameC, nummodC, "QSO"){};
  ~QSOSED(){};

  void writeMag(bool outasc, ofstream &ofsBin, ofstream &ofsDat,
                vector<flt> allFilters, string magtyp) const;

  void readMagBin(ifstream &ins);
};

/// concrete SED implementation for star objects
class StarSED : public SED {
 public:
  StarSED(SED const &p) : SED(p) { nlib = STAR; };
  StarSED(StarSED const &p) : SED(p){};
  StarSED(const string nameC, int nummodC = 0) : SED(nameC, nummodC, "STAR") {
    ;
  }
  ~StarSED() { ; }

  void writeMag(bool outasc, ofstream &ofsBin, ofstream &ofsDat,
                vector<flt> allFilters, string magtyp) const;
  void readMagBin(ifstream &ins);
};

#endif
