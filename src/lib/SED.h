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

//! types of object that LePHARE can treat distinctively
enum object_type {
  GAL, /*!< Galaxy object */
  QSO, /*!< AGN object (QSO naming is historical) */
  STAR /*!< Star object */
};

/*! \brief SED base class
 *
 * The SED class is in charge of representing a template and performing
 * all the necessary computation on it.
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

  //! object type of the SED
  object_type nlib;

  int nummod,    ///< index in the initial list of rest frame SEDs
      index,     ///< index in the full list of SED, redshifted, and modified by
                 ///< extinction, etc...
      index_z0;  ///< index in the full list of SEDs corresponding to the z=0
                 ///< version of the current SED.

  double red,            ///< redshift of this SED
      chi2 = HIGH_CHI2,  ///< best fit chi2 associated with this SED
      dm;                ///< normalization of the SED

  double
      luv,  ///< monochromatic UV luminosity \f$\int_{0.21\,\mu m}^{0.25\,\mu m}
            /// L_{\lambda}\;d\lambda\f$ (in Log unit of erg/s/Hz)
      lopt,  ///< optical luminosity \f$\int_{0.55\,\mu m}^{0.65\,\mu m}
             /// L_{\lambda}\;d\lambda\f$ (in Log unit of erg/s/Hz)
      lnir,  ///< NIR luminosity \f$\int_{2.1\,\mu m}^{2.3\,\mu m}
             /// L_{\lambda}\;d\lambda\f$ (in Log unit of erg/s/Hz)
      ltir;  ///< IR luminosity \f$\int_{8\,\mu m}^{1000\,\mu m}    L_\lambda\;
             ///< d\lambda\f$ in Log unit of \f$L_\odot\f$

  double mass,  ///< mass in \f$M_\odot\f$
      age,      ///< age in year (yr)
      sfr,      ///< Star Formation Rate in \f$M_\odot\f$/yr
      ssfr;     ///< Specific SFR, defined as sfr / mass

  double ebv,  ///< E(B-V) extinction value applied to the SED
      mag0,
      distMod;  ///< Distance modulus of the SED object.

  int extlawId;  ///< index of the extinction law when dust attenuation has been
                 ///< applied

  double qi[4];  ///< Store the number flux (phot/cm\f$^{-2}\f$s\f$^{-1}\f$) of
                 ///< ionizing photons for HeII, HeI, H, and H2. See
                 ///< SED::calc_ph. In practice, qi[2] only is used, and only
                 ///< for the physical modeling of emission lines
                 ///< (EM_LINES="PHYS", see GalMag::read_SED)

  vector<oneElLambda> fac_line;  ///< oneElLambda vector storing emission lines

  /*! Generic constructor
   * \param name Arbitrary name for the SED object
   * \param nummod Identification number of the SED object
   * \param type One of g/G, q/Q or s/S for GAL, QSO, or Star type objects
   */
  SED(const string name, int nummod = 0, string type = "G");
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

  //! Convert string to object_type
  /*!
    \param type String starting with either g, q, or s,
    in either lower or upper case. If it is not the case,
    throw invalid argument exception.

    \return object_type corresponding to input, if valid.
   */
  inline static object_type string_to_object(const string &type) {
    char t = toupper(type[0]);
    if (t == 'S') {
      return STAR;
    } else if (t == 'Q') {
      return QSO;
    } else if (t == 'G') {
      return GAL;
    } else {
      throw invalid_argument("Object type not recognized: " + type);
    }
  }

  //! Return the object type
  object_type get_object_type() const { return nlib; }

  //! Return true if SED is of object_type GAL, false if not
  bool is_gal() { return nlib == GAL ? true : false; }

  //! Return true if SED is of object_type QSO, false if not
  bool is_qso() { return nlib == QSO ? true : false; }

  //! Return true if SED is of object_type STAR, false if not
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
  static void resample(vector<oneElLambda> &lamb_all,
                       vector<oneElLambda> &lamb_new, const int origine,
                       const double lmin, const double lmax);

  ///\brief Generate a calibration SED based on the argument calib
  ///
  ///@param lmin start of the lambda vector
  ///@param lmax end of the lambda vector
  ///@param Nsteps number of intervals between $lambda values (hence there are
  /// Nsteps+1 values of \f$\lambda\f$)
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
      throw invalid_argument("Can't open doc file compiling the SED " +
                             docFile);
    }
    sbinOut.open(binFile.c_str(), ios::binary | ios::out);
    if (!sbinOut) {
      throw invalid_argument("Can't open binary file compiling the SED " +
                             binFile);
    }
    if (nlib == GAL) {
      sphysOut.open(physFile.c_str());
      if (!sphysOut) {
        throw invalid_argument(
            "Can't open physical para file associated to SED " + physFile);
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
      throw invalid_argument("Can't open the binary file compiling the SED " +
                             fname);
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
  virtual void calc_ph() {};

  /*! Compute some integrals to be stored in the object
   * This computes variables SED::luv, SED::lopt, SED::lnir, and SED::ltir
   * luv, lopt, lnir, are monochromatic equivalent luminosities, for a source
   * at 10 parsecs. As the SED unit is taken as erg/cm2/s/Hz, the monochromatic
   * luminosity is obtained by integrating the SED in an interval [lmin, lmax],
   * divided by (lmax-lmin) and multiplied by \f$4\pi(10pc)^2\f$. Note that
   * given the units of an SED, it is defined as dF/dnu for F the corresponding
   * flux. As a result the integral shows a \f$\lambda^2/c\f$ term so that
   * \f$\frac{dF}{d\nu} = \frac{dF}{d\lambda} \frac{\lambda^2}{c}\f$ can be
   * integrated in \f$\lambda\f$. For the variables computed here, in order to
   * speed computation, \f$\lambda^2/c\f$ is evaluated at the center of the
   * interval and taken out of the integral.
   *
   */
  virtual void SEDproperties() {};

  /*! Generate spectrum at given redshift, with given normalization, and
   * adding emission lines and extragalactic extinction
   * \param zin Redshift of the SED
   * \param dmin Scale normalization of the SED
   * \param opaAll Vector of opacities to compute extinction along the line of
   * sight
   */
  void generate_spectra(double zin = 0.0, double dmin = 1.0,
                        vector<opa> opaAll = {});

  ///< clean content of base class
  virtual void clean() {
    lamb_flux.clear();
    mag.clear();
    kcorr.clear();
    fac_line.clear();
  };

  void fit_normalization(const onesource &source, const int imagm);

  /*! Check if the SED is identical to another one
   * @param other: the SED to compare to
   * equality of the following attributes are compared:
   * `nummod`, `ebv`, and `age`
   !*/
  inline bool is_same_model(const SED &other) {
    return ((*this).nummod == other.nummod && (*this).ebv == other.ebv &&
            (*this).age == other.age);
  }

  /*!  Return the pair of vectors [lambdas, vals] of wavelength and spectrum
   * values
   * \param minl: Minimum \f$\lambda\f$ of the vector
   * \param maxl: Maximum \f$\lambda\f$ of the vector
   * \param mag: If true, return magnitudes as vals instead of fluxes
   * \param offset: offset of the mag system to be used in case mag is true.
   */
  pair<vector<double>, vector<double>> get_data_vector(double minl, double maxl,
                                                       bool mag,
                                                       double offset = 0.0);

  /*! Apply transformation of the sed based on the redshift (stored in variable
   * red) \f$val\rightarrow val/(1+z)\f$ and \f$\lambda\rightarrow
   * \lambda(1+z)\f$
   */
  void redshift();

  /*! Apply dust extinction to the SED (GAL and GSO only)
   * \param ebv value of E(B-V)
   * \param obj instance of class ext
   */
  void applyExt(const double ebv, const ext &obj);

  /*! Apply dust extinction to the emission lines (stored in `fac_line`)
   * Only for galaxies and QSO
   * \param obj instance of class `ext`
   */
  void applyExtLines(const ext &obj);

  /*! Apply extinction due to intergalactic medium (only for GAL and QSO)
   * \param opaAll Vector of opacities to compute extinction
   * along the line of sight
   */
  void applyOpa(const vector<opa> &opaAll);

  /// Helper function to append the oneElLambda(lambda, value) object to the sed
  /// vector
  inline void emplace_back(const double lambda, const double value) {
    lamb_flux.emplace_back(lambda, value, 1);
  }

  /*! Helper function to set the sed vector as lambda=x and val = y
   * @param x: vector of lambda value
   * @param y: vector of SED values at each lambda of x
   */
  inline void set_vector(const vector<double> &x, const vector<double> &y) {
    if (x.size() != y.size()) throw runtime_error("vector sizes are different");
    lamb_flux.clear();
    for (size_t k = 0; k < x.size(); k++) {
      emplace_back(x[k], y[k]);
    }
  }
};

/// concrete SED implementation for galaxy objects (object_type GAL)
class GalSED : public SED {
 public:
  vector<double> flEm;
  string format;
  double tau, zmet, d4000,
      fracEm;  //< fraction of the emmission line considered

  /// Copy constructor from base class
  GalSED(SED const &p) : SED(p) { nlib = GAL; };
  /// Copy constructor
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

  /// Standard constructor
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

  /// Compute the k-correction in each filter as :
  /// \f$k = mag(z) - mag(z=0) - \mu\f$
  void kcorrec(const vector<double> &magz0);
  void rescaleEmLines();
  void zdepEmLines(int flag);
  void calc_ph();

  void writeSED(ofstream &ofs, ofstream &ofsPhys, ofstream &ofsDoc);
  void readSEDBin(ifstream &ins);

  void writeMag(bool outasc, ofstream &ofsBin, ofstream &ofsDat,
                vector<flt> allFilters, string magtyp) const;
  void readMagBin(ifstream &ins);

  ///< clean content of class
  void clean() {
    SED::clean();
    flEm.clear();
  }
};

/// concrete SED implementation for AGN/QSO objects (object_type QSO)
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

/// concrete SED implementation for star objects (object_type Star)
class StarSED : public SED {
 public:
  /// copy constructor from `SED` class
  StarSED(SED const &p) : SED(p) { nlib = STAR; };
  /// copy constructor from `StarSED` class
  StarSED(StarSED const &p) : SED(p){};
  /*! constructor
   * @param name: name given to the SED object
   * @param nummod: identity number given to the SED object
   !*/
  StarSED(const string name, int nummod = 0) : SED(name, nummod, "STAR") { ; }
  /// destructor (does nothing)
  ~StarSED() { ; }

  void writeMag(bool outasc, ofstream &ofsBin, ofstream &ofsDat,
                vector<flt> allFilters, string magtyp) const;
  void readMagBin(ifstream &ins);
};

#endif
