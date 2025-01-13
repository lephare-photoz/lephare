/*!
  16/12/2014
  Class to store the flat LCDM cosmology.
  All the functions to measure the age of the Universe, the luminosity distance,
  etc are here
*/

// avoid multiple def of the same class
#ifndef COSMO_H  // check that this keyword has been set already
#define COSMO_H  // define the keyword to be checked

#include <iostream>
#include <vector>

using std::ostream;
using std::vector;

/*! \brief Very simple flat LCDM-type cosmology class
 *
 * LePHARE needs to define a fiducial cosmology class in order to compute
 * distance moduli, absolute magnitudes, etc... This small class provides the
 * necessary computations.
 */
class cosmo {
 private:
  double _h0, _om0, _l0;

 public:
  /// define a cosmology based on the triplet \f$H_0\f$, \f$\Omega_0\f$, and
  /// \f$\Lambda_0\f$
  cosmo(double h0 = 70,    ///< Hubble constant at present time
        double om0 = 0.3,  ///< matter density at present time
        double l0 = 0.7)   ///< cosmological constant
  {
    _h0 = h0;
    _om0 = om0;
    _l0 = l0;
  }

  double distMod(double z) const;  //!< compute the distance modulus at z
  double distMet(double z) const;  //!< compute the metric distance at z
  double time(
      double z) const;  //!< compute the cosmological time from infinity to z
  inline bool operator==(const cosmo &rhs) {
    return _h0 == rhs._h0 && _om0 == rhs._om0 && _l0 == rhs._l0;
  }
  inline bool operator!=(const cosmo &rhs) {
    return _h0 != rhs._h0 || _om0 != rhs._om0 || _l0 != rhs._l0;
  }
  inline friend ostream &operator<<(ostream &output, const cosmo &c) {
    output << c._h0 << "," << c._om0 << "," << c._l0;
    return output;
  }

  //! Return the scaling factor that brings a flux from one value to another,
  //! based only on a change of redshift.
  /*!
    \param z Current redshift
    \param z_t Target redshift

    \return flux scale factor corresponding to a distance change
    from \f$z\f$ to \f$z_t\f$ :
    \f$ scale = 10^{0.4(dm(z_t)-dm(z))}\f$
    where \f$dm\f$ is the distance modulus for a given redshift and is obtained
    calling the function cosmo::distMod.
  */
  double flux_rescaling(double z, double z_t) const;
};

vector<double> zgrid(int gridType, double dz, double zmin, double &zmax);
int indexz(const double red, const vector<double> &gridz);

#endif
