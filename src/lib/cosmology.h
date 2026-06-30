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
  double h0, om0, l0;

 public:
  /// define a cosmology based on the triplet \f$H_0\f$, \f$\Omega_m\f$, and
  /// \f$\Omega_\Lambda\f$
  cosmo(double h0 = 70,    ///< Hubble constant
        double om0 = 0.3,  ///< matter density at present time
        double l0 = 0.7)   ///< cosmological constant
  {
    (*this).h0 = h0;
    (*this).om0 = om0;
    (*this).l0 = l0;
  }

  /*!
   * compute the metric distance to z in Mpc, as
   * \f[ d_M(z) = \frac{c}{H_0}\int_0^z \frac{dz}{\sqrt{\Omega_m(1+z)^3 +
   * (1-\Omega_m-\Omega_\Lambda)(1+z)^2+\Omega_\Lambda}} \quad .\f] The
   * luminosity distance of an object at redshift \f$z\f$ is then defined as \f$
   * d_L(z)=(1+z)\,d_M(z)\f$.
   */
  double distMet(double z) const;

  /*!
   * compute the distance modulus at z, as
   * \f[ \mu(z) = 5\log_{10}d_L(z) + 25 \quad. \f]
   */
  double distMod(double z) const;

  double time(
      double z) const;  //!< compute the cosmological time from infinity to z

  //! Identity operator on a cosmo object
  inline bool operator==(const cosmo &rhs) const {
    return h0 == rhs.h0 && om0 == rhs.om0 && l0 == rhs.l0;
  }

  //! Not equal operator on a cosmo object
  inline bool operator!=(const cosmo &rhs) const {
    return h0 != rhs.h0 || om0 != rhs.om0 || l0 != rhs.l0;
  }

  //! Serializer of the cosmo object
  inline friend ostream &operator<<(ostream &output, const cosmo &c) {
    output << c.h0 << "," << c.om0 << "," << c.l0;
    return output;
  }

  //! Return the scaling factor that brings a flux from one value to another,
  //! based only on a change of redshift.
  /*!
    \param z Current redshift
    \param z_t Target redshift

    \return flux scale factor corresponding to a distance change
    from \f$z\f$ to \f$z_t\f$ :
    \f[ scale = 10^{0.4\large(\mu(z_t)-\mu(z)\large)}\f]
-    where \f$\mu\f$ is the distance modulus for a given redshift and is
obtained calling the function cosmo::distMod.
  */
  double flux_rescaling(double z, double z_t) const;
};

vector<double> zgrid(double dz, double zmin, double &zmax);
int indexz(const double red, const vector<double> &gridz);

#endif
