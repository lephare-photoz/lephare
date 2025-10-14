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

// Class cosmology
class cosmo {
 private:
  double h0, om0, l0;

 public:
  // minimal constructor of the ext class, with its name and the id of the model
  cosmo(double h0C = 70,    ///< Hubble constant
        double om0C = 0.3,  ///< matter density at present time
        double l0C = 0.7)   ///< cosmological constant
  {
    h0 = h0C;
    om0 = om0C;
    l0 = l0C;
  }

  double distMod(double z) const;  //!< compute the distance modulus at z
  double distMet(double z) const;  //!< compute the metric distance at z
  double time(
      double z) const;  //!< compute the cosmological time from infinity to z
  inline bool operator==(const cosmo &rhs) {
    return h0 == rhs.h0 && om0 == rhs.om0 && l0 == rhs.l0;
  }
  inline bool operator!=(const cosmo &rhs) {
    return h0 != rhs.h0 || om0 != rhs.om0 || l0 != rhs.l0;
  }
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
    \f$ scale = 10^{0.4(dm(z_t)-dm(z))}\f$
    where \f$dm\f$ is the distance modulus for a given redshift and is obtained
    calling the function cosmo::distMod.
  */
  double flux_rescaling(double z, double z_t) const;
};

vector<double> zgrid(double dz, double zmin, double &zmax);
int indexz(const double red, const vector<double> &gridz);

#endif
