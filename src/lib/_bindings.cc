#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <string>
#include <vector>

#include "SED.h"
#include "SEDLib.h"
#include "cosmology.h"
#include "ext.h"
#include "flt.h"
#include "globals.h"
#include "keyword.h"
#include "mag.h"
#include "oneElLambda.h"
#include "opa.h"
#include "photoz_lib.h"

template <typename x, typename modT>
void applySEDLibTemplate(modT& m, std::string name) {
  py::class_<SEDLib<x>>(m, name.c_str(), R"pbdoc(
Library class for SED objects of type `{type}`.

This class provides a container for multiple SED objects (e.g., StarSED, QSOSED, GalSED),
with methods to load SEDs, write them, print info, and manage output files.

Parameters
----------
config : str
    Path to configuration file used to define the SED library.
typ : str
    Type identifier for the SEDs (e.g., 'star', 'gal', 'qso').

key_analysed : keymap, optional
    Keymap used for advanced configuration, only for the second constructor.

Methods
-------
print_info()
    Print a summary of the SED library to standard output.

read_model_list()
    Read the list of SED models defined in the library.

write_SED_lib()
    Write the SED library to output files.

print_time_tofile()
    Print timing information to an output file.

close_output_files()
    Close all open output files associated with the SED library.

Usage Example
-------------
>>> import lephare as lp
>>> starlib = lp.StarSEDLib("config.txt", "star")
>>> starlib.print_info()
>>> starlib.read_model_list()
>>> starlib.write_SED_lib()
)pbdoc")
      .def(py::init<std::string, std::string>(), py::arg("config"),
           py::arg("typ"), "Constructor with configuration file and type")
      .def(py::init<keymap&, std::string, std::string>(),
           py::arg("key_analysed"), py::arg("config"), py::arg("typ"),
           "Constructor with keymap, configuration file, and type")
      .def("print_info", &SEDLib<x>::print_info,
           "Print summary of the SED library")
      .def("read_model_list", &SEDLib<x>::read_model_list,
           "Read the list of SED models")
      .def("write_SED_lib", &SEDLib<x>::write_SED_lib,
           "Write the SED library to output files")
      .def("print_time_tofile", &SEDLib<x>::print_time_tofile,
           "Print timing info to a file")
      .def("close_output_files", &SEDLib<x>::close_output_files,
           "Close all output files");
}

PYBIND11_MODULE(_lephare, mod) {
  /* Top-level module docstring. Docs for individual classes are done in
   * __init__.py */
  mod.doc() = R"pbdoc(
    _lephare: Core C++ bindings for the LePHARE Python interface

    This module is the low-level, compiled part of the LePHARE Python package.
    It exposes the C++ classes and functions to Python via pybind11. Users
    generally do not import _lephare directly; instead, use the high-level
    `lephare` package which wraps _lephare and provides convenience functions
    and a cleaner Python interface.

    Key features exposed in this module:
    - Classes: PhotoZ, SED, StarSED, GalSED, QSOSED, flt, cosmo, PDF, onesource
    - Functions: flux/magnitude conversions, extinction laws, library handling
    - Constants: HIGH_CHI2, maptype

    Example usage (high-level):
        import lephare as lp
        lp.data_retrieval.get_auxiliary_data(clone=True)
        from astropy.table import Table
        config = lp.default_cosmos_config.copy()
        lp.prepare(config)
        input_table = Table.read(f"{lp.LEPHAREDIR}/examples/COSMOS.in", format='ascii')
        output, _ = lp.process(config, input_table)

    Documentation:
        Full LePHARE Python documentation is available at:
        https://lephare.readthedocs.io/
    )pbdoc";

  /*object_type enum for python*/
  py::enum_<object_type>(mod, "object_type")
      .value("GAL", object_type::GAL)
      .value("QSO", object_type::QSO)
      .value("STAR", object_type::STAR);

  /******** CLASS ONEELLAMBDA *********/
  py::class_<oneElLambda>(mod, "oneElLambda")
      .def(py::init<double, double, int>(), py::arg("lambin"), py::arg("valin"),
           py::arg("oriin"), "Standard constructor")
      .def(py::init<oneElLambda>(), py::arg("elIn"), "Copy constructor")
      .def_readwrite("lamb", &oneElLambda::lamb, "Wavelength of the element")
      .def_readwrite("val", &oneElLambda::val,
                     "Flux or magnitude value of the element")
      .def_readwrite("ori", &oneElLambda::ori, "Origin or index identifier")
      .def("interp", &oneElLambda::interp, py::arg("previousEl"),
           py::arg("nextEl"),
           R"pbdoc(
Interpolates the value of this element between two neighboring elements.

Parameters
----------
previousEl : oneElLambda
    The previous element to interpolate from.
nextEl : oneElLambda
    The next element to interpolate to.

Returns
-------
float
    The interpolated value for this element.
)pbdoc");

  /******** CLASS COSMOLOGY *********/
  py::class_<cosmo>(mod, "cosmo")
      .def(py::init<double, double, double>(), py::arg("h0") = 70,
           py::arg("om0") = 0.3, py::arg("l0") = 0.7,
           "Standard constructor with optional cosmological parameters.")
      .def("distMod", py::vectorize(&cosmo::distMod),
           R"pbdoc(Compute distance modulus given a redshift.

Parameters
----------
z : float or array_like
    Redshift(s) at which to compute the distance modulus.

Returns
-------
float or ndarray
    Distance modulus corresponding to the input redshift(s).
)pbdoc",
           py::arg("z"))
      .def("distMet", py::vectorize(&cosmo::distMet),
           R"pbdoc(Compute metric (comoving) distance given a redshift.

Parameters
----------
z : float or array_like
    Redshift(s) at which to compute the metric distance.

Returns
-------
float or ndarray
    Metric distance in Mpc.
)pbdoc",
           py::arg("z"))
      .def("time", py::vectorize(&cosmo::time),
           R"pbdoc(Compute cosmic time at a given redshift.

Parameters
----------
z : float or array_like
    Redshift(s) at which to compute the cosmic time.

Returns
-------
float or ndarray
    Cosmic time in Gyr corresponding to redshift(s) `z`.
)pbdoc",
           py::arg("z"))
      .def("flux_rescaling", &cosmo::flux_rescaling,
           R"pbdoc(Compute the flux rescaling factor for this cosmology.

Returns
-------
float
    The multiplicative factor to rescale fluxes according to this cosmology.
)pbdoc")
      .def(py::self == py::self,
           "Equality comparison between two cosmology objects.")
      .def(py::self != py::self,
           "Inequality comparison between two cosmology objects.");

  /******** MODULE-LEVEL FUNCTIONS *********/
  mod.def("zgrid", &zgrid,
          R"pbdoc(Generate a redshift grid.

Returns
-------
ndarray
    Array of redshift values forming the grid.
)pbdoc");

  mod.def("indexz", &indexz,
          R"pbdoc(Get the index corresponding to a redshift value in a grid.

Parameters
----------
z : float
    Redshift value to find the index for.

Returns
-------
int
    Index of `z` in the redshift grid.
)pbdoc");

  /******** CLASS OPA *********/
  py::class_<opa>(mod, "opa")
      .def(py::init<double, std::string>(), py::arg("red"), py::arg("opaFile"),
           "Standard constructor with redshift and opacity file path.")
      .def_readwrite("lamb_opa", &opa::lamb_opa,
                     "Wavelength array from the opacity file.")
      .def_readwrite("opared", &opa::red, "Redshift value of the opacity data.")
      .def("read", &opa::read,
           R"pbdoc(Read the opacity data from the file.

This method populates the `lamb_opa` and any other internal data arrays
from the file specified during object creation.
)pbdoc");
  //   .def("lmin", &opa::lmin, "return smallest wavelength stored")
  //   .def("lmax", &opa::lmax, "return largest wavelength stored")

  /******** CLASS EXT *********/
  py::class_<ext>(mod, "ext")
      .def(py::init<std::string, int>(), py::arg("name"), py::arg("numext"),
           "Standard constructor with name and number of elements.")
      .def_readwrite("lamb_ext", &ext::lamb_ext,
                     "Wavelength array for the extinction curve.")
      .def_readwrite("name", &ext::name, "Name of the extinction curve.")
      .def_readwrite("numext", &ext::numext,
                     "Number of elements in the extinction curve.")
      .def_readonly("lmin", &ext::lmin, "Return smallest wavelength stored.")
      .def_readonly("lmax", &ext::lmax, "Return largest wavelength stored.")
      .def("read", &ext::read, py::arg("extFile"),
           "Read an extinction file and populate the extinction curve.")
      .def("add_element", &ext::add_element,
           "Add a new element to the extinction curve.");

  /******** MODULE-LEVEL FUNCTIONS *********/
  mod.def("compute_filter_extinction", &compute_filter_extinction,
          R"pbdoc(Compute extinction in a filter band.

Parameters
----------
filter_band : object
    Filter object or array specifying the band.

Returns
-------
float or ndarray
    Extinction in the filter band.
)pbdoc");

  mod.def(
      "cardelli_ext", &cardelli_ext,
      R"pbdoc(Compute Galactic extinction in a filter based on Cardelli et al., 1989 (ApJ 345).

Parameters
----------
oneFlt : object
    Filter object to compute extinction for.

Returns
-------
float
    Galactic extinction in the given filter.
)pbdoc",
      py::arg("oneFlt"));

  mod.def(
      "cardelli_law", &cardelli_law,
      R"pbdoc(Compute A_lambda / A_V at a given wavelength using the Cardelli law.

Parameters
----------
lb : float
    Wavelength in Angstroms.

Returns
-------
float
    Extinction ratio A_lambda / A_V at the given wavelength.
)pbdoc",
      py::arg("lb"));

  mod.def(
      "resample", &resample, py::arg("lamb_all"), py::arg("lamb_interp"),
      py::arg("origine"), py::arg("lmin"), py::arg("lmax"),
      R"pbdoc(Resample a wavelength-dependent quantity onto a new wavelength grid.

Parameters
----------
lamb_all : ndarray
    Original wavelength array.
lamb_interp : ndarray
    Wavelength array to interpolate onto.
origine : int
    Origin index for interpolation.
lmin : float
    Minimum wavelength for interpolation.
lmax : float
    Maximum wavelength for interpolation.

Returns
-------
ndarray
    Interpolated values on the new wavelength grid.
)pbdoc");

  mod.def("read_flt", &read_flt, py::arg("sfiltIn"),
          R"pbdoc(Read a filter definition file.

Parameters
----------
sfiltIn : str
    Path to the filter file.

Returns
-------
Filter object
    Object representing the filter.
)pbdoc");

  /******** CLASS KEYWORD *********/
  py::class_<keyword>(mod, "keyword")
      .def_readwrite("name", &keyword::name, "Name of the keyword")
      .def_readwrite("value", &keyword::value,
                     "Value of the keyword as a string")
      .def(py::init(), "Default constructor")
      .def(py::init<std::string, std::string>(), py::arg("n"), py::arg("v"),
           "Constructor with name and value")
      .def("expand_path", &keyword::expand_path,
           "Expand environment variables and user home (`~`) in value")
      .def("split_string", &keyword::split_string,
           "Split value into a list of strings")
      .def("split_int", &keyword::split_int,
           "Split value into a list of integers")
      .def("split_long", &keyword::split_long,
           "Split value into a list of long integers")
      .def("split_double", &keyword::split_double,
           "Split value into a list of doubles/floats")
      .def("split_bool", &keyword::split_bool,
           "Split value into a list of booleans")
      .def(
          "__repr__",
          [](const keyword& a) { return "(" + a.name + ", " + a.value + ")"; },
          "String representation of the keyword as (name, value)");

  /******** MODULE-LEVEL FUNCTIONS *********/
  mod.def(
      "read_command",
      [](std::vector<std::string> args) {
        std::vector<char*> cstrs;
        cstrs.reserve(args.size());
        for (auto& s : args) cstrs.push_back(const_cast<char*>(s.c_str()));
        return read_command(cstrs.size(), cstrs.data());
      },
      R"pbdoc(
Read a command line and parse it into a configuration dictionary.

Parameters
----------
args : list of str
    List of command-line arguments.

Returns
-------
dict
    Dictionary of keywords and their values.
)pbdoc");

  mod.def("read_config", &read_config, R"pbdoc(
Read a configuration file and parse it into a dictionary of keywords.

Parameters
----------
filename : str
    Path to the configuration file.

Returns
-------
dict
    Dictionary of keywords (name-value pairs) from the file.
)pbdoc");

  /******** CLASS FLT *********/
  py::class_<flt>(mod, "flt", py::dynamic_attr())
      .def(py::init<int, std::string, int, int>(), py::arg("id"),
           py::arg("name"), py::arg("trans"), py::arg("calib"),
           "Constructor with ID, name, transmission type, and calibration type")
      .def(py::init<double, double, int>(), py::arg("lmin"), py::arg("lmax"),
           py::arg("nstep"),
           "Top-hat filter from lmin to lmax with nstep points")
      .def("read", static_cast<void (flt::*)(const std::string&)>(&flt::read),
           "Read filter info from file")
      .def("read", static_cast<void (flt::*)(std::ifstream&)>(&flt::read),
           "Read filter info from input stream")
      .def("lambdaMean", &flt::lambdaMean,
           "Compute mean wavelength of the filter")
      .def("lambdaEff", &flt::lambdaEff,
           "Compute effective wavelength of the filter")
      .def("magsun", &flt::magsun,
           "Compute magnitude of the Sun in this filter")
      .def("width", &flt::width, "Compute the width of the filter")
      .def("lmin", &flt::lmin, "Return the minimum wavelength of the filter")
      .def("lmax", &flt::lmax, "Return the maximum wavelength of the filter")
      .def_readonly("name", &flt::name, "Name of the filter")
      .def_readonly("lmean", &flt::lmean, "Mean wavelength of the filter")
      .def_readonly("dwidth", &flt::dwidth, "Width of the filter")
      .def_readwrite("lamb_trans", &flt::lamb_trans,
                     "Transmission data points (wavelength + transmission)")
      .def(
          "data",
          [](const flt& f) {
            int N = f.lamb_trans.size();
            py::array_t<double> result({2, N});
            py::buffer_info buf = result.request();
            double* ptr = static_cast<double*>(buf.ptr);
            for (size_t i = 0; i < N; i++) {
              ptr[i] = f.lamb_trans[i].lamb;     // first row: wavelengths
              ptr[N + i] = f.lamb_trans[i].val;  // second row: transmission
            }
            return result;
          },
          "Return 2xN array of wavelengths and transmission values");

  /******** MODULE-LEVEL FUNCTIONS *********/
  mod.def("write_output_filter", &write_output_filter,
          R"pbdoc(Write filter data to an output file.

Parameters
----------
flt_obj : flt
    Filter object to write.

filename : str
    Path to the output file.
)pbdoc");

  mod.def("read_doc_filters", &read_doc_filters,
          R"pbdoc(Read a document listing filters.

Parameters
----------
filename : str
    Path to the filter documentation file.

Returns
-------
list of flt
    List of filter objects read from the file.
)pbdoc");

  /******** CLASS SED *********/
  py::class_<SED>(mod, "SED")
      .def(py::init<const std::string, int, std::string>(),
           py::arg("name") = "", py::arg("nummod") = 0, py::arg("type") = "G",
           "Create an empty SED with optional name, number of models, and type")
      .def(py::init<const std::string, double, double, int, std::string, int>(),
           py::arg("name"), py::arg("tau"), py::arg("age"), py::arg("nummod"),
           py::arg("type"), py::arg("idAge"),
           "Create SED with physical parameters: tau, age, nummod, type, idAge")
      .def(py::init<const SED>(), "Copy constructor")
      .def_readonly("lamb_flux", &SED::lamb_flux,
                    "Wavelength and flux values of the SED")
      .def_readonly("extlawId", &SED::extlawId,
                    "Extinction law ID applied to this SED")
      .def_readonly("ebv", &SED::ebv, "Color excess E(B-V) for extinction")
      .def_readonly("name", &SED::name, "Name of the SED")
      .def_readonly("nummod", &SED::nummod, "Number of models")
      .def_readonly("mag", &SED::mag, "Magnitudes computed for this SED")
      .def_readwrite("index_z0", &SED::index_z0, "Index for z=0 reference")
      .def("string_to_object", &SED::string_to_object,
           "Convert string to SED object")
      .def("is_gal", &SED::is_gal, "Return True if SED is a galaxy")
      .def("is_star", &SED::is_star, "Return True if SED is a star")
      .def("is_qso", &SED::is_qso, "Return True if SED is a quasar")
      .def("read", &SED::read, "Read SED data from file")
      .def("size", &SED::size, "Return number of wavelength points in SED")
      .def("integrateSED", &SED::integrateSED,
           "Integrate the SED over wavelength")
      .def("resample", &SED::resample,
           "Resample the SED onto a new wavelength grid")
      .def("generateCalib", &SED::generateCalib,
           "Generate calibration vector for the SED")
      .def("rescale", &SED::rescale, "Rescale the SED fluxes by a factor")
      .def("compute_magnitudes", &SED::compute_magnitudes,
           "Compute magnitudes in filters")
      .def("compute_fluxes", &SED::compute_fluxes, "Compute fluxes in filters")
      .def("emplace_back", &SED::emplace_back, "Add a new element to the SED")
      .def("set_vector", &SED::set_vector,
           "Replace the SED vector with new values")
      .def("readSEDBin",
           static_cast<void (SED::*)(const std::string&)>(&SED::readSEDBin),
           "Read a binary SED file")
      .def("writeSED",
           static_cast<void (SED::*)(const std::string&, const std::string&,
                                     const std::string&)>(&SED::writeSED),
           "Write the SED to a file in the specified format and options")
      .def(
          "data",
          [](const SED& f) {
            int N = f.lamb_flux.size();
            py::array_t<double> result({2, N});
            py::buffer_info buf = result.request();
            double* ptr = static_cast<double*>(buf.ptr);
            for (size_t i = 0; i < N; i++) {
              ptr[i] = f.lamb_flux[i].lamb;     // first row: wavelengths
              ptr[N + i] = f.lamb_flux[i].val;  // second row: fluxes
            }
            return result;
          },
          "Return a 2xN numpy array: first row = wavelengths, second row = "
          "fluxes");

  /******** CLASS STARS, QSOS, GALAXIES *********/
  py::class_<StarSED, SED>(mod, "StarSED")
      .def(py::init<const SED&>(), "Create StarSED from existing SED")
      .def(py::init<const StarSED&>(), "Copy constructor")
      .def(py::init<const std::string, int>(), py::arg("name"),
           py::arg("nummod") = 0,
           "Create a new StarSED with optional name and number of models");

  py::class_<QSOSED, SED>(mod, "QSOSED")
      .def(py::init<const SED&>(), "Create QSOSED from existing SED")
      .def(py::init<const QSOSED&>(), "Copy constructor")
      .def(py::init<const std::string, int>(), py::arg("name"),
           py::arg("nummod") = 0,
           "Create a new QSOSED with optional name and number of models");

  py::class_<GalSED, SED>(mod, "GalSED")
      .def(py::init<const SED&>(), "Create GalSED from existing SED")
      .def(py::init<const GalSED&>(), "Copy constructor")
      .def(py::init<const std::string, int>(), py::arg("name"),
           py::arg("nummod") = 0,
           "Create a new GalSED with optional name and number of models")
      .def(py::init<const std::string, double, double, std::string, int,
                    std::string, int>(),
           py::arg("name"), py::arg("tau"), py::arg("age"), py::arg("format"),
           py::arg("nummod"), py::arg("type"), py::arg("idAge"),
           "Create a new GalSED with physical parameters, format, type, and "
           "age ID")
      .def("SEDproperties", &GalSED::SEDproperties,
           "Compute or return SED properties")
      .def("add_neb_cont", &GalSED::add_neb_cont,
           "Add nebular continuum emission")
      .def("generateEmEmpUV", &GalSED::generateEmEmpUV,
           "Generate empirical UV emission")
      .def("generateEmEmpSFR", &GalSED::generateEmEmpSFR,
           "Generate empirical SFR emission")
      .def("generateEmPhys", &GalSED::generateEmPhys,
           "Generate emission lines from physical models")
      .def("sumEmLines", &GalSED::sumEmLines,
           "Sum all emission lines into the SED")
      .def("kcorrec", &GalSED::kcorrec, "Compute K-corrections")
      .def("rescaleEmLines", &GalSED::rescaleEmLines,
           "Rescale emission lines by a factor")
      .def("zdepEmLines", &GalSED::zdepEmLines,
           "Compute redshift-dependent emission lines")
      .def("calc_ph", &GalSED::calc_ph,
           "Compute photometric properties from the SED");

  /******** CLASS SEDLib *********/
  applySEDLibTemplate<StarSED>(mod, "StarSEDLib");
  applySEDLibTemplate<QSOSED>(mod, "QSOSEDLib");
  applySEDLibTemplate<GalSED>(mod, "GalSEDLib");

  /******** CLASS MAG *********/
  // Generic docstring for all magnitude classes
  const char* mag_doc =
      "Class for computing magnitudes for stars, QSOs, or galaxies.\n"
      "This class handles reading SEDs, filters, extinction, and generating "
      "magnitudes.\n"
      "\n"
      "Constructors:\n"
      "Mag(key_analysed: keymap)  -> Initialize with a keymap.\n"
      "Mag()                     -> Default constructor.\n"
      "\n"
      "Methods:\n"
      "open_files()        -> Open all necessary input/output files.\n"
      "close_files()       -> Close all open files.\n"
      "open_opa_files()    -> Open opacity-related files.\n"
      "print_info()        -> Print summary information.\n"
      "read_ext()          -> Read extinction data.\n"
      "read_opa()          -> Read opacity data.\n"
      "read_B12()          -> Read B12 calibration data.\n"
      "read_flt()          -> Read filter data.\n"
      "def_zgrid()         -> Define the redshift grid.\n"
      "set_zgrid()         -> Set the redshift grid.\n"
      "read_SED()          -> Read SEDs for this object type.\n"
      "write_doc()         -> Write documentation files.\n"
      "make_maglib()       -> Generate the magnitude library.\n"
      "write_mag()         -> Write computed magnitudes to output files.\n"
      "\n"
      "Example usage:\n"
      ">>> import lephare as lp\n"
      ">>> mag = lp.StarMag()  # or QSOMag(), GalMag()\n"
      ">>> mag.open_files()\n"
      ">>> mag.read_SED()\n"
      ">>> mag.make_maglib()\n"
      ">>> mag.write_mag()\n";

// Macro to bind magnitude classes using generic docstring
#define MAGDEFS(c, n)                                                     \
  (py::class_<c>(mod, n, mag_doc)                                         \
       .def(py::init<keymap&>(), py::arg("key_analysed"),                 \
            "Constructor with keymap")                                    \
       .def(py::init<>(), "Default constructor")                          \
       .def("open_files", &c::open_files, "Open input/output files")      \
       .def("close_files", &c::close_files, "Close all files")            \
       .def("open_opa_files", &c::open_opa_files, "Open opacity files")   \
       .def("print_info", &c::print_info, "Print summary info")           \
       .def("read_ext", &c::read_ext, "Read extinction data")             \
       .def("read_opa", &c::read_opa, "Read opacity data")                \
       .def("read_B12", &c::read_B12, "Read B12 calibration")             \
       .def("read_flt", &c::read_flt, "Read filter data")                 \
       .def("def_zgrid", &c::def_zgrid, "Define the redshift grid")       \
       .def("set_zgrid", &c::set_zgrid, "Set the redshift grid")          \
       .def("read_SED", &c::read_SED, "Read SEDs for this object type")   \
       .def("write_doc", &c::write_doc, "Write documentation files")      \
       .def("make_maglib", &c::make_maglib, "Generate magnitude library") \
       .def("write_mag", &c::write_mag, "Write computed magnitudes"));

  // Instantiate classes
  MAGDEFS(StarMag, "StarMag");
  MAGDEFS(QSOMag, "QSOMag");
  MAGDEFS(GalMag, "GalMag");

  /******** GLOBAL FUNCTIONS AND CONSTANTS *********/
  // Module-level constants
  mod.attr("HIGH_CHI2") =
      HIGH_CHI2;  // Threshold for chi-squared used in fitting

  // Module-level functions
  mod.def("get_lephare_env", &get_lephare_env,
          R"pbdoc(
Retrieve the LePHARE environment configuration.

Returns
-------
str
    Path or environment variable indicating the LePHARE installation or data location.

Example
-------
>>> lp.get_lephare_env()
'/path/to/lephare/data'
)pbdoc");

  mod.def("check_first_char", &check_first_char,
          R"pbdoc(
Check the first character of a string.

Parameters
----------
s : str
    Input string.

Returns
-------
str
    The first character of the string.

Example
-------
>>> lp.check_first_char('hello')
'h'
)pbdoc");

  mod.def("blackbody", &blackbody,
          R"pbdoc(
Compute the blackbody flux at a given wavelength and temperature.

Parameters
----------
wavelength : float
    Wavelength in meters.
temperature : float
    Temperature in Kelvin.

Returns
-------
float
    Blackbody flux at the given wavelength.

Example
-------
>>> lp.blackbody(500e-9, 5778)
3.1e13
)pbdoc");

  // Utility macros / functions for internal calculations
  mod.def("CHECK_CONTEXT_BIT", &CHECK_CONTEXT_BIT,
          "Check a bit in a context variable (low-level internal utility)");

  mod.def("POW10D", &POW10D, "Compute 10^x (double precision, fast)");

  mod.def("POW10DSLOW", &POW10DSLOW,
          "Compute 10^x (double precision, slower, more accurate)");

  mod.def("LOG10D_SLOW", &LOG10D_SLOW,
          "Compute log10(x) in double precision, slow method");

  mod.def("LOG10D_FAST", &LOG10D_FAST,
          "Compute log10(x) in double precision, fast method");

  // Magnitude / flux conversions
  mod.def("mag2flux", &mag2flux,
          R"pbdoc(
Convert a magnitude to a flux.

Parameters
----------
mag : float
    Magnitude value.
zp : float, optional
    Zero-point (default 0.0).

Returns
-------
float
    Corresponding flux.

Example
-------
>>> lp.mag2flux(25.0)
3.63e-30
)pbdoc");

  mod.def("flux2mag", &flux2mag,
          R"pbdoc(
Convert a flux to a magnitude.

Parameters
----------
flux : float
    Flux value.
zp : float, optional
    Zero-point (default 0.0).

Returns
-------
float
    Corresponding magnitude.

Example
-------
>>> lp.flux2mag(3.63e-30)
25.0
)pbdoc");

  // Vector utilities
  mod.def("indexes_in_vec", &indexes_in_vec,
          R"pbdoc(
Find the indexes of elements in a vector.

Parameters
----------
vec : list or array-like
    Input vector.
value : scalar
    Value to find in vec.

Returns
-------
list of int
    List of indexes where the value occurs in vec.

Example
-------
>>> lp.indexes_in_vec([10, 20, 10], 10)
[0, 2]
)pbdoc");

  /******** CLASS PHOTOZ *********/
  py::class_<PhotoZ>(mod, "PhotoZ")
      .def_readonly("flux", &PhotoZ::flux)
      .def_readonly("fluxIR", &PhotoZ::fluxIR)
      .def_readonly("imagm", &PhotoZ::imagm)
      .def_readonly("fullLib", &PhotoZ::fullLib)
      .def_readonly("zLib", &PhotoZ::zLib)
      .def_readonly("fullLibIR", &PhotoZ::fullLibIR)
      .def_readonly("allFilters", &PhotoZ::allFilters)
      .def_readonly("gridz", &PhotoZ::gridz)
      .def_readonly("outkeywords", &PhotoZ::outkeywords)
      .def_readonly("outpara", &PhotoZ::outpara)
      .def_readonly("pdftype", &PhotoZ::pdftype)
      .def_readwrite("outputHeader", &PhotoZ::outputHeader)
      .def(py::init<keymap&>(), "Constructor with keymap configuration")
      .def("read_autoadapt_sources", &PhotoZ::read_autoadapt_sources,
           "Read sources prepared for auto-adaptive photo-z computation")
      .def("read_photoz_sources", &PhotoZ::read_photoz_sources,
           "Read sources for standard photo-z computation")
      .def(
          "prep_data",
          static_cast<void (PhotoZ::*)(vector<onesource*>)>(&PhotoZ::prep_data),
          "Prepare a vector of sources for processing")
      .def("prep_data",
           static_cast<void (PhotoZ::*)(onesource*)>(&PhotoZ::prep_data),
           "Prepare a single source for processing")
      .def("run_autoadapt", &PhotoZ::run_autoadapt,
           "Run the auto-adaptive photo-z workflow")
      .def("run_photoz", &PhotoZ::run_photoz,
           "Run the standard photo-z computation")
      .def("write_outputs", &PhotoZ::write_outputs,
           "Write computed outputs to files")
      .def("validLib", &PhotoZ::validLib, "Validate the photo-z library")
      .def("compute_offsets", &PhotoZ::compute_offsets,
           "Compute offsets between observed and template fluxes");

  // Module-level utility functions related to PhotoZ
  mod.def("readOutKeywords", &readOutKeywords,
          "Read output keywords for photo-z analysis");

  mod.def("bestFilter", &bestFilter,
          "Identify the best filter for a given source or SED");

  mod.def("maxkcolor", &maxkcolor,
          "Compute maximum k-correction in color space for a source");

  // mod.def("read_lib", [](const string& libName, int ind, vector<int> emMod,
  // int babs) { 			vector<SED*> libFull;
  // int nummodpre[3]; 			string filtname;
  // vector<double> gridz; 			nummodpre[0] = 0;
  // nummodpre[1] = 0; 			nummodpre[2] = 0;
  // read_lib(libFull, ind, nummodpre, libName,
  // filtname, gridz, emMod, babs); 			std::array<int, 3>
  // nummod_tup = {nummodpre[0], nummodpre[1], nummodpre[2]};
  // return std::make_tuple(libFull, ind, nummod_tup, filtname, gridz);
  // 		      }
  // 	  );
  // mod.def("read_doc_filters", [](const string filtFile) {
  // 				bool Fexiste;
  // 				vector<flt> allFilters =
  // read_doc_filters(filtFile, Fexiste); 				return
  // std::make_tuple(allFilters, Fexiste);
  // 			      }
  // );

  /******** CLASS ONESOURCE *********/
  mod.attr("maptype") =
      maptype;  // Global module-level variable indicating PDF map type

  py::class_<onesource>(mod, "onesource")
      .def(py::init<>())
      .def(py::init<const int, vector<double>>())
      .def("setPriors", &onesource::setPriors)
      .def_readonly("priorLib", &onesource::priorLib)
      .def("readsource",
           static_cast<void (onesource::*)(
               const string&, const vector<double>, const vector<double>,
               const long, const double, const string)>(&onesource::readsource))
      .def("set_verbosity", &onesource::set_verbosity)
      .def("get_verbosity", &onesource::get_verbosity)
      .def("fltUsed", &onesource::fltUsed)
      .def("convertFlux", &onesource::convertFlux)
      .def("convertMag", &onesource::convertMag)
      .def("rescale_flux_errors", &onesource::rescale_flux_errors)
      .def("keepOri", &onesource::keepOri)
      .def("adapt_mag", &onesource::adapt_mag)
      .def("fit", &onesource::fit)
      .def("mode", &onesource::mode)
      .def("rm_discrepant", &onesource::rm_discrepant)
      .def("generatePDF", &onesource::generatePDF)
      .def("interp", &onesource::interp)
      .def("uncertaintiesMin", &onesource::uncertaintiesMin)
      .def("uncertaintiesBay", &onesource::uncertaintiesBay)
      .def("secondpeak", &onesource::secondpeak)
      .def("interp_lib", &onesource::interp_lib)
      .def("absmag", &onesource::absmag)
      .def("limits", &onesource::limits)
      .def("computePredAbsMag", &onesource::computePredAbsMag)
      .def("computePredMag", &onesource::computePredMag)
      .def("computeEmFlux", &onesource::computeEmFlux)
      .def("fltUsedIR", &onesource::fltUsedIR)
      .def("substellar", &onesource::substellar)
      .def("generatePDF_IR", &onesource::generatePDF_IR)
      .def("write_out", &onesource::write_out)
      .def("writeSpec", &onesource::writeSpec)
      .def("writeFullChi", &onesource::writeFullChi)
      .def("best_spec_vec", &onesource::best_spec_vec)
      .def_readwrite("spec", &onesource::spec)
      .def_readwrite("consiz", &onesource::consiz)
      .def_readonly("pos", &onesource::pos)
      .def_readonly("cont", &onesource::cont)
      .def_readonly("pdfmap", &onesource::pdfmap)
      .def_readonly("busnorma", &onesource::busnorma)
      .def_readonly("busul", &onesource::busul)
      .def_readonly("nbused", &onesource::nbused)
      .def_readonly("nbul", &onesource::nbul)
      .def_readonly("dm", &onesource::dm)
      .def_readonly("zs", &onesource::zs)
      .def_readonly("ab", &onesource::ab)
      .def_readonly("ab_ori", &onesource::ab_ori)
      .def_readonly("sab", &onesource::sab)
      .def_readonly("mab", &onesource::mab)
      .def_readonly("msab", &onesource::msab)
      .def_readonly("magm", &onesource::magm)
      .def_readonly("magPred", &onesource::magPred)
      .def_readonly("absmagPred", &onesource::absmagPred)
      .def_readonly("zmin", &onesource::zmin)
      .def_readonly("zminIR", &onesource::zminIR)
      .def_readonly("chimin", &onesource::chimin)
      .def_readonly("dmmin", &onesource::dmmin)
      .def_readonly("chiminIR", &onesource::chiminIR)
      .def_readonly("indmin", &onesource::indmin)
      .def_readonly("indminSec", &onesource::indminSec)
      .def_readonly("indminIR", &onesource::indminIR)
      .def_readonly("imasmin", &onesource::imasmin)
      .def_readonly("imasminIR", &onesource::imasminIR)
      .def_readonly("agemed", &onesource::agemed)
      .def_readonly("Ldustmed", &onesource::Ldustmed)
      .def_readonly("LIRmed", &onesource::LIRmed)
      .def_readonly("massmed", &onesource::massmed)
      .def_readonly("SFRmed", &onesource::SFRmed)
      .def_readonly("sSFRmed", &onesource::sSFRmed)
      .def_readonly("col1med", &onesource::col1med)
      .def_readonly("col2med", &onesource::col2med)
      .def_readonly("Mrefmed", &onesource::Mrefmed)
      .def_readonly("limits_zmax", &onesource::limits_zmax)
      .def_readonly("limits_Mfaint", &onesource::limits_Mfaint)
      .def_readonly("results_emission_lines",
                    &onesource::results_emission_lines)
      .def_readonly("fluxEL_SED", &onesource::fluxEL_SED)
      .def_readonly("absfilt", &onesource::absfilt)
      .def_readonly("kap", &onesource::kap)
      .def_readonly("mabs", &onesource::mabs)
      .def_readonly("emabs", &onesource::emabs)
      .def_readonly("str_inp", &onesource::str_inp)
      .def_readonly("results", &onesource::results)
      .def_readonly("zgmin", &onesource::zgmin)
      .def_readonly("zgmode", &onesource::zgmode)
      .def_readonly("zgmed", &onesource::zgmed)
      .def_readonly("zqmin", &onesource::zqmin)
      .def_readonly("zqmode", &onesource::zqmode)
      .def_readonly("zqmed", &onesource::zqmed)
      .def_readonly("zsecMod", &onesource::zsecMod)
      .def_readonly("zsecExtlaw", &onesource::zsecExtlaw)
      .def_readonly("zsec", &onesource::zsec)
      .def_readonly("zsecEbv", &onesource::zsecEbv)
      .def_readonly("zsecProb", &onesource::zsecProb)
      .def_readonly("zsecChi2", &onesource::zsecChi2)
      .def_readonly("zsecScale", &onesource::zsecScale)
      .def_readonly("zsecAge", &onesource::zsecAge);

  /******** CLASS PDF *********/
  py::class_<PDF>(mod, "PDF")
      .def(py::init<double, double, size_t>(), py::arg("min"), py::arg("max"),
           py::arg("size"))
      .def("normalization", &PDF::normalization)
      .def("chi2toPDF", &PDF::chi2toPDF)
      .def("chi2mini", &PDF::chi2mini)
      .def("uncMin", &PDF::uncMin)
      .def("index", &PDF::index)
      .def("get_max", &PDF::get_max)
      .def("get_maxid", &PDF::get_maxid)
      .def("secondMax", &PDF::secondMax)
      .def("size", &PDF::size)
      .def("cumulant", &PDF::cumulant)
      .def("levelCumu2x", &PDF::levelCumu2x)
      .def("credible_interval", &PDF::credible_interval)
      .def("confidence_interval", &PDF::confidence_interval)
      .def("improve_extremum", &PDF::improve_extremum)
      .def_readwrite("vPDF", &PDF::vPDF)
      .def_readwrite("xaxis", &PDF::xaxis)
      .def_readwrite("chi2", &PDF::chi2)
      .def_readwrite("secondX", &PDF::secondX)
      .def_readwrite("secondP", &PDF::secondP)
      .def_readwrite("ind", &PDF::ind)
      .def_readwrite("secondInd", &PDF::secondInd);

  // Module-level utility function
  mod.def("quadratic_extremum", &quadratic_extremum,
          R"pbdoc(
Find the extremum of a discrete function using quadratic interpolation.

Parameters
----------
x : array-like
    The independent variable values.
y : array-like
    The dependent variable values (function values).

Returns
-------
float
    The x-coordinate of the extremum.
)pbdoc");

}  // PYBIND11_MODULE
