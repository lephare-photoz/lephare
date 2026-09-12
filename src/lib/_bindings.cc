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
#include "emission_lines.h"
#include "ext.h"
#include "flt.h"
#include "globals.h"
#include "keyword.h"
#include "mag.h"
#include "oneElLambda.h"
#include "opa.h"
#include "photoz_lib.h"

static_assert(PYBIND11_VERSION_MAJOR >= 2 && PYBIND11_VERSION_MINOR >= 11,
              "pybind11 headers are too old");

template <typename x, typename modT>
void applySEDLibTemplate(modT& m, std::string name) {
  py::class_<SEDLib<x>>(m, name.c_str())
      .def(py::init<string, string>(), py::arg("config"), py::arg("typ"))
      .def(py::init<keymap&, string, string>(), py::arg("key_analysed"),
           py::arg("config"), py::arg("typ"))
      .def("print_info", &SEDLib<x>::print_info)
      .def("read_model_list", &SEDLib<x>::read_model_list)
      .def("readSED", &SEDLib<x>::readSED)
      .def("write_SED_lib", &SEDLib<x>::write_SED_lib)
      .def("close_output_files", &SEDLib<x>::close_output_files);
}

PYBIND11_MODULE(_lephare, mod) {
  /*object_type enum for python*/
  py::enum_<object_type>(mod, "object_type")
      .value("GAL", object_type::GAL)
      .value("QSO", object_type::QSO)
      .value("STAR", object_type::STAR);

  /******** CLASS ONEELLAMBDA *********/
  py::class_<oneElLambda>(mod, "oneElLambda")
      .def(py::init<double, double>(), py::arg("lambin"), py::arg("valin"),
           "standard constructor")
      .def(py::init<oneElLambda>(), py::arg("elIn"), "copy constructor")
      .def_readwrite("lamb", &oneElLambda::lamb)
      .def_readwrite("val", &oneElLambda::val);
  mod.def("make_regular_grid", &make_regular_grid);
  mod.def("make_union_grid", &make_union_grid);
  mod.def("common_interpolate_combined", &common_interpolate_combined);
  mod.def("restricted_resampling", &restricted_resampling);
  mod.def("concatenate_and_sort", &concatenate_and_sort,
          "concatenate and sort two vector of oneElLambda objects. Sorting is "
          "in increasing lambda.");

  /******** CLASS COSMOLOGY*********/
  py::class_<cosmo>(mod, "cosmo")
      .def(py::init<double, double, double>(), "Standard constructor",
           py::arg("h0") = 70, py::arg("om0") = 0.3, py::arg("l0") = 0.7)
      .def("distMod", py::vectorize(&cosmo::distMod),
           "Compute distance modulus given redshift z.", py::arg("z"))
      .def("distMet", py::vectorize(&cosmo::distMet),
           "Compute metric distance given redshift z.", py::arg("z"))
      .def("time", py::vectorize(&cosmo::time), "Compute time at redshift z.",
           py::arg("z"))
      .def("flux_rescaling", &cosmo::flux_rescaling, "Compute flux rescaling.")
      .def(py::self == py::self)
      .def(py::self != py::self);

  mod.def("zgrid", &zgrid, "Generate a redshift grid.");
  mod.def("indexz", &indexz, "Get the index for a redshift value.");

  /******** CLASS OPA *********/
  py::class_<opa>(mod, "opa")
      .def(py::init<double, string>(), py::arg("red"), py::arg("opaFile"),
           "standard constructor")
      .def_readwrite("lamb_opa", &opa::lamb_opa)
      .def_readwrite("opared", &opa::red)
      .def("read", &opa::read);
  //   .def("lmin", &opa::lmin, "return smallest wavelength stored")
  //   .def("lmax", &opa::lmax, "return largest wavelength stored")

  /******** CLASS EXT *********/
  py::class_<ext>(mod, "ext")
      .def(py::init<string, int>(), py::arg("name"), py::arg("numext"),
           "standard constructor")
      .def_readwrite("lamb_ext", &ext::lamb_ext)
      .def_readwrite("name", &ext::name)
      .def_readwrite("numext", &ext::numext)
      .def_readonly("lmin", &ext::lmin, "return smallest wavelength stored")
      .def_readonly("lmax", &ext::lmax, "return largest wavelength stored")
      .def("read", &ext::read, py::arg("extFile"), "read an extinction file")
      .def("add_element", &ext::add_element)
      .def("set_vector", &ext::set_vector);
  mod.def("compute_filter_extinction", &compute_filter_extinction,
          "Compute extinction in a filter band.");
  mod.def("cardelli_ext", &cardelli_ext,
          "Compute galactic extinction in the filter based on Cardelli et "
          "al., 1989, ApJ 345",
          py::arg("oneFlt"));
  mod.def("cardelli_law", &cardelli_law,
          "compute albd/av at a given lambda (A) for the Cardelli law",
          py::arg("lb"));

  /******** CLASS KEYWORD *********/
  py::class_<keyword>(mod, "keyword")
      .def_readwrite("name", &keyword::name)
      .def_readwrite("value", &keyword::value)
      .def(py::init())
      .def(py::init<string, string>(), py::arg("n"), py::arg("v"))
      .def("expand_path", &keyword::expand_path)
      .def("split_string", &keyword::split_string)
      .def("split_int", &keyword::split_int)
      .def("split_long", &keyword::split_long)
      .def("split_double", &keyword::split_double)
      .def("split_bool", &keyword::split_bool)
      .def("__repr__", [](const keyword& a) {
        return "(" + a.name + ", " + a.value + ")";
      });

  mod.def("read_command", [](std::vector<std::string> args) {
    std::vector<char*> cstrs;
    cstrs.reserve(args.size());
    for (auto& s : args) cstrs.push_back(const_cast<char*>(s.c_str()));
    return read_command(cstrs.size(), cstrs.data());
  });
  mod.def("read_config", &read_config);

  /******** CLASS FLT *********/
  py::class_<flt>(mod, "flt", py::dynamic_attr())
      .def(py::init<int, string, int, int>(), py::arg("id"), py::arg("name"),
           py::arg("trans"), py::arg("calib"))
      .def(py::init<double, double, int>(), py::arg("lmin"), py::arg("lmax"),
           py::arg("nstep"),
           "Top hat filter from lmin to lmax with nstep points")
      .def("read", static_cast<void (flt::*)(const string&)>(&flt::read),
           "Read filter info from file")
      .def("read", static_cast<void (flt::*)(ifstream&)>(&flt::read),
           "Read filter info from stream")
      .def("lambdaMean", &flt::lambdaMean, "Mean wavelength of the filter.")
      .def("clean", &flt::clean, "Clear the filter transmission curve.")
      .def("lambdaEff", &flt::lambdaEff,
           "Effective wavelength of the filter for a flat-fnu source.")
      .def("lambdaEff2", &flt::lambdaEff2,
           "Effective wavelength of the filter for a flat-flambda source.")
      .def("vega", &flt::vega, "Vega magnitude of the filter.")
      .def("magsun", &flt::magsun, "Absolute magnitude of the Sun in this filter.")
      .def("abcorr", &flt::abcorr, "AB-to-Vega magnitude offset for this filter.")
      .def("width", &flt::width, "Effective width of the filter.")
      .def("lmin", &flt::lmin, "Lower bound of the filter transmission curve.")
      .def("lmax", &flt::lmax, "Upper bound of the filter transmission curve.")
      .def_readonly("name", &flt::name)
      .def_readonly("lmean", &flt::lmean)
      .def_readonly("fcorr", &flt::fcorr)
      .def_readonly("dwidth", &flt::dwidth)
      .def_readwrite("lamb_trans", &flt::lamb_trans)
      .def("data", [](const flt& f) {
        int N = f.lamb_trans.size();
        // Create a 2D array with shape (2, N) (transposed)
        py::array_t<double> result({2, N});
        py::buffer_info buf = result.request();
        double* ptr = static_cast<double*>(buf.ptr);
        for (size_t i = 0; i < N; i++) {
          ptr[i] = f.lamb_trans[i].lamb;     // First row
          ptr[N + i] = f.lamb_trans[i].val;  // Second row
        }
        return result;
      });
  mod.def("read_filters_from_file", &read_filters_from_file);
  mod.def("write_output_filter", &write_output_filter);
  mod.def("read_doc_filters", &read_doc_filters);

  /******** CLASS SED *********/
  py::class_<SED>(mod, "SED")
      .def(py::init<const string, int, string>(), py::arg("name") = "",
           py::arg("nummod") = 0, py::arg("type") = "G")
      .def(py::init<const string, double, double, int, string, int>(),
           py::arg("name"), py::arg("tau"), py::arg("age"), py::arg("nummod"),
           py::arg("type"), py::arg("idAge"))
      .def(py::init<const SED>())
      .def_readwrite("lamb_flux", &SED::lamb_flux)
      .def_readonly("fac_line", &SED::fac_line)
      .def_readonly("extlawId", &SED::extlawId)
      .def_readonly("ebv", &SED::ebv)
      .def_readonly("luv", &SED::luv)
      .def_readonly("lopt", &SED::lopt)
      .def_readonly("lnir", &SED::lnir)
      .def_readonly("ltir", &SED::ltir)
      .def_readonly("age", &SED::age)
      .def_readonly("qi", &SED::qi)
      .def_readonly("name", &SED::name)
      .def_readonly("nummod", &SED::nummod)
      .def_readonly("mag", &SED::mag)
      .def_readwrite("red", &SED::red)
      .def_readwrite("index_z0", &SED::index_z0)
      .def_readwrite("milky_way_extinction", &SED::milky_way_extinction)
      .def_readwrite("band_pass_correction", &SED::band_pass_correction)
      .def("string_to_object", &SED::string_to_object,
           "Convert a G/Q/S type letter to an object_type.")
      .def("redshift", &SED::redshift,
           "Redshift the SED spectrum using the stored 'red' value.")
      .def("is_gal", &SED::is_gal, "True if this SED is a galaxy template.")
      .def("is_star", &SED::is_star, "True if this SED is a star template.")
      .def("is_qso", &SED::is_qso, "True if this SED is a QSO template.")
      .def("read", &SED::read, "Read the SED flux from an ASCII file.",
           py::arg("sedFile"))
      .def("size", &SED::size, "Number of points in the SED spectrum.")
      .def("sumSpectra", &SED::sumSpectra,
           "Add another SED's flux to this one, with a scaling factor.",
           py::arg("addSED"), py::arg("rescal"))
      .def("integrateSED", &SED::integrateSED,
           "Integrate the SED within a filter bandpass.", py::arg("filter"))
      .def("apply_extinction", &SED::apply_extinction, py::arg("ebv"),
           py::arg("oneext"), py::arg("update_ebv") = true)
      .def("apply_extinction_to_lines", &SED::apply_extinction_to_lines,
           "Apply dust extinction to the emission-line fluxes.")
      .def("applyOpa", &SED::applyOpa,
           "Apply intergalactic-medium opacity along the line of sight.")
      .def("integrate", &SED::integrate)
      .def("generateCalib", &SED::generateCalib)
      .def("rescale", &SED::rescale, "Rescale the SED flux by a factor.")
      .def("compute_magnitudes", &SED::compute_magnitudes,
           "Compute synthetic magnitudes in a set of filters.")
      .def("compute_fluxes", &SED::compute_fluxes,
           "Compute synthetic fluxes in a set of filters.")
      .def("generate_spectra", &SED::generate_spectra,
           "Generate the redshifted, normalized spectrum.", py::arg("zin") = 0.0,
           py::arg("dmin") = 1.0)
      .def("emplace_back", &SED::emplace_back)
      .def("set_vector", &SED::set_vector)
      .def("redshift", &SED::redshift)
      //  .def("applyExt", &SED::applyExt)
      .def("compute_milky_way_extinction", &SED::compute_milky_way_extinction,
           "Compute the Milky Way dust extinction curve for this SED.")
      //  .def("applyExtLines", &SED::applyExtLines)
      .def("applyOpa", &SED::applyOpa)
      .def("get_data_vector", &SED::get_data_vector)
      .def("readSEDBin",
           static_cast<void (SED::*)(const string&)>(&SED::readSEDBin),
           "Read the SED from a binary library file.", py::arg("fname"))
      .def("writeSED",
           static_cast<void (SED::*)(const string&, const string&,
                                     const string&)>(&SED::writeSED),
           "Write the SED to binary/physical-parameters/doc files.",
           py::arg("binFile"), py::arg("physFile"), py::arg("docFile"))
      .def("data", [](const SED& f) {
        int N = f.lamb_flux.size();
        // Create a 2D array with shape (2, N) (transposed)
        py::array_t<double> result({2, N});
        py::buffer_info buf = result.request();
        double* ptr = static_cast<double*>(buf.ptr);
        for (size_t i = 0; i < N; i++) {
          ptr[i] = f.lamb_flux[i].lamb;     // First row
          ptr[N + i] = f.lamb_flux[i].val;  // Second row
        }
        return result;
      });
  mod.attr("_emission_lines") = emission_lines;
  mod.attr("_empirical_ratio") = empirical_ratio;
  mod.attr("_empirical_ratio_ori") = empirical_ratio_ori;
  mod.attr("_ga_total") = ga_total;
  mod.attr("_ga_lamb") = ga_lamb;
  mod.attr("_ga_H_val") = ga_H_val;
  mod.attr("_ga_HeI_val") = ga_HeI_val;
  mod.attr("_ga_2q_val") = ga_2q_val;
  // can add a doc: , mod.attr("_ga_2q_val").doc()="internal use only"
  // other option is getter to ensure that it is readonly:
  // mod.def("ga_HeI_val", [] {return ga_HeI_val;});

  // py::class_<SEDlight>(mod, "SEDlight")
  //     .def(py::init<>());  // Constructeur par défaut

  py::class_<StarSED, SED>(mod, "StarSED")
      .def(py::init<const SED&>())
      .def(py::init<const StarSED&>())
      .def(py::init<const string, int>(), py::arg("name"),
           py::arg("nummod") = 0);

  py::class_<QSOSED, SED>(mod, "QSOSED")
      .def(py::init<const SED&>())
      .def(py::init<const QSOSED&>())
      .def(py::init<const string, int>(), py::arg("name"),
           py::arg("nummod") = 0);

  py::class_<GalSED, SED>(mod, "GalSED")
      .def(py::init<const SED&>())
      .def(py::init<const GalSED&>())
      .def(py::init<const string, int>(), py::arg("name"),
           py::arg("nummod") = 0)
      .def(py::init<const string, double, double, string, int, int>(),
           py::arg("name"), py::arg("tau"), py::arg("age"), py::arg("format"),
           py::arg("nummod"), py::arg("idAge"))
      .def_readonly("tau", &GalSED::tau)
      .def_readonly("d4000", &GalSED::d4000)
      .def_readonly("zmet", &GalSED::zmet)
      .def("compute_luminosities", &GalSED::compute_luminosities,
           "Compute the UV/optical/NIR/IR monochromatic luminosities.")
      .def("add_neb_cont", &GalSED::add_neb_cont,
           "Add the nebular continuum emission from the ionizing photon flux.",
           py::arg("qi"))
      .def("generateEmEmpUV", &GalSED::generateEmEmpUV,
           "Empirical emission-line recipe based on the UV magnitude.",
           py::arg("MNUV_int"), py::arg("NUVR"))
      .def("generateEmEmpSFR", &GalSED::generateEmEmpSFR,
           "Empirical emission-line recipe based on the SFR.",
           py::arg("MNUV_int"), py::arg("NUVR"))
      .def("generateEmPhys", &GalSED::generateEmPhys,
           "Physically-motivated emission-line recipe (photoionization).",
           py::arg("zmet"), py::arg("qi"))
      .def("generateEmSpectra", &GalSED::generateEmSpectra,
           "Resample the emission lines onto a spectrum.", py::arg("nstep"))
      .def("sumEmLines", &GalSED::sumEmLines,
           "Add the emission-line flux to the continuum.")
      .def("kcorrec", &GalSED::kcorrec, "Compute the k-correction.",
           py::arg("magz0"))
      .def("rescaleEmLines", &GalSED::rescaleEmLines,
           "Rescale all emission-line fluxes by the fracEm factor.")
      .def("zdepEmLines", &GalSED::zdepEmLines,
           "Apply a redshift-dependent correction to the [OIII] doublet.",
           py::arg("flag"))
      .def("calc_ph", &GalSED::calc_ph,
           "Compute the number of ionizing photons shortward of the H/He edges.");

  /******** CLASS SEDLib *********/
  applySEDLibTemplate<StarSED>(mod, "StarSEDLib");
  applySEDLibTemplate<QSOSED>(mod, "QSOSEDLib");
  applySEDLibTemplate<GalSED>(mod, "GalSEDLib");
  mod.def("_read_ages_from_file", &read_ages_from_file);
  mod.def("_closeAge", &closeAge);
  mod.def("readBC03", &readBC03);
  mod.def("readPEGASE", &readPEGASE);

  /******** CLASS MAG *********/
#define MAGDEFS(c, n)                                     \
  (py::class_<c>(mod, n)                                  \
       .def(py::init<keymap&>(), py::arg("key_analysed")) \
       .def(py::init<>())                                 \
       .def("open_files", &c::open_files,                 \
            "Open the input/output streams needed for the library build.") \
       .def("close_files", &c::close_files, "Close all opened files.") \
       .def("print_info", &c::print_info,                 \
            "Print a summary of the run configuration.")  \
       .def("read_ext", &c::read_ext,                     \
            "Read the extinction laws into extAll.")       \
       .def("read_B12", &c::read_B12,                     \
            "Read the Bethermin et al. (2012) dust emission templates.") \
       .def("set_zgrid", &c::set_zgrid,                   \
            "Set the redshift grid (dz, zmin, zmax).")     \
       .def("read_SED", &c::read_SED,                     \
            "Read the SED files and apply extinction corrections.") \
       .def("write_doc", &c::write_doc,                   \
            "Write the library documentation file.")       \
       .def("make_maglib", &c::make_maglib,                \
            "Build the synthetic magnitude library for one SED.") \
       .def("write_mag", &c::write_mag,                    \
            "Write the computed magnitudes to the output library.") \
       .def_readonly("extAll", &c::extAll)                \
       .def_readonly("opaAll", &c::opaAll)                \
       .def_readonly("allFlt", &c::allFlt))
  MAGDEFS(StarMag, "StarMag");
  MAGDEFS(QSOMag, "QSOMag");
  MAGDEFS(GalMag, "GalMag");

  //   ;

  //   ;

  /******** FUNCTIONS IN GLOBALS.H *********/
  mod.attr("HIGH_CHI2") = HIGH_CHI2;
  mod.attr("INVALID_VAL") = INVALID_VAL;
  mod.def("get_lephare_env", &get_lephare_env);
  mod.def("check_first_char", &check_first_char);
  mod.def("blackbody", &blackbody);
  mod.def("CHECK_CONTEXT_BIT", &CHECK_CONTEXT_BIT);
  mod.def("POW10D", &POW10D);
  mod.def("LOG10D", &LOG10D);
  mod.def("POW10D_SLOW", &POW10D_SLOW);
  mod.def("POW10D_FAST", &POW10D_FAST);
  mod.def("POW10D_FASTV", &POW10D_FASTV);
  mod.def("POW10D_SLOWV", &POW10D_SLOWV);
  mod.def("LOG10D_SLOW", &LOG10D_SLOW);
  mod.def("LOG10D_FAST", &LOG10D_FAST);
  mod.def("LOG10D_SLOWV", &LOG10D_SLOWV);
  mod.def("LOG10D_FASTV", &LOG10D_FASTV);
  mod.def("mag2flux", &mag2flux);
  mod.def("flux2mag", &flux2mag);
  mod.def("indexes_in_vec", &indexes_in_vec);
  mod.def("fast_interpolate", &fast_interpolate);
  // return a copy to python, only for unit tests
  mod.def("_get_opa_vector", &get_opa_vector);

  /******** FUNCTIONS IN PHOTOZ_LIB.H *********/
  py::class_<PhotoZ>(mod, "PhotoZ")
      .def_readonly("flux", &PhotoZ::flux)
      .def_readonly("fluxIR", &PhotoZ::fluxIR)
      .def_readonly("imagm", &PhotoZ::imagm)
      .def_readonly("fullLib", &PhotoZ::fullLib)
      //.def_readonly("lightLib", &PhotoZ::lightLib)
      .def_readonly("zLib", &PhotoZ::zLib)
      .def_readonly("fullLibIR", &PhotoZ::fullLibIR)
      .def_readonly("allFilters", &PhotoZ::allFilters)
      .def_readonly("gridz", &PhotoZ::gridz)
      .def_readonly("outkeywords", &PhotoZ::outkeywords)
      .def_readonly("outpara", &PhotoZ::outpara)
      .def_readonly("pdftype", &PhotoZ::pdftype)
      .def_readwrite("outputHeader", &PhotoZ::outputHeader)
      .def_readwrite("reddening", &PhotoZ::reddening)
      .def_readwrite("mw_classic_extinction_values",
                     &PhotoZ::mw_classic_extinction_values)
      .def(py::init<keymap&>())
      .def("read_autoadapt_sources", &PhotoZ::read_autoadapt_sources)
      .def("belong_autoadapt", &PhotoZ::belong_autoadapt)
      .def("read_photoz_sources", &PhotoZ::read_photoz_sources)
      .def("read_mw_ebv", &PhotoZ::read_mw_ebv)
      .def("prep_data", static_cast<void (PhotoZ::*)(vector<onesource*>)>(
                            &PhotoZ::prep_data))
      .def("prep_data",
           static_cast<void (PhotoZ::*)(onesource*)>(&PhotoZ::prep_data))
      .def("run_autoadapt", &PhotoZ::run_autoadapt)
      .def("run_photoz", &PhotoZ::run_photoz)
      .def("fit", &PhotoZ::fit)
      .def("fit_uncertainties", &PhotoZ::fit_uncertainties)
      .def("physical_parameters", &PhotoZ::physical_parameters)
      .def("best_template", &PhotoZ::best_template)
      .def("write_spectrum", &PhotoZ::write_spectrum)
      .def("write_outputs", &PhotoZ::write_outputs)
      .def("validLib", &PhotoZ::validLib)
      .def("compute_offsets", &PhotoZ::compute_offsets);
  // mod.def("read_lib", [](const string& libName, int ind, vector<int>
  // emMod, int babs) { 			vector<SED*> libFull; int
  // nummodpre[3]; 			string filtname; vector<double>
  // gridz; 			nummodpre[0] = 0; nummodpre[1] = 0;
  // nummodpre[2] = 0; read_lib(libFull, ind, nummodpre, libName, filtname,
  // gridz, emMod, babs); 			std::array<int, 3> nummod_tup =
  // {nummodpre[0], nummodpre[1], nummodpre[2]}; return
  // std::make_tuple(libFull, ind, nummod_tup, filtname, gridz);
  // 		      }
  // 	  );
  // mod.def("read_doc_filters", [](const string filtFile) {
  // 				bool Fexiste;
  // 				vector<flt> allFilters =
  // read_doc_filters(filtFile, Fexiste); 				return
  // std::make_tuple(allFilters, Fexiste);
  // 			      }
  // );
  mod.def("readOutKeywords", &readOutKeywords);
  mod.def("bestFilter", &bestFilter);
  mod.def("maxkcolor", &maxkcolor);

  mod.attr("maptype") = maptype;
  py::class_<onesource>(mod, "onesource", py::dynamic_attr())
      .def(py::init<>())
      .def(py::init<const int, vector<double>>(), py::arg("pos"), py::arg("gridz"))
      .def("setPriors", &onesource::setPriors,
           "Set the absolute magnitude prior range for the fit.",
           py::arg("magabsB"), py::arg("magabsF"))
      .def_readonly("priorLib", &onesource::priorLib)
      //    .def("readsource", &onesource::readsource)
      .def("readsource",
           static_cast<void (onesource::*)(
               const string&, const vector<double>, const vector<double>,
               const long, const double, const string)>(&onesource::readsource),
           "Set the observed fluxes, errors and metadata for this source.",
           py::arg("identifier"), py::arg("vals"), py::arg("err_vals"),
           py::arg("context"), py::arg("z_spec"), py::arg("additional_input"))
      .def("set_verbosity", &onesource::set_verbosity, py::arg("verbose"))
      .def("get_verbosity", &onesource::get_verbosity)
      .def("computePredAbsMag", &onesource::computePredMag)
      .def("computePredAbsMag", &onesource::computePredAbsMag)
      .def("fltUsed", &onesource::fltUsed,
           "Flag which bands are used in the fit and which are upper limits.",
           py::arg("gbcont"), py::arg("contforb"))
      .def("convertFlux", &onesource::convertFlux,
           "Convert the input magnitudes/fluxes to the internal flux convention.",
           py::arg("catmag"), py::arg("allFilters"))
      .def("convertMag", &onesource::convertMag,
           "Convert fluxes to AB magnitudes and their errors.")
      .def("rescale_flux_errors", &onesource::rescale_flux_errors,
           "Rescale the flux errors (e.g. add a systematic floor).",
           py::arg("min_err"), py::arg("fac_err"))
      .def("keepOri", &onesource::keepOri,
           "Save a copy of the fluxes/magnitudes before correction.")
      .def("adapt_mag", &onesource::adapt_mag,
           "Apply a per-band zero-point offset.", py::arg("a0"))
      .def("fit", &onesource::fit,
           "Fit the source against a SED library by chi2 minimisation.",
           py::arg("lightLib"), py::arg("flux"), py::arg("valid"),
           py::arg("funz0"), py::arg("bp"), py::arg("restrict"))
      .def("nzprior", &onesource::nzprior,
           "Compute the N(z) prior weight applied to the chi2.",
           py::arg("luv"), py::arg("lnir"), py::arg("reds"), py::arg("bp"))
      .def("mode", &onesource::mode,
           "Compute the mode of the marginalized redshift PDFs.")
      .def("rm_discrepant", &onesource::rm_discrepant,
           "Iteratively remove discrepant bands and re-fit.",
           py::arg("lightLib"), py::arg("flux"), py::arg("valid"),
           py::arg("funz0"), py::arg("bp"), py::arg("thresholdChi2"),
           py::arg("restrict"))
      .def("generatePDF", &onesource::generatePDF,
           "Build the marginalized redshift/physical-parameter PDFs.",
           py::arg("lightLib"), py::arg("va"), py::arg("colAnalysis"),
           py::arg("zfix"))
      .def("interp", &onesource::interp,
           "Update the fit solution (fixed redshift or parabolic refinement).",
           py::arg("zfix"), py::arg("zintp"), py::arg("lcdm"))
      .def("uncertaintiesMin", &onesource::uncertaintiesMin,
           "Compute chi2-based confidence intervals around the minimum.")
      .def("uncertaintiesBay", &onesource::uncertaintiesBay,
           "Compute Bayesian confidence intervals from the marginalized PDFs.")
      .def("secondpeak", &onesource::secondpeak,
           "Detect a secondary peak in the marginalized redshift PDF.",
           py::arg("lightLib"), py::arg("dz_win"), py::arg("min_thres"))
      .def("absmag", &onesource::absmag,
           "Compute the absolute magnitude(s) from the best-fit template(s).",
           py::arg("bestFlt"), py::arg("maxkcolor"), py::arg("lcdm"),
           py::arg("gridz"))
      .def("limits", &onesource::limits,
           "Compute the faint-end absolute-magnitude limit for this source.",
           py::arg("fulllib"), py::arg("limits_zbin"), py::arg("limits_ref"),
           py::arg("limits_sel"), py::arg("limits_cut"))
      .def("computeEmFlux", &onesource::computeEmFlux,
           "Compute the emission-line flux/EW of the best-fit template.")
      .def("generatePDF_IR", &onesource::generatePDF_IR,
           "Build the marginalized IR luminosity PDF from the FIR fit.")
      .def("write_out", &onesource::write_out,
           "Write this source's output line in the LePHARE ASCII format.",
           py::arg("stout"), py::arg("outkeywords"))
      .def("redden_flux", &onesource::redden_flux,
           "Apply a per-template reddening correction to a predicted flux array.",
           py::arg("flux"), py::arg("reddening"))
      .def("writeSpec", &onesource::writeSpec,
           "Write the best-fit spectra to a per-source output file.")
      .def("writeFullChi", &onesource::writeFullChi,
           "Write the chi2 of every template to a per-source .chi file.",
           py::arg("lightLib"))
      .def("best_spec_vec", &onesource::best_spec_vec,
           "Flux-averaged best-fit spectrum between two wavelengths.",
           py::arg("sol"), py::arg("fulllib"), py::arg("lcdm"),
           py::arg("minl"), py::arg("maxl"))
      .def_readwrite("spec", &onesource::spec)
      .def_readwrite("consiz", &onesource::consiz)
      .def_readwrite("mw_ebv", &onesource::mw_ebv)
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
      .def_readonly("abIR", &onesource::abIR)
      .def_readonly("ab_ori", &onesource::ab_ori)
      .def_readonly("sab", &onesource::sab)
      .def_readonly("sabIR", &onesource::sabIR)
      .def_readonly("mab", &onesource::mab)
      .def_readonly("msab", &onesource::msab)
      .def_readonly("magm", &onesource::magm)
      .def_readonly("magPred", &onesource::magPred)
      .def_readonly("absmagPred", &onesource::absmagPred)
      .def_readonly("zmin", &onesource::zmin)
      .def_readonly("zminIR", &onesource::zminIR)
      .def_readonly("chimin", &onesource::chimin)
      .def_readonly("chiminIR", &onesource::chiminIR)
      .def_readonly("dmmin", &onesource::dmmin)
      .def_readonly("dmminIR", &onesource::dmminIR)
      .def_readonly("indmin", &onesource::indmin)
      .def_readonly("indminSec", &onesource::indminSec)
      .def_readonly("indminIR", &onesource::indminIR)
      .def_readonly("imasmin", &onesource::imasmin)
      .def_readonly("imasminIR", &onesource::imasminIR)
      .def_readonly("agemed", &onesource::agemed)
      .def_readonly("ebvmed", &onesource::ebvmed)
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
      // output parameters:
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
      .def_readonly("zsecProb", &onesource::zsecMod)
      .def_readonly("zsecChi2", &onesource::zsecMod)
      .def_readonly("zsecScale", &onesource::zsecMod)
      .def_readonly("zsecAge", &onesource::zsecAge);

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
  mod.def("quadratic_extremum", &quadratic_extremum);
}  // PYBIND11_MODULE
