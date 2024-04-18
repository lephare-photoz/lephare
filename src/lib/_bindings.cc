#include <pybind11/numpy.h>
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
void applySEDLibTemplate(modT &m, std::string name) {
  py::class_<SEDLib<x>>(m, name.c_str())
      .def(py::init<string, string>(), py::arg("config"), py::arg("typ"))
      .def(py::init<keymap &, string, string>(), py::arg("key_analysed"),
           py::arg("config"), py::arg("typ"))
      .def("print_info", &SEDLib<x>::print_info)
      .def("read_model_list", &SEDLib<x>::read_model_list)
      .def("write_SED_lib", &SEDLib<x>::write_SED_lib)
      .def("print_time_tofile", &SEDLib<x>::print_time_tofile)
      .def("close_output_files", &SEDLib<x>::close_output_files);
}

PYBIND11_MODULE(_lephare, mod) {
  /******** CLASS ONEELLAMBDA *********/
  py::class_<oneElLambda>(mod, "oneElLambda")
      .def(py::init<double, double, int>(), py::arg("lambin"), py::arg("valin"),
           py::arg("oriin"), "standard constructor")
      .def(py::init<oneElLambda>(), py::arg("elIn"), "copy constructor")
      .def_readwrite("lamb", &oneElLambda::lamb)
      .def_readwrite("val", &oneElLambda::val)
      .def_readwrite("ori", &oneElLambda::ori)
      .def("interp", &oneElLambda::interp, py::arg("previousEl"),
           py::arg("nextEl"));

  /******** CLASS COSMOLOGY*********/
  py::class_<cosmo>(mod, "cosmo")
      .def(py::init<double, double, double>(), py::arg("h0") = 70,
           py::arg("om0") = 0.3, py::arg("l0") = 0.7, "standard constructor")
      .def("distMod", &cosmo::distMod, py::arg("z"))
      .def("distMet", &cosmo::distMet, py::arg("z"))
      .def("time", &cosmo::time, py::arg("z"))
      .def("distMod", py::vectorize(&cosmo::distMod))
      .def("distMet", py::vectorize(&cosmo::distMet))
      .def("time", py::vectorize(&cosmo::time));
  mod.def("zgrid", &zgrid);
  mod.def("indexz", &indexz);

  /******** CLASS OPA *********/
  py::class_<opa>(mod, "opa")
      .def(py::init<double, string>(), py::arg("red"), py::arg("opaFile"),
           "standard constructor")
      .def_readwrite("lamb_opa", &opa::lamb_opa)
      .def_readwrite("opared", &opa::red)
      .def("read", &opa::read)
      //   .def("lmin", &opa::lmin, "return smallest wavelength stored")
      //   .def("lmax", &opa::lmax, "return largest wavelength stored")
      ;

  /******** CLASS EXT *********/
  py::class_<ext>(mod, "ext")
      .def(py::init<string, int>(), py::arg("name"), py::arg("numext"),
           "standard constructor")
      .def_readwrite("lamb_ext", &ext::lamb_ext)
      .def_readwrite("name", &ext::name)
      .def_readwrite("numext", &ext::numext)
      .def_readonly("lmin", &ext::lmin, "return smallest wavelength stored")
      .def_readonly("lmax", &ext::lmax, "return largest wavelength stored")
      .def("read", &ext::read, py::arg("extFile"), "read an extinction file");

  /******** CLASS KEYWORD *********/
  py::class_<keyword>(mod, "keyword")
      .def_readwrite("name", &keyword::name)
      .def_readwrite("value", &keyword::value)
      .def_readwrite("def", &keyword::def)
      .def(py::init())
      .def(py::init<string, string>(), py::arg("n"), py::arg("v"))
      .def("expand_path", &keyword::expand_path)
      .def("split_string", &keyword::split_string)
      .def("split_int", &keyword::split_int)
      .def("split_long", &keyword::split_long)
      .def("split_double", &keyword::split_double)
      .def("split_bool", &keyword::split_bool)
      .def("__repr__", [](const keyword &a) {
        return "(" + a.name + ", " + a.value + ")";
      });

  mod.def("read_command", [](std::vector<std::string> args) {
    std::vector<char *> cstrs;
    cstrs.reserve(args.size());
    for (auto &s : args) cstrs.push_back(const_cast<char *>(s.c_str()));
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
      .def("read", static_cast<void (flt::*)(const string &)>(&flt::read),
           "Read filter info from file")
      .def("read", static_cast<void (flt::*)(ifstream &)>(&flt::read),
           "Read filter info from stream")
      .def("lambdaMean", &flt::lambdaMean)
      .def("lambdaEff", &flt::lambdaEff)
      .def("magsun", &flt::magsun)
      .def("width", &flt::width)
      .def("peak", &flt::peak)
      .def_readonly("name", &flt::name)
      .def_readonly("lmean", &flt::lmean)
      .def_readonly("dwidth", &flt::dwidth)
      .def_readwrite("lamb_trans", &flt::lamb_trans)
      .def("data", [](const flt &f) {
        int N = f.lamb_trans.size();
        py::array_t<double> result = py::array_t<double>(N * 2);
        py::buffer_info buf = result.request();
        double *ptr = (double *)buf.ptr;
        for (size_t i = 0; i < N; i++) {
          ptr[i] = f.lamb_trans[i].lamb;
          ptr[N + i] = f.lamb_trans[i].val;
        }
        result.resize({N, N});
        return result;
      });
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
      .def_readonly("extlawId", &SED::extlawId)
      .def_readonly("ebv", &SED::ebv)
      .def_readonly("name", &SED::name)
      .def_readonly("nummod", &SED::nummod)
      .def_readonly("mag", &SED::mag)
      .def_readwrite("index_z0", &SED::index_z0)
      .def("is_gal", &SED::is_gal)
      .def("is_star", &SED::is_star)
      .def("is_qso", &SED::is_qso)
      .def("read", &SED::read)
      .def("size", &SED::size)
      .def("integrateSED", &SED::integrateSED)
      .def("resample", &SED::resample)
      .def("generateCalib", &SED::generateCalib)
      .def("rescale", &SED::rescale)
      .def("applyShift", &SED::applyShift)
      .def("compute_magnitudes", &SED::compute_magnitudes)
      .def("compute_fluxes", &SED::compute_fluxes)
      .def("set_vector", &SED::set_vector)
      .def("readSEDBin",
           static_cast<void (SED::*)(const string &)>(&SED::readSEDBin))
      .def("writeSED",
           static_cast<void (SED::*)(const string &, const string &,
                                     const string &)>(&SED::writeSED))
      .def("data", [](const SED &f) {
        int N = f.lamb_flux.size();
        py::array_t<double> result = py::array_t<double>(N * 2);
        py::buffer_info buf = result.request();
        double *ptr = (double *)buf.ptr;
        for (size_t i = 0; i < N; i++) {
          ptr[i] = f.lamb_flux[i].lamb;
          ptr[N + i] = f.lamb_flux[i].val;
        }
        result.resize({N, N});
        return result;
      });

  py::class_<StarSED, SED>(mod, "StarSED")
      .def(py::init<const SED &>())
      .def(py::init<const StarSED &>())
      .def(py::init<const string, int>(), py::arg("name"),
           py::arg("nummod") = 0);

  py::class_<QSOSED, SED>(mod, "QSOSED")
      .def(py::init<const SED &>())
      .def(py::init<const QSOSED &>())
      .def(py::init<const string, int>(), py::arg("name"),
           py::arg("nummod") = 0);

  py::class_<GalSED, SED>(mod, "GalSED")
      .def(py::init<const SED &>())
      .def(py::init<const GalSED &>())
      .def(py::init<const string, int>(), py::arg("name"),
           py::arg("nummod") = 0)
      .def(py::init<const string, double, double, string, int, string, int>(),
           py::arg("name"), py::arg("tau"), py::arg("age"), py::arg("format"),
           py::arg("nummod"), py::arg("type"), py::arg("idAge"))
      .def("SEDproperties", &GalSED::SEDproperties)
      .def("add_neb_cont", &GalSED::add_neb_cont)
      .def("generateEmEmpUV", &GalSED::generateEmEmpUV)
      .def("generateEmEmpSFR", &GalSED::generateEmEmpSFR)
      .def("generateEmPhys", &GalSED::generateEmPhys)
      //    .def("generateEmSpectra", &GalSED::generateEmSpectra)
      .def("sumEmLines", &GalSED::sumEmLines)
      .def("kcorrec", &GalSED::kcorrec)
      .def("rescaleEmLines", &GalSED::rescaleEmLines)
      .def("zdepEmLines", &GalSED::zdepEmLines)
      .def("calc_ph", &GalSED::calc_ph);

  /******** CLASS SEDLib *********/
  applySEDLibTemplate<StarSED>(mod, "StarSEDLib");
  applySEDLibTemplate<QSOSED>(mod, "QSOSEDLib");
  applySEDLibTemplate<GalSED>(mod, "GalSEDLib");

  /******** CLASS MAG *********/
#define MAGDEFS(c, n)                                      \
  (py::class_<c>(mod, n)                                   \
       .def(py::init<keymap &>(), py::arg("key_analysed")) \
       .def(py::init<>())                                  \
       .def("open_files", &c::open_files)                  \
       .def("close_files", &c::close_files)                \
       .def("open_opa_files", &c::open_opa_files)          \
       .def("print_info", &c::print_info)                  \
       .def("read_ext", &c::read_ext)                      \
       .def("read_opa", &c::read_opa)                      \
       .def("read_B12", &c::read_B12)                      \
       .def("read_flt", &c::read_flt)                      \
       .def("def_zgrid", &c::def_zgrid)                    \
       .def("read_SED", &c::read_SED)                      \
       .def("write_doc", &c::write_doc)                    \
       .def("make_maglib", &c::make_maglib)                \
       .def("write_mag", &c::write_mag))
  MAGDEFS(StarMag, "StarMag");
  MAGDEFS(QSOMag, "QSOMag");
  MAGDEFS(GalMag, "GalMag");

  //   ;

  //   ;

  /******** FUNCTIONS IN GLOBALS.H *********/
  mod.def("get_lephare_env", &get_lephare_env);
  mod.def("check_first_char", &check_first_char);
  mod.def("blackbody", &blackbody);
  mod.def("CHECK_CONTEXT_BIT", &CHECK_CONTEXT_BIT);
  mod.def("POW10D", &POW10D);
  mod.def("POW10DSLOW", &POW10DSLOW);
  mod.def("LOG10D_SLOW", &LOG10D_SLOW);
  mod.def("LOG10D_FAST", &LOG10D_FAST);
  mod.def("mag2flux", &mag2flux);
  mod.def("flux2mag", &flux2mag);

  /******** FUNCTIONS IN PHOTOZ_LIB.H *********/
  py::class_<PhotoZ>(mod, "PhotoZ")
      .def_readonly("flux", &PhotoZ::flux)
      .def_readonly("fluxIR", &PhotoZ::fluxIR)
      .def_readonly("imagm", &PhotoZ::imagm)
      .def_readonly("fullLib", &PhotoZ::fullLib)
      .def_readonly("fullLibIR", &PhotoZ::fullLibIR)
      .def_readonly("allFilters", &PhotoZ::allFilters)
      .def_readonly("gridz", &PhotoZ::gridz)
      .def_readonly("outkeywords", &PhotoZ::outkeywords)
      .def_readonly("outpara", &PhotoZ::outpara)
      .def_readonly("pdftype", &PhotoZ::pdftype)
      .def_readwrite("outputHeader", &PhotoZ::outputHeader)
      .def(py::init<keymap &>())
      .def("read_autoadapt_sources", &PhotoZ::read_autoadapt_sources)
      .def("read_photoz_sources", &PhotoZ::read_photoz_sources)
      .def("prep_data", static_cast<void (PhotoZ::*)(vector<onesource *>)>(
                            &PhotoZ::prep_data))
      .def("prep_data",
           static_cast<void (PhotoZ::*)(onesource *)>(&PhotoZ::prep_data))
      .def("run_autoadapt", &PhotoZ::run_autoadapt)
      .def("run_photoz", &PhotoZ::run_photoz)
      .def("write_outputs", &PhotoZ::write_outputs);
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
  mod.def("readOutKeywords", &readOutKeywords);
  mod.def("bestFilter", &bestFilter);
  mod.def("maxkcolor", &maxkcolor);

  mod.attr("maptype") = maptype;
  py::class_<onesource>(mod, "onesource")
      .def(py::init<>())
      .def(py::init<const int, vector<double>>())
      .def("setPriors", &onesource::setPriors)
      .def_readonly("priorLib", &onesource::priorLib)
      //    .def("readsource", &onesource::readsource)
      .def("readsource",
           static_cast<void (onesource::*)(
               const string &, const vector<double>, const vector<double>,
               const long, const double, const string)>(&onesource::readsource))
      .def("fltUsed", &onesource::fltUsed)
      .def("convertFlux", &onesource::convertFlux)
      .def("convertMag", &onesource::convertMag)
      .def("rescale_flux_errors", &onesource::rescale_flux_errors)
      .def("keepOri", &onesource::keepOri)
      .def("adapt_mag", &onesource::adapt_mag)
      .def("fit", &onesource::fit)
      .def("validLib", &onesource::validLib)
      .def("rm_discrepant", &onesource::rm_discrepant)
      .def("generatePDF", &onesource::generatePDF)
      .def("interp", &onesource::interp)
      .def("uncertaintiesMin", &onesource::uncertaintiesMin)
      .def("uncertaintiesBay", &onesource::uncertaintiesBay)
      .def("secondpeak", &onesource::secondpeak)
      .def("considered_red", &onesource::considered_red)
      .def("interp_lib", &onesource::interp_lib)
      .def("absmag", &onesource::absmag)
      .def("limits", &onesource::limits)
      .def("computePredAbsMag", &onesource::computePredMag)
      .def("computePredAbsMag", &onesource::computePredAbsMag)
      .def("computeEmFlux", &onesource::computeEmFlux)
      .def("fltUsedIR", &onesource::fltUsedIR)
      .def("substellar", &onesource::substellar)
      .def("generatePDF_IR", &onesource::generatePDF_IR)
      .def("write_out", &onesource::write_out)
      .def("writeSpec", &onesource::writeSpec)
      .def("writeFullChi", &onesource::writeFullChi)
      //    .def("write_pdz", &onesource::write_pdz)
      .def("best_spec_vec", &onesource::best_spec_vec)
      .def_readwrite("spec", &onesource::spec)
      .def_readwrite("consiz", &onesource::consiz)
      .def_readwrite("closest_red", &onesource::closest_red)
      .def_readonly("pos", &onesource::pos)
      .def_readonly("cont", &onesource::cont)
      .def_readonly("pdfmap", &onesource::pdfmap)
      .def_readonly("nbused", &onesource::nbused)
      .def_readonly("nbul", &onesource::nbul)
      .def_readonly("dm", &onesource::dm)
      .def_readonly("zs", &onesource::zs)
      .def_readonly("ab", &onesource::ab)
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
      .def_readonly("imasminIR", &onesource::imasmin)
      .def_readonly("imasmin", &onesource::imasmin)
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
      .def("int_parab", &PDF::int_parab)
      .def("uncMin", &PDF::uncMin)
      .def("index", &PDF::index)
      .def("get_max", &PDF::get_max)
      .def("get_maxid", &PDF::get_maxid)
      .def("secondMax", &PDF::secondMax)
      .def("size", &PDF::size)
      .def("levelCumu2x", &PDF::levelCumu2x)
      .def("credible_interval", &PDF::credible_interval)
      .def_readwrite("vPDF", &PDF::vPDF)
      .def_readwrite("xaxis", &PDF::xaxis)
      .def_readwrite("chi2", &PDF::chi2)
      .def_readwrite("secondX", &PDF::secondX)
      .def_readwrite("secondP", &PDF::secondP)
      .def_readwrite("ind", &PDF::ind)
      .def_readwrite("secondInd", &PDF::secondInd);
}  // PYBIND11_MODULE
