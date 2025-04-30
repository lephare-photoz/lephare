import time

from lephare._lephare import (  # , read_lib, read_doc_filters, readOutKeywords, GalMag, bestFilter, maxkcolor
    GalMag,
    PhotoZ,
    keyword,
)
from lephare.runner import Runner

__all__ = [
    "Zphota",
]

config_keys = {
    "CAT_IN": "input catalog",
    "INP_TYPE": "F if Fluxes, M if magnitudes",
    "CAT_MAG": "AB or VEGA magnitude system of the input catalog (default AB)",
    "CAT_FMT": "MEME or MMEE, to precise the order of value and error columns",
    "CAT_LINES": "min,max range to subselect entries in the input catalog",
    "ZPHOTLIB": "comma-separated list (max 3) of template magnitude libraries,\
    without the suffix (to be found in $LEPHAREWORK/lib_mag)",
    "PARA_OUT": "file listing the output keywords to be written on file (one keyword per line)",
    "CAT_OUT": "output file name",
    "VERBOSE": "add verbosity",
    "CAT_TYPE": "LONG or SHORT, depending on whether zspec, context, and other columns are present\
    after the magnitudes/fluxes and error columns",
    "ERR_SCALE": "values to add to the errors in each band",
    "ERR_FACTOR": "values to scale the errors in each band",
    "GLB_CONTEXT": "global context to define the bands to be used in the fit",
    "FORB_CONTEXT": "forbidden context to define the bands to be discarded from the fit",
    "MAG_ABS": "min,max absolute magnitudes allowed for GAL type",
    "MAG_ABS_QSO": "min,max absolute magnitudes allowed for QSO/AGN type",
    "MAG_REF": "reference band to compute absolute magnitudes",
    "NZ_PRIOR": "band to be used in the n(z) prior; a second band can be given\
    in case the first one is absent in some objects",
    "ZFIX": "fix z to a known redshift, in case one is interested in inferring the best template",
    "Z_INTERP": "perform quadratic interpolation to refine best redshift solution",
    "Z_METHOD": "which redshift solution, between best and pdf median,\
    to use for physical parameter computation",
    "DZ_WIN": "Window search for second peak",
    "MIN_THRES": "lower threshold for second peak",
    "SPEC_OUT": "output the best spectrum of each object",
    "FIR_LIB": "dedicated far infrared synthetic magnitudes files",
    "FIR_LMIN": "minimum value for lambda_mean/(1+z) (in microns) in FIR analysis",
    "FIR_CONT": "context for the FIR bands to be used in the FIR fit",
    "FIR_SCALE": "context for the FIR scale determination in the FIR fit",
    "FIR_FREESCALE": "consider FIR scaling free",
    "FIR_SUBSTELLAR": "substract the stellar component to the observed IR fluxes",
    "MABS_METHOD": "",
    "MABS_CONTEXT": "",
    "MABS_REF": "",
    "MABS_FILT": "",
    "MABS_ZBIN": "",
    "RF_COLORS": "list of 4 band indexes out of which to compute the probability\
    distribution of the 2 corresponding rest frame colors (default=-1 which means\
    to not compute these color distributions)",
    "M_REF": "band index for which to compute the rest frame absolute magnitude\
    probability distribution",
    "APPLY_SYSSHIFT": "list of values (equal to the number of bands) to add to the predicted\
    magnitudes as zerop points  (substracted to the observed magnitudes).\
    These are also the output of the auto adaptation stage",
    "AUTO_ADAPT": "perform zero point auto adaptation based on the provided true redshifts",
    "ADAPT_BAND": "band number for which to check magnitude range for auto adaptation",
    "ADAPT_LIM": "min,max magnitude range for auto adaptation",
    "ADAPT_ZBIN": "min,max value of the true redshifts used for auto adaptation of the zero points",
    "PDZ_OUT": "output the PDF of each object",
    "PDZ_TYPE": "list of the PDFs to output if PDZ_OUT is set",
    "ADD_EMLINES": "add emission lines during the fit",
    "ADDITIONAL_MAG": "file basename (in $LEPHAREWORK/filt/) for the set of additional filters\
    for which to compute magnitudes. Thus a pass to the LePHARE `filter` stage needs to have\
    occurred on these filters.",
    "LIMITS_ZBIN": "list of redshifts on which onesource::limits computes the magnitude limits",
    "LIMITS_MAPP_REF": "",
    "LIMITS_MAPP_SEL": "",
    "LIMITS_MAPP_CUT": "",
}

nonestring = "NONE"


class Zphota(Runner):
    """The specific arguments to the Zphota class are

    CAT_IN:
           input catalog
    INP_TYPE:
           F if Fluxes, M if magnitudes
    CAT_MAG:
           AB or VEGA magnitude system of the input catalog (default AB)
    CAT_FMT:
           MEME or MMEE, to precise the order of value and error columns
    CAT_LINES:
           min,max range to subselect entries in the input catalog
    ZPHOTLIB:
           comma-separated list (max 3) of template magnitude libraries\
without the suffix (to be found in $LEPHAREWORK/lib_mag)
    PARA_OUT:
           file listing the output keywords to be written on file (one keyword per line)
    CAT_OUT:
           output file name
    VERBOSE:
           add verbosity
    CAT_TYPE:
           LONG or SHORT, depending on whether zspec, context, and other columns are present\
after the magnitudes/fluxes and error columns
    ERR_SCALE:
           values to add to the errors in each band
    ERR_FACTOR:
           values to scale the errors in each band
    GLB_CONTEXT:
           global context to define the bands to be used in the fit
    FORB_CONTEXT:
           forbidden context to define the bands to be discarded from the fit
    MAG_ABS:
           min,max absolute magnitudes allowed for GAL type
    MAG_ABS_QSO:
           min,max absolute magnitudes allowed for QSO/AGN type
    MAG_REF:
           reference band to compute absolute magnitudes
    NZ_PRIOR:
           band to be used in the n(z) prior; a second band can be given\
in case the first one is absent in some objects
    ZFIX:
           fix z to a known redshift, in case one is interested in inferring the best template
    Z_INTERP:
           perform quadratic interpolation to refine best redshift solution
    Z_METHOD:
           which redshift solution, between best and pdf median,\
to use for physical parameter computation
    DZ_WIN:
           Window search for second peak
    MIN_THRES:
           lower threshold for second peak
    PROB_INTZ:

    SPEC_OUT:
           output the best spectrum of each object
    FIR_LIB:
           dedicated far infrared synthetic magnitudes files
    FIR_LMIN:
           minimum value for lambda_mean/(1+z) (in microns) in FIR analysis
    FIR_CONT:
           context for the FIR bands to be used in the FIR fit
    FIR_SCALE:
           context for the FIR scale determination in the FIR fit
    FIR_FREESCALE:
           consider FIR scaling free
    FIR_SUBSTELLAR:
           substract the stellar component to the observed IR fluxes
    MABS_METHOD:

    MABS_CONTEXT:

    MABS_REF:

    MABS_FILT:

    MABS_ZBIN:

    RF_COLORS:
           list of 4 band indexes out of which to compute the probability\
distribution of the 2 corresponding rest frame colors (default=-1 which means\
to not compute these color distributions)
    M_REF:
           band index for which to compute the rest frame absolute magnitude\
probability distribution
    APPLY_SYSSHIFT:
           list of values (equal to the number of bands) to add to the magnitudes\
as zerop points. These are also the output of the auto adaptation stage
    AUTO_ADAPT:
           perform zero point auto adaptation based on the provided true redshifts
    ADAPT_BAND:
           band number for which to check magnitude range for auto adaptation
    ADAPT_LIM:
           min,max magnitude range for auto adaptation
    ADAPT_ZBIN:
           min,max value of the true redshifts used for auto adaptation of the zero points
    PDZ_OUT:
           output the PDF of each object
    PDZ_TYPE:
           list of the PDFs to output if PDZ_OUT is set
    ADD_EMLINES:
           add emission lines during the fit
    ADDITIONAL_MAG:
           file basename (in $LEPHAREWORK/filt/) for the set of additional filters\
for which to compute magnitudes. Thus a pass to the LePHARE `filter` stage needs to have\
occurred on these filters.
    LIMITS_ZBIN:
           list of redshifts on which onesource::limits computes the magnitude limits
    LIMITS_MAPP_REF:

    LIMITS_MAPP_SEL:

    LIMITS_MAPP_CUT:
    """

    def add_authorized_keys(self):
        """Add the specific Zphota arguments to the argument parser"""
        for key in config_keys:
            self.parser.add_argument("--%s" % key, type=str, metavar="", help=config_keys[key])
        self.parser.usage = "Execute LePHARE photo-z estimation"

    def __init__(self, config_file=None, config_keymap=None, **kwargs):
        super().__init__(config_keys, config_file, **kwargs)

    def run(self, **kwargs):
        # update keymap and verbosity based on call arguments
        # this is only when the code is called from python session
        self.verbose = kwargs.pop("verbose", self.verbose)
        self.keymap["c"] = keyword("c", self.config)

        photoz = PhotoZ(self.keymap)
        autoadapt = (self.keymap["AUTO_ADAPT"]).split_bool("NO", 1)[0]
        if autoadapt:
            adapt_srcs = photoz.read_autoadapt_sources()
            a0 = photoz.run_autoadapt(adapt_srcs)
        else:
            a0 = []
            for _ in range(photoz.imagm):
                a0.append(0.0)

        opa_out = GalMag.read_opa()  # noqa: F841

        fit_srcs = photoz.read_photoz_sources()
        photoz.run_photoz(fit_srcs, a0)
        photoz.write_outputs(fit_srcs, int(time.time()))

        self.photoz = photoz
        return

    def end(self):
        super().end()


def main():  # pragma no cover
    runner = Zphota()
    runner.run()
    runner.end()


if __name__ == "__main__":  # pragma no cover
    main()
