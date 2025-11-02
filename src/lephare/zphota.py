import time
from contextlib import suppress

from ._lephare import (
    PhotoZ,
    keyword,
)
from .runner import Runner

__all__ = [
    "Zphota",
]

config_keys = {
    "verbose": "increase onscreen verbosity",
    "CAT_IN": "input catalog",
    "INP_TYPE": "F if Fluxes, M if magnitudes",
    "CAT_MAG": "AB or VEGA magnitude system of the input catalog (default AB)",
    "CAT_FMT": "MEME or MMEE, to precise the order of value and error columns",
    "CAT_LINES": "min,max range to subselect entries in the input catalog",
    "ZPHOTLIB": "comma-separated list of template magnitude libraries (GAL, STAR, or QSO), \
    without the suffix (to be found in $LEPHAREWORK/lib_mag)",
    "ADD_EMLINES": "Range of galaxy models for which considering emission lines contribution.",
    "EBV_RANGE": "E(B-V) min and max allowed in the library. Applied to all attenuation laws",
    "Z_RANGE": "Z min and max allowed in the GALAXY library",
    "PARA_OUT": "file listing the output keywords to be written on file (one keyword per line)",
    "CAT_OUT": "output file name",
    "CAT_TYPE": "LONG or SHORT, depending on whether zspec, context, and other columns are present \
    after the magnitudes/fluxes and error columns",
    "ERR_SCALE": "values to be added in quadrature to the errors in magnitude in each band\
    separated by a comma",
    "ERR_FACTOR": "values to scale the errors in flux (a single value for all bands)",
    "GLB_CONTEXT": "global context to define the bands to be used in the fit",
    "FORB_CONTEXT": "forbidden context to define the bands to be discarded from the fit",
    "MAG_ABS": "min,max absolute magnitudes allowed for GAL type (prior)",
    "MAG_ABS_QSO": "min,max absolute magnitudes allowed for QSO/AGN type (prior)",
    "MAG_REF": "reference band to compute absolute magnitudes when using the prior",
    "NZ_PRIOR": "band to be used in the n(z) prior; a second band can be given\
    in case the first one is absent in some objects",
    "ZFIX": "fix z to a known redshift, in case one is interested in inferring the best template",
    "Z_INTERP": "perform quadratic interpolation to refine best redshift solution",
    "Z_METHOD": "which redshift solution, between best and pdf median,\
    to use for physical parameter computation",
    "DZ_WIN": "Window search for second peak",
    "MIN_THRES": "lower threshold for second peak",
    "SPEC_OUT": "output the best spectrum of each object (YES, or the name of the directory\
    to be used",
    "FIR_LIB": "dedicated far infrared synthetic magnitudes files",
    "FIR_LMIN": "minimum value for lambda_mean/(1+z) (in microns) in FIR analysis",
    "FIR_CONT": "context for the FIR bands to be used in the FIR fit",
    "FIR_SCALE": "context for the FIR scale determination in the FIR fit",
    "FIR_FREESCALE": "consider FIR scaling free",
    "FIR_SUBSTELLAR": "substract the stellar component to the observed IR fluxes",
    "MABS_METHOD": "Method to compute absolute magnitudes\
    (if 1 use the most appropriated filter for observed flux, 3 the model)",
    "MABS_CONTEXT": "Context to decide which band can be used to derive the\
    absolute magnitude in method 1",
    "MABS_REF": "Fixed filter used to derived all abs mag (if MABS_METHOD=2)",
    "MABS_FILT": "List of fixed filters chosen to derive abs mag in all bands depending on\
    the redshift bins defined in MABS_ZBIN (if MABS_METHOD=4)",
    "MABS_ZBIN": "List of redshift bins associated with MABS_FILT",
    "RF_COLORS": "list of 4 band indexes out of which to compute the probability\
    distribution of the 2 corresponding rest frame colors (default=-1 which means\
    to not compute these color distributions)",
    "M_REF": "band index for which to compute the rest frame absolute magnitude\
    probability distribution",
    "APPLY_SYSSHIFT": "list of values (equal to the number of bands) to add to the predicted\
    magnitudes as zero-points  (substracted to the observed magnitudes).\
    These are also the output of the auto adaptation stage",
    "AUTO_ADAPT": "perform zero point auto adaptation based on the provided true redshifts",
    "ADAPT_BAND": "band number for which to check magnitude range for auto adaptation",
    "ADAPT_LIM": "min,max magnitude range for auto adaptation",
    "ADAPT_ZBIN": "min,max value of the true redshifts used for auto adaptation of the zero points",
    "PDZ_OUT": "output the PDF of each object",
    "PDZ_TYPE": "list of the PDFs to output if PDZ_OUT is set",
    "ADDITIONAL_MAG": "file basename (in $LEPHAREWORK/filt/) for the set of additional filters\
    for which to compute magnitudes. Thus a pass to the LePHARE `filter` stage needs to have\
    occurred on these filters.",
    "LIMITS_ZBIN": "Used to compute z_max. Redshifts used to split in N bins, separated by a coma.\
    Need N+1 values (start with the minimum redshift).",
    "LIMITS_MAPP_REF": "Used to compute z_max. Band in which the absolute magnitude is computed",
    "LIMITS_MAPP_SEL": "Used to compute z_max. Give the selection band in each redshift bin.\
    Need 1 or N values.",
    "LIMITS_MAPP_CUT": "Used to compute z_max. Apparent mag selection used in each redshift bin.\
    Need 1 or N values.",
    "CHI2_OUT": "Flag to output the chi2 value of each template (one file per source)",
}

nonestring = "NONE"


class Zphota(Runner):
    """The specific arguments to the Zphota class are
    verbose:
           add verbosity
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
           comma-separated list of template magnitude libraries \
without the suffix (to be found in $LEPHAREWORK/lib_mag)
    ADD_EMLINES:
           range of galaxy models for which considering emission lines contribution
    EBV_RANGE:
           E(B-V) min and max allowed in the library. Applied to all attenuation laws
    Z_RANGE:
           Z min and max allowed in the GALAXY library
    PARA_OUT:
           file listing the output keywords to be written on file (one keyword per line)
    CAT_OUT:
           output file name
    CAT_TYPE:
           LONG or SHORT, depending on whether zspec, context, and other columns are present \
after the magnitudes/fluxes and error columns
    ERR_SCALE:
           values to be added in quadrature to the errors in magnitude in each band, \
values separated by a comma
    ERR_FACTOR:
           values to scale the errors in flux (a single value for all bands)
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
           band to be used in the n(z) prior; a second band can be given \
in case the first one is absent in some objects
    ZFIX:
           fix z to a known redshift, in case one is interested in inferring the best template
    Z_INTERP:
           perform quadratic interpolation to refine best redshift solution
    Z_METHOD:
           which redshift solution, between best and pdf median, \
to use for physical parameter computation
    DZ_WIN:
           Window search for second peak
    MIN_THRES:
           lower threshold for second peak
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
           method to compute absolute magnitudes \
(e.g.  1 use the most appropriated filter for observed flux, 3 the model)
    MABS_CONTEXT:
           context to decide which band can be used to derive the \
absolute magnitude in method 1
    MABS_REF:
           fixed filter used to derived all abs mag (if MABS_METHOD=2)
    MABS_FILT:
           list of fixed filters chosen to derive abs mag in all bands depending on \
the redshift bins defined in MABS_ZBIN (if MABS_METHOD=4)
    MABS_ZBIN:
           list of Redshift bins associated with MABS_FILT
    RF_COLORS:
           list of 4 band indexes out of which to compute the probability \
distribution of the 2 corresponding rest frame colors (default=-1 which means \
to not compute these color distributions)
    M_REF:
           band index for which to compute the rest frame absolute magnitude \
probability distribution
    APPLY_SYSSHIFT:
           list of values (equal to the number of bands) to add to the predicted \
magnitudes as zero-points  (substracted to the observed magnitudes)
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
    ADDITIONAL_MAG:
           file basename (in $LEPHAREWORK/filt/) for the set of additional filters \
for which to compute magnitudes. Thus a pass to the LePHARE `filter` stage needs to have \
occurred on these filters.
    LIMITS_ZBIN:
            used to compute z_max. Redshifts used to split in N bins, separated by a coma. \
Need N+1 values (start with the minimum redshift).
    LIMITS_MAPP_REF:
            used to compute z_max. Band in which the absolute magnitude is computed
    LIMITS_MAPP_SEL:
            used to compute z_max. Give the selection band in each redshift bin. \
Need 1 or N values.
    LIMITS_MAPP_CUT:
            used to compute z_max. Apparent mag selection used in each redshift bin. \
Need 1 or N values
    CHI2_OUT
            Flag to output the chi2 value of each template (one file per source)
    """

    def update_help(self):
        """Add the specific Zphota help"""
        doc = "Execute LePHARE photo-z estimation"
        with suppress(Exception):
            self.parser.usage = doc
        self.__doc__ = doc + "\n"  # + inspect.getdoc(Zphota)

    def __init__(self, config_file=None, config_keymap=None, **kwargs):
        super().__init__(config_keys, config_file, **kwargs)

    def run(self, **kwargs):
        super().run(**kwargs)

        self.keymap["c"] = keyword("c", self.config)

        photoz = PhotoZ(self.keymap)

        # Compute offsets depending on the AUTO_ADAPT and APPLY_SYSSHIFT options (0 if none)
        adapt_srcs = photoz.read_autoadapt_sources()
        a0 = photoz.compute_offsets(adapt_srcs)

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
