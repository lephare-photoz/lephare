"""
LePHARE is a code for computing photometric redshifts and physical parameters

For full documentation, see: https://lephare.readthedocs.io/

Since 1999, LePHARE has been a code for computing photometric redshifts and physical parameters
by fitting spectral energy distributions (SEDs) to a dataset of photometric fluxes or apparent
magnitudes. LePHARE was originally written in Fortran (Arnouts et al. 1999; Ilbert et al. 2006).
It has been completely rewritten in C++ with a Python interface, with a larger team of developers
including J. Cohen-Tanugi and R. Shirley.

In order to run LePHARE we need to download auxiliary data such as filters, SEDs,
and attenuation curves which are not shipped with the code. These are explained
in more detail below but for simplicity you can download everything via Python (~1.3 Gb):

    import lephare as lp
    lp.data_retrieval.get_auxiliary_data(clone=True)

The following Python snippet is the most basic example to test that the installation has worked.
This will generate intermediate files and outputs stored in the cache or in user-defined storage
locations. You can also get an example notebook running this code here.

    import lephare as lp
    from astropy.table import Table

    # The following config is highly dependent on your input data and science goals
    # You can change it for your own needs
    config = lp.default_cosmos_config.copy()
    lp.prepare(config)

    # The following example table is in the lephare input format
    input_table = Table.read(f"{lp.LEPHAREDIR}/examples/COSMOS.in", format="ascii")

    # In the next command output is an astropy.table.Table object with the results
    output, _ = lp.process(config, input_table)

    # One can then inspect, for instance, the first 5 lines of output
    output[:5]

This workflow may take over ten minutes to run. To check that everything was successful,
this example should produce a 1-to-1 relationship between the spectroscopic redshift
output['ZSPEC'] and predicted redshift output['Z_BEST'].
"""
# ruff: noqa: E402
# ruff: noqa: F403
# ruff: noqa: F405

# Why is this global?
global LEPHAREDIR

from .data_manager import DataManager

dm = DataManager()
dm.configure_directories()  # noqa: F405
LEPHAREDIR = dm.LEPHAREDIR

from ._lephare import *
from ._version import *

# make LEPHAREDIR and LEPHAREWORK avaliable to the C++ codes
get_lephare_env()  # noqa: F405

from ._flt import *
from ._pdf import *
from ._photoz import *
from ._plot_utils import *
from ._spec import *
from .data_retrieval import *
from .default_cosmos_config import *
from .filter import *
from .filter_extinc import *
from .filterSvc import *
from .mag_gal import *
from .magSvc import *
from .prepare import *
from .process import *
from .runner import *
from .sedtolib import *
from .zphota import *

# Improved documentation for bound classes
oneElLambda.__doc__ = """
Class representing a single element in a spectral energy distribution (SED).

Parameters
----------
lambin : float
    The wavelength value for this element.
valin : float
    The flux or magnitude value associated with `lambin`.
oriin : int
    The origin or index identifier for this element.

Attributes
----------
lamb : float
    Wavelength of the element.
val : float
    Flux or magnitude value of the element.
ori : int
    Origin or index identifier.

Methods
-------
interp(previousEl, nextEl)
    Interpolates the value of this element given two neighboring elements.
"""

cosmo.__doc__ = """
Cosmology class for computing distances, times, and flux rescaling.

This class provides basic cosmological computations such as distance modulus,
metric distance, cosmic time, and flux rescaling, assuming a flat Lambda-CDM universe.

Parameters
----------
h0 : float, optional
    Hubble constant at z=0 in km/s/Mpc (default is 70).
om0 : float, optional
    Matter density parameter (default is 0.3).
l0 : float, optional
    Dark energy density parameter (default is 0.7).

Methods
-------
distMod(z)
    Compute the distance modulus for redshift `z`.
distMet(z)
    Compute the metric (comoving) distance for redshift `z`.
time(z)
    Compute the cosmic time at redshift `z`.
flux_rescaling()
    Compute the flux rescaling factor for this cosmology.
"""

opa.__doc__ = """
Class representing an opacity table for a given redshift.

This class loads and stores opacity data from a file, which can be used in
photometric and SED calculations.

Parameters
----------
red : float
    The redshift at which the opacity is defined.
opaFile : str
    Path to the opacity data file to read.

Attributes
----------
lamb_opa : ndarray
    Wavelength array from the opacity file.
opared : float
    Redshift value of the opacity data.

Methods
-------
read()
    Read the opacity data from the file specified during initialization.
"""

ext.__doc__ = """
Class representing an extinction curve.

This class stores extinction data as a function of wavelength and provides
methods to read extinction files, add elements, and access wavelength limits.

Parameters
----------
name : str
    Name of the extinction curve.
numext : int
    Number of extinction elements stored.

Attributes
----------
lamb_ext : ndarray
    Wavelength array for the extinction curve.
name : str
    Name of the extinction curve.
numext : int
    Number of elements in the extinction curve.
lmin : float, read-only
    Smallest wavelength stored.
lmax : float, read-only
    Largest wavelength stored.

Methods
-------
read(extFile)
    Read an extinction file and populate the internal arrays.
add_element()
    Add a new extinction element to the curve.
"""

keyword.__doc__ = """
Class representing a configuration keyword with a name and value.

This class provides utilities to parse and manipulate keyword values from configuration
files, allowing them to be split into different types (int, long, float, bool, etc.)
and to expand paths.

Attributes
----------
name : str
    Name of the keyword.
value : str
    Value of the keyword as a string.

Methods
-------
expand_path()
    Expand environment variables and user home (`~`) in the keyword value.
split_string()
    Split the keyword value into a list of strings.
split_int()
    Split the keyword value into a list of integers.
split_long()
    Split the keyword value into a list of long integers.
split_double()
    Split the keyword value into a list of doubles/floats.
split_bool()
    Split the keyword value into a list of booleans.
"""

flt.__doc__ = """
Class representing a photometric filter.

This class stores filter transmission data and provides methods to read filter files,
compute effective wavelengths, widths, and magnitudes, and access the transmission
curve.

Constructors
------------
flt(id, name, trans, calib)
    Create a filter with an integer ID, name, transmission type, and calibration type.

flt(lmin, lmax, nstep)
    Create a top-hat filter from lmin to lmax with nstep points.

Attributes
----------
name : str, read-only
    Name of the filter.
lmean : float, read-only
    Mean wavelength of the filter (computed from transmission curve).
dwidth : float, read-only
    Width of the filter.
lamb_trans : list of oneElLambda
    Transmission data points, each containing wavelength and transmission value.

Methods
-------
read(filename)
    Read filter information from a file.

read(stream)
    Read filter information from an input stream.

lambdaMean()
    Compute mean wavelength of the filter.

lambdaEff()
    Compute effective wavelength of the filter.

magsun()
    Compute magnitude of the Sun in this filter.

width()
    Compute the width of the filter.

lmin()
    Return the minimum wavelength of the filter.

lmax()
    Return the maximum wavelength of the filter.

data()
    Return the filter transmission data as a 2xN numpy array,
    with first row = wavelengths, second row = transmission values.
"""

SED.__doc__ = """
Class representing a Spectral Energy Distribution (SED).

This class stores the flux of a model as a function of wavelength and provides
methods to read SEDs, manipulate them, compute magnitudes and fluxes in filters,
resample, integrate, and handle extinction.

Constructors
------------
SED(name="", nummod=0, type="G")
    Create an empty SED with optional name, number of models, and type.

SED(name, tau, age, nummod, type, idAge)
    Create an SED with physical parameters: star formation timescale `tau`,
    age, number of models, type, and age ID.

SED(SED)
    Copy constructor.

Attributes
----------
lamb_flux : list of oneElLambda, read-only
    Wavelength and flux values of the SED.
extlawId : int, read-only
    Extinction law identifier applied to this SED.
ebv : float, read-only
    Color excess E(B-V) used for extinction.
name : str, read-only
    Name of the SED.
nummod : int, read-only
    Number of models.
mag : list of float, read-only
    Magnitudes computed for this SED.
index_z0 : int
    Index for z=0 reference in computations.

Methods
-------
string_to_object()
    Convert a string representation into an SED object.

is_gal()
    Return True if the SED corresponds to a galaxy.

is_star()
    Return True if the SED corresponds to a star.

is_qso()
    Return True if the SED corresponds to a quasar.

read()
    Read SED data from file.

size()
    Return the number of wavelength points in the SED.

integrateSED()
    Integrate the SED over wavelength.

resample()
    Resample the SED onto a new wavelength grid.

generateCalib()
    Generate a calibration vector for the SED.

rescale()
    Rescale the SED fluxes by a given factor.

compute_magnitudes()
    Compute magnitudes in specified filters.

compute_fluxes()
    Compute fluxes in specified filters.

emplace_back()
    Add a new element to the SED.

set_vector()
    Replace the SED internal vector with a new set of values.

readSEDBin(filename)
    Read a binary SED file.

writeSED(filename, fmt, options)
    Write the SED to a file in the specified format and options.

data()
    Return a 2xN numpy array with first row = wavelengths, second row = fluxes.
"""

StarSED.__doc__ = """
Class representing a stellar SED, derived from SED.

This class is used for stars and inherits all methods and attributes from SED.
It provides additional constructors for creating stellar SEDs from existing SEDs or
by specifying a name and number of models.

Constructors
------------
StarSED(sed: SED)
    Create a StarSED from an existing SED object.

StarSED(star_sed: StarSED)
    Copy constructor.

StarSED(name: str, nummod: int = 0)
    Create a new StarSED with the given name and optional number of models.
"""

QSOSED.__doc__ = """
Class representing a quasar SED, derived from SED.

This class is used for QSOs and inherits all methods and attributes from SED.
It provides constructors for creating QSO SEDs from existing SEDs or by name and
number of models.

Constructors
------------
QSOSED(sed: SED)
    Create a QSOSED from an existing SED object.

QSOSED(qso_sed: QSOSED)
    Copy constructor.

QSOSED(name: str, nummod: int = 0)
    Create a new QSOSED with the given name and optional number of models.
"""

GalSED.__doc__ = """
Class representing a galaxy SED, derived from SED.

This class is used for galaxy spectral energy distributions. It inherits all methods
and attributes from SED and provides additional methods to handle emission lines,
nebular continuum, K-corrections, and physical properties.

Constructors
------------
GalSED(sed: SED)
    Create a GalSED from an existing SED object.

GalSED(gal_sed: GalSED)
    Copy constructor.

GalSED(name: str, nummod: int = 0)
    Create a new GalSED with the given name and optional number of models.

GalSED(name: str, tau: float, age: float, format: str, nummod: int, type: str, idAge: int)
    Create a new GalSED with physical parameters, format, type, and age ID.

Methods
-------
SEDproperties()
    Compute or return properties of the SED.

add_neb_cont()
    Add nebular continuum emission to the galaxy SED.

generateEmEmpUV()
    Generate empirical UV emission component.

generateEmEmpSFR()
    Generate empirical SFR-related emission component.

generateEmPhys()
    Generate emission lines based on physical models.

sumEmLines()
    Sum all emission lines into the SED.

kcorrec()
    Compute K-corrections for the SED.

rescaleEmLines()
    Rescale emission lines by a given factor.

zdepEmLines()
    Compute redshift-dependent emission lines.

calc_ph()
    Compute photometric properties from the SED.
"""

PhotoZ.__doc__ = """
Class for computing photometric redshifts and related quantities.

This class handles fluxes, photometric libraries, redshift grids, and
all the processing needed to compute photo-z values from input data.

Attributes
----------
flux : list or array-like
    Fluxes for each source in the main bands.
fluxIR : list or array-like
    Fluxes in the infrared bands.
imagm : list or array-like
    Input magnitudes for sources.
fullLib : list
    Full library of SEDs used for fitting.
fullLibIR : list
    Full IR library of SEDs.
zLib : list
    Redshift grid used for photo-z computation.
allFilters : list
    List of filters used in the analysis.
gridz : list
    Redshift grid values.
outkeywords : list
    Keywords of output parameters.
outpara : list
    Output parameters computed for each source.
pdftype : int
    Type of probability density function used.
outputHeader : str
    Header for output files.

Methods
-------
read_autoadapt_sources()
    Read sources prepared for auto-adaptive photo-z computation.

read_photoz_sources()
    Read sources for standard photo-z computation.

prep_data(sources)
    Prepare data for processing.
    Can be called with a single source or a vector of sources.

run_autoadapt()
    Run the auto-adaptive photo-z workflow.

run_photoz()
    Run the standard photo-z computation.

write_outputs()
    Write computed outputs to files.

validLib()
    Validate the photo-z library.

compute_offsets()
    Compute offsets between observed and template fluxes.
"""

onesource.__doc__ = """
Class representing a single astronomical source for photometric redshift analysis.

The `onesource` class contains all the photometric data, SED fits, fluxes, magnitudes,
physical parameters, and results associated with a single source. It provides methods
for reading input, fitting templates, computing PDFs, and generating outputs.

Constructors
------------
onesource()
    Default constructor.

onesource(id: int, fluxes: list of float)
    Constructor initializing with an ID and a vector of fluxes.

Methods
-------
setPriors(priors)
    Set prior probabilities for this source.

readsource(filename, fluxes, flux_errs, index, zmin, filterType)
    Read input source data from file or arrays.

set_verbosity(level)
    Set verbosity for logging.

get_verbosity() -> int
    Get current verbosity level.

fltUsed()
    Return the list of filters used for fitting.

convertFlux(flux, filter)
    Convert flux into desired units.

convertMag(mag, filter)
    Convert magnitude into flux.

rescale_flux_errors()
    Rescale flux errors to account for systematics.

keepOri()
    Keep original magnitudes for reference.

adapt_mag()
    Adapt magnitudes to current fitting configuration.

fit()
    Fit templates to the source photometry.

mode()
    Return the mode of the redshift PDF.

rm_discrepant()
    Remove discrepant flux points from the data.

generatePDF()
    Generate the redshift probability density function.

interp()
    Interpolate fluxes or magnitudes within templates.

uncertaintiesMin()
    Compute minimum-based uncertainties on parameters.

uncertaintiesBay()
    Compute Bayesian uncertainties on parameters.

secondpeak()
    Identify second peak in PDF if present.

interp_lib()
    Interpolate template library as needed.

absmag()
    Compute absolute magnitudes.

limits()
    Compute observational or physical limits.

computePredMag()
computePredAbsMag()
    Compute predicted magnitudes and absolute magnitudes.

computeEmFlux()
    Compute flux contribution from emission lines.

fltUsedIR()
    Return filters used in the infrared bands.

substellar()
    Determine if source is substellar.

generatePDF_IR()
    Generate PDF for IR data.

write_out()
writeSpec()
writeFullChi()
best_spec_vec()
    Methods to write outputs or return best-fit vectors.
"""

PDF.__doc__ = """
Class representing a probability density function (PDF) for a single parameter.

This class is used to compute, manipulate, and analyze PDFs derived from chi-squared
fits or other likelihood distributions. It provides tools for normalization, extremum
finding, credible intervals, and other statistical analyses.

Constructor
-----------
PDF(min: float, max: float, size: int)
    Create a PDF over the range [min, max] with 'size' bins.

Methods
-------
normalization()
    Normalize the PDF so that its integral equals 1.

chi2toPDF(chi2)
    Convert chi-squared values into a PDF.

chi2mini()
    Find the minimum chi-squared value and corresponding index.

uncMin()
    Compute uncertainties around the minimum of the PDF.

index(value)
    Return the index of a given x-axis value.

get_max()
    Return the maximum PDF value.

get_maxid()
    Return the index of the maximum PDF value.

secondMax()
    Return the second-highest PDF value and its index.

size()
    Return the number of bins in the PDF.

cumulant()
    Return the cumulative distribution of the PDF.

levelCumu2x(level)
    Convert a cumulative probability level to the corresponding x value.

credible_interval(level)
    Compute the credible interval for a given probability level (e.g., 68%, 95%).

confidence_interval(level)
    Compute the confidence interval for a given probability level.

improve_extremum()
    Refine the location of maxima or minima using interpolation.

Attributes
----------
vPDF : array-like
    Values of the PDF at each bin.
xaxis : array-like
    The x-axis values corresponding to each PDF bin.
chi2 : array-like
    Chi-squared values used to generate the PDF.
secondX : float
    Position of the second-highest PDF peak.
secondP : float
    Value of the second-highest PDF peak.
ind : int
    Index of the primary maximum.
secondInd : int
    Index of the secondary maximum.
"""
