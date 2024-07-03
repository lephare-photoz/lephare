Detailed LePHARE user manual
============================

LePHARE was written by Stéphane Arnouts & Olivier Ilbert (Laboratoire 
d’Astrophysique de Marseille) and later by Johann Cohen-Tanugi (Laboratoire 
Univers et Particules de Montpellier). Full details of all contributors can be 
found on the main GitHub `here <https://github.com/lephare-photoz/lephare/>`_.

This section contains more details concerning the 
internal structure of the code and how to run it using the legacy command line arguments.
.

.. _`sect:introduction`:

Introduction
-------------

*LePHARE* is a set of ``C++`` programs to compute photometric
redshifts ( :math:`z_\mathrm{phot}` ) for galaxies and AGN, and galaxy
physical parameters [1]_ by fitting spectral energy distributions (SEDs)
to a dataset of photometric fluxes or apparent magnitudes. Stellar
templates are fitted too. It is based on the previous Fortran version of
*LePHARE*, a code described in Arnouts et al. (1997) and Ilbert et al.
(2006). *LePHARE++* package is composed of four parts:

-  a preliminary step to read the SED and build the SED library;

-  a preliminary step to read the input filters set and build the filter
   library ;

-  a program to compute apparent magnitudes in the filter set for each
   template of the library, along a grid of redshift, adding dust
   attenuation and nebular emission lines and storing the results. This
   step allows the user to extract basic information relative to the
   filters (e.g. :math:`\lambda_{mean}`, AB-corrections, attenuation)
   and SEDs (e.g. k-corrections, color-color diagrams, etc.).

-  The photometric redshift code based on a :math:`\chi^2` fitting
   method. This part can also be used to compute physical parameters.


.. Comment out old download instructions
.. Download and installation
.. ^^^^^^^^^^^^^^^^^^^^^^^^^

.. | The basic package is available from GitHub. Each file is briefly described in Appendix
..   A. Additional SED libraries (in separate tarballs) are also available
..   on the same webpage.
.. | Before starting, you must set two environment variables:

.. -  ``$LEPHAREDIR`` is the root directory of the software (e.g.,
..    ``/home/yourname/LEPHARE/``)

.. -  ``$LEPHAREWORK`` is the path of a new directory which will be created
..    when compiling the code (e.g. ``$LEPHAREDIR/work/``).
..    ``$LEPHAREWORK`` will contain libraries created by *LePHARE++* (it
..    can be placed anywhere). The code will use ``$LEPHAREDIR/work/`` if
..    it doesn’t exist.

.. | These two environment variables could be set definitively in the usual
..   files depending on your shell (e.g., ``.bash_profile`` or
..   ``.csh_envlph``) as any environment variable.
.. | Once downloaded (and unpackaged) *LePHARE++* , enter in the directory
..   ``source``. Then open the file ``Makefile`` and modify its options
..   according to your OS if needed (see notes below).
.. | Once your OS is ready to compile the code, execute the commands:
.. | $ make clean
.. | $ make
.. | if you recompile once again from scratch. Note that the directory
..   ``$LEPHAREWORK`` will be created during this phase.
.. | You can create the directory ``$LEPHAREWORK`` separately with the
..   command:
.. | $ make work
.. | You can use ``Makefile`` to build a gzipped file containing all the
..   source files by the command:
.. | $ make archive
.. | Notes on the compiler options:
.. | With the default set-up, the compiler is GNU ``g++`` with basic
..   optimization flags (``-g -std=c++11 -O3``). If you want to enable
..   parallelization (OPENMP), other options should be included (e.g.,
..   ``-lpthread``). MacOS users must be aware that the ``g++`` installed
..   via XCode is a wrapper of ``clang`` and can raise errors while
..   compiling; the easiest solution is to install the original GNU
..   compiler (in the GCC package by Free Software Foundation) available
..   through Homebrew or MacPorts. Note also that the code has been tested
..   with GCC versions :math:`>6` (6.3.0 and 7.3.0); older versions may
..   result in critical errors.

.. Example with one run
.. ^^^^^^^^^^^^^^^^^^^^

.. | We provide an example of the whole procedure for a test catalog
..   included in the package in
.. | ``$LEPHAREDIR/test/``. For this example, we used the COSMOS2015
..   catalogue (Laigle et al. 2016) but limited to the zCOSMOS bright
..   sample, for which also spectroscopy is provided. You can find more
..   detailed examples in ``$LEPHAREDIR/examples/``. Different
..   configurations are tested, showing how to configure the code for this
..   application. All the key steps are considered.

Structure
^^^^^^^^^

The structure of the package is illustrated in Fig `1 <#fig:skim>`__.
Each step is described in the following sections. First, one has to
define the filters used in the input galaxy catalog (program
``filter``). The other initial step is constructing a library of
rest-frame galaxy SEDs from synthetic models or observed spectra
(program ``sedtolib``). Then, another program builds a grid of apparent
magnitudes at different redshifts, adding dust attenuation and nebular
emission lines to each SED in the library (program ``mag_gal``). The
main *LePHARE++* program (program ``zphota``) will use the file
containing such a grid to fit the photometry from the input catalog.

.. image:: figures/lephare_skim.png
  :width: 700
  :alt: Alternative text
  :name: fig:skim



Syntax
^^^^^^

| All the programs in the suite can be run from a Unix shell with the
  following syntax:
| $ :math:`<`\ program\ :math:`>` -c :math:`<`\ configfile\ :math:`>` [
  -Parameter *value*...]
| where :math:`<`\ program\ :math:`>` is the name of the program,
  followed by a configuration file called with the **-c** option. The
  various code options are defined in this file but can also be given
  through additional instructions in the command line. Using such an
  optional list of parameters, any **-Parameter** *value* statement
  overrides the values in the configuration file.

The parameters associated with the various programs can be included in a
single configuration file (e.g., ``zphot.para``) You can store your
parameter file where you want (e.g., in the directory where you run the
code or in ``$LEPHAREDIR/config/``) to keep configuration files of
different runs. Configuration files must be in ASCII format, compliant
with the following rules:

-  Only one parameter per line, with the syntax: PARAMETER_NAME value(s)

-  Comment line starts with “#”.

-  Depending on the parameter, values can be Float, Integer, or String
   (without quotation marks).

-  When a parameter accepts multiple values, these must be comma
   separated (no space).

-  When a parameter accepts a file location (as a String), the path can
   include environmental variables (``$HOME`` and ``$LEPHAREDIR``).

-  Some parameters are mandatory, *LePHARE++* will print out an error
   message if they are not set (either in the configuration file or via
   the command line)

-  Other parameters can be omitted (*LePHARE++* will assign a default
   value to them)

| In the next sections, we will mark the mandatory parameters with an
  asterisk ("\*").
| You can use the option ``VERBOSE NO`` if you don’t want that
  ``mag_gal`` and ``zphota`` display the template or sources computed.
  It could be helpful if run in batch mode.

.. _models:

Rest-frame SED libraries through ``sedtolib``
---------------------------------------------

Overview
^^^^^^^^

A set of libraries for stars, galaxies, and quasars are available in
$LEPHAREDIR/sed/STAR, $LEPHAREDIR/sed/GAL, $LEPHAREDIR/sed/QSO [2]_
directories and organized in different sub-folders.

Each sub-folder contains a specific collection of SED files, described
in a README (how those SEDs were built, etc.), and a file (usually with
the suffix ``.list``) listing the relative path of the SED files to be
used as input for ``sedtolib``. For STAR and QSO and most of the
galaxies, SEDs are written in ASCII, with :math:`\lambda(\AA)`,
flux[:math:`erg/s/\AA/cm^2`], with increasing :math:`\lambda`\  [3]_.
For Galaxy, in addition to empirical SEDs, output files from stellar
synthesis population models (Pegase and BC03) with a more complex format
can also be used by adding a specific character after the file name in
the SED list file (see end of section 2.2.2).

 ``sedtolib`` program 
^^^^^^^^^^^^^^^^^^^^^

The program **sedtolib** is used to build the different STAR, QSO and
GALAXY libraries from a list of SED files. The goal of this program is
to generate from different kinds of SEDs (star/AGN/galaxy) with various
original formats (ASCII, binary), a unique binary file with direct
access that can be easily read in the following steps. The binary output
file (\*.bin) is saved in the directory $\ *LEPHAREWORK*/lib_bin/ with
an attached doc file (\*.doc) and a file with physical information
(\*.phys) for galaxies. The new SED format is
(:math:`\lambda(\AA)`,flux[:math:`erg/s/\AA/cm^2`]). For models with
input SEDs expressed in luminosity or energy
(:math:`L_{\odot}/\AA`,\ :math:`\nu L_{\nu}`,...), like PEGASE, GISSEL,
or the FIR libraries, the SED are converted in flux
(:math:`erg/s/cm^2/\AA`).

Syntax and parameter values
~~~~~~~~~~~~~~~~~~~~~~~~~~~

| Specific parameters have been duplicated for the STAR, QSO, and
  GAL(axy) categories with different names to simplify this algorithm
  section. The option -t allows you to specify if galaxy (G), star (S),
  or QSO (Q) parameters have to be read.
| The syntax is:
| :math:`\%` **sedtolib** -t G [or Q or S] -c zphot.para

+-------------+--------+---------+----------------------------------+
| parameter   | type   | default | description                      |
+=============+========+=========+==================================+
| XXX_SED(\*) | string | —-      | Full pathname of file with the   |
|             |        |         | list of selected SED files       |
+-------------+--------+---------+----------------------------------+
|             | (n=1)  |         |                                  |
+-------------+--------+---------+----------------------------------+
| XXX_LIB(\*) | string | —-      | Name of the output binary        |
|             |        |         | library (with no extension)      |
+-------------+--------+---------+----------------------------------+
|             | (n=1)  |         | Files *$XXX_LIB*.bin,            |
|             |        |         | *$XXX_LIB*.doc and               |
|             |        |         | *$XXX_LIB*.phys                  |
+-------------+--------+---------+----------------------------------+
|             |        |         | saved in                         |
|             |        |         | $\ *LEPHAREWORK*/lib_bin/        |
+-------------+--------+---------+----------------------------------+
| XXX_FSCALE  | float  | 1.0     | Flux scale to be applied to each |
|             |        |         | SED in the list                  |
+-------------+--------+---------+----------------------------------+
|             | (n=1)  |         |                                  |
+-------------+--------+---------+----------------------------------+
| SEL_AGE     | string | NONE    | Full pathname of file with a     |
|             |        |         | list of ages (Gyr)               |
+-------------+--------+---------+----------------------------------+
|             | (n=1)  |         | to be extracted from GISSEL or   |
|             |        |         | PEGASE SEDs.                     |
+-------------+--------+---------+----------------------------------+
| AGE_RANGE   | float  | —–      | Range of age (Gyr)               |
+-------------+--------+---------+----------------------------------+
|             | (n=2)  |         |                                  |
+-------------+--------+---------+----------------------------------+

| 
| The extracted text from zphota.para, related to the **sedtolib**
  task.The parameter value "XXX" means either GAL or QSO or STAR. Note
  that SEL_AGE and AGE_RANGE are relevant only when using templates
  including an age (e.g. BC03).

Building libraries from a list of SEDs 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| The easiest is to take a predefined list of SED in the existing
  subdirectories and look at the README file.
| *For stars ($LEPHAREDIR/sed/STAR)*, SEDs are available in the
  subdirectories :
| :math:`\bullet` PICKLES/: 131 stellar SEDs from Pickles (1998)
| :math:`\bullet` BD/: Low mass stars library from Chabrier et al.
  (2000)
| :math:`\bullet` BD_NEW/: Brown dwarfs library from Baraffe et al.
  2015, Morley et al. 2012, 2014
| :math:`\bullet` LAGET/: (missing REF)
| :math:`\bullet` WD/: 4 white dwarfs from Bohlin et al. (1995)
| :math:`\bullet` SPEC_PHOT: Spectro-Photometric standards from Hamuy et
  al. (1992, 1994)
| *For QSOs ($LEPHAREDIR/sed/QSO)*, there is a list of observed spectra
  from different authors and some synthetical QSOs listed in the
  subdirectory (synth/). In particular, a list of templates was
  successfully used for computing the photometric redshift of the *XMM*
  and *Chandra* AGN identified in COSMOS. In short, the library includes
  pure QSO and hybrid templates obtained by combining galaxies with
  various AGN and QSO templates with different relative ratios. The
  details of the template construction are outlined in Salvato et al.
  (2009). Note that, unlike for galaxies, the templates to be used in
  QSO depend on the type of AGN and QSO to be fitted (see Salvato et al
  2011, Fotopoulou et al. 2012, Hsu et al. 2014, Ananna et al. 2017)
| *For galaxies ($LEPHAREDIR/sed/GAL)*, SEDs are available in the
  following subdirectories:
| :math:`\bullet` CFHTLS_SED/: 66 SEDs used for CFHTLS photo-z paper
  (Arnouts et al. 2007)
| :math:`\bullet` COSMOS_SED/: 31 SEDs used for COSMOS photo-z paper
  (Ilbert et al. 2009, 2013, Salvato et al. 2011, Dahlen et al. 2013)
| :math:`\bullet` CWW_KINNEY/: original CWW and Kinney spectra
| :math:`\bullet` BC03_CHAB/: SEDs from the BC03 library. These
  templates are derived with exponentially declining Star Formation
  Histories.
| :math:`\bullet` BC03_CHAB_DELAYED/: SEDs from the BC03 library. These
  templates are derived with delayed Star Formation Histories.
| *For Far-Infrared (FIR) SEDs ($LEPHAREDIR/sed/GAL)*, different SEDs
  are available :
| :math:`\bullet` CHARY_ELBAZ/: 105 FIR templates for different
  luminosity
| :math:`\bullet` DALE/ : 64 FIR templates
| :math:`\bullet` LAGACHE/: 46 FIR templates
| :math:`\bullet` SK06/ : different set of starburst models based on
  Siebenmorgen &Krugel (2006)
| Note that for the first 3 libraries (CHARY-ELBAZ, DALE, LAGACHE), we
  have subtracted a stellar component from their SEDs to get only the
  dust contribution at the shortest wavelengths.
| To know the format of the SEDs that are used in your list, an
  additional character must be specified after each SED file, allowing
  you to mix in one list of different types of galaxy SEDs. For example,
  you could prepare a new list which includes:
| BC03_CHAB/bc2003_lr_m52_chab_tau03_dust00.ised_ASCII BC03
| BC03_CHAB/bc2003_lr_m62_chab_tau03_dust00.ised_ASCII BC03
| COSMOS_SED/Ell1_A_0.sed
| COSMOS_SED/Ell2_A_0.sed
| In each list, it is possible to comment a template with #.
| For ASCII SED file, no character is required. The character **BC03**
  is used for the Bruzual and Charlot 2003 models. For the BC03
  templates, the file is in ASCII for the C++ version of LePhare, to
  avoid the problem of portability between various systems.
| For the list with FIR SEDs, the character **LW** (as for Long
  Wavelength) is required.

Physical information for the galaxies 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| For the galaxy templates, an additional file is generated by the
  program ``sedtolib`` with some physical properties (\*.phys). This
  information will be used when running the photo-z code to derive
  physical parameters. It contains the following parameters:
| Model Age :math:`L_{UV}` :math:`L_R` :math:`L_K` :math:`L_{IR}` Mass
  SFR Metallicity Tau :math:`D_{4000}`
| where
| Age is expressed in yr
| :math:`L_{UV}` is NUV monochromatic luminosity (Log([erg/s/Hz]))
  (:math:`\int_{2100}^{2500} L_{\lambda} d\lambda /400 * 2300^2/c` ))
| :math:`L_R` is optical r monochromatic luminosity (Log([erg/s/Hz]))
  (:math:`\int_{5500}^{6500} L_{\lambda} d\lambda /1000 * 6000^2/c` ))
| :math:`L_K` is NIR K monochromatic luminosity (Log([erg/s/Hz]))
  (:math:`\int_{21000}^{23000} L_{\lambda} d\lambda /2000 * 22000^2/c`
  ))
| :math:`L_{IR}` is the IR luminosity (Log([:math:`L_{\odot}`]))
| Mass is the stellar mass (:math:`M_{\odot}`), .i.e. the mass truly in
  stars (not the integral of the SFH)
| SFR is the ongoing star formation rate (:math:`M_{\odot}/yr`)
| Metallicity is the Gas metallicity of the galaxy
| Tau is the e-folding parameter for a star formation history with
  SFH=exp(-t/tau) (yr)
| :math:`D_{4000}` is the 4000A break measured as in Bruzual 1983
  (:math:`D_{4000}= \int_{4050}^{4250} F_{\lambda} d\lambda / \int_{3750}^{3950} F_{\lambda} d\lambda`)
| If not available, the parameters are set to -99.
| The IR luminosity (:math:`L_{IR}`) is derived using LW libraries. For
  the Infra-red libraries ( LW: Dale, Lagache, Chary-Elbaz, Siebenmorgen
  & Krugel) the IR luminosity is measured from 8 to 1000 microns. These
  luminosities may be slightly different then the ones quoted by the
  authors due to the different definitions of the :math:`L_{IR}`
  integration limit and because (at least for Dale, Lagache, and
  Chary-Elbaz) we have subtracted the underlying stellar component from
  the original SEDs.

Adding libraries
~~~~~~~~~~~~~~~~

New SEDs can be easily added to the current ones. They must be located
in the appropriate directory (GAL/STAR/QSO). If they are ASCII files
they must be in :math:`\lambda(\AA)`, flux[:math:`erg/s/\AA/cm^2`], with
increasing :math:`\lambda`.

example
^^^^^^^

| G **-c** zphot.para **-GAL_SED**
  $LEPHAREDIR/sed/GAL/CFHTLS_SED/CFHTLS_MOD.list **-GAL_LIB** LIB_CFHTLS
| This command reads the list of galaxy templates given by the keyword
  **-GAL_SED** (as indicated by **-t** G).
| A binary file LIB_CFHTLS.bin with a LIB_CFHTLS.doc and LIB_CFHTLS.phys
  files are saved in $LEPHAREWORK/lib_bin/.
| S **-c** zphot.para **-STAR_SED** $LEPHAREDIR/sed/STAR/STAR_MOD.list
  **-STAR_LIB** LIB_STAR
| This command reads the list of star templates given by the keyword
  **-STAR_SED** (as indicated by **-t** S).
| A binary file LIB_STAR.bin and a LIB_STAR.doc file are saved in
  $LEPHAREWORK/lib_bin/.
| Q **-c** zphot.para **-QSO_SED** $LEPHAREDIR/sed/QSO/QSO_MOD.list
  **-QSO_LIB** LIB_QSO
| This command reads the list of QSO templates given by the keyword
  **-QSO_SED** (as indicated by **-t** Q).
| A binary file LIB_QSO.bin and a LIB_QSO.doc file are saved in
  $LEPHAREWORK/lib_bin/.
| An example of a misleading command:
| S **-c** zphot.para **-GAL_SED**
  $LEPHAREDIR/sed/GAL/CWW_KINNEY/CWW_MOD.list **-GAL_LIB** LIB_CWW
| This command will work and will read the galaxy templates given by
  **-GAL_SED** but will interprete them as stars rather than galaxies
  because the option is set to S: **-t S** !
| The parameters passed in the command line can also be changed in the
  configuration (zphot.para) file, except -t and -c.

.. _`sec:filter`:

Filters 
-------

Several sets of filters, from different telescopes, are available in the
directory ``$LEPHAREDIR/filt/``. You could find most of the standard
filters (like the Johnson-Kron-Cousins in ``filt/jkc``). New set of
filters can be added there.

By default, the filters are stored in ``$LEPHAREDIR/filt/``, but you
could change this repository using the keyword ``FILTER_REP``. With this
keyword, you could indicate if you prefer to store the filters in
another directory.

Description and outputs
^^^^^^^^^^^^^^^^^^^^^^^

| The program **filter** puts together a list of filter response curves,
  and applies some transformations according to the nature of the
  filters. The resulting file in the directory $\ **LEPHAREWORK/filt/**.

.. _syntax-and-parameter-values-1:

Syntax and parameter values
^^^^^^^^^^^^^^^^^^^^^^^^^^^

| The syntax is : :math:`\%` **filter -c** zphot.para
| The following parameters are considered:

+----------------+----------------+----------------+----------------+
| Parameters     | type           | default        | description    |
+================+================+================+================+
| FILTER_REP     | string         | $LE            | Name of the    |
|                |                | PHAREDIR/filt/ | repository     |
|                |                |                | containing the |
|                |                |                | filters.       |
+----------------+----------------+----------------+----------------+
|                | (n\            |                |                |
|                |  :math:`=`\ 1) |                |                |
+----------------+----------------+----------------+----------------+
| F              | string         | —-             | filter files   |
| ILTER_LIST(\*) |                |                | separated by a |
|                |                |                | comma.         |
+----------------+----------------+----------------+----------------+
|                | Nfilt not      |                |                |
|                | limited        |                |                |
+----------------+----------------+----------------+----------------+
| TRANS_TYPE     | float          | 0              | Filter         |
|                |                |                | transmission   |
|                |                |                | type: 0=       |
|                |                |                | Energy; 1=     |
|                |                |                | Photon         |
+----------------+----------------+----------------+----------------+
|                | n=1 or n=Nfilt |                |                |
+----------------+----------------+----------------+----------------+
| FILTER_CALIB   | integer        | 0              | Filter         |
|                |                |                | calibration    |
|                |                |                | for long       |
|                |                |                | wavelengths    |
|                |                |                | [0-def].       |
+----------------+----------------+----------------+----------------+
|                | n=1 or n=Nfilt |                |                |
+----------------+----------------+----------------+----------------+
| FILTER_FILE    | string         | filter         | Name of the    |
|                |                |                | file with all  |
|                |                |                | combined       |
|                |                |                | filters .      |
+----------------+----------------+----------------+----------------+
|                | (n\            |                | It is saved in |
|                |  :math:`=`\ 1) |                | $\ **LEPHAR    |
|                |                |                | EWORK/filt/**. |
+----------------+----------------+----------------+----------------+

.. _`sec:filter`:

Parameter descriptions
^^^^^^^^^^^^^^^^^^^^^^

| : all the filter names must be separated by a comma. We assume that
  all the filter files are located in the directory
  **$LEPHAREDIR/filt/**, except if the keyword **FILTER_REP** is
  specified. When writing the set of filters to be used, only the
  pathname after the common string **$LEPHAREDIR/filt/** should be
  specified.
| : Type of the transmission curve for each filter, separated by a
  comma. The number of arguments should match the number of filter but
  if only value is given, which will be use for all the filters.
| The transmissions (:math:`T_{\lambda}`) are dimensionless (in % ),
  however they refer either to a transmission in Energy or Photon which
  will slightly modify the magnitude estimates. The magnitude is :

  .. math:: mag(*) = -2.5 \log_{10} \frac{\int F_{\lambda}(*) R_{\lambda} d\lambda}{\int F_{\lambda}(Vega) R_{\lambda} d\lambda}

  If the transmission curve (:math:`T_{\lambda}`) corresponds to energy
  then :math:`R_{\lambda}=T_{\lambda}`,
| If the transmission curve (:math:`T_{\lambda}`) corresponds to number
  of photons (:math:`N_{\varphi}`) then
  :math:`R_{\lambda}= \lambda T_{\lambda}` :

  .. math::

     N_{\varphi} =  \frac{ F_{\lambda} d\lambda }{h\ \nu} = \frac{F_{\lambda} \lambda d\lambda }{h\ c} \rightarrow  
      mag(*)=-2.5 \log_{10} \frac{\int F_{\lambda}(*) \lambda T_{\lambda} d\lambda}{\int F_{\lambda}(Vega) \lambda T_{\lambda} d\lambda}  \rightarrow  R_{\lambda}=\lambda T_{\lambda}

  When building the filter library, the filter shape is changed with
  respect to the original one as follows :

  .. math:: R_{\lambda}=T_{\lambda} ( \frac{\lambda}{< \lambda >})^{tt}

  , where :math:`tt` is the value of TRANS_TYPE parameter and
  :math:`< \lambda >` is the mean wavelength of the filter.
| The modification of filter shape can be significant for long
  wavelength filters and when the filter is broad. Nevertheless it is
  often not the dominant source of errors with respect to other
  uncertainties relative to QE-CCD, telescope transmission, atmospheric
  extinction shape etc...
| In the output filter file specified by the keyword **FILTER_FILE**, we
  save the values (:math:`\lambda (\AA)`,\ :math:`R_{\lambda}`).
| : This keyword allow to consider specific calibrations at long
  wavelengths in order to apply a correction factor to the original flux
  estimated by LEPHARE (see section `3.5 <#sec:filtcalib>`__ for more
  details).
| We define the correction factor as
  fac_corr\ :math:`=\frac{\int  R_{\nu} d\nu}{\int \frac{B_{\nu}}{B_{\nu_0}} R_{\nu} d\nu}= \frac{\int  R_{\lambda} d\lambda/\lambda^2}{1/\lambda_0^2 \int \frac{B_{\lambda}}{B_{\lambda_0}} R_{\lambda} d\lambda}`,
  where :math:`B_{\nu}` is the reference spectrum used to calibrate the
  filters and :math:`\lambda_0` is the effective wavelength defined as
  :math:`\lambda_{0}= \frac{\int R_{\lambda} B_{\lambda} \lambda d\lambda}{\int R_{\lambda}  B_{\lambda}  d\lambda}`.
| The value of **FILTER_CALIB** allows to describe different
  combinations of :math:`\nu_0` and :math:`B_{\nu}`:
| : :math:`\frac{B_{\nu}}{B_{\nu_0}}=1` or :math:`B_{\nu}=ctt`. This is
  the default value used in LEPHARE.
| : :math:`\nu B_{\nu}=ctt`. This describes the SPITZER/IRAC, ISO
  calibrations
| : :math:`B_{\nu}=\nu`. This describes the sub-mm calibrations.
| : :math:`B_{\nu}=`\ black body at T=10,000K.
| : A mix calibration with :math:`\nu_0` defined from
  :math:`\nu B_{\nu}=ctt` and the flux estimated as
  :math:`B_{\nu}=`\ black body at T=10,000K. This appears to be the
  adopted scheme for the SPITZER/MIPS calibration.
| : Similar mix calibration with :math:`\nu_0` defined from
  :math:`\nu B_{\nu}=ctt` and the flux estimated as :math:`B_{\nu}=\nu`.
  This may reflect the SCUBA calibration.

Filter informations
^^^^^^^^^^^^^^^^^^^

Standard filter informations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| As an example, using default values listed in the configuration file
  zphot.para.

+--------------+------------------------------------------------------+
| FILTER_LIST  | tmp/f300.pb,tmp/f450.pb,                             |
|              | tmp/f606.pb,tmp/f814.pb,tmp/Jbb.pb,tmp/H.pb,tmp/K.pb |
+--------------+------------------------------------------------------+
| TRANS_TYPE   | 0                                                    |
+--------------+------------------------------------------------------+
| FILTER_CALIB | 0                                                    |
+--------------+------------------------------------------------------+
| FILTER_FILE  | HDF.filt                                             |
+--------------+------------------------------------------------------+
|              |                                                      |
+--------------+------------------------------------------------------+

| Run the program ::math:`\%` **filter -c** zphot.para.
| It generates the file HDF.filt by combining all filters and saved it
  in $\ *LEPHAREWORK*/filt.
| It returns informations about the filters on screen. Another stand
  alone program allows also to read informations about existing filter
  list (with :math:`\%` **filter_info -f** HDF.filt).
| The following informations are written on the screen :

+-------+----+-------+-------+-------+-------+-------+-------+-------+-----+-------+-------+
| #NAME | ID | :m    | :ma   | FWHM  | ABcor | TGcor | VEGA  | :m    | CAL | :mat  | Fac   |
|       |    | ath:` | th:`\ |       |       |       |       | ath:` |     | h:`\l |       |
|       |    | \lamb | lambd |       |       |       |       | M_{\o |     | ambda |       |
|       |    | da_{m | a_{ef |       |       |       |       | dot}^ |     | _{0}` |       |
|       |    | ean}` | f}^{V |       |       |       |       | {AB}` |     |       |       |
|       |    |       | ega}` |       |       |       |       |       |     |       |       |
+-------+----+-------+-------+-------+-------+-------+-------+-------+-----+-------+-------+
| F300W | 1  | 0     | 0     | 0     | 1.398 | 99.99 | -2    | 7.433 | 0   | 0     | 1.000 |
|       |    | .2999 | .2993 | .0864 |       |       | 1.152 |       |     | .2999 |       |
+-------+----+-------+-------+-------+-------+-------+-------+-------+-----+-------+-------+
| F450W | 2  | 0     | 0     | 0     | -     | -     | -2    | 5.255 | 0   | 0     | 1.000 |
|       |    | .4573 | .4513 | .1077 | 0.074 | 0.339 | 0.609 |       |     | .4573 |       |
+-------+----+-------+-------+-------+-------+-------+-------+-------+-----+-------+-------+
| F606W | 3  | 0     | 0     | 0     | 0.095 | 0.161 | -2    | 4.720 | 0   | 0     | 1.000 |
|       |    | .6028 | .5827 | .2034 |       |       | 1.367 |       |     | .6028 |       |
+-------+----+-------+-------+-------+-------+-------+-------+-------+-----+-------+-------+
| F814W | 4  | 0     | 0     | 0     | 0.417 | 0.641 | -2    | 4.529 | 0   | 0     | 1.000 |
|       |    | .8013 | .7864 | .1373 |       |       | 2.322 |       |     | .8013 |       |
+-------+----+-------+-------+-------+-------+-------+-------+-------+-----+-------+-------+
| Jbb   | 5  | 1     | 1     | 0     | 0.890 | 99.99 | -2    | 4.559 | 0   | 1     | 1.000 |
|       |    | .2370 | .2212 | .2065 |       |       | 3.748 |       |     | .2370 |       |
+-------+----+-------+-------+-------+-------+-------+-------+-------+-----+-------+-------+
| H     | 6  | 1     | 1     | 0     | 1.361 | 99.99 | -2    | 4.702 | 0   | 1     | 1.000 |
|       |    | .6460 | .6252 | .3377 |       |       | 4.839 |       |     | .6460 |       |
+-------+----+-------+-------+-------+-------+-------+-------+-------+-----+-------+-------+
| K     | 7  | 2     | 2     | 0     | 1.881 | 99.99 | -2    | 5.178 | 0   | 2     | 1.000 |
|       |    | .2210 | .1971 | .3967 |       |       | 6.012 |       |     | .2210 |       |
+-------+----+-------+-------+-------+-------+-------+-------+-------+-----+-------+-------+

| 
| where :
| Col 1 : Name put in the first row of the filter file
| Col 2 : incremental number
| Col 3 : Mean wavelength (:math:`\mu m`) :
  :math:`\int R_{\lambda} \lambda d\lambda / \int R_{\lambda} d\lambda`
| Col 4 : Effective wavelength with Vega (:math:`\mu m`) :
  :math:`\int R_{\lambda} F_{\lambda}(Vega)\lambda d\lambda / \int R_{\lambda}F_{\lambda}(Vega) d\lambda`
| Col 5 : Full Width at Half of Maximum (:math:`\mu m`)
| Col 6 : AB Correction where :math:`m_{AB} = m_{VEGA} + ABcor`
| Col 7 : Thuan Gunn correction where :math:`m_{TG} = m_{VEGA} + TGcor`.
  (99.99 if undefined)
| Col 8 : VEGA magnitude :
  :math:`2.5\log_{10}(\int R_{\lambda} F_{\lambda}(Vega) d\lambda / \int R_{\lambda} d\lambda`)
| Col 9 : AB absolute magnitude of the sun (:math:`M^{AB}_{\nu,\odot}`)
   [4]_
| Col 10: value of the calibration used for
  (:math:`B_{\nu}/B_{\nu_0}`,\ :math:`\nu_0`) in **FILTER_CALIB**
| Col 11: Effective wavelength (:math:`\mu m`)
  :math:`\lambda_{0}^{B_{\nu}}= \frac{\int R_{\lambda} B_{\lambda} \lambda d\lambda}{\int R_{\lambda}  B_{\lambda}  d\lambda}`.
| Col 12: Correction factor to be applied to the original flux measured
  by LEPHARE. This correction is included in the programs **mag_gal**
  and **mag_star** as :math:`flux^{cor}= flux^{LePhare}\times`\ fac_cor

Extinction informations
~~~~~~~~~~~~~~~~~~~~~~~

| The stand alone program (**filter_extinc**) returns information about
  atmospheric extinctions and galactic extinctions.
| A set of atmospheric extinction curves and galactic extinction laws
  are available in $LEPHAREDIR/ext/ directory. It includes Calzetti and
  Prevot extinction laws. The Cardelli law is hardcoded in the programs
  and is the default law for the galactic extinction.
| % **filter_extinc** -c COSMOS.para -FILTER_FILE filter_test.dat
| It returns:
| #######################################
| # Computing ATMOSPHERIC AND GALACTIC EXTINCTION
| # with the following options:

=============================== =================
# Filters:                      filter_extinc.dat
# Atmospheric extinction curve: extinc_etc.dat
# Galactic extinction curve:    CARDELLI
# Output file:                  filter_extinc.dat
=============================== =================

| 
| #######################################

====================== ================ ======== ============
Filters                Ext(mag/airmass) Albda/Av Albda/E(B-V)
cosmos/u_cfht          0.486            1.504    4.663
cosmos/B_subaru        0.264            1.297    4.020
cosmos/V_subaru        0.141            1.006    3.118
cosmos/r_subaru        0.096            0.858    2.659
cosmos/i_subaru        0.052            0.643    1.992
cosmos/suprime_FDCCD_z 0.027            0.471    1.461
vista/Y                0.049            0.391    1.211
vista/J                0.096            0.281    0.871
vista/H                0.100            0.181    0.562
vista/K                0.100            0.118    0.364
====================== ================ ======== ============

| 
| Col 2 : Mean atmospheric extinction (mag/airmass) using (EXT_CURVE):
  :math:`A_{\lambda}= \int R_{\lambda} Ext(\lambda) d\lambda / \int R_{\lambda} d\lambda`
| :math:`Ext(\lambda)` comes from any atmospheric extinction curve that
  is put in $\ *LEPHAREDIR*/ext/.
| Col 3 : Mean galactic attenuation (in :math:`A(\lambda)/A_V`) using
  the galactic extinction law (GAL_CURVE). Col 4 : Mean galactic
  attenuation (in :math:`A(\lambda)//E(B-V)`) as a function of color
  excess (E(B-V)) assuming :math:`A_V=R_V\times E(B-V)`.
| For :math:`R_V` coefficients, we assume :math:`R_V=3.1` for most
  extinction laws but Calzetti (:math:`R_V=4.05`) and Prevost
  (:math:`R_V=2.72`).
| Others extinction laws can be added by following the format
  (:math:`\lambda(\AA) , k_{\lambda}`).

.. container:: float
   :name: fig:ext

.. _`sec:filtcalib`:

Application to long wavelengths 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

LEPHARE has been developped for the optical-NIR domain but can be used
at shorter (UV) and longer wavelengths (FIR, submm and radio). In
particular extensive tests have been performed in the long wavelength
domain by E. Le Floc’h to evaluate the photometric accuracy. Some issues
have to be considered :

-  the Vega spectrum is not defined at :math:`\lambda\ge 160\mu m`.
   Thus, AB magnitudes should be used as standard when combining a large
   wavelength domain.

-  The bandpass in radio domain is very narrow and does not require to
   convolve through the filter. However the structure of LEPHARE
   requires to implement a transmission curves for the radio frequencies
   in similar way as in shorter wavelengths.

More important, at long wavelengths the equivalent fluxes are taken as
the monochromatic flux density calculated at the effective wavelength of
the filter and for a reference spectum that would result in the same
energy received on the detector:

.. math:: <F_{\nu}> = \frac{\int F_{\nu} R_{\nu} d\nu}{\int \frac{B_{\nu}}{B_{\nu_0}} R_{\nu} d\nu}

where :math:`B_\nu` is the reference spectrum and :math:`\nu_0` the
effective frequency of the filter. In LEPHARE, the flux estimates are
equivalent to consider :math:`\frac{B_{\nu}}{B_{\nu_0}}=1`
(:math:`B_{\nu}=ctt`). Therefore there is a correction factor to account
for with respect to the original flux estimated by LEPHARE. This
correction is :

.. math:: <F_{\nu}>^{COR} = <F_{\nu}>^{LePhare} \times \frac{\int R_{\nu} d\nu}{\int \frac{B_{\nu}}{B_{\nu_0}} R_{\nu} d\nu}

| At long wavelengths, different conventions have been used for the
  reference spectrum. As an example: SPITZER/IRAC uses a flat spectrum
  (:math:`\nu B_{\nu}=ctt`) as well as ISO; SPITZER/MIPS uses a
  blackbody with temperature T=10000K while SCUBA uses planets which
  have SEDs in submillimeter very close to :math:`B_{\nu}=\nu`. The
  keyword FILTER_CALIB is used to account for these different
  calibration scheme (see section `3.3 <#sec:filter>`__).
| One additional effect is the way the effective wavelength is defined.
  In the case of MIPS, the effective wavelength seems to be defined,
  according to the MIPS handbook, as :math:`\nu B_{\nu}=ctt` while the
  reference spectrum is a black body. This mix definition can be
  described with **FILTER_CALIB=4**.
| In the table below we report the effective wavelengths and the
  correction factors that are applied to LEPHARE fluxes for a set of
  filters spanning from NIR (K band), MIR (SPITZER/IRAC), FIR
  (SPITZER/MIPS), sub-mm (SCUBA) to radio (VLA: 1.4GHz).

+-------+-------+-------+-----+-------+-------+-----+-------+-------+
| #NAME | :m    | :m    | CAL | :mat  | Fac   | CAL | :mat  | Fac   |
|       | ath:` | ath:` |     | h:`\l |       |     | h:`\l |       |
|       | \lamb | M_{\o |     | ambda |       |     | ambda |       |
|       | da_{m | dot}^ |     | _{0}^ |       |     | _{0}^ |       |
|       | ean}` | {AB}` |     | {B_{\ |       |     | {B_{\ |       |
|       |       |       |     | nu}}` |       |     | nu}}` |       |
+=======+=======+=======+=====+=======+=======+=====+=======+=======+
| K     | 2     | 5.178 | 0   | 2     | 1.000 | 0   | 2     | 1.000 |
|       | .2210 |       |     | .2210 |       |     | .2210 |       |
+-------+-------+-------+-----+-------+-------+-----+-------+-------+
| I     | 3     | 6.061 | 1   | 3     | 1.004 | 1   | 3     | 1.004 |
| RAC_1 | .5634 |       |     | .5504 |       |     | .5504 |       |
+-------+-------+-------+-----+-------+-------+-----+-------+-------+
| I     | 4     | 6.559 | 1   | 4     | 1.004 | 1   | 4     | 1.004 |
| RAC_2 | .5110 |       |     | .4930 |       |     | .4930 |       |
+-------+-------+-------+-----+-------+-------+-----+-------+-------+
| I     | 5     | 7.038 | 1   | 5     | 1.005 | 1   | 5     | 1.005 |
| RAC_3 | .7593 |       |     | .7308 |       |     | .7308 |       |
+-------+-------+-------+-----+-------+-------+-----+-------+-------+
| I     | 7     | 7.647 | 1   | 7     | 1.011 | 1   | 7     | 1.011 |
| RAC_4 | .9595 |       |     | .8723 |       |     | .8723 |       |
+-------+-------+-------+-----+-------+-------+-----+-------+-------+
| 24mic | 23    | 9.540 | 4   | 23    | 0.968 | 3   | 23    | 1.006 |
|       | .8437 |       |     | .6750 |       |     | .2129 |       |
+-------+-------+-------+-----+-------+-------+-----+-------+-------+
| 70mic | 72    | 1     | 4   | 71    | 0.932 | 3   | 68    | 1.013 |
|       | .5579 | 2.213 |     | .4211 |       |     | .4725 |       |
+-------+-------+-------+-----+-------+-------+-----+-------+-------+
| 1     | 156   | 1     | 4   | 155   | 0.966 | 3   | 152   | 1.007 |
| 60mic | .9636 | 3.998 |     | .8945 |       |     | .6311 |       |
+-------+-------+-------+-----+-------+-------+-----+-------+-------+
| 8     | 866   | nan   | 5   | 865   | 0.997 | 2   | 862   | 1.000 |
| 50mic | .7652 |       |     | .3377 |       |     | .4710 |       |
+-------+-------+-------+-----+-------+-------+-----+-------+-------+
| VLA_1 | 2     | nan   | 5   | 2     | 1.000 | 2   | 2     | 1.000 |
| .4GHz | 14300 |       |     | 14248 |       |     | 14145 |       |
|       |       |       |     | .3782 |       |     | .1645 |       |
+-------+-------+-------+-----+-------+-------+-----+-------+-------+

| 
| As can be seen from this table :
| :math:`\bullet` For K band, we use FILTER_CALIB=0, so no correcting
  factor is applied.
| :math:`\bullet` For IRAC bands , we adopt :math:`\nu B_{\nu}=ctt`
  (FILTER_CALIB=1). The correction factors are less than 1% and can be
  neglected.
| :math:`\bullet` For MIPS bands (24, 70, 160\ :math:`\mu m`), we adopt
  :math:`B_{\nu}=BB(T=10,000K)` and :math:`\lambda_0` defined as
  :math:`\nu B_ {\nu}=ctt` (FILTER_CALIB=4), which seems to better
  reflect the current MIPS calibration. In this case, correction factors
  between 3% to 7% are applied to the theoretical magnitudes estimated
  with **mag_gal** program. However, we also compare the correction
  factors when both :math:`\lambda_0` and :math:`B_{\nu}` refer to a
  black body at T=10,000K (FILTER_CALIB=3). In this case, the
  corrections become negligeable with :math:`\sim`\ 1%.
| :math:`\bullet` For sub-mm (SCUBA, 850\ :math:`\mu m`) and radio (VLA:
  1.4GHz) wavelengths, no correction is required
| As a general conclusion, the flux measured by LEPHARE appear accurate
  at a level of 1% with respect to most of the calibration scheme
  considered at long wavelength and thus no correction is required. A
  special warning for MIPS calibration, where depending on the
  calibration scheme, a correction up to 7%, may be applied.

Requirement to create a new filter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| Filters are ASCII files with the following format :
| In first row : #   SHORT_NAME_of_FILTER      ADD_COMMENTS
| In next rows : :math:`\lambda (\AA)` Transmission
| Wavelengths must be in increasing order. It is better to put the
  lowest and highest :math:`\lambda` with Transmission=0. The units of
  Transmission are not considered
| In the c++ version, the header, the transmission at 0 on the edges,
  and the transmission sorted in lambda are done internally if not
  prepared by the user.
| As an exemple : I create filter pippo.pb and put it in
  $LEPHAREDIR/filt/pippo.pb :

======= ================================
# PIPPO This is close to window function
5000    0
5001    1
5999    1
6000    0
======= ================================

| 

.. _`sec:mag_gal`:

Predicted magnitudes for galaxy/qso/stars libraries : **mag_gal**
-----------------------------------------------------------------

.. _description-and-outputs-1:

Description and outputs
^^^^^^^^^^^^^^^^^^^^^^^

| The **mag_gal** program predicts the magnitudes expected for
  GALAXY/QSO/STAR templates at various redshifts. It establishes the
  flux library which will be compared later to the data.
| For a set of filters given by **-FILTER_FILE** and an input SED
  library defined by **-GAL_LIB_IN**, the magnitudes are computed at
  different redshifts defined by **-Z_STEP**. Extinctions can be applied
  as specified by the three keywords (**-EXTINC_LAW, -MOD_EXTINC,
  -EB_V**). If evolving stellar population models are used, the
  cosmology (**-COSMOLOGY**) will allow to reject models older than the
  age of the universe. The magnitude in VEGA or AB (defined by
  **-MAGTYPE**) are saved in the binary file defined by **-GAL_LIB_OUT**
  in $LEPHAREWORK/lib_mag/ with an attached doc file.
| An output file (**-LIB_ASCII YES** ) is written to check the
  magnitudes, color tracks with redshift ....

.. _syntax-and-parameter-values-2:

Syntax and parameter values
^^^^^^^^^^^^^^^^^^^^^^^^^^^

| The usual syntax : :math:`\%` **mag_gal -t** G (or Q, or S) **-c**
  zphot.para
| The parameters values :
| (XXX means either GAL/QSO/STAR and are selected with **-t G** / **-t
  Q** ) / **-t S** )

+-----------------+------------------+----------+------------------+
| Parameters      | type             | default  | description      |
+=================+==================+==========+==================+
| FILTER_FILE(\*) | string           | —-       | Name of the      |
|                 |                  |          | filter file      |
+-----------------+------------------+----------+------------------+
|                 | (                |          | file must exist  |
|                 | n\ :math:`=`\ 1) |          | in               |
|                 |                  |          | $\ *LE           |
|                 |                  |          | PHAREWORK*/filt/ |
+-----------------+------------------+----------+------------------+
| XXX_LIB_IN(\*)  | string           | —-       | Name of the      |
|                 |                  |          | GALAXY/QSO/STAR  |
|                 |                  |          | binary library   |
|                 |                  |          | (with no         |
|                 |                  |          | extension)       |
+-----------------+------------------+----------+------------------+
|                 |                  |          | created by       |
|                 |                  |          | **sedtolib**;    |
+-----------------+------------------+----------+------------------+
|                 | (n=1)            |          | Files must exist |
|                 |                  |          | in               |
|                 |                  |          | $\ *LEPHA        |
|                 |                  |          | REWORK*/lib_bin/ |
+-----------------+------------------+----------+------------------+
| XXX_LIB_OUT(\*) | string           | —-       | Name of the      |
|                 |                  |          | magnitude binary |
|                 |                  |          | library (with no |
|                 |                  |          | extension)       |
+-----------------+------------------+----------+------------------+
|                 | (n=1)            |          | files            |
|                 |                  |          | *$GAL[Q          |
|                 |                  |          | SO]_LIB_OUT*.bin |
|                 |                  |          | (.doc)           |
+-----------------+------------------+----------+------------------+
|                 |                  |          | are saved in     |
|                 |                  |          | $\ *LEPHA        |
|                 |                  |          | REWORK*/lib_mag/ |
+-----------------+------------------+----------+------------------+
| MAGTYPE(\*)     | string           | —-       | Magnitude type   |
|                 |                  |          | (AB or VEGA)     |
+-----------------+------------------+----------+------------------+
|                 |                  |          |                  |
+-----------------+------------------+----------+------------------+
| ZGRID_TYPE      | int              | 0        | 0: constant step |
|                 |                  |          | in redshift      |
+-----------------+------------------+----------+------------------+
|                 | (n=1)            |          | 1: evolving step |
|                 |                  |          | in redshift as   |
|                 |                  |          | :math:`          |
|                 |                  |          | dz \times (1+z)` |
+-----------------+------------------+----------+------------------+
| Z_STEP          | float            | 0.04,0,6 | dz,zmin,zmax:    |
|                 |                  |          | redshift step    |
|                 |                  |          | (dz),            |
+-----------------+------------------+----------+------------------+
|                 | (n=3)            |          | the minimum      |
|                 |                  |          | (zmin) and the   |
|                 |                  |          | maximum redshift |
|                 |                  |          | (zmax).          |
+-----------------+------------------+----------+------------------+
| COSMOLOGY(\*)   | float            | —-       | :math:`H_0`,     |
|                 |                  |          | :                |
|                 |                  |          | math:`\Omega_0`, |
|                 |                  |          | :m               |
|                 |                  |          | ath:`\Lambda_0`. |
|                 |                  |          | Used for age     |
|                 |                  |          | constraints.     |
+-----------------+------------------+----------+------------------+
|                 | (n=3)            |          |                  |
+-----------------+------------------+----------+------------------+
| EXTINC_LAW      | string           | NONE     | Extinction laws  |
|                 |                  |          | to be used (in   |
|                 |                  |          | $\ *LEP          |
|                 |                  |          | HAREDIR*/ext/\*) |
+-----------------+------------------+----------+------------------+
|                 | (n\              |          | several files    |
|                 | :math:`\le`\ 10) |          | separated by     |
|                 |                  |          | comma            |
+-----------------+------------------+----------+------------------+
| MOD_EXTINC      | integer          | 0,0      | Range of models  |
|                 |                  |          | for which        |
|                 |                  |          | extinction will  |
|                 |                  |          | be applied       |
+-----------------+------------------+----------+------------------+
|                 | (n\              |          | The numbers      |
|                 | :math:`\le`\ 20) |          | refer to the     |
|                 |                  |          | models in the    |
|                 |                  |          | *$GAL_SED* list  |
+-----------------+------------------+----------+------------------+
|                 |                  |          | Number of values |
|                 |                  |          | must be twice    |
|                 |                  |          | the number of    |
|                 |                  |          | extinction laws. |
+-----------------+------------------+----------+------------------+
| EB_V            | float            | 0.       | Reddening color  |
|                 |                  |          | excess E(B-V)    |
|                 |                  |          | values to be     |
|                 |                  |          | applied          |
+-----------------+------------------+----------+------------------+
|                 | (n\ :            |          | values separated |
|                 | math:`\le`\ 100) |          | by comma.        |
+-----------------+------------------+----------+------------------+
| EM_LINES        | string           | NO       | Add contribution |
|                 |                  |          | of emission      |
|                 |                  |          | lines and        |
|                 |                  |          | specify          |
+-----------------+------------------+----------+------------------+
|                 | (n=1)            |          | how to derive    |
|                 |                  |          | them             |
|                 |                  |          | (``EMP_UV``,     |
|                 |                  |          | ``EMP_SFR``,     |
|                 |                  |          | ``PHYS``)        |
+-----------------+------------------+----------+------------------+
| EM_DISPERSION   | float            | 1        | the emission     |
|                 |                  |          | lines can vary   |
|                 |                  |          | by these         |
|                 |                  |          | fractions from   |
|                 |                  |          | the expected     |
+-----------------+------------------+----------+------------------+
|                 |                  |          | value (example   |
|                 |                  |          | 0.5,1.,1.5)      |
+-----------------+------------------+----------+------------------+
| ADD_DUSTEM      | string           | NO       | Add the dust     |
|                 |                  |          | emission in      |
|                 |                  |          | templates when   |
|                 |                  |          | missing.         |
+-----------------+------------------+----------+------------------+
|                 |                  | (n=1)    | This is based on |
|                 |                  |          | the energy       |
|                 |                  |          | absorbed over    |
|                 |                  |          | the UV-optical   |
|                 |                  |          | range.           |
+-----------------+------------------+----------+------------------+
| LIB_ASCII       | string           | NO       | ASCII file with  |
|                 |                  |          | magnitudes saved |
|                 |                  |          | in               |
|                 |                  |          | $\ *LEPHAREWORK* |
+-----------------+------------------+----------+------------------+
|                 | (n=1)            |          | called           |
|                 |                  |          | *$GAL[Q          |
|                 |                  |          | SO]_LIB_OUT*.dat |
+-----------------+------------------+----------+------------------+

The extinction laws and dust emission
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| A set of extinction laws are available in the directory
  (``$LEPHAREDIR/ext/``). Several extinction laws can be used and set up
  in the keyword **-EXTINC_LAW**. Each extinction law will be applied to
  a range of SED models specified by the keywords **-MOD_EXTINC**. The
  model number corresponds to the rank in the list of SEDs used in
  **-GAL_SED**. The number of models must be twice the number of
  extinction laws. The different values of reddening excess E(B-V) are
  given in the keyword **-EB_V** and will apply to all extinction laws.
  The extinguished flux is :
  :math:`F_{\lambda}^e = F_{\lambda}^0\  10^{-0.4 A_{\lambda}}=  F_{\lambda}^0\  10^{-0.4 k_{\lambda} E(B-V)}`
| If extinction is applied, a new estimate of the IR dust luminosity is
  computed by measuring the amount of light absorbed. Some templates
  don’t include dust emission. We add the possibility of having the dust
  emission by using ADD_DUSTEM YES. In such case, we use the templates
  from Bethermin et al. (2012) and sum their flux contribution to the
  stellar template (e.g. BC03). **Don’t use this option if your
  templates already include dust emission**. The B12 templates are
  different for each redshift. However, a current limitation of the code
  is that an incorrect dust SED is displayed in the .spec file (while
  the fit is correct). Therefore, we use by default only one B12
  template at :math:`z=0`. The fit will be correct if you use all
  templates (but not the final display).

The Emission lines
^^^^^^^^^^^^^^^^^^

The role of nebular emission lines in medium- and even broad-band
filters has been shown to be essential in several cases (Ilbert et 2009,
Schearer et al. 2009, Labbe et al. 2013, Stefanon et al. 2015). Some
templates already include emission lines. In this case, you could use
**EM_LINES NO** to avoid creating additional ones. To include emission
lines in the template SEDs if they don’t exist, one of the available
methods must be selected through the parameter **EM_LINES**. There are
three different options:

-  **EMP_UV** LePHARE accounts for the contribution of emission lines
   with a simple recipe based on the Kennicutt (1998) relations. The SFR
   is estimated from UV luminosity, which in turn defines the
   H\ :math:`\alpha` luminosity. Intensity of other lines
   (:math:`Ly_{\alpha}`, :math:`H_{\alpha}`, :math:`H_{\beta}`, [OII],
   OIII[4959] and OIII[5007]) are defined accordingly by using the flux
   ratios provided in Ilbert et al. (2009) and slightly adjusted since.
   The UV luminosity is derived directly from the SED template. Emission
   lines are not considered in red galaxies with
   :math:`(NUV-r)_{ABS}\ge 4` (rest frame, dust corrected color). This
   option works for any kind of input template.

-  **EMP_SFR** At present, this option can be used only with BC03
   templates. This option can be used with SED templates that have SFR
   already defined (BC03). The SFR is converted in H\ :math:`\alpha`
   according to Kennicutt (1998). It skips the conversion from UV to SFR
   done with the option **EMP_UV**.

-  **PHYS** At present, this option can be used only with BC03
   templates. For each of them, LePhare reads metallicity, fraction of
   photoionizing photons, and other physical quantities needed as input
   in a model (Schearer et al. 2009) that quantifies flux emitted by
   several emission lines. To see details and applications of this
   method in Shun et al. (2019, in prep).

In all the methods, dust attenuation is applied to the emission line
according the continuum value. The MW (Seaton 1979) extinction curve is
considered for the emission lines. A factor :math:`f` is introduce
between the E(B-V) obtained for the stellar content and the E(B-V)
considered for the emission lines. This value is taken as 1.

With the option **EM_DISPERSION**, the emission lines can vary from the
standard value; for example by setting the option to -EM_DISPERSION
0.5,1.,1.5 the code generates three SEDs with identical characteristics,
except the lines will have the standard flux (prescribed by the EMP\_ or
PHY recipe) and :math:`\pm50\%` of that value.

Even if emission lines have been built for the entire library, during
any SED fitting run the user can decide to ignore them for a given
subset of models (see ADD_EMLINES option in Section `5.4 <#fit>`__).

This option is not appropriated for the quasars samples.

ASCII ouput file
^^^^^^^^^^^^^^^^

| An output file is produces in the current directory if **-LIB_ASCII
  YES**. It has the same root name as the binary file with extension
  .dat and contains the following informations :
| Model, Extinc-law, E(B-V), :math:`L_{TIR}(L_{\odot})`, Z, DMod,
  Age(yr), nrec, n , (mag(i),i=1,n),(kcor(i),i=1,n)
| where Model is the number of models based on the original list,
  Extinc-law refers to the number of the extinction laws used,
  :math:`L_{TIR}` the new estimate of the IR luminosity, DMod is the
  distance modulus, nrec is a record (internal use), n the number of
  filters, mag(i) the magnitudes in all filters and kcor(i), the k
  correction in all filters.

Sizing the library
^^^^^^^^^^^^^^^^^^

| You must be aware that the size of the library becomes quickly huge if
  you do not pay attention. You can estimate its size by considering the
  following numbers :
| # of models x # of age x # of z steps x # of extinction law x # of
  EB-V
| For exemple, 10 SEDs with 60 ages, 2 extinction laws and 6 E(B-V) and
  150 z steps will exceed 1,000,000 rows.

.. _example-1:

 Example
^^^^^^^^

| **mag_gal -t Q -c** zphot.para **-FILTER_FILE HDF.filt**
  **-QSO_LIB_IN** LIB_QSO **-QSO_LIB_OUT** QSO_HDF **-EXTINC_LAW** NONE
| It will generate the magnitudes for the QSO library (LIB_QSO.bin)
  through the filters HDF.filt. Those two files have been created during
  the two previous steps. No extinction will be applied. The output
  QSO_HDF.bin and QSO_HDF.doc are saved in $LEPHAREWORK/lib_mag/
| zphot.para **-FILTER_FILE HDF.filt** **-GAL_LIB_IN** LIB_CWW
  **-GAL_LIB_OUT** CWW_HDF **-EXTINC_LAW**
  SMC_prevot.dat,SB_calzetti.dat **-MOD_EXTINC** 3,6,4,8 **-EB_V**
  0.,0.05,0.1,0.2,0.3 **-LIB_ASCII** YES
| It will generate the magnitudes for the galaxy library (LIB_CWW) and
  HDF.filt. The library LIB_CWW is built with the following option in
  sedtolib:
| **sedtolib -t G -c** zphot.para **-GAL_SED**
  $LEPHAREDIR/sed/GAL/CWW_KINNEY/CWW_MOD.list **-GAL_LIB** LIB_CWW
| CWW_MOD.list contains the following SEDs : 1:Ell, 2:Sbc, 3:Scd, 4:Im,
  5:SB1, 6:SB2, 7:SB3, 8:SB4.
| The two extinction laws are applied as follows :
| :math:`\bullet` SMC_prevost is used for models Scd (3), Im(4), SB1(5),
  SB2(6)
| :math:`\bullet` SB_calzetti is used for models Im(4), SB1(5), SB2(6),
  SB3(7), SB4 (8)
| The overlapping models (Im, SB1 and SB2) will be extinguished with the
  2 extinction laws.
| For both extinctions, the same values of E(B-V) are used.
| The files CWW_HDF.bin and CWW_HDF.doc are saved in
  $LEPHAREWORK/lib_mag/ and the ASCII file CWW_HDF.dat is written in the
  current directory.

The photometric redshift program: ``zphota``
--------------------------------------------

| The program ``zphota`` performs a :math:`\chi^2`-based analysis,
  fitting the predicted flux (built in Sect. `4 <#sec:mag_gal>`__) to
  the observed photometry (AB/Vega magnitudes or fluxes). To measure the
  photometric redshift we use a :math:`\chi^2` fitting procedure by
  comparing the observed flux (:math:`F_{obs}`) and its corresponding
  uncertainties (:math:`\sigma`) with the flux from templates
  (:math:`F_{temp}`) defined as:

  .. math:: \chi^2 =   \sum_i [ \frac{F_{obs,i} - s F_{temp,i}}{\sigma_i}]^2

  where i refers to the band used for the analysis and :math:`s` the
  scaling factor that is chosen to minimize the :math:`\chi^2` values
  (:math:`{\it d}\chi^2/{\it d}s=0`):

  .. math:: s =   \sum_j [ \frac{F_{obs,j}  F_{temp,j}}{\sigma_j^2} ]  / \sum_j [ \frac{F_{temp,j}^2}{ \sigma_j^2}]

  where j refers to the band used for the scaling (j can be different
  from i).
| The photometric baseline can span a large wavelength range, as long as
  the templates are established accordingly. Galaxy, star, and QSO
  libraries can be used in the same run, but the :math:`\chi^2`
  minimization process is performed distinctly for each class. For a
  given class (e.g., galaxy SEDs) several libraries can be combined.
| Different options are available to improve the :math:`z_\mathrm{phot}`
  measurement: physical priors, adaptive photometric adjustments,
  addition of nebular emission lines in the synthetic SEDs. If the
  templates include physical information (Sect. `2 <#models>`__, e.g.
  BC03), ``zphota`` gives in output stellar mass, star formation rate,
  etc., for each object.
| As the previous commands, the basic syntax of this program is:
| :math:`\$` zphota -c zphota.para [-PARAM1 VALUE -PARAM2 VALUE ...]
| assuming that zphota.para is the name of the configuration file.

.. _lib:

Input libraries
^^^^^^^^^^^^^^^

| The principle of SED-fitting is to compare observed flux with
  predicted ones. We can extract from this comparison the photometric
  redshift but also physical parameters associated to the galaxies.
  Therefore, a fundamental input of ``zphota`` is a library containing
  predicted flux created with ``mag_gal``.
| The name of this library should be transmitted to ``zphota`` using the
  keyword **ZPHOTLIB**. The name should be a string and points to the
  binary file stored in ``$LEPHAREWORK/lib_mag/``. Indicate only the
  name of the file without extension. For instance, if a file
  ``BC03_LIB.bin`` has been created by ``mag_gal`` and is stored in
  ``$LEPHAREWORK/lib_mag/``, you can simply use the option
  ``-ZPHOTLIB BC03_LIB``.
| Several librairies can be combined, with their name separated with
  coma. You can use as many libraries as you want. Moreover, you can
  combine libraries created with GALAXY/QSO/STAR templates and the code
  will recognize if it corresponds to a GALAXY, QSO, or STAR library.
| Finally, one can modify the properties of the input library by
  applying emission lines to only a sub-sample of the templates and by
  reducing the explored range of E(B-V) and redshift. For instance
  ``ADD_EMLINES`` defines the range of galaxy models (from the .list
  file) in which the code considers the emission lines contribution.
  Similarly ``Z_RANGE`` and ``EBV_RANGE`` could be used to reduce the
  redshift and the E(B-V) coverage allowed in the fit.

.. _input:

Input file
^^^^^^^^^^

| This section describes how to manage the input file.
| **CAT_IN** specifies the location and name of the input file.
| The input catalogue must be an ASCII table including at least for each
  entry:

-  an identification number (Id);

-  the apparent magnitudes (or fluxes);

-  the corresponding errors.

| The format is specified by **CAT_FMT**, whose value must be set to
  **MEME** (“Magnitude-Error-Magnitude-Error”) to use a catalog in the
  format
| *Id mag1 err1 mag2 err2 ... magN errN*...
| while the string **MMEE** (“Magnitude...Magnitude-Error...Error”) is
  used for catalogs written like
| *Id mag1 mag2 ... magN err1 err2 ... errN*...
| Other columns may follow the photometric baseline when the option
  **CAT_TYPE** is set to **LONG** (it is **SHORT** by default). Such
  extended catalog will look like:
| ``Id mag1 err1 mag2 err2 ... magN errN Context``\ :math:`z_\mathrm{spec}`\ ``Extra1 Extra2...``
| The ``Context`` indicates which passbands can be used for the object
  in this row (see below), :math:`z_\mathrm{spec}` is the input redshift
  (can be also equal to -99), and “Extra1”, “Extra2”, etc. are the
  remaining columns (any kind of values) that will be read by the
  program as a single string and propagated in the output if required.
  Only ``Context`` and :math:`z_\mathrm{spec}` are compulsory in the
  LONG format, while Extra1, Extra2, etc. can be left empty.
| The input catalogue could include magnitudes or fluxes. To use fluxes,
  you must specify **F** for the parameter **INP_TYPE** and fluxes must
  be given in :math:`\mathrm{erg}/\mathrm{s}/\mathrm{cm}^2/\mathrm{Hz}`.
  If you use magnitude in input, use ``INP_TYPE M``. In this case, The
  calibration system is declared by the parameter **CAT_MAG**, which can
  be either **VEGA** or **AB**. In any case the filters in the catalog
  must be the same (and in the same order) as in the SED library built
  with ``mag_gal``.
| For a given object, the flux in a given filter could miss (not
  observed or the photometric extraction failed). If the magnitude (or
  flux) and the associated are **both** negative, this filter will be
  ignored. If the measurement is missing because the flux is too faint
  to be detected, one could use an upper limit. In such case, the
  magnitude (or flux) are positive and set to the upper-limit value
  while the error should be negative.
| You can run ``zphota`` on a subsample of sources. **CAT_LINE** gives
  the range of entries which should be considered when running the code.
  For instance, ``CAT_LINE 1,1000`` will run the code only on the first
  1000 lines.
| NOTE: commented lines are NOT considered while reading the catalogue,
  so this range should be intended as the number of entries, not rows.

Context
~~~~~~~

| The Context is an integer value which specifies the filter combination
  to be used. It is defined as the sum of powers of 2 :
  Cont\ :math:`=\sum_{i=1}^{i=N} 2^{i-1}`, where i is the filter number
  as ordered in the input catalog (and in the library), and N is the
  total number of filters.
| As an example, let’s consider a catalog with the following passbands:

================================== = = = = == == == ===
Passband                           U G R I Z  J  H  K
Filter number (i)                  1 2 3 4 5  6  7  8
Filter Context (:math:`2^{(i-1)}`) 1 2 4 8 16 32 64 128
================================== = = = = == == == ===

| 
| :math:`\bullet` If the context is included in the catalog (CAT_TYPE=
  LONG), you can specify a context for each object. One context value
  corresponds to a unique filter combination:
| if an object is observed in all passband but H : Context=191
| if an object is observed in UGRIZ : Context=31
| if an object is observed in GRIZK : Context=158
| :math:`\bullet` If the context is absent in the input catalog
  (CAT_TYPE =SHORT), it is equivalent to use all the passbands for all
  the objects, so Context=255. However, the code checks the error and
  flux values. If both values are negative, the band is not used.
| In practice, the context specified in the input catalog can include
  all the passbands where the object has been observed even the bands
  where it is not detected (upper-limit). Additional options in the
  configuration file will allow to restrict the use of the catalog to
  some specific filter combinations.
| Note 1: if the flux (or mag) and the associated error are negative,
  the filter is ignored in the fit.
| Note 2: In the configuration file, some options refer to a sum of
  filter context:
| GLB_CONTEXT, FORB_CONTEXT, ADAPT_CONTEXT, MABS_CONTEXT, FIR_CONT,
  FIR_SCALE

.. _output:

Output files
^^^^^^^^^^^^

| The name of the output file is given with the ``CAT_OUT`` keyword.
| The format of the output file is flexible. All the columns that the
  user want in output are listed in a parameter file. The name of the
  parameter file should be given in ``PARA_OUT``.
| An example of parameter file with all the existing columns is given in
  ``config/output.para``. Some keywords can be removed or commented with
  #. You can also modify the order of the keywords. The symbol ()
  indicates a vector (with a dimension corresponding to the number of
  filters).
| You can also decide to get the redshift PDF for each source stored in
  a single ascii file. You need to fill the keyword ``PDZ_OUT`` with the
  name of the output file. Don’t put any extension, the code will add it
  for you. You will get the probability measured at each redshift step
  listed in the header. The file will contain one line per object.
| If ``SPEC_OUT YES``, an output file is created for each object. This
  file contains several information on the considered object (like the
  observed magnitudes, the spec-z, the photo-z, etc), but also the PDF
  and the best-fit templates. These files will be named as
  ``IdXXXX.spec`` with XXXX being the Id of the source. The file can be
  read using a python code ``spec.py`` located in the directory
  ``$LEPHAREDIR/tools`` (or the sm macro ``spec.sm`` if you prefer). You
  can create a file containing the figures for several sources using:
| ``python spec.py Id00000*.spec -d pdf``
| It will create a file ``MULTISPEC.pdf`` with all the fit.
| If the user put a name different from YES/NO as argument of SPEC_OUT,
  it will be used as directory to store the .spec files.
| You can also decide to get the full :math:`\chi^2` map (the value of
  the fit for each redshift, template, E(B-V), etc). Be careful that it
  could take a lot of disk space. It could be useful if you have one
  source that you want to study in detail.

.. table:: List of parameters to configure *LePHARE++* . First column is
the keyword to set the parameter, which can be set in the configuration
file or in the command line. The second column is the type of the given
parameter (string, integer, or float) with dimension in square bracket.
For parameters with size :math:`>1` values must be comma-separated
(e.g., :math:`1,2,3`). For parameters having a default value, this is
listed in the third column (a hyphen, —, is shown otherwise). The fourth
column gives a short description of the parameter. Keywords with (\*)
must be defined, all the other keywords are optional.

   +----------------+----------------+----------------+----------------+
   | **Parameters** | **Type**       | **Default      | *              |
   |                |                | val.**         | *Description** |
   +================+================+================+================+
   | LePhare++ SED  |                |                |                |
   | libraries:     |                |                |                |
   +----------------+----------------+----------------+----------------+
   | ZPHOTLIB(\*)   | string         | —-             | Library names  |
   |                |                |                | (with no       |
   |                |                |                | extension)     |
   |                |                |                | like           |
   |                |                |                | XXX_LIB_OUT    |
   +----------------+----------------+----------------+----------------+
   |                | (:ma           |                | Files should   |
   |                | th:`n \geq 1`) |                | exist in       |
   |                |                |                | $\ *LEPHARE    |
   |                |                |                | WORK*/lib_mag/ |
   +----------------+----------------+----------------+----------------+
   | ADD_EMLINES    | int            | 0,0            | Range of       |
   |                |                |                | galaxy models  |
   |                |                |                | in which       |
   +----------------+----------------+----------------+----------------+
   |                | (n\            |                | considering    |
   |                |  :math:`=`\ 2) |                | emission lines |
   |                |                |                | contribution.  |
   +----------------+----------------+----------------+----------------+
   | Z_RANGE        | float          | 0.,99.         | Z min and max  |
   |                |                |                | allowed in the |
   |                |                |                | GALAXY library |
   +----------------+----------------+----------------+----------------+
   |                | (n=2)          |                |                |
   +----------------+----------------+----------------+----------------+
   | EBV_RANGE      | float          | 0,9            | E(B-V) min and |
   |                |                |                | max allowed in |
   |                |                |                | the GALAXY     |
   |                |                |                | library        |
   +----------------+----------------+----------------+----------------+
   |                | (n=2)          |                |                |
   +----------------+----------------+----------------+----------------+
   |                |                |                |                |
   +----------------+----------------+----------------+----------------+
   | Input catalog  |                |                |                |
   | (Sect. `5.2    |                |                |                |
   |  <#input>`__): |                |                |                |
   +----------------+----------------+----------------+----------------+
   |                |                |                |                |
   +----------------+----------------+----------------+----------------+
   | CAT_IN(\*)     | string[1]      | —-             | Name of the    |
   |                |                |                | input          |
   |                |                |                | photometric    |
   |                |                |                | catalogue      |
   |                |                |                | (full path)    |
   +----------------+----------------+----------------+----------------+
   | INP_TYPE(\*)   | string[1]      | —-             | Input values:  |
   |                |                |                | Flux (F) or    |
   |                |                |                | Magnitude (M); |
   |                |                |                | see            |
   |                |                |                | Sect. `5       |
   |                |                |                | .2 <#input>`__ |
   |                |                |                | for units.     |
   +----------------+----------------+----------------+----------------+
   | CAT_MAG(\*)    | string[1]      | —-             | Input          |
   |                |                |                | magnitude type |
   |                |                |                | : AB or VEGA   |
   +----------------+----------------+----------------+----------------+
   | CAT_FMT(\*)    | string[1]      | —-             | Input format   |
   |                |                |                | for photometry |
   |                |                |                | (MEME or MMEE, |
   |                |                |                | see            |
   |                |                |                | Sect. `5.      |
   |                |                |                | 2 <#input>`__) |
   +----------------+----------------+----------------+----------------+
   | CAT_LINES      | integer[2]     | -99,-99        | Min and max    |
   |                |                |                | rows read in   |
   |                |                |                | input catalog  |
   |                |                |                | (starting from |
   |                |                |                | 1)             |
   +----------------+----------------+----------------+----------------+
   | CAT_TYPE       | string[1]      | SHORT          | Input catalog  |
   |                |                |                | format (see    |
   |                |                |                | Sect. `5.2.1   |
   |                |                |                | <#context>`__) |
   +----------------+----------------+----------------+----------------+
   |                |                |                |                |
   +----------------+----------------+----------------+----------------+
   | Output         |                |                |                |
   | catalog:       |                |                |                |
   +----------------+----------------+----------------+----------------+
   | CAT_OUT        | string         | zphot.out      | Name of the    |
   |                |                |                | output file    |
   |                |                |                | (full path)    |
   +----------------+----------------+----------------+----------------+
   |                | (n\            |                | by default     |
   |                |  :math:`=`\ 1) |                | saved in       |
   |                |                |                | working        |
   |                |                |                | directory      |
   +----------------+----------------+----------------+----------------+
   | PARA_OUT(\*)   | string         | —-             | Name of the    |
   |                |                |                | file with      |
   |                |                |                | selected       |
   |                |                |                | output         |
   |                |                |                | parameters     |
   |                |                |                | (full path)    |
   +----------------+----------------+----------------+----------------+
   |                | (n\            |                |                |
   |                |  :math:`=`\ 1) |                |                |
   +----------------+----------------+----------------+----------------+
   | SPEC_OUT       | string         | NO             | Output files   |
   |                |                |                | with           |
   |                |                |                | Gal/Star/QSO   |
   |                |                |                | spectra (one   |
   |                |                |                | file per       |
   |                |                |                | object)        |
   +----------------+----------------+----------------+----------------+
   |                | (n=1)          |                | (if YES: can   |
   |                |                |                | take a lot of  |
   |                |                |                | disk space !)  |
   +----------------+----------------+----------------+----------------+
   |                |                |                | If a string    |
   |                |                |                | different from |
   |                |                |                | NO, save files |
   |                |                |                | in this        |
   |                |                |                | directory.     |
   +----------------+----------------+----------------+----------------+
   | CHI2_OUT       | string         | NO             | Output files   |
   |                |                |                | with all       |
   |                |                |                | :math:`\chi^2` |
   |                |                |                | for galaxy     |
   |                |                |                | library (one   |
   |                |                |                | file per       |
   |                |                |                | object)        |
   +----------------+----------------+----------------+----------------+
   |                | (n=1)          |                | (if YES: can   |
   |                |                |                | take a lot of  |
   |                |                |                | disk space !)  |
   +----------------+----------------+----------------+----------------+

.. _fit:

Managing filters used in the fit
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| The choice of the filters is defined by the context value for each
  object (see `5.2.1 <#context>`__). This context is given in the input
  catalogue. You can also force the analysis to some specific filter
  combination for the whole catalog. If **GLB_CONTEXT** is used, it
  supersedes the individual context. You can also reject some bands with
  **FORB_CONTEXT** keyword. This keyword is useful if you want to
  perform some test without a specific band.
| The empirical and stellar population synthesis libraries of galaxy
  SEDs only account for the stellar light. It is strongly suggested to
  only use filters where the stellar light is dominant. Typically we
  suggest to authorize only the filters with :math:`\lambda\le 5\mu m`.
  Longer wavelength information should be treated separately with the
  FIR libraries.

+----------------+----------------+----------------+----------------+
| **Parameters** | **Type**       | **Default      | *              |
|                |                | val.**         | *Description** |
+================+================+================+================+
| PASSBAND       |                |                |                |
| SELECTION      |                |                |                |
+----------------+----------------+----------------+----------------+
| GLB_CONTEXT    | integer        | -1             | Forces the     |
|                |                |                | context of all |
|                |                |                | objects for    |
|                |                |                | :math:`\chi^2` |
|                |                |                | analysis       |
+----------------+----------------+----------------+----------------+
|                | (n\            |                | defined as :   |
|                |  :math:`=`\ 1) |                | :mat           |
|                |                |                | h:`\sum_{i=0}^ |
|                |                |                | {nbd-1} 2^{i}` |
+----------------+----------------+----------------+----------------+
|                |                |                | 0 means that   |
|                |                |                | all bands are  |
|                |                |                | used           |
+----------------+----------------+----------------+----------------+
|                |                |                | -1 (default)   |
|                |                |                | means that     |
|                |                |                | context per    |
|                |                |                | object is used |
+----------------+----------------+----------------+----------------+
| FORB_CONTEXT   | integer        | -1             | context for    |
|                |                |                | forbidden      |
|                |                |                | bands          |
+----------------+----------------+----------------+----------------+
|                | (n\            |                | defined as :   |
|                |  :math:`=`\ 1) |                | :mat           |
|                |                |                | h:`\sum_{i=0}^ |
|                |                |                | {nbd-1} 2^{i}` |
+----------------+----------------+----------------+----------------+
| RM             | float          | 200            | Threshold in   |
| _DISCREPENT_BD |                |                | chi2 to        |
|                |                |                | consider.      |
+----------------+----------------+----------------+----------------+
|                |                | (n\            | Remove 2 bands |
|                |                |  :math:`=`\ 1) | max, stop when |
|                |                |                | below this     |
|                |                |                | chi2           |
|                |                |                | threshold.     |
+----------------+----------------+----------------+----------------+
| INCREASING     |                |                |                |
| PHOTOMETRIC    |                |                |                |
| ERRORS         |                |                |                |
+----------------+----------------+----------------+----------------+
| ERR_FACTOR     | float          | 1.0            | Scaling factor |
|                |                |                | to the errors  |
|                |                |                | (in flux)      |
+----------------+----------------+----------------+----------------+
|                | (n\            |                |                |
|                |  :math:`=`\ 1) |                |                |
+----------------+----------------+----------------+----------------+
| ERR_SCALE      | float          | -1.            | Systematic     |
|                |                |                | errors (in     |
|                |                |                | mag) add in    |
|                |                |                | quadrature to  |
|                |                |                | the            |
|                |                |                | observations   |
+----------------+----------------+----------------+----------------+
|                | (n\ :ma        |                | must match     |
|                | th:`\le`\ 100) |                | number of      |
|                |                |                | bands, not     |
|                |                |                | used otherwise |
+----------------+----------------+----------------+----------------+
| ANALYSIS OF    |                |                |                |
| THE            |                |                |                |
| :math:`PDF(z)` |                |                |                |
+----------------+----------------+----------------+----------------+
| Z_INTERP       | string         | NO             | Parabolic      |
|                |                |                | interpolation  |
|                |                |                | between        |
|                |                |                | original step  |
|                |                |                | (dz)           |
+----------------+----------------+----------------+----------------+
|                | (n=1)          |                |                |
+----------------+----------------+----------------+----------------+
| DZ_WIN         | float          | 0.25           | “smoothing”    |
|                |                |                | window         |
|                |                |                | function for   |
|                |                |                | 2nd peak       |
|                |                |                | search in L(z) |
+----------------+----------------+----------------+----------------+
|                | (n=1)          |                | (value between |
|                |                |                | 0 to zmax)     |
+----------------+----------------+----------------+----------------+
| MIN_THRES      | float          | 0.1            | threshold for  |
|                |                |                | the detection  |
|                |                |                | of 2nd peak in |
|                |                |                | normalised     |
|                |                |                | L(z)           |
+----------------+----------------+----------------+----------------+
|                | (n=1)          |                | (value between |
|                |                |                | 0 to 1)        |
+----------------+----------------+----------------+----------------+
| PROB_INTZ      | float          | 0.             | redshift       |
|                |                |                | intervalles to |
|                |                |                | compute        |
|                |                |                | probability    |
|                |                |                | from F(z)      |
+----------------+----------------+----------------+----------------+
|                | (              |                | (even number   |
|                | n\ :math:`\le` |                | of values),    |
|                | 100)           |                | output vectors |
|                |                |                | from 0 to 100% |
+----------------+----------------+----------------+----------------+
|                |                |                | 0.-default :   |
|                |                |                | not used       |
+----------------+----------------+----------------+----------------+

| 

Expanding photometric uncertainties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| By definition the :math:`\chi^2` procedure is sensitive to the
  photometric errors, so it is important to provide reliable
  uncertainties. Users must account for a possible underestimation (when
  noise correlation is present in the data) or zero-point calibration
  uncertainties. The keywords ERR_FACTOR and ERR_SCALE allow to tune the
  individual errors. Note that ERR_FACTOR will not change the best
  photo-z solution but just the estimates of the errors, while ERR_SCALE
  can change the relative contribution of the bands and thus the best
  redshift.

Adding prior information
^^^^^^^^^^^^^^^^^^^^^^^^

| Additional constraints can be applied to the :math:`\chi^2` fitting
  procedure with the options below.
| *LePHARE++* could apply a prior on the redshift distribution,
  following a similar procedure than Benitez et al. (2000). This is done
  using the keyword **NZ_PRIOR**. We used the N(z) prior by type
  computed from the VVDS survey in I band and detailed in Ilbert et al.
  (2006).
| A prior could be applied to avoid unrealistically bright galaxies. The
  keyword ``MAG_ABS`` gives the absolute magnitude range allowed in a
  given filter **MAG_REF**. This could be defined by checking the
  luminosity function of the considered population. For field galaxies,
  a common range is -24,8 in the g-band.

+----------------+----------+------------------+------------------+
| **Parameters** | **Type** | **Default val.** | **Description**  |
+================+==========+==================+==================+
|                |          |                  | PRIOR KEYWORDS   |
+----------------+----------+------------------+------------------+
| NZ_PRIOR       | integer  | -1,-1            | N(z) prior as    |
|                |          |                  | function of I    |
|                |          |                  | band.            |
+----------------+----------+------------------+------------------+
|                | (n=2)    |                  | The i-band       |
|                |          |                  | number should be |
|                |          |                  | given in input.  |
+----------------+----------+------------------+------------------+
|                |          |                  | The second       |
|                |          |                  | number indicates |
|                |          |                  | which band to    |
|                |          |                  | use if first     |
|                |          |                  | undefined.       |
+----------------+----------+------------------+------------------+
|                |          |                  | Negative value   |
|                |          |                  | means no prior.  |
+----------------+----------+------------------+------------------+
| MAG_ABS        | float    | 0.,0.            | Absolute         |
|                |          |                  | magnitude range  |
|                |          |                  | acceptable for   |
|                |          |                  | GAL library      |
|                |          |                  | [0,0-def]        |
+----------------+----------+------------------+------------------+
|                | (n=2)    |                  | 0.,0. (default)  |
|                |          |                  | means not used   |
+----------------+----------+------------------+------------------+
| MAG_ABS_QSO    | float    | 0.,0.            | Absolute         |
|                |          |                  | magnitude range  |
|                |          |                  | acceptable for   |
|                |          |                  | QSO library      |
|                |          |                  | [0,0-def]        |
+----------------+----------+------------------+------------------+
|                | (n=2)    |                  | 0.,0. (default)  |
|                |          |                  | means not used   |
+----------------+----------+------------------+------------------+
| MAG_REF        | integer  | 0                | Reference filter |
|                |          |                  | for MAG_ABS (1   |
|                |          |                  | to               |
|                |          |                  | :math:`N_{bd}`)  |
+----------------+----------+------------------+------------------+
|                | (n=1)    |                  | 0 (default)      |
|                |          |                  | means not used   |
+----------------+----------+------------------+------------------+

Adaptive method
^^^^^^^^^^^^^^^

| In the c++ version, we provide the possibility to train the
  zero-points of the photometric catalogue. While this training is less
  sophisticated than the fortran version (which allows for a training of
  the colors and more), this training is sufficient for most of the
  applications.
| In order to turn on this option, use **AUTO_ADAPT YES**. This
  procedure requires to have galaxies with a spec-z within the catalogue
  (format should be LONG with -99 when no spec-z available). This code
  will first fit the best-fit templates to the objects with a spec-z.
  Then, it will measure for each filter the systematic offset which
  minimizes the differences between the predicted and observed
  magnitudes. This procedure is applied iteratively until convergence of
  the systematic offset values (maximum of 10 iterations).
| You can also decide to train the zero-points with a sub-sample of the
  spec-z sample. Galaxies can be selected in a given apparent magnitude
  range (``ADAPT_BAND`` and ``ADAPT_LIM``), in a given redshift range
  (``ADAPT_ZBIN``), in a given model range (``ADAPT_MODBIN``).
| You can decide to train only a specific sub-set of bands which are
  indicated using the keyword ``ADAPT_CONTEXT``.
| If the photometric catalogue contains a large number of objects, you
  can save times by doing the training only on a sub-catalogue with
  spec-z and then apply the offsets by hand to the full catalogue
  (APPLY_SYSSHIFT).
| **Note 1**: for philosophical reason, we decided that these offsets
  are added to the predicted magnitudes (because we don’t know if the
  offsets are due to the imaging, bad knowledge of the filters, bad
  knowledge of the templates). Therefore, if you want to apply them
  directly to the observed magnitude in your catalogue, you need to
  subtract these shifts.
| **Note 2**: when using adaptive mode the redshift, for objects that
  meet the criteria from ADAPT_LIM and ADAPT_ZBIN, is automatically
  fixed to the spectroscopic value during the adaptation, and will be
  let free when adaptation is finished. Do not use the adaption with
  -ZFIX YES.

+----------------+----------------+----------------+----------------+
| **Parameters** | **Type**       | **Default      | *              |
|                |                | val.**         | *Description** |
+================+================+================+================+
| AUTO_ADAPT     | string         | NO             | ZP adaptive    |
|                |                |                | method with    |
|                |                |                | spectro        |
+----------------+----------------+----------------+----------------+
|                | (:math:`n=1`)  |                |                |
+----------------+----------------+----------------+----------------+
| ADAPT_BAND     | integer        | —–             | Reference band |
|                |                |                | for the        |
|                |                |                | selection in   |
|                |                |                | magnitude      |
+----------------+----------------+----------------+----------------+
|                | (:math:`n=1`)  |                |                |
+----------------+----------------+----------------+----------------+
| ADAPT_LIM      | float          | 18.,24.        | Mag range for  |
|                |                |                | spectro in     |
|                |                |                | reference band |
+----------------+----------------+----------------+----------------+
|                | (:math:`n=1`)  |                |                |
+----------------+----------------+----------------+----------------+
| ADAPT_CONTEXT  | integer        | -1             | Context for    |
|                |                |                | bands used for |
|                |                |                | training       |
+----------------+----------------+----------------+----------------+
|                | (:math:`n=1`)  |                | -1 : used      |
|                |                |                | context per    |
|                |                |                | object         |
+----------------+----------------+----------------+----------------+
| ADAPT_ZBIN     | float          | 0.01,6         | Redshift’s     |
|                |                |                | interval used  |
|                |                |                | for training   |
+----------------+----------------+----------------+----------------+
|                | (:math:`n=2`)  |                |                |
+----------------+----------------+----------------+----------------+
| ADAPT_MODBIN   | integer        | 1,1000         | Model’s        |
|                |                |                | interval used  |
|                |                |                | for training   |
+----------------+----------------+----------------+----------------+
|                | (:math:`n=2`)  |                |                |
+----------------+----------------+----------------+----------------+
| APPLY_SYSSHIFT | float          | —–             | Apply          |
|                |                |                | systematic     |
|                |                |                | shifts in each |
|                |                |                | bands          |
+----------------+----------------+----------------+----------------+
|                | (:m            |                | number of      |
|                | ath:`n\le 50`) |                | values must    |
|                |                |                | fit number of  |
|                |                |                | filters        |
+----------------+----------------+----------------+----------------+

Analysing the PDF
^^^^^^^^^^^^^^^^^

We have two methods to extract the information from the fit. Either the
profile likelihood which was the original method in the fortran version
(noted MIN hereafter). Either the Bayesian method (noted BAY hereafter).

The PDF are given in the output file using the keyword ``PDZ_OUT`` for
the root name. You need to indicate the type of PDF you want in output
using the keyword ``PDF_TYPE``:

-  **MIN_ZG** the best :math:`\chi^2` at each redshift step is saved to
   build the function :math:`F(z)=exp[-\chi^2_{min}(z)/2]` (profile
   likelihood);

-  **BAY_ZG** all the probabilities :math:`P=exp(-\chi^2/2)` at each
   redshift step are summed (we marginalize over the redshift).

You can obtain the redshift PDF for the QSO library with similar
keywords **MIN_ZQ** and **MIN_ZQ**. We also propose in output the PDF
for several physical parameters using the BAY approach (sum of
probabilities) with **MASS**, **SFR**, **SSFR**, **AGE**.

The value indicated as ``_BEST`` in the output file are obtained using
the PDF computed with the ``MIN`` method. This PDF is also used to
refine the photo-z ``_BEST`` solution (Z_INTERP YES) with a parabolic
interpolation (Bevington, 1969), and to search for secondary solutions
(``DZ_WIN``, ``MIN_THRES``). The search for a secondary solution is done
by imposing a minimum distance between the two peaks in the PDF
(``DZ_WIN``) and a minimum value with respect to the first peak
(``MIN_THRES``). All the values in the output file indicated as ``_MED``
and ``_MODE`` are derived using the Bayesian method (i.e. summing the
probabilities for a given redshift or physical parameter value). In the
case ``_MED`` , we provide the median of the PDF (the old name \_ML
still works). In the case ``_MODE`` , we provide the main mode of the
PDF.

+----------------+----------------+----------------+----------------+
| **Parameters** | **Type**       | **Default      | *              |
|                |                | val.**         | *Description** |
+================+================+================+================+
| ANALYSIS OF    |                |                |                |
| THE            |                |                |                |
| :math:`PDF(z)` |                |                |                |
+----------------+----------------+----------------+----------------+
| PDZ_TYPE       | string         | NONE           | PDZ in output  |
|                |                |                | [def-BAY].     |
|                |                |                | BAY_ZG sum all |
|                |                |                | probabilities  |
|                |                |                | at a given z.  |
+----------------+----------------+----------------+----------------+
|                |                |                | MIN_ZG takes   |
|                |                |                | ex             |
|                |                |                | p(-chi2_min/2) |
|                |                |                | at a each z.   |
+----------------+----------------+----------------+----------------+
|                |                |                | [BAY_          |
|                |                |                | ZG,BAY_ZQ,MIN_ |
|                |                |                | ZG,MIN_ZQ,MASS |
|                |                |                | ,SFR,SSFR,AGE] |
+----------------+----------------+----------------+----------------+
| PDZ_OUT        | string         | NONE           | Root of the    |
|                |                |                | PDF output     |
|                |                |                | files          |
|                |                |                | [def-NONE]     |
+----------------+----------------+----------------+----------------+
|                |                |                | add            |
|                |                |                | automatically  |
|                |                |                | an extension   |
|                |                |                | [\_z           |
|                |                |                | gbay.prob,...] |
+----------------+----------------+----------------+----------------+
| Z_INTERP       | string         | NO             | Parabolic      |
|                |                |                | interpolation  |
|                |                |                | between        |
|                |                |                | original step  |
|                |                |                | (dz)           |
+----------------+----------------+----------------+----------------+
|                | (n=1)          |                |                |
+----------------+----------------+----------------+----------------+
| DZ_WIN         | float          | 0.25           | “smoothing”    |
|                |                |                | window         |
|                |                |                | function for   |
|                |                |                | 2nd peak       |
|                |                |                | search in F(z) |
+----------------+----------------+----------------+----------------+
|                | (n=1)          |                | (value between |
|                |                |                | 0 to zmax)     |
+----------------+----------------+----------------+----------------+
| MIN_THRES      | float          | 0.1            | threshold for  |
|                |                |                | the detection  |
|                |                |                | of 2nd peak in |
|                |                |                | normalised     |
|                |                |                | F(z)           |
+----------------+----------------+----------------+----------------+
|                | (n=1)          |                | (value between |
|                |                |                | 0 to 1)        |
+----------------+----------------+----------------+----------------+
| PROB_INTZ      | float          | 0.             | redshift       |
|                |                |                | intervalles to |
|                |                |                | compute        |
|                |                |                | probability    |
|                |                |                | from F(z)      |
+----------------+----------------+----------------+----------------+
|                | (              |                | (even number   |
|                | n\ :math:`\le` |                | of values),    |
|                | 100)           |                | output vectors |
|                |                |                | from 0 to 100% |
+----------------+----------------+----------------+----------------+
|                |                |                | 0.-default :   |
|                |                |                | not used       |
+----------------+----------------+----------------+----------------+

| 

Physical parameters
^^^^^^^^^^^^^^^^^^^

| After computing the photometric redshifts, other SED fittings can be
  applied to derive FIR properties, absolute magnitudes or to get
  physical parameters. Often, the photometric redshifts are computed
  first, then the redshift value is fixed with option ``ZFIX YES`` and
  the physical parameters are computed in a second step. The reason for
  this two steps procedure is that the template libraries producing the
  best photo-z are not the same as the ones needed to compute physical
  parameters. However, nothing prevent you for doing the two steps
  together.

Absolute magnitudes
~~~~~~~~~~~~~~~~~~~

| This set of parameters allows the user to specify different methods to
  compute the absolute magnitudes. The absolute magnitudes are computed
  automatically in all the filters of FILTER_LIST. Different methods are
  available :
| :math:`\bullet` MABS_METHOD=0 : A direct method to compute the
  absolute magnitude in a given filter from the apparent magnitude
  measured in the same filter (example:
  :math:`B_{ABS}=B_{obs}-DM(z)-kcor(B)`). This method is extremely
  sensitive to k-correction and to systematic effects in the apparent
  magnitude measurement. This method is likely to be less accurate.
| :math:`\bullet` MABS_METHOD=1 : the goal of this method is to minimize
  the sensitivity to the templates. For example, the absolute magnitude
  in the filter B is computed using the observed apparent magnitude in
  the filter I, which is chosen to be
  :math:`\lambda(I)=\lambda(B)*(1+z)` at :math:`z\sim 0.7` :
  :math:`B_{ABS}= I_{obs} -DM(z=0.7) - (kcor(I) + (B-I)_{ABS})^{template}`.
  This method is described in the appendix of Ilbert et al. (2005). The
  advantage of this method to limit template dependency. Indeed, if
  chosen careful, the term k-correction+color doesn’t depend on the
  template at a given redshift. The drawback of this method is that a
  systematic effect in the observed band will be directly propagated to
  the absolute magnitude (like zero-point calibration, or a band
  systematically with a lower S/N). For this reason, a context
  associated to each filter (MABS_CONTEXT) reduces the filter set used
  for the observed apparent magnitudes (for instance, you don’t want to
  keep in the subset a filter having a large offset between observed and
  predicted magnitude in AUTO_ADAPT).
| :math:`\bullet` MABS_METHOD=2 : used to measure the absolute
  magnitudes in all the rest-frame bands using the observed apparent
  magnitudes always taken in the same observed filter (given by
  MABS_REF). It’s not optimized but you know exactly which filter is
  used to compute the absolute magnitudes. As example if MABS_REF is
  defined as B filter and A could be any filter:
  :math:`A_{ABS}=B_{obs}-DM(z)- kcor(B) + (A-B)_{ABS}^{temp}`
| :math:`\bullet` MABS_METHOD=3 : The absolute magnitudes are directly
  measured from the best-fit template. This method is strongly model
  dependent since you can only derive rest-frame colors which are
  present in your templates. However, a bias affecting the photometry in
  one band could be smooth out.
| :math:`\bullet` MABS_METHOD=4 : imposes the filter depending on the
  redshift. The filters are given in MABS_FILT for the corresponding
  redshift bins listed in MABS_ZBIN.
| The predicted apparent magnitudes and absolute magnitudes can be
  computed in a different set of filters than the standard one. In
  ADDITIONAL_MAG, you can add a different name for the filter file,
  different than the one indicated in FILTER_FILE. New predicted
  apparent and absolute magnitudes (only method 3) will be computed in
  these additional filters.

+----------------+----------------+----------------+----------------+
| **Parameters** | **Type**       | **Default      | *              |
|                |                | val.**         | *Description** |
+================+================+================+================+
| Fixing         |                |                |                |
| redshift       |                |                |                |
+----------------+----------------+----------------+----------------+
| ZFIX           | string         | NO             | Fixed redshift |
|                |                |                | (as defined in |
|                |                |                | CAT_TYPE LONG) |
|                |                |                | and            |
+----------------+----------------+----------------+----------------+
|                | (n=1)          |                | search for     |
|                |                |                | best model     |
+----------------+----------------+----------------+----------------+
| EXTERNALZ_FILE | string         | NONE           | Use the spec-z |
|                |                |                | from an        |
|                |                |                | external file  |
|                |                |                | (format Id,zs) |
+----------------+----------------+----------------+----------------+
|                | (n=1)          |                |                |
+----------------+----------------+----------------+----------------+
| Option to      |                |                |                |
| derive the     |                |                |                |
| absolute       |                |                |                |
| magnitudes     |                |                |                |
+----------------+----------------+----------------+----------------+
| MABS_METHOD    | integer        | 0              | Method used    |
|                |                |                | for absolute   |
|                |                |                | magnitudes in  |
|                |                |                | each filter    |
+----------------+----------------+----------------+----------------+
|                | (n\ :ma        |                | 0 (default):   |
|                | th:`\le`\ 100) |                | mag(filt       |
|                |                |                | er)\ :math:`\r |
|                |                |                | ightarrow M_{A |
|                |                |                | BS}`\ (filter) |
+----------------+----------------+----------------+----------------+
|                |                |                | 1 : mag(best   |
|                |                |                | filt           |
|                |                |                | er)\ :math:`\r |
|                |                |                | ightarrow M_{A |
|                |                |                | BS}`\ (filter) |
+----------------+----------------+----------------+----------------+
|                |                |                | 2 : mag(fixed  |
|                |                |                | filter with    |
|                |                |                | MABS_R         |
|                |                |                | EF)\ :math:`\r |
|                |                |                | ightarrow M_{A |
|                |                |                | BS}`\ (filter) |
+----------------+----------------+----------------+----------------+
|                |                |                | 3 : best SED   |
|                |                |                | :math:`\r      |
|                |                |                | ightarrow M_{A |
|                |                |                | BS}`\ (filter) |
+----------------+----------------+----------------+----------------+
|                |                |                | 4 :            |
|                |                |                | MABS(filter)   |
|                |                |                | derives        |
|                |                |                | according to a |
|                |                |                | fixed filter   |
+----------------+----------------+----------------+----------------+
|                |                |                |    in a fixed  |
|                |                |                | redshift       |
|                |                |                | interval       |
+----------------+----------------+----------------+----------------+
|                |                |                |    as given by |
|                |                |                | MABS_FILT and  |
|                |                |                | MABS_ZBIN      |
+----------------+----------------+----------------+----------------+
| MABS_CONTEXT   | integer        | -1             | Context for    |
|                |                |                | the bands used |
|                |                |                | to derive      |
|                |                |                | :              |
|                |                |                | math:`M_{ABS}` |
+----------------+----------------+----------------+----------------+
|                | (n\ :ma        |                | -1 : used same |
|                | th:`\le`\ 100) |                | context as for |
|                |                |                | photo-z        |
+----------------+----------------+----------------+----------------+
| MABS_REF       | integer        | 0              | Fixed filter   |
|                |                |                | (if            |
|                |                |                | MABS_METHOD=2) |
+----------------+----------------+----------------+----------------+
|                | (n\ :ma        |                | 0 (default)    |
|                | th:`\le`\ 100) |                | means not used |
+----------------+----------------+----------------+----------------+
| MABS_FILT      | integer        | —-             | List of fixed  |
|                |                |                | filters chosen |
|                |                |                | to derive      |
|                |                |                | :              |
|                |                |                | math:`M_{ABS}` |
|                |                |                | in all bands   |
+----------------+----------------+----------------+----------------+
|                | (n\ :ma        |                | according to   |
|                | th:`\le`\ 100) |                | the redshift   |
|                |                |                | bins (if       |
|                |                |                | MABS_METHOD=4) |
+----------------+----------------+----------------+----------------+
| MABS_ZBIN      | float          | —-             | List of        |
|                |                |                | Redshift bins  |
|                |                |                | associated     |
|                |                |                | with MABS_FILT |
+----------------+----------------+----------------+----------------+
|                | (n\ :ma        |                | Even number of |
|                | th:`\le`\ 200) |                | values (if     |
|                |                |                | MABS_METHOD=4) |
+----------------+----------------+----------------+----------------+
| ADDITIONAL_MAG | string         | —-             | name of file   |
|                |                |                | with filters   |
+----------------+----------------+----------------+----------------+

Physical parameters derived from BC03 templates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| Physical parameters are derived as soon as you use a library including
  physical information like the normalisation of the template in stellar
  mass. In *LePHARE++* , such measurement is possible only with the BC03
  templates (but we plan to integrate the PEGASE or MARASTON libraries
  on the long term). You don’t need to turn on any keyword to have these
  measurements. As long as you are using BC03 templates and that the
  corresponding keywords (as MASS_MED, or SFR_MED) appear in the output
  parameter file, you should get the physical parameters in output.
| As for the photo-z, you will find physical parameters measured at the
  minimum :math:`\chi^2` value (indicated with ``_BEST``) and the ones
  obtained by taken the median of the PDF marginalized over the relevant
  parameter.

 FIR libraries
~~~~~~~~~~~~~~

| A set of four FIR libraries are available, and can be used to
  characterize the FIR emission of galaxies assuming that the emission
  is dominated by radiation of dust component heated by star formation
  activity. No implementation of hot dust heated by an AGN component has
  been included yet !
| :math:`\bullet` The user defined the minimal rest-frame wavelength for
  the FIR analysis (FIR_LMIN, default is :math:`\lambda=7\mu m`). The
  global FIR context (FIR_CONT) specifies the set of filters to be used.
  However, the final context will depend on the redshift of the source
  and only filters with :math:`\lambda/(1+z) \ge` FIR_LMIN will be
  considered.
| :math:`\bullet` The contribution from the stellar component can be
  subtracted (FIR_SUBSTELLAR) based on the best galaxy template (used in
  ZPHOTLIB). We arbitrarily add in quadrature the subtracted stellar
  flux in the flux error in a given band, and if the stellar component
  is too large (:math:`F_{obs}-F_{\star}\le 3\sigma_{obs}`) we discard
  the passband in the analysis. When activated, the stellar flux is
  subtracted only if :math:`\lambda/(1+z)\le 25\mu m`, we neglect
  stellar component at longer :math:`\lambda`.
| :math:`\bullet` For each library, we estimate the infrared luminosity
  :math:`L_{IR}=\int_{8\mu m}^{1000\mu m} L_{\lambda} dL_{\lambda}`. In
  most of the case the SED’s distribution is attached to a luminosity.
  However when several FIR bands are available, it can be interesting to
  allow for a free rescaling in order to optimize the SED fitting
  (FIR_FREESCALE , FIR_SCALE).
| :math:`\bullet` the FIR output parameters are described in
  Section `6.1 <#sec:outp>`__. The total IR luminosity :math:`L_{IR}`
  and its uncertainties are derived from the maximum likelihood function
  : :math:`F(L_{IR})=\sum exp(-\chi^2(L_{IR})/2)`.
| If only one passband is available, the FIR parameters luminosity is
  derived from the models with closest predicted flux (no rescaling
  allowed). The median and :math:`\sigma` in :math:`L_{TIR}` is
  estimated from the best models of each library.

+----------------+----------------+----------------+----------------+
| **Parameters** | **Type**       | **Default      | *              |
|                |                | val.**         | *Description** |
+================+================+================+================+
|                |                |                | Additional     |
|                |                |                | Libraries      |
|                |                |                | KEYWORDS       |
+----------------+----------------+----------------+----------------+
| FIR_LIB        | string         | NONE           | Far-IR         |
|                |                |                | libraries      |
|                |                |                | separated by   |
|                |                |                | comma          |
+----------------+----------------+----------------+----------------+
|                | (:             |                |                |
|                | math:`n\le 5`) |                |                |
+----------------+----------------+----------------+----------------+
| FIR_LMIN       | float          | 7.0            | :              |
|                |                |                | math:`\lambda` |
|                |                |                | min for FIR    |
|                |                |                | analysis (in   |
|                |                |                | :math:`\mu m`) |
+----------------+----------------+----------------+----------------+
|                | (:math:`n=1`)  |                |                |
+----------------+----------------+----------------+----------------+
| FIR_CONT       | integer        | -1             | Context for    |
|                |                |                | bands to be    |
|                |                |                | used in Far-IR |
+----------------+----------------+----------------+----------------+
|                | (:math:`n=1`)  |                |                |
+----------------+----------------+----------------+----------------+
| FIR_FREESCALE  | string         | NO             | Allows for     |
|                |                |                | free scaling   |
+----------------+----------------+----------------+----------------+
|                | (:math:`n= 1`) |                |                |
+----------------+----------------+----------------+----------------+
| FIR_SCALE      | integer        | -1             | Context for    |
|                |                |                | bands to be    |
|                |                |                | used for       |
|                |                |                | scaling        |
+----------------+----------------+----------------+----------------+
|                | (:math:`n= 1`) |                |                |
+----------------+----------------+----------------+----------------+
| FIR_SUBSTELLAR | string         | NO             | Removing       |
|                |                |                | stellar        |
|                |                |                | component from |
|                |                |                | best optical   |
|                |                |                | fit            |
+----------------+----------------+----------------+----------------+
|                | (:math:`n =1`) |                |                |
+----------------+----------------+----------------+----------------+

The maximum redshift :math:`z_{max}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We compute the maximum redshift at which a galaxy can be observable
given its SED. This maximum redshift depends on how the sample is
selected in magnitude, at a given redshift. That’s why the user needs to
define the selection criteria of the sample. The :math:`z_{max}` is
computed for each galaxy and is used to compute the :math:`V_{max}` for
the luminosity function in a given reference band.

+----------------+----------------+----------------+----------------+
| **Parameters** | **Type**       | **Default      | *              |
|                |                | val.**         | *Description** |
+================+================+================+================+
|                |                |                | Additional     |
|                |                |                | Libraries      |
|                |                |                | KEYWORDS       |
+----------------+----------------+----------------+----------------+
| LIMITS_ZBIN    | double         | :math:`0,99`   | Redshifts used |
|                |                |                | to split in N  |
|                |                |                | bins,          |
|                |                |                | separated by a |
|                |                |                | coma.          |
+----------------+----------------+----------------+----------------+
|                | (:ma           |                | Need N+1       |
|                | th:`n\le 100`) |                | values (start  |
|                |                |                | with the       |
|                |                |                | minimum        |
|                |                |                | redshift).     |
+----------------+----------------+----------------+----------------+
| L              | integer        | 1              | Band in which  |
| IMITS_MAPP_REF |                |                | the absolute   |
|                |                |                | magnitude is   |
|                |                |                | computed       |
+----------------+----------------+----------------+----------------+
|                | (:math:`n=1`)  |                |                |
+----------------+----------------+----------------+----------------+
| L              | integer        | 1              | Give the       |
| IMITS_MAPP_SEL |                |                | selection band |
|                |                |                | in each        |
|                |                |                | redshift bin.  |
+----------------+----------------+----------------+----------------+
|                | (:mat          |                | Need 1 or N    |
|                | h:`n \le 100`) |                | values.        |
+----------------+----------------+----------------+----------------+
| L              | double         | 90             | Magnitude cut  |
| IMITS_MAPP_CUT |                |                | used in each   |
|                |                |                | redshift bin.  |
+----------------+----------------+----------------+----------------+
|                | (:mat          |                | Need 1 or N    |
|                | h:`n \le 100`) |                | values.        |
+----------------+----------------+----------------+----------------+

The output files and parameters 
-------------------------------

.. _`sec:outp`:

The output parameters
^^^^^^^^^^^^^^^^^^^^^

+---------------+-------------------------------------------+---+---+
| IDENT         | Original IDENT                            |   |   |
+===============+===========================================+===+===+
|               | Best Galaxy model                         |   |   |
+---------------+-------------------------------------------+---+---+
| Z_BEST        | Zphot from minimum chi2 (MIN              |   |   |
|               | distribution)                             |   |   |
+---------------+-------------------------------------------+---+---+
| Z_BEST68_LOW  | Zphot min from :math:`\Delta \chi^2=1.0`  |   |   |
|               | (68%)                                     |   |   |
+---------------+-------------------------------------------+---+---+
| Z_BEST68_HIGH | Zphot max from :math:`\Delta \chi^2=1.0`  |   |   |
|               | (68%)                                     |   |   |
+---------------+-------------------------------------------+---+---+
| Z_BEST90_LOW  | Zphot min from :math:`\Delta \chi^2=2.71` |   |   |
|               | (90%)                                     |   |   |
+---------------+-------------------------------------------+---+---+
| Z_BEST90_HIGH | Zphot max from :math:`\Delta \chi^2=2.71` |   |   |
|               | (90%)                                     |   |   |
+---------------+-------------------------------------------+---+---+
| Z_BEST99_LOW  | Zphot min from :math:`\Delta \chi^2=6.63` |   |   |
|               | (99%)                                     |   |   |
+---------------+-------------------------------------------+---+---+
| Z_BEST99_HIGH | Zphot max from :math:`\Delta \chi^2=6.63` |   |   |
|               | (99%)                                     |   |   |
+---------------+-------------------------------------------+---+---+
| Z_MED         | Zphot from the median of the BAY          |   |   |
|               | distribution                              |   |   |
+---------------+-------------------------------------------+---+---+
| Z_MED68_LOW   | Zphot min at 68% (0.16 quantile of BAY    |   |   |
|               | distribution)                             |   |   |
+---------------+-------------------------------------------+---+---+
| Z_MED68_HIGH  | Zphot max at 68% (0.84 quantile of BAY    |   |   |
|               | distribution)                             |   |   |
+---------------+-------------------------------------------+---+---+
| Z_MED90_LOW   | Zphot min at 90% (0.05 quantile of BAY    |   |   |
|               | distribution)                             |   |   |
+---------------+-------------------------------------------+---+---+
| Z_MED90_HIGH  | Zphot max at 90% (0.95 quantile of BAY    |   |   |
|               | distribution)                             |   |   |
+---------------+-------------------------------------------+---+---+
| Z_MED99_LOW   | Zphot min at 99% (0.01 quantile of BAY    |   |   |
|               | distribution)                             |   |   |
+---------------+-------------------------------------------+---+---+
| Z_MED99_HIGH  | Zphot max at 99% (0.99 quantile of BAY    |   |   |
|               | distribution)                             |   |   |
+---------------+-------------------------------------------+---+---+
| Z_MODE        | Zphot from the mode of the BAY            |   |   |
|               | distribution                              |   |   |
+---------------+-------------------------------------------+---+---+
| Z_MODE68_LOW  | Zphot min at 68% (to encompass 68% of BAY |   |   |
|               | distribution around the mode)             |   |   |
+---------------+-------------------------------------------+---+---+
| Z_MODE68_HIGH | Zphot max at 68% (to encompass 68% of BAY |   |   |
|               | distribution around the mode)             |   |   |
+---------------+-------------------------------------------+---+---+
| Z_MODE90_LOW  | Zphot min at 90% (to encompass 90% of BAY |   |   |
|               | distribution around the mode)             |   |   |
+---------------+-------------------------------------------+---+---+
| Z_MODE90_HIGH | Zphot max at 90% (to encompass 90% of BAY |   |   |
|               | distribution around the mode)             |   |   |
+---------------+-------------------------------------------+---+---+
| Z_MODE99_LOW  | Zphot min at 99% (to encompass 99% of BAY |   |   |
|               | distribution around the mode)             |   |   |
+---------------+-------------------------------------------+---+---+
| Z_MODE99_HIGH | Zphot max at 99% (to encompass 99% of BAY |   |   |
|               | distribution around the mode)             |   |   |
+---------------+-------------------------------------------+---+---+
| CHI_BEST      | lowest galaxy :math:`\chi^2` for galaxy   |   |   |
+---------------+-------------------------------------------+---+---+
| MOD_BEST      | galaxy model for best :math:`\chi^2`      |   |   |
+---------------+-------------------------------------------+---+---+
| EXTLAW_BEST   | Extinction law                            |   |   |
+---------------+-------------------------------------------+---+---+
| EBV_BEST      | E(B-V)                                    |   |   |
+---------------+-------------------------------------------+---+---+
| PDZ_BEST      | Probability distribution :                |   |   |
|               | :math:`\int F(z)dz` between               |   |   |
|               | :math:`z_{best} \pm 0.1(1+z)`             |   |   |
+---------------+-------------------------------------------+---+---+
| SCALE_BEST    | Scaling factor                            |   |   |
+---------------+-------------------------------------------+---+---+
| DIST_MOD_BEST | Distance modulus                          |   |   |
+---------------+-------------------------------------------+---+---+
| NBAND_USED    | Number of band used in the measurement    |   |   |
+---------------+-------------------------------------------+---+---+
| NBAND_ULIM    | Number of band used as upper-limit        |   |   |
+---------------+-------------------------------------------+---+---+
|               | Galaxy Secondary solution from F(z)       |   |   |
|               | function                                  |   |   |
+---------------+-------------------------------------------+---+---+
| Z_SEC         |                                           |   |   |
+---------------+-------------------------------------------+---+---+
| CHI_SEC       |                                           |   |   |
+---------------+-------------------------------------------+---+---+
| MOD_SEC       |                                           |   |   |
+---------------+-------------------------------------------+---+---+
| AGE_SEC       |                                           |   |   |
+---------------+-------------------------------------------+---+---+
| EBV_SEC       |                                           |   |   |
+---------------+-------------------------------------------+---+---+
| ZF_SEC        |                                           |   |   |
+---------------+-------------------------------------------+---+---+
| MAG_ABS_SEC   |                                           |   |   |
+---------------+-------------------------------------------+---+---+
| PDZ_SEC       |                                           |   |   |
+---------------+-------------------------------------------+---+---+
| SCALE_SEC     |                                           |   |   |
+---------------+-------------------------------------------+---+---+

+-----------------+-----------------------------+---+---+
|                 | QSO solutions               |   |   |
+=================+=============================+===+===+
| ZQ_BEST         | Zphot from minimum chi2     |   |   |
|                 | (MIN distribution)          |   |   |
+-----------------+-----------------------------+---+---+
| ZQ_BEST68_LOW   | Zphot min from              |   |   |
|                 | :math:`\Delta \chi^2=1.0`   |   |   |
|                 | (68%)                       |   |   |
+-----------------+-----------------------------+---+---+
| ZQ_BEST68_HIGH  | Zphot max from              |   |   |
|                 | :math:`\Delta \chi^2=1.0`   |   |   |
|                 | (68%)                       |   |   |
+-----------------+-----------------------------+---+---+
| ZQ_BEST90_LOW   | Zphot min from              |   |   |
|                 | :math:`\Delta \chi^2=2.71`  |   |   |
|                 | (90%)                       |   |   |
+-----------------+-----------------------------+---+---+
| ZQ_BEST90_HIGH  | Zphot max from              |   |   |
|                 | :math:`\Delta \chi^2=2.71`  |   |   |
|                 | (90%)                       |   |   |
+-----------------+-----------------------------+---+---+
| ZQ_BEST99_LOW   | Zphot min from              |   |   |
|                 | :math:`\Delta \chi^2=6.63`  |   |   |
|                 | (99%)                       |   |   |
+-----------------+-----------------------------+---+---+
| ZQ_BEST99_HIGH  | Zphot max from              |   |   |
|                 | :math:`\Delta \chi^2=6.63`  |   |   |
|                 | (99%)                       |   |   |
+-----------------+-----------------------------+---+---+
| ZQ_MED          | Zphot from the median of    |   |   |
|                 | the BAY distribution        |   |   |
+-----------------+-----------------------------+---+---+
| ZQ_MED68_LOW    | Zphot min at 68% (0.16      |   |   |
|                 | quantile of BAY             |   |   |
|                 | distribution)               |   |   |
+-----------------+-----------------------------+---+---+
| ZQ_MED68_HIGH   | Zphot max at 68% (0.84      |   |   |
|                 | quantile of BAY             |   |   |
|                 | distribution)               |   |   |
+-----------------+-----------------------------+---+---+
| ZQ_MED90_LOW    | Zphot min at 90% (0.05      |   |   |
|                 | quantile of BAY             |   |   |
|                 | distribution)               |   |   |
+-----------------+-----------------------------+---+---+
| ZQ_MED90_HIGH   | Zphot max at 90% (0.95      |   |   |
|                 | quantile of BAY             |   |   |
|                 | distribution)               |   |   |
+-----------------+-----------------------------+---+---+
| ZQ_MED99_LOW    | Zphot min at 99% (0.01      |   |   |
|                 | quantile of BAY             |   |   |
|                 | distribution)               |   |   |
+-----------------+-----------------------------+---+---+
| ZQ_MED99_HIGH   | Zphot max at 99% (0.99      |   |   |
|                 | quantile of BAY             |   |   |
|                 | distribution)               |   |   |
+-----------------+-----------------------------+---+---+
| ZQ_MODE         | Zphot from the mode of the  |   |   |
|                 | BAY distribution            |   |   |
+-----------------+-----------------------------+---+---+
| ZQ_MODE68_LOW   | Zphot min at 68% (to        |   |   |
|                 | encompass 68% of BAY        |   |   |
|                 | distribution around the     |   |   |
|                 | mode)                       |   |   |
+-----------------+-----------------------------+---+---+
| ZQ_MODE68_HIGH  | Zphot max at 68% (to        |   |   |
|                 | encompass 68% of BAY        |   |   |
|                 | distribution around the     |   |   |
|                 | mode)                       |   |   |
+-----------------+-----------------------------+---+---+
| ZQ_MODE90_LOW   | Zphot min at 90% (to        |   |   |
|                 | encompass 90% of BAY        |   |   |
|                 | distribution around the     |   |   |
|                 | mode)                       |   |   |
+-----------------+-----------------------------+---+---+
| ZQ_MODE90_HIGH  | Zphot max at 90% (to        |   |   |
|                 | encompass 90% of BAY        |   |   |
|                 | distribution around the     |   |   |
|                 | mode)                       |   |   |
+-----------------+-----------------------------+---+---+
| ZQ_MODE99_LOW   | Zphot min at 99% (to        |   |   |
|                 | encompass 99% of BAY        |   |   |
|                 | distribution around the     |   |   |
|                 | mode)                       |   |   |
+-----------------+-----------------------------+---+---+
| ZQ_MODE99_HIGH  | Zphot max at 99% (to        |   |   |
|                 | encompass 99% of BAY        |   |   |
|                 | distribution around the     |   |   |
|                 | mode)                       |   |   |
+-----------------+-----------------------------+---+---+
| CHI_QSO         |                             |   |   |
+-----------------+-----------------------------+---+---+
| MOD_QSO         |                             |   |   |
+-----------------+-----------------------------+---+---+
| Z_ML            | Zphot from Median of ML     |   |   |
|                 | distribution                |   |   |
+-----------------+-----------------------------+---+---+
| Z_ML68_LOW_QSO  | Zphot min at 68% from       |   |   |
|                 | Median of of ML             |   |   |
|                 | distribution                |   |   |
+-----------------+-----------------------------+---+---+
| Z_ML68_HIGH_QSO | Zphot max at 68% from       |   |   |
|                 | Median of of ML             |   |   |
|                 | distribution                |   |   |
+-----------------+-----------------------------+---+---+
|                 | Star solutions              |   |   |
+-----------------+-----------------------------+---+---+
| MOD_STAR        |                             |   |   |
+-----------------+-----------------------------+---+---+
| CHI_STAR        |                             |   |   |
+-----------------+-----------------------------+---+---+

+---------------+-------------------------------------------+---+---+
|               | Best models magnitudes                    |   |   |
+===============+===========================================+===+===+
| MAG_MOD()     | magnitudes from best models               |   |   |
+---------------+-------------------------------------------+---+---+
| K_COR()       | K-corrections                             |   |   |
+---------------+-------------------------------------------+---+---+
| MAG_ABS()     | absolute magnitudes                       |   |   |
+---------------+-------------------------------------------+---+---+
| MABS_FILT()   | Filter used for each absolute magnitude   |   |   |
+---------------+-------------------------------------------+---+---+
| MAG_PRED()    | Apparent magnitudes computed in a spare   |   |   |
|               | filter set                                |   |   |
+---------------+-------------------------------------------+---+---+
| ABSMAG_PRED() | Absolute magnitudes computed in a spare   |   |   |
|               | filter set                                |   |   |
+---------------+-------------------------------------------+---+---+
|               | PHYSICAL PARAMETERS from ZPHOTLIB galaxy  |   |   |
|               | library (from file .phys)                 |   |   |
+---------------+-------------------------------------------+---+---+
|               | with XXX=BEST (minimum chi2), MED (median |   |   |
|               | L(z)), INF (associated to \_MED), SUP     |   |   |
|               | (associated to \_MED)                     |   |   |
+---------------+-------------------------------------------+---+---+
| EBV_XXX       |                                           |   |   |
+---------------+-------------------------------------------+---+---+
| AGE_XXX       | Log10(Age[yr])                            |   |   |
+---------------+-------------------------------------------+---+---+
| MASS_XXX      | log10(Mass (:math:`M_{\odot}`) )          |   |   |
+---------------+-------------------------------------------+---+---+
| SFR_XXX       | log10(SFR ( :math:`M_{\odot}/yr`))        |   |   |
+---------------+-------------------------------------------+---+---+
| SSFR_XXX      | log10(SSFR ( :math:`yr^{-1}`) )           |   |   |
+---------------+-------------------------------------------+---+---+
| LDUST_XXX     | log10(:math:`L_{dust}` (                  |   |   |
|               | :math:`L_{\odot}`) )                      |   |   |
+---------------+-------------------------------------------+---+---+
| LUM_NUV_BEST  | log10(:math:`L_{UV}` ( :math:`L_{\odot}`) |   |   |
|               | )                                         |   |   |
+---------------+-------------------------------------------+---+---+
| LUM_R_BEST    | log10(:math:`L_R` ( :math:`L_{\odot}`))   |   |   |
+---------------+-------------------------------------------+---+---+
| LUM_K_BEST    | log10(:math:`L_K` ( :math:`L_{\odot}`))   |   |   |
+---------------+-------------------------------------------+---+---+
|               | L(TIR) from FIR fit                       |   |   |
+---------------+-------------------------------------------+---+---+
| LTIR_XXX      | log10(:math:`L_{TIR (8->1000 mum}` (      |   |   |
|               | :math:`L_{\odot}`) )                      |   |   |
+---------------+-------------------------------------------+---+---+
|               | Input informations                        |   |   |
+---------------+-------------------------------------------+---+---+
| CONTEXT       |                                           |   |   |
+---------------+-------------------------------------------+---+---+
| ZSPEC         |                                           |   |   |
+---------------+-------------------------------------------+---+---+
| MAG_OBS()     | input observed magnitudes                 |   |   |
+---------------+-------------------------------------------+---+---+
| ERR_MAG_OBS() | input errors                              |   |   |
+---------------+-------------------------------------------+---+---+
| STRING_INPUT  | All the information after zspec in the    |   |   |
|               | input file                                |   |   |
+---------------+-------------------------------------------+---+---+

Appendix A : Content of the main tar file (lephare_main.tar.gz)
----------------------------------------------------------------

+-------------------------+-------------------------------------------+
| $LEPHAREDIR/source/ :   | fortran source files (\*.f)               |
+-------------------------+-------------------------------------------+
| $LEPHAREDIR/vega/ :     | SEDs for reference stars:                 |
+-------------------------+-------------------------------------------+
|                         | VegaLCB.dat from Lejeune et al., (1997)   |
+-------------------------+-------------------------------------------+
|                         | BD\ :math:`+17^{\circ}4708` from Oke &    |
|                         | Gunn, (1983)                              |
+-------------------------+-------------------------------------------+
| $LEPHAREDIR/doc/ :      | documentation files                       |
+-------------------------+-------------------------------------------+
| $LEPHAREDIR/config/ :   | configuration files for photometric       |
|                         | redshifts                                 |
+-------------------------+-------------------------------------------+
| $LEPHAREDIR/ext/ :      | extinction curves for internal galaxy     |
|                         | extinction (see README.extinc)            |
+-------------------------+-------------------------------------------+
| $LEPHAREDIR/opa/ :      | Lyman absorption files produced by        |
|                         | intergalactic medium as function          |
+-------------------------+-------------------------------------------+
|                         | of redshift from Madau (1995).            |
+-------------------------+-------------------------------------------+
| $LEPHAREDIR/filt/ :     | individual filter response curves sorted  |
|                         | by directories:                           |
+-------------------------+-------------------------------------------+
|                         | (ex : ./std/\*.pb, ./sdss/\*.pb,          |
|                         | ./hst/\*.pb, ...)                         |
+-------------------------+-------------------------------------------+
| $LEPHAREDIR/sed/STAR/:  | SED files for STARS sorted by             |
|                         | directories:                              |
+-------------------------+-------------------------------------------+
|                         | ./BD/\*.sed : Brown Dwarves and low mass  |
|                         | stars from Chabrier & Baraffe (2000)      |
+-------------------------+-------------------------------------------+
|                         | ./PICKLES/\*.sed: Atlas from Pickles      |
|                         | (1998)                                    |
+-------------------------+-------------------------------------------+
|                         | ./WD/\*.sed : 4 White dwarves from IRAF   |
|                         | ??                                        |
+-------------------------+-------------------------------------------+
| $LEPHAREDIR/sed/QSO/:   | SED files for QSO                         |
+-------------------------+-------------------------------------------+
|                         | ./qsol.dat : composite template from      |
|                         | Cristiani (?)                             |
+-------------------------+-------------------------------------------+
|                         | ./qso-\*.dat: SEDs based on simulated     |
|                         | slope and strengthes of emission          |
+-------------------------+-------------------------------------------+
|                         | lines (Warren ?) (see README.model_qso)   |
+-------------------------+-------------------------------------------+
|                         | ./a\*.sed : SEDs from Eva Hatziminaoglou  |
|                         | (?)                                       |
+-------------------------+-------------------------------------------+
| $LEPHAREDIR/sed/GAL/:   | SED files for GALAXY sorted by            |
|                         | directories:                              |
+-------------------------+-------------------------------------------+
|                         | CWW_KINNEY/\*.sed: observed SED by        |
|                         | Coleman etal(1980)                        |
+-------------------------+-------------------------------------------+
|                         | + 2 starburst from Kinney (1996) extended |
|                         | in UV and IR by GISSEL                    |
+-------------------------+-------------------------------------------+
|                         | CE/\*.sed : 72 SEDs based on              |
|                         | interpolation between the 4 CWW SEDS      |
+-------------------------+-------------------------------------------+
|                         | + 3 young star-forming galaxies from      |
|                         | GISSEL (Arnouts et al., 1999)             |
+-------------------------+-------------------------------------------+
|                         | 42GISSEL/\*.sed : 42 selected SEDs from   |
|                         | GISSEL                                    |
+-------------------------+-------------------------------------------+
|                         | Additional directories with the other SED |
|                         | tar-files:                                |
+-------------------------+-------------------------------------------+
|                         | GISSEL/\*.sed : 27 theoretical SEDs       |
|                         | including ages from GISSEL96              |
+-------------------------+-------------------------------------------+
|                         | PEGASE2/\*.sed : SEDs including ages from |
|                         | PEGASE2 (Fioc et al., 2000)               |
+-------------------------+-------------------------------------------+
|                         | HYPERZ/\*.sed : SEDs from the HYPERZ      |
|                         | package (Bolzonella et al., 2000)         |
+-------------------------+-------------------------------------------+
| $LEPHAREDIR/tools/ :    | macros for plotting filterset             |
|                         | (filterplot) and extracted SED (spec.sm)  |
+-------------------------+-------------------------------------------+
| $LEPHAREDIR/test/ :     | Quick test with COSMOS2015                |
+-------------------------+-------------------------------------------+
| $LEPHAREDIR/examples/ : | extensive list of tests/commands with     |
|                         | COSMOS.                                   |
+-------------------------+-------------------------------------------+
| $LEPHAREDIR/simul/ :    | configuration files for the simulations   |
|                         | and predefined input files                |
+-------------------------+-------------------------------------------+

Appendix B : keyword differences between the Fortran and the C++ version
-------------------------------------------------------------------------

| We list here the keywords with a new format that you should modify to
  use the c++ version. We don’t list new keywords which correspond to
  new features described in this documentation.
| **Z_STEP format**
| fortran - dz,zmax,dzsup (zmin=0, dzsup step for :math:`z>6`)
| c++ - dz,zmin,zmax
| **ZGRID_TYPE new**
| c++ - 0(def)/1
| **EM_LINES (in mag_gal) format**
| fortran - NO(def)/YES
| c++ - NO(def)/EMP_UV/EMP_SFR/PHYS
| **EM_DISPERSION format**
| fortran - 50,25 (not documented)
| c++ - 0.25,0.5,1,2,4
| **ADD_EMLINES (in zphota) format**
| fortran - NO(def)/YES
| c++ - 0,0(def)

.. [1]
   If BC3-like templates are used, :math:`z_\mathrm{phot}` and galaxy
   physical parameters can be computed at the same time

.. [2]
   A more appropriate name for the directory would be AGN, but QSO is
   kept for compatibility reasons with original work.

.. [3]
   :math:`F_\lambda` or F does not matter, as a normalization will be
   applied.

.. [4]
   convertion from absolute magnitude to luminosity: by combining the
   distance modulus (:math:`m-M=5log D[pc]-5`) and the luminosity
   (:math:`L_{\nu}[erg.s^{-1}.Hz^{-1}]= A\pi D^2[cm^2] f_{\nu}[erg.s^{-1}.cm^{-2}.Hz^{-1}]`)
   you get:
   :math:`L_{\nu,\odot}[erg.s^{-1}.Hz^{-1}]=10^{-0.4(M^{AB}_{\nu,\odot} -51.605)}`
