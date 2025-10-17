
.. _allkeywords-label:

List of keywords available for configuration of the LePhare executables 
=======================================================================

All the input ascii files can have comment lines starting with #. A single configuration file can be used for all the executables.

filter
------

.. list-table:: Keywords
   :widths: 25 25 50
   :header-rows: 1

   * - Keyword name
     - values
     - Comment
   * - FILTER_REP
     - ``$LEPHAREDIR/filt``
     - Repository in which the filters are stored. Absolute path.
   * - FILTER_LIST
     - e.g. ``hsc/gHSC.pb,hsc/rHSC.pb``
     - Comma-separated list of filter file paths relative to ``FILTER_REP``. These files are expected to have two columns: wavelengths and corresponding transmission.
   * - TRANS_TYPE
     - 0[def] or 1 
     - Define the transmission as being of type energy (0), or photon count (1). See ``flt::trans``
   * - FILTER_CALIB
     - 0[def],1,2,3,4,5 
     - comma-separated list of integer corresponding in order exactly to ``FILTER_LIST``, and defining which kind of calibration to apply to the filter transmission curves. 
   * - FILTER_FILE
     - e.g. ``filter_cosmos``
     - output filename written in ``$LEPHAREWORK/filt``

sedtolib
--------

.. list-table:: Keywords
   :widths: 25 25 50
   :header-rows: 1

   * - Keyword name
     - values
     - Comment
   * - GAL_SED
     - e.g. ``examples/COSMOS_MOD.list``
     - Name of file containing the list of Galaxy SED files to include. Absolute or relative to ``$LEPHAREWORK``
   * - GAL_FSCALE
     - float
     - Arbitrary flux scale to be applied to each SED flux in the list
   * - GAL_LIB
     - string
     - basename of the Galaxy SED library output bin and doc files under ``$LEPHAREWORK/lib_bin``
   * - SEL_AGE
     - 
     - filename of the file containing the ages to consider (in Gyr)
   * - AGE_RANGE
     - min,max 
     - age range to be considered (in Gyr)
   * - QSO_SED
     - e.g. ``sed/QSO/SALVATO09/AGN_MOD.list``
     - name of file containing the list of AGN/QSO SED files to include
   * - QSO_FSCALE
     - float
     - Arbitrary flux scale to be applied to each SED flux in the list
   * - QSO_LIB
     - string
     - basename of the AGN/QSO SED library output bin and doc files under ``$LEPHAREWORK/lib_bin``
   * - STAR_SED
     - e.g. ``sed/STAR/STAR_MOD_ALL.list``
     - name of file containing the list of Star SED files to include
   * - STAR_LIB
     - string
     - basename of the Star SED library output bin and doc files under ``$LEPHAREWORK/lib_bin``
   * - STAR_FSCALE
     - float
     - Arbitrary flux scale to be applied to each SED flux in the list

mag_gal
-------

.. list-table:: Keywords
   :widths: 25 25 50
   :header-rows: 1

   * - Keyword name
     - values
     - Comment
   * - COSMOLOGY
     - H0,OmegaM,OmegaL
     - Cosmological parameters
   * - FILTER_FILE
     -  
     - same as above
   * - MAGTYPE
     - AB[def] or VEGA 
     - Magnitude type (AB or VEGA) of the library
   * - EXTINC_LAW
     - calzetti.dat[def] 
     - Extinction laws stored in ``$LEPHAREWORK/ext/``. Could include several files separated by comma.
   * - EB_V
     - 0[def] 
     - Reddening color excess E(B-V) values to be applied. Values separated by comma.
   * - MOD_EXTINC
     - 0,0[def] 
     - Range of models for which extinction will be applied. Two values per attenutation curves, separated by comma.
   * - Z_STEP
     - 0.04,0,6 [def] 
     - dz,zmin,zmax: redshift step (dz), the minimum (zmin) and the maximum redshift (zmax)
   * - GAL_LIB_IN
     - string 
     - Name of the GALAXY library (with no extension). Files must exist in $LEPHAREWORK/lib_bin/ and have been created by a sedtolib run.
   * - QSO_LIB_IN
     - string 
     - Name of the QSO library (with no extension). Files must exist in $LEPHAREWORK/lib_bin/ and have been created by a sedtolib run.
   * - STAR_LIB_IN
     - string 
     - Name of the STAR binary library (with no extension). Files must exist in $LEPHAREWORK/lib_bin/ and have been created by a sedtolib run
   * - GAL_LIB_OUT
     - string 
     - Name of the output library (with no extension) storing predicted fluxes/magnitude/k-corrections. Files .bin and .doc are saved in $LEPHAREWORK/lib_mag/
   * - QSO_LIB_OUT
     - string 
     - Name of the output library (with no extension) storing predicted fluxes/magnitude/k-corrections. Files .bin and .doc are saved in $LEPHAREWORK/lib_mag/
   * - STAR_LIB_OUT
     - string 
     - Name of the output library (with no extension) storing predicted fluxes/magnitude/k-corrections. Files .bin and .doc are saved in $LEPHAREWORK/lib_mag/
   * - LIB_ASCII
     - NO[def] or YES 
     - Writing the magnitude library in ASCII file to $LEPHAREWORK/lib_mag/
   * - EM_LINES
     - NO[def] or EMP_UV or EMP_SFR or PHYS 
     - choice between emission line prescription
   * - EM_DISPERSION
     - 1[def] 
     - Dispersion allowed in the emission line flux factor. Example: 0.5,0.75,1,1.5,2
   * - ADD_DUSTEM
     - NO[def] or YES
     - Add the dust emission in templates when missing (using energy balance).

zphota
------

.. list-table:: Keywords
   :widths: 25 25 50
   :header-rows: 1

   * - Keyword name
     - values
     - Comment
   * - CAT_IN
     - string
     - Absolute path to photometric catalog input. This must have the correct column orders and be a space seperated ascii file with -99 missing values. It is not required when running throught the Python interface.
   * - INP_TYPE
     - F or M 
     - Input  values:  Flux (F) or Magnitude (M)
   * - CAT_TYPE
     - LONG[def] or SHORT 
     - Input catalog format (long requires at minimum spec-z and context)
   * - CAT_MAG
     - AB[def] or VEGA 
     - Input magnitude type
   * - CAT_FMT
     - MEME[def] or MMEE  
     - Input format for photometry (MEME alternate mag-flux with error)
   * - CAT_LINES
     - 0,10000000000[def] 
     - Min and max rows read in input catalog (to run only a subsample)
   * - PARA_OUT
     - e.g. examples/output.para
     - Absolute path of the file with selected output parameters.
   * - CAT_OUT
     - path to output table.
     - Name of the output file
   * - ZPHOTLIB
     - e.g. CE_COSMOS, STAR_COSMOS, QSO_COSMOS
     - Library names (with no extension). Could have several separated by comma. Should be in LEPHAREWORK/lib_mag. Names must correspond to outputs from mag_gal.
   * - ADD_EMLINES
     - e.g. ``0,10000``
     - Range of galaxy models in which considering emission lines contribution.
   * - EBV_RANGE
     - 0,9[def] 
     - E(B-V) min and max allowed in the GAL library
   * - ERR_SCALE
     - 0[def] 
     - Systematic errors (in mag) add in quadrature to the observations. One per filter, separated by comma.
   * - ERR_FACTOR
     - 1.0[def]  
     - Scaling factor to the errors (in flux). Only a single value applied to all filters.
   * - BD_SCALE
     - 0[def] 
     - Band used for scaling  the models to the observations (sum of $2^i$ starting at i=0, as context). 0 means all.
   * - GLB_CONTEXT
     - 0[def] 
     - Forces the context of all objects (sum of 2^i starting at i=0, as context). 0 means all.
   * - FORB_CONTEXT
     - 0[def] 
     - Context for removing some bands from the fit (sum of 2^i starting at i=0, as context). 0 means inactive.
   * - MASS_SCALE
     - 0, 0[def]
     - Prior: allowed range in log10(mass)
   * - MAG_ABS
     - 0, 0[def]
     - Prior: Absolute magnitude range allowed for the GAL library [0,0-def]
   * - MAG_ABS_QSO
     - 0, 0[def]
     - Prior: Absolute magnitude range allowed for the QSO library [0,0-def]
   * - MAG_REF
     - int
     - Reference filter for the prior in abs. mag. (start at 1)
   * - NZ_PRIOR
     - -1[def] 
     - N(z) prior as function of i-band. The i-band number should be given in input (starting filter numbering at 1). The second number indicates which band to use if first undefined (not mandatory). Negative value means no prior.
   * - ZFIX
     - NO[def] or YES 
     - Fixed redshift with the spec-z value (as defined in CAT_TYPE LONG)
   * - EXTERNALZ_FILE
     - string 
     - Name of an external file. Use the spec-z from an external file (format Id,zs) to fix the redshift.
   * - Z_INTERP
     - NO[def] 
     - Parabolic interpolation between original step (dz)
   * - DZ_WIN
     - 0.25[def] 
     - Window function for 2nd peak search in L(z) (minimal distance in dz from the 1st peak)
   * - MIN_THRES
     - 0.1[def]
     - Threshold for the detection of 2nd peak in normalised L(z) (between 0 and 1)
   * - SPEC_OUT
     - NO[def]  
     - Output files with Gal/Star/QSO spectra (one file per object) (if YES: can take a lot of disk space !)
   * - CHI2_OUT
     - NO[def]  
     - Output files with the chi2 for the full library (one file per object) (if YES: can take a lot of disk space !)
   * - PDZ_OUT
     - NONE[.pdz] 
     - Output file name in which PDZ will be stored (full path). The code will add automatically the extension[.pdz]
   * - PDZ_TYPE
     - BAY[def]  or MIN
     - value: BAY\_ZG[def] or/and BAY\_ZQ,MIN\_ZG,MIN\_ZQ,MASS,SFR,SSFR,AGE PDZ in output [def-BAY]. BAY\_ZG sum all probabilities at a given z. MIN_ZG takes ex p(-chi2_min/2) at a each z.
   * - FIR_LIB
     - NONE[def] 
     - Far-IR libraries separated by comma
   * - FIR_LMIN
     - 7[def] 
     - $\lambda$ min for FIR analysis (in $\mu m$)
   * - FIR_CONT
     - -1[def] 
     - Context for bands to be used in Far-IR
   * - FIR_SCALE
     - -1[def]  
     - Context for bands to be used for scaling
   * - FIR_FREESCALE
     - NO[def] 
     - Allows for free scaling
   * - FIR_SUBSTELLAR
     - NO[def] 
     - Removing stellar component from best optical fit
   * - MABS_METHOD
     - 0[def], 1,2,3, or 4 
     - Method used for absolute magnitudes in each filter. 0 (default): mag(filter) -> Mabs(filter). 1 : mag(best filter) -> Mabs(filter). 2 : mag(fixed filter with MABS_REF)-> Mabs(filter). 3 : best SED -> Mabs(filter). 4 : MABS(filter) derives according to a fixed filter in a fixed redshift interval as given by MABS_FILT and MABS_ZBIN
   * - Z_METHOD
     - BEST[def] or ML 
     - Compute the absolute magnitude at a given redshift solution maginalised over all models (ML) or for the best model which minimizes the chi2 (BEST)
   * - MABS_CONTEXT
     - 0[def]  
     - Context for the bands used to derive  Mabs.
   * - MABS_REF
     - 1[def]  
     - Filter in observed frame used to derive all the Mabs if method=2
   * - MABS_FILT
     - 1[def]   
     - For method 4: list of  fixed filters chosen to derive Mabs in all bands according to the redshift bins
   * - MABS_ZBIN
     - 0,6[def] 
     - For method 4: list of Redshift bins associated with  MABS_FILT. Even number of values.
   * - ADDITIONAL_MAG
     - string
     - Name of file compiling several filters (in $LEPHAREWORK/filt, created by filter)to derive Mabs in additional filters
   * - APPLY_SYSSHIFT
     - 0[def] 
     - Apply systematic shifts in each bands (convention: values in magnitude to be substracted to the observed magnitudes). Number of values must correspond to the number of filters.
   * - AUTO_ADAPT
     - NO[def] 
     - Optimize zero-points with spec-z 
   * - ADAPT_BAND
     - 1[def] 
     - Reference band for the selection in magnitude (start at 1)
   * - ADAPT_LIM
     - 15,35[def] 
     - Mag range for spectro in reference band
   * - ADAPT_ZBIN
     - 0.01,6[def] 
     - Redshift's interval used for training
   * - LIMITS_ZBIN
     - 0.0,90.[def] 
     - Redshift limits used to split in N bins, separated by a coma.
   * - LIMITS_MAPP_REF
     - 1[def] 
     - Compute z-max. Band in which the absolute magnitude is computed.
   * - LIMITS_MAPP_SEL
     - 1[def] 
     - Compute z-max. Give the selection band in each redshift bin.  Need 1 or N values.
   * - LIMITS_MAPP_CUT
     - 90[def] 
     - Compute z-max. Magnitude cut used in each redshift bin. Need 1 or N values.
   * - RM_DISCREPENT_BD
     - 200[def] 
     - Threshold in chi2 to stop removing bands. Remove 2 bands max, stop when below this chi2 threshold.
   * - Z_RANGE
     - 0.,99.[def] 
     - Z min and max allowed in the GAL library
