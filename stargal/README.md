StarGal - LePHARE for star-galaxy separation.

Initially designed for estimating galaxy redshifts & intrinsic properties, we want to explore the uses and efficiency of LePHARE as a star-galaxy classifier (stargal). StarGal is an extension of LePHARE mainly functioning in Python, where only a few modifications to the original source code in C++ have been made (only to retrieve the fit parameters of each source to the star library). Two aspects stand out in StarGal:
 - A series of statistical tests to analyze the behavior of LePHARE over a star SED library.
 - The development of SPLIT (Source-type Probability with Likelihood Identification Tool), which is the main code to determine whether a source is a star or a galaxy (or more, like QSO) using likelihood data.


The stargal folder contains:
 - **`notebooks{.ipynb}`**:
    - We advise to do your tests here.
    - `cat_to_lephare_input`: Convert catalogs to `.dat` format to make them suitable for LePHARE.
    - `chi2_stargal_separation`: Example of χ² comparison for stargal (first part of SPLIT).
    - `chi2_stats`: Look at the z and χ² distributions.
    - `plot_filters`: Plot the filter libraries.
    - `plot_seds`: Plot the SED libraries.
    - `SPLIT_dev`: Implementation of SPLIT in a notebook.
    - `star_pdf_flagger_part1`: Analyzing the behavior of source classes over physical parameters of stars.
    - `star_pdf_flagger_part2`: Use the results of "star_pdf_flagger_part1" to classify sources.
    - `lephare_full_run`: Run LePHARE with the Python interface.
    - `statistic_tests_for_sedlib_improvement`: To verify if a richer star SED library improves the separation.
    - `imports.py`: The usual packages to import at the beginning of a notebook. 


 - **`scripts{.py}`**:
    - The scripts holding functions to gain writing time.
    - `format_to_dat`: Holds the functions used in "cat_to_lephare_input.ipynb".
    - `statsplot`: Plot functions.
    - `utils`: Utilities for data manipulation.
    - `split_src`: Functions to run SPLIT.
    
 - **`simulated_sed`**:
    - Expected folders: `bt_spectra`, `phoenix_spectra`, `buzzard_base`.
    - `download_bt_spectra`: Script to download `BTNextGen_spectra`.
    - `download_phoenix_spectra`: Script to download `phoenix_spectra`.

 - **`catalogs`**: 
    - We advise to put the catalogs an the result of your tests here.
    - `DC1`: Simulated galaxies from LSST Data Challenge 1.
    - `DC2`: Simulated stars and galaxies from LSST Data Challenge 2.
    - `DES`: Reference star catalogs from DES.
    - `DP1`: Observed sources from LSST Data Preview 1.

 - **`others`**:
    - `usual_commands.txt`: Usual bash commands to run in the terminal to execute LePHARE.
    - `stargal_config.para`: Config file for stargal with LePHARE.
    - `stargal_output.para`: Output configuration for stargal with LePHARE.