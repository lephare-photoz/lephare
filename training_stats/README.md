Lephare statistics by @github/Tototime

A collection of notebooks and scripts to perform statistics on the output catalog from lephare. 4 folders added in training stats:
 - simulation_catalogs which is suppose to holds filters, SEDs and catalogs on which Lephare ran.
 - output_data where belong output data from lephare.
 - scripts
 - notebooks and scripts with codes on simulation_catalogs and output_data.
 

Simulation_catalogs:
 - not on github 
 - originally to perform training algorithm on lephare.
 - see in-folder README.

notebooks:
 - fits_to_dat: few manipulation on catalogs in .fits file format. 
 - DC1_results_stats: main notebook to perform statistics on photoz results on DC1 simulated catalog
 - plot_SEDs
 - plot_filters
 
scripts:
 - apply_context.py: apply values to context column in the input.dat catalog for zphota. Use threshold in magnitudes.
 - output_catalog_flagger.py: use zphota output catalog and pdz_out to add flag to PDZ (associated to z_best or z_mode values) depending on the quality of the PDZ.
 - stats_utils.py: future useful functions for stats.
 - spec.py: from lephare
 - make_SEDlib_DES_stars: run SEDtoLib with python interface for PICKLES SEDs and filter
 - make_SEDlib_DC1

data (examples):
 - catalogs
 - pdzs
 - specs
 - spec_pdfs
 - output.para
 - config_file.para
 
Any question about the data used, contact : theotime.hallouin@uca.fr - @github/Tototime
