Lephare statistics by @github/Tototime

A collection of notebooks and scripts to perform statistics on the output catalog from lephare. 4 folders added in training stats:
 - simulation_catalogs holding filters, SEDs and catalogs on which Lephare ran.
 - output_data where belong output data from lephare.
 - notebooks and scripts with codes on simulation_catalogs and output_data.
 

Simulation_catalogs:
 - not on github
 - originally to perform training algorithm on lephare.
 - see in-folder README.

notebooks:
 - fits_to_dat.ipynb: few manipulation on catalogs in .fits file format. 
 - lephare_results_tats.ipynb: main notebook to perform statistics.
 - lephare_training.ipynb: to make filter ans SED libraries via python interface.
 
scripts:
 - apply_context.py: apply values to context column in the input.dat catalog for zphota. Use threshold in magnitudes.
 - filter_plot.py: plot filters from simulation_catalogs.
 - SED_plot.py: plot some SEDs from simulation_catalogs.
 - output_catalog_flagger.py: use zphota output catalog and pdz_out to add flag to PDZ (associated to z_best or z_mode values) depending on the quality of the PDZ.
 - stats_utils.py: useful functions used in lephare_results_tats.
 - spec.py: makes plots on given .spec files from zphota (after specified SPEC_OUT keyword in zphota). Ouput is a .pdf. Originally from Lephare-data.

data:
 - catalogs
 - pdzs
 - specs
 - spec_pdfs
 - output.para
 - config_file.para
 
