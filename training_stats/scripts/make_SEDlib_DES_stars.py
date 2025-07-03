import os
import lephare as lp
import numpy as np
from matplotlib import pylab as plt

##### Create filter library #####
base_dir = os.path.abspath(os.path.join(os.getcwd(), '..')) #change to your lephare base_dir
filter_path = os.path.join(base_dir, 'lephare/training_stats/simulation_catalogs/DES/DES_STARCAT/WORK_COMPLETE2/filt') #paste your relative filter path
print(filter_path)
filter_rep = lp.keyword("FILTER_REP", filter_path)
filter_list = lp.keyword("FILTER_LIST",
                         "DES_filter_g.res,DES_filter_r.res,\
                         DES_filter_i.res,DES_filter_z.res,DES_filter_Y.res") #edit your filter list
#crash if add comma

filterLib = lp.Filter(config_keymap={"FILTER_REP":filter_rep,
                                   "FILTER_LIST":filter_list,
                                   "TRANS_TYPE":lp.keyword("TRANS_TYPE", "0"),
                                   "FILTER_CALIB":lp.keyword("FILTER_CALIB", "0"),
                                   "FILTER_FILE":lp.keyword("FILTER_FILE", "photozDES")})


filterLib.run()



#### Create SED library #####
#ready in lephare data

SED_list_path = f"{base_dir}/lephare/training_stats/simulation_catalogs/DES/DES_STARCAT/WORK_COMPLETE2/lib_bin/STAR_MOD.list"
sedLib = lp.Sedtolib(config_keymap={
    "STAR_SED": lp.keyword("STAR_SED", SED_list_path),
    "STAR_FSCALE": lp.keyword("STAR_FSCALE", "3.0e-9"),
    "STAR_LIB": lp.keyword("STAR_LIB", 'PICKLES_SEDs')},)
sedLib.run(typ="S")

# #if you use gal seds
# SED_list_path = f"{base_dir}/lephare/training_stats/simulation_catalogs/buzzard_base/SEDS/updated_Buzzard_SEDs/updated_Buzzard_SEDs.list" 
# sedLib = lp.Sedtolib(config_keymap={
#     "GAL_SED": lp.keyword("GAL_SED", SED_list_path),
#     "GAL_FSCALE": lp.keyword("GAL_FSCALE", "1.0"),
#     "GAL_LIB": lp.keyword("GAL_LIB", 'buzzard_SEDs_DES')},)
# sedLib.run(typ="G")