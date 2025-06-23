import os
import lephare as lp
import numpy as np
from matplotlib import pylab as plt

##### Create filter library #####
base_dir = os.path.abspath(os.path.join(os.getcwd(), '..')) #change to your lephare base_dir
filter_path = os.path.join(base_dir, 'lephare/training_stats/simulation_catalogs/buzzard_base/FILTERS') 
print(filter_path)
filter_rep = lp.keyword("FILTER_REP", filter_path)
filter_list = lp.keyword("FILTER_LIST",
"BuzzardLSSTu.res,BuzzardLSSTg.res,BuzzardLSSTr.res,\
     BuzzardLSSTi.res,BuzzardLSSTz.res,BuzzardLSSTy4.res")



filterLib = lp.Filter(config_keymap={"FILTER_REP":filter_rep,
                                      "FILTER_LIST":filter_list,
"TRANS_TYPE":lp.keyword("TRANS_TYPE", "0"),
"FILTER_CALIB":lp.keyword("FILTER_CALIB", "0"),
"FILTER_FILE":lp.keyword("FILTER_FILE", "photozdc1")
                                     })


filterLib.run()



##### Create SED library #####

SED_path = os.path.join(base_dir, 'lephare/training_stats/simulation_catalogs/buzzard_base/SEDS/updated_Buzzard_SEDs')

# SED_list = sorted([
#     f for f in os.listdir(SED_path)
#     if os.path.isfile(os.path.join(SED_path, f)) and f.endswith(".sed")
# ])

SED_list = [f'kmeansbuzzard_{i}.sed' for i in range(0,100)]

#save to file.list because we must do it like this
SED_list_path = os.path.join(SED_path, "updated_Buzzard_SEDs.list")
with open(SED_list_path, "w") as f:
    for sed_file in SED_list:
        f.write("Buzzard/"+sed_file + "\n")

SED_rep = lp.keyword("SED_REP", SED_path)
print(SED_rep)

###add SED to cahe file because ...###
import shutil

#Lephare cache file
LEPHARE_SED_GAL_PATH = os.path.expanduser("~/.cache/lephare/data/sed/GAL/Buzzard")

#Create cache file in case
os.makedirs(LEPHARE_SED_GAL_PATH, exist_ok=True)

# Cpypast all file to cache
for sed_file in SED_list + ["updated_Buzzard_SEDs.list"]:
    src = os.path.join(SED_path, sed_file)
    dst = os.path.join(LEPHARE_SED_GAL_PATH, sed_file)
    shutil.copyfile(src, dst)

print(f"past {len(SED_list)} file.sed and file.list in {LEPHARE_SED_GAL_PATH}")



SED_list_path = f"{base_dir}/lephare/training_stats/simulation_catalogs/buzzard_base/SEDS/updated_Buzzard_SEDs/updated_Buzzard_SEDs.list"
sedLib = lp.Sedtolib(config_keymap={
    "SED_REP": SED_rep,
    "GAL_SED": lp.keyword("GAL_SED", SED_list_path),
    "GAL_FSCALE": lp.keyword("GAL_FSCALE", "1.0"),
    "GAL_LIB": lp.keyword("GAL_LIB", 'buzzard_SEDs')},)
sedLib.run(typ="G")
