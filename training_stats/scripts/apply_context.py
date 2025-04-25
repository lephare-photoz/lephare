### Tests on rm_discerpant function

import numpy as np
import os

###load .dat
base_dir = '/home/hallouin/Documents/thall_2025/lephare/'
file_dat_path = os.path.join(base_dir, 'simulation_catalogs/buzzard_base/Final_Buzzard_training_file.dat')
badmags = np.loadtxt(file_dat_path)


corr_badmags = np.copy(badmags)
source_indices = np.arange(len(corr_badmags), dtype='int')

#loop over all sources
for k, source in zip(source_indices,corr_badmags):
    mag_source = source[1:7]
    badmag_count = len(mag_source)
    indices = np.arange(len(mag_source), dtype='int')
    #loop over filters
    for i, mag in zip(indices, mag_source):
        if mag < 98:
            badmag_count -= 1
            #apply context
            corr_badmags[k][13] += 2**i
    # print(corr_badmags[k][13])

CAT_IN_PATH = os.path.join(base_dir, 'simulation_catalogs/buzzard_base/context_Final_Buzzard_training_file.dat')
formats = ['%d'] + ['%.5f'] * 12 + ['%d'] + ['%.5f']
np.savetxt(CAT_IN_PATH, corr_badmags, fmt=formats)




# ###
# CAT_IN_PATH = os.path.join(base_dir, 'simulation_catalogs/buzzard_base/context_aberrantmags_Buzzard_training_file.dat')
# formats = ['%d'] + ['%.5f'] * 12 + ['%d'] + ['%.5f']
# np.savetxt(CAT_IN_PATH, file_dat_shuffled, fmt=formats)