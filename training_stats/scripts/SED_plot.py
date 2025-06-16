import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib 
matplotlib.rc('ytick', labelsize=14)
matplotlib.rc('xtick', labelsize=14)
matplotlib.rcParams.update({'font.size': 18})

base_dir = '/home/hallouin/Documents/thall_2025/photoz/lephare/training_stats/simulation_catalogs/buzzard_base/SEDS'
SED_folder = os.path.join(base_dir, 'updated_Buzzard_SEDs')


plt.figure(figsize=(10, 6))

for i, SED in enumerate(os.listdir(SED_folder)):
    print(i)
    if i%19!=0:  # Adjust if your SED files use a different extension
        continue
    file = os.path.join(SED_folder, f'{SED}')
    print(file)
    df = np.loadtxt(file)
    plt.loglog(df.T[0], df.T[1], label=f'{SED}')

# plt.title("SEDs")
plt.xlabel("Wavelength (A)")
plt.axvspan(2500, 12000, color='gray', alpha=0.1, label='SDSS filter region')
# plt.xlim(2500,12000)
# plt.ylim(1e-3,12000)
plt.ylabel("Transmission (Normed at 800A)")
# plt.legend()
plt.show()