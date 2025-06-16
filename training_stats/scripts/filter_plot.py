import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib 
matplotlib.rc('ytick', labelsize=14)
matplotlib.rc('xtick', labelsize=14)
matplotlib.rcParams.update({'font.size': 18})
# plt.rcParams["font.family"] = "DejaVu Sans"


base_dir = os.path.abspath(os.path.join(os.getcwd(), '..'))
filter_folder = '/home/hallouin/Documents/thall_2025/photoz/lephare/training_stats/simulation_catalogs/buzzard_base/FILTERS'
print()
plt.figure(figsize=(10, 6))

filters = ['BuzzardLSSTu.res', 'BuzzardLSSTg.res', 'BuzzardLSSTr.res', 'BuzzardLSSTi.res', 'BuzzardLSSTz.res', 'BuzzardLSSTy4.res',]

cmap = cm.get_cmap('turbo')
colors = [cmap(i / (len(filters) - 1)) for i in range(len(filters))]

for i, filt in enumerate(filters):
    if not filt.endswith(".res"):  # Adjust if your filter files use a different extension
        continue
    filtre = os.path.join(filter_folder, f'{filt}')
    df = np.loadtxt(filtre)
    plt.plot(df.T[0], df.T[1], label=f'{filt}'[11], color=colors[i])

# plt.title("Filters transmission as a function of wavelength")
plt.xlabel("Wavelength (A)")
plt.ylabel("Transmission (Normed at maximum)")
plt.xlim(2500,12000)
plt.legend()
plt.show()