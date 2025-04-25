import os
import matplotlib.pyplot as plt
import numpy as np

base_dir = '/home/hallouin/Documents/thall_2025/lephare/simulation_catalogs/buzzard_base/SEDS'
SED_folder = os.path.join(base_dir, 'updated_Buzzard_SEDs')

plt.figure(figsize=(10, 6))

for SED in os.listdir(SED_folder):
    if not SED.endswith("9.sed"):  # Adjust if your SED files use a different extension
        continue
    file = os.path.join(SED_folder, f'{SED}')
    print(file)
    df = np.loadtxt(file)
    plt.semilogy(df.T[1], label=f'{SED}')

plt.title("SEDs")
plt.xlabel("Wavelength (A)")
plt.ylabel("Transmission")
plt.legend()
plt.show()