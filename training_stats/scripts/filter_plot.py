import os
import matplotlib.pyplot as plt
import numpy as np

base_dir = os.path.abspath(os.path.join(os.getcwd(), '..'))
filter_folder = os.path.join(base_dir, 'FILTERS')

plt.figure(figsize=(10, 6))

for filt in os.listdir(filter_folder):
    if not filt.endswith(".res"):  # Adjust if your filter files use a different extension
        continue
    filtre = os.path.join(filter_folder, f'{filt}')
    print(filtre)
    df = np.loadtxt(filtre)
    plt.plot(df.T[1], label=f'{filt}'[11])

plt.title("Filters")
plt.xlabel("Wavelength (A)")
plt.ylabel("Transmission")
plt.legend()
plt.show()