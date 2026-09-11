import numpy as np
import matplotlib.pyplot as plt
import sys
import os


# Check if filename is provided
if len(sys.argv) < 2:
    print("Usage: python3 plot_bands.py <datafile1> <datafile2> ... <datafileN>")
    sys.exit(1)

# plot energy window [-EnergyWindow, EnergyWindow] (EnergyWindow in eV)
EnergyWindow = 0.1

# OffsetZero = 0: plot energies as in output file. OffsetZero = 1: fix Dirac point (if appropriate) at zero. 
OffsetZero = 0

plt.figure(figsize=(6, 8))


for nn in range(1,len(sys.argv)):
    filename = sys.argv[nn]

    # Load the data
    data = np.loadtxt(filename)

    kpoints = data[:, 0]
    bands = data[:, 1:]

    # Plot with automatic color cycling
    for i in range(bands.shape[1]):
        plt.plot(kpoints, bands[:, i] - OffsetZero*bands[0,bands.shape[1]//2], linewidth=1,marker='o',markersize=2, c='C' + str(nn))


plt.plot(kpoints,0.0*kpoints,linewidth=1)

plt.ylabel("Energy (eV)", fontsize=14)
plt.title(os.path.basename(filename), fontsize=14)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

# Custom ticks
xnames = ["$K'_M$", "$\\Gamma_M$", "$K_M$", "$M_M$", "$\\Gamma_M$", "$K_M$"]
xticklist = [0., 0.22904127, 0.45808254, 0.57260317, 0.77095873, 1.]
plt.xticks(xticklist, xnames)


plt.xlim(1, 0.45808254)     # force x-axis exactly from 0 to 1
plt.ylim(-EnergyWindow, EnergyWindow)   

#plt.tight_layout()
plt.show()

