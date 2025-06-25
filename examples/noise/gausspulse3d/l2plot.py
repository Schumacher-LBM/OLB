import matplotlib.pyplot as plt
import numpy as np
from glob import glob

plt.rcParams['text.usetex'] = True

fig, ax = plt.subplots(figsize=(6,3), dpi=300)
filename = glob("tmp_eternal*/gnuplotData/data/l2_absolute.dat")[0]
data = np.loadtxt(filename)
ax.plot(data[:,0], data[:,1], label="ideal case (larger domain)")
filename = glob("tmp_periodic*/gnuplotData/data/l2_absolute.dat")[0]
data = np.loadtxt(filename)
ax.plot(data[:,0], data[:,1], label="periodic BC")
filename = glob("tmp_local*/gnuplotData/data/l2_absolute.dat")[0]
data = np.loadtxt(filename)
ax.plot(data[:,0], data[:,1], label="local BCs")
filename = glob("tmp_damping*/gnuplotData/data/l2_absolute.dat")[0]
data = np.loadtxt(filename)
ax.plot(data[:,0], data[:,1], label=f"sponge layer")
ax.set_yscale('log')
ax.set_ylim((1e-3,1.1))
ax.set_ylabel(r"$\frac{L_p}{L_{p0}}$", rotation=0, ha='right', fontsize=15)
ax.set_xlabel(r"$t_\mathrm{LU}$", fontsize=11, va='top')
ax.legend()
plt.tight_layout()
plt.show()
fig.savefig("l2ComparisonPlot.png")
