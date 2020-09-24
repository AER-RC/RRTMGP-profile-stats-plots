from netCDF4 import Dataset
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
import sys

gwt = [
    0.1527534276, 0.1491729617, 0.1420961469, 0.1316886544,
    0.1181945205, 0.1019300893, 0.0832767040, 0.0626720116,
    0.0424925   , 0.0046269894, 0.0038279891, 0.0030260086,
    0.0022199750, 0.0014140010, 0.000533    , 0.000075
    ]

gas = Dataset('gas_components.nc','r')
# major_tau(ncols, nlays, ngpts)
ncols = gas.variables['major_tau'].shape[0]
nlays = gas.variables['major_tau'].shape[1]
ngpts = gas.variables['major_tau'].shape[2]
tau_major = gas.variables['major_tau']
# col_gas(0:ngas, nlays, ncols)
col_gas = gas.variables['col_gas']
# p/t_lay(nlays, ncols)
t_lay = gas.variables['t_lay']
p_lay = gas.variables['p_lay']
etainterp = gas.variables['interpolated_eta']

wt_od = np.ndarray([ncols, nlays])
abs_co = np.ndarray([ncols, nlays])
# water
igas = 1
tmp_wt = 0.0
for icol in range(ncols):
    for ilay in range(nlays):
        tmp_wt = 0.0
        for igpt in range(16):
            tmp_wt = tmp_wt + tau_major[icol, ilay, igpt] * gwt[igpt]
            wt_od[icol, ilay] = np.sum(tau_major[icol, ilay, 0:16] * gwt[0:16])
            abs_co[icol, ilay] = wt_od[icol, ilay]/col_gas[igas, ilay, icol]
        #     print 'HERE', icol, ilay, igpt, tmp_wt, wt_od[icol, ilay]
        # sys.exit()

# Write to CSV file
columns = ['Column', 'Layer', 'Pressure', 'Temperature', 'Col_Gas(H2O)', 'Wt_OD', 'Abs_Coefficient']
f = open('garand_absco.txt', 'wt')

writer = csv.writer(f)
n = 0
for icol in range(ncols):
    for ilay in range(nlays):
        writer.writerow((n, icol, ilay, p_lay[icol, ilay], t_lay[icol, ilay],
                         col_gas[1, ilay, icol]))
        n = n + 1

# Plot T v. Abs Co
npts = nlays * ncols
t = np.ndarray([npts])
k = np.ndarray([npts])
p = np.ndarray([npts])
n = 0
for icol in range(ncols):
    for ilay in range(nlays):
        t[n] = t_lay[ilay, icol]
        p[n] = p_lay[ilay, icol]
        k[n] = abs_co[ilay, icol]
        n = n + 1
        print 'ETAINTERP', icol, ilay, etainterp[ilay, icol, 0:10]
    sys.exit()

plt.scatter(t, k)
plt.xlabel('Layer Temperature [K]')
plt.ylabel('Absorption Coefficient')
plt.ylim(min(k), max(k))
plt.savefig('garand_t_absco.png')
