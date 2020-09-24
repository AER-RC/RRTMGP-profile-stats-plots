from matplotlib import pyplot as plt
import pandas as pd
import netCDF4
import numpy as np
import sys

wt = [
    0.1527534276, 0.1491729617, 0.1420961469, 0.1316886544,
    0.1181945205, 0.1019300893, 0.0832767040, 0.0626720116,
    0.0424925   , 0.0046269894, 0.0038279891, 0.0030260086,
    0.0022199750, 0.0014140010, 0.000533    , 0.000075
    ]

fname_lbl = 'TAUSORT_LBL_RRTMGP.csv'
lbl = pd.read_csv(fname_lbl, header=2)

fname_rrtmgp = 'rrtmgp-inputs-outputs.nc'
nc = netCDF4.Dataset(fname_rrtmgp, 'r')

m_tau = nc.variables['tau_minor']
datatmp = nc.variables['tau_minor'][9][:]
data1 = np.empty([42,224,1])
for j in range(224):
    for k in range(42):
        data1[k,j,0] = datatmp[j,k,0]

data1 = data1[0:,80:96].flatten()

# Calculate integrated OD
lbl['INTEG-OD-RRTMGP'] = np.zeros([672])
lbl['AVE-OD-RRTMGP'] = data1

# Calculate integrated OD
ktmp = 0
for j in range(42):
    b = 0.
    for k in range(16):
        ktmp = j*16 + k
        b = b + wt[k] * lbl.loc[ ktmp, ['AVE-OD-RRTMGP']].values
        lbl.loc[ktmp, ['INTEG-OD-RRTMGP']] = b

fname_lbl = 'compare_taus_self_atm17.csv'
lbl.rename(columns = {'INTEG-OD': 'INTEG-OD-LBL', 'AVE-OD': 'AVE-OD-LBL'}, inplace=True)
lbl['INTEG-OD-(RRTMGP-LBL)'] = lbl['INTEG-OD-RRTMGP'] - lbl['INTEG-OD-LBL']
lbl['AVE-OD-(RRTMGP-LBL)'] = lbl['AVE-OD-RRTMGP'] - lbl['AVE-OD-LBL']
headers = ['LAYER', 'GPT', 'TAVE', 'PAVE',
           'AVE-OD-RRTMGP', 'AVE-OD-LBL', 'AVE-OD-(RRTMGP-LBL)',
           'INTEG-OD-RRTMGP', 'INTEG-OD-LBL',  'INTEG-OD-(RRTMGP-LBL)']
lbl.to_csv(fname_lbl, columns = headers, index_label = 'INDEX')

