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

fname_rrtmgp = 'gas_components.nc'
nc = netCDF4.Dataset(fname_rrtmgp, 'r')

fname_rrtmgp = 'planck_components.nc'
pl = netCDF4.Dataset(fname_rrtmgp, 'r')
pl_data = pl.variables['plank_fraction'][:][0]
pl_data = pl_data[0:, 176:192].flatten()

# m_tau = nc.variables['minor_tau']
# data1 = nc.variables['minor_tau'][:][0]
# data1 = data1[0:, 0:16, 8].flatten()
m_tau = nc.variables['major_tau']
data1 = nc.variables['major_tau'][:][0]
data1 = data1[0:, 176:192].flatten()

# Calculate integrated OD

#layer = np.zeros([672])
# gpt = np.zeros([672])
# for i in range(42):
#     #print i, type(layer), layer.shape
#     layer[i*16:i*16+16] = i + 1
#     #print layer[i*16:i*16+16]
#     gpt[i*16:i*16+16] = range(1, 17, 1)
#
# layer = layer.astype('int')
# gpt = gpt.astype('int')

#index = range(672)
#df = pd.DataFrame(index = index, columns = ['m_tau_rrtmgp'])
# df['INTEG-OD-LBL'] = data1
# df['LAYER'] = layer
# df['IG'] = gpt

lbl['INTEG-OD-RRTMGP'] = np.zeros([672])
lbl['AVE-OD-RRTMGP'] = data1

lbl['INTEG-PL-RRTMGP'] = np.zeros([672])
lbl['AVE-PL-RRTMGP'] = pl_data

# Calculate integrated OD
ktmp = 0
for j in range(42):
    b = 0.
    c = 0.
    for k in range(16):
        ktmp = j*16 + k
        b = b + wt[k] * lbl.loc[ ktmp, ['AVE-OD-RRTMGP']].values
        lbl.loc[ktmp, ['INTEG-OD-RRTMGP']] = b
        c = c + lbl.loc[ ktmp, ['AVE-PL-RRTMGP']].values
        lbl.loc[ktmp, ['INTEG-PL-RRTMGP']] = c

#lbl['AVE-OD-LBL-RRTMGP'] = lbl['AVE-OD-LBL'] - lbl['AVE-OD']
fname_lbl = 'compare_taus_preind-key-minor-quad_co2_42_12_r1464.csv'
lbl.rename(columns = {'INTEG-OD': 'INTEG-OD-LBL', 'AVE-OD': 'AVE-OD-LBL'}, inplace=True)
lbl['INTEG-OD-(RRTMGP-LBL)'] = lbl['INTEG-OD-RRTMGP'] - lbl['INTEG-OD-LBL']
lbl['AVE-OD-(RRTMGP-LBL)'] = lbl['AVE-OD-RRTMGP'] - lbl['AVE-OD-LBL']
headers = ['LAYER', 'GPT', 'TAVE', 'PAVE',
           'AVE-OD-RRTMGP', 'AVE-OD-LBL', 'AVE-OD-(RRTMGP-LBL)',
           'INTEG-OD-RRTMGP', 'INTEG-OD-LBL',  'INTEG-OD-(RRTMGP-LBL)']
lbl.to_csv(fname_lbl, columns = headers, index_label = 'INDEX')

fname_lbl = 'compare_planck_preind-key-minor-quad_co2_42_12_r1464.csv'
lbl.rename(columns = {'INTEG-PL': 'INTEG-PL-LBL', 'AVE-PL': 'AVE-PL-LBL'}, inplace=True)
lbl['INTEG-PL-(RRTMGP-LBL)'] = lbl['INTEG-PL-RRTMGP'] - lbl['INTEG-PL-LBL']
lbl['AVE-PL-(RRTMGP-LBL)'] = lbl['AVE-PL-RRTMGP'] - lbl['AVE-PL-LBL']
headers = ['LAYER', 'GPT', 'TAVE', 'PAVE',
           'AVE-PL-RRTMGP', 'AVE-PL-LBL', 'AVE-PL-(RRTMGP-LBL)',
           'INTEG-PL-RRTMGP', 'INTEG-PL-LBL',  'INTEG-PL-(RRTMGP-LBL)']
lbl.to_csv(fname_lbl, columns = headers, index_label = 'INDEX')

# lbl.PLot(y='AVE-OD-(RRTMGP-LBL)',x='AVE-OD-LBL', kind='scatter')
# PLt.show()
