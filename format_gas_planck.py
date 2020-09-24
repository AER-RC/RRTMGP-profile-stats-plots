from matplotlib import pyplot as plt
import pandas as pd
import netCDF4
import numpy as np
import sys
import csv

def get_gpt_bnd(b, nc):
    # Use "REAL" band number, not indexed version
    # key_species(bnd, atmos_layer, pair)
    bnd = len(nc.dimensions['bnd'])
    gpt = len(nc.dimensions['gpt'])
    gpt_bnd = gpt/bnd
    g1 = (b-1)*gpt_bnd
    g2 = (b-1)*gpt_bnd + (gpt_bnd-1)

    return g1, g2

# Using python indices, which flavor for this band?
iflav = 1
icol = 0
bndnum = 12 # This is the REAL band number, not the python index



wt = [
    0.1527534276, 0.1491729617, 0.1420961469, 0.1316886544,
    0.1181945205, 0.1019300893, 0.0832767040, 0.0626720116,
    0.0424925   , 0.0046269894, 0.0038279891, 0.0030260086,
    0.0022199750, 0.0014140010, 0.000533    , 0.000075
    ]

fname_rrtmgp = 'coefficients.nc'
coeffs = netCDF4.Dataset(fname_rrtmgp, 'r')

fname_rrtmgp = 'gas_components.nc'
gas = netCDF4.Dataset(fname_rrtmgp, 'r')

fname_rrtmgp = 'planck_components.nc'
planck = netCDF4.Dataset(fname_rrtmgp, 'r')

g1, g2 = get_gpt_bnd(bndnum, coeffs)
g2 = g2 + 1
# Using python indices, which flavor for this band?
iflav = 1
icol = 0
# fname_rrtmgp = 'planck_components.nc'
# planck = netCDF4.Dataset(fname_rrtmgp, 'r')

ngpt = len(gas.dimensions['var1_dim1'])
nlay = len(gas.dimensions['var1_dim2'])
ncol = len(gas.dimensions['var1_dim3'])
nflavor = len(gas.dimensions['var13_dim1'])

p_lay = gas.variables['p_lay']
t_lay = gas.variables['t_lay']
interpolated_eta = gas.variables['interpolated_eta']
etaval = gas.variables['etaval']

fname_out = 'rrtmgp_components_preind-key-minor-quad_co2_42_12.csv'
header = ['gpt', 'lay', 'p_lay', 't_lay',
          'tau_major', 'integ_tau_major',
          'pf', 'integ_pf',
          'eta_interpolated',
          'f_press', 'j_press', 'f_temp', 'j_temp',
          'etaval_1', 'etaval_2', 'f_etaval_1', 'f_etaval_2',
          'j_eta_1', 'j_eta_2',
          'w_comb_1', 'w_comb_2',
          'fmajor_1', 'w_kmajor_1', 'pf_1',
          'fmajor_2', 'w_kmajor_2', 'pf_2',
          'fmajor_3', 'w_kmajor_3', 'pf_3',
          'fmajor_4', 'w_kmajor_4', 'pf_4',
          'fmajor_5', 'w_kmajor_5', 'pf_5',
          'fmajor_6', 'w_kmajor_6', 'pf_6',
          'fmajor_7', 'w_kmajor_7', 'pf_7',
          'fmajor_8', 'w_kmajor_8', 'pf_8'
          ]

with open (fname_out, 'wb') as csvFile:
    outputwriter = csv.writer(csvFile, delimiter=',')
    outputwriter.writerow(header)
    # Account for fact that output in gas_components is in Fortran indexing (start: 1),
    # subtract 1 for python
    for ilay in range(nlay):
        j_press_1 = gas.variables['jpress_for_k'][ilay][icol] - 1 - 1
        j_press_2 = gas.variables['jpress_for_k'][ilay][icol] - 1
        j_temp_1 = gas.variables['j_temp'][ilay][icol] - 1
        j_temp_2 = gas.variables['j_temp'][ilay][icol] + 1 - 1
        j_eta_1 = gas.variables['j_eta'][ilay][icol][iflav][0] - 1
        j_eta_2 = gas.variables['j_eta'][ilay][icol][iflav][1] - 1
        col_mix_1 = gas.variables['col_mix'][ilay][icol][iflav][0]
        col_mix_2 = gas.variables['col_mix'][ilay][icol][iflav][1]
        tau_integ = 0
        pf_integ = 0
        ict = 0
        for igpt in range(g1, g2, 1):
            tau_integ = tau_integ + gas.variables['major_tau'][icol][ilay][igpt]
            pf_integ = pf_integ + planck.variables['plank_fraction'][icol][ilay][igpt]
            ict = ict + 1
            outputwriter.writerow(
                [
                    igpt + 1, ilay + 1,
                    gas.variables['p_lay'][ilay][0],
                    gas.variables['t_lay'][ilay][0],
                    gas.variables['major_tau'][icol][ilay][igpt],
                    tau_integ,
                    planck.variables['plank_fraction'][icol][ilay][igpt],
                    pf_integ,
                    gas.variables['interpolated_eta'][ilay][icol][iflav],
                    gas.variables['f_press'][ilay][icol],
                    gas.variables['j_press'][ilay][icol],
                    gas.variables['f_temp'][ilay][icol],
                    j_temp_1 + 1,
                    gas.variables['etaval'][ilay][icol][iflav][0],
                    gas.variables['etaval'][ilay][icol][iflav][1],
                    gas.variables['f_etaval'][ilay][icol][iflav][0],
                    gas.variables['f_etaval'][ilay][icol][iflav][1],
                    j_eta_1 + 1, j_eta_2 + 1,
                    col_mix_1, col_mix_2,
                    gas.variables['f_major'][ilay][icol][iflav][0][0][0],
                    col_mix_1 * coeffs.variables['kmajor'][j_temp_1][j_press_1][j_eta_1][igpt],
                    coeffs.variables['plank_fraction'][j_temp_1][j_press_1][j_eta_1][igpt],
                    gas.variables['f_major'][ilay][icol][iflav][0][0][1],
                    col_mix_1 * coeffs.variables['kmajor'][j_temp_1][j_press_1][j_eta_1 + 1][igpt],
                    coeffs.variables['plank_fraction'][j_temp_1][j_press_1][j_eta_1 + 1][igpt],
                    gas.variables['f_major'][ilay][icol][iflav][0][1][0],
                    col_mix_1 * coeffs.variables['kmajor'][j_temp_1][j_press_2][j_eta_1][igpt],
                    coeffs.variables['plank_fraction'][j_temp_1][j_press_2][j_eta_1][igpt],
                    gas.variables['f_major'][ilay][icol][iflav][0][1][1],
                    col_mix_1 * coeffs.variables['kmajor'][j_temp_1][j_press_2][j_eta_1 + 1][igpt],
                    coeffs.variables['plank_fraction'][j_temp_1][j_press_2][j_eta_1 + 1][igpt],
                    gas.variables['f_major'][ilay][icol][iflav][1][0][0],
                    col_mix_2 * coeffs.variables['kmajor'][j_temp_1][j_press_1][j_eta_2][igpt],
                    coeffs.variables['plank_fraction'][j_temp_1][j_press_1][j_eta_2][igpt],
                    gas.variables['f_major'][ilay][icol][iflav][1][0][1],
                    col_mix_2 * coeffs.variables['kmajor'][j_temp_1][j_press_1][j_eta_2 + 1][igpt],
                    coeffs.variables['plank_fraction'][j_temp_1][j_press_1][j_eta_2 + 1][igpt],
                    gas.variables['f_major'][ilay][icol][iflav][1][1][0],
                    col_mix_2 * coeffs.variables['kmajor'][j_temp_1][j_press_2][j_eta_2][igpt],
                    coeffs.variables['plank_fraction'][j_temp_1][j_press_2][j_eta_2][igpt],
                    gas.variables['f_major'][ilay][icol][iflav][1][1][1],
                    col_mix_2 * coeffs.variables['kmajor'][j_temp_1][j_press_2][j_eta_2 + 1][igpt],
                    coeffs.variables['plank_fraction'][j_temp_1][j_press_2][j_eta_2 + 1][igpt],
                    ]
            )

csvFile.close()
gas.close()
coeffs.close()
