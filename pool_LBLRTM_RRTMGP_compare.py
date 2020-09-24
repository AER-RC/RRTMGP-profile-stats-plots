#!/usr/bin/env python

# on Mac OS, this was preventing this library from running:
# mpl.rcParams['backend'] = MacOSX
import matplotlib as mpl
mpl.use('TkAgg')

import os, sys, argparse
import numpy as np
import matplotlib.pyplot as plot
import ConfigParser

sys.path.append('/home/rpernak/revised_python_libraries')
import utils

sys.path.append('/rd47/scratch/RRTMGP/pyscripts')
import LBLRTM_RRTMGP_compare as rrtmgp

# for multi-page PDF files
from matplotlib.backends.backend_pdf import PdfPages

# for RP local environment
import socket
host = socket.gethostname()
if host == 'rd47':
  ncDir = '%s/band_generation/model_output' % os.getenv('RRTMGP')
  ncRef = '%s/lblrtm-lw-flux-inputs-outputs-garand-all.nc' % ncDir
  ncTest = '%s/rrtmgp-lw-flux-inputs-outputs-garand-all.nc' % ncDir
else:
  ncRef = 'lblrtm-lw-flux-inputs-outputs-garand-all.nc'
  ncTest = 'rrtmgp-lw-flux-inputs-outputs-garand-all.nc'
# end host

# trial and error spacing for subplots
hspace = 0.5

def prepForPool(refFile, testFile, outDir='.', shortWave=False, \
  tPauseP=100.0, plotMean=False, yLog=False, inBand=None, \
  atmType='Garand Atmospheres', \
  prefix='PROFS_sens_key_cnt_add_min_04_dbl_r472_trynewn2opres'):

  """
  Call
    broadbandList = prepForPool(ref, test)

  Input
    refFile -- string, path to netCDF for reference model
    testFile -- string, path to netCDF for test model

  Output
    List of dictionaries that will be used in poolProfPDFs() -- each 
    list element is a dictionary for a single band

  Keywords
    outDir -- string, top level directory for output figures
    shortWave -- boolean, specifies that short-wave radiation flux 
      units should be specified in axis labels
    tPauseP -- float, tropopause pressure threshold (mbar)
    prefix -- string that will be placed in front of the band  
      number in the name of the output PDF file
    plotMean -- boolean, plot mean of variable for a given pressure 
      and band over all columns
    yLog -- boolean, plot ordinate axis on log scale instead of linear
    atmType -- string, atmosphere type for profiles
  """

  plotVars = ['band_flux_up', 'band_flux_dn', 'band_heating_rate', \
    'p_lay', 'p_lev', 'band_lims_wvn']

  refDict = rrtmgp.getVars(refFile, attrList=plotVars)
  testDict = rrtmgp.getVars(testFile, attrList=plotVars)

  # some quality control (consistency check)
  dum = plotVars[0]
  if refDict[dum].shape != testDict[dum].shape:
    sys.exit('%s and %s do not have equal dimensions, returning' % \
      refFile, testFile)

  # grab dimensions of variables
  varShape = refDict[dum].shape
  nBand = varShape[2]; nCol = varShape[1]

  broadbandList = []

  bandList = range(nBand) if inBand is None else [inBand]
  refWN = refDict['band_lims_wvn']
  for band in bandList:
    outFile = '%s/%s_%02d.pdf' % (outDir, prefix, band+1)
    wnRange = refWN[:, band]
    bandDict = {'plotVars': plotVars, 'wnRange': wnRange, \
      'pdf': outFile, 'band': band, 'numberCol': nCol, \
      'sw': shortWave, 'pLay': refDict['p_lay'], \
      'pLev': refDict['p_lev'],  \
      'band_flux_up_ref': refDict['band_flux_up'][:, :, band], \
      'band_flux_dn_ref': refDict['band_flux_dn'][:, :, band], \
      'band_heating_rate_ref': \
        refDict['band_heating_rate'][:, :, band], \
      'band_flux_up_test': testDict['band_flux_up'][:, :, band], \
      'band_flux_dn_test': testDict['band_flux_dn'][:, :, band], \
      'band_heating_rate_test': \
        testDict['band_heating_rate'][:, :, band], \
      'tPauseP': tPauseP, 'plot_mean': plotMean, 'log_y': yLog, \
      'atmType': atmType, 'doBroadband': False}
    broadbandList.append(bandDict)
  # end bands

  # for summing arrays:
  axNum = 2

  # and now the broadband dictionary
  outFile = '%s/%s_broadband.pdf' % (outDir, prefix)
  bandDict = {'plotVars': plotVars, \
    'wnRange': [refWN[0, 0], refWN[1, -1]], \
    'pdf': outFile, 'band': 'Broad', 'numberCol': nCol, \
    'sw': shortWave, 'pLay': refDict['p_lay'], \
    'pLev': refDict['p_lev'],  \
    'band_flux_up_ref': np.sum(refDict['band_flux_up'], axis=axNum), \
    'band_flux_dn_ref': np.sum(refDict['band_flux_dn'], axis=axNum), \
    'band_heating_rate_ref': \
      np.sum(refDict['band_heating_rate'], axis=axNum), \
    'band_flux_up_test': \
      np.sum(testDict['band_flux_up'], axis=axNum), \
    'band_flux_dn_test': \
       np.sum(testDict['band_flux_dn'], axis=axNum), \
    'band_heating_rate_test': \
       np.sum(testDict['band_heating_rate'], axis=axNum), \
    'tPauseP': tPauseP, 'plot_mean': plotMean, 'log_y': yLog, \
    'atmType': atmType, 'doBroadband': True}

  broadbandList.append(bandDict)
  return broadbandList

# end prepForPool()

def poolProfPDFs(paramDict):
  """
  Plot upward flux, downward flux, and heating rate as well as the 
  test-ref difference for each of those three variables on a single 
  page (3x2 array of figures) for each column and a single band 
  specified in ref and test. Fluxes and HR are plotted as functions 
  of pressure, so we are plotting profiles. This should be used with 
  the map() method in the multiprocessing library (Pool object)

  Call
    from multiprocessing import Pool, cpu_count, Process
    p = Pool(nCores)
    compareDat = prepForPool(ref, test)
    p.map(poolProfPDFs, compareDat)

  Input
    paramDict -- dictionary, output from prepForPool

  Output
    PDF for a single waveband as specified in paramDict

  """

  # unpack dictionary
  plotVars = paramDict['plotVars']
  outFile = paramDict['pdf']
  wnRange = paramDict['wnRange']
  band = paramDict['band']
  nCol = paramDict['numberCol']
  shortWave = paramDict['sw']
  atmType = paramDict['atmType']
  doBroad = paramDict['doBroadband']

  plotTitle = ['Upward Flux', 'Downward Flux', 'Heating Rate', \
    'Pressure (mbar)', 'Pressure (mbar)', 'Wavenumber Range']
  pdf = PdfPages(outFile)

  for aCol in range(nCol):
    # plot 1 atm column per page, 3x2 array of subfigures
    print 'Column %d of %d' % (aCol+1, nCol)

    figTitle = '%s Column %d, %.0f-%.0f cm$^{-1}$' % \
      (atmType, aCol+1, wnRange[0], wnRange[1])
    fig = plot.figure()
    fig.set_size_inches(8.5, 11)
    fig.suptitle(figTitle, fontweight='bold')

    ctr = 1

    # space between plots
    plot.subplots_adjust(hspace=hspace, wspace=hspace)

    # now loop over subplot rows and columns
    for iRow in range(3):
      pVar = plotVars[iRow]
      pTitle = plotTitle[iRow]
      for iCol in range(2):
        plot.subplot(3, 2, ctr)

        # fluxes are by level, heat rate are by layer
        pStr = 'pLay' if iRow == 2 else 'pLev'
        pressure = paramDict[pStr]

        # some other plotting params
        delta = True if iCol % 2 == 1 else False
        if 'Flux' in pTitle:
          units = '[W m$^{-2}$]' if shortWave else '[W m$^{-2}$]'
          tempVar = 'Flux'
        else:
          units = '[K day$^{-1}$]'
          tempVar = 'HR'
        # end units and tempVar

        if delta:
          xt = 'RRTMGP-LBLRTM %s Difference %s' % (tempVar, units)
          yt = ''
        else:
          xt = '%s %s' % (tempVar, units)
          yt = 'Pressure [mbar]'
        # end delta

        varRef = paramDict[pVar + '_ref']
        varTest = paramDict[pVar + '_test']
        rrtmgp.plotProfiles(varRef[:, aCol], varTest[:, aCol], \
          pressure[:, aCol], pTitle=pTitle, xTitle=xt, yTitle=yt, \
          plotDelta=delta, tPauseP=paramDict['tPauseP'], \
          plotMean=paramDict['plot_mean'], yLog=paramDict['log_y'])
        ctr += 1
      # end iCol
    # end iRow

    # write page for column
    pdf.savefig()
    plot.close()
  # end column

  pdf.close()
  bandStr = 'broad' if doBroad else '%d' % (band+1) 
  print 'Processed band %s' % bandStr
  return
# end poolProfPDFs()

def poolStatPDFs(paramDict):
  """
  ...This should be used with 
  the map() method in the multiprocessing library (Pool object)

  Call
    from multiprocessing import Pool, cpu_count, Process
    p = Pool(nCores)
    compareDat = prepForPool(ref, test)
    p.map(poolStatPDFs, compareDat)

  Input
    paramDict -- dictionary, output from prepForPool

  Output
    PDF...

  """

  # unpack dictionary
  plotVars = paramDict['plotVars']
  outFile = paramDict['pdf']
  wnRange = paramDict['wnRange']
  band = paramDict['band']
  nCol = paramDict['numberCol']
  shortWave = paramDict['sw']

  sys.exit('This needs a lot of work')
  rrtmgp.statPDF(refFile, testFile, singlePDF=True, \
    tPauseP=paramDict['tPauseP'], xTitle=xt, yTitle=yt, \
    prefix=statPrefix, atmType=aType, statCSV=statCSV)

  return
# end poolProfPDFs()

if __name__ == '__main__':

  parser = argparse.ArgumentParser(\
    description='Attempt to multithread the plotting process ' + \
    'as designed in LBLRTM_RRTMGP_compare.py.')
  parser.add_argument('--reference_file', type=str, default=ncRef, \
    help='Full path to netCDF file that contains reference ' + \
    'model output (probably LBLRTM output).')
  parser.add_argument('--test_file', type=str, default=ncTest, \
    help='Full path to netCDF file that contains test model ' + \
    'output (probably RRTMGP output).')
  parser.add_argument('--cores', type=int, default=4, \
    help='Number of cores to use.')
  parser.add_argument('--profiles', action='store_true', \
    help='Plot RRTM-LBLRTM profiles (flux and heating rate).')
  parser.add_argument('--stats', action='store_true', \
    help='Plot statistics for RRTM-LBLRTM differences (flux and ' + \
    'heating rate). This should be set to string to where the ' + \
    'combined PDF files will be saved.')
  parser.add_argument('--band', type=int, \
    help='Number of band to plot. Default (band=None) is all bands.')
  parser.add_argument('--log_y', action='store_true', \
    help='Generate a semilog-y plot.')
  parser.add_argument('--mean', action='store_true', \
    help='Plot column-averaged parameters in a separate page.')
  parser.add_argument('--atm_type', type=str, \
    default='Garand Atmospheres', \
    help='Atmospheric type for all of the profiles.')
  parser.add_argument('--config_file', type=str, \
    help='Path to configuration file that contains the values ' + \
    'to keyword arguments (so this method can be used instead of ' + \
    'providing the keywords; config_file data supercede any ' + \
    'keyword argument input).')
  args = parser.parse_args()

  from multiprocessing import Pool, cpu_count, Process

  inBand = None if args.band is None else args.band-1

  conFile = args.config_file
  if conFile:
    conFile = args.config_file
    utils.file_check(conFile)
    cParse = ConfigParser.ConfigParser()
    cParse.read(conFile)
    cRefName = cParse.get('Plot Params', 'reference_model')
    cRefMD = cParse.get('Plot Params', 'reference_description')
    cTestMD = cParse.get('Plot Params', 'test_description')
    cTestName = cParse.get('Plot Params', 'test_model')
    aType = cParse.get('Plot Params', 'atmosphere')

    refFile = cParse.get('Filename Params', 'reference_path')
    testFile = cParse.get('Filename Params', 'test_path')
    profPrefix = cParse.get('Filename Params', 'profiles_prefix')
    statPrefix = cParse.get('Filename Params', 'stats_prefix')

    xt = cRefName
    yt = '%s - %s' % (cTestName, cRefName)
  else:
    refFile = args.reference_file
    testFile = args.test_file
    xt = 'LBLRTM'
    yt = 'RRTMGP - LBLRTM'
    statPrefix = 'stats_lblrtm_rrtmgp'
    profPrefix = 'PROFS_sens_key_cnt_add_min_04_dbl_r472_trynewn2opres'
    fluxPrefix = 'flux_compare_LBLRTM_RRTMGP'
    aType = args.atm_type
  # end config_file

  # plot profile statistics
  if args.stats:
    rrtmgp.statPDF(refFile, testFile, forcing=False, singlePDF=True, \
      xTitle=xt, yTitle=yt)
  else:
    compareDat = prepForPool(refFile, testFile, plotMean=args.mean, \
      yLog=args.log_y, inBand=inBand, atmType=aType, \
      prefix=profPrefix)

    nCores = args.cores
    totCores = cpu_count()
    nCores = nCores if nCores < totCores else totCores-1

    p = Pool(nCores)

    if args.profiles: p.map(poolProfPDFs, compareDat)
  # end plot_stats
# end main()

