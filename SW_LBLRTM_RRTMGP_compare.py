#!/usr/bin/env python

import os, sys, argparse
import numpy as np

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

import matplotlib.pyplot as plot
from matplotlib import rc
import matplotlib.font_manager as font
import netCDF4 as nc
import datetime as DT
import ConfigParser
import time # for clocking processes

# for multi-page PDF files
from matplotlib.backends.backend_pdf import PdfPages

# should be in working dir
import utils
import LBLRTM_RRTMGP_compare as compare

# some plotting parameters
FUNITS = '[W m$^{-2}$]'; HRUNITS = '[K day$^{-1}$]'
FUPSTR = r'$F^{\uparrow}$'; HRSTR = r'$\dot{T}_{max}$'
DIRSTR = r'$F^{\downarrow}_{Direct}$'
DIFSTR = r'$F^{\downarrow}_{Diffuse}$'
font_prop = font.FontProperties()
font_prop.set_size(8)

# for flux to HR conversion, which i believe is described in
# RRTM and in SW_direct_beam_validations.txt
# https://rrtmgp.slack.com/files/jdelamere/F6NRG4PNX/sw_direct_beam_validations.txt
# flux needs to be in W/m2 and P in mbar
HEATFAC = 8.4391

def parseConfig(inFile):
  """
  Read configuration file input into main() and extract necessary 
  information from it

  Input
    inFile -- string, full path to input configuration file

  Output
    dictionary with the following keys:
      ref -- string, full path to reference model netCDF (LBLRTM)
      test -- string, full path to test model netCDF (RRTMGP)
      stat -- string, prefix for stats PDF file
      prof -- string, prefix for profiles PDF file
      x -- string, x-axis label
      y -- string, y-axis label
      atm -- string, atmosphere type (e.g., Garand)
      forcing -- boolean, are forcing files being compared?
  """

  # extract metadata from configuration file
  cParse = ConfigParser.ConfigParser()
  cParse.read(inFile)
  cRefName = cParse.get('Plot Params', 'reference_model')
  cRefMD = cParse.get('Plot Params', 'reference_description')
  cTestMD = cParse.get('Plot Params', 'test_description')
  cTestName = cParse.get('Plot Params', 'test_model')
  aType = cParse.get('Plot Params', 'atmosphere')

  # flux files
  refFile = cParse.get('Filename Params', 'reference_path')
  testFile = cParse.get('Filename Params', 'test_path')
  utils.file_check(refFile); utils.file_check(testFile)

  # forcing files
  rForceFile = \
    cParse.get('Filename Params', 'reference_force_path')
  tForceFile = cParse.get('Filename Params', 'test_force_path')

  if len(rForceFile) > 0 and len(tForceFile) > 0:
    forcing = True
  else:
    forcing = False

  if forcing:
    utils.file_check(rForceFile); utils.file_check(tForceFile)
    refFile, testFile = compare.forcingDiff(refFile, testFile, \
      rForceFile, tForceFile, \
      repVars=['flux_dir_dn', 'band_flux_dir_dn'])

    cRefForceName = \
      cParse.get('Plot Params', 'reference_forcing_model')
    cRefForceMD = cParse.get('Plot Params', 'reference_description')
    cTestForceName = cParse.get('Plot Params', 'test_forcing_model')
    cTestForceMD = cParse.get('Plot Params', 'test_description')

    xt = cRefForceName
    yt = '%s - %s' % (cTestForceName, cRefForceName)
  else:
    xt = cRefName
    yt = '%s - %s' % (cTestName, cRefName)
  # endif forcing

  statPrefix = cParse.get('Filename Params', 'stats_prefix')
  profPrefix = cParse.get('Filename Params', 'profiles_prefix')

  return {'ref': refFile, 'test': testFile, 'x': xt, 'y': yt, \
    'atm': aType, 'stat': statPrefix, 'prof': profPrefix, \
    'forcing': forcing}
# end parseConfig()

def varSetup(ref, test, diffuse=False, broadband=False, charts=False):
  """
  Extract the variables we need for further calculation and plotting

  Input
    ref -- string, path to reference netCDF file
    test -- string, path to test netCDF file

  Output
    refDict -- dictionary with reference values for flux (direct or
      diffuse), heating rate, pressure (level and layer), and TSI
      (total solar irradiance)
    testDict -- like refDict, but with test values
    numBands -- int, number of SW bands
    numLev -- int, number of levels per profile
    numProf -- int, number of columns per band
    plotVars -- string list of ordinate variables to be plotted
    fluxPlotLab -- string, used for flux (direct or diffuse) plot
      labels; one of either [DIRSTR, DIFSTR] global variables
    fluxFileStr -- string, either "direct" or "diffuse"; used in
      output PDF files
  Keywords
    diffuse -- boolean, calculate differences in diffuse flux rather
      than direct beam
    broadband -- boolean, extract and plot broadband parameters
    charts -- boolean, extract parameters for CHARTS plots
  """

  if broadband:
    dnParam = 'flux_dif_dn' if diffuse else 'flux_dir_dn'
  else:
    dnParam = 'band_flux_dif_dn' if diffuse else 'band_flux_dir_dn'
  # endif broadband

  # get the output for plotting
  if broadband:
    if charts:
      plotVars = ['flux_dn', 'flux_dir_dn', 'flux_up', 'flux_dif_dn']
    else:
      plotVars = [dnParam, 'flux_up', 'heating_rate']
    # endif charts
  else:
    if charts:
      plotVars = ['band_flux_dn', 'band_flux_dir_dn', \
        'band_flux_up', 'band_flux_dif_dn']
    else:
      plotVars = [dnParam, 'band_flux_up', 'band_heating_rate']
    # endif charts
  # endif broadband

  plotVars += ['']
  plotVars += ['p_lay', 'p_lev', 'band_lims_wvn', \
    'total_solar_irradiance']

  refDict = compare.getVars(ref, attrList=plotVars)
  testDict = compare.getVars(test, attrList=plotVars)
  #refDict['p_lev'] *= -1

  # just for now, until we correct LBLRTM pressure levels
  #refDict['p_lev'] = np.array(testDict['p_lev'])
  iNegRow, iNegCol = np.where(refDict['p_lev'] < 0)
  if iNegRow.size != 0:
    refDict['p_lev'][iNegRow, iNegCol] = \
      testDict['p_lev'][iNegRow, iNegCol]

    print 'Replacing LBLRTM pressures with RRTMGP pressures ' + \
      'for (lev, profile): '
    for row, col in zip(iNegRow, iNegCol): print row, col
  # endif iNeg

  # for SW, should only be plotting fluxes and heating rate
  # WARNING: for the rest of the module, it is assumed that this
  # order/convention (heating rate, down flux, up flux) is followed
  nVar = 4 if charts else 3
  plotVars = plotVars[:nVar]

  # only because i had insane RRTMGP band_flux_dif_dn values (inf)
  if charts: 
    testDict[plotVars[3]] = testDict[plotVars[0]] - testDict[plotVars[1]]
  # endif charts

  # some quality control (consistency check)
  # also some perturbations to test file for script verification
  for dum in plotVars:
    if refDict[dum].shape != testDict[dum].shape:
      errMsg = 'Returning: '
      errMsg += '%s and %s do not have equal dimensions for %s' % \
        (ref, test, dum)
      sys.exit(errMsg)
    # endif

    # for SW, had to introduce some pertubations to the test files so
    # we could see test-ref differences
    # doing a Gaussian noise distribution (mu=0, sigma=1)
    # testDict[dum] += np.random.normal(size=testDict[dum].shape)
  # end QC loop

  # grab dimensions of flux variables (HR is on layers, not levels)
  varShape = refDict[plotVars[0]].shape
  if broadband:
    numLev = varShape[0]; numProf = varShape[1]
    numBands = 1
  else:
    numLev = varShape[0]; numProf = varShape[1]
    numBands = varShape[2]
  # endif broadband

  fluxPlotLab = str(DIFSTR) if diffuse else str(DIRSTR)
  fluxFileStr = 'diffuse' if diffuse else 'direct'
  if charts: fluxFileStr = 'both'

  return refDict, testDict, numBands, numLev, numProf, plotVars, \
    fluxPlotLab, fluxFileStr
# end varSetup()

def profPDFs(ref, test, deltaStr, outDir='.', \
  tPauseP=100.0, broadOnly=False, diffuse=False, \
  prefix='PROFS_sens_key_cnt_add_min_04_dbl_r472_trynewn2opres', \
  atmType='Garand Atmospheres', inBand=None, logy=False):

  """
  Plot flux heating rate as well as the test-ref difference for each
  of the two variables on a single page (2x2 array of figures) for
  each column and band specified in ref and test. Fluxes and HR are
  plotted as functions of pressure, so we are plotting profiles.

  Call
    profPDFs(ref, test, outDir='.', \
      prefix='PROFS_sens_key_cnt_add_min_04_dbl_r472_trynewn2opres')

  Input
    ref -- string, path to netCDF for reference model
    test -- string, path to netCDF for test model
    deltaStr -- string, something like 'test-reference'

  Output
    prefix_band.pdf in outDir for all bands provided
      in ref and test

  Keywords
    outDir -- string, top level directory for output figures
    tPauseP -- float, tropopause pressure threshold (mbar)
    broadOnly -- boolean, only generate a broadband plot
    diffuse -- boolean, calculate differences in diffuse flux rather
      than direct beam
    prefix -- string that will be placed in front of the band
      number in the name of the output PDF file
    atmType -- float, atmosphere type (e.g., Garand)
    inBand -- int, band number (default is all)
    logy -- boolean, plots with a log-y axis rather than linear
  """

  refDict, testDict, nBands, nLev, nCol, pVars, \
    fluxStr, fluxOutFile = varSetup(ref, test, diffuse=diffuse, \
    broadband=broadOnly)

  pTitles = [''] * 6
  xTitles = ['%s %s' % (fluxStr, FUNITS), \
    '%s %s Difference %s' % (deltaStr, fluxStr, FUNITS), \
    '%s %s' % (FUPSTR, FUNITS), \
    '%s %s Difference %s' % (deltaStr, FUPSTR, FUNITS), \
    '%s %s' % (HRSTR, HRUNITS), \
    '%s %s Difference %s' % (deltaStr, HRSTR, HRUNITS)]

  # if inBand is not specified, cover all bands
  bandList = range(nBands) if inBand is None else [inBand]
  if broadOnly: bandList = [nBands]
  for iBand in bandList:
    # flux-to-heating rate conversion DIRECT DOWN!
    tempRefHR = HEATFAC * \
      np.diff(refDict[pVars[0]], axis=0)[:, :, iBand] / \
      np.diff(refDict['p_lev'], axis=0)
    tempTestHR = HEATFAC * \
      np.diff(refDict[pVars[0]], axis=0)[:, :, iBand] / \
      np.diff(testDict['p_lev'], axis=0)

    # make 1 PDF per band
    if broadOnly:
      bandStr = 'Broadband'
      outFile = '%s/%s_%s_broadband.pdf' % \
        (outDir, prefix, fluxOutFile)
      wnRange1 = refDict['band_lims_wvn'][0,:][0]
      wnRange2 = refDict['band_lims_wvn'][-1,:][1]
      wnRange = [wnRange1, wnRange2]
    else:
      bandStr = 'Band %d' % (iBand+1)
      outFile = '%s/%s_%s_%02d.pdf' % \
        (outDir, prefix, fluxOutFile, iBand+1)
      wnRange = refDict['band_lims_wvn'][iBand, :]
    # endif broadOnly

    print bandStr

    pdf = PdfPages(outFile)

    t1 = time.clock()

    for iCol in range(nCol+1):
      isLast = (iCol == nCol)
      if isLast: axMean = 1

      # plot 1 atm column per page, 3x2 array of subfigures
      print 'Column %d of %d' % (iCol+1, nCol)

      tsi = np.nan if isLast else \
        refDict['total_solar_irradiance'][iCol]

      colStr = 'Mean' if isLast else '%d' % (iCol+1)
      colStr = 'Column %s' % colStr
      figTitle = '%s %s, %s (%.2f-%.2f cm$^{-1}$)' % \
        (atmType, colStr, bandStr, wnRange[0], wnRange[1])

      fig = plot.figure()
      fig.set_size_inches(8.5, 11)
      fig.suptitle(figTitle, fontweight='bold')

			# try to center the timestamp on the bottom
      fig.text(0.45, 0.01, \
        DT.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), \
	      multialignment='center')

      for ctr in range(1, 7):
        if ctr in [1, 2]:
          # down fluxes in panels 1-2
          keyAbs = pVars[0]
        elif ctr in [3, 4]:
          # up fluxes in panels 3-4
          keyAbs = pVars[1]
        else:
          # heating rates in panels 5-6
          keyAbs = pVars[2]
        # endif keyAbs

        # heating rate (ctr = [5, 6]) is on layers
        keyOrd = 'p_lev' if ctr < 5 else 'p_lay'

        # grab plotting arrays for given band and column
        if broadOnly:
          if isLast:
            tAbscissa = np.mean(testDict[keyAbs], axis=axMean)
            rAbscissa = np.mean(refDict[keyAbs], axis=axMean)
          else:
            tAbscissa = testDict[keyAbs][:, iCol]
            rAbscissa = refDict[keyAbs][:, iCol]
          # endif isLast
        else:
          if isLast:
            tAbscissa = np.mean(testDict[keyAbs][:, :, iBand], \
              axis=axMean)
            rAbscissa = np.mean(refDict[keyAbs][:, :, iBand], \
              axis=axMean)
          else:
            tAbscissa = testDict[keyAbs][:, iCol, iBand]
            rAbscissa = refDict[keyAbs][:, iCol, iBand]
          # endif isLast
        # endif broadOnly

        dAbscissa = tAbscissa - rAbscissa
        if isLast:
          rOrdinate = np.mean(refDict[keyOrd], axis=axMean)
          tOrdinate = np.mean(testDict[keyOrd], axis=axMean)
        else:
          rOrdinate = refDict[keyOrd][:, iCol]
          tOrdinate = testDict[keyOrd][:, iCol]
        # endif isLast

        if ctr % 2 == 0:
          # profile difference
          plot.subplot(3, 2, ctr)
          plot.plot(dAbscissa, rOrdinate, 'k')

          # plot horizontal line at tropopause
          ax = plot.gca()
          ax.axhline(tPauseP, color='k', linestyle=':')
          ax.axvline(0, color='k', linestyle=':')
        else:
          # test and reference profiles
          plot.subplot(3, 2, ctr)
          plot.plot(tAbscissa, tOrdinate, 'b')
          plot.plot(rAbscissa, rOrdinate, 'r')
          plot.legend(['Test', 'Reference'], loc='upper right', \
            numpoints=1, prop=font_prop)
          plot.ylabel('Pressure [mbar]')
        # end % 2

        yRange = [max(rOrdinate), min(rOrdinate)]
        plot.ylim(yRange)
        plot.xlabel(xTitles[ctr-1])
        plot.title(pTitles[ctr-1])
        plot.xticks(rotation=15)
        if logy: plot.semilogy()

      # end ctr loop

      # write page for column
      pdf.savefig()
      plot.close()
      #break
    # end iCol loop

    pdf.close()

    print time.clock() - t1
    if broadOnly:
      print 'Processed broadband'
    else:
      print 'Processed band %d' % (iBand+1)
  # end band loop
# end profPDFs()

def profCHARTS(ref, test, deltaStr, outDir='.', \
  tPauseP=100.0, broadOnly=False, diffuse=False, \
  prefix='PROFS_sens_key_cnt_add_min_04_dbl_r472_trynewn2opres', \
  atmType='Garand Atmospheres', inBand=None, logy=False):

  """
  Plot flux heating rate as well as the test-ref difference for each
  of the two variables on a single page (2x2 array of figures) for
  each column and band specified in ref and test. Fluxes and HR are
  plotted as functions of pressure, so we are plotting profiles.

  Call
    profPDFs(ref, test, outDir='.', \
      prefix='PROFS_sens_key_cnt_add_min_04_dbl_r472_trynewn2opres')

  Input
    ref -- string, path to netCDF for reference model
    test -- string, path to netCDF for test model
    deltaStr -- string, something like 'test-reference'

  Output
    prefix_band.pdf in outDir for all bands provided
      in ref and test

  Keywords
    outDir -- string, top level directory for output figures
    tPauseP -- float, tropopause pressure threshold (mbar)
    broadOnly -- boolean, only generate a broadband plot
    diffuse -- boolean, calculate differences in diffuse flux rather
      than direct beam
    prefix -- string that will be placed in front of the band
      number in the name of the output PDF file
    atmType -- float, atmosphere type (e.g., Garand)
    inBand -- int, band number (default is all)
    logy -- boolean, plots with a log-y axis rather than linear
  """

  refDict, testDict, nBands, nLev, nCol, pVars, \
    fluxStr, fluxOutFile = varSetup(ref, test, diffuse=diffuse, \
    broadband=broadOnly, charts=True)

  downStr = r'$F^{\downarrow}$'
  pTitles = [''] * 6
  xTitles = ['%s %s' % (downStr, FUNITS), \
    '%s %s Differences %s' % (deltaStr, downStr, FUNITS), \
    '%s %s' % (FUPSTR, FUNITS), \
    '%s %s Difference %s' % (deltaStr, FUPSTR, FUNITS), \
    '%s %s' % (DIFSTR, FUNITS), \
    '%s %s Difference %s' % (deltaStr, DIFSTR, FUNITS)]

  # if inBand is not specified, cover all bands
  bandList = range(nBands) if inBand is None else [inBand]
  if broadOnly: bandList = [nBands]
  for iBand in bandList:
    # make 1 PDF per band
    if broadOnly:
      bandStr = 'Broadband'
      outFile = '%s/%s_%s_broadband.pdf' % \
        (outDir, prefix, fluxOutFile)
      wnRange1 = refDict['band_lims_wvn'][0,:][0]
      wnRange2 = refDict['band_lims_wvn'][-1,:][1]
      wnRange = [wnRange1, wnRange2]
    else:
      bandStr = 'Band %d' % (iBand+1)
      outFile = '%s/%s_%s_%02d.pdf' % \
        (outDir, prefix, fluxOutFile, iBand+1)
      wnRange = refDict['band_lims_wvn'][iBand, :]
    # endif broadOnly

    print bandStr

    pdf = PdfPages(outFile)

    t1 = time.clock()

    for iCol in range(nCol+1):
      isLast = (iCol == nCol)
      if isLast: axMean = 1

      # plot 1 atm column per page, 3x2 array of subfigures
      print 'Column %d of %d' % (iCol+1, nCol)

      tsi = np.nan if isLast else \
        refDict['total_solar_irradiance'][iCol]

      colStr = 'Mean' if isLast else '%d' % (iCol+1)
      colStr = 'Column %s' % colStr
      figTitle = '%s %s, %s (%.2f-%.2f cm$^{-1}$)' % \
        (atmType, colStr, bandStr, wnRange[0], wnRange[1])

      fig = plot.figure()
      fig.set_size_inches(8.5, 11)
      fig.suptitle(figTitle, fontweight='bold')

			# try to center the timestamp on the bottom
      fig.text(0.45, 0.01, \
        DT.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), \
	      multialignment='center')

      for ctr in range(1, 7):
        if ctr in [1, 2]:
          # down fluxes in panels 1-2
          keyAbs = pVars[0]
          keyAbsDif = pVars[1]
        elif ctr in [3, 4]:
          # up fluxes in panels 3-4
          keyAbs = pVars[2]
        else:
          # heating rates in panels 5-6
          keyAbs = pVars[3]
        # endif keyAbs

        # heating rate (ctr = [5, 6]) is on layers
        keyOrd = 'p_lev' #if ctr < 5 else 'p_lay'

        # grab plotting arrays for given band and column
        if broadOnly:
          if isLast:
            tAbscissa = np.mean(testDict[keyAbs], axis=axMean)
            rAbscissa = np.mean(refDict[keyAbs], axis=axMean)

            if ctr in [1, 2]:
              tAbscissaDif = np.mean(testDict[keyAbsDif], axis=axMean)
              rAbscissaDif = np.mean(refDict[keyAbsDif], axis=axMean)
            # endif ctr
          else:
            tAbscissa = testDict[keyAbs][:, iCol]
            rAbscissa = refDict[keyAbs][:, iCol]

            if ctr in [1, 2]:
              tAbscissaDif = testDict[keyAbsDif][:, iCol]
              rAbscissaDif = refDict[keyAbsDif][:, iCol]
            # endif ctr
          # endif isLast
        else:
          if isLast:
            tAbscissa = np.mean(testDict[keyAbs][:, :, iBand], \
              axis=axMean)
            rAbscissa = np.mean(refDict[keyAbs][:, :, iBand], \
              axis=axMean)

            if ctr in [1, 2]:
              tAbscissaDif = np.mean(\
                testDict[keyAbsDif][:, :, iBand], axis=axMean)
              rAbscissaDif = np.mean(\
                refDict[keyAbsDif][:, :, iBand], axis=axMean)
            # endif ctr
          else:
            tAbscissa = testDict[keyAbs][:, iCol, iBand]
            rAbscissa = refDict[keyAbs][:, iCol, iBand]

            if ctr in [1, 2]:
              tAbscissaDif = testDict[keyAbsDif][:, iCol, iBand]
              rAbscissaDif = refDict[keyAbsDif][:, iCol, iBand]
            # endif ctr
          # endif isLast
        # endif broadOnly

        dAbscissa = tAbscissa - rAbscissa
        if ctr in [1, 2]: dAbscissaDif = tAbscissaDif - rAbscissaDif

        if isLast:
          rOrdinate = np.mean(refDict[keyOrd], axis=axMean)
          tOrdinate = np.mean(testDict[keyOrd], axis=axMean)
        else:
          rOrdinate = refDict[keyOrd][:, iCol]
          tOrdinate = testDict[keyOrd][:, iCol]
        # endif isLast

        if ctr % 2 == 0:
          # profile difference
          plot.subplot(3, 2, ctr)
          plot.plot(dAbscissa, rOrdinate, 'k')

          # plot horizontal line at tropopause
          ax = plot.gca()
          ax.axhline(tPauseP, color='k', linestyle=':')
          ax.axvline(0, color='k', linestyle=':')

          # plot dir down and total down together
          if ctr in [1, 2]: plot.plot(dAbscissaDif, tOrdinate, 'k--')
        else:
          # test and reference profiles
          plot.subplot(3, 2, ctr)
          plot.plot(tAbscissa, tOrdinate, 'b')
          plot.plot(rAbscissa, rOrdinate, 'r')
          plot.legend(['Test', 'Reference'], loc='upper right', \
            numpoints=1, prop=font_prop)
          plot.ylabel('Pressure [mbar]')

          # plot dir down and total down together
          if ctr in [1, 2]:
            plot.plot(tAbscissaDif, tOrdinate, 'b--')
            plot.plot(rAbscissaDif, rOrdinate, 'r--')
          # endif ctr

        # end % 2

        yRange = [max(rOrdinate), min(rOrdinate)]
        plot.ylim(yRange)
        plot.xlabel(xTitles[ctr-1])
        plot.title(pTitles[ctr-1])
        plot.xticks(rotation=15)
        if logy: plot.semilogy()

      # end ctr loop

      # write page for column
      pdf.savefig()
      plot.close()
      #break
    # end iCol loop

    pdf.close()

    print time.clock() - t1
    if broadOnly:
      print 'Processed broadband'
    else:
      print 'Processed band %d' % (iBand+1)
  # end band loop
# end profCHARTS()

def calcDiffs(refArr, testArr, tPauseIdx, heatrate=False):
  """
  Given a reference and test array and an index corresponding to the
  tropopause, calculate the differences in absorption over the entire
  troposphere and stratosphere

  Input
    refArr -- float array, reference model direct downwelling
      flux array
    testArr -- float array, test model direct downwelling flux array
    tPauseIdx -- int, index corresponding to tropopause

  Output
    outDict -- dictionary with the following keys:

      ref_surface -- float, reference surface flux value
      ref_trop -- float, reference total tropospheric flux value
      ref_strat -- float, reference total stratospheric flux value
      diff_surface -- float, test-ref difference in surface
        flux
      diff_trop -- float, test-ref difference in total tropospheric
        flux
      diff_strat -- float, test-ref difference in total stratospheric
        flux

  Keywords
    heatrate -- boolean, calculate maximum heating rate differences
     instead of total flux for a given layer
  """

  if heatrate:
    refSfc = np.nan; diffSfc = np.nan
    refTrop = max(refArr[:tPauseIdx])
    testTrop = max(testArr[:tPauseIdx])
    refStrat = max(refArr[tPauseIdx:])
    testStrat = max(testArr[tPauseIdx:])
  else:
    # surface is the first index of the array, so last index is TOA
    refSfc = refArr[0]
    diffSfc = testArr[0] - refArr[0]
    refTrop = refArr[tPauseIdx] - refArr[0]
    testTrop = testArr[tPauseIdx] - testArr[0]
    refStrat = refArr[-1] - refArr[tPauseIdx]
    testStrat = testArr[-1] - testArr[tPauseIdx]
  # endif heatrate

  diffTrop = testTrop - refTrop
  diffStrat = testStrat - refStrat

  outDict = {'ref_surface': refSfc, 'diff_surface': diffSfc, \
    'ref_trop': refTrop, 'diff_trop': diffTrop, \
    'ref_strat': refStrat, 'diff_strat': diffStrat}
  """
  outDict = {'ref_surface': refSfc, 'diff_surface': diffSfc, \
    'ref_trop': refTrop, 'diff_trop': diffTrop, \
    'ref_strat': refStrat, 'diff_strat': diffStrat, \
    'dTropMean': np.mean(diffTrop), \
    'dTropSD': np.std(diffTrop, ddof=1), \
    'dStratMean': np.mean(diffStrat), \
    'dStratSD': np.std(diffStrat, ddof=1), \
    'dTropAbsMean': np.mean(np.abs(diffTrop)), \
    'dTropAbsSD': np.std(np.abs(diffTrop), ddof=1), \
    'dStratAbsMean': np.mean(np.abs(diffStrat)), \
    'dStratAbsSD': np.std(np.abs(diffStrat), ddof=1)}
  """

  return outDict
# end calcFluxDiffs()

def statPDF(ref, test, tPauseP=100.0, atmType='Garand', \
  xTitle='Reference', yTitle='Test-Reference', forcing=False,
  prefix='column_SW_flux_HR_stats_', diffuse=False, outDir='.'):
  """
  Generate column RRTMGP-LBLRTM residual arrays as a function of
  associated LBLRTM measurements (for SW flux and heating rate) at the
  surface and tropopause, then plot the arrays

  Call
    statPDF(ref, test)

  Input
    ref -- string, path to netCDF for reference model
    test -- string, path to netCDF for test model

  Output
    prefix_all_bands.pdf in outDir for all columns provided
      in ref and test

  Keywords
    prefix -- string that will be placed in front of the band
      number in the name of the output PDF file
    tPauseP -- float, tropopause pressure threshold (mbar)
    atmType -- float, atmosphere type (e.g., Garand)
    diffuse -- boolean, calculate differences in diffuse flux rather
      than direct beam
    xTitle -- string, label for x-axis (unitless)
    yTitle -- string, label for y-axis (unitless)
    forcing -- boolean, if set, plot the 5% and 10% error lines on
      all plots
    outDir -- string, top level directory for output figures
  """
  refDict, testDict, nBands, nLev, nCol, pVars, \
    fluxStr, fluxOutFile = varSetup(ref, test, diffuse=diffuse)
  # broadband page will be the last page of output PDF, and we only
  # need the broadband dictionaries and pVar list
  refDictBB, testDictBB, dum, dum, dum, pVarsBB, dum, dum = \
    varSetup(ref, test, diffuse=diffuse, broadband=True)
  outFile = '%s/%s_%s_all_bands.pdf' % (outDir, prefix, fluxOutFile)
  pdfObj = PdfPages(outFile)

  xTitles = ['', '', '', '', \
    #'%s %s %s' % (xTitle, fluxStr, FUNITS), \
    #'%s %s %s' % (xTitle, fluxStr, FUNITS), \
    #'%s %s %s' % (xTitle, fluxStr, FUNITS), \
    #'%s %s %s' % (xTitle, HRSTR, HRUNITS), \
    '%s %s %s' % (xTitle, fluxStr, FUNITS), \
    '%s %s %s' % (xTitle, HRSTR, HRUNITS)]
  yTitles = ['%s %s' % (yTitle, FUNITS), '', \
    '%s %s' % (yTitle, FUNITS), '', \
    '%s %s' % (yTitle, HRUNITS), '']

  # for error lines in forcing plots
  mSize = 0.5; oPlotSym = 'k:'

  # for statistics and easier string formatting
  muStr = r'$\mu$'; muAbsStr = r'$|\mu|$'
  sdStr = r'$\sigma$'; sdAbsStr = r'$|\sigma|$'

  # universal font size for this page
  font = {'size': 8}
  rc('font', **font)

  for iBand in range(nBands+1):
    isBB = (iBand == nBands)
    if isBB:
      print 'Processing broadband'

      wnRange1 = refDict['band_lims_wvn'][:,0][0]
      wnRange2 = refDict['band_lims_wvn'][:,-1][1]
      wnRange = [wnRange1, wnRange2]
      bandStr = 'Broadband'
      pVars = list(pVarsBB)
    else:
      bandStr = 'Band %d' % (iBand+1)
      print 'Processing %s' % bandStr

      wnRange = refDict['band_lims_wvn'][iBand, :]
    # endif iBB

    # test-ref differences in flux and HR at surface and tropopause
    dFluxSfc, dFluxTrop, dFluxStrat, dHRTrop, dHRStrat = \
      [], [], [], [], []

    # reference values for 3 parameters
    rFluxSfc, rFluxTrop, rFluxStrat, rHRTrop, rHRStrat = \
      [], [], [], [], []

    for iCol in range(nCol):
      # grab total solar irradiance for figure caption/title
      tsi = refDict['total_solar_irradiance'][iCol]

      # find tropopause index
      pTemp = refDict['p_lay'][:, iCol]
      iPause = np.argmin(np.abs(pTemp-tPauseP))

      # grab the values for a given band and column to be used in
      # calcDiffs() and plotted
      if isBB:
        tempRefFlux = refDictBB[pVars[0]][:, iCol]
        tempTestFlux = testDictBB[pVars[0]][:, iCol]
      else:
        tempRefFlux = refDict[pVars[0]][:, iCol, iBand]
        tempTestFlux = testDict[pVars[0]][:, iCol, iBand]
      # endif isBB

      # flux-to-heating rate conversion
      tempRefHR = HEATFAC * \
        np.diff(tempRefFlux) / np.diff(refDict['p_lev'][:, iCol])
      tempTestHR = HEATFAC * \
        np.diff(tempTestFlux) / np.diff(testDict['p_lev'][:, iCol])

      # for y values (test-ref differences)
      dFlux = calcDiffs(tempRefFlux, tempTestFlux, iPause)
      dHR = calcDiffs(tempRefHR, tempTestHR, iPause, heatrate=True)

      # append single x (reference) values for a given column to lists
      rFluxSfc.append(dFlux['ref_surface'])
      rFluxTrop.append(dFlux['ref_trop'])
      rFluxStrat.append(dFlux['ref_strat'])
      rHRTrop.append(dHR['ref_trop'])
      rHRStrat.append(dHR['ref_strat'])

      # append single y values (test-ref differences)
      # for a given column to lists
      dFluxSfc.append(dFlux['diff_surface'])
      dFluxTrop.append(dFlux['diff_trop'])
      dFluxStrat.append(dFlux['diff_strat'])
      dHRTrop.append(dHR['diff_trop'])
      dHRStrat.append(dHR['diff_strat'])
    # end column loop

    ordinates = [dFluxSfc, dFluxTrop, dHRTrop, \
      dFluxStrat, dHRStrat]

    # statistics for panel titles
    meanOrd = [np.mean(o) for o in ordinates]
    meanAbsOrd = [np.mean(np.abs(o)) for o in ordinates]
    sdOrd = [np.std(o, ddof=1) for o in ordinates]
    sdAbsOrd = [np.std(np.abs(o), ddof=1) for o in ordinates]

    # skipping over second plot, so have to put an empty list at
    # second list element
    ordinates.insert(1, [])
    meanOrd.insert(1, [])
    meanAbsOrd.insert(1, [])
    sdOrd.insert(1, [])
    sdAbsOrd.insert(1, [])

    abscissae = [rFluxSfc, [], rFluxTrop, rHRTrop, \
      rFluxStrat, rHRStrat]

    # figure (i.e. page) settings
    figTitle = '%s %s (%.2f-%.2f cm$^{-1}$)\n' % \
      (atmType, bandStr, wnRange[0], wnRange[1])
    figTitle += 'Reference TSI = %.2f' % tsi
    fig = plot.figure()
    fig.set_size_inches(8.5, 11)
    fig.suptitle(figTitle, fontweight='bold')

    # start making plots; 3x2 panels of flux at surface and
    # tropopause (downwelling, then upwelling), then heating rate at
    # surface and tropopause
    t1 = time.clock()
    for ctr in range(1, 7):
      if ctr == 2: continue
      abscissa = np.array(abscissae[ctr-1])
      ordinate = np.array(ordinates[ctr-1])

      plot.subplot(3, 2, ctr)
      for x, y, i in zip(abscissa, ordinate, range(nCol)):
        plot.plot(x, y, 'k', marker='$%d$' % (i+1), markersize=10)

      plot.xlabel(xTitles[ctr-1])
      plot.ylabel(yTitles[ctr-1])
      plot.xticks(rotation=15)
      plot.title('%s = %.3f, %s = %.3f\n%s = %.3f, %s = %.3f' % \
        (muStr, meanOrd[ctr-1], sdStr, sdOrd[ctr-1], \
         muAbsStr, meanAbsOrd[ctr-1], sdAbsStr, sdAbsOrd[ctr-1]) )

      # overplot 5% and 10% error lines (+/-)
      if forcing:
        xZeros = np.array([min(abscissa), max(abscissa)])

        # make a y = ax line for each threshold
        err5 = np.poly1d([0.05, 0]); err10 = np.poly1d([0.10, 0])
        plot.plot(xZeros, err5(xZeros), oPlotSym, markersize=mSize)
        plot.plot(xZeros, err10(xZeros), oPlotSym, markersize=mSize)

        err5 = np.poly1d([-0.05, 0]); err10 = np.poly1d([-0.10, 0])
        plot.plot(xZeros, err5(xZeros), oPlotSym, markersize=mSize)
        plot.plot(xZeros, err10(xZeros), oPlotSym, markersize=mSize)
      # end error lines
    # end ctr loop

		# try to center the timestamp on the bottom of the page
		# all of the positioning is trial and error given my setup
    fig.text(0.45, 0.01, \
      DT.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), \
      multialignment='center')

		# now try to center common x and y labels (layer labels)
    xPos = 0.02
    fig.text(xPos, 0.25, 'Stratosphere', multialignment='center', \
      rotation='vertical')
    fig.text(xPos, 0.54, 'Troposphere', multialignment='center', \
      rotation='vertical')
    fig.text(xPos, 0.8, 'Surface', multialignment='center', \
      rotation='vertical')

    pdfObj.savefig()
  # end band loop

  pdfObj.close()
  return True
# end statPDF()

if __name__ == '__main__':
  parser = argparse.ArgumentParser(\
    description='Generate comparative plots for calculations ' + \
    'from two radiative transfer models.')
  parser.add_argument('--config_file', type=str, \
    default='LBLRTM_RRTMGP_SW_config.ini', \
    help='Path to configuration file that contains the values ' + \
    'to keyword arguments.')
  parser.add_argument('--plot_profiles', action='store_true', \
    help='Generates an N-band page PDF, each with a 3x2 array of ' + \
    'plots that compare the upward and downwward fluxes from ' + \
    'the reference and test models. Also plots heating rate for ' + \
    'each model.')
  parser.add_argument('--plot_stats', action='store_true', \
    help='Generates a single page with RRTMGP-LBLRTM vs. LBLRTM ' + \
    'statistics for TOA upwelling, SFC downwelling, and Max Net ' + \
    'flux cases as well as Max Tropospheric and Stratospheric HR.')
  parser.add_argument('--tropopause_pressure', type=float, \
    default=np.exp(4.6), \
    help='Pressure threshold that separates troposphere and ' + \
    'stratosphere.')
  parser.add_argument('--band', type=int, nargs='+', \
    help='Number of band to plot. Default (band=None) is all bands.')
  parser.add_argument('--log_y', action='store_true', \
    help='Generate a semilog-y plot.')
  parser.add_argument('--diffuse', action='store_true', \
    help='Plot diffuse downwelling flux array instead of flux ' + \
    'from direct beam.')
  parser.add_argument('--broad_only', action='store_true', \
    help='Only generate a broadband plot.')
  parser.add_argument('--charts', action='store_true', \
    help='Plot panels relevant for CHARTS study (PROFILES ONLY).')
  parser.add_argument('--out_dir', type=str, default='.', \
    help='Path to output directory in which PDF files are saved.')
  args = parser.parse_args()

  conFile = args.config_file; utils.file_check(conFile)
  pTrop = args.tropopause_pressure
  doDif = args.diffuse
  figOutDir = args.out_dir
  utils.file_check(figOutDir)

  configDict = parseConfig(conFile)

  # plot the test and reference upward and downward flux together
  # with the heating rates (and differences for each)
  inBand = None if args.band is None else np.array(args.band)-1
  if args.plot_profiles:
    if inBand is None:
      if args.charts:
        profCHARTS(configDict['ref'], configDict['test'], \
          configDict['y'], tPauseP=pTrop, \
          prefix=configDict['prof'], atmType=configDict['atm'], \
          inBand=inBand, logy=args.log_y, broadOnly=args.broad_only, \
          diffuse=doDif, outDir=figOutDir)
      else:  
        profPDFs(configDict['ref'], configDict['test'], \
          configDict['y'], tPauseP=pTrop, \
          prefix=configDict['prof'], atmType=configDict['atm'], \
          inBand=inBand, logy=args.log_y, broadOnly=args.broad_only, \
          diffuse=doDif, outDir=figOutDir)
      # endif charts
    else:
      for iBand in inBand:
        if args.charts:
          profCHARTS(configDict['ref'], configDict['test'], \
            configDict['y'], tPauseP=pTrop, inBand=iBand, \
            prefix=configDict['prof'], atmType=configDict['atm'], \
            logy=args.log_y, diffuse=doDif, outDir=figOutDir)
        else:  
          profPDFs(configDict['ref'], configDict['test'], \
            configDict['y'], tPauseP=pTrop, inBand=iBand, \
            prefix=configDict['prof'], atmType=configDict['atm'], \
            logy=args.log_y, diffuse=doDif, outDir=figOutDir)
        # endif charts
      # end iBand loop

      # for specified bands AND broadband
      if args.broad_only:
        if args.charts:
          profCHARTS(configDict['ref'], configDict['test'], \
            configDict['y'],  tPauseP=pTrop, \
            prefix=configDict['prof'], atmType=configDict['atm'], \
            logy=args.log_y, broadOnly=args.broad_only, \
            diffuse=doDif, outDir=figOutDir)
        else:
          profPDFs(configDict['ref'], configDict['test'], \
            configDict['y'],  tPauseP=pTrop, \
            prefix=configDict['prof'], atmType=configDict['atm'], \
            logy=args.log_y, broadOnly=args.broad_only, \
            diffuse=doDif, outDir=figOutDir)
        # endif charts
      # end broadband plot
    # end inBand
  # end plot_profiles

  # plot profile statistics
  if args.plot_stats:
    statPDF(configDict['ref'], configDict['test'], tPauseP=pTrop, \
      atmType=configDict['atm'], \
      xTitle=configDict['x'], yTitle=configDict['y'], \
      forcing=configDict['forcing'], prefix=configDict['stat'], \
      diffuse=args.diffuse, outDir=figOutDir)

# end main()
