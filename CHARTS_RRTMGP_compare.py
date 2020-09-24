#!/usr/bin/env python

import os, sys, argparse
import numpy as np
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
import SW_LBLRTM_RRTMGP_compare as sw_compare

# like in sw_compare, strings for plotting axes
DOWNDIRSTR = r'$F^{\downarrow}_{Direct}$'
NETSTR = r'$F_{net}^{max}$'
DOWNTOTSTR = r'$F^{\downarrow}$'

class refTestFluxCompare():
  def __init__(self, inDict, outDir=os.getcwd(), tPauseP=100.0, \
    bands=None, yLog=False, broadband=False):
    """
    Extract the variables we need for further calculation and plotting

    It starts with just fluxes, then heating rates are derived from 
    those fluxes.

    Input
      inDict -- dictionary generated with sw_compare.parseConfig()

    Keywords
      outDir -- string, top level directory for output figures
      tPauseP -- float, tropopause pressure threshold (mbar)
      bands -- int array (or list) of bands to plot
      yLog -- boolean, plot the test-model differences on a log scale
      broadband -- boolean, plot differences integrated over all 
        bands
    """

    ref = inDict['ref']; test = inDict['test']
    self.plotVars = \
      ['flux_dn', 'flux_dir_dn', 'flux_dif_dn', 'flux_up', \
      'flux_net', 'heating_rate', \
      'band_flux_dn', 'band_flux_dir_dn', 'band_flux_dif_dn', \
      'band_flux_up', 'band_flux_net', 'band_heating_rate', \
      'p_lay', 'p_lev', 'band_lims_wvn', 'total_solar_irradiance']

    # generate dictionaries with arrays corresponding to all plotVars
    self.refDict = compare.getVars(ref, attrList=self.plotVars)
    self.testDict = compare.getVars(test, attrList=self.plotVars)

    # only because i had insane RRTMGP band_flux_dif_dn values (inf)
    self.testDict['flux_dif_dn'] = \
      self.testDict['flux_dn'] - self.testDict['flux_dir_dn']
    self.testDict['band_flux_dif_dn'] = \
      self.testDict['band_flux_dn'] - self.testDict['band_flux_dir_dn']

    # some quality control (consistency check)
    # also some perturbations to test file for script verification
    for dum in self.plotVars:
      if self.refDict[dum].shape != self.testDict[dum].shape:
        errMsg = 'Returning: '
        errMsg += '%s and %s do not have equal dimensions for %s' % \
          (ref, test, dum)
        sys.exit(errMsg)
      # endif
    # end QC loop

    # grab dimensions of flux variables (using a by-band variable 
    # instead of broadband so we can deduce number of bands)
    varShape = self.refDict['band_flux_dn'].shape
    self.numLev = varShape[0]
    self.numLay = varShape[0] - 1
    self.numProf = varShape[1]
    self.numBands = varShape[2]

    self.tPauseP = float(tPauseP)

    # other parameters extracted from the configuration file
    self.statPrefix = inDict['stat']
    self.profPrefix = inDict['prof']
    self.xTitle = inDict['x']
    self.yTitle = inDict['y']
    self.atmType = inDict['atm']
    self.forcing = inDict['forcing']
    self.outDir = str(outDir)

    self.broadBand = bool(broadband)
    self.yLog = bool(yLog)
    self.bands = None if bands is None else np.array(bands)
  # end constructor

  def calcDiffs(self, inRefDict, inTestDict, tPauseIdx, diffVars):
    """
    Compute differences in TOA up, direct down, total down, and 
    net flux. This is for stats only

    Input
      inRefDict -- dictionary with the following keys:
        band_flux_up, band_flux_dn, band_flux_dir_dn, band_flux_net 
          band_heating_rate if by-band
        flux_up, flux_dn, flux_dir_dn, flux_net, heating_rate 
          if broadband

      inTestDict -- dictionary, same as inRefDict but with test arrays
      tPauseIdx -- int, index corresponding to tropopause
      diffVars -- string list, names of variables to diff

    Output
      outDict -- dictionary with the following keys (with associated 
        float array values; references and test-ref differences):

        ref_flux_TOA_up, diff_flux_TOA_up
        ref_flux_Sfc_dir_dn, diff_flux_Sfc_dir_dn
        ref_flux_Sfc_dif_dn, diff_flux_Sfc_dif_dn
        ref_flux_Total_dn, diff_flux_Total_dn
        ref_flux_Max_net, diff_flux_Max_net
        ref_HR_Max_strat, diff_HR_Max_strat
        ref_HR_Max_tropo, diff_HR_Max_tropo

    Keywords
      broadband -- boolean, do the broadband (instead of by-band) 
        calculations (so broadband in*Dict keys are expected)
    """

    # surface is the first index of the array, so last index is TOA
    dVar = diffVars[0]
    refDnSfc = inRefDict[dVar][0]
    diffDnSfc = (inTestDict[dVar]-inRefDict[dVar])[0]

    dVar = diffVars[1]
    refDirDnSfc = inRefDict[dVar][0]
    diffDirDnSfc = (inTestDict[dVar]-inRefDict[dVar])[0]

    dVar = diffVars[2]
    refDifDnSfc = inRefDict[dVar][0]
    diffDifDnSfc = (inTestDict[dVar]-inRefDict[dVar])[0]

    dVar = diffVars[3]
    refUpTOA = inRefDict[dVar][-1]
    diffUpTOA = (inTestDict[dVar]-inRefDict[dVar])[-1]

    # find index of max diff
    dVar = diffVars[4]
    diffNet = inTestDict[dVar]-inRefDict[dVar]
    iMax = np.argmax(np.abs(diffNet))
    refNetMax = inRefDict[dVar][iMax]
    diffNetMax = diffNet[iMax]

    # find index of max HR diff in troposphere and stratosphere
    dVar = diffVars[5]
    diffTropoHR = \
      (inTestDict[dVar]-inRefDict[dVar])[:tPauseIdx]
    iTropo = np.argmax(np.abs(diffTropoHR))
    refHrTropoMax = inRefDict[dVar][iTropo]
    diffHrTropoMax = diffTropoHR[iTropo]

    diffStratHR = \
      (inTestDict[dVar]-inRefDict[dVar])[tPauseIdx:]
    iStrat = np.argmax(np.abs(diffStratHR))
    refHrStratMax = inRefDict[dVar][iStrat]
    diffHrStratMax = diffStratHR[iStrat]

    outDict = {'ref_flux_toa_up': refUpTOA, \
      'diff_flux_toa_up': diffUpTOA, \
      'ref_flux_sfc_dir_dn': refDirDnSfc, \
      'diff_flux_sfc_dir_dn': diffDirDnSfc, \
      'ref_flux_tot_dn': refDnSfc, \
      'diff_flux_tot_dn': diffDnSfc, \
      'ref_flux_sfc_dif_dn': refDifDnSfc, \
      'diff_flux_sfc_dif_dn': diffDifDnSfc, \
      'ref_flux_max_net': refNetMax, \
      'diff_flux_max_net': diffNetMax, \
      'ref_HR_max_strat': refHrStratMax, \
      'diff_HR_max_strat': diffHrStratMax, \
      'ref_HR_max_tropo': refHrTropoMax, \
      'diff_HR_max_tropo': diffHrTropoMax}

    return outDict
  # end calcDiffs()

  def statPDF(self):
    """
    Generate column RRTMGP-CHARTS residual arrays as a function of
    associated CHARTS output (for SW flux and heating rate) at 
    the surface and TOA, then plot the arrays
    """

    outFile = '%s/%s_all_bands.pdf' % (self.outDir, self.statPrefix)
    pdfObj = PdfPages(outFile)

    xTitle = self.xTitle; yTitle = self.yTitle
    xTitles = ['%s TOA %s %s' % \
        (self.xTitle, sw_compare.FUPSTR, sw_compare.FUNITS), \
      '%s Surface %s %s' % \
        (self.xTitle, sw_compare.DIRSTR, sw_compare.FUNITS), \
      '%s Surface %s %s' % \
        (self.xTitle, DOWNTOTSTR, sw_compare.FUNITS), \
      '%s %s %s' % (self.xTitle, NETSTR, sw_compare.FUNITS), \
      '%s Trop %s %s' % \
        (self.xTitle, sw_compare.HRSTR, sw_compare.HRUNITS), \
      '%s Strat %s %s' % \
        (self.xTitle, sw_compare.HRSTR, sw_compare.HRUNITS)]
    yTitles = ['%s TOA %s %s' % \
        (self.yTitle, sw_compare.FUPSTR, sw_compare.FUNITS), \
      '%s Surface %s %s' % \
        (self.yTitle, sw_compare.DIRSTR, sw_compare.FUNITS), \
      '%s Surface %s %s' % \
        (self.yTitle, DOWNTOTSTR, sw_compare.FUNITS), \
      '%s %s %s' % (self.yTitle, NETSTR, sw_compare.FUNITS), \
      '%s Trop %s %s' % \
        (self.yTitle, sw_compare.HRSTR, sw_compare.HRUNITS), \
      '%s Strat %s %s' % \
        (self.yTitle, sw_compare.HRSTR, sw_compare.HRUNITS)]

    # for error lines in forcing plots
    mSize = 0.5; oPlotSym = 'k:'

    # for statistics and easier string formatting
    muStr = r'$\mu$'; muDifStr = r'$\mu_{Dif}$'
    sdStr = r'$\sigma$'; sdDifStr = r'$\sigma_{Dif}$'

    # universal font size for this page
    font = {'size': 8}
    rc('font', **font)

    for iBand in range(self.numBands+1):
      isBB = (iBand == self.numBands)
      if isBB:
        print 'Processing broadband'

        wnRange1 = self.refDict['band_lims_wvn'][:,0][0]
        wnRange2 = self.refDict['band_lims_wvn'][:,-1][1]
        wnRange = [wnRange1, wnRange2]
        bandStr = 'Broadband'
        pVars = list(self.plotVars[:6])
      else:
        bandStr = 'Band %d' % (iBand+1)
        print 'Processing %s' % bandStr

        wnRange = self.refDict['band_lims_wvn'][iBand, :]
        pVars = list(self.plotVars[6:12])
      # endif iBB

      # test-ref differences in flux and HR at surface and tropopause
      dUpTOA, dDirSfc, dDownSfc, dDifSfc, dNetMax, dHRTropMax, \
        dHRStratMax = [], [], [], [], [], [], []

      # reference values for 3 parameters
      rUpTOA, rDirSfc, rDownSfc, rDifSfc, rNetMax, rHRTropMax, \
        rHRStratMax = [], [], [], [], [], [], []

      for iCol in range(self.numProf):
        # grab total solar irradiance for figure caption/title
        tsi = self.refDict['total_solar_irradiance'][iCol]

        # find tropopause index
        pLay = self.refDict['p_lay'][:, iCol]
        iPause = np.argmin(np.abs(pLay-self.tPauseP))

        # grab the values for a given band and column to be used in
        # calcDiffs() and plotted
        tempRef, tempTest = {}, {}
        for pVar in pVars:
          if isBB:
            tempRef[pVar] = self.refDict[pVar][:, iCol]
            tempTest[pVar] = self.testDict[pVar][:, iCol]
          else:
            tempRef[pVar] = self.refDict[pVar][:, iCol, iBand]
            tempTest[pVar] = self.testDict[pVar][:, iCol, iBand]
          # endif isBB
        # end pVars loop

        # for x (reference) values
        diffDict = self.calcDiffs(tempRef, tempTest, iPause, pVars)
        self.statsDict = dict(diffDict)
        rUpTOA.append(diffDict['ref_flux_toa_up'])
        rDirSfc.append(diffDict['ref_flux_sfc_dir_dn'])
        rDifSfc.append(diffDict['ref_flux_sfc_dif_dn'])
        rDownSfc.append(diffDict['ref_flux_tot_dn'])
        rNetMax.append(diffDict['ref_flux_max_net'])
        rHRTropMax.append(diffDict['ref_HR_max_tropo'])
        rHRStratMax.append(diffDict['ref_HR_max_strat'])

        # append single y values (test-ref differences)
        # for a given column to lists
        dUpTOA.append(diffDict['diff_flux_toa_up'])
        dDirSfc.append(diffDict['diff_flux_sfc_dir_dn'])
        dDifSfc.append(diffDict['diff_flux_sfc_dif_dn'])
        dDownSfc.append(diffDict['diff_flux_tot_dn'])
        dNetMax.append(diffDict['diff_flux_max_net'])
        dHRTropMax.append(diffDict['diff_HR_max_tropo'])
        dHRStratMax.append(diffDict['diff_HR_max_strat'])
      # end column loop

      ordinates = [dUpTOA, dDirSfc, dDifSfc, dDownSfc, dNetMax, \
        dHRTropMax, dHRStratMax]

      # statistics for panel titles
      meanOrd = [np.mean(o) for o in ordinates]
      sdOrd = [np.std(o, ddof=1) for o in ordinates]

      # remove diffuse from ordinates
      del ordinates[2]
      avgDiffuse = meanOrd.pop(2)
      sdDiffuse = sdOrd.pop(2)

      abscissae = [rUpTOA, rDirSfc, rDownSfc, rNetMax, rHRTropMax, \
        rHRStratMax]

      # figure (i.e. page) settings
      figTitle = '%s %s (%.2f-%.2f cm$^{-1}$)\n' % \
        (self.atmType, bandStr, wnRange[0], wnRange[1])
      figTitle += 'Reference TSI = %.2f' % tsi
      fig = plot.figure()
      fig.set_size_inches(8.5, 11)
      fig.suptitle(figTitle, fontweight='bold')

      # start making plots; 3x2 panels of flux at surface and
      # tropopause (downwelling, then upwelling), then heating rate at
      # surface and tropopause
      t1 = time.clock()
      for ctr in range(1, 7):
        abscissa = np.array(abscissae[ctr-1])
        ordinate = np.array(ordinates[ctr-1])

        plot.subplot(3, 2, ctr)
        for x, y, i in zip(abscissa, ordinate, range(self.numProf)):
          plot.plot(x, y, 'k', marker='$%d$' % (i+1), markersize=10)

        plot.xlabel(xTitles[ctr-1])
        plot.ylabel(yTitles[ctr-1])
        if ctr == 3:
          mom1 = '%s = %.3f, %s = %.3f' % \
            (muStr, meanOrd[ctr-1], sdStr, sdOrd[ctr-1])
          mom2 = '%s = %.3f, %s = %.3f' % \
            (muDifStr, avgDiffuse, sdDifStr, sdDiffuse)
          plot.title('%s; %s' % (mom1, mom2))
        else:
          plot.title('%s = %.3f, %s = %.3f' % \
            (muStr, meanOrd[ctr-1], sdStr, sdOrd[ctr-1]) )
        # end ctr 3

        # overplot 5% and 10% error lines (+/-)
        if self.forcing:
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

      pdfObj.savefig()
    # end band loop

    pdfObj.close()
    print '%s completed' % outFile
    print 'Stats plotting completed'

  # end statPDF()

  def profPDF(self):

    """
    Test and reference fluxes are plotted as functions of pressure, 
    then test-model differences are plotted as functions of pressure 
    as well. This is done for 3 sets of variables 
    (total down, up, down diffuse; down direct is then overplotted 
     onto total down)
    """

    # this method was copied from SW_LBLRTM_RRTMGP_compare.py and 
    # originally extracted these guys with varSetup() in that module,
    # but we already do this in the constructor, and for now it's 
    # just easier to assign
    refDict = self.refDict
    testDict = self.testDict
    deltaStr = self.yTitle

    nBands, nLev, nCol = self.numBands, self.numLev, self.numProf
    pVars = self.plotVars[:4] if self.broadBand else \
      self.plotVars[6:10]

    # and since this was copied from another library, i exposed a flaw
    # that the pVar order could be different. here i will conform to
    # the convention expected in SW_LBLRTM_RRTMGP_compare.py
    iSort = [0, 1, 3, 2]
    pVars = [pVars[i] for i in iSort]
 
    pTitles = [''] * 6
    xTitles = ['%s %s' % (DOWNTOTSTR, sw_compare.FUNITS), \
      '%s %s Differences %s' % \
        (deltaStr, DOWNTOTSTR, sw_compare.FUNITS), \
      '%s %s' % (sw_compare.FUPSTR, sw_compare.FUNITS), \
      '%s %s Difference %s' % \
        (deltaStr, sw_compare.FUPSTR, sw_compare.FUNITS), \
      '%s %s' % (sw_compare.DIFSTR, sw_compare.FUNITS), \
      '%s %s Difference %s' % \
        (deltaStr, sw_compare.DIFSTR, sw_compare.FUNITS)]

    # if bands is not specified, cover all bands
    bandList = range(nBands) if self.bands is None else self.bands-1
    if self.broadBand: bandList = [nBands]
    for iBand in bandList:
      # make 1 PDF per band
      if self.broadBand:
        bandStr = 'Broadband'
        outFile = '%s/%s_broadband.pdf' % \
          (self.outDir, self.profPrefix)
        wnRange1 = refDict['band_lims_wvn'][0,:][0]
        wnRange2 = refDict['band_lims_wvn'][-1,:][1]
        wnRange = [wnRange1, wnRange2]
      else:
        bandStr = 'Band %d' % (iBand+1)
        outFile = '%s/%s_%02d.pdf' % \
          (self.outDir, self.profPrefix, iBand+1)
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
          (self.atmType, colStr, bandStr, wnRange[0], wnRange[1])

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
          if self.broadBand:
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
            ax.axhline(self.tPauseP, color='k', linestyle=':')
            ax.axvline(0, color='k', linestyle=':')

            # plot dir down and total down together
            if ctr in [1, 2]: plot.plot(dAbscissaDif, tOrdinate, 'k--')
          else:
            # test and reference profiles
            plot.subplot(3, 2, ctr)
            plot.plot(tAbscissa, tOrdinate, 'b')
            plot.plot(rAbscissa, rOrdinate, 'r')
            plot.ylabel('Pressure [mbar]')

            # plot dir down and total down together
            if ctr in [1, 2]:
              plot.plot(tAbscissaDif, tOrdinate, 'b--')
              plot.plot(rAbscissaDif, rOrdinate, 'r--')
              plot.legend(\
                ['Test %s' % DOWNTOTSTR, 'Ref %s' % DOWNTOTSTR, \
                 'Test %s' % DOWNDIRSTR, 'Ref %s' %DOWNDIRSTR], \
                loc='upper left', numpoints=1, \
                prop=sw_compare.font_prop)
            else:
              plot.legend(['Test', 'Reference'], loc='best', \
                numpoints=1, prop=sw_compare.font_prop)
            # endif ctr

          # end % 2

          yRange = [max(rOrdinate), min(rOrdinate)]
          plot.ylim(yRange)
          plot.xlabel(xTitles[ctr-1])
          plot.title(pTitles[ctr-1])
          plot.xticks(rotation=15)
          if self.yLog: plot.semilogy()

        # end ctr loop

        # write page for column
        pdf.savefig()
        plot.close()
      # end iCol loop

      pdf.close()
      print '%s completed' % outFile

      print time.clock() - t1
      if self.broadBand:
        print 'Processed broadband'
      else:
        print 'Processed band %d' % (iBand+1)
    # end band loop

    print 'Profile plotting completed'

  # end profCHARTS()

# end refTestFluxCompare

if __name__ == '__main__':
  parser = argparse.ArgumentParser(\
    description='Generate comparative plots for calculations ' + \
    'from two radiative transfer models.')
  parser.add_argument('--config_file', type=str, \
    default='rrtmgp_charts_config_garand-key.ini', \
    help='Path to configuration file that contains the values ' + \
    'to keyword arguments.')
  parser.add_argument('-prof', '--plot_profiles', \
    action='store_true', \
    help='Generates an N-band page PDF, each with a 3x2 array of ' + \
    'plots that compare the upward and downwward fluxes from ' + \
    'the reference and test models. Also plots heating rate for ' + \
    'each model.')
  parser.add_argument('-stats', '--plot_stats', action='store_true', \
    help='Generates a single page with RRTMGP-CHARTS vs. CHARTS ' + \
    'statistics for TOA upwelling, SFC direct downwelling, ' + \
    'SFC total downwelling, Max Net flux cases as well as ' + \
    'Max Tropospheric and Stratospheric HR.')
  parser.add_argument('--tropopause_pressure', type=float, \
    default=np.exp(4.6), \
    help='Pressure threshold that separates troposphere and ' + \
    'stratosphere.')
  parser.add_argument('--band', type=int, nargs='+', \
    help='Number of band to plot. Default (band=None) is all bands.')
  parser.add_argument('--log_y', action='store_true', \
    help='Generate a semilog-y plot.')
  parser.add_argument('--broadband', action='store_true', \
    help='Generate a broadband plot for profiles (stats process ' + \
    'every band and broadband relatively quickly). WARNING: ' + \
    'broadband supercedes band, so if this is set, the only ' + \
    'result will be a broadband figure.')
  parser.add_argument('--out_dir', type=str, default=os.getcwd(), \
    help='Path to output directory in which PDF files are saved.')
  args = parser.parse_args()

  conFile = args.config_file; utils.file_check(conFile)
  figOutDir = args.out_dir
  utils.file_check(figOutDir)

  configDict = sw_compare.parseConfig(conFile)
  plotObj = refTestFluxCompare(configDict, outDir=figOutDir, \
    tPauseP=args.tropopause_pressure, bands=args.band, \
    yLog=args.log_y, broadband=args.broadband)

  if args.plot_stats: plotObj.statPDF()
  if args.plot_profiles: plotObj.profPDF()

# end main()

