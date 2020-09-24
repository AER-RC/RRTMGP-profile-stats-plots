#!/usr/bin/env python

import os, sys, argparse, glob
import ConfigParser
import numpy as np

from multiprocessing import Pool, cpu_count

# should be in working directory
import utils
import SW_LBLRTM_RRTMGP_compare as sw_compare

def configToDict(inFile, tPauseP=np.exp(4.6), diffuse=False, \
  bands=None, broadband=False, yLog=False):
  """
  Convert configuration file into a dictionary for input into either 
  poolStats() or poolProfs()

  Input
    inFile -- string, .ini file to parse and save into outDict

  Output
    outDict -- dictionary with the following keys:
      ref_name -- string, reference model name
      test_name -- string, test model name
      ref_desc -- string, reference model description
      test_desc -- string, test model description
      ref_file -- string, path to reference model netCDF
      test_file -- string, path to test model netCDF
      ref_force_file -- string, path to reference model forcing netCDF
      test_force_file -- string, path to test model forcing netCDF
      stat_prefix -- string, prefix for statistics PDF filename 
      profile_prefix --string, prefix for profile PDF filenames
      x_label -- string, label for x axes in plots (stats only)
      y_label -- string, label for y axes in profile and stats plots
      atm_type -- string, type of atmosphere provided by user into 
        LBLRTM to produce ref_file and test_file
      tpause_pres -- float, pressure at tropopause (tPauseP value)
      ini_file -- string, inFile
      diffuse -- boolean, generate plots for diffuse flux rather than 
        of direct beam flux (diffuse keyword value)
      forcing -- boolean, plot differences in ref and test forcing

      The dictionary will also contain the following fields, but they
      will probably only be used with the profile plots (poolProfs())

      bands -- int array, zero-offset list of bands to process (only 
        included if profile plotting is being done)
      broadband -- boolean, for profiles, make a page with broadband
        parameter differences
      log_y -- boolean, plot flux and HR profiles on a semilog plot
        (pressure is the log axis)

  Keywords
    tPauseP -- float, pressure at tropopause
    diffuse -- boolean, generate plots for diffuse flux rather than 
      of direct beam flux
    bands -- int array, zero-offset list of bands to process (only 
      necessary with profile plotting)
    broadband -- plot profiles for SW broadband
    yLog -- plot profiles on a log scale on the y axis (pressure)
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
  utils.file_check(rForceFile); utils.file_check(tForceFile)

  # forcing keyword is True if both ref and test forcing files are 
  # provided (still a lot of work to do here)
  forcing = True if rForceFile and tForceFile else False
#  if forcing:
#    cRefName, cTestName = sw_compare.compare.forcingDiff(\
#      cRefName, cTestName, \
#      inDict['ref_force_file'], inDict['test_force_file'])
#  # endif forcing

  statPrefix = cParse.get('Filename Params', 'stats_prefix')
  profPrefix = cParse.get('Filename Params', 'profiles_prefix')

  xt = cRefName
  yt = '%s - %s' % (cTestName, cRefName)

  outDict = {'ref_name': cRefName, 'test_name': cTestName, \
    'ref_desc': cRefMD, 'test_desc': cTestMD, \
    'ref_file': refFile, 'test_file': testFile, \
    'ref_force_file': rForceFile, 'test_force_file': tForceFile, \
    'stat_prefix': statPrefix, 'profile_prefix': profPrefix, \
    'x_label': xt, 'y_label': yt, 'atm_type': aType, \
    'tpause_pres': tPauseP, 'ini_file': inFile, \
    'diffuse': diffuse, 'forcing': forcing}

  outDict['bands'] = bands
  outDict['broadband'] = broadband
  outDict['log_y'] = yLog

  return outDict
# end configToDict()

def poolProfs(inDict):
  """
  Single function call to profile plotting routine profPDF() in 
  sw_compare module, designed to work with Pool.map(), which is 
  called in main().

  Input
    inDict -- dictionary generated with configToDict()

  Output
    None, plots are generated with sw_compare.profPDF():

  """

  # placeholder code
  """
  # plot the test and reference upward and downward flux together
  # with the heating rates (and differences for each)
  inBand = None if args.band is None else np.array(args.band)-1

  if inBand is None:
    profPDFs(refFile, testFile, yt, tPauseP=pTrop, \
      prefix=profPrefix, atmType=aType, inBand=inBand, \
      logy=args.log_y, broadOnly=args.broad_only, \
      diffuse=doDif, outDir=figOutDir)
  else:
    for iBand in inBand:
      profPDFs(refFile, testFile, yt, tPauseP=pTrop, \
        prefix=profPrefix, atmType=aType, inBand=iBand, \
        logy=args.log_y, diffuse=doDif, outDir=figOutDir)
    # end iBand loop

    # for specified bands AND broadband
    if args.broad_only:
      profPDFs(refFile, testFile, yt, tPauseP=pTrop, \
        prefix=profPrefix, atmType=aType, \
        logy=args.log_y, broadOnly=args.broad_only, \
        diffuse=doDif, outDir=figOutDir)
    # end broadband plot
  # end inBand
  """

  if not os.path.exists('figures'): os.mkdir('figures')
  figOutDir = 'figures/%s' % \
    os.path.basename(inDict['ini_file'])[:-4]
  if not os.path.exists(figOutDir): os.mkdir(figOutDir)

  # single-band profiles
  inBand = inDict['bands']
  if inBand is None:
    sw_compare.profPDFs(inDict['ref_file'], inDict['test_file'], \
      inDict['y_label'], tPauseP=inDict['tpause_pres'], \
      prefix=inDict['profile_prefix'], atmType=inDict['atm_type'], \
      inBand=inBand, logy=inDict['log_y'], \
      diffuse=inDict['diffuse'], outDir=figOutDir)
  else:
    for iBand in inBand:
      sw_compare.profPDFs(inDict['ref_file'], inDict['test_file'], \
        inDict['y_label'], tPauseP=inDict['tpause_pres'], \
        prefix=inDict['profile_prefix'], atmType=inDict['atm_type'], \
        inBand=iBand, logy=inDict['log_y'], \
        diffuse=inDict['diffuse'], outDir=figOutDir)
    # end iBand loop
  # endif inBand

  # broadband profiles (if specified)
  if inDict['broadband']:
    sw_compare.profPDFs(inDict['ref_file'], inDict['test_file'], \
      inDict['y_label'], tPauseP=inDict['tpause_pres'], \
      prefix=inDict['profile_prefix'], atmType=inDict['atm_type'], \
      logy=inDict['log_y'], broadOnly=inDict['broadband'], \
      diffuse=inDict['diffuse'], outDir=figOutDir)

  return True
# end poolProfs()

def poolStats(inDict):
  """
  Single function call to profile plotting routine statPDF() in 
  sw_compare module, designed to work with Pool.map(), which is 
  called in main().

  Input
    inDict -- dictionary generated with configToDict()

  Output
    None, plots are generated with sw_compare.statPDF() and saved into
      the figures/ directory (which is generated if it does not 
      already exist in the working directory) under the appropriate 
      subdirectory (which is a string extracted from the .ini file)

  """

  if not os.path.exists('figures'): os.mkdir('figures')
  figOutDir = 'figures/%s' % \
    os.path.basename(inDict['ini_file'])[:-4]
  if not os.path.exists(figOutDir): os.mkdir(figOutDir)

  sw_compare.statPDF(inDict['ref_file'], inDict['test_file'], \
    tPauseP=inDict['tpause_pres'], atmType=inDict['atm_type'], \
    xTitle=inDict['x_label'], yTitle=inDict['y_label'], \
    forcing=inDict['forcing'], prefix=inDict['stat_prefix'], \
    diffuse=inDict['diffuse'], outDir=figOutDir)

  return True
# end poolStats()

if __name__ == '__main__':
  parser = argparse.ArgumentParser(\
    description='Multithread the SW LBLRTM-RRTMGP plot (profiles ' + \
    'or stats) process. This should not be used with the LW.')
  parser.add_argument('--stats', action='store_true', \
    help='Plot test-reference difference stats.')
  parser.add_argument('--profiles', action='store_true', \
    help='Plot test-reference flux and difference profiles.')
  parser.add_argument('--indir', type=str, default='config_files/', \
    help='Directory with .ini configuration files to use as input ' + \
    'into plotting functions.')
  parser.add_argument('--tropopause_pressure', type=float, \
    default=np.exp(4.6), \
    help='Pressure threshold that separates troposphere and ' + \
    'stratosphere.')
  parser.add_argument('--diffuse', action='store_true', \
    help='Plot diffuse downwelling flux array instead of flux ' + \
    'from direct beam.')
  parser.add_argument('--bands', type=int, nargs='+', \
    help='Band numbers to plot. Default (band=None) is all bands.')
  parser.add_argument('--broadband', action='store_true', \
    help='Generate a broadband PDF.')
  parser.add_argument('--log_y', action='store_true', \
    help='Generate a semilog-y plot.')
  parser.add_argument('--cores', type=int, default=4, \
    help='Number of cores to use in parallelization.')
  args = parser.parse_args()

  stats = args.stats; profs = args.profiles
  if not stats and not profs:
    sys.exit('Please use --stats or --profs in arguments.')

  # grab all configuration (.ini) files
  inDir = args.indir; utils.file_check(inDir)
  iniFiles = sorted(glob.glob('%s/*.ini' % inDir))
  pTrop = args.tropopause_pressure; difFlux = args.diffuse

  # for each ini file, generate a dictionary to be used in 
  # pool functions
  inDicts = []
  for iniFile in iniFiles:
    if stats:
      cDict = configToDict(iniFile, tPauseP=pTrop, diffuse=difFlux)

    if profs:
      doBands = True if args.bands else False
      if doBands:
        # generate a zero-offset band array from the band list 
        # specified by the user
        bandArr = np.array(args.bands)-1
      else:
        bandArr = None
      # endif doBands

      cDict = configToDict(iniFile, tPauseP=pTrop, diffuse=difFlux, \
        bands=bandArr, broadband=args.broadband, yLog=args.log_y)
    # endif profs

    inDicts.append(cDict)
  # end loop over ini files

  # always leave at least one extra core for the user
  nCores = args.cores
  totCores = cpu_count()
  nCores = nCores if nCores < totCores else totCores-1

  p = Pool(nCores)
  if args.stats: p.map(poolStats, inDicts)
  if args.profiles: p.map(poolProfs, inDicts)

# end main()

