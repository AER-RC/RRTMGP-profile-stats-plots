#!/usr/bin/env python

# for Python 3 compatibility
from __future__ import print_function

import os, sys, argparse
import subprocess as sub
import shutil

# NOTE: this *should* be flexible enough for LW and SW, LBLRTM and 
# CHARTS, but the immediate need is CHARTS (we have other, Bash, 
# scripts for LBLRTM-RRTMGP validation). should also handle forcing, 
# but does not do that yet

class refTestPlot():
  def __init__(self, inFile, charts=False):
    """
    Starting with a template .ini file, replace fields as necessary 
    (plot labels, file name, revision numbers, etc.), then run 
    comparison plotting routine to generate profiles and/or stats 
    PDF.

    Input
      inFile -- string, .cfg configuration file that contains 
        parameters for the .ini configuration file

    Keywords
      charts -- boolean, if set, uses CHARTS-LBLRTM fluxes as a 
        reference instead of standalone LBLRTM fluxes
    """

    self.file_check(inFile)
    self.inFile = inFile
    self.cfgDict = self.cfgRead()

    # sanity check
    #for key in self.cfgDict.keys(): print(key, self.cfgDict[key])

    self.topDir = os.getcwd()
    self.charts = charts
    self.iniPrefix = 'rrtmgp_charts_config' if charts else \
      'rrtmgp_lblrtm_config'

    # now do stuff
    self.stdout()
    self.iniList, self.plotDirs = self.generateIni()
    self.makePlot()

  # end constructor

  def file_check(self, path):
    """
    Quick check if path exists.  Use before reading a file.
    """
    if not os.path.exists(path):
      sys.exit('Could not find %s, returning' % path)
  # end file_check()

  def cfgRead(self):
    """
    Read in .cfg file and convert to a dictionary (whose fields depend
    on the .cfg convention)
    """

    cfgDat = open(self.inFile).read().splitlines()
    cfgDict = {}

    # first remove comments
    for iLine, line in enumerate(cfgDat):
      if '#' in line: cfgDat[iLine] = line.split('#')[0]
    # end line loop

    for line in cfgDat:
      # we also assume "=" is in the ASCII line and that the 
      # left side is the dictionary key, right side is value
      split = line.split('=')

      # empty lines or ones without 2 fields
      if (len(line) == 0) or (len(split) <= 1): continue

      # remove end quotes and white spaces if necessary
      split1 = split[1].strip()
      if '"' in split1: split1 = split1[1:-1]
      cfgDict[split[0]] = str(split1)

      # scenario strings (store as a list)
      if split[0] == 'scenarios': cfgDict[split[0]] = split1.split()
      if split[0] == 'forcing': cfgDict[split[0]] = split1.split()

    # end line loop

    return cfgDict
  # end cfgRead()

  def stdout(self):
    """
    Print out what scenarios will being plotted
    """

    print()
    print('Scenarios:')
    for scenario in self.cfgDict['scenarios']:
      print('\t%s' % scenario)

    print('Forcing Scenarios:')
    if len(self.cfgDict['forcing']) == 0:
      self.forcing = False
      print('\tNone')
    else:
      self.forcing = True
      for scenario in self.cfgDict['forcing']:
        print('\t%s' % scenario)
    # endif forcing
  # endif stdout

  def generateIni(self):
    """
    Using the dictionary from cfgRead(), generate a .ini file that 
    uses the cfg-specified parameters

    Output
      newIni -- string list of new .ini configuration files that can 
        be used in the plotting scripts
    """

    # some shortcuts (less typing)
    inDict = self.cfgDict
    pDir = inDict['pdir']
    vRRTMGP = inDict['rrtmgpvernum']
    vLBL = inDict['lblvernum']
    vAtm = inDict['atmvernum']
    domain = inDict['lwsw']
    ang = inDict['angles']
    catFluxDir = 'concatenate_fluxes_%s' % domain
    dirLBL = '%s/output_lbl_%s_atm_%s_%s' % \
      (inDict['lbldir'], vLBL, vAtm, domain)
    testModel = 'RRTMGP_%s' % vRRTMGP
    refModel = 'CHARTS' if self.charts else \
      'LBLRTM_%s' % vLBL

    # hard set CHARTS netCDF (for now)
    chartsPath = '/project/rc/rc2/rpernak/RRTMGP/' + \
      'Garand_fluxes/netCDF_alb_0.0/' + \
      'charts-sw-inputs-outputs-clear_alb_0.0.nc'

    # generate .ini file for each scenario
    newIni, plotDirList = [], []
    for run in inDict['scenarios']:
      print('RUNNING...%s' % run)
      testNC = 'rrtmgp-%s-flux-inputs-outputs-%s-all.nc' % \
        (domain, run)
      lblNC = 'lblrtm-%s-flux-inputs-outputs-%s-all.nc' % (ang, run)

      # new .ini file name
      iniFile = '%s_%s.ini' % (self.iniPrefix, run)
      if os.path.exists(iniFile): os.remove(iniFile)

      # .ini template
      iniTemp = '%s/rrtmgp_validations/%s' % \
        (inDict['pythondir'], inDict['ptemplate'])
      self.file_check(iniTemp)
      shutil.copyfile(iniTemp, iniFile)

      # where we'll run the plotting script
      plotDir = '%s/output_rrtmgp_%s_atm_%s_%s/%s/%s' % \
        (pDir, vRRTMGP, vAtm, domain, run, ang)
      self.file_check(plotDir)
      plotDirList.append(plotDir)

      # now loop through fields and replace placeholders (XYZ) with 
      # trial strings
      xyz = 'XYZ'
      iniDat = open(iniFile).read().splitlines()
      outFP = open(iniFile, 'w')
      for iLine, line in enumerate(iniDat):
        # not particularly robust, but for now we only need to 
        # replace a few lines, so this should suffice
        # output file strings
        if 'atmosphere = XYZ' in line:
          line = line.replace(xyz, '%s_%s' % (run, vAtm))
        if 'profiles_prefix = XYZ' in line:
          line = line.replace(xyz, 'profs_%s_%s_%s_%s' % \
            (inDict['y_descrip'], refModel, testModel, run))
        if 'stats_prefix = XYZ' in line:
          line = line.replace(xyz, 'stats_%s_%s_%s' % \
            (refModel, testModel, run))
        if 'stats_csv = XYZ' in line:
          line = line.replace(xyz, '%s_%s_%s_diff.csv' % \
            (refModel, testModel, run))

        # model strings for plot labels
        if 'test_model = XYZ' in line:
          line = line.replace(xyz, testModel)
        if 'reference_model = XYZ' in line:
          line = line.replace(xyz, refModel)
        if 'test_forcing_model = XYZ' in line:
          line = line.replace(xyz, '')
        if 'reference_forcing_model = XYZ' in line:
          line = line.replace(xyz, '')
        
        # model paths
        if 'test_path = XYZ' in line:
          testPath = '%s/%s/%s' % (plotDir, catFluxDir, testNC)
          line = line.replace(xyz, testPath)
        # endif test_path

        if 'reference_path = XYZ' in line:
          if self.charts:
            refPath = str(chartsPath)
          else:
            refPath = '%s/%s/%s/%s/%s' % \
              (dirLBL, run, ang, catFluxDir, lblNC)
          # endif charts

          line = line.replace(xyz, refPath)
        # end reference_path

        if 'test_force_path = XYZ' in line:
          line = line.replace(xyz, '')
        if 'reference_force_path = XYZ' in line:
          line = line.replace(xyz, '')
        
        # replace template line with modified line
        outFP.write('%s\n' % line)
      # end line loop

      # close up shop for run and save .ini file for future reference
      outFP.close()
      newIni.append(iniFile)

    # end run loop

    return newIni, plotDirList

  # end generateIni()

  def makePlot(self):
    """
    Call the validation plotting scripts with new .ini files
    """

    inDict = self.cfgDict
    domain = inDict['lwsw']
    scriptDir = '%s/rrtmgp_validations' % inDict['pythondir']
    bandSplit = inDict['pl_bands'].split()

    for ini, pDir in zip(self.iniList, self.plotDirs):
      print('Plotting for %s' % ini)
      os.rename(ini, '%s/%s' % (pDir, ini) )
      os.chdir(pDir)

      if self.charts:
        # couldn't get subprocess.call() to work, so we're hacking
        # also not doing stats since CHARTS stats are not yet done
        plotScript = '%s/CHARTS_RRTMGP_compare.py' % scriptDir
        cmd = '%s --config_file %s %s %s %s -stats' % \
          (plotScript, ini, inDict['profiles'], inDict['pl_bands'], \
           inDict['y_axis'])
      else:
        plotPrefix = 'SW_' if domain == 'sw' else ''
        plotScript = '%s/%sLBLRTM_RRTMGP_compare.py' % \
          (scriptDir, plotPrefix)
        cmd = '%s --config_file %s %s %s %s' % \
          (plotScript, ini, inDict['profiles'], inDict['pl_bands'], \
           inDict['y_axis'])
        if domain == 'lw': cmd += ' --single_stat'
      # end charts

      call = sub.Popen(cmd, shell=True)
      os.chdir(self.topDir)
    # end ini loop
  # end makePlot()
# end refTestPlot

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('cfg_file', type=str, \
    help='Configuration file (.cfg extenstion) that contains ' + \
    'parameters that will be used in the .ini configuration file.')
  parser.add_argument('--charts', action='store_true', \
    help='Use CHARTS instead of LBLRTM as reference')
  args = parser.parse_args()

  cfgFile = args.cfg_file
  refTest = refTestPlot(cfgFile, charts=args.charts)
  
# endif main()

