#!/usr/bin/env python

import os, sys, argparse

# should be in working directory
import utils

parser = argparse.ArgumentParser(\
  description='Given a FORTRAN 90 file with ' + \
  'BND_LIMITS_WAVENUMBER_LW and BND_LIMITS_WAVENUMBER_SW arrays, ' + \
  'generate an LNFL TAPE5 file for each. The output files follow ' + \
  'the convention: TAPE5.ww_bb where "ww" = [lw, sw] and bb is ' + \
  'the 2-digit band number.')
parser.add_argument('--infile', type=str, \
  default='mo_band_wavenumbers.f90', \
  help='Fortran 90 file with wavenumber limit arrays.')
parser.add_argument('--nmolecules', type=int, default=7, \
  help='The number of molecules to use in LNFL.')
parser.add_argument('--wn_pad', type=float, default=50.0, \
  help='The number of wavenumbers with which to pad on both ' + \
  'sides of the band limits. This is for the TAPE5 band range.')
args = parser.parse_args()

inFile = args.infile; utils.file_check(inFile)
inDat = open(inFile).read().splitlines()

lwLims, swLims = [], []
for iLine, line in enumerate(inDat):
  # change the lims ID based on band domain
  # lims.append() will append to lwLims
  if 'BND_LIMITS_WAVENUMBER_LW' in line: lims = lwLims
  if 'BND_LIMITS_WAVENUMBER_SW' in line: lims = swLims

  # sanity check: does this method work?
  """
  try:
    print iLine, id(lims), id(lwLims), id(swLims)
  except:
    continue
  """

  # append each band's limits as a tuple to the appropriate list
  if '_wp' in line:
    # first remove unnecessary characters
    split = line.split(',')
    wn1 = split[0].replace('_wp', '')
    wn2 = split[1].replace('_wp', '').replace('/', '').replace(')', '')
    lims.append( (float(wn1), float(wn2)) )
  # endif _wp
# print line

# sanity check: did i append correctly?
"""
for lw in lwLims: print lw
print
for sw in swLims: print sw
"""

# common ASCII lines for LW and SW
t5head = '$ '; t5head += '\n'
t5mol = '1' * args.nmolecules; t5mol += '\n'
t5tail = '%' * 20; t5tail += '\n'

wnPad = args.wn_pad

# loop over LW bands and write associated TAPE5
for iBand, lw in enumerate(lwLims):
  outFile = 'TAPE5.lw_%02d' % (iBand+1)
  outFP = open(outFile, 'w')
  outFP.write(t5head)

  outLine = '%10.3f%10.3f\n' % (lw[0]-wnPad, lw[1]+wnPad)
  outFP.write(outLine)
  outFP.write(t5mol)
  outFP.write(t5tail)

  outFP.close()
  print outFile
# end lw loop

# loop over LW bands and write associated TAPE5
for iBand, sw in enumerate(swLims):
  outFile = 'TAPE5.sw_%02d' % (iBand+1)
  outFP = open(outFile, 'w')
  outFP.write(t5head)

  outLine = '%10.3f%10.3f\n' % (sw[0]-wnPad, sw[1]+wnPad)
  outFP.write(outLine)
  outFP.write(t5mol)
  outFP.write(t5tail)

  outFP.close()
  print outFile
# end lw loop

