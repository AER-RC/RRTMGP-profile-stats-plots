#!/usr/bin/env python
import os
import sys
import argparse
import netCDF4 as nc
import numpy as np


def getVars(ncFile):
    print 'Reading %s' % ncFile
    ncObj = nc.Dataset(ncFile, 'r')

    outDict = {}
    for var in ncObj.variables:
        outDict[var] = np.array(ncObj.variables[var])

    return outDict


def main():
    parser = argparse.ArgumentParser(\
        description='Read and parse gas_components file')
    parser.add_argument('--reference_file', type=str,
                        default='gas_components.nc',
                        help='Full path to tau gas components file')
    parser.add_argument('--test_file', type=str,
                        default='gas_components.nc',
                        help='Full path to tau gas components file')
    args = parser.parse_args()

    refDict = getVars(args.reference_file)
    testDict = getVars(args.test_file)

if __name__ == '__main__':

    main()
