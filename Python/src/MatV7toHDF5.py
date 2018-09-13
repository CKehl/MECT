#!/usr/bin/python
'''
Created on Apr 16, 2018

@author: christian
'''
import os
import numpy
from scipy import io
import h5py
import string
import itertools
from optparse import OptionParser
from volumeFileInterface import ini

if __name__ == '__main__':
    optionParser = OptionParser("usage: %prog -i <input matrix file> -o <output matrix file> [-n <field name inside mat file>]")
    optionParser.add_option("-i","--input",action="store",dest="input",default="",help="input path to matrix file")
    optionParser.add_option("-o","--output",action="store",dest="output",default="",help="output path to target matrix file")
    optionParser.add_option("-f","--fieldName",action="store",dest="fieldName",default="data",help="optional: field name to be extracted from .mat file")
    (options,args) = optionParser.parse_args()
    
    fieldName=options.fieldName
    matNumPyArray = numpy.zeros((1,1,1))
    try:
        # for matlab files version <=7.2
        matNumPyArray = io.loadmat(options.input)[fieldName]
    except ValueError as e:
        # for hdf5 matlab files
        e.printTraceback()

    f = h5py.File(options.output,'w')
    data_group = f.create_group('data')
    data_hdf5 = matNumPyArray.transpose();
    data_group.create_dataset('value', data=data_hdf5);
    f.close()
    del data_hdf5