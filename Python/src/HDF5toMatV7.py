'''
Created on Jan 09, 2018

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
    except (NotImplementedError,ValueError) as e:
        # for hdf5 matlab files
        f = h5py.File(options.input,'r')
        #matNumPyArray = numpy.array(f[fieldName]['value'], order='F').transpose()
        matNumPyArray = numpy.array(f[fieldName], order='F').transpose()
        f.close()

    io.savemat(options.output, {"data":matNumPyArray}, do_compression=True)