#!/usr/bin/python
'''
Created on Mar 8, 2018

@author: christian
'''
import os
import numpy
from scipy import io
import h5py
import string
import itertools
from optparse import OptionParser
from tools import loadFAMSoutput


if __name__ == '__main__':
    optionParser = OptionParser("usage: %prog -i <input matrix file> -o <output file/folder for MHDs [-n <field name inside mat file>]")
    optionParser.add_option("-i","--input",action="store",dest="input",default="",help="input path to matrix file")
    optionParser.add_option("-o","--output",action="store",dest="output",default="",help="output file (2D/3D) or folder (4D+) for MHD files")
    optionParser.add_option("-f","--fieldName",action="store",dest="fieldName",default="data",help="optional: field name to be extracted from .mat file")
    optionParser.add_option("-X","--dimX",action="store",dest="dimX",type=int,default=-1,help="number of bins in dimension 0")
    optionParser.add_option("-Y","--dimY",action="store",dest="dimY",type=int,default=-1,help="number of bins in dimension 1")
    optionParser.add_option("-Z","--dimZ",action="store",dest="dimZ",type=int,default=-1,help="number of bins in dimension 2")
    (options,args) = optionParser.parse_args()
    
    fieldName=options.fieldName
    ndims = 3
    if(options.dimZ==-1):
        ndims = 2
    if (options.dimX==-1) or (options.dimY==-1):
        print "Error - less than 2 dimensions are specified. Exiting ..."
        exit(1)

    dims = []
    dims.append(options.dimX)
    dims.append(options.dimY)
    ashape = (dims[0], dims[1])
    if(ndims > 2):
        dims.append(options.dimZ)
        ashape = (dims[0], dims[1], dims[2])
    
    print dims
    print ashape
    
    matNumPyArray = numpy.zeros((1,1,1), dtype = numpy.ushort)
    
    if ndims==2:
        reader = loadFAMSoutput.loadFAMSoutput()
        reader.load(options.input)
        print "# features: %d" % (reader._nfeatures)
        segments = reader._segments
        segments = numpy.reshape(segments,newshape=ashape,order='F')
        
        f = h5py.File(options.output,'w')
        data_group = f.create_group(options.fieldName)
        segment_data_hdf5 = segments.transpose();
        data_group.create_dataset('value', data=segment_data_hdf5);
        f.close()
    elif ndims==3:
        reader = loadFAMSoutput.loadFAMSoutput()
        reader.load(options.input)
        print "# features: %d" % (reader._nfeatures)
        segments = reader._segments
        segments = numpy.reshape(segments,newshape=ashape,order='F')
        
        f = h5py.File(options.output,'w')
        data_group = f.create_group(options.fieldName)
        segment_data_hdf5 = segments.transpose();
        data_group.create_dataset('value', data=segment_data_hdf5);
        f.close()
    