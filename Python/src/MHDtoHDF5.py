#!/usr/bin/python
'''
Created on Jan 10, 2018

@author: christian
'''
import os
import numpy
from scipy import io
import h5py
import string
import itertools
from optparse import OptionParser
from volumeFileInterface import mhd
import std_typedefs

if __name__ == '__main__':
    optionParser = OptionParser("usage: %prog -i <input matrix file> -o <output file/folder for MHDs> [-n <field name inside mat file>]")
    optionParser.add_option("-i","--input",action="store",dest="input",default="",help="input path to matrix file")
    optionParser.add_option("-o","--output",action="store",dest="output",default="",help="output file (2D/3D) or folder (4D+) for MHD files")
    optionParser.add_option("-f","--fieldName",action="store",dest="fieldName",default="data",help="optional: field name to be extracted from .mat file")
    (options,args) = optionParser.parse_args()
    
    fieldName=options.fieldName
    matNumPyArray = numpy.zeros((1,1,1))
    
    fileReader = mhd.MhdFile()
    fileReader.setFilename(options.input)
    fileReader.readFile()
    nDims = fileReader.getNumberOfDimensions()
    #print nDims
    dim = fileReader.GetDimensions(nDims)
    #print dim
    datatype = fileReader.getDatatype()
    #print datatype
    fsize = 1
    for i in range(0,nDims):
        fsize*=dim[i]
    #print fsize

    dArray = numpy.zeros(fsize)
    if (datatype==std_typedefs.CHAR):
        dArray=fileReader.getDataAsChar(fsize)
    elif (datatype==std_typedefs.UCHAR):
        dArray=fileReader.getDataAsUChar(fsize)
    elif (datatype==std_typedefs.SHORT):
        dArray=fileReader.getDataAsShort(fsize)
    elif (datatype==std_typedefs.USHORT):
        dArray=fileReader.getDataAsUShort(fsize)
    elif (datatype==std_typedefs.INT):
        dArray=fileReader.getDataAsInt(fsize)
    elif (datatype==std_typedefs.UINT):
        dArray=fileReader.getDataAsUInt(fsize)
    elif (datatype==std_typedefs.LONG):
        dArray=fileReader.getDataAsLong(fsize)
    elif (datatype==std_typedefs.ULONG):
        dArray=fileReader.getDataAsULong(fsize)
    elif (datatype==std_typedefs.FLOAT):
        dArray=fileReader.getDataAsFloat(fsize)
    elif (datatype==std_typedefs.DOUBLE):
        dArray=fileReader.getDataAsDouble(fsize)

    if(nDims==2):
        matNumPyArray = numpy.reshape(dArray, (dim[0], dim[1]), 'F') #
    elif (nDims==3):
        matNumPyArray = numpy.reshape(dArray, (dim[0], dim[1], dim[2]), 'F') #
    elif (nDims==4):
        matNumPyArray = numpy.reshape(dArray, (dim[0], dim[1], dim[2], dim[3]), 'F') #
    
    #io.savemat(raw_mat_path, {'data':raw_data})
    f = h5py.File(options.output,'w')
    data_group = f.create_group(fieldName)
    data_hdf5 = matNumPyArray.transpose();
    data_group.create_dataset('value', data=data_hdf5);
    f.close()
    del data_hdf5