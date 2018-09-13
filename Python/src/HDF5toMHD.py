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

if __name__ == '__main__':
    optionParser = OptionParser("usage: %prog -i <input matrix file> -o <output file/folder for MHDs> [-n <field name inside mat file>]")
    optionParser.add_option("-i","--input",action="store",dest="input",default="",help="input path to matrix file")
    optionParser.add_option("-o","--output",action="store",dest="output",default="",help="output file (2D/3D) or folder (4D+) for MHD files")
    optionParser.add_option("-f","--fieldName",action="store",dest="fieldName",default="data",help="optional: field name to be extracted from .mat file")
    optionParser.add_option("-s","--squeeze",action="store_true",dest="squeeze",help="tell the tool to squeeze the input data, removing ghost dimensions")
    (options,args) = optionParser.parse_args()
    
    fieldName=options.fieldName
    matNumPyArray = numpy.zeros((1,1,1))
    #try:
    #    # for matlab files version <=7.2
    #    matNumPyArray = io.loadmat(options.input)[fieldName]
    #except:
    #    # for hdf5 matlab files
    #    f = h5py.File(options.input,'r')
    #    matNumPyArray = numpy.array(f[fieldName]['value'], order='F').transpose()
    #    f.close()

    f = h5py.File(options.input,'r')
    matNumPyArray = numpy.array(f[fieldName]['value'], order='F').transpose()
    f.close()

    if options.squeeze:
        matNumPyArray = numpy.squeeze(matNumPyArray)
    numDims = len(matNumPyArray.shape)
    arrayType = matNumPyArray.dtype
    #print "Matrix type: %s." % (str(arrayType))

    if(numDims==2):
        #fname = options.output + os.path.sep + "image.mhd"
        fname = options.output
        dimArray = numpy.array([matNumPyArray.shape[0], matNumPyArray.shape[1]], dtype=numpy.uintc)
        spaceArray = numpy.ones(3, dtype=numpy.float32)
        arrayLen = matNumPyArray.shape[0]*matNumPyArray.shape[1]
        selectedArray = numpy.reshape(matNumPyArray[:,:], arrayLen, 'F')
        sliceFile = mhd.MhdFile()
        sliceFile.SetDimensions(dimArray)
        sliceFile.SetSpacing(spaceArray)
        #sliceFile.setDataAsDouble(selectedArray.astype(numpy.double))
        if (arrayType == numpy.byte) or (arrayType == numpy.int8):
            sliceFile.setDataAsChar(selectedArray)
        elif (arrayType == numpy.ubyte) or (arrayType == numpy.uint8):
            sliceFile.setDataAsUChar(selectedArray)
        elif (arrayType == numpy.short) or (arrayType == numpy.int16):
            sliceFile.setDataAsShort(selectedArray)
        elif (arrayType == numpy.ushort) or (arrayType == numpy.uint16):
            sliceFile.setDataAsUShort(selectedArray)
        elif (arrayType == numpy.intc) or (arrayType == numpy.int32):
            sliceFile.setDataAsInt(selectedArray)
        elif (arrayType == numpy.uintc) or (arrayType == numpy.uint32):
            sliceFile.setDataAsUInt(selectedArray)
        elif (arrayType == numpy.int_) or (arrayType == numpy.int64):   # this is "long" python-style
            sliceFile.setDataAsLong(selectedArray)
        elif (arrayType == numpy.uint) or (arrayType == numpy.uint64):   # this is "ulong" python-style
            sliceFile.setDataAsULong(selectedArray)
        elif (arrayType == numpy.single) or (arrayType == numpy.float32):
            sliceFile.setDataAsFloat(selectedArray.astype(numpy.float32))
        elif (arrayType == numpy.float_) or (arrayType==numpy.float) or (arrayType==numpy.double) or (arrayType == numpy.float64):
            sliceFile.setDataAsDouble(selectedArray.astype(numpy.double))
        sliceFile.setFilename(fname)
        sliceFile.writeFile()
    elif(numDims==3):
        #fname = options.output + os.path.sep + "image.mhd"
        fname = options.output
        dimArray = numpy.array([matNumPyArray.shape[0], matNumPyArray.shape[1], matNumPyArray.shape[2]], dtype=numpy.uintc)
        spaceArray = numpy.ones(3, dtype=numpy.float32)
        arrayLen = matNumPyArray.shape[0]*matNumPyArray.shape[1]*matNumPyArray.shape[2]
        selectedArray = numpy.reshape(matNumPyArray[:,:,:], arrayLen, 'F')
        sliceFile = mhd.MhdFile()
        sliceFile.SetDimensions(dimArray)
        sliceFile.SetSpacing(spaceArray)
        #sliceFile.setDataAsDouble(selectedArray.astype(numpy.double))
        #sliceFile.setDataAsInt(selectedArray.astype(numpy.int32))
        if (arrayType == numpy.byte) or (arrayType == numpy.int8):
            sliceFile.setDataAsChar(selectedArray.astype(numpy.byte))
        elif (arrayType == numpy.ubyte) or (arrayType == numpy.uint8):
            sliceFile.setDataAsUChar(selectedArray.astype(numpy.ubyte))
        elif (arrayType == numpy.short) or (arrayType == numpy.int16):
            sliceFile.setDataAsShort(selectedArray.astype(numpy.short))
        elif (arrayType == numpy.ushort) or (arrayType == numpy.uint16):
            sliceFile.setDataAsUShort(selectedArray.astype(numpy.ushort))
        elif (arrayType == numpy.intc) or (arrayType == numpy.int32):
            sliceFile.setDataAsInt(selectedArray.astype(numpy.intc))
        elif (arrayType == numpy.uintc) or (arrayType == numpy.uint32):
            sliceFile.setDataAsUInt(selectedArray.astype(numpy.uintc))
        elif (arrayType == numpy.int_) or (arrayType == numpy.int64):   # this is "long" python-style
            sliceFile.setDataAsLong(selectedArray)
        elif (arrayType == numpy.uint) or (arrayType == numpy.uint64):   # this is "ulong" python-style
            sliceFile.setDataAsULong(selectedArray)
        elif (arrayType == numpy.single) or (arrayType == numpy.float32):
            sliceFile.setDataAsFloat(selectedArray.astype(numpy.float32))
        elif (arrayType == numpy.float_) or (arrayType==numpy.float) or (arrayType==numpy.double) or (arrayType == numpy.float64):
            sliceFile.setDataAsDouble(selectedArray.astype(numpy.double))
        sliceFile.setFilename(fname)
        sliceFile.writeFile()
    elif(numDims==4):
        if(os.path.exists(options.output)==False):
            os.mkdir(options.output)
        
        numChannels = matNumPyArray.shape[3]
        for channelIndex in itertools.islice(itertools.count(),0,numChannels):
            fname = options.output + os.path.sep + ("image_channel%04d.mhd" % (channelIndex))
            dimArray = numpy.array([matNumPyArray.shape[0], matNumPyArray.shape[1], matNumPyArray.shape[2]], dtype=numpy.uintc)
            spaceArray = numpy.ones(3, dtype=numpy.float32)
            arrayLen = matNumPyArray.shape[0]*matNumPyArray.shape[1]*matNumPyArray.shape[2]
            selectedArray = numpy.reshape(matNumPyArray[:,:,:,channelIndex], arrayLen, 'F')
            sliceFile = mhd.MhdFile()
            sliceFile.SetDimensions(dimArray)
            sliceFile.SetSpacing(spaceArray)
            #sliceFile.setDataAsDouble(selectedArray.astype(numpy.double))
            #sliceFile.setDataAsShort(selectedArray.astype(numpy.int16))
            if (arrayType == numpy.byte) or (arrayType == numpy.int8):
                sliceFile.setDataAsChar(selectedArray.astype(numpy.byte))
            elif (arrayType == numpy.ubyte) or (arrayType == numpy.uint8):
                sliceFile.setDataAsUChar(selectedArray.astype(numpy.ubyte))
            elif (arrayType == numpy.short) or (arrayType == numpy.int16):
                sliceFile.setDataAsShort(selectedArray.astype(numpy.short))
            elif (arrayType == numpy.ushort) or (arrayType == numpy.uint16):
                sliceFile.setDataAsUShort(selectedArray.astype(numpy.ushort))
            elif (arrayType == numpy.intc) or (arrayType == numpy.int32):
                sliceFile.setDataAsInt(selectedArray.astype(numpy.intc))
            elif (arrayType == numpy.uintc) or (arrayType == numpy.uint32):
                sliceFile.setDataAsUInt(selectedArray.astype(numpy.uintc))
            elif (arrayType == numpy.int_) or (arrayType == numpy.int64):   # this is "long" python-style
                sliceFile.setDataAsLong(selectedArray)
            elif (arrayType == numpy.uint) or (arrayType == numpy.uint64):   # this is "ulong" python-style
                sliceFile.setDataAsULong(selectedArray)
            elif (arrayType == numpy.single) or (arrayType == numpy.float32):
                sliceFile.setDataAsFloat(selectedArray.astype(numpy.float32))
            elif (arrayType == numpy.float_) or (arrayType==numpy.float) or (arrayType==numpy.double) or (arrayType == numpy.float64):
                sliceFile.setDataAsDouble(selectedArray.astype(numpy.double))
            sliceFile.setFilename(fname)
            sliceFile.writeFile()

        fname = options.output + os.path.sep + "image_4D.mhd"
        dimArray = numpy.array([matNumPyArray.shape[0], matNumPyArray.shape[1], matNumPyArray.shape[2], matNumPyArray.shape[3]], dtype=numpy.uintc)
        spaceArray = numpy.ones(4, dtype=numpy.float32)
        arrayLen = matNumPyArray.shape[0]*matNumPyArray.shape[1]*matNumPyArray.shape[2]*matNumPyArray.shape[3]
        selectedArray = numpy.reshape(matNumPyArray, arrayLen, 'F')
        sliceFile = mhd.MhdFile()
        sliceFile.SetDimensions(dimArray)
        sliceFile.SetSpacing(spaceArray)
        sliceFile.setDataAsDouble(selectedArray.astype(numpy.double))
        sliceFile.setFilename(fname)
        sliceFile.writeFile()
    else:
        exit()