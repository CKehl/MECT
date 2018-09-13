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
    optionParser = OptionParser("usage: %prog -i <input matrix file> -o <output file/folder for MHDs [-n <field name inside mat file>]")
    optionParser.add_option("-i","--input",action="store",dest="input",default="",help="input path to matrix file")
    optionParser.add_option("-o","--output",action="store",dest="output",default="",help="output file (2D/3D) or folder (4D+) for MHD files")
    optionParser.add_option("-f","--fieldName",action="store",dest="fieldName",default="data",help="optional: field name to be extracted from .mat file")
    (options,args) = optionParser.parse_args()
    
    fieldName=options.fieldName
    matNumPyArray = numpy.zeros((1,1,1))
    try:
        # for matlab files version <=7.2
        matNumPyArray = io.loadmat(options.input)[fieldName]
    except:
        try:
            # for hdf5 matlab files
            f = h5py.File(options.input,'r')
            matNumPyArray = numpy.array(f[fieldName]['value'], order='F').transpose()
            f.close()
        except:
            print "Error loading mat file - as v4-v7, and as hdf5. EXITING."
            exit(1)

    #f = h5py.File(options.input,'r')
    #matNumPyArray = numpy.array(f[fieldName]['value'], order='F').transpose()
    #f.close()

    matNumPyArray = numpy.squeeze(matNumPyArray)
    numDims = len(matNumPyArray.shape)
    
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
        sliceFile.setDataAsDouble(selectedArray.astype(numpy.double))
        #sliceFile.setDataAsInt(selectedArray.astype(numpy.int32))
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
            sliceFile.setDataAsDouble(selectedArray.astype(numpy.double))
            #sliceFile.setDataAsShort(selectedArray.astype(numpy.int16))
            sliceFile.setFilename(fname)
            sliceFile.writeFile()

        #fname = options.output + os.path.sep + "image_4D.mhd"
        #dimArray = numpy.array([matNumPyArray.shape[0], matNumPyArray.shape[1], matNumPyArray.shape[2], matNumPyArray.shape[3]], dtype=numpy.uintc)
        #spaceArray = numpy.ones(4, dtype=numpy.float32)
        #arrayLen = matNumPyArray.shape[0]*matNumPyArray.shape[1]*matNumPyArray.shape[2]*matNumPyArray.shape[3]
        #selectedArray = numpy.reshape(matNumPyArray, arrayLen, 'F')
        #sliceFile = mhd.MhdFile()
        #sliceFile.SetDimensions(dimArray)
        #sliceFile.SetSpacing(spaceArray)
        #sliceFile.setDataAsDouble(selectedArray.astype(numpy.double))
        #sliceFile.setFilename(fname)
        #sliceFile.writeFile()
    else:
        exit()