#!/usr/bin/python
'''
Created on Mar 1, 2018

@author: christian
'''
import os
import numpy
from scipy import io
import h5py
import string
import itertools
from optparse import OptionParser
from volumeFileInterface import fams_ascii
import math
import time
import sys
import traceback

if __name__ == '__main__':
    optionParser = OptionParser("usage: %prog -i <input matrix file> -o <output file/folder for MHDs [-n <field name inside mat file>]")
    optionParser.add_option("-i","--input",action="store",dest="input",default="",help="input path to matrix file")
    optionParser.add_option("-o","--output",action="store",dest="output",default="",help="output file (2D/3D) or folder (4D+) for MHD files")
    optionParser.add_option("-f","--fieldName",action="store",dest="fieldName",default="data",help="optional: field name to be extracted from .mat file")

    optionParser.add_option("-S","--spectral_stride",action="store",dest="spectralStride",type="int",default=-1,help="stride between spectral bins selected for evaluation")
    optionParser.add_option("-s","--train_stride",action="store",dest="trainStride",type="int",default=-1,help="(spatial) stride between samples for the training dataset")
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

    matNumPyArray = numpy.squeeze(matNumPyArray)
    numDims = len(matNumPyArray.shape)
    arrayType = matNumPyArray.dtype
    
    aMax, aMin = matNumPyArray.max(), matNumPyArray.min()
    matNumPyArray = (matNumPyArray - aMin)/(aMax - aMin)
    matNumPyArray = matNumPyArray*65535
    
    spectralStride = options.spectralStride
    trainStride = options.trainStride
    
    selected_channels = [19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75,77]
    if(spectralStride>=0):
        selected_channels=[]
        for binIdx in range(0,60,spectralStride):
            selected_channels.append(19+binIdx)

    multispectral_train = numpy.zeros(matNumPyArray.shape, dtype=matNumPyArray.dtype, order='F')
    if (trainStride>0) and (numDims==4):
        oshape  = matNumPyArray.shape
        nfeatures = len(selected_channels)
        nshape = (int(math.ceil(float(oshape[0])/trainStride)), int(math.ceil(float(oshape[1])/trainStride)), int(math.ceil(float(oshape[2])/trainStride)), nfeatures)
        numSamples = int(math.ceil(float(oshape[0])/trainStride)) * int(math.ceil(float(oshape[1])/trainStride)) * int(math.ceil(float(oshape[2])/trainStride))
        multispectral_train = numpy.zeros( nshape, dtype=numpy.double, order='F' )
        print "multispectral_train shape:" + str(multispectral_train.shape)
        start_ordering = time.time()
        D_index=0
        T_index=0
        #for D_index in itertools.islice(itertools.count(),0,nshape[0]*nshape[1]*nshape[2]):
        for X_index in itertools.islice(itertools.count(),0,oshape[0]):
            for Y_index in itertools.islice(itertools.count(),0,oshape[1]):
                for Z_index in itertools.islice(itertools.count(),0,oshape[2]):
                    for E_index in itertools.islice(itertools.count(),0,len(selected_channels)):
                        E_data_index = selected_channels[E_index]
                        X_train_index = int(math.floor(X_index/trainStride))
                        Y_train_index = int(math.floor(Y_index/trainStride))
                        Z_train_index = int(math.floor(Z_index/trainStride))
                        
                        try:
                            if(math.isnan(matNumPyArray[X_index,Y_index,Z_index,E_index])):
                                matNumPyArray[X_index,Y_index,Z_index,E_data_index]=0
                            if(math.isinf(matNumPyArray[X_index,Y_index,Z_index,E_index])):
                                matNumPyArray[X_index,Y_index,Z_index,E_data_index]=0

                            if (X_index%trainStride == 0) and (Y_index%trainStride == 0) and (Z_index%trainStride == 0):
                                multispectral_train[X_train_index,Y_train_index,Z_train_index,E_index] = matNumPyArray[X_index,Y_index,Z_index,E_data_index]
                        except:
                            print "error with %i, %i, %i\n" % (X_index, Y_index, Z_index)
                            exc_type, exc_value, exc_traceback = sys.exc_info()
                            print repr(traceback.format_exception(exc_type, exc_value, exc_traceback))

        end_ordering = time.time()
        print "Set up the training data in %d seconds" % (end_ordering-start_ordering)
    elif(numDims==4):
        multispectral_train = matNumPyArray[:,:,:,selected_channels]
    else:
        ##multispectral_train = matNumPyArray
        pass

    
    if(numDims==3):
        fname = options.output + os.path.sep + "fams.txt"
        #dimArray = numpy.array([matNumPyArray.shape[0], matNumPyArray.shape[1], matNumPyArray.shape[2]], dtype=numpy.uintc)
        dimArray = numpy.array([matNumPyArray.shape[0], matNumPyArray.shape[1], len(selected_channels)], dtype=numpy.uintc)
        #arrayLen = matNumPyArray.shape[0]*matNumPyArray.shape[1]*matNumPyArray.shape[2]
        arrayLen = matNumPyArray.shape[0]*matNumPyArray.shape[1]*len(selected_channels)
        selectedArray = matNumPyArray[:,:,selected_channels]
        print selectedArray.shape
        selectedArray = numpy.reshape(selectedArray, arrayLen, 'F')
        sliceFile = fams_ascii.FamsFile_ASCII()
        sliceFile.SetDimensions(dimArray)
        sliceFile.setDataorder(fams_ascii.COLUMN_MAJOR)
        sliceFile.setDataByUShort(selectedArray.astype(numpy.ushort))
        sliceFile.setFilename(fname)
        sliceFile.writeFile()
    elif(numDims==4):
        fname = options.output + os.path.sep + "fams.txt"
        dimArray = numpy.array([multispectral_train.shape[0], multispectral_train.shape[1], multispectral_train.shape[2], multispectral_train.shape[3]], dtype=numpy.uintc)
        arrayLen = multispectral_train.shape[0]*multispectral_train.shape[1]*multispectral_train.shape[2]*multispectral_train.shape[3]

        multispectral_train = numpy.reshape(multispectral_train, arrayLen, 'F')
        sliceFile = fams_ascii.FamsFile_ASCII()
        sliceFile.SetDimensions(dimArray)
        sliceFile.setDataorder(fams_ascii.COLUMN_MAJOR)
        sliceFile.setDataByUShort(multispectral_train.astype(numpy.ushort))
        sliceFile.setFilename(fname)
        sliceFile.writeFile()
    else:
        exit()