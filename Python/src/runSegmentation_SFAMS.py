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
from volumeFileInterface import fams_ascii
import math
import time
import sys
import traceback
from tools import loadFAMSmodes
import Image4D
import ValueRangeReduce_4D
import std_typedefs
from volumeFileInterface import mhd

def is_locked(filepath):
    """Checks if a file is locked by opening it in append mode.
    If no exception thrown, then the file is not locked.
    """
    locked = None
    file_object = None
    if os.path.exists(filepath):
        try:
            #print "Trying to open %s." % filepath
            buffer_size = 8
            # Opening file in append mode and read the first 8 characters.
            file_object = open(filepath, 'a', buffer_size)
            if file_object:
            #    print "%s is not locked." % filepath
                locked = False
        except IOError, message:
            #print "File is locked (unable to open in append mode). %s." % (message)
            locked = True
        finally:
            if file_object:
                file_object.close()
            #    print "%s closed." % filepath
    else:
        pass
        #print "%s not found." % filepath
    return locked



if __name__ == '__main__':
    
    optionParser = OptionParser("usage: %prog -i <input matrix file> -o <output file/folder for MHDs> [-n <field name inside mat file>]")
    optionParser.add_option("-i","--input",action="store",dest="input",default="",help="input path to matrix file")
    optionParser.add_option("-d","--dir_fams",action="store",dest="famsDir",default="",help="local dir for fams data")
    optionParser.add_option("-o","--output",action="store",dest="output",default="",help="output file (2D/3D) or folder (4D+) for MHD files")
    optionParser.add_option("-f","--fieldName",action="store",dest="fieldName",default="data",help="optional: field name to be extracted from .mat file")

    optionParser.add_option("-S","--spectral_stride",action="store",dest="spectralStride",type="int",default=-1,help="stride between spectral bins selected for evaluation")
    optionParser.add_option("-s","--train_stride",action="store",dest="trainStride",type="int",default=-1,help="(spatial) stride between samples for the training dataset")

    optionParser.add_option("-K","--famsK",action="store",dest="famsK",type="int",default=24,help="K parameter for fast adaptive mean shift")
    optionParser.add_option("-L","--famsL",action="store",dest="famsL",type="int",default=35,help="L parameter for fast adaptive mean shift")
    optionParser.add_option("-k","--famsk",action="store",dest="famsk",type="int",default=100,help="kNN parameter for max number of samples per bin")
    optionParser.add_option("-D","--mhdDimension",action="store",dest="numDims",type="int",default=4,help="number of dims in data; optional param, only needed for MHD input")
    (options,args) = optionParser.parse_args()
    
    fieldName=options.fieldName
    numDims=options.numDims
    matNumPyArray = numpy.zeros((1,1,1))
    if('.h5' in options.input):
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
        numDims = len(matNumPyArray.shape)
        arrayType = matNumPyArray.dtype
        oshape  = matNumPyArray.shape
    elif('.mhd' in options.input):
        sliceFile = mhd.MhdFile()
        sliceFile.setFilename(options.input)
        sliceFile.readFile()
        oshape = numpy.asarray(sliceFile.GetDimensions(numDims),dtype=numpy.long)
        #print oshape
        arrayLen = 1
        for i in range(0, numDims):
            arrayLen *= oshape[i]
        matNumPyArray = sliceFile.getDataAsDouble(arrayLen)
        matNumPyArray = numpy.reshape(matNumPyArray, oshape, order='F')
        #matNumPyArray = numpy.reshape(matNumPyArray, oshape, order='C')
    
    matNumPyArray = numpy.squeeze(matNumPyArray)
    numDims = len(matNumPyArray.shape)
    oshape  = matNumPyArray.shape
    #print numDims
    #**************************************************************************************************#
    #**************************************************************************************************#
    #                           P R E P R O C E S S I N G   B L O C K                                  #
    #**************************************************************************************************#
    #**************************************************************************************************#
    start_preprocessing = time.time()
    # ==== ==== NEW NAN & INF REPLACEMENT ==== ==== #
    flatIterator = matNumPyArray.flat
    for entry in flatIterator:
        try:
            if(math.isnan(entry)):
                entry=0
            if(math.isinf(entry)):
                entry=0
        except:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            print repr(traceback.format_exception(exc_type, exc_value, exc_traceback))

    # ==== ==== OLD NAN & INF REPLACEMENT ==== ==== #
    #for X_index in itertools.islice(itertools.count(),0,oshape[0]):
    #    for Y_index in itertools.islice(itertools.count(),0,oshape[1]):
    #        for Z_index in itertools.islice(itertools.count(),0,oshape[2]):
    #            for E_index in itertools.islice(itertools.count(),0,oshape[3]):
    #                try:
    #                    if(math.isnan(matNumPyArray[X_index,Y_index,Z_index,E_index])):
    #                        matNumPyArray[X_index,Y_index,Z_index,E_index]=0
    #                    if(math.isinf(matNumPyArray[X_index,Y_index,Z_index,E_index])):
    #                        matNumPyArray[X_index,Y_index,Z_index,E_index]=0
    #                except:
    #                    print "error with %i, %i, %i\n" % (X_index, Y_index, Z_index)
    #                    exc_type, exc_value, exc_traceback = sys.exc_info()
    #                    print repr(traceback.format_exception(exc_type, exc_value, exc_traceback))

    # ==== ==== COMPUTE THE (NORMALIZED) SPECTRAL GRADIENT BASED ON THE ORIGINAL VALUES ==== ==== #
    aMax, aMin = numpy.nanmax(matNumPyArray), numpy.nanmin(matNumPyArray)
    #print "Min: %f, Max: %f" % (aMin, aMax)
    dLambda_data = numpy.gradient((matNumPyArray - aMin)/(aMax - aMin)) # normalized gradient
    #dLambda_data = numpy.gradient(matNumPyArray)    # non-normalized gradient
    #print dLambda_data.shape
    dLambda_data=(dLambda_data[numDims-1]*1023)+1023        #remember gradient range [-1,1]
    dLambda_data=dLambda_data.astype(numpy.ushort)
    
    # ==== ==== OLD DATA SCALING / VALUE RANGE REDUCTION ==== ==== #
    matNumPyArray = (matNumPyArray - aMin)/(aMax - aMin)
    matNumPyArray = matNumPyArray*65535
    matNumPyArray = matNumPyArray.astype(numpy.ushort)
    ushrtArray = numpy.reshape(matNumPyArray, oshape, order='C')
    dLambda_data = numpy.reshape(dLambda_data, oshape, order='C')
    # ======== ======== NEW VALUE RANGE REDUCTION ====== ====== #
    #img_in = Image4D.Image4Ddouble()
    #img_in.setDataorder(std_typedefs.ROW_MAJOR)
    #img_in.setDatatype(std_typedefs.DOUBLE)
    #npDims = numpy.array([oshape[0], oshape[1], oshape[2], oshape[3]], dtype=numpy.uint32)
    #img_in.SetDimensions(npDims)
    #img_in.SetData(matNumPyArray.flatten())
    #
    #img_out = Image4D.Image4Dushort()
    #img_out.setDataorder(std_typedefs.ROW_MAJOR)
    #img_out.setDatatype(std_typedefs.USHORT)
    #img_out.SetDimensions(npDims)
    #img_out.createImage()
    #
    #reducer = ValueRangeReduce_4D.ValueRangeReduce_4D()
    #reducer.reduce(img_in, img_out)
    #
    #ushrtArray = img_out.GetData(oshape[0]*oshape[1]*oshape[2]*oshape[3])
    #ushrtArray = numpy.reshape(ushrtArray, oshape, order='C') # 'C' together with mhd reshape 'F' working for vis, but not for segmentation
    # ======== ======== VALUE RANGE REDUCTION END ====== ====== #
    
    spectralStride = options.spectralStride
    trainStride = options.trainStride
    
    selected_channels=[]
    #selected_channels = [0,1,2,3,4,5]
    #selected_channels = [19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75,77]
    if(spectralStride>=0):
        selected_channels=[]
    #    for binIdx in range(0,60,spectralStride):
    #        selected_channels.append(19+binIdx)
        offset = math.floor(math.floor(89/spectralStride)/2)
        for binIdx in range(0,89,spectralStride):
            if((16+offset+binIdx) < 105):
                selected_channels.append(int(16+offset+binIdx))
    else:
        if(oshape[numDims-1]<128):
            #print oshape[numDims-1]
            for i in range(0,oshape[numDims-1]):
                selected_channels.append(i)
        else:
            selected_channels = [19,27,35,43,51,59,67,75,83,91,99]
    #print selected_channels

    end_preprocessing = time.time()
    print "Preprocessing dataset in %d seconds" % (end_preprocessing-start_preprocessing)



    #***************************************************************************************************************#
    #***************************************************************************************************************#
    #                                S U B S A M P L I N G   S T E P                                                #
    #***************************************************************************************************************#
    #***************************************************************************************************************#    
    nspectra=len(selected_channels)
    nfeatures = 2*nspectra
    multispectral_train = numpy.zeros((1,1), dtype=ushrtArray.dtype)
    multispectral_data = numpy.zeros((1,1), dtype=ushrtArray.dtype)
    if numDims==4:
        multispectral_train = numpy.zeros( (oshape[0],oshape[1],oshape[2], nfeatures), dtype=ushrtArray.dtype )
        multispectral_data = numpy.zeros( (oshape[0],oshape[1],oshape[2], nfeatures), dtype=ushrtArray.dtype )
    elif numDims==3:
        multispectral_train = numpy.zeros( (oshape[0],oshape[1], nfeatures), dtype=ushrtArray.dtype )
        multispectral_data = numpy.zeros( (oshape[0],oshape[1], nfeatures), dtype=ushrtArray.dtype )
    if (trainStride>0) and (numDims==4):
        #oshape  = matNumPyArray.shape
        #nspectra=len(selected_channels)
        #nfeatures = 2*nspectra+3
        nshape = (int(math.ceil(float(oshape[0])/trainStride)), int(math.ceil(float(oshape[1])/trainStride)), int(math.ceil(float(oshape[2])/trainStride)), nfeatures)
        numSamples = int(math.ceil(float(oshape[0])/trainStride)) * int(math.ceil(float(oshape[1])/trainStride)) * int(math.ceil(float(oshape[2])/trainStride))
        multispectral_train = numpy.zeros( nshape, dtype=numpy.double)
        print "multispectral_train shape:" + str(multispectral_train.shape)
        print "multispectral_data shape:" + str(multispectral_data.shape)
        start_ordering = time.time()
        D_index=0
        #for D_index in itertools.islice(itertools.count(),0,nshape[0]*nshape[1]*nshape[2]):
        for X_index in itertools.islice(itertools.count(),0,oshape[0]):
            for Y_index in itertools.islice(itertools.count(),0,oshape[1]):
                for Z_index in itertools.islice(itertools.count(),0,oshape[2]):
                    X_train_index = int(math.floor(X_index/trainStride))
                    Y_train_index = int(math.floor(Y_index/trainStride))
                    Z_train_index = int(math.floor(Z_index/trainStride))

                    for E_index in itertools.islice(itertools.count(),0,len(selected_channels)):
                        E_data_index = selected_channels[E_index]

                        multispectral_data[X_index,Y_index,Z_index,E_index*2+0] = ushrtArray[X_index,Y_index,Z_index,E_data_index]
                        multispectral_data[X_index,Y_index,Z_index,E_index*2+1] = dLambda_data[X_index,Y_index,Z_index,E_data_index]
                        if (X_index%trainStride == 0) and (Y_index%trainStride == 0) and (Z_index%trainStride == 0):
                            multispectral_train[X_train_index,Y_train_index,Z_train_index,E_index*2+0] = ushrtArray[X_index,Y_index,Z_index,E_data_index]
                            multispectral_train[X_train_index,Y_train_index,Z_train_index,E_index*2+1] = dLambda_data[X_index,Y_index,Z_index,E_data_index]


        end_ordering = time.time()
        print "Set up the training data in %d seconds" % (end_ordering-start_ordering)
    elif(numDims==4):
        for X_index in itertools.islice(itertools.count(),0,oshape[0]):
            for Y_index in itertools.islice(itertools.count(),0,oshape[1]):
                for Z_index in itertools.islice(itertools.count(),0,oshape[2]):
                    multispectral_data[X_index,Y_index,Z_index,:] = numpy.concatenate((ushrtArray[X_index,Y_index,Z_index,selected_channels],dLambda_data[X_index,Y_index,Z_index,selected_channels]))
                    multispectral_train[X_index,Y_index,Z_index,:] = numpy.concatenate((ushrtArray[X_index,Y_index,Z_index,selected_channels],dLambda_data[X_index,Y_index,Z_index,selected_channels]))
    elif(numDims==3):
        for X_index in itertools.islice(itertools.count(),0,oshape[0]):
            for Y_index in itertools.islice(itertools.count(),0,oshape[1]):
                multispectral_data[X_index,Y_index,:] = numpy.concatenate((ushrtArray[X_index,Y_index,selected_channels],dLambda_data[X_index,Y_index,selected_channels]))
                multispectral_train[X_index,Y_index,:] = numpy.concatenate((ushrtArray[X_index,Y_index,selected_channels],dLambda_data[X_index,Y_index,selected_channels]))
    else:
        pass

    ## == NOT NOW - TO BE DONE LATER == ##
    program_file = ""
    if(numDims==3):
        fname = options.famsDir + os.path.sep + "fams.txt"
        dimArray = numpy.array([multispectral_train.shape[0], multispectral_train.shape[1], multispectral_train.shape[2]], dtype=numpy.uint32)
        arrayLen = multispectral_train.shape[0]*multispectral_train.shape[1]*multispectral_train.shape[2]
        multispectral_train = numpy.reshape(multispectral_train, arrayLen) #, 'F'
        sliceFile = fams_ascii.FamsFile_ASCII()
        sliceFile.SetDimensions(dimArray)
        #sliceFile.setDataorder(fams_ascii.COLUMN_MAJOR)    # switched after address correction
        sliceFile.setDataorder(fams_ascii.ROW_MAJOR)
        sliceFile.setDataByUShort(multispectral_train.astype(numpy.ushort))
        sliceFile.setFilename(fname)
        sliceFile.writeFile()
        program_file = "fams2D"
    elif(numDims==4):
        fname = options.famsDir + os.path.sep + "fams.txt"
        dimArray = numpy.array([multispectral_train.shape[0], multispectral_train.shape[1], multispectral_train.shape[2], multispectral_train.shape[3]], dtype=numpy.uint32)
        arrayLen = multispectral_train.shape[0]*multispectral_train.shape[1]*multispectral_train.shape[2]*multispectral_train.shape[3]
        multispectral_train = numpy.reshape(multispectral_train, arrayLen) #, 'F'
        sliceFile = fams_ascii.FamsFile_ASCII()
        sliceFile.SetDimensions(dimArray)
        #sliceFile.setDataorder(fams_ascii.COLUMN_MAJOR)    # switched after address correction
        sliceFile.setDataorder(fams_ascii.ROW_MAJOR)
        sliceFile.setDataByUShort(multispectral_train.astype(numpy.ushort))
        sliceFile.setFilename(fname)
        sliceFile.writeFile()
        program_file = "fams3D"
    else:
        exit()


    #*************************************************************************************************************#
    #*************************************************************************************************************#
    #                            E X E C U T I N G   F A M S   S E G M E N T A T I O N                            #
    #*************************************************************************************************************#
    #*************************************************************************************************************#
    start_famsproc = time.time()
    # calling convention: fams<2D/3D> <K> <L> <kNN> <filename> <working dir>
    #cmmd = r'"%s" %i %i %i %s %s' % ("fams3D", 24, 35, 100, "fams", options.famsDir)
    #cmmd = r'"%s" %i %i %i %s %s' % ("fams3D", 21, 8, 100, "fams", options.famsDir)
    cmmd = r'"%s" %i %i %i %s %s' % (program_file, options.famsK, options.famsL, options.famsk, "fams", options.famsDir)
    print "calling %s" %(cmmd)
    status = os.system(cmmd)

    if(status==0):
        famsFilePath = os.path.join(options.famsDir, "modes_fams.txt")
        while (not os.path.exists(famsFilePath)) or (is_locked(famsFilePath)):
            print "'%s' not ready for reading ..." % (famsFilePath)
            time.sleep(0.25)
        print "'%s' ready for processing ..." %(famsFilePath)
        
    else:
        if(status==1):
            print "Error: %s" % ("cpp reconstruction failed")
            exit(status)
        elif(status==2):
            print "Error: %s" % ("invalid parameters for cpp reconstruction - FAIL.")
            exit(status)
        elif(status==126):
            print "Error: %s" % ("cannot invoke cpp reconstruction - FAIL.")
            exit(status)
        elif(status==127):
            print "Error: %s" % ("executable file for cpp reconstruction could not be found - FAIL.")
            exit(status)
        elif(status==128):
            print "Error: %s" % ("cpp reconstruction with invalid exit code - FAIL.")
            exit(status)
        elif(status==130):
            print "Error: %s" % ("manual abort of cpp reconstruction - FAIL.")
            exit(status)

    f_out_name = options.famsDir + os.path.sep + "modes_fams.txt"
    reader = loadFAMSmodes.loadFAMSmodes()
    reader.load(f_out_name)
    modes = reader._modes
    dist = numpy.zeros(modes.shape[0], dtype=numpy.float32)
    #print "# segments: %i" % (modes.shape[0])
    if numDims==4:
        segments = numpy.zeros((oshape[0],oshape[1],oshape[2]), dtype=numpy.ushort) #, order='F'
    elif numDims==3:
        segments = numpy.zeros((oshape[0],oshape[1]), dtype=numpy.ushort) #, order='F'
    end_famsproc = time.time()
    print "Ran FAMS segmentation in %d seconds" % (end_famsproc-start_famsproc)
    multispectral_data = numpy.asarray(multispectral_data, dtype=numpy.float32)


    #*************************************************************************************************************#
    #*************************************************************************************************************#
    #                            A P P L Y I N G   F A M S   S E G M E N T A T I O N                              #
    #*************************************************************************************************************#
    #*************************************************************************************************************#
    final_segments = []
    start_segmentData = time.time()
    if(numDims==3):
        segments = numpy.zeros((oshape[0],oshape[1]), dtype=numpy.ushort)
        for X_index in itertools.islice(itertools.count(),0,oshape[0]):
            for Y_index in itertools.islice(itertools.count(),0,oshape[1]):
                for M_index in itertools.islice(itertools.count(), modes.shape[0]):
                    dv = modes[M_index,:] - multispectral_data[X_index,Y_index,:]
                    dist[M_index] = numpy.linalg.norm(dv)
                seg = numpy.argmin(dist)
                segments[X_index,Y_index]=seg
                if seg not in final_segments:
                    final_segments.append(seg)
                    #print "new seg at (%d,%d,%d): %d" (X_index, Y_index, Z_index, seg)
    elif(numDims==4):
        for X_index in itertools.islice(itertools.count(),0,oshape[0]):
            for Y_index in itertools.islice(itertools.count(),0,oshape[1]):
                for Z_index in itertools.islice(itertools.count(),0,oshape[2]):
                    for M_index in itertools.islice(itertools.count(),0,modes.shape[0]):
                        dv = modes[M_index,:] - multispectral_data[X_index,Y_index,Z_index,:]
                        dist[M_index] = numpy.linalg.norm(dv)
                    seg = numpy.argmin(dist)
                    segments[X_index,Y_index,Z_index]=seg
                    if seg not in final_segments:
                        final_segments.append(seg)
                        #print "new seg at (%d,%d,%d): %d" (X_index, Y_index, Z_index, seg)
    end_segmentData = time.time()
    print "Applied Segmentation to the rest of the data in %d seconds (final # segments: %d)" % (end_segmentData-start_segmentData, len(final_segments))

    #os.remove(os.path.join(options.famsDir, "fams.txt"))
    os.remove(os.path.join(options.famsDir, "modes_fams.txt"))
    os.remove(os.path.join(options.famsDir, "out_fams.txt"))
    os.remove(os.path.join(options.famsDir, "pilot_%d_fams.txt" % (options.famsk)))

    f = h5py.File(options.output,'w')
    data_group = f.create_group(options.fieldName)
    segment_data_hdf5 = segments.transpose();
    data_group.create_dataset('value', data=segment_data_hdf5);
    f.close()
