'''
Created on Mar 7, 2017

@author: christian
'''
import dataStructures
import os
import numpy
import struct
import sys
import math
import traceback
import itertools

class loadFAMSoutput():
    _modeCenters = numpy.zeros((1,1), dtype=numpy.float32, order='F')
    _segments = numpy.zeros(1, dtype=numpy.ushort, order='F')
    _modes = numpy.zeros((1,1), dtype=numpy.float32, order='F')
    _nfeatures = 0
    
    def __init__(self):
        self._modeCenters = numpy.zeros((1,1), dtype=numpy.float32, order='F')
        self._segments = numpy.zeros(1, dtype=numpy.ushort, order='F')
        self._modes = numpy.zeros((1,1), dtype=numpy.float32, order='F')
        self._nfeatures=0

    def load(self, filename):
        if(os.path.exists(filename)==False):
            return False
        try:
            modeCenters = []
            dataString=""
            with open(filename, 'r') as text_file:
                for line in text_file:
                    dataStringArray = line.split(" ")
                    modeData = []
                    for entryString in dataStringArray:
                        entryString=entryString.strip()
                        #print entryString
                        if(len(entryString)==0):
                            continue
                        modeData.append(float(entryString))
                    modeCenters.append(modeData)
            self._modeCenters = numpy.array(modeCenters, dtype=numpy.float32, order='F')
            self._nfeatures=self._modeCenters.shape[1]
            modes = []
            for Dindex in itertools.islice(itertools.count(),0,self._modeCenters.shape[0]):
                atuple = tuple(self._modeCenters[Dindex,:])
                if atuple not in modes:
                    print atuple
                    modes.append(atuple)
            segLen = len(modes)
            self._modes = numpy.array(modes, dtype=numpy.float32, order='F')
            segments = []
            for Dindex in itertools.islice(itertools.count(),0,self._modeCenters.shape[0]):
                atuple = tuple(self._modeCenters[Dindex,:])
                indicator = modes.index(atuple)
                segments.append(indicator)
            self._segments = numpy.array(segments, dtype=numpy.ushort, order='F')
        except (TypeError, ValueError, IOError):
            exc_type, exc_value, exc_traceback = sys.exc_info()
            print repr(traceback.format_exception(exc_type, exc_value, exc_traceback))
            return False
        except:
            e = sys.exc_info()
            print e
            return False
        
        return True

