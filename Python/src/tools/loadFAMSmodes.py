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

class loadFAMSmodes():
    _modes = numpy.zeros((1,1), dtype=numpy.float32, order='C')
    _numEntries = numpy.zeros(1, dtype=numpy.uintc, order='C')
    
    def __init__(self):
        self._modes = numpy.zeros(1, dtype=numpy.float32, order='C')
        self._numEntries = numpy.zeros(1, dtype=numpy.uintc, order='C')

    def load(self, filename):
        if(os.path.exists(filename)==False):
            return False
        try:
            modes = []
            numEntries = []
            dataString=""
            with open(filename, 'r') as text_file:
                for line in text_file:
                    dataStringArray = line.split(" ")
                    numEntries.append(int(dataStringArray[0]))
                    modeData = []
                    for entryString in dataStringArray[1:len(dataStringArray)]:
                        entryString=entryString.strip()
                        #print entryString
                        if(len(entryString)==0):
                            continue
                        modeData.append(float(entryString))
                    modes.append(modeData)
                self._numEntries = numpy.array(numEntries, dtype=numpy.uintc, order='C')
                self._modes = numpy.array(modes, dtype=numpy.float32, order='C')
        except (TypeError, ValueError, IOError):
            exc_type, exc_value, exc_traceback = sys.exc_info()
            print repr(traceback.format_exception(exc_type, exc_value, exc_traceback))
            return False
        except:
            e = sys.exc_info()
            print e
            return False
        
        return True