'''
Created on 2012-07-24

@author: jwilso01
'''
import os
import csv
import sys
import time
import numpy
import scipy
import argparse
import nipy as nip
import nibabel as nib
from scipy import stats
import matplotlib.pyplot as pyplot

#===============================================================================
# PARSE INPUT
#===============================================================================
parser = argparse.ArgumentParser(description="Test program to do b0 statistics on a nifit image...")
parser.add_argument("-D", "--ImageDir", dest="niftiDir", default=os.getcwd(), type=str, help="specify nifti directory...")
parser.add_argument("-S", "--Subject", dest="Subject", type=str, help="specify subject for prepending to output file...")
parser.add_argument("-N", "--DWIImage", dest="dwiFile", type=str, help="specify nifti DWI file...")
parser.add_argument("-B", "--DWIBvals", dest="bvalsFile", type=str, help="specify nifti B Vals file...")
parser.add_argument("-O", "--OutputDir", dest="outputDir", type=str, help="specify where to write the output files...")

InputArgs = parser.parse_args()
niftiDir = InputArgs.niftiDir
dwiFile = InputArgs.dwiFile
bvalsFile = InputArgs.bvalsFile
outputDir = InputArgs.outputDir
Subject = InputArgs.Subject
#===============================================================================

sTime = time.time()
printData = True
plotData = False

#===============================================================================
def fPrintData( inputMeanVec, inputStdVec, Subject, outputDir ):
    
    StatsListLen = len(inputMeanVec) 
    headerStr = ['Slice', 'b0Mean', 'b0Std'];
    fileName = Subject + '_b0Statistics.txt'

    fileID = open(outputDir +os.sep+ fileName, 'wb')
    fileID.write(headerStr[0] + "\t")
    fileID.write(headerStr[1] + "\t")
    fileID.write(headerStr[2] + "\n")

    for i in xrange(0, StatsListLen):
        fileID.write('%s' % str(i+1) + "\t")
        fileID.write('%.8f' % inputMeanVec[i] + "\t")
        fileID.write('%.8f' % inputStdVec[i] + "\n")
#===============================================================================
def fStripExtension( inputName ):
    inputNameNoExtension, inputNameExtension = os.path.splitext(inputName)
    
    if inputNameExtension == '.gz':
        inputNameNoExtension, inputNameExtension = os.path.splitext(inputNameNoExtension)
        return inputNameNoExtension
    else:
        return inputNameNoExtension
#===============================================================================

#===============================================================================
# read the nifti file...
#===============================================================================
inputFileName = dwiFile
dwiFile = niftiDir +os.sep+ dwiFile
print "Nifti File: " + dwiFile

nimBabel = nib.load(dwiFile)
numpyNibData = nimBabel.get_data()
#numpyNibData = numpy.float32(numpy.asarray(numpyNibData))
numpyNibDataSz = numpy.asarray(numpyNibData.shape)
print numpyNibDataSz

#===============================================================================
# read the bval file...
#===============================================================================
inputBvalsName = bvalsFile
bvalsFile = niftiDir +os.sep+ bvalsFile
print "Bvals File: " + bvalsFile

bvalsFileReader = csv.reader(open(bvalsFile), delimiter=" ")
bvalsList = bvalsFileReader.next()

b0List = list()
for i in xrange( 0, numpyNibDataSz[3] ):
    if (bvalsList[i] == '5'):
        b0List.append(i)

b0Volume = numpyNibData[:,:,:,b0List]

b0Mean = numpy.mean(b0Volume, axis=3)
b0Std = numpy.std(b0Volume, axis=3)

b0MeanVec = numpy.mean(b0Mean, axis=1)
b0StdVec = numpy.std(b0Std, axis=1)

b0MeanVec = numpy.mean(b0MeanVec, axis=0)
b0StdVec = numpy.std(b0StdVec, axis=0)

if printData:
    
    outputFileName = fStripExtension( inputFileName )
    fPrintData( b0MeanVec, b0StdVec, outputFileName, outputDir  )
    
print("Duration: %s" % (time.time() - sTime))
            








