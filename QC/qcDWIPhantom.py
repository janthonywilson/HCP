'''
Created on Apr 30, 2013

@author: Tony
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
# USAGE: --ImageDir C:\Users\Tony\workspace\data\qcPhantom --DWIImage AGAR_EB_2013430_DWI_MB3_RL.nii.gz --DWIBvals AGAR_EB_2013430_DWI_MB3_RL.bval --DWIBvecs AGAR_EB_2013430_DWI_MB3_RL.bvec --OutputDir C:\tmp\
#===============================================================================
parser = argparse.ArgumentParser(description="Test program to do b0 statistics on a nifit image...")
parser.add_argument("-D", "--ImageDir", dest="niftiDir", default=os.getcwd(), type=str, help="specify nifti directory...")
parser.add_argument("-S", "--Subject", dest="Subject", type=str, help="specify subject for prepending to output file...")
parser.add_argument("-N", "--DWIImage", dest="dwiFile", type=str, help="specify nifti DWI file...")
parser.add_argument("-B", "--DWIBvals", dest="bvalsFile", type=str, help="specify nifti B Vals file...")
parser.add_argument("-V", "--DWIBvecs", dest="bvecsFile", type=str, help="specify nifti B Vecs file...")
parser.add_argument("-O", "--OutputDir", dest="outputDir", type=str, help="specify where to write the output files...")

InputArgs = parser.parse_args()
niftiDir = InputArgs.niftiDir
dwiFile = InputArgs.dwiFile
bvalsFile = InputArgs.bvalsFile
bvecsFile = InputArgs.bvecsFile
outputDir = InputArgs.outputDir
Subject = InputArgs.Subject
#===============================================================================

sTime = time.time()
printData = False
plotData = True

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
# read the nifti file...
#===============================================================================
if (niftiDir[-1] != os.sep):
    niftiDir = niftiDir + os.sep
        
if (outputDir[-1] != os.sep):
    outputDir = outputDir + os.sep
    
InputFileBase = dwiFile.split('.')[0]

dwiDirFile = os.path.normpath(niftiDir +os.sep+ dwiFile)
print "Nifti File: " + dwiDirFile

nibDWI = nib.load(dwiDirFile)
nibDWIData = nibDWI.get_data()
nibDWIInfo = nib.get_info()
nibDWIAffine = nibDWI.get_affine()
nibDWIHeader = nibDWI.get_header().structarr
nibDWIPixDim = nibDWIHeader['pixdim']
nibDWIDataSz = numpy.asarray(nibDWIData.shape)
print nibDWIDataSz

#===============================================================================
# read the bval file...
#===============================================================================
if bvalsFile:
    bvalsFile = os.path.normpath(niftiDir +os.sep+ bvalsFile)
else:
    bvalsFile = os.path.normpath(niftiDir +os.sep+ InputFileBase + '.bval')
print 'Bvals File: %s' % bvalsFile

with open(bvalsFile, 'rU') as bvalsFileObj:
    bvalsFileReader = csv.reader(bvalsFileObj, delimiter=' ')
    for row in bvalsFileReader:
        # this is a hack to get rid of trailing space on input file...
        bvalsList = ' '.join(row).strip().split()
            

b0List = list()
for i in xrange( 0, len(bvalsList) ):
    if (int(bvalsList[i]) <= 5):
        b0List.append(i)
#===============================================================================
# read the becs file...
#===============================================================================
bvecsList = list()
if bvecsFile:
    bvecsFile = os.path.normpath(niftiDir +os.sep+ bvecsFile)
else:
    bvecsFile = os.path.normpath(niftiDir +os.sep+ InputFileBase + '.bvec')
print "Bvecs File: " + bvecsFile

with open(bvecsFile, 'rU') as bvecsFileObj:
    bvecsFileReader = csv.reader(open(bvecsFile), delimiter=' ')
    for row in bvecsFileReader:
        # this is a hack to get rid of trailing space on input file...
        bvecsList.append(' '.join(row).strip().split())
#===============================================================================
# do B0 stuff...
#===============================================================================
b0Volume = nibDWIData[:,:,:,b0List]
if (len(b0List) == 1):
    b0Volume = numpy.reshape(b0Volume, [nibDWIDataSz[0], nibDWIDataSz[1], nibDWIDataSz[2]])
    

b0Mean = numpy.mean(b0Volume, axis=3)
b0Std = numpy.std(b0Volume, axis=3)

b0MeanVec = numpy.mean(b0Mean, axis=1)
b0StdVec = numpy.std(b0Std, axis=1)

b0MeanVec = numpy.mean(b0MeanVec, axis=0)
b0StdVec = numpy.std(b0StdVec, axis=0)

if printData:
    fPrintData( b0MeanVec, b0StdVec, outputFileName, outputDir  )
    
print("Duration: %s" % (time.time() - sTime))
            

