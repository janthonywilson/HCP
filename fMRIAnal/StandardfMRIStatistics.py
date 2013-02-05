'''
Created on Aug 20, 2012

@author: Tony
'''
import os
import time
import math
import numpy 
import scipy
import socket
import argparse
#import nipy as nip
import nibabel as nib
#from scipy import stats

#import matplotlib.cm as cm
#import matplotlib.pyplot as pyplot
#import matplotlib.animation as animation
#import matplotlib.image as pyimg
#import pymorph as pym
print "Running on " + socket.gethostname()

#===============================================================================
# PARSE INPUT
#===============================================================================
parser = argparse.ArgumentParser(description="Test program to do SNR on a nifit image...")
parser.add_argument("-I", "--Input", dest="niftiFileDir", type=str, help="specify nifti input, full path and file...")
parser.add_argument("-O", "--OutputDir", dest="outputDir", type=str, help="specify location to write...")
parser.add_argument("-P", "--PrintData", dest="printData", type=bool, help="specify whether to write output...")
InputArgs = parser.parse_args()
niftiFileDir = InputArgs.niftiFileDir
outputDir = InputArgs.outputDir
printData = InputArgs.printData
#===============================================================================
def fStripExtension( inputName ):
    inputNameNoExtension, inputNameExtension = os.path.splitext(inputName)
    
    if inputNameExtension == '.gz':
        inputNameNoExtension, inputNameExtension = os.path.splitext(inputNameNoExtension)
        return inputNameNoExtension
    else:
        return inputNameNoExtension
#===============================================================================
def fPrintData( inputData, inputFileRoot, outputDir ):
    
    inputListSz = numpy.asarray(inputData.shape) 
    headerStr = ['Volume', 'RMS/mm', 'RMS', 'SD']
    fileName = inputFileRoot + '_' + 'fMRIStatistics.txt'

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    fileID = open(outputDir +os.sep+ fileName, 'wb')
    fileID.write(headerStr[0]+'\t')
    fileID.write(headerStr[1]+'\t')
    fileID.write(headerStr[2]+'\t')
    fileID.write(headerStr[3]+'\n')
    
    for i in xrange(0, inputListSz[1]):
            fileID.write('%i' % i + "\t")
            fileID.write('%.8f' % inputData[0,i] + "\t")
            fileID.write('%.8f' % inputData[1,i] + "\t")
            fileID.write('%.8f' % inputData[2,i] + "\n")
#===============================================================================

sTime = time.time()

InputFileName = os.path.basename(niftiFileDir)
(InputFileDir, InputFileName) = os.path.split(niftiFileDir)
InputfileBaseName = fStripExtension( InputFileName )

#print InputFileName, InputfileBaseName, InputFileDir

nimBabel = nib.load(niftiFileDir)
nimBabelData = nimBabel.get_data()
nimBabelInfo = nib.get_info()
nimBabelAffine = nimBabel.get_affine()
nimBabelHeader = nimBabel.get_header().structarr
nimBabelPixDim = nimBabelHeader['pixdim']
numpyBabelDataSz = numpy.asarray(nimBabelData.shape)

if len(numpyBabelDataSz) <= 3:
    print 'ERROR: Wrong dim size input...'
    exit -1


# Temporal SNR in brain mask    
# Simple SNR (mean/SNR)        
# of frames greater than .5 from motion outliers in FLS    
# of frames with DVARs above 0.5    
# of frames with DVARs above 0.7    
# of frames with rms/mm frame above .2    
# of frames with rms/mm frame above .3    
# rms mm    
# rms mm/frame    
# voxelwise SD
 
#VoxelSD = scipy.std( nimBabelData, axis=3 )
frameRMS = list()
frameRMS_MM = list()
frameSD = list()

frameRMSVol = numpy.zeros([numpyBabelDataSz[0], numpyBabelDataSz[1], numpyBabelDataSz[2]])

for i in xrange(0, numpyBabelDataSz[0]):
    for j in xrange(0, numpyBabelDataSz[1]):
        for k in xrange(0, numpyBabelDataSz[2]):
            currVec = nimBabelData[i,j,k,:]
            currRMS = math.sqrt(  numpy.sum( (scipy.mean(currVec) - currVec)**2,  ) / (numpyBabelDataSz[0] * numpyBabelDataSz[1] * numpyBabelDataSz[2])  )
        
            frameRMSVol[i, j, k] = currRMS
#            frameRMS.append( currRMS )
#            frameRMS_MM.append( currRMS / (nimBabelPixDim[1] * nimBabelPixDim[2] * nimBabelPixDim[3]))
#            frameSD.append( scipy.std(currFrame.ravel()) ) 
            
nimCVImage = nib.Nifti1Image(frameRMSVol, nimBabelAffine)
nimCVImage.to_filename(outputDir +os.sep+ InputfileBaseName + '_RMS.nii.gz')
#print frameRMS 
#print frameRMS_MM
print frameRMSVol

#numpyData = numpy.vstack( (numpy.asarray(frameRMS_MM), numpy.asarray(frameRMS)) )
#numpyData = numpy.vstack( (numpyData, numpy.asarray(frameSD)) )

#if printData:
#    fPrintData( numpyData, InputfileBaseName, outputDir )







print("Duration: %s" % (time.time() - sTime))








