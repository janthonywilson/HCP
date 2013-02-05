'''
Created on Aug 26, 2012

@author: Tony
'''
import os
#import csv
import time
import math
import pywt
import numpy 
import scipy
import socket
import argparse
import nibabel as nib
from scipy import stats

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
parser.add_argument("-P", "--PrintData", dest="printData", type=bool, help="specify wheather to write output...")
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
    
    inputSz = numpy.asarray(inputData.shape) 
    print inputSz
    headerStr = ['Voxel', 'Scale', 'Min', 'Max', 'Mean', 'Variance', 'Skew', 'Kurtosis']
    fileName = inputFileRoot + '_' + 'fMRIWavelet.txt'

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    fileID = open(outputDir +os.sep+ fileName, 'wb')
    fileID.write(headerStr[0]+'\t')
    fileID.write(headerStr[1]+'\t')
    fileID.write(headerStr[2]+'\t')
    fileID.write(headerStr[3]+'\t')
    fileID.write(headerStr[4]+'\t')
    fileID.write(headerStr[5]+'\t')
    fileID.write(headerStr[6]+'\t')
    fileID.write(headerStr[7]+'\n')
    
    for h in xrange(0, inputSz[2]):
        for i in xrange(0, inputSz[1]):
            fileID.write('%i' % h + "\t")
            fileID.write('%i' % i + "\t")
            for j in xrange(0, inputSz[0]):
                
                if (j < inputSz[0]-1):
                    fileID.write('%.8f' % inputData[i,j,h] + "\t")
                else:
                    fileID.write('%.8f' % inputData[i,j,h] + "\n")
#===============================================================================

sTime = time.time()

waveName = 'coif1';
waveLevel = 5

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
 
#VoxelSD = scipy.std( nimBabelData, axis=3 )
#frameRMS = list()
#frameRMS_MM = list()
#frameSD = list()

nStats = 6

statsMat = numpy.ones([ waveLevel+1, nStats, numpyBabelDataSz[0] * numpyBabelDataSz[1] * numpyBabelDataSz[2]])
voxelIdx = 0

for i in xrange(0, numpyBabelDataSz[0]):
    for j in xrange(0, numpyBabelDataSz[1]):
        for k in xrange(0, numpyBabelDataSz[2]):
            
            currVec = nimBabelData[i,j,k,:]
    
            coeffs = pywt.wavedec(currVec, waveName, level=waveLevel)
            
            for m in xrange(0, waveLevel+1):
                currCoeffs = numpy.asarray(coeffs[m])
                
#                if (voxelIdx == 7492):
#                    print scipy.var(currCoeffs), scipy.std(currCoeffs)
                if (sum(currCoeffs) <= 0):
                    print currCoeffs

    
                if (sum(abs(currCoeffs)) > 0):
                    n, min_max, mean, var, skew, kurt = stats.describe(currCoeffs.ravel())
                    currMin = min_max[0]
                    currMax = min_max[1]
                else:
                    n = len(currCoeffs)
                    currMin = 0
                    currMax = 0 
                    mean = 0 
                    var = 0 
                    skew = 0 
                    kurt = 0
                    
                statsVec = [currMin, currMax, mean, var, skew, kurt]
                
                statsMat[m, :, voxelIdx] = statsVec
            
            voxelIdx += 1
            
#print statsMat




#                    kurtList.append(kurt)
#                    kurtMat[i,0,j] = kurt

#            print i, j, k, numpy.std(coeffs[0]), numpy.std(coeffs[1]), numpy.std(coeffs[2]), numpy.std(coeffs[3]), numpy.std(coeffs[4]), numpy.std(coeffs[5])
            
#    currFrame = nimBabelData[:,:,:,i]
#    currRMS = math.sqrt(  numpy.sum( (scipy.mean(currFrame.ravel()) - currFrame.ravel())**2,  ) / (numpyBabelDataSz[0] * numpyBabelDataSz[1] * numpyBabelDataSz[2])  )
#
#    frameRMS.append( currRMS )
#    frameRMS_MM.append( currRMS / (nimBabelPixDim[1] * nimBabelPixDim[2] * nimBabelPixDim[3]))
#    frameSD.append( scipy.std(currFrame.ravel()) ) 
    
#print frameRMS 
#print frameRMS_MM

#numpyData = numpy.vstack( (numpy.asarray(frameRMS_MM), numpy.asarray(frameRMS)) )
#numpyData = numpy.vstack( (numpyData, numpy.asarray(frameSD)) )

if printData:
    fPrintData( statsMat, InputfileBaseName, outputDir )







print("Temporal Wavelet Duration: %s" % (time.time() - sTime))









