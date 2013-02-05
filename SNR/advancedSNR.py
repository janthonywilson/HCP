'''
Created on 2012-06-14

@author: jwilso01
'''
import io
import pickle
import os
import csv
import math
import time
import numpy 
import scipy
import argparse
import nipy as nip
import nibabel as nib
from scipy import stats


#===============================================================================
# PARSE INPUT
#===============================================================================
parser = argparse.ArgumentParser(description="Program to do brain extraction as in PreFreeSurferPipeline...")
parser.add_argument("-D", "--ImageDir", dest="imageDir", default=os.getcwd(), type=str, help="specify nifti directory...")
parser.add_argument("-N", "--ImageName", dest="imageFileName", type=str, help="specify nifti file...")
parser.add_argument("-M", "--MaskDir", dest="maskDir", type=str, help="specify nifti mask location from BrainExtractionHCP...")
parser.add_argument("-n", "--MaskName", dest="maskFileName", type=str, help="specify nifti mask name from BrainExtractionHCP...")
parser.add_argument("-O", "--OutputDir", dest="outputDir", default=os.getcwd(), type=str, help="specify where to write the output files...")
parser.add_argument("-S", "--Subject", dest="Subject", type=str, help="specify subject for prepending to output file...")

InputArgs = parser.parse_args()
imageDir = InputArgs.imageDir
imageFileName = InputArgs.imageFileName
maskDir = InputArgs.maskDir
maskFileName = InputArgs.maskFileName
outputDir = InputArgs.outputDir
Subject = InputArgs.Subject

if not os.path.exists(outputDir):
    os.makedirs(outputDir)

sTime = time.time()
printMasks = False
saveNiftiMasks = False
printStats = True

#===============================================================================
# FUNCTION DEFINITIONS
#===============================================================================
def printVolume( Volume, FileName ):
    XYZ = numpy.asarray(Volume.shape, dtype='int64')
    
    HeaderStr = ['X', 'Y', 'Z', 'Val']
    fileId = csv.writer(open(FileName, 'wb'), delimiter='\t')
    fileId.writerow(HeaderStr)
    
    for i in range(0, XYZ[0]):
        for j in range(0, XYZ[1]):
            for k in range(0, XYZ[2]):
                fileId.writerow([i, j, k, Volume[i,j,k]])
#===============================================================================
def fPrintStats( StatsList, Subject, outputDir ):
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    
    StatsListSz = StatsList.shape
    headerStr = ['SNR', 'sigMean', 'sigMedian', 'sigStd', 'bgMean', 'bgMedian', 'bgStd', 'KSstat', 'KSpval']
    fileName = Subject + '_Statistics_3' + '.txt'
    fileStatsId = csv.writer(open(outputDir +os.sep+ fileName, 'wb'), delimiter='\t')
    fileStatsId.writerow(headerStr)
        
    for j in xrange(0, StatsListSz[0]):
        # NOTE: there is also "writerows" to get rid of the for loop...
        fileStatsId.writerow(StatsList[j,:])
#===============================================================================

structElement = scipy.ndimage.generate_binary_structure(3, 3)
#ndimage.iterate_structure(structElement, 2).astype(int)
structElement5 = scipy.ndimage.iterate_structure(structElement, 5)
structElement9 = scipy.ndimage.iterate_structure(structElement, 9)
structElement10 = scipy.ndimage.iterate_structure(structElement, 10)


#===============================================================================
# load input nifti file...
#===============================================================================
nimImage = nib.load(imageDir +os.sep+ imageFileName)
nimImageData = nimImage.get_data()
#nimImageData = numpy.asarray(nimImage.get_data(), dtype='bool')
nimImageSize = numpy.asarray(nimImage.shape)

#===============================================================================
# load brain mask...
#===============================================================================
nimMask = nib.load(maskDir +os.sep+ maskFileName)
BrainMaskData = numpy.asarray(nimMask.get_data(), dtype='int')
BrainMaskDilate = scipy.ndimage.binary_dilation(BrainMaskData, structElement10).astype(BrainMaskData.dtype)
#BrainMaskDataDilate.
print BrainMaskData.shape

#===============================================================================
# make background mask...
#===============================================================================
BackgroundMask = numpy.zeros([nimImageSize[0], nimImageSize[1], nimImageSize[2]], dtype=numpy.int8)
BackgroundIdx = numpy.nonzero( BrainMaskDilate == 0 )
BackgroundMask[BackgroundIdx] = 1
BackgroundMask.astype(bool)



SignalMaskEro = scipy.ndimage.binary_erosion(BrainMaskData, structElement5).astype(BrainMaskData.dtype)
BackgroundMaskEro = scipy.ndimage.binary_erosion(BackgroundMask, structElement5).astype(BackgroundMask.dtype)

nStats = 9
StatsList = numpy.zeros([1, nStats])

sigValsIdx = numpy.nonzero( (nimImageData * SignalMaskEro) > 0 )
bgValsIdx = numpy.nonzero( (nimImageData * BackgroundMaskEro) > 0 )

sigVals = nimImageData[sigValsIdx]
bgVals = nimImageData[bgValsIdx]

sigMean = sigVals.mean()
sigMed =  numpy.median(sigVals)
sigStd = numpy.std(sigVals)

bgMean =  bgVals.mean()
bgMed =  numpy.median(bgVals)
bgStd =  numpy.std(bgVals)

sigSS = sigVals**2
bgSS = bgVals**2

SNR = 10 * math.log10( math.sqrt(sigSS.mean()) / math.sqrt(bgSS.mean()) )
print("SNR: %s" % SNR)

#===========================================================================
# Lets do ks test...
#===========================================================================
sortSigIdx = sigVals.argsort()
sortBgIdx = bgVals.argsort()
kstestSigBg = stats.ks_2samp(sigVals[sortSigIdx], bgVals[sortBgIdx])
print kstestSigBg

StatsList[0, :] = [SNR, sigMean, sigMed, sigStd, bgMean, bgMed, bgStd, kstestSigBg[0], kstestSigBg[1]]

if printStats:
    fPrintStats( StatsList, Subject, outputDir )
    
if saveNiftiMasks:
    #===============================================================================
    # save masks to nifti...
    #===============================================================================
    maskSignalImg = nib.Nifti1Image(SignalMaskEro, None, nimMask.get_header())
    nib.save(maskSignalImg, outputDir +os.sep+ 'SignalMaskEro.nii.gz')
    
    maskBackgroundImg = nib.Nifti1Image(BackgroundMaskEro, None, nimMask.get_header())
    nib.save(maskBackgroundImg, outputDir +os.sep+ 'BackgroundMaskEro.nii.gz')
    
if printMasks:
    print "Printing volumes to text..."
    printVolume( BrainMaskDilate, outputDir +os.sep+ 'DilatedBrainMask.txt' )
    printVolume( BackgroundMaskEro, outputDir +os.sep+ 'BackgroundMaskEro.txt' )
    printVolume( SignalMaskEro, outputDir +os.sep+ 'SignalMaskEro.txt' )
    printVolume( nimImageData, outputDir +os.sep+ 'InputData.txt' )
#    printVolume( BrainMaskData, 'BrainMask.txt' )


tTime = time.time() - sTime
print("Duration: %s" % tTime)




