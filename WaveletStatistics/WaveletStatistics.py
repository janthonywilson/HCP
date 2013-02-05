'''
Created on Jun 19, 2012

@author: Tony
'''

import os
import csv
import time
import math
import numpy 
import pywt
import socket
import argparse
import nipy as nip
import nibabel as nib
from scipy import stats

import matplotlib.cm as cm
import matplotlib.pyplot as pyplot
import matplotlib.animation as animation
#import matplotlib.image as pyimg
#import pymorph as pym
print "Running on " + socket.gethostname()

#===============================================================================
# PARSE INPUT
#===============================================================================
parser = argparse.ArgumentParser(description="Test program to do wavelet statistics on a nifit image...")
parser.add_argument("-D", "--ImageDir", dest="niftiDir", default=os.getcwd(), type=str, help="specify nifti directory...")
parser.add_argument("-N", "--Image", dest="niftiFile", type=str, help="specify nifti file...")
parser.add_argument("-O", "--OutputDir", dest="outputDir", type=str, help="specify where to write the output files...")
parser.add_argument("-F", "--FractionElim", dest="fracElimSlices", type=float, default=0.19, help="specify what fraction of top and bottom slices to eliminate...")
parser.add_argument("-U", "--UseFractionElim", dest="fracElimSlicesUse", type=str, default="y", help="specify if fraction of slices to be eliminated...")
InputArgs = parser.parse_args()
niftiDir = InputArgs.niftiDir
niftiFile = InputArgs.niftiFile
outputDir = InputArgs.outputDir
fracElimSlices = InputArgs.fracElimSlices
fracElimSlicesUse = InputArgs.fracElimSlicesUse
#===============================================================================

sTime = time.time()
printMasks = False
showAni = False
showSNR = False
printStats = True
osName = os.name
hostName = socket.gethostname()
numpy.seterr(divide = 'ignore')

waveName = 'coif1';
waveLevel = 5


#===============================================================================
# FUNCTIONS
#===============================================================================
def fPrintWaveletStats( inputMat, inputType, inputSubject, outputDir ):
    
    StatsListSz = inputMat.shape
    headerStr = ['Volume', 'cA', 'cH Scale 1', 'cV Scale 1', 'cD Scale 1', 'cH Scale 2', 'cV Scale 2', 'cD Scale 2', 'cH Scale 3', 'cV Scale 3', 'cD Scale 3', 'cH Scale 4', 'cV Scale 4', 'cD Scale 4', 'cH Scale 5', 'cV Scale 5', 'cD Scale 5']
    fileName = inputSubject + '_WaveletKurt' +inputType+ '_' +str(len(numpyNipyDataSz))+ '.txt'

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    fileStatsId = csv.writer(open(outputDir +os.sep+ fileName, 'wb'), delimiter='\t')
    fileStatsId.writerow(headerStr)
        
    if len(StatsListSz) == 2:
        fileStatsId.writerow([1, inputMat[0,0], inputMat[1,1], inputMat[2,1], inputMat[3,1], inputMat[1,2], inputMat[2,2], inputMat[3,2], inputMat[1,3], inputMat[2,3], inputMat[3,3], inputMat[1,4], inputMat[2,4], inputMat[3,4], inputMat[1,5], inputMat[2,5], inputMat[3,5]])
    else:
        for j in xrange(0, StatsListSz[2]):
                fileStatsId.writerow([j, inputMat[0,0,j], inputMat[1,1,j], inputMat[2,1,j], inputMat[3,1,j], inputMat[1,2,j], inputMat[2,2,j], inputMat[3,2,j], inputMat[1,3,j], inputMat[2,3,j], inputMat[3,3,j], inputMat[1,4,j], inputMat[2,4,j], inputMat[3,4,j], inputMat[1,5,j], inputMat[2,5,j], inputMat[3,5,j]])
#===============================================================================
def fStripSession( inputName ):
    # check for session on input subject string...
    if (inputName.find("_strc") != -1) or (inputName.find("_diff") != -1) or (inputName.find("_fnc") != -1) or (inputName.find("_xtr") != -1):
        # strip out the session stuff.  Total hack with the index < stuff...
        sessionIdx = inputName.index("_")
        inputSubject = inputName[0:sessionIdx]
        try:
            fileIdx = inputName.index(".")
            # ACK!  Hard coding...
            if (sessionIdx < 8):
                outputName = inputSubject +'.'+ inputName[fileIdx:]
        except:
            sessionIdxEnd = inputName[sessionIdx+1:].index("_")
            inputName = inputName[sessionIdxEnd+sessionIdx+2:]
            outputName = inputSubject +'_'+ inputName
        
        else:
            outputName = inputName
        
    return outputName
#===============================================================================
def fStripExtension( inputName ):
    inputNameNoExtension, inputNameExtension = os.path.splitext(inputName)
    
    if inputNameExtension == '.gz':
        inputNameNoExtension, inputNameExtension = os.path.splitext(inputNameNoExtension)
        return inputNameNoExtension
    else:
        return inputNameNoExtension
#===============================================================================
def fPrintVec( inputVec, nSlices, nFrames, inputSubject, outputDir ):
    
    StatsListLen = len(inputVec) 
    headerStr = ['Volume', 'Slice', 'contrastRMS']
    fileName = inputSubject + '_contrastRMS_' +str(len(numpyNipyDataSz))+ '.txt'

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    fileStatsId = csv.writer(open(outputDir +os.sep+ fileName, 'wb'), delimiter='\t')
    fileStatsId.writerow(headerStr)
    
    linIdx = 0
    for i in xrange(0, nFrames):
        for j in xrange(0, nSlices):
            fileStatsId.writerow([(i+1), (j+1), inputVec[linIdx]])
            linIdx += 1
#===============================================================================


    
inputFileName = niftiFile
niftiFile = niftiDir +os.sep+ niftiFile
print "File: " + niftiFile 
    

nimBabel = nib.load(niftiFile)
nimBabelAffine = nimBabel.get_affine()
nimBabelScale = [nimBabelAffine[0,0], nimBabelAffine[1,1], nimBabelAffine[2,2]]
numpyNipyData = nip.load_image(niftiFile)
numpyNipyData = numpy.float32(numpy.asarray(numpyNipyData))
numpyNipyDataSz = numpy.asarray(numpyNipyData.shape)


#set up bool for keeping middle slices...
nElimSlices = round(numpyNipyDataSz[2] * fracElimSlices)
boolElimSlices = numpy.zeros(numpyNipyDataSz[2], dtype=numpy.int)
if (fracElimSlicesUse == 'y'):
    boolElimSlices[nElimSlices:(numpyNipyDataSz[2] - nElimSlices)] = 1
else:
    boolElimSlices[0:numpyNipyDataSz[2]] = 1
ElimSlicesIdx = numpy.nonzero(boolElimSlices)

print "Eliminating " +str( round(numpyNipyDataSz[2] - sum(boolElimSlices)) )+ " of " +str(numpyNipyDataSz[2])+ " slices ..."

kurtList = list()
kurtMeanList = list()
kurtStdList = list()
contrastList = list()
if len(numpyNipyDataSz) == 4:
    kurtMatFunctionalMean = numpy.zeros([4, waveLevel+1, numpyNipyDataSz[3]])
    kurtMatFunctionalStd = numpy.zeros([4, waveLevel+1, numpyNipyDataSz[3]])
    kurtMatFunctionalMax = numpy.zeros([4, waveLevel+1, numpyNipyDataSz[3]])
    kurtMatFunctionalMin = numpy.zeros([4, waveLevel+1, numpyNipyDataSz[3]])
    for h in xrange(0, numpyNipyDataSz[3]):
        currVol = numpyNipyData[:,:,:,h]
        
        kurtMat = numpy.zeros([numpyNipyDataSz[2], 4, waveLevel+1])
        for i in xrange(0, numpyNipyDataSz[2]):
            numpyData = currVol[:,:,i]
            
            # frame by frame variation...
            currRMS = math.sqrt(  numpy.sum((numpy.mean(numpyData) - numpyData.ravel())**2 ) / (numpyNipyDataSz[0] * numpyNipyDataSz[1]) )
            contrastList.append(currRMS)

            #[phi, psi, x] = pywt.Wavelet('db2').wavefun(level=4)
            #cA, (cH, cV, cD) = pywt.dwt2(numpyData, 'db3')
            coeffs = pywt.wavedec2(numpyData, waveName, level=waveLevel)
            
        
            for j in xrange(0, waveLevel+1):
                currCoeffs = numpy.asarray(coeffs[j])
    
                #===================================================================
                # low pass coeffs
                #===================================================================
                if len(currCoeffs.shape) == 2:
                    n, min_max, mean, var, skew, kurt = stats.describe(currCoeffs.ravel())
                    kurtList.append(kurt)
                    kurtMat[i,0,j] = kurt

                #===================================================================
                # HVD coeffs 
                #===================================================================
                elif len(currCoeffs.shape) == 3:
                    for k in xrange(0, 3):
                        hvdCoeffs = currCoeffs[k,:]
                        n, min_max, mean, var, skew, kurt = stats.describe(hvdCoeffs.ravel())
                        kurtList.append(kurt)
                        kurtMat[i,k+1,j] = kurt
                        
        kurtMeanList.append(numpy.mean(kurtList))
        kurtStdList.append(numpy.std(kurtList))
        
#        kurtMatMean = numpy.mean(kurtMat, axis=0)
#        kurtMatStd = numpy.std(kurtMat, axis=0)
#        kurtMatMax = numpy.max(kurtMat, axis=0)
#        kurtMatMin = numpy.min(kurtMat, axis=0)
        
        kurtMatMean = numpy.mean(kurtMat[ElimSlicesIdx[0],:,:], axis=0)
        kurtMatStd = numpy.std(kurtMat[ElimSlicesIdx[0],:,:], axis=0)
        kurtMatMax = numpy.max(kurtMat[ElimSlicesIdx[0],:,:], axis=0)
        kurtMatMin = numpy.min(kurtMat[ElimSlicesIdx[0],:,:], axis=0)
        
        kurtMatFunctionalMean[:,:,h] = kurtMatMean
        kurtMatFunctionalStd[:,:,h] = kurtMatStd
        kurtMatFunctionalMax[:,:,h] = kurtMatMax
        kurtMatFunctionalMin[:,:,h] = kurtMatMin
        
        
    kurtMatMean = kurtMatFunctionalMean
    kurtMatStd = kurtMatFunctionalStd
    kurtMatMax = kurtMatFunctionalMax
    kurtMatMin = kurtMatFunctionalMin
    
elif len(numpyNipyDataSz) == 3:
    currVol = numpyNipyData
    
    kurtMat = numpy.zeros([numpyNipyDataSz[2], 4, waveLevel+1])
    for i in xrange(0, numpyNipyDataSz[2]):
        
        
        numpyData = currVol[:,:,i]
        
        currRMS = math.sqrt(  numpy.sum((numpy.mean(numpyData) - numpyData.ravel())**2 ) / (numpyNipyDataSz[0] * numpyNipyDataSz[1]) )
        contrastList.append(currRMS)
        #print numpyData.shape()
        #cA, (cH, cV, cD) = pywt.dwt2(numpyData, 'db3')
        coeffs = pywt.wavedec2(numpyData, waveName, level=waveLevel)
#        print [coeffs[1,0]]
        
        
        for j in xrange(0, waveLevel+1):
            currCoeffs = numpy.asarray(coeffs[j])
#            tmpCoeffs = [[[currCoeffs[0]], [currCoeffs[1,0]], [currCoeffs[1,1]], [currCoeffs[1,2]]]]

            #===================================================================
            # low pass coeffs
            #===================================================================
            if len(currCoeffs.shape) == 2:
                n, min_max, mean, var, skew, kurt = stats.describe(currCoeffs.ravel())
                kurtList.append(kurt)
                kurtMat[i,0,j] = kurt
#                pyplot.imshow(currCoeffs)
#                pyplot.show()
            #===================================================================
            # HVD coeffs 
            #===================================================================
            elif len(currCoeffs.shape) == 3:
                for k in xrange(0, 3):
                    hvdCoeffs = currCoeffs[k,:]
                    n, min_max, mean, var, skew, kurt = stats.describe(hvdCoeffs.ravel())
                    kurtList.append(kurt)
                    kurtMat[i,k+1,j] = kurt
                    
#                    pyplot.imshow(hvdCoeffs, interpolation='none', cmap=cm.gray)
#                    pyplot.show()
                    
        
        kurtMeanList.append(numpy.mean(kurtList))
        kurtStdList.append(numpy.std(kurtList))

#    kurtMatMeanNS = numpy.mean(kurtMat, axis=0) 
#    kurtMatStd = numpy.std(kurtMat, axis=0)
#    kurtMatMax = numpy.max(kurtMat, axis=0)
#    kurtMatMin = numpy.min(kurtMat, axis=0)
    
    kurtMatMean = numpy.mean(kurtMat[ElimSlicesIdx[0],:,:], axis=0)
    kurtMatStd = numpy.std(kurtMat[ElimSlicesIdx[0],:,:], axis=0)
    kurtMatMax = numpy.max(kurtMat[ElimSlicesIdx[0],:,:], axis=0)
    kurtMatMin = numpy.min(kurtMat[ElimSlicesIdx[0],:,:], axis=0)

#pyplot.imshow(kurtMat[:,:,0], interpolation='none', cmap=cm.gray)
#pyplot.show()

tTime = time.time() - sTime
print("Duration: %s" % tTime)

if printStats:
    niftiFileName = fStripExtension( inputFileName )
    outputName = fStripSession( niftiFileName )
    fPrintWaveletStats( kurtMatMean, 'Mean', outputName, outputDir )
    fPrintWaveletStats( kurtMatStd, 'STD', outputName, outputDir )
    fPrintWaveletStats( kurtMatMax, 'Max', outputName, outputDir )
    fPrintWaveletStats( kurtMatMin, 'Min', outputName, outputDir )
    
    if len(numpyNipyDataSz) == 3:
        fPrintVec( contrastList, numpyNipyDataSz[2], 1, outputName, outputDir )
    else:
        fPrintVec( contrastList, numpyNipyDataSz[2], numpyNipyDataSz[3], outputName, outputDir )

#pyplot.figure()
#if len(numpyNipyDataSz) == 4:
##    pyplot.plot(kurtMeanList, linestyle='--', marker='o', color='r', markersize=5)
#    pyplot.errorbar(xrange(0,len(kurtMeanList)), kurtMeanList, yerr=kurtStdList, linestyle='--', marker='o', color='r', markersize=5)
#else:
#    pyplot.hist(kurtList, bins=100)
##    pyplot.plot(kurtMeanList, linestyle='--', marker='o', color='r', markersize=5)
#    pyplot.figure()
#    pyplot.plot(contrastList, linestyle='--', marker='o', color='r', markersize=5)
#pyplot.show()

