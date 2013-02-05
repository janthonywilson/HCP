'''
Created on 2012-05-15

@author: jwilso01
'''

import os
import csv
import sys
import math
import time
import numpy 
import scipy
from scipy import stats
import nipy as nip
import nibabel as nib
import matplotlib.pyplot as pyplot
import matplotlib.animation as animation
import argparse

#===============================================================================
# PARSE INPUT
#===============================================================================
parser = argparse.ArgumentParser(description="Test program to do SNR on a nifit image...")
parser.add_argument("-D", "--ImageDir", dest="niftiDir", default=os.getcwd(), type=str, help="specify nifti directory...")
parser.add_argument("-N", "--Image", dest="niftiFile", type=str, help="specify nifti file...")
parser.add_argument("-O", "--OutputDir", dest="outputDir", type=str, help="specify location to write...")
InputArgs = parser.parse_args()
niftiDir = InputArgs.niftiDir
niftiFile = InputArgs.niftiFile
outputDir = InputArgs.outputDir

sTime = time.time()
printStats = True
printMethod = 'Char'
printMasks = False
showAni = False
showSNR = False

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

#subjectDir = 'CP10104_v1'
#niftiType = 'T1w'
inputFileName = niftiFile
niftiFile = niftiDir + os.sep + niftiFile

if not os.path.exists(outputDir):
        os.makedirs(outputDir)
        
inputFileName = fStripExtension( inputFileName )
outputName = fStripSession( inputFileName )
#print niftiFile

#===============================================================================
# totally redundant, but keeping here as a reminder
#===============================================================================
nimBabel = nib.load(niftiFile)
#nimBabelData = nimBabel.get_data()
#nimBabelData = nimBabel.get_data()[:,:,:,4]
nimBabelInfo = nib.get_info()
nimBabelAffine = nimBabel.get_affine()
nimBabelHeader = nimBabel.get_header().structarr
nimBabelPixDim = nimBabelHeader['pixdim']

nipNipyData = nip.load_image(niftiFile)
nipNipyDataSz = nipNipyData.shape
print nipNipyData.shape

numpyNipyData = numpy.float32(numpy.asarray(nipNipyData))
numpyNipyDataSz = numpy.asarray(numpyNipyData.shape)
if len(numpyNipyDataSz) > 3:
    numpyMeanImage = scipy.mean(numpyNipyData, axis=3)
    numpyStdImage = scipy.std(numpyNipyData, axis=3)
    numpyStdImageIdx = numpy.nonzero( numpyStdImage == 0 )
    
    # fix the div by zero error...
    if (len(numpyStdImageIdx[0]) > 0):    
        noiseData = numpy.random.rand(nipNipyDataSz[0], nipNipyDataSz[1], nipNipyDataSz[2]) * 1.0e-16
        numpyStdImage[numpyStdImageIdx] = noiseData[numpyStdImageIdx]
        
    CVImage = numpyMeanImage / numpyStdImage
    nimCVImage = nib.Nifti1Image(CVImage, nimBabelAffine)
    nimCVImage.to_filename(outputDir +os.sep+ outputName + '_CV.nii.gz')
    
    #===========================================================================
    # create mean image for bet...
    #===========================================================================
    nimMeanImage = nib.Nifti1Image(numpyMeanImage, nimBabelAffine)
    nimMeanImage.to_filename(outputDir +os.sep+ outputName + '_Mean.nii.gz')


    




#===============================================================================
# Make noise and add/mult to sig
#===============================================================================
#noiseSigma = 16
#noiseData = numpy.random.randn(nipNipyDataSz[0], nipNipyDataSz[1], nipNipyDataSz[2]) * noiseSigma
#numpyNipyData = numpyNipyData + noiseData


#===============================================================================
# make the mask with sphere at center and points on corners...
#===============================================================================
CenterPoint = numpyNipyDataSz / 2
CornerPoint = numpy.array([[0, 0, 0], [0, 0, numpyNipyDataSz[2]], [0, numpyNipyDataSz[1], numpyNipyDataSz[2]], [numpyNipyDataSz[0], numpyNipyDataSz[1], numpyNipyDataSz[2]], [numpyNipyDataSz[0], 0, 0], [numpyNipyDataSz[0], numpyNipyDataSz[1], 0], [numpyNipyDataSz[0], 0, numpyNipyDataSz[2]], [0, numpyNipyDataSz[1], 0]])
CornerPointSz = CornerPoint.shape
Radius = 32

#===============================================================================
# setup edge mask...
#===============================================================================
#numpyNipyEdgeMask = numpy.zeros([numpyNipyDataSz[0],numpyNipyDataSz[1],numpyNipyDataSz[2]])
#numpyNipyEdgeMask[:,:,0] = 1
#numpyNipyEdgeMask[:,0,:] = 1
#numpyNipyEdgeMask[0,:,:] = 1
#numpyNipyEdgeMask[-1,:,:] = 1
#numpyNipyEdgeMask[:,-1,:] = 1
#numpyNipyEdgeMask[:,:,-1] = 1
#
#EdgeIdx = numpy.nonzero( numpyNipyEdgeMask > 0 )
#print numpy.sum(numpyNipyEdgeMask)
    
#===============================================================================
# for i in xrange(0, numpyNipyDataSz[0]):
#    for j in xrange(0, numpyNipyDataSz[1]):
#        for k in xrange(0, numpyNipyDataSz[2]):
#       
#       
#           #===================================================================
#           # CENTER
#           #===================================================================
#           currPoint = ([[i, j, k]] - CenterPoint)**2
#           currDist = math.sqrt( currPoint.sum() )
#           if currDist <= Radius:
#               numpyNipyCenterMask[i, j, k] = 1
#               
#           #===================================================================
#           # CORNERS (Just one)
#           #===================================================================
#           cornerIdx = 0
#           currCorner = ([[i, j, k]] - CornerPoint[cornerIdx])**2
#           currDist = math.sqrt( currCorner.sum() )
#           if currDist <= Radius:
#               numpyNipyCornerMask[i, j, k] = 1
#===============================================================================

ndMesh =  numpy.asarray(numpy.mgrid[0:numpyNipyDataSz[0],0:numpyNipyDataSz[1],0:numpyNipyDataSz[2]])

numpyNipyCenterMask = numpy.zeros([numpyNipyDataSz[0],numpyNipyDataSz[1],numpyNipyDataSz[2]])
numpyNipyCornerMask = numpy.zeros([numpyNipyDataSz[0],numpyNipyDataSz[1],numpyNipyDataSz[2]])

#===============================================================================
# do the center sphere
#===============================================================================
xGrid = (ndMesh[0,:] - CenterPoint[0])**2
yGrid = (ndMesh[1,:] - CenterPoint[1])**2
zGrid = (ndMesh[2,:] - CenterPoint[2])**2  

xyzSum = (xGrid + yGrid + zGrid)
Dist = numpy.sqrt( xyzSum )

CenterIdx = numpy.nonzero( Dist <= Radius )
numpyNipyCenterMask[CenterIdx] = 1
#print numpy.sum(numpyNipyCenterMask)
#===============================================================================
# do corner spheres....
#===============================================================================
for i in range(CornerPointSz[0]):
#    currCorner = CornerPoint[i]
    
    xGrid = (ndMesh[0,:] - CornerPoint[i,0])**2
    yGrid = (ndMesh[1,:] - CornerPoint[i,1])**2
    zGrid = (ndMesh[2,:] - CornerPoint[i,2])**2
    
    xyzSum = (xGrid + yGrid + zGrid)
    Dist = numpy.sqrt( xyzSum )
    CornerIdx = numpy.nonzero( Dist <= Radius )
    numpyNipyCornerMask[CornerIdx] = 1

CornerIdx = numpy.nonzero( numpyNipyCornerMask > 0 )


#===============================================================================
# collect data...
#===============================================================================
nStats = 9
if len(numpyNipyDataSz) == 3:
    StatsList = numpy.zeros([1, nStats])
    
    sigVals = numpyNipyData[CenterIdx]
    bgVals = numpyNipyData[CornerIdx]
#    edgeVals = numpyNipyData[EdgeIdx]
    
    sigMean = sigVals.mean()
    sigMed =  numpy.median(sigVals)
    sigStd = numpy.std(sigVals)
    
    bgMean =  bgVals.mean()
    bgMed =  numpy.median(bgVals)
    bgStd =  numpy.std(bgVals)
    
#    edgeMean =  edgeVals.mean()
#    edgeMed =  numpy.median(edgeVals)    
#    edgeStd = numpy.std(sigVals)
    
    print sigMean, sigMed, sigStd, bgMean, bgMed, bgStd
    
    #===========================================================================
    # Lets do KS...
    #===========================================================================
    sortSigIdx = sigVals.argsort()
    sortBgIdx = bgVals.argsort()
    kstestSigBg = stats.ks_2samp(numpy.float64(sigVals[sortSigIdx]), numpy.float64(bgVals[sortBgIdx]))

    #MSE...
    #sigSS = (sigVals - sigMean)**2
    #bgSS = (bgVals - bgMean)**2
    #RMS...
    sigSS = sigVals**2
    bgSS = bgVals**2
#    edgeSS = edgeVals**2
    
    SNRsph = 10 * math.log10( math.sqrt(sigSS.mean()) / math.sqrt(bgSS.mean()) )
    print("SNR: %s" % SNRsph)
    
#    SNRedge = 10 * math.log10( math.sqrt(sigSS.mean()) / math.sqrt(edgeSS.mean()) ) 
    
    StatsList[0, :] = [SNRsph, sigMean, sigMed, sigStd, bgMean, bgMed, bgStd, kstestSigBg[0], kstestSigBg[1]]
    
elif len(numpyNipyDataSz) == 4:
    
    #===========================================================================
    # read brain mask...
    #===========================================================================
    BET = '/nrgpackages/tools/fsl-nrg/bin/bet2 ' + outputDir +os.sep+ outputName + '_Mean.nii.gz' +' '+  outputDir +os.sep+ outputName + ' -f 0.3 -n -m'

    if sys.platform == 'win32':
        print BET
    else:
        print BET
        os.system(BET)
        
    nipNipyMask = nip.load_image(outputDir +os.sep+ outputName + '_mask.nii.gz')
    nipNipyDataMask = nipNipyMask.get_data()
    MaskIdx = numpy.nonzero( nipNipyDataMask == 1 )
    CVMask = numpyMeanImage[MaskIdx] / numpyStdImage[MaskIdx]

    
    StatsList = numpy.zeros([numpyNipyDataSz[3], nStats])
    SNRsph = numpy.zeros(numpyNipyDataSz[3])
    SNRedge = numpy.zeros(numpyNipyDataSz[3])
    sigVoxelMat = numpy.zeros( [ numpyNipyDataSz[3], len(CenterIdx[0]) ] )
    
    for i in xrange(0, numpyNipyDataSz[3]):
        volVals = numpyNipyData[:,:,:,i]
        
        sigVals = volVals[MaskIdx]
        sigValsSph = volVals[CenterIdx]
        bgVals = volVals[CornerIdx]
#        edgeVals = volVals[EdgeIdx]
#        if edgeVals.mean() < 1e-16:
#            print 'ERROR:', sigVals.mean(), bgVals.mean(), edgeVals.mean()
        
        sigVoxelMat[i,:] = sigVals
        sigMean = sigVals.mean()
        sigMed =  numpy.median(sigVals)
        sigStd = numpy.std(sigVals)
        
        bgMean =  bgVals.mean()
        bgMed =  numpy.median(bgVals)
        bgStd =  numpy.std(bgVals)
        
#        edgeMean =  edgeVals.mean()
#        edgeMed =  numpy.median(edgeVals)    
#        edgeStd = numpy.std(sigVals)
        
        # KS Test...
        sortSigIdx = sigVals.argsort()
        sortBgIdx = bgVals.argsort()
        kstestSigBg = scipy.stats.ks_2samp(numpy.float64(sigVals[sortSigIdx]), numpy.float64(bgVals[sortBgIdx]))

        
        sigSS = sigVals**2
        bgSS = bgVals**2
#        edgeSS = edgeVals**2
        
        SNRsph[i] = 10 * math.log10( math.sqrt(sigSS.mean()) / math.sqrt(bgSS.mean()) )

        StatsList[i, :] = [SNRsph[i], sigMean, sigMed, sigStd, bgMean, bgMed, bgStd, kstestSigBg[0], kstestSigBg[1]]
        
    sigVoxelMatMean = scipy.mean(sigVoxelMat, axis=1)
    sigVoxelMatStd = scipy.std(sigVoxelMat, axis=1)
    sigVoxelMatRatio = sigVoxelMatMean / sigVoxelMatStd
    tSNR = scipy.mean(sigVoxelMat.ravel()) / scipy.std(sigVoxelMat.ravel())
    tSNRVec = numpy.zeros([1, nStats])
    tSNRVec[0,1] = tSNR
    tSNRVec[0,0] = numpy.mean(CVMask)

    # complete hack here, overwriting SNRs with CVs...
    StatsList[:, 0] = sigVoxelMatRatio
    StatsList = numpy.vstack([StatsList, tSNRVec])

        
    if showSNR:
        fig = pyplot.figure()
        pyplot.plot(sigVals[sortSigIdx], color='red', hold='on')
        pyplot.plot(bgVals[sortBgIdx], color='blue')
        pyplot.show()

if printStats:
    
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    
    StatsListSz = StatsList.shape
    headerStr = ['Volume', 'SNR', 'sigMean', 'sigMedian', 'sigStd', 'bgMean', 'bgMedian', 'bgStd', 'KSstat', 'KSpval']
    fileName = outputName + '_simpleSNR_' +str(len(numpyNipyDataSz))+ '.txt'

    if (printMethod == 'Row'):
        fileStatsId = csv.writer(open(outputDir +os.sep+ fileName, 'wb'), delimiter='\t')
        fileStatsId.writerow(headerStr)
            
        for j in xrange(0, StatsListSz[0]):
            currVol = j
            currStatsList = StatsList[j,:]
            print numpy.hstack([currVol, currStatsList])
            currVol = numpy.hstack([currVol, currStatsList])
            fileStatsId.writerow( currVol )
            
            
    elif (printMethod == 'Char'):
        fileID = open(outputDir +os.sep+ fileName, 'wb')
        for j in range(0, len(headerStr)):
            if j == len(headerStr)-1:
                fileID.write(headerStr[j] + "\n")
            else:
                fileID.write(headerStr[j] + "\t")
    
        for i in xrange(0, StatsListSz[0]):
            fileID.write('%i' % (i+1) + "\t")
            fileID.write('%.8f' % StatsList[i,0] + "\t")
            fileID.write('%.8f' % StatsList[i,1] + "\t")
            fileID.write('%.8f' % StatsList[i,2] + "\t")
            fileID.write('%.8f' % StatsList[i,3] + "\t")
            fileID.write('%.8f' % StatsList[i,4] + "\t")
            fileID.write('%.8f' % StatsList[i,5] + "\t")
            fileID.write('%.8f' % StatsList[i,6] + "\t")
            fileID.write('%.8f' % StatsList[i,7] + "\t")
            fileID.write('%.8e' % StatsList[i,8] + "\n")
        
        

tTime = time.time() - sTime
print("Duration: %s" % tTime)



if showAni:
    numpyNipyMask = numpyNipyCornerMask + numpyNipyCenterMask
    fig = pyplot.figure()
    ims = []
#    print nipNipyDataSz
    if len(nipNipyDataSz) > 3:
        for j in range(0, nipNipyDataSz[3]):
            for k in range(0, nipNipyDataSz[2]):
        #    for j in range(0, 1):
                currImg = pyplot.imshow(nipNipyData[:,:,k,j], animated=True)
                ims.append([currImg])
    else:
        for j in range(0, nipNipyDataSz[2]):
        #    for j in range(0, 1):
                currImg = pyplot.imshow(nipNipyData[:,:,j], animated=True)
                ims.append([currImg])
        
    ani = animation.ArtistAnimation(fig, ims, interval=5, blit=True, repeat_delay=2500)
    pyplot.show()

if printMasks:
    headerStr = ['X', 'Y', 'Z', 'Val']
    fileCenterId = csv.writer(open('CenterMask.txt', 'wb'), delimiter='\t')
    fileCenterId.writerow(headerStr)
    
    fileCornerId = csv.writer(open('CornerMask.txt', 'wb'), delimiter='\t')
    fileCornerId.writerow(headerStr)
    
#    fileEdgeId = csv.writer(open('EdgeMask.txt', 'wb'), delimiter='\t')
#    fileEdgeId.writerow(headerStr)
        
    for i in xrange(0, nipNipyDataSz[0]):
        for j in xrange(0, nipNipyDataSz[1]):
            for k in range(0, nipNipyDataSz[2]):
                fileCenterId.writerow([i, j, k, numpyNipyCenterMask[i, j, k]])
                fileCornerId.writerow([i, j, k, numpyNipyCornerMask[i, j, k]])
#                fileEdgeId.writerow([i, j, k, numpyNipyEdgeMask[i, j, k]])

    
    
    
    

#===============================================================================
# fig = pyplot.figure()
# imsMask = []
# for i in range(0, nipNipyDataSz[2]):
#    currImg = pyplot.imshow(numpyNipyEdgeMask[:,:,i], animated=True)
#    pyplot.colorbar()
#    imsMask.append([currImg])
#    
# aniMaks = animation.ArtistAnimation(fig, imsMask, interval=5, blit=True, repeat_delay=2500)
# pyplot.show()
#===============================================================================

#===============================================================================
# pyplot.figure()
# pyplot.subplot(211)
# pyplot.imshow(numpyNipyData[:,:,64])
# 
# pyplot.subplot(212)
# pyplot.imshow(numpyNipyDataNoisy[:,:,64])
# pyplot.show()
#===============================================================================

#===============================================================================
# pyplot.figure()
# pyplot.subplot(211)
# pyplot.hist(numpyNipyData.reshape(-1), 100)
# pyplot.subplot(212)
# pyplot.hist(numpyNipyDataNoisy.reshape(-1), 100)
# pyplot.show()
#===============================================================================



#=======================================================================
# sigMean = sigVals.mean()
# bgMean =  bgVals.mean()
# sigStd = numpy.std(sigVals)
# bgStd =  numpy.std(bgVals)
# sigMin = numpy.min(sigVals)
# bgMin =  numpy.min(bgVals)
# sigMax = numpy.max(sigVals)
# bgMax =  numpy.max(bgVals)
# print sigMean, sigStd, sigMin, sigMax, bgMean, bgStd, bgMin, bgMax
#=======================================================================





