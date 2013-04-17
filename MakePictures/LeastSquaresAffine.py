'''
Created on Dec 15, 2012

@author: Tony
'''
import os
import sys
import time
import math
import numpy
#import scipy
import argparse
import subprocess
import nibabel as nib
from scipy import spatial
#from PIL import Image
#from PIL import ImageFont
#from PIL import ImageDraw
import cv2
#===============================================================================
# PLOT
#===============================================================================
from mpl_toolkits.mplot3d import Axes3D, art3d
import matplotlib
import matplotlib.pyplot as plt

#===============================================================================
# PARSE INPUT
# -I R:\nifti\CP10104_v1\Biomotion\CP10104_v3_BOLD_BIOMOTION1.nii.gz -O R:\tmp\pretty -S 6 -E 64 -N 36
#===============================================================================
parser = argparse.ArgumentParser(description="Script to convert nifit image slices into PNG...")
parser.add_argument("-d", "--ImageDir", dest="niftiDir", default=os.getcwd(), type=str, help="specify nifti directory...")
parser.add_argument("-D", "--ImageFileDir", dest="niftiFileDir", default=os.getcwd(), type=str, help="specify nifti directory...")
parser.add_argument("-N", "--ImageFile", dest="niftiFile", default=os.getcwd(), type=str, help="specify nifti directory...")
parser.add_argument("-O", "--OutputDir", dest="outputDir", type=str, help="specify where to write the output files...")
parser.add_argument("-B", "--BetDir", dest="betDir", type=str, help="specify where to write the output files...")
parser.add_argument("-MXR", "--maxRange", dest="MaxRange", type=float, default=-1., help="specify max range...")
parser.add_argument("-MNR", "--minRange", dest="MinRange", type=float, default=-1., help="specify min range...")

InputArgs = parser.parse_args()
niftiFileDir = InputArgs.niftiFileDir
niftiDir = InputArgs.niftiDir
niftiFile = InputArgs.niftiFile
outputDir = InputArgs.outputDir
betDir = InputArgs.betDir
MaxRange = InputArgs.MaxRange
MinRange = InputArgs.MinRange
##===============================================================================

##===============================================================================
sTime = time.time()
if not os.path.exists(outputDir):
    os.makedirs(outputDir)
if (niftiDir[-1] != os.sep):
    niftiDir += os.sep 
if (betDir[-1] != os.sep):
    betDir += os.sep
    
Vis = False
Write = False
#===============================================================================
# functions...
#===============================================================================
def fStripExtension( inputName ):
    inputNameNoExtension, inputNameExtension = os.path.splitext(inputName)
    
    if inputNameExtension == '.gz':
        inputNameNoExtension, inputNameExtension = os.path.splitext(inputNameNoExtension)
        return inputNameNoExtension
    else:
        return inputNameNoExtension
#===============================================================================
def fNormImage( inputImage, Gamma, MaxRange, MinRange ):
    
    inputImage = numpy.power( inputImage, Gamma )
    ImgMin = min(inputImage.ravel())
    ImgStretch = (inputImage - ImgMin)
    ImgMax = max(ImgStretch.ravel())
    ImgNorm = numpy.uint8( (ImgStretch / ImgMax) * 255 )
    return ImgNorm 
#===============================================================================

#===============================================================================
# read the nifti file...
#===============================================================================
InputFileName = os.path.basename(niftiFileDir)
(InputFileDir, InputFileName) = os.path.split(niftiFileDir)
niftiFileBaseName = fStripExtension( InputFileName )

nimBabel = nib.load(niftiFileDir)
nimBabelData = nimBabel.get_data()
nimBabelInfo = nib.get_info()
nimBabelAffine = nimBabel.get_affine()
nimBabelHeader = nimBabel.get_header().structarr
nimBabelPixDim = nimBabelHeader['pixdim']
numpyBabelDataSz = numpy.asarray(nimBabelData.shape)
print "Input size: " + str(numpyBabelDataSz) + ' Max: ' + str(numpy.max(nimBabelData.ravel())) + ' Min: ' + str(numpy.min(nimBabelData.ravel()))
#===============================================================================
# read the MNI file...
#===============================================================================
nibMNI = nib.load(niftiDir + 'MNI152_T1_0.7mm.nii.gz')
nibMNIBrain = nib.load(niftiDir + 'MNI152_T1_0.7mm_brain.nii.gz')
nibMNIBrainMask = nib.load(niftiDir + 'MNI152_T1_0.7mm_brain_mask.nii.gz')
nibMNIData = nibMNI.get_data()

nibNMIBrainMaskIdx = numpy.nonzero( nibMNIBrainMask.get_data() > 0 )
nibNMIBrainMaskSize = numpy.asarray(nibNMIBrainMaskIdx).shape
standardCentroid = numpy.sum(numpy.asarray(nibNMIBrainMaskIdx, dtype=numpy.float32), axis=1) / nibNMIBrainMaskIdx[0].shape
standardXYZ = numpy.subtract(numpy.asarray(nibNMIBrainMaskIdx, dtype=numpy.float32), standardCentroid.reshape(3,1))
print "Standard size: " + str(numpy.asarray(nibMNIData.shape)) + ' Max: ' + str(numpy.max(nibMNIData.ravel())) + ' Min: ' + str(numpy.min(nibMNIData.ravel()))
#===============================================================================
# Do BET, save mask to nifti...
#===============================================================================
img = nib.AnalyzeImage( nimBabelData, nimBabelAffine )
img.to_filename(niftiDir + niftiFileBaseName)

if not os.path.exists(niftiDir + niftiFileBaseName + '_mask.hdr'):
    commandBET = '%sbet.exe %s%s %s%s -m -s -n' % (betDir, niftiDir, niftiFileBaseName, niftiDir, niftiFileBaseName)
    with open(os.devnull, "w") as fnull:
        subprocBET = subprocess.call( commandBET, stdout = fnull, stderr = fnull )
    
analBabel = nib.load(niftiDir + niftiFileBaseName + '_mask.hdr')
analBabelData = analBabel.get_data()
analMaskDataTmp = numpy.asarray(analBabelData, dtype='int')
analMaskData = numpy.reshape(analMaskDataTmp, [analMaskDataTmp.shape[0], analMaskDataTmp.shape[1], analMaskDataTmp.shape[2]])
maskImg = nib.Nifti1Image(analMaskData, None, nimBabel.get_header())
nib.save(maskImg, niftiDir +niftiFileBaseName + '_mask.nii.gz')
nibInputBrainMask = nib.load(niftiDir +niftiFileBaseName + '_mask.nii.gz')
inputBrainMaskIdx = numpy.nonzero( analMaskData > 0 )
inputBrainMaskSize = numpy.asarray(inputBrainMaskIdx).shape
inputBrainMaskCentroid = numpy.sum(numpy.asarray(inputBrainMaskIdx, dtype=numpy.float32), axis=1) / inputBrainMaskIdx[0].shape
inputXYZ = numpy.subtract(numpy.asarray(inputBrainMaskIdx, dtype=numpy.float32), inputBrainMaskCentroid.reshape(3,1))

#print numpy.max(nibNMIBrainMaskIdx, axis = 1), numpy.max(inputBrainMaskIdx, axis = 1)
#===============================================================================
# Do convex hull and scale input to standard volume...
#===============================================================================
inputHull = spatial.ConvexHull(inputXYZ.T)
standardHull = spatial.ConvexHull(standardXYZ.T)

inputVolume = 0
for i in xrange(0, inputHull.nsimplex):
    inputVolume += (numpy.abs(numpy.linalg.det([inputHull.points[inputHull.simplices[i,0],:], inputHull.points[inputHull.simplices[i,1],:], inputHull.points[inputHull.simplices[i,2],:]])) / 6)
    
standardVolume = 0
for i in xrange(0, standardHull.nsimplex):
    standardVolume += (numpy.abs(numpy.linalg.det([standardHull.points[standardHull.simplices[i,0],:], standardHull.points[standardHull.simplices[i,1],:], standardHull.points[standardHull.simplices[i,2],:]])) / 6)

scaleFactor = (standardVolume / inputVolume) ** (1./3.)
scaleMat = [[scaleFactor, 0, 0], [0, scaleFactor, 0], [0, 0, scaleFactor]]
inputScaleXYZ = numpy.dot(scaleMat, inputXYZ)
#===============================================================================
# Check scaling with new convex hull...
#===============================================================================
inputScaleHull = spatial.ConvexHull(inputScaleXYZ.T)
inputScaleVolume = 0
for i in xrange(0, inputScaleHull.nsimplex):
    inputScaleVolume += (numpy.abs(numpy.linalg.det([inputScaleHull.points[inputScaleHull.simplices[i,0],:], inputScaleHull.points[inputScaleHull.simplices[i,1],:], inputScaleHull.points[inputScaleHull.simplices[i,2],:]])) / 6)
rescaleFactor = (standardVolume / inputScaleVolume)
#===============================================================================
# Do some subsampling to look at stuff...
#===============================================================================
subsampleStep = 2500
standardSubsampleVec = numpy.arange(0, standardHull.npoints, subsampleStep)
inputSubsampleVec = numpy.arange(0, inputHull.npoints, subsampleStep)

if len(standardSubsampleVec) > len(inputSubsampleVec):
    usableSubsampleVec = inputSubsampleVec
else:
    usableSubsampleVec = standardSubsampleVec

print inputScaleHull.npoints, inputScaleHull.nsimplex, standardHull.npoints, standardHull.nsimplex
#inputScaleSubsampleHull = spatial.ConvexHull(inputScaleHull.points[inputSubsampleVec,:])
#print inputScaleSubsampleHull.npoints, inputScaleSubsampleHull.nsimplex 

if Vis:
    fig = plt.figure()
    ax = Axes3D(fig)
    
    ax.scatter(inputScaleHull.points[inputScaleHull.simplices,0], inputScaleHull.points[inputScaleHull.simplices,1], inputScaleHull.points[inputScaleHull.simplices,2], s=20, c='r', marker='o')
    for simplex in inputScaleHull.simplices:
        simplexLooped = numpy.append(simplex, simplex[0])
        ax.plot3D(inputScaleHull.points[simplexLooped,0], inputScaleHull.points[simplexLooped,1], inputScaleHull.points[simplexLooped,2], 'm-')
#    ax.plot_surface(inputScaleHull.points[:,0], inputScaleHull.points[:,1], inputScaleHull.points[:,2])

#    fig = plt.figure()
#    ax = Axes3D(fig)
    ax.scatter(standardHull.points[standardHull.simplices,0], standardHull.points[standardHull.simplices,1], standardHull.points[standardHull.simplices,2], s=10, c='b', marker='o')
    for simplex in standardHull.simplices:
        simplexLooped = numpy.append(simplex, simplex[0])
        ax.plot3D(standardHull.points[simplexLooped,0], standardHull.points[simplexLooped,1], standardHull.points[simplexLooped,2], 'c-')
    plt.show() 
#===============================================================================
# Do some matching on voxel locations...
#===============================================================================
inputRoundXYZ = numpy.asarray(numpy.round(inputScaleXYZ, decimals=0), dtype=numpy.int16)
standardRoundXYZ = numpy.asarray(numpy.round(standardXYZ, decimals=0), dtype=numpy.int16)

nibMNIData = nibMNI.get_data()
nibMNIMaskData = nibMNIData[standardRoundXYZ[0], standardRoundXYZ[1], standardRoundXYZ[2]]
nimBabelMaskData = nimBabelData[inputRoundXYZ[0], inputRoundXYZ[1], inputRoundXYZ[2]]

standardHullIdxs = numpy.unique(standardHull.simplices.ravel())
inputScaleHullIdxs = numpy.unique(inputScaleHull.simplices.ravel())

#standardHullIdxs = standardHull.simplices.ravel()
#inputScaleHullIdxs = inputScaleHull.simplices.ravel()

standardHullData = standardHull.points[standardHullIdxs,:]
inputScaleHullData = inputScaleHull.points[inputScaleHullIdxs,:]

dInput = spatial.distance.pdist(inputScaleHullData, 'euclidean')
dStandard = spatial.distance.pdist(standardHullData, 'euclidean')
dInputStandard = spatial.distance.cdist(inputScaleHullData, standardHullData, 'euclidean')
#print numpy.sum((dInputStandard == 0)), numpy.sum((dStandard == 0)), numpy.sum((dInput == 0)),

distanceCrit = math.sqrt(2)
matchIdx = numpy.nonzero(dInputStandard < distanceCrit) 
if Vis:
    print 'foo...'
    #fig = plt.figure()
    #ax = fig.add_subplot(131)
    #im = ax.imshow(dInputStandard)
    #fig.colorbar(im)
    #
    #ax = fig.add_subplot(132)
    #im = ax.imshow(spatial.distance.squareform(dInput))
    #fig.colorbar(im)
    #
    #ax = fig.add_subplot(133)
    #im = ax.imshow(spatial.distance.squareform(dStandard))
    #fig.colorbar(im)
    #plt.show()
#===============================================================================
# Covariance matrix...[V,S,W] = svd(C), d = sign(det(W * V.T)) 
#===============================================================================
covInputStandard = numpy.asmatrix(inputScaleHullData[matchIdx[0],:], dtype=numpy.float32).T * numpy.asmatrix(standardHullData[matchIdx[1],:], dtype=numpy.float32)



U, S, V = numpy.linalg.svd(covInputStandard)
d = numpy.sign(numpy.linalg.det(V * U.T)) 
EyeMat = numpy.eye(3)
EyeMat[2,2] = d
At = V * EyeMat * U.T

transMat = [[1, 0, 0, -inputBrainMaskCentroid[0]], [0, 1, 0, -inputBrainMaskCentroid[1]], [0, 0, 1, -inputBrainMaskCentroid[2]], [0, 0, 0, 1]]

lsRot = numpy.linalg.lstsq(inputScaleHullData[matchIdx[0],:], standardHullData[matchIdx[1],:])
#lsRot = numpy.linalg.solve(inputScaleHullData[matchIdx[0],:].T, standardHullData[matchIdx[1],:])
lsMat = numpy.identity(4, dtype=numpy.float32)
lsMat[:3, :3] = lsRot[0]

eyeMat = numpy.identity(4, dtype=numpy.float32)
rotMat = numpy.linalg.inv(At)
eyeMat[:3, :3] = rotMat
rotMat = eyeMat

eyeScale = numpy.identity(4, dtype=numpy.float32)
eyeScale[:3, :3] = scaleMat
scaleMat = eyeScale

scaleTransMat = numpy.dot(scaleMat, transMat)
#M = (numpy.dot(scaleTransMat, rotMat))
M = (numpy.dot(scaleTransMat, lsMat))
invM = numpy.linalg.inv(M)

numpy.set_printoptions(precision=3, suppress=True)
print covInputStandard
print scaleTransMat
print lsRot[0]
print rotMat
print M
#print invM
#print numpy.linalg.inv(M)


#nimBabelData, 256 x 320 x 320
inputXYZIdx = (numpy.nonzero( nimBabelData > -1 ))
onesArray = numpy.ones(inputXYZIdx[0].shape)
tmpX = inputXYZIdx[0]
tmpY = inputXYZIdx[1]
tmpZ = inputXYZIdx[2]
inputXYZIdx = numpy.array([tmpX, tmpY, tmpZ, onesArray])

newXYZ = numpy.dot(inputXYZIdx.T, M)
newRoundXYZ = numpy.asarray(numpy.round(newXYZ, decimals=0), dtype=numpy.int16)


newRoundXYZ[numpy.nonzero(newRoundXYZ[:,0] > nimBabelData.shape[0]-1), 0] = nimBabelData.shape[0]-1
newRoundXYZ[numpy.nonzero(newRoundXYZ[:,1] > nimBabelData.shape[1]-1), 1] = nimBabelData.shape[1]-1
newRoundXYZ[numpy.nonzero(newRoundXYZ[:,2] > nimBabelData.shape[2]-1), 2] = nimBabelData.shape[2]-1

newRoundXYZ[numpy.nonzero(newRoundXYZ[:,0] < 0), 0] = 0
newRoundXYZ[numpy.nonzero(newRoundXYZ[:,1] < 0), 1] = 0
newRoundXYZ[numpy.nonzero(newRoundXYZ[:,2] < 0), 2] = 0

nimBabelDataVec = nimBabelData.ravel()
newInputData = numpy.zeros(nimBabelData.shape, dtype=numpy.int16)
newInputData[newRoundXYZ[:,0], newRoundXYZ[:,1], newRoundXYZ[:,2]] = nimBabelDataVec
           
#print nimBabel.get_header()
maskSignalImg = nib.Nifti1Image(newInputData, numpy.eye(4))
nib.save(maskSignalImg, niftiDir + niftiFileBaseName + '_Final.nii.gz')

    
#fig = plt.figure()          
#hist, bins = numpy.histogram(numpy.nonzero(newInputData.ravel()), bins = 100)
#plt.bar(bins, numpy.append(hist, [0]), align = 'edge', width = 0.8)
#plt.show()

if Write:
    with open(outputDir +os.sep+ 'OutputData.txt', 'wb') as outputFileObj:
        for i in xrange(0, len(nimBabelDataVec)):
            writeStr = '%d\t%d\t%d\t%d\n' % (newRoundXYZ[i,0], newRoundXYZ[i,1], newRoundXYZ[i,2], nimBabelDataVec[i])
            outputFileObj.write(writeStr)
    #        outputFileObj.writelines('\t'.join(i) + '\n' for i in newInputData)
    #        writeCode = outputFileObj.write(newInputData)
    
#print newInputData[:,:,0] 
#fig = plt.figure()
#ax = fig.add_subplot(131)
#im = ax.imshow(newInputData[:,:,0])
#fig.colorbar(im)
#
#ax = fig.add_subplot(132)
#im = ax.imshow(newInputData[:,:,200])
#fig.colorbar(im)
#
#ax = fig.add_subplot(133)
#im = ax.imshow(newInputData[:,:,280])
#fig.colorbar(im)
#plt.show()

#matchList = list()
#nLoop = int(numpy.round(inputXYZ.shape[1] * 0.0001))
#halfAlert = True
#for i in xrange(nLoop):
#    if (i > (nLoop * 0.5)) and halfAlert:
#        sys.stdout.write("\r%d of %d" % (i, nLoop))
#        sys.stdout.flush()
#        halfAlert = False
#        
#    currXYZ = inputXYZ[:,i]
#    dXYZ = numpy.sqrt(numpy.sum(numpy.square(numpy.subtract(standardXYZ, currXYZ.reshape(3,1))), axis=0))
#    
#    currVoxel = nimBabelMaskData[i]
#    dVoxel = numpy.sqrt(numpy.square(numpy.subtract(nibMNIMaskData, currVoxel)))
#    
#    dXYZVoxelSum = numpy.add(dXYZ, dVoxel)
#    matchList.append((i, dXYZVoxelSum.argmin()))
    
#    if Vis:
#        dXYZRangeBool = numpy.nonzero(dXYZ < 32.)
#            
#        plt.figure()
#        plt.subplot(211)
#        hist, bins = numpy.histogram(dXYZ, bins = 100)
#        plt.bar(bins, numpy.append(hist, [0]), align = 'edge', width = 0.8)
#        plt.subplot(212)
#        hist, bins = numpy.histogram(dVoxel, bins = 100)
#        plt.bar(bins, numpy.append(hist, [0]), align = 'edge', width = 0.8)
#        plt.show()
#        
#        plt.figure()
#        plt.plot(dXYZ[dXYZRangeBool], dVoxel[dXYZRangeBool], 'ko')
#        plt.show()
#        
#        plt.figure()
#        plt.plot(numpy.arange(0,len(dXYZVoxelSum)), dXYZVoxelSum, 'ko')
#        plt.show()


#===============================================================================
# fill new volume with rotated data....
#===============================================================================
#newSize = numpy.max(ijkT, axis=1)+1
#nimBabelVol = numpy.zeros((newSize[0], newSize[1], newSize[2]), dtype=numpy.int32)
#numpyBabelDataSz = numpy.asarray(nimBabelData.shape)
## Speed up of 6.1526x v for loop...
#nimBabelVol[ijkT[0,:], ijkT[1,:], ijkT[2,:]] = dataVec 

print("Duration: %s" % (time.time() - sTime))
            









