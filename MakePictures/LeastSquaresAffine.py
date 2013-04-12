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
parser.add_argument("-MXR", "--maxRange", dest="MaxRange", type=float, default=-1, help="specify max range...")
parser.add_argument("-MNR", "--minRange", dest="MinRange", type=float, default=-1, help="specify min range...")

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

nibNMIBrainMaskIdx = numpy.nonzero( nibMNIBrainMask.get_data() > 0 )
nibNMIBrainMaskSize = numpy.asarray(nibNMIBrainMaskIdx).shape
standardCentroid = numpy.sum(numpy.asarray(nibNMIBrainMaskIdx, dtype=numpy.float32), axis=1) / nibNMIBrainMaskIdx[0].shape
standardXYZ = numpy.subtract(numpy.asarray(nibNMIBrainMaskIdx, dtype=numpy.float32), standardCentroid.reshape(3,1))
#===============================================================================
# Do BET, save mask to nifti...
#===============================================================================
img = nib.AnalyzeImage( nimBabelData, nimBabelAffine )
img.to_filename(niftiDir + niftiFileBaseName)

if not os.path.exists(niftiDir + niftiFileBaseName + '.hdr'):
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
# Do some matching on voxel locations and inensities...
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
covInputStandard = numpy.asmatrix(inputScaleHullData[matchIdx[0],:], dtype=numpy.float32).T * numpy.asmatrix(standardHullData[matchIdx[1],:], dtype=numpy.float32)

U, S, V = numpy.linalg.svd(covInputStandard)
d = numpy.sign(numpy.linalg.det(V * U.T)) 
EyeMat = numpy.eye(3)
EyeMat[2,2] = d
At = V * EyeMat * U.T

print numpy.linalg.inv(At)
#rotXYZ = numpy.asarray(numpy.transpose(inputXYZ.T * numpy.linalg.inv(At)), dtype=numpy.float32)




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

plt.plot(numpy.min(dInputStandard, axis=1), 'r.-')
plt.plot(numpy.min(dInputStandard, axis=0), 'bo-')
plt.show()


covInputStandard = numpy.asmatrix(inputScaleHullData, dtype=numpy.float32).T * numpy.asmatrix(standardHullData, dtype=numpy.float32)

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
# Covariance matrix...[V,S,W] = svd(C), d = sign(det(W * V.T)) 
#===============================================================================
#inputXYZ = numpy.random.randn(3,32)
#standardXYZ = numpy.random.randn(3,32)
#covMat = numpy.asmatrix(inputXYZ, dtype=numpy.float32) * numpy.asmatrix(standardXYZ, dtype=numpy.float32).T
#match input INTO standard...
inputXYMatchIdx = dXYZVoxelSum[0]
standardXYZMatchIdx = dXYZVoxelSum[1]

covMat = numpy.asmatrix(inputXYZ[:,usableSubsampleVec], dtype=numpy.float32) * numpy.asmatrix(standardXYZ[:,usableSubsampleVec], dtype=numpy.float32).T
U, S, V = numpy.linalg.svd(covMat)

d = numpy.sign(numpy.linalg.det(V * U.T)) 

EyeMat = numpy.eye(3)
EyeMat[2,2] = d
At = V * EyeMat * U.T


rotXYZ = numpy.asarray(numpy.transpose(inputXYZ.T * numpy.linalg.inv(At)), dtype=numpy.float32)
#Dist = sqrt( sum( (XYZoneMatch(:) - XYZrot(:)) .^ 2) )

if Vis:
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(rotXYZ[0,usableSubsampleVec], rotXYZ[1,usableSubsampleVec], rotXYZ[2,usableSubsampleVec], s=10, c='k', marker='o')
    ax.scatter(inputXYZ[0,usableSubsampleVec], inputXYZ[1,usableSubsampleVec], inputXYZ[2,usableSubsampleVec], s=20, c='r', marker='o')
    ax.scatter(standardXYZ[0,usableSubsampleVec], standardXYZ[1,usableSubsampleVec], standardXYZ[2,usableSubsampleVec], s=20, c='b', marker='o')
    plt.show()


#===============================================================================
# rotate...could also try scipy.ndarray.rotate
#===============================================================================
#Theta = numpy.deg2rad(90)
#Rx = [[1, 0, 0, 0], [0, cos(Theta), -sin(Theta), 0], [0, sin(Theta), cos(Theta), 0], [0, 0, 0, 1]]
#Ry = [[cos(Theta), 0, sin(Theta), 0], [0, 1, 0, 0], [-sin(Theta), 0, cos(Theta), 0], [0, 0, 0, 1]]
#Rz = [[cos(Theta), -sin(Theta), 0, 0], [sin(Theta), cos(Theta), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]

R = [[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
T = [[1, 0, 0, -numpyBabelDataSz[0]/2], [0, 1, 0, -numpyBabelDataSz[1]/2], [0, 0, 1, -numpyBabelDataSz[2]/2], [0, 0, 0, 1]]
M = numpy.dot(R, T)
# need to update this with new code...
invM = numpy.linalg.inv(M)

            
# Speed up of 5.3997x versus for looping...
ndMesh = numpy.asarray(numpy.mgrid[0:numpyBabelDataSz[0], 0:numpyBabelDataSz[1], 0:numpyBabelDataSz[2]])
Is = numpy.ravel(ndMesh[0, :])
Js = numpy.ravel(ndMesh[1, :])
Ks = numpy.ravel(ndMesh[2, :])

IJKs = numpy.vstack((Is, Js))
IJKs = numpy.vstack((IJKs, Ks))
IJKs = numpy.vstack((IJKs, numpy.ones((numpyBabelDataSz[0] * numpyBabelDataSz[1] * numpyBabelDataSz[2]), dtype=numpy.int32)))

dataVec = numpy.ravel(nimBabelVol)
xyzR = numpy.dot(M, IJKs)
xyzMin = numpy.min(xyzR, 1)
invT = [[1, 0, 0, -(xyzMin[0])], [0, 1, 0, -(xyzMin[1])], [0, 0, 1, -(xyzMin[2])], [0, 0, 0, 1]]
ijkT = numpy.dot(invT, xyzR)
#dataVecT = numpy.dot(invT, numpy.dot(M, IJKs))

#===============================================================================
# fill new volume with rotated data....
#===============================================================================
newSize = numpy.max(ijkT, axis=1)+1
nimBabelVol = numpy.zeros((newSize[0], newSize[1], newSize[2]), dtype=numpy.int32)
numpyBabelDataSz = numpy.asarray(nimBabelData.shape)
# Speed up of 6.1526x v for loop...
nimBabelVol[ijkT[0,:], ijkT[1,:], ijkT[2,:]] = dataVec 

print("Duration: %s" % (time.time() - sTime))
            









