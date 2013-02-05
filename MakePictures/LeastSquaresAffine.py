'''
Created on Dec 15, 2012

@author: Tony
'''
import os
import sys
import time
import numpy
import scipy
import argparse
import nibabel as nib
from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw
import cv2



#===============================================================================
# PARSE INPUT
# -I R:\nifti\CP10104_v1\Biomotion\CP10104_v3_BOLD_BIOMOTION1.nii.gz -O R:\tmp\pretty -S 6 -E 64 -N 36
#===============================================================================
parser = argparse.ArgumentParser(description="Script to convert nifit image slices into PNG...")
parser.add_argument("-I", "--ImageDir", dest="niftiFileDir", default=os.getcwd(), type=str, help="specify nifti directory...")
parser.add_argument("-O", "--OutputDir", dest="outputDir", type=str, help="specify where to write the output files...")
parser.add_argument("-D", "--SliceDim", dest="sliceDim", type=int, default=2, help="specify which dimension [0, 1, or 2] to slice...")
#parser.add_argument("-S", "--StartingSlice", dest="startSlice", type=int, help="specify which slice to start at...")
#parser.add_argument("-E", "--EndingSlice", dest="endSlice", type=int, help="specify which slice to end at...")
#parser.add_argument("-N", "--NumberSlice", dest="nSlice", type=int, default=36, help="specify how many slices to include...")
parser.add_argument("-G", "--GammaCorrect", dest="Gamma", type=float, default=0.75, help="specify how much gamma to apply to image...")
parser.add_argument("-T", "--maxRange", dest="MaxRange", type=float, default=-1, help="specify max range...")
parser.add_argument("-B", "--minRange", dest="MinRange", type=float, default=-1, help="specify min range...")

InputArgs = parser.parse_args()
niftiFileDir = InputArgs.niftiFileDir
outputDir = InputArgs.outputDir
sliceDim = InputArgs.sliceDim
#startSlice = InputArgs.startSlice
#endSlice = InputArgs.endSlice
#nSlice = InputArgs.nSlice
Gamma = InputArgs.Gamma
MaxRange = InputArgs.MaxRange
MinRange = InputArgs.MinRange
#===============================================================================

sTime = time.time()
if not os.path.exists(outputDir):
    os.makedirs(outputDir)

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
#    ImgStretch = (inputImage - ImgMin) + MinRange
    ImgStretch = (inputImage - ImgMin)
#    ImgStretch = (inputImage - MinRange)
#    ImgNorm = ImgStretch / MaxRange
    ImgMax = max(ImgStretch.ravel())
    ImgNorm = numpy.uint8( (ImgStretch / ImgMax) * 255 )
#    ImgNorm[numpy.nonzero( ImgNorm > 1 )] = 1
#    ImgNorm[numpy.nonzero( ImgNorm < 0 )] = 0
    
#    print ImgMax, MaxRange, max(ImgNorm.ravel())
    return ImgNorm 
#===============================================================================

#===============================================================================
# read the nifti file...
#===============================================================================
InputFileName = os.path.basename(niftiFileDir)
(InputFileDir, InputFileName) = os.path.split(niftiFileDir)
InputfileBaseName = fStripExtension( InputFileName )

nimBabel = nib.load(niftiFileDir)
nimBabelData = nimBabel.get_data()
nimBabelInfo = nib.get_info()
nimBabelAffine = nimBabel.get_affine()
nimBabelHeader = nimBabel.get_header().structarr
nimBabelPixDim = nimBabelHeader['pixdim']
numpyBabelDataSz = numpy.asarray(nimBabelData.shape)
print "Input size: " + str(numpyBabelDataSz) + ' Max: ' + str(numpy.max(nimBabelData.ravel())) + ' Min: ' + str(numpy.min(nimBabelData.ravel()))

#===============================================================================
# average, if functional input...
#===============================================================================
if len(numpyBabelDataSz) > 3:
    nimBabelVol = numpy.mean(nimBabelData, axis=3)
else:
    nimBabelVol = nimBabelData

if (MaxRange == -1):
    MaxRange = max(nimBabelVol.ravel())
    
if (MinRange == -1):
    MinRange = min(nimBabelVol.ravel())
    
ColorbarWidth = 24.0
ColorbarN = 16.0
nimBabelVolMax = MaxRange
nimBabelVolMin = MinRange
ColorbarVals = scipy.linspace(nimBabelVolMax, nimBabelVolMin, ColorbarN)
ColorbarValsNorm = scipy.linspace(255, 0, ColorbarN)
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
#===============================================================================
# make slices and montage...
#===============================================================================
slices = numpy.round( scipy.linspace(256, 8, num=64) )
for i in xrange(0, len(slices)):
    
    if (sliceDim == 0):
        currSlice = fNormImage( nimBabelVol[slices[i],:,:], Gamma, MaxRange, MinRange )
    if (sliceDim == 1):
        currSlice = fNormImage( nimBabelVol[:,slices[i],:], Gamma, MaxRange, MinRange )
    if (sliceDim == 2):
        currSlice = nimBabelVol[:,:,slices[i]]
        print numpy.min(currSlice.ravel()), numpy.max(currSlice.ravel())
        currSlice = fNormImage( nimBabelVol[:,:,slices[i]], Gamma, MaxRange, MinRange )
        
    print numpy.min(currSlice.ravel()), numpy.max(currSlice.ravel())
    cv2Surf = cv2.SURF()
    mask = numpy.uint8(numpy.ones(currSlice.shape))
#    keypoints, descritors = cv2Surf.detect(currSlice, numpy.uint8(numpy.ones(currSlice.shape)))
#    keypoints = cv2Surf.detect(currSlice, mask)
    keypoints = cv2Surf.detect(currSlice, None)
    
    for k in keypoints[::1000]:
        cv2.circle(currSlice,(int(k.pt[0]),int(k.pt[1])),2,(255,255,255),-1)
        cv2.circle(currSlice,(int(k.pt[0]),int(k.pt[1])),int(k.size),(255,255,255),2)
    cv2.imshow('local descriptors', currSlice)
    cv2.waitKey()
        
    currSlice = numpy.asarray(numpy.round( currSlice * 255 ), dtype=numpy.uint8)
    currImg = Image.fromarray( currSlice )
    currDraw = ImageDraw.Draw( currImg )
    currDraw.text((1, 1), str( numpy.asarray(slicesWrite[i], dtype=numpy.int) ), fill="white")
    currDraw = ImageDraw.Draw( currImg )
    
    if (i < 10):
        currImg.save(outputDir +os.sep+ InputfileBaseName +'_00'+ str(i) + '.png')
    elif (i >= 10) and (i < 100):
        currImg.save(outputDir +os.sep+ InputfileBaseName +'_0'+ str(i) + '.png')
    else:
        currImg.save(outputDir +os.sep+ InputfileBaseName +'_'+ str(i) + '.png')


    MontageImg.paste(currImg, MontageList[i])

MontageImg.save(outputDir +os.sep+ InputfileBaseName +'_Montage.png')

#===============================================================================
# make colorbar and save...
#===============================================================================
if sys.platform == 'win32':
    resourceRoot = ''
else:
    resourceRoot = '/nrgpackages/tools/hcp_qc/MakePictures/'
    
fontSize = 16
nRows = round(ColorbarImgSz[1] / ColorbarN)
rowOffset = numpy.asarray(numpy.arange(0, ColorbarImgSz[1], nRows), dtype=numpy.uint16)
ColorbarFont = ImageFont.truetype(resourceRoot + 'resources/arial.ttf', fontSize)
TextOffset = (numpy.int(nRows) / 2) - (fontSize / 2)
for i in xrange(0, len(ColorbarVals)):
    currMat = numpy.ones([nRows, ColorbarWidth], dtype=numpy.uint8) * numpy.asarray(ColorbarValsNorm[i], dtype=numpy.uint8)
    currImg = Image.fromarray( currMat )
    
    currImg = Image.new('RGBA', [numpy.int(ColorbarWidth)*2, numpy.int(nRows)])
    currDraw = ImageDraw.Draw( currImg )
    currDraw.text((1, TextOffset), str( numpy.asarray(ColorbarVals[i], dtype=numpy.int) ), fill=(253, 253, 253), font=ColorbarFont )
    currDraw = ImageDraw.Draw( currImg )
    ColorbarImg.paste(currImg, ( numpy.int(ColorbarWidth), rowOffset[i] ))

    currColor = numpy.asarray(ColorbarValsNorm[i], dtype=numpy.int)
    if (currColor == 0): currColor += 1
    currImg = Image.new('RGBA', [numpy.int(ColorbarWidth), numpy.int(nRows)], color=(currColor, currColor, currColor))
    ColorbarImg.paste(currImg, ( 0, rowOffset[i] ))
    
ColorbarImg.save(outputDir +os.sep+ InputfileBaseName +'_Colorbar.png', transparency=0)

print("Duration: %s" % (time.time() - sTime))
            









