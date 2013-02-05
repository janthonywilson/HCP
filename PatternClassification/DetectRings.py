'''
Created on 2012-07-30

@author: jwilso01
'''
import os
import cv2
import cv 
import sys
import math
import time
import Image
import numpy 
import scipy
import zipfile
import nipy as nip
import nibabel as nib
from scipy import stats
import matplotlib.pyplot as pyplot

#===============================================================================
def fPrintData( inputData, Subject, outputDir ):
    
    inputDataSz = numpy.asarray( inputData.shape )
    headerStr = ['X', 'Y', 'Rho']
    fileName = Subject + '_DetectRingsOutput.txt'

    fileID = open(outputDir +os.sep+ fileName, 'wb')
    fileID.write(headerStr[0] + "\t")
    fileID.write(headerStr[1] + "\t")
    fileID.write(headerStr[2] + "\n")

    for i in xrange(0, inputDataSz[1]):
        fileID.write('%.8f' % inputData[0, i, 0] + "\t")
        fileID.write('%.8f' % inputData[0, i, 1] + "\t")
        fileID.write('%.8f' % inputData[0, i, 2] + "\n")
#===============================================================================
def fStripExtension( inputName ):
    inputNameNoExtension, inputNameExtension = os.path.splitext(inputName)
    
    if inputNameExtension == '.gz':
        inputNameNoExtension, inputNameExtension = os.path.splitext(inputNameNoExtension)
        return inputNameNoExtension
    else:
        return inputNameNoExtension
#===============================================================================
def fCircle( x0, y0, rho ):
    angVec = numpy.arange(0, 2 * math.pi, 0.01)
    xp = x0 + (rho * numpy.cos(angVec))
    yp = y0 + (rho * numpy.sin(angVec))

    return xp, yp
#===============================================================================

sTime = time.time()
imgName = 'CP10180_v2_BOLD_balanced_MB8_PE4_NII.nii'
#imgName = 'CP10181_v1_BOLD_balanced_MB6_PE3_Mtx100.nii.gz'
imgNameStrip = fStripExtension( imgName )
printData = False

if sys.platform == 'win32':
    imgRoot = 'R:\\nifti\\QC_Tests\\'
    outputDir = imgRoot
else:
    imgRoot = '/home/NRG/jwilso01/nifti/QC_Tests'
    outputDir = imgRoot
    
niftiFile = os.path.normpath(imgRoot +os.sep+ imgName)
nimBabel = nib.load(niftiFile)
#nimBabelData = nimBabel.get_data()
nimBabelData = numpy.asarray(nimBabel.get_data()[:,:,:,0], dtype=numpy.float32) 
nimBabelData = ( nimBabelData / numpy.max(nimBabelData.ravel())) * 255 
nimBabelDataSz = numpy.asarray(nimBabelData.shape)
print numpy.max(nimBabelData.ravel())
    
#ringImg = cv2.imread(imgRoot + imgName, cv.CV_LOAD_IMAGE_GRAYSCALE)
#ringImgSz = numpy.asarray(ringImg.shape)
#print type(ringImg), ringImgSz, cv.CV_HOUGH_GRADIENT


#ringImg = cv.LoadImageM(imgRoot + imgName, cv.CV_LOAD_IMAGE_GRAYSCALE)
#ringStorage = cv.CreateMat(ringImgSz[0], 1, cv.CV_32FC3)
ringStorage = numpy.zeros([nimBabelDataSz[0], 1])

#for i in xrange(0, nimBabelDataSz[2]):
for i in xrange(0, 35):

    ringImg = numpy.asarray(nimBabelData[:,:,i], dtype=numpy.uint8)
    print numpy.min(ringImg.ravel()), numpy.max(ringImg.ravel()), numpy.min(nimBabelData.ravel()), numpy.max(nimBabelData.ravel())
    
    
    
    ringCoords = cv2.HoughCircles(ringImg, cv.CV_HOUGH_GRADIENT, 1, 0.0001, ringStorage, 100, 40);
    
    if (ringCoords is not None):
        
        if printData:
            fPrintData( ringCoords, imgNameStrip + '_' + str(i), outputDir )

        implot = pyplot.imshow(ringImg)
        

        ringCoordsSz = numpy.asarray(ringCoords.shape)
        for j in xrange(ringCoordsSz[1]):
            
            pX, pY = fCircle( ringCoords[0, j, 0], ringCoords[0, j, 1], ringCoords[0, j, 2] )
            pyplot.plot(pX, pY, 'r-')
            
            
#            cv2.circle(ringImg, (ringCoords[0, j, 0], ringCoords[0, j, 1]), ringCoords[0, j, 2], cv.Scalar(0, 0, 255))
#            
#            radius = ringCoords[0, j, 2]
#            center = (ringCoords[0, j, 0], ringCoords[0, j, 1])
#            
#            print (radius, center)
#        
#        cv2.imshow('Image', ringImg)
#        cv2.waitKey(0)
        pyplot.show()
    print 'looping...' + str(i)
#print ringImg.height, ringImg.width 
#print ringImg

#pyplot.imshow( ringImg, holdon=True )
#pyplot.show()


print("Duration: %s" % (time.time() - sTime))