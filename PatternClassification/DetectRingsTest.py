'''
Created on 2012-07-30

@author: jwilso01
'''
import os
import cv2
import cv 
import sys
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

sTime = time.time()
imgName = 'Rings3.png'

if sys.platform == 'win32':
    imgRoot = 'R:\\tmp\\'
    outputDir = imgRoot
else:
    imgRoot = '/home/NRG/jwilso01/tmp/'
    outputDir = imgRoot
    
ringImg = cv2.imread(imgRoot + imgName, cv.CV_LOAD_IMAGE_GRAYSCALE)
ringImgSz = numpy.asarray(ringImg.shape)
print type(ringImg), ringImgSz, cv.CV_HOUGH_GRADIENT


#ringImg = cv.LoadImageM(imgRoot + imgName, cv.CV_LOAD_IMAGE_GRAYSCALE)
#ringStorage = cv.CreateMat(ringImgSz[0], 1, cv.CV_32FC3)
ringStorage = numpy.zeros([ringImgSz[0], 1])

ringCoords = cv2.HoughCircles(ringImg, cv.CV_HOUGH_GRADIENT, 1, 0.0001, ringStorage, 100, 250);

fPrintData( ringCoords, imgName, outputDir )





ringCoordsSz = numpy.asarray(ringCoords.shape)
for i in xrange(ringCoordsSz[1]):
    
    cv2.circle(ringImg, (ringCoords[0, i, 0], ringCoords[0, i, 1]), ringCoords[0, i, 2], cv.Scalar(0, 0, 255))
    
    radius = ringCoords[0, i, 2]
    center = (ringCoords[0, i, 0], ringCoords[0, i, 1])
    
    print (radius, center)

cv2.imshow('Image', ringImg)
cv2.waitKey(0)
#print ringImg.height, ringImg.width 
#print ringImg

#pyplot.imshow( ringImg, holdon=True )
#pyplot.show()


print("Duration: %s" % (time.time() - sTime))