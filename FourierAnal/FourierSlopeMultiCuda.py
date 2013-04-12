'''
Created on 2013-03-22

@author: jwilso01
'''


import scipy.misc
#import numpy.fft as nfft
import multiprocessing

from pyfft.cuda import Plan
from pycuda.tools import make_default_context
import pycuda.tools as pytools
import pycuda.gpuarray as garray
import pycuda.driver as drv

import os
import csv
import time
import numpy
import scipy
import socket
import ctypes
import argparse
import nipy as nip
import nibabel as nib
from scipy import stats
import matplotlib.pyplot as pyplot

import multiprocessing as mp
from multiprocessing import Pool, Process, Queue
from multiprocessing.sharedctypes import Value, Array
#import pydevd;pydevd.settrace()

#===============================================================================
# -D R:\nifti\QC_Tests  -N  CP10086_v3_T2w2_0_7mm.nii.gz -O C:\tmp\bar
# -D C:\Users\jwilso01\Downloads  -N 111312_strc_T1w_MPR1.nii.gz -O C:\tmp\bar
#===============================================================================

#===============================================================================
# PARSE INPUT
#===============================================================================
parser = argparse.ArgumentParser(description="Test program to do Fourier slope statistics on a nifit image...")
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


print "Running on " + socket.gethostname()
freqCutoffLow = 2;
freqCutoffHigh = 60;
plotSlope = False
plotResults = False
printData = True





inputFileName = niftiFile
niftiFile = niftiDir +os.sep+ niftiFile
print "File: " + niftiFile

numpyNipyData = nip.load_image(niftiFile)
numpyNipyData = numpy.float32(numpy.asarray(numpyNipyData))
numpyNipyDataSz = numpy.asarray(numpyNipyData.shape)

CenterPoint = numpyNipyDataSz / 2
ndMesh =  numpy.asarray(numpy.mgrid[0:numpyNipyDataSz[0],0:numpyNipyDataSz[1]])
xGrid = (ndMesh[0,:] - CenterPoint[0])**2
yGrid = (ndMesh[1,:] - CenterPoint[1])**2 


#Volume = numpy.empty(numpyNipyData.shape)
#xGrid = numpy.empty([numpyNipyData.shape[0], numpyNipyData.shape[1]])
#yGrid = numpy.empty([numpyNipyData.shape[0], numpyNipyData.shape[1]])

coeffsArray = Array('f', numpyNipyData.shape[2], lock=False)
SliceIdxs = Array('i', range(numpyNipyData.shape[2]), lock=False)

Volume = numpyNipyData

#===============================================================================
# FUNCTIONS
#===============================================================================
def fFourierRegressCoeffs( SliceIdxs, coeffsArray ):
    global Volume, xGrid, yGrid
    for i in xrange(0, len(SliceIdxs)):
        Slice = Volume[:,:,SliceIdxs[i]]
#        print i, SliceIdxs[i], numpy.sum(Slice.ravel())
        fftImg = numpy.fft.fft2(Slice)
    #    fftImgVec = numpy.abs(numpy.real(fftImg.ravel()))

        SliceSz = Slice.shape
        xySum = (xGrid + yGrid)
        radFreqs = numpy.sqrt( xySum )
        freqIdx = numpy.nonzero( (radFreqs > freqCutoffLow) & (radFreqs < freqCutoffHigh) )
     
        plotFreqs = numpy.log10(radFreqs[freqIdx[0], freqIdx[1]])
        imageFFT = numpy.fft.fftshift(fftImg)
        imagePow = abs(imageFFT)**2
        imageZeroIdx = numpy.where(imagePow <= 0)
        if len(imageZeroIdx[0]) > 0:
            imagePow += (numpy.random.rand(SliceSz[0], SliceSz[1]) * 1e-16)
     
     
        plotPows = numpy.log10(imagePow[freqIdx[0], freqIdx[1]])
        
        #slope, intercept, r-value, p-value, stderr
        regressCoeffs = stats.linregress( plotFreqs, plotPows )
                
        coeffsArray[SliceIdxs[i]] = regressCoeffs[0]
#        print regressCoeffs
    return coeffsArray
#===============================================================================

if __name__ == "__main__":

#    print mp.cpu_count()
    sTime = time.time()
    p = Process(target=fFourierRegressCoeffs, args=(SliceIdxs, coeffsArray))
    p.start()
    p.join()
    print("Duration: %s" % str(time.time() - sTime))
    
    print p, coeffsArray[:]
    
#===============================================================================
# http://www.sqlservercentral.com/blogs/sqlwise/2012/12/21/a-really-simple-multiprocessing-python-example/
# http://www.scipy.org/ParallelProgramming
#===============================================================================

#===============================================================================
# from multiprocessing import Pool 
# import sys
# import os
# import os.path
# import Image
# resize_factor = 0.5
# dest = os.getcwd()
# 
# def resize(x):
# 
# try:
# # Attempt to open an image file
# filepath = x
# image = Image.open(filepath)
# except IOError, e:
# # Report error, and then skip to the next argument
# print "Problem opening", filepath, ":", e
# return
# h,w = image.size
# h,w = (int(h * resize_factor), int(w * resize_factor))
# # Resize the image
# image = image.resize((h,w), Image.ANTIALIAS)
# fname = os.path.basename(filepath)
# # Split our original filename into name and extension
# (name, extension) = os.path.splitext(fname)
# # Save the thumbnail as "(original_name)_thumb.png"
# image.save(os.path.join(dest,name + '.jpg'),quality=80)
# image = None
# 
# 
# 
# if __name__ == '__main__':
# core_ct = os.sysconf('SC_NPROCESSORS_ONLN')
# pool = Pool(processes=core_ct)
# pool.map(resize,sys.argv[1:])
# pool.close()
# pool.join()
#===============================================================================

#===============================================================================
# from multiprocessing import Pool
# import time
# import numpy
# 
# def takeuptime(ntrials):
#    for ii in xrange(ntrials):
#        junk = numpy.std(numpy.random.randn(1e5))
#    return junk
# 
# if __name__ == "__main__":
#    start = time.time()
#    map(takeuptime, [500, 500])
#    print "Serial time: %f" % (time.time() - start)
# 
#    start = time.time()
#    pool = Pool(processes=2)
#    pool.map(takeuptime, [500, 500])
#    print "Parallel time: %f" % (time.time() - start)
#===============================================================================
    








