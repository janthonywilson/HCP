'''
Created on Feb 21, 2013

@author: Tony
'''
import os
import sys
#import csv
#import math
import time
import numpy 
import scipy
import subprocess
import nipy as nip
import nibabel as nib
import nipy.labs.utils.mask as mask
import argparse

sTime = time.time()
#===============================================================================
# PARSE INPUT
#===============================================================================
# -D C:\Python27\lib\site-packages\nibabel\tests\data -N example4d.nii.gz
# -D C:\Users\Tony\workspace\data\nifti -N 100307_T1w_MPR2.nii.gz
#===============================================================================
parser = argparse.ArgumentParser(description="Program to do brain extraction as in PreFreeSurferPipeline...")
parser.add_argument("-D", "--ImageDir", dest="niftiFileDir", default=os.getcwd(), type=str, help="specify nifti directory...")
parser.add_argument("-N", "--ImageName", dest="niftiFileName", type=str, help="specify nifti file...")
parser.add_argument("-Y", "--ImageType", dest="niftiFileType", default="T1w", type=str, help="specify nifti file...")
parser.add_argument("-AD", "--AnalyzeDir", dest="analFileDir", default=os.getcwd(), type=str, help="specify analyze directory...")
parser.add_argument("-AN", "--AnalyseName", dest="analFileName", type=str, help="specify nifti file...")
parser.add_argument("-S", "--SEIter", dest="seIter", type=int, default=10, help="specify iterations for SE...")
parser.add_argument("-SC", "--SEConn", dest="seConn", type=int, default=2, help="specify connectivity for SE...")

InputArgs = parser.parse_args()
niftiDir = InputArgs.niftiFileDir
niftiName = InputArgs.niftiFileName
niftiType = InputArgs.niftiFileType
seIter = InputArgs.seIter
seConn = InputArgs.seConn

analDir = InputArgs.analFileDir
analName = InputArgs.analFileName
#===============================================================================
# read the nifti file...
#===============================================================================
if (niftiDir[-1] != os.sep):
    niftiDir = niftiDir + os.sep
        
if (analDir[-1] != os.sep):
    analDir = analDir + os.sep
    
InputNameList = niftiName.split('.')
InputNameBase = InputNameList[0]

nimBabel = nib.load(niftiDir + niftiName)
nimBabelData = nimBabel.get_data()
nimBabelInfo = nib.get_info()
nimBabelAffine = nimBabel.get_affine()
nimBabelHeader = nimBabel.get_header().structarr
nimBabelPixDim = nimBabelHeader['pixdim']
numpyBabelDataSz = numpy.asarray(nimBabelData.shape)
print "Input size: " + str(numpyBabelDataSz) + ' Max: ' + str(numpy.max(nimBabelData.ravel())) + ' Min: ' + str(numpy.min(nimBabelData.ravel()))

nipMask = mask.compute_mask(nimBabel.get_data(), reference_volume=None, m=0.4, M=0.8, cc=True, opening=16, exclude_zeros=False)

img = nib.AnalyzeImage( nimBabel.get_data(), nimBabel.get_affine() )
img.to_filename(analDir + analName)



commandBET = 'resources' +os.sep+ 'bet.exe '+analDir+analName+' '+analDir+analName+ ' -m -s -n'
with open(os.devnull, "w") as fnull:
    subprocBET = subprocess.call( commandBET, stdout = fnull, stderr = fnull )



analBabel = nib.load(analDir + analName + '_mask.hdr')
analBabelData = analBabel.get_data()
analMaskDataTmp = numpy.asarray(analBabelData, dtype='int')
analMaskData = numpy.reshape(analMaskDataTmp, [analMaskDataTmp.shape[0], analMaskDataTmp.shape[1], analMaskDataTmp.shape[2]])


structElement = scipy.ndimage.generate_binary_structure(3, seConn)
structElement = scipy.ndimage.iterate_structure(structElement, seIter)
BrainMaskDilate = scipy.ndimage.binary_dilation(analMaskData, structure=structElement).astype(analMaskData.dtype)
BrainMaskOpen = scipy.ndimage.binary_opening(analMaskData, structure=structElement).astype(analMaskData.dtype)

#===========================================================================
# maskSignalImg = nib.Nifti1Image(SignalMaskEro, None, nimMask.get_header())
# nib.save(maskSignalImg, outputDir +os.sep+ 'SignalMaskEro.nii.gz')
#===========================================================================
dilImg = nib.Nifti1Image(BrainMaskDilate, None, nimBabel.get_header())
nib.save(dilImg, niftiDir +niftiName+'MaskDilate'+ str(seIter) +'_'+ str(seConn) +'.nii.gz')

openImg = nib.Nifti1Image(BrainMaskOpen, None, nimBabel.get_header())
nib.save(openImg, niftiDir +niftiName+'MaskOpen'+ str(seIter) +'_'+ str(seConn) +'.nii.gz')

nipMaskImg = nib.Nifti1Image(nipMask, None, nimBabel.get_header())
nib.save(nipMaskImg, niftiDir +niftiName+'MaskNipy.nii.gz')








print("Duration: %s" % (time.time() - sTime))