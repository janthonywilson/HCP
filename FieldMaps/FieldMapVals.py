'''
Created on 2012-07-21

@author: jwilso01
'''
import os
import csv
import sys
import time
import numpy
import scipy
import argparse
import nipy as nip
import nibabel as nib
from scipy import stats
import matplotlib.pyplot as pyplot

#===============================================================================
# PARSE INPUT
#===============================================================================
parser = argparse.ArgumentParser(description="Test program to do wavelet statistics on a nifit image...")
parser.add_argument("-D", "--ImageDir", dest="niftiDir", default=os.getcwd(), type=str, help="specify nifti directory...")
parser.add_argument("-S", "--Subject", dest="Subject", type=str, help="specify subject for prepending to output file...")
parser.add_argument("-M", "--MagnitudeImage", dest="magFile", type=str, help="specify nifti magnitudue file...")
parser.add_argument("-P", "--PhaseImage", dest="phaFile", type=str, help="specify nifti phase file...")
parser.add_argument("-O", "--OutputDir", dest="outputDir", type=str, help="specify where to write the output files...")
parser.add_argument("-TE", "--TE", dest="TE", default=2.46, type=float, help="specify TE for fieldmaps...")
InputArgs = parser.parse_args()
niftiDir = InputArgs.niftiDir
magFile = InputArgs.magFile
phaFile = InputArgs.phaFile
outputDir = InputArgs.outputDir
TE = InputArgs.TE
Subject = InputArgs.Subject
#===============================================================================

sTime = time.time()
printData = True
plotData = True

print  sys.platform
if sys.platform == 'win32':
    PipelineScripts = 'R:\\dev\\Pipelines\\PreFreeSurfer\\scripts\\'
    GlobalScripts = 'R:\\dev\\Pipelines\\global\\scripts\\'
    TemplateDir = 'R:\\dev\\Pipelines\\global\\templates\\'
    ConfigDir = 'R:\\dev\\Pipelines\\global\\config\\'
    FSLDir = 'NoWinFSL\\'
else:
    PipelineScripts = '/home/NRG/jwilso01/dev/Pipelines/PreFreeSurfer/scripts/'
    GlobalScripts = '/nrgpackages/tools/HCP/scripts/'
    TemplateDir = '/nrgpackages/atlas/HCP/'
    ConfigDir = '/nrgpackages/tools/HCP/conf/'
    FSLDir = '/nrgpackages/tools/fsl-nrg/'
    
if not os.path.exists(outputDir):
    os.makedirs(outputDir)
        
#===============================================================================
def fPrintData( inputVec, Subject, outputDir ):
    
    StatsListLen = len(inputVec) 
    headerStr = 'FieldMapValues'
    fileName = Subject + '_FieldMapValues.txt'

    fileID = open(outputDir +os.sep+ fileName, 'wb')
    fileID.write(headerStr + "\n")

    for i in xrange(0, StatsListLen):
        fileID.write('%.8f' % inputVec[i] + "\n")
#===============================================================================


inputFileName = 'Fieldmap'
magFile = niftiDir +os.sep+ magFile
phaFile = niftiDir +os.sep+ phaFile
print "Input Files: " + magFile +'  '+ phaFile

nimBabel = nib.load(magFile)
nimBabelAffine = nimBabel.get_affine()
numpyMagData = nip.load_image(magFile)
numpyMagData = numpy.mean(numpy.float32(numpy.asarray(numpyMagData)), axis=3)
numpyMagDataSz = numpy.asarray(numpyMagData.shape)

#===============================================================================
# save mag matrix to nifti...
#===============================================================================
nimMagImage = nib.Nifti1Image(numpyMagData, nimBabelAffine)
nimMagImage.to_filename(outputDir +os.sep+ 'Magnitude.nii.gz')

#===============================================================================
# do BET on magnitude image...
#===============================================================================
BETstr = FSLDir + 'bin' +os.sep+ 'bet ' +outputDir +os.sep+ 'Magnitude.nii.gz ' +outputDir +os.sep+ 'Magnitude_brain.nii.gz -f .35 -m'
if sys.platform == 'win32':
    print BETstr
else: 
    print BETstr
    os.system(BETstr)

#===============================================================================
# run FMRIB script..."$GlobalScripts"/fmrib_prepare_fieldmap.sh SIEMENS "$WorkingDirectory"/Phase.nii.gz "$WorkingDirectory"/Magnitude_brain.nii.gz "$WorkingDirectory"/FieldMap.nii.gz "$TE"
#===============================================================================   
magBrain = outputDir +os.sep+ 'Magnitude_brain.nii.gz'
outFieldMap = outputDir +os.sep+ 'FieldMap.nii.gz'
FMRIB_prepare_fieldmap = GlobalScripts + 'fmrib_prepare_fieldmap.sh SIEMENS ' + phaFile +' '+ magBrain +' '+ outFieldMap +' '+ str(TE)
if sys.platform == 'win32':
    print FMRIB_prepare_fieldmap
else: 
    print FMRIB_prepare_fieldmap
    os.system(FMRIB_prepare_fieldmap)
    
numpyFieldMapData = nip.load_image(outFieldMap)
numpyFieldMapData = numpy.float32(numpy.asarray(numpyFieldMapData))
numpyFieldMapDataSz = numpy.asarray(numpyFieldMapData.shape)
print 'Fieldmap Size: ' + str(numpyMagDataSz) +'   '+ str(numpyFieldMapDataSz)

numpyFieldMapDataFlat = numpyFieldMapData.ravel()
numpyNonZsIdx = numpy.flatnonzero(numpyFieldMapData)
#print numpy.min(numpyFieldMapDataFlat[numpyNonZsIdx]), numpy.max(numpyFieldMapDataFlat[numpyNonZsIdx])

print("Duration: %s" % (time.time() - sTime))



if plotData:
    pyplot.figure()
    pyplot.hist(numpyFieldMapDataFlat[numpyNonZsIdx], bins=100)
    pyplot.figure()
    pyplot.plot( numpy.sort(numpyFieldMapDataFlat[numpyNonZsIdx]), xrange(0, len(numpyFieldMapDataFlat[numpyNonZsIdx])) )
    pyplot.show()

if printData:
    fPrintData( numpy.sort(numpyFieldMapDataFlat[numpyNonZsIdx]), Subject, outputDir )
    
#numpyPhaDataNorm = (numpyPhaData / numpy.power(2,12)) * numpy.pi
#print numpy.min(numpyPhaDataNorm[-1]), numpy.max(numpyPhaDataNorm[-1])
#numpyPhaUnwrap = numpy.unwrap(numpyPhaDataNorm) # , discont=numpy.pi
#print numpy.min(numpyPhaUnwrap[-1]), numpy.max(numpyPhaUnwrap[-1])
#pyplot.figure()
#pyplot.plot(numpy.mean(numpyNipyData, axis=2))
#pyplot.show







