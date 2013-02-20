'''
Created on 2012-06-11

@author: jwilso01
'''

import os
import sys
#import csv
#import math
import time
import numpy 
import nipy as nip
import nibabel as nib
#import matplotlib.pyplot as pyplot
#import matplotlib.animation as animation
import argparse

#===============================================================================
# PARSE INPUT
#===============================================================================
parser = argparse.ArgumentParser(description="Program to do brain extraction as in PreFreeSurferPipeline...")
parser.add_argument("-D", "--ImageDir", dest="niftiDir", default=os.getcwd(), type=str, help="specify nifti directory...")
parser.add_argument("-N", "--ImageName", dest="niftiFileName", type=str, help="specify nifti file...")
parser.add_argument("-Y", "--ImageType", dest="niftiFileType", default="T1w", type=str, help="specify nifti file...")
parser.add_argument("-S", "--Subject", dest="Subject", type=str, help="specify subject file...")
parser.add_argument("-B", "--BrainSize", dest="BrainSize", type=str, help="specify brain size...")
parser.add_argument("-t", "--Template", dest="TemplateFile", type=str, help="specify template file name...")
#additions for FNIRT...
parser.add_argument("-m", "--MaskFile", dest="MaskFile", type=str, help="specify mask file file name...")
parser.add_argument("-T", "--Template2mm", dest="Template2mm", type=str, help="specify template file file name...")
parser.add_argument("-M", "--MaskFile2mm", dest="MaskFile2mm", type=str, help="specify mask file file name...")
parser.add_argument("-C", "--FNIRTconfig", dest="FNIRTconfig", type=str, help="specify FNIRT config file name...")

InputArgs = parser.parse_args()
niftiDir = InputArgs.niftiDir
niftiFileName = InputArgs.niftiFileName
niftiFileType = InputArgs.niftiFileType
Subject = InputArgs.Subject
BrainSize = InputArgs.BrainSize
TemplateFile = InputArgs.TemplateFile
MaskFile = InputArgs.MaskFile
Template2mmFile = InputArgs.Template2mm
MaskFile2mm = InputArgs.MaskFile2mm
FNIRTconfig = InputArgs.FNIRTconfig

#===============================================================================
def fPrintData( inputData, Subject, outputDir ):
    
    inputDataSz = numpy.asarray( inputData.shape )
    headerStr = 'BrainVolume'
    fileName = Subject + '_BrainVolume.txt'

    fileID = open(outputDir +os.sep+ fileName, 'wb')
    fileID.write(headerStr + "\n")
    fileID.write('%.8f' % inputData + "\n")
#===============================================================================

sTime = time.time()

niftiNoExtension, niftiExtension = os.path.splitext(niftiFileName)
if niftiExtension == '.gz':
    niftiNoExtension, niftiExtension = os.path.splitext(niftiNoExtension)
    


print  sys.platform
if sys.platform == 'win32':
    PipelineScripts = 'R:\\dev\\Pipelines\\PreFreeSurfer\\scripts\\'
    GlobalScripts = 'R:\\dev\\Pipelines\\global\\scripts\\'
    TemplateDir = 'R:\\dev\\Pipelines\\global\\templates\\'
    ConfigDir = 'R:\\dev\\Pipelines\\global\\config\\'
else:
    PipelineScripts = '/home/NRG/jwilso01/dev/Pipelines/PreFreeSurfer/scripts/'
    GlobalScripts = '/nrgpackages/tools/HCP/scripts/'
    TemplateDir = '/nrgpackages/atlas/HCP/'
    ConfigDir = '/nrgpackages/tools/HCP/conf/'
    
#===============================================================================
# build input nifti path and file...
#===============================================================================
niftiPath = niftiDir + os.sep + Subject + os.sep + niftiFileType
niftiPathName = niftiPath + os.sep + niftiFileName
nimBabel = nib.load(niftiPathName)
#nimData = nimBabel.get_data()
    
#print nimBabel.header_class.binaryblock, nimBabel.file_map, 
print nimBabel.get_sform() , nimBabel.get_header(), nimBabel.to_file_map()
img = nib.AnalyzeImage( nimBabel, nimBabel.get_sform(), header=nimBabel.get_header(), extra=None, file_map=nimBabel.to_file_map() )
nib.loadsave.save(img, 'foo')

#nipNipyData = nip.load_image(niftiFile)
#nipNipyDataSz = nipNipyData.shape

#===============================================================================
# set up directories...
#===============================================================================
xfmsDir = niftiDir + os.sep + Subject + os.sep + niftiFileType + os.sep + 'xfms'
if not os.path.exists(xfmsDir):
    os.makedirs(xfmsDir)
    
ACPCdir = niftiDir + os.sep + Subject + os.sep + niftiFileType + os.sep + 'ACPCAlignment' 
if not os.path.exists(ACPCdir):
    os.makedirs(ACPCdir)
    
BrainExtractFNIRTdir = niftiDir + os.sep + Subject + os.sep + niftiFileType + os.sep + 'BrainExtraction_FNIRTbased' 
if not os.path.exists(BrainExtractFNIRTdir):
    os.makedirs(BrainExtractFNIRTdir)


    
#===============================================================================
# Start of ACPC input variables...
#===============================================================================
WorkingDirectory = niftiPath +os.sep+ 'ACPCAlignment'
Input = niftiPathName
Reference = TemplateDir +os.sep+ TemplateFile
Output = niftiPath +os.sep+ niftiNoExtension +'_acpc'
OutputMatrix = niftiDir +os.sep+ Subject +os.sep+ niftiFileType + '/xfms/acpc.mat'

#acpc align T1w image to 0.8mm MNI T1wTemplate to create native volume space
#"$PipelineScripts"/ACPCAlignment.sh "$T1wFolder"/ACPCAlignment "$T1wFolder"/"$T1wImage" "$T1wTemplate" "$T1wFolder"/"$T1wImage"_acpc "$T1wFolder"/xfms/acpc.mat "$GlobalScripts" "$BrainSize"
#/home/NRG/jwilso01/dev/Pipelines/PreFreeSurfer/scripts//ACPCAlignment.sh /home/NRG/jwilso01/nifti//CP10104_v1/T1w/ACPCAlignment /home/NRG/jwilso01/nifti//CP10104_v1/T1w/T1w /nrgpackages/atlas/HCP//MNI152_T1_0.8mm.nii.gz /home/NRG/jwilso01/nifti//CP10104_v1/T1w/T1w_acpc /home/NRG/jwilso01/nifti//CP10104_v1/T1w/xfms/acpc.mat /home/NRG/jwilso01/dev/Pipelines/global/scripts/ 150
ACPCAlignment = PipelineScripts + 'ACPCAlignment.sh ' + WorkingDirectory +' '+ Input +' '+ Reference +' '+ Output +' '+ OutputMatrix +' '+ GlobalScripts +' '+ BrainSize

if sys.platform == 'win32':
    print ACPCAlignment
else:
    print ACPCAlignment
    os.system(ACPCAlignment)

#===============================================================================
# Start of BrainExtraction_FNIRTbased input variables...
#===============================================================================
#"$PipelineScripts"/BrainExtraction_FNIRTbased.sh "$T1wFolder"/BrainExtraction_FNIRTbased "$T1wFolder"/"$T1wImage"_acpc "$T1wTemplate" #"$TemplateMask" "$T1wTemplate2mm" "$Template2mmMask" "$T1wFolder"/"$T1wImage"_acpc_brain "$T1wFolder"/"$T1wImage"_acpc_brain_mask "$FNIRTConfig"
WorkingDirectory = niftiPath +os.sep+ 'BrainExtraction_FNIRTbased'
Input = niftiPath +os.sep+ niftiNoExtension +'_acpc'
Reference = TemplateDir +os.sep+ TemplateFile
ReferenceMask = TemplateDir +os.sep+ MaskFile
Reference2mm = TemplateDir +os.sep+ Template2mmFile
Reference2mmMask = TemplateDir +os.sep+ MaskFile2mm
OutputBrainExtractedImage = niftiPath +os.sep+ niftiNoExtension + '_acpc_brain.nii.gz' 
OutputBrainMask = niftiPath +os.sep+ niftiNoExtension + '_acpc_brain_mask.nii.gz'
FNIRTConfig = ConfigDir + FNIRTconfig


#BrainExtraction_FNIRTbased.sh /home/NRG/jwilso01/nifti//CP10104_v1/T1w/BrainExtraction_FNIRTbased  /home/NRG/jwilso01/nifti//CP10104_v1/T1w/T1w_acpc /nrgpackages/atlas/HCP//MNI152_T1_0.8mm.nii.gz /nrgpackages/atlas/HCP//MNI152_T1_0.8mm_brain_mask.nii.gz /nrgpackages/atlas/HCP//MNI152_T1_2mm.nii.gz /nrgpackages/atlas/HCP//MNI152_T1_2mm_brain_mask_dil.nii.gz /home/NRG/jwilso01/nifti//CP10104_v1/T1w/T1w_acpc_brain /home/NRG/jwilso01/nifti//CP10104_v1/T1w/T1w_acpc_brain_mask /nrgpackages/tools/HCP/conf//T1_2_MNI152_2mm.cnf
BrainExtraction_FNIRTbased = PipelineScripts + 'BrainExtraction_FNIRTbased.sh ' + WorkingDirectory +' '+ Input +' '+ Reference +' '+ ReferenceMask +' '+ Reference2mm +' '+ Reference2mmMask +' '+ OutputBrainExtractedImage +' '+ OutputBrainMask +' '+ FNIRTConfig        

if sys.platform == 'win32':
    print BrainExtraction_FNIRTbased
else: 
    print BrainExtraction_FNIRTbased
    os.system(BrainExtraction_FNIRTbased)


nimBabel = nib.load(OutputBrainMask)
nimBabelHeader = nimBabel.get_header().structarr
nimBabelPixdim = nimBabelHeader['pixdim']
#nimBabelAffine = nimBabel.get_affine()
#nimBabelScale = numpy.abs([nimBabelAffine[0,0], nimBabelAffine[1,1], nimBabelAffine[2,2]])
nimBabelScale = numpy.abs([nimBabelPixdim[1], nimBabelPixdim[2], nimBabelPixdim[3]])
numpyNipyData = nip.load_image(OutputBrainMask)

fPrintData( numpy.sum(numpyNipyData) * (nimBabelScale[0] * nimBabelScale[1] * nimBabelScale[2]), Subject, niftiPath )
print 'Brain Volume: ' + str(numpy.sum(numpyNipyData) * (nimBabelScale[0] * nimBabelScale[1] * nimBabelScale[2]))

tTime = time.time() - sTime
print("Duration: %s" % tTime)

#21/07/2012 NOTE:
#1848449.10663
#Duration: 4053.55621982
#an average brain volume of 1273.6cc for men, ranging from 1052.9 to 1498.5cc, and 1131.1cc for women, ranging from 974.9 to 1398.1cc.
#from FSL: import 3610252 1848449.250000

 

#===============================================================================
# Example input...
#===============================================================================
#-D R:\nifti -N CP10104_v1_T1w2.nii.gz -Y T1w -S CP10104_v1 -B 150 -t MNI152_T1_0.8mm.nii.gz -m MNI152_T1_0.8mm_brain_mask.nii.gz -T MNI152_T1_2mm.nii.gz -M MNI152_T1_2mm_brain_mask_dil.nii.gz -C T1_2_MNI152_2mm.cnf
#python BrainExtractionHCP.py -D /home/NRG/jwilso01/nifti -N CP10104_v1_T1w2.nii.gz -Y T1w -S CP10104_v1 -B 150 -t MNI152_T1_0.8mm.nii.gz -m MNI152_T1_0.8mm_brain_mask.nii.gz -T MNI152_T1_2mm.nii.gz -M MNI152_T1_2mm_brain_mask_dil.nii.gz -C T1_2_MNI152_2mm.cnf

