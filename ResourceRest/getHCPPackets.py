'''
Created on 2013-02-19

@author: jwilso01
'''

import os
import sys
import socket
import operator
import argparse
import itertools
import subprocess

# Time manipulation...
import time
# HCP interface class...
from pyHCP import pyHCP, getHCP, writeHCP

sTime = time.time()

#===============================================================================
# PARSE INPUT
#===============================================================================
# Examples:
# -User tony -Password passfoo -Server hcpx-dev-cuda00.nrg.mir -Project HCP_Q1 -Subject 100307 -Session 100307_3T -Scan 108 -FileType NIFTI -OutputDir C:\tmp\HCP -Strip strc
# -User tony -Password passfoo -Server intradb.humanconnectome.org -Project HCP_Phase2 -Subject 100307 -Session 100307_strc -Scan 19 -FileType NIFTI -OutputDir C:\tmp\HCP -Strip None
# -User tony -Password passfoo -Server hcpx-dev-cuda00.nrg.mir -Project HCP_Q1 -Subject 100307 -Session 100307_3T -Scan 108 -FileType NIFTI -OutputDir /home/NRG/jwilso01/tmp/HCP -Strip 3T
# -User tony -Password passfoo -Server hcpx-demo3.humanconnectome.org -Project HCP_Q1 -Subject 100307 -Session 100307_3T -Resource tfMRI_WM_RL_unproc -FileType NIFTI -OutputDir C:\tmp\HCP -Strip None
# -User tony -Password passfoo -Server hcpx-dev-cuda00.nrg.mir -Project HCP_Q1 -Subject 100307 -Session 100307_3T -Resource tfMRI_WM_RL_unproc -ResourceRoot LINKED_DATA -FileType NIFTI,TXT -OutputDir C:\tmp\HCP -Strip None
# -User tony -Password passfoo -Server hcpx-demo.humanconnectome.org -Project HCP_Q1 -Subject 100307 -Session 100307_3T -Resource tfMRI_WM_RL_unproc -ResourceRoot tfMRI_WM_RL_unproc -OutputDir C:\tmp\HCP -Strip None
#===============================================================================
parser = argparse.ArgumentParser(description="Program to get files from resources for HCP pipelines...")
# input...
parser.add_argument("-User", "--User", dest="User", default='tony', type=str)
parser.add_argument("-Password", "--Password", dest="Password", type=str)
parser.add_argument("-Project", "--Project", dest="Project", type=str, default=None, help="pick project")
parser.add_argument("-Subject", "--Subject", dest="Subject", type=str, default=None, help="pick subject")
parser.add_argument("-Package", "--Package", dest="Package", type=str, default=None, help="pick package")
parser.add_argument("-FunctionalSeries", "--FunctionalSeries", dest="FunctSeries", type=str, default=None, help="pick functional series")
#parser.add_argument("-PackageRoot", "--PackageRoot", dest="PackageRoot", type=str, default=None, help="pick package root")
parser.add_argument("-Strip", "--Strip", dest="StripSession", type=str, default=None, help="strip session information")
# remote...
parser.add_argument("-Server", "--Server", dest="Server", type=str, default="https://intradb.humanconnectome.org", help="pick server")
# output...
parser.add_argument("-OutputDir", "--OutputDir", dest="OutputDir", type=str, default=None, help="specify location to write output")
parser.add_argument("-OutputFile", "--OutputFile", dest="OutputFile", type=str, default=None, help="specify filename to write output")
parser.add_argument("-OutputDirFile", "--OutputDirFile", dest="OutputDirFile", type=str, default=None, help="specify filename and location to write output")
# timeout...
parser.add_argument("-t", "--time_out", dest="Timeout", type=float, default=16.0, help="change timeout")

# version...
parser.add_argument('--version', action='version', version='%(prog)s 0.0.1')

args = parser.parse_args()
Server = args.Server
User = args.User
Password = args.Password

Project = args.Project
Subject = args.Subject
Package = args.Package
#PackageRoot = args.PackageRoot

FunctSeries = args.FunctSeries



StripSession = args.StripSession
if (StripSession == 'None'): StripSession = None

OutputFile = args.OutputFile
OutputDir = args.OutputDir
OutputDirFile = args.OutputDirFile

if (OutputDir is not None):
    OutputDir = os.path.normpath(OutputDir)+os.sep
else:
    print 'ERROR: Output directory not specified.'
    sys.exit()
    
if (OutputDirFile is not None):
    OutputDirFile = os.path.normpath(args.OutputDirFile)
elif (OutputDir is None) and (OutputDirFile is None):
    print 'ERROR: No output directory location specified.'
    sys.exit()
    
    
if not os.path.exists(OutputDir):
    os.makedirs(OutputDir)
#===============================================================================

print "Running %s on %s with IP %s" % (os.path.split(sys.argv[0])[1], socket.gethostname(), socket.gethostbyname(socket.gethostname()))

pyHCP = pyHCP(User, Password, Server)
getHCP = getHCP(pyHCP)
writeHCP = writeHCP(getHCP, OutputDir)
writeHCP.Flatten = False

getHCP.Project = Project
getHCP.Subject = Subject

SubjectResourcesMeta = getHCP.getSubjectResourcesMeta()

FileName = list()
FileURIName = list()
FilePathName = list()
FilePathNameReadable = list()
#===============================================================================
# STRUCTURAL...
#===============================================================================
if (Package == 'Structural'):
    
    if (Server.find('intradb') != -1) or (Server.find('hcpi') != -1):
        ResourceRoot = 'Details'
    else:
        ResourceRoot = 'Structural_preproc'
        
    getHCP.Resource = ResourceRoot
    getHCP.Session = SubjectResourcesMeta.get('Session')[SubjectResourcesMeta.get('Label').index(ResourceRoot)]
    ResourceMeta = getHCP.getSubjectResourceMeta()
    
    SourceDir = 'T1w'
    T1wList = list()
    T1wList.append('%s/%s/T1w_acpc_dc.nii.gz' % (ResourceRoot, SourceDir))
    T1wList.append('%s/%s/T1w_acpc_dc_restore.nii.gz' % (ResourceRoot, SourceDir))
    T1wList.append('%s/%s/T1w_acpc_dc_restore_brain.nii.gz' % (ResourceRoot, SourceDir))
    T1wList.append('%s/%s/T2w_acpc_dc.nii.gz' % (ResourceRoot, SourceDir))
    T1wList.append('%s/%s/T2w_acpc_dc_restore.nii.gz' % (ResourceRoot, SourceDir))
    T1wList.append('%s/%s/T2w_acpc_dc_restore_brain.nii.gz' % (ResourceRoot, SourceDir))
    T1wList.append('%s/%s/T1wDividedByT2w.nii.gz' % (ResourceRoot, SourceDir))
    T1wList.append('%s/%s/T1wDividedByT2w_ribbon.nii.gz' % (ResourceRoot, SourceDir))
    T1wList.append('%s/%s/brainmask_fs.nii.gz' % (ResourceRoot, SourceDir))
    T1wList.append('%s/%s/wmparc.nii.gz' % (ResourceRoot, SourceDir))
    T1wList.append('%s/%s/ribbon.nii.gz' % (ResourceRoot, SourceDir))
    T1wList.append('%s/%s/BiasField_acpc_dc.nii.gz' % (ResourceRoot, SourceDir))
    
    SourceDir = 'MNINonLinear/xfms'
    xfmList = list()
    xfmList.append('%s/%s/acpc_dc2standard.nii.gz' % (ResourceRoot, SourceDir))
    xfmList.append('%s/%s/NonlinearRegJacobians.nii.gz' % (ResourceRoot, SourceDir))
    xfmList.append('%s/%s/standard2acpc_dc.nii.gz' % (ResourceRoot, SourceDir))
    
    # SourceFile = 'ALL'
    SourceList = list()
    SourceList.append('%s/%s/' % (ResourceRoot, 'MNINonLinear'))
    SourceList.append('%s/%s/' % (ResourceRoot, 'MNINonLinear/Native'))
    SourceList.append('%s/%s/' % (ResourceRoot, 'MNINonLinear/fsaverage_LR32k'))
    SourceList.append('%s/%s/' % (ResourceRoot, 'T1w/Native'))

    
    for i in xrange(0, len(ResourceMeta.get('Path'))):
        for j in xrange(0, len(T1wList)):
            if (ResourceMeta.get('Path')[i].find(T1wList[j]) != -1) and (ResourceMeta.get('Path')[i] not in FilePathName):
                FileName.append(ResourceMeta.get('Name')[i])
                FilePathName.append(ResourceMeta.get('Path')[i])
                FilePathNameReadable.append(ResourceMeta.get('Readable')[i])
                
                FileURIName.append(ResourceMeta.get('URI')[i])
                
        for j in xrange(0, len(xfmList)):
            if (ResourceMeta.get('Path')[i].find(xfmList[j]) != -1) and (ResourceMeta.get('Path')[i] not in FilePathName):
                FileName.append(ResourceMeta.get('Name')[i])
                FilePathName.append(ResourceMeta.get('Path')[i])
                FilePathNameReadable.append(ResourceMeta.get('Readable')[i])
                
                FileURIName.append(ResourceMeta.get('URI')[i])
                
        # SourceList/* has stronger constraint, thus .find(SourceList[j]+ResourceMeta.get('Name')[i])
        for j in xrange(0, len(SourceList)):
            if (ResourceMeta.get('Path')[i].find(SourceList[j]+ResourceMeta.get('Name')[i]) != -1) and (ResourceMeta.get('Path')[i] not in FilePathName):
                ResourceMetaPathSplit = ResourceMeta.get('Path')[i].split('/')
                FileName.append(ResourceMeta.get('Name')[i])
                FilePathName.append(ResourceMeta.get('Path')[i])
                FilePathNameReadable.append(ResourceMeta.get('Readable')[i])
                
                FileURIName.append(ResourceMeta.get('URI')[i])
#===============================================================================
# DIFFUSION...
#===============================================================================
elif (Package == 'Diffusion'):
    
    if (Server.find('intradb') != -1) or (Server.find('hcpi') != -1):
        ResourceRoot = 'Diffusion'
    else:
        ResourceRoot = 'Diffusion_preproc'
        
    getHCP.Resource = ResourceRoot
    getHCP.Session = SubjectResourcesMeta.get('Session')[SubjectResourcesMeta.get('Label').index(ResourceRoot)]
    ResourceMeta = getHCP.getSubjectResourceMeta()
        
    SourceDir = 'Diffusion/data'
    DataList = list()
    DataList.append('%s/%s/bvals' % (ResourceRoot, SourceDir))
    DataList.append('%s/%s/bvecs' % (ResourceRoot, SourceDir))
    DataList.append('%s/%s/data.nii.gz' % (ResourceRoot, SourceDir))
    DataList.append('%s/%s/nodif_brain_mask.nii.gz' % (ResourceRoot, SourceDir))
    DataList.append('%s/%s/grad_dev.nii.gz' % (ResourceRoot, SourceDir))

    
    SourceDir = 'T1w/xfms'
    T1wXfmList = list()
    T1wXfmList.append('%s/%s/diff2str.mat' % (ResourceRoot, SourceDir))
    T1wXfmList.append('%s/%s/str2diff.mat' % (ResourceRoot, SourceDir))
    
    SourceDir = 'MNINonLinear/xfms'
    MNIXfmList = list()
    MNIXfmList.append('%s/%s/diff2standard.nii.gz' % (ResourceRoot, SourceDir))
    MNIXfmList.append('%s/%s/standard2diff.nii.gz' % (ResourceRoot, SourceDir))
    
    AllList = DataList + T1wXfmList + MNIXfmList
  
    for i in xrange(0, len(ResourceMeta.get('Path'))):
        for j in xrange(0, len(AllList)):
            if (ResourceMeta.get('Path')[i].find(AllList[j]) != -1) and (ResourceMeta.get('Path')[i] not in FilePathName):
                FileName.append(ResourceMeta.get('Name')[i])
                FilePathName.append(ResourceMeta.get('Path')[i])
                FilePathNameReadable.append(ResourceMeta.get('Readable')[i])
                
                FileURIName.append(ResourceMeta.get('URI')[i])

#===============================================================================
# FUNCTIONAL...
#===============================================================================
elif (Package == 'Functional'):
    
    if (Server.find('intradb') != -1):
        ResourceRoot = 'Functional'
    else:
        ResourceRoot = 'Functional_preproc'
#===============================================================================
# FIX...
#===============================================================================
#${SubjectID}/MNINonLinear/Results/${fMRIName}/${fMRIName}_hp2000_clean.nii.gz
#${SubjectID}/MNINonLinear/Results/${fMRIName}/${fMRIName}_Atlas_hp2000_clean.dtseries.nii
#
#${SubjectID}/MNINonLinear/Results/${fMRIName}/${fMRIName}_hp2000.ica/.fix
#${SubjectID}/MNINonLinear/Results/${fMRIName}/${fMRIName}_hp2000.ica/fix4melview_HCP_hp2000_thr5.txt
#${SubjectID}/MNINonLinear/Results/${fMRIName}/${fMRIName}_hp2000.ica/filtered_func_data.ica/*
#===============================================================================
elif (Package == 'FIX'):

    if (Server.find('intradb') != -1) or (Server.find('hcpi') != -1):
        ResourceRoot = 'FIX'
    else:
        ResourceRoot = 'FIX_preproc'
        
    getHCP.Resource = ResourceRoot
    getHCP.Session = SubjectResourcesMeta.get('Session')[SubjectResourcesMeta.get('Label').index(ResourceRoot)]
    ResourceMeta = getHCP.getSubjectResourceMeta()
#    print '\n'.join(ResourceMeta.get('Path'))
    
    SourceDir = '%s/%s' % (ResourceRoot, FunctSeries)
    FIXList = list()
    FIXList.append('%s/%s_hp2000_clean.nii.gz' % (SourceDir, FunctSeries))
    FIXList.append('%s/%s_Atlas_hp2000_clean.dtseries.nii' % (SourceDir, FunctSeries))
    
    SourceDir = '%s/%s/%s_hp2000.ica' % (ResourceRoot, FunctSeries, FunctSeries)
    FIXList.append('%s/.fix' % (SourceDir))
    FIXList.append('%s/fix4melview_HCP_hp2000_thr5.txt' % (SourceDir))
    
    # SourceFile = 'ALL'
    SourceDir = '%s/%s/%s_hp2000.ica/filtered_func_data.ica' % (ResourceRoot, FunctSeries, FunctSeries)
    FIXAllList = list()
    FIXAllList.append('%s/' % (SourceDir))
    SourceDir = '%s/%s/%s_hp2000.ica/filtered_func_data.ica/stats' % (ResourceRoot, FunctSeries, FunctSeries)
    FIXAllList.append('%s/' % (SourceDir))
    SourceDir = '%s/%s/%s_hp2000.ica/filtered_func_data.ica/report' % (ResourceRoot, FunctSeries, FunctSeries)
    FIXAllList.append('%s/' % (SourceDir))
    
    
    for i in xrange(0, len(ResourceMeta.get('Path'))):
        for j in xrange(0, len(FIXList)):
            if (ResourceMeta.get('Path')[i].find(FIXList[j]) != -1) and (ResourceMeta.get('Path')[i] not in FilePathName):
                FileName.append(ResourceMeta.get('Name')[i])
                FilePathName.append(ResourceMeta.get('Path')[i])
                FilePathNameReadable.append(ResourceMeta.get('Readable')[i])
                
                FileURIName.append(ResourceMeta.get('URI')[i])
        # SourceList/* has stronger constraint, thus .find(SourceList[j]+ResourceMeta.get('Name')[i]), without this a list all subdirs deep would be generated.
        for j in xrange(0, len(FIXAllList)):
            if (ResourceMeta.get('Path')[i].find(FIXAllList[j]+ResourceMeta.get('Name')[i]) != -1) and (ResourceMeta.get('Path')[i] not in FilePathName):
                ResourceMetaPathSplit = ResourceMeta.get('Path')[i].split('/')
                FileName.append(ResourceMeta.get('Name')[i])
                FilePathName.append(ResourceMeta.get('Path')[i])
                FilePathNameReadable.append(ResourceMeta.get('Readable')[i])
                
                FileURIName.append(ResourceMeta.get('URI')[i])

#===============================================================================
# END...NOW PRINT
#===============================================================================        
        
        
        
#for i in xrange(0, len(FilePathName)):
#    print FilePathName[i]
    
if all(FilePathNameReadable):
    print ('cp %s %s' % (OutputDir, FilePathName))
    writeHCP.writeFileFromPath(FilePathName, FileName)
else:
    print ('Destination: %s  URI: %s' % (OutputDir, FileURIName))
    writeHCP.writeFileFromURL(getHCP, FileURIName, FileName)
    print 'Streamed bytes: ' + ', '.join(map(str, writeHCP.BytesStream))
    

print 'Written bytes: ' + ', '.join(map(str, writeHCP.BytesWrite))
print 'Delta bytes: ' + ', '.join(map(str, list(itertools.imap(operator.sub, writeHCP.BytesStream, writeHCP.BytesWrite))))
print("Duration: %s" % (time.time() - sTime))

    
#===============================================================================

    


    
    
    
    
    
    