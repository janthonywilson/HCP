'''
Created on 2013-02-19

@author: jwilso01
'''

import os
import sys
import socket
import argparse
import subprocess
# Time manipulation...
import time
# HCP interface class...
from pyHCP import pyHCP, getHCP, writeHCP

sTime = time.time()
#===============================================================================
# TO DO:
#===============================================================================
# -Scan should accept a list...
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
# -Server https://hcpx-dev-cuda00.nrg.mir/ -User tony -Password passfoo -Project HCP_Q1 -Subject 100307 -Session 100307_3T -Scan 109 -FileType NIFTI -OutputFile 100307_T1w_MPR1.nii.gz -OutputDir C:\tmp
# DIFFUSION: -Files T1w_acpc_dc.nii.gz,T1w_acpc_dc_restore.nii.gz,T1w_acpc_dc_restore_brain.nii.gz,BiasField_acpc_dc.nii.gz,brainmask_fs.nii.gz,
#===============================================================================
parser = argparse.ArgumentParser(description="Program to get files from resources for HCP pipelines...")
# input...
parser.add_argument("-User", "--User", dest="User", default='tony', type=str)
parser.add_argument("-Password", "--Password", dest="Password", type=str)
parser.add_argument("-Project", "--Project", dest="Project", type=str, default=None, help="pick project")
parser.add_argument("-Subject", "--Subject", dest="Subject", type=str, default=None, help="pick subject")
parser.add_argument("-Session", "--Session", dest="Session", type=str, default=None, help="pick session")
parser.add_argument("-Scan", "--Scan", dest="Scan", type=str, default=None, help="pick scan")
parser.add_argument("-Resource", "--Resource", dest="Resource", type=str, default=None, help="pick resource")
parser.add_argument("-ResourcePath", "--ResourcePath", dest="ResourcePath", type=str, default=None, help="pick resource path off of resource(root)")
parser.add_argument("-Files", "--Files", dest="Files", type=str, default=None, help="pick resource file(s)")
parser.add_argument("-Flatten", "--Flatten", dest="Flatten", default=True, help="bool to flatten write")
parser.add_argument("-FileType", "--FileType", dest="FileType", type=str, default=None, help="pick filetype for copy")
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
User = args.User
Password = args.Password
Project = args.Project
Subject = args.Subject
Session = args.Session
Scan = args.Scan
Files = args.Files
FileType = args.FileType
Flatten = args.Flatten
Server = args.Server
Resource = args.Resource
ResourcePath = args.ResourcePath
StripSession = args.StripSession
if (StripSession == 'None'): StripSession = None

# OUTPUT...
OutputDir = args.OutputDir
OutputFile = args.OutputFile
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
    
if (FileType is not None) and (FileType.find(',') != -1):
    FileType = FileType.split(',')
    
if (Files is not None) and (Files.find(',') != -1):
    Files = Files.split(',')
    
if (ResourcePath is None):
    ResourcePath = Resource
    
if (Flatten == 'False'):
    Flatten = False
else:
    Flatten = True
    
if not os.path.exists(OutputDir):
    os.makedirs(OutputDir)
#===============================================================================

print "Running %s on %s with IP %s" % (os.path.split(sys.argv[0])[1], socket.gethostname(), socket.gethostbyname(socket.gethostname()))

pyHCP = pyHCP(User, Password, Server)
getHCP = getHCP(pyHCP)
writeHCP = writeHCP(getHCP, OutputDir)
writeHCP.Flatten = Flatten

getHCP.Project = Project
getHCP.Subject = Subject
getHCP.Session = Session
getHCP.Scan = Scan

FileName = list()
FileURIName = list()
FilePathName = list()
FilePathNameReadable = list()

OutputFileName = list()
OutputFileURIName = list()
OutputFilePathName = list()

if (Resource is not None):
    
    #===========================================================================
    # check for session on input, if not present do the full Resources call, this is a resource hog...
    #===========================================================================
    if (getHCP.Session == None):
        SubjectResourcesMeta = getHCP.getSubjectResourcesMeta()
        getHCP.Session = SubjectResourcesMeta.get('Session')[SubjectResourcesMeta.get('Label').index(Resource)]
        if (Resource in SubjectResourcesMeta.get('Label')):
            getHCP.Resource = Resource
            ResourceMeta = getHCP.getSubjectResourceMeta()
        else:
            print 'ERROR: %s not in %s resources.' % (Resource, Subject)
            sys.exit()
                
    else:
        getHCP.Resource = Resource
        ResourceMeta = getHCP.getSubjectResourceMeta()
    
    for i in xrange(0, len(ResourceMeta.get('Path'))):
        
        #=======================================================================
        # if entire RESOURCE is requested...
        #=======================================================================
        if (Files == None):
            if (ResourceMeta.get('Path')[i].find(Resource +'/'+ ResourcePath +'/') > 0) and (ResourcePath[-1] == '/*'):
                ResourceMetaList = ResourceMeta.get('Path')[i].split('/')
    
                    #===========================================================
                    # <!--For diffusion...-->
                    # <!--T1w/*, T1w/xfms/*, T1w/<subject>/surf/*, T1w/<subject>/label/*, T1w/<subject>/mri/*, T1w/<subject>/mri/transforms/*, T1w/<subject>/stats/*, T1w/<subject>/touch/*-->
                    # <!--MNINonLinear/*, MNINonLinear/xfms/*-->
                    #===========================================================
                if (Resource in ResourceMetaList):
                    FilePathNameReadable.append(ResourceMeta.get('Readable')[i])
#                    ResourceRestIdx = ResourceMetaList.index(ResourcePath)
                    FilePathName.append(ResourceMeta.get('Path')[i])
                    FileURIName.append(ResourceMeta.get('URI')[i])
                else:
                    print 'ERROR: %s is not in %s.' % (ResourcePath, ResourceMeta.get('Path')[i])
    #                sys.exit()
        #=======================================================================
        # else grab just some files from RESOURCE...
        #=======================================================================
        else:
            for j in xrange(0, len(Files)):
                if (ResourceMeta.get('Path')[i].find(Files[j]) != -1) and (ResourceMeta.get('Path')[i] not in FilePathName):
                    FileName.append(ResourceMeta.get('Name')[i])
                    FilePathName.append(ResourceMeta.get('Path')[i])
                    FilePathNameReadable.append(ResourceMeta.get('Readable')[i])
                    
                    FileURIName.append(ResourceMeta.get('URI')[i])
                
#        CollectionsIdx = [i for i, x in enumerate(ResourceMeta.get('Format')) if x == FileType]

    
else:
    ScanMeta = getHCP.getScanMeta()
    CollectionsIdx = ScanMeta.get('Format').index(FileType)
    FilePathNameReadable.append(ScanMeta.get('Readable')[CollectionsIdx])
    
    currPath = ScanMeta.get('Path')[CollectionsIdx]
    currURI = ScanMeta.get('URI')[CollectionsIdx]
    
    if (OutputFile is not None):
        currPathBase = os.path.dirname(currPath)
        currURIBase = os.path.dirname(currURI)
        OutputFileName.append(OutputFile)
        OutputFilePathName.append(os.path.normpath(OutputDir +os.sep+ OutputFile))
        OutputFileURIName.append(os.path.normpath(OutputDir +os.sep+ OutputFile))
    else:
        OutputFileName.append(ScanMeta.get('Name')[CollectionsIdx])
        FileName.append(ScanMeta.get('Name')[CollectionsIdx])

    FilePathName.append(currPath)
    FileURIName.append(currURI)
    
                
if all(FilePathNameReadable):
    print ('cp %s %s' % (FilePathName, OutputFileName))
    writeHCP.writeFileFromPath(FilePathName, OutputFileName)
else:
    print ('Destination: %s  URI: %s' % (OutputDir, FileName))
    writeHCP.writeFileFromURL(getHCP, FileURIName, FileName)
    print 'Streamed bytes: ' + ', '.join(map(str, writeHCP.BytesStream))
    
print 'Written bytes: ' + ', '.join(map(str, writeHCP.BytesWrite))
print("Duration: %s" % (time.time() - sTime))
    
    
    
    
    
    