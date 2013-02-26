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
parser.add_argument("-Session", "--Session", dest="Session", type=str, default=None, help="pick session")
parser.add_argument("-Scan", "--Scan", dest="Scan", type=str, default=None, help="pick scan")
parser.add_argument("-Resource", "--Resource", dest="Resource", type=str, default=None, help="pick resource")
parser.add_argument("-ResourceRoot", "--ResourceRoot", dest="ResourceRoot", type=str, default=None, help="pick resource root")
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
FileType = args.FileType
Server = args.Server
Resource = args.Resource
ResourceRoot = args.ResourceRoot
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
    
if (FileType is not None) and (FileType.find(',') != -1):
    FileType = FileType.split(',')
    
if (ResourceRoot is None):
    ResourceRoot = Resource
    
if not os.path.exists(OutputDir):
    os.makedirs(OutputDir)
#===============================================================================

print "Running %s on %s with IP %s" % (os.path.split(sys.argv[0])[1], socket.gethostname(), socket.gethostbyname(socket.gethostname()))

pyHCP = pyHCP(User, Password, Server)
getHCP = getHCP(pyHCP)
writeHCP = writeHCP(getHCP, OutputDir)

getHCP.Project = Project
getHCP.Subject = Subject
getHCP.Session = Session
getHCP.Scan = Scan

FilePathNameReadable = list()

if (Resource is not None):
    #===========================================================================
    # not complete, may finish later if management bitches about "resources"
    #===========================================================================
    FileName = None
    FilePathsNames = list()
    
    
    SubjectResources = getHCP.getSubjectResourcesMeta()
    if (Resource in SubjectResources.get('Label')):
        getHCP.Resource = Resource
        ResourceMeta = getHCP.getResourceMeta()
        for i in xrange(0, len(ResourceMeta.get('Path'))):
            ResourceMetaList = ResourceMeta.get('Path')[i].split('/')

            
            if (ResourceRoot in ResourceMetaList):
                FilePathNameReadable.append(ResourceMeta.get('Readable')[i])
                ResourceRestIdx = ResourceMetaList.index(ResourceRoot)
#                print '/'.join(ResourceMetaList[0:ResourceRestIdx]), '/'.join(ResourceMetaList[ResourceRestIdx:])
                if (ResourceMeta.get('Readable')[i]):
                    FilePathsNames.append(ResourceMeta.get('Path')[i])
                else:
                    FilePathsNames.append(ResourceMeta.get('URI')[i])
            else:
                print 'ERROR: %s is not in %s.' % (ResourceRoot, ResourceMeta.get('Path')[i])
                sys.exit()
                    
#        CollectionsIdx = [i for i, x in enumerate(ResourceMeta.get('Format')) if x == FileType]
    else:
        print 'ERROR: %s not in %s resources.' % (Resource, Subject)
        sys.exit()
    
else:
    ScanMeta = getHCP.getScanMeta()
    CollectionsIdx = ScanMeta.get('Format').index(FileType)
    FilePathNameReadable.append(ScanMeta.get('Readable')[CollectionsIdx])
    if ScanMeta.get('Readable')[CollectionsIdx]:
        FilePathsNames = ScanMeta.get('Path')[CollectionsIdx]
    else:
        FilePathsNames = ScanMeta.get('URI')[CollectionsIdx]
    
    #for i in xrange(0, len(ResourcePaths)): print 'Path: %s  \nURI: %s ' % (ResourcePaths[i], ResourceURIs[i])
    FileName = os.path.basename(FilePathsNames)
    FilePath = os.path.dirname(FilePathsNames)
    
    if (StripSession is not None):
        FileNameSplit = FileName.split('_')
        SessionIdx = FileNameSplit.index(StripSession)
        FileNameSans = FileNameSplit[0]
        for i in xrange(1, len(FileNameSplit)):
            if (i != SessionIdx):
                FileNameSans = FileNameSans + '_' + FileNameSplit[i]
        FileName = FileNameSans
        
    print FilePathsNames, FileName, FilePath
    
    if (OutputFile is not None):
        OutputDirFile = OutputDir + OutputFile
    else:
        OutputDirFile = OutputDir + FileName
                
if all(FilePathNameReadable):
    print ('cp %s %s' % (FilePathsNames, OutputDir))
    writeHCP.writeFileFromPath(FilePathsNames, FileName)
else:
    print ('Destination: %s  URI: %s' % (OutputDir, FilePathsNames))
    writeHCP.writeFileFromURL(getHCP, FilePathsNames, FileName)
    

print 'Streamed bytes: ' + ', '.join(map(str, writeHCP.BytesStream))
print 'Written bytes: ' + ', '.join(map(str, writeHCP.BytesWrite))
print("Duration: %s" % (time.time() - sTime))
    
    
    
    
    
    