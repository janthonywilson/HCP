'''
Created on Dec 23, 2012

@author: Tony
'''

#import glob
#import base64
# multiplatform stuff...
import os
import sys
import socket
import argparse
# Time manipulation...
import time
# XML parsing...
import xml.etree.ElementTree as ET

#sys.path.append('..' +os.sep) 
from pyHCP import pyHCP, getHCP, writeHCP

sTime = time.time()

#-U hodgem3 -P Michael333



#===============================================================================
# PARSE INPUT
# -D C:\tmp\QC -U tony -P passfoo -W hcpx-demo.humanconnectome.org
#===============================================================================
parser = argparse.ArgumentParser(description="Program to figure out pipelines status via WORKFLOW XML...")
# input...
parser.add_argument("-U", "--User", dest="User", default='tony', type=str)
parser.add_argument("-P", "--Password", dest="Password", type=str)
parser.add_argument("-S", "--Subjects", dest="Subjects", default=None, help="pick subject, or a list of subjects")
# output...
parser.add_argument("-D", "--outputDir", dest="outputDir", type=str, help="where do you want to write output tab-text")
# timeout...
parser.add_argument("-t", "--time_out", dest="Timeout", type=float, default=16.0, help="change timeout")
# remote...
parser.add_argument("-W", "--Web", dest="Server", type=str, default="https://intradb.humanconnectome.org", help="pick server")
parser.add_argument("-Proj", "--Project", dest="Project", type=str, default="HCP_Phase2", help="pick project")
# version...
parser.add_argument('--version', action='version', version='%(prog)s 0.1.1')

args = parser.parse_args()
User = args.User
Password = args.Password
Subjects = args.Subjects
    
outputDir = os.path.normpath(args.outputDir)

Project = args.Project
Server = args.Server
Server.strip()
if (Server[-1] != '/'):
    Server = Server + '/'
if (Server.find('http') == -1):
    Server = 'https://' + Server
    
showClassUsage = True
printLists = True
fromFileXML = False
Verbose = True
TimeoutStep = 8
TimeoutMax = 1024
Timeout = args.Timeout
TimeoutDefault = args.Timeout
PipelineProcList = ['FunctionalHCP', 'StructuralHCP', 'DiffusionHCP']

#===============================================================================
# FUNCTIONS
#===============================================================================
def fPrint( outputDirFile, headerStr, *args ):

    outputFile = os.path.basename(outputDirFile)
    ouputFileBase, outputFileExt = os.path.splitext(outputFile)
    outputDir = os.path.dirname(os.path.normpath(outputDirFile)) + os.sep

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
        
    with open(outputDirFile, 'wb') as fileID:
        
        for i in xrange(0, len(headerStr)):
            if (i < len(headerStr)-1):
                fileID.write(headerStr[i]+'\t')
            else:
                fileID.write(headerStr[i]+'\n')
                
        for i in xrange(0, len(args[0])):
            for j in xrange(0, len(args)):
                
                if (j < len(args)-1):
                    fileID.write('%s' % args[j][i] + "\t")
                else:
                    fileID.write('%s' % args[j][i] + "\n")
#===============================================================================

ET.register_namespace('wrk', 'http://nrg.wustl.edu/workflow')
ET.register_namespace('pip', 'http://nrg.wustl.edu/pipeline')

#===============================================================================
# Init lists
# headerStr = ['SubjectID', 'DataType', 'Pipeline', 'Time [RFC822]', 'EpochTime']
# New Lists
# headerStr = ['SubjectID', 'DataType', 'Pipeline', 'Status', 'PercentComplete', 'Session', 'TimeLaunch', 'TimeCompletion', 'TimeLaunchEpoch']
#===============================================================================
AltSubjectList = list()
CompleteTimeList = list()
DataTypeList = list()
FunctionalSeriesList = list()
JobIdList = list()
LaunchEpochTimeList = list()
LaunchReadTimeList = list()
LaunchTimeList = list()
PercentCompleteList = list()
PipelineList = list()
SessionList = list()
StatusList = list()
SubjectList = list()
#===============================================================================
# init interface to server and get subjects if none input...
#===============================================================================
print "Running %s on %s with IP %s" % (os.path.split(sys.argv[0])[1], socket.gethostname(), socket.gethostbyname(socket.gethostname()))
print sys.path[0]

pyHCP = pyHCP(User, Password, Server)
getHCP = getHCP(pyHCP)
writeHCP = writeHCP(getHCP, outputDir)
getHCP.Verbose = True
if (getHCP.Server.find('intradb') != -1): 
    getHCP.Scan = '2'
    getHCP.Project = 'HCP_Phase2'
else:
    getHCP.Scan = '105'
    getHCP.Project = 'HCP_Q1'

DestinationDir = outputDir + os.sep
writeHCP = writeHCP(getHCP, DestinationDir)

#===============================================================================
# Setup subjects...
#===============================================================================
if (Subjects == None):
    getHCP.Subjects = getHCP.getSubjects()
elif (Subjects != None):
    getHCP.Subjects = Subjects.split(',')

    
#outputFileExt = '.txt'
#outputDirFile = '%s\\%s%s%s' % (outputDir, outputFileBase, outputFileAppend, outputFileExt)
#===============================================================================
# show usage of class...
#===============================================================================
if (showClassUsage):
    
#==============================================================================
# getAssessorIDs
# getAssessorOutputFile
# getFileInfo
# getProjects
# getScanMeta
# getScanParms
# getServer
# getSessionId
# getSessionMeta
# getSessionQuality
# getSubjectResourceMeta
# getSubjectResourcesMeta
# getSubjects
# getSubjectSessions
# getSubjectsSessions
# getURLString
#==============================================================================

    Projects = getHCP.getProjects()
    getHCP.Subject = getHCP.Subjects[0]
    getHCP.SubjectSessions = getHCP.getSubjectSessions()
    getHCP.Session = getHCP.SubjectSessions.get('Sessions')[0]
    getHCP.SessionMeta = getHCP.getSessionMeta()
    

    ScanParms = getHCP.getScanParms()
    ScanMeta = getHCP.getScanMeta()
    
    print ScanMeta.get('URIs')[ScanMeta.get('Collections').index('NIFTI')]
    writeHCP.writeFileFromURL(getHCP, ScanMeta.get('URIs')[ScanMeta.get('Collections').index('NIFTI')])
    
    if (getHCP.Server.find('intradb') != -1):
        AssessorIDs = getHCP.getAssessorIDs( )
        AssessorFileURIList = getHCP.getAssessorOutputFile( AssessorIDs )
        
        
    print getHCP.getProjects()
    print getHCP.SessionId, getHCP.getSessionId()
    print getHCP.Subjects
    print getHCP.getSubjectSessions()
    
#    print getHCP.getFileInfo('https://intradb.humanconnectome.org/data/projects/HCP_Phase2/subjects/192439/experiments/192439_strc/resources/Details/files/StructuralHCP.log')
#    print 'Bytes: ' +getHCP.getFileInfo('https://intradb.humanconnectome.org/data/projects/HCP_Phase2/subjects/197550/experiments/197550_diff/resources/Diffusion/files/Diffusion/data/bvals').get('Bytes')
    print getHCP.getProjects()
    print getHCP.Session
    
#    Quality, ScanId, Series, Sessions, ScanType = getHCP.fGetSessionQuality()
#    print Quality
#    print ScanType[2]
#    print getHCP.fGetScanParms(ScanId[2])

#    Projects = getHCP.fGetProjects()
#    Sessions, Quality, Series, ScanIds, ScanType = getHCP.fGetSessionQuality()
#    fPrint( outputDir + 'Quality.txt', ['Sessions', 'Quality', 'Series', 'ScanIds', 'ScanType'], Sessions, Quality, Series, ScanIds, ScanType)

#    getHCP.AssessorDataType = 'qcAssessmentData'
#    AssessorIDs = getHCP.fGetAssessorIDs( )
#    AssessorFileURIList = getHCP.fGetAssessorOutputFile( AssessorIDs )
#    getHCP.fWriteFileFromURL(AssessorFileURIList, ('Fourier', 'txt'))



#===============================================================================



print("Duration: %s" % (time.time() - sTime))

