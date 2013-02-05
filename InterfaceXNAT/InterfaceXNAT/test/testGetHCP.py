'''
Created on Dec 23, 2012

@author: Tony
'''

#import glob
#import base64
# multiplatform stuff...
import os
import sys
import argparse
# Time manipulation...
import time
#import pytz
#from datetime import datetime
# XML parsing...
import xml.etree.ElementTree as ET
# Web stuff...
import socket
#import urllib
#import urllib2
#import urlnorm
#from ssl import SSLError
#from urllib2 import URLError, HTTPError

#sys.path.append('..' +os.sep) 
from pyHCP import getHCP

sTime = time.time()

#===============================================================================
# PARSE INPUT
#===============================================================================
parser = argparse.ArgumentParser(description="Program to figure out pipelines status via WORKFLOW XML...")
# input...
parser.add_argument("-U", "--User", dest="iUser", default='tony', type=str)
parser.add_argument("-P", "--Password", dest="iPassword", type=str)
parser.add_argument("-S", "--inputSubjects", dest="inputSubjects", default=None, help="pick subject, or a list of subjects")
# output...
parser.add_argument("-D", "--outputDir", dest="outputDir", type=str, help="where do you want to write output tab-text")
# timeout...
parser.add_argument("-t", "--time_out", dest="Timeout", type=float, default=16.0, help="change timeout")
# remote...
parser.add_argument("-Web", "--WS", dest="WebServer", type=str, default="https://intradb.humanconnectome.org", help="pick server")
parser.add_argument("-Proj", "--Project", dest="inputProject", type=str, default="HCP_Phase2", help="pick project")
# version...
parser.add_argument('--version', action='version', version='%(prog)s 0.1.1')

args = parser.parse_args()
User = args.iUser
Password = args.iPassword
inputSubjects = args.inputSubjects
    
outputDir = os.path.normpath(args.outputDir)

inputProject = args.inputProject
Server = args.WebServer
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
print "Running %s on %s" % (os.path.split(sys.argv[0])[1], socket.gethostname())
print sys.path
getHCP.Verbose = True
getHCP.DestinationDir = outputDir + os.sep

getHCP = getHCP(User, Password, Server)

#===============================================================================
# Setup output...
#===============================================================================
if (inputSubjects == None):
    inputSubjectsList = getHCP.Subjects
elif (inputSubjects != None):
    inputSubjectsList = inputSubjects.split(',')

outputFileBase = 'PipelineStatus'
if (len(inputSubjectsList) == 1):
    outputFileAppend = inputSubjectsList[0]
else:
    outputFileAppend = inputProject
    
outputFileExt = '.txt'
outputDirFile = '%s\\%s%s%s' % (outputDir, outputFileBase, outputFileAppend, outputFileExt)
#===============================================================================
# show usage of class...
#===============================================================================
if (showClassUsage):
    
    getHCP.Project = 'HCP_Phase2'
    getHCP.Subject = '103515'
    getHCP.Session = getHCP.Subject + '_strc'
    getHCP.SessionType = 'Structural'
    
    getHCP.Scan = '5'
    ScanIds, ScanTypes, ScanSeries, ScanQualty, ScanXnatId = getHCP.fGetSessionMeta()
    Names, Sizes, URIs, Collections = getHCP.fGetScanMeta()
    getHCP.fWriteFileFromURL(URIs, ['nii'])
    
    print getHCP.fGetProjects()
    print getHCP.SessionId
    print getHCP.Subjects
    print getHCP.fGetSubjectSessions()
    
    print getHCP.fGetFileInfo('https://intradb.humanconnectome.org/data/projects/HCP_Phase2/subjects/192439/experiments/192439_strc/resources/Details/files/StructuralHCP.log')
    print 'Bytes: ' +getHCP.fGetFileInfo('https://intradb.humanconnectome.org/data/projects/HCP_Phase2/subjects/197550/experiments/197550_diff/resources/Diffusion/files/Diffusion/data/bvals').get('Bytes')
    print getHCP.fGetProjects()
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

