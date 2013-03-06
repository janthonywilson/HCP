'''
Created on Dec 24, 2012

@author: Tony
'''


# multiplatform stuff...
import os
import sys
import argparse
# Time manipulation...
import time
# Web stuff...
import socket

from pyHCP import pyHCP, getHCP, writeHCP




sTime = time.time()

#===============================================================================
# PARSE INPUT
#===============================================================================
parser = argparse.ArgumentParser(description="Get outputs of structural QC using getHCP...")
# input...
parser.add_argument("-U", "--User", dest="iUser", default='tony', type=str)
parser.add_argument("-P", "--Password", dest="iPassword", type=str)
parser.add_argument("-S", "--inputSubjects", dest="inputSubjects", default=None, help="pick subject, or a list of subjects")
# output...
parser.add_argument("-D", "--outputDir", dest="outputDir", type=str, help="where do you want to write output tab-text")
# timeout...
parser.add_argument("-t", "--time_out", dest="Timeout", type=float, default=16.0, help="change timeout")
# remote...
parser.add_argument("-iWeb", "--WS", dest="WebServer", type=str, default="https://intradb.humanconnectome.org", help="pick server")
parser.add_argument("-iProj", "--Project", dest="inputProject", type=str, default="HCP_Phase2", help="pick project")
# version...
parser.add_argument('--version', action='version', version='%(prog)s 0.1.2')

args = parser.parse_args()
User = args.iUser
Password = args.iPassword
inputSubjects = args.inputSubjects
    
OutputDir = os.path.normpath(args.outputDir)

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

#===============================================================================
# FUNCTIONS
#===============================================================================
def fPrintTabList( outputDirFile, headerStr, *args ):

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

#===============================================================================
# init interface to server and get subjects if none input...
#===============================================================================
print "Running %s on %s" % (os.path.split(sys.argv[0])[1], socket.gethostname())

pyHCP = pyHCP(User, Password, Server)
getHCP = getHCP(pyHCP)
writeHCP = writeHCP(getHCP, OutputDir)

getHCP.Verbose = True
getHCP.Project = 'HCP_Phase2'



#===============================================================================
# Setup output...
#===============================================================================
if (inputSubjects == None):
    inputSubjectsList = getHCP.getSubjects()
elif (inputSubjects != None):
    inputSubjectsList = inputSubjects.split(',')


if (len(inputSubjectsList) == 1):
    outputFileAppend = inputSubjectsList[0]
else:
    outputFileAppend = inputProject
    
outputFileExt = '.txt'

#===============================================================================
# grab output files...
# NOTES: need to deal with xtr QC...
#===============================================================================

for i in xrange(0, len(inputSubjectsList)):
    
    getHCP.Subject = inputSubjectsList[i]
    print i, getHCP.Subject
    getHCP.Session = getHCP.Subject + '_strc'
    getHCP.SessionType = 'Structural'
#    getHCP.AssessorDataType = 'qcAssessmentData'
    AssessorIDs = getHCP.getAssessorIDs( )

    currSessions, currSessionType = getHCP.getSubjectSessions()
    StrcIdx = [k for k, element in enumerate(currSessionType) if element == 'strc']

    for j in xrange(0, len(StrcIdx)):
        writeScanQuality = list()
        writeScanSeries = list()
        writeAssessorFileURIList = list()
        getHCP.Session = currSessions[StrcIdx[j]]
        
        # get quality and scan type to match up with QC output...
        Quality, ScanId, Series, Sessions, ScanType = getHCP.getSessionQuality()
        
        #===============================================================================
        # make a list of quality and series for printing...
        ScanTypeIdx = [k for k, element in enumerate(ScanType) if (element == 'T1w') or (element == 'T2w')]
        for k in xrange(0, len(ScanTypeIdx)):
            writeScanQuality.append(Quality[ScanTypeIdx[k]])
            writeScanSeries.append(Series[ScanTypeIdx[k]])
        outputFileAppend = 'Quality'
        outputDirPrepend = 'Quality'
        outputFileBase = getHCP.Session
        outputDirFile = '%s\\%s\\%s%s' % (outputDir, outputDirPrepend, outputFileBase, outputFileExt)
        fPrintTabList( outputDirFile, ['Quality', 'Series'], writeScanQuality, writeScanSeries)
        #===============================================================================
        
        AssessorIDs = getHCP.getAssessorIDs( )
        AssessorFileURIList = getHCP.getAssessorOutputFile( AssessorIDs )
        
        for k in xrange(0, len(AssessorFileURIList)):
            currURI = AssessorFileURIList[k]
            if ( (currURI.find('Fourier') != -1) and (currURI.find('txt') != -1) and not (currURI.find('BOLD') != -1) ):
                writeAssessorFileURIList.append(currURI)
                
        getHCP.writeFileFromURL(writeAssessorFileURIList)



#===============================================================================



print("Duration: %s" % (time.time() - sTime))
