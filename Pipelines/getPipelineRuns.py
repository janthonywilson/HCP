'''
Created on Sep 16, 2012

@author: Tony
'''

import urllib2
import urllib
import base64
import socket
import os
import sys
import time
import argparse
from ssl import SSLError
from urllib2 import URLError, HTTPError


sTime = time.time()
print "Running on " + socket.gethostname()

#===============================================================================
# PARSE INPUT
#===============================================================================
parser = argparse.ArgumentParser(description="Alpha program to figure out what pipelines have been run by their logs in XNAT...")

parser.add_argument("-W", "--server", dest="Server", default="https://intradb.humanconnectome.org", type=str, help="specify which server to connect to")
parser.add_argument("-u", "--username", dest="User", type=str, help="username must be specified")
parser.add_argument("-p", "--password", dest="Password", type=str, help="password must be specified")
parser.add_argument("-s", "--subject", dest="Subject", default="ALL", type=str, help="specify subject of interest")
parser.add_argument("-P", "--project", dest="Project", default="HCP_Phase2", type=str, help="specify project")


parser.add_argument("-D", "--destination_dir", dest="destDir", default='tmp', type=str, help="specify the directory for output")
parser.add_argument("-M", "--print_csv", dest="printLists", default=False, help="print the lists to a csv file for looking at")

parser.add_argument('--version', action='version', version='%(prog)s 0.1')

args = parser.parse_args()

Project = args.Project
Subject = args.Subject

User = args.User
Password = args.Password
Server = args.Server
destDir = args.destDir

printLists = args.printLists

if (Server.find('http') == -1):
    Server = 'https://' + Server

restURI = Server
restRoot = Server + '/REST'

#===============================================================================
# FUNCTIONS
#===============================================================================
def fPrintData( inputSubject, inputType, inputPPL, inputTime, inputEpoch, inputStem, outputDir ):
    
    inputListLen = len(inputSubject) 
    headerStr = ['SubjectID', 'DataType', 'Pipeline', 'Time [RFC822]', 'EpochTime']
    fileName = inputStem + '_PipelineStatus.txt'

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    fileID = open(outputDir +os.sep+ fileName, 'wb')
    fileID.write(headerStr[0]+'\t')
    fileID.write(headerStr[1]+'\t')
    fileID.write(headerStr[2]+'\t')
    fileID.write(headerStr[3]+'\t')
    fileID.write(headerStr[4]+'\n')
    
    for i in xrange(0, inputListLen):
            fileID.write('%s' % inputSubject[i] + "\t")
            fileID.write('%s' % inputType[i] + "\t")
            fileID.write('%s' % inputPPL[i] + "\t")
            fileID.write('%s' % inputTime[i] + "\t")
            fileID.write('%s' % inputEpoch[i] + "\n")
#===============================================================================
def fReadURL( URL, SessionID, ModDate ):
    
    restRequest = urllib2.Request(URL)
    restRequest.add_header("Cookie", "JSESSIONID=" + SessionID);
    
    if ModDate:    
        restRequest.get_method = lambda : 'HEAD'
        try:
            restConnHandle = urllib2.urlopen(restRequest)
        except HTTPError, e:
            print 'HTTPError code: ' +str(e.code)+ ' for ' +URL
        restResults = restConnHandle.info().getheader('Last-Modified')
#        restResultsSplit = restResults.split(',')
#        restResults = restResultsSplit[1]
    else:

        try:
            restConnHandle = urllib2.urlopen(restRequest)
        except HTTPError, e:
            print 'HTTPError code: ' +str(e.code)+ ' for ' +URL
            sys.exit()
        except URLError, e:
            print 'URLError code: ' +str(e.reason)+ ' for ' +URL
            sys.exit()
        except SSLError, e:
            print 'SSLError code: ' +str(e.message)+ ' for ' +URL
            sys.exit()
        except socket.timeout:
            print 'Socket timed out for ' +URL
            sys.exit()
        else:
            restResults = restConnHandle.read()
    return restResults
#===============================================================================
# Get session ID...
#===============================================================================
restURL = restURI + '/data/JSESSION'
restPost = urllib.urlencode({'foo' : 'bar'})
restRequest = urllib2.Request(restURL, restPost)
restAuthHeader = "Basic %s" % base64.encodestring('%s:%s' % (User, Password))[:-1]
restRequest.add_header("Authorization", restAuthHeader)
Timeout = 15.0
try:
    restConnHandle = urllib2.urlopen(restRequest, None, Timeout)
except URLError, e:
    print 'URLError code: ' +str(e.reason)+ '. Timeout was ' +str(Timeout)+' seconds for JSESSION cookie...'
    
try:        
    restSessionID = restConnHandle.read()
except HTTPError, e:
    print 'HTTPError code: ' +str(e.code)+ '. File read timeout for ' +str(Timeout)+ ' seconds for ' +restURL
    
#===============================================================================
# Get subject list if input empty...
#===============================================================================
idList = list()
subjectList = list()
    
restURL = restRoot + '/projects/' + Project + '/subjects?format=csv'
#print restURL
#    "ID","project","label","insert_date","insert_user","URI"
restResults = fReadURL(restURL, restSessionID, False)

restResultsSplit = restResults.split('\n')
restEndCount = restResults.count('\n')
restSessionHeader = restResultsSplit[0]
restSessionHeaderSplit = restSessionHeader.split(',')

idIdx = restSessionHeaderSplit.index('"ID"')
subjectIdx = restSessionHeaderSplit.index('"label"')

for i in xrange(1, restEndCount):
    currRow = restResultsSplit[i]
    
    currRowSplit = currRow.split(',')

    idList.append(currRowSplit[idIdx].replace('"', ''))
    subjectList.append(currRowSplit[subjectIdx].replace('"', ''),)
    
if (Subject != 'ALL'):
    singleSubjectIdx = subjectList.index(Subject)
    subjectList = subjectList[singleSubjectIdx]
    idList = idList[singleSubjectIdx]

#===============================================================================
# Loop over subjects to get details...
#===============================================================================
sessionList = list()
subjectFullList = list()
pipelineList = list()
logNameList = list() 
FileModTimeList = list()
timeEpochList = list()
totalRestCalls = 0
if isinstance(subjectList, str): loopLen = 1
else: loopLen = len(subjectList)
for i in xrange(0, loopLen):
    
    if isinstance(subjectList, str): 
        currSubject = subjectList
    else: currSubject = subjectList[i]
    
    restURL = restRoot + '/projects/' + Project + '/subjects/' + currSubject + '/experiments?format=csv'
    
    restResults = fReadURL(restURL, restSessionID, False)
    totalRestCalls += 1
    
    restResultsSplit = restResults.split('\n')
    restEndCount = restResults.count('\n')
    restResultsHeader = restResultsSplit[0]
    restResultsHeaderSplit = restResultsHeader.split(',')
    restLabelIdx = restResultsHeaderSplit.index('"label"')
    
    for j in xrange(1,restEndCount):
        currRow = restResultsSplit[j]
        if (len(currRow) == 0): 
            j =+ 1 
            break
        currRowSplit = currRow.split(',')
        currRowCount = currRow.count(',')
        currLabel = currRowSplit[restLabelIdx].replace('"', '')
        
        if (currLabel.find("_strc") != -1) or (currLabel.find("_diff") != -1) or (currLabel.find("_fnc") != -1) or (currLabel.find("_xtr") != -1):
            
            restURL = restRoot +'/projects/'+ Project +'/subjects/'+ currSubject +'/experiments/'+ currLabel + '/resources?format=csv'
            
            restResults = fReadURL(restURL, restSessionID, False)
            totalRestCalls += 1
            
            
            restSessionResultsSplit = restResults.split('\n')
            restSessionCount = restResults.count('\n')
            restSessionHeader = restSessionResultsSplit[0]
            restSessionHeaderSplit = restSessionHeader.split(',')
            restExpLabelIdx = restSessionHeaderSplit.index('"label"')

            for k in xrange(1,restSessionCount):
    
                currRow = restSessionResultsSplit[k]
                currRowSplit = currRow.split(',')
                currRowCount = currRow.count(',')
                currRowClean = currRowSplit[restExpLabelIdx].replace('"', '')

                restURL = restRoot +'/projects/'+ Project +'/subjects/'+ currSubject +'/experiments/'+ currLabel +'/resources/'+ currRowClean +'/files?format=csv'
                
                currRestResults = fReadURL(restURL, restSessionID, False)
                totalRestCalls += 1

                currRestResultsSplit = currRestResults.split('\n')
                restFileCount = currRestResults.count('\n')
                restResultsHeader = currRestResultsSplit[0]
                restResultsHeaderSplit = restResultsHeader.split(',')
                restFileNameIdx = restResultsHeaderSplit.index('"Name"')
                
                for m in xrange(1, restFileCount):
                    
                    currRow = currRestResultsSplit[m]
                    currRowSplit = currRow.split(',')
                    currRowName = (currRowSplit[restFileNameIdx].replace('"', ''))
                    
                    if (currRowName.find('FunctionalHCP.log') != -1) or (currRowName.find('StructuralHCP.log') != -1) or (currRowName.find('StructHCP.log') != -1) or \
                     (currRowName.find('DiffusionHCP.log') != -1):

                        restURL = restRoot +'/projects/'+ Project +'/subjects/'+ currSubject +'/experiments/'+ currLabel +'/resources/'+ currRowClean +'/files/'+ currRowName
                        
                        currFileModTime = fReadURL(restURL, restSessionID, True)
                        totalRestCalls += 1
                        
                        timeStruct = time.strptime(currFileModTime, "%a, %d %b %Y %H:%M:%S GMT")
                        timeEpoch = time.mktime(timeStruct)
                        logNameNoExtension, logNameExtension = os.path.splitext(currRowName)
                        if (currRowClean == 'Details'): currPipelineName = 'STRUCTURAL'
                        else: currPipelineName = currRowClean
                        
                        if not printLists:
                            print currSubject, currPipelineName, logNameNoExtension, currFileModTime, timeEpoch
                        
                        sessionList.append(currSubject)
#                        subjectFullList.append(currSubject)
                        pipelineList.append(currPipelineName)
                        logNameList.append(logNameNoExtension) 
                        FileModTimeList.append(currFileModTime)
                        timeEpochList.append(timeEpoch)
        

    
    
if printLists:
        fPrintData( sessionList, pipelineList, logNameList, FileModTimeList, timeEpochList, Subject, destDir )

    
    
print "Rest Calls: " + str(totalRestCalls)
print("Duration: %s" % ( time.time() - sTime ))




