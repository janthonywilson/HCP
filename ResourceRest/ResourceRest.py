'''
Created on 2012-05-17

@author: jwilso01
'''
import urllib
import urllib2
from urllib2 import URLError, HTTPError
from ssl import SSLError
import base64
import socket
import hashlib
import csv
import zipfile
import os
import sys
import time
import stat
import shutil
import argparse
import datetime
import subprocess

sTime = time.time()


#=NOTES:========================================================================
# -u tony -p ###### -W https://intradb.humanconnectome.org -P HCP_Phase2 -S 585862 -d MNINonLinear/Native -f ALL  -D C:\tmp\packets\MNINonLinear -x true
# Duration: 4993.046
#===============================================================================

#===============================================================================
# PARSE INPUT
#===============================================================================
parser = argparse.ArgumentParser(description="Alpha program to pull NIFTI data from XNAT and put it somewhere...")

parser.add_argument("-W", "--server", dest="restRoot", default="https://intradb.humanconnectome.org", type=str, help="specify which server to connect to")
parser.add_argument("-u", "--username", dest="restUser", type=str, help="username must be specified")
parser.add_argument("-p", "--password", dest="restPass", type=str, help="password must be specified")

parser.add_argument("-P", "--project", dest="inputProject", default="PARCEL_PILOT", type=str, help="specify project")
parser.add_argument("-S", "--subject", dest="inputSubject", default="CP10104", type=str, help="specify subject of interest")
parser.add_argument("-T", "--type", dest="inputDataType", default="structural", type=str, help="specify datatype")

parser.add_argument("-d", "--source_dir", dest="SourceDir", type=str, help="specify the XNAT source directory under the parent")
parser.add_argument("-f", "--source_file", dest="SourceFile", type=str, help="specify the XNAT source file")
parser.add_argument("-x", "--strip_file", dest="SourceStrip", type=str, default="false", help="strip session info from results files")
parser.add_argument("-t", "--time_out", dest="Timeout", type=float, default=30.0, help="change timeout")

parser.add_argument("-D", "--destination_dir", dest="destDir", default='tmp', type=str, help="specify the directory for output")
parser.add_argument("-M", "--print_csv", dest="printLists", default=False, help="print the lists to a csv file for looking at")
parser.add_argument("-V", "--Verbose", dest="Verbose", type=str, default='false', help="show more verbose output")

parser.add_argument('--version', action='version', version='%(prog)s 0.3.2')

args = parser.parse_args()

inputProject = args.inputProject
inputSubject = args.inputSubject
inputDataType = args.inputDataType.lower()
SourceDir = args.SourceDir
SourceFile = args.SourceFile
SourceStrip = args.SourceStrip
Verbose = str.lower(args.Verbose)

restUser = args.restUser
restPass = args.restPass
restRoot = args.restRoot
destDir = os.path.normpath(args.destDir)

printLists = args.printLists

DownloadRetryMax = 8
TimeoutMax = 1024.0
TimeoutStep = 8.0
Timeout = args.Timeout
TimeoutDefault = args.Timeout

restURI = restRoot
restRoot = args.restRoot + '/REST'

if (destDir[-1] != os.sep):
    destDir = destDir + os.sep


inputDir = SourceDir

# hack for using command line string input...
if (Verbose == 'false'): Verbose = False
elif (Verbose == 'true'): Verbose = True
#===============================================================================

print "Running on " + socket.gethostname()

#===============================================================================
# FUNCTIONS
#===============================================================================
def fReadURL( URL, SessionID, Timeout ):
    restRequest = urllib2.Request(URL)
    restRequest.add_header("Cookie", "JSESSIONID=" + SessionID);
    

    while (Timeout <= TimeoutMax):
        try:
            restConnHandle = urllib2.urlopen(restRequest, None, Timeout)
        except HTTPError, e:
            Timeout += TimeoutStep
            print 'HTTPError code: ' +str(e.code)+ '. Timeout increased to ' +str(Timeout)+' seconds for ' +URL
        except URLError, e:
            Timeout += TimeoutStep
            print 'URLError code: ' +str(e.reason)+ '. Timeout increased to ' +str(Timeout)+' seconds for ' +URL
        except SSLError, e:
            Timeout += TimeoutStep
            print 'SSLError code: ' +str(e.message)+ '. Timeout increased to ' +str(Timeout)+' seconds for ' +URL
        except socket.timeout:
            Timeout += TimeoutStep
            print 'Socket timed out. Timeout increased to ' +str(Timeout)+ ' seconds for ' +URL
            
        else:
            try:
#                socket.setdefaulttimeout(1.0)
                restResults = restConnHandle.read()
                return restResults
            except HTTPError, e:
                print 'HTTPError code: ' +str(e.code)+ '. File read timeout for ' +str(Timeout)+ ' seconds for ' +URL
                
    print 'ERROR: No reasonable timeout limit could be found for ' + URL
    sys.exit()
#===============================================================================
def fStripSession( inputName ):
    # check for session on input subject string...
    if (inputName.find("_strc") != -1) or (inputName.find("_diff") != -1) or (inputName.find("_fnc") != -1) or (inputName.find("_xtr") != -1):
        # strip out the session stuff.  Total hack with the index < stuff...
        sessionIdx = inputName.index("_")
        inputSubject = inputName[0:sessionIdx]
        try:
            fileIdx = inputName.index(".")
            try:
                underscoreIdx = inputName[sessionIdx+1:-1].index("_")
            except:
                underscoreIdx = float('inf')

            # ACK!  Hard coding...
            if (sessionIdx < underscoreIdx):
                outputName = inputSubject + inputName[fileIdx:]
            else:
                outputName = inputSubject + inputName[underscoreIdx+sessionIdx+1:]
        except:
            sessionIdxEnd = inputName[sessionIdx+1:].index("_")
            inputName = inputName[sessionIdxEnd+sessionIdx+2:]
            outputName = inputSubject +'_'+ inputName
        
    else:
        outputName = inputName
        
    return outputName
#===============================================================================
def fGetFileInfo( URL, SessionID, Timeout ):
    
    restRequest = urllib2.Request(URL)
    restRequest.add_header("Cookie", "JSESSIONID=" + SessionID);
     
    restRequest.get_method = lambda : 'HEAD'
    while (Timeout <= TimeoutMax):
            try:
                restConnHandle = urllib2.urlopen(restRequest, None)
                FileInfo = { 'ModDate': restConnHandle.info().getheader('Last-Modified'), 'Bytes': restConnHandle.info().getheader('Content-Length'), 'URL': URL }
            except HTTPError, e:
                if (e.code != 404):
                    Timeout += TimeoutStep
                    print 'HTTPError code: ' +str(e.code)+ '. Timeout increased to ' +str(Timeout)+' seconds for ' +URL
                else:
                    return '404 Error'
        
            if (FileInfo.get( 'Bytes' ) == None): 
                FileInfo[ 'Bytes' ] = '0'

            return FileInfo
#===============================================================================
# Get session ID...
#===============================================================================
restURL = restURI + '/data/JSESSION'
restPost = urllib.urlencode({'foo' : 'bar'})
restRequest = urllib2.Request(restURL, restPost)
restAuthHeader = "Basic %s" % base64.encodestring('%s:%s' % (restUser, restPass))[:-1]
restRequest.add_header("Authorization", restAuthHeader)
#Timeout = 0.1
while (Timeout <= TimeoutMax):
    try:
        restConnHandle = urllib2.urlopen(restRequest, None, Timeout)
#        print restConnHandle.code
        Timeout = TimeoutDefault
        break
    except URLError, e:
        Timeout += TimeoutStep
        print 'URLError code: ' +str(e.reason)+ '. Timeout increased to ' +str(Timeout)+' seconds for JSESSION cookie...'
        
        
restSessionID = restConnHandle.read()
#===============================================================================

restURL = restRoot + '/projects?format=csv'

restResults = fReadURL(restURL, restSessionID, Timeout)

restResultsSplit = restResults.split('\n')
restEndCount = restResults.count('\n')

ProjectMatch = False
for i in range(0,restEndCount):
    # check for match...
    if ProjectMatch: break

    currRow = restResultsSplit[i]
    currRowSplit = currRow.split(',')
    currRowCount = currRow.count(',')
    for j in range(0, currRowCount):
        currRowClean = currRowSplit[j].replace('"', '')
        if currRowClean == inputProject:
            ProjectMatch = True
            break


#===============================================================================
# PROJECTS
#===============================================================================
restURL = restRoot + '/projects/' + inputProject + '/subjects?format=csv'
if Verbose: print restURL

restResults = fReadURL(restURL, restSessionID, Timeout)

restResultsSplit = restResults.split('\n')
restEndCount = restResults.count('\n')
SubjectMatch = False
for i in range(0,restEndCount):
    # check for match...
    if SubjectMatch: break

    currRow = restResultsSplit[i]
    currRowSplit = currRow.split(',')
    currRowCount = currRow.count(',')
    for j in range(0, currRowCount):
        currRowClean = currRowSplit[j].replace('"', '')
        if currRowClean == inputSubject:
            SubjectMatch = True
            break

#===============================================================================
# SESSIONS
#===============================================================================
restURL = restRoot + '/projects/' + inputProject + '/subjects/' + inputSubject + '/experiments?format=csv'
if Verbose: print restURL

restResults = fReadURL(restURL, restSessionID, Timeout)

restResultsSplit = restResults.split('\n')
restEndCount = restResults.count('\n')
sessionList = list()
sessionListIdx = 0
for i in xrange(0,restEndCount):
    currRow = restResultsSplit[i]
    currRowSplit = currRow.split(',')
    currRowCount = currRow.count(',')
    for j in range(0, currRowCount):
        currRowClean = currRowSplit[j].replace('"', '')
        if i == 0:
            if currRowClean == 'label':
                labelIdx = j
                break
        elif j == labelIdx:
            sessionList.append(currRowClean)
            sessionListIdx += 1

#===============================================================================
# lets get the resource header info...
#===============================================================================
restURL = restRoot + '/projects/' + inputProject + '/subjects/' + inputSubject + '/experiments/' + sessionList[0] + '/resources?format=csv'
if Verbose: print restURL

restResults = fReadURL(restURL, restSessionID, Timeout)

restResultsSplit = restResults.split('\n')
restSessionHeader = restResultsSplit[0]
restSessionHeaderSplit = restSessionHeader.split(',')
restSessionHeaderCount = restSessionHeader.count(',')
restHeaderList = list()
for i in range(0, restSessionHeaderCount + 1):
    restHeaderList.append(restSessionHeaderSplit[i].replace('"', ''))


#==============================================================================
# relevant fields: label...
#==============================================================================
labelIdx = restHeaderList.index('label')
labelList = list()
fileNameList = list()
fileURIList = list()
fileSessionList = list()
fileSessionTypeList = list()
fileLabelList = list()
for i in xrange(0,len(sessionList)):
    currSession = sessionList[i]
    
    matchSession = False
    if (inputDataType == 'structural') and ( (currSession.find('strc') != -1) or (currSession.find('xtr') != -1) ):
        matchSession = True
    elif (inputDataType == 'functional') and ( (currSession.find('fnc') != -1) or (currSession.find('xtr') != -1) ):
        matchSession = True
    elif (inputDataType == 'diffusion') and ( (currSession.find('diff') != -1) or (currSession.find('xtr') != -1) ):
        matchSession = True
        
    if matchSession:
        restURL = restRoot + '/projects/' + inputProject + '/subjects/' + inputSubject + '/experiments/' + currSession + '/resources?format=csv'
    
        restResults = fReadURL(restURL, restSessionID, Timeout)
    
        restResultsSplit = restResults.split('\n')
        restEndCount = restResults.count('\n')
    
        #    print restResultsSplit[labelIdx]
        if (restEndCount > 1):
            for j in xrange(1,restEndCount):
    
                currRow = restResultsSplit[j]
                currRowSplit = currRow.split(',')
                currRowCount = currRow.count(',')
                currRowClean = currRowSplit[labelIdx].replace('"', '')
                labelList.append(currRowClean)
    
                restTypeURL = restRoot +'/projects/' + inputProject + '/subjects/' + inputSubject + '/experiments/' + currSession + '/scans?format=csv&columns=ID,type'
                currRestTypeResults = fReadURL(restTypeURL, restSessionID, Timeout)
                currRestTypeResultsSplit = currRestTypeResults.split('\n')
                currRestTypeEndCount = currRestTypeResults.count('\n')
                
                for k in xrange(1,currRestTypeEndCount):
                    currRowSplit = currRestTypeResultsSplit[k].split(',')
                    if ('"T1w"' in currRowSplit):
                        SubjectSessionType = 'structural'
                    elif ('"dMRI"' in currRowSplit):
                        SubjectSessionType = 'diffusion'
                    elif ('"tfMRI"' in currRowSplit):
                        SubjectSessionType = 'functional'
                
                restURL = restRoot +'/projects/'+ inputProject +'/subjects/'+ inputSubject +'/experiments/'+ currSession +'/resources/'+ currRowClean +'/files?format=csv'
                currRestResults = fReadURL(restURL, restSessionID, Timeout)
                currRestResultsSplit = currRestResults.split('\n')
                currRestEndCount = currRestResults.count('\n')
                
                for k in xrange(1,currRestEndCount):
                    newRow = currRestResultsSplit[k]
                    currRowSplit = newRow.split(',')
                    fileNameList.append(currRowSplit[0])
                    fileURIList.append(currRowSplit[2])
                    fileSessionList.append(currSession)
                    fileLabelList.append(currRowClean)
                    
                    fileSessionTypeList.append(SubjectSessionType)

#===============================================================================
# print lists if you want...
#===============================================================================
if printLists:
    headerStr = ['FileName', 'URI', 'Session']
    fileResultId = csv.writer(open('Results.txt', 'wb'), delimiter='\t')
    fileResultId.writerow(headerStr)

    for i in xrange(0, len(fileNameList)):
        fileResultId.writerow([fileNameList[i], fileURIList[i], fileSessionList[i]])

#===============================================================================
# make working directory...
#===============================================================================
#print  sys.platform
#    destDir = destDir + os.sep + inputSubject + os.sep + SeriesName
if not os.path.exists(destDir):
    os.makedirs(destDir)



#===============================================================================
# Now loop across SourceDir and SourceFile save and move...
#===============================================================================
fileSessionListUniq = set(fileSessionList)
if (len(fileSessionListUniq) > 1):
    for i in xrange(len(fileNameList)-1, -1, -1):
        currSession = fileSessionList[i]
        currSessionType = fileSessionTypeList[i]
        if (currSession.find('xtr') == -1) and (currSessionType != inputDataType):
            del fileNameList[i]
            del fileURIList[i]
            del fileSessionList[i]
            del fileLabelList[i]
    
    
restResultsTot = 0
for i in xrange(0,len(fileNameList)):
    currFileName = fileNameList[i].replace('"', '')
    currFileURI = fileURIList[i].replace('"', '')
    currSession = fileSessionList[i]
    currLabel = fileLabelList[i].replace('"', '')

    if (SourceDir == "mri") or (SourceDir == "surf") or (SourceDir == "mri/transforms"):
        SourceDir = 'T1w/' + inputSubject +'/'+ SourceDir

    Match = False
    if (SourceFile == "ALL"):
        if (currFileURI.find('/'+ SourceDir +'/'+ currFileName) != -1):

            if (inputDataType == 'structural') and (currLabel == 'Details'):
                if (currSession.find('strc') != -1) or (currSession.find('xtr') != -1):
                    restURL = restRoot +'/projects/'+ inputProject +'/subjects/'+ inputSubject +'/experiments/'+ currSession +'/resources/'+ currLabel +'/files/'+ SourceDir +'/'+ currFileName
                    if Verbose: print restURL
                    Match = True
            elif (inputDataType == 'functional') and (currLabel.find('BOLD') != -1):
                if (currSession.find('fnc') != -1) or (currSession.find('xtr') != -1):
                    restURL = restRoot +'/projects/'+ inputProject +'/subjects/'+ inputSubject +'/experiments/'+ currSession +'/resources/'+ currLabel +'/files/'+ SourceDir +'/'+ currFileName
                    if Verbose: print restURL
                    Match = True
            elif (inputDataType == 'diffusion') and (currLabel.find('Diffusion') != -1):
                if (currSession.find('diff') != -1) or (currSession.find('xtr') != -1):
                    restURL = restRoot +'/projects/'+ inputProject +'/subjects/'+ inputSubject +'/experiments/'+ currSession +'/resources/'+ currLabel +'/files/'+ SourceDir +'/'+ currFileName
                    if Verbose: print restURL
                    Match = True
                    
    else:
        if (currFileURI.find(SourceDir +'/'+ SourceFile) != -1):
            if (inputDataType == 'structural') and (currLabel == 'Details'):
                if (currSession.find('strc') != -1) or (currSession.find('xtr') != -1):
                    restURL = restURI + currFileURI
                    if Verbose: print restURL
                    Match = True
            elif (inputDataType == 'functional') and (currLabel.find('BOLD') != -1):
                if (currSession.find('fnc') != -1) or (currSession.find('xtr') != -1):
                    restURL = restURI + currFileURI
                    if Verbose: print restURL
                    Match = True
            elif (inputDataType == 'diffusion') and (currLabel == 'Diffusion'):
                if (currSession.find('diff') != -1) or (currSession.find('xtr') != -1):
                    restURL = restURI + currFileURI
                    if Verbose: print restURL
                    Match = True
                    

    if (Match == True):
        
        # rename the file...
        if (SourceStrip == 'true'):
            newFileName = fStripSession( currFileName )
        else:
            newFileName = currFileName

        if not os.path.isfile(destDir + newFileName):

            currFileInfo = fGetFileInfo(restURL, restSessionID, Timeout)
            restResults = fReadURL(restURL, restSessionID, Timeout)
            restResultsTot = restResultsTot + len(restResults)
            if (currFileInfo.get('Bytes') != str(len(restResults))):
                print 'WARNING: Expected ' +currFileInfo.get('Bytes')+ ' bytes and downloaded ' +str(len(restResults))+ ' bytes for file ' +newFileName
                
                # date format: 2012/10/08 10:35:27
                downloadRetries = 1
                while (currFileInfo.get('Bytes') != str(len(restResults))):
                    if (DownloadRetryMax >= downloadRetries):
                        print 'NOTE: Download retry ' +str(downloadRetries)+ ' (' +str(datetime.datetime.now())+ ')'
                        restResults = fReadURL(restURL, restSessionID, Timeout)
                        downloadRetries += 1
                    else:
                        print 'WARNING: Aborted download retries at ' + str(downloadRetries-1)
                        break

            with open(destDir + newFileName, 'wb') as outputFileId:
                writeCode = outputFileId.write(restResults)
                print 'File Name ' +newFileName+ ' and write code ' +str(writeCode)
                if sys.platform != 'win32':
                    subprocess.call('sync')
                outputFileId.close()
            outputFileSize = os.path.getsize(destDir + newFileName)
            if (currFileInfo.get('Bytes') != str(outputFileSize)):
                print 'WARNING: WROTE ' +str(len(restResults))+ ' bytes but expected ' +str(outputFileSize)+ ' bytes for file ' +newFileName

            print 'Dest dir and file: ' + destDir + newFileName


            if not os.path.isfile(destDir +os.sep+ newFileName):
                # move the file....
                if sys.platform != 'win32':
                    moveString = 'mv ' +newFileName+ ' ' +destDir
                    os.system(moveString)
                else:
                    shutil.move(newFileName, destDir)
                    if Verbose:
                        print 'Moving ' +newFileName+ ' to ' +destDir
                
                if os.path.isfile(os.getcwd() +os.sep+ newFileName):
                    try:
                        os.remove(os.getcwd() +os.sep+ newFileName)
                    except Exception, err:
                        if sys.platform == 'win32':
                            os.chmod(os.getcwd() + os.sep + newFileName, stat.S_IWRITE)
                        else:
                            os.chmod(os.getcwd() + os.sep + newFileName, 0766)
                        os.remove(os.getcwd() +os.sep+ newFileName)
                        print 'WARNING: Changed mod to delete ' + os.getcwd() + os.sep + newFileName + '.  This may or may not have been successful...'
                        pass
        else:
            # delete downloaded file...
            print "File already exists in destination directory..."           

    if (SourceFile != "ALL") and (Match == True):
        break


tTime = time.time() - sTime
#print outputNameList[niftiIdx]
print("Bytes written: %s" % str(restResultsTot))
print("Duration: %s" % tTime)
