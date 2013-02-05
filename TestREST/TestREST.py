'''
Created on 2012-05-17

@author: jwilso01
'''
import urllib2
import base64
import csv
import zipfile
import os
import sys
import shutil
#import shptree
#import optparse
import argparse


#===============================================================================
# PARSE INPUT
#===============================================================================
parser = argparse.ArgumentParser(description="Test program to pull NIFTI data from XNAT and put it somewhere...")

parser.add_argument("-P", "--project", dest="inputProject", default="PARCEL_PILOT", type=str, help="specify project")
parser.add_argument("-S", "--subject", dest="inputSubject", default="CP10104", type=str, help="specify subject of interest")
parser.add_argument("-s", "--series", dest="SeriesName", default='BOLD_MOTOR1', type=str, help="specify the series name of interest")
parser.add_argument("-t", "--series_type", dest="TypeStr", default='tfMRI', type=str, help="specify the series type")
parser.add_argument("-u", "--username", dest="restUser", type=str, help="username must be specified")
parser.add_argument("-p", "--password", dest="restPass", type=str, help="password must be specified")
parser.add_argument("-D", "--destination_dir", dest="destDir", default='/tmp', type=str, help="specify the directory for output")
parser.add_argument('--version', action='version', version='%(prog)s 0.1')

args = parser.parse_args()

inputProject = args.inputProject
if inputProject == 'PARCEL_PILOT':
    altProject = 'Phase1Reg'
else: 
    altProject = ''
inputSubject = args.inputSubject
SeriesName = args.SeriesName
TypeStr = args.TypeStr
restUser = args.restUser
restPass = args.restPass
destDir = args.destDir
#===============================================================================


#===============================================================================
# FUNCTIONS
#===============================================================================
def fReadURL( URL, User, Password ):
    restRequest = urllib2.Request(restURL)

    restAuthHeader = "Basic %s" % base64.encodestring('%s:%s' % (restUser, restPass))[:-1]
    restRequest.add_header("Authorization", restAuthHeader)

    restConnHandle = urllib2.urlopen(restRequest)
    restResults = restConnHandle.read()
    return restResults
#===============================================================================

printLists = False

#out = "<html>%s%s%s%s</html>" % (head, prologue, query, tail)
#restRoot = 'https://hcpi-dev-cuda00.nrg.mir/data/archive/'
restRoot = 'https://hcpi-dev-cuda00.nrg.mir/REST/'
restURL = restRoot + 'projects?format=csv'

restResults = fReadURL(restURL, restUser, restPass)

#restEndIdx = restResults.find('\n',0,-1)

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
        
        
# PROJECTS
restURL = restRoot + 'projects/' + inputProject + '/subjects?format=csv'
print restURL

restResults = fReadURL(restURL, restUser, restPass)

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

# SESSIONS
restURL = restRoot + 'projects/' + inputProject + '/subjects/' + inputSubject + '/experiments?format=csv'
print restURL

restResults = fReadURL(restURL, restUser, restPass)

restResultsSplit = restResults.split('\n')
restEndCount = restResults.count('\n')
sessionList = list()
sessionListIdx = 0
for i in range(0,restEndCount):
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
# lets get the header info...
#===============================================================================
restURL = restRoot + 'projects/' + inputProject + '/subjects/' + inputSubject + '/experiments/' + sessionList[0] + '/scans?format=csv'
print restURL

restResults = fReadURL(restURL, restUser, restPass)

restResultsSplit = restResults.split('\n')
restSessionHeader = restResultsSplit[0]
restSessionHeaderSplit = restSessionHeader.split(',')
restSessionHeaderCount = restSessionHeader.count(',')
restHeaderList = list()
for i in range(0, restSessionHeaderCount + 1):
    restHeaderList.append(restSessionHeaderSplit[i].replace('"', ''))
    

#==============================================================================
# relevant fields: session, ID, series_description...
#==============================================================================
idIdx = restHeaderList.index('ID')
seriesIdx = restHeaderList.index('series_description')
typeIdx = restHeaderList.index('type')

sessionRepList = list()
idList = list()
typeList = list()
seriesList = list()
visitList = list()

#===============================================================================
# Now loop across sessions with header...
#===============================================================================
for i in range(0,len(sessionList)):
    currSession = sessionList[i]
    restURL = restRoot + 'projects/' + inputProject + '/subjects/' + inputSubject + '/experiments/' + currSession + '/scans?format=csv'

    restResults = fReadURL(restURL, restUser, restPass)
    
    restResultsSplit = restResults.split('\n')
    restEndCount = restResults.count('\n')
    
    for j in range(1,restEndCount):
        currRow = restResultsSplit[j]
        currRowSplit = currRow.split(',')
        currRowCount = currRow.count(',')
        
        sessionRepList.append(currSession)
        idList.append(currRowSplit[idIdx].replace('"', ''))
        typeList.append(currRowSplit[typeIdx].replace('"', ''))
        seriesList.append(currRowSplit[seriesIdx].replace('"', ''))
        currSessionSplit = currSession.split('_')
        visitList.append(currSessionSplit[1])



if printLists:
    headerStr = ['Session', 'ID', 'Type', 'Series']
    fileResultId = csv.writer(open('Results.txt', 'wb'), delimiter='\t')
    fileResultId.writerow(headerStr)
    
    for i in range(0, len(idList)):
        fileResultId.writerow([sessionRepList[i], idList[i], typeList[i], seriesList[i]])
    

SeriesNameListIdx = list()
for i in range(0, len(seriesList)):
    if seriesList[i] == SeriesName and typeList[i] == TypeStr:
        SeriesNameListIdx.append(i)
        restURL = restRoot + 'projects/' + inputProject + '/subjects/' + inputSubject + '/experiments/' + sessionRepList[i] + '/scans/' + idList[i] + '/resources/NIFTI/files?format=zip'
        
        print restURL

#===============================================================================
# Now get the file and save it, extract, move...   
#===============================================================================
restResults = fReadURL(restURL, restUser, restPass)
    
# save data to disk
outputFilename = inputSubject + '_' + SeriesName + '.zip'
print "Saving to " + os.getcwd() + os.sep + outputFilename
outputFileId = open(outputFilename, 'wb')
outputFileId.write(restResults)
outputFileId.close()

outputZipObj = zipfile.ZipFile(outputFilename)
outputZipInfo = outputZipObj.infolist()
outputNameList = outputZipObj.namelist()
#outputZipObj.printdir()

#===============================================================================
# find nifti in list...may or may not use this...
#===============================================================================
niftiIdx = -1
for i in range(0, len(outputNameList)):
    currName = outputNameList[i]
    if currName.find('nii.gz') != -1:
        niftiIdx = i
    
#===============================================================================
# make working directory...
#===============================================================================
#print  sys.platform
if os.name == 'nt':
#    destDir = 'R:\\tmp\\' + SeriesName
    destDir = 'R:\\tmp\\' + inputSubject + '\\' + SeriesName
else:
    destDir = '/home/NRG/jwilso01/tmp/' + inputSubject + '/' + SeriesName
    
if not os.path.exists(destDir):
    os.makedirs(destDir)

#===============================================================================
# put contents into the destDir...
#===============================================================================
print "Unpacking files to " + destDir + os.sep
for member in outputZipObj.namelist():
    currFilename = os.path.basename(member)

    # copy file (taken from zipfile's extract)
    source = outputZipObj.open(member)
    target = file(os.path.join(destDir, currFilename), "wb")
    shutil.copyfileobj(source, target)
    source.close()
    target.close()

outputZipObj.close()


#print outputNameList[niftiIdx]
print 'Bytes written: ' + str(len(restResults))

     



    


    
    
    
    
    
    
    
