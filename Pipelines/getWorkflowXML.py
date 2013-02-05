'''
Created on 2012-12-11

@author: jwilso01
'''


import glob
import base64
# multiplatform stuff...
import os
import sys
import argparse
# Time manipulation...
import time
import pytz
from datetime import datetime
# XML parsing...
import xml.etree.ElementTree as ET
# Web stuff...
import socket
import urllib
import urllib2
from ssl import SSLError
from urllib2 import URLError, HTTPError

#from pytz import timezone
#from dateutil.parser import parse
#from dateutil import tz

time



sTime = time.time()
#===============================================================================
# -iD C:\Users\jwilso01\workspace\data\xml\workflow\ -iP asdf -oD C:\tmp\work -W http://hcpi-dev-cuda00.nrg.mir/
#===============================================================================

#===============================================================================
# PARSE INPUT
#===============================================================================
parser = argparse.ArgumentParser(description="Program to figure out pipelines status via WORKFLOW XML...")
# input...
parser.add_argument("-iU", "--User", dest="iUser", default='tony', type=str)
parser.add_argument("-iP", "--Password", dest="iPassword", type=str)
parser.add_argument("-iS", "--inputSubjects", dest="inputSubjects", default=None, help="pick subject, or a list of subjects")
# output...
parser.add_argument("-oD", "--outputDir", dest="outputDir", type=str, help="where do you want to write output tab-text")
# timeout...
parser.add_argument("-t", "--time_out", dest="Timeout", type=float, default=16.0, help="change timeout")
# remote...
parser.add_argument("-W", "--WS", dest="WebServer", type=str, default="https://intradb.humanconnectome.org", help="pick server")
parser.add_argument("-P", "--Project", dest="inputProject", type=str, default="HCP_Phase2", help="pick project")
# version...
parser.add_argument('--version', action='version', version='%(prog)s 0.1')

args = parser.parse_args()
User = args.iUser
Password = args.iPassword
inputSubjects = args.inputSubjects
if (inputSubjects == None):
    print 'get all subjects...'
    
outputDir = args.outputDir

inputProject = args.inputProject
Server = args.WebServer
if (Server[-1] != '/'):
    Server = Server + '/'

showUsage = False
printLists = True
fromFileXML = False
TimeoutStep = 8
TimeoutMax = 1024
Timeout = args.Timeout
TimeoutDefault = args.Timeout
PipelineProcList = ['FunctionalHCP', 'StructuralHCP', 'DiffusionHCP']
#===============================================================================
# CLASSES
#===============================================================================
class intradb:
    """intraDB Interfacing Class"""
    def __init__( self, User, Password, Timeout, TimeoutMax, TimeoutStep ):
        self.User = User
        self.Password = Password
        self.Timeout = Timeout
        self.TimeoutMax = TimeoutMax
        self.TimeoutStep = TimeoutStep
        self.SessionId = self.fGetSessionId()
        self.ReadResult = []
        self.FileInfo = {}
        self.SubjectSessions = ()
        self.SubjectSessionsUniq = ()
        self.SessionParms = {}
    #===============================================================================
    def fGetSessionId( self ):
        """Get session id for session spawn"""
        restURL = Server + '/data/JSESSION'
        restPost = urllib.urlencode({'foo' : 'bar'})
        restRequest = urllib2.Request(restURL, restPost)
        restAuthHeader = "Basic %s" % base64.encodestring('%s:%s' % (self.User, self.Password))[:-1]
        restRequest.add_header("Authorization", restAuthHeader)

        while (self.Timeout <= self.TimeoutMax):
            try:
                restConnHandle = urllib2.urlopen(restRequest, None, self.Timeout)
                break
            except URLError, e:
                self.Timeout += self.TimeoutStep
                print 'URLError code: ' +str(e.reason)+ '. Timeout increased to ' +str(self.Timeout)+' seconds for JSESSION cookie...'
                
                
        self.SessionId = restConnHandle.read()
        return self.SessionId
    #===============================================================================
    def fGetURLString( self, URL ):
        restRequest = urllib2.Request(URL)
        restRequest.add_header("Cookie", "JSESSIONID=" + self.SessionId);
    
        while (self.Timeout <= self.TimeoutMax):
            try:
                restConnHandle = urllib2.urlopen(restRequest, None, self.Timeout)
            except HTTPError, e:
                self.Timeout += self.TimeoutStep
                print 'HTTPError code: ' +str(e.code)+ '. Timeout increased to ' +str(self.Timeout)+' seconds for ' +URL
            except URLError, e:
                self.Timeout += self.TimeoutStep
                print 'URLError code: ' +str(e.reason)+ '. Timeout increased to ' +str(self.Timeout)+' seconds for ' +URL
            except SSLError, e:
                self.Timeout += self.TimeoutStep
                print 'SSLError code: ' +str(e.message)+ '. Timeout increased to ' +str(self.Timeout)+' seconds for ' +URL
            except socket.timeout:
                self.Timeout += self.TimeoutStep
                print 'Socket timed out. Timeout increased to ' +str(self.Timeout)+ ' seconds for ' +URL
                
            else:
                try:
                    self.ReadResults = restConnHandle.read()
                    return self.ReadResults
                except HTTPError, e:
                    print 'HTTPError code: ' +str(e.code)+ '. File read timeout for ' +str(self.Timeout)+ ' seconds for ' +URL
                    
        print 'ERROR: No reasonable timeout limit could be found for ' + URL
        sys.exit()
    #===============================================================================
    def fGetFileInfo( self, URL ):
        
        restRequest = urllib2.Request(URL)
        restRequest.add_header("Cookie", "JSESSIONID=" + self.SessionId);
         
        restRequest.get_method = lambda : 'HEAD'
        try:
            restConnHandle = urllib2.urlopen(restRequest, None)
            self.FileInfo = { 'ModDate': restConnHandle.info().getheader('Last-Modified'), 'Bytes': restConnHandle.info().getheader('Content-Length') }
        except HTTPError, e:
                print 'HTTPError code: ' +str(e.code)+ '. Timeout was ' +str(self.Timeout)+' seconds for ' +URL
    
        if (self.FileInfo.get( 'Bytes' ) == None): 
            self.FileInfo[ 'Bytes' ] = '0'

        return self.FileInfo
    #===============================================================================
    def fGetSubjectSessions(self, inputSubject):
        SubjectSessions = list()
        SubjectSessionUniq = list()
        
        restURL = Server + 'data/projects/' + inputProject + '/subjects/' + inputSubject + '/experiments?format=csv&columns=ID,label'
        restResults = self.fGetURLString(restURL)
    
        restResultsSplit = restResults.split('\n')
        restEndCount = restResults.count('\n')
        restSessionHeader = restResultsSplit[0]
        restSessionHeaderSplit = restSessionHeader.split(',')
        
        labelIdx = restSessionHeaderSplit.index('"label"')
        uniqInternalId = restSessionHeaderSplit.index('"ID"')
        
        
        for j in xrange(1, restEndCount):
            currRow = restResultsSplit[j]
            
            currRowSplit = currRow.split(',')
            currSession = currRowSplit[labelIdx].replace('"', '')
            currSessionUniq = currRowSplit[uniqInternalId].replace('"', '')
            
            if (currSession.find('fnc') != -1) or (currSession.find('str') != -1) or (currSession.find('diff') != -1) or (currSession.find('xtr') != -1):
                SubjectSessions.append(currSession)
                SubjectSessionUniq.append(currSessionUniq)
#                idList, typeList, seriesList, qualityList, xnatidList = self.fGetSessionParms(inputSubject, currSession)
                
                #===============================================================
                # I think im getting off course here...
                #===============================================================
                
#                try: 
#                    seriesMatch = seriesList[seriesList.index(inputSeries)]
#                    seriesResultList.append(currSeries)
#                except ValueError:
#                    pass
    
#        return SubjectSessionUniq
        return SubjectSessions
    #===============================================================================
    def fGetSessionParms(self, inputSubject, inputSession):
        SessionParms = {}
        idList = list()
        typeList = list()
        seriesList = list()
        qualityList = list()
        xnatidList = list()

    
#        Server = 'https://intradb.humanconnectome.org/data/'
        restURL = Server + 'data/projects/' + inputProject + '/subjects/' + inputSubject + '/experiments/' + inputSession + '/scans?format=csv&columns=ID,type,series_description,quality,xnat:mrSessionData/id'
        
        restResults = self.fGetURLString(restURL)
        
        restResultsSplit = restResults.split('\n')
        restEndCount = restResults.count('\n')
        restSessionHeader = restResultsSplit[0]
        restSessionHeaderSplit = restSessionHeader.split(',')
        
        idIdx = restSessionHeaderSplit.index('"ID"')
        seriesIdx = restSessionHeaderSplit.index('"series_description"')
        typeIdx = restSessionHeaderSplit.index('"type"')
        qualityIdx = restSessionHeaderSplit.index('"quality"')
        xnatidIdx = restSessionHeaderSplit.index('"xnat:mrsessiondata/id"')
        
        for j in xrange(1, restEndCount):
            currRow = restResultsSplit[j]
            
            currRowSplit = currRow.split(',')
    
            idList.append(currRowSplit[idIdx].replace('"', ''))
            typeList.append(currRowSplit[typeIdx].replace('"', ''))
            seriesList.append(currRowSplit[seriesIdx].replace('"', ''))
            qualityList.append(currRowSplit[qualityIdx].replace('"', ''))
            xnatidList.append(currRowSplit[xnatidIdx].replace('"', ''))
            
#        return SessionParms
        return idList, typeList, seriesList, qualityList, xnatidList
#===============================================================================
# FUNCTIONS
#===============================================================================
def fPrintData( inputSubject, inputType, inputPPL, inputStatus, inputPercent, inputSession, inputLaunchTime, inputCompletionTime, inputEpoch, inputStem, outputDir ):
    
    inputListLen = len(inputSubject) 
    headerStr = ['SubjectID', 'DataType', 'Pipeline', 'Status', 'PercentComplete', 'Session', 'TimeLaunch', 'TimeComplete', 'TimeLaunchEpoch']
    fileName = inputStem + '_PipelineStatus.txt'

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    fileID = open(outputDir +os.sep+ fileName, 'wb')
    fileID.write(headerStr[0]+'\t')
    fileID.write(headerStr[1]+'\t')
    fileID.write(headerStr[2]+'\t')
    fileID.write(headerStr[3]+'\t')
    fileID.write(headerStr[4]+'\t')
    fileID.write(headerStr[5]+'\t')
    fileID.write(headerStr[6]+'\t')
    fileID.write(headerStr[7]+'\t')
    fileID.write(headerStr[8]+'\n')
    
    for i in xrange(0, inputListLen):
            fileID.write('%s' % inputSubject[i] + "\t")
            fileID.write('%s' % inputType[i] + "\t")
            fileID.write('%s' % inputPPL[i] + "\t")
            fileID.write('%s' % inputStatus[i] + "\t")
            fileID.write('%s' % inputPercent[i] + "\t")
            fileID.write('%s' % inputSession[i] + "\t")
            fileID.write('%s' % inputLaunchTime[i] + "\t")
            fileID.write('%s' % inputCompletionTime[i] + "\t")
            fileID.write('%s' % inputEpoch[i] + "\n")
#===============================================================================

ET.register_namespace('wrk', 'http://nrg.wustl.edu/workflow')
ET.register_namespace('pip', 'http://nrg.wustl.edu/pipeline')

#===============================================================================
# Init lists
# headerStr = ['SubjectID', 'DataType', 'Pipeline', 'Time [RFC822]', 'EpochTime']
# New Lists
# headerStr = ['SubjectID', 'DataType', 'Pipeline', 'Status', 'PercentComplete', 'Session', 'TimeLaunch', 'TimeCompletion', 'TimeLaunchEpoch']
#===============================================================================
SubjectList = list()
AltSubjectList = list()
DataTypeList = list()
PipelineList = list()
StatusList = list()
PercentCompleteList = list()
LaunchTimeList = list()
CompleteTimeList = list()
LaunchEpochTimeList = list()
LaunchReadTimeList = list()
FunctionalSeriesList = list()
SessionList = list()

#===============================================================================
# init interface to server
#===============================================================================
intradb = intradb(User, Password, Timeout, TimeoutMax, TimeoutStep)
if (showUsage):
    print intradb.SessionId
    print intradb.fGetURLData('https://intradb.humanconnectome.org/data/projects?format=csv')
    print intradb.fGetFileInfo('https://intradb.humanconnectome.org/data/projects/HCP_Phase2/subjects/HCPIntradb_S01642/experiments/HCPIntradb_E04497/resources/107687/files/Diffusion/data/bvals')
    print intradb.fGetFileInfo('https://intradb.humanconnectome.org/data/projects/HCP_Phase2/subjects/HCPIntradb_S01642/experiments/HCPIntradb_E04497/resources/107687/files/Diffusion/data/bvals').get('Bytes')

#===============================================================================
# Other vars...
#===============================================================================
#DirFiles = glob.glob(inputDir +os.sep+ '*.xml')

if (inputSubjects != None):
    inputSubjectsList = inputSubjects.split(',')

inputFileBase, inputFileExt = os.path.splitext('foo_bar.xml')
inputFileBase = 'foo_bar'
#inputFilePPL = inputFileBase.split('_')[0]
#inputFileBase = inputFileBase.split('_')[0]
outputFileAppend = '_Status'
outputFileExt = '.txt'
outputDirFile = '%s\\%s%s%s' % (outputDir, inputFileBase, outputFileAppend, outputFileExt)

#===============================================================================
# Looper...
#===============================================================================
for i in xrange(0, len(inputSubjectsList)):
    
    currSubject = inputSubjectsList[i]

    # get subject session...
    currSessions = intradb.fGetSubjectSessions(currSubject)
    
    # get workflow XML...
    currURL = Server + 'data/services/workflows/StructuralHCP.xml?project=HCP_Phase2&experiment=' + currSessions[0]
    currXML = intradb.fGetURLString(currURL)
    
    #===========================================================================
    # XNATRestClient -host https://hcpi-dev-mohana.nrg.mir -u XNAT_LOGIN_HERE -p PWD_HERE -remote "/data/services/workflows/StructuralHCP.xml?experiment=HCPIntradb_E05412" -m GET
    # /data/services/workflows?status=Running - gets you ALL pipelines across all projects which are in Running state
    # /data/services/workflows/StructuralHCP.xml?project=HCP_Phase2&experiment=156637_strc
    #===========================================================================
    
    if (fromFileXML):
        inputDataET = ET.parse(inputDirFile)
        inputDataRoot = inputDataET.getroot().attrib
    else:
        inputDataET = ET.fromstring(currXML)
        inputDataRoot = inputDataET.attrib
        
    currStepLaunchTime = inputDataRoot.get('current_step_launch_time')
    currLaunchTime = inputDataRoot.get('launch_time')
    timeStruct = time.strptime(currLaunchTime, '%Y-%m-%dT%H:%M:%S')
    timeLaunch = time.strftime('%d %b %Y %H:%M:%S', timeStruct)
    timeEpoch = time.mktime(timeStruct)
    timeRead = time.strftime('%d %b %Y %H:%M:%S', timeStruct)
    
    currProject = inputDataRoot.get('ExternalID')
    uniqSessionId = inputDataRoot.get('ID')
    

    
        
    for ExecutionEnvironment in inputDataET.findall('{http://nrg.wustl.edu/workflow}executionEnvironment'):
        for Pipeline in ExecutionEnvironment.findall('{http://nrg.wustl.edu/workflow}pipeline'):
            currPipeline = os.path.splitext(os.path.basename(Pipeline.text))[0]
            
            # check list of pipelines...limit because other pipelines have odd parms...
            if ( ''.join(PipelineProcList).find(currPipeline) != -1):
                PipelineList.append( currPipeline )
                PercentCompleteList.append(inputDataRoot.get('percentageComplete'))
                StatusList.append( str.lower(inputDataRoot.get('status')) )
                
                StepLaunchTime = inputDataRoot.get('current_step_launch_time')
                StepLaunchTimeStruct = time.strptime(StepLaunchTime, '%Y-%m-%dT%H:%M:%S')
                timeStep = time.strftime('%d %b %Y %H:%M:%S', StepLaunchTimeStruct)
                timeStepEpoch = time.mktime(StepLaunchTimeStruct) - (time.timezone)
                
                LaunchTime = inputDataRoot.get('launch_time')
                LaunchTimeStruct = time.strptime(LaunchTime, '%Y-%m-%dT%H:%M:%S')
                timeLaunch = time.strftime('%d %b %Y %H:%M:%S', LaunchTimeStruct)
                timeLaunchEpoch = time.mktime(LaunchTimeStruct) - (time.timezone)
                
                LaunchTimeList.append(timeLaunch)
                LaunchEpochTimeList.append(timeEpoch)
                LaunchReadTimeList.append(timeRead)
                
                for Parms in ExecutionEnvironment.findall('{http://nrg.wustl.edu/workflow}parameters'):
                    for Parm in Parms.findall('{http://nrg.wustl.edu/workflow}parameter'):
                        
                        #=======================================================
                        # Structural
                        #=======================================================
                        if (currPipeline == 'StructuralHCP'):
                            
                            if (str.lower(Parm.attrib.get('name')) == 'subjects'):
                                currSubject = Parm.text
                                if ('\n' in currSubject):
                                    currSubject = currSubject.strip()
                                SubjectList.append(currSubject)
                                
                            if (str.lower(Parm.attrib.get('name')) == 'sessionid'):
                                currSession = Parm.text
                                if ('\n' in currSession):
                                    currSession = currSession.strip()
                                SessionList.append(Parm.text.split('_')[1])
                                
                            if (str.lower(Parm.attrib.get('name')) == 'project'):
                                currProject = Parm.text
                                if ('\n' in currProject):
                                    currProject = currProject.strip()
                                
                            if (Parm.attrib.get('name') == 't1seriesdesc_1'):
                                DataTypeList.append('Structural')
                                
                            ResourceSubdir = 'Details'
                            try:
                                fileURL = Server +'data/projects/'+ currProject +'/subjects/'+ currSubject +'/experiments/'+ currSession +'/resources/'+ ResourceSubdir +'/files/'+ currPipeline +'.log'
                            except:
                                pass
                        #=======================================================
                        # Functional
                        #=======================================================
                        elif (currPipeline == 'FunctionalHCP'):
                            
                            if (Parm.attrib.get('name') == 'subjects'):
                                currSubject = Parm.text
                                SubjectList.append(Parm.text)
                                
                            if (Parm.attrib.get('name') == 'sessionid'):
                                currSession = Parm.text
                                SessionList.append(Parm.text.split('_')[1])
                                
                            if (Parm.attrib.get('name') == 'functionalseries'):
                                currFunctSeries = Parm.text
                                DataTypeList.append(Parm.text)
                                
                            try:
                                fileURL = Server +'data/projects/'+ currProject +'/subjects/'+ currSubject +'/experiments/'+ currSession +'/resources/'+ currFunctSeries +'/files/'+ currPipeline +'.log'
                            except:
                                pass                      
                        #=======================================================
                        # Diffusion      
                        #=======================================================
                        elif (currPipeline == 'DiffusionHCP'):
                            
                            if (Parm.attrib.get('name') == 'subjects'):
                                currSubject = Parm.text
                                SubjectList.append(Parm.text)
                                
                            if (Parm.attrib.get('name') == 'sessionid'):
                                currSession = Parm.text
                                SessionList.append(Parm.text.split('_')[1])
                                
                            if (str.lower(Parm.attrib.get('name')) == 'project'):
                                currProject = Parm.text
                                
                            if (Parm.attrib.get('name') == 'LR_1ScanId'):
                                DataTypeList.append('Diffusion')
                            
                            ResourceSubdir = 'Diffusion'
                            try:
                                fileURL = Server +'data/projects/'+ currProject +'/subjects/'+ currSubject +'/experiments/'+ currSession +'/resources/'+ ResourceSubdir +'/files/'+ currPipeline +'.log'
                            except:
                                pass
                            

                            
                        else:
                            
                            if (Parm.attrib.get('name') == 'subjects'):
                                AltSubjectList.append(Parm.text)
                            
#                        print Parm.attrib, Parm.text
                        
                    # Fri, 07 Dec 2012 20:59:17 GMT
                    EpochDatetimeOrig = pytz.timezone('US/Central').localize(datetime(year=1970,month=1,day=1))
                    ModDate = intradb.fGetFileInfo(fileURL).get('ModDate')
                    
                    
#                    print datetime.strptime(ModDate, '%a, %d %b %Y %H:%M:%S GMT')
#                    print pytz.timezone('US/Central').localize(datetime.strptime(ModDate, '%a, %d %b %Y %H:%M:%S GMT'))
#                    print pytz.timezone('US/Central').localize(datetime.strptime(ModDate, '%a, %d %b %Y %H:%M:%S GMT')).astimezone(pytz.utc)
                    
#                    2012-12-07 20:59:17-06:00
#                    print pytz.timezone('UTC').localize(datetime.strptime(ModDate, '%a, %d %b %Y %H:%M:%S GMT')).astimezone(pytz.timezone('US/Central'))
                    ModDatetimeStruct = pytz.timezone('UTC').localize(datetime.strptime(ModDate, '%a, %d %b %Y %H:%M:%S GMT')).astimezone(pytz.timezone('US/Central'))

                    
                    
#                    print time.timezone

                    timeModEpoch = time.mktime(ModDatetimeStruct.timetuple())
                    timeModLocalStruct = time.localtime(timeModEpoch)
                    timeModLocal = time.strftime('%d %b %Y %H:%M:%S', timeModLocalStruct)
                    timeModLocalEpoch = time.mktime(timeModLocalStruct) - (time.timezone)
                    
                    CompleteTimeList.append(timeModLocal)
                    
                    print currPipeline, (timeModEpoch - timeLaunchEpoch) / 60, (timeModLocalEpoch - timeStepEpoch) / 60
                    


if printLists:
    # headerStr = ['SubjectID', 'DataType', 'Pipeline', 'Session', 'LaunchTime', 'CompletionTime', 'EpochTime']
    # fPrintData( inputSubject, inputType, inputPPL, inputSession, inputLaunchTime, inputCompletionTime, inputEpoch, inputStem, outputDir ):
    fPrintData( SubjectList, DataTypeList, PipelineList, StatusList, PercentCompleteList, SessionList, LaunchTimeList, CompleteTimeList, LaunchEpochTimeList, 'foo', outputDir )



print("Duration: %s" % (time.time() - sTime))




