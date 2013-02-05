'''
Created on 2012-12-11

@author: jwilso01
'''


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

from getHCP import getHCP

sTime = time.time()

#===============================================================================
# -iD C:\Users\jwilso01\workspace\data\xml\workflow\ -iP asdf -oD C:\tmp\work -W http://hcpi-dev-cuda00.nrg.mir/
#===============================================================================

#===============================================================================
# PARSE INPUT
#===============================================================================
parser = argparse.ArgumentParser(description="Program to figure out pipelines status via WORKFLOW XML...")
# input...
parser.add_argument("-U", "--User", dest="User", default='tony', type=str)
parser.add_argument("-P", "--Password", dest="Password", type=str)
parser.add_argument("-S", "--Subjects", dest="Subjects", default=None, help="pick subject, or a list of subjects")
parser.add_argument("-Proj", "--Project", dest="Project", type=str, default="HCP_Phase2", help="pick project")
parser.add_argument("-W", "--WS", dest="WebServer", type=str, default="https://intradb.humanconnectome.org", help="pick server")
# output...
parser.add_argument("-D", "--outputDir", dest="outputDir", type=str, help="where do you want to write output tab-text")
# timeout...
parser.add_argument("-t", "--time_out", dest="Timeout", type=float, default=16.0, help="change timeout")
# version...
parser.add_argument('--version', action='version', version='%(prog)s 0.4.1')

args = parser.parse_args()
User = args.User
Password = args.Password
Subjects = args.Subjects
Server = args.WebServer
Project = args.Project
    
outputDir = os.path.normpath(args.outputDir) + os.sep
#===============================================================================
# CLEAN UP INPUT
#===============================================================================
Server.strip()
if (Server[-1] != '/'):
    Server = Server + '/'
if (Server.find('http') == -1):
    Server = 'https://' + Server
    
if (outputDir[-1] != os.sep):
    outputDir = outputDir + os.sep
#===============================================================================
# GLOBALS
#===============================================================================
printLists = True
fromFileXML = False
TimeoutStep = 8.0
TimeoutMax = 1024.0
Timeout = args.Timeout
TimeoutDefault = args.Timeout
SessionProcList = ['strc', 'fnc', 'diff', 'xtr']
PipelineProcList = ['FunctionalHCP', 'StructuralHCP', 'DiffusionHCP']
#===============================================================================
# FUNCTIONS
#===============================================================================  
def fParseWorkflowData( CSV ):
    """Get Workflow Data"""

    # "label","label","id","externalid","pipeline_name","launch_time","jobid","status","current_step_launch_time","current_step_id","builddir"
    Subject = list()
    Session = list()
    ID = list()
    ExternalID = list()
    PipelinePathName = list()
    LaunchTime = list()
    JobID = list()
    JobStatus = list()
    StepLaunchTime = list()
    StepID = list()
    BuildDir = list()

    
    restResultsSplit = CSV.split('\n')
    restEndCount = CSV.count('\n')
    restSessionHeader = restResultsSplit[0]
    restSessionHeaderSplit = restSessionHeader.split(',')
    
    labelIdx = restSessionHeaderSplit.index('"label"')
    idIdx = restSessionHeaderSplit.index('"id"')
    pipelineIdx = restSessionHeaderSplit.index('"pipeline_name"')
    launchTimeIdx = restSessionHeaderSplit.index('"launch_time"')
    jobIdx = restSessionHeaderSplit.index('"jobid"')
    statusIdx = restSessionHeaderSplit.index('"status"')
    stepLaunchIdx = restSessionHeaderSplit.index('"current_step_launch_time"')
    stepIdx = restSessionHeaderSplit.index('"current_step_id"')
    buildDirIdx = restSessionHeaderSplit.index('"builddir"')
    
    for j in xrange(1, restEndCount):
        currRow = restResultsSplit[j]
        
        currRowSplit = currRow.split(',')

        Subject.append(currRowSplit[labelIdx].replace('"', ''))
        ID.append(currRowSplit[idIdx].replace('"', ''))
        PipelinePathName.append(currRowSplit[pipelineIdx].replace('"', ''))
        LaunchTime.append(currRowSplit[launchTimeIdx].replace('"', ''))
        JobID.append(currRowSplit[jobIdx].replace('"', ''))
        JobStatus.append(currRowSplit[statusIdx].replace('"', ''))
        StepLaunchTime.append(currRowSplit[stepLaunchIdx].replace('"', ''))
        StepID.append(currRowSplit[stepIdx].replace('"', ''))
        BuildDir.append(currRowSplit[buildDirIdx].replace('"', ''))
        
    ParsedData = {'Subject': Subject, 'ID': ID, 'PipelinePathName': PipelinePathName, 'LaunchTime': LaunchTime, 'JobID': JobID, 'JobStatus': JobStatus, 'StepLaunchTime': StepLaunchTime, 'StepID': StepID, 'BuildDir': BuildDir}
    return ParsedData
#===============================================================================    
def fPrint( outputDirFile, headerStr, *args ):

    outputFile = os.path.basename(outputDirFile)
    ouputFileBase, outputFileExt = os.path.splitext(outputFile)
    outputDir = os.path.dirname(os.path.normpath(outputDirFile)) + os.sep

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
        
    fileID = open(outputDirFile, 'wb')
        
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

#ET.register_namespace('wrk', 'http://nrg.wustl.edu/workflow')
#ET.register_namespace('pip', 'http://nrg.wustl.edu/pipeline')

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
getHCP = getHCP(User, Password, Server)
getHCP.Project = Project

if (Subjects == None):
    SubjectsList = getHCP.getSubjects()
elif (Subjects != None):
    SubjectsList = Subjects.split(',')
#===============================================================================
# Setup output...
#===============================================================================
#DirFiles = glob.glob(inputDir +os.sep+ '*.xml')
#inputFileBase, inputFileExt = os.path.splitext('foo_bar.xml')

outputFileBase = 'PipelineStatus'
if (len(SubjectsList) == 1):
    outputFileAppend = SubjectsList[0]
else:
    outputFileAppend = Project
    
outputFileExt = '.txt'
outputDirFile = '%s\\%s%s%s' % (outputDir, outputFileBase, outputFileAppend, outputFileExt)

#===============================================================================
# Looper...
#===============================================================================
for i in xrange(0, len(SubjectsList)):
    
    getHCP.Subject = SubjectsList[i]

    # get subject sessions...
    Sessions = getHCP.getSubjectSessions()
    currSessions = Sessions.get('Sessions')
    currSessionsType = Sessions.get('SessionTypes')


    for j in xrange(0, len(currSessions)):
        getHCP.Session = currSessions[j]
        
        # idList, typeList, seriesList, qualityList, xnatidList
        currSessionMeta = getHCP.getSessionMeta()
        
        if ( ''.join(SessionProcList).find(getHCP.Session.split('_')[1]) != -1):
                if (getHCP.Session.find('fnc') != -1):
                    currPipeline = 'FunctionalHCP'
                elif (getHCP.Session.find('strc') != -1):
                    currPipeline = 'StructuralHCP'
                elif (getHCP.Session.find('diff') != -1):
                    currPipeline = 'DiffusionHCP'
                elif (getHCP.Session.find('xtr') != -1):
                    # now find the session type...
                    if 'tfMRI' in currSessionMeta.get('Types'):
                        currPipeline = 'FunctionalHCP'
                    elif 'dMRI' in currSessionMeta.get('Types'):
                        currPipeline = 'DiffusionHCP'
                    elif 'T1w' in currSessionMeta.get('Types'):
                        currPipeline = 'StructuralHCP'
        
        #===============================================================
        # XNATRestClient -u mohanar -p admin -host https://hcpi-dev-mohana.nrg.mir -m GET -remote "/data/services/workflows/StructuralHCP?display=LATEST&columns=builddir&format=csv" > wrk_param_latest_structural.dat
        # XNATRestClient -u mohanar -p admin -host https://hcpi-dev-mohana.nrg.mir -m GET -remote "/data/services/workflows/Diffusion?display=LATEST&columns=builddir&format=csv" > wrk_param_latest_diffusion.dat
        # XNATRestClient -u mohanar -p admin -host https://hcpi-dev-mohana.nrg.mir -m GET -remote "/data/services/workflows/FunctionalHCP?columns=functionalseries,builddir&format=csv&latest_by_param=functionalseries" > wrk_param_latest_functional.dat
        #===============================================================
        
        # get workflow XML...
        if (currPipeline == 'FunctionalHCP'):
            currURL = Server + 'data/services/workflows/' + currPipeline + '?columns=functionalseries,builddir&format=csv&latest_by_param=functionalseries'
        else:
            currURL = Server + 'data/services/workflows/' + currPipeline + '?display=LATEST&columns=builddir&format=csv'
            
        currCSV = getHCP.getURLString(currURL)
        
        ParsedWorkflow = fParseWorkflowData(currCSV)
            
#        if (currCSV[0:9] == '<ResultSet>'):
#                break
        

        
        
            
        if (currCSV != '404 Error'):
            if (fromFileXML):
                inputDataET = ET.parse(inputDirFile)
                inputDataRoot = inputDataET.getroot().attrib
            else:
                inputDataET = ET.fromstring(currXML)
                inputDataRoot = inputDataET.attrib
                
#                currStepLaunchTime = inputDataRoot.get('current_step_launch_time')
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
                        JobIdList.append(inputDataRoot.get('jobID'))
                        
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
                                        SessionList.append(currSession.split('_')[1])
                                        
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
                                        if ('\n' in currSubject):
                                            currSubject = currSubject.strip()
                                        SubjectList.append(currSubject)
                                        
                                    if (Parm.attrib.get('name') == 'sessionid'):
                                        currSession = Parm.text
                                        if ('\n' in currSession):
                                            currSession = currSession.strip()
                                        SessionList.append(currSession.split('_')[1])
                                        
                                    if (str.lower(Parm.attrib.get('name')) == 'project'):
                                        currProject = Parm.text
                                        if ('\n' in currProject):
                                            currProject = currProject.strip()
                                            
                                        
                                    if (Parm.attrib.get('name') == 'functionalseries'):
                                        currFunctSeries = Parm.text
                                        if ('\n' in currFunctSeries):
                                            currFunctSeries = currFunctSeries.strip()
                                        DataTypeList.append(currFunctSeries)
                                        
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
                                        if ('\n' in currSubject):
                                            currSubject = currSubject.strip()
                                        SubjectList.append(currSubject)
                                        
                                    if (Parm.attrib.get('name') == 'sessionid'):
                                        currSession = Parm.text
                                        if ('\n' in currSession):
                                            currSession = currSession.strip()
                                        SessionList.append(currSession.split('_')[1])
                                        
                                    if (str.lower(Parm.attrib.get('name')) == 'project'):
                                        currProject = Parm.text
                                        if ('\n' in currProject):
                                            currProject = currProject.strip()
                                        
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
                            fileInfo = getHCP.getFileInfo(fileURL)
                            if (fileInfo != '404 Error'):
                                ModDate = fileInfo.get('ModDate')
                            else:
                                ModDate = 'Fri, 13 Dec 2013 19:13:13 GMT'
                            
        #                    print datetime.strptime(ModDate, '%a, %d %b %Y %H:%M:%S GMT')
        #                    print pytz.timezone('US/Central').localize(datetime.strptime(ModDate, '%a, %d %b %Y %H:%M:%S GMT'))
        #                    print pytz.timezone('US/Central').localize(datetime.strptime(ModDate, '%a, %d %b %Y %H:%M:%S GMT')).astimezone(pytz.utc)
                            
        #                    2012-12-07 20:59:17-06:00
                            ModDatetimeStruct = pytz.timezone('UTC').localize(datetime.strptime(ModDate, '%a, %d %b %Y %H:%M:%S GMT')).astimezone(pytz.timezone('US/Central'))
        
        
                            timeModEpoch = time.mktime(ModDatetimeStruct.timetuple())
                            timeModLocalStruct = time.localtime(timeModEpoch)
                            timeModLocal = time.strftime('%d %b %Y %H:%M:%S', timeModLocalStruct)
                            timeModLocalEpoch = time.mktime(timeModLocalStruct) - (time.timezone)
                            
                            CompleteTimeList.append(timeModLocal)
                            
#                                print currPipeline, (timeModEpoch - timeLaunchEpoch) / 60, (timeModLocalEpoch - timeStepEpoch) / 60, (timeModLocalEpoch - timeLaunchEpoch) / 60
                        


if printLists:
    fPrint( outputDirFile, ['SubjectID', 'DataType', 'Pipeline', 'Status', 'PercentComplete', 'JobID', 'Session', 'TimeLaunch', 'TimeComplete', 'TimeLaunchEpoch'], SubjectList, DataTypeList, PipelineList, StatusList, PercentCompleteList, JobIdList, SessionList, LaunchTimeList, CompleteTimeList, LaunchEpochTimeList )
#    fPrintData( SubjectList, DataTypeList, PipelineList, StatusList, PercentCompleteList, SessionList, LaunchTimeList, CompleteTimeList, LaunchEpochTimeList, outputDirFile )



print("Duration: %s" % (time.time() - sTime))




