'''
Created on 2012-12-11

@author: jwilso01
'''
# for reasonable sorting...
#import numpy
from collections import OrderedDict
# multiplatform stuff...
import os
import sys
import argparse
# Time manipulation...
import time
from datetime import datetime
# XML parsing...
#import xml.etree.ElementTree as ET

from pyHCP import pyHCP, getHCP

sTime = time.time()

#===============================================================================
# -U tony -P passfoo -WS db.humanconnectome.org -PPL FunctionalHCP,StructrualHCP,DiffusionHCP -D C:\tmp\wfHCP -F workflowsHCP.txt
#===============================================================================

#===============================================================================
# PARSE INPUT
#===============================================================================
parser = argparse.ArgumentParser(description="Program to figure out pipelines status via WORKFLOW XML...")
# input...
parser.add_argument("-U", "--User", dest="User", default='tony', type=str)
parser.add_argument("-P", "--Password", dest="Password", type=str)
parser.add_argument("-PPL", "--Pipelines", dest="Pipelines", type=str)
parser.add_argument("-S", "--Subjects", dest="Subjects", default=None, help="pick subject, or a list of subjects")
parser.add_argument("-Prj", "--Project", dest="Project", type=str, default="HCP_Q2", help="pick project")
parser.add_argument("-WS", "--Server", dest="WebServer", type=str, default="https://intradb.humanconnectome.org", help="pick server")
# output...
parser.add_argument("-D", "--OutputDir", dest="OutputDir", type=str, help="output dir")
parser.add_argument("-F", "--OutputFile", dest="OutputFile", type=str, help="output file name")
parser.add_argument("-DF", "--OutputDirFile", dest="OutputDirFile", type=str, help="output dir and filename")
# timeout...
parser.add_argument("-t", "--time_out", dest="Timeout", type=float, default=256.0, help="change timeout")
# version...
parser.add_argument('--version', action='version', version='%(prog)s 0.0.1')

args = parser.parse_args()
User = args.User
Password = args.Password
Subjects = args.Subjects
Server = args.WebServer
Project = args.Project
Pipelines = args.Pipelines
OutputDir = args.OutputDir
OutputFile = args.OutputFile
OutputDirFile = args.OutputDirFile
#===============================================================================
# CHECK SOME INPUTS
#===============================================================================
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
#===============================================================================
# GLOBALS
#===============================================================================
Print = True
TimeoutStep = 8.0
TimeoutMax = 1024.0
Timeout = args.Timeout
TimeoutDefault = args.Timeout
PipelineProcList = ['FunctionalHCP', 'StructuralHCP', 'DiffusionHCP']
if (',' in Pipelines):
    PipelineList = Pipelines.split(',')
else:
    PipelineList = [Pipelines]
#===============================================================================
# SET UP OUTPUT
#===============================================================================
if (OutputDir[-1] != os.sep):
    OutputDir = OutputDir + os.sep
    
if not os.path.exists(OutputDir):
    os.makedirs(OutputDir)
            
# headerStr = ['SubjectID', 'DataType', 'Pipeline', 'Status', 'PercentComplete', 'Session', 'TimeLaunch', 'TimeCompletion', 'TimeLaunchEpoch']
HeaderStr = ['SubjectID', 'Series', 'Status', 'JobID', 'LaunchTime', 'StepTime', 'LaunchEpochTime', 'StepEpochTime']
if Print:
    OutputFileObj = open(OutputDir + OutputFile, 'wb')
    writeCode = OutputFileObj.write('\t'.join(HeaderStr))
#===============================================================================
# INTERFACE...
# https://db.humanconnectome.org/data/services/workflows/FunctionalHCP?columns=functionalseries,builddir&format=csv&latest_by_param=functionalseries
#===============================================================================
pyHCP = pyHCP(User, Password, Server)
getHCP = getHCP(pyHCP)
getHCP.Project = Project
getHCP.Timeout = Timeout
#===============================================================================
AllFunctionalSeries = ['rfMRI_REST1_RL', 'rfMRI_REST1_LR', 'rfMRI_REST2_RL', 'rfMRI_REST2_LR', \
             'tfMRI_WM_RL', 'tfMRI_WM_LR', 'tfMRI_GAMBLING_RL', 'tfMRI_GAMBLING_LR', 'tfMRI_MOTOR_RL', 'tfMRI_MOTOR_LR', \
             'tfMRI_LANGUAGE_RL', 'tfMRI_LANGUAGE_LR', 'tfMRI_SOCIAL_RL', 'tfMRI_SOCIAL_LR', 'tfMRI_RELATIONAL_LR', 'tfMRI_RELATIONAL_RL', 'tfMRI_EMOTION_RL', 'tfMRI_EMOTION_LR']
AllStatus = ['Failed', 'Complete', 'Running']
#===============================================================================
# Looper...
#unique sorted....
#100307,103414,103515,103818,105115,110411,111312,113619,114924,115320,116120,117122,118730,118932,119833,120212,123117,124422,125525,128632,129028,130013,133827,133928,134324,135932,136833,137128,138231,139637,140420,142828,143325,144226,149337,149539,150423,151223,151627,156637,158035,161731,162329,163432,167743,169343,172332,175439,177746,182739,182840,185139,191437,192439,192540,193239,194140,195647,196144,197550,199150,199251,200614,201111,205119,205725,209733,210617,212318,214019,214221,214423,217429,221319,249947,293748,298051,304020,307127,329440,397760,414229,448347,485757,499566,528446,530635,552544,559053,579665,581349,585862,598568,627549,638049,645551,654754,665254,672756,677968,685058,688569,702133,729557,732243,734045,748258,753251,788876,792564,826353,856766,857263,859671,861456,865363,872158,877168,885975,887373,889579,894673,896778,896879,901139,917255,932554,937160,984472,
#===============================================================================
for i in xrange(0, len(PipelineList)):
    
    getHCP.Pipeline = PipelineList[i]
    
    workflow = getHCP.getParsedWorkflow()
    # for filtering (Generator Expressions vs. List Comprehension): data_line = (data_line[i] for i in good_cols)
    Subjects = workflow.get('Subject') 
    FunctionalSeries = workflow.get('FunctionalSeries')
    JobIDs = workflow.get('JobID') 
    StepIDs = workflow.get('StepID')
    UniqSessionIds = workflow.get('ID')
    Projects = workflow.get('ExternalID')
    JobStatus = workflow.get('JobStatus')

    StepLaunchTimes = workflow.get('StepLaunchTime')
    LaunchTimes = workflow.get('LaunchTime')
    
    
    TimeLaunch = list()
    TimeLaunchEpoch = list()
    
    TimeStep = list()
    TimeStepEpoch = list()
    
    TimeComplete = list()
    
    for j in xrange(0, len(LaunchTimes)):
        timeLaunchStruct = time.strptime(LaunchTimes[j], '%Y-%m-%d %H:%M:%S.%f')
        timeLaunch = time.strftime('%d %b %Y %H:%M:%S', timeLaunchStruct)
        TimeLaunch.append(timeLaunch)
        timeLaunchEpoch = time.mktime(timeLaunchStruct)
        TimeLaunchEpoch.append(timeLaunchEpoch)
        
        timeStepStruct = time.strptime(StepLaunchTimes[j], '%Y-%m-%d %H:%M:%S.%f')
        timeStep = time.strftime('%d %b %Y %H:%M:%S', timeStepStruct)
        TimeStep.append(timeStep)
        timeStepEpoch = time.mktime(timeStepStruct)
        TimeStepEpoch.append(timeStepEpoch)
        
    # uniqify a list keeping it in order...
    uniqSubjects = list(OrderedDict.fromkeys(Subjects))
    
    for j in xrange(0, len(uniqSubjects)):
#        getHCP.Subject = '581349'
        getHCP.Subject = uniqSubjects[j]
        getHCP.Session = '%s_3T' % getHCP.Subject
        sessionMeta = getHCP.getSessionMeta()
        subjectIdx = [item for item in range(len(Subjects)) if Subjects[item] == getHCP.Subject]
        subjectWorkflowSeries = [FunctionalSeries[i] for i in subjectIdx]
        subjectWorkflowStatus = [JobStatus[i] for i in subjectIdx]
        subjectTimeStepEpoch = [TimeStepEpoch[i] for i in subjectIdx]
        subjectTimeStep = [StepLaunchTimes[i] for i in subjectIdx]
        subjectTimeLaunchEpoch = [TimeLaunchEpoch[i] for i in subjectIdx]
        subjectTimeLaunch = [LaunchTimes[i] for i in subjectIdx]
        subjectJobId = [JobIDs[i] for i in subjectIdx]

        
        if (getHCP.Pipeline == 'FunctionalHCP'):
            subjectSessionSeries = list()
            for k in xrange(0, len(sessionMeta.get('Types'))):
                if (sessionMeta.get('Types')[k] == 'tfMRI') or (sessionMeta.get('Types')[k] == 'rfMRI'):
                    subjectSessionSeries.append(sessionMeta.get('Series')[k])
                    
                    # Too network/filesystem intensive....
#                    getHCP.Resource = '%s_preproc' % (sessionMeta.get('Series')[k])
#                    subjectMeta = getHCP.getSubjectResourceMeta()
              
            subjectResourcesSeries = list()
            subjectResources = getHCP.getSubjectResources()
            for k in xrange(0, len(subjectResources.get('Names'))):
                if ('preproc' in subjectResources.get('Names')[k]) and not ('Structural' in subjectResources.get('Names')[k]):
                    subjectResourcesSeries.append(subjectResources.get('Names')[k])
            
            subjectSessionSeriesAlt = list()
            subjectWorkflowSeriesAlt = list()
            subjectWorkflowStatusAlt = list()
            subjectTimeStepEpochAlt = list()
            subjectTimeStepAlt = list()
            subjectTimeLaunchEpochAlt = list()
            subjectTimeLaunchAlt = list()
            subjectJobIdAlt = list()

            
            for k in xrange(0, len(AllFunctionalSeries)):
                
                if (AllFunctionalSeries[k] not in subjectWorkflowSeries):
                    subjectWorkflowSeriesAlt.append(AllFunctionalSeries[k])
                    # Check Resource List...substring matching so looping.
                    ResourcePresent = False
                    for m in xrange(0, len(subjectResourcesSeries)):
                        if (AllFunctionalSeries[k] in subjectResourcesSeries[m]):
                            ResourcePresent = True
                            print 'Absent Workflow: Resource for series %s found for subject %s ' % (AllFunctionalSeries[k], getHCP.Subject)
                            
                    if ResourcePresent:
                        subjectWorkflowStatusAlt.append('Complete')
                        subjectJobIdAlt.append('PUT')
                        
                        # This should be replaced with a call to getResourceMeta() for time...
                        #=======================================================
                        # getHCP.Resource = '%s_preproc' % (AllFunctionalSeries[k])
                        # subjectMeta = getHCP.getSubjectResourceMeta()
                        # logIdx = subjectMeta.get('Name').index('FunctionalHCP.log')
                        # logMeta = getHCP.getFileInfo(subjectMeta.get('URI')[logIdx])
                        #=======================================================
                        
                        subjectTimeStepEpochAlt.append('ResourceTimeStepEpoch')
                        subjectTimeStepAlt.append('ResourceTimeStepReal')
                        subjectTimeLaunchEpochAlt.append('ResourceTimeLaunchEpoch')
                        subjectTimeLaunchAlt.append('ResourceTimeLaunch')
                        
                    elif (AllFunctionalSeries[k] not in subjectSessionSeries):
                        subjectWorkflowStatusAlt.append('Absent')
                        subjectTimeStepEpochAlt.append('0000000001')
                        subjectTimeStepAlt.append('1970-01-01 00:00:01.0')
                        subjectTimeLaunchEpochAlt.append('0000000001')
                        subjectTimeLaunchAlt.append('1970-01-01 00:00:01.0')
                        subjectJobIdAlt.append('AbsentJobID')
                        
                    else:
                        subjectWorkflowStatusAlt.append('NoLaunch')
                        subjectTimeStepEpochAlt.append('0000000001')
                        subjectTimeStepAlt.append('1970-01-01 00:00:01.0')
                        subjectTimeLaunchEpochAlt.append('0000000001')
                        subjectTimeLaunchAlt.append('1970-01-01 00:00:01.0')
                        subjectJobIdAlt.append('NoLaunchJobID')
                   
                elif (AllFunctionalSeries[k] in subjectWorkflowSeries): 
                    # if failed, check resource....
                    ResourcePresent = False
                    if (subjectWorkflowStatus[subjectWorkflowSeries.index(AllFunctionalSeries[k])] == 'Failed') or (subjectWorkflowStatus[subjectWorkflowSeries.index(AllFunctionalSeries[k])] == 'Running'):
                        for m in xrange(0, len(subjectResourcesSeries)):
                            if (AllFunctionalSeries[k] in subjectResourcesSeries[m]):
                                ResourcePresent = True
                                print 'Failed/Running Workflow: Resource for series %s found for subject %s ' % (AllFunctionalSeries[k], getHCP.Subject)

                        
                    if ResourcePresent:
                        subjectWorkflowSeriesAlt.append(AllFunctionalSeries[k])
                        subjectWorkflowStatusAlt.append('Complete')
                        subjectJobIdAlt.append('PUT')
                        # This should be replaced with a call to getResourceMeta() for time...
                        subjectTimeStepEpochAlt.append('ResourceTimeStepEpoch')
                        subjectTimeStepAlt.append('ResourceTimeStepReal')
                        subjectTimeLaunchEpochAlt.append('ResourceTimeLaunchEpoch')
                        subjectTimeLaunchAlt.append('ResourceTimeLaunch')
                
                    else:
                        subjectWorkflowSeriesAlt.append(AllFunctionalSeries[k])
                        # Catch the queuing issues before 09/05/2013... 
                        if (subjectTimeLaunchEpoch[subjectWorkflowSeries.index(AllFunctionalSeries[k])] < 1368136800) and (subjectWorkflowStatus[subjectWorkflowSeries.index(AllFunctionalSeries[k])] == 'Running'):
                            subjectWorkflowStatusAlt.append('NotTrusted')
                            # LAUNCH TIMES...
                            subjectTimeLaunchEpochAlt.append('2971215073')
                            subjectTimeLaunchAlt.append('2064-02-26 01:13:13.0')
                        elif (subjectTimeLaunchEpoch[subjectWorkflowSeries.index(AllFunctionalSeries[k])] < 1368136800) and (subjectWorkflowStatus[subjectWorkflowSeries.index(AllFunctionalSeries[k])] == 'Complete'):
                            subjectWorkflowStatusAlt.append(subjectWorkflowStatus[subjectWorkflowSeries.index(AllFunctionalSeries[k])])
                            # LAUNCH TIMES...
                            subjectTimeLaunchEpochAlt.append('2971215073')
                            subjectTimeLaunchAlt.append('2064-02-26 01:13:13.0')
                        else:
                            subjectWorkflowStatusAlt.append(subjectWorkflowStatus[subjectWorkflowSeries.index(AllFunctionalSeries[k])])
                            # LAUNCH TIMES...
                            subjectTimeLaunchEpochAlt.append(subjectTimeLaunchEpoch[subjectWorkflowSeries.index(AllFunctionalSeries[k])])
                            subjectTimeLaunchAlt.append(subjectTimeLaunch[subjectWorkflowSeries.index(AllFunctionalSeries[k])])
                        
                        # STEP TIMES...
                        subjectTimeStepEpochAlt.append(subjectTimeStepEpoch[subjectWorkflowSeries.index(AllFunctionalSeries[k])])
                        subjectTimeStepAlt.append(subjectTimeStep[subjectWorkflowSeries.index(AllFunctionalSeries[k])])
                        
                        subjectJobIdAlt.append(subjectJobId[subjectWorkflowSeries.index(AllFunctionalSeries[k])])
                        
                else:
                    subjectWorkflowSeriesAlt.append('Fault')
                    subjectWorkflowStatusAlt.append('Fault')
                    subjectJobIdAlt.append('Fault')
                    # START TIMES...
                    subjectTimeStepEpochAlt.append('Fault')
                    subjectTimeStepAlt.append('Fault')
                    # LAUNCH TIMES...
                    subjectTimeLaunchEpochAlt.append('Fault')
                    subjectTimeLaunchAlt.append('Fault')
                    
                    

            if Print:
                # print len(subjectWorkflowSeriesAlt), len(subjectWorkflowStatusAlt), len(subjectTimeLaunchAlt), len(subjectTimeStepAlt), len(subjectTimeLaunchEpochAlt), len(subjectTimeStepEpochAlt)
                # HeaderStr = ['SubjectID', 'Series', 'Status', 'LaunchTime', 'CompletionTime', 'LaunchEpochTime', 'CompleteEpochTime']
                for k in xrange(0, len(subjectWorkflowStatusAlt)):
                    OutputFileObj.write('\n' )
                    OutputFileObj.write('\t'.join( [getHCP.Subject, subjectWorkflowSeriesAlt[k], subjectWorkflowStatusAlt[k], subjectJobIdAlt[k], subjectTimeLaunchAlt[k], subjectTimeStepAlt[k], str(subjectTimeLaunchEpochAlt[k]), str(subjectTimeStepEpochAlt[k])]) )
                    
    
if Print:
    OutputFileObj.close()
    
print("Duration: %s" % (time.time() - sTime))




