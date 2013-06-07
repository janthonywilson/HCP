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
# 100307,103414,103515,103818,105115,110411,111312,113619,114924,115320,116120,117122,118730,118932,119833,120212,123117,124422,125525,128632,129028,130013,133827,133928,134324,135932,136833,137128,138231,139637,140420,142828,143325,144226,149337,149539,150423,151223,151627,156637,158035,159239,161731,162329,163432,167743,169343,172332,175439,177746,182739,182840,185139,191437,192439,192540,193239,194140,195647,196144,197550,199150,199251,200614,201111,205119,205725,209733,210617,212318,214019,214221,214423,217429,221319,249947,293748,298051,304020,307127,329440,397760,414229,448347,485757,499566,528446,530635,552544,559053,579665,581349,585862,598568,627549,638049,645551,654754,665254,672756,677968,685058,702133,729557,732243,734045,748258,753251,788876,792564,826353,856766,857263,859671,861456,865363,872158,877168,885975,887373,889579,894673,896778,896879,901139,917255,932554,937160,984472
#===============================================================================

#===============================================================================
# PARSE INPUT
#===============================================================================
parser = argparse.ArgumentParser(description="Program to figure out pipelines status via WORKFLOW XML...")
# input...
parser.add_argument("-U", "--User", dest="User", default='tony', type=str)
parser.add_argument("-P", "--Password", dest="Password", type=str)
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
#if (',' in Pipelines):
#    PipelineList = Pipelines.split(',')
#else:
#    PipelineList = [Pipelines]
#===============================================================================
# SET UP OUTPUT
#===============================================================================
if (OutputDir[-1] != os.sep):
    OutputDir = OutputDir + os.sep
    
if not os.path.exists(OutputDir):
    os.makedirs(OutputDir)
            
# headerStr = ['SubjectID', 'DataType', 'Pipeline', 'Status', 'PercentComplete', 'Session', 'TimeLaunch', 'TimeCompletion', 'TimeLaunchEpoch']

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
PackageTail = ['_Atlas.dtseries.nii', '_Jacobian.nii.gz', '_SBRef.nii.gz', '.nii.gz']
PackageNames = ['Movement_Regressors_dt.txt', 'Movement_Regressors.txt', 'goodvoxels.nii.gz']
#===============================================================================
# Looper...
#unique sorted....
#100307,103414,103515,103818,105115,110411,111312,113619,114924,115320,116120,117122,118730,118932,119833,120212,123117,124422,125525,128632,129028,130013,133827,133928,134324,135932,136833,137128,138231,139637,140420,142828,143325,144226,149337,149539,150423,151223,151627,156637,158035,161731,162329,163432,167743,169343,172332,175439,177746,182739,182840,185139,191437,192439,192540,193239,194140,195647,196144,197550,199150,199251,200614,201111,205119,205725,209733,210617,212318,214019,214221,214423,217429,221319,249947,293748,298051,304020,307127,329440,397760,414229,448347,485757,499566,528446,530635,552544,559053,579665,581349,585862,598568,627549,638049,645551,654754,665254,672756,677968,685058,688569,702133,729557,732243,734045,748258,753251,788876,792564,826353,856766,857263,859671,861456,865363,872158,877168,885975,887373,889579,894673,896778,896879,901139,917255,932554,937160,984472,
#===============================================================================
SubjectList = list()
SubjectResource = list()
Subjects = Subjects.split(',')
SuccessBool = list()
#Subjects = ['174437']

for i in xrange(0, len(Subjects)):
    print Subjects[i]
    getHCP.Subject = Subjects[i]
    getHCP.Session = '%s_3T' % getHCP.Subject
    subjectResources = getHCP.getSubjectResources()
    subjectSessionMeta = getHCP.getSessionMeta()
    subjectSeries = subjectSessionMeta.get('Series')
    subjectType = subjectSessionMeta.get('Types')
    
    

    
    FunctionalList = list()
    for j in xrange(0, len(subjectSessionMeta.get('Types'))):
        if (subjectSessionMeta.get('Types')[j] == 'tfMRI') or (subjectSessionMeta.get('Types')[j] == 'rfMRI'):
            FunctionalList.append(subjectSessionMeta.get('Series')[j])
                
                
    if ('404 Error' in subjectResources):
        pass
    else:
        subjectResources = subjectResources.get('Names')
        SubjectListTmp = list()
        SubjectResourceTmp = list()

        for j in xrange(0, len(subjectResources)):
            if ('preproc' in subjectResources[j]) and ('Structural' not in subjectResources[j]):
                currResource = subjectResources[j]
                preprocIdx = currResource.find('_preproc')
                currSeries = currResource[0:preprocIdx]

                getHCP.Resource = subjectResources[j]
                preprocMeta = getHCP.getSubjectResourceMeta()
                tailCount = 0
                for k in xrange(0, len(PackageTail)):
                    if (currSeries + PackageTail[k] in preprocMeta.get('Name') ):
#                        print '%s%s found...' % (currSeries, PackageTail[k])
                        tailCount += 1
                
                namesCount = 0
                for k in xrange(0, len(PackageNames)):
                    if (PackageNames[k] in preprocMeta.get('Name') ):
                        print '%s found...' % (PackageNames[k])
                        namesCount += 1
                        

                if (tailCount == len(PackageTail)) and (namesCount == len(PackageNames)):
                    SubjectResourceTmp.append(subjectResources[j])
                    SubjectListTmp.append(getHCP.Subject)
    
    
    for j in xrange(0, len(FunctionalList)):
        SubjectResource.append(FunctionalList[j])
        SubjectList.append(getHCP.Subject)
        foundResource = False
        for k in xrange(0, len(SubjectResourceTmp)):
            if (FunctionalList[j] in SubjectResourceTmp[k]):
                foundResource = True
                
        if foundResource:
            SuccessBool.append('True')
        else:
            SuccessBool.append('False')
                

if Print:
    # print len(subjectWorkflowSeriesAlt), len(subjectWorkflowStatusAlt), len(subjectTimeLaunchAlt), len(subjectTimeStepAlt), len(subjectTimeLaunchEpochAlt), len(subjectTimeStepEpochAlt)
    # HeaderStr = ['SubjectID', 'Series', 'Status', 'LaunchTime', 'CompletionTime', 'LaunchEpochTime', 'CompleteEpochTime']
    with open(OutputDir + OutputFile, 'wb') as OutputFileObj:
        HeaderStr = ['SubjectID', 'Resource', 'Completion']
        OutputFileObj.write('\t'.join(HeaderStr))
        OutputFileObj.write('\n')
        for i in xrange(0, len(SubjectResource)):
            OutputFileObj.write('\t'.join( [SubjectList[i], SubjectResource[i], SuccessBool[i]] ) )
            OutputFileObj.write('\n')
                    
    
print("Duration: %s" % (time.time() - sTime))




