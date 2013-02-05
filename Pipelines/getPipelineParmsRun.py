'''
Created on Sep 3, 2012

@author: jwilso01
'''
import base64
import sys
import xml
import os
import time
import numpy
#import array
import socket
import argparse
import datetime
import xml.etree.ElementTree as ET

from pyHCP import getHCP


sTime = time.time()

#===============================================================================
# PARSE INPUT
#===============================================================================
parser = argparse.ArgumentParser(description="Script to generate proper command for XNAT functional pipeline lauching ...")

#MANDATORY....
parser.add_argument("-User", "--User", dest="User", default='tony', type=str)
parser.add_argument("-Password", "--Password", dest="Password", default='none', type=str)
parser.add_argument("-Pipeline", "--Pipeline", dest="Pipeline", default='fMRIVolume', type=str)
parser.add_argument("-Subjects", "--Subjects", dest="Subjects", default='00', type=str)
parser.add_argument("-Server", "--Server", dest="Server", default='http://hcpi-dev-cuda00.nrg.mir/', type=str)
parser.add_argument('--version', action='version', version='%(prog)s 0.0.1')
# IF Pipeline == Functional...
parser.add_argument("-FunctSeries", "--FunctSeries", dest="FunctSeries", default='BOLD_00', type=str)
#END MANDATORY....
parser.add_argument("-Project", "--Project", dest="Project", default='HCP_Phase2', type=str)
parser.add_argument("-Shadow", "--Shadow", dest="Shadow", default=None, type=str)
parser.add_argument("-Build", "--Build", dest="Build", default=None, type=str)


parser.add_argument("-Production", "--Production", dest="iProduction", default=False, type=bool)


parser.add_argument("-HCPid", "--HCPid", dest="iHCPid", default='HCPIntradb_E00000', type=str)
parser.add_argument("-Label", "--Label", dest="iLabel", default='00', type=str)
parser.add_argument("-FuncScanId", "--FuncScanId", dest="iFuncScanId", default='00', type=str)
parser.add_argument("-ScoutScanId", "--ScoutScanId", dest="iScoutScanId", default='00', type=str)
parser.add_argument("-UnwarpDir", "--UnwarpDir", dest="iUnwarpDir", default='x', type=str)
# structural...
parser.add_argument("-StructSeries", "--SturctSeries", dest="iStructSeries", default='T1w_MPR1', type=str)
parser.add_argument("-T1wScanId_1", "--T1wScanId_1", dest="iT1wScanId_1", default='00', type=str)
parser.add_argument("-T1wScanId_2", "--T1wScanId_2", dest="iT1wScanId_2", default='00', type=str)
parser.add_argument("-T2wScanId_1", "--T2wScanId_1", dest="iT2wScanId_1", default='00', type=str)
parser.add_argument("-T2wScanId_2", "--T2wScanId_2", dest="iT2wScanId_2", default='00', type=str)
parser.add_argument("-T1wSeriesDesc_1", "--T1wSeriesDesc_1", dest="iT1wSeriesDesc_1", default='T1w_MPR1', type=str)
parser.add_argument("-T1wSeriesDesc_2", "--T1wSeriesDesc_2", dest="iT1wSeriesDesc_2", default='T1w_MPR2', type=str)
parser.add_argument("-T2wSeriesDesc_1", "--T2wSeriesDesc_1", dest="iT2wSeriesDesc_1", default='T2w_SPC1', type=str)
parser.add_argument("-T2wSeriesDesc_2", "--T2wSeriesDesc_2", dest="iT2wSeriesDesc_2", default='T2w_SPC2', type=str)
# diffusion...
parser.add_argument("-DiffSeries", "--DiffSeries", dest="iDiffSeries", default='DWI_RL_dir95', type=str)

args = parser.parse_args()

#MANDATORY....
User = args.User
Password = args.Password
Server = args.Server
Pipeline = args.Pipeline
Subjects = args.Subjects
FunctSeries = args.FunctSeries
#END MANDATORY....START ALT
Project = args.Project
Shadow = args.Shadow
Build = args.Build
#===============================================================================
# DEPRICATED...
iHCPid = args.iHCPid
iLabel = args.iLabel
iFuncScanId = args.iFuncScanId
iScoutScanId = args.iScoutScanId
iUnwarpDir = args.iUnwarpDir

# structural...
iStructSeries = args.iStructSeries
iT1wScanId_1 = args.iT1wScanId_1
iT1wScanId_2 = args.iT1wScanId_2
iT2wScanId_1 = args.iT2wScanId_1
iT2wScanId_2 = args.iT2wScanId_2
iT1wSeriesDesc_1 = args.iT1wSeriesDesc_1
iT1wSeriesDesc_2 = args.iT1wSeriesDesc_2 
iT2wSeriesDesc_1 = args.iT2wSeriesDesc_1
iT2wSeriesDesc_2 = args.iT2wSeriesDesc_2

#diffusion...
iDiffSeries = args.iDiffSeries
#===============================================================================


#if (iPipeline.find('Funct') != -1):
#    iStructFunctDiffSeries = iFunctSeries
#elif (iPipeline.find('Struct') != -1):
#    iStructFunctDiffSeries = iStructSeries
#elif (iPipeline.find('Diff') != -1):
#    iStructFunctDiffSeries = iDiffSeries
#    
#tmpStructFunctDiffSeries = iStructFunctDiffSeries

#===============================================================================
# STATIC PARMS...
#===============================================================================
JobSubmitter = '/data/intradb/pipeline/bin/PipelineJobSubmitter '
PipelineLauncher = '/data/intradb/pipeline/bin/XnatPipelineLauncher '
DataType = '-dataType xnat:mrSessionData '
SupressNotify = '-supressNotification '
NotifyUser = '-notify wilsont@mir.wustl.edu ' 
NotifyAdmin = '-notify db-admin@humanconnectome.org '
MailHost = '-parameter mailhost=mail.nrg.wustl.edu '
UserFullName = '-parameter userfullname=T.Wilson ' 
XnatServer = '-parameter xnatserver=HCPIntradb '
AdminEamil = '-parameter adminemail=db-admin@humanconnectome.org '
UserEmail = '-parameter useremail=wilsont@mir.wustl.edu '
#===============================================================================
# HACK for REST RL/LR...
#===============================================================================
preChangeRLList = ('100307','111312','114924','119833','125525','138231','144266','150423','159239','162329',\
                   '167743','174437','185139','192439','197550','199251','217429','249947','255639','329440',\
                   '355542','499566','585862','611231','665254','672756','792564','826353','877168','896778')
#===============================================================================
# INTERFACE...
#===============================================================================
getHCP = getHCP(User, Password, Server)
getHCP.Project = Project
#===============================================================================
# FUNCTIONS...
#===============================================================================
def fGetAllIndices(inputVal, inputList):
    allIndicesList = list()
    currIdx = -1
    while True:
        try:
            currIdx = inputList.index(inputVal, currIdx + 1)
            allIndicesList.append(currIdx)
        except ValueError:
            break
    return numpy.asarray(allIndicesList, dtype=numpy.int)
#===============================================================================

SleepTime = 3
SubjectsList = Subjects.split(',')
if (Shadow != None):
    ShadowList = Shadow.split(',')
else:
    ShadowList = ('')
    
if (Build != None):
    BuildList = Build.split(',')
else:
    BuildList = ('')
    
#SeriesList = iStructFunctDiffSeries.split(',')

if (len(SubjectsList) > len(ShadowList)):
#    ShadowArray = numpy.tile(ShadowList, (numpy.ceil(len(SubjectsList) / len(ShadowList))))
    ShadowArray = numpy.tile(ShadowList, (numpy.ceil(len(SubjectsList))))
else: 
    ShadowArray = numpy.tile(ShadowList, (1))

if (len(SubjectsList) > len(BuildList)):
    BuildArray = numpy.tile(BuildList, (numpy.ceil(len(SubjectsList))))
else: 
    BuildArray = numpy.tile(BuildList, 1)
    
if (Pipeline == 'FunctionalHCP'): PipelineSubString = 'fnc'
elif (Pipeline == 'StructuralHCP'): PipelineSubString = 'strc'
elif (Pipeline == 'DiffusionHCP'): PipelineSubString = 'diff'
else: PipelineSubString = None
    
UsableList = ['good', 'excellent', 'usable', 'undetermined']

SeriesList = FunctSeries.split(',')

linIdx = 0
for h in xrange(0, len(SubjectsList)): 
    getHCP.Subject = SubjectsList[h]
    SubjectSessions = getHCP.getSubjectSessions()
    
    if (PipelineSubString not in SubjectSessions.get('Types')):
        print 'ERROR: No ' +Pipeline+ ' session could be found for subject ' +getHCP.Subject
    
    for i in xrange(0, len(SeriesList)):
        linIdx += 1
        
        if (Pipeline.find('Funct') != -1):
            currFunctSeries = SeriesList[i]
            tmpFunctSeries = currFunctSeries
        
            if ( (tmpFunctSeries.find('REST') != -1) and (preChangeRLList.count(getHCP.Subject) > 0) ):
                if (tmpFunctSeries.find('REST3') != -1) or (tmpFunctSeries.find('REST1') != -1):
                    currFunctSeries = 'BOLD_'+ tmpFunctSeries +'_RL'
                elif (tmpFunctSeries.find('REST4') != -1) or (tmpFunctSeries.find('REST2') != -1):
                    currFunctSeries = 'BOLD_'+ tmpFunctSeries +'_LR'
            elif ( (tmpFunctSeries.find('REST') != -1) and (preChangeRLList.count(getHCP.Subject) == 0) ):
                if (tmpFunctSeries.find('REST3') != -1) or (tmpFunctSeries.find('REST2') != -1):
                    currFunctSeries = 'BOLD_'+ tmpFunctSeries +'_LR'
                elif (tmpFunctSeries.find('REST4') != -1) or (tmpFunctSeries.find('REST1') != -1):
                    currFunctSeries = 'BOLD_'+ tmpFunctSeries +'_RL'

            # Check that the functional scan is usable, return the session with usable scan...
            for j in xrange(0, len(SubjectSessions.get('Sessions'))):
                getHCP.Session = SubjectSessions.get('Sessions')[j]
                sessionMeta = getHCP.getSessionMeta()
                if (currFunctSeries in sessionMeta.get('Series')):
                    currUsability = sessionMeta.get('Quality')[sessionMeta.get('Series').index(currFunctSeries)]
                    if (currUsability in UsableList): break
                
            SessionId = getHCP.Session

        elif (Pipeline.find('Diff') != -1):
            diffSessionIdx = 0
#            print SubjectSessions.get('Types').count('diff')
            for j in xrange(0, SubjectSessions.get('Types').count('diff')):
                diffSessionIdx =+ SubjectSessions.get('Types').index('diff')
                getHCP.Session = SubjectSessions.get('Sessions')[diffSessionIdx]
                sessionMeta = getHCP.getSessionMeta()
                
        elif (Pipeline.find('Struct') != -1):
            print 'Looking at StructuralHCP...'
        elif (Pipeline.find('Funct') != -1):
            print 'Looking at FunctionalHCP...'
        else:
            print 'ERROR: Pipline not found...'
            
        
        
        if (BuildArray.size >= len(SubjectsList)):
            BuildDirRoot = '/data/intradb/build' + str(BuildArray[h]) + '/' + Project + '/'
        else:
            BuildDirRoot = '/data/intradb/build/' + Project + '/'
#            print iBuildDirRoot
        
        launcherProject = '-parameter project=%s ' % Project 
        launcherPipeline = '-pipeline /data/intradb/pipeline/catalog/%s/%s.xml ' % (Pipeline, Pipeline)
        launcherUser = '-u %s ' % User 
        launcherPassword = '-pwd %s ' % Password 
        launcherLabel = '-label %s ' % getHCP.Session
        launcherHCPid = '-id %s ' % sessionMeta.get('XNATID')[0]
        launcherXnatId = '-parameter xnat_id=%s ' % sessionMeta.get('XNATID')[0] 
        launcherSession = '-parameter sessionid=%s ' % getHCP.Session 
        launcherSubject = '-parameter subjects=%s ' % getHCP.Subject

        
        if (len(SeriesList) > 1):
            currBuildDir = BuildDirRoot + str(numpy.asarray(round(time.time()), dtype=numpy.uint64))
        else:
            currBuildDir = BuildDirRoot + str(numpy.asarray(round(sTime), dtype=numpy.uint64))
        BuildDir = '-parameter builddir=' + currBuildDir + ' '

        
        if not os.path.exists(currBuildDir + os.sep + getHCP.Subject) and sys.platform != 'win32':
            os.makedirs(currBuildDir + os.sep + getHCP.Subject)
            
        RedirectionStr = ' > ' + currBuildDir.replace(' ', '') + os.sep + getHCP.Subject + os.sep + Pipeline + 'LaunchSTDOUT.txt'
        
        if (socket.gethostname() == 'intradb.humanconnectome.org') and (Shadow != None):
                Host = '-host https://intradb-shadow' + ShadowArray[h] + '.nrg.mir '
        else: 
            Host = '-host https://' + socket.gethostname() + ' '

                
            
        #===============================================================================
        # START BUILDING STRING....
        #===============================================================================
        if (Pipeline == 'DiffusionHCP'):
            #===============================================================================
            # /data/intradb/pipeline/bin/PipelineJobSubmitter /data/intradb/pipeline/bin/XnatPipelineLauncher -pipeline /data/intradb/pipeline/catalog/DiffusionHCP/DiffusionHCP.xml 
            # -id HCPIntradb_E07951 -host https://intradb.humanconnectome.org -u tony -pwd asdf -dataType xnat:mrSessionData -label 142828_diff -supressNotification -parameter project=HCP_Phase2 -notify wilsont@mir.wustl.edu 
            # -parameter xnat_id=HCPIntradb_E07951 -parameter sessionid=142828_diff -parameter subjects=142828 -parameter LR_Dir1=95 -parameter LR_Dir2=96 -parameter LR_Dir3=EMPTY -parameter RL_Dir1=95 -parameter RL_Dir2=96 -parameter RL_Dir3=EMPTY 
            # -parameter EchoSpacing=0.7800117313764398 -parameter PhaseEncodingDir=1 -parameter RL_1ScanId=14 -parameter RL_2ScanId=18 -parameter RL_3ScanId=20 -parameter LR_1ScanId=16 -parameter LR_2ScanId=20 -parameter LR_3ScanId=20 
            # -notify db-admin@humanconnectome.org -parameter mailhost=mail.nrg.wustl.edu -parameter userfullname=T.Wilson -parameter builddir=/data/intradb/build1/HCP_Phase2/1357674525 -parameter xnatserver=HCPIntradb 
            # -parameter adminemail=db-admin@humanconnectome.org -parameter useremail=wilsont@mir.wustl.edu -project HCP_Phase2 > /data/intradb/build1/HCP_Phase2/1357674525/142828/DiffusionHCPLaunchSTDOUT.txt
            #===============================================================================
            
            #===================================================================
            # grad a dummy scan id to feed to XML if scan does not exist.  XML must have scan id, else it will break...
            #===================================================================
            DummyScanId = sessionMeta.get('IDs')[0]
            EchoSpacing = '-parameter EchoSpacing=0.7800117313764398 '
            PhaseEncodingDir = '-parameter PhaseEncodingDir=1 '
            
            
            DiffusionSeriesList = ['DWI_RL_dir95','DWI_RL_dir96','DWI_RL_dir97','DWI_LR_dir95','DWI_LR_dir96','DWI_LR_dir97']
            
            DiffusionScanIdList = ['RL_1ScanId', 'RL_2ScanId', 'RL_3ScanId', 'LR_1ScanId', 'LR_2ScanId', 'LR_3ScanId']
            DiffusionScanIdDict = {'RL_1ScanId' : None, 'RL_2ScanId' : None, 'RL_3ScanId' : None, 'LR_1ScanId' : None, 'LR_2ScanId' : None, 'LR_3ScanId' : None}
            DiffusionDirList = ['RL_Dir1', 'RL_Dir2', 'RL_Dir3', 'LR_Dir1', 'LR_Dir2', 'LR_Dir3']
            DiffusionDirDict = {'RL_Dir1' : '95', 'RL_Dir2' : '96', 'RL_Dir3' : '97', 'LR_Dir1' : '95', 'LR_Dir2' : '96', 'LR_Dir3' : '97' }
            
#            DiffusionSeriesIntersectList = list(set(DiffusionSeriesList) & set(seriesList))

            #===============================================================================
            # HERE BE DRAGONS...
            #===============================================================================

            for j in xrange(0, len(DiffusionSeriesList)):
                currDiffDesc = DiffusionSeriesList[j]
                if (sessionMeta.get('Series').count(currDiffDesc) > 0):
                    currDiffIdx = sessionMeta.get('Series').index(currDiffDesc)
                    currScanId = sessionMeta.get('IDs')[currDiffIdx]
                    currQuality = sessionMeta.get('Quality')[currDiffIdx]
                    getHCP.Scan = currScanId
                    scanParms = getHCP.getScanParms()
                    scanMeta = getHCP.getScanMeta()
                    
                    # ScanIdDict['LR_2ScanId'] = '-parameter LR_2ScanId=%s ' % str(currScanId)
                    if (currQuality in UsableList):
                        DiffusionScanIdDict[DiffusionScanIdList[DiffusionSeriesList.index(currDiffDesc)]] = '-parameter %s=%s ' % (DiffusionScanIdList[DiffusionSeriesList.index(currDiffDesc)], currScanId)
                    else:
                        DiffusionScanIdDict[DiffusionScanIdList[DiffusionSeriesList.index(currDiffDesc)]] = '-parameter %s=%s ' % (DiffusionScanIdList[DiffusionSeriesList.index(currDiffDesc)], DummyScanId)
                        DiffusionDirDict[DiffusionDirList[DiffusionSeriesList.index(currDiffDesc)]] = 'EMPTY'
                    
                else:
                    DiffusionScanIdDict[DiffusionScanIdList[DiffusionSeriesList.index(currDiffDesc)]] = '-parameter %s=%s ' % (DiffusionScanIdList[DiffusionSeriesList.index(currDiffDesc)], DummyScanId)
                    DiffusionDirDict[DiffusionDirList[DiffusionSeriesList.index(currDiffDesc)]] = 'EMPTY'
            
            
            LR_Dir1 = '-parameter LR_Dir1=%s ' % DiffusionDirDict['LR_Dir1']
            LR_Dir2 = '-parameter LR_Dir2=%s ' % DiffusionDirDict['LR_Dir2']
            LR_Dir3 = '-parameter LR_Dir3=%s ' % DiffusionDirDict['LR_Dir3']
            RL_Dir1 = '-parameter RL_Dir1=%s ' % DiffusionDirDict['RL_Dir1']
            RL_Dir2 = '-parameter RL_Dir2=%s ' % DiffusionDirDict['RL_Dir2']
            RL_Dir3 = '-parameter RL_Dir3=%s ' % DiffusionDirDict['RL_Dir3']
            
            
            SubmitStr = JobSubmitter + PipelineLauncher + launcherPipeline + launcherHCPid + Host + launcherUser + launcherPassword + DataType + launcherLabel + SupressNotify + launcherProject + NotifyUser + launcherXnatId + launcherSession + \
            launcherSubject + LR_Dir1 + LR_Dir2 + LR_Dir3 + RL_Dir1 + RL_Dir2 + RL_Dir3 + EchoSpacing + PhaseEncodingDir + \
            DiffusionScanIdDict['RL_1ScanId'] + DiffusionScanIdDict['RL_2ScanId'] + DiffusionScanIdDict['RL_3ScanId'] + DiffusionScanIdDict['LR_1ScanId'] + DiffusionScanIdDict['LR_2ScanId'] + DiffusionScanIdDict['LR_3ScanId'] + \
            NotifyAdmin + MailHost + UserFullName + BuildDir + XnatServer + AdminEamil + UserEmail + RedirectionStr
            
            if sys.platform == 'win32':
                print 'Index of EMPTY ' + str(SubmitStr.find('EMPTY'))
#                    if (SubmitStr.find('EMPTY') == -1):
                print SubmitStr
#                    else: 
#                        print 'ERROR: EMPTY found for diffusion scan...'
            else:
                print SubmitStr
                os.system(SubmitStr)
                
        elif (Pipeline == 'FunctionalHCP'):
                
                
            # NOTE: not liking all the +/- 1 stuff here for ints v strings...
            if (seriesList.count(iFunctSeries) == 1):
                iFuncScanId = idList[seriesList.index(iFunctSeries)]
                funcQuality = qualityList[seriesList.index(iFunctSeries)]
            else:
                funcIndicesArray = fGetAllIndices(iFunctSeries, seriesList)
                funcQualityList = list()
                for j in xrange(0, len(funcIndicesArray)):
                    currQuality = qualityList[funcIndicesArray[j]]
                    if (currQuality == 'undetermined') or (currQuality == 'usable'):
                        iFuncScanId = funcIndicesArray[j] + 1
                        funcQuality = currQuality
                    else:
                        print 'WARNING: Functional scan ' +str(funcIndicesArray[j]+1)+ ' is neither undetermined or usable...'
                
            if (seriesList.count(iFunctSeries + '_SBRef') == 1):
                iScoutScanId = idList[seriesList.index(iFunctSeries + '_SBRef')]
                scoutQuality = qualityList[seriesList.index(iFunctSeries + '_SBRef')]
            else:
                
                scoutIndicesArray = fGetAllIndices(iFunctSeries + '_SBRef', seriesList)
                for j in xrange(0, len(scoutIndicesArray)):
                    seriesEnumerate = enumerate(seriesList)
                    currQuality = qualityList[scoutIndicesArray[j]]
                    if (currQuality == 'undetermined') or (currQuality == 'usable'):
                        iScoutScanId = scoutIndicesArray[j] + 1
                        iScoutQuality = currQuality
                    else:
                        print 'WARNING: Functional Scout scan ' +str(scoutIndicesArray[j]+1)+ ' is neither undetermined or usable...'
            
            funcScanAcqTime = fGetAcquisitionTime(iUser, iPassword, iProject, iSubject, iSessionId, str(iFuncScanId))
            
            magScanCount = seriesList.count('BOLD_LR_SB_SE')
            magScanIdList = list()
            phaScanIdList = list()
            magScanTimeList = list()
            magScanDiffList = list()
            currMagScanId = -1
            currPhaScanId = -1
            for j in xrange(0, magScanCount):
                currMagScanId = idList[seriesList.index('BOLD_LR_SB_SE', int(currMagScanId) + 1)]
                magScanIdList.append(currMagScanId)
                currPhaScanId = idList[seriesList.index('BOLD_RL_SB_SE', int(currPhaScanId) + 1)]
                phaScanIdList.append(currPhaScanId)
                magScanAcqTime = fGetAcquisitionTime(iUser, iPassword, iProject, iSubject, iSessionId, currMagScanId)
                magScanTimeList.append(magScanAcqTime)
                magScanDelta = datetime.datetime.strptime(funcScanAcqTime, '%H:%M:%S') - datetime.datetime.strptime(magScanAcqTime, '%H:%M:%S')
                magScanDiffList.append(magScanDelta.seconds)
        
            #///////////////////////////////////////////////////////////////
            scanParms = fGetScanParms( iUser, iPassword, iProject, iSubject, iSessionId, str(iFuncScanId) )
            shimCurrent = scanParms.get('alShimCurrent')
            linearOffset = scanParms.get('LinearOffset')
            #///////////////////////////////////////////////////////////////
            
            minIdx = magScanDiffList.index(min(magScanDiffList)) 
            iMagScanId = magScanIdList[minIdx]
            iPhaScanId = phaScanIdList[minIdx]
            #------------------------------------------
            MagScanId = '-parameter magscanid=' + iMagScanId + ' '
            PhaScanId = '-parameter phascanid=' + iPhaScanId + ' '
            #------------------------------------------
            FuncScanId = '-parameter functionalscanid=' + str(iFuncScanId) + ' '
            ScoutScanId = '-parameter scoutscanid=' + str(iScoutScanId) + ' '
            FunctSeries = '-parameter functionalseries=' + iFunctSeries + ' '
            LR_Fieldmap = '-parameter lr_fieldmapseries=BOLD_LR_SB_SE '
            RL_Fieldmap = '-parameter rl_fieldmapseries=BOLD_RL_SB_SE '
            DwellTime = '-parameter DwellTime=0.00058 '
            TE = '-parameter TE=2.46 '
            
            if (iFunctSeries.find('RL') != -1):
                iUnwarpDir = 'x'
            elif (iFunctSeries.find('LR') != -1):
                iUnwarpDir = 'x-'
            UnwarpDir = '-parameter UnwarpDir=' + iUnwarpDir + ' '
            DistortionCorrect = '-parameter DistortionCorrection=TOPUP '
            #-------------------------------------------
            
            SubmitStr = JobSubmitter + PipelineLauncher + Pipeline + HCPid + Host + User + Password + DataType + Label + SupressNotify + Project + NotifyUser + NotifyAdmin + \
            MailHost + UserFullName + BuildDir + XnatServer + AdminEamil + UserEmail + XnatId + SessionId + Subjects + MagScanId + PhaScanId + FuncScanId + ScoutScanId + \
            FunctSeries + LR_Fieldmap + RL_Fieldmap + DwellTime + TE + UnwarpDir + DistortionCorrect + RedirectionStr
            
            if sys.platform == 'win32':
                print SubmitStr
            else:
                print SubmitStr
                os.system(SubmitStr)
            

        elif (iPipeline == 'StructuralHCP'):
            
            iT1wSeriesDesc_1 = 'T1w_MPR1'
            iT1wSeriesDesc_2 = 'T1w_MPR2' 
            iT2wSeriesDesc_1 = 'T2w_SPC1'
            iT2wSeriesDesc_2 = 'T2w_SPC2'
            
            for j in xrange(0, len(seriesList)):
                currSeriesDesc = seriesList[j]
                currTypeList = typeList[j]
                if (currSeriesDesc.find(iT1wSeriesDesc_1) != -1): 
                    iT1wSeriesDesc_1 = currSeriesDesc
                    iT1wScanId_1 = idList[j]
                if (currSeriesDesc.find(iT1wSeriesDesc_2) != -1): 
                    iT1wSeriesDesc_2 = currSeriesDesc
                    iT1wScanId_2 = idList[j]
                if (currSeriesDesc.find(iT2wSeriesDesc_1) != -1): 
                    iT2wSeriesDesc_1 = currSeriesDesc
                    iT2wScanId_1 = idList[j]
                if (currSeriesDesc.find(iT2wSeriesDesc_2) != -1): 
                    iT2wSeriesDesc_2 = currSeriesDesc
                    iT2wScanId_2 = idList[j]
                    
                # grab the fieldmap ids
                if (currSeriesDesc.find('FieldMap_Magnitude') != -1) and (currTypeList.find('FieldMap') != -1): 
                    iMagScanId = idList[j]
                if (currSeriesDesc.find('FieldMap_Phase') != -1) and (currTypeList.find('FieldMap') != -1): 
                    iPhaScanId = idList[j]
            
            #===========================================================================
            # NOW HANDLES TERNARY CASES...
            #===========================================================================
            if (seriesList.count(iT1wSeriesDesc_1) == 1):
                if (qualityList[seriesList.index(iT1wSeriesDesc_1)] == 'usable') or (qualityList[seriesList.index(iT1wSeriesDesc_1)] == 'excellent') or (qualityList[seriesList.index(iT1wSeriesDesc_1)] == 'good'):
                    iT1wScanId_1 = idList[seriesList.index(iT1wSeriesDesc_1)]
                else:
                    iT1wScanId_1 = '00'
                    iT1wSeriesDesc_1 = 'XXX'
                    curriT1wSeriesIdx = -1
                    for j in xrange(0, seriesList.count(iT1wSeriesDesc_2)):
                        curriT1wSeriesIdx = seriesList.index(iT1wSeriesDesc_2, curriT1wSeriesIdx + 1)
                        if (qualityList[curriT1wSeriesIdx] == 'usable') or (qualityList[curriT1wSeriesIdx] == 'excellent') or (qualityList[curriT1wSeriesIdx] == 'good'):
                            iT1wScanId_1 = idList[curriT1wSeriesIdx]
                            iT1wSeriesDesc_1 = iT1wSeriesDesc_2
                    if (qualityList[seriesList.index(iT1wSeriesDesc_2)] == 'usable') or (qualityList[seriesList.index(iT1wSeriesDesc_2)] == 'excellent') or (qualityList[seriesList.index(iT1wSeriesDesc_2)] == 'good'):
                        iT1wScanId_1 = idList[seriesList.index(iT1wSeriesDesc_2)]
                        iT1wSeriesDesc_1 = iT1wSeriesDesc_2
            if (seriesList.count(iT1wSeriesDesc_2) == 1):
                if (qualityList[seriesList.index(iT1wSeriesDesc_2)] == 'usable') or (qualityList[seriesList.index(iT1wSeriesDesc_2)] == 'excellent') or (qualityList[seriesList.index(iT1wSeriesDesc_2)] == 'good'):
                    iT1wScanId_2 = idList[seriesList.index(iT1wSeriesDesc_2)]
                else:
                    iT1wScanId_2 = '00'
                    iT1wSeriesDesc_2 = 'XXX'
                    if (qualityList[seriesList.index(iT1wSeriesDesc_1)] == 'usable') or (qualityList[seriesList.index(iT1wSeriesDesc_1)] == 'excellent') or (qualityList[seriesList.index(iT1wSeriesDesc_1)] == 'good'):
                        iT1wScanId_2 = idList[seriesList.index(iT1wSeriesDesc_1)]
                        iT1wSeriesDesc_2 = iT1wSeriesDesc_1
            if (seriesList.count(iT2wSeriesDesc_1) == 1):
                if (qualityList[seriesList.index(iT2wSeriesDesc_1)] == 'usable') or (qualityList[seriesList.index(iT2wSeriesDesc_1)] == 'excellent') or (qualityList[seriesList.index(iT2wSeriesDesc_1)] == 'good'):
                    iT2wScanId_1 = idList[seriesList.index(iT2wSeriesDesc_1)]
                else:
                    iT2wScanId_1 = '00'
                    iT2wSeriesDesc_1 = 'XXX'
                    curriT2wSeriesIdx = -1
                    for j in xrange(0, seriesList.count(iT2wSeriesDesc_2)):
                        curriT2wSeriesIdx = seriesList.index(iT2wSeriesDesc_2, curriT2wSeriesIdx + 1)
                        if (qualityList[curriT2wSeriesIdx] == 'usable') or (qualityList[curriT2wSeriesIdx] == 'excellent') or (qualityList[curriT2wSeriesIdx] == 'good'):
                            iT2wScanId_1 = idList[curriT2wSeriesIdx]
                            iT2wSeriesDesc_1 = iT2wSeriesDesc_2
            if (seriesList.count(iT2wSeriesDesc_2) == 1):
                if (qualityList[seriesList.index(iT2wSeriesDesc_2)] == 'usable') or (qualityList[seriesList.index(iT2wSeriesDesc_2)] == 'excellent') or (qualityList[seriesList.index(iT2wSeriesDesc_2)] == 'good'):
                    iT2wScanId_2 = idList[seriesList.index(iT2wSeriesDesc_2)]
                else:
                    iT2wScanId_2 = '00'
                    iT2wSeriesDesc_2 = 'XXX'
                    if (qualityList[seriesList.index(iT2wSeriesDesc_1)] == 'usable') or (qualityList[seriesList.index(iT2wSeriesDesc_1)] == 'excellent') or (qualityList[seriesList.index(iT2wSeriesDesc_1)] == 'good'):
                        iT2wScanId_2 = idList[seriesList.index(iT2wSeriesDesc_1)]
                        iT2wSeriesDesc_2 = iT2wSeriesDesc_1
                
        
            MagScanId = '-parameter magscanid=%s ' % (iMagScanId)
            PhaScanId = '-parameter phascanid=%s ' % (iPhaScanId)
            
            T1wScanId_1 = '-parameter t1scanid_1=%s ' % (iT1wScanId_1)
            T1wScanId_2 = '-parameter t1scanid_2=%s ' % (iT1wScanId_2)
            T2wScanId_1 = '-parameter t2scanid_1=%s ' % (iT2wScanId_1)
            T2wScanId_2 = '-parameter t2scanid_2=%s ' % (iT2wScanId_2)
            
            T1wSeriesDesc_1 = '-parameter t1seriesdesc_1=%s ' % (iT1wSeriesDesc_1)
            T1wSeriesDesc_2 = '-parameter t1seriesdesc_2=%s ' % (iT1wSeriesDesc_2)
            T2wSeriesDesc_1 = '-parameter t2seriesdesc_1=%s ' % (iT2wSeriesDesc_1)
            T2wSeriesDesc_2 = '-parameter t2seriesdesc_2=%s ' % (iT2wSeriesDesc_2)
        
            
            TE = '-parameter TE=2.46 '
            sampleSpacingT1w = fGetScanParms( iUser, iPassword, iProject, iSubject, iSessionId, iT1wScanId_1 ).get('SampleSpacing')
            sampleSpacingT2w = fGetScanParms( iUser, iPassword, iProject, iSubject, iSessionId, iT2wScanId_1 ).get('SampleSpacing')
            T1wSampleSpacing = "-parameter T1wSampleSpacing=%1.9f " % (float(sampleSpacingT1w)/1.0e+9)
            T2wSampleSpacing = "-parameter T2wSampleSpacing=%1.9f " % (float(sampleSpacingT2w)/1.0e+9)
            T1wTemplate = '-parameter T1wTemplate=MNI152_T1_0.7mm.nii.gz '
            T1wTemplateBrain = '-parameter T1wTemplateBrain=MNI152_T1_0.7mm_brain.nii.gz '
            T2wTemplate = '-parameter T2wTemplate=MNI152_T2_0.7mm.nii.gz '
            T2wTemplateBrain = '-parameter T2wTemplateBrain=MNI152_T2_0.7mm_brain.nii.gz '
            TemplateMask = '-parameter TemplateMask=MNI152_T1_0.7mm_brain_mask.nii.gz '
            
            # for PostFS...
            FinalTemplateSpace = '-parameter FinalTemplateSpace=MNI152_T1_0.7mm.nii.gz'
            
            SubmitStr = JobSubmitter + PipelineLauncher + Pipeline + HCPid + Host + User + Password + DataType + Label + SupressNotify + Project + NotifyUser + NotifyAdmin + \
            MailHost + UserFullName + BuildDir + XnatServer + AdminEamil + UserEmail + XnatId + SessionId + Subjects + MagScanId + PhaScanId + T1wScanId_1 + T1wScanId_2 + \
            T2wScanId_1 + T2wScanId_2 + T1wSeriesDesc_1 + T1wSeriesDesc_2 + T2wSeriesDesc_1 + T2wSeriesDesc_2 + TE + T1wSampleSpacing + T2wSampleSpacing + T1wTemplate + \
            T1wTemplateBrain + T2wTemplate + T2wTemplateBrain + TemplateMask + FinalTemplateSpace + RedirectionStr
            
            if sys.platform == 'win32':
                print SubmitStr
            else:
                print SubmitStr
                os.system(SubmitStr)
            
        
            
            
            
        if (linIdx < ( len(SubjectsList) * len(SeriesList) )):
            print 'Sleeping for ' + str(SleepTime) + ' seconds...'         
            time.sleep(SleepTime)
        else:
            print 'Done...total launch time was %s seconds for %s jobs with a sleep time of %s seconds per job...' % ( (time.time() - sTime), ( len(SubjectsList) * len(SeriesList) ), str(SleepTime) ) 
    
if __name__ == '__main__':
    pass
