'''
Created on Sep 3, 2012

@author: jwilso01
'''
import os
import sys
import time
import numpy
#import array
import socket
import argparse
import datetime

from pyHCP import pyHCP, getHCP, writeHCP


sTime = time.time()

#===============================================================================
# PARSE INPUT
#===============================================================================
# Examples:
# -User tony -Password passfoo -Server hcpx-demo.humanconnectome.org -Pipeline DiffusionHCP -Subject 100307 -Project ReleaseTest 
# -User tony -Password passfoo -Pipeline FunctionalHCP -Server hcpx-dev-cuda00.nrg.mir -Subjects 100307 -Project HCP_Q1
# -User tony -Password passfoo -Pipeline FunctionalHCP -Server https://hcpx-demo.humanconnectome.org -Subject 100307 -Project ReleaseTest
# -User tony -Password passfoo -Server hcpi-dev-cuda00.nrg.mir -Pipeline StructuralHCP -Subject 100307 -Project HCP_Phase2
#===============================================================================
parser = argparse.ArgumentParser(description="Script to generate proper command for XNAT functional pipeline lauching ...")

#MANDATORY....
parser.add_argument("-User", "--User", dest="User", default='tony', type=str)
parser.add_argument("-Password", "--Password", dest="Password", default='none', type=str)
parser.add_argument("-Pipeline", "--Pipeline", dest="Pipeline", default='fMRIVolume', type=str)
parser.add_argument("-Subjects", "--Subjects", dest="Subjects", default='00', type=str)
parser.add_argument("-Server", "--Server", dest="Server", default='http://hcpi-dev-cuda00.nrg.mir/', type=str)
parser.add_argument('--version', action='version', version='%(prog)s 2.0.0')
# IF Pipeline == Functional...
parser.add_argument("-FunctSeries", "--FunctSeries", dest="FunctSeries", default=None, type=str)
#END MANDATORY....
parser.add_argument("-Project", "--Project", dest="Project", default='HCP_Phase2', type=str)
parser.add_argument("-Shadow", "--Shadow", dest="Shadow", default=None, type=str)
parser.add_argument("-Build", "--Build", dest="Build", default=None, type=str)
parser.add_argument("-Compute", "--Compute", dest="Compute", default='NRG', type=str)

parser.add_argument("-Production", "--Production", dest="iProduction", default=False, type=bool)


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
Compute = args.Compute

if (Server.find('intra') != -1) or (Server.find('hcpi') != -1):
    dbStr = 'intradb'
    if (Pipeline == 'FunctionalHCP'): 
        print 'Series Descriptions Changed: IntraDB FunctionalHCP broken...'
        sys.exit()
        
elif (Server.find('hcpx') != -1) or (Server.find('db.humanconnectome') != -1):
    dbStr = 'hcpdb'

#===============================================================================
# STATIC PARMS...
#===============================================================================
JobSubmitter = '/data/%s/pipeline/bin/PipelineJobSubmitter ' % (dbStr)
PipelineLauncher = '/data/%s/pipeline/bin/XnatPipelineLauncher ' % (dbStr)
DataType = '-dataType xnat:mrSessionData '
SupressNotify = '-supressNotification '
NotifyUser = '-notify wilsont@mir.wustl.edu ' 
NotifyAdmin = '-notify db-admin@humanconnectome.org '
MailHost = '-parameter mailhost=mail.nrg.wustl.edu '
UserFullName = '-parameter userfullname=T.Wilson '
if (dbStr == 'intradb'): 
    XnatServer = '-parameter xnatserver=HCPIntradb '
else:
    XnatServer = '-parameter xnatserver=ConnectomeDB '
AdminEmail = '-parameter adminemail=db-admin@humanconnectome.org '
UserEmail = '-parameter useremail=wilsont@mir.wustl.edu '
if (Compute == 'NRG'):
    TemplatesDir = '/nrgpackages/atlas/HCP/'
    ConfigDir = '/nrgpackages/tools/HCP/conf/'
    CaretAtlasDir = '/nrgpackages/atlas/HCP/standard_mesh_atlases/'
elif (Compute == 'CHPC'):
    TempaltesDir = '/NRG/BlueArc/nrgpackages/atlas/HCP/'
    ConfigDir = '/NRG/BlueArc/nrgpackages/tools/HCP/conf/'
    CaretAtlasDir = '/NRG/BlueArc/nrgpackages/atlas/HCP/standard_mesh_atlases/'
     
#===============================================================================
# HACK for REST RL/LR...
#===============================================================================
preChangeRLList = ('100307','111312','114924','119833','125525','138231','144266','150423','159239','162329',\
                   '167743','174437','185139','192439','197550','199251','217429','249947','255639','329440',\
                   '355542','499566','585862','611231','665254','672756','792564','826353','877168','896778')
#===============================================================================
# INTERFACE...
#===============================================================================
pyHCP = pyHCP(User, Password, Server)
getHCP = getHCP(pyHCP)
getHCP.Project = Project
#===============================================================================
# FUNCTIONS...
#===============================================================================
#def fGetAllIndices(inputVal, inputList):
#    allIndicesList = list()
#    currIdx = -1
#    while True:
#        try:
#            currIdx = inputList.index(inputVal, currIdx + 1)
#            allIndicesList.append(currIdx)
#        except ValueError:
#            break
#    return numpy.asarray(allIndicesList, dtype=numpy.int)
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
    

if (len(SubjectsList) > len(ShadowList)):
#    ShadowArray = numpy.tile(ShadowList, (numpy.ceil(len(SubjectsList) / len(ShadowList))))
    ShadowArray = numpy.tile(ShadowList, (numpy.ceil(len(SubjectsList))))
else: 
    ShadowArray = numpy.tile(ShadowList, (1))

if (len(SubjectsList) > len(BuildList)):
    BuildArray = numpy.tile(BuildList, (numpy.ceil(len(SubjectsList))))
else: 
    BuildArray = numpy.tile(BuildList, 1)
    
if (Pipeline == 'FunctionalHCP'): PipelineSubString = ['fnc', 'task', 'rest']
elif (Pipeline == 'StructuralHCP'): PipelineSubString = ['strc']
elif (Pipeline == 'DiffusionHCP'): PipelineSubString = ['diff']
elif (Pipeline == 'FIX_HCP'): PipelineSubString = ['fnc', 'rest']
else: PipelineSubString = None
    
UsableList = ['good', 'excellent', 'usable', 'undetermined']



linIdx = 0
for h in xrange(0, len(SubjectsList)): 
    getHCP.Subject = SubjectsList[h]
    SubjectSessions = getHCP.getSubjectSessions()
    
    # find correct session...
    for i in xrange(0, len(SubjectSessions.get('Sessions'))):
        getHCP.Session = SubjectSessions.get('Sessions')[i]
        sessionMeta = getHCP.getSessionMeta()
        if (FunctSeries != None) and (FunctSeries in sessionMeta.get('Series')):
            break
        else: 
            getHCP.Session = SubjectSessions.get('Sessions')[SubjectSessions.get('Types').index(PipelineSubString[0])]
            
    sessionMeta = getHCP.getSessionMeta()
    seriesList = sessionMeta.get('Series')
    typeList = sessionMeta.get('Types')
    idList = sessionMeta.get('IDs')
    qualityList = sessionMeta.get('Quality')
    
    if (FunctSeries == None) and (Pipeline == 'FunctionalHCP'):
        FunctionalList = list()
        for i in xrange(0, len(sessionMeta.get('Types'))):
            if (sessionMeta.get('Types')[i] == 'tfMRI') or (sessionMeta.get('Types')[i] == 'rfMRI'):
                FunctionalList.append(sessionMeta.get('Series')[i])
    elif (Pipeline == 'FunctionalHCP') or (Pipeline == 'FIX_HCP'):
        FunctionalList = FunctSeries.split(',')
    else:
        FunctionalList = ['Other']
        
    
#    print set(PipelineSubString) & set(SubjectSessions.get('Types'))
    if (len(set(PipelineSubString) & set(SubjectSessions.get('Types'))) == 0):
        print 'ERROR: No ' +Pipeline+ ' session could be found for subject ' +getHCP.Subject
        
    
    for i in xrange(0, len(FunctionalList)):
        linIdx += 1
        

            
        if (Pipeline.find('Funct') != -1) or (Pipeline.find('FIX') != -1):
            currSeries = FunctionalList[i]
            
            if (currSeries in sessionMeta.get('Series')):
                currUsability = sessionMeta.get('Quality')[sessionMeta.get('Series').index(currSeries)]
            else:
                print 'Current usability could not be determined for subject %s, session %s...' % (getHCP.Subject, getHCP.Session)

        elif (Pipeline.find('Diff') != -1):
            diffSessionIdx = 0
#            print SubjectSessions.get('Types').count('diff')
            for j in xrange(0, SubjectSessions.get('Types').count('diff')):
                diffSessionIdx =+ SubjectSessions.get('Types').index('diff')
                getHCP.Session = SubjectSessions.get('Sessions')[diffSessionIdx]
                sessionMeta = getHCP.getSessionMeta()
                
        elif (Pipeline.find('Struct') != -1):
#            print 'Looking at StructuralHCP...'
            structSessionIdx = SubjectSessions.get('Types').index('strc')
            getHCP.Session = SubjectSessions.get('Sessions')[structSessionIdx]
            sessionMeta = getHCP.getSessionMeta()
            
        else:
            print 'ERROR: Pipline not found...'
            
        
        
        if (BuildArray.size >= len(SubjectsList)):
            BuildDirRoot = '/data/%s/build%s/%s/' % (dbStr, str(BuildArray[h]), Project)
        else:
            BuildDirRoot = '/data/%s/build/%s/' % (dbStr, Project)
#            print iBuildDirRoot
        
        launcherProject = '-parameter project=%s ' % Project 
        if (Pipeline == 'FIX_HCP'):
            launcherPipeline = '-pipeline /data/%s/pipeline/catalog/FIX/%s.xml ' % (dbStr, Pipeline)
        else:
            launcherPipeline = '-pipeline /data/%s/pipeline/catalog/%s/%s.xml ' % (dbStr, Pipeline, Pipeline)
        launcherUser = '-u %s ' % User 
        launcherPassword = '-pwd %s ' % Password 
        launcherLabel = '-label %s ' % getHCP.Session
        launcherHCPid = '-id %s ' % sessionMeta.get('XNATID')[0]
        launcherXnatId = '-parameter xnat_id=%s ' % sessionMeta.get('XNATID')[0] 
        launcherSession = '-parameter sessionid=%s ' % getHCP.Session 
        launcherSubject = '-parameter subjects=%s ' % getHCP.Subject
        
        launcherTemplatesDir = '-parameter templatesdir=%s ' % TemplatesDir
        launcherConfigDir = '-parameter configdir=%s ' % ConfigDir
        launcherCaretAtlasDir = '-parameter CaretAtlasDir=%s ' % CaretAtlasDir

        
        if (len(FunctionalList) > 1):
            currBuildDir = BuildDirRoot + str(numpy.asarray(round(time.time()), dtype=numpy.uint64))
        else:
            currBuildDir = BuildDirRoot + str(numpy.asarray(round(sTime), dtype=numpy.uint64))
        BuildDir = '-parameter builddir=%s ' % (currBuildDir)

        
#        if not os.path.exists(currBuildDir + os.sep + getHCP.Subject) and sys.platform != 'win32':
#            os.makedirs(currBuildDir + os.sep + getHCP.Subject)
            
        RedirectionStr = ' > ' + currBuildDir.replace(' ', '') + os.sep + getHCP.Subject + os.sep + Pipeline + 'LaunchSTDOUT.txt'
        
        if (socket.gethostname() == 'intradb.humanconnectome.org') and (Shadow != None):
            Host = '-host https://intradb-shadow' + ShadowArray[h] + '.nrg.mir '
        else: 
            Host = '-host https://%s ' % (socket.gethostname())

                
            
        #===============================================================================
        # DiffusionHCP....
        #===============================================================================
        if (Pipeline == 'DiffusionHCP'):
            #===================================================================
            # grab a dummy scan id to feed to XML if scan does not exist.  XML must have scan id, else it will break...
            #===================================================================
            DummyScanId = sessionMeta.get('IDs')[0]
            EchoSpacingList = list()
            PhaseEncodingDirList = list() 

            
            # if intradb...
            if (dbStr == 'intradb'):
                DiffusionSeriesList = ['DWI_RL_dir95','DWI_RL_dir96','DWI_RL_dir97','DWI_LR_dir95','DWI_LR_dir96','DWI_LR_dir97']
            elif (dbStr == 'hcpdb'):
                DiffusionSeriesList = ['DWI_dir95_RL','DWI_dir96_RL','DWI_dir97_RL','DWI_dir95_LR','DWI_dir96_LR','DWI_dir97_LR']
            
            DiffusionScanIdList = ['RL_1ScanId', 'RL_2ScanId', 'RL_3ScanId', 'LR_1ScanId', 'LR_2ScanId', 'LR_3ScanId']
            DiffusionScanIdDict = {'RL_1ScanId' : None, 'RL_2ScanId' : None, 'RL_3ScanId' : None, 'LR_1ScanId' : None, 'LR_2ScanId' : None, 'LR_3ScanId' : None}
            DiffusionDirList = ['RL_Dir1', 'RL_Dir2', 'RL_Dir3', 'LR_Dir1', 'LR_Dir2', 'LR_Dir3']
            DiffusionDirDict = {'RL_Dir1' : '95', 'RL_Dir2' : '96', 'RL_Dir3' : '97', 'LR_Dir1' : '95', 'LR_Dir2' : '96', 'LR_Dir3' : '97' }
            
#            DiffusionSeriesIntersectList = list(set(DiffusionSeriesList) & set(SeriesList))


            for j in xrange(0, len(DiffusionSeriesList)):
                currDiffDesc = DiffusionSeriesList[j]
                if (sessionMeta.get('Series').count(currDiffDesc) > 0):
                    currDiffIdx = sessionMeta.get('Series').index(currDiffDesc)
                    currScanId = sessionMeta.get('IDs')[currDiffIdx]
                    currQuality = sessionMeta.get('Quality')[currDiffIdx]
                    getHCP.Scan = currScanId
                    scanParms = getHCP.getScanParms()
                    scanMeta = getHCP.getScanMeta()
                    
                    EchoSpacingList.append(float(scanParms.get('EchoSpacing')) * 1.0e+3)
                    PhaseEncodingDirList.append(scanParms.get('PhaseEncodingDir'))
                    
                    # ScanIdDict['LR_2ScanId'] = '-parameter LR_2ScanId=%s ' % str(currScanId)
                    if (currQuality in UsableList):
                        DiffusionScanIdDict[DiffusionScanIdList[DiffusionSeriesList.index(currDiffDesc)]] = '-parameter %s=%s ' % (DiffusionScanIdList[DiffusionSeriesList.index(currDiffDesc)], currScanId)
                    else:
                        DiffusionScanIdDict[DiffusionScanIdList[DiffusionSeriesList.index(currDiffDesc)]] = '-parameter %s=%s ' % (DiffusionScanIdList[DiffusionSeriesList.index(currDiffDesc)], DummyScanId)
                        DiffusionDirDict[DiffusionDirList[DiffusionSeriesList.index(currDiffDesc)]] = 'EMPTY'
                    
                else:
                    DiffusionScanIdDict[DiffusionScanIdList[DiffusionSeriesList.index(currDiffDesc)]] = '-parameter %s=%s ' % (DiffusionScanIdList[DiffusionSeriesList.index(currDiffDesc)], DummyScanId)
                    DiffusionDirDict[DiffusionDirList[DiffusionSeriesList.index(currDiffDesc)]] = 'EMPTY'
            
            launcherEchoSpacing = '-parameter EchoSpacing=%s ' % (sum(EchoSpacingList) / float(len(EchoSpacingList)))
            # 1 for RL/LR phase encoding and 2 for AP/PA phase encoding
            launcherPhaseEncodingDir = '-parameter PhaseEncodingDir=1 '
            
            launcherLR_Dir1 = '-parameter LR_Dir1=%s ' % DiffusionDirDict['LR_Dir1']
            launcherLR_Dir2 = '-parameter LR_Dir2=%s ' % DiffusionDirDict['LR_Dir2']
            launcherLR_Dir3 = '-parameter LR_Dir3=%s ' % DiffusionDirDict['LR_Dir3']
            launcherRL_Dir1 = '-parameter RL_Dir1=%s ' % DiffusionDirDict['RL_Dir1']
            launcherRL_Dir2 = '-parameter RL_Dir2=%s ' % DiffusionDirDict['RL_Dir2']
            launcherRL_Dir3 = '-parameter RL_Dir3=%s ' % DiffusionDirDict['RL_Dir3']
            
            SubmitStr = JobSubmitter + PipelineLauncher + launcherPipeline + launcherHCPid + DataType + Host + XnatServer + launcherProject + launcherXnatId + launcherSession + launcherLabel + launcherUser + launcherPassword +  SupressNotify + NotifyUser + NotifyAdmin + AdminEmail + UserEmail + MailHost + UserFullName +\
            launcherEchoSpacing + launcherPhaseEncodingDir + launcherSubject + BuildDir + launcherLR_Dir1 + launcherLR_Dir2 + launcherLR_Dir3 + launcherRL_Dir1 + launcherRL_Dir2 + launcherRL_Dir3 + \
            DiffusionScanIdDict['RL_1ScanId'] + DiffusionScanIdDict['RL_2ScanId'] + DiffusionScanIdDict['RL_3ScanId'] + DiffusionScanIdDict['LR_1ScanId'] + DiffusionScanIdDict['LR_2ScanId'] + DiffusionScanIdDict['LR_3ScanId'] + RedirectionStr
            
            if sys.platform == 'win32':
                print 'Index of EMPTY ' + str(SubmitStr.find('EMPTY'))
                print SubmitStr
            else:
                print SubmitStr
                os.system(SubmitStr)
                
        #=======================================================================
        # StructuralHCP
        #=======================================================================
        elif (Pipeline == 'StructuralHCP'):
            
            PathMatch = list()
            ScanIdList = list()
            StructResources = ['T1w_MPR1_unproc', 'T1w_MPR2_unproc', 'T2w_SPC1_unproc', 'T2w_SPC2_unproc']
            getHCP.Resource = StructResources[0]
            resourceMeta = getHCP.getSubjectResourceMeta()
            

            StructuralSeriesDescDict = {'T1w_MPR1' : 'T1w_MPR1', 'T1w_MPR2' : 'T1w_MPR2', 'T2w_SPC1' : 'T2w_SPC1', 'T2w_SPC2' : 'T2w_SPC2'}
            StructuralSeriesDescScanIdDict = {'T1w_MPR1' : None, 'T1w_MPR2' : None, 'T2w_SPC1' : None, 'T2w_SPC2' : None}
            StructuralSeriesQualityDict = {'T1w_MPR1' : None, 'T1w_MPR2' : None, 'T2w_SPC1' : None, 'T2w_SPC2' : None}
            StructuralSeriesList = ['T1w_MPR1', 'T1w_MPR2', 'T2w_SPC1', 'T2w_SPC2']
            

            # grab the fieldmap ids...
            fieldmapMagIdx = seriesList.index('FieldMap_Magnitude')
            filedmapPhaIdx = seriesList.index('FieldMap_Phase')
            if (typeList[fieldmapMagIdx] == 'FieldMap') and (qualityList[fieldmapMagIdx] in UsableList): 
                MagScanId = idList[fieldmapMagIdx]
                getHCP.Scan = MagScanId
                magScanParms = getHCP.getScanParms()
                
                
            if (typeList[filedmapPhaIdx] == 'FieldMap') and (qualityList[filedmapPhaIdx] in UsableList):
                PhaScanId = idList[filedmapPhaIdx]
                getHCP.Scan = PhaScanId
                phaScanParms = getHCP.getScanParms()
                
                    
            # collect quality, series descriptions, and scan ids...
            for j in xrange(0, len(seriesList)):
                currSeriesDesc = seriesList[j]
                currTypeList = typeList[j]
                currQuality = qualityList[j]
                if (currSeriesDesc in StructuralSeriesList) and (qualityList[j] in UsableList): 
                    StructuralSeriesDescScanIdDict[currSeriesDesc] = idList[j]
                    StructuralSeriesQualityDict[currSeriesDesc] = qualityList[j]
                    

                
            # this should check for absence and quality of scans and swap if bad or absent...
            for j in xrange(0, len(StructuralSeriesList)):
                currSeries = StructuralSeriesList[j]
                if (StructuralSeriesDescScanIdDict.get(currSeries) == None):
                    if (currSeries == 'T1w_MPR1'):
                        if (StructuralSeriesDescScanIdDict.get('T1w_MPR2') != None): # or (StructuralSeriesQualityDict.get('T1w_MPR2') not in UsableList):
                            StructuralSeriesDescScanIdDict['T1w_MPR1'] = StructuralSeriesDescScanIdDict.get('T1w_MPR2')
                            StructuralSeriesDescDict['T1w_MPR1'] = 'T1w_MPR2'
                    elif (currSeries == 'T1w_MPR2'):
                        if (StructuralSeriesDescScanIdDict.get('T1w_MPR1') != None): # or (StructuralSeriesQualityDict.get('T1w_MPR1') not in UsableList):
                            StructuralSeriesDescScanIdDict['T1w_MPR2'] = StructuralSeriesDescScanIdDict.get('T1w_MPR1')
                            StructuralSeriesDescDict['T1w_MPR2'] = 'T1w_MPR1'
                    elif (currSeries == 'T2w_SPC1'):
                        if (StructuralSeriesDescScanIdDict.get('T2w_SPC2') != None): # or (StructuralSeriesQualityDict.get('T2w_SPC2') not in UsableList):
                            StructuralSeriesDescScanIdDict['T2w_SPC1'] = StructuralSeriesDescScanIdDict.get('T2w_SPC2')
                            StructuralSeriesDescDict['T2w_SPC1'] = 'T2w_SPC2'
                    elif (currSeries == 'T2w_SPC2'):
                        if (StructuralSeriesDescScanIdDict.get('T2w_SPC1') != None): # or (StructuralSeriesQualityDict.get('T2w_SPC1') not in UsableList):
                            StructuralSeriesDescScanIdDict['T2w_SPC2'] = StructuralSeriesDescScanIdDict.get('T2w_SPC1')
                            StructuralSeriesDescDict['T2w_SPC2'] = 'T2w_SPC1'
                            
            launcherMagScanId = '-parameter magscanid=%s ' % (MagScanId)
            launcherPhaScanId = '-parameter phascanid=%s ' % (PhaScanId)
            
            launcherT1wScanId_1 = '-parameter t1scanid_1=%s ' % StructuralSeriesDescScanIdDict.get('T1w_MPR1')
            launcherT1wScanId_2 = '-parameter t1scanid_2=%s ' % StructuralSeriesDescScanIdDict.get('T1w_MPR2')
            launcherT2wScanId_1 = '-parameter t2scanid_1=%s ' % StructuralSeriesDescScanIdDict.get('T2w_SPC1')
            launcherT2wScanId_2 = '-parameter t2scanid_2=%s ' % StructuralSeriesDescScanIdDict.get('T2w_SPC2')
            
            launcherT1wSeriesDesc_1 = '-parameter t1seriesdesc_1=%s ' % StructuralSeriesDescDict.get('T1w_MPR1')
            launcherT1wSeriesDesc_2 = '-parameter t1seriesdesc_2=%s ' % StructuralSeriesDescDict.get('T1w_MPR2')
            launcherT2wSeriesDesc_1 = '-parameter t2seriesdesc_1=%s ' % StructuralSeriesDescDict.get('T2w_SPC1')
            launcherT2wSeriesDesc_2 = '-parameter t2seriesdesc_2=%s ' % StructuralSeriesDescDict.get('T2w_SPC2')
        
            # Collect scan ids for later testing...
            ScanIdList.append(MagScanId)
            ScanIdList.append(PhaScanId)
            ScanIdList.append(StructuralSeriesDescScanIdDict.get('T1w_MPR1'))
            ScanIdList.append(StructuralSeriesDescScanIdDict.get('T1w_MPR2'))
            ScanIdList.append(StructuralSeriesDescScanIdDict.get('T2w_SPC1'))
            ScanIdList.append(StructuralSeriesDescScanIdDict.get('T2w_SPC2'))
            ScanIdList = list(set(ScanIdList))
            
            getHCP.Scan = StructuralSeriesDescScanIdDict.get('T1w_MPR1')
            scanParms = getHCP.getScanParms( )
            sampleSpacingT1w = scanParms.get('SampleSpacing')
            
            getHCP.Scan = StructuralSeriesDescScanIdDict.get('T2w_SPC1')
            sampleSpacingT2w = getHCP.getScanParms( ).get('SampleSpacing')
            
            TE = magScanParms.get('DeltaTE')
            launcherTE = '-parameter TE=%s ' % (TE)
            launcherT1wSampleSpacing = '-parameter T1wSampleSpacing=%1.9f ' % (float(sampleSpacingT1w)/1.0e+9)
            launcherT2wSampleSpacing = '-parameter T2wSampleSpacing=%1.9f ' % (float(sampleSpacingT2w)/1.0e+9)
            
            launcherT1wTemplate = '-parameter T1wTemplate=MNI152_T1_0.7mm.nii.gz '
            launcherT1wTemplateBrain = '-parameter T1wTemplateBrain=MNI152_T1_0.7mm_brain.nii.gz '
            launcherT2wTemplate = '-parameter T2wTemplate=MNI152_T2_0.7mm.nii.gz '
            launcherT2wTemplateBrain = '-parameter T2wTemplateBrain=MNI152_T2_0.7mm_brain.nii.gz '
            launcherTemplateMask = '-parameter TemplateMask=MNI152_T1_0.7mm_brain_mask.nii.gz '
            
            # for PostFS...
            launcherFinalTemplateSpace = '-parameter FinalTemplateSpace=MNI152_T1_0.7mm.nii.gz '
            
            SubmitStr = JobSubmitter + PipelineLauncher + launcherPipeline + launcherHCPid + DataType + Host + XnatServer + launcherProject + launcherXnatId + launcherLabel + launcherUser + launcherPassword +  SupressNotify + NotifyUser + NotifyAdmin + AdminEmail + UserEmail + MailHost + UserFullName +\
            BuildDir + launcherSession + launcherSubject + launcherMagScanId + launcherPhaScanId + launcherT1wScanId_1 + launcherT1wScanId_2 + \
            launcherT2wScanId_1 + launcherT2wScanId_2 + launcherT1wSeriesDesc_1 + launcherT1wSeriesDesc_2 + launcherT2wSeriesDesc_1 + launcherT2wSeriesDesc_2 + launcherTE + launcherT1wSampleSpacing + launcherT2wSampleSpacing + launcherT1wTemplate + \
            launcherT1wTemplateBrain + launcherT2wTemplate + launcherT2wTemplateBrain + launcherTemplateMask + launcherFinalTemplateSpace + launcherTemplatesDir + launcherConfigDir + launcherCaretAtlasDir + RedirectionStr
            
            # print scanParms.get('GEFieldMapGroup'), magScanParms.get('GEFieldMapGroup'), phaScanParms.get('GEFieldMapGroup')
            if (scanParms.get('GEFieldMapGroup') == magScanParms.get('GEFieldMapGroup') == phaScanParms.get('GEFieldMapGroup')):
                if sys.platform == 'win32':
                    print SubmitStr
                else:
                    # do T1w and T2w path test...
                    for j in xrange(0, len(StructResources)):
                        getHCP.Resource = StructResources[j]
                        resourcePath = getHCP.getSubjectResourceMeta().get('RealPath')
                        getHCP.Scan = StructuralSeriesDescScanIdDict.get(StructuralSeriesList[j])
                        # this is a hack to account for archive1,2,3 and the discrepancy in the DB...
                        scanPathSplit = getHCP.getScanMeta().get('Path')[0].split('/')
                        scanPathSub = '/'.join(scanPathSplit[scanPathSplit.index(Project):])
                        if resourcePath:
                            if (' '.join(resourcePath).index(scanPathSub) != -1): PathMatch.append(True)
                            else: PathMatch.append(False)
                    
                    if all(PathMatch):
                        print SubmitStr
                        #os.system(SubmitStr)
                    else:
                        print 'ERROR: file paths mismatch for subject %s, session %s, pipeline %s, on server %s.' % (getHCP.Subject, getHCP.Session, Pipeline, getHCP.Server)
            else:
                print 'ERROR: GEFieldMapGroup mismatch for subject %s, session %s, pipeline %s, on server %s.' % (getHCP.Subject, getHCP.Session, Pipeline, getHCP.Server) 
                
        #=======================================================================
        # FunctionalHCP
        #=======================================================================
        elif (Pipeline == 'FunctionalHCP'):
                
            if (FunctionalList.count(currSeries) == 1):
                FuncScanId = idList[sessionMeta.get('Series').index(currSeries)]
                FuncQuality = qualityList[sessionMeta.get('Series').index(currSeries)]
                getHCP.Scan = FuncScanId
                FuncScanParms = getHCP.getScanParms()
            else:
                print 'OOPS, FunctionalHCP mismatch with FunctionalList and FunctSeries'
                
            if (seriesList.count(currSeries + '_SBRef') == 1):
                ScoutScanId = idList[seriesList.index(currSeries + '_SBRef')]
                ScoutQuality = qualityList[seriesList.index(currSeries + '_SBRef')]
                getHCP.Scan = ScoutScanId
                ScoutScanParms = getHCP.getScanParms()
            else:
                print 'OOPS, FunctionalHCP mismatch with FunctionalList SBRef and FunctSeries'
            
            getHCP.Scan = FuncScanId
            fucntScanMeta = getHCP.getScanMeta()
            functScanParms = getHCP.getScanParms()
            
            
            #===================================================================
            # Here be dragons...series descriptions changed for CDB, broke IntraDB
            #===================================================================
            magScanCount = seriesList.count('SpinEchoFieldMap_LR')
            
            magScanIdList = list()
            phaScanIdList = list()
            magScanTimeList = list()
            phaScanTimeList = list()
            magShimGroupList = list()
            phaShimGroupList = list()
            magScanDiffList = list()
            phaScanDiffList = list()
            
            currMagScanIdx = 0
            currPhaScanIdx = 0
            for j in xrange(0, magScanCount):
                currMagScanId = idList[seriesList.index('SpinEchoFieldMap_LR', int(currMagScanIdx))]
                currMagScanIdx = seriesList.index('SpinEchoFieldMap_LR', int(currMagScanIdx)) + 1
                magScanIdList.append(currMagScanId)
                getHCP.Scan = currMagScanId
                magScanMeta = getHCP.getScanMeta()
                magScanParms = getHCP.getScanParms()
                magScanTimeList.append(magScanParms.get('AcquisitionTime'))
                magShimGroupList.append(magScanParms.get('ShimGroup'))
                
                currPhaScanId = idList[seriesList.index('SpinEchoFieldMap_RL', int(currPhaScanIdx))]
                currPhaScanIdx = seriesList.index('SpinEchoFieldMap_RL', int(currPhaScanIdx)) + 1
                phaScanIdList.append(currPhaScanId)
                getHCP.Scan = currPhaScanId
                phaScanMeta = getHCP.getScanMeta()
                phaScanParms = getHCP.getScanParms()
                phaScanTimeList.append(magScanParms.get('AcquisitionTime'))
                phaShimGroupList.append(magScanParms.get('ShimGroup'))
                

#                magScanTimeList.append(magScanAcqTime)
                magScanDelta = datetime.datetime.strptime(FuncScanParms.get('AcquisitionTime'), '%H:%M:%S') - datetime.datetime.strptime(magScanParms.get('AcquisitionTime'), '%H:%M:%S')
                magScanDiffList.append(magScanDelta.seconds)
                
                phaScanDelta = datetime.datetime.strptime(FuncScanParms.get('AcquisitionTime'), '%H:%M:%S') - datetime.datetime.strptime(phaScanParms.get('AcquisitionTime'), '%H:%M:%S')
                phaScanDiffList.append(phaScanDelta.seconds)
        
            
            minMagIdx = magScanDiffList.index(min(magScanDiffList)) 
            minPhaIdx = phaScanDiffList.index(min(phaScanDiffList)) 
            functMagGroupIdx = magShimGroupList.index(functScanParms.get('ShimGroup'))
            functPhaGroupIdx = phaShimGroupList.index(functScanParms.get('ShimGroup'))
            

                
            MagScanId = magScanIdList[functMagGroupIdx]
            PhaScanId = phaScanIdList[functPhaGroupIdx]
            #------------------------------------------
            launcherMagScanId = '-parameter magscanid=%s ' % (MagScanId)
            launcherPhaScanId = '-parameter phascanid=%s ' % (PhaScanId)
            #------------------------------------------
            launcherFuncScanId = '-parameter functionalscanid=%s ' % (FuncScanId)
            launcherScoutScanId = '-parameter scoutscanid=%s ' % (ScoutScanId)
            launcherFunctSeries = '-parameter functionalseries=%s ' % (currSeries)
#            100307_3T_SpinEchoFieldMap_LR.nii.gz
            launcherLR_Fieldmap = '-parameter lr_fieldmapseries=SpinEchoFieldMap_LR '
            launcherRL_Fieldmap = '-parameter rl_fieldmapseries=SpinEchoFieldMap_RL '
            launcherDwellTime = '-parameter DwellTime=%s ' % (str( float(FuncScanParms.get('EchoSpacing')) ))
            launcherUnwarpDir = '-parameter UnwarpDir=%s ' % (FuncScanParms.get('PhaseEncodingDir'))
            launcherDistortionCorrect = '-parameter DistortionCorrection=TOPUP '
            TE = '-parameter TE=2.46 '
            # MG: The TE parameter is actually Delta TE and is for the field map, not the T1w or T2w scans.  If you look for delta TE under the field map in the DB you will find it: 2.46ms.  This will change for 7T vs 3T but will otherwise always be the same.
            # NOTE: also important for functionalHCP is distortion correction is TOPUP, so fieldmap distortion correction is not even used.  TE could be anything and it would not matter.
            #-------------------------------------------
            
            SubmitStr = JobSubmitter + PipelineLauncher + launcherPipeline + launcherHCPid + DataType + Host + XnatServer + launcherProject + launcherXnatId + launcherLabel + launcherUser + launcherPassword +  SupressNotify + NotifyUser + NotifyAdmin + AdminEmail + UserEmail + MailHost + UserFullName +\
            BuildDir + launcherSession + launcherSubject + launcherMagScanId + launcherPhaScanId + launcherFuncScanId + launcherScoutScanId + \
            launcherFunctSeries + launcherLR_Fieldmap + launcherRL_Fieldmap + launcherDwellTime + TE + launcherUnwarpDir + launcherDistortionCorrect + launcherTemplatesDir + launcherConfigDir + launcherCaretAtlasDir + RedirectionStr
            
            if (magShimGroupList[minMagIdx] == phaShimGroupList[minPhaIdx] == functScanParms.get('ShimGroup')) and (functScanParms.get('SEFieldMapGroup') == magScanParms.get('SEFieldMapGroup') == phaScanParms.get('SEFieldMapGroup')):
                print 'ShimGroup and SEFieldMapGroup match successful...'
                if sys.platform == 'win32':
                    print SubmitStr
                else:
                    print SubmitStr
                    os.system(SubmitStr)
            else:
                print 'WARNING: ShimGroup or SEFieldMapGroup mismatch for subject %s, session %s, series %s, on server %s.' % (getHCP.Subject, getHCP.Session, currSeries, getHCP.Server)
                

                
        #===============================================================================
        # FIX HCP....
        #===============================================================================
        elif (Pipeline == 'FIX_HCP'):

            launcherBP = '-parameter BP=%s ' % (str(2000))
            launcherFunctSeries = '-parameter functseries=%s ' % (currSeries)
            
            #===================================================================
            # -label 100307_fnca 
            # -host https://hcpi-dev-cuda00.nrg.mir
            # -parameter xnatserver=ConnectomeDB 
            # -parameter xnat_id=HCPIntradb_E04466
            # -parameter project=HCP_Phase2
            # -parameter sessionid=100307_fnca 
            # -parameter subjects=100307 
            # -parameter BP=2000 
            # -parameter functionalseries=BOLD_REST1_RL
            # -parameter builddir=/data/hcpdb/build/HCP_Phase2/1362341272 
            #===================================================================
            SubmitStr = JobSubmitter + PipelineLauncher + launcherPipeline + launcherHCPid + DataType + Host + XnatServer + launcherProject + launcherXnatId + launcherLabel + launcherUser + launcherPassword +\
            SupressNotify + NotifyUser + NotifyAdmin + AdminEmail + UserEmail + MailHost + UserFullName + launcherSubject + launcherSession + BuildDir + launcherBP + launcherFunctSeries + RedirectionStr
            
            if sys.platform == 'win32':
                print SubmitStr
            else:
                print SubmitStr
                os.system(SubmitStr)

        
            
        
            
            
            
        if (linIdx < ( len(SubjectsList) * len(FunctionalList) )):
            print 'Sleeping for ' + str(SleepTime) + ' seconds...'         
            time.sleep(SleepTime)
        else:
            print 'Done...total launch time was %s seconds for %s jobs with a sleep time of %s seconds per job...' % ( (time.time() - sTime), ( len(SubjectsList) * len(FunctionalList) ), str(SleepTime) ) 
    
if __name__ == '__main__':
    pass
