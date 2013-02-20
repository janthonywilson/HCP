'''
Created on Sep 3, 2012

@author: Tony
'''
import base64
import sys
import xml
import os
import time
import numpy
#import array
import urllib2
import argparse
import datetime
import xml.etree.ElementTree as ET


sTime = time.time()

#===============================================================================
# PARSE INPUT
# -User tony -Password ###### -Pipeline FunctionalHCP -Production True -Shadow 3 -Subjects 255639 -FunctSeries BOLD_REST3_RL
# for testing diffusion...
# -User tony -Password asdf -Pipeline DiffusionHCP -Subjects 197550 -Build 1,2,3
# structural failures...
# -User tony -Password asdf -Pipeline StructuralHCP -Build 1,2,3 -Shadow 1,2,3,4 -Subjects 167743,174437,182739,199251,688569
#===============================================================================
parser = argparse.ArgumentParser(description="Script to generate proper command for XNAT functional pipeline lauching ...")

parser.add_argument("-User", "--User", dest="iUser", default='tony', type=str)
parser.add_argument("-Password", "--Password", dest="iPassword", default='none', type=str)
parser.add_argument("-Project", "--Project", dest="iProject", default='HCP_Phase2', type=str)
parser.add_argument("-Production", "--Production", dest="iProduction", default=False, type=bool)
parser.add_argument("-Shadow", "--Shadow", dest="iShadowInt", default=None, type=str)
parser.add_argument("-Build", "--Build", dest="iBuildInt", default=None, type=str)

parser.add_argument("-Pipeline", "--Pipeline", dest="iPipeline", default='fMRIVolume', type=str)
parser.add_argument("-HCPid", "--HCPid", dest="iHCPid", default='HCPIntradb_E00000', type=str)
parser.add_argument("-Label", "--Label", dest="iLabel", default='00', type=str)
parser.add_argument("-Subjects", "--Subjects", dest="iSubjects", default='00', type=str)
parser.add_argument("-FuncScanId", "--FuncScanId", dest="iFuncScanId", default='00', type=str)
parser.add_argument("-ScoutScanId", "--ScoutScanId", dest="iScoutScanId", default='00', type=str)
parser.add_argument("-FunctSeries", "--FunctSeries", dest="iFunctSeries", default='BOLD_00', type=str)
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

parser.add_argument('--version', action='version', version='%(prog)s 0.2')

args = parser.parse_args()

iUser = args.iUser
iPassword = args.iPassword
iProject = args.iProject
iProduction = args.iProduction
iShadowInt = args.iShadowInt
iBuildInt = args.iBuildInt

iPipeline = args.iPipeline
iHCPid = args.iHCPid
iLabel = args.iLabel
iSubjects = args.iSubjects
iFuncScanId = args.iFuncScanId
iScoutScanId = args.iScoutScanId
iFunctSeries = args.iFunctSeries
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

if (iPipeline.find('Funct') != -1):
    iStructFunctDiffSeries = iFunctSeries
elif (iPipeline.find('Struct') != -1):
    iStructFunctDiffSeries = iStructSeries
elif (iPipeline.find('Diff') != -1):
    iStructFunctDiffSeries = iDiffSeries
    
tmpStructFunctDiffSeries = iStructFunctDiffSeries

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
preChangeRLList = ('100307','111312','114924','119833','125525','138231','144266','150423','159239','162329','167743','174437','185139','192439','197550','199251','217429','249947','255639','329440','355542','499566','585862','611231','665254','672756','792564','826353','877168','896778')


#===============================================================================
# FUNCTIONS
#===============================================================================
def fReadURL(URL, User, Password):
    restRequest = urllib2.Request(URL)

    restAuthHeader = "Basic %s" % base64.encodestring('%s:%s' % (User, Password))[:-1]
    restRequest.add_header("Authorization", restAuthHeader)

    restConnHandle = urllib2.urlopen(restRequest)
    restResults = restConnHandle.read()
    return restResults
#===============================================================================
def fGetSessionParms(inputUser, inputPassword, inputProject, inputSubject, inputSession):
    idList = list()
    typeList = list()
    seriesList = list()
    qualityList = list()
    xnatidList = list()

#    restRoot = 'https://intradb.humanconnectome.org/REST/'
    restRoot = 'https://hcpi-dev-cuda00.nrg.mir/data/'
    restURL = restRoot + 'projects/' + inputProject + '/subjects/' + inputSubject + '/experiments/' + inputSession + '/scans?format=csv&columns=ID,type,series_description,quality,xnat:mrSessionData/id'
    
    restResults = fReadURL(restURL, inputUser, inputPassword)
    
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
        
    return idList, typeList, seriesList, qualityList, xnatidList
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
def fGetSession(inputUser, inputPassword, inputProject, inputSubject, inputSeries):

#    restRoot = 'https://intradb.humanconnectome.org/REST/'
    restRoot = 'https://hcpi-dev-cuda00.nrg.mir/data/'
    restURL = restRoot + 'projects/' + inputProject + '/subjects/' + inputSubject + '/experiments?format=csv&columns=ID,label'
    
    restResults = fReadURL(restURL, inputUser, inputPassword)
    
    restResultsSplit = restResults.split('\n')
    restEndCount = restResults.count('\n')
    restSessionHeader = restResultsSplit[0]
    restSessionHeaderSplit = restSessionHeader.split(',')
    
#    idIdx = restSessionHeaderSplit.index('"ID"')
    labelIdx = restSessionHeaderSplit.index('"label"')
    
    seriesResultList = list()
    for j in xrange(1, restEndCount):
        currRow = restResultsSplit[j]
        
        currRowSplit = currRow.split(',')

#        currId = currRowSplit[idIdx].replace('"', '')
        currSeries = currRowSplit[labelIdx].replace('"', '')
        
        if (currSeries.find('fnc') != -1) or (currSeries.find('str') != -1) or (currSeries.find('diff') != -1) or (currSeries.find('xtr') != -1):
            idList, typeList, seriesList, qualityList, xnatidList = fGetSessionParms(inputUser, inputPassword, inputProject, inputSubject, currSeries)
            
            
            try: 
                seriesMatch = seriesList[seriesList.index(inputSeries)]
                seriesResultList.append(currSeries)
            except ValueError:
                pass

    return seriesResultList
#===============================================================================
def fGetAcquisitionTime(inputUser, inputPassword, inputProject, inputSubject, inputSession, inputScan):
#    https://intradb.humanconnectome.org/data/projects/HCP_Phase2/subjects/792564/experiments/792564_strc/scans/2
#    restRoot = 'https://intradb.humanconnectome.org/data/'
    restRoot = 'https://hcpi-dev-cuda00.nrg.mir/data/'
    restURL = restRoot + 'projects/' + inputProject + '/subjects/' + inputSubject + '/experiments/' + inputSession + '/scans/' + inputScan
    xmlData = fReadURL(restURL, inputUser, inputPassword)
    acqTimeET = ET.fromstring(xmlData)
#    print acqTime.tag.
    acqTime = acqTimeET.find('{http://nrg.wustl.edu/xnat}startTime').text
    return acqTime
#===============================================================================
def fGetScanParms(inputUser, inputPassword, inputProject, inputSubject, inputSession, inputScan):
    resultsParms = list()
#    restRoot = 'https://intradb.humanconnectome.org/data/'
    restRoot = 'https://hcpi-dev-cuda00.nrg.mir/data/'
    restURL = restRoot + 'projects/' +inputProject+ '/subjects/' +inputSubject+ '/experiments/' +inputSession+ '/scans/' +inputScan
    xmlData = fReadURL( restURL, inputUser, inputPassword )
    parmsET = ET.fromstring(xmlData)
    scanParms = parmsET.find('{http://nrg.wustl.edu/xnat}parameters')
    sampleSpacing = scanParms.find('{http://nrg.wustl.edu/xnat}readoutSampleSpacing').text
    resultsParms.append(sampleSpacing)
    
    for addParms in scanParms.findall('{http://nrg.wustl.edu/xnat}addParam'):
        addParmsAttrib = addParms.attrib
        
        if (addParmsAttrib.get('name') == 'Siemens GRADSPEC alShimCurrent'):
#            print addParms.text
            alShimCurrent = addParms.text
            resultsParms.append(addParms.text)
            
        if (addParmsAttrib.get('name') == 'Siemens GRADSPEC lOffset'):
#            print addParms.text
            LinOffset = addParms.text
            resultsParms.append(addParms.text)
    
    resultsParmsDict = { 'SampleSpacing': sampleSpacing, 'alShimCurrent': alShimCurrent, 'LinearOffset':  LinOffset }
    return resultsParmsDict
#===============================================================================
# Volume Parms...
# -User tony -Password asdf -Pipeline FunctionalHCP -Production True -Shadow 1,2,3,4 -Subjects 100307,103515,255639,156637,125525,161731 -FunctSeries BOLD_SOCIAL1_RL,BOLD_SOCIAL2_LR,BOLD_GAMBLING1_RL,BOLD_WM2_LR
#===============================================================================

SleepTime = 3
SubjectsList = iSubjects.split(',')
if (iShadowInt != None):
    ShadowList = iShadowInt.split(',')
else:
    ShadowList = ('inf')
    
if (iBuildInt != None):
    BuildList = iBuildInt.split(',')
else:
    BuildList = ('')
    
SeriesList = iStructFunctDiffSeries.split(',')

if (len(SubjectsList) > len(ShadowList)):
#    ShadowArray = numpy.tile(ShadowList, (numpy.ceil(len(SubjectsList) / len(ShadowList))))
    ShadowArray = numpy.tile(ShadowList, (numpy.ceil(len(SubjectsList))))
else: 
    ShadowArray = numpy.tile(ShadowList, (1))

if (len(SubjectsList) > len(BuildList)):
    BuildArray = numpy.tile(BuildList, (numpy.ceil(len(SubjectsList))))
else: 
    BuildArray = numpy.tile(BuildList, 1)
    
linIdx = 0
for h in xrange(0, len(SubjectsList)): 
    iSubject = SubjectsList[h]
    
    for i in xrange(0, len(SeriesList)):
        linIdx += 1
        
        if (iPipeline.find('Funct') != -1):
            iFunctSeries = SeriesList[i]
            iStructFunctDiffSeries = iFunctSeries
        elif (iPipeline.find('Struct') != -1):
            iStructFunctDiffSeries = iStructSeries
        elif (iPipeline.find('Diffusion') != -1):
            iStructFunctDiffSeries = iDiffSeries
        
        tmpFunctSeries = tmpStructFunctDiffSeries
    
        if ( (tmpFunctSeries.find('REST') != -1) and (preChangeRLList.count(iSubject) > 0) ):
            if (tmpFunctSeries.find('REST3') != -1) or (tmpFunctSeries.find('REST1') != -1):
                iFunctSeries = 'BOLD_'+ tmpFunctSeries +'_RL'
                iStructFunctDiffSeries = iFunctSeries
            elif (tmpFunctSeries.find('REST4') != -1) or (tmpFunctSeries.find('REST2') != -1):
                iFunctSeries = 'BOLD_'+ tmpFunctSeries +'_LR'
                iStructFunctDiffSeries = iFunctSeries
        elif ( (tmpFunctSeries.find('REST') != -1) and (preChangeRLList.count(iSubject) == 0) ):
            if (tmpFunctSeries.find('REST3') != -1) or (tmpFunctSeries.find('REST2') != -1):
                iFunctSeries = 'BOLD_'+ tmpFunctSeries +'_LR'
                iStructFunctDiffSeries = iFunctSeries
            elif (tmpFunctSeries.find('REST4') != -1) or (tmpFunctSeries.find('REST1') != -1):
                iFunctSeries = 'BOLD_'+ tmpFunctSeries +'_RL'
                iStructFunctDiffSeries = iFunctSeries
                
#        print iUser, iPassword, iProject, iSubject, iStructFunctDiffSeries
        SessionId = fGetSession( iUser, iPassword, iProject, iSubject, iStructFunctDiffSeries  )
        
        if (len(SessionId) < 1):
            print 'ERROR: No session id could be found for subject ' + iSubject + ' with series ' +iStructFunctDiffSeries

        elif (len(SessionId) >= 1):
            if (len(SessionId) == 1):
                iSessionId = SessionId[0]
            elif (len(SessionId) > 1):
                if ((iSubject +'_xtra') in SessionId):
                    iSessionId = SessionId[SessionId.index(iSubject +'_xtra')]
                elif ((iSubject +'_xtrb') in SessionId):
                    iSessionId = SessionId[SessionId.index(iSubject +'_xtrb')]
            
            idList, typeList, seriesList, qualityList, xnatidList = fGetSessionParms( iUser, iPassword, iProject, iSubject, iSessionId  )
            
            iXNATid = xnatidList[0]
            if (BuildArray.size >= len(SubjectsList)):
                iBuildDirRoot = '/data/intradb/build' + str(BuildArray[h]) + '/' + iProject + '/'
            else:
                iBuildDirRoot = '/data/intradb/build/' + iProject + '/'
#            print iBuildDirRoot
            
            Project = '-parameter project=%s ' % iProject 
            Pipeline = '-pipeline /data/intradb/pipeline/catalog/%s/%s.xml ' % (iPipeline, iPipeline)
            User = '-u %s ' % iUser 
            Password = '-pwd %s ' % iPassword 
            Label = '-label %s ' % iSessionId
            HCPid = '-id %s ' % iXNATid
            XnatId = '-parameter xnat_id=%s ' % iXNATid 
            SessionId = '-parameter sessionid=%s ' % iSessionId 
            Subjects = '-parameter subjects=%s ' % iSubject

            
            if (len(SeriesList) > 1):
                currBuildDir = iBuildDirRoot + str(numpy.asarray(round(time.time()), dtype=numpy.uint64))
            else:
                currBuildDir = iBuildDirRoot + str(numpy.asarray(round(sTime), dtype=numpy.uint64))
            BuildDir = '-parameter builddir=' + currBuildDir + ' '

            
            if not os.path.exists(currBuildDir + os.sep + iSubject) and sys.platform != 'win32':
                os.makedirs(currBuildDir + os.sep + iSubject)
                
            RedirectionStr = ' > ' + currBuildDir.replace(' ', '') + os.sep + iSubject + os.sep + iPipeline + 'LaunchSTDOUT.txt'
            
            if (iProduction):
                if (iShadowInt != None):
                    Host = '-host https://intradb-shadow' + ShadowArray[h] + '.nrg.mir '
                else:
                    Host = '-host https://intradb.humanconnectome.org '
            else:
                Host = '-host https://hcpi-dev-cuda00.nrg.mir ' 
                
            
            #===============================================================================
            # START BUILDING STRING....
            #===============================================================================
            if (iPipeline == 'fMRIVolume') or (iPipeline == 'FunctionalHCP'):
                
                    
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
                
            elif (iPipeline == 'DiffusionHCP'):
                #===============================================================================
                # -User tony -Password asdf -Pipeline DiffusionHCP -Subjects 103515,100307,111312,114924,119833,125525,138231,144266,150423,159239,162329,167743,174437,185139,192439,197550,199251,217429,249947,255639,329440,355542,499566,585862,611231,665254,672756,792564,826353,877168,896778 -DiffSeries DWI_RL_dir95
                #
                # /data/intradb/pipeline/bin/PipelineJobSubmitter /data/intradb/pipeline/bin/XnatPipelineLauncher -pipeline /data/intradb/pipeline/catalog/DiffusionHCP/DiffusionHCP.xml -id HCPIntradb_E05324 -host https://hcpi-dev-cuda00.nrg.mir -u tony -pwd asdf -dataType xnat:mrSessionData -label 111312_diff -supressNotification -project HCP_Phase2 -parameter xnat_id=HCPIntradb_E05324 -parameter sessionid=111312_diff -parameter subjects=111312 -parameter LR_1ScanId=12 -parameter LR_2ScanId=16 -parameter LR_3ScanId=20 -parameter RL_1ScanId=10 -parameter RL_2ScanId=14 -parameter RL_3ScanId=18 -parameter EchoSpacing=0.7800117313764398 -parameter PhaseEncodingDir=1 -notify wilsont@mir.wustl.edu -notify db-admin@humanconnectome.org -parameter mailhost=mail.nrg.wustl.edu -parameter userfullname=T.Wilson -parameter builddir=/data/intradb/build/HCP_Phase2/20121128_00004 -parameter xnatserver=HCPIntradb -parameter adminemail=db-admin@humanconnectome.org -parameter useremail=wilsont@mir.wustl.edu -parameter Dir1=95 -parameter Dir2=96 -parameter Dir3=97 -parameter project=HCP_Phase2
                # /data/intradb/pipeline/bin/PipelineJobSubmitter /data/intradb/pipeline/bin/XnatPipelineLauncher -pipeline /data/intradb/pipeline/catalog/DiffusionHCP/DiffusionHCP.xml -id HCPIntradb_E04787 -host https://hcpi-dev-cuda00.nrg.mir -u tony -pwd asdf -dataType xnat:mrSessionData -label 114924_diff -supressNotification -project HCP_Phase2 -parameter xnat_id=HCPIntradb_E04787 -parameter sessionid=114924_diff -parameter subjects=114924 -parameter LR_1ScanId=12 -parameter LR_2ScanId=16 -parameter LR_3ScanId=20 -parameter RL_1ScanId=10 -parameter RL_2ScanId=14 -parameter RL_3ScanId=18 -parameter EchoSpacing=0.7800117313764398 -parameter PhaseEncodingDir=1 -notify wilsont@mir.wustl.edu -notify db-admin@humanconnectome.org -parameter mailhost=mail.nrg.wustl.edu -parameter userfullname=T.Wilson -parameter builddir=/data/intradb/build/HCP_Phase2/20121128_00005 -parameter xnatserver=HCPIntradb -parameter adminemail=db-admin@humanconnectome.org -parameter useremail=wilsont@mir.wustl.edu -parameter Dir1=95 -parameter Dir2=96 -parameter Dir3=97 -parameter project=HCP_Phase2
                # /data/intradb/pipeline/bin/PipelineJobSubmitter /data/intradb/pipeline/bin/XnatPipelineLauncher -pipeline /data/intradb/pipeline/catalog/DiffusionHCP/DiffusionHCP.xml -id HCPIntradb_E07771 -host https://intradb-shadow1.nrg.mir -u tony -pwd asdf -dataType xnat:mrSessionData -label 103515_strc -supressNotification -parameter project=HCP_Phase2 -notify wilsont@mir.wustl.edu -parameter xnat_id=HCPIntradb_E07771 -parameter sessionid=103515_strc -parameter subjects=103515 -parameter Dir1=95 -parameter Dir2=96 -parameter Dir3=97 -parameter EchoSpacing=0.7800117313764398 -parameter PhaseEncodingDir=1 -parameter LR_1ScanId=12 -parameter LR_2ScanId=16 -parameter LR_3ScanId=20 -parameter RL_1ScanId=10 -parameter RL_2ScanId=14 -parameter RL_3ScanId=18 -notify db-admin@humanconnectome.org -parameter mailhost=mail.nrg.wustl.edu -parameter userfullname=T.Wilson -parameter builddir=/data/intradb/build/HCP_Phase2/1354406629 -parameter xnatserver=HCPIntradb -parameter adminemail=db-admin@humanconnectome.org -parameter useremail=wilsont@mir.wustl.edu  > /data/intradb/build/HCP_Phase2/1354406629\103515\DiffusionHCPLaunchSTDOUT.txt
                #===============================================================================
                
                EchoSpacing = '-parameter EchoSpacing=0.7800117313764398 '
                PhaseEncodingDir = '-parameter PhaseEncodingDir=1 '
                
                DiffusionSeriesList = ('DWI_RL_dir95','DWI_LR_dir95','DWI_RL_dir96','DWI_LR_dir96','DWI_RL_dir97','DWI_LR_dir97')
                ScanIdDict = {'RL_1ScanId' : None, 'RL_2ScanId' : None, 'RL_3ScanId' : None, 'LR_1ScanId' : None, 'LR_2ScanId' : None, 'LR_3ScanId' : None}
                
                DiffusionSeriesIntersectList = list(set(DiffusionSeriesList) & set(seriesList))
                
                for j in xrange(0, len(DiffusionSeriesList)):
                    currDiffDesc = DiffusionSeriesList[j]
                    if (seriesList.count(currDiffDesc) > 0):
                        currDiffIdx = seriesList.index(currDiffDesc)
                        currTypeList = typeList[currDiffIdx]
                        currScanId = idList[currDiffIdx]
                        currQuality = qualityList[currDiffIdx]
                        
                    if ((currDiffDesc.find('RL_dir95') != -1) and ( (currQuality == 'usable') or (currQuality == 'excellent') or (currQuality == 'good') or (currQuality == 'undetermined') )):
                        ScanIdDict['RL_1ScanId'] = '-parameter RL_1ScanId=%s ' % str(currScanId)
                    elif ((currDiffDesc.find('RL_dir96') != -1) and ( (currQuality == 'usable') or (currQuality == 'excellent') or (currQuality == 'good') or (currQuality == 'undetermined') )):
                        ScanIdDict['RL_2ScanId'] = '-parameter RL_2ScanId=%s ' % str(currScanId)
                    elif ((currDiffDesc.find('RL_dir97') != -1) and ( (currQuality == 'usable') or (currQuality == 'excellent') or (currQuality == 'good') or (currQuality == 'undetermined') )):
                        ScanIdDict['RL_3ScanId'] = '-parameter RL_3ScanId=%s ' % str(currScanId)
                    elif ((currDiffDesc.find('LR_dir95') != -1) and ( (currQuality == 'usable') or (currQuality == 'excellent') or (currQuality == 'good') or (currQuality == 'undetermined') )):
                        ScanIdDict['LR_1ScanId'] = '-parameter LR_1ScanId=%s ' % str(currScanId)
                    elif ((currDiffDesc.find('LR_dir96') != -1) and ( (currQuality == 'usable') or (currQuality == 'excellent') or (currQuality == 'good') or (currQuality == 'undetermined') )):
                        ScanIdDict['LR_2ScanId'] = '-parameter LR_2ScanId=%s ' % str(currScanId)
                    elif ((currDiffDesc.find('LR_dir97') != -1) and ( (currQuality == 'usable') or (currQuality == 'excellent') or (currQuality == 'good') or (currQuality == 'undetermined') )):
                        ScanIdDict['LR_3ScanId'] = '-parameter LR_3ScanId=%s ' % str(currScanId)
                    else:
                        if (currDiffDesc.find('LR_dir95') > 0): 
                            ScanIdDict['LR_1ScanId'] = '-parameter LR_1ScanId=EMPTY '
                        elif (currDiffDesc.find('LR_dir96') > 0): 
                            #===================================================
                            # Oi! Look here nitwit...
                            #===================================================
                            print [i for i, x in enumerate(seriesList) if x == currDiffDesc]
                            ScanIdDict['LR_2ScanId'] = '-parameter LR_2ScanId=EMPTY '
                        elif (currDiffDesc.find('LR_dir97') > 0): 
                            ScanIdDict['LR_3ScanId'] = '-parameter LR_3ScanId=EMPTY '
                        elif (currDiffDesc.find('RL_dir95') > 0): 
                            ScanIdDict['RL_1ScanId'] = '-parameter RL_1ScanId=EMPTY '
                        elif (currDiffDesc.find('RL_dir96') > 0): 
                            ScanIdDict['RL_2ScanId'] = '-parameter RL_2ScanId=EMPTY '
                        elif (currDiffDesc.find('RL_dir97') > 0): 
                            ScanIdDict['RL_3ScanId'] = '-parameter RL_3ScanId=EMPTY '
                
                
                Dir1 = '-parameter Dir1=95 '
                Dir2 = '-parameter Dir2=96 '
                Dir3 = '-parameter Dir3=97 '
                
                
                SubmitStr = JobSubmitter + PipelineLauncher + Pipeline + HCPid + Host + User + Password + DataType + Label + SupressNotify + Project + NotifyUser + XnatId + SessionId + Subjects + Dir1 + Dir2 + Dir3 + \
                EchoSpacing + PhaseEncodingDir + ScanIdDict['RL_1ScanId'] + ScanIdDict['RL_2ScanId'] + ScanIdDict['RL_3ScanId'] + ScanIdDict['LR_1ScanId'] + ScanIdDict['LR_2ScanId'] + ScanIdDict['LR_3ScanId'] + \
                NotifyAdmin + MailHost + UserFullName + BuildDir + XnatServer + AdminEamil + UserEmail + RedirectionStr
                
                if sys.platform == 'win32':
                    print SubmitStr.find('EMPTY')
#                    if (SubmitStr.find('EMPTY') == -1):
                    print SubmitStr
#                    else: 
#                        print 'ERROR: EMPTY found for diffusion scan...'
                else:
                    print SubmitStr
                    if (SubmitStr.find('EMPTY') == -1):
                        os.system(SubmitStr)
                    else:
                        continue
                
                
                
            if (linIdx < ( len(SubjectsList) * len(SeriesList) )):
                print 'Sleeping for ' + str(SleepTime) + ' seconds...'         
                time.sleep(SleepTime)
            else:
                print 'Done...total launch time was %s seconds for %s jobs with a sleep time of %s seconds per job...' % ( (time.time() - sTime), ( len(SubjectsList) * len(SeriesList) ), str(SleepTime) ) 
        
        
            
            
            
            
            
            
        
        
        
