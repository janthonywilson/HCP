'''
Created on 2012-12-19

@author: jwilso01
'''

# multiplatform system stuff...
import os
import sys
# Web stuff...
import socket
import base64
import urllib
import urllib2
#import requests
#import httplib2
from ssl import SSLError
import xml.etree.ElementTree as ET
from urllib2 import URLError, HTTPError

#===============================================================================
# CLASSES
#===============================================================================
class getHCP:
    """HCP Interfacing Class for GETs"""
    def __init__( self, User, Password, Server  ):
        self.User = User
        self.Password = Password
        
        Server.strip()
        if (Server[-1] != '/'):
            Server = Server + '/'
        if (Server.find('http') == -1):
            Server = 'https://' + Server
        self.Server = Server
        
        self.Verbose = False
        self.Timeout = 8
        self.TimeoutMax = 256
        self.TimeoutStep = 8
        
        self.Project = ''
        self.Projects = []
        self.ProjectNames = []
        self.ProjectsSecondary = []
        
        self.Session = ''
        self.Sessions = []
        self.SessionType = ''
        self.SessionTypes = []
        self.SessionParms = {}
        
        self.Subject = ''
        self.Subjects = []
        self.SubjectSessions = []
        self.SubjectSessionsUniq = []
        
        self.Scan = ''
        self.Scans = []
        
        self.FileInfo = {}
        
#        self.SessionId = 'BE07716A225958F2FDD51D26E0D26449'
        self.SessionId = self.getSessionId()
    #===============================================================================
    def getSessionId( self ):
        """Get session id for getHCP session spawn"""
        URL = self.Server + 'data/JSESSION'

        #=======================================================================
        # # httplib2
        # h = httplib2.Http(".cache")
        # h.add_credentials(self.User, self.Password)
        # r, content = h.request(URL, "GET")
        #=======================================================================

        # URLLIB2
        Request = urllib2.Request(URL)
        basicPasswordManager = urllib2.HTTPPasswordMgrWithDefaultRealm()
        basicPasswordManager.add_password(None, URL, self.User, self.Password)
        basicAuthHandler = urllib2.HTTPBasicAuthHandler(basicPasswordManager)
        openerURL = urllib2.build_opener(basicAuthHandler)
        urllib2.install_opener(openerURL)
        
#        restPost = urllib.urlencode({'foo' : 'bar'})
#        restRequest = urllib2.Request(URL, restPost)
#        restAuthHeader = "Basic %s" % base64.encodestring('%s:%s' % (self.User, self.Password))[:-1]
#        restRequest.add_header("Authorization", restAuthHeader)

        while (self.Timeout <= self.TimeoutMax):
            try:
                connHandle = urllib2.urlopen(Request, None, self.Timeout)
                break
            except URLError, e:
                self.Timeout += self.TimeoutStep
                print 'URLError code: ' +str(e.reason)+ '. Timeout increased to ' +str(self.Timeout)+' seconds for JSESSION cookie...'
                
                
        self.SessionId = connHandle.read()
        return self.SessionId
    #===============================================================================
    def getURLString( self, URL ):
        """Get URL results as string"""
        restRequest = urllib2.Request(URL)
        restRequest.add_header("Cookie", "JSESSIONID=" + self.SessionId);
    
        while (self.Timeout <= self.TimeoutMax):
            try:
                restConnHandle = urllib2.urlopen(restRequest, None, self.Timeout)
            except HTTPError, e:
                if (e.code != 404):
                    self.Timeout += self.TimeoutStep
                    print 'HTTPError code: ' +str(e.code)+ '. Timeout increased to ' +str(self.Timeout)+' seconds for ' +URL
                else:
                    return '404 Error'
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
                    ReadResults = restConnHandle.read()
                    return ReadResults
                except HTTPError, e:
                    print 'HTTPError code: ' +str(e.code)+ '. File read timeout for ' +str(self.Timeout)+ ' seconds for ' +URL
                    
        print 'ERROR: No reasonable timeout limit could be found for ' + URL
        sys.exit()
    #===============================================================================
    def getProjects( self ):
        """Get a list of project names from XNAT instance"""
        
        restURL = self.Server + 'data/projects?format=csv'
        restResults = self.getURLString(restURL)
        
        restResultsSplit = restResults.split('\n')
        restEndCount = restResults.count('\n')
        restProjectHeader = restResultsSplit[0]
        restProjectHeaderSplit = restProjectHeader.split(',')
        
        projectNameIdx = restProjectHeaderSplit.index('"name"')
        projectIdx = restProjectHeaderSplit.index('"ID"')
        projectSecondaryIdx = restProjectHeaderSplit.index('"secondary_ID"')
        
        for i in xrange(1, restEndCount):
            currRow = restResultsSplit[i]
            
            currRowSplit = currRow.split(',')
            currProjectName = currRowSplit[projectNameIdx].replace('"', '')
            currProjectId = currRowSplit[projectIdx].replace('"', '')
            currProjectSecondary = currRowSplit[projectSecondaryIdx].replace('"', '')
            
            self.Projects.append(currProjectId)
            self.ProjectNames.append(currProjectName)
            self.ProjectsSecondary.append(currProjectSecondary)
            
        return self.Projects, self.ProjectNames, self.ProjectsSecondary
    #===============================================================================
    def getSubjects( self ):
        """Get all subjects for a given project"""
        Subjects = list()
        
        if (not self.Project):
            print 'Project is empty.  Must assign a Project before getting subjects.  Try this one ''' +self.getProjects()[0][0]+ ''
            return Subjects
#            sys.exit()
        else:
            
            restURL = self.Server + 'data/projects/' + self.Project + '/subjects?format=csv'
            restResults = self.getURLString(restURL)
            
            restResultsSplit = restResults.split('\n')
            restEndCount = restResults.count('\n')
            restSessionHeader = restResultsSplit[0]
            restSessionHeaderSplit = restSessionHeader.split(',')
            labelIdx = restSessionHeaderSplit.index('"label"')
            
            for i in range(1,restEndCount):
                currRow = restResultsSplit[i]
                
                currRowSplit = currRow.split(',')
                currSubject = currRowSplit[labelIdx].replace('"', '')
                
                Subjects.append(currSubject)
                
            #=======================================================================
            # if list is not unique, use set(Subjects)
            #=======================================================================
            return Subjects
    #===============================================================================    
    def getSubjectSessions( self ):
        """Get all sessions and session types for a given subject"""
        SubjectSessionsID = list()
        SubjectSessionUniq = list()
        SubjectSessionsType = list()
        
        if (not self.Project):
            print 'ERROR: No project specified...'
        
        restURL = self.Server + 'data/projects/' + self.Project + '/subjects/' + self.Subject + '/experiments?format=csv&columns=ID,label'
        restResults = self.getURLString(restURL)
        
        if (restResults != '404 Error'):
    
            restResultsSplit = restResults.split('\n')
            restEndCount = restResults.count('\n')
            restSessionHeader = restResultsSplit[0]
            restSessionHeaderSplit = restSessionHeader.split(',')
            
            labelIdx = restSessionHeaderSplit.index('"label"')
            uniqInternalId = restSessionHeaderSplit.index('"ID"')
            
            
            for i in xrange(1, restEndCount):
                currRow = restResultsSplit[i]
                
                currRowSplit = currRow.split(',')
                currSession = currRowSplit[labelIdx].replace('"', '')
                currSessionUniq = currRowSplit[uniqInternalId].replace('"', '')
                
                if (currSession.find('fnc') != -1) or (currSession.find('str') != -1) or (currSession.find('diff') != -1) or (currSession.find('xtr') != -1) or (currSession.find('3T') != -1):
                    if (currSession.find('xtr') != -1):
                        SubjectSessionsID.append(currSession)
                        self.Session = currSession
                        SessionTypeList = self.getSessionMeta( ).get('Types')
                        if ('T1w' in SessionTypeList):
                            SubjectSessionsType.append('strc')
                        elif ('dMRI' in SessionTypeList):
                            SubjectSessionsType.append('diff')
                        elif ('tfMRI' in SessionTypeList):
                            SubjectSessionsType.append('fnc')
                        else:
                            SubjectSessionsType.append('unknown')
                        
                    elif (currSession.find('3T') != -1):
                        SubjectSessionsID.append(currSession)
                        self.Session = currSession
                        SessionTypeList = self.getSessionMeta( ).get('Types')
                        if ('T1w' in SessionTypeList) and ('T2w' in SessionTypeList):
                            SubjectSessionsType.append('strc')
                        if ('dMRI' in SessionTypeList):
                            SubjectSessionsType.append('diff')
                        if ('tfMRI' in SessionTypeList):
                            SubjectSessionsType.append('task')
                        if ('rfMRI' in SessionTypeList):
                            SubjectSessionsType.append('rest')
                            
                    else:
                        if (currSession.find('fnc') != -1): SubjectSessionsType.append('fnc')
                        elif (currSession.find('strc') != -1): SubjectSessionsType.append('strc')
                        elif (currSession.find('diff') != -1): SubjectSessionsType.append('diff')
                        else: SubjectSessionsType.append('unknown')
                        SubjectSessionsID.append(currSession)
                
            SubjectSessionUniq.append(currSessionUniq)

        SubjectSessions = {'Sessions': SubjectSessionsID, 'Types': SubjectSessionsType}
        return SubjectSessions
    #===============================================================================    
    def getSubjectsSessions( self ):
        """Get all sessions and session types for all subjects"""
        
        if (not self.Subjects):
            print 'ERROR: No must specify a list of subjects...'
            print 'Correcting...'
            self.Subjects = self.getSubjects()
            print '...Please try again.'
        else:
            for i in xrange(0, len(self.Subjects)):
                self.Subject = self.Subjects[i]
                print self.getSubjectSessions()
                
            
    #===========================================================================
    # Session Level Meta Data
    #===========================================================================
    def getSessionMeta( self ):
        """Get ID, Type, Series, Quality, and XNAT ID for a given subject and session"""

        ScanIds = list()
        ScanTypes = list()
        ScanSeries = list()
        ScanQualty = list()
        ScanXnatId = list()

        if not self.Session:
            print 'No session for Session Meta Data...'
            
        restURL = self.Server + 'data/projects/' + self.Project + '/subjects/' + self.Subject + '/experiments/' + self.Session + '/scans?format=csv&columns=ID,type,series_description,quality,xnat:mrSessionData/id'
        
        restResults = self.getURLString(restURL)
        
        restResultsSplit = restResults.split('\n')
        restEndCount = restResults.count('\n')
        restSessionHeader = restResultsSplit[0]
        restSessionHeaderSplit = restSessionHeader.split(',')
        

        # ['"xnat_imagescandata_id"', '"ID"', '"type"', '"series_description"', '"quality"', '"xnat:mrsessiondata/id"', '"URI"']
        idIdx = restSessionHeaderSplit.index('"ID"')
        seriesIdx = restSessionHeaderSplit.index('"series_description"')
        typeIdx = restSessionHeaderSplit.index('"type"')
        qualityIdx = restSessionHeaderSplit.index('"quality"')
        xnatidIdx = restSessionHeaderSplit.index('"xnat:mrsessiondata/id"')
        
        for j in xrange(1, restEndCount):
            currRow = restResultsSplit[j]
            
            currRowSplit = currRow.split(',')
    
            ScanIds.append(currRowSplit[idIdx].replace('"', ''))
            ScanTypes.append(currRowSplit[typeIdx].replace('"', ''))
            ScanSeries.append(currRowSplit[seriesIdx].replace('"', ''))
            ScanQualty.append(currRowSplit[qualityIdx].replace('"', ''))
            ScanXnatId.append(currRowSplit[xnatidIdx].replace('"', ''))
            
        SessionMeta = {'IDs':ScanIds, 'Types':ScanTypes, 'Series':ScanSeries, 'Quality':ScanQualty, 'XNATID':ScanXnatId }
        return SessionMeta
    #===============================================================================
    def getSessionQuality( self ):
        """QC: Get Session, Subject, Scan Type, Quality of all data"""

        Quality = list()
        Series = list()
        ScanIds = list()
        Sessions = list()
        ScanType = list()
        
        if (not self.Subjects):
            self.Subjects = self.getSubjects( )
        else:
            self.Subjects = self.Subject.split()
            
        for i in xrange(0, len(self.Subjects)):
            self.Subject = self.Subjects[i]
            
            if (not self.Session) and (not self.Sessions):
                self.Sessions = self.getSubjectSessions()[0]
            else:
                self.Sessions = self.Session.split()
            
            for j in xrange(0, len(self.Sessions)):
                
                self.Session = self.Sessions[j]
                
                SessionMeta = self.getSessionMeta()
                
                Quality.extend(SessionMeta.get('Quality'))
                Series.extend(SessionMeta.get('Series'))
                ScanIds.extend(SessionMeta.get('IDs'))
                ScanType.extend(SessionMeta.get('Types'))
                
                Sessions.extend( self.Session.split(',') * len(SessionMeta.get('IDs')))
                
        return Quality, ScanIds, Series, Sessions, ScanType
    #===============================================================================
    def getSubjectResources(self):
        
        ResourceHeader = list()
        FileNames = list()
        FileURIs = list()
        FileSessions = list()
        FileLabels = list()
        
        SubjectSessions = self.getSubjectSessions()[0]
    
        for i in range(0, len(SubjectSessions)):
    #        restURL = self.Server + 'data/projects/' + self.Project + '/subjects/' + self.Subject + '/experiments/' + self.Session + '/scans?format=csv&columns=ID,type,series_description,quality,xnat:mrSessionData/id'
            restURL = self.Server + 'data/projects/' + self.Project + '/subjects/' + self.Subject + '/experiments/' + SubjectSessions[i] + '/resources?format=csv'
            if self.Verbose: print restURL
            
            restResults = self.getURLString(restURL)
            
            restResultsSplit = restResults.split('\n')
            restEndCount = restResults.count('\n')
        
            restSessionHeader = restResultsSplit[0]
            restSessionHeaderSplit = restSessionHeader.split(',')
            restSessionHeaderCount = restSessionHeader.count(',')
            
            for j in range(0, restSessionHeaderCount + 1):
                ResourceHeader.append(restSessionHeaderSplit[j].replace('"', ''))
            
            labelIdx = ResourceHeader.index('label')
            if (restEndCount > 1):
                for j in xrange(1,restEndCount):
        
                    currRow = restResultsSplit[j]
                    currRowSplit = currRow.split(',')
#                    currRowCount = currRow.count(',')
                    currLabel = currRowSplit[labelIdx].replace('"', '')
        
                    restURL = self.Server +'data/projects/'+ self.Project +'/subjects/'+ self.Subject +'/experiments/'+ SubjectSessions[i] +'/resources/'+ currLabel +'/files?format=csv'
        
                    restResults = self.getURLString(restURL)
        
                    currRestResultsSplit = restResults.split('\n')
                    currRestEndCount = restResults.count('\n')
        
                    for k in xrange(1,currRestEndCount):
                        newRow = currRestResultsSplit[k]
                        currRowSplit = newRow.split(',')
                        FileNames.append(currRowSplit[0].replace('"', ''))
                        FileURIs.append(currRowSplit[2].replace('"', ''))
                        FileSessions.append(SubjectSessions[i].replace('"', ''))
                        FileLabels.append(currLabel.replace('"', ''))
        
        SubjectResources = { 'FileNames': FileNames, 'FileURIs': FileURIs, 'FileSessions': FileSessions, 'FileLabels': FileLabels }
        return SubjectResources
    #===============================================================================
    def getFileInfo( self, URL ):
        """Get info about a file on the server"""
        restRequest = urllib2.Request(URL)
        restRequest.add_header("Cookie", "JSESSIONID=" + self.SessionId);
         
        restRequest.get_method = lambda : 'HEAD'
        while (self.Timeout <= self.TimeoutMax):
            try:
                restConnHandle = urllib2.urlopen(restRequest, None)
                self.FileInfo = { 'ModDate': restConnHandle.info().getheader('Last-Modified'), 'Bytes': restConnHandle.info().getheader('Content-Length'), 'URL': URL }
            
            except HTTPError, e:
                if (e.code != 404):
                    self.Timeout += self.TimeoutStep
                    print 'HTTPError code: ' +str(e.code)+ '. Timeout increased to ' +str(self.Timeout)+' seconds for ' +URL
                else:
                    return '404 Error'
        
            if (self.FileInfo.get( 'Bytes' ) == None): 
                self.FileInfo[ 'Bytes' ] = '0'
            
            return self.FileInfo        

    #===============================================================================
    def getAssessorIDs( self ):
        """QC: Get assessor for subject and session"""

        IDs = list()
        SessionIDs = list()
        SessionLabels = list()
        Labels = list()
        XnatIDs = list()
        URIs = list()
        XsiType = list()
            
        restURL = self.Server + 'data/projects/' + self.Project + '/subjects/' + self.Subject + '/experiments/' + self.Session + '/assessors?format=csv'
        restResults = self.getURLString(restURL)
        
        if (restResults != '404 Error'):
        
    #        AssessorSubjectSessionsList, AssessorDataTypeList = self.getSubjectSessions()
            
            restResultsSplit = restResults.split('\n')
            restEndCount = restResults.count('\n')
            restSessionHeader = restResultsSplit[0]
            restSessionHeaderSplit = restSessionHeader.split(',')
            
            idIdx = restSessionHeaderSplit.index('"ID"')
            sessionIdIdx = restSessionHeaderSplit.index('"session_ID"')
            sessionLabelIdx = restSessionHeaderSplit.index('"session_label"')
            xnatidIdx = restSessionHeaderSplit.index('"xnat:imageassessordata/id"')
            labelIdx = restSessionHeaderSplit.index('"label"')
            uriIdx = restSessionHeaderSplit.index('"URI"')
            xsiTypeIdx = restSessionHeaderSplit.index('"xsiType"')
            
            for i in xrange(1, restEndCount):
                currRow = restResultsSplit[i]
                currRowSplit = currRow.split(',')
                
                # need to call getSubjectSessions before this....
                if (currRowSplit[xsiTypeIdx].replace('"', '').find(self.AssessorDataType) != -1):
                
                    IDs.append(currRowSplit[idIdx].replace('"', ''))
                    SessionIDs.append(currRowSplit[sessionIdIdx].replace('"', ''))
                    SessionLabels.append(currRowSplit[sessionLabelIdx].replace('"', ''))
                    Labels.append(currRowSplit[labelIdx].replace('"', ''))
                    XnatIDs.append(currRowSplit[xnatidIdx].replace('"', ''))
                    URIs.append(currRowSplit[uriIdx].replace('"', ''))
                    XsiType.append(currRowSplit[xsiTypeIdx].replace('"', ''))
                
        return Labels
    #===============================================================================
    def getAssessorOutputFile( self, AssessorIDs ):
        """QC: Get assessor output files as a list"""

        AssessorOutputFileURI = list()
        AssessorOutputFileSize = list()
        
        for i in xrange(0, len(AssessorIDs)):
            currID = AssessorIDs[i]
            
            restURL = self.Server + 'data/projects/' + self.Project + '/subjects/' + self.Subject + '/experiments/' + self.Session + '/assessors/' + currID + '/files?format=csv'
#            restURL = Server + 'data/projects/' + self.Project + '/subjects/' + self.Subject + '/experiments/' + self.Session + '/assessors/' + currID + '/files?format=csv'
            restResults = self.getURLString(restURL)
            
            #===================================================================
            # restResultsET = ET.fromstring(restResults)
            # for assessorOutput in restResultsET.findall('{http://nrg.wustl.edu/xnat}out'):
            #        for assessorFile in assessorOutput.findall('{http://nrg.wustl.edu/xnat}file'):
            #            tmpURL = self.Server + assessorFile.attrib.get('URI')
            #            tmpResults = self.getURLString(tmpURL)
            #            AssessorOutputFileURI.append(os.path.normpath(assessorFile.attrib.get('URI')))
            #===================================================================
            
            restResultsSplit = restResults.split('\n')
            restEndCount = restResults.count('\n')
            restSessionHeader = restResultsSplit[0]
            restSessionHeaderSplit = restSessionHeader.split(',')
            
            #===================================================================
            # Hopefully, digest will come out here?
            # digestIdx = restSessionHeaderSplit.index('"Digest"')
            # "Name","Size","URI","collection","file_tags","file_format","file_content","cat_ID"
            #===================================================================
            # Maybe useful later...
            #===================================================================
            # nameIdx = restSessionHeaderSplit.index('"Name"')
            # collectionIdx = restSessionHeaderSplit.index('"collection"')
            # fileTagsIdx = restSessionHeaderSplit.index('"file_tags"')
            # fileFormatIdx = restSessionHeaderSplit.index('"file_format"')
            # fileContentIdx = restSessionHeaderSplit.index('"file_content"')
            # catIdIdx = restSessionHeaderSplit.index('"cat_ID"')     
            #===================================================================
            sizeIdx = restSessionHeaderSplit.index('"Size"')
            uriIdx = restSessionHeaderSplit.index('"URI"')      
           
            for i in xrange(1, restEndCount):
                currRow = restResultsSplit[i]
                currRowSplit = currRow.split(',')
                AssessorOutputFileURI.append(currRowSplit[uriIdx].replace('"', ''))
                AssessorOutputFileSize.append(currRowSplit[sizeIdx].replace('"', ''))
                        
        return AssessorOutputFileURI
    #===============================================================================
    def getScanParms( self ):
        """HCP: Get scan parms from a scan, duh."""
        # def getScanParms(inputUser, inputPassword, inputProject, inputSubject, inputSession, inputScan):

        restURL = self.Server + 'data/projects/' +self.Project+ '/subjects/' +self.Subject+ '/experiments/' +self.Session+ '/scans/' +self.Scan
        xmlData = self.getURLString( restURL )
        parmsET = ET.fromstring(xmlData)
        
        acquisitionTime = parmsET.find('{http://nrg.wustl.edu/xnat}startTime').text
        
        scanParms = parmsET.find('{http://nrg.wustl.edu/xnat}parameters')
        sampleSpacing = scanParms.find('{http://nrg.wustl.edu/xnat}readoutSampleSpacing').text
        voxelResolution = scanParms.find('{http://nrg.wustl.edu/xnat}voxelRes').attrib
        orientation = scanParms.find('{http://nrg.wustl.edu/xnat}orientation').text
        FOV = scanParms.find('{http://nrg.wustl.edu/xnat}readoutSampleSpacing').text
        TR = scanParms.find('{http://nrg.wustl.edu/xnat}tr').text
        TE = scanParms.find('{http://nrg.wustl.edu/xnat}te').text
        flipAngle = scanParms.find('{http://nrg.wustl.edu/xnat}flip').text
        scanSequence = scanParms.find('{http://nrg.wustl.edu/xnat}scanSequence').text
        pixelBandwidth = scanParms.find('{http://nrg.wustl.edu/xnat}pixelBandwidth').text
        echoSpacing = scanParms.find('{http://nrg.wustl.edu/xnat}echoSpacing').text
        
        
        for addParms in scanParms.findall('{http://nrg.wustl.edu/xnat}addParam'):
            addParmsAttrib = addParms.attrib
            
            if (addParmsAttrib.get('name') == 'Siemens GRADSPEC alShimCurrent'):
                alShimCurrent = addParms.text
                
            if (addParmsAttrib.get('name') == 'Siemens GRADSPEC lOffset'):
                LinOffset = addParms.text
        
        scanParms = { 'SampleSpacing': sampleSpacing, 'alShimCurrent': alShimCurrent, 'LinearOffset':  LinOffset, 'AcquisitionTime': acquisitionTime, 'VoxelResolution': voxelResolution, \
                            'Orientation': orientation, 'FOV': FOV, 'TR': TR, 'TE': TE, 'FlipAngle': flipAngle, 'ScanSequence': scanSequence, 'PixelBandwidth': pixelBandwidth, 'EchoSpacing': echoSpacing }
        return scanParms
    #===============================================================================    
    def getScanMeta( self ):
        """Get Scan ID, Type, Series, Quality, and XNAT ID for a given subject and session"""

        Names = list()
        Sizes = list()
        URIs = list()
        Collections = list()
        FileTags = list()
        FileFormats = list()
        FileContents = list()

        if (not self.Scan):
            print 'Scan not defined...'
        else:
            restURL = self.Server + 'data/projects/' + self.Project + '/subjects/' + self.Subject + '/experiments/' + self.Session + '/scans/' + self.Scan + '/files?format=csv&columns=ID,type,series_description,quality,xnat:mrSessionData/id'
        
        restResults = self.getURLString(restURL)
        
        restResultsSplit = restResults.split('\n')
        restEndCount = restResults.count('\n')
        restSessionHeader = restResultsSplit[0]
        restSessionHeaderSplit = restSessionHeader.split(',')
        
        # ['"Name"', '"Size"', '"URI"', '"collection"', '"file_tags"', '"file_format"', '"file_content"', '"cat_ID"']
        nameIdx = restSessionHeaderSplit.index('"Name"')
        sizeIdx = restSessionHeaderSplit.index('"Size"')
        uriIdx = restSessionHeaderSplit.index('"URI"')
        collectionIdx = restSessionHeaderSplit.index('"collection"')
        fileTagsIdx = restSessionHeaderSplit.index('"file_tags"')
        fileFormatIdx = restSessionHeaderSplit.index('"file_format"')
        fileContentIdx = restSessionHeaderSplit.index('"file_content"')
        
        for j in xrange(1, restEndCount):
            currRow = restResultsSplit[j]
            
            currRowSplit = currRow.split(',')
    
            Names.append(currRowSplit[nameIdx].replace('"', ''))
            Sizes.append(currRowSplit[sizeIdx].replace('"', ''))
            URIs.append(currRowSplit[uriIdx].replace('"', ''))
            Collections.append(currRowSplit[collectionIdx].replace('"', ''))
            FileTags.append(currRowSplit[fileTagsIdx].replace('"', ''))
            FileFormats.append(currRowSplit[fileFormatIdx].replace('"', ''))
            FileContents.append(currRowSplit[fileContentIdx].replace('"', ''))
            
        ScanMeta = {'Names': Names, 'Bytes': Sizes, 'URIs': URIs, 'Collections': Collections}
        return ScanMeta
    
#===============================================================================
# WRITE
#===============================================================================
class writeHCP:
    """HCP Write Class"""
    def __init__( self, DestinationDir  ):
        self.DestinationDir = DestinationDir
        
        
    def writeFileFromURL( self, FileURI ):
    
        FileURIList = FileURI.split(',')
        WriteCode = True
        if (self.DestinationDir[-1] != os.sep):
            self.DestinationDir = self.DestinationDir + os.sep

        if not os.path.exists(self.DestinationDir):
            os.makedirs(self.DestinationDir)
            
        for i in xrange(len(FileURIList)):
            currURI = FileURIList[i]
            currFileName = os.path.basename(currURI)
                
            #===================================================================
            # need getHCP class here...
            #===================================================================
            fileURL = self.Server + currURI
            fileInfo = self.getFileInfo(fileURL)
            fileResults = self.getURLString(fileURL)
            
#                WriteTotal = fileInfo.get('Bytes') + WriteTotal
            
            if (fileInfo.get('Bytes') != str(len(fileResults))):
                print 'WARNING: Expected ' +fileInfo.get('Bytes')+ ' bytes and downloaded ' +str(len(fileResults))+ ' bytes for file ' +currFileName
                WriteCode = False
            else:
                with open(self.DestinationDir + currFileName, 'wb') as outputFileId:
                    writeCode = outputFileId.write(fileResults)
                    if (self.Verbose):
                        print 'File: ' +self.DestinationDir+currFileName+ '  Write Code: ' +str(writeCode)
                    outputFileId.flush()
                    os.fsync(outputFileId)
                    outputFileId.close()
                    
                # check file size after write...
                writeFileSize = os.path.getsize(self.DestinationDir + currFileName)
                if (fileInfo.get('Bytes') != str(writeFileSize)):
                    print 'WARNING: WROTE ' +str(len(fileResults))+ ' bytes but expected ' +str(writeFileSize)+ ' bytes for file ' +self.DestinationDir+currFileName
                    WriteCode = False
                    
        return WriteCode
#===============================================================================
# END CLASS DEFs
#===============================================================================

if __name__ == "__main__":
    getHCP(sys.argv[1], sys.argv[2], sys.argv[3])





