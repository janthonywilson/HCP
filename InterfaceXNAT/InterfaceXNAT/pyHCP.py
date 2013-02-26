'''
Created on 2012-12-19

@author: jwilso01
'''

# multiplatform system stuff...
import os
import sys
import errno
import subprocess
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
class pyHCP(object):
    """Main HCP Interfacing Class"""
    def __init__( self, User, Password, Server ):
        super(pyHCP, self).__init__()
        self.User = User
        self.Password = Password
        self.Server = self.cleanServer(Server)
        
    def cleanServer(self, Server):
        Server.strip()
        if (Server[-1] != '/'):
            Server = Server + '/'
        if (Server.find('http') == -1):
            Server = 'https://' + Server
        self.Server = Server
        return self.Server
    
    def getServer(self):
        return self.Server
        
class getHCP(pyHCP):
    """HCP Interfacing Class for GETs"""

    def __init__( self, pyHCP ):
        
        #=======================================================================
        # need to add a check for project existence.../
        #=======================================================================
        self.User = pyHCP.User
        self.Password = pyHCP.Password
        self.Server = pyHCP.Server
        
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
        self.SubjectResourcesMeta = {}
        
        self.Scan = ''
        self.Scans = []
        
        self.Resource = ''
        self.Resources = []
        
        self.FileInfo = {}
        self.ScanMeta = {}
        self.ResourceMeta = {}
        
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

        while (self.Timeout <= self.TimeoutMax):
            try:
                connHandle = urllib2.urlopen(Request, None, self.Timeout)
                break
            except URLError, e:
                try:
                    code = e.code
                    if  (code != 401):
                        self.Timeout += self.TimeoutStep
                        print 'URLError code: ' +str(e.reason)+ '. Timeout increased to ' +str(self.Timeout)+' seconds for JSESSION cookie...'
                    else:
                        print 'URLError code: ' +str(e.reason)+ '. getSessionId Failed with wrong password.'
                        sys.exit(401)
                except:
                    print 'URL: %s failed with code: %s ' % (URL, e.code)
                    sys.exit()
            except SSLError, e:
                self.Timeout += self.TimeoutStep
                print 'SSLError code: ' +str(e.message)+ '. Timeout increased to ' +str(self.Timeout)+' seconds for ' +URL
            except socket.timeout:
                self.Timeout += self.TimeoutStep
                print 'Socket timed out. Timeout increased to ' +str(self.Timeout)+ ' seconds for ' +URL
            
                
                
        self.SessionId = connHandle.read()
        return self.SessionId
    #===============================================================================
    def getURLString( self, URL ):
        """Get URL results as a string"""
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
                    print 'READ HTTPError code: ' +str(e.code)+ '. File read timeout for ' +str(self.Timeout)+ ' seconds for ' +URL
                except URLError, e:
                    print 'READ URLError code: ' +str(e.reason)+ '. File read timeout for ' +str(self.Timeout)+' seconds for ' +URL
                except SSLError, e:
                    print 'READ SSLError code: ' +str(e.message)+ '. File read timeout for ' +str(self.Timeout)+' seconds for ' +URL
                except socket.timeout:
                    print 'READ Socket timed out. File read timeout for ' +str(self.Timeout)+ ' seconds for ' +URL
                    
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
        
        AllProject = self.getProjects()
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
                        
                        self.Session = currSession
                        SessionTypeList = self.getSessionMeta( ).get('Types')
                        if ('T1w' in SessionTypeList) and ('T2w' in SessionTypeList):
                            SubjectSessionsType.append('strc')
                            SubjectSessionsID.append(currSession)
                        if ('dMRI' in SessionTypeList):
                            SubjectSessionsType.append('diff')
                            SubjectSessionsID.append(currSession)
                        if ('tfMRI' in SessionTypeList):
                            SubjectSessionsType.append('task')
                            SubjectSessionsID.append(currSession)
                        if ('rfMRI' in SessionTypeList):
                            SubjectSessionsType.append('rest')
                            SubjectSessionsID.append(currSession)
                            
                    else:
                        if (currSession.find('fnc') != -1): SubjectSessionsType.append('fnc')
                        elif (currSession.find('strc') != -1): SubjectSessionsType.append('strc')
                        elif (currSession.find('diff') != -1): SubjectSessionsType.append('diff')
                        else: SubjectSessionsType.append('unknown')
                        SubjectSessionsID.append(currSession)
                
            SubjectSessionUniq.append(currSessionUniq)

            SubjectSessions = {'Sessions': SubjectSessionsID, 'Types': SubjectSessionsType}
            return SubjectSessions
        else:
            print 'ERROR(getSubjectSessions()): No subject sessions found for subject %s under project %s' % (self.Subject, self.Project)
            sys.exit(-1)
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
            print 'No session for getSessionMeta()...'
            sys.exit(-1)
            
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
    def getSubjectResourcesMeta(self):
        
        ResourceHeader = list()
        FileNames = list()
        FileURIs = list()
        FileSessions = list()
        FileLabels = list()
        FileTags = list()
        FileFormats = list()
        FileContents = list()
        FileReadable = list()
        FilePath = list()
        
        SubjectSessions = list(set(self.getSubjectSessions().get('Sessions')))
    
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
        
                    #===========================================================
                    # all this nonsense should be replaced with a call to getResourceMeta()
                    #===========================================================
                    restURL = self.Server +'data/projects/'+ self.Project +'/subjects/'+ self.Subject +'/experiments/'+ SubjectSessions[i] +'/resources/'+ currLabel +'/files?format=csv'
                    restResults = self.getURLString(restURL)
        
                    currRestResultsSplit = restResults.split('\n')
                    currRestEndCount = restResults.count('\n')
                    currRestSessionHeader = currRestResultsSplit[0]
                    currRestSessionHeaderSplit = currRestSessionHeader.split(',')
                    
                    nameIdx = currRestSessionHeaderSplit.index('"Name"')
                    sizeIdx = currRestSessionHeaderSplit.index('"Size"')
                    uriIdx = currRestSessionHeaderSplit.index('"URI"')
                    collectionIdx = currRestSessionHeaderSplit.index('"collection"')
                    fileTagsIdx = currRestSessionHeaderSplit.index('"file_tags"')
                    fileFormatIdx = currRestSessionHeaderSplit.index('"file_format"')
                    fileContentIdx = currRestSessionHeaderSplit.index('"file_content"')
        
                    for k in xrange(1,currRestEndCount):
                        newRow = currRestResultsSplit[k]
                        currRowSplit = newRow.split(',')
                        FileNames.append(currRowSplit[nameIdx].replace('"', ''))
                        FileURIs.append(currRowSplit[uriIdx].replace('"', ''))
                        FileTags.append(currRowSplit[fileTagsIdx].replace('"', ''))
                        FileFormats.append(currRowSplit[fileFormatIdx].replace('"', ''))
                        FileContents.append(currRowSplit[fileContentIdx].replace('"', ''))
                        
                        FileSessions.append(SubjectSessions[i].replace('"', ''))
                        FileLabels.append(currLabel.replace('"', ''))
                        
                    # do the path query...    
                    restURL = self.Server + 'data/projects/' + self.Project + '/subjects/' + self.Subject + '/experiments/' + SubjectSessions[i] + '/resources/' + currLabel + '/files?format=csv&locator=absolutePath'
                    newRestResults = self.getURLString(restURL)
        
                    newRestResultsSplit = newRestResults.split('\n')
                    newRestEndCount = newRestResults.count('\n')
                    newRestSessionHeader = newRestResultsSplit[0]
                    newRestSessionHeaderSplit = newRestSessionHeader.split(',')
                    
                    # ['"Name"', '"Size"', '"URI"', '"collection"', '"file_tags"', '"file_format"', '"file_content"', '"cat_ID"']
                    pathIdx = newRestSessionHeaderSplit.index('"absolutePath"')
                    
                    for k in xrange(1, newRestEndCount):
                        currRow = newRestResultsSplit[k]
                        currRowSplit = currRow.split(',')
                        FilePath.append(currRowSplit[pathIdx].replace('"', ''))
                        try:
                            FileObj = open(FilePath[-1], 'r')
                            # if readable:
                            FileReadable.append(True)
                        except IOError, e:
                            if self.Verbose:
                                print 'getSubjectResources(): File read error number: %s, error code: %s, and error message: %s' % (e.errno, errno.errorcode[e.errno], os.strerror(e.errno))
                            FileReadable.append(False)
        
        SubjectResources = { 'Name': FileNames, 'URI': FileURIs, 'Session': FileSessions, 'Label': FileLabels, 'Content': FileContents, 'Format': FileFormats, 'Path': FilePath, 'Readable': FileReadable }
        return SubjectResources
    #===============================================================================  
    def getResourceMeta(self):
        """Get file info about a given resource"""
        
        Names = list()
        Sizes = list()
        URIs = list()
        Collections = list()
        FilePath = list()
        FileTags = list()
        FileFormats = list()
        FileContents = list()
        FileReadable = list()
        
        # https://hcpx-dev-cuda00.nrg.mir/data/projects/HCP_Q1/subjects/100307/experiments/100307_3T/resources/tfMRI_WM_RL_unproc/files?format=csv&columns=ID,type,series_description,quality,xnat:mrSessionData/id
        restURL = self.Server + 'data/projects/' + self.Project + '/subjects/' + self.Subject + '/experiments/' + self.Session + '/resources/' + self.Resource + '/files?format=csv&columns=ID,type,series_description,quality,xnat:mrSessionData/id,file_tags,file_format,file_content'
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
        
        # do the path query...    
        restURL = self.Server + 'data/projects/' + self.Project + '/subjects/' + self.Subject + '/experiments/' + self.Session + '/resources/' + self.Resource + '/files?format=csv&locator=absolutePath'
        restResults = self.getURLString(restURL)
        
        restResultsSplit = restResults.split('\n')
        restEndCount = restResults.count('\n')
        restSessionHeader = restResultsSplit[0]
        restSessionHeaderSplit = restSessionHeader.split(',')
        
        # ['"Name"', '"Size"', '"URI"', '"collection"', '"file_tags"', '"file_format"', '"file_content"', '"cat_ID"']
        pathIdx = restSessionHeaderSplit.index('"absolutePath"')
        
        for j in xrange(1, restEndCount):
            currRow = restResultsSplit[j]
            currRowSplit = currRow.split(',')
            FilePath.append(currRowSplit[pathIdx].replace('"', ''))
            try:
                FileObj = open(FilePath[-1], 'r')
                # if readable:
                FileReadable.append(True)
            except IOError, e:
                if self.Verbose:
                    print 'getScanMeta(): File read error number: %s, error code: %s, and error message: %s' % (e.errno, errno.errorcode[e.errno], os.strerror(e.errno))
                FileReadable.append(False)
                
        ResourceMeta = {'Name': Names, 'Bytes': Sizes, 'URI': URIs, 'Path': FilePath, 'Readable': FileReadable, 'Label': Collections, 'Format': FileFormats, 'Contents': FileContents}
        return ResourceMeta
    #===============================================================================
    def getFileInfo( self, URL ):
        """Get mod-date, size, and URL for a file on the server"""

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

        if (self.Server.find('intradb') > 0):
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
        else:
            print 'ERROR: Assessor data only on intradb.humanconnectome.org.'
            return -1
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
        
        flipAngle = scanParms.find('{http://nrg.wustl.edu/xnat}flip').text
        scanSequence = scanParms.find('{http://nrg.wustl.edu/xnat}scanSequence').text
        pixelBandwidth = scanParms.find('{http://nrg.wustl.edu/xnat}pixelBandwidth').text
        
        # alt parms that are not present in all scan types...
        try:
            readoutDirection = scanParms.find('{http://nrg.wustl.edu/xnat}readoutDirection').text
        except:
            readoutDirection = 'NA'
        try:
            echoSpacing = scanParms.find('{http://nrg.wustl.edu/xnat}echoSpacing').text
        except:
            echoSpacing = 'NA'
        try:
            peDirection = scanParms.find('{http://nrg.wustl.edu/xnat}peDirection').text
        except:
            peDirection = 'NA'
        try:
            shimGroup = scanParms.find('{http://nrg.wustl.edu/xnat}shimGroup').text
        except:
            shimGroup = 'NA'
        try:
            seFieldMapGroup = scanParms.find('{http://nrg.wustl.edu/xnat}seFieldMapGroup').text
        except:
            seFieldMapGroup = 'NA'
        try:
            deltaTE = scanParms.find('{http://nrg.wustl.edu/xnat}deltaTE').text
        except:
            deltaTE = 'NA'
        try:
            TE = scanParms.find('{http://nrg.wustl.edu/xnat}te').text
        except:
            TE = 'NA'
        
        
        for addParms in scanParms.findall('{http://nrg.wustl.edu/xnat}addParam'):
            addParmsAttrib = addParms.attrib
            
            if (addParmsAttrib.get('name') == 'Siemens GRADSPEC alShimCurrent'):
                alShimCurrent = addParms.text
                
            if (addParmsAttrib.get('name') == 'Siemens GRADSPEC lOffset'):
                LinOffset = addParms.text
        
        scanParms = { 'SampleSpacing': sampleSpacing, 'alShimCurrent': alShimCurrent, 'LinearOffset':  LinOffset, 'AcquisitionTime': acquisitionTime, 'VoxelResolution': voxelResolution, 'Orientation': orientation, \
                            'FOV': FOV, 'TR': TR, 'TE': TE, 'FlipAngle': flipAngle, 'ScanSequence': scanSequence, 'PixelBandwidth': pixelBandwidth, 'ReadoutDirection': readoutDirection, 'EchoSpacing': echoSpacing, \
                            'PhaseEncodingDir': peDirection, 'ShimGroup': shimGroup, 'SEFieldmapGroup': seFieldMapGroup, 'DeltaTE': deltaTE  }
        return scanParms
    #===============================================================================    
    def getScanMeta( self ):
        """Get Scan ID, Type, Series, Quality, and XNAT ID for a given subject and session"""

        Names = list()
        Sizes = list()
        URIs = list()
        Collections = list()
        FilePath = list()
        FileTags = list()
        FileFormats = list()
        FileContents = list()
        FileReadable = list()

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
        
        # do the path query...    
        restURL = self.Server + 'data/projects/' + self.Project + '/subjects/' + self.Subject + '/experiments/' + self.Session + '/scans/' + self.Scan + '/files?format=csv&locator=absolutePath'
        restResults = self.getURLString(restURL)
        
        restResultsSplit = restResults.split('\n')
        restEndCount = restResults.count('\n')
        restSessionHeader = restResultsSplit[0]
        restSessionHeaderSplit = restSessionHeader.split(',')
        
        # ['"Name"', '"Size"', '"URI"', '"collection"', '"file_tags"', '"file_format"', '"file_content"', '"cat_ID"']
        pathIdx = restSessionHeaderSplit.index('"absolutePath"')
        
        for j in xrange(1, restEndCount):
            currRow = restResultsSplit[j]
            currRowSplit = currRow.split(',')
            FilePath.append(currRowSplit[pathIdx].replace('"', ''))
            try:
                FileObj = open(FilePath[-1], 'r')
                # if readable:
                FileReadable.append(True)
            except IOError, e:
                if self.Verbose:
                    print 'getScanMeta(): File read error number: %s, error code: %s, and error message: %s' % (e.errno, errno.errorcode[e.errno], os.strerror(e.errno))
                FileReadable.append(False)
            
        ScanMeta = {'Names': Names, 'Bytes': Sizes, 'URI': URIs, 'Path': FilePath, 'Readable': FileReadable, 'Collections': Collections, 'Format': FileFormats, 'Content': FileContents }
        return ScanMeta

#===============================================================================
# WRITE
#===============================================================================
class writeHCP(getHCP):
    """HCP Write Class"""
    def __init__( self, getHCP, DestinationDir  ):
        self.DestinationDir = DestinationDir
        self.Server = getHCP.Server
        self.SessionId = getHCP.SessionId
        self.Timeout = getHCP.Timeout
        self.TimeoutMax = getHCP.TimeoutMax
        self.TimeoutStep = getHCP.TimeoutStep
        self.FileInfo = getHCP.FileInfo
        self.Verbose = getHCP.Verbose
        
        self.BytesStream = list()
        self.BytesWrite = list()
    #===============================================================================
    def getURLString(self, fileURL):
        return super(writeHCP, self).getURLString(fileURL)
    #===============================================================================
    def getFileInfo(self, fileURL):
        return super(writeHCP, self).getFileInfo(fileURL)
    #===============================================================================
    def writeFileFromURL( self, getHCP, FileURI, FileName ):
    
        try:
            FileURIList = FileURI.split(',')
        except:
            FileURIList = FileURI
        try:
            FileNameList = FileName.split(',')
        except:
            FileNameList = FileName

        
        WriteCode = True
        if (self.DestinationDir[-1] != os.sep):
            self.DestinationDir = self.DestinationDir + os.sep

        if not os.path.exists(self.DestinationDir):
            os.makedirs(self.DestinationDir)
            
        for i in xrange(len(FileURIList)):
            currURI = FileURIList[i]
            currURISplit = currURI.split('/')
            currFileNameIdx = currURISplit.index(os.path.basename(currURI))
            currResrouceRootIdx = currURISplit.index('files') + 1
            
            if (currFileNameIdx > currResrouceRootIdx):
#                print self.DestinationDir +os.sep.join(currURISplit[currResrouceRootIdx:currFileNameIdx])+os.sep
                newDestinationDir = self.DestinationDir +os.sep.join(currURISplit[currResrouceRootIdx:currFileNameIdx])+os.sep
                
                if not os.path.exists(newDestinationDir):
                    os.makedirs(newDestinationDir)
            else:
                newDestinationDir = self.DestinationDir

                
            if (FileName == None):
                currFileName = os.path.basename(currURI)
            else:
                currFileName = FileNameList[i]
                
            fileURL = self.Server + currURI
            fileInfo = self.getFileInfo(fileURL)
            fileResults = getHCP.getURLString(fileURL)
            self.BytesStream.append(len(fileResults))
            
#                WriteTotal = fileInfo.get('Bytes') + WriteTotal
            
            if (fileInfo.get('Bytes') != str(len(fileResults))):
                print 'WARNING: Expected ' +fileInfo.get('Bytes')+ ' bytes and downloaded ' +str(len(fileResults))+ ' bytes for file ' +currFileName
                WriteCode = False
            else:
                with open(newDestinationDir + currFileName, 'wb') as outputFileObj:
                    writeCode = outputFileObj.write(fileResults)
                    if (self.Verbose):
                        print 'File: ' +newDestinationDir+currFileName+ '  Write Code: ' +str(writeCode)
                    outputFileObj.flush()
                    os.fsync(outputFileObj)
                    outputFileObj.close()
                    
                # check file size after write...
                writeFileSize = os.path.getsize(newDestinationDir + currFileName)
                self.BytesWrite.append(writeFileSize)
                
                if (fileInfo.get('Bytes') != str(writeFileSize)):
                    print 'WARNING: WROTE ' +str(len(fileResults))+ ' bytes but expected ' +str(writeFileSize)+ ' bytes for file ' +newDestinationDir+currFileName
                    WriteCode = False
                    
        return WriteCode
    #===============================================================================
    def writeFileFromPath( self, FilePathName, FileName ):
        
        try:
            FilePathNameList = FilePathName.split(',')
        except:
            FilePathNameList = FilePathName
            
        try:
            FileNameList = FileName.split(',')
        except:
            FileNameList = FileName
       
        WriteCode = True
        if (self.DestinationDir[-1] != os.sep):
            self.DestinationDir = self.DestinationDir + os.sep

        if not os.path.exists(self.DestinationDir):
            os.makedirs(self.DestinationDir)
        
        for i in xrange(len(FilePathNameList)):
            currFilePathName = FilePathNameList[i]
            
            currFilePathNameSplit = currFilePathName.split('/')
            currFileNameIdx = currFilePathNameSplit.index(os.path.basename(currFilePathName))
            currResrouceRootIdx = currFilePathNameSplit.index('RESOURCES') + 1
            
            if (currFileNameIdx > currResrouceRootIdx):
                print self.DestinationDir +os.sep.join(currFilePathNameSplit[currResrouceRootIdx:currFileNameIdx])+os.sep
                newDestinationDir = self.DestinationDir +os.sep.join(currFilePathNameSplit[currResrouceRootIdx:currFileNameIdx])+os.sep
                
                if not os.path.exists(newDestinationDir):
                    os.makedirs(newDestinationDir)
            else:
                newDestinationDir = self.DestinationDir
            
            if (FileName == None):
                currFileName = os.path.basename(currFilePathName)
            else:
                currFileName = FileNameList[i]
            
            WriteCode = subprocess.call(('cp %s %s' % (currFilePathName, newDestinationDir + currFileName)), shell=True)
        
        return WriteCode
#===============================================================================
# END CLASS DEFs
#===============================================================================
    
if __name__ == "__main__":
    pyHCP(sys.argv[1], sys.argv[2], sys.argv[3])





