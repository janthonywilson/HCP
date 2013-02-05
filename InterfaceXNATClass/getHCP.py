'''
Created on 2012-12-19

@author: jwilso01
'''


import base64
# multiplatform stuff...
import os
import sys
# Web stuff...
import socket
import urllib
import urllib2
from ssl import SSLError
from urllib2 import URLError, HTTPError


#===============================================================================
# CLASSES
#===============================================================================
class getHCP:
    """intraDB Interfacing Class"""
    def __init__( self, User, Password, Timeout, TimeoutMax, TimeoutStep ):
        self.User = User
        self.Password = Password
        self.Timeout = Timeout
        self.TimeoutMax = TimeoutMax
        self.TimeoutStep = TimeoutStep
        
        
        self.Session = ''
        self.Subject = ''
        self.ProjectIDs = []
        self.ProjectNames = []
        self.Subjects = []
        self.FileInfo = {}
        self.SubjectSessions = []
        self.SubjectSessionsUniq = []
        self.SessionParms = {}
        self.SessionType = []
        self.Sessions = []
        
        self.ProjectNames = []
        self.ProjectIDs = []
        self.Projects = ()
        
        
        self.SessionId = self.fGetSessionId()
#        self.ProjectNames, self.ProjectIDs = self.fGetProjects()
#        self.Projects = self.fGetProjects()
        self.Subjects = self.fGetSubjects()
        
        self.fGetAssessorIDs
    #===============================================================================
    def fGetSessionId( self ):
        """Get session id for getHCP session spawn"""
        restURL = self.Server + '/data/JSESSION'
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
    def fGetFileInfo( self, URL ):
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
    def fGetProjects( self ):
        """Get a list of project names from XNAT instance"""
        
        restURL = self.Server + 'data/projects?format=csv'
        restResults = self.fGetURLString(restURL)
        
        restResultsSplit = restResults.split('\n')
        restEndCount = restResults.count('\n')
        restSessionHeader = restResultsSplit[0]
        restSessionHeaderSplit = restSessionHeader.split(',')
        
        projectNameIdx = restSessionHeaderSplit.index('"name"')
        projectIdx = restSessionHeaderSplit.index('"ID"')
        
        for i in xrange(1, restEndCount):
            currRow = restResultsSplit[i]
            
            currRowSplit = currRow.split(',')
            currProjectName = currRowSplit[projectNameIdx].replace('"', '')
            currProjectId = currRowSplit[projectIdx].replace('"', '')
            
            self.ProjectNames.append(currProjectName)
            self.ProjectIDs.append(currProjectId)
            
        return self.ProjectNames, self.ProjectIDs
        
        
    #===============================================================================    
    def fGetSubjectSessions( self ):
        """Get all sessions and session types for a given subject"""
        SubjectSessions = list()
        SubjectSessionUniq = list()
        SubjectSessionsType = list()
        
        restURL = self.Server + 'data/projects/' + self.Project + '/subjects/' + self.Subject + '/experiments?format=csv&columns=ID,label'
        restResults = self.fGetURLString(restURL)
        
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
                
                if (currSession.find('fnc') != -1) or (currSession.find('str') != -1) or (currSession.find('diff') != -1) or (currSession.find('xtr') != -1):
                    if (currSession.find('xtr') != -1):
                        SubjectSessions.append(currSession)
                        SessionTypeList = self.fGetSessionMeta( )[1]
                        if ('T1w' in SessionTypeList):
                            SubjectSessionsType.append('strc')
                        elif ('dMRI' in SessionTypeList):
                            SubjectSessionsType.append('diff')
                        elif ('tfMRI' in SessionTypeList):
                            SubjectSessionsType.append('fnc')
                        else:
                            SubjectSessionsType.append('unknown')
                        
                    else:
                        SubjectSessionsType.append(currSession[currSession.find('_')+1:])
                        SubjectSessions.append(currSession)
                
            SubjectSessionUniq.append(currSessionUniq)

        return SubjectSessions, SubjectSessionsType
    #===============================================================================
    def fGetSubjects( self ):
        """Get all subjects for a given project"""
        Subjects = list()
        
        restURL = self.Server + 'data/projects/' + self.Project + '/subjects?format=csv'
        restResults = self.fGetURLString(restURL)
        
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
    def fGetSessionMeta( self ):
        """Get ID, Type, Series, Quality, and XNAT ID for a given subject and session"""
#        SessionParms = {}
        idList = list()
        typeList = list()
        seriesList = list()
        qualityList = list()
        xnatidList = list()

        restURL = self.Server + 'data/projects/' + self.Project + '/subjects/' + self.Subject + '/experiments/' + self.Session + '/scans?format=csv&columns=ID,type,series_description,quality,xnat:mrSessionData/id'
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
    def fGetSessionQuality( self ):
        """QC: Get Session, Subject, Scan Type, Quality of all data"""

        Quality = list()
        Series = list()
        ScanIds = list()
        Sessions = list()
        ScanType = list()
        
        if (not self.Subjects):
            self.Subjects = self.fGetSubjects( )
        else:
            self.Subjects = self.Subject.split()
            
        for i in xrange(0, len(self.Subjects)):
            self.Subject = self.Subjects[i]
            
            if (not self.Session) and (not self.Sessions):
                self.Sessions = self.fGetSubjectSessions( )[0]
            else:
                self.Sessions = self.Session.split()
            
            for j in xrange(0, len(self.Sessions)):
                
                self.Session = self.Sessions[j]
                
                idList, typeList, seriesList, qualityList, xnatidList = self.fGetSessionMeta( )
                
                Quality.extend(qualityList)
                Series.extend(seriesList)
                ScanIds.extend(idList)
                ScanType.extend(typeList)
                
                Sessions.extend( self.Session.split(',') * len(idList))
                
        return Quality, ScanIds, Series, Sessions, ScanType
    #===============================================================================
    def fGetAssessorIDs( self ):
        """QC: Get assessor for subject and session"""

        IDs = list()
        SessionIDs = list()
        SessionLabels = list()
        Labels = list()
        XnatIDs = list()
        URIs = list()
        XsiType = list()
        
#        if (str.lower(self.SessionType).find('str') != -1):
#            Session = '%s_%s' % (self.Subject, 'strc')
#        elif (str.lower(self.SessionType).find('fnc') != -1):
#            Session = '%s_%s' % (self.Subject, 'fncXXX')
#            print 'WARNING: Need to code this for fnca or fncb...'
#        elif (str.lower(self.SessionType).find('dif') != -1):
#            Session = '%s_%s' % (self.Subject, 'diff')
            
        restURL = self.Server + 'data/projects/' + self.Project + '/subjects/' + self.Subject + '/experiments/' + self.Session + '/assessors?format=csv'
        restResults = self.fGetURLString(restURL)
        
        if (restResults != '404 Error'):
        
    #        AssessorSubjectSessionsList, AssessorDataTypeList = self.fGetSubjectSessions()
            
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
                
                # need to call fGetSubjectSessions before this....
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
    def fGetAssessorOutputFile( self, AssessorIDs ):
        """QC: Get assessor output files as a list"""

        AssessorOutputFileURI = list()
        AssessorOutputFileSize = list()
        
        for i in xrange(0, len(AssessorIDs)):
            currID = AssessorIDs[i]
            
            restURL = self.Server + 'data/projects/' + self.Project + '/subjects/' + self.Subject + '/experiments/' + self.Session + '/assessors/' + currID + '/files?format=csv'
#            restURL = Server + 'data/projects/' + self.Project + '/subjects/' + self.Subject + '/experiments/' + self.Session + '/assessors/' + currID + '/files?format=csv'
            restResults = self.fGetURLString(restURL)
            
            #===================================================================
            # restResultsET = ET.fromstring(restResults)
            # for assessorOutput in restResultsET.findall('{http://nrg.wustl.edu/xnat}out'):
            #        for assessorFile in assessorOutput.findall('{http://nrg.wustl.edu/xnat}file'):
            #            tmpURL = self.Server + assessorFile.attrib.get('URI')
            #            tmpResults = self.fGetURLString(tmpURL)
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
            nameIdx = restSessionHeaderSplit.index('"Name"')
            sizeIdx = restSessionHeaderSplit.index('"Size"')
            uriIdx = restSessionHeaderSplit.index('"URI"')
            collectionIdx = restSessionHeaderSplit.index('"collection"')
            fileTagsIdx = restSessionHeaderSplit.index('"file_tags"')
            fileFormatIdx = restSessionHeaderSplit.index('"file_format"')
            fileContentIdx = restSessionHeaderSplit.index('"file_content"')
            catIdIdx = restSessionHeaderSplit.index('"cat_ID"')           
           
            for i in xrange(1, restEndCount):
                currRow = restResultsSplit[i]
                currRowSplit = currRow.split(',')
                AssessorOutputFileURI.append(currRowSplit[uriIdx].replace('"', ''))
                AssessorOutputFileSize.append(currRowSplit[sizeIdx].replace('"', ''))
                        
        return AssessorOutputFileURI
    
    #===============================================================================
    def fWriteFileFromURL( self, FileURIList, IncludeList ):
    
            WriteCode = True
            for i in xrange(len(FileURIList)):
                currURI = FileURIList[i]
                currFileName = os.path.basename(currURI)
                
                # build list to test for desirable IncludeList values
                IncludeBool = list()
                for j in xrange(0, len(IncludeList)):
                    
                    if ( currFileName.find( IncludeList[j] ) != -1 ):
                        IncludeBool.append(True)
                    else:
                        IncludeBool.append(False)
                        
                if (IncludeBool.count(False) == 0):
                    
                    fileURL = self.Server + currURI
                    fileInfo = self.fGetFileInfo(fileURL)
                    fileResults = self.fGetURLString(fileURL)
                    
    #                WriteTotal = fileInfo.get('Bytes') + WriteTotal
                    
                    if (fileInfo.get('Bytes') != str(len(fileResults))):
                        print 'WARNING: Expected ' +fileInfo.get('Bytes')+ ' bytes and downloaded ' +str(len(fileResults))+ ' bytes for file ' +currFileName
                        WriteCode = False
                    else:
                        with open(self.DestinationDir + currFileName, 'wb') as outputFileObj:
                            writeCode = outputFileObj.write(fileResults)
                            if (self.Verbose):
                                print 'File: ' +self.DestinationDir+currFileName+ '  Write Code: ' +str(writeCode)
                            outputFileObj.flush()
                            os.fsync(outputFileObj)
                            outputFileObj.close()
                            
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
    getHCP(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])





