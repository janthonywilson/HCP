'''
Created on Jan 4, 2013

@author: Tony
'''
import os
import sys
import time
import socket
import argparse

from pyHCP import getHCP

sTime = time.time()

#===============================================================================
# PARSE INPUT
#===============================================================================
parser = argparse.ArgumentParser(description="Alpha program to pull NIFTI data from XNAT and put it somewhere...")

parser.add_argument("-W", "--server", dest="Server", default="https://intradb.humanconnectome.org", type=str, help="specify which server to connect to")
parser.add_argument("-u", "--username", dest="User", type=str, help="username must be specified")
parser.add_argument("-p", "--password", dest="Password", type=str, help="password must be specified")

parser.add_argument("-P", "--project", dest="Project", default="HCP_Phase2", type=str, help="specify project")
parser.add_argument("-S", "--subject", dest="Subject", default="100307", type=str, help="specify subject of interest")
parser.add_argument("-T", "--type", dest="DataType", default="Structural", type=str, help="specify datatype")

parser.add_argument("-d", "--source_dir", dest="SourceDir", type=str, help="specify the XNAT source directory under the parent")
parser.add_argument("-f", "--source_file", dest="SourceFile", type=str, help="specify the XNAT source file")
parser.add_argument("-x", "--strip_file", dest="SourceStrip", type=str, default="false", help="strip session info from results files")

parser.add_argument("-D", "--destination_dir", dest="DestDir", default='tmp', type=str, help="specify the directory for output")
parser.add_argument("-M", "--print_csv", dest="printLists", default=False, help="print the lists to a csv file for looking at")
parser.add_argument("-V", "--Verbose", dest="Verbose", type=str, default='false', help="show more verbose output")

parser.add_argument('--version', action='version', version='%(prog)s 0.1')

args = parser.parse_args()

User = args.User
Password = args.Password
Project = args.Project
Subject = args.Subject
Server = args.Server
DestDir = args.DestDir
DataType = args.DataType
SourceDir = args.SourceDir
SourceFile = args.SourceFile
SourceStrip = args.SourceStrip
Verbose = str.lower(args.Verbose)

DestDir = os.path.normpath(args.DestDir)+os.sep

#===============================================================================
# init interface to server and get subjects if none input...
#===============================================================================
print "Running %s on %s" % (os.path.split(sys.argv[0])[1], socket.gethostname())
getHCP = getHCP(User, Password, Server)
#===============================================================================

#===============================================================================
# SourceFile = 'ALL'
# SourceDir = 'MNINonLinear'
# SourceFile = 'ALL'
# SourceDir = 'MNINonLinear/Native'
# SourceFile = 'ALL'
# SourceDir = 'MNINonLinear/fsaverage_LR32k'
# SourceFile = 'ALL'
# SourceDir = 'T1w/Native'
# SourceFile = 'ribbon.nii.gz'
#
# SourceDir = 'T1w' 
T1wList = list()
T1wList.append('T1w_acpc_dc.nii.gz')
T1wList.append('T1w_acpc_dc_restore.nii.gz')
T1wList.append('T1w_acpc_dc_restore_brain.nii.gz')
T1wList.append('T2w_acpc_dc.nii.gz')
T1wList.append('T2w_acpc_dc_restore.nii.gz')
T1wList.append('T2w_acpc_dc_restore_brain.nii.gz')
T1wList.append('T1wDividedByT2w.nii.gz')
T1wList.append('T1wDividedByT2w_ribbon.nii.gz')
T1wList.append('brainmask_fs.nii.gz')
T1wList.append('wmparc.nii.gz')
T1wList.append('BiasField_acpc_dc.nii.gz')
# 
# SourceDir = '/MNINonLinear/xfms'
xfmsList = list()
xfmsList.append('acpc_dc2standard.nii.gz')
xfmsList.append('NonlinearRegJacobians.nii.gz')
xfmsList.append('standard2acpc_dc.nii.gz')
#===============================================================================
getHCP.DestinationDir = DestDir
getHCP.Project = Project
getHCP.Subject = Subject
SubjectSessions, SubjectSessionsType = getHCP.getSubjectSessions()
SubjectResources = getHCP.getSubjectResources()
SubjectSessionsUniq = set(SubjectResources.get('FileSessions'))

print 'Length of unique File Names, URIs, Sessions, and Labels: %s %s %s %s' % ( len(set(SubjectResources.get('FileNames'))), len(set(SubjectResources.get('FileURIs'))), len(set(SubjectResources.get('FileSessions'))), len(set(SubjectResources.get('FileLabels'))) )

Keywords = ['MNINonLinear', 'MNINonLinear/Native', 'MNINonLinear/fsaverage_LR32k', 'T1w/Native']
FileURIs = SubjectResources.get('FileURIs')
FileSessions = SubjectResources.get('FileSessions')
FileNames = SubjectResources.get('FileNames')
#matching = [s for s in FileURIs if Keyword in s]
#print all(x in Keywords for x in FileURIs)

if (DataType == 'Structural') and (SubjectSessionsType.count('strc') > 1):
    for i in xrange(0, len(SubjectSessions)):
        if (SubjectSessions[i].find('xtr') != -1):
            relevantSession = SubjectSessions[i]
            
PrintList = list()
#===============================================================================
# Do ALL for a given set of directories...
#===============================================================================
for i in xrange(0, len(FileURIs)):
    for j in xrange(0, len(Keywords)):
        if (FileURIs[i].find(Keywords[j] +'/'+ FileNames[i].replace('"', '')) != -1) and (FileSessions[i] == relevantSession):
            getHCP.DestinationDir = DestDir + Keywords[j] 
            getHCP.writeFileFromURL([FileURIs[i]])
            PrintList.append(FileURIs[i])

#currDestDir = os.path.normpath(destDir +os.sep+ SourceDir)
#getHCP.fWriteFileFromURL(PrintList)

PrintList = list()
#===============================================================================
# Do T1w file specific matching...
#===============================================================================
SourceDir = 'T1w/'
for i in xrange(0, len(FileURIs)):
    for j in xrange(0, len(T1wList)):
        if (FileURIs[i].find(SourceDir + T1wList[j]) != -1) and (FileSessions[i] == relevantSession):
            getHCP.DestinationDir = DestDir + SourceDir 
            getHCP.writeFileFromURL([FileURIs[i]])
            PrintList.append(FileURIs[i])
            
PrintList = list()
#===============================================================================
# Do xfms file specific matching....
#===============================================================================
SourceDir = '/MNINonLinear/xfms/'
for i in xrange(0, len(FileURIs)):
    for j in xrange(0, len(xfmsList)):
        if (FileURIs[i].find(SourceDir + xfmsList[j]) != -1) and (FileSessions[i] == relevantSession):
            getHCP.DestinationDir = DestDir + SourceDir
            getHCP.writeFileFromURL([FileURIs[i]])
            PrintList.append(FileURIs[i])
            
print("Duration: %s" % (time.time() - sTime))
