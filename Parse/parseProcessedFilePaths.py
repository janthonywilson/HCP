'''
Created on Jan 25, 2013

@author: Tony
'''
import base64
import sys
import xml
import os
import csv
import time
import filecmp
import urllib2
import difflib
import argparse
import xml.etree.ElementTree as ET

sTime = time.time()

#===============================================================================
# PARSE INPUT
#===============================================================================
parser = argparse.ArgumentParser(description="Script to parse provenance XNAT file ...")

parser.add_argument("-User", "--User", dest="User", default='tony', type=str)
parser.add_argument("-Password", "--Password", dest="Password", default='none', type=str)
parser.add_argument("-Subject", "--Subject", dest="Subject", default='none', type=str)
parser.add_argument("-InputDir", "--InputDir", dest="inputDir", default='none', type=str)
parser.add_argument("-InputDirFile", "--InputDirFile", dest="inputDirFile", default='none', type=str)
parser.add_argument("-OutputDir", "--OutputDir", dest="outputDir", default='none', type=str)
parser.add_argument("-OutputDirFile", "--OutputDirFile", dest="outputDirFile", default='none', type=str)

parser.add_argument('--version', action='version', version='%(prog)s 0.1')

args = parser.parse_args()

User = args.User
Password = args.Password
Subject = args.Subject
inputDir = os.path.normpath(args.inputDir) + os.sep
inputDirFile = args.inputDirFile

outputDir = os.path.normpath(args.outputDir) + os.sep
outputDirFile = args.outputDirFile
#===============================================================================
if not os.path.exists(outputDir):
    os.makedirs(outputDir)

headerStr = ['FilePathStruct', 'FilePathFunct', 'Val']
fileDiffObj = csv.writer(open(outputDirFile, 'wb'), delimiter='\t')
fileDiffObj.writerow(headerStr)
    
#for j in xrange(0, StatsListSz[0]):
#    # NOTE: there is also "writerows" to get rid of the for loop...
#    fileStatsId.writerow(StatsList[j,:])
#===============================================================================
refFileObj = open(inputDir + Subject + '_strc_list.txt', 'r')
refFileData = refFileObj.read().split('\n')
refFileObj.close()

compareFileObj = open(inputDir + Subject + '_fnca_list.txt', 'r')
compareFileData = compareFileObj.read().split('\n')
compareFileObj.close()

compareFilesExt = list()
compareFiles = list()
compareDirs = list()
for i in xrange(0, len(compareFileData)):
    compareDir, compareFileName = os.path.split(compareFileData[i])
    compareFileBaseName, compareFileExt = os.path.splitext(compareFileName)
    if ( len(compareFileName) != 0 ):
        compareFilesExt.append(compareFileExt)
        compareFiles.append(compareFileName)
        compareDirs.append(compareDir)
    
compareFilesUniq = set(compareFiles)
compareFilesExtUniq = set(compareFilesExt)
    
for i in xrange(0, len(refFileData)):
    currRefFileName = refFileData[i]
    if (currRefFileName[-1] != '/') and (currRefFileName.find('.') != -1):
        refDir, refFileName = os.path.split(currRefFileName)
        refFileBaseName, refFileExt = os.path.splitext(refFileName)
        if (len(refDir) > 0):
            refDirSplit = refDir.split('/')
#            refIdx = refDir.index('RESOURCES/Details')
#            refDirTail = refDir[refIdx + len('RESOURCES/Details') :]
    
            #print compareFilesUniq
            if (refFileName in compareFilesUniq):
                #print refFileName      
                compareFileIdx = compareFiles.index(refFileName)
                compareDir = compareDirs[compareFileIdx]
                compareDirSplit = compareDirs[compareFileIdx].split('/')
                compareFileName = compareFiles[compareFileIdx]
                countFilesep = compareDir.count('/')
                
                refResourcesIdx = refDirSplit.index('RESOURCES')
                compareResourcesIdx = compareDirSplit.index('RESOURCES')
                    
                    
                #print  os.sep.join(refDirSplit[refResourcesIdx+2:-1]), os.sep.join(compareDirSplit[compareResourcesIdx+2:-1])
                if ( os.sep.join(refDirSplit[refResourcesIdx+2:-1]) +os.sep+ refFileName ) == ( os.sep.join(compareDirSplit[compareResourcesIdx+2:-1]) +os.sep+ compareFileName ):
#                    if (refFileName.find('brainmask_fs.nii.gz') != -1) or (refFileName.find('orig.mgz') != -1):
#                        print  ( os.sep.join(refDirSplit[refResourcesIdx+2:-1]) +os.sep+ refFileName ), ( os.sep.join(compareDirSplit[compareResourcesIdx+2:-1]) +os.sep+ compareFileName )
                
                    shellCompareCmd = 'diff -s --brief %s %s > /dev/null' % (refDir +os.sep+ refFileName, compareDir +os.sep+ compareFileName)
                    
                    if sys.platform == 'win32':
    #                        print shellCompareCmd
                        try:
                            cmp = filecmp.cmp(refDir +os.sep+ refFileName, compareDir +os.sep+ compareFileName, shallow=False)
                        except:
                            pass
                    else:
    #                        print shellCompareCmd
                        try:
                            sysDiff = os.system(shellCompareCmd)
#                            print sysDiff
                            
                            if sysDiff.find('identical'):
                                compareValue = True
                            elif sysDiff.find('differ'):
                                compareValue = False
                            else:
                                compareValue = 'Unknown'
                            
                            print [refDir +os.sep+ refFileName, compareDir +os.sep+ compareFileName, compareValue]
                            
                        except:
                            pass
                        
                        try: 
                            cmp = filecmp.cmp(refDir +os.sep+ refFileName, compareDir +os.sep+ compareFileName, shallow=False)
#                            print cmp, refDir +os.sep+ refFileName, compareDir +os.sep+ compareFileName
                            fileDiffObj.writerow( [refDir +os.sep+ refFileName, compareDir +os.sep+ compareFileName, cmp] )
                            
                        except:
                            pass
    
                        
    #                        os.system(shellCompareCmd)
        
    #        print i, compareFileIdx, compareDir, refFileName, refDirName






