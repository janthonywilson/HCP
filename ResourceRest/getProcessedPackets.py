'''
Created on Aug 25, 2012

@author: Tony
'''


import os
import time
#import base64
import hashlib
import argparse
#import ResourceRest

sTime = time.time()

#===============================================================================
# PARSE INPUT
#===============================================================================
parser = argparse.ArgumentParser(description="Alpha program to pull NIFTI data from XNAT and put it somewhere...")

parser.add_argument("-W", "--server", dest="restRoot", default="https://intradb.humanconnectome.org", type=str, help="specify which server to connect to")
parser.add_argument("-u", "--username", dest="restUser", type=str, help="username must be specified")
parser.add_argument("-p", "--password", dest="restPass", type=str, help="password must be specified")
parser.add_argument("-P", "--project", dest="inputProject", default="HCP_Phase2", type=str, help="specify project")
parser.add_argument("-S", "--subject", dest="inputSubject", default="792564", type=str, help="specify subject of interest")
parser.add_argument("-K", "--package", dest="inputPacket", default="Structural", type=str, help="specify packet of interest")
parser.add_argument("-fn", "--fMRIname", dest="inputfMRIName", default="BOLD_EMOTION1_RL", type=str, help="specify name of fMRI packet of interest")
parser.add_argument("-D", "--destination_dir", dest="destDir", default='tmp', type=str, help="specify the directory for output")
parser.add_argument("-x", "--strip_file", dest="SourceStrip", type=str, default="true", help="strip session info from results files")
parser.add_argument("-t", "--time_out", dest="Timeout", type=float, default=15.0, help="change timeout")
parser.add_argument("-V", "--Verbose", dest="Verbose", type=str, default="false", help="show more verbose output")

parser.add_argument("-FS", "--structualFS", dest="structuralFS", default="false", type=str, help="specify if FreeSurfer data is desired")
parser.add_argument('--version', action='version', version='%(prog)s 0.3.1')

args = parser.parse_args()

inputProject = args.inputProject
inputSubject = args.inputSubject
inputPacket = args.inputPacket
SourceStrip = args.SourceStrip.lower()
structuralFS = args.structuralFS.lower()
inputfMRIName = args.inputfMRIName
timeoutURL = str(args.Timeout)
Verbose = args.Verbose.lower()

restUser = args.restUser
restPass = args.restPass
restRoot = args.restRoot
destDir = args.destDir

hashUser = hashlib.sha1(restUser).hexdigest()
hashPass = hashlib.sha1(restPass).hexdigest()

#hashUser = hashlib.md5(restUser).hexdigest()
#hashPass = hashlib.md5(restPass).hexdigest()

# hack for using command line string input...
if (Verbose == 'false'): VerboseBool = False
elif (Verbose == 'true'): VerboseBool = True


#==============================================================================
# MNINonLinear/*
# T1w/Native/*
# T1w/T1w_acpc_dc_restore.nii.gz
# T1w/T2w_acpc_dc_restore.nii.gz
# T1w/T1wDividedByT2w.nii.gz
# T1w/T1wDividedByT2w_ribbon.nii.gz
# T1w/brainmask_fs.nii.gz
# T1w/wmparc.nii.gz
# T1w/ribbon.nii.gz
# T1w/BiasField_acpc_dc.nii.gz
# 
# Some people might also want the FreeSurfer outputs:
# T1w/$SubjectID/*
#
# Changes on 29/09/2012
#DELETE: ${Subject}/T1w/lh.ribbon.nii.gz (wrong file)
#ADD: ${Subject}/T1w/ribbon.nii.gz (correct file)
#ADD: ${Subject}/T1w/T1w_acpc_dc.nii.gz
#ADD: ${Subject}/T1w/T1w_acpc_dc_restore_brain.nii.gz
#ADD: ${Subject}/T1w/T1w_acpc_dc_restore_brain.nii.gz
#ADD: ${Subject}/T1w/T1w_acpc_dc_restore_brain.nii.gz
#ADD: ${Subject}/MNINonLinear/Native/*
#ADD: ${Subject}/MNINonLinear/fsaverage_LR32k/*
#ADD: ${Subject}/MNINonLinear/xfms/acpc_dc2standard.nii.gz
#ADD: ${Subject}/MNINonLinear/xfms/NonlinearRegJacobians.nii.gz
#ADD: ${Subject}/MNINonLinear/xfms/standard2acpc_dc.nii.gz

#==============================================================================
if (inputPacket == 'Structural'):
    
    SourceFile = 'ALL'
    SourceDir = 'MNINonLinear'
    currDestDir = os.path.normpath(destDir +os.sep+ SourceDir)
    structCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -f ' +SourceFile+ ' -D ' +currDestDir+ ' -x ' +SourceStrip+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -f ' +SourceFile+ ' -D ' +currDestDir+ ' -x ' +SourceStrip+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    os.system(structCall)
    
    SourceFile = 'ALL'
    SourceDir = 'MNINonLinear/Native'
    currDestDir = os.path.normpath(destDir +os.sep+ SourceDir)
    structCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -f ' +SourceFile+ ' -D ' +currDestDir+ ' -x ' +SourceStrip+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -f ' +SourceFile+ ' -D ' +currDestDir+ ' -x ' +SourceStrip+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    os.system(structCall)
    
    SourceFile = 'ALL'
    SourceDir = 'MNINonLinear/fsaverage_LR32k'
    currDestDir = os.path.normpath(destDir +os.sep+ SourceDir)
    structCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -f ' +SourceFile+ ' -D ' +currDestDir+ ' -x ' +SourceStrip+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -f ' +SourceFile+ ' -D ' +currDestDir+ ' -x ' +SourceStrip+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    os.system(structCall)
    
    SourceFile = 'ALL'
    SourceDir = 'T1w/Native'
    currDestDir = os.path.normpath(destDir +os.sep+ SourceDir)
    structCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -f ' +SourceFile+ ' -D ' +currDestDir+ ' -x ' +SourceStrip+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -f ' +SourceFile+ ' -D ' +currDestDir+ ' -x ' +SourceStrip+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    os.system(structCall)
    
    SourceFile = 'ribbon.nii.gz'
    SourceDir = 'T1w'
    currDestDir = os.path.normpath(destDir +os.sep+ SourceDir)
    structCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -f ' +SourceFile+ ' -D ' +currDestDir+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -f ' +SourceFile+ ' -D ' +currDestDir+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    os.system(structCall)
        
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
#    T1wList.append('ribbon.nii.gz')
    T1wList.append('BiasField_acpc_dc.nii.gz')
    
    SourceDir = 'T1w'
    currDestDir = os.path.normpath(destDir +os.sep+ SourceDir)
    for i in xrange(0, len(T1wList)):
        SourceFile = T1wList[i]
        
        structCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -f ' +SourceFile+ ' -D ' +currDestDir+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
        if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -f ' +SourceFile+ ' -D ' +currDestDir+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
        os.system(structCall)
            
    try:
        badRibbonFileId = open(currDestDir + '/T1w/lh.ribbon.nii.gz', 'r')
        if (badRibbonFileId != -1):
            badRibbonFileId.close()
            os.remove(currDestDir + '/T1w/lh.ribbon.nii.gz')
    except:
        print 'bad ribbon file is not present...'
        
    xfmList = list()
    xfmList.append('acpc_dc2standard.nii.gz')
    xfmList.append('NonlinearRegJacobians.nii.gz')
    xfmList.append('standard2acpc_dc.nii.gz')

    SourceDir = '/MNINonLinear/xfms'
    currDestDir = os.path.normpath(destDir +os.sep+ SourceDir)
    for i in xrange(0, len(xfmList)):
        SourceFile = xfmList[i]
        
        structCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -f ' +SourceFile+ ' -D ' +currDestDir+ ' -x ' +SourceStrip+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
        if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -f ' +SourceFile+ ' -D ' +currDestDir+ ' -x ' +SourceStrip+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
        os.system(structCall)
                
    if (structuralFS == 'true'):
        print "Warning: No Structual FS will be downloaded.  This feature is not currently supported..."

#        SourceDir = 'T1w/' +inputSubject
#        SourceFile = 'ALL'
#        currDestDir = os.path.normpath(destDir +os.sep+ SourceDir)
#        structCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -f ' +SourceFile+ ' -D ' +currDestDir+ ' -x ' +SourceStrip
#        print structCall
#        os.system(structCall)

    
    
#===============================================================================
# MNINonLinear/Results/$fMRIName/$fMRIName .nii.gz
# MNINonLinear/Results/$fMRIName/$fMRIName _Atlas.dtseries.nii
# MNINonLinear/Results/$fMRIName/Movement_Regressors.txt
# MNINonLinear/Results/$fMRIName/Movement_Regressors_dt.txt
# 
# Other files currently needed for task fMRI analysis (tfMRI_{Task}): 
# MNINonLinear/Results/$fMRIName/$fMRIName _s2.atlasroi.L.32k_fs_LR.func.gii
# MNINonLinear/Results/$fMRIName/$fMRIName _s2.atlasroi.R.32k_fs_LR.func.gii
# MNINonLinear/Results/$fMRIName/$fMRIName _AtlasSubcortical_s2.nii.gz
# 
# Changes on 29/09/2012...
# Functional Packet:
# Lets always include (currently your script has an option inside of it that turns this on, but I've decided it shouldn't be optional):
# ${Subject}/MNINonLinear/Results/${fMRIName}/${fMRIName}_AtlasSubcortical_s2.nii.gz
# ${Subject}/MNINonLinear/Results/${fMRIName}/${fMRIName}_s2.atlasroi.L.32k_fs_LR.func.gii
# ${Subject}/MNINonLinear/Results/${fMRIName}/${fMRIName}_s2.atlasroi.R.32k_fs_LR.func.gii
# 
# ADD: ${Subject}/MNINonLinear/Results/${fMRIName}/${fMRIName}_SBRef.nii.gz
# ADD: ${Subject}/MNINonLinear/Results/${fMRIName}/RibbonVolumeToSurfaceMapping/goodvoxels.nii.gz
#===============================================================================

elif (inputPacket == 'Functional'):
    
    SourceDir = 'MNINonLinear/Results/' + inputfMRIName
    currDestDir = os.path.normpath(destDir +os.sep+ SourceDir)
    
    currFileName = inputfMRIName+ '.nii.gz'
    functCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    os.system(functCall)
    
    currFileName = inputfMRIName+ '_SBRef.nii.gz' 
    functCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    os.system(functCall)

    currFileName = inputfMRIName+ '_Atlas.dtseries.nii' 
    functCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    os.system(functCall)
    
    currFileName = 'Movement_Regressors.txt'
    functCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    os.system(functCall)
    
    currFileName = 'Movement_Regressors_dt.txt'
    functCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    os.system(functCall)
    
#    ${SubjectID}/MNINonLinear/Results/${fMRIName}/${fMRIName}_Jacobian.nii.gz
    currFileName = inputfMRIName+ '_Jacobian.nii.gz' 
    functCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    os.system(functCall)
    
#    if TaskData:
    currFileName = inputfMRIName+ '_s2.atlasroi.L.32k_fs_LR.func.gii'
    functCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    os.system(functCall)
    
    currFileName = inputfMRIName+ '_s2.atlasroi.R.32k_fs_LR.func.gii'
    functCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    os.system(functCall)

    currFileName = inputfMRIName+ '_AtlasSubcortical_s2.nii.gz'
    functCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    os.system(functCall)
    
    SourceDir = 'MNINonLinear/Results/' + inputfMRIName + '/RibbonVolumeToSurfaceMapping'
    currDestDir = os.path.normpath(destDir +os.sep+ SourceDir)
    currFileName = 'goodvoxels.nii.gz'
    functCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    os.system(functCall)

    
elif (inputPacket == 'Diffusion'):
    #===============================================================================
    # ${SubjectID}/Diffusion/data/bvals
    # ${SubjectID}/Diffusion/data/bvecs
    # ${SubjectID}/Diffusion/data/data.nii.gz
    # ${SubjectID}/Diffusion/data/nodif_brain_mask.nii.gz
    # ${SubjectID}/Diffusion/data/grad_dev.nii.gz
    # ${SubjectID}/T1w/xfms/diff2str.mat
    # ${SubjectID}/T1w/xfms/str2diff.mat
    # ${SubjectID}/MNINonLinear/xfms/diff2standard.nii.gz
    # ${SubjectID}/MNINonLinear/xfms/standard2diff.nii.gz
    #===============================================================================

    SourceDir = 'Diffusion/data'
    currDestDir = os.path.normpath(destDir +os.sep+ SourceDir)
    
    currFileName = 'bvals'
    functCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    os.system(functCall)
    
    currFileName = 'bvecs'
    functCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    os.system(functCall)
    
    currFileName = 'data.nii.gz'
    functCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    os.system(functCall)
    
    currFileName = 'nodif_brain_mask.nii.gz'
    functCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    os.system(functCall)
    
    currFileName = 'grad_dev.nii.gz'
    functCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    os.system(functCall)
    
    SourceDir = 'T1w/xfms'
    currDestDir = os.path.normpath(destDir +os.sep+ SourceDir)
    
    currFileName = 'diff2str.mat'
    functCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    os.system(functCall)
    
    currFileName = 'str2diff.mat'
    functCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    os.system(functCall)
    
    SourceDir = 'MNINonLinear/xfms'
    currDestDir = os.path.normpath(destDir +os.sep+ SourceDir)
    
    currFileName = 'diff2standard.nii.gz'
    functCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    os.system(functCall)
    
    currFileName = 'standard2diff.nii.gz'
    functCall = 'python ResourceRest.py -u ' +restUser+ ' -p ' +restPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    if VerboseBool: print 'ResourceRest.py -u ' +hashUser+ ' -p ' +hashPass+ ' -W ' +restRoot+ ' -P ' +inputProject+ ' -S ' +inputSubject+ ' -d ' +SourceDir+ ' -D ' +currDestDir+ ' -f ' +currFileName+ ' -T ' +inputPacket+ ' -t ' +timeoutURL+ ' -V ' +Verbose
    os.system(functCall)
    

print("Processed Packets Duration: %s" % (time.time() - sTime))

