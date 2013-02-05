import os
import csv
import time
import numpy
import scipy
import socket
import argparse
import nipy as nip
import nibabel as nib
from scipy import stats
import matplotlib.pyplot as pyplot

#===============================================================================
# -D R:\nifti\QC_Tests  -N  CP10086_v3_T2w2_0_7mm.nii.gz -O C:\tmp\bar
#===============================================================================

#===============================================================================
# PARSE INPUT
#===============================================================================
parser = argparse.ArgumentParser(description="Test program to do Fourier slope statistics on a nifit image...")
parser.add_argument("-D", "--ImageDir", dest="niftiDir", default=os.getcwd(), type=str, help="specify nifti directory...")
parser.add_argument("-N", "--Image", dest="niftiFile", type=str, help="specify nifti file...")
parser.add_argument("-O", "--OutputDir", dest="outputDir", type=str, help="specify where to write the output files...")
parser.add_argument("-F", "--FractionElim", dest="fracElimSlices", type=float, default=0.19, help="specify what fraction of top and bottom slices to eliminate...")
parser.add_argument("-U", "--UseFractionElim", dest="fracElimSlicesUse", type=str, default="y", help="specify if fraction of slices to be eliminated...")
InputArgs = parser.parse_args()
niftiDir = InputArgs.niftiDir
niftiFile = InputArgs.niftiFile
outputDir = InputArgs.outputDir
fracElimSlices = InputArgs.fracElimSlices
fracElimSlicesUse = InputArgs.fracElimSlicesUse
#===============================================================================

sTime = time.time()
print "Running on " + socket.gethostname()
freqCutoffLow = 2;
freqCutoffHigh = 60;
plotSlope = False
plotResults = False
printData = True

#===============================================================================
# FUNCTIONS
#===============================================================================
def fFourierRegressCoeffs( Slice, freqCutoffLow, freqCutoffHigh, plotSlope ):
    fftImg = numpy.fft.fft2(Slice)
#    fftImgVec = numpy.abs(numpy.real(fftImg.ravel()))

    SliceSz = Slice.shape
    xySum = (xGrid + yGrid)
    radFreqs = numpy.sqrt( xySum )
    freqIdx = numpy.nonzero( (radFreqs > freqCutoffLow) & (radFreqs < freqCutoffHigh) )

    plotFreqs = numpy.log10(radFreqs[freqIdx[0], freqIdx[1]])
    imageFFT = numpy.fft.fftshift(fftImg)
    imagePow = abs(imageFFT)**2
    imageZeroIdx = numpy.where(imagePow <= 0)
    if len(imageZeroIdx[0]) > 0:
        imagePow += (numpy.random.rand(SliceSz[0], SliceSz[1]) * 1e-16)


    plotPows = numpy.log10(imagePow[freqIdx[0], freqIdx[1]])
    
    #slope, intercept, r-value, p-value, stderr
    regressCoeffs = stats.linregress( plotFreqs, plotPows )
    
    if plotSlope:
        xPlot = numpy.arange(numpy.min(plotFreqs), numpy.max(plotFreqs), 0.1)
        pyplot.figure()
        pyplot.scatter( plotFreqs, plotPows, hold=True )
        pyplot.plot( xPlot, regressCoeffs[0] * xPlot + regressCoeffs[1], color='k', linewidth=3 )
        pyplot.show()
            
    return regressCoeffs
#===============================================================================
def fPrintData( inputMat, inputSubject, outputDir ):
    
    StatsListLen = len(inputMat) 
    headerStr = ['Volume', 'Slice', 'Slope', 'Intercept', 'r-value', 'p-value', 'stderr']
    fileName = inputSubject + '_FourierSlopeStatistics_' +str(len(numpyNipyDataSz))+ '.txt'

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    fileStatsId = csv.writer(open(outputDir +os.sep+ fileName, 'wb'), delimiter='\t')
    fileStatsId.writerow(headerStr)
    
    for i in xrange(0, StatsListLen):
            fileStatsId.writerow(inputMat[i,:])
#===============================================================================
def fPrintMeanStdData( inputMat, inputSubject, outputDir ):
    
    StatsListLen = len(inputMat) 
    headerStr = ['Volume', 'Mean','Std']
    fileName = inputSubject + '_FourierSlopeStatistics_' +str(len(numpyNipyDataSz))+ '_MeanStd.txt'

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    fileStatsId = csv.writer(open(outputDir +os.sep+ fileName, 'wb'), delimiter='\t')
    fileStatsId.writerow(headerStr)
    
    for i in xrange(0, StatsListLen):
            fileStatsId.writerow(inputMat[i,:])

#===============================================================================
def fStripSession( inputName ):
    # check for session on input subject string...
    if (inputName.find("_v") != -1):
        # strip out the session stuff.
        underscoreStartIdx = inputName.index("_")
        underscoreEndIdx = inputName[underscoreStartIdx+1:len(inputName)].index("_")
        outputName = inputName[0:underscoreStartIdx] + inputName[underscoreStartIdx+underscoreEndIdx+1:len(inputName)]
    else:
        outputName = inputName
        
    return outputName
#===============================================================================
def fStripExtension( inputName ):
    inputNameNoExtension, inputNameExtension = os.path.splitext(inputName)
    
    if inputNameExtension == '.gz':
        inputNameNoExtension, inputNameExtension = os.path.splitext(inputNameNoExtension)
        return inputNameNoExtension
    else:
        return inputNameNoExtension
#===============================================================================

inputFileName = niftiFile
niftiFile = niftiDir +os.sep+ niftiFile
print "File: " + niftiFile

numpyNipyData = nip.load_image(niftiFile)
numpyNipyData = numpy.float32(numpy.asarray(numpyNipyData))
numpyNipyDataSz = numpy.asarray(numpyNipyData.shape)
CenterPoint = numpyNipyDataSz / 2

ndMesh =  numpy.asarray(numpy.mgrid[0:numpyNipyDataSz[0],0:numpyNipyDataSz[1]])
xGrid = (ndMesh[0,:] - CenterPoint[0])**2
yGrid = (ndMesh[1,:] - CenterPoint[1])**2 


#set up bool for keeping middle slices...
nElimSlices = round(numpyNipyDataSz[2] * fracElimSlices)
boolElimSlices = numpy.zeros(numpyNipyDataSz[2], dtype=numpy.int)
if (fracElimSlicesUse == 'y'):
    boolElimSlices[nElimSlices:(numpyNipyDataSz[2] - nElimSlices)] = 1
else:
    boolElimSlices[0:numpyNipyDataSz[2]] = 1
ElimSlicesIdx = numpy.nonzero(boolElimSlices)
print "Eliminating " +str( round(numpyNipyDataSz[2] - sum(boolElimSlices)) )+ " of " +str(numpyNipyDataSz[2])+ " slices ..."


if len(numpyNipyDataSz) == 3:
    currVol = numpyNipyData    
    regressCoeffsResults = numpy.zeros([numpyNipyDataSz[2], 7])
    for i in xrange(0, numpyNipyDataSz[2]):
        currSlice = currVol[:,:,i]

        regressCoeffs = fFourierRegressCoeffs( currSlice, freqCutoffLow, freqCutoffHigh, plotSlope )
        regressCoeffsResults[i,:] = numpy.concatenate([ [1, (i+1)], regressCoeffs ])
        
    if plotResults:
        pyplot.hist(regressCoeffsResults[:,2], bins=100)
        pyplot.figure()
        pyplot.scatter(xrange(0, numpyNipyDataSz[2]), regressCoeffsResults[:,2])
#        pyplot.figure()
#        pyplot.imshow(regressCoeffsResults, interpolation='none')
        pyplot.show()

elif len(numpyNipyDataSz) == 4:
    linIdx = 0
    regressCoeffsResults = numpy.zeros([numpyNipyDataSz[2] * numpyNipyDataSz[3], 7])
    volumeMeansStd = numpy.zeros([numpyNipyDataSz[3], 3])
    for h in xrange(0, numpyNipyDataSz[3]):
        currVol = numpyNipyData[:,:,:,h]
        sliceSlopes = numpy.zeros(numpyNipyDataSz[2])
        for i in xrange(0, numpyNipyDataSz[2]):
            currSlice = currVol[:,:,i]
            
            regressCoeffs = fFourierRegressCoeffs( currSlice, freqCutoffLow, freqCutoffHigh, plotSlope )
            sliceSlopes[i] = regressCoeffs[0]
            regressCoeffsResults[linIdx,:] = numpy.concatenate([ [(h+1), (i+1)], regressCoeffs ])
            linIdx += 1
        volumeMeansStd[h,0] = h+1
        volumeMeansStd[h,1] = numpy.average(sliceSlopes) 
        volumeMeansStd[h,2] = numpy.std(sliceSlopes)
        
if printData:
    niftiFileName = fStripExtension( inputFileName )
    outputName = fStripSession( niftiFileName )
    fPrintData( regressCoeffsResults, outputName, outputDir )
    if len(numpyNipyDataSz) == 4:
       fPrintMeanStdData( volumeMeansStd, outputName, outputDir )
    

tTime = time.time() - sTime
print("Duration: %s" % tTime)


#pngImg = Image.open("C:\Users\Tony\Pictures\Lena.jpg")
#pngImgSz = numpy.asarray(pngImg.size);
#
#print pngImg.size, pngImg.format
#pngImg.show()
#pngArray = numpy.asarray(pngImg)

#fftConj = numpy.fft.fftshift(numpy.real(fftImg)+numpy.imag(fftImg))
#fftImgConj = numpy.log10( fftConj )

#print 'Max/Min: ', [numpy.max(fftImgVec), numpy.min(fftImgVec)]
##pyplot.hist(numpy.real(fftImg.ravel()), bins=10000)
#print numpy.real(fftImg.ravel())
#print scipy.arange(len(numpy.real(fftImg.ravel()))+1)
#pyplot.semilogy(numpy.linspace(0, len(fftImg.ravel()), len(fftImg.ravel())), numpy.fft.fftshift(numpy.real(fftImg.ravel())))
#pyplot.figure()
#pyplot.imshow( fftImgConj )
#pyplot.colorbar()
#
#pyplot.figure()
#pyplot.imshow(radFreqs)
#pyplot.colorbar()
#
#pyplot.figure()
#pyplot.plot( plotFreqs, plotPows )
#pyplot.colorbar()




#x = numpy.arange(0, 2*math.pi, 0.1, dtype='float32')
#y = numpy.cos(x)
#z = y * numpy.sin(y)
#nGauss = numpy.random.randn(len(y)) * 0.15

#xLen = scipy.alen(x)
#xCorr = scipy.corrcoef(x, y)
#yFFT = numpy.abs(scipy.fft(y + nGauss))
#fft.fft(y, x, axis)


#printStr = 'X: ' + x + ' Length: ' + xCorr.astype(string)
#print(printStr)

#pyplot.scatter(x, y + nGauss)
#pyplot.figure()
#pyplot.hist(nGauss, 100)
#pyplot.figure()
#pyplot.stem(x, yFFT)

#fig = plt.figure()
#ax = Axes3D(fig)
#ax.scatter(x,y,z)
#pyplot.draw()
#pyplot.show()
