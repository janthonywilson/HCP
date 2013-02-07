local_path = os.path.dirname(sys.argv[0])
os.chdir(local_path)

dataNames = {}
dataNames['imageFolder'] = local_path + '\\MyInputFolder\\'
dataNames['workFolder'] = local_path + '\\MyOutputFolder\\'
dataNames['imageFile'] = 'MyImage.img'

parameters = {}
# Firts deformation
parameters['def1_MaxIter'] = 1
parameters['def1_gamma'] = 0.65
parameters['def1_nu'] = 0.4 
parameters['def1_lambda'] = 0.15
parameters['def1_delta'] = 0.5  
parameters['def1_l'] = 15
parameters['def1_Df'] = 20.
        
# Second deformation
parameters['def2_MaxIter'] = 1
parameters['def2_gamma'] = 0.3
parameters['def2_nu'] = 0.4   
parameters['def2_lambda'] = 0.4
parameters['def2_delta'] = 0.5
parameters['def2_l'] = 8
parameters['def2_Df'] = 1.
parameters['def2_S'] = 2
parameters['def2_dmin'] = 4.
parameters['def2_dmax'] = 5.
parameters['def2_dmean'] = 2.
parameters['def2_D'] = 0.3
parameters['def2_dp'] = 0.5

# Third deformation
parameters['def3_MaxIter'] = 1 
parameters['def3_gamma'] = 0.3
parameters['def3_nu'] = 0.4    
parameters['def3_lambda'] = 0.4
parameters['def3_delta'] = 0.5
parameters['def3_l'] = 8
parameters['def3_Df'] = 1.
parameters['def3_S'] = 2
parameters['def3_dmin'] = 4.
parameters['def3_dmax'] = 5.
parameters['def3_dmean'] = 2.
parameters['def3_D'] = 0.3
parameters['def3_dp'] = 0.5 


from ToolsSMHASS import *
import scipy.ndimage.morphology as ndimageMorphology
import scipy.ndimage.filters as ndimageFilters
import scipy.ndimage.measurements as ndimageMeasurements
import scipy.ndimage.interpolation as ndimageInterpolation
from scipy.optimize.minpack import leastsq
from scipy.signal.bsplines import cubic
##from pylab import *
import copy
import sys
import timeit

# ######### Names #############3
dataNames['preSegmentedImageFile'] = dataNames['imageFile'][:-4] + '_PreSeg.img'
dataNames['mesh3deformedFilename'] = dataNames['imageFile'][:-4] + '_3def'
dataNames['dataPreSeg'] = dataNames['imageFile'][:-4] + '_DataPreSeg'
dataNames['cortexMeshModel'] = 'mallaModeloCortezaST_LONI'
dataNames['bainBinModel'] = 'LoniBinAtlas.img'
dataNames['cortexRegisterFileName'] =  dataNames['imageFile'][:-4] + '_Reg'
dataNames['cortexDef1'] = dataNames['imageFile'][:-4] + '_Def1'
dataNames['cortexDef2'] = dataNames['imageFile'][:-4] + '_Def2'
dataNames['cortexDef3'] = dataNames['imageFile'][:-4] + '_Def3'
dataNames['binSegImage'] = dataNames['imageFile'][:-4] + '_SegBin'
dataNames['SegImage'] = dataNames['imageFile'][:-4] + '_Seg'


##
##os.chdir('C:\\Chubo\\Doctorado\\Python\\Programas\\Para_publicar')

def SimplexDeform3DASM(simplexMesh, alfa, beta, gamma, kappa, forces, ITER, S, sigma, recomputeParameters, thresholdITER, Scontour, AngleGammaObjetivo, angleSimplexObjetivo, modo, extraData = None, iterCorrection = Inf):

    """    
    sigma:   parametro para la acomodacion de la malla
    S:       tamaño de la vecindad para calcular el angulo objetivo
    alpha:   parametro de fuerza interna
    gamma:   parametro de viscosidad
    kappa:   parametro de fuerza externa
    fx,fy:   campos de fuerza externa
    thresholdITER: umbral de desplazamiento promedio en la malla para detener las iteraciones
    """

    dibujos_tag = 0
    calculoZonaEncerrada = 1
    flag_desplegarDatos = 0
    if calculoZonaEncerrada:
        import copy

    def calculateFext1erDef(dSample, gmax, points, normals, imageData):

        dx = imageData['dx']
        dy = imageData['dy']
        dz = imageData['dz']
        imagen_bin = ['imagen_bin']
        
        dSampleIn = 40
        dSampleOut = parameters['def1_l'] #10
        
        steps = arange(-dSampleIn,dSampleOut+1)
        
        gradient_aux = zeros((points.shape[0], 3), 'float32')
        xTar = zeros((points.shape[0], 3), 'float32')
        samplesImage = zeros((points.shape[0], len(steps)), 'float32')
        samplesImage_Borrar = zeros((simplex.shape[0], len(steps)), 'float32') # BORRAR

        
        for i in range(len(steps)):
            gradient_aux[:,1] = map_coordinates(dx, (points + (normals * steps[i])).T, output_type = None, output = None, order = 1, mode = 'constant', cval = 0.0, prefilter = False)
            gradient_aux[:,0] = map_coordinates(dy, (points + (normals * steps[i])).T, output_type = None, output = None, order = 1, mode = 'constant', cval = 0.0, prefilter = False)
            gradient_aux[:,2] = map_coordinates(dz, (points + (normals * steps[i])).T, output_type = None, output = None, order = 1, mode = 'constant', cval = 0.0, prefilter = False)
            
            normGradient_aux = (gradient_aux * gradient_aux).sum(1)
            
##                G = (gradient_aux * normals).sum(1) * ((gmax * (gmax + normGradient_aux)) / (gmax * gmax + normGradient_aux * normGradient_aux))
            
            G = (gradient_aux * normals).sum(1)
            
            # criterios de eliminacion
            # imponer una regla para que solo quede el borde mas lejano (para que se segmente la envolvente)
##                F = G.copy()
            F = G.copy()
            F[ G < 0.] = 0.
            # ####
            
            samplesImage[:,i] =  F
            
        indice = range(steps.shape[0])
        indice.reverse()
        for i in range(samplesImage.shape[0]):
            if samplesImage[i,:].max() > 0.000000001:
                aux = samplesImage[i,:].max() * 0.5
                for j in indice:
                    if samplesImage[i,j] > aux:
                        index_aux = j - dSampleIn
                        break
            else:
                index_aux = 0                
            xTar[i,:] = points[i] + normals[i] * index_aux

##        # ############### para Visualizar
##        PtoVTK = PythonToPolyData()
##        PtoVTK.SetInputPoints(xTar)
##        PtoVTK.WritePolyData('..\\Borrar0.vtk')
##        # ###############
        
        gradient_aux[:,1] = map_coordinates(dx, xTar.T, output_type = None, output = None, order = 1, mode = 'constant', cval = 0.0, prefilter = False)
        gradient_aux[:,0] = map_coordinates(dy, xTar.T, output_type = None, output = None, order = 1, mode = 'constant', cval = 0.0, prefilter = False)
        gradient_aux[:,2] = map_coordinates(dz, xTar.T, output_type = None, output = None, order = 1, mode = 'constant', cval = 0.0, prefilter = False)
            
        normGradient_aux = sqrt((gradient_aux * gradient_aux).sum(1)).reshape(-1,1)
            
        Eext = ((gradient_aux / normGradient_aux) * (xTar - points)).sum(1).reshape(-1,1)
        
        Eext[normGradient_aux == 0.] = 0.
        
##            Fext[:,:] = xTar - points
        Fext[:,:] = normals * Eext



        desp = (xTar - points) # MODIFICADO 07-02-2011
        despNorm = sqrt((desp * desp).sum(1)).reshape(-1,1)
        # aplicacion de funcion de stiffness
        D = parameters['def1_Df'] # distancia sobre la cual la fueza decrece exponencialmente
        despNormD = despNorm - D
        resp = zeros(despNorm.shape, despNorm.dtype)
        cond1 = (despNormD<0.)
##        cond1 = less(despNormD, 1.)
        resp[cond1] = 1.
        resp[~cond1] = 1. / exp(despNormD[~cond1])

        Fext[:,:] = Fext[:,:] * resp



        return Fext
    
    def calculateFext2aDef(dSample, gmax, points, normals, imageData):
        dSample = 15
        dCalculo = 5
        D = 0.5
        

        imagen = imageData['imagen']
        dx = imageData['dx']
        dy = imageData['dy']
        dz = imageData['dz']
        
        gradient_aux = zeros((points.shape[0], 3), 'float32')
        xTar = zeros((points.shape[0], 3), 'float32')
        
        steps = arange(-(dSample+dCalculo),(dSample+dCalculo)+1)
        steps2 = arange(-dSample,dSample+1)
        samplesImage = zeros((points.shape[0], len(steps)), 'float32')
        for i in range(len(steps)):
            samplesImage[:,i] = map_coordinates(imagen, (points + (normals * steps[i])).T, output_type = None, output = None, order = 1, mode = 'constant', cval = 0.0, prefilter = False)
            
##            normGradient_aux = (gradient_aux * gradient_aux).sum(1)
##            
####                G = (gradient_aux * normals).sum(1) * ((gmax * (gmax + normGradient_aux)) / (gmax * gmax + normGradient_aux * normGradient_aux))
##            
##            G = (gradient_aux * normals).sum(1)
##            
##            # criterios de eliminacion
##            # imponer una regla para que solo quede el borde mas lejano (para que se segmente la envolvente)
####                F = G.copy()
##            F = G.copy()
##            F[ G < 0.] = 0.
##            # ####
##            
##            samplesImage[:,i] =  F


        
##        indice = range(1, dSample*2+1)
        indice = range(dCalculo, dCalculo+dSample*2+1)
        indice2 = range(dSample*2+1)
        data = zeros(31, 'float32')
        xTarIndex = zeros(points.shape[0], 'float32')
        G = zeros((samplesImage.shape[0], dSample*2+1), 'float32')
##        val_Borrar = zeros(samplesImage.shape, 'float32')
##        val_Borrar1 = zeros(samplesImage.shape, 'float32')
##        val_Borrar2 = zeros(samplesImage.shape, 'float32')
##        val_Borrar3 = zeros(samplesImage.shape, 'float32')
        for i in range(samplesImage.shape[0]):
            for j, j2 in zip(indice, indice2):
                
                std_aux = samplesImage[i,j-dCalculo:j].std()
                prom_aux = samplesImage[i,j-dCalculo:j].mean()
##                val_Borrar1[i,j] = ((prom_aux - samplesImage[i,j]) /prom_aux) / std_aux
                G[i,j2] = ((prom_aux - samplesImage[i,j:j+3].mean())) / std_aux
##                val_Borrar1[i,j] = ((prom_aux - samplesImage[i,j]))
##                val_Borrar2[i,j] = std_aux
##                val_Borrar3[i,j] = prom_aux
                
                
##                prom_aux = samplesImage[i,:j].mean()
##                val_Borrar[i,j] = ((samplesImage[i,j] - prom_aux) / prom_aux)
##                val_Borrar1[i,j] = samplesImage[i,:j].std()
##                val_Borrar2[i,j] = samplesImage[i,j:].std()
##                val_Borrar3[i,j] = val_Borrar1[i,j]- val_Borrar2[i,j]
                
                
            xTarIndex[i] = (G[i,:] - D * (steps2 * steps2)).argmax()

        if 1: #dibujo
            from pylab import imshow, figure, plot, colorbar, gray, hold, show
            rango = [9600,9660]
            figure()
            imshow(samplesImage[rango[0]:rango[1]], interpolation='nearest', origin='lower')
            colorbar()
            gray()
            hold('on')
            plot(xTarIndex[rango[0]:rango[1]] + dCalculo , range(rango[1]-rango[0]))
            show()
            
            figure()
            imshow(G[rango[0]:rango[1]], interpolation='nearest', origin='lower')
            colorbar()
            gray()
            show()
            
            
        xTar[:,:] = points + normals * (xTarIndex - dSample).reshape(-1,1)
            

                    

        # ############### para Visualizar
        PtoVTK = PythonToPolyData()
        PtoVTK.SetInputPoints(xTar)
        PtoVTK.WritePolyData('..\\Borrar0.vtk')
        # ###############
        
##        gradient_aux[:,1] = map_coordinates(dx, xTar.T, output_type = None, output = None, order = 1, mode = 'nearest', cval = 0.0, prefilter = False)
##        gradient_aux[:,0] = map_coordinates(dy, xTar.T, output_type = None, output = None, order = 1, mode = 'nearest', cval = 0.0, prefilter = False)
##        gradient_aux[:,2] = map_coordinates(dz, xTar.T, output_type = None, output = None, order = 1, mode = 'nearest', cval = 0.0, prefilter = False)
##            
##        normGradient_aux = sqrt((gradient_aux * gradient_aux).sum(1)).reshape(-1,1)
##            
##        Eext = ((gradient_aux / normGradient_aux) * (xTar - points)).sum(1).reshape(-1,1)
        
##        Eext[normGradient_aux == 0.] = 0.
        
        Fext[:,:] = xTar - points
##        Fext[:,:] = normals * Eext

        return Fext
    
    def calculateFext2aDef2(dSample, gmax, points, normals, neighbors, imageData, extraData):
        calculoGradiente = 'fuerzasDeEntrada' # 'fuerzasDeEntrada', 'MuestreoDeLineas'
        dSample = parameters['def2_l']#8 # 8
        sampleStep = parameters['def2_delta'] #0.5 # 1.
        dCalculo = 0
        porce = 0.8 # porcentaje (suponiendo que el maximo F es 1) que se restara al valor del gradiente al final de la linea de muestreo para encontrar el Xtarget
##        D = porce/(8.*8.) #0.5 8.
        D = (porce*22.)/(8.*8.) #0.5 8. # el 22 esta agregado porque ese es el valor maximo que puede tomar el gradiante al normalizar el mapa de bordes en [0,1] y aplicar el filtro sobel
        D = parameters['def2_D']
        gmax = (dSample*dSample * D) * 2.
        
        tipoMuestreoImagen = 'muestreoConCruces' # 'muestreoConCruces', 'muestreoEnLinea'
        
        imagen = imageData['imagen']
        dx = imageData['dx']
        dy = imageData['dy']
        dz = imageData['dz']
        
        gradient_aux = zeros((points.shape[0], 3), 'float32')
        xTar = zeros((points.shape[0], 3), 'float32')
        
        steps = arange(-dSample,dSample+sampleStep, sampleStep)
        steps2 = arange(-(dSample+dCalculo),(dSample+dCalculo)+ sampleStep, sampleStep)
        centerIndex = int(dSample/sampleStep)
        samplesImage = zeros((points.shape[0], len(steps2)), 'float32')
        F = zeros((points.shape[0], len(steps)), 'float32')

        if calculoGradiente == 'fuerzasDeEntrada':
            for i in range(len(steps2)):
                gradient_aux[:,1] = map_coordinates(dx, (points + (normals * steps2[i])).T, output_type = None, output = None, order = 2, mode = 'constant', cval = 0.0, prefilter = False)
                gradient_aux[:,0] = map_coordinates(dy, (points + (normals * steps2[i])).T, output_type = None, output = None, order = 2, mode = 'constant', cval = 0.0, prefilter = False)
                gradient_aux[:,2] = map_coordinates(dz, (points + (normals * steps2[i])).T, output_type = None, output = None, order = 2, mode = 'constant', cval = 0.0, prefilter = False)


    ##            normGradient_aux = (gradient_aux * gradient_aux).sum(1)
    ##
##                G = (gradient_aux * -normals).sum(1)
    ##            G *= -1.
    ##
    ##            # criterios de eliminacion
    ##            # imponer una regla para que solo quede el borde mas lejano (para que se segmente la envolvente)
    ####                F = G.copy()
##                F_aux = G.copy()
    ##            F_aux[ G < 0.] = 0.
    ##            # ####
                F[:,i] = (gradient_aux * -normals).sum(1)
                

        # ### muestreo de datos en imagen, valorex de voxeles en la imagen original ###
        for i in range(len(steps2)):
            samplesImage[:,i] = map_coordinates(imagen, (points + (normals * steps2[i])).T, output_type = None, output = None, order = 1, mode = 'constant', cval = 0.0, prefilter = False)
        if tipoMuestreoImagen == 'muestreoConCruces': # muestreo con cruces   # CORREGIDO 3-2-2011
            dc = 4. #largo de cruz mm
            xc = -dSample #posicion de cruz mm
            stepsC = arange(sampleStep,dc + sampleStep, sampleStep)
            samplesImageC = zeros((points.shape[0], 4, (dc/sampleStep)), 'float32')
            direcciones = zeros((points.shape[0],4,3), 'float32')
            direcciones[:,0,:] = points[neighbors[:,0]]- points
            direcciones[:,0,:] = direcciones[:,0,:] / sqrt((direcciones[:,0,:] * direcciones[:,0,:]).sum(1)).reshape(-1,1)
            direcciones[:,1,:] = -direcciones[:,0,:]
            direcciones[:,2,:] = cross(normals, direcciones[:,0,:])
            direcciones[:,3,:] = -direcciones[:,2,:]
            for direccion_index in range(4): # cada una de las 4 partes de la cruz
                for i in range(len(stepsC)):
                    samplesImageC[:,direccion_index,i] = map_coordinates(imagen, (points + (normals * xc) + (direcciones[:,direccion_index,:] * stepsC[i])).T, output_type = None, output = None, order = 1, mode = 'constant', cval = 0.0, prefilter = False)

##        # ######## filtro gaussiano 1D##########
##        radio = 3
##        s = sqrt(-((radio)**2)/(log(0.01)*2))
##        samplesImage_S = zeros(samplesImage.shape, 'float32')
##        ndimageFilters.gaussian_filter1d(samplesImage, sigma = s, axis = 1, order = 0, output = samplesImage_S, mode = 'constant', cval = 0.0)
##        samplesImage = samplesImage_S
##        del samplesImage_S

        if calculoGradiente == 'MuestreoDeLineas':
            # ############ 
            # hacia atras
            dex = array([-1,0,1], 'float32')
            dex = dex.reshape(1,3)        
            F = ndimageFilters.convolve(samplesImage, dex, output = None, mode = 'nearest', cval = 0.0, origin = 0)
        
##        F = contrast3D(F, -1., 1.)

##        F[((samplesImage < extraData['Totsu_grisOrig']*0.1) + (samplesImage > (extraData['u_nivel_grisOrig'][1] + 3.*extraData['s_nivel_grisOrig'][1]))).nonzero()] = 0
##        F[( (samplesImage > (extraData['u_nivel_grisOrig'][1] + 3.*extraData['s_nivel_grisOrig'][1]))).nonzero()] = 0
##        F[((samplesImage < extraData['Totsu_grisOrig']*0.1)).nonzero()] = 0        
        
##        MBx0 = abs(MBx0)
                   
##        indice = range(1, dSample*2+1)
        indice2 = range(int(dCalculo/sampleStep), int(dCalculo/sampleStep)+int(dSample/sampleStep)*2+1)
        indice = range(int(dSample/sampleStep)*2+1)
        data = zeros(31, 'float32')
        xTarIndex = ones(points.shape[0], 'int') * int(dSample/sampleStep)
        G = zeros((samplesImage.shape[0], len(steps)), 'float32')

##        for i in range(samplesImage.shape[0]):
####            if i == 9638:
####                print 'alto'
##            for j, j2 in zip(indice, indice2):
##
##                    if samplesImage[i,j2] < 200: # del liquido o del craneo
##                            if j > 2:
##                                xTarIndex[i] = F[i,j-2:j+3].argmax() + j - 1
##                            break
        aux = samplesImage.shape[1]
        aux_index = range(aux)
        listPointsForward = []
        listPointsBackward = []
        
        list1=[]
        list2=[]
        list3=[]
        list4=[]
        for i in range(samplesImage.shape[0]):
            
##            xTarIndex[i] = (F[i,:] - D * (steps2 * steps2)).argmax()
##            # ##### agregado 17-11-2010
####            if samplesImage[i][dSample-1:dSample+2].mean() > extraData['u_nivel_grisOrig'][0]:
####            if samplesImage[i][dSample-1:dSample+2].mean() > extraData['Totsu_grisOrig']:
##            mb = samplesImage[i][:dSample+1].max()
##            if samplesImage[i][dSample-4:dSample+1].mean() < (mb*0.66):
##                xTarIndex[i] = dSample + 6
##            # #########################

            Imini = samplesImage[i][int(centerIndex-(parameters['def2_dmin']/sampleStep)):centerIndex + 1 + (1./sampleStep)].min() # CORREGIDO 2-2-2011
            if tipoMuestreoImagen == 'muestreoEnLinea':
                mb = samplesImage[i][:centerIndex+1].max()  # muestreo en linea
            elif tipoMuestreoImagen == 'muestreoConCruces':
                mb = max(samplesImageC[i,:,:].max(), samplesImage[i][:centerIndex+1].max())  # muestreo con cruces   # CORREGIDO 3-2-2011

##            if mb < extraData['u_nivel_grisOrig'][2] - 2. * extraData['s_nivel_grisOrig'][2]:
##                mb = extraData['u_nivel_grisOrig'][2] - 2. * extraData['s_nivel_grisOrig'][2]

            if Imini > (mb*0.66): # se llego a una parte de intensidad alta
##                if samplesImage[i][:centerIndex + int(5./sampleStep) + 1].max() > (extraData['u_nivel_grisOrig'][2]+2.5*abs(extraData['s_nivel_grisOrig'][2])): # hay una intensidad muy alta adelante -> ojos
                if samplesImage[i][:centerIndex + int(parameters['def2_dmax']/sampleStep) + 1].max() > (mb*1.3) : # hay una intensidad muy alta adelante -> ojos
                    xTarIndex[i] = (F[i,:] - D * (steps2 * steps2)).argmax() #+ 1
                    list2.append(i)
                else:
##                    xTarIndex[i] = centerIndex + int(4/sampleStep)# 4
                    xTarIndex[i] = centerIndex + 1# 4
                    listPointsForward.append(i)
                    list3.append(i)
            elif Imini < extraData['u_nivel_grisOrig'][1]*0.1: # no hay duda que se llego al LCR
##                if samplesImage[i][int(centerIndex-(2./sampleStep)):centerIndex + 1].mean() < (mb*0.3): # esta en LCR asi que se hace retroceder
##                if samplesImage[i][int(centerIndex-(2./sampleStep)):centerIndex + 1].mean() < (extraData['u_nivel_grisOrig'][1] - 7*extraData['s_nivel_grisOrig'][1]): # esta en LCR asi que se hace retroceder
                if samplesImage[i][int(centerIndex-(parameters['def2_dmean']/sampleStep)):centerIndex + 1].mean() < (extraData['u_nivel_grisOrig'][1] - 8 * extraData['s_nivel_grisOrig'][1]): # esta en LCR asi que se hace retroceder
##                if samplesImage[i][int(centerIndex-(2./sampleStep)):centerIndex + 1].mean() < -1:
                    xTarIndex[i] = centerIndex - 1# 4
                    listPointsBackward.append(i)
                    list4.append(i)
                else:
                    xTarIndex[i] = (F[i,:] - D * (steps2 * steps2)).argmax() #+ 1
                    list1.append(i)
            else:
##                if samplesImage[i][int(centerIndex-(2./sampleStep)):centerIndex + 1].mean() < (mb*0.3): # esta en LCR asi que se hace retroceder
##                if samplesImage[i][int(centerIndex-(2./sampleStep)):centerIndex + 1].mean() < (extraData['u_nivel_grisOrig'][1] - 7*extraData['s_nivel_grisOrig'][1]): # esta en LCR asi que se hace retroceder
                if samplesImage[i][int(centerIndex-(2./sampleStep)):centerIndex + 1].mean() < (extraData['u_nivel_grisOrig'][1] - 8 * extraData['s_nivel_grisOrig'][1]): # esta en LCR asi que se hace retroceder
##                if samplesImage[i][int(centerIndex-(2./sampleStep)):centerIndex + 1].mean() < -1:
                    xTarIndex[i] = centerIndex - 1# 4
                    listPointsBackward.append(i)
                    list4.append(i)
                else:
                    xTarIndex[i] = (F[i,:] - D * (steps2 * steps2)).argmax() #+ 1
                    list1.append(i)                    
                
                
        if 0: #dibujo
            from pylab import imshow, figure, plot, colorbar, gray, hold, show
            rango = [6600,6660]
            figure()
            imshow(samplesImage[rango[0]:rango[1]], interpolation='nearest', origin='lower')
            colorbar()
            gray()
            hold('on')
            plot(xTarIndex[rango[0]:rango[1]] + dCalculo , range(rango[1]-rango[0]))
            show()
            
            figure()
            imshow(F[rango[0]:rango[1]], interpolation='nearest', origin='lower')
            colorbar()
            gray()  
            show()
            
            
        listPointsBackward = set(listPointsBackward)
        listPointsBackward = list(listPointsBackward)        
        xTar[:,:] = points + normals * (steps2[xTarIndex]).reshape(-1,1)
        xTar[listPointsForward,:] = points[listPointsForward,:] + normals[listPointsForward,:] * parameters['def2_dp'] #0.8
        xTar[listPointsBackward,:] = points[listPointsBackward,:] - normals[listPointsBackward,:] * parameters['def2_dp']
        
##        # ############### para Visualizar
##        PtoVTK = PythonToPolyData()
##        PtoVTK.SetInputPoints(xTar)
##        PtoVTK.WritePolyData('..\\xTar.vtk')
##        # ###############
##        
##        # ############### para Visualizar
##        PtoVTK = PythonToPolyData()
##        PtoVTK.SetInputPoints(points[list1])
##        PtoVTK.WritePolyData('..\\L1.vtk')
##        # ###############
##                # ############### para Visualizar
##        PtoVTK = PythonToPolyData()
##        PtoVTK.SetInputPoints(points[list2])
##        PtoVTK.WritePolyData('..\\L2.vtk')
##        # ###############
##        # ############### para Visualizar
##        PtoVTK = PythonToPolyData()
##        PtoVTK.SetInputPoints(points[list3])
##        PtoVTK.WritePolyData('..\\L3.vtk')
##        # ###############
##        # ############### para Visualizar
##        PtoVTK = PythonToPolyData()
##        PtoVTK.SetInputPoints(points[list4])
##        PtoVTK.WritePolyData('..\\L4.vtk')
##        # ###############

        
##        gradient_aux[:,1] = map_coordinates(dx, xTar.T, output_type = None, output = None, order = 1, mode = 'nearest', cval = 0.0, prefilter = False)
##        gradient_aux[:,0] = map_coordinates(dy, xTar.T, output_type = None, output = None, order = 1, mode = 'nearest', cval = 0.0, prefilter = False)
##        gradient_aux[:,2] = map_coordinates(dz, xTar.T, output_type = None, output = None, order = 1, mode = 'nearest', cval = 0.0, prefilter = False)
##            
##        normGradient_aux = sqrt((gradient_aux * gradient_aux).sum(1)).reshape(-1,1)
##            
##        Eext = ((gradient_aux / normGradient_aux) * (xTar - points)).sum(1).reshape(-1,1)
        
##        Eext[normGradient_aux == 0.] = 0.

        desp = (xTar - points)
        despNorm = sqrt((desp * desp).sum(1)).reshape(-1,1)

        # aplicacion de funcion de stiffness
        D = parameters['def2_Df'] # distancia sobre la cual la fueza decrece exponencialmente
        despNormD = despNorm - D
        resp = zeros(despNorm.shape, despNorm.dtype)
        cond1 = (despNormD<0.)
##        cond1 = less(despNormD, 1.)
        resp[cond1] = 1.
        resp[~cond1] = 1. / exp(despNormD[~cond1])

        Fext[:,:] = resp * desp

        
##        Fext[:,:] = xTar - points
##        Fext[:,:] = normals * Eext

        return Fext
    
    def calculateFext2aDef3(dSample, gmax, points, normals, neighbors, imageData, extraData):
        calculoGradiente = 'fuerzasDeEntrada' # 'fuerzasDeEntrada', 'MuestreoDeLineas'
        dSample = parameters['def3_l'] # 8
        sampleStep = parameters['def3_delta'] # 1.
        dCalculo = 0
        porce = 0.8 # porcentaje (suponiendo que el maximo F es 1) que se restara al valor del gradiente al final de la linea de muestreo para encontrar el Xtarget
##        D = porce/(8.*8.) #0.5 8.
##        D = (porce*22.)/(8.*8.) #0.5 8. # el 22 esta agregado porque ese es el valor maximo que puede tomar el gradiante al normalizar el mapa de bordes en [0,1] y aplicar el filtro sobel
        D = parameters['def3_D']
        gmax = (dSample*dSample * D) * 2.
        
        tipoMuestreoImagen = 'muestreoConCruces' # 'muestreoConCruces', 'muestreoEnLinea'
        
        imagen = imageData['imagen']
        dx = imageData['dx']
        dy = imageData['dy']
        dz = imageData['dz']
        
        gradient_aux = zeros((points.shape[0], 3), 'float32')
        xTar = zeros((points.shape[0], 3), 'float32')
        
        steps = arange(-dSample,dSample+sampleStep, sampleStep)
        steps2 = arange(-(dSample+dCalculo),(dSample+dCalculo)+ sampleStep, sampleStep)
        centerIndex = int(dSample/sampleStep)
        samplesImage = zeros((points.shape[0], len(steps2)), 'float32')
        F = zeros((points.shape[0], len(steps)), 'float32')

        if calculoGradiente == 'fuerzasDeEntrada':
            for i in range(len(steps2)):
                gradient_aux[:,1] = map_coordinates(dx, (points + (normals * steps2[i])).T, output_type = None, output = None, order = 2, mode = 'constant', cval = 0.0, prefilter = False)
                gradient_aux[:,0] = map_coordinates(dy, (points + (normals * steps2[i])).T, output_type = None, output = None, order = 2, mode = 'constant', cval = 0.0, prefilter = False)
                gradient_aux[:,2] = map_coordinates(dz, (points + (normals * steps2[i])).T, output_type = None, output = None, order = 2, mode = 'constant', cval = 0.0, prefilter = False)


    ##            normGradient_aux = (gradient_aux * gradient_aux).sum(1)
    ##
##                G = (gradient_aux * -normals).sum(1)
    ##            G *= -1.
    ##
    ##            # criterios de eliminacion
    ##            # imponer una regla para que solo quede el borde mas lejano (para que se segmente la envolvente)
    ####                F = G.copy()
##                F_aux = G.copy()
    ##            F_aux[ G < 0.] = 0.
    ##            # ####
                F[:,i] = (gradient_aux * -normals).sum(1)
                

        # ### muestreo de datos en imagen, valorex de voxeles en la imagen original ###
        for i in range(len(steps2)):
            samplesImage[:,i] = map_coordinates(imagen, (points + (normals * steps2[i])).T, output_type = None, output = None, order = 1, mode = 'constant', cval = 0.0, prefilter = False)
        if tipoMuestreoImagen == 'muestreoConCruces': # muestreo con cruces   # CORREGIDO 3-2-2011
            dc = 4. #largo de cruz mm
            xc = -dSample #posicion de cruz mm
            stepsC = arange(sampleStep,dc + sampleStep, sampleStep)
            samplesImageC = zeros((points.shape[0], 4, (dc/sampleStep)), 'float32')
            direcciones = zeros((points.shape[0],4,3), 'float32')
            direcciones[:,0,:] = points[neighbors[:,0]]- points
            direcciones[:,0,:] = direcciones[:,0,:] / sqrt((direcciones[:,0,:] * direcciones[:,0,:]).sum(1)).reshape(-1,1)
            direcciones[:,1,:] = -direcciones[:,0,:]
            direcciones[:,2,:] = cross(normals, direcciones[:,0,:])
            direcciones[:,3,:] = -direcciones[:,2,:]
            for direccion_index in range(4): # cada una de las 4 partes de la cruz
                for i in range(len(stepsC)):
                    samplesImageC[:,direccion_index,i] = map_coordinates(imagen, (points + (normals * xc) + (direcciones[:,direccion_index,:] * stepsC[i])).T, output_type = None, output = None, order = 1, mode = 'constant', cval = 0.0, prefilter = False)

##        # ######## filtro gaussiano 1D##########
##        radio = 3
##        s = sqrt(-((radio)**2)/(log(0.01)*2))
##        samplesImage_S = zeros(samplesImage.shape, 'float32')
##        ndimageFilters.gaussian_filter1d(samplesImage, sigma = s, axis = 1, order = 0, output = samplesImage_S, mode = 'constant', cval = 0.0)
##        samplesImage = samplesImage_S
##        del samplesImage_S

        if calculoGradiente == 'MuestreoDeLineas':
            # ############ 
            # hacia atras
            dex = array([-1,0,1], 'float32')
            dex = dex.reshape(1,3)        
            F = ndimageFilters.convolve(samplesImage, dex, output = None, mode = 'nearest', cval = 0.0, origin = 0)
        
##        F = contrast3D(F, -1., 1.)

##        F[((samplesImage < extraData['Totsu_grisOrig']*0.1) + (samplesImage > (extraData['u_nivel_grisOrig'][1] + 3.*extraData['s_nivel_grisOrig'][1]))).nonzero()] = 0
##        F[( (samplesImage > (extraData['u_nivel_grisOrig'][1] + 3.*extraData['s_nivel_grisOrig'][1]))).nonzero()] = 0
##        F[((samplesImage < extraData['Totsu_grisOrig']*0.1)).nonzero()] = 0        
        
##        MBx0 = abs(MBx0)
                   
##        indice = range(1, dSample*2+1)
        indice2 = range(int(dCalculo/sampleStep), int(dCalculo/sampleStep)+int(dSample/sampleStep)*2+1)
        indice = range(int(dSample/sampleStep)*2+1)
        data = zeros(31, 'float32')
        xTarIndex = ones(points.shape[0], 'int') * int(dSample/sampleStep)
        G = zeros((samplesImage.shape[0], len(steps)), 'float32')

##        for i in range(samplesImage.shape[0]):
####            if i == 9638:
####                print 'alto'
##            for j, j2 in zip(indice, indice2):
##
##                    if samplesImage[i,j2] < 200: # del liquido o del craneo
##                            if j > 2:
##                                xTarIndex[i] = F[i,j-2:j+3].argmax() + j - 1
##                            break
        aux = samplesImage.shape[1]
        aux_index = range(aux)
        listPointsForward = []
        listPointsBackward = []

        
        list1=[]
        list2=[]
        list3=[]
        list4=[]
        for i in range(samplesImage.shape[0]):
            
##            xTarIndex[i] = (F[i,:] - D * (steps2 * steps2)).argmax()
##            # ##### agregado 17-11-2010
####            if samplesImage[i][dSample-1:dSample+2].mean() > extraData['u_nivel_grisOrig'][0]:
####            if samplesImage[i][dSample-1:dSample+2].mean() > extraData['Totsu_grisOrig']:
##            mb = samplesImage[i][:dSample+1].max()
##            if samplesImage[i][dSample-4:dSample+1].mean() < (mb*0.66):
##                xTarIndex[i] = dSample + 6
##            # #########################

            Imini = samplesImage[i][int(centerIndex-(parameters['def3_dmin']/sampleStep)):centerIndex + 1 + (1./sampleStep)].min() # CORREGIDO 2-2-2011
            if tipoMuestreoImagen == 'muestreoEnLinea':
                mb = samplesImage[i][:centerIndex+1].max()  # muestreo en linea
            elif tipoMuestreoImagen == 'muestreoConCruces':
                mb = max(samplesImageC[i,:,:].max(), samplesImage[i][:centerIndex+1].max())  # muestreo con cruces   # CORREGIDO 3-2-2011

            if mb < extraData['u_nivel_grisOrig'][2] - 2. * extraData['s_nivel_grisOrig'][2]:
                mb = extraData['u_nivel_grisOrig'][2] - 2. * extraData['s_nivel_grisOrig'][2]

            if Imini > (mb*0.66): # se llego a una parte de intensidad alta
##                if samplesImage[i][:centerIndex + int(5./sampleStep) + 1].max() > (extraData['u_nivel_grisOrig'][2]+2.5*abs(extraData['s_nivel_grisOrig'][2])): # hay una intensidad muy alta adelante -> ojos
                if samplesImage[i][:centerIndex + int(parameters['def3_dmax']/sampleStep) + 1].max() > (mb*1.3) : # hay una intensidad muy alta adelante -> ojos
                    xTarIndex[i] = (F[i,:] - D * (steps2 * steps2)).argmax() #+ 1
                    list2.append(i)
                else:
##                    xTarIndex[i] = centerIndex + int(4/sampleStep)# 4
                    xTarIndex[i] = centerIndex + 1# 4
                    listPointsForward.append(i)
                    list3.append(i)
            elif Imini < extraData['u_nivel_grisOrig'][1]*0.1: # no hay duda que se llego al LCR
##                if samplesImage[i][int(centerIndex-(2./sampleStep)):centerIndex + 1].mean() < (mb*0.5): # esta en LCR asi que se hace retroceder
                if samplesImage[i][int(centerIndex-(parameters['def3_dmean']/sampleStep)):centerIndex + 1].mean() < (extraData['u_nivel_grisOrig'][1] - 8 * extraData['s_nivel_grisOrig'][1]): # esta en LCR asi que se hace retroceder
                    xTarIndex[i] = centerIndex - 1# 4
                    listPointsBackward.append(i)
                else:
                    xTarIndex[i] = (F[i,:] - D * (steps2 * steps2)).argmax()# + 1
                    list1.append(i)
            else:
##                if samplesImage[i][int(centerIndex-(2./sampleStep)):centerIndex + 1].mean() < (mb*0.5): # esta en LCR asi que se hace retroceder
                if samplesImage[i][int(centerIndex-(2./sampleStep)):centerIndex + 1].mean() < (extraData['u_nivel_grisOrig'][1] - 8*extraData['s_nivel_grisOrig'][1]): # esta en LCR asi que se hace retroceder
                    xTarIndex[i] = centerIndex - 1# 4
                    listPointsBackward.append(i)
                else:
                    xTarIndex[i] = (F[i,:] - D * (steps2 * steps2)).argmax() #+ 1
                    list4.append(i)
                
                
        if 0: #dibujo
            from pylab import imshow, figure, plot, colorbar, gray, hold, show
            rango = [6600,6660]
            figure()
            imshow(samplesImage[rango[0]:rango[1]], interpolation='nearest', origin='lower')
            colorbar()
            gray()
            hold('on')
            plot(xTarIndex[rango[0]:rango[1]] + dCalculo , range(rango[1]-rango[0]))
            show()
            
            figure()
            imshow(F[rango[0]:rango[1]], interpolation='nearest', origin='lower')
            colorbar()
            gray()  
            show()
            
            
        listPointsBackward = set(listPointsBackward)
        listPointsBackward = list(listPointsBackward)        
        xTar[:,:] = points + normals * (steps2[xTarIndex]).reshape(-1,1)
        xTar[listPointsForward,:] = points[listPointsForward,:] + normals[listPointsForward,:] * parameters['def3_dp'] # 0.1
        xTar[listPointsBackward,:] = points[listPointsBackward,:] - normals[listPointsBackward,:] * parameters['def3_dp']


##        # ############### para Visualizar
##        PtoVTK = PythonToPolyData()
##        PtoVTK.SetInputPoints(xTar)
##        PtoVTK.WritePolyData('..\\xTar.vtk')
##        # ###############
##        
##        # ############### para Visualizar
##        PtoVTK = PythonToPolyData()
##        PtoVTK.SetInputPoints(points[list1])
##        PtoVTK.WritePolyData('..\\L1.vtk')
##        # ###############
##                # ############### para Visualizar
##        PtoVTK = PythonToPolyData()
##        PtoVTK.SetInputPoints(points[list2])
##        PtoVTK.WritePolyData('..\\L2.vtk')
##        # ###############
##        # ############### para Visualizar
##        PtoVTK = PythonToPolyData()
##        PtoVTK.SetInputPoints(points[list3])
##        PtoVTK.WritePolyData('..\\L3.vtk')
##        # ###############
##        # ############### para Visualizar
##        PtoVTK = PythonToPolyData()
##        PtoVTK.SetInputPoints(points[list4])
##        PtoVTK.WritePolyData('..\\L4.vtk')
##        # ###############

        
##        gradient_aux[:,1] = map_coordinates(dx, xTar.T, output_type = None, output = None, order = 1, mode = 'nearest', cval = 0.0, prefilter = False)
##        gradient_aux[:,0] = map_coordinates(dy, xTar.T, output_type = None, output = None, order = 1, mode = 'nearest', cval = 0.0, prefilter = False)
##        gradient_aux[:,2] = map_coordinates(dz, xTar.T, output_type = None, output = None, order = 1, mode = 'nearest', cval = 0.0, prefilter = False)
##            
##        normGradient_aux = sqrt((gradient_aux * gradient_aux).sum(1)).reshape(-1,1)
##            
##        Eext = ((gradient_aux / normGradient_aux) * (xTar - points)).sum(1).reshape(-1,1)
        
##        Eext[normGradient_aux == 0.] = 0.

        desp = (xTar - points)
        despNorm = sqrt((desp * desp).sum(1)).reshape(-1,1)

        # aplicacion de funcion de stiffness
        D = parameters['def3_Df'] # distancia sobre la cual la fueza decrece exponencialmente
        despNormD = despNorm - D
        resp = zeros(despNorm.shape, despNorm.dtype)
        cond1 = (despNormD<0.)
##        cond1 = less(despNormD, 1.)
        resp[cond1] = 1.
        resp[~cond1] = 1. / exp(despNormD[~cond1])

        Fext[:,:] = resp * desp

        
##        Fext[:,:] = xTar - points
##        Fext[:,:] = normals * Eext

        return Fext
    
    def calculateFextVentri1erDef(dSample, gmax, points, normals, imageData):

        dx = imageData['dx']
        dy = imageData['dy']
        dz = imageData['dz']
        imagen_bin = ['imagen_bin']
        
        dSampleIn = 8
        dSampleOut = 8
        
        steps = arange(-dSampleIn,dSampleOut+1)
        
        gradient_aux = zeros((points.shape[0], 3), 'float32')
        xTar = zeros((points.shape[0], 3), 'float32')
        samplesImage = zeros((points.shape[0], len(steps)), 'float32')
        samplesImage_Borrar = zeros((simplex.shape[0], len(steps)), 'float32') # BORRAR

        
        for i in range(len(steps)):
            gradient_aux[:,1] = map_coordinates(dx, (points + (normals * steps[i])).T, output_type = None, output = None, order = 1, mode = 'constant', cval = 0.0, prefilter = False)
            gradient_aux[:,0] = map_coordinates(dy, (points + (normals * steps[i])).T, output_type = None, output = None, order = 1, mode = 'constant', cval = 0.0, prefilter = False)
            gradient_aux[:,2] = map_coordinates(dz, (points + (normals * steps[i])).T, output_type = None, output = None, order = 1, mode = 'constant', cval = 0.0, prefilter = False)
            
            normGradient_aux = (gradient_aux * gradient_aux).sum(1)
            
##                G = (gradient_aux * normals).sum(1) * ((gmax * (gmax + normGradient_aux)) / (gmax * gmax + normGradient_aux * normGradient_aux))
            
            G = (gradient_aux * normals).sum(1)
            
            # criterios de eliminacion
            # imponer una regla para que solo quede el borde mas lejano (para que se segmente la envolvente)
##                F = G.copy()
            F = G.copy()
            F[ G < 0.] = 0.
            # ####
            
            samplesImage[:,i] =  F


        if 0:
            # Dibujar normales Simplex
            points_aux = concatenate((points, points + (normals * steps[0])), 0)
            line_aux = []
            N_points_aux = points.shape[0]
            for iii in range(N_points_aux):
                line_aux.append([iii, iii + N_points_aux])
            PtoVTK = PythonToPolyData()
            PtoVTK.SetInputLines(points_aux, line_aux)
            PtoVTK.WritePolyData('..\\per_ini.vtk')
            # ####

        if 0:
            # Dibujar normales Simplex
            points_aux = concatenate((points, points + (normals * steps[-1])), 0)
            line_aux = []
            N_points_aux = points.shape[0]
            for iii in range(N_points_aux):
                line_aux.append([iii, iii + N_points_aux])
            PtoVTK = PythonToPolyData()
            PtoVTK.SetInputLines(points_aux, line_aux)
            PtoVTK.WritePolyData('..\\per_fin.vtk')
            # ####


        indice = range(steps.shape[0])
        indice.reverse()
        for i in range(samplesImage.shape[0]):
            if samplesImage[i,:].max() > 0.000000001:
                aux = samplesImage[i,:].max() * 0.5
                for j in indice:
                    if samplesImage[i,j] > aux:
                        index_aux = j - dSampleIn
                        break
            else:
                index_aux = 0                
            xTar[i,:] = points[i] + normals[i] * index_aux

##        # ############### para Visualizar
##        PtoVTK = PythonToPolyData()
##        PtoVTK.SetInputPoints(xTar)
##        PtoVTK.WritePolyData('..\\Borrar0.vtk')
##        # ###############

        # ############### para Visualizar
        PtoVTK = PythonToPolyData()
        PtoVTK.SetInputPoints(xTar)
        PtoVTK.WritePolyData('..\\xTar.vtk')
        # ###############

        
        gradient_aux[:,1] = map_coordinates(dx, xTar.T, output_type = None, output = None, order = 1, mode = 'constant', cval = 0.0, prefilter = False)
        gradient_aux[:,0] = map_coordinates(dy, xTar.T, output_type = None, output = None, order = 1, mode = 'constant', cval = 0.0, prefilter = False)
        gradient_aux[:,2] = map_coordinates(dz, xTar.T, output_type = None, output = None, order = 1, mode = 'constant', cval = 0.0, prefilter = False)
            
        normGradient_aux = sqrt((gradient_aux * gradient_aux).sum(1)).reshape(-1,1)
            
        Eext = ((gradient_aux / normGradient_aux) * (xTar - points)).sum(1).reshape(-1,1)
        
        Eext[normGradient_aux == 0.] = 0.
        
##            Fext[:,:] = xTar - points
        Fext[:,:] = normals * Eext



        desp = (xTar - points) # MODIFICADO 07-02-2011
        despNorm = sqrt((desp * desp).sum(1)).reshape(-1,1)
        # aplicacion de funcion de stiffness
        D = 20. # distancia sobre la cual la fueza decrece exponencialmente
        despNormD = despNorm - D
        resp = zeros(despNorm.shape, despNorm.dtype)
        cond1 = (despNormD<0.)
##        cond1 = less(despNormD, 1.)
        resp[cond1] = 1.
        resp[~cond1] = 1. / exp(despNormD[~cond1])

        Fext[:,:] = Fext[:,:] * resp


        return Fext


    def calculateFextVentri2aDef(dSample, gmax, points, normals, imageData):
        calculoGradiente = 'fuerzasDeEntrada' # 'fuerzasDeEntrada', 'MuestreoDeLineas'
        dSample = 8
        sampleStep = 0.5 # 1.
        dCalculo = 0
        porce = 0.8 # 0.8 porcentaje (suponiendo que el maximo F es 1) que se restara al valor del gradiente al final de la linea de muestreo para encontrar el Xtarget
##        D = porce/(8.*8.) #0.5 8.
        D = (porce*22.)/(8.*8.) #0.5 8. # el 22 esta agregado porque ese es el valor maximo que puede tomar el gradiante al normalizar el mapa de bordes en [0,1] y aplicar el filtro sobel
        gmax = (dSample*dSample * D) * 2.

        imagen = imageData['imagen']
        dx = imageData['dx']
        dy = imageData['dy']
        dz = imageData['dz']
        
        gradient_aux = zeros((points.shape[0], 3), 'float32')
        xTar = zeros((points.shape[0], 3), 'float32')
        
        steps = arange(-dSample,dSample+sampleStep, sampleStep)
        steps2 = arange(-(dSample+dCalculo),(dSample+dCalculo)+ sampleStep, sampleStep)
        centerIndex = int(dSample/sampleStep)
        samplesImage = zeros((points.shape[0], len(steps2)), 'float32')
        F = zeros((points.shape[0], len(steps)), 'float32')

        if calculoGradiente == 'fuerzasDeEntrada':
            for i in range(len(steps)):
                gradient_aux[:,1] = map_coordinates(dx, (points + (normals * steps[i])).T, output_type = None, output = None, order = 2, mode = 'nearest', cval = 0.0, prefilter = False)
                gradient_aux[:,0] = map_coordinates(dy, (points + (normals * steps[i])).T, output_type = None, output = None, order = 2, mode = 'nearest', cval = 0.0, prefilter = False)
                gradient_aux[:,2] = map_coordinates(dz, (points + (normals * steps[i])).T, output_type = None, output = None, order = 2, mode = 'nearest', cval = 0.0, prefilter = False)

##                normGradient_aux = sqrt((gradient_aux * gradient_aux).sum(1))
##
##    ##            G = (gradient_aux * normals).sum(1)
##                G = (gradient_aux * normals).sum(1) * ((gmax * (gmax + normGradient_aux)) / (gmax * gmax + normGradient_aux * normGradient_aux))
##
##                # criterios de eliminacion
##                # imponer una regla para que solo quede el borde mas lejano (para que se segmente la envolvente)
##    ##                F = G.copy()
##                F_aux = G.copy()
##                F_aux[ G < 0.] = 0.
##                # ####
##                F[:,i] = F_aux
                F[:,i] = (gradient_aux * normals).sum(1)
                
        for i in range(len(steps)):
            samplesImage[:,i] = map_coordinates(imagen, (points + (normals * steps2[i])).T, output_type = None, output = None, order = 1, mode = 'nearest', cval = 0.0, prefilter = False)

        if calculoGradiente == 'MuestreoDeLineas':
            # ######## filtro gaussiano 1D##########
            radio = 3
            s = sqrt(-((radio)**2)/(log(0.01)*2))
            samplesImage_S = zeros(samplesImage.shape, 'float32')
            ndimageFilters.gaussian_filter1d(samplesImage, sigma = s, axis = 1, order = 0, output = samplesImage_S, mode = 'constant', cval = 0.0)
            samplesImage = samplesImage_S
            del samplesImage_S

            # ############ para detectar los pics
            # hacia atras
            dex = array([1,0,-1], 'float32')
            dex = dex.reshape(1,3)
            F = ndimageFilters.convolve(samplesImage, dex, output = None, mode = 'nearest', cval = 0.0, origin = 0)
            
                   
##        indice = range(1, dSample*2+1)
        indice2 = range(int(dCalculo/sampleStep), int(dCalculo/sampleStep)+int(dSample/sampleStep)*2+1)
        indice = range(int(dSample/sampleStep)*2+1)
        data = zeros(31, 'float32')
        xTarIndex = ones(points.shape[0], 'int') * int(dSample/sampleStep)
        G = zeros((samplesImage.shape[0], len(steps)), 'float32')

        listPointsForward = []
        listPointsBackward = []


        for i in range(samplesImage.shape[0]):
            
            Imean = samplesImage[i][centerIndex:centerIndex + 1 + (2./sampleStep)].mean()
    
            if Imean < (extraData['u_nivel_grisOrig'][1] - 8 * extraData['s_nivel_grisOrig'][1]): # esta en LCR asi que se hace avanzar
                listPointsForward.append(i)
##            elif Imean > (extraData['u_nivel_grisOrig'][2] - 4 * extraData['s_nivel_grisOrig'][1]): # esta en la materia blanca asi que se hace retroceder
##                listPointsBackward.append(i)
            else:
                xTarIndex[i] = (F[i,:] - D * (steps2 * steps2)).argmax()# + 1


        if 0: #dibujo
            from pylab import imshow, figure, plot, colorbar, gray, hold, show
            rango = [6600,6660]
            figure()
            imshow(samplesImage[rango[0]:rango[1]], interpolation='nearest', origin='lower')
            colorbar()
            gray()
            hold('on')
            plot(xTarIndex[rango[0]:rango[1]] + dCalculo , range(rango[1]-rango[0]))
            show()
            
            figure()
            imshow(F[rango[0]:rango[1]], interpolation='nearest', origin='lower')
            colorbar()
            gray()
            show()
            
            
        xTar[:,:] = points + normals * (steps2[xTarIndex]).reshape(-1,1)
        xTar[listPointsForward,:] = points[listPointsForward,:] + normals[listPointsForward,:] * 0.2
        xTar[listPointsBackward,:] = points[listPointsBackward,:] - normals[listPointsBackward,:] * 0.2

        # ############### para Visualizar
        PtoVTK = PythonToPolyData()
        PtoVTK.SetInputPoints(xTar)
        PtoVTK.WritePolyData('..\\Borrar0.vtk')
        # ###############
        
##        gradient_aux[:,1] = map_coordinates(dx, xTar.T, output_type = None, output = None, order = 1, mode = 'nearest', cval = 0.0, prefilter = False)
##        gradient_aux[:,0] = map_coordinates(dy, xTar.T, output_type = None, output = None, order = 1, mode = 'nearest', cval = 0.0, prefilter = False)
##        gradient_aux[:,2] = map_coordinates(dz, xTar.T, output_type = None, output = None, order = 1, mode = 'nearest', cval = 0.0, prefilter = False)
##            
##        normGradient_aux = sqrt((gradient_aux * gradient_aux).sum(1)).reshape(-1,1)
##            
##        Eext = ((gradient_aux / normGradient_aux) * (xTar - points)).sum(1).reshape(-1,1)
        
##        Eext[normGradient_aux == 0.] = 0.

        desp = (xTar - points)
        despNorm = sqrt((desp * desp).sum(1)).reshape(-1,1)

        # aplicacion de funcion de stiffness
        D = 1. # distancia sobre la cual la fueza decrece exponencielmente
        despNormD = despNorm - D
        resp = zeros(despNorm.shape, despNorm.dtype)
        cond1 = (despNormD<0.)
        resp[cond1] = 1.
        resp[~cond1] = 1. / exp(despNormD[~cond1])

        Fext[:,:] = resp * desp
                
        
##        Fext[:,:] = xTar - points
##        Fext[:,:] = normals * Eext

        return Fext

    def calculateFextVentri3aDef(dSample, gmax, points, normals, imageData):
        calculoGradiente = 'fuerzasDeEntrada' # 'fuerzasDeEntrada', 'MuestreoDeLineas'
        dSample = 8
        sampleStep = 0.5 # 1.
        dCalculo = 0
        porce = 0.8 # 0.8 porcentaje (suponiendo que el maximo F es 1) que se restara al valor del gradiente al final de la linea de muestreo para encontrar el Xtarget
##        D = porce/(8.*8.) #0.5 8.
        D = (porce*22.)/(8.*8.) #0.5 8. # el 22 esta agregado porque ese es el valor maximo que puede tomar el gradiante al normalizar el mapa de bordes en [0,1] y aplicar el filtro sobel
        gmax = (dSample*dSample * D) * 2.

        imagen = imageData['imagen']
        dx = imageData['dx']
        dy = imageData['dy']
        dz = imageData['dz']
        
        gradient_aux = zeros((points.shape[0], 3), 'float32')
        xTar = zeros((points.shape[0], 3), 'float32')
        
        steps = arange(-dSample,dSample+sampleStep, sampleStep)
        steps2 = arange(-(dSample+dCalculo),(dSample+dCalculo)+ sampleStep, sampleStep)
        centerIndex = int(dSample/sampleStep)
        samplesImage = zeros((points.shape[0], len(steps2)), 'float32')
        F = zeros((points.shape[0], len(steps)), 'float32')

        if calculoGradiente == 'fuerzasDeEntrada':
            for i in range(len(steps)):
                gradient_aux[:,1] = map_coordinates(dx, (points + (normals * steps[i])).T, output_type = None, output = None, order = 2, mode = 'nearest', cval = 0.0, prefilter = False)
                gradient_aux[:,0] = map_coordinates(dy, (points + (normals * steps[i])).T, output_type = None, output = None, order = 2, mode = 'nearest', cval = 0.0, prefilter = False)
                gradient_aux[:,2] = map_coordinates(dz, (points + (normals * steps[i])).T, output_type = None, output = None, order = 2, mode = 'nearest', cval = 0.0, prefilter = False)

##                normGradient_aux = sqrt((gradient_aux * gradient_aux).sum(1))
##
##    ##            G = (gradient_aux * normals).sum(1)
##                G = (gradient_aux * normals).sum(1) * ((gmax * (gmax + normGradient_aux)) / (gmax * gmax + normGradient_aux * normGradient_aux))
##
##                # criterios de eliminacion
##                # imponer una regla para que solo quede el borde mas lejano (para que se segmente la envolvente)
##    ##                F = G.copy()
##                F_aux = G.copy()
##                F_aux[ G < 0.] = 0.
##                # ####
##                F[:,i] = F_aux
                F[:,i] = (gradient_aux * normals).sum(1)
                
        for i in range(len(steps)):
            samplesImage[:,i] = map_coordinates(imagen, (points + (normals * steps2[i])).T, output_type = None, output = None, order = 1, mode = 'nearest', cval = 0.0, prefilter = False)

        if calculoGradiente == 'MuestreoDeLineas':
            # ######## filtro gaussiano 1D##########
            radio = 3
            s = sqrt(-((radio)**2)/(log(0.01)*2))
            samplesImage_S = zeros(samplesImage.shape, 'float32')
            ndimageFilters.gaussian_filter1d(samplesImage, sigma = s, axis = 1, order = 0, output = samplesImage_S, mode = 'constant', cval = 0.0)
            samplesImage = samplesImage_S
            del samplesImage_S

            # ############ para detectar los pics
            # hacia atras
            dex = array([1,0,-1], 'float32')
            dex = dex.reshape(1,3)
            F = ndimageFilters.convolve(samplesImage, dex, output = None, mode = 'nearest', cval = 0.0, origin = 0)
            
                   
##        indice = range(1, dSample*2+1)
        indice2 = range(int(dCalculo/sampleStep), int(dCalculo/sampleStep)+int(dSample/sampleStep)*2+1)
        indice = range(int(dSample/sampleStep)*2+1)
        data = zeros(31, 'float32')
        xTarIndex = ones(points.shape[0], 'int') * int(dSample/sampleStep)
        G = zeros((samplesImage.shape[0], len(steps)), 'float32')

        listPointsForward = []
        listPointsBackward = []


        for i in range(samplesImage.shape[0]):
##            Imean = samplesImage[i][centerIndex:centerIndex + 1 + (2./sampleStep)].mean()
##    
##            if Imean < (extraData['u_nivel_grisOrig'][1] - 8 * extraData['s_nivel_grisOrig'][1]): # esta en LCR asi que se hace avanzar
##                listPointsForward.append(i)
####            elif Imean > (extraData['u_nivel_grisOrig'][2] - 4 * extraData['s_nivel_grisOrig'][1]): # esta en la materia blanca asi que se hace retroceder
####                listPointsBackward.append(i)
##            else:
##                xTarIndex[i] = (F[i,:] - D * (steps2 * steps2)).argmax()# + 1
##                
            xTarIndex[i] = (F[i,:] - D * (steps2 * steps2)).argmax()# + 1

        if 0: #dibujo
            from pylab import imshow, figure, plot, colorbar, gray, hold, show
            rango = [6600,6660]
            figure()
            imshow(samplesImage[rango[0]:rango[1]], interpolation='nearest', origin='lower')
            colorbar()
            gray()
            hold('on')
            plot(xTarIndex[rango[0]:rango[1]] + dCalculo , range(rango[1]-rango[0]))
            show()
            
            figure()
            imshow(F[rango[0]:rango[1]], interpolation='nearest', origin='lower')
            colorbar()
            gray()
            show()
            
            
        xTar[:,:] = points + normals * (steps2[xTarIndex]).reshape(-1,1)
##        xTar[listPointsForward,:] = points[listPointsForward,:] + normals[listPointsForward,:] * 0.2
##        xTar[listPointsBackward,:] = points[listPointsBackward,:] - normals[listPointsBackward,:] * 0.2

        # ############### para Visualizar
        PtoVTK = PythonToPolyData()
        PtoVTK.SetInputPoints(xTar)
        PtoVTK.WritePolyData('..\\Borrar0.vtk')
        # ###############
        
##        gradient_aux[:,1] = map_coordinates(dx, xTar.T, output_type = None, output = None, order = 1, mode = 'nearest', cval = 0.0, prefilter = False)
##        gradient_aux[:,0] = map_coordinates(dy, xTar.T, output_type = None, output = None, order = 1, mode = 'nearest', cval = 0.0, prefilter = False)
##        gradient_aux[:,2] = map_coordinates(dz, xTar.T, output_type = None, output = None, order = 1, mode = 'nearest', cval = 0.0, prefilter = False)
##            
##        normGradient_aux = sqrt((gradient_aux * gradient_aux).sum(1)).reshape(-1,1)
##            
##        Eext = ((gradient_aux / normGradient_aux) * (xTar - points)).sum(1).reshape(-1,1)
        
##        Eext[normGradient_aux == 0.] = 0.

        desp = (xTar - points)
        despNorm = sqrt((desp * desp).sum(1)).reshape(-1,1)

        # aplicacion de funcion de stiffness
        D = 1. # distancia sobre la cual la fueza decrece exponencielmente
        despNormD = despNorm - D
        resp = zeros(despNorm.shape, despNorm.dtype)
        cond1 = (despNormD<0.)
        resp[cond1] = 1.
        resp[~cond1] = 1. / exp(despNormD[~cond1])

        Fext[:,:] = resp * desp
                
        
##        Fext[:,:] = xTar - points
##        Fext[:,:] = normals * Eext

        return Fext

    simplexMesh.ComputeNormalsOfPoints()
    
    simplex = simplexMesh.points
    neighbors = simplexMesh.neighbors
    facesSimplex = simplexMesh.faces
    meshesContour = simplexMesh.contoursPoints
    normals = simplexMesh.normalsOfPoints
    
       
    simplexA = zeros(simplex.shape, 'float32')
    simplexB = zeros(simplex.shape, 'float32')
    simplexA[:,:] = simplex

    count_verifMov = 0
    iterToVerifMov = 5
    simplexAA = []
    for i in range(iterToVerifMov): simplexAA.append(zeros(simplex.shape, 'float32'))

    N = simplex.shape[0];
    desp_auxAnt = zeros((N,3), 'float32') + 1.
##    axisInterpx = arange(Fx.shape[0])
##    axisInterpy = arange(Fy.shape[1])
##    axisInterpz = arange(Fz.shape[2])
##
##    Ifx = InterpolatingFunction((axisInterpx,axisInterpy,axisInterpz), Fx)
##    Ify = InterpolatingFunction((axisInterpx,axisInterpy,axisInterpz), Fy)
##    Ifz = InterpolatingFunction((axisInterpx,axisInterpy,axisInterpz), Fz)

    #invAI = inv(gamma * diag(ones(1,N)));
    ep = ones(neighbors.shape) * (1./3.)
    epContours = []
    if len(meshesContour)>0:
        epContours = []
        for i in range(len(meshesContour)):
            epContours.append(ones( (meshesContour[i].shape[0],2), 'float32' ) * 0.5)
    activateComputeParameters = 0
    j = recomputeParameters

    Fext = zeros(simplex.shape, 'float32')
    Fint = zeros(simplex.shape, 'float32')

##    # ############### para Visualizar
##    PtoVTK = PythonToPolyData()
##    PtoVTK.SetInputSimplex(simplexMesh)
##    PtoVTK.WritePolyData('..\\Borrar1.vtk')
##    # ###############

    simplexMesh.ComputeNormalsOfPoints()
    normals = simplexMesh.normalsOfPoints
    normals_A = simplexMesh.normalsOfPoints.copy()

    
    count_iterCorrection = 0
    corregir_anomalia = 0
    pointsToModify = []
    pointsToModify_N = []

##    trianglesMeshToVerifIntersec = simplexMesh.SimplexToTrianglesNotDualConservacion()
    trianglesMeshToVerifIntersec = simplexMesh.SimplexToTriangles4()
    if calculoZonaEncerrada:
        trianglesMeshToVerifIntersec.SearchNeighborsUp3()
        trianglesMeshToVerifIntersec_neighbors = copy.deepcopy(trianglesMeshToVerifIntersec.neighbors)
        trianglesMeshToVerifIntersec.SearchMeshes()
        trianglesMeshToVerifIntersec_meshes = copy.deepcopy(trianglesMeshToVerifIntersec.meshes)
    
    dibujar_desp = 0
    if dibujar_desp:
        depNormalList = []
        
    TT = timeit.Timer()
    for count in range(ITER):
        T1 = TT.timer()
        
        normals_A[:] = normals[:]        

        if flag_desplegarDatos:    
            print count
        if 0: # calculo normal
            coordenadas = simplex.T
            Fext[:,1] = map_coordinates(Fx, coordenadas, output_type = None, output = None, order = 1, mode = 'constant', cval = 0.0, prefilter = False)
            Fext[:,0] = map_coordinates(Fy, coordenadas, output_type = None, output = None, order = 1, mode = 'constant', cval = 0.0, prefilter = False)
            Fext[:,2] = map_coordinates(Fz, coordenadas, output_type = None, output = None, order = 1, mode = 'constant', cval = 0.0, prefilter = False)
            
            # conservando solo las fuerzas normales a la malla
            normals[:,:] = cross(simplex[neighbors[:,1],:]-simplex[neighbors[:,0],:], simplex[neighbors[:,2],:]-simplex[neighbors[:,0],:])
            normals_norm = reshape(sqrt(sum(normals*normals,1)), (-1, 1))
            normals[:,:] = normals / normals_norm
            ##normals[(normals_norm==0).nonzero()] = 1/sqrt(3)
            
            Fext[:,:] = normals * reshape(sum(Fext * normals, 1), (-1, 1))
        else:

            dSample = 15 # largo del perfil de muestreo (solo hacia un lado)
            gmax = 1.

##            if ((count/200.) - floor(count/200.))==0:                
##                print 'alto'

            if modo == '1erDef': 
                Fext = calculateFext1erDef(dSample, gmax, points = simplex, normals = normals, imageData = forces)
            elif modo == '2aDef':
                Fext = calculateFext2aDef(dSample, gmax, points = simplex, normals = normals, imageData = forces)
            elif modo == '2aDef2':
                Fext = calculateFext2aDef2(dSample, gmax, points = simplex, normals = normals, neighbors = neighbors, imageData = forces, extraData = extraData)
            elif modo == '2aDef3':
                Fext = calculateFext2aDef3(dSample, gmax, points = simplex, normals = normals, neighbors = neighbors, imageData = forces, extraData = extraData)                
            elif modo == '1erDefVentri':
                Fext = calculateFextVentri1erDef(dSample, gmax, points = simplex, normals = normals, imageData = forces)
            elif modo == '2aDefVentri':
                Fext = calculateFextVentri2aDef(dSample, gmax, points = simplex, normals = normals, imageData = forces)
            elif modo == '3aDefVentri':
                Fext = calculateFextVentri3aDef(dSample, gmax, points = simplex, normals = normals, imageData = forces)

                
        j -= 1;
        if j<0:
            j=recomputeParameters
            activateComputeParameters=1

        Fint[:,:] = FintSimplex(simplex,neighbors,normals,S,sigma,ep,activateComputeParameters,facesSimplex, meshesContour,epContours, Scontour, AngleGammaObjetivo, angleSimplexObjetivo) 
               
        activateComputeParameters = 0       
        
        simplexB[:,:] = simplex
        
        simplex[:,:] = simplex + (1-gamma) * (simplex-simplexA) + kappa * Fext + alfa * Fint
        
        simplexA[:,:] = simplexB

        simplexMesh.ComputeNormalsOfPoints()
        normals = simplexMesh.normalsOfPoints

    ##    if count>-1:
    ##        vtkData = MakePolyDataFromSimplex(simplex, facesSimplex)
    ##        VisualizePolyData(vtkData)
        if dibujos_tag:
            # ############### para Visualizar
            PtoVTK = PythonToPolyData()
            PtoVTK.SetInputSimplex(simplexMesh)
            PtoVTK.WritePolyData('..\\Borrar1.vtk')
            # ##############


        # revisar si la normal se invirtio en algun punto (significa cruse de puntos)    
        rev_normals = (normals_A * normals).sum(1)
        rev_normals = (rev_normals < 0)
        if (rev_normals < 0).any():
            if flag_desplegarDatos:
                print 'Normal Inversa'
            pointsToModify_N = rev_normals.nonzero()[0]
##            pointsToModify_N = list(set(neighbors[pointsToModify_N,:].flatten())) # expando en uno los puntos
            corregir_anomalia = 1
            

        # revicion de auto-intersecciones
        count_iterCorrection += 1
        if (count_iterCorrection > iterCorrection) or (count == (ITER-1)):            
##            trianglesMesh = simplexMesh.SimplexToTrianglesNotDualConservacion()
            trianglesMesh = simplexMesh.SimplexToTriangles4()
            trianglesMeshToVerifIntersec.points[:,:] = trianglesMesh.points[:,:] #para no perder las otras cosas que se han calculado en la malla
            autoIntersecciones = AutoIntersectMesh(trianglesMeshToVerifIntersec)
            autoIntersecciones.PrepareIntersection()
            autoIntersecciones.ComputeIntersections()
            
            if dibujos_tag:
                # ############### para Visualizar
                PtoVTK = PythonToPolyData()
                PtoVTK.SetInputTriangles(trianglesMeshToVerifIntersec)
                PtoVTK.WritePolyData('..\\Borrar1_triangleMesh.vtk')
                # ##########################

            if len(autoIntersecciones.pointsOfIntersection) > 0:
                if flag_desplegarDatos:
                    print 'Intersecciones'
                # encontrar todos los triangulos que tienen problemas de intersecciones.
                trianglesWithIntersections = []
                for point_index in arange(len(autoIntersecciones.pointsOfIntersectionToEdges)):
                    for edge_index in autoIntersecciones.pointsOfIntersectionToEdges[point_index]:
                        trianglesWithIntersections.extend(autoIntersecciones.mesh.edgesToTriangles[edge_index])
                for triangle_index in arange(len(autoIntersecciones.trianglesMeshIntersectionPointsInSurface)):
                    if autoIntersecciones.trianglesMeshIntersectionPointsInSurface[triangle_index]:
                        trianglesWithIntersections.append(triangle_index)
                trianglesWithIntersections = list(set(trianglesWithIntersections))
                
##                if len(trianglesWithIntersections) > 0:                                
##                    pointsToModify = set(trianglesMeshToVerifIntersec.triangles[trianglesWithIntersections,:].flatten())
##                    pointsToModify = list(pointsToModify.intersection(range(simplexMesh.points.shape[0]))) # solo los puntos debo adjuntar los puntos que correspoden a puntos de la malla simplex
####                    pointsToModify = list(set(simplexMesh.neighbors[pointsToModify,:].flatten())) # expando en uno los puntos
##                    corregir_anomalia = 1

                if len(trianglesWithIntersections) > 0:
                    pointsToModify = list(trianglesWithIntersections)
                    pointsToModify = list(set(simplexMesh.neighbors[pointsToModify,:].flatten())) # expando en uno los puntos
                    corregir_anomalia = 1
                    
                    if calculoZonaEncerrada:
                     # Buscar areas encerradas
                        trianglesWithIntersections_aux = pointsToModify[:] # para incluir el alargue que hice
                     
                        trianglesMeshToVerifIntersec_neighbors_aux = copy.deepcopy(trianglesMeshToVerifIntersec_neighbors)
                        for one_triangle in trianglesWithIntersections_aux:
                            trianglesMeshToVerifIntersec_neighbors_aux[one_triangle] = array([], 'int')
                        trianglesMeshToVerifIntersec.neighbors = trianglesMeshToVerifIntersec_neighbors_aux
                        trianglesMeshToVerifIntersec.SearchMeshesUp3()

                        new_INtrianglesWithIntersections = set([])
                        for index_oneMesh,oneMesh in enumerate(trianglesMeshToVerifIntersec.meshes):
                            for index_oneIniMesh,oneIniMesh in enumerate(trianglesMeshToVerifIntersec_meshes):
                                if (oneIniMesh == oneMesh[0]).any():
                                    sizeMeshOfZone = len(oneIniMesh)
                                    break                            
                            if len(oneMesh) < sizeMeshOfZone * 0.1:
                                new_INtrianglesWithIntersections.update(list(oneMesh))
                        new_INtrianglesWithIntersections = list(new_INtrianglesWithIntersections)
                        
                        if dibujos_tag:
                            PtoVTK = PythonToPolyData()
                            PtoVTK.SetInputTriangles(TrianglesMesh(trianglesMeshToVerifIntersec.points,trianglesMeshToVerifIntersec.triangles[new_INtrianglesWithIntersections,:]))
                            PtoVTK.WritePolyData('..\\Borrar1_trianglesEncerrados.vtk')

                        # agregar areas encerradas a los puntos
                        if len(new_INtrianglesWithIntersections) > 0:
##                            if tipo_triangulacion == 'dual':
                            pointsToModify = set(pointsToModify)
                            pointsToModify.update(new_INtrianglesWithIntersections)
                            pointsToModify = list(pointsToModify)
                                
##                            elif tipo_triangulacion == 'detallada':
##                                pointsToModify_aux = set(trianglesMeshToVerifIntersec.triangles[list(new_INtrianglesWithIntersections)].flatten())
##                                pointsToModify_aux = list(pointsToModify_aux.intersection(range(simplexMesh.points.shape[0]))) # solo los puntos debo adjuntar los puntos que correspoden a puntos de la malla simplex
##                                pointsToModify = list(set(pointsToModify).update(pointsToModify_aux))        

                    if dibujos_tag:
                        # ############### para Visualizar
                        PtoVTK = PythonToPolyData()
                        PtoVTK.SetInputTriangles(TrianglesMesh(trianglesMeshToVerifIntersec.points,trianglesMeshToVerifIntersec.triangles[trianglesWithIntersections,:]))
                        PtoVTK.WritePolyData('..\\Borrar1_trianglesWithIntersections.vtk')
                        
                        PtoVTK = PythonToPolyData()
                        PtoVTK.SetInputPoints(simplexMesh.points[pointsToModify])
                        PtoVTK.WritePolyData('..\\Borrar1_intersectionPoints.vtk')
                        # ##########################


            count_iterCorrection = 0


        if corregir_anomalia:
            softFint_flag = 0 # suavizado por fuerza interna o por suavizado laplaciano
            
            pointsToModify = set(pointsToModify)
            pointsToModify.union(pointsToModify_N)
            pointsToModify = list(pointsToModify)
            if flag_desplegarDatos:
                print 'puntos a modificar',len(pointsToModify)
            
##            pointsToModify = list(set(simplexMesh.neighbors[pointsToModify,:].flatten())) # expando en uno los puntos
            pointsModified = set(pointsToModify)

            # iteraciones para corregir la malla            
            for inter_correction in range(0):
                # un primer suavizado basado en el centro de gravedad de las caras
                trianglesMesh = simplexMesh.SimplexToTriangles()
                simplexMesh_aux_int = trianglesMesh.TrianglesToSimplex()
                simplex[pointsToModify,:] = simplexMesh_aux_int.points[pointsToModify,:]
                

                angleSimplexObjetivo_int = -1.
                S_int = 2
                recomputeParameters_int = 100  #cada cuantas iteraciones se recalculan los parametros de acomodacion
                sigma_int = 0.1  #parametro para la acomodacion de la malla
                gamma_int = 0.1 #0.8 #biscoso 

                alfa_int = 0.5 # 0.5 (10-5-2011) # 0.7 # 0.5 #0.7 ponderacion fuerza interna
                kappa_int = 0.4 # 0.4  (10-5-2011)  # 0.3 # 0.5 #0.3 #ponderacion fuerza externa
                ##        alfa = 0.3 #0.2 ponderacion fuerza interna
                ##        kappa = 0.5 #0.15 #ponderacion fuerza externa        
                beta_int = 0 #; presion (-infla)

                Scontour_int = 2 # si < 0 el angulo simplex objetivo de cada punto del contorno sera 0
                AngleGammaObjetivo_int = 0 # para contorno

                activateComputeParameters_int = 0
                ep_int = ones(neighbors.shape) * (1./3.)
                ITER_int = 50                
                epContours_int = []
                
                if len(meshesContour)>0:
                    epContours_int = []
                    for i in range(len(meshesContour)):
                        epContours_int.append(ones( (meshesContour[i].shape[0],2), 'float32' ) * 0.5)

                simplexB_int = simplex.copy()
                simplexA_int = simplex.copy()

        ##            TT = timeit.Timer()
                for count_int in range(ITER_int):
                    simplexMesh.ComputeNormalsOfPoints()
                    normals = simplexMesh.normalsOfPoints

                #   [Fint, F_tangencial, F_normal]
        ##                T1 = TT.timer()
                    Fint_int = FintSimplex(simplex,neighbors,normals,S_int,sigma_int,ep_int,activateComputeParameters_int,facesSimplex, meshesContour,epContours_int, Scontour_int, AngleGammaObjetivo_int, angleSimplexObjetivo_int, pointsToModify)
        ##                T2 = TT.timer()
        ##                print 'Tiempo'
        ##                print T2-T1
                    
                    simplexB_int[pointsToModify,:]= simplex[pointsToModify,:]

                    simplex[pointsToModify,:] = simplex[pointsToModify,:] + (1-gamma_int) * (simplex[pointsToModify,:]-simplexA_int[pointsToModify,:]) + alfa_int * Fint_int[pointsToModify,:]

                    simplexA_int[pointsToModify,:]= simplexB_int[pointsToModify,:]

                    desplazamiento_int = sqrt(((simplex-simplexA_int)**2).sum(1))

                    
                    if desplazamiento_int.sum()/N < 0.001:
                        break
                    

                # verificar si ya no hay intersecciones
                
##                trianglesMesh = simplexMesh.SimplexToTrianglesNotDualConservacion() 
                trianglesMesh = simplexMesh.SimplexToTriangles4()
                trianglesMeshToVerifIntersec.points[:,:] = trianglesMesh.points[:,:]
                autoIntersecciones = AutoIntersectMesh(trianglesMeshToVerifIntersec)
                autoIntersecciones.PrepareIntersection()
                autoIntersecciones.ComputeIntersections()     
                if len(autoIntersecciones.pointsOfIntersection) > 0:
                    if flag_desplegarDatos:
                        print 'Intersecciones'
                    # encontrar todos los triangulos que tienen problemas de intersecciones.
                    trianglesWithIntersections = []
                    for point_index in arange(len(autoIntersecciones.pointsOfIntersectionToEdges)):
                        for edge_index in autoIntersecciones.pointsOfIntersectionToEdges[point_index]:
                            trianglesWithIntersections.extend(autoIntersecciones.mesh.edgesToTriangles[edge_index])
                    for triangle_index in arange(len(autoIntersecciones.trianglesMeshIntersectionPointsInSurface)):
                        if autoIntersecciones.trianglesMeshIntersectionPointsInSurface[triangle_index]:
                            trianglesWithIntersections.append(triangle_index)
                    trianglesWithIntersections = list(set(trianglesWithIntersections))
                    
##                    if len(trianglesWithIntersections) > 0:                                
##                        pointsToModify_aux = set(trianglesMeshToVerifIntersec.triangles[trianglesWithIntersections,:].flatten())
##                        pointsToModify_aux = list(pointsToModify_aux.intersection(range(simplexMesh.points.shape[0]))) # solo los puntos debo adjuntar los puntos que correspoden a puntos de la malla simplex
##                        pointsToModify = list(set(pointsToModify + pointsToModify_aux)) # sumo mas puntos por si aparecioeron nuevas intersecciones

                    if len(trianglesWithIntersections) > 0:
                        pointsToModify_aux = trianglesWithIntersections.copy()
                        pointsToModify = list(set(pointsToModify + pointsToModify_aux)) # sumo mas puntos por si aparecioeron nuevas intersecciones
                        
                else:                    
                    pointsToModify = []
                    pointsToModify_N = []
                    corregir_anomalia = 0
                    break



            for inter_correction in range(10):
                if flag_desplegarDatos:
                    print 'cicloCorrec:',inter_correction
                pointsModified.update(pointsToModify)
                if dibujos_tag:
                    # ############### para Visualizar
                    PtoVTK = PythonToPolyData()
                    PtoVTK.SetInputPoints(simplexMesh.points[pointsToModify,:].reshape(-1,3))
                    PtoVTK.WritePolyData('..\\Borrar1_PuntosIntersec.vtk')
                    # ##########################
                if softFint_flag:
                    angleSimplexObjetivo_int = -1.
                    S_int = 2
                    recomputeParameters_int = 100  #cada cuantas iteraciones se recalculan los parametros de acomodacion
                    sigma_int = 0.1  #parametro para la acomodacion de la malla
                    gamma_int = 1. #0.8 #biscoso 

                    alfa_int = 0.5 # 0.5 (10-5-2011) # 0.7 # 0.5 #0.7 ponderacion fuerza interna
                    kappa_int = 0.4 # 0.4  (10-5-2011)  # 0.3 # 0.5 #0.3 #ponderacion fuerza externa
                    ##        alfa = 0.3 #0.2 ponderacion fuerza interna
                    ##        kappa = 0.5 #0.15 #ponderacion fuerza externa        
                    beta_int = 0 #; presion (-infla)

                    Scontour_int = 2 # si < 0 el angulo simplex objetivo de cada punto del contorno sera 0
                    AngleGammaObjetivo_int = 0 # para contorno

                    activateComputeParameters_int = 0
                    ep_int = ones(neighbors.shape) * (1./3.)
                    ITER_int = 50                
                    epContours_int = []
                    
                    if len(meshesContour)>0:
                        epContours_int = []
                        for i in range(len(meshesContour)):
                            epContours_int.append(ones( (meshesContour[i].shape[0],2), 'float32' ) * 0.5)

                    simplexB_int = simplex.copy()
                    simplexA_int = simplex.copy()

            ##            TT = timeit.Timer()
                    for count_int in range(ITER_int):
                        simplexMesh.ComputeNormalsOfPoints()
                        normals = simplexMesh.normalsOfPoints

                    #   [Fint, F_tangencial, F_normal]
            ##                T1 = TT.timer()
                        Fint_int = FintSimplex(simplex,neighbors,normals,S_int,sigma_int,ep_int,activateComputeParameters_int,facesSimplex, meshesContour,epContours_int, Scontour_int, AngleGammaObjetivo_int, angleSimplexObjetivo_int, pointsToModify)
            ##                T2 = TT.timer()
            ##                print 'Tiempo'
            ##                print T2-T1
                        
                        simplexB_int[pointsToModify,:]= simplex[pointsToModify,:]

                        simplex[pointsToModify,:] = simplex[pointsToModify,:] + (1-gamma_int) * (simplex[pointsToModify,:]- simplexA_int[pointsToModify,:]) + alfa_int * Fint_int[pointsToModify,:]

                        simplexA_int[pointsToModify,:]= simplexB_int[pointsToModify,:]

                        desplazamiento_int = sqrt(((simplex-simplexA_int)**2).sum(1))

                        
                        if desplazamiento_int.sum()/N < 0.001:
                            break
                    
                else:                        
                    for count_int in range(50):
                        simplexMesh.points[pointsToModify,:] = simplexMesh.points[simplexMesh.neighbors[pointsToModify,:]].mean(1)

                if dibujos_tag:
                    # ############### para Visualizar
                    PtoVTK = PythonToPolyData()
                    PtoVTK.SetInputSimplex(simplexMesh)
                    PtoVTK.WritePolyData('..\\Borrar1_' + str(inter_correction) + '_IntermediateCorrection.vtk')
                    # ##########################

                trianglesMesh = simplexMesh.SimplexToTriangles4()
                trianglesMeshToVerifIntersec.points[:,:] = trianglesMesh.points[:,:]
                autoIntersecciones = AutoIntersectMesh(trianglesMeshToVerifIntersec)
                autoIntersecciones.PrepareIntersection()
                autoIntersecciones.ComputeIntersections()     

                if len(autoIntersecciones.pointsOfIntersection) > 0:
                    if flag_desplegarDatos:
                        print 'Intersecciones'
                    # encontrar todos los triangulos que tienen problemas de intersecciones.
                    trianglesWithIntersections = []
                    for point_index in arange(len(autoIntersecciones.pointsOfIntersectionToEdges)):
                        for edge_index in autoIntersecciones.pointsOfIntersectionToEdges[point_index]:
                            trianglesWithIntersections.extend(autoIntersecciones.mesh.edgesToTriangles[edge_index])
                    for triangle_index in arange(len(autoIntersecciones.trianglesMeshIntersectionPointsInSurface)):
                        if autoIntersecciones.trianglesMeshIntersectionPointsInSurface[triangle_index]:
                            trianglesWithIntersections.append(triangle_index)
                    trianglesWithIntersections = list(set(trianglesWithIntersections))
                        
                    if len(trianglesWithIntersections) > 0:
                        # para ocupar solo los puntos nuevos expandidos                
                        pointsToModify = list(trianglesWithIntersections)
                        for alargue in range(inter_correction + 1):
                            pointsToModify = list(set(simplexMesh.neighbors[pointsToModify,:].flatten())) # expando en uno los puntos

                        # para ocuopar los puntos nuevos mas los antiguos extendidos                    
        ##                pointsToModify_aux = list(trianglesWithIntersections)                  
        ##                pointsToModify = list(set(pointsToModify + pointsToModify_aux)) # sumo mas puntos por si aparecioeron nuevas intersecciones
                        if dibujos_tag:
                            # ############### para Visualizar
                            PtoVTK = PythonToPolyData()
                            PtoVTK.SetInputPoints(simplexMesh.points[pointsToModify])
                            PtoVTK.WritePolyData('..\\Borrar1_' + str(inter_correction) + '_IntermediateIntersectionPoints.vtk')
                            # ##########################
                            
                        if calculoZonaEncerrada:
                         # Buscar areas encerradas
                            trianglesWithIntersections_aux = pointsToModify[:] # para incluir el alargue que hice
                         
                            trianglesMeshToVerifIntersec_neighbors_aux = copy.deepcopy(trianglesMeshToVerifIntersec_neighbors)
                            for one_triangle in trianglesWithIntersections_aux:
                                trianglesMeshToVerifIntersec_neighbors_aux[one_triangle] = array([], 'int')
                            trianglesMeshToVerifIntersec.neighbors = trianglesMeshToVerifIntersec_neighbors_aux
                            trianglesMeshToVerifIntersec.SearchMeshesUp3()

                            new_INtrianglesWithIntersections = set([])
                            for index_oneMesh,oneMesh in enumerate(trianglesMeshToVerifIntersec.meshes):
                                for index_oneIniMesh,oneIniMesh in enumerate(trianglesMeshToVerifIntersec_meshes):
                                    if (oneIniMesh == oneMesh[0]).any():
                                        sizeMeshOfZone = len(oneIniMesh)
                                        break                            
                                if len(oneMesh) < sizeMeshOfZone * 0.1:
                                    new_INtrianglesWithIntersections.update(list(oneMesh))
                            new_INtrianglesWithIntersections = list(new_INtrianglesWithIntersections)
                            
                            if dibujos_tag:
                                PtoVTK = PythonToPolyData()
                                PtoVTK.SetInputTriangles(TrianglesMesh(trianglesMeshToVerifIntersec.points,trianglesMeshToVerifIntersec.triangles[new_INtrianglesWithIntersections,:]))
                                PtoVTK.WritePolyData('..\\Borrar1_trianglesEncerrados.vtk')

                            # agregar areas encerradas a los puntos
                            if len(new_INtrianglesWithIntersections) > 0:
    ##                            if tipo_triangulacion == 'dual':
                                pointsToModify = set(pointsToModify)
                                pointsToModify.update(new_INtrianglesWithIntersections)
                                pointsToModify = list(pointsToModify)
                                    
    ##                            elif tipo_triangulacion == 'detallada':
    ##                                pointsToModify_aux = set(trianglesMeshToVerifIntersec.triangles[list(new_INtrianglesWithIntersections)].flatten())
    ##                                pointsToModify_aux = list(pointsToModify_aux.intersection(range(simplexMesh.points.shape[0]))) # solo los puntos debo adjuntar los puntos que correspoden a puntos de la malla simplex
    ##                                pointsToModify = list(set(pointsToModify).update(pointsToModify_aux))        

                else:
                    if dibujos_tag:
                        # ############### para Visualizar
                        PtoVTK = PythonToPolyData()
                        PtoVTK.SetInputSimplex(simplexMesh)
                        PtoVTK.WritePolyData('..\\Borrar1_finalCorrection.vtk')
                        # ##########################
                    pointsModified = list(pointsModified)
                    simplexA[pointsModified,:] = simplex[pointsModified,:]
                    corregir_anomalia = 0
                    break


        
        if dibujos_tag:
            # ############### para Visualizar
            PtoVTK = PythonToPolyData()
            PtoVTK.SetInputSimplex(simplexMesh)
            PtoVTK.WritePolyData('..\\Borrar1_OneIterationFinalDeformed.vtk')
            # ##############
        

        if count >= iterToVerifMov:
            desplazamientoRetard = sqrt(((simplex-simplexAA[count_verifMov])**2).sum(1))
        else:
            desplazamientoRetard = array([Inf])
        simplexAA[count_verifMov][:,:] = simplex[:,:]
        
        count_verifMov += 1
        if count_verifMov == iterToVerifMov:
            count_verifMov = 0

        desp_aux = simplex-simplexA
        desplazamiento = sqrt(((simplex-simplexA)**2).sum(1))
        
        desplazamiento_index = (((desp_aux * desp_auxAnt).sum(1)) > 0.).nonzero()[0]
        desplazamiento2 = desplazamiento[desplazamiento_index]
        desp_auxAnt = desp_aux

        if flag_desplegarDatos:
            print 'desplazamiento2 prom =',desplazamiento2.sum()/N        
            print 'desplazamiento prom =',desplazamiento.sum()/N
            print 'desplazamiento retardado prom =',desplazamientoRetard.sum()/N
            print 'desplazamiento max =',desplazamiento.max()

        T2 = TT.timer()
        if flag_desplegarDatos:
            print 'Tiempo:',T2-T1

        if dibujar_desp:
            desp_normal = (normals * desp_aux).sum(1)            
####            width = 1.
##            figure()
##            plot(desp_normal)
####            bar(arange(len(desp_normal)) - (width), desp_normal, width = width, bottom = 0, color = 'g' , yerr = None, xerr = None, ecolor = 'k', capsize = 3)
##            title('Desplazamiento Normal')
##            show()

            # ############### para Visualizar
            PtoVTK = PythonToPolyData()
            PtoVTK.SetInputSimplex(simplexMesh)
            PtoVTK.WritePolyData('..\\' 'iter' + str(count)+ '.vtk')
            # ###############
            depNormalList.append(desp_normal)
            SaveData(workFolder + 'datosDeformacionOs', ('depNormalList',), (depNormalList,))

        
##        if ((desplazamiento.sum() / N ) < thresholdITER ) or count == ITER or isnan(simplex).any():
        if ((desplazamientoRetard.sum() / N) < thresholdITER ) or count == ITER or isnan(simplex).any():

            # revicion de auto-intersecciones
##            trianglesMesh = simplexMesh.SimplexToTrianglesNotDualConservacion()
            trianglesMesh = simplexMesh.SimplexToTriangles4()
            trianglesMeshToVerifIntersec.points[:,:] = trianglesMesh.points[:,:] #para no perder las otras cosas que se han calculado en la malla
            autoIntersecciones = AutoIntersectMesh(trianglesMeshToVerifIntersec)
            autoIntersecciones.PrepareIntersection()
            autoIntersecciones.ComputeIntersections()
            
            if len(autoIntersecciones.pointsOfIntersection) > 0:
                if flag_desplegarDatos:
                    print 'Intersecciones'
                # encontrar todos los triangulos que tienen problemas de intersecciones.
                trianglesWithIntersections = []
                for point_index in arange(len(autoIntersecciones.pointsOfIntersectionToEdges)):
                    for edge_index in autoIntersecciones.pointsOfIntersectionToEdges[point_index]:
                        trianglesWithIntersections.extend(autoIntersecciones.mesh.edgesToTriangles[edge_index])
                for triangle_index in arange(len(autoIntersecciones.trianglesMeshIntersectionPointsInSurface)):
                    if autoIntersecciones.trianglesMeshIntersectionPointsInSurface[triangle_index]:
                        trianglesWithIntersections.append(triangle_index)
                trianglesWithIntersections = list(set(trianglesWithIntersections))
                
##                if len(trianglesWithIntersections) > 0:                                
##                    pointsToModify = set(trianglesMeshToVerifIntersec.triangles[trianglesWithIntersections,:].flatten())
##                    pointsToModify = list(pointsToModify.intersection(range(simplexMesh.points.shape[0]))) # solo los puntos debo adjuntar los puntos que correspoden a puntos de la malla simplex
####                    pointsToModify = list(set(simplexMesh.neighbors[pointsToModify,:].flatten())) # expando en uno los puntos
##                    corregir_anomalia = 1

                if len(trianglesWithIntersections) > 0:
                    pointsToModify = list(trianglesWithIntersections)
                    pointsToModify = list(set(simplexMesh.neighbors[pointsToModify,:].flatten())) # expando en uno los puntos
                    corregir_anomalia = 1
                    
                    if calculoZonaEncerrada:
                     # Buscar areas encerradas
                        trianglesWithIntersections_aux = pointsToModify[:] # para incluir el alargue que hice
                        trianglesMeshToVerifIntersec_neighbors_aux = copy.deepcopy(trianglesMeshToVerifIntersec_neighbors)
                        for one_triangle in trianglesWithIntersections_aux:
                            trianglesMeshToVerifIntersec_neighbors_aux
                            trianglesMeshToVerifIntersec_neighbors_aux[one_triangle] = array([], 'int')
                        trianglesMeshToVerifIntersec.neighbors = trianglesMeshToVerifIntersec_neighbors_aux
                        trianglesMeshToVerifIntersec.SearchMeshesUp3()

                        new_INtrianglesWithIntersections = set([])
                        for index_oneMesh,oneMesh in enumerate(trianglesMeshToVerifIntersec.meshes):
                            for index_oneIniMesh,oneIniMesh in enumerate(trianglesMeshToVerifIntersec_meshes):
                                if (oneIniMesh == oneMesh[0]).any():
                                    sizeMeshOfZone = len(oneIniMesh)
                                    break                            
                            if len(oneMesh) < sizeMeshOfZone * 0.1:
                                new_INtrianglesWithIntersections.update(list(oneMesh))
                        new_INtrianglesWithIntersections = list(new_INtrianglesWithIntersections)
                        
                        if dibujos_tag:
                            PtoVTK = PythonToPolyData()
                            PtoVTK.SetInputTriangles(TrianglesMesh(trianglesMeshToVerifIntersec.points,trianglesMeshToVerifIntersec.triangles[new_INtrianglesWithIntersections,:]))
                            PtoVTK.WritePolyData('..\\Borrar1_trianglesEncerrados.vtk')

                        # agregar areas encerradas a los puntos
                        if len(new_INtrianglesWithIntersections) > 0:
##                            if tipo_triangulacion == 'dual':
                            pointsToModify = set(pointsToModify)
                            pointsToModify.update(new_INtrianglesWithIntersections)
                            pointsToModify = list(pointsToModify)
                                

            count_iterCorrection = 0


            if corregir_anomalia:
                
                softFint_flag = 1 # suavizado por fuerza interna o por suavizado laplaciano
                
                pointsToModify = set(pointsToModify)
                pointsToModify.union(pointsToModify_N)
                pointsToModify = list(pointsToModify)
                if flag_desplegarDatos:
                    print 'puntos a modificar',len(pointsToModify)
                
##                pointsToModify = list(set(simplexMesh.neighbors[pointsToModify,:].flatten())) # expando en uno los puntos


                # iteraciones para corregir la malla
                for inter_correction in range(0):
                    
                    pointsToModify = list(set(simplexMesh.neighbors[pointsToModify,:].flatten())) # expando en uno los puntos
                    # un primer suavizado basado en el centro de gravedad de las caras
                    trianglesMesh = simplexMesh.SimplexToTriangles()
                    simplexMesh_aux_int = trianglesMesh.TrianglesToSimplex()
                    simplex[pointsToModify,:] = simplexMesh_aux_int.points[pointsToModify,:]
                    

                    angleSimplexObjetivo_int = -1.
                    S_int = 2
                    recomputeParameters_int = 100  #cada cuantas iteraciones se recalculan los parametros de acomodacion
                    sigma_int = 0.1  #parametro para la acomodacion de la malla
                    gamma_int = 0.1 #0.8 #biscoso 

                    alfa_int = 0.5 # 0.5 (10-5-2011) # 0.7 # 0.5 #0.7 ponderacion fuerza interna
                    kappa_int = 0.4 # 0.4  (10-5-2011)  # 0.3 # 0.5 #0.3 #ponderacion fuerza externa
                    ##        alfa = 0.3 #0.2 ponderacion fuerza interna
                    ##        kappa = 0.5 #0.15 #ponderacion fuerza externa        
                    beta_int = 0 #; presion (-infla)

                    Scontour_int = 2 # si < 0 el angulo simplex objetivo de cada punto del contorno sera 0
                    AngleGammaObjetivo_int = 0 # para contorno

                    activateComputeParameters_int = 0
                    ep_int = ones(neighbors.shape) * (1./3.)
                    ITER_int = 50                
                    epContours_int = []
                    
                    if len(meshesContour)>0:
                        epContours_int = []
                        for i in range(len(meshesContour)):
                            epContours_int.append(ones( (meshesContour[i].shape[0],2), 'float32' ) * 0.5)

                    simplexB_int = simplex.copy()
                    simplexA_int = simplex.copy()

            ##            TT = timeit.Timer()
                    for count_int in range(ITER_int):
                        simplexMesh.ComputeNormalsOfPoints()
                        normals = simplexMesh.normalsOfPoints

                    #   [Fint, F_tangencial, F_normal]
            ##                T1 = TT.timer()
                        Fint_int = FintSimplex(simplex,neighbors,normals,S_int,sigma_int,ep_int,activateComputeParameters_int,facesSimplex, meshesContour,epContours_int, Scontour_int, AngleGammaObjetivo_int, angleSimplexObjetivo_int, pointsToModify)
            ##                T2 = TT.timer()
            ##                print 'Tiempo'
            ##                print T2-T1
                        
                        simplexB_int[pointsToModify,:]= simplex[pointsToModify,:]

                        simplex[pointsToModify,:] = simplex[pointsToModify,:] + (1-gamma_int) * (simplex[pointsToModify,:]-simplexA_int[pointsToModify,:]) + alfa_int * Fint_int[pointsToModify,:]

                        simplexA_int[pointsToModify,:]= simplexB_int[pointsToModify,:]

                        desplazamiento_int = sqrt(((simplex-simplexA_int)**2).sum(1))

                        
                        if desplazamiento_int.sum()/N < 0.001:
                            break
                        

                    # verificar si ya no hay intersecciones
                    
    ##                trianglesMesh = simplexMesh.SimplexToTrianglesNotDualConservacion() 
                    trianglesMesh = simplexMesh.SimplexToTriangles4()
                    trianglesMeshToVerifIntersec.points[:,:] = trianglesMesh.points[:,:]
                    autoIntersecciones = AutoIntersectMesh(trianglesMeshToVerifIntersec)
                    autoIntersecciones.PrepareIntersection()
                    autoIntersecciones.ComputeIntersections()     
                    if len(autoIntersecciones.pointsOfIntersection) > 0:
                        print 'Intersecciones'
                        # encontrar todos los triangulos que tienen problemas de intersecciones.
                        trianglesWithIntersections = []
                        for point_index in arange(len(autoIntersecciones.pointsOfIntersectionToEdges)):
                            for edge_index in autoIntersecciones.pointsOfIntersectionToEdges[point_index]:
                                trianglesWithIntersections.extend(autoIntersecciones.mesh.edgesToTriangles[edge_index])
                        for triangle_index in arange(len(autoIntersecciones.trianglesMeshIntersectionPointsInSurface)):
                            if autoIntersecciones.trianglesMeshIntersectionPointsInSurface[triangle_index]:
                                trianglesWithIntersections.append(triangle_index)
                        trianglesWithIntersections = list(set(trianglesWithIntersections))
                        
    ##                    if len(trianglesWithIntersections) > 0:                                
    ##                        pointsToModify_aux = set(trianglesMeshToVerifIntersec.triangles[trianglesWithIntersections,:].flatten())
    ##                        pointsToModify_aux = list(pointsToModify_aux.intersection(range(simplexMesh.points.shape[0]))) # solo los puntos debo adjuntar los puntos que correspoden a puntos de la malla simplex
    ##                        pointsToModify = list(set(pointsToModify + pointsToModify_aux)) # sumo mas puntos por si aparecioeron nuevas intersecciones

                        if len(trianglesWithIntersections) > 0:
                            pointsToModify_aux = trianglesWithIntersections.copy()
                            pointsToModify = list(set(pointsToModify + pointsToModify_aux)) # sumo mas puntos por si aparecioeron nuevas intersecciones
                            
                    else:                    
                        pointsToModify = []
                        pointsToModify_N = []
                        corregir_anomalia = 0
                        break



                for inter_correction in range(10):
                    print 'ciclocorrec:',inter_correction

                    if softFint_flag:
                        angleSimplexObjetivo_int = -1.
                        S_int = 2
                        recomputeParameters_int = 100  #cada cuantas iteraciones se recalculan los parametros de acomodacion
                        sigma_int = 0.1  #parametro para la acomodacion de la malla
                        gamma_int = 1. #0.8 #biscoso 

                        alfa_int = 0.5 # 0.5 (10-5-2011) # 0.7 # 0.5 #0.7 ponderacion fuerza interna
                        kappa_int = 0.4 # 0.4  (10-5-2011)  # 0.3 # 0.5 #0.3 #ponderacion fuerza externa
                        ##        alfa = 0.3 #0.2 ponderacion fuerza interna
                        ##        kappa = 0.5 #0.15 #ponderacion fuerza externa        
                        beta_int = 0 #; presion (-infla)

                        Scontour_int = 2 # si < 0 el angulo simplex objetivo de cada punto del contorno sera 0
                        AngleGammaObjetivo_int = 0 # para contorno

                        activateComputeParameters_int = 0
                        ep_int = ones(neighbors.shape) * (1./3.)
                        ITER_int = 50                
                        epContours_int = []
                        
                        if len(meshesContour)>0:
                            epContours_int = []
                            for i in range(len(meshesContour)):
                                epContours_int.append(ones( (meshesContour[i].shape[0],2), 'float32' ) * 0.5)

                        simplexB_int = simplex.copy()
                        simplexA_int = simplex.copy()

                ##            TT = timeit.Timer()
                        for count_int in range(ITER_int):
                            simplexMesh.ComputeNormalsOfPoints()
                            normals = simplexMesh.normalsOfPoints

                        #   [Fint, F_tangencial, F_normal]
                ##                T1 = TT.timer()
                            Fint_int = FintSimplex(simplex,neighbors,normals,S_int,sigma_int,ep_int,activateComputeParameters_int,facesSimplex, meshesContour,epContours_int, Scontour_int, AngleGammaObjetivo_int, angleSimplexObjetivo_int, pointsToModify)
                ##                T2 = TT.timer()
                ##                print 'Tiempo'
                ##                print T2-T1
                            
                            simplexB_int[pointsToModify,:]= simplex[pointsToModify,:]

                            simplex[pointsToModify,:] = simplex[pointsToModify,:] + (1-gamma_int) * (simplex[pointsToModify,:]- simplexA_int[pointsToModify,:]) + alfa_int * Fint_int[pointsToModify,:]

                            simplexA_int[pointsToModify,:]= simplexB_int[pointsToModify,:]

                            desplazamiento_int = sqrt(((simplex-simplexA_int)**2).sum(1))

                            
                            if desplazamiento_int.sum()/N < 0.001:
                                break
                        
                    else:                        
                        for count_int in range(50):
                            simplexMesh.points[pointsToModify,:] = simplexMesh.points[simplexMesh.neighbors[pointsToModify,:]].mean(1)


                    trianglesMesh = simplexMesh.SimplexToTriangles4()
                    trianglesMeshToVerifIntersec.points[:,:] = trianglesMesh.points[:,:]
                    autoIntersecciones = AutoIntersectMesh(trianglesMeshToVerifIntersec)
                    autoIntersecciones.PrepareIntersection()
                    autoIntersecciones.ComputeIntersections()     

                    if len(autoIntersecciones.pointsOfIntersection) > 0:
                        print 'Intersecciones'
                        # encontrar todos los triangulos que tienen problemas de intersecciones.
                        trianglesWithIntersections = []
                        for point_index in arange(len(autoIntersecciones.pointsOfIntersectionToEdges)):
                            for edge_index in autoIntersecciones.pointsOfIntersectionToEdges[point_index]:
                                trianglesWithIntersections.extend(autoIntersecciones.mesh.edgesToTriangles[edge_index])
                        for triangle_index in arange(len(autoIntersecciones.trianglesMeshIntersectionPointsInSurface)):
                            if autoIntersecciones.trianglesMeshIntersectionPointsInSurface[triangle_index]:
                                trianglesWithIntersections.append(triangle_index)
                        trianglesWithIntersections = list(set(trianglesWithIntersections))
                            
                        if len(trianglesWithIntersections) > 0:
                            # para ocupar solo los puntos nuevos expandidos                
                            pointsToModify = list(trianglesWithIntersections)
                            for alargue in range(inter_correction + 1):
                                pointsToModify = list(set(simplexMesh.neighbors[pointsToModify,:].flatten())) # expando en uno los puntos

                            # para ocuopar los puntos nuevos mas los antiguos extendidos                    
            ##                pointsToModify_aux = list(trianglesWithIntersections)                  
            ##                pointsToModify = list(set(pointsToModify + pointsToModify_aux)) # sumo mas puntos por si aparecioeron nuevas intersecciones

                            if calculoZonaEncerrada:
                             # Buscar areas encerradas
                                trianglesWithIntersections_aux = pointsToModify[:] # para incluir el alargue que hice
                             
                                trianglesMeshToVerifIntersec_neighbors_aux = copy.deepcopy(trianglesMeshToVerifIntersec_neighbors)
                                for one_triangle in trianglesWithIntersections_aux:
                                    trianglesMeshToVerifIntersec_neighbors_aux[one_triangle] = array([], 'int')
                                trianglesMeshToVerifIntersec.neighbors = trianglesMeshToVerifIntersec_neighbors_aux
                                trianglesMeshToVerifIntersec.SearchMeshesUp3()

                                new_INtrianglesWithIntersections = set([])
                                for index_oneMesh,oneMesh in enumerate(trianglesMeshToVerifIntersec.meshes):
                                    if len(oneMesh) < trianglesMeshToVerifIntersec.triangles.shape[0] * 0.1:
                                        new_INtrianglesWithIntersections.update(list(oneMesh))
                                new_INtrianglesWithIntersections = list(new_INtrianglesWithIntersections)
                                
                                if dibujos_tag:
                                    PtoVTK = PythonToPolyData()
                                    PtoVTK.SetInputTriangles(TrianglesMesh(trianglesMeshToVerifIntersec.points,trianglesMeshToVerifIntersec.triangles[new_INtrianglesWithIntersections,:]))
                                    PtoVTK.WritePolyData('..\\Borrar1_trianglesEncerrados.vtk')

                                # agregar areas encerradas a los puntos
                                if len(new_INtrianglesWithIntersections) > 0:
        ##                            if tipo_triangulacion == 'dual':
                                    pointsToModify = set(pointsToModify)
                                    pointsToModify.update(new_INtrianglesWithIntersections)
                                    pointsToModify = list(pointsToModify)
                                        
        ##                            elif tipo_triangulacion == 'detallada':
        ##                                pointsToModify_aux = set(trianglesMeshToVerifIntersec.triangles[list(new_INtrianglesWithIntersections)].flatten())
        ##                                pointsToModify_aux = list(pointsToModify_aux.intersection(range(simplexMesh.points.shape[0]))) # solo los puntos debo adjuntar los puntos que correspoden a puntos de la malla simplex
        ##                                pointsToModify = list(set(pointsToModify).update(pointsToModify_aux))        

                            
                    else:
                        break

            break

    return simplex



def SegmentationSteps(dataNames,segmentation_Flags):
    imageFolder = dataNames['imageFolder']
    workFolder = dataNames['workFolder']
    imageFile = dataNames['imageFile']
    prePegmentedImageFile = dataNames['preSegmentedImageFile']   
    
    if segmentation_Flags['flag_segStep_craneo']: # eliminar el craneo
        extras = {}
##        extras['filDataEstimationTs'] = filData
        EliminarCraneo(dataNames, extras)
    if segmentation_Flags['flag_segStep_registro']: # ajustar mallas con transformaciones geometricas       
        dibujos_flag = 0
        desplegarDatos = 0
        transformationAffinFlag = 1

        NameFileImagen = workFolder + prePegmentedImageFile
        NameFileMeshCortex = dataNames['cortexMeshModel']            
        
        # leer malla corteza
        
        tupla=ReadData(NameFileMeshCortex)
        script=''
        for i in range(len(tupla[0])):
            exec(script.join((tupla[0][i],'=tupla[1][i]')))
        if tupla[0] == ('trianglesMesh',):
            trianglesMesh.__init__() # inicializar por si cambie algo en la clase
            cortexSimplexMesh = trianglesMesh.TrianglesToSimplex4()
        else:
            simplexMesh.__init__() # inicializar por si cambie algo en la clase
            cortexSimplexMesh = simplexMesh
        del tupla

        if dibujos_flag:
            # ############### para Visualizar
            PtoVTK = PythonToPolyData()
            PtoVTK.SetInputSimplex(cortexSimplexMesh)
            PtoVTK.WritePolyData('..\\Borrar0.vtk')
            # ###############

        # leer imagen segmentada de cerebro
        hdr = medimages.read_ANALYZE_tags(NameFileImagen)
        imagen = medimages.read_ANALYZE(hdr)


##        # ######################################################
##        cortexSimplexMesh.points =  cortexSimplexMesh.points / array(hdr['Spacing'])
##        cortexSimplexMesh.points = cortexSimplexMesh.points.astype('float32')
##        # ######################################################

        ##vtkData = MakePolyDataFromSimplex(simplexCortex, facesSimplexCortex)
        ##VisualizePolyData(vtkData['vtkData'], 'Wireframe') # 'Wireframe', 'Surface'

        imagen_bin = imagen>0
        indexImagen_bin = imagen_bin.nonzero()
        if 0: # estimacion de primer ajuste simple
            maxMesh = cortexSimplexMesh.points.max(0)
            minMesh = cortexSimplexMesh.points.min(0)
            sz = float(indexImagen_bin[2].max() - indexImagen_bin[2].min()) / (maxMesh[2] - minMesh[2])
            sy = float(indexImagen_bin[1].max() - indexImagen_bin[1].min()) / (maxMesh[1] - minMesh[1])
##            sx = float(indexImagen_bin[0].max() - indexImagen_bin[0].min()) / (maxMesh[0] - minMesh[0])
            sx = (sz + sy) / 2.

            cortexSimplexMesh.points = cortexSimplexMesh.points * array([sx, sy, sz], 'float32') #ajustar tamaño
            ##simplexCortex2 = simplexCortex * array([1., sy, sz], 'float32') #ajustar tamaño

            skullSimplexMesh.points = skullSimplexMesh.points * array([sx, sy, sz], 'float32') #ajustar tamaño

            ventriclesSimplexMesh.points = ventriclesSimplexMesh.points * array([sx, sy, sz], 'float32') #ajustar tamaño
            ##simplexVentricle2 = simplexVentricle * array([1., sy, sz], 'float32') #ajustar tamaño

            tentoriumSimplexMesh.points = tentoriumSimplexMesh.points * array([sx, sy, sz], 'float32') #ajustar tamaño
            hozSimplexMesh.points = hozSimplexMesh.points * array([sx, sy, sz], 'float32') #ajustar tamaño

            c_mesh = cortexSimplexMesh.points.mean(0)
            c_brain = array([indexImagen_bin[0].mean(), indexImagen_bin[1].mean(), indexImagen_bin[2].mean()])
            T = c_brain - c_mesh # ajustar posicion de centro de gravedad

            if desplegarDatos:
                print ' '
                print 'Traslacion central'
                print 'T =',T

            cortexSimplexMesh.points += T
            skullSimplexMesh.points += T
            ventriclesSimplexMesh.points += T
            tentoriumSimplexMesh.points += T
            hozSimplexMesh.points += T
            
        else: # estimacion de primer ajuste con perfiles de la mascara
            
        # #### estimacion de transformacion con perfiles de mascara
##            z_eyeMesh = 40. # calculado a mano desde la malla
            z_eyeMesh = 43. # calculado a mano desde la malla, con modelo de loni
            bbImagen = array([[indexImagen_bin[0].min(), indexImagen_bin[0].max()], [indexImagen_bin[1].min(), indexImagen_bin[1].max()], [indexImagen_bin[2].min(), indexImagen_bin[2].max()]])
            bbMesh = concatenate((cortexSimplexMesh.points.min(0).reshape(-1,1),cortexSimplexMesh.points.max(0).reshape(-1,1)),1)
            lMeshY = float(bbMesh[1,1] - bbMesh[1,0])
            lMeshX = float(bbMesh[2,1] - bbMesh[2,0])

            centX = (bbImagen[2,0] + bbImagen[2,1])/2
            lcent = (bbImagen[2,1] - bbImagen[2,0])/8
            proyection = imagen_bin[:,:,centX-lcent:centX+lcent].sum(2)
            
    ##                    proyection = M1.sum(2)
            proyection_index = proyection.nonzero()
            proyection_profile = zeros(((bbImagen[0,1] - bbImagen[0,0])+1))
            for index_p in range(((bbImagen[0,1] - bbImagen[0,0])+1)):
                try:
                    proyection_profile[index_p] = proyection_index[1][proyection_index[0] == (index_p + bbImagen[0,0])].max()
                except:
                    None
    ##                    limProfile = bbM1[1,1] - (bbM1[1,1] - bbM1[1,0])*0.2
            max_profile = 0
            maxImagenY = 0
            limProfile = -1
            for index_p in arange(len(proyection_profile)-1,-1,-1):
                val_profile = proyection_profile[index_p]
                if max_profile < val_profile:
                    max_profile = val_profile
                    maxImagenY_profileIndex = index_p + bbImagen[0,0]
                    limProfile = max_profile - (max_profile - bbImagen[1,0])*0.2
                if proyection_profile[index_p] < limProfile:
                    maxImagenY = max_profile
                    minImagenY = imagen_bin.sum(2)[maxImagenY_profileIndex,:].nonzero()[0].min()
                    break
                
            z_eyeImagen = index_p + bbImagen[0,0]

            if maxImagenY > 0:
                sy = (maxImagenY - minImagenY) / lMeshY
            else:
                sy = (bbImagen[1,1] - bbImagen[1,0]) / lMeshY
                
            sx = (bbImagen[2,1] - bbImagen[2,0]) / lMeshX
            sz = float(bbImagen[0,1] - z_eyeImagen) / (bbMesh[0,1] - z_eyeMesh)
##            if sz > 1.25:
##                sz = 1.25
##            elif sz < 0.75:
##                sz = 0.75

            # #################

            matrix = array([[sz,0,0],[0,sy,0],[0,0,sx]])
            offset = array([bbImagen[0,1]-sz*bbMesh[0,1],bbImagen[1,1]-sy*bbMesh[1,1],bbImagen[2,1]-sx*bbMesh[2,1]]).astype('float32')

            cortexSimplexMesh.points = cortexSimplexMesh.points * array([sz, sy, sx], 'float32') + offset #ajustar tamaño


        if dibujos_flag:
            # ############### para Visualizar
            PtoVTK = PythonToPolyData()
            PtoVTK.SetInputSimplex(cortexSimplexMesh)
            PtoVTK.WritePolyData('..\\Borrar3.vtk')
            # ###############
            # ############### para Visualizar
            PtoVTK = PythonToPolyData()
            PtoVTK.SetInputSimplex(skullSimplexMesh)
            PtoVTK.WritePolyData('..\\Borrar4.vtk')
            # ###############
            # ############### para Visualizar
            PtoVTK = PythonToPolyData()
            PtoVTK.SetInputSimplex(ventriclesSimplexMesh)
            PtoVTK.WritePolyData('..\\Borrar5.vtk')
            # ###############

        ##simplexCortex2 = simplexCortex2.astype('float32')
        ##simplexVentricle2 = simplexVentricle2.astype('float32')


        # ### calcular la distancia a los bordes ####
        if 1:
            dx = zeros((1,2,1), 'float32'); dy = zeros((2,1,1), 'float32'); dz = zeros((1,1,2), 'float32')
            dy[0,:,:]=1; dy[1,:,:]=-1
            dx[:,0,:]=1; dx[:,1,:]=-1
            dz[:,:,0]=1; dz[:,:,1]=-1

            Dx1 = ndimageFilters.convolve(imagen_bin, dx, output = None, mode = 'reflect', cval = 0.0, origin = 0)
            Dy1 = ndimageFilters.convolve(imagen_bin, dy, output = None, mode = 'reflect', cval = 0.0, origin = 0)
            Dz1 = ndimageFilters.convolve(imagen_bin, dz, output = None, mode = 'reflect', cval = 0.0, origin = 0)

            dx = zeros((1,3,1), 'float32'); dy = zeros((3,1,1), 'float32'); dz = zeros((1,1,3), 'float32')
            dy[0,:,:]=-1; dy[1,:,:]=1
            dx[:,0,:]=-1; dx[:,1,:]=1
            dz[:,:,0]=-1; dz[:,:,1]=1

            Dx2 = ndimageFilters.convolve(imagen_bin, dx, output = None, mode = 'reflect', cval = 0.0, origin = 0)
            Dy2 = ndimageFilters.convolve(imagen_bin, dy, output = None, mode = 'reflect', cval = 0.0, origin = 0)
            Dz2 = ndimageFilters.convolve(imagen_bin, dz, output = None, mode = 'reflect', cval = 0.0, origin = 0)

            Dx = (Dx1>0) + (Dx2>0)
            Dy = (Dy1>0) + (Dy2>0)
            Dz = (Dz1>0) + (Dz2>0)

            D = Dx | Dy | Dz # bordes
            D = D==0
            del Dx, Dy, Dz, Dx1, Dy1, Dz1, Dx2, Dy2, Dz2
            distance = ndimageMorphology.distance_transform_edt(D, sampling = None, return_distances = True, return_indices = False, distances = None, indices = None)
            SaveData('..\DistanciaCabeza1', ('distance',), (distance,))
        # ####
        else:
            tupla=ReadData('..\DistanciaCabeza1')
            script=''
            for i in range(len(tupla[0])):
                exec(script.join((tupla[0][i],'=tupla[1][i]')))
            del tupla
        if transformationAffinFlag:
            parametrosIni = concatenate((eye(3,3).flatten(), zeros(3, 'float32')), 0)
            paramatrosAjust = array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.], 'float32')

        else:    
            paramatrosAjust = array([1., 1., 1., 1., 1., 1., 1., 1., 1.], 'float32')
            parametrosIni = array([0., 0., 0., 1., 1., 1., 0., 0., 0.], 'float32') #[tetax, tetay, tetaz, sx, sy, sz, tx, ty, tz]
            parametrosIni *= concatenate((zeros(3, 'float32'), paramatrosAjust[3:6]**(-1), zeros(3,  'float32')),1)

        # ajustar corteza y con ella todo el resto
        if transformationAffinFlag:

##            D = [0.004, 0.004, 0.004, 0.0005, 0.0005, 0.0005, 1, 1, 1]
            FF = 0.1 #100
            parametrosOpt = leastsq(TransformationAdjustmentAffin, parametrosIni.copy(), args = (cortexSimplexMesh.points.copy(), distance), Dfun = None, full_output = 1, col_deriv = 0, ftol = 1.49012e-008, xtol = 1.49012e-008, gtol = 0.0, maxfev = 0, epsfcn = 0.0, factor = FF, diag = paramatrosAjust)
            parametrosOptimos = parametrosOpt[0]

            if desplegarDatos:            
                print parametrosOptimos 
            
            cortexSimplexMesh.points = AppplyTransformationAffin(cortexSimplexMesh.points.copy(), parametrosOptimos)
            
        else:

            parametrosOpt = leastsq(TransformationAdjustment, parametrosIni.copy(), args = (cortexSimplexMesh.points.copy(), distance, paramatrosAjust), Dfun = None, full_output = 1, col_deriv = 0, ftol = 1.49012e-008, xtol = 1.49012e-008, gtol = 0.0, maxfev = 0, epsfcn = 0.0, factor = 100, diag = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
            
            parametrosOptimos = parametrosOpt[0] * paramatrosAjust
            tetax = parametrosOptimos[0]
            tetay = parametrosOptimos[1]
            tetaz = parametrosOptimos[2]
            sx = parametrosOptimos[3]
            sy = parametrosOptimos[4]
            sz = parametrosOptimos[5]
            tx = parametrosOptimos[6]
            ty = parametrosOptimos[7]
            tz = parametrosOptimos[8]
            R=array([[cos(tetay)*cos(tetaz), -cos(tetay)*sin(tetaz), -sin(tetay), 0.],
                     [-sin(tetax)*sin(tetay)*cos(tetaz)+cos(tetax)*sin(tetaz), sin(tetax)*sin(tetay)*sin(tetaz)+cos(tetax)*cos(tetaz), -sin(tetax)*cos(tetay), 0.],
                     [cos(tetax)*sin(tetay)*cos(tetaz)+sin(tetax)*sin(tetaz), -cos(tetax)*sin(tetay)*sin(tetaz)+sin(tetax)*cos(tetaz), cos(tetax)*cos(tetay), 0.],
                     [0., 0., 0., 1.]], 'float32')
            ST = array([[sx, 0., 0., tx],
                       [0., sy, 0., ty],
                       [0., 0., sz, tz],
                       [0., 0., 0., 1.]], 'float32')
            transformation = dot(ST, R)


            cortexSimplexMesh.points = AppplyTransformation(cortexSimplexMesh.points.copy(), transformation)

        if desplegarDatos:
            print ' '
            print 'Ajustar corteza y con ella ventriculos'
            print 'tetax =', tetax*180./pi
            print 'tetay =', tetay*180./pi
            print 'tetaz =', tetaz*180./pi
            print 'sx =', sx
            print 'sy =', sy
            print 'sz =', sz
            print 'tx =', tx
            print 'ty =', ty
            print 'tz =', tz
            print 'Distancia Antes =', TransformationAdjustment(parametrosIni.copy(), cortexSimplexMesh.points.copy(), distance, paramatrosAjust).sum()
            print 'Distancia Despues =', TransformationAdjustment(parametrosOpt[0], cortexSimplexMesh.points.copy(), distance, paramatrosAjust).sum()


        SaveData(workFolder + dataNames['cortexRegisterFileName'] , ('simplexMesh',), (cortexSimplexMesh,))
        PtoVTK = PythonToPolyData()
        PtoVTK.SetInputSimplex(cortexSimplexMesh)
        PtoVTK.WritePolyData(workFolder+ dataNames['cortexRegisterFileName'] + '.vtk')

        if dibujos_flag:
            # ############### para Visualizar
            PtoVTK = PythonToPolyData()
            PtoVTK.SetInputSimplex(cortexSimplexMesh)
            PtoVTK.WritePolyData('..\\Borrar6.vtk')
            PtoVTK.WritePolyData(workFolder+ dataNames['cortexRegisterFileName'] + '.vtk')
            # ###############     
        
    if segmentation_Flags['flag_segStep_deform1']: # deformacion de corteza usando segmentacion por pixeles

        dibujos_Flag = 1      
        NameFileImagen_aux = workFolder + prePegmentedImageFile
        NameFileMeshCortex = workFolder + dataNames['cortexRegisterFileName']
        flag_desplegarDatos = 0
        
        # leer malla corteza
        tupla=ReadData(NameFileMeshCortex)
        script=''
        for i in range(len(tupla[0])):
            exec(script.join((tupla[0][i],'=tupla[1][i]')))
        if tupla[0] == ('trianglesMesh',):
            cortexSimplexMesh = trianglesMesh.TrianglesToSimplex4()
        else:
            cortexSimplexMesh = simplexMesh
        del tupla


        # leer imagen segmentada de cerebro
        hdr = medimages.read_ANALYZE_tags(NameFileImagen_aux)
        imagen = medimages.read_ANALYZE(hdr)

        imagen_bin = imagen>0
        imagen_bin = imagen_bin.astype('float32')
        if 0:
            # closing (optativo)
            radio = 3 # 3
            iteraciones = 2 #1
            SE = zeros(((radio*2)+1,(radio*2)+1,(radio*2)+1))
            for i in range(-radio,radio+1):
                for j in range(-radio,radio+1):
                    for z in range(-radio,radio+1):
                        if sqrt(i**2+j**2+z**2)<=radio:
                            SE[i+radio,j+radio,z+radio] = 1
            imagen_binClose = ndimageMorphology.binary_closing(imagen_bin, structure = SE, iterations = iteraciones, output = None, origin = 0)    
            imagen_bin = imagen_binClose
            del imagen_binClose
            imagen_bin = imagen_bin.astype('float32')
            # ##

        imagen_bin *= -1

        # sobel
        a = array([[1, 3, 1], [3, 6, 3], [1, 3, 1]])
        dx = zeros((3,3,3), 'float32'); dy = zeros((3,3,3), 'float32'); dz = zeros((3,3,3), 'float32')
        dy[0,:,:]=a; dy[2,:,:]=-a
        dx[:,0,:]=a; dx[:,2,:]=-a;
        dz[:,:,0]=a; dz[:,:,2]=-a;

        MBx = ndimageFilters.convolve(imagen_bin, dx, output = None, mode = 'reflect', cval = 0.0, origin = 0)
        MBy = ndimageFilters.convolve(imagen_bin, dy, output = None, mode = 'reflect', cval = 0.0, origin = 0)
        MBz = ndimageFilters.convolve(imagen_bin, dz, output = None, mode = 'reflect', cval = 0.0, origin = 0)

        forces = {'dx' : MBx, 'dy' : MBy, 'dz' : MBz, 'imagen_bin' :imagen_bin} 


        ITER = parameters['def1_MaxIter'] #100#250
##        thresholdITER = 0.0001 #umbral de desplazamiento promedio en la malla
        thresholdITER = 0.3 #umbral de desplazamiento promedio en la malla
        S = 3 #tamaño de la vecindad tomada en cuenta en el angulo objetivo
##        angleSimplexObjetivo = -1. # si < 0 se toma en cuenta la vecindad y se ponderan los vecinos
        cortexSimplexMesh.ComputeAngleSimplex()
        angleSimplexObjetivo = cortexSimplexMesh.angleSimplex.copy()
        recomputeParameters = 100  #cada cuantas iteraciones se recalculan los parametros de acomodacion
        sigma= 0.1  #parametro para la acomodacion de la malla
        gamma= parameters['def1_gamma'] #0.65  #biscoso
        iterCorrection = 20 # cada cuantas iteraciones se corrigen problemas de intersecciones 

        alfa = parameters['def1_nu'] # 0.8 #0.2 ponderacion fuerza interna
        kappa = parameters['def1_lambda'] # 0.3 #0.15 #ponderacion fuerza externa
        beta = 0#; presion (-infla)

        Scontour = 2 # si < 0 el angulo simplex objetivo de cada punto del contorno sera 0
        AngleGammaObjetivo = 0 # para contorno
        modo = '1erDef'

        pointsCortex_nuevo = SimplexDeform3DASM(cortexSimplexMesh, alfa, beta, gamma, kappa, forces, ITER, S, sigma, recomputeParameters, thresholdITER, Scontour, AngleGammaObjetivo, angleSimplexObjetivo, modo = modo, extraData = None, iterCorrection = iterCorrection)

        cortexSimplexMesh.points = pointsCortex_nuevo

        SaveData(workFolder + dataNames['cortexDef1'], ('simplexMesh',), (cortexSimplexMesh,))
        PtoVTK = PythonToPolyData()
        PtoVTK.SetInputSimplex(cortexSimplexMesh)
        PtoVTK.WritePolyData(workFolder+ dataNames['cortexDef1'] + '.vtk')

        
        if dibujos_Flag:
            # ############### para Visualizar
            PtoVTK = PythonToPolyData()
            PtoVTK.SetInputSimplex(cortexSimplexMesh)
            PtoVTK.WritePolyData('..\\Borrar7.vtk')
            PtoVTK.WritePolyData(workFolder+ dataNames['cortexDef1'] + '.vtk')
            # ###############

        
        SaveData(workFolder + 'datos_aux', ('angleSimplexObjetivo',), (angleSimplexObjetivo,))

    if segmentation_Flags['flag_segStep_deform2']: # deformacion de corteza con imagen entera
        dibujos_Flag = 0
        NameFileImagen_aux = imageFolder + imageFile
        NameFileMeshCortex = workFolder + dataNames['cortexDef1']
        flag_desplegarDatos = 0
        
        RESHAPE111 = 1 # Para crear imagenes con pixeles de 1x1x1
        
        # leer malla corteza
        tupla=ReadData(NameFileMeshCortex)
        script=''
        for i in range(len(tupla[0])):
            exec(script.join((tupla[0][i],'=tupla[1][i]')))
        if tupla[0] == ('trianglesMesh',):
            cortexSimplexMesh = trianglesMesh.TrianglesToSimplex4()
        else:
            cortexSimplexMesh = simplexMesh
        del tupla

        tupla = ReadData(workFolder + dataNames['dataPreSeg']) #'Totsu' 'u_nivel_grisOrig','s_nivel_grisOrig'
        script=''
        for i in range(len(tupla[0])):
            exec(script.join((tupla[0][i],'=tupla[1][i]')))
        del tupla

        tupla = ReadData(workFolder + 'datos_aux')
        script=''
        for i in range(len(tupla[0])):
            exec(script.join((tupla[0][i],'=tupla[1][i]')))
        if flag_desplegarDatos:
            print tupla[0]
        del tupla

        if 0: # Aumentar resolucion
            # ###################### Prueba # ######################
            cortexTrianglesMesh = cortexSimplexMesh.SimplexToTriangles4()
            PtoVTK = PythonToPolyData()
            PtoVTK.SetInputTriangles(cortexTrianglesMesh)
            VTK_CortexTrianglesMesh = PtoVTK.GetVTKPolyData()
            VTK_Butterfly  = vtkButterflySubdivisionFilter()
            VTK_Butterfly.SetInputConnection(VTK_CortexTrianglesMesh['vtkData'].GetProducerPort())
            VTK_Butterfly.Update()
            VTK_CortexTrianglesMesh_sub = VTK_Butterfly.GetOutput()
            VTKtoP = PolyDataToPython()
            VTKtoP.SetInput(VTK_CortexTrianglesMesh_sub)
            cortexTrianglesMesh_sub = VTKtoP.GetOutput()
            cortexSimplexMesh = cortexTrianglesMesh_sub.TrianglesToSimplex4()
            # ###################### Prueba # ######################
        

        if dibujos_Flag:
            # ############### para Visualizar
            PtoVTK = PythonToPolyData()
            PtoVTK.SetInputSimplex(cortexSimplexMesh)
            PtoVTK.WritePolyData('..\\Borrar6.vtk')
            # ###############

        # leer imagen
        hdr = medimages.read_ANALYZE_tags(NameFileImagen_aux)
        imagen = medimages.read_ANALYZE(hdr)

        if RESHAPE111 and (hdr['Spacing'] != [1,1,1]):
            # Cambiar tamano de imagen
            transformation = eye(3)
            transformation[0,0] =  1./hdr['Spacing'][0]
            transformation[1,1] =  1./hdr['Spacing'][1]
            transformation[2,2] =  1./hdr['Spacing'][2]
            output_shape_aux = array(imagen.shape) * hdr['Spacing']
            imagen = affine_transform(imagen, transformation, offset = 0.0, output_shape = output_shape_aux, output_type = None, output = None, order = 1, mode = 'constant', cval = 0.0, prefilter = False)


        # sobel
        imagenF = contrast3D(imagen.copy(), 0., 1.)
        a = array([[1, 3, 1], [3, 6, 3], [1, 3, 1]])
        dx = zeros((3,3,3), 'float32'); dy = zeros((3,3,3), 'float32'); dz = zeros((3,3,3), 'float32')
        dy[0,:,:]=a; dy[2,:,:]=-a
        dx[:,0,:]=a; dx[:,2,:]=-a;
        dz[:,:,0]=a; dz[:,:,2]=-a;

        MBx = ndimageFilters.convolve(imagenF, dx, output = None, mode = 'reflect', cval = 0.0, origin = 0)
        MBy = ndimageFilters.convolve(imagenF, dy, output = None, mode = 'reflect', cval = 0.0, origin = 0)
        MBz = ndimageFilters.convolve(imagenF, dz, output = None, mode = 'reflect', cval = 0.0, origin = 0)
        del imagenF

        forces = {'dx' : MBx, 'dy' : MBy, 'dz' : MBz, 'imagen' :imagen}
        
        ITER = parameters['def2_MaxIter']#300#200
##        thresholdITER = 0.03 #umbral de desplazamiento promedio en la malla
        thresholdITER = 0.01 #umbral de desplazamiento promedio en la malla
        S = 2 #3 amaño de la vecindad tomada en cuenta en el angulo objetivo
        try:
            angleSimplexObjetivo
        except:
            angleSimplexObjetivo = -1. # si < 0 se toma en cuenta la vecindad y se ponderan los vecinos
        angleSimplexObjetivo = -1.
        recomputeParameters = 100  #cada cuantas iteraciones se recalculan los parametros de acomodacion
        sigma = 0.1  #parametro para la acomodacion de la malla
        gamma= parameters['def2_gamma'] #0.1(final) #0.8 #biscoso
        iterCorrection = 10 # cada cuantas iteraciones se corrigen problemas de intersecciones
        
        alfa = parameters['def2_nu'] #0.5 #0.5 # 0.5 (10-5-2011) # 0.7 # 0.5 #0.7 ponderacion fuerza interna
        kappa = parameters['def2_lambda'] #0.7 # 0.7 # 0.4  (10-5-2011)  # 0.3 # 0.5 #0.3 #ponderacion fuerza externa
##        alfa = 0.3 #0.2 ponderacion fuerza interna
##        kappa = 0.5 #0.15 #ponderacion fuerza externa        
        beta=0 #; presion (-infla)

        Scontour = parameters['def2_S'] #2 # si < 0 el angulo simplex objetivo de cada punto del contorno sera 0
        AngleGammaObjetivo = 0 # para contorno
        modo = '2aDef2'

        pointsCortex_nuevo = SimplexDeform3DASM(cortexSimplexMesh, alfa, beta, gamma, kappa, forces, ITER, S, sigma, recomputeParameters, thresholdITER, Scontour, AngleGammaObjetivo, angleSimplexObjetivo, modo = modo, extraData = skullStrippingData, iterCorrection = iterCorrection)


        if RESHAPE111:
            # ############### para Visualizar
            PtoVTK = PythonToPolyData()
            PtoVTK.SetInputSimplex(cortexSimplexMesh)
            PtoVTK.WritePolyData(workFolder + dataNames['cortexDef2'] + '.vtk')
            # ###############            
            # Cambiar tamano de imagen
            pointsCortex_nuevo = (pointsCortex_nuevo / array(hdr['Spacing'])).astype('float32')

        cortexSimplexMesh.points = pointsCortex_nuevo

        SaveData(workFolder + dataNames['cortexDef2'] , ('simplexMesh',), (cortexSimplexMesh,))
        if dibujos_Flag:
            # ############### para Visualizar
            PtoVTK = PythonToPolyData()
            PtoVTK.SetInputSimplex(cortexSimplexMesh)
            PtoVTK.WritePolyData('..\\Borrar7.vtk')
            PtoVTK.WritePolyData(workFolder + dataNames['cortexDef2'] + '.vtk')
            # ###############

    if segmentation_Flags['flag_segStep_deform3']: # deformacion de corteza con imagen entera
        dibujos_Flag = 0
        flag_desplegarDatos = 0
        NameFileImagen_aux = dataNames['imageFolder'] + dataNames['imageFile']
        NameFileMeshCortex = dataNames['workFolder'] + dataNames['cortexDef2']
        
        RESHAPE111 = 1 # Para crear imagenes con pixeles de 1x1x1
        
        # leer malla corteza
        tupla=ReadData(NameFileMeshCortex)
        script=''
        for i in range(len(tupla[0])):
            exec(script.join((tupla[0][i],'=tupla[1][i]')))
        if tupla[0] == ('trianglesMesh',):
            cortexSimplexMesh = trianglesMesh.TrianglesToSimplex4()
        else:
            cortexSimplexMesh = simplexMesh
        del tupla

        # leer imagen
        hdr = medimages.read_ANALYZE_tags(NameFileImagen_aux)
        imagen = medimages.read_ANALYZE(hdr)

        # ############ solo por PRUEBA ########################
        if RESHAPE111:
            cortexSimplexMesh.points  = (cortexSimplexMesh.points  * array(hdr['Spacing'])).astype('float32')
        # ########################################################        

        if 1: # Aumentar resolucion
            # ###################### Prueba # ######################
            cortexTrianglesMesh = cortexSimplexMesh.SimplexToTriangles4()
            PtoVTK = PythonToPolyData()
            PtoVTK.SetInputTriangles(cortexTrianglesMesh)
            VTK_CortexTrianglesMesh = PtoVTK.GetVTKPolyData()
            VTK_Butterfly  = vtkButterflySubdivisionFilter()
            VTK_Butterfly.SetInputConnection(VTK_CortexTrianglesMesh['vtkData'].GetProducerPort())
            VTK_Butterfly.Update()
            VTK_CortexTrianglesMesh_sub = VTK_Butterfly.GetOutput()
            VTKtoP = PolyDataToPython()
            VTKtoP.SetInput(VTK_CortexTrianglesMesh_sub)
            cortexTrianglesMesh_sub = VTKtoP.GetOutput()
            cortexSimplexMesh = cortexTrianglesMesh_sub.TrianglesToSimplex4()
            # ###################### Prueba # ######################
        
        if dibujos_Flag:
            # ############### para Visualizar
            PtoVTK = PythonToPolyData()
            PtoVTK.SetInputSimplex(cortexSimplexMesh)
            PtoVTK.WritePolyData('..\\Borrar6.vtk')
            # ###############

        tupla = ReadData(workFolder + dataNames['dataPreSeg']) #'Totsu' 'u_nivel_grisOrig','s_nivel_grisOrig'
        script=''
        for i in range(len(tupla[0])):
            exec(script.join((tupla[0][i],'=tupla[1][i]')))
        del tupla

        tupla = ReadData(workFolder + 'datos_aux')
        script=''
        for i in range(len(tupla[0])):
            exec(script.join((tupla[0][i],'=tupla[1][i]')))
        if flag_desplegarDatos:
            print tupla[0]
        del tupla


        if RESHAPE111 and (hdr['Spacing'] != [1,1,1]):
            # Cambiar tamano de imagen
            transformation = eye(3)
            transformation[0,0] =  1./hdr['Spacing'][0]
            transformation[1,1] =  1./hdr['Spacing'][1]
            transformation[2,2] =  1./hdr['Spacing'][2]
            output_shape_aux = array(imagen.shape) * hdr['Spacing']
            imagen = affine_transform(imagen, transformation, offset = 0.0, output_shape = output_shape_aux, output_type = None, output = None, order = 1, mode = 'constant', cval = 0.0, prefilter = False)


        # sobel
        imagenF = contrast3D(imagen.copy(), 0., 1.)
        a = array([[1, 3, 1], [3, 6, 3], [1, 3, 1]])
        dx = zeros((3,3,3), 'float32'); dy = zeros((3,3,3), 'float32'); dz = zeros((3,3,3), 'float32')
        dy[0,:,:]=a; dy[2,:,:]=-a
        dx[:,0,:]=a; dx[:,2,:]=-a;
        dz[:,:,0]=a; dz[:,:,2]=-a;

        MBx = ndimageFilters.convolve(imagenF, dx, output = None, mode = 'reflect', cval = 0.0, origin = 0)
        MBy = ndimageFilters.convolve(imagenF, dy, output = None, mode = 'reflect', cval = 0.0, origin = 0)
        MBz = ndimageFilters.convolve(imagenF, dz, output = None, mode = 'reflect', cval = 0.0, origin = 0)
        del imagenF

        forces = {'dx' : MBx, 'dy' : MBy, 'dz' : MBz, 'imagen' :imagen}

        ITER = parameters['def3_MaxIter'] #100 
        thresholdITER = 0.01 #umbral de desplazamiento promedio en la malla
        S = 2 #3 amaño de la vecindad tomada en cuenta en el angulo objetivo
        try:
            angleSimplexObjetivo
        except:
            angleSimplexObjetivo = -1. # si < 0 se toma en cuenta la vecindad y se ponderan los vecinos
        angleSimplexObjetivo = -1.
        recomputeParameters = 100  #cada cuantas iteraciones se recalculan los parametros de acomodacion
        sigma = 0.1  #parametro para la acomodacion de la malla
        gamma= parameters['def3_gamma'] #0.1 #0.8 #biscoso
        iterCorrection = 10 # cada cuantas iteraciones se corrigen problemas de intersecciones 

        alfa = parameters['def3_nu'] # 0.4 # 0.5 # 0.7 # 0.5 #0.7 ponderacion fuerza interna
        kappa = parameters['def3_lambda'] #0.32 # 0.4 # 0.3 # 0.5 #0.3 #ponderacion fuerza externa
##        alfa = 0.3 #0.2 ponderacion fuerza interna
##        kappa = 0.5 #0.15 #ponderacion fuerza externa        
        beta=0 #; presion (-infla)

        Scontour = parameters['def3_S'] #2 # si < 0 el angulo simplex objetivo de cada punto del contorno sera 0
        AngleGammaObjetivo = 0 # para contorno
        modo = '2aDef3'

        pointsCortex_nuevo = SimplexDeform3DASM(cortexSimplexMesh, alfa, beta, gamma, kappa, forces, ITER, S, sigma, recomputeParameters, thresholdITER, Scontour, AngleGammaObjetivo, angleSimplexObjetivo, modo, extraData = skullStrippingData, iterCorrection = iterCorrection)


        if RESHAPE111:
            # ############### para Visualizar
            PtoVTK = PythonToPolyData()
            PtoVTK.SetInputSimplex(cortexSimplexMesh)
            PtoVTK.WritePolyData(workFolder + dataNames['cortexDef3'] + '.vtk')
            # ###############            
            # Cambiar tamano de imagen
            pointsCortex_nuevo = (pointsCortex_nuevo / array(hdr['Spacing'])).astype('float32') # a tamano de imagen original

        cortexSimplexMesh.points = pointsCortex_nuevo

        SaveData(dataNames['workFolder'] + dataNames['cortexDef3'] , ('simplexMesh',), (cortexSimplexMesh,))
        
        if dibujos_Flag:
            # ############### para Visualizar
            PtoVTK = PythonToPolyData()
            PtoVTK.SetInputSimplex(cortexSimplexMesh)
            PtoVTK.WritePolyData('..\\Borrar7.vtk')
            PtoVTK.WritePolyData(dataNames['workFolder'] + dataNames['mesh3deformedFilename'] + '.vtk')
            # ###############
        

    if segmentation_Flags['flag_createBinImage']: # crear la imagen binaria
        dibujos_Flag = 0
        # Datos de la imagen que se segmento
        filePathIn = dataNames['imageFolder'] + dataNames['imageFile']
        hdrOrig = medimages.read_ANALYZE_tags(filePathIn) # esta debe tener las dimenciones correctas del sitio web
        fh = open(filePathIn[:-4]+'.hdr','r')
        ss = fh.read()
        fh.close()
        spacing_original = array(hdrOrig['Spacing'], 'float32')
        
        NameFileMeshCortex = dataNames['workFolder'] + dataNames['cortexDef3']
        
        # leer malla corteza
        tupla=ReadData(NameFileMeshCortex)
        script=''
        for i in range(len(tupla[0])):
            exec(script.join((tupla[0][i],'=tupla[1][i]')))
        if tupla[0] == ('simplexMesh',):
            if segmentation_Options['flag_SimplexTrianMethodBinImag'] == 'tangentesDual':
                trianglesMesh = simplexMesh.SimplexToTriangles4()
            elif segmentation_Options['flag_SimplexTrianMethodBinImag'] == 'tangentesNotDual':
                trianglesMesh = simplexMesh.SimplexToTrianglesNotDualConservacion()
        elif tupla[0] == ('trianglesMesh',):
            None
        else:
            print 'La malla no se cargo completamente'

        trianglesMesh.ComputeNormalsOfTriangles()

        # Datos para el muestreo ###
        origen_muestreo = [0,0,0]
        fin_muestreo = [hdrOrig['Slices'], hdrOrig['Rows'], hdrOrig['Columns']]
        spacing_muestreo = [1.,1.,1.] #hdr_image['Spacing'][:]

        [pointsBin, imagenBin]  = MakeImageBinAndPointsFromMesh(trianglesMesh, spacing_muestreo, array([origen_muestreo, fin_muestreo]))

        if dibujos_Flag:
            # ############### para Visualizar
            PtoVTK = PythonToPolyData()
            PtoVTK.SetInputPoints(pointsBin)
            PtoVTK.WritePolyData('..\\points.vtk')
            # ##############
        del pointsBin

        # ########### Guardar imagen Binaria ##########
        if dibujos_Flag:
            fileNameImageBinVTK = dataNames['binSegImage'] + '_malla' + '.vtk'
            # en VTK
            SaveImage3DasVTK(imagenBin, workFolder + fileNameImageBinVTK, spacing_muestreo)
        fileNameImageBinAnalize = dataNames['binSegImage'] + '.img'
        folderImageAnalize = dataNames['workFolder']

        # en Analize
##        medimages.save_ANALYZE(imagenBin,fileNameImageBinAnalize, spacing=spacing_muestreo, folder = folderImageAnalize)
        medimages.save_ANALYZE(imagenBin,fileNameImageBinAnalize, spacing = spacing_original, folder = folderImageAnalize)        
        # #############################################
        ss_new = ['\x00']*348 # datos para encabezado nuevo
        ss_new[:] = ss[:] # copiar los datos del encabezado de la imagen original
        echar='<'
        datatype = imagenBin.dtype # formato numerico de la imagen
        if datatype=='int16':
            dt = 4
            bitpix = 16
        elif datatype=='int8':
            dt = 2
            bitpix = 8
        elif datatype=='float32':
            dt = 16
            bitpix = 32
        import struct
        ss_new[72:74] = struct.pack(echar+'H',bitpix) # correccion del formato numerico de la imagen, las mascaras son binarias entonces se pueden guardar en 'int8'
        ss_new[70:72] = struct.pack(echar+'H',dt)
        ss_new = "".join(ss_new)
        fh = open((folderImageAnalize + fileNameImageBinAnalize)[:-4]+'.hdr','wb') # escribir encabezado de la imagen binaria con datos correctos
        fh.write(ss_new)
        fh.close()
        
        del imagenBin
        
    if segmentation_Flags['flag_segStep_CondDilatation']:
        
        RESHAPE111 = 1 # Para crear imagenes con pixeles de 1x1x1
        SE = ones((3,3,3))
        
        # cargar la imagen original ######
        hdrOrig = medimages.read_ANALYZE_tags(dataNames['imageFolder'] + dataNames['imageFile'])
        imagenOrig = medimages.read_ANALYZE(hdrOrig)
        shapeOrig = imagenOrig.shape
        # ################
        
        # Cargar imagen binaria de la segmentacion por simplex
        hdrBin = medimages.read_ANALYZE_tags(dataNames['workFolder'] + dataNames['binSegImage'] + '.img')
        imagenBin = medimages.read_ANALYZE(hdrBin)
        # ################

        # Cargar datos de la pre-segmentacion
        tupla = ReadData(dataNames['workFolder'] + dataNames['dataPreSeg']) #'Totsu' 'u_nivel_grisOrig','s_nivel_grisOrig'
        script=''
        for i in range(len(tupla[0])):
            exec(script.join((tupla[0][i],'=tupla[1][i]')))
        del tupla        


        if RESHAPE111 and (hdrOrig['Spacing'] != [1,1,1]):
            # Cambiar tamano de imagen
            transformation = eye(3)
            transformation[0,0] =  1./hdrOrig['Spacing'][0]
            transformation[1,1] =  1./hdrOrig['Spacing'][1]
            transformation[2,2] =  1./hdrOrig['Spacing'][2]
            output_shape_aux = array(imagenOrig.shape) * hdrOrig['Spacing']
            imagenOrig = affine_transform(imagenOrig, transformation, offset = 0.0, output_shape = output_shape_aux, output_type = None, output = None, order = 3, mode = 'constant', cval = 0.0, prefilter = False)        
            imagenBin = affine_transform(imagenBin, transformation, offset = 0.0, output_shape = output_shape_aux, output_type = None, output = None, order = 0, mode = 'constant', cval = 0.0, prefilter = False)
            
        imagenBin = imagenBin > 0.5

        # ##################### PROCESAMIENTO ##################################

        # ##### Eliminacion de voxeles de borde
        for iter_Ero in range(2):
            # erocion de imagen binaria
            imagenBinE = imagenBin.copy()
            ndimageMorphology.binary_erosion(imagenBin.copy(), structure = SE, iterations = 1, mask = None, output = imagenBinE, border_value = 0, origin = 0, brute_force = False)
            # dejar solo la parte erocionada
            imagenBinE =  imagenBin - imagenBinE

            imagenOrigMaskD = imagenOrig * imagenBinE
            # ver que voxeles se erodan
            umbral = (skullStrippingData['u_nivel_grisOrig'][1] - 8 * skullStrippingData['s_nivel_grisOrig'][1])
            if umbral <= 0:
                print 'umbral menor que 0'
                umbral = 0.001
            imagenBinE_aux = (imagenOrigMaskD >= umbral) # permanecen
            imagenBinE = imagenBinE - imagenBinE_aux # los que se deben eliminar

            # nueva mascara
            imagenBinE = imagenBin - imagenBinE
            
            imagenBin = imagenBinE.copy()

        

        # ##### Dilatacion de voxeles de borde
        for iter_Dil in range(1):
            # enmascara imagen Original
            imagenOrigMask = imagenOrig * imagenBin

            # aplicar filtro maximo a imagen original enmascarada
            imagenOrigMax = ndimageFilters.maximum_filter(imagenOrigMask, size = None, footprint = SE, output = None, mode = 'constant', cval = 0.0, origin = 0)
            del imagenOrigMask
            
            # dilatar imagen binaria
            iterationsD = 1
            imagenBinD = imagenBin.copy()
            ndimageMorphology.binary_dilation(imagenBin.copy(), structure = SE, iterations = iterationsD, mask = None, output = imagenBinD, border_value = 0, origin = 0, brute_force = False)
            # dejar solo la parte dilatada
            imagenBinD =  (imagenBinD.astype('uint8') + imagenBin.astype('uint8'))
            imagenBinD = (imagenBinD > 0.5) * (imagenBinD < 1.5)

            imagenOrigMaskD = imagenOrig * imagenBinD
            
            # ver que voxeles se dilataran
        ##    imagenBinD = (imagenOrigMaskD >= (imagenOrigMax - 5 * skullStrippingData['s_nivel_grisOrig'][1])) * (imagenOrigMax >= (skullStrippingData['u_nivel_grisOrig'][1] - 8 * skullStrippingData['s_nivel_grisOrig'][1])) * imagenBinD
            imagenBinD = (imagenOrigMaskD >= (imagenOrigMax - 5 * skullStrippingData['s_nivel_grisOrig'][1])) * (imagenOrigMaskD >= (skullStrippingData['u_nivel_grisOrig'][1] - 8 * skullStrippingData['s_nivel_grisOrig'][1])) * imagenBinD
            
        ##    imagenBinD = (imagenOrigMaskD >= 0) * imagenBinD

            # nueva mascara
            imagenBinD = (imagenBin + imagenBinD) > 0.5            

            imagenBin = imagenBinD.copy() # para las iteraciones de dilatacion



        # #### GUARDAR LA NUEVA SEGMENTACION #######

##        ndimageMorphology.binary_closing(imagenBin.copy(), structure = SE, iterations = 1, output = imagenBin, origin = 0)

        if RESHAPE111 and (hdrOrig['Spacing'] != [1,1,1]):
            # Cambiar tamano de imagen
            transformation = eye(3)
            transformation[0,0] =  float(hdrOrig['Spacing'][0])
            transformation[1,1] =  float(hdrOrig['Spacing'][1])
            transformation[2,2] =  float(hdrOrig['Spacing'][2])
            output_shape_aux = shapeOrig
            imagenBin = affine_transform(imagenBin.astype('uint8'), transformation, offset = 0.0, output_shape = output_shape_aux, output_type = None, output = None, order = 3, mode = 'constant', cval = 0.0, prefilter = False)
        
        imagenBin = imagenBin > 0.5
        imagenBin = imagenBin.astype('uint8')

        # ######## Leer encabezado de imagen original ######
        fh = open(dataNames['imageFolder'] + dataNames['imageFile'][:-4]+'.hdr','r')
        ss = fh.read()
        fh.close()    
        
        spacing_muestreo_guardar = array(hdrOrig['Spacing'], 'float32')
       
        # en Analize
        res = medimages.save_ANALYZE(imagenBin, dataNames['binSegImage'] + '.img', spacing=spacing_muestreo_guardar, folder = dataNames['workFolder'] )

        ss_new = ['\x00']*348 # datos para encabezado nuevo
        ss_new[:] = ss[:] # copiar los datos del encabezado de la imagen original
        echar='<'
        datatype = imagenBin.dtype # formato numerico de la imagen
        if datatype=='int16':
            dt = 4
            bitpix = 16
        elif datatype=='int8':
            dt = 2
            bitpix = 8
        elif datatype=='uint8':
            dt = 2
            bitpix = 8        
        elif datatype=='float32':
            dt = 16
            bitpix = 32
            
        import struct
        ss_new[72:74] = struct.pack(echar+'H',bitpix) # correccion del formato numerico de la imagen, las mascaras son binarias entonces se pueden guardar en 'int8'
        ss_new[70:72] = struct.pack(echar+'H',dt)
        ss_new = "".join(ss_new)
        fh = open(dataNames['workFolder'] +  dataNames['binSegImage'] + '.hdr','wb') # escribir encabezado de la imagen binaria con datos correctos
        fh.write(ss_new)
        fh.close()

        del imagenBin
        del imagenOrig

    if segmentation_Flags['flag_masking']:
        # cargar la imagen original ######
        hdrOrig = medimages.read_ANALYZE_tags(dataNames['imageFolder'] + dataNames['imageFile'])
        imagenOrig = medimages.read_ANALYZE(hdrOrig)
        # ################
        
        # cargar imagen binaria
        hdrBin = medimages.read_ANALYZE_tags(dataNames['workFolder'] + dataNames['binSegImage'] + '.img')
        imagenBin = medimages.read_ANALYZE(hdrBin)
        # ################        # en Analize

        # aplicar la mascara
        imagenOrig = imagenOrig * imagenBin
        del imagenBin

        # en Analize        
        res = medimages.save_ANALYZE(imagenOrig, dataNames['SegImage'] + '.img', spacing = hdrOrig['Spacing'], folder = dataNames['workFolder'] )

        # ######## Leer encabezado de imagen original ######
        fh = open(dataNames['imageFolder'] + dataNames['imageFile'][:-4]+'.hdr','r')
        ss = fh.read()
        fh.close()
        
        ss_new = ['\x00']*348 # datos para encabezado nuevo
        ss_new[:] = ss[:] # copiar los datos del encabezado de la imagen original
        echar='<'
        datatype = imagenOrig.dtype # formato numerico de la imagen
        if datatype=='int16':
            dt = 4
            bitpix = 16
        elif datatype=='int8':
            dt = 2
            bitpix = 8
        elif datatype=='uint8':
            dt = 2
            bitpix = 8        
        elif datatype=='float32':
            dt = 16
            bitpix = 32
            
        import struct
        ss_new[72:74] = struct.pack(echar+'H',bitpix) # correccion del formato numerico de la imagen, las mascaras son binarias entonces se pueden guardar en 'int8'
        ss_new[70:72] = struct.pack(echar+'H',dt)
        ss_new = "".join(ss_new)
        fh = open(dataNames['workFolder'] +  dataNames['SegImage'] + '.hdr','wb') # escribir encabezado de la imagen binaria con datos correctos
        fh.write(ss_new)
        fh.close()

        
print 'Funciones para "SegmentationSteps" en la parte de: eliminar el craneo'

def HistogramAdjustment(parametros, histograma, parametrosIni):
##        print parametros
    u=parametros[0:3] #punto medio de la función
    s=parametros[3:6] #desviación estandar
    Np=parametros[6:9] #altura en punto medio
    
    histograma_parametros = zeros(len(histograma), 'float32')
    for i in range(len(histograma)):
        histograma_parametros[i]=0
        for k in range(3):
            histograma_parametros[i] += Np[k]*exp((-1./2.)*(((i-u[k])/s[k])**2.))
    diferencia = (histograma_parametros-histograma) # el cuadrado lo calcula la funcion leastsq al minimizar

##    restriccion1 = parametrosIni[7]/10.
##    restriccion2 = abs(parametrosIni[1] - parametrosIni[2]) / 2.
##    restriccion3 = abs(parametrosIni[1] - parametrosIni[2]) / 3.
##    if (Np[0] + restriccion1 - Np[1]) > 0:
##        diferencia += abs((Np[0] + restriccion1 - Np[1]))
##    if (Np[2] + restriccion1 - Np[1]) > 0:
##        diferencia += abs((Np[2] + restriccion1 - Np[1]))
##    for i in [0,1,2]:
##        if Np[i] < restriccion1:
##            diferencia += abs(Np[i]) + restriccion1
##    if (u[0] + restriccion2 - u[1]) > 0.:
##        diferencia += (u[0] + restriccion2 - u[1]) * 1000
##    if (u[1] + restriccion2 - u[2]) > 0:
##        diferencia += (u[1] + restriccion2 - u[2]) * 1000
##    if (u[0] + restriccion2*2 - u[2]) > 0.:
##        diferencia += (u[0] + restriccion2*2 - u[2]) * 1000
##    for i in [0,1,2]:
##        if s[i] < restriccion3:
##            diferencia += (restriccion3 - s[i]) * 10000
##

    if 0: # dibujo
        print 'diferencia histogramas = ' + str(sqrt(diferencia*diferencia).sum())
        figure()
        plot(histograma_parametros,'r')
        hold(True)
        plot(histograma,'b')
        print u
        print s
        print Np


##    a = diferencia.sum()
##        a = diferencia * diferencia
##    print  a.sum()
    return diferencia



def EliminarCraneo(dataNames, extras=None):

    flag_printData = 0
##    filData = extras['filDataEstimationTs']
##    filData.write('Imagen' + imageNumber + '\n')
    
    imageFolder = dataNames['imageFolder']
    workFolder = dataNames['workFolder']
    imageFile = dataNames['imageFile']
    prePegmentedImageFile = dataNames['preSegmentedImageFile']

    DIBUJO_TESIS = 0 # Para crear imagenes para la tesis
    GRAFICOS_Hist = 0 # Para crear graficos de los histogramas
    RESHAPE111 = 1 # Para crear imagenes con pixeles de 1x1x1

    
    # ######## Direccion de archivos ####################
    filePathIn = imageFolder + imageFile            
    fileOut = prePegmentedImageFile        
    folderOut = workFolder
    # ###################################################
    
    min_image = 1.;
    max_image = 257.; #257
    nbins = 256.; #256


    FactTs = 0.7 # se modifica mas adelante si esta encendido flag_FactTsEstimation #0.6 (paper 4./5.)0.7  #Ts = Totsu + FactTs * (ugm-Totsu)
    flag_FactTsEstimation = 1
    RadioOpeningM1= 3 #2 (paper 3)
    IteracionesOpeningM1 = 2 #1(paper 1)
    IteracionesDilatationM1= 2 #(paper 3) antes puse 2

    FiltroHistograma = True

    RadioOpeningM2= 5 #2(Brain) (paper 5)
    IteracionesOpeningM2=1


    # ######## Leer imagen ##########################
    hdr = medimages.read_ANALYZE_tags(filePathIn)
    imagen = medimages.read_ANALYZE(hdr)
    # ###############################################
    if RESHAPE111 and (hdr['Spacing'] != [1,1,1]):
        # Cambiar tamano de imagen
        transformation = eye(3)
        transformation[0,0] =  1./hdr['Spacing'][2]
        transformation[1,1] =  1./hdr['Spacing'][1]
        transformation[2,2] =  1./hdr['Spacing'][0]
        output_shape_aux = array(imagen.shape) * [hdr['Spacing'][2], hdr['Spacing'][1], hdr['Spacing'][0]]
        imagen = affine_transform(imagen, transformation, offset = 0.0, output_shape = output_shape_aux, output_type = None, output = None, order = 3, mode = 'constant', cval = 0.0, prefilter = False)

    # ########## modelo para comparaciones para estimar segundo umbral #################
    if flag_FactTsEstimation:
        hdr = medimages.read_ANALYZE_tags(dataNames['bainBinModel'])
        maskModelo_aux = medimages.read_ANALYZE(hdr)
        shapeMask = array([imagen.shape,maskModelo_aux.shape]).max(0)
        maskModelo = zeros(shapeMask, 'int8')            
        maskModelo[0:maskModelo_aux.shape[0],0:maskModelo_aux.shape[1],0:maskModelo_aux.shape[2]] = maskModelo_aux
        del maskModelo_aux
        maskModelo_Zmax = 139
        cMMCarlota = [88, 71] # [y,x] centro de masa de z = 73:
        # si ocupo transfor5macion affin
        bbMaskModelo = (maskModelo>0).nonzero()
        bbMaskModelo = array([[bbMaskModelo[0].min(), bbMaskModelo[0].max()], [bbMaskModelo[1].min(), bbMaskModelo[1].max()], [bbMaskModelo[2].min(), bbMaskModelo[2].max()]])
        if 1: # transformacion affin
            proyection = maskModelo.sum(2)
            proyection_index = proyection.nonzero()
            proyection_profile = zeros(((bbMaskModelo[0,1] - bbMaskModelo[0,0])+1))
            for index_p in range(((bbMaskModelo[0,1] - bbMaskModelo[0,0])+1)):
                try:
                    proyection_profile[index_p] = proyection_index[1][proyection_index[0]== (index_p + bbMaskModelo[0,0])].max()
                except:
                    None
            limProfile = bbMaskModelo[1,1] - (bbMaskModelo[1,1] - bbMaskModelo[1,0])*0.2
            for index_p in arange(argmax(proyection_profile),-1,-1):
                if proyection_profile[index_p] < limProfile:
                    break
            z_eyeMaskModelo = index_p + bbMaskModelo[0,0]

    
    # #########################################################

    max_imagenIni = double(imagen.max())
    min_imagenIni = double(imagen.min())

    ##min_image = 1.;
    ##max_image = 100.;
    ##nbins = 256.;
    p =((max_image-min_image)/nbins)
    edges_hist = arange(min_image,max_image+(p/2), p)
    centros_hist = arange(min_image+(((max_image-min_image)/nbins)/2), max_image, ((max_image-min_image)/nbins))
    imagenOrig = imagen.copy()
    imagen = contrast3D(imagen, min_image, max_image)
    tamano = imagen.shape

    # ############## umbral de Otsu ###########################################

    histograma = ndimageMeasurements.histogram(imagen, min_image, max_image, int(nbins), labels=None, index=None)
    histograma[-1] += (imagen==max_image).sum() # porque no toma en cuenta <= del ultimo bin

    if GRAFICOS_Hist: # dibujar
        width = abs(centros_hist[0] - centros_hist[1])
        figure()
        bar(centros_hist - (width), histograma, width = width, bottom = 0, color = 'g' , yerr = None, xerr = None, ecolor = 'k', capsize = 3)
##                bar(centros_hist - (histograma), histograma/histograma.sum(), width = width, bottom = 0, color = 'g' , yerr = None, xerr = None, ecolor = 'k', capsize = 3)
##                plot(centros_hist,histograma1)
        title('Histograma inicial ' + imageFile)
        show()

    N=float32(imagen.size)
    #u=imagen.mean()
    u=sum(centros_hist*histograma)/N # promedio ponderado
        
    Nb=0.;
    No=N;
    ub=0.;
    uo=u;
    var_between_n=zeros(nbins-1, 'float32')
    for T in range(int(nbins)-1):
        Nt = histograma[T]
        Nb_a = Nb
        No_a = No
        Nb += Nt
        No -= Nt
        if Nb == 0.:
            ub = 0.
        else:
            ub = ((ub*Nb_a)+(Nt*centros_hist[T]))/Nb # aproximacion
        uo = ((uo*No_a)-(Nt*centros_hist[T]))/No
        var_between = Nb*No*((ub-uo)**2)
        var_between_n[T] = var_between

    T = ndimageMeasurements.maximum_position(var_between_n, labels = None, index = None)[0]

    Totsu = edges_hist[T+1]
    
##    alfa = (max_image-min_image)/(max_imagenIni-min_imagenIni)
##    beta = max_image-(alfa*max_imagenIni)    
##    Totsu = (15. * alfa) + beta
##    
    M0 = imagen>Totsu
    histograma1 = histograma.copy(); histograma1[:T+1]=0
    B1 = imagen * M0  #eliminar el fondo
    del M0

    # #### para pasar informacion de los niveles de gris al paso ########
    # #### de segmentracion con modelos deformables              ########
    alfa = (max_image-min_image)/(max_imagenIni-min_imagenIni)
    beta = max_image-(alfa*max_imagenIni)
    Totsu_grisOrig = (Totsu - beta)/alfa
    # ################################            
    
    if DIBUJO_TESIS:
        if RESHAPE111:
            res = medimages.save_ANALYZE(B1.astype('int16'),'Debug_TOtsu',spacing=[1,1,1], folder = folderOut)
        else:
            res = medimages.save_ANALYZE(B1.astype('int16'),'Debug_TOtsu1',spacing=hdr['Spacing'], folder = folderOut)    
    # ########################################################

    if GRAFICOS_Hist: # dibujar
        width = abs(centros_hist[0] - centros_hist[1])
        figure()
        bar(centros_hist - (width), histograma1, width = width, bottom = 0, color = 'g' , yerr = None, xerr = None, ecolor = 'k', capsize = 3)
##                bar(centros_hist - (width), histograma1/histograma1.sum(), width = width, bottom = 0, color = 'g' , yerr = None, xerr = None, ecolor = 'k', capsize = 3)
##                plot(centros_hist,histograma1)
        title('Histograma 1 ' + imageFile)                
        show()

    I = ndimageMeasurements.maximum_position(histograma1, labels = None, index = None)[0]

    ugm = centros_hist[I]
    
# ##########Estimacion de segundo umbral            
    if flag_FactTsEstimation: # con estimacion por medio de mascara
        for iter_Ts in range(10):
            FactTs = 0.0 + iter_Ts * 0.1
            
            Ts = Totsu + FactTs * (ugm-Totsu)

            M1 = B1>Ts
            
            # opening of M1
            radio = RadioOpeningM1 
            SE = zeros(((radio*2)+1,(radio*2)+1,(radio*2)+1))
            for i in range(-radio,radio+1):
                for j in range(-radio,radio+1):
                    for z in range(-radio,radio+1):
                        if sqrt(i**2+j**2+z**2)<=radio:
                            SE[i+radio,j+radio,z+radio] = 1

            ndimageMorphology.binary_opening(M1, structure = SE, iterations = IteracionesOpeningM1, output = M1, origin = 0)

            L = ndimageMeasurements.label(M1, structure = None, output = None)

            nL = ndimageMeasurements.histogram(L[0], -0.5, L[1]+0.5, L[1]+1, labels=None, index=None)
            nL[0] = -1 #para no tomar en cuenta los voxeles con 0, (fondo)

            I1 = ndimageMeasurements.maximum_position(nL, labels = None, index = None)[0]
            nL[I1] = -1
            if max(nL) > 0: # si hay mas de una estructura hay que verificar cual corresponde al cerebro
                I2 = ndimageMeasurements.maximum_position(nL, labels = None, index = None)[0]

                M_aux1 = L[0]==I1;
                M_aux2 = L[0]==I2;

                volumen1 = M_aux1.sum()
                volumen2 = M_aux2.sum()

                bb1 = (M_aux1>0).nonzero()
                bb1 = array([[bb1[0].min(), bb1[0].max()], [bb1[1].min(), bb1[1].max()], [bb1[2].min(), bb1[2].max()]])
                bb2 = (M_aux2>0).nonzero()
                bb2 = array([[bb2[0].min(), bb2[0].max()], [bb2[1].min(), bb2[1].max()], [bb2[2].min(), bb2[2].max()]])

                if (float(volumen1)/volumen2)>4 or (float(volumen1)/volumen2)<0.25: # son muy diferentes
                    M1[:] = M_aux1 # elijo la mayor
                elif (bb1[:,0]>=bb2[:,0]).all() and (bb1[:,1]<=bb2[:,1]).all(): # M_aux1 esta contenido en M_aux2
                    M1[:] = M_aux1
                elif (bb1[:,0]<=bb2[:,0]).all() and (bb1[:,1]>=bb2[:,1]).all(): # M_aux2 esta contenido en M_aux1
                    M1[:] = M_aux2
                elif (float(volumen1)/volumen2)<2 and (float(volumen1)/volumen2)>0.5: # no son muy diferentes en tamaño
                    CM1 = ndimageMeasurements.center_of_mass(M_aux1, labels = None, index = None)
                    CM2 = ndimageMeasurements.center_of_mass(M_aux2, labels = None, index = None)    
                    if CM1[0]>CM2[0]: # asumo que el cuello es el otro y que está mas abajo
                        M1[:] = M_aux1
                    else:
                        M1[:] = M_aux2
                    del CM1, CM2
                else:
                    if volumen1>volumen2:
                        M1[:] = M_aux1
                    else:
                        M1[:] = M_aux2
                
##                    if (bb1[:,0]>=bb2[:,0]).all() and (bb1[:,1]<=bb2[:,1]).all(): # M_aux1 esta contenido en M_aux2
##                        M1[:] = M_aux1
##                    elif (bb1[:,0]<=bb2[:,0]).all() and (bb1[:,1]>=bb2[:,1]).all(): # M_aux2 esta contenido en M_aux1
##                        M1[:] = M_aux2
##                    elif (float(volumen1)/volumen2)<2 and (float(volumen1)/volumen2)>0.5: # no son muy diferentes en tamaño
##                        CM1 = ndimageMeasurements.center_of_mass(M_aux1, labels = None, index = None)
##                        CM2 = ndimageMeasurements.center_of_mass(M_aux2, labels = None, index = None)    
##                        if CM1[0]>CM2[0]: # asumo que el cuello es el otro y que está mas abajo
##                            M1[:] = M_aux1
##                        else:
##                            M1[:] = M_aux2
##                        del CM1, CM2
##                    else:
##                        if volumen1>volumen2:
##                            M1[:] = M_aux1
##                        else:
##                            M1[:] = M_aux2

                del M_aux1,M_aux2, volumen1, volumen2
                
            else:
                M1[:] = L[0]==I1;

            #cambie las iteraciones (eran 3)
            ndimageMorphology.binary_dilation(M1.copy(), structure = SE, iterations = IteracionesDilatationM1, mask = None, output = M1, border_value = 0, origin = 0, brute_force = False)

            B2 = B1 * M1 #dejar solo el cerebro y algunos tejidos alrededor

            # comparacion con modelo
            bbM1 = M1.nonzero()
            if 0: # transformacion rigida
                M1_Zmax = bbM1[0].max()
                if (M1_Zmax-73) > 0:
                    cCarlota = M1[M1_Zmax-73:,:,:].sum(0)
                else:
                    cCarlota = M1[0:,:,:].sum(0)                    
                cCarlota = ndimageMeasurements.center_of_mass(cCarlota, labels = None, index = None)
                shift = array([M1_Zmax-maskModelo_Zmax, cCarlota[0]-cMMCarlota[0], cCarlota[1]-cMMCarlota[1]])                
                maskModelo_reg = zeros(maskModelo.shape, 'uint8')
                maskModelo_reg = ndimageInterpolation.shift(maskModelo, shift, output_type = None, output = None, order = 0, mode = 'constant', cval = 0.0, prefilter = False)

            else: # transformacion affin
                bbM1 = array([[bbM1[0].min(), bbM1[0].max()], [bbM1[1].min(), bbM1[1].max()], [bbM1[2].min(), bbM1[2].max()]])
                lMaskModeloY = float(bbMaskModelo[1,1] - bbMaskModelo[1,0])
                lMaskModeloX = float(bbMaskModelo[2,1] - bbMaskModelo[2,0])

                centX = (bbM1[2,0] + bbM1[2,1])/2

##                lcent = (bbM1[2,1] - bbM1[2,0])/16
##                lcent = (bbM1[2,1] - bbM1[2,0])/8
                lcent = (bbM1[2,1] - bbM1[2,0])/30 #26-05-2011
                proyection = M1[:,:,centX-lcent:centX+lcent].sum(2)
                
##                    proyection = M1.sum(2)
                proyection_index = proyection.nonzero()
                proyection_profile = zeros(((bbM1[0,1] - bbM1[0,0])+1))
                for index_p in range(((bbM1[0,1] - bbM1[0,0])+1)):
                    try:
                        proyection_profile[index_p] = proyection_index[1][proyection_index[0] == (index_p + bbM1[0,0])].max()
                    except:
                        None
##                    limProfile = bbM1[1,1] - (bbM1[1,1] - bbM1[1,0])*0.2
                max_profile = 0
                maxM1Y = 0
                limProfile = -1
                for index_p in arange(len(proyection_profile)-1,-1,-1):
                    val_profile = proyection_profile[index_p]
                    if max_profile < val_profile:
                        max_profile = val_profile
                        maxM1Y_profileIndex = index_p + bbM1[0,0]
                        limProfile = max_profile - (max_profile - bbM1[1,0])*0.2
                    if proyection_profile[index_p] < limProfile:
                        z_eyeM1 = index_p + bbM1[0,0]
##                        maxM1Y = max_profile
##                        minM1Y = M1.sum(2)[maxM1Y_profileIndex,:].nonzero()[0].min()
                        minM1Y = M1.sum(2)[z_eyeM1:,:].nonzero()[1].min() # 16-06-2011 (el minimo dorsal desde los ojos hacia arriba)
                        maxM1Y = M1.sum(2)[(z_eyeM1 + maxM1Y_profileIndex)/2.:,:].nonzero()[1].max() # 17-06-2011
                        break
                        
                index_aux = M1.sum(1)[z_eyeM1:,:].nonzero()[1]
                minM1X = index_aux.min() # 16-06-2011 (calculando el maximo lateral desde los ojos hacia arriba)
                maxM1X = index_aux.max() # 16-06-2011 (calculando el maximo lateral desde los ojos hacia arriba)
                del index_aux

                if maxM1Y > 0:
                    sy = lMaskModeloY / (maxM1Y - minM1Y)
                else:
                    sy = lMaskModeloY / (bbM1[1,1] - bbM1[1,0])
                    maxM1Y = bbM1[1,1]
##                sx = lMaskModeloX / (bbM1[2,1] - bbM1[2,0])
                sx = lMaskModeloX / (maxM1X - minM1X) # 16-06-2011
                sz = (bbMaskModelo[0,1] - z_eyeMaskModelo) / float(bbM1[0,1] - z_eyeM1)
                if sz > 1.25:
                    sz = 1.25
                elif sz < 0.75:
                    sz = 0.75
                
##                    sz = (sx + sy) / 2.
                matrix = array([[sz,0,0],[0,sy,0],[0,0,sx]])
##                offset = array([bbMaskModelo[0,1]-sz*bbM1[0,1],bbMaskModelo[1,1]-sy*bbM1[1,1],bbMaskModelo[2,1]-sx*bbM1[2,1]])
                offset = array([bbMaskModelo[0,1]-sz*bbM1[0,1],bbMaskModelo[1,1]-sy*maxM1Y,bbMaskModelo[2,1]-sx*maxM1X])
                maskModelo_reg = ndimageInterpolation.affine_transform(maskModelo, matrix, offset, output_shape = None, output_type = None, output = None, order = 0, mode = 'constant', cval = 0.0, prefilter = False)
        
            imageComp = (M1 * 2) + maskModelo_reg[:M1.shape[0],:M1.shape[1],:M1.shape[2]]
            TP = float((imageComp==3).sum())
            FP = float((imageComp==2).sum())
            FN = float((imageComp==1).sum())
            TN = float((imageComp==0).sum())
            jaccard = TP/(TP+FP+FN)
            jaccard2 = TP/(TP+FP)
            porcent1 = FP/(TP+FN)
            porcent2 = FP/(TP)

            # BORRAR                
            aux = float((B1>0).sum())                
            Ts_test = (aux - M1.sum()) / aux
            if flag_printData:
                print 'FactTs: ',FactTs
                print 'porcent1: ',porcent1
                print 'porcent2: ',porcent2
##            filData.write('Ts %.7f\n' %(Ts))
##            filData.write('FactTs %.7f\n' %(FactTs))                
##            filData.write('Ts_Test %.7f\n' %(Ts_test))
##            filData.write('Jaccard %.7f\n' %(jaccard))
##            filData.write('Jaccard_sin_FN %.7f\n' %(jaccard2))
##            filData.write('porcent1 %.7f\n' %(porcent1))
##            filData.write('porcent2 %.7f\n' %(porcent2))
            #
            
            if porcent1 <= 0.08:
                if flag_printData:
                    print 'nn'
                break
                
        if DIBUJO_TESIS:
            if RESHAPE111:
                res = medimages.save_ANALYZE(B2.astype('int16'),'Debug_T_stripping',spacing=[1,1,1], folder = folderOut)
                res = medimages.save_ANALYZE(imageComp.astype('int16'),'imageComp',spacing=[1,1,1], folder = folderOut)
            else:
                res = medimages.save_ANALYZE(B2.astype('int16'),'Debug_T_stripping',spacing=hdr['Spacing'], folder = folderOut)

        if 0: # dibujar
            figure()
            imshow(B2[:,:,100], interpolation='nearest', origin='lower')
            show()


        del M1, B1, maskModelo_reg, maskModelo, proyection_index, proyection, imageComp

    else: # sin estimacion
        Ts = Totsu + FactTs * (ugm-Totsu)
        M1 = B1>Ts

        # opening of M1
        radio = RadioOpeningM1 
        SE = zeros(((radio*2)+1,(radio*2)+1,(radio*2)+1))
        for i in range(-radio,radio+1):
            for j in range(-radio,radio+1):
                for z in range(-radio,radio+1):
                    if sqrt(i**2+j**2+z**2)<=radio:
                        SE[i+radio,j+radio,z+radio] = 1

        ndimageMorphology.binary_opening(M1, structure = SE, iterations = IteracionesOpeningM1, output = M1, origin = 0)

        L = ndimageMeasurements.label(M1, structure = None, output = None)

        nL = ndimageMeasurements.histogram(L[0], -0.5, L[1]+0.5, L[1]+1, labels=None, index=None)
        nL[0] = -1 #para no tomar en cuenta los voxeles con 0, (fondo)

        I1 = ndimageMeasurements.maximum_position(nL, labels = None, index = None)[0]
        nL[I1] = -1
        if max(nL) > 0: # si hay mas de una estructura hay que verificar cual corresponde al cerebro
            I2 = ndimageMeasurements.maximum_position(nL, labels = None, index = None)[0]

            M_aux1 = L[0]==I1;
            M_aux2 = L[0]==I2;

            volumen1 = M_aux1.sum()
            volumen2 = M_aux2.sum()

            bb1 = (M_aux1>0).nonzero()
            bb1 = array([[bb1[0].min(), bb1[0].max()], [bb1[1].min(), bb1[1].max()], [bb1[2].min(), bb1[2].max()]])
            bb2 = (M_aux2>0).nonzero()
            bb2 = array([[bb2[0].min(), bb2[0].max()], [bb2[1].min(), bb2[1].max()], [bb2[2].min(), bb2[2].max()]])

            if (float(volumen1)/volumen2)>4 or (float(volumen1)/volumen2)<0.25: # son muy diferentes
                M1[:] = M_aux1 # elijo la mayor
            elif (bb1[:,0]>=bb2[:,0]).all() and (bb1[:,1]<=bb2[:,1]).all(): # M_aux1 esta contenido en M_aux2
                M1[:] = M_aux1
            elif (bb1[:,0]<=bb2[:,0]).all() and (bb1[:,1]>=bb2[:,1]).all(): # M_aux2 esta contenido en M_aux1
                M1[:] = M_aux2
            elif (float(volumen1)/volumen2)<2 and (float(volumen1)/volumen2)>0.5: # no son muy diferentes en tamaño
                CM1 = ndimageMeasurements.center_of_mass(M_aux1, labels = None, index = None)
                CM2 = ndimageMeasurements.center_of_mass(M_aux2, labels = None, index = None)    
                if CM1[0]>CM2[0]: # asumo que el cuello es el otro y que está mas abajo
                    M1[:] = M_aux1
                else:
                    M1[:] = M_aux2
                del CM1, CM2
            else:
                if volumen1>volumen2:
                    M1[:] = M_aux1
                else:
                    M1[:] = M_aux2

            del M_aux1,M_aux2, volumen1, volumen2
            
        else:
            M1[:] = L[0]==I1;

        #cambie las iteraciones (eran 3)
        ndimageMorphology.binary_dilation(M1.copy(), structure = SE, iterations = IteracionesDilatationM1, mask = None, output = M1, border_value = 0, origin = 0, brute_force = False)

        B2 = B1 * M1 #dejar solo el cerebro y algunos tejidos alrededor

        if DIBUJO_TESIS:
            if RESHAPE111:
                res = medimages.save_ANALYZE(B2.astype('int16'),'Debug_T_stripping',spacing=[1,1,1], folder = folderOut)
            else:
                res = medimages.save_ANALYZE(B2.astype('int16'),'Debug_T_stripping',spacing=hdr['Spacing'], folder = folderOut)    

        if 0: # dibujar
            figure()
            imshow(B2[:,:,100], interpolation='nearest', origin='lower')
            show()

        del M1, B1




    # ###############################################################

    histograma2 = ndimageMeasurements.histogram(B2, min_image, max_image, int(nbins), labels=None, index=None)
    histograma2[-1] += (B2==max_image).sum() # porque no toma en cuenta <= del ultimo bin

    
    if GRAFICOS_Hist: #dibujar
        width = abs(centros_hist[0] - centros_hist[1])
        figure()
        bar(centros_hist - (width), histograma2, width = width, bottom = 0, color = 'g' , yerr = None, xerr = None, ecolor = 'k', capsize = 3)
##                bar(centros_hist - (width), histograma2/histograma2.sum(), width = width, bottom = 0, color = 'g' , yerr = None, xerr = None, ecolor = 'k', capsize = 3)
##                plot(centros_hist,histograma2)
        title('Histograma 2 ' + imageFile)
        show()


##             # ######## filtro histograma - optativo #########
##             # para tener una curba pareja
##            if FiltroHistograma:
##                for i in range(5):
##                    umb=histograma2.max() * 0.01
##                    for i in range(1,int(nbins)-1):    
##                        if (histograma2[i-1]-histograma2[i] > umb) and (histograma2[i+1] - histograma2[i] > umb):
##                            histograma2[i] = (histograma2[i-1] + histograma2[i+1]) / 2
##             # ###############################################

##             # ######## filtro histograma - optativo #########
##             # para tener una curba pareja
##            if FiltroHistograma:
##                umb=histograma2.max() * 0.1
##                for i in range(1,int(nbins)-1):    
##                    if (histograma2[i-1]-histograma2[i] > umb) and (histograma2[i+1] - histograma2[i] > umb):
##                        histograma2[i] = (histograma2[i-1] + histograma2[i+1]) / 2
##             # ###############################################

    if GRAFICOS_Hist: #dibujar
        width = abs(centros_hist[0] - centros_hist[1])
        figure()
        bar(centros_hist - (width), histograma2, width = width, bottom = 0, color = 'g' , yerr = None, xerr = None, ecolor = 'k', capsize = 3)
##                bar(centros_hist - (width), histograma2/histograma2.sum(), width = width, bottom = 0, color = 'g' , yerr = None, xerr = None, ecolor = 'k', capsize = 3)
##                plot(centros_hist,histograma2)
        title('Histograma 2 filtrado ' + imageFile)
        show()


    # ########### parametros iniciales
    # ###################### kernel density estimates

    histograma2 = histograma2.astype('float32')
    N = histograma2.sum()
    
    histograma_aux = histograma2.copy()
    histograma_aux_index = (histograma_aux > (histograma_aux.max()*0.1)).nonzero()[0] # para borrar las colas producen problemas al encontrar las modas
    histograma_aux[:histograma_aux_index[0]] = 0.
    histograma_aux[histograma_aux_index[-1]:] = 0.

    K=2. #numero de modos que se quieren encontrar (centroides)
    h_alto = nbins/10.
    h_bajo = 0.1;
    h = (h_alto+h_bajo)/2.
    dif1D = array([1, -1], 'float32')

    contt = 0;
    pp = arange(nbins)
    
    while 1:
        p = zeros(nbins, 'float32')    
        cont = -1
        for i in centros_hist:    
            cont += 1
            aux = 0
            for ii in range(int(nbins)):
                aux += ((1./sqrt(2.*pi)) * exp((-( ((i - centros_hist[ii]) * (1./h) )**2. ))/2.)) * histograma_aux[ii]
            p[cont] = (1./(N*h)) * aux;         

        if 0: # dibujos
            figure()
        ##    plot(pp,p)
            plot(range(len(p)),p)
        ##    h
            show()
            
        p_prim_der = convolve(p, dif1D, mode = 'valid')
    ##    figure()
    ##    plot(p_prim_der)
    ##    show()
        
        p_sec_der = convolve(p_prim_der, dif1D, mode = 'valid')
    ##    figure()
    ##    plot(p_sec_der)
    ##    show()
        
        modes=[]
        for i in range(len(p_prim_der)-1):
            if (p_prim_der[i]*p_prim_der[i+1]) < 0:
                if p_sec_der[i]<0:
                    ## modes.append(pp[i+1])
                    modes.append(i+1)
                    
        if flag_printData:
            print modes
        
        if len(modes) <= K:
            h_alto = h
        elif len(modes) > K:
            h_bajo = h    
        
        h_ant = h
        h = (h_alto+h_bajo)/2.
        
        contt += 1
       
        
        if ((len(modes)==K) and (abs(h-h_ant)<0.01)):
            if K == 2:
                if p[modes[0]] > p[modes[1]]: # para evitar que elija uno pequeño al de LCR al comienzo y el de MG al medio
                    break
                else:
                    K = 3
                    contt = 0
                    h_alto = nbins/10.
                    h_bajo = 0.05;
                    h = (h_alto+h_bajo)/2.
            else: # K==3:
                modes = [modes[1], modes[2]]
                break
                
##                elif (contt > 34) and (len(modes)==K): # si da muchas iteraciones elijo los dos modos mas grandes
        elif (contt > 34): # si da muchas iteraciones elijo los dos modos mas grandes 
            modes_aux = [0, 0]
            i_aux = p[modes].argmax()
            modes_aux[0] = modes[i_aux]
            modes[i_aux] = -1
            modes_aux[1] = modes[p[modes].argmax()]
            modes_aux.sort()
            modes = modes_aux
            break


    if GRAFICOS_Hist: #dibujar
            figure()
            plot(centros_hist,p)
            title('Calculo kernel density estimates ' + imageFile)
            show()
            

##    ini=array([ round((modes[0]*0.75)), modes[0], modes[1]])
##    parametrosIni=array([ini[0], ini[1], ini[2], 1, 1, 1, histograma2[ini[0]], histograma2[ini[1]], histograma2[ini[2]] ])
##            modes = [101, 135]

    rangDesv = (histograma2> (histograma2.max() * 0.01)).nonzero()
    rangDesv = rangDesv[0][-1] - rangDesv[0][0]
    desvIni = rangDesv/3./3.
        
    ini=array([ round((modes[0]*0.75)), modes[0], modes[1]])
    ini_h = (histograma_aux.sum()/p.sum()) * p[ini.astype('int')]
##    parametrosIni=array([ini[0], ini[1], ini[2], desvIni, desvIni, desvIni, histograma2[ini[0]], histograma2[ini[1]], histograma2[ini[2]] ])
    parametrosIni=array([ini[0], ini[1], ini[2], desvIni, desvIni, desvIni, ini_h[0], ini_h[1], ini_h[2] ])
    
    
    # ####fin parametros iniciales####


    # ##### Calculo de parametros optimo
    D = [0.004, 0.004, 0.004, 0.0005, 0.0005, 0.0005, 1, 1, 1]
    FF = 0.1
    parametrosOpt = leastsq(HistogramAdjustment, parametrosIni.copy(), args = (histograma2, parametrosIni.copy()), Dfun = None, full_output = 1, col_deriv = 0, ftol = 1.49012e-008, xtol = 1.49012e-008, gtol = 0.0, maxfev = 0, epsfcn = 0.0, factor = FF, diag = D)
    ##    from scipy.optimize.optimize import *
    ##parametrosOpt = fmin_powell(HistogramAdjustment, parametrosIni.copy(), args = (histograma2, ), xtol = 0.0001, ftol = 0.0001, maxiter = None, maxfun = None, full_output = 0, disp = 1, retall = 0, callback = None)
    ##parametrosOpt = fmin(HistogramAdjustment, parametrosIni.copy(), args = (histograma2, ), xtol = 0.0001, ftol = 0.0001, maxiter = None, maxfun = None, full_output = 1, disp = 1, retall = 0, callback = None)

##    print 'Error de aporoximacion:' + str(sqrt(parametrosOpt[2]['fvec']*parametrosOpt[2]['fvec']).sum())
    

    # dibujar los histogramas
    if GRAFICOS_Hist:
        u = parametrosOpt[0][0:3] #punto medio de la función
        s = parametrosOpt[0][3:6] #desviación estandar
        Np = parametrosOpt[0][6:9] #altura en punto medio
        
        diferencia = 0
        histograma_parametros = zeros(len(histograma), 'float32')
        G = zeros((3,len(histograma)), 'float32');
        for i in range(len(histograma)):
            histograma_parametros[i]=0
            for k in range(3):
                histograma_parametros[i] += Np[k]*exp((-1./2.)*(((i-u[k])/s[k])**2.))
                G[k,i] += Np[k]*exp((-1./2.)*(((i-u[k])/s[k])**2.))
        diferencia = (histograma_parametros-histograma2)
        
        if flag_printData:
            print 'diferencia histogramas = ' + str(sqrt(diferencia*diferencia).sum())

        figure()
        plot(histograma_parametros,'r')
        hold(True)
        plot(histograma2,'k--')
    ##    figure()
        plot(G[0,:],'g')
    ##    figure()
        plot(G[1,:],'y')
    ##    figure()
        plot(G[2,:],'b')
        title('Ajuste Gaussianas' + imageFile)
        # ######################
##                figure()
##                plot(centros_hist,histograma2)
##
##                figure()
##                plot(histograma2)
        show()

    # ##### Fin - Calculo de parametros optimo
    
    u = parametrosOpt[0][0:3]
    s = parametrosOpt[0][3:6]
    Np = parametrosOpt[0][6:9]

##    u_nivel_gris=(min_image-(((max_image-min_image)/nbins)/2.))+((max_image-min_image)/nbins)*u
##    s_nivel_gris=abs(((max_image-min_image)/nbins)*s)
    u_nivel_gris = min_image+((max_image-min_image)/nbins)*u
    s_nivel_gris=abs(((max_image-min_image)/nbins)*s)

    # #### para pasar informacion de los niveles de gris al paso ########
    # #### de segmentracion con modelos deformables              ########
    alfa = (max_image-min_image)/(max_imagenIni-min_imagenIni)
    beta = max_image-(alfa*max_imagenIni)
    u_nivel_grisOrig = (u_nivel_gris - beta)/alfa
    s_nivel_grisOrig = (s_nivel_gris - beta)/alfa
    # ################################

    M2 = ((u_nivel_gris[1]-2.5*s_nivel_gris[1]<=B2) & (B2<=u_nivel_gris[2]+2.5*s_nivel_gris[2]))

    # ## opening of M2
    radio = RadioOpeningM2
    SE = zeros(((radio*2)+1,(radio*2)+1,(radio*2)+1))
    for i in range(-radio,radio+1):
        for j in range(-radio,radio+1):
            for z in range(-radio,radio+1):
                if sqrt(i**2+j**2+z**2)<=radio:
                    SE[i+radio,j+radio,z+radio] = 1
    ndimageMorphology.binary_opening(M2.copy(), structure = SE, iterations = IteracionesOpeningM2, output = M2, origin = 0)                

    L = ndimageMeasurements.label(M2, structure = None, output = None)
    nL = ndimageMeasurements.histogram(L[0], -0.5, L[1]+0.5, L[1]+1, labels=None, index=None)
    nL[0] = -1 #para no tomar en cuenta los voxeles con 0, (fondo)
    I = ndimageMeasurements.maximum_position(nL, labels = None, index = None)[0]
    M2 = L[0]==I

    B3 = B2 * M2 # eliminar los ultimos tejidos que quedan    

    if 0: # dibujar
        figure()
        imshow(B3[:,:,100], interpolation='nearest', origin='lower')
        colorbar()
        show()

    if GRAFICOS_Hist: # dibujar
        histograma3 = ndimageMeasurements.histogram(B3, min_image, max_image, int(nbins), labels=None, index=None)
        histograma3[-1] += (B3==max_image).sum() # porque no toma en cuenta <= del ultimo bin                
        width = abs(centros_hist[0] - centros_hist[1])
        figure()
        bar(centros_hist - (width), histograma3, width = width, bottom = 0, color = 'g' , yerr = None, xerr = None, ecolor = 'k', capsize = 3)
##                bar(centros_hist - (width), histograma3/histograma3.sum(), width = width, bottom = 0, color = 'g' , yerr = None, xerr = None, ecolor = 'k', capsize = 3)
##                plot(centros_hist,histograma3)
        title('Histograma 3 ' + imageFile)                
        show()

    # #######################################################################
    
    B3 = imagenOrig * M2 # eliminar los ultimos tejidos que quedan
    if RESHAPE111:
        spacing_sux = hdr['Spacing']
        spacing_sux[0] = 1
        spacing_sux[1] = 1
        spacing_sux[2] = 1
    else:        
        spacing_sux = hdr['Spacing']
    res = medimages.save_ANALYZE((B3).astype('int16'),prePegmentedImageFile,spacing=hdr['Spacing'], folder = folderOut)


    del M2, B2, B3

    # datos para la segmentacion posterior
    skullStrippingData = {}
    skullStrippingData['Totsu_grisOrig'] = Totsu_grisOrig
    skullStrippingData['u_nivel_grisOrig'] = u_nivel_grisOrig
    skullStrippingData['s_nivel_grisOrig'] = s_nivel_grisOrig
    SaveData(workFolder + dataNames['dataPreSeg'], ('skullStrippingData',), (skullStrippingData,))          

def TransformationAdjustmentAffin(parametros, pointsIn, distance):
##    PointsIn = PointsIn_distance[0]
##    distance = PointsIn_distance[1]    
    
   ## print parametros1
    transformation = zeros((4,4), 'float32')
    transformation[:3,:3] = parametros[:9].reshape(3,3)
    transformation[:3,3] = parametros[9:]
    transformation[3,3] = 1

    points = concatenate((pointsIn.T, ones((1,pointsIn.shape[0]))))    
    
    ##Points = ((dot(Transformation, Points)).round().astype('int'))[:3,:]
    points = dot(transformation, points)[:3,:]

##    L = array(distance.shape).reshape(-1,1)
##    
##    PointsOut=((Points>L) | (Points<0)).any(0).nonzero() # puntos fuera de la matriz de la imagen
##    
##    Points[:,PointsOut] = 0
##    Values = zeros((Points.shape[1]))
    
    values = map_coordinates(distance, points, output_type = None, output = None, order = 1, mode = 'nearest', cval = 0.0, prefilter = False)
  
##    for i in range(Points.shape[1]):
##        Values[i] = distanceInterpolarion(Points[0,i], Points[1,i], Points[2,i])
        
##    print Values.sum()
    
    return values   

def TransformationAdjustment(parametros1, pointsIn, distance, paramatrosAjust):
##    PointsIn = PointsIn_distance[0]
##    distance = PointsIn_distance[1]    
    
    parametros = parametros1 * paramatrosAjust

   ## print parametros1
            
    tetax = parametros[0]
    tetay = parametros[1]
    tetaz = parametros[2]
    sx = parametros[3]
    sy = parametros[4]
    sz = parametros[5]
    tx = parametros[6]
    ty = parametros[7]
    tz = parametros[8]
    R = array([[cos(tetay)*cos(tetaz), -cos(tetay)*sin(tetaz), -sin(tetay), 0.],
             [-sin(tetax)*sin(tetay)*cos(tetaz)+cos(tetax)*sin(tetaz), sin(tetax)*sin(tetay)*sin(tetaz)+cos(tetax)*cos(tetaz), -sin(tetax)*cos(tetay), 0.],
             [cos(tetax)*sin(tetay)*cos(tetaz)+sin(tetax)*sin(tetaz), -cos(tetax)*sin(tetay)*sin(tetaz)+sin(tetax)*cos(tetaz), cos(tetax)*cos(tetay), 0.],
             [0., 0., 0., 1.]], 'float32')
    ST = array([[sx, 0., 0., tx],
               [0., sy, 0., ty],
               [0., 0., sz, tz],
               [0., 0., 0., 1.]], 'float32')
    SS = array([[1, tan(Syx), tan(Szx), 0],
                [tan(Sxy), 1, tan(Szy), 0],
                [tan(Sxz), tan(Syz), 1, 0],
                [0, 0, 0, 0]], 'float32')
    transformation = dot(ST, R)

    points = concatenate((pointsIn.T, ones((1,pointsIn.shape[0]))))    
    
    ##Points = ((dot(Transformation, Points)).round().astype('int'))[:3,:]
    points = dot(transformation, points)[:3,:]

##    L = array(distance.shape).reshape(-1,1)
##    
##    PointsOut=((Points>L) | (Points<0)).any(0).nonzero() # puntos fuera de la matriz de la imagen
##    
##    Points[:,PointsOut] = 0
##    Values = zeros((Points.shape[1]))
    
    values = map_coordinates(distance, points, output_type = None, output = None, order = 1, mode = 'nearest', cval = 0.0, prefilter = False)
  
##    for i in range(Points.shape[1]):
##        Values[i] = distanceInterpolarion(Points[0,i], Points[1,i], Points[2,i])
        
##    print Values.sum()
    
    return values   

def AppplyTransformation(points_aux, transformation):
    points_aux = concatenate((points_aux.T, ones((1,points_aux.shape[0]), 'float32')))
    points_aux = dot(transformation, points_aux)[:3,:]
    points_aux = points_aux.T
    return points_aux

def AppplyTransformationAffin(points_aux, parametros):
    transformation = zeros((4,4), 'float32')
    transformation[:3,:3] = parametros[:9].reshape(3,3)
    transformation[:3,3] = parametros[9:]
    transformation[3,3] = 1
    
    points_aux = concatenate((points_aux.T, ones((1,points_aux.shape[0]), 'float32')))
    points_aux = dot(transformation, points_aux)[:3,:]
    points_aux = points_aux.T
    return points_aux








def MakeImageBinAndPointsFromMesh(mesh, spacing, box):
    # Datos para el muestreo ###
    origen_muestreo = box[0,:]
    fin_muestreo = box[1,:]
    spacing_muestreo = spacing

    xp = arange(origen_muestreo[0], fin_muestreo[0], spacing_muestreo[0])
    yp = arange(origen_muestreo[1], fin_muestreo[1], spacing_muestreo[1])
    zp = arange(origen_muestreo[2], fin_muestreo[2], spacing_muestreo[2])

    xn = len(xp)
    yn = len(yp)
    zn = len(zp)
    # #####
    [xp,yp,zp] = meshgrid2(xp,yp,zp)           
    qp = concatenate((xp.reshape(-1,1),yp.reshape(-1,1), zp.reshape(-1,1)), 1)
    del xp,yp,zp

    # puntos dentro de la malla
##    TT = timeit.Timer()
##    T1 = TT.timer()   
    inn = InPolyedron(mesh.points ,mesh.triangles , mesh.normalsOfTriangles ,qp)
##    T2 = TT.timer()
##    ttt = T2-T1
##    print ttt

    indexBin = inn.nonzero()[0]
    pointsBin = qp[indexBin,:]
    
    # creacion de imagen binaria
    imagenBin = zeros((xn * yn * zn), 'int8')
    imagenBin[indexBin] = 1
    imagenBin = imagenBin.reshape(zn, yn, xn) # se crea con x,z invertida para que quede en la posicion correcta al pasar a vtk o analize

    return [pointsBin, imagenBin] 

##rootFolder = 'E:\\Users\\FG\\Doctorado\\'
##rootWorkFolder = rootFolder + 'Validacion de Segmentacion\\'
##rootImagesFolder = rootFolder + 'Imagenes\\'

    
flag_steps = {}
flag_steps['segmentation'] = 1

segmentation_Flags = {}
segmentation_Flags['flag_segStep_craneo'] = 0
segmentation_Flags['flag_segStep_registro'] = 0
segmentation_Flags['flag_segStep_deform1'] = 0
segmentation_Flags['flag_segStep_deform2'] = 0
segmentation_Flags['flag_segStep_deform3'] = 0
segmentation_Flags['flag_createBinImage'] = 0
segmentation_Flags['flag_segStep_CondDilatation'] = 0
segmentation_Flags['flag_masking'] = 1

segmentation_Options = {}
segmentation_Options['flag_SimplexTrianMethodBinImag'] = 'tangentesNotDual' # 'tangentesNotDual', 'tangentesDual'



TT = timeit.Timer()
T1 = TT.timer()
if 1: #esto lo cambio por el "try" anterior
    if flag_steps['segmentation']: #segmentacion
        SegmentationSteps(dataNames, segmentation_Flags)

T2 = TT.timer()
print 'Tiempo:',T2-T1