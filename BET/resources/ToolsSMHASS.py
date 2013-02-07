from numpy import *
##from Scientific.Functions.Interpolation import InterpolatingFunction
from scipy.ndimage.interpolation import affine_transform, map_coordinates
from scipy import sparse
from numpy.linalg import det, norm
from vtk import vtkPolyData, vtkFloatArray, vtkPoints, vtkIdTypeArray, vtkCellArray #solo para MakePolyDataFromSimplex

from vtkImageExportToArray import *
from vtkImageImportFromArray import *
from vtk import *
##from numpy import *
import medimages

import scipy.ndimage.filters as ndimageFilters
import time
import timeit
import os
from scipy.optimize.optimize import *
from scipy.optimize.cobyla import *
import copy
import glob
import cPickle
import gzip

import scipy.ndimage.measurements as ndimageMeasurements

##os.chdir('C:\\Chubo\\Doctorado\\Python\\Programas\\Para_publicar')

##try:
##    os.chdir('C:\\Chubo\\Doctorado\\Python\\Programas')
##except:
##    os.chdir('E:\\Users\\FG\\Doctorado\\Python\\Programas')


class SimplexMesh:
    def __init__(self, points = None , neighbors = None, faces = None, pointToFace = None, meshes = None, contoursPoints = None, facesSimplexContour = None): # constructor

        if points!=None:
            self.points = points.astype('float32')
        if neighbors!=None:
            self.neighbors = neighbors.astype('int')
        if faces!=None:
            self.faces = faces[:]
        if pointToFace!=None:
            self.pointToFace = pointToFace
        if meshes!=None:
            self.meshes = meshes
        if contoursPoints!=None:
            self.contoursPoints = contoursPoints
        if facesSimplexContour!=None:
            self.facesContour = facesSimplexContour
        self.flag_printData = 0
    def ComputeNormalsOfPoints(self):
        """ Calcula las normales de una malla simplex"""
        self.normalsOfPoints = cross(self.points[self.neighbors[:,1],:]-self.points[self.neighbors[:,0],:], self.points[self.neighbors[:,2],:]-self.points[self.neighbors[:,0],:])
        normals_norm = reshape(sqrt(sum(self.normalsOfPoints*self.normalsOfPoints,1)), (-1, 1))
        self.normalsOfPoints[:,:] = self.normalsOfPoints / normals_norm
    def ComputeAreaOfBaseTriangles(self):
        w = self.points[self.neighbors[:,1]] - self.points[self.neighbors[:,0]]
        v = self.points[self.neighbors[:,2]] - self.points[self.neighbors[:,0]]
        c = cross(w, v)
        self.areaOfBaseTriangles = sqrt((c * c).sum(1)) * 0.5
    def ComputeMeanCurvatureInPoints(self, method = 'simplex'):
        if method == 'vtk':
            "Calcula la curvatura media en cada punto"
            PtoVTK = PythonToPolyData()
            PtoVTK.SetInputSimplex(self)
            vtkSimplexMeshDic = PtoVTK.GetVTKPolyData()
            vtkSimplexMesh = vtkSimplexMeshDic['vtkData']
            vtk_Curvatures = vtkCurvatures()
            vtk_Curvatures.SetCurvatureTypeToMean()
            vtk_Curvatures.SetInput(vtkSimplexMesh)
            vtk_Curvatures.Update()
            vtk_CurvaturesPolyData = vtk_Curvatures.GetOutput()

            vtkFloatArrayScalars = vtk_CurvaturesPolyData.GetPointData().GetScalars()
            arrayScalars = zeros((vtkFloatArrayScalars.GetNumberOfTuples(),vtkFloatArrayScalars.GetNumberOfComponents()), 'float64')
            vtkFloatArrayScalars.ExportToVoidPointer(arrayScalars)

            self.meanCurvatureInPoints = arrayScalars
            
        if method == 'simplex':
            if not hasattr(self, 'angleSimplex'):
                self.ComputeAngleSimplex()            
            self.meanCurvatureInPoints = sin(self.angleSimplex) / self.radSimplexCircle
        
    def ComputeAngleSimplex(self):
        """
        Calcula los angulos simplex de la malla
        """

        if not hasattr(self, 'normals'):
            self.ComputeNormals()  

        neighbors1 = self.points[self.neighbors[:,0],:]
        neighbors2 = self.points[self.neighbors[:,1],:]
        neighbors3 = self.points[self.neighbors[:,2],:]


        #calculo del centro de la esfera que incluye los 4 puntos
        angulos_cero = zeros(self.points.shape[0], 'int8')

        Csphere = zeros(self.points.shape, 'float32')

        A = zeros((4,4), 'float32')
        A[:,3] = 1;

        for i in range(self.points.shape[0]):
            A[0,:3] = self.points[i,:]
            A[1,:3] = neighbors1[i,:]
            A[2,:3] = neighbors2[i,:]
            A[3,:3] = neighbors3[i,:]
            bajo = det(A)
            
            if abs(bajo) > 0.00001:

                A[0,0] = sum(self.points[i,:]**2); A[0,1:3] = self.points[i,1:]
                A[1,0] = sum(neighbors1[i,:]**2); A[1,1:3] = neighbors1[i,1:]
                A[2,0] = sum(neighbors2[i,:]**2); A[2,1:3] = neighbors2[i,1:]
                A[3,0] = sum(neighbors3[i,:]**2); A[3,1:3] = neighbors3[i,1:]
                Csphere[i,0] = det(A) / (2 * bajo)

                A[0,1] = self.points[i,0]; A[0,2] = self.points[i,2]
                A[1,1] = neighbors1[i,0]; A[1,2] = neighbors1[i,2]
                A[2,1] = neighbors2[i,0]; A[2,2] = neighbors2[i,2]
                A[3,1] = neighbors3[i,0]; A[3,2] = neighbors3[i,2]
                Csphere[i,1] = -det(A) / (2 * bajo)

                A[0,1:3] = self.points[i,:2]
                A[1,1:3] = neighbors1[i,:2]
                A[2,1:3] = neighbors2[i,:2]
                A[3,1:3] = neighbors3[i,:2]
                Csphere[i,2] = det(A) / (2 * bajo)
                
            else:
                Csphere[i,:] = 0
                angulos_cero[i]=1


        # calculo del centro del circulo que incluye a los vecinos

        v12 = neighbors2 - neighbors1
        v13 = neighbors3 - neighbors1
        v23 = neighbors3 - neighbors2
        a = sqrt(sum(v23**2,1)) #1
        b = sqrt(sum(v13**2,1)) #2
        c = sqrt(sum(v12**2,1)) #3
        C1 = (a**2)*(b**2 + c**2 - a**2)
        C2 = (b**2)*(a**2 + c**2 - b**2)
        C3 = (c**2)*(a**2 + b**2 - c**2) 
        ww = C1 + C2 + C3
        ww[ww == 0] = 0.00001 # ajuste num
        C1 = C1 / ww
        C2 = C2 / ww
        C3 = C3 / ww
        Ccircle = neighbors1 * C1.reshape(-1,1) + neighbors2 * C2.reshape(-1,1) + neighbors3 * C3.reshape(-1,1)


        # calculo del radio del circulo y del radio de la esfera

        Rsphere = sqrt(sum((Csphere - neighbors1)**2, 1))
        Rcircle = sqrt(sum((Ccircle - neighbors1)**2, 1))


        # calculo de angulos simplex    
        aux = ( Rcircle/Rsphere ) * ((sum((neighbors1 - self.points) * self.normals,1)<0)*2-1)
        aux[aux > 1] = 1 # ajuste num
        aux[aux > -1] = -1 # ajuste num
        angleSimplex = arcsin(  ( Rcircle/Rsphere ) * ((sum((neighbors1 - self.points) * self.normals,1)<0)*2-1) )
        angleSimplex[(Rcircle>Rsphere).nonzero()]=pi/2
        angleSimplex[angulos_cero.nonzero()] = 0

        self.angleSimplex = angleSimplex
        self.radSimplexCircle = Rcircle
        self.radSimplexSphere = Rsphere

    def ComputeRadCircle(self):
        if not hasattr(self, 'normals'):
            self.ComputeNormals()  

        neighbors1 = self.points[self.neighbors[:,0],:]
        neighbors2 = self.points[self.neighbors[:,1],:]
        neighbors3 = self.points[self.neighbors[:,2],:]

        # calculo del centro del circulo que incluye a los vecinos
        # por calculo de Circumcenter

        v12 = neighbors2 - neighbors1
        v13 = neighbors3 - neighbors1
        v23 = neighbors3 - neighbors2

        # version 1
        dca = (v13 * v12).sum(1).reshape(-1,1)
        dba = (v23 * -v12).sum(1).reshape(-1,1)
        dcb = (-v13 * -v23).sum(1).reshape(-1,1)
        
        n1 = dba * dcb
        n2 = dcb * dca
        n3 = dca * dba

        Rcircle = sqrt(((dca + dba) * (dba + dcb) * (dcb + dca)) / (n1 + n2 + n3)) / 2.
        return Rcircle
##        Ccircle = ((n2 + n3) * neighbors1 + (n3 + n1) * neighbors2 + (n1 + n2) * neighbors3) / (2 * (n1 + n2 + n3))
##
##        # calculo del radio del circulo y del radio de la esfera
##
##        Rsphere = sqrt(sum((Csphere - neighbors1)**2, 1))
##        Rcircle = sqrt(sum((Ccircle - neighbors1)**2, 1))

        
    def ComputeNormals(self):
        
        normals = cross(self.points[self.neighbors[:,1],:] - self.points[self.neighbors[:,0],:], self.points[self.neighbors[:,2],:] - self.points[self.neighbors[:,0],:])
        normals_norm = reshape(sqrt(sum(normals*normals,1)), (-1, 1))
        normals[:,:] = normals / normals_norm
        self.normals = normals

    def SimplexToTriangles4(self):
        """
        v.2
        Es la triangulacion Dual    
        Ahora incluye manejo de bodes
        Ocupa una estimacion del punto central (estimacion por tangentes) para calcular los pesos de las distancias a los puntos de la cara
        """

        if not hasattr(self, 'areaOfBaseTriangles'):    
            self.ComputeAreaOfBaseTriangles()
            
        TT = timeit.Timer()  
        T1 = TT.timer()
        if not hasattr(self, 'normalsOfPoints'):
            self.ComputeNormalsOfPoints()
    
        rCircle2 = self.ComputeRadCircle()  #  no deberia ser el cuadrado?
                                            # esta parte puede traer problemas si los puntos son colineales => rad -> Inf
##        rCircle2 = rCircle2 * rCircle2
        
        Lfaces = zeros(len(self.faces))
        for i in range(len(self.faces)): # cuantas caras con i numero de puntos hay
            Lfaces[i] = self.faces[i].shape[0]
            
        pointsContours = array([], 'int')    
        for j in range(len(self.contoursPoints)):
            pointsContours = concatenate((pointsContours, self.contoursPoints[j]), 1) # poner todos los puntos pertenecientes a todos los contornos de la malla en una sola lista
                
        Nfaces = Lfaces.sum()    # numero total de caras
        points = zeros((Nfaces,3), 'float32') # los puntos que habran en la triangulacion
        triangles = zeros((self.points.shape[0] - len(pointsContours), 3), 'int') # los triangulos que habran en la triangulacion
        ubicacion_puntos = zeros((len(self.faces), Lfaces.max()), 'int') # referencia para encontrar a que cara pertecece cada punto simplex

        pointToTriangle = [] #lista del largo de numero de puntos donde cada elemento es una lita con los triangulos a los que pertecece el punto
        for i in range(int(Nfaces)): pointToTriangle.append([])
        
        
        n = 0
        planes = concatenate((self.normalsOfPoints, - (self.points * self.normalsOfPoints).sum(1).reshape(-1,1)), 1) # calculo de los planos correspondientes a cada punto
        T2 = TT.timer()
        if self.flag_printData:
            print 'Tiempo1 StoTTang'
            print T2-T1
        
##        list_beta = [] #BORRAR
        for i in range(1, len(self.faces)): # creacion de los puntos de la triangulacion, un punto por cara simplex (los del bode se corregiran luego)
            w = 0.5 / i
            for j in range(self.faces[i].shape[0]):
                ubicacion_puntos[i,j] = n

                pointLinst = self.faces[i][j]
                points_aux = self.points[pointLinst,:]
                
                # calculo de pesos
##                areas = self.areaOfBaseTriangles[pointLinst]
##                sumAreas = areas.sum()
##                alfa = areas / sumAreas

                areas = rCircle2[pointLinst]
                sumAreas = areas.sum()
                alfa = areas / sumAreas


                

                pCentral = points_aux.sum(0) / i
                q = zeros(3, 'float32')
                q = ((self.normalsOfPoints[pointLinst] * (points_aux - pCentral)).sum(1).reshape(-1,1) * self.normalsOfPoints[pointLinst]).sum(0)
                pointAprox = pCentral + q * w
                pointAprox = concatenate((pointAprox, [1]))
                pointAprox = pointAprox.reshape(-1,1)                
        
                A = zeros((1,3), dtype='float32')
                for one_alfa, t in zip(alfa, pointLinst):
                    A  += one_alfa * dot(planes[t,:].reshape(1,-1), pointAprox) * self.normalsOfPoints[t]
                B = (points_aux.sum(0).reshape(-1,1) - i * pointAprox[:3])
                    
                dotB = dot(B.T,B)
                if dotB == 0.:
                    beta = 2.
                else:
                    beta = (dot(A,B)/dot(B.T,B))[0,0]
                    if beta < 0.4:
                        beta = 0.4
                    elif beta > 2.:
                        beta = 2.

##                list_beta.append(beta) #BORRAR
                
##                if beta < 0.4:
##                    beta = 0.4
##                elif beta > 2.:
##                    beta = 2.
                # aporte de los puntos de la cara
                Q = zeros((4,4), dtype='float32')
                Q_aux = eye(4, dtype='float32')
                Q_aux[3,3] = 0.
                aux = zeros((3), 'float32')

                for k in points_aux:
                    Q_aux[3,3] = dot(k, k)
                    Q_aux[:3,3] = -k
                    Q_aux[3,:3] = -k
                    Q += beta * Q_aux

                # aporte de los planos de los puntos de la cara
                for one_alfa, t in zip(alfa, pointLinst):                    
                    Q  += one_alfa * dot(planes[t,:].reshape(-1,1), planes[t,:].reshape(1,-1))

                        
                
                points[n,:] = (dot(linalg.inv(Q[:3,:3]),-Q[:3,3].reshape(-1,1))).flatten()

                n += 1

        T2 = TT.timer()
        if self.flag_printData:
            print 'Tiempo2 StoTTang'
            print T2-T1


##        from pylab import * #BORRAR
##        figure()
##        plot(list_beta)
##        show()
                
        # ######################################################################################   
        # ## ACA DEBO AJUSTAR LA UBICACION DE LOS PUNTOS DEL CONTORNO DE LA MALLA TRIANGULAR
        # ####################################################################################
        for i in range(self.facesContour.shape[0]): # para ubicar bien los puntos que pertenecen al borde de la malla triangular, sino queda al medio de la cara
            aux = intersect1d(pointsContours, self.faces[self.facesContour[i,0]][self.facesContour[i,1],:])
            points[ubicacion_puntos[self.facesContour[i,0], self.facesContour[i,1]]] = self.points[aux].sum(0) / 2.
            
        # creacion de los triangulos   
        n = 0
        indicesSelection = ones(self.pointToFace.shape[0]); indicesSelection[pointsContours] = 0  #para eliminar los puntos del contorno, como pertenecen solo a 2 caras no corresponde crear un triangulo
        indices = arange(self.pointToFace.shape[0]); indices = indices[indicesSelection.nonzero()]
        for i in indices :
            P1 = ubicacion_puntos[self.pointToFace[i,0,0], self.pointToFace[i,0,1]]
            P2 = ubicacion_puntos[self.pointToFace[i,1,0], self.pointToFace[i,1,1]]
            P3 = ubicacion_puntos[self.pointToFace[i,2,0], self.pointToFace[i,2,1]]
            triangles[n,:] =[P1, P2, P3]
            pointToTriangle[P1].append(n)
            pointToTriangle[P2].append(n)
            pointToTriangle[P3].append(n)    
            n += 1
        
        OUT = TrianglesMesh(points, triangles) 
        OUT.pointToTriangle = pointToTriangle
        return OUT

    def SimplexToTrianglesNotDualConservacion(self):
        """
        v.2
        Es la triangulacion Dual    
        Ahora incluye manejo de bodes
        Ocupa una estimacion del punto central (estimacion por tangentes) para calcular los pesos de las distancias a los puntos de la cara
        """

        if not hasattr(self, 'areaOfBaseTriangles'):    
            self.ComputeAreaOfBaseTriangles()
            
        if not hasattr(self, 'normalsOfPoints'):
            self.ComputeNormalsOfPoints()

        rCircle2 = self.ComputeRadCircle()
##        rCircle2 = rCircle2 * rCircle2
        
        Lfaces = zeros(len(self.faces))
        for i in range(len(self.faces)): # cuantas caras con i numero de puntos hay
            Lfaces[i] = self.faces[i].shape[0]
            
        pointsContours = array([], 'int')    
        for j in range(len(self.contoursPoints)):
            pointsContours = concatenate((pointsContours, self.contoursPoints[j]), 1) # poner todos los puntos pertenecientes a todos los contornos de la malla en una sola lista

        Ntriangles = 0
        for i in range(len(Lfaces)):
            Ntriangles = Ntriangles + (i * Lfaces[i])
        Nfaces = Lfaces.sum()    # numero total de caras
        Npoints = self.points.shape[0] + Nfaces        
        NpointsC = Nfaces
        
        points = zeros((Npoints,3), 'float32') # los puntos que habran en la triangulacion
        pointsC = zeros((NpointsC,3), 'float32') # los puntos que se crean al centro de las caras simplex
        triangles = zeros((Ntriangles, 3), 'int') # los triangulos que habran en la triangulacion
        ubicacion_puntos = zeros((len(self.faces), Lfaces.max()), 'int') # referencia para encontrar a que cara pertecece cada punto simplex

        pointCToTriangle = [] #lista del largo de numero de puntos donde cada elemento es una lita con los triangulos a los que pertecece el punto
        for i in range(int(NpointsC)): pointCToTriangle.append([])        
        
        n = 0
        planes = concatenate((self.normalsOfPoints, - (self.points * self.normalsOfPoints).sum(1).reshape(-1,1)), 1) # calculo de los planos correspondientes a cada punto
        
##        list_beta = [] #BORRAR
        for i in range(1, len(self.faces)): # creacion de los puntos de la triangulacion, un punto por cara simplex (los del bode se corregiran luego)
            w = 0.5 / i
            for j in range(self.faces[i].shape[0]):
                ubicacion_puntos[i,j] = n

                pointLinst = self.faces[i][j]
                points_aux = self.points[pointLinst,:]
                
                # calculo de pesos
##                areas = self.areaOfBaseTriangles[pointLinst]
##                sumAreas = areas.sum()
##                alfa = areas / sumAreas

                areas = rCircle2[pointLinst]
                sumAreas = areas.sum()
                alfa = areas / sumAreas


                

                pCentral = points_aux.sum(0) / i
                q = zeros(3, 'float32')
                q = ((self.normalsOfPoints[pointLinst] * (points_aux - pCentral)).sum(1).reshape(-1,1) * self.normalsOfPoints[pointLinst]).sum(0)
                pointAprox = pCentral + q * w
                pointAprox = concatenate((pointAprox, [1]))
                pointAprox = pointAprox.reshape(-1,1)                
        
                A = zeros((1,3), dtype='float32')
                for one_alfa, t in zip(alfa, pointLinst):
                    A  += one_alfa * dot(planes[t,:].reshape(1,-1), pointAprox) * self.normalsOfPoints[t]
                B = (points_aux.sum(0).reshape(-1,1) - i * pointAprox[:3])
                    
                dotB = dot(B.T,B)
                if dotB == 0.:
                    beta = 2.
                else:
                    beta = (dot(A,B)/dot(B.T,B))[0,0]
                    if beta < 0.4:
                        beta = 0.4
                    elif beta > 2.:
                        beta = 2.

##                list_beta.append(beta) #BORRAR
                
##                if beta < 0.4:
##                    beta = 0.4
##                elif beta > 2.:
##                    beta = 2.
                # aporte de los puntos de la cara
                Q = zeros((4,4), dtype='float32')
                Q_aux = eye(4, dtype='float32')
                Q_aux[3,3] = 0.
                aux = zeros((3), 'float32')

                for k in points_aux:
                    Q_aux[3,3] = dot(k, k)
                    Q_aux[:3,3] = -k
                    Q_aux[3,:3] = -k
                    Q += beta * Q_aux

                # aporte de los planos de los puntos de la cara
                for one_alfa, t in zip(alfa, pointLinst):
                    Q  += one_alfa * dot(planes[t,:].reshape(-1,1), planes[t,:].reshape(1,-1))
                    
                pointsC[n,:] = (dot(linalg.inv(Q[:3,:3]),-Q[:3,3].reshape(-1,1))).flatten()
                n += 1

##        from pylab import * #BORRAR
##        figure()
##        plot(list_beta)
##        show()
            
        # creacion de los triangulos   
        n = 0
        m = self.points.shape[0]
        k = 0
        for i in range(1, len(self.faces)):
            index_aux = range(1,i) + [0]
            for j in range(self.faces[i].shape[0]):
                triangles[k:k+i,0] = (m+n) * ones(i)
                triangles[k:k+i,1] = self.faces[i][j,:]
                triangles[k:k+i,2] = self.faces[i][j,index_aux]
                k += i
                n += 1
                
        points = concatenate((self.points.copy(),pointsC), 0)
        OUT = TrianglesMesh(points, triangles) 
##        OUT.pointToTriangle = pointToTriangle
        return OUT

class TrianglesMesh:
    def __init__(self, points, triangles): # constructor
        self.points = points.astype('float32')
        self.triangles = triangles.astype('int')
        self.flag_printData = 0
    
    def SearchNeighbors(self):
        self.neighbors = zeros((self.triangles.shape[0],4), 'int') - 1
        self.neighbors[:,3] = 1

    ##	cont=0
    ##	cont1=0
        for i in range(self.triangles.shape[0]):
    ##		cont = cont+1
    ##		if cont==100:
    ##			cont = 0
    ##			cont1 += 1
    ##			print i
            if self.neighbors[i,3]<4:
                revisar=self.triangles[(i+1):,:]
                VFneighbors=self.neighbors[(i+1):,3]<4
                if any(VFneighbors):
                    VF0 = any(self.triangles[i,0]==revisar,1) * any(self.triangles[i,1]==revisar,1)
                    VF1 = any(self.triangles[i,1]==revisar,1) * any(self.triangles[i,2]==revisar,1)
                    VF2 = any(self.triangles[i,2]==revisar,1) * any(self.triangles[i,0]==revisar,1)	
                    if any(VF0):
                        j = VF0.nonzero()[0][0]+i+1
                        self.neighbors[i,0] = j
                        if all( self.triangles[j,[0,1]] == self.triangles[i,[0,1]] ) or all( self.triangles[j,[0,1]]==self.triangles[i,[1,0]] ):
                            self.neighbors[j,0] = i
                        elif all( self.triangles[j,[1,2]] == self.triangles[i,[0,1]] ) or all( self.triangles[j,[1,2]]==self.triangles[i,[1,0]] ):
                            self.neighbors[j,1] = i
                        elif all( self.triangles[j,[2,0]] == self.triangles[i,[0,1]] ) or all( self.triangles[j,[2,0]]==self.triangles[i,[1,0]] ):
                            self.neighbors[j,2] = i
                        self.neighbors[j,3] += 1
                        self.neighbors[i,3] += 1
                    if any(VF1):
                        j=VF1.nonzero()[0][0]+i+1
                        self.neighbors[i,1]=j
                        if all( self.triangles[j,[0,1]]==self.triangles[i,[1,2]] ) or all( self.triangles[j,[0,1]]==self.triangles[i,[2,1]] ):
                            self.neighbors[j,0]=i
                        elif all( self.triangles[j,[1,2]]==self.triangles[i,[1,2]] ) or all( self.triangles[j,[1,2]]==self.triangles[i,[2,1]] ):
                            self.neighbors[j,1]=i
                        elif all( self.triangles[j,[2,0]]==self.triangles[i,[1,2]] ) or all( self.triangles[j,[2,0]]==self.triangles[i,[2,1]] ):
                            self.neighbors[j,2]=i
                        self.neighbors[j,3] += 1
                        self.neighbors[i,3] += 1
                    if any(VF2):
                        j=VF2.nonzero()[0][0]+i+1
                        self.neighbors[i,2]=j
                        if all( self.triangles[j,[0,1]]==self.triangles[i,[2,0]] ) or all( self.triangles[j,[0,1]]==self.triangles[i,[0,2]] ):
                            self.neighbors[j,0]=i
                        elif all( self.triangles[j,[1,2]]==self.triangles[i,[2,0]] ) or all( self.triangles[j,[1,2]]==self.triangles[i,[0,2]] ):
                            self.neighbors[j,1]=i
                        elif all( self.triangles[j,[2,0]]==self.triangles[i,[2,0]] ) or all( self.triangles[j,[2,0]]==self.triangles[i,[0,2]] ):
                            self.neighbors[j,2]=i
                        self.neighbors[j,3] += 1
                        self.neighbors[i,3] += 1
    ##	print i
    ##	if any(self.neighbors[:,3]==4):
    ##		print 'malla no cerrada'
        self.neighbors=self.neighbors[:,0:3]

    ##	print 'ordenados'
        
        ##        return self.neighbors
    def SearchNeighborsUp3(self):
        "para mallas con mas de 3 vecinos por triangulo"
        listTrianglesInVertex = []
        self.neighbors = []
        for i in range(self.points.shape[0]):
            listTrianglesInVertex.append( (self.triangles == i).any(1).nonzero()[0])
            
        for i in range(self.triangles.shape[0]):
            aux = concatenate( (intersect1d( listTrianglesInVertex[self.triangles[i,0]], listTrianglesInVertex[self.triangles[i,1]]), intersect1d( listTrianglesInVertex[self.triangles[i,1]], listTrianglesInVertex[self.triangles[i,2]]), intersect1d( listTrianglesInVertex[self.triangles[i,2]], listTrianglesInVertex[self.triangles[i,0]])), 0)
            aux = aux[(aux!=i).nonzero()[0]]
            self.neighbors.append(aux)
        
    def CalculateMeshes(self, neighbors):
        puntos_revisados = ones(len(neighbors), 'int')
        n = 0
        meshes = []
        Nmeshes = 0
        while 1:
            faltan = puntos_revisados.nonzero()[0]
            if faltan.shape[0] == 0:
                break
            else:
                n += 1
                #puntos_para_revisar=array(faltan[0])
                puntos_para_revisar = ones(faltan.shape, 'int') * -1        
                puntos_para_revisar[0] = faltan[0]
                revisar_N = 0;
                agregados_revision = ones(len(neighbors), 'int')
                agregar_revision_N = 1
                meshes.append(zeros(faltan.shape))
                Nmeshes += 1
                n_meshesN = 0;
                while 1:
                    if revisar_N == agregar_revision_N:
                        break
                    if puntos_revisados[puntos_para_revisar[revisar_N]] == 1:
                        meshes[-1][n_meshesN] = puntos_para_revisar[revisar_N]
                        n_meshesN += 1
                        puntos_revisados[puntos_para_revisar[revisar_N]]=0
                        for i in range(len(neighbors[puntos_para_revisar[revisar_N]])):
                            indice_rev = neighbors[puntos_para_revisar[revisar_N]][i]
                            if indice_rev != -1:
                                if puntos_revisados[indice_rev] == 1:
                                    if agregados_revision[indice_rev] ==1:
                                        puntos_para_revisar[agregar_revision_N] = indice_rev;
                                        agregar_revision_N += 1
                                        agregados_revision[indice_rev]=0
                    revisar_N += 1
                meshes[-1] = meshes[-1][:n_meshesN].astype('int')
        return meshes
                
    def SearchMeshes(self):
        if not hasattr(self, 'neighbors'):
            self.SearchNeighbors()
        self.meshes = self.CalculateMeshes(self.neighbors)
    def SearchMeshesUp3(self):
        if not hasattr(self, 'neighbors'):
            self.SearchNeighborsUp3()
        self.meshes = self.CalculateMeshesUp3(self.neighbors)        
        
    def CalculateMeshesUp3(self, neighbors):
        if not hasattr(self, 'neighbors'):
            self.SearchNeighborsUp3()
        
        puntos_revisados = ones(len(neighbors), 'int')
        n = 0
        meshes = []
        Nmeshes = 0
        while 1:
            faltan = puntos_revisados.nonzero()[0]
            if faltan.shape[0] == 0:
                break
            else:
                n += 1
                #puntos_para_revisar=array(faltan[0])
                puntos_para_revisar = ones(faltan.shape, 'int') * -1        
                puntos_para_revisar[0] = faltan[0]
                revisar_N = 0;
                agregados_revision = ones(len(neighbors), 'int')
                agregar_revision_N = 1
                meshes.append(zeros(faltan.shape))
                Nmeshes += 1
##                print Nmeshes
                n_meshesN = 0;
                while 1:
                    if revisar_N == agregar_revision_N:
                        break
                    if puntos_revisados[puntos_para_revisar[revisar_N]] == 1:
                        meshes[-1][n_meshesN] = puntos_para_revisar[revisar_N]
                        n_meshesN += 1
                        puntos_revisados[puntos_para_revisar[revisar_N]]=0
                        for i in range(len(neighbors[puntos_para_revisar[revisar_N]])):
                            indice_rev = neighbors[puntos_para_revisar[revisar_N]][i]
                            if indice_rev != -1:
                                if puntos_revisados[indice_rev] == 1:
                                    if agregados_revision[indice_rev] ==1:
                                        puntos_para_revisar[agregar_revision_N] = indice_rev;
                                        agregar_revision_N += 1
                                        agregados_revision[indice_rev]=0
                    revisar_N += 1
                meshes[-1] = meshes[-1][:n_meshesN].astype('int')
        return meshes

    def SearchEdges(self):
         # buscar todas las aristas de la malla    
        self.edges = zeros((self.triangles.shape[0]*3,2), 'int')
        self.edgesToTriangles = []
        self.trianglesToEdges = zeros((self.triangles.shape[0],3), 'int')
        edgesAdd_aux = zeros((self.points.shape[0],self.points.shape[0]), 'int') - 1
        edges_count = 0;
        
        for i in range(self.triangles.shape[0]):
            if edgesAdd_aux[self.triangles[i,0], self.triangles[i,1]]==-1:
                self.edges[edges_count] = self.triangles[i,0:2]
                self.edgesToTriangles.append([i])                
                self.trianglesToEdges[i,0] = edges_count
                edgesAdd_aux[self.triangles[i,1], self.triangles[i,0]] = edges_count
                edges_count += 1
            else:
                self.edgesToTriangles[edgesAdd_aux[self.triangles[i,0], self.triangles[i,1]]].extend([i])
                self.trianglesToEdges[i,0] = edgesAdd_aux[self.triangles[i,0], self.triangles[i,1]]
            if edgesAdd_aux[self.triangles[i,1], self.triangles[i,2]] == -1:
                self.edges[edges_count] = self.triangles[i,1:3]
                self.edgesToTriangles.append([i])
                self.trianglesToEdges[i,1] = edges_count
                edgesAdd_aux[self.triangles[i,2], self.triangles[i,1]] = edges_count
                edges_count += 1
            else:
                self.edgesToTriangles[edgesAdd_aux[self.triangles[i,1], self.triangles[i,2]]].extend([i])
                self.trianglesToEdges[i,1] = edgesAdd_aux[self.triangles[i,1], self.triangles[i,2]]
            if edgesAdd_aux[self.triangles[i,2], self.triangles[i,0]]==-1:
                self.edges[edges_count] = self.triangles[i,[2,0]]
                self.edgesToTriangles.append([i])
                self.trianglesToEdges[i,2] = edges_count
                edgesAdd_aux[self.triangles[i,0], self.triangles[i,2]] = edges_count
                edges_count += 1
            else:
                self.edgesToTriangles[edgesAdd_aux[self.triangles[i,2], self.triangles[i,0]]].extend([i])
                self.trianglesToEdges[i,2] = edgesAdd_aux[self.triangles[i,2], self.triangles[i,0]]
        self.edges = self.edges[:edges_count,:]        
        del edgesAdd_aux
        del edges_count        

    def SearchEdgesUp3(self):
         # buscar todas las aristas de la malla, para mallas con mas de 3 vecinos
        self.edges = zeros((self.triangles.shape[0]*3,2), 'int')
        self.edgesToTriangles = []
        self.trianglesToEdges = zeros((self.triangles.shape[0],3), 'int')
##        edgesAdd_aux = zeros((self.points.shape[0],self.points.shape[0]), 'int') - 1
        edgesAdd_aux = sparse.lil_matrix((self.points.shape[0], self.points.shape[0]), dtype='float32') # el menor tipo que acepta lil_matrix es float32
        edges_count = 0;
        
        for i in range(self.triangles.shape[0]):
##            if edgesAdd_aux[self.triangles[i,0], self.triangles[i,1]]==-1:
            if edgesAdd_aux[self.triangles[i,0], self.triangles[i,1]] == 0:
                self.edges[edges_count] = self.triangles[i,0:2]
                self.edgesToTriangles.append([i])                
                self.trianglesToEdges[i,0] = edges_count 
##                edgesAdd_aux[self.triangles[i,1], self.triangles[i,0]] = edges_count
##                edgesAdd_aux[self.triangles[i,0], self.triangles[i,1]] = edges_count
                edges_count_aux = edges_count + 1
                edgesAdd_aux[self.triangles[i,1], self.triangles[i,0]] = edges_count_aux
                edgesAdd_aux[self.triangles[i,0], self.triangles[i,1]] = edges_count_aux

                edges_count += 1
            else:
                edges_count_aux = int(edgesAdd_aux[self.triangles[i,0], self.triangles[i,1]] - 1)
                self.edgesToTriangles[ edges_count_aux ].extend([i])
                self.trianglesToEdges[i,0] = edges_count_aux
            if edgesAdd_aux[self.triangles[i,1], self.triangles[i,2]] == 0:
                self.edges[edges_count] = self.triangles[i,1:3]
                self.edgesToTriangles.append([i])
                self.trianglesToEdges[i,1] = edges_count                
##                edgesAdd_aux[self.triangles[i,2], self.triangles[i,1]] = edges_count
##                edgesAdd_aux[self.triangles[i,1], self.triangles[i,2]] = edges_count
                edges_count_aux = edges_count + 1
                edgesAdd_aux[self.triangles[i,2], self.triangles[i,1]] = edges_count_aux
                edgesAdd_aux[self.triangles[i,1], self.triangles[i,2]] = edges_count_aux
                edges_count += 1
            else:
                edges_count_aux = int(edgesAdd_aux[self.triangles[i,1], self.triangles[i,2]] - 1)
                self.edgesToTriangles[edges_count_aux].extend([i])
                self.trianglesToEdges[i,1] = edges_count_aux
            if edgesAdd_aux[self.triangles[i,2], self.triangles[i,0]] == 0:
                self.edges[edges_count] = self.triangles[i,[2,0]]
                self.edgesToTriangles.append([i])
                self.trianglesToEdges[i,2] = edges_count
##                edgesAdd_aux[self.triangles[i,0], self.triangles[i,2]] = edges_count
##                edgesAdd_aux[self.triangles[i,2], self.triangles[i,0]] = edges_count
                edges_count_aux = edges_count + 1
                edgesAdd_aux[self.triangles[i,0], self.triangles[i,2]] = edges_count_aux
                edgesAdd_aux[self.triangles[i,2], self.triangles[i,0]] = edges_count_aux
                edges_count += 1
            else:
                edges_count_aux = int(edgesAdd_aux[self.triangles[i,2], self.triangles[i,0]] - 1)
                self.edgesToTriangles[edges_count_aux].extend([i])
                self.trianglesToEdges[i,2] = edges_count_aux
        self.edges = self.edges[:edges_count,:]        
        del edgesAdd_aux
        del edges_count         

    def ComputeNormalsOfTriangles(self):
        self.normalsOfTriangles = cross((self.points[self.triangles[:,1],:].astype('float64') - self.points[self.triangles[:,0],:].astype('float64')), (self.points[self.triangles[:,2],:].astype('float64') - self.points[self.triangles[:,0],:].astype('float64')))
##        self.normalsOfTriangles = cross(self.points[self.triangles[:,0],:].astype('float64'), self.points[self.triangles[:,1],:].astype('float64')) + cross(self.points[self.triangles[:,1],:].astype('float64'), self.points[self.triangles[:,2],:].astype('float64')) + cross(self.points[self.triangles[:,2],:].astype('float64'), self.points[self.triangles[:,0],:].astype('float64'))
        self.normalsOfTriangles /= sqrt((self.normalsOfTriangles*self.normalsOfTriangles).sum(1)).reshape(-1,1)  

    def ComputeNormalsOfPoints(self):
        PtoVTK = PythonToPolyData()
        PtoVTK.SetInputTriangles(self)
        PtoVTKOut = PtoVTK.GetVTKPolyData()
        trianglesMeshPolyData = PtoVTKOut['vtkData']
        getNormalsPolyData = vtkPolyDataNormals()
        getNormalsPolyData.SetInput(trianglesMeshPolyData)
        getNormalsPolyData.SplittingOff()
        
        normalsPolyData = getNormalsPolyData.GetOutput()
        normalsPolyData.Update()
        
        vtkFloatArrayNormals = normalsPolyData.GetPointData().GetNormals()
        arrayNormals = zeros( ( vtkFloatArrayNormals.GetNumberOfTuples(), vtkFloatArrayNormals.GetNumberOfComponents() ),  'f')
        vtkFloatArrayNormals.ExportToVoidPointer(arrayNormals)
        arrayNormals_bak = arrayNormals.copy()
        arrayNormals[:,0] = arrayNormals_bak[:,2]
        arrayNormals[:,2] = arrayNormals_bak[:,0]
        
        self.normalsOfPoints = arrayNormals

    def ComputeMeanCurvatureInPoints(self):
        "Calcula la curvatura media en cada punto"
        PtoVTK = PythonToPolyData()
        PtoVTK.SetInputTriangles(self)
        vtkTrianglesMeshDic = PtoVTK.GetVTKPolyData()
        vtkTrianglesMesh = vtkTrianglesMeshDic['vtkData']
        vtk_Curvatures = vtkCurvatures()
        vtk_Curvatures.SetCurvatureTypeToMean()
        vtk_Curvatures.SetInput(vtkTrianglesMesh)
        vtk_Curvatures.Update()
        vtk_CurvaturesPolyData = vtk_Curvatures.GetOutput()

        vtkFloatArrayScalars = vtk_CurvaturesPolyData.GetPointData().GetScalars()
        arrayScalars = zeros((vtkFloatArrayScalars.GetNumberOfTuples(),vtkFloatArrayScalars.GetNumberOfComponents()), 'float64')
        vtkFloatArrayScalars.ExportToVoidPointer(arrayScalars)

        self.meanCurvatureInPoints = arrayScalars
        
        
    def ComputeAreaOfTriangles(self):
        w = self.points[self.triangles[:,1]] - self.points[self.triangles[:,0]]
        v = self.points[self.triangles[:,2]] - self.points[self.triangles[:,0]]
        c = cross(w, v)
        self.areaOfTriangles = sqrt((c * c).sum(1)) * 0.5
        
    def TrianglesToSimplex4(self, tangentPlane = 'points'):
        """
        ENTRDAS
        pointsTriangles, triangles, neighbors, meshes
        
        SALIDAS
        simplex = array(N x 3) ; N:numero de puntos 
        facesSimplex = list(N) ; N:numero de tipos de caras, osea con 1,2,3 ... o N puntos
        facesSimplex[i] = array(N x i) ; N:numero de caras con i puntos, osea es un arreglo con todas las caras con i puntos
        pointFace = array(N x 3 x 2) ; N:numero de puntos; 
        pointFace[punto_simplex, Nº puntos en la face(para buscarlo en facesSimplex), posicion en la lista (dentro de la lista con N puntos)]
        meshes = list(N) , N:numero de mallas
        meshes[i] = array(Np); Np:numero de puntos en la malla i ; se crea denuevo porque se agregan los puntos de borde de la malla simplex
        meshesContour = list(N); N:numero de mallas
        meshesContour[i] = list(Nc); Nc:numero de contornos en la malla i
        meshesContour[i][j] = array(Np) ; Np:numero de puntos en el contorno j de la malla i
        facesSimplexContour = array(Ncx2); Nc: numero de caras que pertenecen al contorno
        facesSimplexContour[numero de puntos i de la cara, lugar dentro de la lista (lista facesSimplex) con i puntos]
        neighbors = array(N x 3); N:numero puntos; se crea denuevo porque se agregan los puntos del borde de la malla simplex

        Ocupa una estimacion del punto central (estimacion por tangentes) para calcular los pesos de las distancias a los puntos del triangulo       
        """
        
        #pointsTriangles, triangles, neighbors, meshes
        w = 0.5 / 3.        
        
        if not hasattr(self, 'neighbors'):
            self.SearchNeighbors()
        if not hasattr(self, 'meshes'):
            self.SearchMeshes()        
        if not hasattr(self, 'normalsOfTriangles'):
            self.ComputeNormalsOfTriangles()
            
        TT = timeit.Timer()
        T1 = TT.timer()           
        if not hasattr(self, 'normalsOfPoints'):
            self.ComputeNormalsOfPoints()
        if not hasattr(self, 'areaOfTriangles'):    
            self.ComputeAreaOfTriangles()
        if not hasattr(self, 'pointToTriangle'):    
            self.SearchPointToTriangle()
            
        from numpy import linalg as linalg
        import copy
        meshesSimplex = copy.deepcopy(self.meshes)
        neighborsSimplex = copy.deepcopy(self.neighbors)

        simplex = zeros(self.triangles.shape, 'float32')
        
        # calculo de los planos correspondientes a cada triangulo o punto
        
        if tangentPlane == 'faces':
            planes = concatenate((self.normalsOfTriangles, -(self.points[self.triangles[:,0],:] * self.normalsOfTriangles).sum(1).reshape(-1,1)), 1)
        elif tangentPlane == 'points':
            planes = concatenate((self.normalsOfPoints, -(self.points * self.normalsOfPoints).sum(1).reshape(-1,1)), 1)
        
##        list_beta = [] #BORRAR
        
        T2 = TT.timer()
        if self.flag_printData:
            print 'Tiempo1 TtoSTang'
            print T2-T1
        
        areasGroup = zeros(3, 'float32')
        for i in range(simplex.shape[0]):

##            T_iter = time.time()
#            beta = [1.,1.,1.]
             # calculo de pesos
            if tangentPlane == 'faces':
                trianglesList = set(((self.triangles == self.triangles[i,0]) + (self.triangles == self.triangles[i,1]) + (self.triangles == self.triangles[i,2])).any(1).nonzero()[0])
                trianglesList.discard(i)
                trianglesList = list(trianglesList) # lista de los triangulos que comparten alguno de los puntos del triangulo en cuestion
                
                sumAreas = self.areaOfTriangles[trianglesList].sum()
                alfa = self.areaOfTriangles[trianglesList] / sumAreas # Calculo de los pesos alfa solo por area de triangulo pero sin emplear angulo
            elif tangentPlane == 'points':
                areasGroup[0] = self.areaOfTriangles[self.pointToTriangle[self.triangles[i,0]]].sum()
                areasGroup[1] = self.areaOfTriangles[self.pointToTriangle[self.triangles[i,1]]].sum()
                areasGroup[2] = self.areaOfTriangles[self.pointToTriangle[self.triangles[i,2]]].sum()
                    
##                    areasGroup[trit] = self.areaOfTriangles[trianglesList].sum()
                alfa = areasGroup / areasGroup.sum()
                
##            print 'Tiempo seleccion triangulos'
##            print time.time() - T_iter
##            T_iter = time.time()
             
            triangle_aux = self.triangles[i,:]
            triangle_points = self.points[triangle_aux,:]
            pCentral = triangle_points.sum(0) / 3.
            q = (dot(self.normalsOfPoints[triangle_aux[0]], triangle_points[0] - pCentral) * self.normalsOfPoints[triangle_aux[0]]) + (dot(self.normalsOfPoints[triangle_aux[1]], triangle_points[1] - pCentral) * self.normalsOfPoints[triangle_aux[1]]) + (dot(self.normalsOfPoints[triangle_aux[2]], triangle_points[2] - pCentral) * self.normalsOfPoints[triangle_aux[2]])            
            pointAprox = pCentral + q * w # aproximacion por tangentes al nuevo punto simplex
            pointAprox = concatenate((pointAprox, [1]))
            pointAprox = pointAprox.reshape(-1,1)
            
            A = zeros((1,3), dtype='float32')
            if tangentPlane == 'faces':
                for one_alfa, t in zip(alfa, trianglesList):
                    A  += one_alfa * dot(planes[t,:].reshape(1,-1), pointAprox) * self.normalsOfTriangles[t]
            elif tangentPlane == 'points':
                for one_alfa, t in zip(alfa, triangle_aux):
                    A  += one_alfa * dot(planes[t,:].reshape(1,-1), pointAprox) * self.normalsOfPoints[t]
                
            B = (triangle_points.sum(0).reshape(-1,1) - 3. * pointAprox[:3])
                
            dotB = dot(B.T,B)
            if dotB == 0.:
                beta = 2.
            else:
                beta = (dot(A,B)/dotB)[0,0]
                if beta < 0.1:
                    beta = 0.1
                elif beta > 2.:
                    beta = 2.
        
        
            # aporte de los puntos del triangulo
            Q = zeros((4,4), dtype='float32')
            Q_aux = eye(4, dtype='float32')
            Q_aux[3,3] = 0.
            
            for k in triangle_points:
                Q_aux[3,3] = dot(k, k)
                Q_aux[:3,3] = -k
                Q_aux[3,:3] = -k
                Q += beta * Q_aux
            
            # aporte de los planos de los triangulos vecinos
            if tangentPlane == 'faces':
                for one_alfa, t in zip(alfa, trianglesList):
                    Q  += one_alfa * dot(planes[t,:].reshape(-1,1), planes[t,:].reshape(1,-1))
            elif tangentPlane == 'points':
                for one_alfa, t in zip(alfa, triangle_aux):
                    Q  += one_alfa * dot(planes[t,:].reshape(-1,1), planes[t,:].reshape(1,-1))                
            
            simplex[i,:]=(dot(linalg.inv(Q[:3,:3]),-Q[:3,3].reshape(-1,1))).flatten()
            
##            print 'Tiempo resto'
##            print time.time() - T_iter
##            T_iter = time.time()
            
        T2 = TT.timer()
        if self.flag_printData:
            print 'Tiempo2 TtoSTang'
            print T2-T1

                
        sCont = simplex.shape[0]
        contoursPointsSimplex = [] # solo se llena si se encuentran contornos, en caso contrario queda como una lista bacia
        
        # encontrar si hay bordes
        nMod3=[0,1,2,0,1,2,0,1,2]
        pointTriangle_pointEdge_pointTriangle = []
        triangle2pointEdge = []
        trianglesEdge = (neighborsSimplex == -1).any(1).nonzero()[0] # lista de triangulos que estan en borde (los vecinos -1 son que no tiene => hay un borde)
        if len(trianglesEdge)>0: # hay bordes        
            simplexEdges = zeros((2*len(trianglesEdge), 3), 'float32')
            neighborsEdges = zeros((2*len(trianglesEdge), 3), 'int')
            neighborsEdgesAux = zeros((2*len(trianglesEdge), 1), 'int')
            
            triangle2pointEdge = zeros((2*len(trianglesEdge), 2), 'int')
            pointTriangle_pointEdge_pointTriangle = zeros((2*len(trianglesEdge), 3), 'int')
            triangleCent = zeros((3), 'int')
            edge = zeros((2), 'int')
            pCont = 0
            for i in trianglesEdge:
                triangleCent[:] = self.triangles[i,:]
                triangulosComp = neighborsSimplex[i,:].copy()
                triangulosComp = triangulosComp[(triangulosComp!=-1).nonzero()] # dejo los triangulos vecinos eliminando el -1 (borde)
                if len(triangulosComp) == 2:
                    p1 = intersect1d_nu(triangleCent, self.triangles[triangulosComp[0],:]) # los 2 puntos comunes del primer triangulo
                    p2 = intersect1d_nu(triangleCent, self.triangles[triangulosComp[1],:]) # los 2 puntos comunes del segundo triangulo
                    pC = intersect1d_nu(p1, p2) # el punto comun de los 2 triangulos vecinos al triangulo del borde
                    edge[:] = triangleCent[(triangleCent!=pC).nonzero()] # los 2 puntos del triangulo central que dan al borde
                    simplexEdges[pCont,:] = (self.points[edge,:].sum(0) / 2.)
                    aux = triangleCent; aux[(aux==pC).nonzero()[0][0]] = -1  # ocupare estos puntos para, al final, unirlos en el orden correcto con los vecinos
                    neighborsEdges[pCont,:] = aux
                    neighborsEdgesAux[pCont] = i 
                    triangle2pointEdge[pCont,:] = [i,sCont] # referencia para MakeFace
                    pointTriangle_pointEdge_pointTriangle[pCont,:] = [edge[0], sCont, edge[1]]  # referencia para MakeFace (no esta ordenado con respecto a la normal)
                    pCont += 1
                    neighborsSimplex[i,(neighborsSimplex[i,:] == -1).nonzero()[0]] = sCont # le pongo el punto vecino que faltaba                
                    sCont += 1
                else:
                    p1 = intersect1d_nu(triangleCent, self.triangles[triangulosComp,:]) # los 2 puntos comunes con el traingulo vecino 
                    ipC = ((triangleCent!=p1[0]) * (triangleCent!=p1[1])).nonzero()[0][0]
                    pC = triangleCent[ipC]
                    iNeighbor = (neighborsSimplex[i,:]!=-1).nonzero()[0][0]
                    for k in [1, -1]:
                        simplexEdges[pCont,:] = (self.points[[pC,triangleCent[nMod3[ipC+k]]],:].sum(0) / 2.) #
                        
                        aux = triangleCent.copy(); aux[nMod3[ipC-k]] = -1  # ocupare estos puntos para, al final, unirlos en el orden correcto con los vecinos
                        neighborsEdges[pCont,:] = aux
                        
                        neighborsEdgesAux[pCont] = i
                        triangle2pointEdge[pCont,:] = [i,sCont] # referencia para MakeFace
                        pointTriangle_pointEdge_pointTriangle[pCont,:] = [pC, sCont, triangleCent[nMod3[ipC+k]]]  # referencia para MakeFace 
                        pCont += 1
                        neighborsSimplex[i, nMod3[iNeighbor-k] ] = sCont
                        sCont += 1
                        
            simplexEdges = simplexEdges[:pCont,:]
            neighborsEdges = neighborsEdges[:pCont,:]
            neighborsEdgesAux = neighborsEdgesAux[:pCont,:]
            pointTriangle_pointEdge_pointTriangle = pointTriangle_pointEdge_pointTriangle[:pCont,:]
            triangle2pointEdge = triangle2pointEdge[:pCont,:]
               
            aux = neighborsEdges[(neighborsEdges != -1).nonzero()].reshape(-1,2)
            aux2 = zeros(aux.shape, 'int')
            for i in range(pCont): # busco los puntos del borde asociados con los mismo puntos de la triangulacion, para asociarlos como vecinos
                for k in range(2):
                    indice = (aux == aux[i,k]).nonzero()
                    aux2[indice] = i; aux2[i,k] = indice[0][0]
                
            neighborsEdges[(neighborsEdges!=-1).nonzero()] = aux2.flatten() + neighborsSimplex.shape[0]  # los puntos simplex de borde vecinos
                
            # crear las listas con los bodes (contornos)
            EdgesAdded = ones(pCont, 'int') # para marcar los puntos que ya se agregaron a un contorno
            contour = []
            while 1:
                try:
                    pIni = EdgesAdded.nonzero()[0][0] # ver si queda algun punto sin ser agregado
                except:
                    break
                p = pIni
                contour.append([])
                contour[-1] = zeros(pCont, 'int')
                contour[-1][0] = p
                EdgesAdded[p] = 0
                contPContour = 1
                p = neighborsEdges[p, nMod3[(neighborsEdges[p,:]==-1).nonzero()[0][0] + 1]] - neighborsSimplex.shape[0]
                while p != pIni :
                    contour[-1][contPContour] = p
                    EdgesAdded[p] = 0
                    contPContour += 1
                    p = neighborsEdges[p, nMod3[(neighborsEdges[p,:]==-1).nonzero()[0][0] + 1]] - neighborsSimplex.shape[0]
                contour[-1] = contour[-1][:contPContour]
                contour[-1] = contour[-1] + neighborsSimplex.shape[0]
                
            # agregar los contornos a la malla (mesh) correspondiente
            for i in range(len(contour)):
                contoursPointsSimplex.append(contour[i][::-1])
                for j in range(len(meshesSimplex)):
                    if len((meshesSimplex[j] == neighborsEdgesAux[contour[i][0]- neighborsSimplex.shape[0]]).nonzero()[0]) > 0:
                        meshesSimplex[j] = concatenate((meshesSimplex[j], contour[i]), 0)
                        break
            # ## 
            neighborsEdges[(neighborsEdges==-1).nonzero()] = neighborsEdgesAux.flatten() # los puntos vecinos dentro de la malla
            
            
            simplex = concatenate((simplex,simplexEdges),0)
            neighborsSimplex = concatenate((neighborsSimplex, neighborsEdges),0)
        
        pointFace = ones((simplex.shape[0],3,2), 'int') * -1 # pointFace[punto_simplex, Nº puntos en la face(para buscarlo en facesSimplex), posicion en la lista (dentro de la lista con N puntos)]

        N_parches_inicial=20
        N_parches_por_parche_inicial=100;
        facesSimplex=[]
        for i in range(N_parches_inicial): facesSimplex.append(zeros((N_parches_por_parche_inicial,i), 'int'))
        facesSimplex_cont=zeros(N_parches_inicial)

        # ## para prueba
        facesSimplexB=[]
        for i in range(N_parches_inicial): facesSimplexB.append(zeros((N_parches_por_parche_inicial,i), 'int'))
        facesSimplexB_cont=zeros(N_parches_inicial)
        # #######

        facesSimplexContour = zeros((N_parches_por_parche_inicial, 2), 'int') # [numero de puntos i de la cara, lugar dentro de la lista con i puntos]
        facesSimplexContour_cont = 0
        
        cont=0

        for i in range(self.points.shape[0]):

            [face, FaceBorde] = self.__MakeFace__(i, neighborsSimplex, pointTriangle_pointEdge_pointTriangle, triangle2pointEdge)            
            L = len(face)
            
            facesSimplex[L][facesSimplex_cont[L],:]=face
            if FaceBorde:
                # ## para prueba
                facesSimplexB[L][facesSimplexB_cont[L],:]=face
                # #######
                facesSimplexContour[facesSimplexContour_cont,:] = [L, facesSimplexB_cont[L]]
                facesSimplexContour_cont += 1
            
            for j in range(L):
                try:
                    if self.triangles[face[j],0]==i:
                        pointFace[face[j],0,:]=[L, facesSimplex_cont[L]]
                    elif self.triangles[face[j],1]==i:
                        pointFace[face[j],1,:]=[L, facesSimplex_cont[L]]
                    elif self.triangles[face[j],2]==i:
                        pointFace[face[j],2,:]=[L, facesSimplex_cont[L]]
                except: # para los puntos del borde
                    aux = pointTriangle_pointEdge_pointTriangle[(pointTriangle_pointEdge_pointTriangle[:,1]==face[j]).nonzero()[0], [0,2]]
                    if aux[0] == i:
                        pointFace[face[j],0,:]=[L, facesSimplex_cont[L]]
                    elif aux[1] == i:
                        pointFace[face[j],1,:]=[L, facesSimplex_cont[L]]
                    
            facesSimplex_cont[L] += 1
            if facesSimplex_cont[L] == facesSimplex[L].shape[0]:
                facesSimplex[L] = concatenate((facesSimplex[L],zeros((N_parches_por_parche_inicial,L), 'int')),0)

            if facesSimplexContour_cont == facesSimplexContour.shape[0]:
                facesSimplexContour = concatenate((facesSimplexContour,zeros((N_parches_por_parche_inicial, 2), 'int')),0)

            # ## para prueba
            facesSimplexB_cont[L] += 1
            if facesSimplexB_cont[L] == facesSimplexB[L].shape[0]:
                facesSimplexB[L] = concatenate((facesSimplexB[L],zeros((N_parches_por_parche_inicial,L), 'int')),0)        
            # ######

            
        L=0
        for i in range(N_parches_inicial):
            if facesSimplex_cont[i]>0:
                facesSimplex[i]=facesSimplex[i][0:facesSimplex_cont[i],:]
                L=i
            else:
                facesSimplex[i]=array([])
        facesSimplex = facesSimplex[0:L+1]

        facesSimplexContour = facesSimplexContour[:facesSimplexContour_cont,:]

        # ## para prueba
        L=0
        for i in range(N_parches_inicial):
            if facesSimplexB_cont[i]>0:
                facesSimplexB[i]=facesSimplexB[i][0:facesSimplexB_cont[i],:]
                L=i
            else:
                facesSimplexB[i]=array([])
        facesSimplexB = facesSimplexB[0:L+1]
        # ##########
        
        OUT = SimplexMesh(simplex, neighborsSimplex, facesSimplex, pointFace, meshesSimplex, contoursPointsSimplex, facesSimplexContour)
    
##        out=[simplex, self.neighbors, facesSimplex, pointFace, self.meshes, meshesContour ,facesSimplexB, facesSimplexContour]  
        return OUT
    
    def __MakeFace__(self, punto_central, neighbors, pointTriangle_pointEdge_pointTriangle, triangle2pointEdge):
        """
        Version 2
        """
        
        triangulo_inicial = (any(self.triangles==punto_central,1)).nonzero()[0][0]    # elige el primer triangulo que encuentre al que pertenesca el punto
        ipunto_central = (self.triangles[triangulo_inicial,:]==punto_central).nonzero()[0][0] #indice del punto dentral dentro del triangulo
        
        face = ones(20, 'int')*(-1)        
        face[0] = triangulo_inicial
        contP=1

        FaceBorde = 0    

        # cada triangulo corresponde a un punto simplex    
        termino = 0
        triangulo_anterior = triangulo_inicial
        while termino==0:       
            punto_borde=self.triangles[triangulo_inicial,ipunto_central-1]
            encontrado = False
            for triangulos_vecinos in range(3):
                triangulo_en_comprobacion = neighbors[triangulo_inicial,triangulos_vecinos]
                try:
                    puntosTrianguloComp = self.triangles[triangulo_en_comprobacion,:]
                    if any(puntosTrianguloComp == punto_central) and any(puntosTrianguloComp==punto_borde):
                        face[contP] = triangulo_en_comprobacion
                        contP += 1
                        encontrado = True
                        break 
                except: # el vecino es mas grande que la lista de triangulos porque es un punto simplex de borde agregado
                    puntosTrianguloComp = triangulo_inicial # no hacer nada 
            triangulo_anterior = triangulo_inicial
            if encontrado: # encontro un triangulo vecino que corresponde
                triangulo_inicial = triangulo_en_comprobacion
                ipunto_central = (self.triangles[triangulo_inicial,:]==punto_central).nonzero()[0][0]
            else: # no lo encontro porque es un borde
                FaceBorde = 1
                aux = triangle2pointEdge[(triangle2pointEdge[:,0]==triangulo_inicial).nonzero(),1]
                if aux.shape[1]>1: # si el triangulo tiene 2 caras que son bordes, osea hay 2 puntos de borde(simplex) asociados al triangulo
                    l1 = pointTriangle_pointEdge_pointTriangle[(pointTriangle_pointEdge_pointTriangle[:,1]==aux[0,0]).nonzero(), [0, 2]]
                    l2 = pointTriangle_pointEdge_pointTriangle[(pointTriangle_pointEdge_pointTriangle[:,1]==aux[0,1]).nonzero(), [0, 2]]
                    v1 = (l1==punto_central).any()
                    v2 = (l2==punto_central).any()
                    if v1*v2: # 
                        ##pm = intersect1d_nu(l1,l2) # vertice del triangulo entre las 2 caras que son bordes
                        pm = punto_central
                        pm0 = self.triangles[triangulo_inicial,:][(self.triangles[triangulo_inicial,:]==pm).nonzero()[0]-1] #punto antes de pm en el triangulo
                        if (l1==pm0).any():
                            aux = aux[0,0]
                        else:
                            aux = aux[0,1]
                    elif v1:
                        aux = aux[0,0]
                    elif v2:
                        aux = aux[0,1]
                face[contP] = aux # punto simplex de borde
                aux = pointTriangle_pointEdge_pointTriangle[concatenate(((pointTriangle_pointEdge_pointTriangle[:,0] == punto_central).reshape(-1,1), (pointTriangle_pointEdge_pointTriangle[:,2] == punto_central).reshape(-1,1)),1).any(1).nonzero()[0],:]
                aux = aux[ (aux[:,1]!=face[contP]).nonzero()[0],: ][0]
                contP += 1
                face[contP] = aux[1] # el siguiente punto simplex de borde
                contP += 1
                face[contP] = triangle2pointEdge[(triangle2pointEdge[:,1]==aux[1]).nonzero()[0][0],0] # el el punto en el centro del triangulo que es borde
                triangulo_inicial = face[contP]
                contP += 1 
                ipunto_central = (self.triangles[triangulo_inicial,:]==punto_central).nonzero()[0][0]
                
            if triangulo_inicial == face[0]: termino = 1
        face = face[0:contP-1]
        return face, FaceBorde

    def DelMeshes(self, meshesForDelete, flagIndex = False):
        noDelTriangles = ones(self.triangles.shape[0], 'int')
        for delTriangles in meshesForDelete:
            noDelTriangles[self.meshes[delTriangles]] = 0
            
        # borrar los triangulos
        tianglesNew = self.triangles[noDelTriangles.nonzero()[0]]

        # encontrar los puntos que ya no estan referenciados
        pointsIndex = arange(self.points.shape[0])
        noDelPoints = setmember1d(pointsIndex, tianglesNew.flatten())
        self.points = self.points[noDelPoints.nonzero()[0]]
        pointsIndex = pointsIndex[noDelPoints.nonzero()[0]]
                             
        # actualizar los indices a los puntos en la lista de triangulos
        trianglesIndex = arange(self.triangles.shape[0]) # tomar los indices antes de borrarlo
        self.triangles = zeros(tianglesNew.shape, 'int')
        for i in range(pointsIndex.shape[0]):
            self.triangles[(tianglesNew==pointsIndex[i]).nonzero()] = i
            
        # actualizar las mallas y sus indices a los triangulos en las mallas        
        trianglesIndex = trianglesIndex[noDelTriangles.nonzero()[0]]
        meshesNew = []
        meshesNewIndex = range(len(self.meshes))
        for i in meshesForDelete:
            meshesNewIndex.remove(i)
        for i in meshesNewIndex:
            oneMeshNew = zeros(self.meshes[i].shape[0], 'int')
            for j in range(trianglesIndex.shape[0]):
                oneMeshNew[(self.meshes[i]==trianglesIndex[j]).nonzero()] = j
            meshesNew.append(oneMeshNew)
                
        self.meshes = meshesNew
        if flagIndex:
            self.indexDelPoints = pointsIndex
            self.indexDeltriangles = trianglesIndex
    def DeleteNotUsedPoints(self, flagIndex=False):

        # encontrar los puntos que ya no estan referenciados
        pointsIndex = arange(self.points.shape[0])
        noDelPoints = setmember1d(pointsIndex, self.triangles.flatten())
        self.points = self.points[noDelPoints.nonzero()[0]]
        pointsIndex = pointsIndex[noDelPoints.nonzero()[0]]
                             
        # actualizar los indices a los puntos en la lista de triangulos
        trianglesIndex = arange(self.triangles.shape[0]) # tomar los indices antes de borrarlo
        tianglesNew = zeros(self.triangles.shape, 'int')
        for i in range(pointsIndex.shape[0]):
            tianglesNew[(self.triangles==pointsIndex[i]).nonzero()] = i
        self.triangles = tianglesNew
        del tianglesNew
                        
        if flagIndex:
            self.indexDelPoints = pointsIndex
            
    def DeleteTrianglesOf2Edges(self, method_flag = 'juntos'):
        """
        Elimina los triangulos que tiene 2 aristas en el borde de la malla (sin vecinos )
        dividiendolos en 2 triangulos, para esto tambien divide el unico vecino de ese triangulo
        y debe agregar un punto en la arista comun de ambos. Por esto entrega una nueva lista de
        puntos y vecinos. Las entradas se modifican y se entregan en la salida. Para que no cambien
        hay que entrar una copia *.copy()
        """
        
        if method_flag == 'deAuno': # calculo de a uno, MUY INEFICIENTE PERO PRESENTA MENOS PROBLEMAS
            self.SearchNeighbors()
            for i in range(self.triangles.shape[0]):
                TrianglesWith2Edges = ((self.neighbors == -1).sum(1) == 2).nonzero()[0]
                if len(TrianglesWith2Edges) > 0:
                    one_triangle = TrianglesWith2Edges[0]
                    indexNeighborsBack = (self.neighbors[one_triangle,:] != -1).nonzero()[0][0]
                    neighborsBack = self.neighbors[one_triangle,:][indexNeighborsBack]
                    indexPointBack = (self.neighbors[neighborsBack,:] == one_triangle).nonzero()[0][0] - 1
                    pointBack = self.triangles[neighborsBack, indexPointBack] # el punto del tr1angulo vecino al con 2 bordes, que no es un vertice del triangulo con 2 bordes
                    
                    pointsT2 = self.triangles[one_triangle, [indexNeighborsBack -2, indexNeighborsBack -1, indexNeighborsBack]].T
                    pointCentNew = (self.points[pointsT2[0],:] + self.points[pointsT2[2],:]) / 2.
                    indexpointCentNew = self.points.shape[0] # nuevo punto en la arista entre el triangulo con 2 bordes y su vecino
                    self.points = concatenate((self.points, pointCentNew.reshape(-1,3)) ,0)
                    
                    TrianglesNew1 = array([indexpointCentNew, pointsT2[2:3][0], pointBack])
                    TrianglesNew2 = array([indexpointCentNew, pointBack, pointsT2[0:1][0]])
                    TrianglesNew3 = array([indexpointCentNew, pointsT2[:2][0],pointsT2[:2][1]]).reshape(-1,3)
                    TrianglesNew4 = array([indexpointCentNew, pointsT2[1:][0],pointsT2[1:][1]]).reshape(-1,3)

                    # agregar los nuevos triangulos
                    self.triangles[one_triangle,:] = TrianglesNew1
                    self.triangles[neighborsBack,:] = TrianglesNew2
                    self.triangles = concatenate((self.triangles, TrianglesNew3, TrianglesNew4), 0)

                    self.SearchNeighbors()

                else:
                    return 1
        elif method_flag == 'juntos': #todo junto
            
            if not hasattr(self, 'neighbors'):
                self.SearchNeighbors()
            
            TrianglesWith2Edges = ((self.neighbors == -1).sum(1) == 2).nonzero()[0]
            Ntriangles2Edges = TrianglesWith2Edges.shape[0]
            indexNeighborsBack = (self.neighbors[TrianglesWith2Edges,:] != -1).nonzero()
            neighborsBack = self.neighbors[TrianglesWith2Edges,:][indexNeighborsBack] #triangulos vecinos del que tiene 2 aristas de borde

            indexPointBack = (self.neighbors[neighborsBack,:] == TrianglesWith2Edges.reshape(-1,1)).nonzero()[1] - 1
            
            pointBack = self.triangles[neighborsBack, indexPointBack ].reshape(-1,1) # el punto del truangulo vecino al con 2 bordes, que no es un vertice del triangulo con 2 bordes
            pointsT2 = self.triangles[TrianglesWith2Edges, [indexNeighborsBack[1]-2, indexNeighborsBack[1]-1, indexNeighborsBack[1]]].T
            pointCentNew = (self.points[pointsT2[:,0],:] + self.points[pointsT2[:,2],:]) / 2.
            indexpointCentNew = arange(self.points.shape[0], self.points.shape[0]+Ntriangles2Edges).reshape(-1,1) # nuevo punto en la arista entre el triangulo con 3 bordes y su vecino
            self.points = concatenate((self.points, pointCentNew) ,0) # agregar en nuevo punto

            # nuevos triangulos
            TrianglesNew1 = concatenate((indexpointCentNew, pointsT2[:,2:3], pointBack) ,1)
            TrianglesNew2 = concatenate((indexpointCentNew, pointBack, pointsT2[:,0:1]) ,1)
            TrianglesNew3 = concatenate((indexpointCentNew, pointsT2[:,:2]) ,1)
            TrianglesNew4 = concatenate((indexpointCentNew, pointsT2[:,1:]) ,1)

            NtrianglesOld = self.triangles.shape[0]

            # agregar los nuevos triangulos
            self.triangles[TrianglesWith2Edges,:] = TrianglesNew1
            self.triangles[neighborsBack,:] = TrianglesNew2
            self.triangles = concatenate((self.triangles, TrianglesNew3, TrianglesNew4), 0)

            indexTrianglesNew3 = arange(NtrianglesOld, NtrianglesOld+Ntriangles2Edges)
            indexTrianglesNew4 = arange(NtrianglesOld+Ntriangles2Edges, NtrianglesOld+Ntriangles2Edges*2)

            neighbors_aux = self.neighbors[neighborsBack, [indexPointBack-1, indexPointBack]].T # los vecinos del triangulo "neighborsBack", sin con contar los triangulos que queremos eliminar

            index_aux = (neighbors_aux[:,0] != -1).nonzero()[0]
            index_aux2 = (self.neighbors[neighbors_aux[index_aux, 0], :] == neighborsBack[index_aux].reshape(-1,1)).nonzero()[1]    
            self.neighbors[neighbors_aux[index_aux, 0], index_aux2] = TrianglesWith2Edges[index_aux] # cambiar el vecino de los triangulos que ahora seran vecinos de TrianglesNew1
           
            # vecinos de los nuevos triangulos
            neighborsNew1 = zeros(TrianglesNew1.shape, 'int') - 1
            neighborsNew1[:,0] = indexTrianglesNew4
            neighborsNew1[:,1] = neighbors_aux[:,0]
            neighborsNew1[:,2] = neighborsBack

            neighborsNew2 = zeros(TrianglesNew1.shape, 'int') - 1
            neighborsNew2[:,0] = TrianglesWith2Edges
            neighborsNew2[:,1] = neighbors_aux[:,1]
            neighborsNew2[:,2] = indexTrianglesNew3

            neighborsNew3 = zeros(TrianglesNew1.shape, 'int') - 1
            neighborsNew3[:,0] = neighborsBack
            neighborsNew3[:,2] = indexTrianglesNew4

            neighborsNew4 = zeros(TrianglesNew1.shape, 'int') - 1
            neighborsNew4[:,0] = indexTrianglesNew3
            neighborsNew4[:,2] = TrianglesWith2Edges

            # agregar los nuevos vecinos
            self.neighbors[TrianglesWith2Edges,:] = neighborsNew1
            self.neighbors[neighborsBack,:] = neighborsNew2
            self.neighbors = concatenate((self.neighbors, neighborsNew3, neighborsNew4), 0)
    def SearchExternalContours(self):
        """
        v.1
        Encuentra los bordes de una malla de triangulos
        SALIDAS
        ContoursTriangles = list(Nm) ; Nm:numero de mallas
        ContoursTriangles[i] = list(Nc); Nc:numero de contornos (bordes) dentro de la malla i
        ContoursTriangles[i][j] = array(Nt); Nt:numero de triangulos que forman el contorno j de la mala i. Es un arreglo con
                                            los triangulos que forman el contorno. No se consideran formando contorno los triangulos
                                            que solo tiene un punto en el contorno, es decir solo se consideran los que no tienen 3 vecinos.
        
        ContoursPointsTriangles = list(Nm) ; Nm:numero de mallas
        ContoursPointsTriangles[i] = list(Nc); Nc:numero de contornos (bordes) dentro de la malla i
        ContoursPointsTriangles[i][j] = array(Np); Np:numero de puntos que forman el contorno j de la mala i. Es un arreglo con los vertices que forman el contorno
        
        ContoursPointsTriangles = 
        """
        
        self.contoursTriangles = []
        self.contoursPoints = []
            
        EdgeTriangles = (self.neighbors == -1).any(1).nonzero()[0] # triangulos que forman bordes
        TrianglesAdded = ones(EdgeTriangles.shape, 'int') # para marcar los triangulos que ya se agregaron a un contorno
        
        TrianglesEdgePoints  = [];
        menosUno = array([-1])
        # se crea la lista "TrianglesEdgePoints" compuesta de len(EdgeTriangles) vectores,
        # en donde cada vector i son los puntos del triangulo EdgeTriangles[i] que forman el contorno
        for i in range(len(EdgeTriangles)):
            Nbordes = (self.neighbors[EdgeTriangles[i],:]==-1).sum()
            if Nbordes == 1: #triangulos que tienen 2 puntos en el borde
                # para que los segmentos queden con la normal al lado correcto (la misma orientacion que los triangulos)
                neighborsReal = setdiff1d(self.neighbors[EdgeTriangles[i],:],menosUno)
                p1 = intersect1d_nu(self.triangles[EdgeTriangles[i],:], self.triangles[neighborsReal[0],:]) # los 2 puntos comunes del primer triangulo
                p2 = intersect1d_nu(self.triangles[EdgeTriangles[i],:], self.triangles[neighborsReal[1],:]) # los 2 puntos comunes del segundo triangulo
                pC = intersect1d_nu(p1, p2)
                pCi = (self.triangles[EdgeTriangles[i],:] == pC).nonzero()[0]
                if pCi==1: #para que los segmentos queden con la normal al lado correcto
                    TrianglesEdgePoints.append(self.triangles[EdgeTriangles[i],[2,0]])
                else:
                    TrianglesEdgePoints.append(setdiff1d(self.triangles[EdgeTriangles[i],:],pC))
            else: #triangulos que tienen 3 puntos en el borde
                #para que los segmentos queden con la normal al lado correcto
                neighborsReal = setdiff1d(self.neighbors[EdgeTriangles[i],:],menosUno)
                p1 = intersect1d_nu(self.triangles[EdgeTriangles[i],:], self.triangles[neighborsReal[0],:]) # los 2 puntos comunes del primer triangulo
                pC = setdiff1d(self.triangles[EdgeTriangles[i],:], p1)
                pCi = (self.triangles[EdgeTriangles[i],:] == pC).nonzero()[0]
                TrianglesEdgePoints.append(self.triangles[EdgeTriangles[i],[pCi-1, pCi, pCi-2]].flatten())
        Contours_cont = -1
        # se agrupan los segmentos contenidos en "TrianglesEdgePoints" para crear los contornos
        while 1: # loop para los contornos. N iteraciones => N contornos
            TrianglesFaltantes = TrianglesAdded.nonzero()[0] # lista con los triangulos que falta agregar a un contorno
            if len(TrianglesFaltantes)>0: # significa que quedan triangulos sin agregar
                self.contoursTriangles.append(zeros(len(TrianglesFaltantes), 'int'))
                self.contoursPoints.append(zeros(len(TrianglesFaltantes) * 3, 'int'))
                Contours_cont += 1
                Points_count = 0
                Triangles_count = 1
                # se crea el primer punto del contorno
                self.contoursTriangles[Contours_cont][Triangles_count] = EdgeTriangles[TrianglesFaltantes[0]]
                self.contoursPoints[Contours_cont][Points_count:Points_count+len(TrianglesEdgePoints[TrianglesFaltantes[0]])] = TrianglesEdgePoints[TrianglesFaltantes[0]].copy()
                Points_count += len(TrianglesEdgePoints[TrianglesFaltantes[0]])
                TrianglesAdded[TrianglesFaltantes[0]] = 0
                TrianglesFaltantesList = TrianglesFaltantes.tolist()
                TrianglesFaltantesList.remove(TrianglesFaltantes[0])
                while 1: # loop para formar el contorno agregando ordenadamente todos los triangulos que tengan puntos de contorno comunes
                    NOencontro = 1 # maraca para ver si se encontro un triangulo que tenga punto de contorno comunes con el ultimo agregado
                    for i in range(len(TrianglesFaltantesList)): # buscar en la lista "TrianglesFaltantesList" a ver si hay un triangulo que tenga puntos de contorno comunes con el ultimo agregado al contorno
                        if (TrianglesEdgePoints[TrianglesFaltantesList[i]] == self.contoursPoints[Contours_cont][Points_count-1]).any():
                            self.contoursTriangles[Contours_cont][Triangles_count] = EdgeTriangles[TrianglesFaltantesList[i]]
                            self.contoursPoints[Contours_cont][Points_count:Points_count+len(TrianglesEdgePoints[TrianglesFaltantesList[i]])-1] = TrianglesEdgePoints[TrianglesFaltantesList[i]][1:]
                            Triangles_count += 1
                            Points_count += len(TrianglesEdgePoints[TrianglesFaltantesList[i]]) - 1
                            TrianglesAdded[TrianglesFaltantesList[i]] = 0
                            TrianglesFaltantesList.remove(TrianglesFaltantesList[i])
                            NOencontro = 0 # se marca que se encontro otro triangulo para el contorno
                            break
                    if NOencontro: # no se encontro otro triangulo por lo tanto el contorno se debe haber cerrado. No se confirma que el triangulo que sigue al ultimo agregado es el primero que se puso en el contorno, pero se asume que si no se encontro ninguno debe ser el final
                        break
                self.contoursTriangles[Contours_cont] = self.contoursTriangles[Contours_cont][:Triangles_count] # se eliminar del arreglo el espacio que no se utilizo
                self.contoursPoints[Contours_cont] = self.contoursPoints[Contours_cont][:Points_count-1] # se eliminar del arreglo el espacio que no se utilizo
            else: # no quedan triangulos por agregar. Ya se encontraron todos los contornos de la malla M
                break

##        return [ContoursTriangles, ContoursPointsTriangles]

    def SearchExternalContoursUp3(self):
        '''
        Busca los contornos externos de la malla, pero solo funciona con los vecinos calculados con SearchNeighborsOrderUp3
        '''
        boolTrianglesNeighbors = zeros(self.neighbors.shape)
        for i in arange(self.neighbors.shape[0]):
            for j in [0,1,2]:
                try:
                    if self.neighbors[i,j] == -1:
                        boolTrianglesNeighbors[i,j] = 1
                        break
                except:
                    None
##        boolTrianglesNeighbors = negative((self.neighbors != -1))
        edgeTriangles = boolTrianglesNeighbors.any(1).nonzero()[0] # triangulos que forman bordes
        edgesList = []
        trianglesList = []
        aux_edges = [[0,1],[1,2],[2,0]]
        for one_EdgeTriangles in edgeTriangles:
            for edge in boolTrianglesNeighbors[one_EdgeTriangles].nonzero()[0]:
                edgesList.append(self.triangles[one_EdgeTriangles,aux_edges[edge]])
                trianglesList.append(one_EdgeTriangles)
            
        contour_count = 0
        contoursTriangles = []
        contoursPoints = []
        while len(trianglesList) > 0:
            contoursTriangles.append([])
            contoursPoints.append([])
            contoursTriangles[-1].append(trianglesList.pop())
            contoursPoints[-1].extend(edgesList.pop())

            while 1:
                if len(trianglesList) == 0:
                    break
                edgesListArray = concatenate(edgesList).reshape(-1,2)
                candidats = (edgesListArray[:,0] == contoursPoints[-1][-1]).nonzero()[0]
                if len(candidats) > 0:
                    contoursPoints[-1].extend([edgesList[candidats[0]][1]])
                    edgesList.pop(candidats[0])
                    contoursTriangles[-1].extend([trianglesList.pop(candidats[0])])
                else:
                    break
                
            contoursPoints[-1].reverse()
            contoursTriangles[-1].reverse()
            while 1:
                if len(trianglesList) == 0:
                    break
                edgesListArray = concatenate(edgesList).reshape(-1,2)
                candidats = (edgesListArray[:,1] == contoursPoints[-1][-1]).nonzero()[0]
                if len(candidats) > 0:
                    contoursPoints[-1].extend([edgesList[candidats[0]][0]])
                    edgesList.pop(candidats[0])
                    contoursTriangles[-1].extend([trianglesList.pop(candidats[0])])
                else:
                    break
            contoursPoints[-1].reverse()
            contoursTriangles[-1].reverse()
            
        self.contoursPoints = contoursPoints
        self.contoursTriangles = contoursTriangles        

            
    def SearchMeshesSimpleSurface(self):
        self.SearchNeighborsSimpleSurface()
        self.meshesSimpleSurface = self.CalculateMeshesUp3(self.neighborsSimpleSurface)
            
    def SearchNeighborsSimpleSurface(self):
        listTrianglesInVertex = []
        self.neighborsSimpleSurface = []
        for i in range(self.points.shape[0]):
            listTrianglesInVertex.append( (self.triangles == i).any(1).nonzero()[0])
            
        for i in range(self.triangles.shape[0]):
            aux = []
            for i0, i1 in ((0,1),(1,2),(2,0)):
                neighbors_aux = intersect1d( listTrianglesInVertex[self.triangles[i,i0]], listTrianglesInVertex[self.triangles[i,i1]])
                if len(neighbors_aux) == 2:
                    aux.append(neighbors_aux)
            if aux: # este if solo esta puesto para el caso de un triangulo sin vecinos (Entonces seria una malla de un triangulo)
                aux = concatenate(aux,0)
                aux = aux[(aux!=i).nonzero()[0]]
                self.neighborsSimpleSurface.append(aux)
            else:
                self.neighborsSimpleSurface.append(array([], 'int'))

    def SearchNeighborsOrderUp3(self):
        listTrianglesInVertex = []
        self.neighbors = zeros((self.triangles.shape[0],3), dtype=object) - 1

        # crear lista de triangulos en vertices    
        for i in range(self.points.shape[0]):
            listTrianglesInVertex.append( (self.triangles == i).any(1).nonzero()[0])
            
        for i in range(self.triangles.shape[0]):
            # para cada una de los lados de cada triangulo
            neighbors_aux = intersect1d( listTrianglesInVertex[self.triangles[i,0]], listTrianglesInVertex[self.triangles[i,1]]) # triangulos comunes a los 2 puntos del triangulo
            if len(neighbors_aux) > 1:
                self.neighbors[i,0] = neighbors_aux[(neighbors_aux!=i).nonzero()[0]]

            neighbors_aux = intersect1d( listTrianglesInVertex[self.triangles[i,1]], listTrianglesInVertex[self.triangles[i,2]]) # triangulos comunes a los 2 puntos del triangulo
            if len(neighbors_aux) > 1:
                self.neighbors[i,1] = neighbors_aux[(neighbors_aux!=i).nonzero()[0]]

            neighbors_aux = intersect1d( listTrianglesInVertex[self.triangles[i,2]], listTrianglesInVertex[self.triangles[i,0]]) # triangulos comunes a los 2 puntos del triangulo
            if len(neighbors_aux) > 1:
                self.neighbors[i,2] = neighbors_aux[(neighbors_aux!=i).nonzero()[0]]

    def SearchPointToTriangle(self):
        self.pointToTriangle = [[] for i in range(self.points.shape[0])]
        
        for i in range(self.triangles.shape[0]):
            self.pointToTriangle[self.triangles[i,0]].append(i)
            self.pointToTriangle[self.triangles[i,1]].append(i)
            self.pointToTriangle[self.triangles[i,2]].append(i)
        
        
    def SearchExternalEdgesAndNeighborsSimpleSurface(self):
        '''
        Busca los bordes externos de la malla y los vecinos considerando una superficie simple, es decir, si hay mas de 2 triangulos compartiendo una arista no se consideran como vecinos
        '''
        listTrianglesInVertex = []
        self.neighborsSimpleSurface = []
        for i in range(self.points.shape[0]):
            listTrianglesInVertex.append( (self.triangles == i).any(1).nonzero()[0])

        segmentsListExternalContours = []
        trianglesListExternalContours = []
        for i in range(self.triangles.shape[0]):
            aux = []
            for i0, i1 in ((0,1),(1,2),(2,0)):
                neighbors_aux = intersect1d( listTrianglesInVertex[self.triangles[i,i0]], listTrianglesInVertex[self.triangles[i,i1]])
                if len(neighbors_aux) == 2: # solo son vecinos los triangulos de una superficie simple
                    aux.append(neighbors_aux)
                elif len(neighbors_aux) == 1: # el triangulo no tiene vecinos en esta arista
                    segmentsListExternalContours.append(self.triangles[i,[i0, i1]])
                    trianglesListExternalContours.append(i)
                    
            if aux: # este if solo esta puesto para el caso de un triangulo sin vecinos (Entonces seria una malla de un triangulo)
                aux = concatenate(aux,0)
                aux = aux[(aux!=i).nonzero()[0]]
                self.neighborsSimpleSurface.append(aux)
            else:
                self.neighborsSimpleSurface.append(array([], 'int'))


        self.contoursTriangles = []
        self.contoursPoints = []
        while trianglesListExternalContours: # loop para los contornos. N iteraciones => N contornos:
            self.contoursTriangles.append([])
            self.contoursPoints.append([])
            self.contoursTriangles[-1].append(trianglesListExternalContours.pop(0))
            self.contoursPoints[-1].extend(segmentsListExternalContours.pop(0))
            
            if trianglesListExternalContours:
                while len(segmentsListExternalContours): # alargar la linea hacia un lado
                    index = (array(segmentsListExternalContours)[:,0] ==  self.contoursPoints[-1][-1]).nonzero()[0]
                    if len(index):
                        index = index[0]
                        self.contoursTriangles[-1].append(trianglesListExternalContours.pop(index))
                        self.contoursPoints[-1].append(segmentsListExternalContours.pop(index)[1])
                    else:
                        if self.contoursPoints[-1][-1] != self.contoursPoints[-1][0]:
                            self.contoursTriangles[-1].reverse()
                            self.contoursPoints[-1].reverse()
                        break
                        
                if self.contoursPoints[-1][-1] != self.contoursPoints[-1][0]:
                    while len(segmentsListExternalContours): # alargar la linea hacia el otro lado
                        index = (array(segmentsListExternalContours)[:,1] ==  self.contoursPoints[-1][-1]).nonzero()[0]
                        if len(index):
                            index = index[0]
                            self.contoursTriangles[-1].append(trianglesListExternalContours.pop(index))
                            self.contoursPoints[-1].append(segmentsListExternalContours.pop(index)[0])
                        else:
                            break
                    self.contoursTriangles[-1].reverse()
                    self.contoursPoints[-1].reverse()


    def FilterDegenerateTriangles(self, edge_limit, area_limit, print_iter_flag = 0):
        """
        Elimina los triangulos que tiene aristas menores que edge_limit y superficie menor a area_limit
        """
##        edge_limit = 0.5
##        area_limit = 0.1
##        print_iter_flag = 0

        for iter_2 in range(20):
            if print_iter_flag:
                print 'iter total:',iter_2
            for iter in range(10):
                if print_iter_flag:
                    print 'iter largo:',iter
                triangle_list_exim = set([])
                salir_largo = 1
                index_rand = random.rand(self.triangles.shape[0])
        ##        for T_n, one_triangle in enumerate(self.triangles):
                for iii in range(self.triangles.shape[0]):
                    T_n = index_rand.argmax()
                    index_rand[T_n] = -1
                    if ~triangle_list_exim.issuperset(set([T_n])):
                        one_triangle = self.triangles[T_n]
                        L = self.points[one_triangle] - self.points[[one_triangle[1], one_triangle[2], one_triangle[0]]]
                        L = sqrt((L*L).sum(1))
                        index = (L < edge_limit)
                        if index.any():
                            triangle_list_exim.add(T_n)
                            salir_largo = 0
                            if index[0]:
                                p0 = one_triangle[0]
                                p1 = one_triangle[1]
                            elif index[1]:
                                p0 = one_triangle[1]
                                p1 = one_triangle[2]                
                            elif index[2]:
                                p0 = one_triangle[2]
                                p1 = one_triangle[0]
                                
                            e0 = (self.triangles == p0).nonzero()
                            e1 = (self.triangles == p1).nonzero()
                            list_T2 = intersect1d_nu(e0[0], e1[0]) # tienen los dos puntos
                            triangle_list_exim.update(set(list_T2))
                            self.triangles[list_T2,:] = -1


                            self.points[p0] = (self.points[p0] + self.points[p1]) / 2.
        ##                    self.points[p0] = self.points[p1]
                            self.triangles[e1] = p0
                if salir_largo:
                    break
                index = ones(self.triangles.shape[0],'Bool')
                index[list(triangle_list_exim)] = 0
                self.triangles = self.triangles[index]

        # #####################
            for iter in range(1):
        ##        self.ComputeAreaOfTriangles()
                triangle_list_exim = set([])
                if print_iter_flag:
                    print 'iter corto:',iter
                salir_corto = 1
                index_rand = random.rand(self.triangles.shape[0])
                
        ##        for T_n, one_triangle in enumerate(self.triangles):
                triangles_new = array([])
                triangles_new = triangles_new.reshape(-1,3)
                triangles_new = triangles_new.astype('int')        
                for iii in range(self.triangles.shape[0]):            
                    T_n = index_rand.argmax()
                    index_rand[T_n] = -1            
                    if ~triangle_list_exim.issuperset(set([T_n])):
                        one_triangle = self.triangles[T_n]
                        w = self.points[one_triangle[1]] - self.points[one_triangle[0]]
                        v = self.points[one_triangle[2]] - self.points[one_triangle[0]]
                        c = cross(w, v)
                        area = sqrt((c * c).sum()) * 0.5
                        if area < area_limit:
                            salir_corto = 0
                            triangle_list_exim.add(T_n)
                            L = self.points[one_triangle] - self.points[[one_triangle[1], one_triangle[2], one_triangle[0]]]
                            L = sqrt((L*L).sum(1))
                            Lm = argmax(L)
                            
                            if Lm == 0:
                                p0 = one_triangle[0]
                                p1 = one_triangle[1]
                                pL = one_triangle[2]
                            elif Lm == 1:
                                p0 = one_triangle[1]
                                p1 = one_triangle[2]
                                pL = one_triangle[0]
                            elif Lm == 2:
                                p0 = one_triangle[2]
                                p1 = one_triangle[0]
                                pL = one_triangle[1]
                                                    
                            e0 = (self.triangles == p0).nonzero()
                            e1 = (self.triangles == p1).nonzero()
                            eL = (self.triangles == pL).nonzero()
                            

                            self.points[pL] = self.points[one_triangle].sum(0) / 3.
                            self.triangles[T_n,:] = -1
                            
                            list_T2 = set(intersect1d_nu(e0[0], e1[0])) # tienen los dos puntos
                            list_T2.remove(T_n)
                            list_T2 = list(list_T2)
                            T_aux = self.triangles[list_T2]

                                                
        ##                    self.triangles[list_T2][(self.triangles[list_T2] == p1)] = pL
                            TT_aux = self.triangles[list_T2]
                            TT_aux[TT_aux==p1] = pL
                            self.triangles[list_T2] = TT_aux
                            
                            T_aux[(T_aux == p0)] = pL
                            self.triangles = concatenate((self.triangles, T_aux.reshape(-1,3)), 0)
    ##                        triangles_new = concatenate((triangles_new,T_aux.reshape(-1,3)),0)
                            
                if salir_corto:
                    break
                index = ones(self.triangles.shape[0],'Bool')
                index[list(triangle_list_exim)] = 0
                self.triangles = self.triangles[index]
    ##            self.triangles = concatenate((self.triangles,triangles_new.reshape(-1,3)),0)
                            
            if salir_largo and salir_corto:
                self.DeleteNotUsedPoints()
                break
            
            

        
        
        
class AutoIntersectMesh:
    def __init__(self, mesh):
        self.mesh = mesh
        self.edgeEdgeIntersectionTolerance = 0.0001 # esta al cuadrado par  a no tener que calcular la raiz dentro del algoritmo
        self.nearPointsTolerance = 0.0001 
        self.edgeModifyTolerance = 1.
        
        self.edgeIntersectionTolerance = 0.0001 # umbral de cercania entre aristas para calificarla como una interseccion de aristas
        self.vertexIntersectionTolerance = 0.1 # umbral de cercania entre aristas-vertice para calificarla como una interseccion
        
        self.vertexVertexIntersectionTolerance = 0.01 # umbral para unir intersecciones de vertice en ambas mallas
        self.vertexEdgeIntersectionTolerance = 0.01 # umbral para unir intersecciones de vertice e interseccion de arist en ambas mallas        

        self.pointsSelectedIncluded = []

        # ############### para Visualizar
        PtoVTK = PythonToPolyData()
        PtoVTK.SetInputTriangles(mesh)
        PtoVTK.WritePolyData('..\\BorrarVtkDataTrianglesMesh0.vtk')
        # ###############
    def PrepareIntersection(self):
        # Calcular los datos si no estan

##        if not hasattr(self.mesh, 'neighbors'):
##            self.mesh.SearchNeighborsUp3()
        if not hasattr(self.mesh, 'edges'):
            self.mesh.SearchEdgesUp3()
    
        self.pointsOfIntersection = []
        self.pointsOfIntersectionToEdges = []
        self.pointsOfIntersection = array(self.pointsOfIntersection, 'float32')
        self.correspondenceIntersectionPointToVertex = []
        self.correspondenceIntersectionPointToVertex = array(self.correspondenceIntersectionPointToVertex, 'int').reshape(-1,3)

        self.trianglesMeshIntersectionPointsInSurface = [] # Mesh
        for i in range(self.mesh.triangles.shape[0]): self.trianglesMeshIntersectionPointsInSurface.append([])

        self.edgesMeshToIntersectionTriangles = [] # Mesh0ToMesh0 Mesh0ToMesh1
        self.edgesMeshToIntersectionPoints = [] # Mesh0ToMesh0 Mesh0ToMesh1

        for i in range(self.mesh.edges.shape[0]):
            self.edgesMeshToIntersectionTriangles.append([])
            self.edgesMeshToIntersectionPoints.append([])
            
    def ComputeTrianglesGrid(self):

        # definir cuadricula
##        N0 = 50
##        N1 = 50
##        N2 = 50
        
        if not hasattr(self.mesh, 'edges'):
            self.mesh.SearchEdgesUp3()
            
        boundingBoxMin = [self.mesh.points[:,0].min(), self.mesh.points[:,1].min(), self.mesh.points[:,2].min()]
        boundingBoxMax = [self.mesh.points[:,0].max(), self.mesh.points[:,1].max(), self.mesh.points[:,2].max()]

        l0 = boundingBoxMax[0] - boundingBoxMin[0]
        l1 = boundingBoxMax[1] - boundingBoxMin[1]
        l2 = boundingBoxMax[2] - boundingBoxMin[2]

        # definir tamano de la cuadricula en base en lago promedio de las aristas de los triangulos
        lengEdge = self.mesh.points[self.mesh.edges[:,0]] - self.mesh.points[self.mesh.edges[:,1]]
        lengEdge = sqrt((lengEdge * lengEdge).sum(1))
        lengEdgeProm = lengEdge.sum() / len(lengEdge)

        d0 =  lengEdgeProm
        d1 =  lengEdgeProm
        d2 =  lengEdgeProm
        dd = array([d0,d1,d2])
        
        N0 = int(l0 / d0) + 1
        N1 = int(l1 / d1) + 1
        N2 = int(l2 / d2) + 1
        
##        d0 = l0 / N0
##        d1 = l1 / N1
##        d2 = l2 / N2
        # informacion de cuales traingulos tienen su bounding box chocando con cada celda
        gridOfTriangles = empty((N0,N1,N2), dtype=object)
        gridListOfTriangles = empty((N0*N1*N2), dtype=object)
        trianglesListOfGrid = empty((self.mesh.triangles.shape[0]), dtype=object)
        for i in range(self.mesh.triangles.shape[0]): trianglesListOfGrid[i] = []

        index_range = zeros((2,3), 'int')
        for triangleIndex,triangleCurrent in enumerate(self.mesh.triangles):
            minT = self.mesh.points[triangleCurrent].min(0)
            maxT = self.mesh.points[triangleCurrent].max(0)
            
            index_range[0,:] = ((minT - boundingBoxMin) / dd).astype('int')
            index_range[1,:] = ((maxT - boundingBoxMin) / dd).astype('int') + 1
            
            for i in range(index_range[0,0], index_range[1,0]):
                for j in range(index_range[0,1], index_range[1,1]):
                    for k in range(index_range[0,2], index_range[1,2]):
                        index_aux = i + j * N0 + k * N0 * N1
                        if gridOfTriangles[i,j,k] == None:
                            gridOfTriangles[i,j,k] = []
                            gridListOfTriangles[index_aux] = []
                        gridOfTriangles[i,j,k].append(triangleIndex)
                        boxIndex = index_aux
                        gridListOfTriangles[boxIndex].append(triangleIndex)
                        trianglesListOfGrid[triangleIndex].append(boxIndex)
                        
                
        self.gridOfTriangles = gridOfTriangles
        self.gridListOfTriangles = gridListOfTriangles
        self.trianglesListOfGrid = trianglesListOfGrid



        # informacion de cuales aristas tienen su bounding box chocando con cada celda
        gridOfEdges = empty((N0,N1,N2), dtype=object)
        gridListOfEdges = empty((N0*N1*N2), dtype=object)
        edgesListOfGrid = empty((self.mesh.edges.shape[0]), dtype=object)
        for i in range(self.mesh.edges.shape[0]): edgesListOfGrid[i] = []

        for edgeIndex,edgeCurrent in enumerate(self.mesh.edges):
            minT = self.mesh.points[edgeCurrent].min(0)
            maxT = self.mesh.points[edgeCurrent].max(0)


            index_range[0,:] = ((minT - boundingBoxMin) / dd).astype('int')
            index_range[1,:] = ((maxT - boundingBoxMin) / dd).astype('int') + 1
            
            for i in range(index_range[0,0], index_range[1,0]):
                for j in range(index_range[0,1], index_range[1,1]):
                    for k in range(index_range[0,2], index_range[1,2]):
                        index_aux = i + j * N0 + k * N0 * N1
                        if gridOfEdges[i,j,k] == None:
                            gridOfEdges[i,j,k] = []
                            gridListOfEdges[index_aux] = []
                        gridOfEdges[i,j,k].append(edgeIndex)
                        boxIndex = index_aux
                        gridListOfEdges[boxIndex].append(edgeIndex)
                        edgesListOfGrid[edgeIndex].append(boxIndex)
                        
                
        self.gridOfEdges = gridOfEdges
        self.gridListOfEdges = gridListOfEdges
        self.edgesListOfGrid = edgesListOfGrid



    def ComputeIntersections(self):
        """
        Optimizado para calcular solo un area, no todos los triangulos 
        Le agrege deteccion de intersecciones en vertices y calculo inicial de distancias entre puntos y planos de los triangulos
        """        
            
        # crear una grilla en donde cada celda tiene los triangulos que la intersectan. Cada triangulo puede estar en mas de una celda
        if not hasattr(self, 'gridOfTriangles'):
##            TT = timeit.Timer()
##            T1 = TT.timer()  
            gridOfTriangles = self.ComputeTrianglesGrid() # crea self.gridOfTriangles
##            T2 = TT.timer()
##            print 'Time Grid:',T2-T1

        trianglesListIn = arange(self.mesh.triangles.shape[0])
        edgesListIn = arange(self.mesh.edges.shape[0])
        
##        T1 = TT.timer()  
        self.ComputeIntersectionsZone(edgesListIn, trianglesListIn)
##        T2 = TT.timer()
##        print 'Time inter cal:',T2-T1
        
        self.correspondenceIntersectionPointToVertex = array(self.correspondenceIntersectionPointToVertex, 'int').reshape(-1,3)

##        # ############### para Visualizar
##        if len(self.pointsOfIntersection) > 0:
##            PtoVTK = PythonToPolyData()
##            PtoVTK.SetInputPoints(self.pointsOfIntersection)
##            PtoVTK.WritePolyData('..\\BorrarVtkDataIntersectionsPoints.vtk')
##        # ###############



    def ComputeIntersectionsZone(self, edgesListIn, trianglesListIn):
        """
        Le agrege deteccion de intersecciones en vertices y calculo inicial de distancias entre puntos y planos de los triangulos
        """
        # Iteracion para las 2 mallas intersectando con la otra y entre si
        self.pointsOfIntersection = self.pointsOfIntersection.tolist()
        PI_count = len(self.pointsOfIntersection)
        
        pointListMeshEdges = self.mesh.edges[edgesListIn,:].flatten().tolist() # lista de los puntos que forman las aristas
        pointListMeshEdges = set(pointListMeshEdges) # para eliminar puntos repetidos
        pointListMeshEdges = list(pointListMeshEdges)
        
        # puntos de cada triangulo
        T0_aux = self.mesh.points[self.mesh.triangles[trianglesListIn,0],:].astype('float64')
        T1_aux = self.mesh.points[self.mesh.triangles[trianglesListIn,1],:].astype('float64')
        T2_aux = self.mesh.points[self.mesh.triangles[trianglesListIn,2],:].astype('float64')                  
            
        # calculo de las normales de los triangulos de la malla Triangles
##                if not hasattr(meshTriangles, 'normalsOfTriangles'):
        self.mesh.ComputeNormalsOfTriangles()            
            
        for i in edgesListIn:
            # cajas de las aristas
            boxCurrentEdge = self.edgesListOfGrid[i]
            trianglesOfCurrentEdge = self.mesh.edgesToTriangles[i]
            edgePoints = self.mesh.edges[i,:]
            
            # todos los triangulos que intersectan con las cajas
            trianglesOfCurrentBoxes = []
            for box_index in boxCurrentEdge:
                trianglesOfCurrentBoxes.extend(self.gridListOfTriangles[box_index])
            trianglesOfCurrentBoxes = array(list(set(trianglesOfCurrentBoxes)))
##            trianglesOfCurrentBoxes.difference_update(set(trianglesOfCurrentEdge)) # eliminar los triangulos a los cuales pertenece la arista que estoy analizando
            # eliminar los planos de los triangulos vecinos se incluyen los de la arista
            trianglesOfCurrentBoxes = trianglesOfCurrentBoxes[(((self.mesh.triangles[trianglesOfCurrentBoxes] != edgePoints[0]) * (self.mesh.triangles[trianglesOfCurrentBoxes] != edgePoints[1])).all(1)).nonzero()[0]] # con respecto a planesWithIntersection

            # distancia de los puntos de la arista al los planos de los triangulo
            d0_aux = ((T1_aux[trianglesOfCurrentBoxes] - self.mesh.points[edgePoints[0],:].astype('float64')) * self.mesh.normalsOfTriangles[trianglesOfCurrentBoxes,: ]).sum(1)
            d1_aux = ((T1_aux[trianglesOfCurrentBoxes] - self.mesh.points[edgePoints[1],:].astype('float64')) * self.mesh.normalsOfTriangles[trianglesOfCurrentBoxes,: ]).sum(1)

            
            indexPlanesWithIntersection = ((d0_aux * d1_aux)<=0).nonzero()[0] # indice con respecto a trianglesOfCurrentBoxes
            if indexPlanesWithIntersection.shape[0]:
                planesWithIntersection = trianglesListIn[trianglesOfCurrentBoxes[indexPlanesWithIntersection]] # planos en la lista comopleta de la malla o trianglesListIn

##                    if meshEdges==meshTriangles:trianglesListIn[meshTrianglesNo][planesWithIntersection_aux]
##                        for mismoPlano in meshEdges.edgesToTriangles[i]:
##                            if len((planesWithIntersection_aux == mismoPlano).nonzero()[0])>0:
##                                print 'alto'
##                            planesWithIntersection_aux = planesWithIntersection_aux[(planesWithIntersection_aux != mismoPlano).nonzero()[0]]
                
                # eliminar de la lista los planos de los triangulos con los cuales ya tiene una interseccion, para no pasar sobre el mismo triangulo, si ya lo marque por interseccion de aristas o vertice
                planesWithIntersection_aux = []
                indexPlanesWithIntersection_aux = []
                for ii in range(len(planesWithIntersection)):
                    plane_aux = planesWithIntersection[ii]
                    indexPlane_aux = indexPlanesWithIntersection[ii]
                    if not self.edgesMeshToIntersectionTriangles[i].count(plane_aux):
                        planesWithIntersection_aux.append(plane_aux)
                        indexPlanesWithIntersection_aux.append(indexPlane_aux)
                planesWithIntersection = planesWithIntersection_aux[:]
                indexPlanesWithIntersection = indexPlanesWithIntersection_aux[:]
                
                if len(planesWithIntersection)>0:
                    Vsegment_aux = self.mesh.points[edgePoints[1],:].astype('float64') - self.mesh.points[edgePoints[0],:].astype('float64')
                    for iPlane in range(len(planesWithIntersection)):
                        plane = planesWithIntersection[iPlane]
                        
                        if not self.edgesMeshToIntersectionTriangles[i].count(plane): # por si se agrego antes con una interseccion de vertice
                            indexPlane = indexPlanesWithIntersection[iPlane] # indice del triangulo en trianglesOfCurrentBoxes
                            
                            # ###### precision multiple
##                                    import mpmath
##                                    mpmath.mp.dps = 60
##                                    Vsegment_aux_mpf = empty((1,3), 'object')
##                                    for jj in range(3):
##                                        Vsegment_aux_mpf[0,jj] = mpmath.mpf(str(meshEdges.points[meshEdges.edges[i,1],jj])) - mpmath.mpf(str(meshEdges.points[meshEdges.edges[i,0],jj]))
##                                        
##                                    aux0 = (mpmath.mpf(str(T1_aux[plane][0])) - mpmath.mpf(str(meshEdges.points[meshEdges.edges[i,0],0]))) * meshTriangles.normalsOfTriangles_mpf[trianglesListIn[meshTrianglesNo][plane],0 ]
##                                    aux1 = (mpmath.mpf(str(T1_aux[plane][1])) - mpmath.mpf(str(meshEdges.points[meshEdges.edges[i,0],1]))) * meshTriangles.normalsOfTriangles_mpf[trianglesListIn[meshTrianglesNo][plane],1 ]
##                                    aux2 = (mpmath.mpf(str(T1_aux[plane][2])) - mpmath.mpf(str(meshEdges.points[meshEdges.edges[i,0],2]))) * meshTriangles.normalsOfTriangles_mpf[trianglesListIn[meshTrianglesNo][plane],2 ]
##                                    d0_aux_mpf = aux0 + aux1 + aux2
##                                    
##                                    r_mpf = d0_aux_mpf / ( Vsegment_aux_mpf * meshTriangles.normalsOfTriangles_mpf[plane,:]).sum()
##                                    
##                                    vertexEdge_aux = empty((1,3), 'object')
##                                    for jj in range(3):
##                                        vertexEdge_aux[0,jj] = mpmath.mpf(str(meshEdges.points[meshEdges.edges[i,0],jj].astype('float64')))
##                                    PI_mpf = vertexEdge_aux + (r_mpf * Vsegment_aux_mpf)
##                                    
##                                    PI_mpf_aux = zeros(3, 'float32')
##                                    PI_mpf_aux[0] = float(PI_mpf[0,0])
##                                    PI_mpf_aux[1] = float(PI_mpf[0,1])
##                                    PI_mpf_aux[2] = float(PI_mpf[0,2])
##                                    PI_mpf_aux = PI_mpf_aux.reshape(1,3)
##                                    PI_mpf = PI_mpf_aux
                            # ###### 
                            
                           
                            r = d0_aux[indexPlane] / ( Vsegment_aux * self.mesh.normalsOfTriangles[plane]).sum()
                            PI = self.mesh.points[edgePoints[0],:].astype('float64') + (r * Vsegment_aux)
                            # comprobar si la interseccion cae en el triangulo.
                            # La ecuacion parametrica del punto en el plano del triangulo seria: T(s,t) = T0_aux + s(T1_aux-T0_aux) + t(T2_aux-T0_aux) = T0_aux + su + tv
                            u = T1_aux[plane] - T0_aux[plane]
                            v = T2_aux[plane] - T0_aux[plane]
                            w = PI - T0_aux[plane]
                            
                            uv = dot(u,v)
                            wu = dot(w,u)
                            wv = dot(w,v)
                            vv = dot(v,v)
                            uu = dot(u,u)
                            denom = (uv * uv) - uu * vv
                            s = (uv * wv - vv * wu) / denom
                            t = (uv * wu - uu * wv) / denom

                            if  abs(denom) < 1e-016:
                                print 'Posible triangulo degenerado!!!'
##                                s = 10 # para que no tome en cuenta la interseccion
##                                t = 10

##                            vp = cross(meshTriangles.normalsOfTriangles[plane], v)                            
##                            up = cross(meshTriangles.normalsOfTriangles[plane], u)
##                            s = (w*vp).sum() / (u*vp).sum()
##                            t = (w*up).sum() / (v*up).sum()
                            
                            if (s>=-0.00001 and t>=-0.00001 and (s+t)<=1.00001):


##                                if exec("edgesEdgeIntersectionEdgeMesh" + meshTrianglesStr + "ToTriangles[k]")
##                                
##                                if (s+t) > self.edgeEdgeIntersectionTolerance:
##                                    # encuentro las aristas del triangulo de la interseccion
##                                    edgesCheck = meshTriangles.trianglesToEdges[plane]
##                                    trianglesCheck = array(meshEdges.edgesToTriangles[i])
##                                    for k in edgesCheck:
##                                        exec("edgesEdgeIntersectionEdgeMesh" + meshTrianglesStr + "ToTriangles[k] = [trianglesCheck," + meshEdgesStr + "]")
####                                    edgesEdgeIntersectionEdgeMesh'meshTrianglesStr'ToTriangles[edgesCheck] = [trianglesCheck, 'meshEdgesStr']
##                                        exec("edgesEdgeIntersectionEdgeMesh" + meshTrianglesStr + "ToPoints[k].append(PI_count)")
####                                    edgesEdgeIntersectionEdgeMesh'meshTrianglesStr'ToPoints[edgesCheck].append(PI_count)
##

##                                if len(self.pointsOfIntersection)>0: # ver si por tolerancia coincide con otro punto
##                                    distance = array(self.pointsOfIntersection) - PI
##                                    distance = sqrt(distance*distance).sum(1)
##                                    rep_index = (distance < self.nearPointsTolerance).nonzero()[0]
##                                else:
##                                    rep_index = []
                                rep_index = [] # lo anterior no lo estoy tomando en cuanta. BORRAR si no funciona bien
                                    

                                if abs(s)<0.0001 or abs(t)<0.0001 or abs(1-(s+t))<0.0001: # si hay posibilidad de interseccion de vertice o de arista
                                    # distancia desde el punto a las aristas
                                    edgesOfTriangle = self.mesh.trianglesToEdges[plane,:]
                                    p0 = self.mesh.points[self.mesh.edges[edgesOfTriangle][:,0]]
                                    p1 = self.mesh.points[self.mesh.edges[edgesOfTriangle][:,1]]                                    
                                    d = p1 - p0
                                    d2 = (d*d).sum(1)
                                    
                                    tt = ((d * (PI - p0)).sum(1) / d2).reshape(-1,1)
                                    p = p0 + tt * d
                                    p[(tt <= 0.).nonzero()[0],:] = p0[(tt <= 0.).nonzero()[0],:]
                                    p[(tt >= 1.).nonzero()[0],:] = p1[(tt >= 1.).nonzero()[0],:]
                                    dist = (PI - p)
                                    dist = sqrt((dist * dist).sum(1))
                                    indexMinEdgeDistance = dist.argmin()
                                    valMinEdgeDistance = dist[indexMinEdgeDistance]

                                    pointsOfIntersection = self.mesh.triangles[plane,:]
                                    distance = self.mesh.points[pointsOfIntersection,:] - PI
                                    distance = sqrt((distance*distance).sum(1))
                                    indexMinVertexDistance = distance.argmin()
                                    valMinVertexDistance = distance[indexMinVertexDistance]
                                    

                                if (abs(s)<0.0001 or abs(t)<0.0001 or abs(1-(s+t))<0.0001) and (valMinEdgeDistance < self.edgeIntersectionTolerance): # posible interseccion de aristas
                                    # actualizar intersecciones
                                    edgeMinDist = edgesOfTriangle[indexMinEdgeDistance]

                                    if len(rep_index)==0: 
                                        self.pointsOfIntersection.append(PI)
                                        self.pointsOfIntersectionToEdges.append([i])                                          
                                        PI_index = PI_count
                                        PI_count += 1
                                        self.pointsOfIntersectionToEdges[PI_index].append(edgeMinDist)
                                    else: 
                                        PI_index = rep_index[0]                                            
                                        
                                    intersectedTrianglesList = self.mesh.edgesToTriangles[edgeMinDist]
                                    for oneTriangle in intersectedTrianglesList: # agregar el punto a los triangulos intersectados de la malla Edges
                                        if not self.trianglesMeshIntersectionPointsInSurface[oneTriangle].count(PI_index): # agregar el punto al triangulo intersectado
                                            self.trianglesMeshIntersectionPointsInSurface[oneTriangle].append(PI_index)
                                            
                                    for oneTriangle in intersectedTrianglesList:
                                        if not self.edgesMeshToIntersectionTriangles[i].count(oneTriangle):
                                            self.edgesMeshToIntersectionTriangles[i].append(oneTriangle)
                                            self.edgesMeshToIntersectionPoints[i].append(PI_index)
                                        else: # si la arista ya tenia una interseccion con este triangulo la cambio por el nuevo punto
                                            index = self.edgesMeshToIntersectionTriangles[i].index(oneTriangle)
                                            self.edgesMeshToIntersectionPoints[i][index] = PI_index
                                            
                                            
                                    intersectedTrianglesList = self.mesh.edgesToTriangles[i]
                                    for oneTriangle in intersectedTrianglesList: # agregar el punto a los triangulos intersectados de la malla Edges
                                        if not self.trianglesMeshIntersectionPointsInSurface[oneTriangle].count(PI_index): # agregar el punto al triangulo intersectado
                                            self.trianglesMeshIntersectionPointsInSurface[oneTriangle].append(PI_index)                                        
                                    for oneTriangle in intersectedTrianglesList:
##                                            if not self.edgesMeshToIntersectionPoints[meshTrianglesNo][meshEdgesNo][edgeMinDist].count(PI_index):
                                        self.edgesMeshToIntersectionTriangles[edgeMinDist].append(oneTriangle)
                                        self.edgesMeshToIntersectionPoints[edgeMinDist].append(PI_index)
##                                        else: # si la arista ya tenia una interseccion con este triangulo la cambio por el nuevo punto
##                                            index = self.edgesMeshToIntersectionTriangles[meshTrianglesNo][meshEdgesNo][edgeMinDist].index(oneTriangle)
##                                            self.edgesMeshToIntersectionPoints[meshTrianglesNo][meshEdgesNo][edgeMinDist][index] = PI_index
                                                
                                else: # interseccion en una cara
                                    
                                    if not self.edgesMeshToIntersectionTriangles[i].count(plane): # por si se agrego antes con una interseccion de vertice
                                        if len(rep_index)==0: # interseccion en una cara, por lo tanto se agrega en punto
                                            self.pointsOfIntersection.append(PI)
                                            self.pointsOfIntersectionToEdges.append([i])
                                            PI_index = PI_count
                                            PI_count += 1
                                        else: # hay interseccion de aristas, por lo tanto el punto es el mismo que otro anterior
                                            PI_index = rep_index[0]
                                            if not self.pointsOfIntersectionToEdges[PI_index].count(i): # para no repetir la arista al intersectar con caras contiguas
                                                self.pointsOfIntersectionToEdges[PI_index].append(i)
                                        self.trianglesMeshIntersectionPointsInSurface[plane].append(PI_index)
                                        if not self.edgesMeshToIntersectionTriangles[i].count(plane): # para no agregarlo si ya estaba a causa de una interseccion de arista o de vertices
                                            self.edgesMeshToIntersectionTriangles[i].append(plane)
                                            self.edgesMeshToIntersectionPoints[i].append(PI_index)                                

##        # ############### para Visualizar
##        PtoVTK = PythonToPolyData()
##        PtoVTK.SetInputPoints(array(self.pointsOfIntersection, 'float32'))
##        PtoVTK.WritePolyData('..\\BorrarVtkDataIntersectionsPoints.vtk')
##        # ###############
        
        self.pointsOfIntersection = array(self.pointsOfIntersection, 'float32')
        del T0_aux, T1_aux, T2_aux




    
class PolyDataToPython:
    def __init__(self):
        pass     
    def ExtractPoints(self):
        """
        vtkData es un vtkPolyData
        """    
        vtkFloatArrayPoints = self.vtkPolyDataIn.GetPoints().GetData()
        arrayPoints=zeros( ( vtkFloatArrayPoints.GetNumberOfTuples(), vtkFloatArrayPoints.GetNumberOfComponents() ),  'f')
        vtkFloatArrayPoints.ExportToVoidPointer(arrayPoints)
        arrayPoints_bak=arrayPoints.copy()   # para dejar los puntos de la misma manera que la matriz numpy con la imagen
        arrayPoints[:,0] = arrayPoints_bak[:,2]
        arrayPoints[:,2] = arrayPoints_bak[:,0]
        self.arrayPoints = arrayPoints

    def ExtractTriangles(self):
        vtkIdTypeArrayTriangles = self.vtkPolyDataIn.GetPolys().GetData()
        arrayTriangles = zeros((vtkIdTypeArrayTriangles.GetNumberOfTuples(),1), 'i')
        vtkIdTypeArrayTriangles.ExportToVoidPointer(arrayTriangles)
        arrayTriangles = arrayTriangles.reshape(-1,4)[:,1:4]
        self.arrayTriangles = arrayTriangles

        arrayTriangles_bak = arrayTriangles.copy()   # para dejar los triangulos con la normal hacia afuera
        arrayTriangles[:,0] = arrayTriangles_bak[:,2]
        arrayTriangles[:,2] = arrayTriangles_bak[:,0]
        
    def ExtractFaces(self):
        vtkIdTypeArrayFaces = self.vtkPolyDataIn.GetPolys().GetData()
        arrayFaces = zeros((vtkIdTypeArrayFaces.GetNumberOfTuples(),1), 'i')
        vtkIdTypeArrayFaces.ExportToVoidPointer(arrayFaces)

        faces = []
        for i in range(30):
            faces.append([])
        n = 0
        while n < arrayFaces.shape[0]:
            largo = arrayFaces[n,0]
            faces[largo].append(arrayFaces[n+1:n+largo+1,0])
            n += largo+1
            
        index = range(len(faces))
        index.reverse()
        for i in index:
            if len(faces[i]) == 0:
                faces.pop()
            else:
                break
        for i in range(len(faces)):
            faces[i] = array(faces[i])

        self.faces = faces    

    def ExtractCellScalars(self):

        # escalares de las caras
        vtkFloatArrayScalars = self.vtkPolyDataIn.GetCellData().GetScalars()
        if vtkFloatArrayScalars != None:
            arrayScalars = zeros((vtkFloatArrayScalars.GetNumberOfTuples(),vtkFloatArrayScalars.GetNumberOfComponents()), 'float32')
            vtkFloatArrayScalars.ExportToVoidPointer(arrayScalars)
            self.arrayCellScalars = arrayScalars
        else:
            self.arrayCellScalars = array([],'float32')
        
    def SetInput(self, vtkPolyDataIn):
        self.vtkPolyDataIn = vtkPolyDataIn
    def ReadFile(self, fileName):
        polyDataReader = vtkPolyDataReader()
        polyDataReader.SetFileName(fileName)
        polyDataReader.Update()
        self.vtkPolyDataIn = polyDataReader.GetOutput()
        
    def GetPolyDataIn(self):
        return self.vtkPolyDataIn
        
    def GetOutput(self):
        """
        Solo para mallas de triangulos
        """
        self.ExtractPoints()
        if (self.vtkPolyDataIn.GetPolys().GetData().GetNumberOfTuples()) == self.vtkPolyDataIn.GetNumberOfCells()*4:
            self.ExtractTriangles()
            self.Output = TrianglesMesh(self.arrayPoints, self.arrayTriangles)
        else:
            self.ExtractFaces()
            self.Output = SimplexMesh(points = self.arrayPoints, faces = self.faces)
    
        self.ExtractCellScalars()
        scalarsCheck = ones(self.arrayCellScalars.shape[0], 'int')
        if self.arrayCellScalars.shape[0] > 1:
            meshes = []            
            for i in range(self.arrayCellScalars.shape[0]):
                if scalarsCheck[i]:
                    index = (self.arrayCellScalars == self.arrayCellScalars[i]).nonzero()[0]
                    meshes.append(index)
                    scalarsCheck[index] = 0
                    if not scalarsCheck.any():
                        break
            self.Output.meshes = meshes
        return self.Output
   
class PythonToPolyData:
    def __init__(self):
        pass
    def SetInput(self, mesh):
        if isinstance(mesh, SimplexMesh):
            self.SetInputSimplex(mesh)
        elif isinstance(mesh, TrianglesMesh):
            self.SetInputTriangles(mesh)
        elif isinstance(mesh, TetrahedraMesh):
            self.SetInputTetrahedra(mesh)
        else:
            print 'formato desconocido'
        
    def SetInputSimplex(self, simplexMesh, curvaturePoints = False, scalarsInPoints = False):
        
        self.__pointsPolyData = simplexMesh.points.copy() #para que queden en el orden que ocupa VTK
        self.__pointsPolyData[:,0] = simplexMesh.points[:,2]
        self.__pointsPolyData[:,2] = simplexMesh.points[:,0]

        self.__vtkData = vtkPolyData()

        # Points
        vtkFloatArrayPoints=vtkFloatArray()
        vtkFloatArrayPoints.SetNumberOfComponents(3)
        vtkFloatArrayPoints.SetVoidArray(self.__pointsPolyData, self.__pointsPolyData.shape[0]*3 ,1)

        vtkMyPoints=vtkPoints()
        vtkMyPoints.SetData(vtkFloatArrayPoints)

        self.__vtkData.SetPoints(vtkMyPoints)

        # Faces
        # Ordenarlos
        NDatos = 0
        for i in range(len(simplexMesh.faces)):
            NDatos += (i+1)*len(simplexMesh.faces[i])
        self.__faces=zeros((NDatos,1), 'int')
        count=0
        NCells=0
        for i in range(len(simplexMesh.faces)):
            if len(simplexMesh.faces[i]):
                NCells += len(simplexMesh.faces[i])
                facesAux = zeros((simplexMesh.faces[i].shape[0],simplexMesh.faces[i].shape[1]+1), 'int')
                facesAux[:,1:] = simplexMesh.faces[i]
                facesAux[:,0]=i
                facesAux = facesAux.reshape((-1,1))
                self.__faces[count:count + facesAux.shape[0]] = facesAux
                count = count + facesAux.shape[0] 

        vtkFacesArray=vtkIdTypeArray()    
        vtkFacesArray.SetVoidArray(self.__faces, len(self.__faces), 1)

        vtkCellArrayFaces=vtkCellArray()
        vtkCellArrayFaces.SetCells(NCells,vtkFacesArray)

        self.__vtkData.SetPolys(vtkCellArrayFaces)    

        if curvaturePoints:
            if simplexMesh.meanCurvatureInPoints.dtype == 'float32':
                vtkArrayScalars = vtkFloatArray()
            elif simplexMesh.meanCurvatureInPoints.dtype == 'float64':
                vtkArrayScalars = vtkDoubleArray()
            vtkArrayScalars.SetNumberOfComponents(1)
            vtkArrayScalars.SetVoidArray(simplexMesh.meanCurvatureInPoints, simplexMesh.meanCurvatureInPoints.shape[0] ,1)
            self.__vtkData.GetPointData().SetScalars(vtkArrayScalars)

        if scalarsInPoints:
            if simplexMesh.scalarsInPoints.dtype == 'float32':
                vtkArrayScalars = vtkFloatArray()
            elif simplexMesh.scalarsInPoints.dtype == 'float64':
                vtkArrayScalars = vtkDoubleArray()
            vtkArrayScalars.SetNumberOfComponents(1)
            vtkArrayScalars.SetVoidArray(simplexMesh.scalarsInPoints, simplexMesh.scalarsInPoints.shape[0] ,1)
            self.__vtkData.GetPointData().SetScalars(vtkArrayScalars)



##        dict={'vtkData': vtkData, 'Faces':faces, 'PointsVTK':simplexVTK}
    
##        return dict
        
    def SetInputTriangles(self, trianglesMesh, colorMesh = False, curvaturePoints = False):
        
        # Points            
        self.__vtkData = vtkPolyData()
        if trianglesMesh.points.dtype == 'float32':
            vtkFloatArrayPoints = vtkFloatArray()
        elif trianglesMesh.points.dtype == 'float64':
            vtkFloatArrayPoints = vtkDoubleArray()
        vtkFloatArrayPoints.SetNumberOfComponents(3)
        self.__pointsPolyData = trianglesMesh.points.copy()
        ##pointsNewPolyData[:,0] = pointsNew[:,2]; pointsNewPolyData[:,2] = pointsNew[:,0]+10
        self.__pointsPolyData[:,0] = trianglesMesh.points[:,2]; self.__pointsPolyData[:,2] = trianglesMesh.points[:,0]
        vtkFloatArrayPoints.SetVoidArray(self.__pointsPolyData, trianglesMesh.points.shape[0]*3 ,1)

        vtkMyPoints=vtkPoints()
        vtkMyPoints.SetData(vtkFloatArrayPoints)

        self.__vtkData.SetPoints(vtkMyPoints)

        # Faces
        # Ordenarlos
        triangles_aux = trianglesMesh.triangles.copy()
        triangles_aux[:,0] = trianglesMesh.triangles[:,2]; triangles_aux[:,2] = trianglesMesh.triangles[:,0]
        self.__faces = concatenate((ones((trianglesMesh.triangles.shape[0],1),'int')*3, triangles_aux),1).flatten()

        vtkFacesArray = vtkIdTypeArray()    
        vtkFacesArray.SetVoidArray(self.__faces, len(self.__faces), 1)

        vtkCellArrayFaces=vtkCellArray()
        vtkCellArrayFaces.SetCells(trianglesMesh.triangles.shape[0],vtkFacesArray)
        
        self.__vtkData.SetPolys(vtkCellArrayFaces)

     
        # escalares de las caras
        self.scalarFaces = zeros(trianglesMesh.triangles.shape[0], 'float32')
        if colorMesh:            
            for i in range(len(trianglesMesh.meshes)):
                self.scalarFaces[trianglesMesh.meshes[i]] = i
        vtkFloatArrayScalars = vtkFloatArray()
        vtkFloatArrayScalars.SetNumberOfComponents(1)
        vtkFloatArrayScalars.SetVoidArray(self.scalarFaces, self.scalarFaces.shape[0] ,1)

##        vtkIntArrayScalars=vtkIntArray()
##        vtkIntArrayScalars.SetNumberOfComponents(1)
##        vtkIntArrayScalars.SetVoidArray(self.__scalarFaces, self.__scalarFaces.shape[0] ,1)


        self.__vtkData.GetCellData().SetScalars(vtkFloatArrayScalars)


##        # escalares de los puntos
##        self.scalarPoints = zeros(trianglesMesh.points.shape[0], 'float32')
##        if colorMesh:            
##            for i in range(len(trianglesMesh.meshes)):
##                self.scalarPoints[ trianglesMesh.triangles[trianglesMesh.meshes[i]].flatten() ] = i
####        vtkFloatArrayScalars = vtkFloatArray()
####        vtkFloatArrayScalars.SetNumberOfComponents(1)
####        vtkFloatArrayScalars.SetVoidArray(self.scalarPoints, self.scalarPoints.shape[0] ,1)
##
##        vtkIntArrayScalars=vtkIntArray()
##        vtkIntArrayScalars.SetNumberOfComponents(1)
##        vtkIntArrayScalars.SetVoidArray(self.scalarPoints, self.scalarPoints.shape[0] ,1)
##
##        self.__vtkData.GetPointData().SetScalars(vtkIntArrayScalars)

        if curvaturePoints:
            if trianglesMesh.meanCurvatureInPoints.dtype == 'float32':
                vtkArrayScalars = vtkFloatArray()
            elif trianglesMesh.meanCurvatureInPoints.dtype == 'float64':
                vtkArrayScalars = vtkDoubleArray()
            vtkArrayScalars.SetNumberOfComponents(1)
            vtkArrayScalars.SetVoidArray(trianglesMesh.meanCurvatureInPoints, trianglesMesh.meanCurvatureInPoints.shape[0] ,1)
            self.__vtkData.GetPointData().SetScalars(vtkArrayScalars)

        
        ##        dict={'vtkData': vtkData, 'Faces':Faces, 'PointsVTK':pointsPolyData}
        ##        return dict

    def SetInputTetrahedra(self, tetrahedraMesh):
        # Points            
        self.__vtkData = vtkPolyData()
        if tetrahedraMesh.points.dtype == 'float32':
            vtkFloatArrayPoints = vtkFloatArray()
        elif tetrahedraMesh.points.dtype == 'float64':
            vtkFloatArrayPoints = vtkDoubleArray()
        vtkFloatArrayPoints.SetNumberOfComponents(3)
        self.__pointsPolyData = tetrahedraMesh.points.copy()
        ##pointsNewPolyData[:,0] = pointsNew[:,2]; pointsNewPolyData[:,2] = pointsNew[:,0]+10
        self.__pointsPolyData[:,0] = tetrahedraMesh.points[:,2]; self.__pointsPolyData[:,2] = tetrahedraMesh.points[:,0]
        vtkFloatArrayPoints.SetVoidArray(self.__pointsPolyData, tetrahedraMesh.points.shape[0]*3 ,1)

        vtkMyPoints=vtkPoints()
        vtkMyPoints.SetData(vtkFloatArrayPoints)

        self.__vtkData.SetPoints(vtkMyPoints)

        # Faces
        # Ordenarlos
        N = tetrahedraMesh.tetrahedra.shape[0]
        triangles_aux = zeros((N * 4, 3), 'int')
        
        triangles_aux[:N,:] = tetrahedraMesh.tetrahedra[:,[2,1,0]]
        triangles_aux[N:N*2,:] = tetrahedraMesh.tetrahedra[:,[1,3,0]]
        triangles_aux[N*2:N*3,:] = tetrahedraMesh.tetrahedra[:,[2,3,1]]
        triangles_aux[N*3:N*4,:] = tetrahedraMesh.tetrahedra[:,[0,3,2]]
        
        self.__tetra = concatenate((ones((triangles_aux.shape[0],1),'int')*3, triangles_aux),1).flatten()

        vtkFacesArray = vtkIdTypeArray()    
        vtkFacesArray.SetVoidArray(self.__tetra, len(self.__tetra), 1)

        vtkCellArrayFaces=vtkCellArray()
        vtkCellArrayFaces.SetCells(triangles_aux.shape[0],vtkFacesArray)
        
        self.__vtkData.SetPolys(vtkCellArrayFaces)

     
##        # NO IMPLEMENTADO
##        self.scalarFaces = zeros(tetrahedraMesh.tetrahedra.shape[0], 'float32')
##        if colorMesh:            
##            for i in range(len(tetrahedraMesh.meshes)):
##                self.scalarFaces[tetrahedraMesh.meshes[i]] = i
##        vtkFloatArrayScalars = vtkFloatArray()
##        vtkFloatArrayScalars.SetNumberOfComponents(1)
##        vtkFloatArrayScalars.SetVoidArray(self.scalarFaces, self.scalarFaces.shape[0] ,1)
##
##        self.__vtkData.GetCellData().SetScalars(vtkFloatArrayScalars)
        
    def SetInputPoints(self, points):
        
        self.__vtkData=vtkPolyData()
        
        if points.dtype == 'float32':
            vtkFloatArrayPoints = vtkFloatArray()
        elif points.dtype == 'float64':
            vtkFloatArrayPoints = vtkDoubleArray()            
        vtkFloatArrayPoints.SetNumberOfComponents(3)
        self.__pointsPolyData = points.copy()
        ##pointsNewPolyData[:,0] = pointsNew[:,2]; pointsNewPolyData[:,2] = pointsNew[:,0]+10
        if len(points):
            if points.shape[1]>2:
                self.__pointsPolyData[:,0] = points[:,2]; self.__pointsPolyData[:,2] = points[:,0]
            vtkFloatArrayPoints.SetVoidArray(self.__pointsPolyData, self.__pointsPolyData.shape[0] * points.shape[1] ,1)

        vtkMyPoints = vtkPoints()
        vtkMyPoints.SetData(vtkFloatArrayPoints)

        self.__vtkData.SetPoints(vtkMyPoints)
        
        # Faces
        # Ordenarlos
        self.__faces = concatenate((ones((points.shape[0],1),'int'), arange(points.shape[0]).reshape(-1,1)),1).flatten()    
        vtkFacesArray=vtkIdTypeArray()
        vtkFacesArray.SetVoidArray(self.__faces, len(self.__faces), 1)

        vtkCellArrayFaces=vtkCellArray()
        vtkCellArrayFaces.SetCells(points.shape[0],vtkFacesArray)

        self.__vtkData.SetVerts(vtkCellArrayFaces)
        
##        dict={'vtkData': vtkData, 'Faces':Faces, 'PointsVTK':pointsPolyData}
##        return dict
        
    def SetInputLines(self, points, lines):
        
        if not points == None:            
            # Points
            self.__vtkData = vtkPolyData()

            vtkFloatArrayPoints=vtkFloatArray()
            vtkFloatArrayPoints.SetNumberOfComponents(3)
            self.__pointsPolyData = points.copy()
            ##pointsNewPolyData[:,0] = pointsNew[:,2]; pointsNewPolyData[:,2] = pointsNew[:,0]+10
            if points.shape[1]>2:
                self.__pointsPolyData[:,0] = points[:,2]; self.__pointsPolyData[:,2] = points[:,0]
            vtkFloatArrayPoints.SetVoidArray(self.__pointsPolyData, self.__pointsPolyData.shape[0]*points.shape[1] ,1)

            vtkMyPoints=vtkPoints()
            vtkMyPoints.SetData(vtkFloatArrayPoints)

            self.__vtkData.SetPoints(vtkMyPoints)
            
            # Faces
            self.__faces = []
            for oneLine in lines:
                self.__faces.extend([len(oneLine)])
                self.__faces.extend(oneLine)
            self.__faces = array(self.__faces, 'int')

            vtkFacesArray=vtkIdTypeArray()    
            vtkFacesArray.SetVoidArray(self.__faces, len(self.__faces), 1)

            vtkCellArrayFaces=vtkCellArray()
            vtkCellArrayFaces.SetCells(1,vtkFacesArray)

            self.__vtkData.SetLines(vtkCellArrayFaces)
                        
        else:
            # Faces
            faces = []
            for oneLine in lines:
                faces.extend([len(oneLine)])
                faces.extend(oneLine)
            faces = array(faces, 'int')

            self.__faces = [self.__faces]
            self.__faces.append(faces)
            
            vtkFacesArray=vtkIdTypeArray()
            vtkFacesArray.SetVoidArray(faces, len(faces), 1)

            vtkCellArrayFaces=vtkCellArray()
            vtkCellArrayFaces.SetCells(1,vtkFacesArray)

            self.__vtkData.SetLines(vtkCellArrayFaces)                
        
        ##        dict={'vtkData': vtkData, 'Faces':Faces, 'PointsVTK':pointsPolyData}
        ##        return dict
    def SetInputPolygon(self, points, polygon):
        # Points
        self.__vtkData = vtkPolyData()

        vtkFloatArrayPoints=vtkFloatArray()
        vtkFloatArrayPoints.SetNumberOfComponents(3)
        self.__pointsPolyData = points.copy()
        ##pointsNewPolyData[:,0] = pointsNew[:,2]; pointsNewPolyData[:,2] = pointsNew[:,0]+10
        if points.shape[1]>2:
            self.__pointsPolyData[:,0] = points[:,2]; self.__pointsPolyData[:,2] = points[:,0]
        vtkFloatArrayPoints.SetVoidArray(self.__pointsPolyData, self.__pointsPolyData.shape[0]*points.shape[1] ,1)

        vtkMyPoints=vtkPoints()
        vtkMyPoints.SetData(vtkFloatArrayPoints)

        self.__vtkData.SetPoints(vtkMyPoints)

        # Faces
        self.__faces = zeros(len(polygon)+1, 'int')
        self.__faces[1:] = polygon[:]
        self.__faces[0] = len(polygon)

        vtkFacesArray=vtkIdTypeArray()    
        vtkFacesArray.SetVoidArray(self.__faces, len(self.__faces), 1)

        vtkCellArrayFaces=vtkCellArray()
        vtkCellArrayFaces.SetCells(1,vtkFacesArray)

        self.__vtkData.SetPolys(vtkCellArrayFaces)
        
        ##        dict={'vtkData': vtkData, 'Faces':Faces, 'PointsVTK':pointsPolyData}
        ##        return dict

    def GetVTKPolyData(self):
        if hasattr(self, 'scalarFaces'):
##            dict={'vtkData': self.__vtkData, 'Faces':self.__faces, 'pointsVTK':self.__pointsPolyData, 'scalarFaces':self.scalarFaces, 'scalarPoints':self.__scalarPoints}
            dict={'vtkData': self.__vtkData, 'Faces':self.__faces, 'pointsVTK':self.__pointsPolyData, 'scalarFaces':self.scalarFaces}
        else:
            dict={'vtkData': self.__vtkData, 'Faces':self.__faces, 'pointsVTK':self.__pointsPolyData}
        return dict

    def WritePolyData(self, FileName, format = 'Binary'):
        vtksave = vtkPolyDataWriter ()
        vtksave.SetFileName(FileName)
        vtksave.SetInput(self.__vtkData)
        if format == 'Binary':
            vtksave.SetFileTypeToBinary()
        else:
            vtksave.SetFileTypeToASCII()
        vtksave.Update()
        vtksave.Write()


def WriteMeshToPLY(fileName, mesh):

    PtoVTK = PythonToPolyData()
    PtoVTK.SetInput(mesh)
    vtkMesh = PtoVTK.GetVTKPolyData()

    plyWriter = vtkPLYWriter()
    plyWriter.SetFileName(fileName)
    ##plyWriter.SetFileTypeToASCII()
    plyWriter.SetFileTypeToBinary()
    plyWriter.SetInput(vtkMesh['vtkData'])
    plyWriter.Write()

def WriteTrianglesMeshToOFF(fileName, mesh):
    Points = mesh.points
    Triangles = mesh.triangles

    if fileName[-3:] == 'off':
        fileName = fileName[:-4]
    
    file = open(fileName + '.off','w')

    i=0;

    file.write("OFF\n")
    file.write("%s" %(len(Points)))
    file.write(" %s" %(len(Triangles)))
    file.write(" 0\n")

    while i<len(Points) :
        file.write("%.8f " %(Points[i,2]))
        file.write("%.8f " %(Points[i,1]))
        file.write("%.8f\n" %(Points[i,0]))
        i=i+1	

    i=0
    while i<len(Triangles) :
        file.write("3 ")
        file.write("%i " %(Triangles[i,2]))
        file.write("%i " %(Triangles[i,1]))
        file.write("%i\n" %(Triangles[i,0]))
        i=i+1	

    file.close()
def ReadTrianglesMeshFromOFF(fileName):

    if fileName[-3:] == 'off':
        fileName = fileName[:-4]
        
    # leer los puntos
    fil = file(fileName + '.off', 'r')

    aux = fil.readline().split()
    if aux != 'OFF':
        print 'error de lectura malla OFF'
        return None
    
    aux = fil.readline().split()
    n_points = int(aux[0])
    n_triangles = int(aux[1])
    
    points = zeros((n_points, 3), 'float32')

    for i in range(n_points):
        aux = fil.readline().split()
        points[i,0] = float(aux[1])
        points[i,1] = float(aux[2])
        points[i,2] = float(aux[3])
        
    triangles = zeros((n_triangles, 3), 'int')
    
    for i in range(n_triangles):
        aux = fil.readline().split()
        triangles[i,0] = int(aux[1])
        triangles[i,1] = int(aux[2])
        triangles[i,2] = int(aux[3])
        
    trianglesMesh = TrianglesMesh(points, triangles)
    
    fil.close()
    return trianglesMesh
    
  
def ReadPLY(fileName):
    # solo para mallas de triangulos
    plyReader = vtkPLYReader()
    plyReader.SetFileName(fileName)
    plyReader.Update()
    vtkMesh = plyReader.GetOutput()
    trianglesMesh = TrianglesMesh(ExtractPoints(vtkMesh), ExtractTriangles(vtkMesh))
    return trianglesMesh
    
def FintSimplex(simplex, neighbors,normals,S,sigma,ep,activateComputeParameters,facesSimplex, meshesContour,epContours,Scontour,AngleGammaObjetivo, angleSimplexObjetivo, pointsSelected = None):
    """
    S: tamaño de la vecindad
    sigma: parametro para la acomodacion de la malla
    Scontour < 0 indica que en algulo simplex objetivo es 0
    """

    if pointsSelected != None:
        index_compute_simplex = pointsSelected
    else:
        index_compute_simplex = range(simplex.shape[0])

    neighbors1=simplex[neighbors[:,0],:]
    neighbors2=simplex[neighbors[:,1],:]
    neighbors3=simplex[neighbors[:,2],:]


    #calculo del centro de la esfera que incluye los 4 puntos
    angulos_cero = zeros(simplex.shape[0], 'int8')
    Csphere = zeros(simplex.shape, 'float32')

    A = zeros((4,4), 'float32')
    A[:,3] = 1;

    for i in index_compute_simplex:
        A[0,:3] = simplex[i,:]
        A[1,:3] = neighbors1[i,:]
        A[2,:3] = neighbors2[i,:]
        A[3,:3] = neighbors3[i,:]
        bajo = det(A)
        
        if abs(bajo) > 0.00001:

            A[0,0] = sum(simplex[i,:]**2); A[0,1:3] = simplex[i,1:]
            A[1,0] = sum(neighbors1[i,:]**2); A[1,1:3] = neighbors1[i,1:]
            A[2,0] = sum(neighbors2[i,:]**2); A[2,1:3] = neighbors2[i,1:]
            A[3,0] = sum(neighbors3[i,:]**2); A[3,1:3] = neighbors3[i,1:]
            Csphere[i,0] = det(A) / (2 * bajo)

            A[0,1] = simplex[i,0]; A[0,2] = simplex[i,2]
            A[1,1] = neighbors1[i,0]; A[1,2] = neighbors1[i,2]
            A[2,1] = neighbors2[i,0]; A[2,2] = neighbors2[i,2]
            A[3,1] = neighbors3[i,0]; A[3,2] = neighbors3[i,2]
            Csphere[i,1] = -det(A) / (2 * bajo)

            A[0,1:3] = simplex[i,:2]
            A[1,1:3] = neighbors1[i,:2]
            A[2,1:3] = neighbors2[i,:2]
            A[3,1:3] = neighbors3[i,:2]
            Csphere[i,2] = det(A) / (2 * bajo)
            
        else:
            Csphere[i,:] = 0
            angulos_cero[i]=1


    # calculo del centro del circulo que incluye a los vecinos

    v12=neighbors2-neighbors1
    v13=neighbors3-neighbors1
    v23=neighbors3-neighbors2
    a=sqrt(sum(v23**2,1)) #1
    b=sqrt(sum(v13**2,1)) #2
    c=sqrt(sum(v12**2,1)) #3
    C1=(a**2)*(b**2 + c**2 - a**2)
    C2=(b**2)*(a**2 + c**2 - b**2)
    C3=(c**2)*(a**2 + b**2 - c**2) 
    ww=C1+C2+C3
    ww[ww == 0] = 0.00001 # ajuste por problemas numericos
    C1=C1/ww
    C2=C2/ww
    C3=C3/ww
    Ccircle=neighbors1*C1.reshape(-1,1)+neighbors2*C2.reshape(-1,1)+neighbors3*C3.reshape(-1,1)


    # calculo del radio del circulo y del radio de la esfera

    Rsphere = sqrt(sum((Csphere - neighbors1)**2, 1))
    Rcircle = sqrt(sum((Ccircle - neighbors1)**2, 1))


    # calculo de angulos simplex    
    aux = ( Rcircle/Rsphere ) * ((sum((neighbors1-simplex)*normals,1)<0)*2-1)
    aux[aux > 1] = 1 # ajuste por problemas numericos
    aux[aux > -1] = -1 # ajuste por problemas numericos
    angleSimplex = arcsin(  ( Rcircle/Rsphere ) * ((sum((neighbors1-simplex)*normals,1)<0)*2-1) )
    angleSimplex[(Rcircle>Rsphere).nonzero()]=pi/2
    angleSimplex[angulos_cero.nonzero()] = 0


    # calculo de la proyeccion de los puntos en el plano
    L = sum((simplex-neighbors1)*normals,1) # producto punto
    proyeccion = simplex-(normals*L.reshape(-1,1));

    # #########
    #calculo de angulo para contraccion continua
    ponderado=0  #1: para ponderar la contribucion de los vecinos por su distancia 0: para solo promediar
    vecindad = zeros((2,1+sum((ones(S)*2)**arange(0,S)*3)), 'int') #el maximo tamaño que podria ocurrir    
    count_vecindad=0

    ##vecindad_aux = zeros((2**(S-1))*3)
    ##Ivecindad_aux = ones((2**(S-1))*3)
    ##count_vecindad_aux=0

    angulo_objetivo = zeros(simplex.shape[0], 'float32')

    if all(angleSimplexObjetivo < 0):
        for i in index_compute_simplex:
        # #################### calculo de vecindad #########################
            count_vecindad = 0
            vecindad[:] = -1; vecindad[0,0] = i; vecindad[1,0] = 0
            count_vecindad = 1

            vecindad_aux = i
            ##count_vecindad_aux = 3
            for numVecindad in range(S):
                vecindad_aux = neighbors[vecindad_aux,:].reshape((-1))
                vecindad_aux = set(vecindad_aux)
                vecindad_aux.difference_update(vecindad[0,:])
                vecindad_aux = list(vecindad_aux)
                
                Lvecindad = len(vecindad_aux)
                vecindad[0,count_vecindad:count_vecindad+Lvecindad] = vecindad_aux
                vecindad[1,count_vecindad:count_vecindad+Lvecindad] = numVecindad+1
                count_vecindad += Lvecindad 
                
            if ponderado:
                #falta
                angulo_objetivo[i]=  arcsin(Rcircle[i] * sum( sin(angleSimplex[vecindad[0,:count_vecindad]]) / Rcircle[vecindad[0,:count_vecindad]] ) * (1. / count_vecindad)   ) 
            else:
                ponderacion = Rcircle[i] / Rcircle[vecindad[0,:count_vecindad]]
                ponderacion = ponderacion / ponderacion.sum()           
                angulo_objetivo[i]=  arcsin( sum( sin(angleSimplex[vecindad[0,:count_vecindad]]) * ponderacion ) )
    else:
        angulo_objetivo[:] = angleSimplexObjetivo


    proyeccion_objetivo = (ep[:,0].reshape(-1,1) * neighbors1) + ( ep[:,1].reshape(-1,1) * neighbors2) + (ep[:,2].reshape(-1,1) * neighbors3)
    d=sqrt(sum((Ccircle - proyeccion_objetivo)**2,1))
    Rcircle2=Rcircle**2

    L_objetivo= ((Rcircle2-(d**2))*tan(angulo_objetivo)) / (  sqrt( Rcircle2 + ( (Rcircle2-(d**2)) * (tan(angulo_objetivo)**2) ) )+Rcircle)

    # calculo de la fuerza interna
    F_tangencial = proyeccion_objetivo - proyeccion
    F_normal = normals * (L_objetivo-L).reshape(-1,1)
    Fint = F_tangencial + F_normal
    

    # ##### calculo para los contornos #########
    Ncontours = len(meshesContour)
    if len(meshesContour)>0:
        contours = []
        epContoursT = []
        for i in range(len(meshesContour)):
            contours.append(meshesContour[i])
            epContoursT.append(epContours[i])
        
        T = []; B = []; Q = []; R = []; RT = []; angleSimplexContours = []; AngleGamma = []
        Ncontours = len(contours)
        for i in range(Ncontours):
            concatenate((contours[i][1:], contours[i][:1]), 0)
            Pi = simplex[contours[i]].copy()
            Pimas1 = concatenate((simplex[contours[i][1:]],simplex[contours[i][:1]]),0)
            Pimenos1 = concatenate((simplex[contours[i][-1:]],simplex[contours[i][:-1]]),0)
            Pimas2 = concatenate((simplex[contours[i][2:]],simplex[contours[i][:2]]),0)
            Pimenos2 = concatenate((simplex[contours[i][-2:]],simplex[contours[i][:-2]]),0)
            
            Pimas1MenPimenos1 = Pimas1 - Pimenos1
            Dcontours = sqrt(sum((Pimas1MenPimenos1*Pimas1MenPimenos1),1)).reshape(-1,1)
            T.append(Pimas1MenPimenos1/Dcontours)
            
            PiMenPimenos1 = Pi - Pimenos1
            B.append(cross(PiMenPimenos1, Pimas1-Pi))
            Bnorm = sqrt(sum((B[-1]*B[-1]),1)).reshape(-1,1)
            B[-1] = B[-1]/Bnorm
            B0 = (Bnorm==0).nonzero()
            
            # si B es 0 fijarlo de manera que la curvatura se produzca en el mismo plano de la malla, para eso fijo B utlizando la diraccion del vecino al punto pero que no pertenece al contorno, de manara que sea ortogonal a esa direccion y a T
            lContour = len(contours[i])
            for val0 in B0[0]:
                neighborsVal0 = neighbors[contours[i][val0]]
                neighborsInMesh = neighborsVal0[((neighborsVal0!=contours[i][val0-1]) * (neighborsVal0!=contours[i][val0+1-lContour])).nonzero()][0] # encontrar el vecino del punto de contorno, pero que se encuentra dentro de la malla, no en el contorno 
                B[-1][val0] = cross(simplex[neighborsInMesh] - Pi[val0], T[-1][val0])
                B[-1][val0] /= norm(B[-1][val0])
            
##            tetay = 1.; tetax = 0.; tetaz = 0.     
##            Rotp=array([[cos(tetay)*cos(tetaz), -cos(tetay)*sin(tetaz), -sin(tetay)],
##            [-sin(tetax)*sin(tetay)*cos(tetaz)+cos(tetax)*sin(tetaz), sin(tetax)*sin(tetay)*sin(tetaz)+cos(tetax)*cos(tetaz), -sin(tetax)*cos(tetay)],
##            [cos(tetax)*sin(tetay)*cos(tetaz)+sin(tetax)*sin(tetaz), -cos(tetax)*sin(tetay)*sin(tetaz)+sin(tetax)*cos(tetaz), cos(tetax)*cos(tetay)]], 'float32')
##            Rotp = Rotp.T
##                
##            if len(B0[0]) > 0: # por si hay 2 trozos de contorno iguales, entonces el producto crus es 0, hay que elegir un vector perpendicular a ellos.
##                B[-1][B0[0],:] = cross(PiMenPimenos1[B0[0],:], dot(PiMenPimenos1[B0[0],:],Rotp) ) # a lo mejor a una manera mejor
##                B[-1][B0[0],:] = B[-1][B0[0],:]/sqrt(sum((B[-1][B0[0],:]*B[-1][B0[0],:]),1)).reshape(-1,1)
                
##            Q.append(cross(B[-1],T[-1])) # Normal
            Q.append(cross(T[-1],B[-1])) # Normal
            
            Pimenos1MenPimenos2 = Pimenos1-Pimenos2
            R.append(cross(T[-1], cross(Pimenos1MenPimenos2,Pimas2-Pimas1)))
##            R.append(cross(cross(Pimenos1MenPimenos2,Pimas2-Pimas1),T[-1]))
            Rnorm = sqrt(sum((R[-1]*R[-1]),1)).reshape(-1,1)
            R[-1] = R[-1]/Rnorm
            R0 = (Rnorm==0).nonzero()


            # si R es 0 fijarlo de manera que siga las mismas caracteristicas que B 
            for val0 in R0[0]:
                R[-1][val0] = cross(T[-1][val0], B[-1][val0])
                R[-1][val0] /= norm(R[-1][val0])
                
##            if len(R0[0]) > 0: # por si hay 2 trozos de contorno iguales, entonces el producto crus es 0, hay que elegir un vector perpendicular a ellos.
##                R[-1][R0[0],:] = cross(cross(Pimenos1MenPimenos2[R0[0],:], dot(Pimenos1MenPimenos2[R0[0],:],Rotp)) , T[-1][R0[0],:]) # a lo mejor a una manera mejor
##                R[-1][R0[0],:] = R[-1][R0[0],:]/sqrt(sum((R[-1][R0[0],:]*R[-1][R0[0],:]),1)).reshape(-1,1)
             


            
##            aux1 = Pimenos1 - Pi; aux1 = aux1/sqrt(sum((aux1*aux1),1)).reshape(-1,1) # por arriba
##            aux2 = Pimas1 - Pi; aux2 = aux2/sqrt(sum((aux2*aux2),1)).reshape(-1,1)
##            angleSimplexContours.append(arccos(sum(aux1*aux2,1).reshape(-1,1)))
            
            aux1 = Pi - Pimenos1; aux1 = aux1/sqrt(sum((aux1*aux1),1)).reshape(-1,1) # por el lado
            aux2 = Pimas1 - Pi; aux2 = aux2/sqrt(sum((aux2*aux2),1)).reshape(-1,1)
            aux1Dotaux2 = (aux1*aux2).sum(1)
            aux1Dotaux2[(aux1Dotaux2 > 1.).nonzero()] = 1.
            aux1Dotaux2[(aux1Dotaux2 < -1.).nonzero()] = -1.
            angleSimplexContours.append(arccos(aux1Dotaux2.reshape(-1,1)))
            
            QdotR = sum(Q[-1]*R[-1],1).reshape(-1,1)
            QdotR[(QdotR>1).nonzero()] = 1.
            QdotR[(QdotR<-1).nonzero()] = -1.
            RT.append(cross(T[-1],R[-1]))
            AngleGamma.append(arccos(QdotR))
            QdotRT = (Q[-1] * RT[-1]).sum(1)
            index_aux = (QdotRT < 0).nonzero()[0]
            AngleGamma[-1][index_aux] = - AngleGamma[-1][index_aux]
            
##            AngleGammaObjetivo = AngleGamma[-1].copy() #cuidado
            
            proyeccionContours = Pimenos1 + T[-1] * sum(T[-1] * PiMenPimenos1,1).reshape(-1,1)
            proyeccion_objetivoContours = epContoursT[i][:,:1]*Pimenos1 + epContoursT[i][:,1:2]*Pimas1
            
            # calculo del angulo promedio
            angulo_objetivoContours = angleSimplexContours[-1].copy()
            if Scontour > 0:
                for k in range(Scontour+1)[1:]:
                    angulo_objetivoContours += concatenate((angleSimplexContours[-1][k:],angleSimplexContours[-1][:k]),0)
                    angulo_objetivoContours += concatenate((angleSimplexContours[-1][-k:],angleSimplexContours[-1][:-k]),0)
                angulo_objetivoContours = angulo_objetivoContours / (Scontour*2.+1)
            elif Scontour < 0:
                angulo_objetivoContours = 0
            
            Rcontours = Dcontours/2.; Rcontours2 = Rcontours * Rcontours
##            dcontours = Rcontours * (epContoursT[i][:,:1] - 0.5); dcontours2 = dcontours * dcontours
##            dcontours = proyeccionContours - ((Pimas1 + Pimenos1) / 2.)  ; dcontours2 = sum(dcontours * dcontours,1).reshape(-1,1)
            dcontours = proyeccion_objetivoContours - ((Pimas1 + Pimenos1) / 2.)  ; dcontours2 = sum(dcontours * dcontours,1).reshape(-1,1)
            L_objetivoContours = ((Rcontours2-dcontours2)*tan(angulo_objetivoContours)) / (  sqrt( Rcontours2 + ( (Rcontours2-dcontours2) * (tan(angulo_objetivoContours)**2) ) )+Rcontours)
            
##            L_Contours = ((Rcontours2-dcontours2)*tan(angleSimplexContours[-1])) / (  sqrt( Rcontours2 + ( (Rcontours2-dcontours2) * (tan(angleSimplexContours[-1])**2) ) )+Rcontours)
            L_Contours = (((Pi - proyeccionContours) * Q[-1]).sum(1)).reshape(-1,1)
##            L_Contours2 = sqrt(sum(L_Contours2*L_Contours2,1).reshape(-1,1))
##            L_Contours = L_Contours2
            
            # calculo de la fuerza interna
            F_tangencialContours = proyeccion_objetivoContours - proyeccionContours
##            F_tangencialContours = 0
##            F_tangencialContours = (Pimenos1 - Pi) * epContoursT[i][:,:1] +   (Pimas1 - Pi) * epContoursT[i][:,1:2]
##            F_normalContours = (L_objetivoContours * cos(AngleGammaObjetivo) ) * R[-1] + (L_objetivoContours * sin(AngleGammaObjetivo) ) * cross(T[-1],R[-1])
            
            F_normalContours = (L_objetivoContours * cos(AngleGammaObjetivo) - L_Contours * cos(AngleGamma[-1]) ) * R[-1] + (L_objetivoContours * sin(AngleGammaObjetivo) - L_Contours * sin(AngleGamma[-1]) ) * RT[-1]
##            F_normalContours = 0.
            
            
            FintContours = F_tangencialContours + F_normalContours
            Fint[contours[i],:] = FintContours

    if 0: #dibujos        
        points_dibujo = array([], 'float32')
        points_dibujo = points_dibujo.reshape(0,3)
        points_dibujo_R = array([], 'float32')
        points_dibujo_R = points_dibujo_R.reshape(0,3)
        points_dibujo_TR = array([], 'float32')
        points_dibujo_TR = points_dibujo_TR.reshape(0,3)
        points_dibujo_Q = array([], 'float32')
        points_dibujo_Q = points_dibujo_Q.reshape(0,3)
        
        lineas_dibujo = []
        for i, aux in  enumerate(meshesContour):
            LL = len(aux)
            points_dibujo = concatenate((points_dibujo, simplex[aux]), 0)
            points_dibujo_R = concatenate((points_dibujo_R, simplex[aux] + R[i]), 0)
            points_dibujo_TR = concatenate((points_dibujo_TR, simplex[aux] + cross(T[i],R[i])), 0)
            points_dibujo_Q = concatenate((points_dibujo_Q, simplex[aux] + Q[i]), 0)
        LL = points_dibujo.shape[0]
        lineas_dibujo_R = []
        lineas_dibujo_RT = []
        lineas_dibujo_Q = []
        for j in range(points_dibujo.shape[0]):
            lineas_dibujo_R.append([j, j+LL])
            lineas_dibujo_RT.append([j, j+2*LL])
            lineas_dibujo_Q.append([j, j+3*LL])

        points_dibujo_all = concatenate((points_dibujo, points_dibujo_R, points_dibujo_TR, points_dibujo_Q), 0)
        # ############### para Visualizar
        PtoVTK = PythonToPolyData()
        PtoVTK.SetInputLines(points_dibujo_all, lineas_dibujo_R)
        PtoVTK.WritePolyData('..\\R.vtk')

        PtoVTK = PythonToPolyData()
        PtoVTK.SetInputLines(points_dibujo_all, lineas_dibujo_RT)
        PtoVTK.WritePolyData('..\\RT.vtk')

        PtoVTK = PythonToPolyData()
        PtoVTK.SetInputLines(points_dibujo_all, lineas_dibujo_Q)
        PtoVTK.WritePolyData('..\\Q.vtk')
        
        # ###############
        

    return Fint
                    

def NearPointPointTriangle(point, triangle):
    B = triangle[0]
    e0 = triangle[1] - B
    e1 = triangle[2] - B
    a = dot(e0, e0)
    b = dot(e0, e1)
    c = dot(e1, e1)
    B_P = B - point # d
    d = dot(e0, B_P)
    e = dot(e1, B_P)
    f = dot(B_P, B_P)

    bb = a * c - b * b
    s = (b * e - c * d) / bb
    t = (b * d - a * e) / bb
    
    if s>=0 and t>=0 and (s+t)<=1: # el punto mas cercano es la proyeccion ortogonal sobre el triangulo
        None
    else:
        dett = a * c - b * b
        s = (b * e - c * d)
        t = (b * d - a * e)
        if (s + t) <= dett:
            if s < 0:
                if t < 0:
                    region = 4
                else:
                    region = 3
            elif t < 0:
                region = 5
            else:
                region = 0
        else:
            if  s < 0:
                region = 2
            elif t < 0:
                region = 6
            else:
                region = 1
                
            
        if region == 0:
            invdett = 1./dett
            s *= invdett
            t *= invdett
        elif region == 1:
            numer = c + e - b - d
            if numer <= 0:
                s = 0
            else:
                denom = a - 2 * b + c
                if numer >= denom:
                    s = 1
                else:
                    s = numer / denom
            t =  1 - s

        elif (region == 5):
            t = 0
            if d >= 0:
                s = 0
            else:
                if -d >= a:
                    s = 1
                else:
                    s = -d / a
                    
        elif (region == 3):
            s = 0
            if e >= 0:
                t = 0
            else:
                if -e >= c:
                    t = 1
                else:
                    t = -e / c
                  
        elif (region == 2):
            tmp0 = b + d
            tmp1 = c + e
            if tmp1 > tmp0:
                numer = tmp1 - tmp0
                denom = a - 2 * b + c
                if numer >= denom:
                    s = 1
                else:
                    s = numer / denom
                t = 1 - s
            else:
                s = 0
                if tmp1 <= 0:
                    t = 1
                else:
                    if e >= 0:
                        t = 0
                    else:
                        t = -e / c

        elif (region == 6):
            tmp0 = b + e
            tmp1 = a + d
            if tmp1 > tmp0:
                numer = tmp1 - tmp0
                denom = - a - 2 * b + c
                if numer >= denom:
                    t = 1
                else:
                    t = numer / denom
                s = 1 - t
            else:
                t = 0
                if tmp1 <= 0:
                    s = 1
                else:
                    if d >= 0:
                        s = 0
                    else:
                        s = -d / a

        elif (region == 4):
            if d < 0.:
                t = 0
                if d > 0:
                    s = 0
                else:
                    if a + d <= 0:
                        s = 1
                    else:
                        s = -d/a
            else:
                s = 0
                if e > 0:
                    t = 0
                else:
                    if c + e <= 0:
                        t = 1
                    else:
                        t = -e/c
            
                        
    return B + s * e0 + t * e1


def contrast3D (Din, min, max):
    minIn = Din.min()
    maxIn = Din.max()

    alfa=(max-min)/(maxIn-minIn)
    beta=max-(alfa*maxIn)

    Dout=(Din*alfa)+beta
    return Dout


def SaveData(FileName, DataName, Data ):
    if len(DataName) != len(Data):
        print 'El numero de nombres no coincide con el numero de datos'
        return 0
    f = gzip.GzipFile(FileName, 'wb')
    f.write(cPickle.dumps((DataName,Data),1))    
    f.close()
    return 1
    

def ReadData(FileName):
    f = gzip.GzipFile(FileName, 'rb')
    buffer = ""
    while 1:
        data = f.read()
        if data == "":
            break
        buffer += data
    tupla = cPickle.loads(buffer)
    f.close()    
    return tupla
    # para asignarlos nombres de ArchiveName a los datos se debe agregar esto en el codigo despues de utilizar ReadData   
##    script=''
def meshgrid2(*arrs):
    arrs = tuple(reversed(arrs))  #edit
    lens = map(len, arrs)
    dim = len(arrs)

    sz = 1
    for s in lens:
        sz*=s

    ans = []    
    for i, arr in enumerate(arrs):
        slc = [1]*dim
        slc[i] = lens[i]
        arr2 = asarray(arr).reshape(slc)
        for j, sz in enumerate(lens):
            if j!=i:
                arr2 = arr2.repeat(sz, axis=j) 
        ans.append(arr2)

    return tuple(ans)

def InPolyedron(p,t,tnorm,qp):
    ##% InPolyedron detects points inside a manifold closed surface. The normal
    ##% orientation is assumed to be outward. Warning there is no check for
    ##% closed surface or normal orientation. The v3 version, differently from v1
    ##% and v2 supports points lying on edges or point. In the limit of numerical
    ##% accuracy, points lying on the surface will be considered in.
    ##% 
    ##% 
    ##% Input:
    ##% 
    ##% p: POints of the surface npx3 array
    ##% 
    ##% t:triangles indexes, first points flagged as one, ntx3 array
    ##% 
    ##% tnorm: outward normals of traingles, ntx3 array
    ##% 
    ##% qp: points to be queried, nqx3 array
    ##% 
    ##% 
    ##% Output:
    ##% 
    ##% in: nqx1 logical vector, true for in, false for out
    ##% 
    ##% 
    ##% Author: Giaccari Luigi 
    ##% Created: 15/05/2009
    ##% e-mail: giaccariluigi@msn.com
    ##% 
    ##% Visti my website: http://giaccariluigi.altervista.org/blog/
    ##% 
    ##% This work is free thanks to users gratitude. If you find it usefull
    ##% consider making a donation on my website.
    ##
    ##
    ##
    ##
    ##%possibili migliorie:
    ##% - Effettuare il test con il numero di intersioni pari o dispari per la
    ##% maggiorparte dei punti e passare all'attuale test per i casi dubbi.
    ##
    ##
    def GetBox2tMap(p,t,tnorm,rect,k):
    ##    %rect contiene gli estrameli del rettagnolo di divisione.
    ##    %k è il fattore di divisione delle boxes, per k=1 il numero di scatole sarà
    ##    %uguale al numero di punti. k=2 qudruplicherà il numero delle scatole

        firsttoll = 1e-9;
        secondtoll = 1e-10;


        # size parameters
        np = p.shape[0]
        nt = t.shape[0]



        # bounding box analysis
        minx = rect[0] #min(p(:,0));
        miny = rect[1] #min(p(:,1));

        maxx = rect[2] #max(p(:,1));
        maxy = rect[3] #max(p(:,2));

        A = (maxx-minx)*(maxy-miny)
        step = sqrt(A/(np*k)) #step quasi square

        nx = int(floor((maxx-minx)/step)) #nota: non cambiare floor con ceil
        ny = int(floor((maxy-miny)/step))


        if nx==0: #check thin mapping
            px = firsttoll+maxx-minx
            nx=1
        else:
            px=(maxx-minx+firsttoll)/nx #eps per aumentare il passo

        if ny==0: #check thin mapping
            py=firsttoll+(maxy-miny)
            ny=1
        else:
            py=(maxy-miny+firsttoll)/ny


    ##    %nota 1e-9 per aumentare il passo in modo da evitare nel passo successivo u
    ##    %indice più grande del numero delle scatole


        N = nx*ny #real Boxes number

        Box2tMap = list([] for i in range(N))
        BoxesId = zeros((np,2), 'int32')

    ##    %first loop to count points and get ceiling round
        for i in range(np):
            idx = ceil((p[i,0]-minx+secondtoll)/px)
            idy = ceil((p[i,1]-miny+secondtoll)/py)
    ##    %     id=ny*(idx-0)+idy
            BoxesId[i,0]=idx
            BoxesId[i,1]=idy
             
        c = zeros((N,1),'int32') - 1

        #now loop trough all triangles a build the map

        for i in range(nt):
             
            if abs(tnorm[i,2])<0.00000000001: #jump vertical triangles
                continue
            
            p1 = t[i,0]
            p2 = t[i,1]
            p3 = t[i,2]
            
            #get extreme values 
            idxmax = max(BoxesId[[p1,p2,p3],0]);
            idxmin = min(BoxesId[[p1,p2,p3],0]);
            idymax = max(BoxesId[[p1,p2,p3],1]);
            idymin = min(BoxesId[[p1,p2,p3],1]);
            
            # loop trough all boxes that may contains the triangle
            for idx in range(idxmin,idxmax + 1):
                for idy in range(idymin,idymax + 1):
                   
                    #get boxes id and increase counter
                    id = ny*(idx-1) + idy - 1
                    c[id] = c[id]+1
                    Box2tMap[id].append(i) #insert traingle into map

        return Box2tMap


    def SortCounterclockwise(t,x,y):
    # SortCounterclockwise
        #get points coordinate vectors
        x1 = x[t[:,0]]
        x2 = x[t[:,1]]
        x3 = x[t[:,2]]
        y1 = y[t[:,0]]
        y2 = y[t[:,1]]
        y3 = y[t[:,2]]


        cx = (x1+x2+x3)/3
        cy=(y1+y2+y3)/3 #centroid
        del x3,y3

        v1x = x1-cx
        v1y = y1-cy

        v2x = x2-cx
        v2y = y2-cy


        cp = (v1x*v2y - v1y*v2x) < 0 #fails cross product criterion (magnitud de coordenada z del producto cruz)

        t[cp,:] = t[cp,:][:,[1,0,2]] #%get counterclockwise orientation


        return t






    # errors check

    [m,n] = p.shape
    if n !=3:
        print 'Wrong points dimension'

    [m1,n] = t.shape
    if n!=3:
        print 'Wrong t dimension'

    [m2,n] = tnorm.shape
    if n!=3:
        print 'Wrong tnorm dimension'

    if m1 != m2:
        print 't dismatch tnorm dimensions'

    [m,n]= qp.shape
    if n != 3:
        print 'Wrong qp dimension'


    # internal parameters

    k = 1 #boxes subdivison factor
    #future improvement find anoptimal value depending on t and qp size



    #rectangle for box2triangle map
    rect = zeros(4)
    rect[0] = p[:,0].min()
    rect[1] = p[:,1].min()

    rect[2] = p[:,0].max()
    rect[3] = p[:,1].max()

    firsttoll = 1e-9;
    secondtoll = 1e-10;

    #get size data

    np = p.shape[0]
    ## nt=size(t,1);
    nq = qp.shape[0]


    #Sort counterclockwise
    t = SortCounterclockwise(t,p[:,0],p[:,1])


    # Box2tMap
    Box2tMap = GetBox2tMap(p,t,tnorm,rect,k)

    inn = zeros((nq,1), 'int8')

    #get box reference
    minx = rect[0] #min(p(:,1));
    miny = rect[1] #min(p(:,2));

    maxx = rect[2] #max(p(:,1));
    maxy = rect[3] #%max(p(:,2));

    A = (maxx-minx)*(maxy-miny)
    step = sqrt(A/(np*k)) #step quasi square
    nx = int(floor((maxx-minx)/step))
    ny = int(floor((maxy-miny)/step))

    if nx == 0: #check thin mapping
        px = firsttoll+maxx-minx
        nx = 1
    else:
        px = (maxx-minx+firsttoll)/nx #eps per aumentare il passo

    if ny==0:  #check thin mapping
        py = firsttoll+(maxy-miny)
        ny=1
    else:
        py = (maxy-miny+firsttoll)/ny

    #loop trough all query points

    n = zeros((2,1))  #edge normal
    for i in range(nq):
        #make temp scalar
        x = qp[i,0]
        y = qp[i,1]
        z = qp[i,2]
        
        #get box coordinates
        idx = int(ceil((x-rect[0]+secondtoll)/px))
        if idx<1 or idx>nx:
            continue #points is outside
            
        idy = int(ceil((y-rect[1]+secondtoll)/py))
        if idy<1 or idy>ny:
            continue #points is outside
        
        id = ny*(idx-1)+idy - 1
        
        #get mapped triagnles
        ttemp = Box2tMap[id]

        #loop trough all triangles   
        mindist = inf
        N = 1
        for j in range(len(ttemp)):
            idt = ttemp[j]
            p1 = t[idt,0]
            p2 = t[idt,1]
            p3 = t[idt,2]

            #run inside triangle test


            #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # edge1
            #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            n[0] = -p[p2,1] + p[p1,1]
            n[1]=p[p2,0]-p[p1,0] #normals to triangle edge
            test = n[0]*x+n[1]*y-n[0]*p[p1,0] - n[1]*p[p1,1];

            #debug
    ##%        close(figure(1));
    ##%        figure(1)
    ##%        hold on
    ##%        plot(p([p2,p1],1),p([p2,p1],2),'r-')
    ##%        plot(qp(i,1),qp(i,2),'g*')
    ##%        

            if test < 0: #mettendo l'uguale si escludono i punti sulla superficie
                continue #test failed pints is outside of the triangle

            #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # edge2
            #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            n[0] = -p[p3,1]+p[p2,1]
            n[1] = p[p3,0]-p[p2,0] #normals to triangle edge
            test = n[0]*x+n[1]*y-n[0]*p[p2,0]-n[1]*p[p2,1]

    ##       %debug
    ##%        figure(1)
    ##%        hold on
    ##%        plot(p([p3,p2],1),p([p3,p2],2),'r-')
    ##%      
           
            if test < 0: #mettendo l'uguale si escludono i punti sulla superficie
                continue #%test failed pints is outside of the triangle

            #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            #% edge3
            #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           
            n[0] = -p[p1,1]+p[p3,1]
            n[1]=p[p1,0]-p[p3,0] #normals to triangle edge
            test = n[0]*x+n[1]*y-n[0]*p[p1,0]-n[1]*p[p1,1];

    ##       %debug
    ##%        figure(1)
    ##%        hold on
    ##%        plot(p([p1,p3],1),p([p1,p3],2),'r-')
    ##%      
           
            if test < 0: #mettendo l'uguale si escludono i punti sulla superficie
                continue #%test failed pints is outside of the triangle

           
    ##       %debug
    ##%        close(figure(2));
    ##%        figure(2);
    ##%        hold on
    ##%        axis equal
    ##%        trisurf(t(idt,:),p(:,1),p(:,2),p(:,3),'facecolor','c');
    ##%        cc=(p(t(idt,1),:)+p(t(idt,2),:)+p(t(idt,3),:))/3;
    ##%        quiver3(cc(1),cc(2),cc(3),tnorm(idt,1),tnorm(idt,2),tnorm(idt,3))
    ##%        plot3(x,y,z,'g*');
    ##%        
    ##%    end
       
       
       
       #run semispace test with closest triangle

            n1 = tnorm[idt,0]
            n2 = tnorm[idt,1]
            n3 = tnorm[idt,2]
            p1 = t[idt,0]
            d = -n1*p[p1,0]-n2*p[p1,1]-n3*p[p1,2]

            #run distance test triangle test

            dist = z-(-n1*x-n2*y-d)/n3 #distance along ray
            if abs(dist) < abs(mindist):
                mindist = dist
                N = n3

        inn[i] = mindist*N <= 0
    return inn
def SaveImage3DasVTK(image, fileName, spacing = [1.,1.,1.]):
    vtkImage3DImport = vtkImageImportFromArray()
    vtkImage3DImport.SetArray(image)    
    vtkImage3D = vtkImage3DImport.GetOutput()
    vtkImage3D.Update()
    vtkImage3D.SetSpacing(spacing)

    vtkImageDataWrite = vtkStructuredPointsWriter()
    vtkImageDataWrite.SetFileName(fileName)
    vtkImageDataWrite.SetInput(vtkImage3D)
    vtkImageDataWrite.SetFileTypeToBinary()
    vtkImageDataWrite.Write()