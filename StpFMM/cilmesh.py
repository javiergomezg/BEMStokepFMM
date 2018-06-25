import numpy as np 
from numpy import pi, sin, cos
from math import atan2
from matplotlib import pyplot, cm
from matplotlib import rcParams



def geometria_cil(R, n, h, cx, cy):
    #### preambulos, definicion de arreglos para los nodos
    ### R: radio circunferencia
    ### n: numero de paneles en el que esta dividido el cilindro
    ### h: diferencial de angulo
    ### cx, cy: centro de cilindro
    
    ### X, Y_ pan: nodos del panel
    X_pan = np.zeros(n+1, dtype = np.float64)
    Y_pan = np.zeros(n+1, dtype = np.float64)
    
    ### X, Y_col: puntos de colocacion del panel
    X_col = np.zeros(n, dtype = np.float64)
    Y_col = np.zeros(n, dtype = np.float64)
    ### R_pan: longitud del panel
    R_pan = np.zeros(n, dtype = np.float64)
    #theta_pan = np.zeros(n)
    
    ### se definen los nodos de cada panel de la circunferencia
    nn = np.arange(0, n, dtype = np.float64)
    X_pan[:-1] = R*cos(nn[:]*h) + cx
    Y_pan[:-1] = R*sin(nn[:]*h) + cy
    X_pan[-1] = X_pan[0]
    Y_pan[-1] = Y_pan[0]
    del nn
    
    ### El ciclo for ubica los nodos en la circunferencia
        
    #for i in range(n+1):
        #X_pan[i] = R*cos(i*h) + cx
        #Y_pan[i] = R*sin(i*h) + cy
        #print X_pan[i], Y_pan[i]
    
    ### Las siguiente lineas son para ubicar el punto de colocacion
    ### Al centro del panel que corresponda
    X_col[0], Y_col[0] = 0.5*(X_pan[1]+X_pan[0]), 0.5*(Y_pan[1]+Y_pan[0])
    X_col[1:n-1] = 0.5*(X_pan[2:n]+X_pan[1:n-1])  
    Y_col[1:n-1] = 0.5*(Y_pan[2:n]+Y_pan[1:n-1])
    X_col[n-1], Y_col[n-1] = 0.5*(X_pan[0]+X_pan[n-1]), 0.5*(Y_pan[0]+Y_pan[n-1])
    

    
    #for i in range(len(theta_pan)):
    #    theta_pan[i] = atan2((Y_pan[i+1]-Y_pan[i]),(X_pan[i+1]-X_pan[i]))
    
    ### se define la longitud de cada panel
    R_pan[:] = ((X_pan[1:]-X_pan[:-1])**2 + (Y_pan[1:]-Y_pan[:-1])**2)**0.5 

    
    #for i in range(n_pan-1):
    #    print X_pan[i], Y_pan[i], X_col[i], Y_col[i]
    
    #XX, YY = np.empty(5), np.empty(5)
    
    #XX[0], YY[0] = -1.2*R + cx , 1.2*R + cy
    #XX[1], YY[1] = 1.2*R + cx , 1.2*R + cy
    #XX[2], YY[2] = 1.2*R + cx , -1.2*R + cy
    #XX[3], YY[3] = -1.2*R + cx , -1.2*R + cy
    #XX[4], YY[4] = XX[0], YY[0]
        
    
#    fig2 = pyplot.figure(num = None, figsize = (8, 6), dpi = 80, 
#                     facecolor = 'w', edgecolor = 'k')
#    pyplot.plot(X_pan,Y_pan, 'go--', linewidth = 3)
#    pyplot.plot(X_col,Y_col, 'ro', linewidth = 3)
#    pyplot.plot(XX,YY, 'bo-', linewidth = 3)
#    pyplot.grid(True)
#    pyplot.axis([-1.3*(R+abs(cx)),1.3*(R+abs(cx)),-1.3*(R+abs(cy)),1.3*(R+abs(cy))])
#    pyplot.grid(color = '0.3', linestyle = '--', linewidth = 0.3)
#    pyplot.savefig('circulito.png')
#    pyplot.show()
    
    return X_col, Y_col, X_pan, Y_pan, R_pan
    

def normal_vector2(X, Y, h):
    rx = np.zeros(len(X)-1, dtype = np.float64)
    ry = np.zeros(len(Y)-1, dtype = np.float64)
    
    nx = np.zeros(len(X)-1, dtype = np.float64)
    ny = np.zeros(len(Y)-1, dtype = np.float64)    
    
    rx[:] = (X[1:]-X[:-1])/h[:]
    ry[:] = (Y[1:]-Y[:-1])/h[:]
    
    nx[:] = -ry[:]
    ny[:] =  rx[:]
    
    return nx, ny

def unit_vector1(X, Y, h):
    tx = np.zeros(len(X)-1, dtype = np.float64)
    ty = np.zeros(len(Y)-1, dtype = np.float64)
    
    nx = np.zeros(len(X)-1, dtype = np.float64)
    ny = np.zeros(len(Y)-1, dtype = np.float64)    
    
    tx[:] = (X[:-1]-X[1:])/h[:]
    ty[:] = (Y[:-1]-Y[1:])/h[:]
    
    nx[:] = -ty[:]
    ny[:] =  tx[:]
    
    return tx, ty, nx, ny

def tan_vetor(nx,ny):
    tx = np.zeros(len(nx), dtype = np.float64)
    ty = np.zeros(len(nx), dtype = np.float64)
    
    ty[:] = -ny[:]
    tx[:] =  nx[:]
    
    return tx, ty

