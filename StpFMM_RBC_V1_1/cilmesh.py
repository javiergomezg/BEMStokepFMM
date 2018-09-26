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


#X_pan, Y_pan, X_col, Y_col, R_pan, n_pan_tot = RBC_mesh(d,h,b,r, n_pan_fc, X_c, Y_c, theta_rot)
def RBC_mesh(d,h,b,r, n_pan_fc, X_c, Y_c, theta_rot):

    n_pan_r = n_pan_fc/2
    n_pan_d = n_pan_fc/2
    
    step_r = r/(n_pan_r - 1)
    
    X_r = np.zeros(n_pan_r, dtype = np.float64)
    X_r[0] = r
    
    for i in range(1,n_pan_r):
        X_r[i] = X_r[i-1] - step_r 
    X_r[-1] = 0.0
       
    #print X_r, len(X_r), n_pan_r, n_pan_d, step_r
    #print
    A1 = 0.5*b
    A2 = 0.5*(h-b)/r**2
    A3 = (h-b)/r**3
    
    
    Y_r = np.zeros(n_pan_r, dtype = np.float64)
    
    Y_r[:] = A1 + A2*(X_r[:])**2 - A3*(X_r[:]**2)*(X_r[:]-r)
    
    del A1, A2, A3
    
    #print X_r
    #print
    #print Y_r
    
    ####SI X > r
    
    c = 0.5*d - r
    
    XX_i = -0.25*h*np.sqrt(2)
    XX_f = 0.5*c*np.sqrt(2)
    
    #print XX_i
    #print XX_f
    #print
    
    step_XX = (XX_f - XX_i)/(n_pan_d)
    
    XX = np.zeros(n_pan_d, dtype = np.float64)
    XX[0] = XX_f
    for i in range(1,n_pan_r):
        XX[i] = XX[i-1] - step_XX 
       
    
    #XX = np.arange(XX_f, XX_i, -step_XX)
    
    #print
    #print XX, len(XX)
    
    AA1 = -XX_i
    AA2 = -4.*np.sqrt(2)*h/((2.*c+h)**2)
    AA3 = -16.*(2*c-h)/((2.*c+h)**3)
    BB1 = XX_f
    
    YY = np.zeros(n_pan_d, dtype = np.float64)
    
    YY[:] = (AA1 + (XX[:] + AA1)
                 + AA2*(XX[:]+AA1)**2 
                 + AA3*(XX[:]-BB1)*(XX[:]+AA1)**2)
                 
    del AA1, AA2, AA3, BB1, XX_i, XX_f
                 
    X_d, Y_d = np.zeros(n_pan_d, dtype=np.float64), np.zeros(n_pan_d, dtype=np.float64)
    
    X_d[:] = 0.5*(XX[:] + YY[:])*np.sqrt(2) + r
    Y_d[:] = 0.5*(-XX[:] + YY[:])*np.sqrt(2)
    
    #print X_d
    
    del XX, YY
    
    X_pc = np.zeros(n_pan_fc, dtype=np.float64)
    Y_pc = np.zeros(n_pan_fc, dtype=np.float64)
    
    X_pc[:n_pan_d] = X_d[:]
    Y_pc[:n_pan_d] = Y_d[:]
    
    X_pc[n_pan_d:] = X_r[:]
    Y_pc[n_pan_d:] = Y_r[:]
    
    del X_d, Y_d, X_r, Y_r
    
    X_pan = np.zeros(4*n_pan_fc - 2, dtype = np.float64)
    Y_pan = np.zeros(4*n_pan_fc - 2, dtype = np.float64)
    
    X_pan[:n_pan_fc] = X_pc[:]
    Y_pan[:n_pan_fc] = Y_pc[:]
    
    del X_pc, Y_pc
    
    for i in range(n_pan_fc, 2*n_pan_fc-1):
        #print i, 2*n_pan -i -2
        X_pan[i] = -X_pan[2*n_pan_fc -i -2]
        Y_pan[i] = Y_pan[2*n_pan_fc -i -2]
        
    for i in range(2*n_pan_fc - 1, 4*n_pan_fc - 3):
        #print i, 4*n_pan - i - 4
        X_pan[i] = X_pan[4*n_pan_fc -i - 4]
        Y_pan[i] = -Y_pan[4*n_pan_fc -i - 4]
        
    n_pan_tot = 4*n_pan_fc - 3
    
    X_pan[n_pan_tot] = X_pan[0]
    Y_pan[n_pan_tot] = Y_pan[0]
    
    X_aux = np.zeros(n_pan_tot+1, dtype = np.float64)
    Y_aux = np.zeros(n_pan_tot+1, dtype = np.float64)
    
    #X_aux[:] = X_pan[:]*cos(theta_rot)
    #Y_aux[:] = Y_pan[:]*sin(theta_rot)
    
    X_aux[:] = X_pan[:]*cos(theta_rot) - Y_pan[:]*sin(theta_rot)
    Y_aux[:] = X_pan[:]*sin(theta_rot) + Y_pan[:]*cos(theta_rot)
    
    X_pan[:] = X_aux[:] 
    Y_pan[:] = Y_aux[:]
    
    del X_aux, Y_aux
    
    X_pan[:] += X_c
    Y_pan[:] += Y_c
    
    R_pan = np.zeros(n_pan_tot, dtype = np.float64)
    X_col = np.zeros(n_pan_tot, dtype = np.float64)
    Y_col = np.zeros(n_pan_tot, dtype = np.float64)

    R_pan[:] = ((X_pan[1:] - X_pan[:n_pan_tot])**2 + 
                (Y_pan[1:] - Y_pan[:n_pan_tot])**2)**0.5
                
    X_col[:] = (X_pan[1:]+X_pan[:n_pan_tot])*0.5
    Y_col[:] = (Y_pan[1:]+Y_pan[:n_pan_tot])*0.5
    
    fig2 = pyplot.figure(num = None, figsize = (8, 6), dpi = 80, 
                     facecolor = 'w', edgecolor = 'k')
    #pyplot.plot(X_r,Y_r, 'go--', linewidth = 3)
    pyplot.plot(X_pan,Y_pan, 'ro-', linewidth = 3)
    pyplot.plot(X_col,Y_col, 'bo', linewidth = 3)
    #pyplot.plot(X_c,Y_c, 'go', linewidth = 3)
    #pyplot.grid(True)
    pyplot.axis([-1.3*d/2 + X_c,1.3*d/2 + X_c,-1.3*d/2 + Y_c,1.3*d/2 + Y_c])
    #    pyplot.grid(color = '0.3', linestyle = '--', linewidth = 0.3)
    pyplot.savefig('cell.png')
    pyplot.show()

    return X_pan, Y_pan, X_col, Y_col, R_pan, n_pan_tot
    

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

