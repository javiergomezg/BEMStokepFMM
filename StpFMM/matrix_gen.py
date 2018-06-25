import numpy as np
from numpy import pi, log, sqrt
import cilmesh
import cuadratura

"""""""""""""""""""""""""""""""""
Amatrix =  components of the single layer potential
Bmatrix =  components of the double layer potential

"""""""""""""""""""""""""""""""""


### Usando Gauss Chebyshev

def Bmatrix(X_pan, Y_pan, X_col, Y_col, R_pan, t, nx_pan, ny_pan):
    n = len(X_col)
    bxx = np.zeros([n,n], dtype = np.float64)
    byx = np.zeros([n,n], dtype = np.float64)
    byy = np.zeros([n,n], dtype = np.float64)
    BB = np.zeros([2*n,2*n], dtype = np.float64)
         
    for q in range(n):
        for p in range(n):
            Xg, Yg = nodos_gauss_panel(t, X_pan[p], Y_pan[p], 
                                           X_pan[p+1], Y_pan[p+1])
            b1,b2,b3,b4 = beta_comp_gs(t,Xg,Yg, X_col[q], Y_col[q], 
                                       nx_pan[p], ny_pan[p], R_pan[p])                               
                                           
            bxx[q,p] = b1
            byx[q,p] = b3
            byy[q,p] = b4
            
    for i in range(n):
        bxx[i,i], byy[i,i], byx[i,i] = -2.*pi, -2.*pi, 0.0

    BB[:n,:n], BB[:n,n:] = bxx[:,:], byx[:,:]
    BB[n:,:n], BB[n:,n:] = byx[:,:], byy[:,:]
    
    BB = (-1./(2.*pi)) * BB
                
    return BB
    
    
def Amatrix(X_pan, Y_pan, X_col, Y_col, R_pan, t, nx_pan, ny_pan, mu):
    n = len(X_col)
    axx = np.zeros([n,n], dtype = np.float64)
    axy = np.zeros([n,n], dtype = np.float64)
    ayy = np.zeros([n,n], dtype = np.float64)
    AA = np.zeros([2*n,2*n], dtype = np.float64)
                                     
    
    for q in range(n):
        for p in range(n):
            Xg, Yg = nodos_gauss_panel(t, X_pan[p], Y_pan[p], 
                                           X_pan[p+1], Y_pan[p+1])
                                           
            a1, a2, a3, a4 = alpha_comp_gs(t,Xg,Yg, X_col[q], 
                                           Y_col[q], R_pan[p])
            
            axx[q,p] = a1
            axy[q,p] = a2
            ayy[q,p] = a4


    for i in range(n):
        axx[i,i] = int_alpha_xx(X_pan[i],X_pan[i+1],Y_pan[i],Y_pan[i+1])
        axy[i,i] = int_alpha_xy(X_pan[i],X_pan[i+1],Y_pan[i],Y_pan[i+1]) 
        ayy[i,i] = int_alpha_yy(X_pan[i],X_pan[i+1],Y_pan[i],Y_pan[i+1]) 
        
    AA[:n,:n], AA[:n,n:] = axx[:,:], axy[:,:]
    AA[n:,:n], AA[n:,n:] = axy[:,:], ayy[:,:]
    
    AA = (1.0/(2.0*pi*mu)) * AA
        
    return AA

    
### Usando Gauss Legendre

def Bmatrix2(X_pan, Y_pan, X_col, Y_col, R_pan, w, t, nx_pan, ny_pan):
    n = len(X_col)
    bxx = np.zeros([n,n], dtype = np.float64)
    byx = np.zeros([n,n], dtype = np.float64)
    byy = np.zeros([n,n], dtype = np.float64)
    BB = np.zeros([2*n,2*n], dtype = np.float64)
         
    for q in range(n):
        for p in range(n):
            Xg, Yg = nodos_gauss_panel(t, X_pan[p], Y_pan[p], 
                                           X_pan[p+1], Y_pan[p+1])
            b1,b2,b3,b4 = beta_comp_gl(w,t,Xg,Yg, X_col[q], Y_col[q], 
                                       nx_pan[p], ny_pan[p], R_pan[p])                               
                                           
            bxx[q,p] = b1
            byx[q,p] = b3
            byy[q,p] = b4
            
    for i in range(n):
        bxx[i,i], byy[i,i], byx[i,i] = -2.*pi, -2.*pi, 0.0

    BB[:n,:n], BB[:n,n:] = bxx[:,:], byx[:,:]
    BB[n:,:n], BB[n:,n:] = byx[:,:], byy[:,:]
    
    BB = (-1./(2.*pi)) * BB
                
    return BB
    
    
def Amatrix2(X_pan, Y_pan, X_col, Y_col, R_pan, w, t, nx_pan, ny_pan, mu):
    n = len(X_col)
    axx = np.zeros([n,n], dtype = np.float64)
    axy = np.zeros([n,n], dtype = np.float64)
    ayy = np.zeros([n,n], dtype = np.float64)
    AA = np.zeros([2*n,2*n], dtype = np.float64)
                                     
    
    for q in range(n):
        for p in range(n):
            Xg, Yg = nodos_gauss_panel(t, X_pan[p], Y_pan[p], 
                                           X_pan[p+1], Y_pan[p+1])
                                           
            a1, a2, a3, a4 = alpha_comp_gl(w, t,Xg,Yg, X_col[q], 
                                           Y_col[q], R_pan[p])
            
            axx[q,p] = a1
            axy[q,p] = a2
            ayy[q,p] = a4


    for i in range(n):
        axx[i,i] = int_alpha_xx(X_pan[i],X_pan[i+1],Y_pan[i],Y_pan[i+1])
        axy[i,i] = int_alpha_xy(X_pan[i],X_pan[i+1],Y_pan[i],Y_pan[i+1]) 
        ayy[i,i] = int_alpha_yy(X_pan[i],X_pan[i+1],Y_pan[i],Y_pan[i+1]) 
        
    AA[:n,:n], AA[:n,n:] = axx[:,:], axy[:,:]
    AA[n:,:n], AA[n:,n:] = axy[:,:], ayy[:,:]
    
    AA = (1.0/(2.0*pi*mu)) * AA
        
    return AA

    