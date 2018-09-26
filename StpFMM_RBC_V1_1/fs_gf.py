import numpy as np
from numpy import log, pi, sqrt

# Funcion de Green para la velocidad
def gf_vel_fs(X,Y,X0,Y0):
    dx, dy = np.zeros(len(X), dtype = np.float64), np.zeros(len(X), dtype = np.float64)
    r = np.zeros(len(X), dtype = np.float64)
    
    dx[:], dy[:] = X0 - X[:], Y0 - Y[:]
    r[:] = sqrt(dx[:]**2 + dy[:]**2)
    
    Gxx = np.zeros(len(X), dtype = np.float64)
    Gxy = np.zeros(len(X), dtype = np.float64)
    Gyx = np.zeros(len(X), dtype = np.float64)
    Gyy = np.zeros(len(X), dtype = np.float64)
    
    Gxx[:] = -log(r[:]) + (dx[:]**2)/(r[:]**2)
    Gxy[:] = dx[:]*dy[:]/(r[:]**2)
    Gyx[:] = Gxy[:]
    Gyy[:] = -log(r[:]) + (dy[:]**2)/(r[:]**2)
    
    del dx, dy, r    
    
    return Gxx, Gxy, Gyx, Gyy
    

# Funcion de Green para la presion
def gf_pres_fs(X,Y,X0,Y0):
    dx, dy = np.zeros(len(X), dtype = np.float64), np.zeros(len(X), dtype = np.float64)
    r = np.zeros(len(X), dtype = np.float64)
    
    dx[:], dy[:] = X[:] - X0, Y[:] - Y0
    r[:] = sqrt(dx[:]**2 + dy[:]**2)
    
    px = np.zeros(len(X), dtype = np.float64)
    py = np.zeros(len(X), dtype = np.float64)
    
    px[:] = 2*dx[:]/(r[:]**2)
    py[:] = 2*dy[:]/(r[:]**2)
    
    del dx, dy, r
    
    return px, py
    
def gf_stress_fs(X,Y,X0,Y0):
    dx, dy = np.zeros(len(X), dtype = np.float64), np.zeros(len(X), dtype = np.float64)
    r = np.zeros(len(X), dtype = np.float64)
    
    dx[:], dy[:] = X[:] - X0, Y[:] - Y0
    r[:] = sqrt(dx[:]**2 + dy[:]**2)
    
    Txxx = np.zeros(len(X), dtype = np.float64)
    Txxy = np.zeros(len(X), dtype = np.float64)
    
    Tyxx = np.zeros(len(X), dtype = np.float64)
    Tyxy = np.zeros(len(X), dtype = np.float64)
    
    Txyx = np.zeros(len(X), dtype = np.float64)
    Txyy = np.zeros(len(X), dtype = np.float64)
    
    Tyyx = np.zeros(len(X), dtype = np.float64)
    Tyyy = np.zeros(len(X), dtype = np.float64)
    
    Txxx[:] = -4.*(dx[:]**3)/(r[:]**4)
    Txxy[:] = -4.*(dy[:]*dx[:]**2)/(r[:]**4)
    
    Tyxx[:] = Txxy[:]
    Tyxy[:] = -4.*(dx[:]*dy[:]**2)/(r[:]**4)
    
    Txyx[:] = Txxy[:]
    Txyy[:] = Tyxy[:]
    
    Tyyx[:] = Txyy[:]
    Tyyy[:] = -4*(dy[:]**3)/(r[:]**4)
    
    return Txxx, Txxy, Tyxx, Tyxy, Txyx, Txyy, Tyyx, Tyyy
    
    