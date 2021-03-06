import numpy as np
import gmres
import matrixfreefunction
from numpy import pi, cos, sin, arctan2


def f_traction(A,B, ux, uy):
    n = len(ux)
    u = np.zeros([2*n,1], dtype = np.float64)
    f = np.zeros([2*n,1], dtype = np.float64)
    
    u[:n,0], u[n:,0] = ux[:], uy[:]
    
    b = np.matmul(B,u)
    print b
    f = np.linalg.solve(A, b)
    
    #print f
    
    fx = np.zeros(n, dtype = np.float64) 
    fy = np.zeros(n, dtype = np.float64)
    
    fx[:], fy[:] = f[:n,0], f[n:,0]
    
    #print fx
    #print
    #print fy
    
    return fx, fy
    
    
def f_traction_freeM(A,T, S, nc, ux, uy,Rp,nx,ny):
    n = len(ux)
    #u = np.zeros([2*n,1], dtype = np.float64)
    f = np.zeros([2*n,1], dtype = np.float64)
    
    #u[:n,0], u[n:,0] = ux[:], uy[:]
    
    b = vector_b_p2p(T,S,nc, ux, uy, Rp, nx, ny)    
    
    #b = np.matmul(B,u)
    #print b
    f = np.linalg.solve(A, b)
    
    #print f
    
    fx = np.zeros(n, dtype = np.float64) 
    fy = np.zeros(n, dtype = np.float64)
    
    fx[:], fy[:] = f[:n,0], f[n:,0]
    
    #print fx
    #print
    #print fy
    
    return fx, fy
    
def f_traction_gmres(A,B, ux, uy):
    n = len(ux)
    nn = 2*n
    u = np.zeros([2*n,1], dtype = np.float64)
    f = np.zeros([2*n,1], dtype = np.float64)
    
    u[:n,0], u[n:,0] = ux[:], uy[:]
    f = np.zeros(nn, dtype = np.float64)
    f[:] = 1.0E-16
    
    #b = prueba_cont(T,S,nc, ux, uy, Rp,nx, ny)
    b = np.matmul(B,u)
   
    
    #D = np.zeros([nn,nn], dtype=np.float64)
    #for i in range(nn):
    #    D[i,i] = A[i,i]**(-1)


    #AA = np.matmul(D,A)
    #bb = np.matmul(D,b)
    b2 = np.transpose(b)
    
    #print b
    #print b2
    f,_ = gmres_mgs(A,f,b2, nn, 500, 1.0E-4)
    
    fx = np.zeros(n, dtype = np.float64) 
    fy = np.zeros(n, dtype = np.float64)
    
    fx[:], fy[:] = f[:n], f[n:]
    
    del f
    
    #print fx
    #print
    #print fy
    
    return fx, fy

def f_traction_gmres_Mfree(T, S, nc,Xp, Yp, ux, uy, Rp, nx, ny, mu):
    n = len(ux)
    nn = 2*n
    #u = np.zeros([2*n,1], dtype = np.float64)
    f = np.zeros([2*n,1], dtype = np.float64)
    
    #u[:n,0], u[n:,0] = ux[:], uy[:]
    f = np.zeros(nn, dtype = np.float64)
    f[:] = 1.0E-16
    
    #b = prueba_cont(T,S,nc, ux, uy, Rp,nx, ny)
    #b = np.matmul(B,u)
    b = vector_b_p2p(T,S,nc, ux, uy, Rp, nx, ny)
   
    #D = np.zeros([nn,nn], dtype=np.float64)
    #for i in range(nn):
    #    D[i,i] = A[i,i]**(-1)


    #AA = np.matmul(D,A)
    #bb = np.matmul(D,b)
    b2 = np.transpose(b)

    
    #print b

    #f,_ = gmres_mgs(A,f,b2, nn, 500, 1.0E-14)
    f,_ = gmres_mgs(T,S,nc,f,b2,Rp,Xp,Yp,mu, nn, nn, 1.0E-4)
    
    fx = np.zeros(n, dtype = np.float64) 
    fy = np.zeros(n, dtype = np.float64)
    
    fx[:], fy[:] = f[:n], f[n:]
    
    del f
    
    #print fx
    #print
    #print fy
    
    return fx, fy
    
def f_traction_gmres_Mfree_fmm(T, S, nc,Xp, Yp, ux, uy, Rp, nx, ny, mu, n_cheb):
    n = len(ux)
    nn = 2*n
    #u = np.zeros([2*n,1], dtype = np.float64)
    f = np.zeros([nn,1], dtype = np.float64)
    
    #u[:n,0], u[n:,0] = ux[:], uy[:]
    f = np.zeros(nn, dtype = np.float64)
    f[:] = 1.0E-16
    
    #b = prueba_cont(T,S,nc, ux, uy, Rp,nx, ny)
    #b = np.matmul(B,u)
    b = vector_b_fmm(T,S,nc, ux, uy, Rp, nx, ny, n_cheb)
   
    #D = np.zeros([nn,nn], dtype=np.float64)
    #for i in range(nn):
    #    D[i,i] = A[i,i]**(-1)


    #AA = np.matmul(D,A)
    #bb = np.matmul(D,b)
    b2 = np.transpose(b)

    
    #print b

    #f,_ = gmres_mgs(A,f,b2, nn, 500, 1.0E-14)
    print "Entrando al GMRES"
    print
    f,_ = gmres_mgs_fmm(T,S,nc,f,b2,Rp,Xp,Yp,mu, n_cheb, nn, nn, 1.0E-4)
    
    fx = np.zeros(n, dtype = np.float64) 
    fy = np.zeros(n, dtype = np.float64)
    
    fx[:], fy[:] = f[:n], f[n:]
    
    del f
    
    #print fx
    #print
    #print fy
    
    return fx, fy
    
    
    
def tangential_of_ft(fx,fy, tx, ty, h):
    ftn = np.zeros(len(fx), dtype = np.float64)
    dDrag = np.zeros(len(fx), dtype = np.float64)    
    
    ftn[:] = fx[:]*tx[:] + fy[:]*ty[:]
#    print fx
#    print 
#    print fy
#    print 
#    print ftn    
    
    dDrag[:] = fx[:]*h[:]
    
    return sum(dDrag)
    
def calc_ndrag(fx,h):
    dDrag = np.zeros(len(fx), dtype = np.float64)
    dDrag[:] = fx[:]*h[:]
    
    return sum(dDrag)
    
def TORQUE(fx, fy, X_col, Y_col, X_c, Y_c, Rp):
    theta = np.zeros(len(X_col), dtype = np.float64)
    longToc = np.zeros(len(X_col), dtype = np.float64)
    
    for i in range(len(X_col)):
        theta[i] = arctan2(Y_col[i]-Y_c,X_col[i]-X_c)
        longToc[i] = ((X_col[i]-X_c)**2 + (Y_col[i] - Y_c)**2)**0.5
        
    torqueV = np.zeros(len(X_col), dtype = np.float64)
    
    torqueV[:] = (fx[:]*Rp[:]*sin(theta[:]) + fy[:]*Rp[:]*cos(theta[:]))*longToc[:]
    
    torque = sum(torqueV)

    del theta, torqueV, longToc 
    
    return torque
    