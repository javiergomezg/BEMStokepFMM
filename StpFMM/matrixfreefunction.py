import numpy as np
import Classnode
import cuadratura
import fs_gf
from numpy import pi


def targetvector(Ncol,Xcol,Ycol):
    target = Target(Ncol)
    target.z[:].real = Xcol[:]    
    target.z[:].imag = Ycol[:]
    return target
    
def sourcevector(Xpan, Ypan, t):
    N_sources = (len(Xpan)-1)*len(t)
    sources = Sources(N_sources)
    #print len(Xpan)
    #print len(t)    
    
    #print sources.z
    #print len(sources.z)
    
    i = 0
    for p in range(len(Xpan)-1):
        Xg, Yg = nodos_gauss_panel(t, Xpan[p], Ypan[p], Xpan[p+1], Ypan[p+1])
        #print Xpan[p], Ypan[p]
        #print Xg
        #print Yg
        #print Xpan[p+1], Ypan[p+1]
        
        for j in range(len(t)):
            sources.z[i] = Xg[len(t)-1-j] + 1j*Yg[len(t)-1-j]
            sources.w[i] = (pi/len(t))*(1.0-(t[len(t)-1-j])**2)**(0.5)
            #sources.z[i] = Xg[j] + 1j*Yg[j]
            #sources.w[i] = (1.0-t[j]**2)**(0.5)
            i+=1
            
    #print sources.z
    #print
    #print t
    #print
    #print sources.w
    
    return sources
    
def vector_b_p2p(T,S,nc, ux, uy, Rp,nx, ny):
    Nq = len(T.z)
    Ns = len(S.z)
    
    #print Nq, Ns
    
    bx,by = np.zeros(Nq, dtype=np.float64), np.zeros(Nq, dtype=np.float64)
    bx[:],by[:] = ux[:], uy[:]
    
    
    #print bx, ux
    
    X = np.zeros(Ns)
    Y = np.zeros(Ns)

    X[:] = S.z[:].real
    Y[:] = S.z[:].imag    
    
    #print X
    #print Y    
    
    for q in range(Nq):
        X0 = T.z[q].real
        Y0 = T.z[q].imag
        Txxx, Txxy, Tyxx, Tyxy, Txyx, Txyy, Tyyx, Tyyy =gf_stress_fs(X,Y,X0,Y0)
        for ic in range(Ns):
            if (ic//nc != q):
                bx[q] = bx[q] + (-1./(2.*pi))*(Rp[ic//nc]/2.)*S.w[ic]*(ux[ic//nc]*
                (Txxx[ic]*nx[ic//nc]+Txxy[ic]*ny[ic//nc]) + uy[ic//nc]*
                (Tyxx[ic]*nx[ic//nc]+Tyxy[ic]*ny[ic//nc]))
                
                by[q] = by[q] + (-1./(2.*pi))*(Rp[ic//nc]/2.)*S.w[ic]*(ux[ic//nc]*
                (Txyx[ic]*nx[ic//nc]+Txyy[ic]*ny[ic//nc]) + uy[ic//nc]*
                (Tyyx[ic]*nx[ic//nc]+Tyyy[ic]*ny[ic//nc]))
                
    b = np.zeros([2*Nq,1], dtype = np.float64)
    
    b[:Nq,0], b[Nq:,0] = bx[:],by[:]
    
    #print b
    #print
    
    return b
    
def vector_Ax_p2p(T,S,nc,f,Rp,Xp,Yp,mu):
    Nq, Ns = len(T.z), len(S.z)

    fx, fy = np.zeros(Nq, dtype = np.float64), np.zeros(Nq, dtype = np.float64)
    fx[:], fy[:] = f[:Nq], f[Nq:]
    
    Ax = np.zeros(2*Nq, dtype = np.float64)
    
        
    for q in range(Nq):
        Ax[q] = (1./(2*pi*mu))*(  # (Rp[q]/2.)*(
                int_alpha_xx(Xp[q],Xp[q+1],Yp[q],Yp[q+1])*fx[q] + 
                int_alpha_xy(Xp[q],Xp[q+1],Yp[q],Yp[q+1])*fy[q])
        Ax[Nq + q] = (1./(2*pi*mu))*(  # (Rp[q]/2.)*(
                int_alpha_xy(Xp[q],Xp[q+1],Yp[q],Yp[q+1])*fx[q] + 
                int_alpha_yy(Xp[q],Xp[q+1],Yp[q],Yp[q+1])*fy[q])
                
    #print Ax
    
    X = np.zeros(Ns)
    Y = np.zeros(Ns)

    X[:] = S.z[:].real
    Y[:] = S.z[:].imag    
    
    for q in range(Nq):
        X0 = T.z[q].real
        Y0 = T.z[q].imag        
        Gxx, Gxy, Gyx, Gyy = gf_vel_fs(X,Y,X0,Y0)
        for ic in range(Ns):
            if (ic//nc != q):
                Ax[q] = Ax[q] + (1./(2*pi*mu))*(Rp[ic//nc]/2.)*S.w[ic]*(
                        Gxx[ic]*fx[ic//nc] + Gxy[ic]*fy[ic//nc])
                        
                Ax[q + Nq] = Ax[q + Nq] + (1./(2.*pi*mu))*(
                             Rp[ic//nc]/2.)*S.w[ic]*(
                             Gyx[ic]*fx[ic//nc] + Gyy[ic]*fy[ic//nc])
    
    return Ax
    
    
    
                
    
    