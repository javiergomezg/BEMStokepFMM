import cilmesh
import cuadratura
import matrix_gen
import matrixfreefunction
import Classnode
import solver
import numpy as np
import time
from time import time
from numpy import pi
from numpy import log, sqrt

R_cil = 0.5   # radio del cilindo [m]
u_0 = 1.0      # velocidad [m/s]

L = 2.*R_cil        # Diametro cilindro [m] (longitud caracteristica)

n_pan = [100, 200, 400, 1000, 2000, 4000, 10000]#, 20000, 40000] 

nq = 8         # nodos cuadratura gaussiana
mu = 1.0          
rho = 1.0

ni = mu/rho

Re = u_0*(L)/ni

nc_pol = 8 #cantidad de nodos que se utilizan para aproximar el kernel

#R_cil = 5.E-7         # radio e-6 [m]
#u_0 = 1.E-6            # velocidad e-6 [m/s] 
print
print "###########################################################"
print "############PRUEBA: ECUACION DIMENSIONAL##################"
print "###########################################################"
print
print 'Re = ', Re

X_c, Y_c = 0.0 , 0.0 # Centro del cilindro

print

################################################################
########## Formulas Drag de Lamb################################
################################################################

eul_const = 0.577215664901532860606
eps = (0.5-eul_const - log(Re/4.))**(-1)
DRAGL = 4.*pi*mu*u_0*eps

Del1 = (log(4./Re)-eul_const-0.5)**(-1)
Cd_teo = (4.*pi/Re)*(Del1 -0.87*(Del1**3))

################################################################
########## Formulas Drag de Lamb################################
################################################################


def test_DRAG(R_cil, n_pan, nq, nc_pol, DRAGL, u_0, X_c, Y_c, mu, rho):
    
    h_theta = 2.0*pi/(n_pan) 
    X_col, Y_col, X_pan, Y_pan, R_pan = geometria_cil(R_cil, n_pan, h_theta, X_c, Y_c)    
    t = nodos_gauss_chebyshev(nq)   
    tx_pan, ty_pan, nx_pan, ny_pan = unit_vector1(X_pan, Y_pan, R_pan)
    
    target = targetvector(n_pan,X_col,Y_col)
    
    sources = sourcevector(X_pan, Y_pan, t)
    
    ux = np.zeros(n_pan, dtype = np.float64)
    uy = np.zeros(n_pan, dtype = np.float64)
        
    ux[:] = u_0
    
##################################################################################    
    ## Calculo directo del DRAG (Sin speed up)
            
    dirT_i = time()
    fx, fy = f_traction_gmres_Mfree(target, sources, 
                nq,X_pan, Y_pan, ux, uy, R_pan, nx_pan, ny_pan, mu)
                
    NDRAG_direct = calc_ndrag(fx, R_pan)
    dirT_f = time()
     
##################################################################################    
    
##################################################################################    
    ## Calculo del DRAG usando FMM
    
    spupT_i = time()
    
    fx_fmm, fy_fmm = f_traction_gmres_Mfree_fmm(target, sources, 
                nq,X_pan, Y_pan, ux, uy, R_pan, nx_pan, ny_pan, mu, nc_pol)
    
    NDRAG_fmm = calc_ndrag(fx_fmm, R_pan)
    spupT_f = time()
##################################################################################  
    
    dirT = dirT_f - dirT_i
    
    fmmT = spupT_f - spupT_i
    
    err = abs(NDRAG_fmm - NDRAG_direct)
    
    print
    print "NDRAG_direct = ", NDRAG_direct
    print "tiempo calculo directo = ", dirT_f - dirT_i, "[s]"     
    print    
    print "NDRAG_fmm = ", NDRAG_fmm
    print "tiempo calculo con bbfmm = ", spupT_f - spupT_i, "[s]"
    print
    
    print "error = ", abs(NDRAG_fmm - NDRAG_direct)
    print
   
    print 'Numero de paneles: ', n_pan
    print 'numero de particulas: ', n_pan*nq + n_pan
    
    return NDRAG_direct, NDRAG_fmm, dirT, fmmT, err



#test_DRAG(R_cil, n_pan, nq, nc_pol, DRAGL, u_0, X_c, Y_c, mu, rho)

NDRAG_dir = np.zeros(len(n_pan), dtype = np.float64)
NDRAG_fmm = np.zeros(len(n_pan), dtype = np.float64)
T_dir = np.zeros(len(n_pan), dtype = np.float64)
T_fmm = np.zeros(len(n_pan), dtype = np.float64)
error = np.zeros(len(n_pan), dtype = np.float64)

for i in range(len(n_pan)):
    NDRAG_dir[i], NDRAG_fmm[i], T_dir[i], T_fmm[i], error[i] = test_DRAG(R_cil, 
                                                                n_pan[i], 
                                                                nq, nc_pol, 
                                                                DRAGL, u_0, 
                                                                X_c, Y_c, 
                                                                mu, rho)

for i in range(len(n_pan)):
    print
    print "Numero de paneles: ", n_pan[i]
    print "Particulas totales: ", n_pan[i]*nq + n_pan[i]
    print "DRAG calculo directo: ", NDRAG_dir[i], "[N/m]"
    print "DRAG calculo BBFMM: ", NDRAG_fmm[i], "[N/m]"
    print "error: ", error[i]
    print "tiempo calculo directo: ", T_dir[i], "[s]"
    print "tiempo calculo con BBFMM: ", T_fmm[i], "[s]"

#print (NDRAG_dir[0]*NDRAG_dir[2]-NDRAG_dir[1]**2)/(NDRAG_dir[0]-2*NDRAG_dir[1]+NDRAG_dir[2])