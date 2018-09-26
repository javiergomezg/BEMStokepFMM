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

""""
Parametros:

d: Diametro principal del globulo rojo
h: Altura mayor
b: altura menor
r: radio interno
"""""""""

d = 7.82E-6
h = 2.58E-6
b = 0.81E-6
r = d/3

#d = 10.
#h = 4.
#b = 2.
#r = d/3

n_pan_fc = [40]#, 80, 160]

#centro del eritrocito

X_c = 0.0
Y_c = 0.0

theta_rot = -10.0
theta_rot = theta_rot*pi/180

u_0 = 1.0

#n_pan = [100, 200, 400, 1000, 2000, 4000, 10000]#, 20000, 40000] 

nq = 8         # nodos cuadratura gaussiana
mu = 1.0          
rho = 1.0

ni = mu/rho

Re = u_0*(d)/ni

nc_pol = 8 #cantidad de nodos que se utilizan para aproximar el kernel

#R_cil = 5.E-7         # radio e-6 [m]
#u_0 = 1.E-6            # velocidad e-6 [m/s] 
print
print "###########################################################"
print "############PRUEBA: ECUACION DIMENSIONAL##################"
print "###########################################################"
print
print 'Re = ', Re


################################################################
########## Formulas Drag de Lamb################################
################################################################

def test_DRAG_RBC(d,h,b,r, theta_rot, n_pan_fc, nq, nc_pol, u_0, X_c, Y_c, mu, rho):
#def test_DRAG(R_cil, n_pan, nq, nc_pol, DRAGL, u_0, X_c, Y_c, mu, rho):
    
    #h_theta = 2.0*pi/(n_pan)
    X_pan, Y_pan, X_col, Y_col, R_pan, n_pan = RBC_mesh(d,h,b,r, 
                                                        n_pan_fc, 
                                                        X_c, Y_c, 
                                                        theta_rot)
    #X_col, Y_col, X_pan, Y_pan, R_pan = geometria_cil(R_cil, n_pan, h_theta, X_c, Y_c)    
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
    torque_direct = TORQUE(fx, fy, X_col, Y_col, X_c, Y_c, R_pan)
    dirT_f = time()
     
##################################################################################    
    
##################################################################################    
    ## Calculo del DRAG usando FMM
    
    spupT_i = time()
    
    fx_fmm, fy_fmm = f_traction_gmres_Mfree_fmm(target, sources, 
                nq,X_pan, Y_pan, ux, uy, R_pan, nx_pan, ny_pan, mu, nc_pol)
    torque_fmm = TORQUE(fx_fmm, fy_fmm, X_col, Y_col, X_c, Y_c, R_pan)
    NDRAG_fmm = calc_ndrag(fx_fmm, R_pan)
    spupT_f = time()
##################################################################################  
    
    dirT = dirT_f - dirT_i
    
    fmmT = spupT_f - spupT_i
    
    err = abs(NDRAG_fmm - NDRAG_direct)
    
    print
    #print "NDRAG_direct = ", NDRAG_direct
    #print "tiempo calculo directo = ", dirT_f - dirT_i, "[s]"     
    #print    
    print "NDRAG_fmm = ", NDRAG_fmm
    print "tiempo calculo con bbfmm = ", spupT_f - spupT_i, "[s]"
    #print
    
    #print "error = ", abs(NDRAG_fmm - NDRAG_direct)
    print
   
    print 'Numero de paneles: ', n_pan
    print 'numero de particulas: ', n_pan*nq + n_pan
    
    return NDRAG_direct, torque_direct, NDRAG_fmm, torque_fmm, dirT, fmmT, err, n_pan
    #return NDRAG_fmm, torque_fmm, fmmT, n_pan



#test_DRAG(R_cil, n_pan, nq, nc_pol, DRAGL, u_0, X_c, Y_c, mu, rho)

NDRAG_dir = np.zeros(len(n_pan_fc), dtype = np.float64)
tq_dir = np.zeros(len(n_pan_fc), dtype = np.float64)
NDRAG_fmm = np.zeros(len(n_pan_fc), dtype = np.float64)
tq_fmm = np.zeros(len(n_pan_fc), dtype = np.float64)
T_dir = np.zeros(len(n_pan_fc), dtype = np.float64)
T_fmm = np.zeros(len(n_pan_fc), dtype = np.float64)
error = np.zeros(len(n_pan_fc), dtype = np.float64)
n_pan = np.zeros(len(n_pan_fc), dtype = np.int64)

for i in range(len(n_pan_fc)):
    NDRAG_dir[i], tq_dir[i], NDRAG_fmm[i], tq_fmm[i], T_dir[i], T_fmm[i], error[i], n_pan[i] = test_DRAG_RBC(d,h,b,r, 
                                                                          theta_rot, 
                                                                          n_pan_fc[i], 
                                                                          nq, nc_pol,
                                                                          u_0, 
                                                                          X_c, Y_c, 
                                                                          mu, rho)                                                

#    NDRAG_fmm[i], tq_fmm[i], T_fmm[i], n_pan[i] = test_DRAG_RBC(d,h,b,r, 
#                                                                theta_rot, 
#                                                                n_pan_fc[i], 
#                                                                 nq, nc_pol,
#                                                                 u_0, 
#                                                                 X_c, Y_c, 
#                                                                 mu, rho)                                                


for i in range(len(n_pan_fc)):
    print
    print "Numero de paneles: ", n_pan[i]
    print "Particulas totales: ", n_pan[i]*nq + n_pan[i]
    print "DRAG calculo directo: ", NDRAG_dir[i], "[N/m]"
    print "Torque calculo directo: ", tq_dir[i], "[N]"
    print "DRAG calculo BBFMM: ", NDRAG_fmm[i], "[N]"
    print "Torque calculo BBFMM: ", tq_fmm[i], "[N]"
    print "error: ", error[i]
    print "tiempo calculo directo: ", T_dir[i], "[s]"
    print "tiempo calculo con BBFMM: ", T_fmm[i], "[s]"

#print (NDRAG_dir[0]*NDRAG_dir[2]-NDRAG_dir[1]**2)/(NDRAG_dir[0]-2*NDRAG_dir[1]+NDRAG_dir[2])