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
#n_pan = 80        # numero de paneles discretizados         
u_0 = 1.0      # velocidad [m/s]

#R_cil = [0.5, 0.05, 0.1, 0.25]
#u_0 = [1., 10., 5., 2.]

L = 2.*R_cil        # Diametro cilindro [m] (longitud caracteristica)
#L = 3.0

#n_pan = 160
n_pan = [40, 64, 80, 128, 160] 


nq = 6         # nodos cuadratura gaussiana
mu = 1.0          
rho = 1.0

ni = mu/rho

Re = u_0*(L)/ni

#R_cil = 5.E-7         # radio e-6 [m]
#u_0 = 1.E-6            # velocidad e-6 [m/s] 
print
print "###########################################################"
print "############PRUEBA: ECUACION DIMENSIONAL##################"
print "###########################################################"
print
print 'Re = ', Re

X_c, Y_c = 0.0 , 0.0
#X_cad, Y_cad = X_c/L, Y_c/L
#u_0ad = u_0/u_0
#R_ad = R_cil/L 
#factor = (mu*u_0)/L
#print 
#print factor
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


def test_DRAG(R_cil, n_pan, nq, DRAGL, u_0, X_c, Y_c, mu, rho):
    h_theta = 2.0*pi/(n_pan) 
    X_col, Y_col, X_pan, Y_pan, R_pan = geometria_cil(R_cil, n_pan, h_theta, X_c, Y_c)    
    t = nodos_gauss_chebyshev(nq)   
    tx_pan, ty_pan, nx_pan, ny_pan = unit_vector1(X_pan, Y_pan, R_pan)
    
    print
    
    #print X_col
    #print Y_col
    
    target = targetvector(n_pan,X_col,Y_col)
    #print target.z
    
    sources = sourcevector(X_pan, Y_pan, t)
    
    
    #print X_col
    #print Y_col
    #print 
    #print nx_pan
    #print ny_pan
    #print
    #print tx_pan
    #print ty_pan
    
    B = Bmatrix(X_pan, Y_pan, X_col, Y_col, R_pan, t, nx_pan, ny_pan)
    A = Amatrix(X_pan, Y_pan, X_col, Y_col, R_pan, t, nx_pan, ny_pan, mu)
    
    #B = Bmatrix2(X_pan, Y_pan, X_col, Y_col, R_pan, w2, t2, nx_pan, ny_pan)
    #A = Amatrix2(X_pan, Y_pan, X_col, Y_col, R_pan, w2, t2, nx_pan, ny_pan, mu)
        
    ux = np.zeros(n_pan, dtype = np.float64)
    uy = np.zeros(n_pan, dtype = np.float64)
        
    ux[:] = u_0
    
    #vector_b_p2p(target,sources, nq, ux, uy, R_pan,nx_pan,ny_pan)
    
    
    #####linea 111 a 128 es debugueo, descomentar si persiste el bug###
    
    #f = np.zeros(2*n_pan, dtype = np.float64)
    #r = np.zeros(2*n_pan, dtype = np.float64)
    #f[:] = 2.0
    
    #r = np.matmul(A,f)
    
    #print r
    
    #print
    #print
    #print A
    #print
    
    #for i in range(2*n_pan):
    #    print A[i,i]
    
    #print 
    
    #vector_Ax_p2p(target,sources,nq,f,R_pan,X_pan,Y_pan,mu)
    
    
    #####linea 111 a 128 es debugueo, descomentar si persiste el bug###
    
    
    #print Y_col
    #print ux
        
    #fx, fy = f_traction(A,B, ux, uy)
    #fx, fy = f_traction_freeM(A,target, sources, nq, ux, uy,R_pan,nx_pan,ny_pan)
    
    #fx,fy = f_traction_gmres_Mfree(A,target,
    #        sources, nq, ux, uy, R_pan, nx_pan, ny_pan)
            
    fx, fy = f_traction_gmres_Mfree(target, sources, 
                nq,X_pan, Y_pan, ux, uy, R_pan, nx_pan, ny_pan, mu)
    #fx, fy = f_traction_gmres(A,B, ux, uy)
    
    
    #NDRAG = tangential_of_ft(fx,fy, tx_pan, ty_pan, L*R_pan)
    NDRAG = calc_ndrag(fx, R_pan)
    print NDRAG 

    #print fx
    #print 
    #print fy
    
    #norm = np.zeros(len(fx), dtype = np.float64)
    
    #norm[:] = sqrt(fx[:]**2 + fy[:]**2)
    
    #    dR = R_pan[0]
    
    
    
    #    DRAG = 4.*pi*mu*u_0/(0.5-eul_const-log(u_0*R_cil/(4.*ni)))
    #    DRAGE = 6.*pi*mu*R_cil*u_0
    
    print 'Numero de paneles: ', n_pan
    #print
    print 'Calculo del Drag formula cilindro = ', DRAGL
    #print DRAGL
    #print 
    #print 'Calculo del Drag formula esfera = ', DRAGE
    #    print
    print 'Calculo Drag (fx * dR)= ', NDRAG
    #print 'Calculo Cd = ', NDRAG/(rho*u_0**2)
    
    #    print
    #   print 'error absoluto = ', abs(NDRAG - DRAGL)
    #    print 
    print 'error porcentual relativo = ', 100.*abs(NDRAG - DRAGL)/abs(DRAGL), '%'
    #print 'Calculo Drag (fx * dR) = ', sum(fx)*dR
    #print sum(fy)*dR
    #print B
    
    return abs(NDRAG - DRAGL), 100.*abs(NDRAG - DRAGL)/abs(DRAGL), NDRAG

#print A

t_calc = np.zeros(len(n_pan), dtype = np.float64)
e_abs = np.zeros(len(n_pan), dtype = np.float64)
e_rel = np.zeros(len(n_pan), dtype = np.float64)
ndr = np.zeros(len(n_pan), dtype = np.float64)
#cd = np.zeros(len(n_pan), dtype = np.float64)
#ecd = np.zeros(len(n_pan), dtype = np.float64)

#print 'numero de paneles = ', n_pan
#ti = time()
#e1,e2, nd = test_DRAG(R_cil, n_pan, DRAGL, u_0ad, X_cad, Y_cad, L, factor)
#tf = time()
#t_cal = tf - ti
#print 'timepo de calculo [s]= ', t_cal

for i in range(len(n_pan)):
    print
    print 'numero de paneles = ', n_pan[i]
#    ti = time()
    e_abs[i], e_rel[i], ndr[i] = test_DRAG(R_cil, n_pan[i], nq,
                                    DRAGL, u_0, X_c, Y_c, mu, rho)
#    tf = time()
#    t_calc[i] = tf - ti
#    print 'timepo de calculo [s]= ', t_calc[i]
#    print 

#for i in range(4):
#    e_abs[i], e_rel[i], ndr[i] = test_DRAG(R_cil[i], n_pan[0], nq,
#                                    DRAGL, u_0[i], X_c, Y_c, mu, rho)
    

   
#cd[:] = ndr[:]/(rho*L*(u_0**2))
#ecd[:] = 100*abs(Cd_teo - cd[:])/abs(Cd_teo)    
    
#ap_conv = log((ndr[0]-ndr[1])/(ndr[1]-ndr[2]))/log(2.)

#print
#print ndr

#print 
#print 'La cantidad de nodos de cuadratura: ', nq
#print 'convergencia aparente= ',ap_conv


#print     
#print Cd_teo    
#print cd
#print ecd[:]


#fig2 = pyplot.figure(num = None, figsize = (8, 6), dpi = 80, 
#                 facecolor = 'w', edgecolor = 'k')
#pyplot.loglog(n_pan,t_calc, 'go--', linewidth = 3)
#pyplot.title("Tiempo de calculo Vs numero de paneles", fontsize = 14)
#pyplot.xlabel("Numero de paneles")
#pyplot.ylabel("Tiempo de calculo [s]")
#pyplot.grid(True)
#pyplot.grid(color = '0.3', linestyle = '--', linewidth = 0.3)
#pyplot.savefig('time.png')
#pyplot.show()

#fig3 = pyplot.figure(num = None, figsize = (8, 6), dpi = 80, 
#                 facecolor = 'w', edgecolor = 'k')
#pyplot.plot(n_pan,e_abs, 'go--', linewidth = 3)
#pyplot.title("error absoluto Vs numero de paneles", fontsize = 14)
#pyplot.xlabel("Numero de paneles")
#pyplot.ylabel("error")
#pyplot.grid(True)
#pyplot.grid(color = '0.3', linestyle = '--', linewidth = 0.3)
#pyplot.savefig('error_abs.png')
#    pyplot.show()

#fig3 = pyplot.figure(num = None, figsize = (8, 6), dpi = 80, 
#                 facecolor = 'w', edgecolor = 'k')
#pyplot.loglog(n_pan,e_rel, 'go--', linewidth = 3)
#pyplot.title("error relativo Vs numero de paneles", fontsize = 14)
#pyplot.xlabel("Numero de paneles")
#pyplot.ylabel("error relativo %")
#pyplot.grid(True)
#pyplot.grid(color = '0.3', linestyle = '--', linewidth = 0.3)
#pyplot.savefig('error_rel.png')
#pyplot.show()